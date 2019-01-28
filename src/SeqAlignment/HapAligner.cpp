#include <algorithm>
#include <climits>
#include <string>

#include "AlignmentState.h"
#include "AlignmentTraceback.h"
#include "HapAligner.h"

// Minimum distance of a seed base from an indel, mismatch or a repetitive region
const int32_t MIN_SEED_DIST = 5;

/*
   Computes the position in the provided region that is the furthest from region boundaries and repetitive haplotype blocks.
   Stores the resulting position and distance using the references. If no such position exists, both values are set to -1
 */
void HapAligner::calc_best_seed_position(int32_t region_start, int32_t region_end,
					 int32_t& best_dist, int32_t& best_pos){
  best_dist = best_pos = -1;
  int32_t pos = region_start;
  int repeat_index = 0;
  while (repeat_index < repeat_starts_.size() && pos <= region_end){
    if (pos < repeat_starts_[repeat_index]){
      int32_t dist = 1 + (std::min(region_end, repeat_starts_[repeat_index]-1)-pos)/2;
      if (dist >= best_dist){
	best_dist = dist;
	best_pos  = dist - 1 + pos;
      }
      pos = repeat_ends_[repeat_index++];
    }
    else if (pos < repeat_ends_[repeat_index])
      pos = repeat_ends_[repeat_index++];
    else
      repeat_index++;
  }
  if (pos <= region_end){
    int32_t dist = 1 + (region_end-pos)/2;
    if (dist >= best_dist){
      best_dist = dist;
      best_pos  = dist - 1 + pos;
    }
  }
}

/* 
 * Identify the base with the largest minimum distance from an insertion, a deletion and a stutter block
 * as defined by its alignment to the reference genome
 */
int HapAligner::calc_seed_base(const Alignment& aln){
  int32_t pos          = aln.get_start();
  int best_seed = -1, cur_base = 0, max_dist = MIN_SEED_DIST;
  for (auto cigar_iter = aln.get_cigar_list().begin(); cigar_iter != aln.get_cigar_list().end(); cigar_iter++){
    switch(cigar_iter->get_type()){
    case '=': {
      int32_t min_region = pos, max_region = pos + cigar_iter->get_num() - 1; 

      // Clip region so that the seed doesn't extend beyond the maximal sequence flanking the repeats
      // Only this region will be used to construct haplotypes
      min_region = std::max(min_region, fw_haplotype_->get_first_block()->start());
      max_region = std::min(max_region, fw_haplotype_->get_last_block()->end()-1);

      if (min_region <= max_region){
	int32_t distance, dist_pos;
	calc_best_seed_position(min_region, max_region, distance, dist_pos);
	if (distance >= max_dist){
	  max_dist  = distance;
	  best_seed = cur_base + (dist_pos - pos);
	}
      }
      pos      += cigar_iter->get_num();
      cur_base += cigar_iter->get_num();
      break;
    }
    case 'I': {
      cur_base += cigar_iter->get_num();
      break;
    }
    case 'X': {
      pos      += cigar_iter->get_num();
      cur_base += cigar_iter->get_num();
      break;
    }
    case 'D': {
      pos += cigar_iter->get_num();
      break;
    }
    default: {
      printErrorAndDie("Unrecognized CIGAR char in calc_seed_base()");
    }
    }
  }

  // Verify seed validity
  if (best_seed < -1 || best_seed == 0 || best_seed >= ((int)aln.get_sequence().size())-1)
    printErrorAndDie("Invalid alignment seed");
  return best_seed;
}

void HapAligner::process_reads(const std::vector<Alignment>& alignments,
			       std::vector<AlignmentMatrixCache*>& fw_matrix_caches,
			       std::vector<AlignmentMatrixCache*>& rv_matrix_caches,
			       const BaseQuality* base_quality, const std::vector<bool>& realign_read,
			       double* aln_probs, int* seed_positions){
  assert(alignments.size() == realign_read.size());
  assert((alignments.size() == fw_matrix_caches.size()) && (alignments.size() == rv_matrix_caches.size()));

  AlignmentTrace trace(fw_haplotype_->num_blocks());
  double* prob_ptr = aln_probs;
  for (unsigned int i = 0; i < alignments.size(); ++i){
    if (!realign_read[i]){
      prob_ptr += fw_haplotype_->num_combs();
      continue;
    }

    int seed_base = calc_seed_base(alignments[i]);
    seed_positions[i] = seed_base;
    if (seed_base == -1){
      // Assign all haplotypes the same zero LL
      for (unsigned int i = 0; i < fw_haplotype_->num_combs(); ++i, ++prob_ptr)
	*prob_ptr = 0;
    }
    else {
      process_read(alignments[i], seed_base, base_quality, false, fw_matrix_caches[i], rv_matrix_caches[i],
		   prob_ptr, trace, false);
      prob_ptr += fw_haplotype_->num_combs();
    }
  }
}

void HapAligner::process_read(const Alignment& aln, int seed, const BaseQuality* base_quality, bool retrace_aln,
			      AlignmentMatrixCache* fw_cache, AlignmentMatrixCache* rv_cache,
			      double* prob_ptr, AlignmentTrace& trace, bool debug){
  assert(seed != -1);
  assert(aln.get_sequence().size() > 0);
  assert(aln.get_sequence().size() == aln.get_base_qualities().size());
  const int seq_len = (int)aln.get_sequence().size();

  // Extract probabilities related to base quality scores for fw and rv configurations
  double* fw_log_wrong = new double[seq_len]; // log10(Prob(error))
  double* fw_log_right = new double[seq_len]; // log10(Prob(correct))
  double* rv_log_wrong = new double[seq_len];
  double* rv_log_right = new double[seq_len];
  const char* fw_seq   = aln.get_sequence().c_str();
  char* rv_seq         = new char[seq_len];
  const std::string& qual_string = aln.get_base_qualities();
  int k = seq_len-1;
  for (int j = 0; j < seq_len; ++j, --k){
    fw_log_wrong[j] = base_quality->log_prob_error(qual_string[j]);
    fw_log_right[j] = base_quality->log_prob_correct(qual_string[j]);
    rv_log_wrong[k] = fw_log_wrong[j];
    rv_log_right[k] = fw_log_right[j];
    rv_seq[k]       = fw_seq[j];
  }

  // True iff we should reuse alignment information from the previous haplotype to accelerate computations
  bool reuse_alns = false;

  // Generate data structures we'll use to perform and keep track of the Fw and Rv alignments
  AlignmentState fw_state(fw_haplotype_, fw_seq, seq_len, fw_log_wrong, fw_log_right, seed);
  AlignmentState rv_state(rv_haplotype_, rv_seq, seq_len, rv_log_wrong, rv_log_right, seq_len-1-seed);

  double max_LL = -100000000;
  do {
    if (!realign_to_hap_[fw_haplotype_->cur_index()]){
      prob_ptr++;
      reuse_alns = false;
      continue;
    }

    // Perform alignment to current haplotype and compute the LL of the optimal alignment
    fw_state.align_seq_to_haplotype(reuse_alns, fw_cache);
    rv_state.align_seq_to_haplotype(reuse_alns, rv_cache);
    double LL = calc_maximum_likelihood_alignment(fw_state, rv_state);
    
    *prob_ptr = LL;
    prob_ptr++;
    reuse_alns = true;

    if (LL > max_LL){
      max_LL = LL;
      if (retrace_aln)
	stitch(fw_state, rv_state, aln, trace, debug);
    }
  } while (fw_haplotype_->next() && rv_haplotype_->next());
  fw_haplotype_->reset(); rv_haplotype_->reset();

  // Deallocate arrays
  delete [] fw_log_wrong;
  delete [] fw_log_right;
  delete [] rv_log_wrong;
  delete [] rv_log_right;
  delete [] rv_seq;

  // Clear the matrix caches if appropriate
  // TO DO: Add clause to control when this occurs
  fw_cache->clear();
  rv_cache->clear();
}

AlignmentTrace* HapAligner::trace_optimal_aln(const Alignment& orig_aln, int seed, int best_haplotype, const BaseQuality* base_quality,
					      AlignmentMatrixCache* fw_cache, AlignmentMatrixCache* rv_cache,
					      bool debug){
  fw_haplotype_->go_to(best_haplotype);
  fw_haplotype_->fix();
  rv_haplotype_->go_to(best_haplotype);
  fw_haplotype_->fix();
  double prob;
  AlignmentTrace* trace = new AlignmentTrace(fw_haplotype_->num_blocks());
  process_read(orig_aln, seed, base_quality, true, fw_cache, rv_cache, &prob, *trace, debug);
  fw_haplotype_->unfix();
  rv_haplotype_->unfix();
  return trace;
}
