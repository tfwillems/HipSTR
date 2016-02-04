#include <algorithm>
#include <climits>
#include <map>
#include <set>
#include <sstream>

#include "AlignmentModel.h"
#include "AlignmentTraceback.h"
#include "HapAligner.h"
#include "HapBlock.h"
#include "../mathops.h"
#include "RepeatBlock.h"
#include "StutterAlignerClass.h"


// Minimum distance of a seed base from an indel, mismatch or a repetitive region
const int32_t MIN_SEED_DIST = 5;

// Large negative value to prevent impossible or undesirable configurations 
const double IMPOSSIBLE = -1000000000;

void HapAligner::align_seq_to_hap(Haplotype* haplotype,
				  const char* seq_0, int seq_len, const double* base_log_wrong, const double* base_log_correct,
				  double* match_matrix, double* insert_matrix, double* deletion_matrix,
				  int* best_artifact_size, int* best_artifact_pos, double& left_prob){
  // NOTE: Input matrix structure: Row = Haplotype position, Column = Read index
  double* L_log_probs = new double[seq_len];
 
  // Initialize first row of matrix (each base position matched with leftmost haplotype base)
  left_prob = 0.0;
  char first_hap_base = haplotype->get_first_char();
  for (int j = 0; j < seq_len; ++j){
    match_matrix[j]    = (seq_0[j] == first_hap_base ? base_log_correct[j] : base_log_wrong[j]) + left_prob;
    insert_matrix[j]   = base_log_correct[j] + left_prob;
    deletion_matrix[j] = IMPOSSIBLE;
    left_prob         += base_log_correct[j];
    L_log_probs[j]     = left_prob;
  }

  int haplotype_index = 1;
  int matrix_index    = seq_len;
  int stutter_R       = -1; // Haplotype index for right boundary of most recent stutter block

  // Fill in matrix row by row, iterating through each haplotype block
  for (int block_index = 0; block_index < haplotype->num_blocks(); block_index++){
    const std::string& block_seq = haplotype->get_seq(block_index);
    bool stutter_block           = (haplotype->get_block(block_index)->get_repeat_info()) != NULL;

    // Skip any blocks to the left of the last changed block (as we can reuse the alignments)
    if (haplotype->last_changed() != -1 && block_index < haplotype->last_changed()){
      haplotype_index += block_seq.size() + (block_index == 0 ? -1 : 0);
      matrix_index     = seq_len*haplotype_index;
      if (stutter_block)
	stutter_R = haplotype_index - 1;
      continue;
    }

    if (stutter_block){
      RepeatStutterInfo* rep_info   = haplotype->get_block(block_index)->get_repeat_info();
      int period                    = rep_info->get_period();
      int block_option              = haplotype->cur_index(block_index);
      int block_len                 = block_seq.size();
      int prev_row_index            = seq_len*(haplotype_index-1);            // Index into matrix for haplotype character preceding stutter block (column = 0) 
      matrix_index                  = seq_len*(haplotype_index+block_len-1);  // Index into matrix for rightmost character in stutter block (column = 0)
      int num_stutter_artifacts     = (rep_info->max_insertion()-rep_info->max_deletion())/period + 1;
      const char* end_block_seq_arr = block_seq.c_str() + (block_seq.size()-1);
      StutterAlignerClass* stutter_aligner = haplotype->get_block(block_index)->get_stutter_aligner(block_option);

      /*
      // Precompute all match probabilites
      std::vector< std::vector<int> > suffix_match_probs(seq_len, std::vector<int>());
      for (int j = 0; j < seq_len; ++j){
	int min_index = std::max(-1, j-block_len);
	suffix_match_probs[j].reserve(j-min_index);
	double total  = 0;
	int ref_index = block_len-1;
	for (int k = j; k > min_index; --k, --ref_index){
	  total += (seq_0[k] == block_seq[ref_index] ? base_log_correct[k] : base_log_wrong[k]);
	  suffix_match_probs[j].push_back(total);
	}
      }
      */

      std::vector<double> block_probs(num_stutter_artifacts); // Reuse in each iteration to avoid reallocation penalty
      int j = 0;

      // If this haplotype and its predecessor have a suffix match that exceeds the maximum
      // haplotype bases used for a subset of the read position, we can reuse the match probabilities
      // for those positions as they're identical
      if (haplotype->last_changed() != -1){
	int suffix_match_length = haplotype->get_block(block_index)->suffix_match_len(block_option);
	int old_matrix_index    = seq_len*(haplotype_index+haplotype->get_block(block_index)->get_seq(block_option-1).size()-1);
	int num_copies          = std::min(seq_len, suffix_match_length + rep_info->max_deletion());
	for (; j < num_copies; ++j, ++matrix_index, ++old_matrix_index){
	  match_matrix[matrix_index]    = match_matrix[old_matrix_index];
	  insert_matrix[matrix_index]   = IMPOSSIBLE;
	  deletion_matrix[matrix_index] = IMPOSSIBLE;
	  // NOTE: No need to update artifact size and position as they're unchanged from last iteration
	}
      }

      for (; j < seq_len; ++j, ++matrix_index){
	// Consider valid range of insertions and deletions, including no stutter artifact
	int art_idx    = 0;
	double best_LL = IMPOSSIBLE;
	best_artifact_size[j] = -10000;
	for (int artifact_size = rep_info->max_deletion(); artifact_size <= rep_info->max_insertion(); artifact_size += period){
	  int art_pos          = -1;
	  int base_len         = std::min(block_len+artifact_size, j+1);
	  double prob          = stutter_aligner->align_stutter_region_reverse(base_len, seq_0+j, base_log_wrong+j, base_log_correct+j, artifact_size, art_pos);
	  double pre_prob      = (j-base_len < 0 ? 0 : match_matrix[j-base_len + prev_row_index]);
	  block_probs[art_idx] = rep_info->log_prob_pcr_artifact(block_option, artifact_size) + prob + pre_prob;
	  if (block_probs[art_idx] > best_LL){
	    best_artifact_size[j] = artifact_size;
	    best_artifact_pos[j]  = art_pos;
	    best_LL               = block_probs[art_idx];
	  }
	  art_idx++;
	}

	match_matrix[matrix_index]    = fast_log_sum_exp(block_probs);
	insert_matrix[matrix_index]   = IMPOSSIBLE;
	deletion_matrix[matrix_index] = IMPOSSIBLE;
      }
      
      // Adjust indices appropriately
      stutter_R         = haplotype_index + block_len - 1;
      haplotype_index  += block_len;
    }
    else {
      // Handle normal n -> n-1 transitions while preventing sequencing indels from extending into preceding stutter blocks
      int coord_index      = (block_index == 0 ? 1 : 0);
      int homopolymer_len  = haplotype->homopolymer_length(block_index, std::max(0, coord_index-1));

      for (; coord_index < block_seq.size(); ++coord_index, ++haplotype_index){
	assert(matrix_index == seq_len*haplotype_index);
	char hap_char = block_seq[coord_index];
	
	// Update the homopolymer tract length
	homopolymer_len = std::min(MAX_HOMOP_LEN, std::max(haplotype->homopolymer_length(block_index, coord_index),
							   haplotype->homopolymer_length(block_index, std::max(0, coord_index-1))));

	// Boundary conditions for leftmost base in read
	match_matrix[matrix_index]    = (seq_0[0] == hap_char ? base_log_correct[0] : base_log_wrong[0]);
	insert_matrix[matrix_index]   = (haplotype_index == stutter_R+1 ? IMPOSSIBLE : base_log_correct[0]);
	deletion_matrix[matrix_index] = (haplotype_index == stutter_R+1 ? IMPOSSIBLE :
					 std::max(deletion_matrix[matrix_index-seq_len]+LOG_DEL_TO_DEL, match_matrix[matrix_index-seq_len]+LOG_DEL_TO_MATCH));
	matrix_index++;
	
	// Stutter block must be followed by a match
	if (haplotype_index == stutter_R + 1){
	  int prev_match_index = matrix_index - seq_len - 1;
	  for (int j = 1; j < seq_len; ++j, ++matrix_index, ++prev_match_index){
	    double match_emit             = (seq_0[j] == hap_char ? base_log_correct[j] : base_log_wrong[j]);
	    match_matrix[matrix_index]    = match_emit + match_matrix[prev_match_index];
	    insert_matrix[matrix_index]   = IMPOSSIBLE;
	    deletion_matrix[matrix_index] = IMPOSSIBLE;
	  }
	  continue;
	}

	std::vector<double> match_probs; match_probs.reserve(3); // Reuse for each iteration to avoid reallocation penalty
	for (int j = 1; j < seq_len; ++j, ++matrix_index){
	  // Compute all match-related deletion probabilities (including normal read extension, where k = 1)
	  match_probs.push_back(insert_matrix[matrix_index-1]           + LOG_MATCH_TO_INS[homopolymer_len]);
	  match_probs.push_back(match_matrix[matrix_index-seq_len-1]    + LOG_MATCH_TO_MATCH[homopolymer_len]);
	  match_probs.push_back(deletion_matrix[matrix_index-seq_len-1] + LOG_MATCH_TO_DEL[homopolymer_len]);

	  double match_emit             = (seq_0[j] == hap_char ? base_log_correct[j] : base_log_wrong[j]);
	  match_matrix[matrix_index]    = match_emit          + std::max(match_probs[0], std::max(match_probs[1], match_probs[2]));
	  insert_matrix[matrix_index]   = base_log_correct[j] + std::max(match_matrix[matrix_index-seq_len-1] + LOG_INS_TO_MATCH,
									insert_matrix[matrix_index-1]         + LOG_INS_TO_INS);
	  deletion_matrix[matrix_index] = std::max(match_matrix[matrix_index-seq_len]    + LOG_DEL_TO_MATCH,
						   deletion_matrix[matrix_index-seq_len] + LOG_DEL_TO_DEL);
	  match_probs.clear();
	}	
      }
    }
  }
  delete [] L_log_probs;
  assert(haplotype_index == haplotype->cur_size());
}

double HapAligner::compute_aln_logprob(int base_seq_len, int seed_base,
				       char seed_char, double log_seed_wrong, double log_seed_correct,
				       double* l_match_matrix, double* l_insert_matrix, double* l_deletion_matrix, double l_prob,
				       double* r_match_matrix, double* r_insert_matrix, double* r_deletion_matrix, double r_prob,
				       int& max_index){
  int lflank_len = seed_base;
  int rflank_len = base_seq_len-seed_base-1;
  int hapsize    = fw_haplotype_->cur_size();

  // Compute number of viable seed positions and the corresponding uniform prior
  int num_seeds = 0;
  for (int block_index = 0; block_index < fw_haplotype_->num_blocks(); block_index++)
    if (fw_haplotype_->get_block(block_index)->get_repeat_info() == NULL)
      num_seeds += fw_haplotype_->get_seq(block_index).size();
  double SEED_LOG_MATCH_PRIOR = -int_log(num_seeds);
  
  double max_LL;
  std::vector<double> log_probs;
  // Left flank entirely outside of haplotype window, seed aligned with 0   
  log_probs.push_back(SEED_LOG_MATCH_PRIOR + (seed_char == fw_haplotype_->get_first_char() ? log_seed_correct: log_seed_wrong)
		      + l_prob + r_match_matrix[rflank_len*(hapsize-1)-1]);
  max_index = 0;
  max_LL    = log_probs[0];

  // Right flank entirely outside of haplotype window, seed aligned with n-1
  int last_block = fw_haplotype_->num_blocks()-1;
  int char_index = fw_haplotype_->get_seq(last_block).size()-1;
  log_probs.push_back(SEED_LOG_MATCH_PRIOR + (seed_char == fw_haplotype_->get_last_char() ? log_seed_correct: log_seed_wrong)
		      + r_prob + l_match_matrix[lflank_len*(hapsize-1)-1]);
  if (log_probs[1] > max_LL){
    max_index = fw_haplotype_->cur_size()-1;
    max_LL    = log_probs[1];
  }

  // NOTE: Rationale for matrix indices:
  // lflank_len-1 with i-1 = lflank_len-1 with i-1 = lflank_len*(i-1) + lflank_len-1  = lflank_len*i - 1 ;
  // rflank_len-1 with i+1 = rflank_len-1 with hap_size-1-(i+1)       = rflank_len-1 with (hap_size-i-2) 
  //                       = rflank_len*(hap_size-i-2) + rflank_len-1 = rflank_len*(hap_size-i-1) -1;

  // Seed base aligned with each haplotype base
  double* l_match_ptr  = l_match_matrix  + (lflank_len - 1);
  double* r_match_ptr  = r_match_matrix  + (rflank_len*(hapsize-2) - 1);
  int hap_index = 1;
  for (int block_index = 0; block_index < fw_haplotype_->num_blocks(); ++block_index){
    const std::string& block_seq = fw_haplotype_->get_seq(block_index);
    bool stutter_block           = fw_haplotype_->get_block(block_index)->get_repeat_info() != NULL;
    if (stutter_block){
      // Update matrix pointers
      l_match_ptr += lflank_len*block_seq.size();
      r_match_ptr -= rflank_len*block_seq.size();
      hap_index   += block_seq.size();
      continue;
    }
    else {
      int coord_index     = (block_index == 0 ? 1 : 0); // Avoid situation where seed is aligned with first base
      int end_coord_index = (block_index == fw_haplotype_->num_blocks()-1 ? block_seq.size()-1 : block_seq.size()); // Avoid situation where seed is aligned with last base
      for (; coord_index < end_coord_index; ++coord_index, ++hap_index){
	log_probs.push_back(SEED_LOG_MATCH_PRIOR + (seed_char == block_seq[coord_index] ? log_seed_correct : log_seed_wrong) + *l_match_ptr + *r_match_ptr);
	if (log_probs.back() > max_LL){
	  max_index = hap_index;
	  max_LL    = log_probs.back();
	}
	l_match_ptr += lflank_len;
	r_match_ptr -= rflank_len;
      }
    }
  }
  double total_LL = fast_log_sum_exp(log_probs);
  assert(total_LL < TOLERANCE);
  return total_LL;
}

/* 
 * Identify the base with the largest minimum distance from an insertion, a deletion and a stutter block
 * as defined by its alignment to the reference genome
 */
int HapAligner::calc_seed_base(Alignment& aln){
  assert(fw_haplotype_->num_blocks() == 3);
  int32_t pos          = aln.get_start();
  int32_t repeat_start = fw_haplotype_->get_block(1)->start();
  int32_t repeat_stop  = fw_haplotype_->get_block(1)->end();
  int best_seed = -1, cur_base = 0, max_dist = MIN_SEED_DIST;
  for (auto cigar_iter = aln.get_cigar_list().begin(); cigar_iter != aln.get_cigar_list().end(); cigar_iter++){
    switch(cigar_iter->get_type()){
    case '=': {
      int32_t min_region = pos, max_region = pos + cigar_iter->get_num() - 1; 
      if (min_region >= repeat_start)
	min_region = std::max(min_region, repeat_stop);
      if (max_region < repeat_stop)
	max_region = std::min(max_region, repeat_start-1);

      // Clip region so that the seed doesn't extend beyond the maximal sequence flanking the repeats
      // Only this region will be used to construct putative haplotypes
      min_region = std::max(min_region, fw_haplotype_->get_block(0)->start());
      max_region = std::min(max_region, fw_haplotype_->get_block(2)->end());

      if (min_region <= max_region){
	// Choose larger of valid two regions
	if (min_region < repeat_start && max_region >= repeat_stop){
	  if (repeat_start-1-min_region > max_region-repeat_stop){
	    //min_region = min_region;
	    max_region = repeat_start-1;
	  }
	  else {
	    min_region = repeat_stop;
	    //max_region = max_region;
	  }
	}
	
	int dist = 1 + (max_region-min_region)/2;
	if (dist >= max_dist){
	  max_dist  = dist;
	  best_seed = cur_base + (min_region-pos) + (max_region-min_region)/2;	  
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
      break;  
    }
    }
  }

  // Verify seed validity
  if (best_seed < -1 || best_seed == 0 || best_seed >= ((int)aln.get_sequence().size())-1)
    printErrorAndDie("Invalid alignment seed " + std::to_string(best_seed));
  return best_seed;
}

void HapAligner::process_reads(std::vector<Alignment>& alignments, int init_read_index, BaseQuality* base_quality,
			       double* aln_probs, int* seed_positions){
  AlignmentTrace trace(0);
  double* prob_ptr = aln_probs + (init_read_index*fw_haplotype_->num_combs());
  for (unsigned int i = 0; i < alignments.size(); i++){
    int seed_base = calc_seed_base(alignments[i]);
    seed_positions[init_read_index+i] = seed_base;
    if (seed_base == -1){
      // Assign all haplotypes the same zero LL
      for (unsigned int i = 0; i < fw_haplotype_->num_combs(); ++i, ++prob_ptr)
	*prob_ptr = 0;
    }
    else {
      process_read(alignments[i], seed_base, base_quality, false, prob_ptr, trace);
      prob_ptr += fw_haplotype_->num_combs();
    }
  }
}

const double TRACE_LL_TOL = 0.001;
inline int triple_min_index(double v1, double v2, double v3){
  if (v1 > v2+TRACE_LL_TOL)
    return (v1 > v3+TRACE_LL_TOL ? 0 : 2);
  else
    return (v2 > v3+TRACE_LL_TOL ? 1 : 2);
}

inline int rev_triple_min_index(double v1, double v2, double v3){
  if (v3 > v2+TRACE_LL_TOL)
    return (v3 > v1+TRACE_LL_TOL ? 2 : 0);
  else
    return (v2 > v1+TRACE_LL_TOL ? 1 : 0);
}

inline int     pair_min_index(double v1, double v2){ return (v1 > v2+TRACE_LL_TOL ? 0 : 1); }
inline int rev_pair_min_index(double v1, double v2){ return (v2 > v1+TRACE_LL_TOL ? 1 : 0); }

std::string HapAligner::retrace(Haplotype* haplotype, const char* read_seq,
				  int seq_len, int block_index, int base_index, int matrix_inde
				double* match_matrix, double* insert_matrix, double* deletion_matrix, int* best_artifact_size, int* best_artifact_pos,
				int& flank_ins_size, int& flank_del_size, int& stutter_size, std::string& str_seq,
				std::string& full_str_seq, std::vector< std::pair<int,int> >& flank_indel_data, bool right_to_left){
  const int MATCH = 0, DEL = 1, INS = 2, NONE = -1; // Types of matrices
  int seq_index   = seq_len-1;
  int matrix_type = MATCH;
  flank_ins_size  = 0;
  flank_del_size  = 0;
  std::stringstream aln_ss, str_ss, full_str_ss;

  int (*pair_index_fn)(double, double);
  int (*triple_index_fn)(double, double, double);

  if (right_to_left){
    pair_index_fn    = &pair_min_index;
    triple_index_fn  = &triple_min_index;
  }
  else {
    pair_index_fn   = &rev_pair_min_index;
    triple_index_fn = &rev_triple_min_index;
  }

  while (block_index >= 0){
    bool stutter_block = haplotype->get_block(block_index)->get_repeat_info() != NULL;
    if (stutter_block){
      const std::string& block_seq = haplotype->get_seq(block_index);
      int block_len = block_seq.size();
      stutter_size  = best_artifact_size[seq_index];
      assert(matrix_type == MATCH && base_index+1 == block_len);
      int i = 0;
      for (; i < std::min(seq_index+1, best_artifact_pos[seq_index]); i++){
	aln_ss      << "M";
	str_ss      << read_seq[seq_index-i];
	full_str_ss << read_seq[seq_index-i];
      }
      if (best_artifact_size[seq_index] < 0)
	aln_ss << std::string(-best_artifact_size[seq_index], 'D');
      else
	for (; i < std::min(seq_index+1, best_artifact_pos[seq_index] + best_artifact_size[seq_index]); i++){
	  aln_ss      << "I";
	  str_ss      << read_seq[seq_index-i];
	  full_str_ss << read_seq[seq_index-i];
	}
      for (; i < std::min(block_len + best_artifact_size[seq_index], seq_index+1); i++){
	aln_ss      << "M";
	str_ss      << read_seq[seq_index-i];
	full_str_ss << read_seq[seq_index-i];
      }
      str_seq = str_ss.str();

      // Add the non-spanned stutter block bases to generate what would be the full STR sequence if the read were longer
      // NOTE: Some weird edge case behavior can arise here if there's a stutter indel that's not spanned by the read.
      // In these very rare instances, the extracted string won't necessarily be correct
      int block_seq_index = std::min(block_len-1, block_len-1+best_artifact_size[seq_index]-i);
      while (block_seq_index >= 0)
	full_str_ss << block_seq[block_seq_index--];
      full_str_seq = full_str_ss.str();

      if (block_len + best_artifact_size[seq_index] >= seq_index+1)
	return aln_ss.str(); // Sequence doesn't span stutter block
      else {
	matrix_index -= (block_len + best_artifact_size[seq_index] + seq_len*block_len);
	matrix_type   = MATCH;
	seq_index    -= (block_len + best_artifact_size[seq_index]);
      }
    }
    else {
      int homopolymer_len     = haplotype->homopolymer_length(block_index, std::max(0, base_index-1));
      int prev_matrix_type    = NONE;
      std::string block_seq   = haplotype->get_seq(block_index);
      int32_t pos             = haplotype->get_block(block_index)->start() + (haplotype->reversed() ? -base_index : base_index);
      const int32_t increment = (haplotype->reversed() ? 1 : -1);
      int32_t indel_seq_index, indel_position;

      // Retrace flanks while tracking any indels that occur
      // Indels are ultimately reported as (position, size) tuples, where position is the left-most
      // start coordinate and sizes are +/- for insertions and deletions, respectively
      while (base_index >= 0 && seq_index >= 0){
	// Update the homopolymer tract length
	homopolymer_len = std::min(MAX_HOMOP_LEN, std::max(haplotype->homopolymer_length(block_index, base_index),
							   haplotype->homopolymer_length(block_index, std::max(0, base_index-1))));

	if (matrix_type != prev_matrix_type){
	  // Record any processed indels
	  if (prev_matrix_type == DEL){
	    int del_size = (haplotype->reversed() ? pos - indel_position : indel_position - pos);
	    if (haplotype->reversed())
	      flank_indel_data.push_back(std::pair<int,int>(indel_position, indel_position - pos));
	    else
	      flank_indel_data.push_back(std::pair<int,int>(pos+1, pos - indel_position));
	  }
	  else if (prev_matrix_type == INS)
	    flank_indel_data.push_back(std::pair<int,int>(indel_position + (haplotype->reversed() ? 0 : 1), indel_seq_index - seq_index));

	  // Initialize fields for indels that just began
	  if (matrix_type == DEL || matrix_type == INS){
	    indel_seq_index = seq_index;
	    indel_position  = pos;
	  }
	  prev_matrix_type = matrix_type;
        }

	// Extract alignment character for current base and update indices
	switch (matrix_type){
	case MATCH:
	  aln_ss << "M";
	  seq_index--;
	  base_index--;
	  pos += increment;
	  break;
	case DEL:
	  flank_del_size++;
	  aln_ss << "D";
	  base_index--;
	  pos += increment;
	  break;
	case INS:
	  flank_ins_size++;
	  aln_ss << "I";
	  seq_index--;
	  break;
	default:
	  printErrorAndDie("Invalid matrix type when retracing alignments");
	  break;
	}

	if (seq_index == -1)
	  return aln_ss.str();
	if (base_index == -1 && block_index == 0){
	  while(seq_index != -1){
	    aln_ss << "S";
	    seq_index--;
	  }
	  return aln_ss.str();
	}

	int best_opt;
	switch (matrix_type){
	case MATCH:
	  best_opt = triple_index_fn(insert_matrix[matrix_index-1]           + LOG_MATCH_TO_INS[homopolymer_len],
				     deletion_matrix[matrix_index-seq_len-1] + LOG_MATCH_TO_DEL[homopolymer_len],
				     match_matrix[matrix_index-seq_len-1]    + LOG_MATCH_TO_MATCH[homopolymer_len]);
	  if (best_opt == 0){
	    matrix_type   = INS;
	    matrix_index -= 1;
	  }
	  else if (best_opt == 1){
	    matrix_type   = DEL;
	    matrix_index -= (seq_len + 1);
	  }
	  else {
	    matrix_type   = MATCH;
	    matrix_index -= (seq_len + 1);
	  }
	  break;
	case DEL:
	  best_opt = pair_index_fn(deletion_matrix[matrix_index-seq_len] + LOG_DEL_TO_DEL,
				   match_matrix[matrix_index-seq_len]    + LOG_DEL_TO_MATCH);
	  if (best_opt == 0){
	    matrix_type   = DEL;
	    matrix_index -= seq_len;
	  }
	  else {
	    matrix_type   = MATCH;
	    matrix_index -= seq_len;
	  }
	  break;
	case INS:
	  best_opt = pair_index_fn(insert_matrix[matrix_index-1]        + LOG_INS_TO_INS,
				   match_matrix[matrix_index-seq_len-1] + LOG_INS_TO_MATCH);
	  if (best_opt == 0){
	    matrix_type   = INS;
	    matrix_index -= 1;
	  }
	  else {
	    matrix_type   = MATCH;
	    matrix_index -= (seq_len + 1);
	  }
	  break;
	default:
	  printErrorAndDie("Invalid matrix type when retracing alignments");
	  break;
	}
      }
    }
    base_index = haplotype->get_seq(--block_index).size()-1;
  }
  return aln_ss.str();
}

void HapAligner::process_read(Alignment& aln, int seed_base, BaseQuality* base_quality, bool retrace_aln,
			      double* prob_ptr, AlignmentTrace& trace){
  assert(seed_base != -1);
  assert(aln.get_sequence().size() == aln.get_base_qualities().size());

  // Extract probabilites related to base quality scores
  double base_log_wrong[aln.get_sequence().size()];   // log10(Prob(error))
  double base_log_correct[aln.get_sequence().size()]; // log10(Prob(correct))
  const std::string& qual_string = aln.get_base_qualities();
  for (unsigned int j = 0; j < qual_string.size(); j++){
    base_log_wrong[j]   = base_quality->log_prob_error(qual_string[j]);
    base_log_correct[j] = base_quality->log_prob_correct(qual_string[j]);
  }

  const char* base_seq = aln.get_sequence().c_str();
  int base_seq_len     = (int)aln.get_sequence().size();

  // Allocate scoring matrices based on the maximum haplotype size
  int max_hap_size          = fw_haplotype_->max_size();
  double* l_match_matrix    = new double [seed_base*max_hap_size];
  double* l_insert_matrix   = new double [seed_base*max_hap_size];
  double* l_deletion_matrix = new double [seed_base*max_hap_size];
  int* l_best_artifact_size = new int    [seed_base];
  int* l_best_artifact_pos  = new int    [seed_base];
  double* r_match_matrix    = new double [(base_seq_len-seed_base-1)*max_hap_size];
  double* r_insert_matrix   = new double [(base_seq_len-seed_base-1)*max_hap_size];
  double* r_deletion_matrix = new double [(base_seq_len-seed_base-1)*max_hap_size];
  int* r_best_artifact_size = new int    [(base_seq_len-seed_base-1)];
  int* r_best_artifact_pos  = new int    [(base_seq_len-seed_base-1)];
  double max_LL             = -100000000;

  // Reverse bases and quality scores for the right flank
  std::string rev_rseq = aln.get_sequence().substr(seed_base+1);
  std::reverse(rev_rseq.begin(), rev_rseq.end());
  std::reverse(base_log_wrong+seed_base+1,   base_log_wrong+base_seq_len);
  std::reverse(base_log_correct+seed_base+1, base_log_correct+base_seq_len);

  do {
    // Perform alignment to current haplotype
    double l_prob, r_prob;
    int max_index;
    align_seq_to_hap(fw_haplotype_, base_seq, seed_base, base_log_wrong, base_log_correct,
		     l_match_matrix, l_insert_matrix, l_deletion_matrix, l_best_artifact_size, l_best_artifact_pos, l_prob);

    align_seq_to_hap(rev_haplotype_, rev_rseq.c_str(), rev_rseq.size(), base_log_wrong+seed_base+1, base_log_correct+seed_base+1,
		     r_match_matrix, r_insert_matrix, r_deletion_matrix, r_best_artifact_size, r_best_artifact_pos, r_prob);
    
    double LL = compute_aln_logprob(base_seq_len, seed_base, base_seq[seed_base], base_log_wrong[seed_base], base_log_correct[seed_base],
				    l_match_matrix, l_insert_matrix, l_deletion_matrix, l_prob, r_match_matrix, r_insert_matrix, r_deletion_matrix, r_prob, max_index);
    *prob_ptr = LL;
    prob_ptr++;

    if (LL > max_LL){
      max_LL = LL;
      if (retrace_aln){
	// Retrace sequence to left of seed (if appropriate)
	assert(max_index >= 0 && max_index < fw_haplotype_->cur_size());
	int fw_seed_block, fw_seed_coord, l_flank_ins = -1, l_flank_del = -1, l_stutter_size = 0;
	std::vector< std::pair<int,int> > flank_indel_data;
	fw_haplotype_->get_coordinates(max_index, fw_seed_block, fw_seed_coord);
	assert(fw_seed_block != 1);
	int l_matrix_index = seed_base*max_index - 1;
	std::string left_aln, l_full_str_seq = "", l_str_seq = "";
	if (max_index == 0){
	  left_aln    = std::string(seed_base, 'S'); // Soft clip read to left of seed as it extends beyond haplotype. Don't retrace
	  l_flank_ins = l_flank_del = 0;
	}
	else {
	  if (fw_seed_coord == 0){
	    int prev_block_size = fw_haplotype_->get_seq(fw_seed_block-1).size();
	    left_aln = retrace(fw_haplotype_, base_seq, seed_base, fw_seed_block-1, prev_block_size-1, l_matrix_index, l_match_matrix, l_insert_matrix, l_deletion_matrix,
			       l_best_artifact_size, l_best_artifact_pos, l_flank_ins, l_flank_del, l_stutter_size, l_str_seq, l_full_str_seq, flank_indel_data, true);
	  }
	  else
	    left_aln = retrace(fw_haplotype_, base_seq, seed_base, fw_seed_block, fw_seed_coord-1, l_matrix_index, l_match_matrix, l_insert_matrix, l_deletion_matrix,
			       l_best_artifact_size, l_best_artifact_pos, l_flank_ins, l_flank_del, l_stutter_size, l_str_seq, l_full_str_seq, flank_indel_data, true);
	  assert(l_flank_ins != -1 && l_flank_del != -1);
	}
	std::reverse(left_aln.begin(), left_aln.end());   // Alignment is backwards for left flank
	std::reverse(l_str_seq.begin(), l_str_seq.end());
	std::reverse(l_full_str_seq.begin(), l_full_str_seq.end());
	assert(left_aln.size() - std::count(left_aln.begin(), left_aln.end(), 'D') == seed_base);

	// Retrace sequence to right of seed (if appropriate)
	int rev_max_index = fw_haplotype_->cur_size()-1-max_index;
	assert(rev_max_index >= 0 && rev_max_index < rev_haplotype_->cur_size());
	int rev_seed_block, rev_seed_coord, r_flank_ins = -1, r_flank_del = -1, r_stutter_size = 0;
	rev_haplotype_->get_coordinates(rev_max_index, rev_seed_block, rev_seed_coord);
	assert(rev_seed_block != -1);
	int r_matrix_index = (base_seq_len-1-seed_base)*rev_max_index - 1;
	std::string right_aln, r_full_str_seq = "", r_str_seq = "";
	if (rev_max_index == 0){
	  right_aln   = std::string(base_seq_len-1-seed_base, 'S'); // Soft clip read to right of seed as it extends beyond haplotype. Don't retrace
	  r_flank_ins = r_flank_del = 0;
	}
	else {
	  if (rev_seed_coord == 0){
	    int prev_block_size = rev_haplotype_->get_seq(rev_seed_block-1).size();
	    right_aln = retrace(rev_haplotype_, rev_rseq.c_str(), base_seq_len-1-seed_base, rev_seed_block-1, prev_block_size-1, r_matrix_index, r_match_matrix,
				r_insert_matrix, r_deletion_matrix, r_best_artifact_size, r_best_artifact_pos, r_flank_ins, r_flank_del,
				r_stutter_size, r_str_seq, r_full_str_seq, flank_indel_data, false);
	  }
	  else
	    right_aln = retrace(rev_haplotype_, rev_rseq.c_str(), base_seq_len-1-seed_base, rev_seed_block, rev_seed_coord-1, r_matrix_index, r_match_matrix,
				r_insert_matrix, r_deletion_matrix, r_best_artifact_size, r_best_artifact_pos, r_flank_ins, r_flank_del,
				r_stutter_size, r_str_seq, r_full_str_seq, flank_indel_data, false);
	  assert(r_flank_ins != -1 && r_flank_del != -1);
	}
	assert(right_aln.size() - std::count(right_aln.begin(), right_aln.end(), 'D') == base_seq_len-1-seed_base);
	assert(l_str_seq.empty() || r_str_seq.empty());

	trace.set_flank_ins_size(l_flank_ins + r_flank_ins);
	trace.set_flank_del_size(l_flank_del + r_flank_del);
	std::string read_aln_to_hap = left_aln + "M" + right_aln;
	trace.set_hap_aln(read_aln_to_hap);
	trace.set_str_seq(!l_str_seq.empty() ? l_str_seq : r_str_seq);
	trace.set_full_str_seq(!l_full_str_seq.empty() ? l_full_str_seq : r_full_str_seq);

	// Only one of the flanks enters the stutter block, so only one can have a non-zero size
	assert(l_stutter_size == 0 || r_stutter_size == 0);
	trace.set_stutter_size(l_stutter_size+r_stutter_size);
	trace.set_flank_indel_data(flank_indel_data);

	stitch_alignment_trace(fw_haplotype_->get_block(0)->start(), fw_haplotype_->get_aln_info(),
			       read_aln_to_hap, max_index, seed_base, aln, trace.traced_aln());
      }
    }
  } while (fw_haplotype_->next() && rev_haplotype_->next());
  fw_haplotype_->reset();
  rev_haplotype_->reset();

  // Deallocate scoring matrices
  delete [] l_match_matrix;
  delete [] l_insert_matrix;
  delete [] l_deletion_matrix;
  delete [] l_best_artifact_size;
  delete [] l_best_artifact_pos;
  delete [] r_match_matrix;
  delete [] r_insert_matrix;
  delete [] r_deletion_matrix;
  delete [] r_best_artifact_size;
  delete [] r_best_artifact_pos;
}

AlignmentTrace* HapAligner::trace_optimal_aln(Alignment& orig_aln, int seed_base, int best_haplotype, BaseQuality* base_quality){
  fw_haplotype_->go_to(best_haplotype);
  fw_haplotype_->fix();
  rev_haplotype_->go_to(best_haplotype);
  fw_haplotype_->fix();
  double prob;
  AlignmentTrace* trace = new AlignmentTrace(best_haplotype);
  process_read(orig_aln, seed_base, base_quality, true, &prob, *trace);
  fw_haplotype_->unfix();
  rev_haplotype_->unfix();
  return trace;
}
