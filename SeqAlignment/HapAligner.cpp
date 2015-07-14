#include <algorithm>
#include <climits>
#include <map>
#include <set>

#include "AlignmentModel.h"
#include "HapAligner.h"
#include "HapBlock.h"
#include "../mathops.h"
#include "RepeatBlock.h"
#include "StutterAligner.h"

// Minimum distance of a seed base from an indel, mismatch or a repetitive region
const int32_t MIN_SEED_DIST = 5;

// Large negative value to prevent impossible or undesirable configurations 
const double IMPOSSIBLE = -1000000000;

void HapAligner::align_seq_to_hap(Haplotype* haplotype,
				  const char* seq_0, int seq_len, const double* base_log_wrong, const double* base_log_correct,
				  double* match_matrix, double* insert_matrix, double* deletion_matrix, int* best_artifact, double& left_prob){
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

      std::vector<double> block_probs(num_stutter_artifacts); // Reuse in each iteration to avoid reallocation penalty
      for (int j = 0; j < seq_len; ++j, ++matrix_index){
	// Consider valid range of insertions and deletions, including no stutter artifact
	int art_idx      = 0;
	best_artifact[j] = -10000;
	double best_LL   = IMPOSSIBLE;
	for (int artifact_size = rep_info->max_deletion(); artifact_size <= rep_info->max_insertion(); artifact_size += period){
	  int base_len         = std::min(block_len+artifact_size, j+1);
	  double prob          = align_stutter_region_reverse(block_len, end_block_seq_arr, base_len, seq_0+j, base_log_wrong+j, base_log_correct+j, artifact_size, period);
	  double pre_prob      = (j-base_len < 0 ? 0 : match_matrix[j-base_len + prev_row_index]);
	  block_probs[art_idx] = rep_info->log_prob_pcr_artifact(block_option, artifact_size) + prob + pre_prob;
	  if (block_probs[art_idx] > best_LL){
	    best_artifact[j] = artifact_size;
	    best_LL          = block_probs[art_idx];
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
      // Handle normal  n-> n-1 transitions while preventing sequencing indels from extending into preceding stutter blocks
      int coord_index = (block_index == 0 ? 1 : 0);
      char prev_char  = ' ';
      int homopolymer_len;

      for (; coord_index < block_seq.size(); ++coord_index, ++haplotype_index){
	assert(matrix_index == seq_len*haplotype_index);
	char hap_char = block_seq[coord_index];
	
	// Update the homopolymer tract length whenever we encounter a new character
	if (hap_char != prev_char){
	  prev_char       = hap_char;
	  homopolymer_len = std::min(MAX_HOMOP_LEN, haplotype->homopolymer_length(block_index, coord_index));
	}
		  
	// Boundary conditions for leftmost base in read
	match_matrix[matrix_index]    = (seq_0[0] == hap_char ? base_log_correct[0] : base_log_wrong[0]);
	insert_matrix[matrix_index]   = base_log_correct[0];
	deletion_matrix[matrix_index] = (haplotype_index == stutter_R+1 ? IMPOSSIBLE :
					 fast_log_sum_exp(deletion_matrix[matrix_index-seq_len]+LOG_DEL_TO_DEL, match_matrix[matrix_index-seq_len]+LOG_DEL_TO_MATCH));
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
	  match_matrix[matrix_index]    = match_emit          + fast_log_sum_exp(match_probs); 
	  insert_matrix[matrix_index]   = base_log_correct[j] + fast_log_sum_exp(match_matrix[matrix_index-seq_len-1] + LOG_INS_TO_MATCH,
										 insert_matrix[matrix_index-1]        + LOG_INS_TO_INS);
	  deletion_matrix[matrix_index] = fast_log_sum_exp(match_matrix[matrix_index-seq_len]    + LOG_DEL_TO_MATCH,
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
    if (fw_haplotype_->get_block(block_index)->get_repeat_info() != NULL)
      num_seeds += fw_haplotype_->get_seq(block_index).size();
  double SEED_LOG_MATCH_PRIOR = -log(num_seeds);
  
  double max_LL;
  std::vector<double> log_probs;
  // Left flank entirely outside of haplotype window, seed aligned with 0   
  log_probs.push_back(SEED_LOG_MATCH_PRIOR + (seed_char == fw_haplotype_->get_first_char() ? log_seed_correct: log_seed_wrong)
		      + l_prob
		      + LOG_MATCH_TO_MATCH[fw_haplotype_->homopolymer_length(0, 0)]
		      + r_match_matrix[rflank_len*(hapsize-1)-1]);
  max_index = 0;
  max_LL    = log_probs[0];

  // Right flank entirely outside of haplotype window, seed aligned with n-1
  int last_block = fw_haplotype_->num_blocks()-1;
  int char_index = fw_haplotype_->get_seq(last_block).size()-1;
  log_probs.push_back(SEED_LOG_MATCH_PRIOR + (seed_char == fw_haplotype_->get_last_char() ? log_seed_correct: log_seed_wrong)
		      + r_prob
		      + LOG_MATCH_TO_MATCH[fw_haplotype_->homopolymer_length(last_block, char_index)]
		      + l_match_matrix[lflank_len*(hapsize-1)-1]);
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
      char prev_char  = ' ';
      int homopolymer_len;
      int coord_index     = (block_index == 0 ? 1 : 0); // Avoid situation where seed is aligned with first base
      int end_coord_index = (block_index == fw_haplotype_->num_blocks()-1 ? block_seq.size()-1 : block_seq.size()); // Avoid situation where seed is aligned with last base
      for (; coord_index < end_coord_index; ++coord_index, ++hap_index){
	if (block_seq[coord_index] != prev_char){
          prev_char       = block_seq[coord_index];
          homopolymer_len = std::min(MAX_HOMOP_LEN, fw_haplotype_->homopolymer_length(block_index, coord_index));
        }

	log_probs.push_back(SEED_LOG_MATCH_PRIOR + (seed_char == block_seq[coord_index] ? log_seed_correct : log_seed_wrong)
			    + LOG_MATCH_TO_MATCH[homopolymer_len] + *l_match_ptr
			    + LOG_MATCH_TO_MATCH[homopolymer_len] + *r_match_ptr);
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
	    min_region = min_region; max_region = repeat_start-1;
	  }
	  else {
	    min_region = repeat_stop; max_region = max_region;
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
  if (best_seed < -1 || best_seed == 0 || best_seed >= ((int)aln.get_sequence().size())-1){
    std::cerr << "'" << best_seed << "' " << aln.get_sequence().size() << " " << aln.get_sequence().size()-1 << std::endl;
    printErrorAndDie("Invalid alignment seed " + std::to_string(best_seed));
  }
  return best_seed;
}

void HapAligner::process_reads(std::vector<Alignment>& alignments, int init_read_index, BaseQuality* base_quality,
			       double* aln_probs, int* seed_positions){
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
      process_read(alignments[i], seed_base, base_quality, false, prob_ptr);
      prob_ptr += fw_haplotype_->num_combs();
    }
  }
}


void HapAligner::retrace(Haplotype* haplotype,
			 int seq_len, int block_index, int base_index, int matrix_index,
			 double* match_matrix, double* insert_matrix, double* deletion_matrix, int* best_artifact){
  const int MATCH = 0, DEL = 1, INS = 2; // Types of matrices
  int seq_index   = seq_len-1;
  int matrix_type = MATCH;
  std::stringstream aln_ss;

  while (block_index >= 0){
    bool stutter_block = haplotype->get_block(block_index)->get_repeat_info() != NULL;
    if (stutter_block){
      int block_len = haplotype->get_seq(block_index).size();
      assert(matrix_type == MATCH && base_index+1 == block_len);
      if (block_len + best_artifact[seq_index] >= seq_index+1){
	// Flank doesn't span stutter block
	// TO DO: Extract stutter block alignment

	return;
      }
      else {
	// TO DO: Extract stutter block alignment
	// For now, just left-align all indels
	if (best_artifact[seq_index] < 0){
	  for (int i = 0; i < block_len + best_artifact[seq_index]; i++)
	    aln_ss << "M";
	  for (int i = 0; i > best_artifact[seq_index]; i--)
	    aln_ss << "D";
	}
	else {
	  for (int i = 0; i < block_len; i++)
	    aln_ss << "M";
	  for (int i = 0; i < best_artifact[seq_index]; i++)
	    aln_ss << "I";
	}

	matrix_index -= (block_len + best_artifact[seq_index] + seq_len*block_len);
	matrix_type   = MATCH;
	seq_index    -= (block_len + best_artifact[seq_index]);
      }
    }
    else {
      int homopolymer_len;
      char prev_char        = ' ';
      std::string block_seq = haplotype->get_seq(block_index);
      while (base_index >= 0 && seq_index >= 0){
	//std::cerr << "BASE=" << base_index << " SEQ=" << seq_index << std::endl;

	// Update the homopolymer tract length whenever we encounter a new character
	char hap_char = block_seq[base_index];
        if (hap_char != prev_char){
          prev_char       = hap_char;
          homopolymer_len = std::min(MAX_HOMOP_LEN, haplotype->homopolymer_length(block_index, base_index));
        }

	// Extract alignment character for current base and update indices
	switch (matrix_type){
	case MATCH:
	  aln_ss << "M";
	  seq_index--;
	  base_index--;
	  break;
	case DEL:
	  aln_ss << "D";
	  base_index--;
	  break;
	case INS:
	  aln_ss << "I";
	  seq_index--;
	  break;
	default:
	  printErrorAndDie("Invalid matrix type when retracing alignments");
	  break;
	}

	if (seq_index == -1)
	  return;
	if (base_index == -1 && block_index == 0){
	  while(seq_index != -1){
	    aln_ss << "S";
	    seq_index--;
	  }
	  return;
	}

	// Determine index for next element based on maximum likelihood transitions
	switch (matrix_type){
	case MATCH:
	  if (match_matrix[matrix_index-seq_len-1]    + LOG_MATCH_TO_MATCH[homopolymer_len] >
	      deletion_matrix[matrix_index-seq_len-1] + LOG_MATCH_TO_DEL[homopolymer_len]){
	    if (match_matrix[matrix_index-seq_len-1]  + LOG_MATCH_TO_MATCH[homopolymer_len] >
		insert_matrix[matrix_index-1]         + LOG_MATCH_TO_INS[homopolymer_len]){
	      matrix_type   = MATCH;
	      matrix_index -= (seq_len + 1);
	    }
	    else {
	      matrix_type   = INS;
	      matrix_index -= 1;
	    }
	  }
	  else {
	    if (deletion_matrix[matrix_index-seq_len-1] + LOG_MATCH_TO_DEL[homopolymer_len] >
		insert_matrix[matrix_index-1]           + LOG_MATCH_TO_INS[homopolymer_len]){
	      matrix_type   = DEL;
	      matrix_index -= (seq_len + 1);
	    }
	    else {
	      matrix_type   = INS;
	      matrix_index -= 1;
	    }
	  }
	  break;
	case DEL:
	  if (match_matrix[matrix_index-seq_len] + LOG_DEL_TO_MATCH >
	      deletion_matrix[matrix_index-seq_len] + LOG_DEL_TO_DEL){
	    matrix_type   = MATCH;
	    matrix_index -= seq_len;
	  }
	  else {
	    matrix_type   = DEL;
	    matrix_index -= seq_len;
	  }
	  break;
	case INS:
	  if (match_matrix[matrix_index-seq_len-1] + LOG_INS_TO_MATCH >
	      insert_matrix[matrix_index-1]        + LOG_INS_TO_INS){
	    matrix_type   = MATCH;
	    matrix_index -= (seq_len + 1);
	  }
	  else {
	    matrix_type   = INS;
	    matrix_index -= 1;
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
}

void HapAligner::process_read(Alignment& aln, int seed_base, BaseQuality* base_quality, bool retrace_aln,
			      double* prob_ptr){
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
  int offset           = base_seq_len-1;

  // Allocate scoring matrices based on the maximum haplotype size
  int max_hap_size          = fw_haplotype_->max_size();
  double* l_match_matrix    = new double [seed_base*max_hap_size];
  double* l_insert_matrix   = new double [seed_base*max_hap_size];
  double* l_deletion_matrix = new double [seed_base*max_hap_size];
  int* l_best_artifact      = new int    [seed_base];
  double* r_match_matrix    = new double [(base_seq_len-seed_base-1)*max_hap_size];
  double* r_insert_matrix   = new double [(base_seq_len-seed_base-1)*max_hap_size];
  double* r_deletion_matrix = new double [(base_seq_len-seed_base-1)*max_hap_size];
  int* r_best_artifact      = new int    [(base_seq_len-seed_base-1)];
  double max_LL = -100000000;


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
		     l_match_matrix, l_insert_matrix, l_deletion_matrix, l_best_artifact, l_prob);

    align_seq_to_hap(rev_haplotype_, rev_rseq.c_str(), rev_rseq.size(), base_log_wrong+seed_base+1, base_log_correct+seed_base+1,
		     r_match_matrix, r_insert_matrix, r_deletion_matrix, r_best_artifact, r_prob);
    
    double LL = compute_aln_logprob(base_seq_len, seed_base, base_seq[seed_base], base_log_wrong[seed_base], base_log_correct[seed_base],
				    l_match_matrix, l_insert_matrix, l_deletion_matrix, l_prob, r_match_matrix, r_insert_matrix, r_deletion_matrix, r_prob, max_index);
    *prob_ptr = LL;
    prob_ptr++;

    if (LL > max_LL){
      max_LL = LL;
      if (retrace_aln){
	  // lflank_len-1 with i-1 = lflank_len-1 with i-1 = lflank_len*(i-1) + lflank_len-1  = lflank_len*i - 1 ;
  // rflank_len-1 with i+1 = rflank_len-1 with hap_size-1-(i+1)       = rflank_len-1 with (hap_size-i-2) 
  //                       = rflank_len*(hap_size-i-2) + rflank_len-1 = rflank_len*(hap_size-i-1) -1;

	assert(max_index >= 0 && max_index < fw_haplotype_->cur_size());
	int seed_block, seed_coord;
	fw_haplotype_->get_coordinates(max_index, seed_block, seed_coord);
	//std::cerr << "Max index = " << max_index << std::endl;
	assert(seed_block != 1);
	int l_matrix_index = seed_base*max_index - 1;
	//int r_matrix_index = ()*() - 1;

	// TO DO: Retrace sequence to left of seed (if appropriate)
	if (max_index == 0){
	  // Soft clip read to left of seed as it extends beyond haplotype. Don't retrace

	}
	else {
	  if (seed_coord == 0){
	    int prev_block_size = fw_haplotype_->get_seq(seed_block-1).size();
	    retrace(fw_haplotype_, seed_base, seed_block-1, prev_block_size-1, l_matrix_index, l_match_matrix, l_insert_matrix, l_deletion_matrix, l_best_artifact);
	  }
	  else {
	    retrace(fw_haplotype_, seed_base, seed_block, seed_coord-1, l_matrix_index, l_match_matrix, l_insert_matrix, l_deletion_matrix, l_best_artifact);
	  }
	}

	/*

	  // NOTE: Have to use different get_coordinates function for right flank I believe

	// TO DO: Retrace sequence to right of seed (if appropriate)
	if (max_index == haplotype_->cur_size()-1){
	  // Soft clip read to right of seed as it extends beyond haplotype. Don't retrace

	}
	else {
	int matrix_index;
	  int cur_block_size = haplotype_->get_seq(seed_block).size();
	  retrace_right(base_seq_len-1-seed_base, (seed_coord == cur_block_size-1 ? seed_block+1 : seed_block), (seed_coord == cur_block_size-1 ? 0 : seed_coord+1),
			matrix_index, r_match_matrix, r_insert_matrix, r_deletion_matrix, r_best_artifact);
	}
	*/
      }
    }
  } while (fw_haplotype_->next() && rev_haplotype_->next());
  fw_haplotype_->reset();
  rev_haplotype_->reset();

  // Deallocate scoring matrices
  delete [] l_match_matrix;
  delete [] l_insert_matrix;
  delete [] l_deletion_matrix;
  delete [] l_best_artifact;
  delete [] r_match_matrix;
  delete [] r_insert_matrix;
  delete [] r_deletion_matrix;
  delete [] r_best_artifact;
}


Alignment HapAligner::trace_optimal_aln(Alignment& alignment, int seed_base, int best_haplotype, BaseQuality* base_quality){
  fw_haplotype_->go_to(best_haplotype);
  fw_haplotype_->fix();
  rev_haplotype_->go_to(best_haplotype);
  fw_haplotype_->fix();
  double probs[1];
  //process_read(alignment, seed_base, base_quality, true, probs);
  fw_haplotype_->unfix();
  rev_haplotype_->unfix();

  // TO DO: Report traced value

  return alignment;
}

