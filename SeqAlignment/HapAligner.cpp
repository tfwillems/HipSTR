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

void HapAligner::align_left_flank(const char* seq_0, int seq_len, 
				  const double* base_log_wrong, const double* base_log_correct,
				  double* match_matrix, double* insert_matrix, double& left_prob){
  // NOTE: Input matrix structure: Row = Haplotype position, Column = Read index
  double L_log_probs[seq_len];
 
  // Initialize first row of matrix (each base position matched with leftmost haplotype base)
  left_prob = 0.0;
  char first_hap_base = haplotype_->get_first_char();
  for (int j = 0; j < seq_len; ++j){
    match_matrix[j]  = (seq_0[j] == first_hap_base ? base_log_correct[j] : base_log_wrong[j]) + left_prob;
    insert_matrix[j] = base_log_correct[j] + left_prob;
    left_prob       += base_log_correct[j];
    L_log_probs[j]   = left_prob;
  }

  int haplotype_index = 1;
  int matrix_index    = seq_len;
  int stutter_R       = -1; // Haplotype index for right boundary of most recent stutter block

  // Fill in matrix row by row, iterating through each haplotype block
  for (int block_index = 0; block_index < haplotype_->num_blocks(); block_index++){
    const std::string& block_seq = haplotype_->get_seq(block_index);
    bool stutter_block           = (haplotype_->get_block(block_index)->get_repeat_info()) != NULL;

    // Skip any blocks to the left of the last changed block (as we can reuse the alignments)
    if (haplotype_->last_changed() != -1 && block_index < haplotype_->last_changed()){
      haplotype_index += block_seq.size() + (block_index == 0 ? -1 : 0);
      matrix_index     = seq_len*haplotype_index;
      if (stutter_block)
	stutter_R = haplotype_index - 1;
      continue;
    }

    if (stutter_block){
      // NOTE: Stutter block emission probability is invariant to direction of alignment, so we can align
      // left->right instead right->left as would be typically required

      RepeatStutterInfo* rep_info = haplotype_->get_block(block_index)->get_repeat_info();
      int period                  = rep_info->get_period();
      int block_option            = haplotype_->cur_index(block_index);
      int block_len               = block_seq.size();
      int prev_row_index          = seq_len*(haplotype_index-1);            // Index into matrix for haplotype character preceding stutter block (column = 0) 
      matrix_index                = seq_len*(haplotype_index+block_len-1);  // Index into matrix for rightmost character in stutter block (column = 0)
      int num_stutter_artifacts   = (rep_info->max_insertion()-rep_info->max_deletion())/period + 1;

      std::vector<double> block_probs(num_stutter_artifacts); // Reuse in each iteration to avoid reallocation penalty
      for (int j = 0; j < seq_len; ++j, ++matrix_index){
	// Consider valid range of insertions and deletions, including no stutter artifact
	int art_idx = 0;
	for (int artifact_size = rep_info->max_deletion(); artifact_size <= rep_info->max_insertion(); artifact_size += period){
	  int base_len    = std::min(block_len+artifact_size, j+1);
	  int ptr_offset  = j - base_len + 1;
	  double prob     = align_stutter_region(block_len, block_seq, base_len, seq_0+ptr_offset, base_log_wrong+ptr_offset, base_log_correct+ptr_offset, artifact_size);
	  double pre_prob = (j-base_len < 0 ? 0 : match_matrix[j-base_len + prev_row_index]);
	  block_probs[art_idx++] = rep_info->log_prob_pcr_artifact(block_option, artifact_size) + prob + pre_prob;
	}
	match_matrix[matrix_index] = log_sum_exp(block_probs);
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
	  homopolymer_len = std::min(MAX_HOMOP_LEN, haplotype_->homopolymer_length(block_index, coord_index));
	}
		  
	// Boundary conditions for leftmost base in read
	match_matrix[matrix_index]  = (seq_0[0] == hap_char ? base_log_correct[0] : base_log_wrong[0]);
	insert_matrix[matrix_index] = base_log_correct[0];
	matrix_index++;
	
	// Stutter block must be followed by a match
	if (haplotype_index == stutter_R + 1){
	  int prev_match_index = matrix_index - seq_len - 1;
	  for (int j = 1; j < seq_len; ++j, ++matrix_index, ++prev_match_index){
	    double match_emit           = (seq_0[j] == hap_char ? base_log_correct[j] : base_log_wrong[j]);
	    match_matrix[matrix_index]  = match_emit + match_matrix[prev_match_index];
	    insert_matrix[matrix_index] = IMPOSSIBLE;
	  }
	  continue;
	}

	std::vector<double> match_probs; match_probs.reserve(MAX_SEQ_DEL+1); // Reuse for each iteration to avoid reallocation penalty
	for (int j = 1; j < seq_len; ++j, ++matrix_index){
	  // Compute all match-related deletion probabilities (including normal read extension, where k = 1)
	  int del_index = matrix_index - seq_len - 1;
	  if (stutter_R == -1){
	    for (int k = 1; k <= std::min(haplotype_index, MAX_SEQ_DEL); k++){
	      match_probs.push_back(match_matrix[del_index]+LOG_DEL_N[homopolymer_len][k]);
	      del_index -= seq_len;
	    }
	    // Add deletion transitions to L flank state
	    for (int k = haplotype_index+1; k <= MAX_SEQ_DEL; k++)
	      match_probs.push_back(L_log_probs[j-1]+LOG_DEL_N[homopolymer_len][k]);
	  }
	  else {
	    // Add deletion transitions up until stutter_R+1
	    for (int k = 1; k <= std::min(haplotype_index-stutter_R-1, MAX_SEQ_DEL); k++){
	      match_probs.push_back(match_matrix[del_index]+LOG_DEL_N[homopolymer_len][k]);
	      del_index -= seq_len;
	    }
	  }

	  match_probs.push_back(insert_matrix[matrix_index-1]+LOG_MATCH_TO_INS[homopolymer_len]); // Add insertion-related probability
	  double match_emit           = (seq_0[j] == hap_char ? base_log_correct[j] : base_log_wrong[j]);
	  match_matrix[matrix_index]  = match_emit          + log_sum_exp(match_probs); 
	  insert_matrix[matrix_index] = base_log_correct[j] + log_sum_exp(match_matrix[matrix_index-seq_len-1]+LOG_INS_TO_MATCH, 
									  insert_matrix[matrix_index-1]+LOG_INS_TO_INS);
	  match_probs.clear();
	}	
      }
    }
  }
  assert(haplotype_index == haplotype_->cur_size());
}


void HapAligner::align_right_flank(const char* seq_n, int seq_len, 
				   const double* base_log_wrong, const double* base_log_correct,
				   double* match_matrix, double* insert_matrix, double& right_prob){
  // Note: Input matrix structure: Row = Haplotype position, Column = Read index
  double R_log_probs[seq_len];
 
  // Initialize first row of matrix (each base position matched with rightmost haplotype base)
  right_prob = 0.0;
  char last_hap_base = haplotype_->get_last_char();
  for (int j = 0; j < seq_len; ++j){
    match_matrix[j]  = (seq_n[-j] == last_hap_base ? base_log_correct[-j] : base_log_wrong[-j]) + right_prob;
    insert_matrix[j] = base_log_correct[-j] + right_prob;
    right_prob      += base_log_correct[-j];
    R_log_probs[j]   = right_prob;
  }

  int haplotype_index = 1;
  int matrix_index    = seq_len;
  int stutter_L       = -1; // Haplotype index for left boundary of most recent stutter block

  // Fill in matrix row by row, iterating through each haplotype block
  for (int block_index = haplotype_->num_blocks()-1; block_index >= 0; block_index--){
    const std::string& block_seq = haplotype_->get_seq(block_index);
    bool stutter_block           = haplotype_->get_block(block_index)->get_repeat_info() != NULL;

    // Skip any blocks to the right of the last changed block (as we can reuse the alignments)
    if (haplotype_->last_changed() != -1 && block_index > haplotype_->last_changed()){
      haplotype_index += block_seq.size() + (block_index == haplotype_->num_blocks()-1 ? -1 : 0);
      matrix_index     = seq_len*haplotype_index;
      if (stutter_block)
	stutter_L = haplotype_index - 1;
      continue;
    }

    if (stutter_block){
      RepeatStutterInfo* rep_info = haplotype_->get_block(block_index)->get_repeat_info();
      int period                  = rep_info->get_period();
      int block_option            = haplotype_->cur_index(block_index);
      int block_len               = block_seq.size();
      int prev_row_index          = seq_len*(haplotype_index-1);            // Index into matrix for haplotype character preceding stutter block (column = 0) 
      matrix_index                = seq_len*(haplotype_index+block_len-1);  // Index into matrix for rightmost character in stutter block (column = 0)
      int num_stutter_artifacts   = (rep_info->max_insertion()-rep_info->max_deletion())/period + 1;

      std::vector<double> block_probs(num_stutter_artifacts); // Reuse in each iteration to avoid reallocation penalty
      for (int j = 0; j < seq_len; ++j, ++matrix_index){
	// Consider valid range of insertions and deletions, including no stutter artifact
	int art_idx = 0;
	for (int artifact_size = rep_info->max_deletion(); artifact_size <= rep_info->max_insertion(); artifact_size += period){
	  int base_len    = std::min(block_len+artifact_size, j+1);
	  double prob     = align_stutter_region(block_len, block_seq, base_len, seq_n-j, base_log_wrong-j, base_log_correct-j, artifact_size);
	  double pre_prob = (j-base_len < 0 ? 0 : match_matrix[j-base_len + prev_row_index]);
	  block_probs[art_idx++] = rep_info->log_prob_pcr_artifact(block_option, artifact_size) + prob + pre_prob;
	}
	match_matrix[matrix_index] = log_sum_exp(block_probs);
      }
      
      // Adjust indices appropriately
      stutter_L         = haplotype_index + block_len - 1;
      haplotype_index  += block_len;
    }
    else {
      // Handle normal  n-> n-1 transitions while preventing sequencing indels from extending into preceding stutter blocks
      int coord_index = (block_index == haplotype_->num_blocks()-1 ? block_seq.size()-2 : block_seq.size()-1);
      char prev_char  = ' ';
      int homopolymer_len;

      for (; coord_index >= 0; --coord_index, ++haplotype_index){
	assert(matrix_index == seq_len*haplotype_index);
	char hap_char = block_seq[coord_index];
	
	// Update the homopolymer tract length whenever we encounter a new character
	if (hap_char != prev_char){
	  prev_char       = hap_char;
	  homopolymer_len = std::min(MAX_HOMOP_LEN, haplotype_->homopolymer_length(block_index, coord_index));
	}

	// Boundary conditions for rightmost base in read
	match_matrix[matrix_index]  = (seq_n[0] == hap_char ? base_log_correct[0] : base_log_wrong[0]);
	insert_matrix[matrix_index] = base_log_correct[0];
	matrix_index++;

	// Stutter block must be followed by a match
	if (haplotype_index == stutter_L + 1){
	  int prev_match_index = matrix_index - seq_len - 1;
	  for (int j = 1; j < seq_len; ++j, ++matrix_index, ++prev_match_index){
	    double match_emit           = (seq_n[-j] == hap_char ? base_log_correct[-j] : base_log_wrong[-j]);
	    match_matrix[matrix_index]  = match_emit + match_matrix[prev_match_index];
	    insert_matrix[matrix_index] = IMPOSSIBLE;
	  }
	  continue;
	}

	std::vector<double> match_probs; match_probs.reserve(MAX_SEQ_DEL+1); // Reuse for each iteration to avoid reallocation penalty
	for (int j = 1; j < seq_len; ++j, ++matrix_index){
	  // Compute all match-related deletion probabilities (including normal read extension, where k = 1)
	  int del_index = matrix_index - 1 - seq_len;
	  if (stutter_L == -1){
	    for (int k = 1; k <= std::min(haplotype_index, MAX_SEQ_DEL); k++){
	      match_probs.push_back(match_matrix[del_index]+LOG_DEL_N[homopolymer_len][k]);
	      del_index -= seq_len;
	    }
	    // Add deletion transitions to R flank state
	    for (int k = haplotype_index+1; k <= MAX_SEQ_DEL; k++)
	      match_probs.push_back(R_log_probs[j-1]+LOG_DEL_N[homopolymer_len][k]);
	  }
	  else {
	    // Add deletion transitions up until stutter_L+1
	    for (int k = 1; k <= std::min(haplotype_index-stutter_L-1, MAX_SEQ_DEL); k++){
	      match_probs.push_back(match_matrix[del_index]+LOG_DEL_N[homopolymer_len][k]);
	      del_index -= seq_len;
	    }
	  }

	  match_probs.push_back(insert_matrix[matrix_index-1]+LOG_MATCH_TO_INS[homopolymer_len]); // Add insertion-related probability
	  double match_emit           = (seq_n[-j] == hap_char ? base_log_correct[-j] : base_log_wrong[-j]);
	  match_matrix[matrix_index]  = match_emit           + log_sum_exp(match_probs); 
	  insert_matrix[matrix_index] = base_log_correct[-j] + log_sum_exp(match_matrix[matrix_index-seq_len-1]+LOG_INS_TO_MATCH, 
									   insert_matrix[matrix_index-1]+LOG_INS_TO_INS);
	  match_probs.clear();
	}
      }
    }
  }
  assert(haplotype_index == haplotype_->cur_size());
}

double HapAligner::compute_aln_logprob(int base_seq_len, int seed_base, 
				       char seed_char, double log_seed_wrong, double log_seed_correct,
				       double* l_match_matrix, double* l_insert_matrix, double l_prob,
				       double* r_match_matrix, double* r_insert_matrix, double r_prob){
  int lflank_len = seed_base;
  int rflank_len = base_seq_len-seed_base-1;
  int hapsize    = haplotype_->cur_size();

  // TO DO: Set this based on number of viable positions
  double SEED_LOG_MATCH_PRIOR = 0.0;

  std::vector<double> log_probs; 
  // Left flank entirely outside of haplotype window, seed aligned with 0   
  log_probs.push_back(SEED_LOG_MATCH_PRIOR + (seed_char == haplotype_->get_first_char() ? log_seed_correct: log_seed_wrong)
		      + l_prob
		      + LOG_DEL_N[haplotype_->homopolymer_length(0, 0)][1] + r_match_matrix[rflank_len*(hapsize-1)-1]); 

  // Right flank entirely outside of haplotype window, seed aligned with n-1
  int last_block = haplotype_->num_blocks()-1;
  int char_index = haplotype_->get_seq(last_block).size()-1;
  log_probs.push_back(SEED_LOG_MATCH_PRIOR + (seed_char == haplotype_->get_last_char() ? log_seed_correct: log_seed_wrong) 
		      + r_prob 
		      + LOG_DEL_N[haplotype_->homopolymer_length(last_block, char_index)][1] + l_match_matrix[lflank_len*(hapsize-1)-1]);
  
  // TO DO: Seed base outside of haplotype window


  // NOTE: Rationale for matrix indices:
  // lflank_len-1 with i-1 = lflank_len-1 with i-1 = lflank_len*(i-1) + lflank_len-1  = lflank_len*i - 1 ;
  // rflank_len-1 with i+1 = rflank_len-1 with hap_size-1-(i+1)       = rflank_len-1 with (hap_size-i-2) 
  //                       = rflank_len*(hap_size-i-2) + rflank_len-1 = rflank_len*(hap_size-i-1) -1;

  // Seed base aligned with each haplotype base
  double* l_match_ptr  = l_match_matrix  + (lflank_len - 1);
  double* r_match_ptr  = r_match_matrix  + (rflank_len*(hapsize-2) - 1);
  for (int block_index = 0; block_index < haplotype_->num_blocks(); block_index++){
    const std::string& block_seq = haplotype_->get_seq(block_index);
    bool stutter_block           = haplotype_->get_block(block_index)->get_repeat_info() != NULL;
    if (stutter_block){
      // Update matrix pointers
      l_match_ptr += lflank_len*block_seq.size();
      r_match_ptr -= rflank_len*block_seq.size();
      continue;
    }
    else {
      char prev_char  = ' ';
      int homopolymer_len;
      int coord_index     = (block_index == 0 ? 1 : 0); // Avoid situation where seed is aligned with first base
      int end_coord_index = (block_index == haplotype_->num_blocks()-1 ? block_seq.size()-1 : block_seq.size()); // Avoid situation where seed is aligned with last base
      for (; coord_index < end_coord_index; ++coord_index){
	if (block_seq[coord_index] != prev_char){
          prev_char       = block_seq[coord_index];
          homopolymer_len = std::min(MAX_HOMOP_LEN, haplotype_->homopolymer_length(block_index, coord_index));
        }

	log_probs.push_back(SEED_LOG_MATCH_PRIOR + (seed_char == block_seq[coord_index] ? log_seed_correct : log_seed_wrong)
			    + LOG_DEL_N[homopolymer_len][1] + *l_match_ptr
			    + LOG_DEL_N[homopolymer_len][1] + *r_match_ptr);
	l_match_ptr += lflank_len;
	r_match_ptr -= rflank_len;
      }
    }
  }
  double total_LL = log_sum_exp(log_probs);
  return total_LL;
}

/* 
 * Identify the base with the largest minimum distance from an insertion, a deletion and a stutter block
 * as defined by its alignment to the reference genome
 */
int HapAligner::calc_seed_base(Alignment& aln){  
  int32_t pos          = aln.get_start();
  int32_t repeat_start = region_->start() > 5 ? region_->start() - 5 : 0;
  int32_t repeat_stop  = region_->stop() + 5;
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
      if (region_->start() > max_ref_flank_len_)
	min_region = std::max(min_region, (int32_t)region_->start() - max_ref_flank_len_);
      max_region = std::min(max_region, (int32_t)region_->stop() + max_ref_flank_len_);

      if (min_region <= max_region){
	// Choose larger of valid two regions
	if (min_region < repeat_start && max_region >= repeat_stop)
	  if (repeat_start-1-min_region > max_region-repeat_stop){
	    min_region = min_region; max_region = repeat_start-1;
	  }
	  else {
	    min_region = repeat_stop; max_region = max_region;
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

void HapAligner::process_reads(std::vector<Alignment>& alignments, int init_read_index){
  for (unsigned int i = 0; i < alignments.size(); i++){
    int seed_base = calc_seed_base(alignments[i]);
    if (seed_base != -1){
      assert(alignments[i].get_sequence().size() == alignments[i].get_base_qualities().size());

      // Extract probabilites related to base quality scores
      double base_log_wrong[alignments[i].get_sequence().size()];   // log10(Prob(error))
      double base_log_correct[alignments[i].get_sequence().size()]; // log10(Prob(correct))
      const std::string& qual_string = alignments[i].get_base_qualities();
      for (int j = 0; j < qual_string.size(); j++){
	base_log_wrong[j]   = base_quality_->log_prob_error(qual_string[j]);
	base_log_correct[j] = base_quality_->log_prob_correct(qual_string[j]);
      }                                                                         

      const char* base_seq = alignments[i].get_sequence().c_str();
      int base_seq_len     = (int)alignments[i].get_sequence().size();
      int offset           = base_seq_len-1;

      // Allocate scoring matrices based on maximum haplotype size      
      int max_hap_size        = haplotype_->max_size();
      double* l_match_matrix  = new double [seed_base*max_hap_size];
      double* l_insert_matrix = new double [seed_base*max_hap_size];
      double* r_match_matrix  = new double [(base_seq_len-seed_base-1)*max_hap_size];
      double* r_insert_matrix = new double [(base_seq_len-seed_base-1)*max_hap_size];

      std::cerr << "Aligning read " << i+init_read_index << " " << alignments[i].get_sequence() << std::endl;
      //		<<  alignments[i].get_sequence().substr(0, seed_base) << " " <<  alignments[i].get_sequence().substr(seed_base) << std::endl;

      do {
	// Perform alignment to current haplotype
	double l_prob, r_prob;
	align_left_flank(base_seq, seed_base, base_log_wrong, base_log_correct, l_match_matrix, l_insert_matrix, l_prob);
	align_right_flank(base_seq+offset, base_seq_len-1-seed_base, base_log_wrong+offset, base_log_correct+offset, r_match_matrix, r_insert_matrix, r_prob);
	double LL = compute_aln_logprob(base_seq_len, seed_base, base_seq[seed_base], base_log_wrong[seed_base], base_log_correct[seed_base],
					l_match_matrix, l_insert_matrix, l_prob, r_match_matrix, r_insert_matrix, r_prob);
	//std::cerr << "\t" << LL << "\t";
	//haplotype_->print(std::cerr);
      } while (haplotype_->next());
      haplotype_->reset();
      //std::cerr << std::endl;

      // Deallocate scoring matrices
      delete [] l_match_matrix;
      delete [] l_insert_matrix;
      delete [] r_match_matrix;
      delete [] r_insert_matrix; 
    }
  }
  init_read_index++;
}
