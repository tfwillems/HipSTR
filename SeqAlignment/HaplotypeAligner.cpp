#include <algorithm>
#include <climits>
#include <map>
#include <set>

#include "AlignmentModel.h"
#include "HaplotypeAligner.h"
#include "HapBlock.h"
#include "../mathops.h"
#include "RepeatBlock.h"
#include "StutterAligner.h"

#include "constants.h"

// Large negative value used to prevent certain configurations
const double IMPOSSIBLE  = -1000000000;

// Minimum distance of a seed base from an indel or a stutter block
const int32_t MIN_SEED_DIST = 10;

void HaplotypeAligner::align_left_flank(const char* seq_0, 
					int seq_len, 
					const double* base_log_wrong,
					const double* base_log_correct,
					double* match_matrix, 
					double* insert_matrix,
					double& left_prob){
  
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

  int  haplotype_index = 1;
  int  matrix_index    = seq_len;
  int  stutter_R       = -1; // Haplotype index for right boundary of most recent stutter block

  // Fill in matrix row by row, iterating through each haplotype block
  for (int block_index = 0; block_index < haplotype_->num_blocks(); block_index++){
    const std::string& block_seq = haplotype_->get_seq(block_index);
    bool stutter_block           = haplotype_->get_block(block_index)->get_repeat_info() != NULL;

    if (stutter_block){
      // NOTE: Stutter block emission probability is invariant to direction of alignment, so we can align
      // left->right instead right->left as would be typically required

      RepeatStutterInfo* rep_info = haplotype_->get_block(block_index)->get_repeat_info();
      int period                  = rep_info->get_period();
      int block_option            = haplotype_->cur_index(block_index);
      int block_len               = block_seq.size();
      int prev_row_index          = seq_len*(haplotype_index-1);            // Index into matrix for haplotype character preceding stutter block (column = 0) 
      matrix_index                = seq_len*(haplotype_index+block_len-1);  // Index into matrix for rightmost character in stutter block (column = 0)

      for (int j = 0; j < seq_len; ++j, ++matrix_index){
	std::vector<double> block_probs;

	// Consider valid range of insertions and deletions, including no stutter artifact
	for (int artifact_size = rep_info->max_deletion(); artifact_size <= rep_info->max_insertion(); artifact_size += period){
	  int base_len    = std::min(block_len+artifact_size, j+1);
	  int ptr_offset  = j - base_len + 1;
	  double prob     = align_stutter_region(block_len, block_seq, base_len, seq_0+ptr_offset, base_log_wrong+ptr_offset, base_log_correct+ptr_offset, artifact_size);
	  double pre_prob = (j-base_len < 0 ? 0 : match_matrix[j-base_len + prev_row_index]);
	  block_probs.push_back(rep_info->log_prob_pcr_artifact(block_option, artifact_size) + prob + pre_prob);
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

	for (int j = 1; j < seq_len; ++j, ++matrix_index){
	  // Compute all match-related deletion probabilities (including normal read extension, where k = 1)
	  std::vector<double> match_probs;	
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
	    for (int k = 1; k  <= std::min(haplotype_index-stutter_R-1, MAX_SEQ_DEL); k++){
	      match_probs.push_back(match_matrix[del_index]+LOG_DEL_N[homopolymer_len][k]);
	      del_index -= seq_len;
	    }
	  }

	  match_probs.push_back(insert_matrix[matrix_index-1]+LOG_MATCH_TO_INS[homopolymer_len]); // Add insertion-related probability
	  double match_emit           = (seq_0[j] == hap_char ? base_log_correct[j] : base_log_wrong[j]);
	  match_matrix[matrix_index]  = match_emit          + log_sum_exp(match_probs); 
	  insert_matrix[matrix_index] = base_log_correct[j] + log_sum_exp(match_matrix[matrix_index-seq_len-1]+LOG_INS_TO_MATCH, 
									  insert_matrix[matrix_index-1]+LOG_INS_TO_INS);
	}
      }
    }
  }
}
void HaplotypeAligner::align_right_flank(const char* seq_n, 
					 int seq_len, 
					 const double* base_log_wrong, 
					 const double* base_log_correct,
					 double* match_matrix, 
					 double* insert_matrix,
					 double& right_prob){
  // Note: Input matrix structure: Row = Haplotype position, Column = Read index

  double R_log_probs[seq_len];
 
  // Initialize first row of matrix (each base position matched with rightmost haplotype base)
  right_prob = 0.0;
  char last_hap_base = haplotype_->get_last_char();
  for (int j = 0; j < seq_len; j++){
    match_matrix[j]  = (seq_n[-j] == last_hap_base ? base_log_correct[-j] : base_log_wrong[-j]) + right_prob;
    insert_matrix[j] = base_log_correct[-j] + right_prob;
    right_prob      += base_log_correct[-j];
    R_log_probs[j]   = right_prob;
  }

  int  haplotype_index = 1;
  int  matrix_index    = seq_len;
  int  stutter_L       = -1; // Haplotype index for left boundary of most recent stutter block

  // Fill in matrix row by row, iterating through each haplotype block
  for (int block_index = haplotype_->num_blocks()-1; block_index >= 0; block_index--){
    const std::string& block_seq = haplotype_->get_seq(block_index);
    bool stutter_block           = haplotype_->get_block(block_index)->get_repeat_info() != NULL;

    if (stutter_block){
      RepeatStutterInfo* rep_info = haplotype_->get_block(block_index)->get_repeat_info();
      int period                  = rep_info->get_period();
      int block_option            = haplotype_->cur_index(block_index);
      int block_len               = block_seq.size();
      int prev_row_index          = seq_len*(haplotype_index-1);            // Index into matrix for haplotype character preceding stutter block (column = 0) 
      matrix_index                = seq_len*(haplotype_index+block_len-1);  // Index into matrix for rightmost character in stutter block (column = 0)

      for (int j = 0; j < seq_len; ++j, ++matrix_index){
	std::vector<double> block_probs;

	// Consider valid range of insertions and deletions, including no stutter artifact
	for (int artifact_size = rep_info->max_deletion(); artifact_size <= rep_info->max_insertion(); artifact_size += period){
	  int base_len    = std::min(block_len+artifact_size, j+1);
	  double prob     = align_stutter_region(block_len, block_seq, base_len, seq_n-j, base_log_wrong-j, base_log_correct-j, artifact_size);
	  double pre_prob = (j-base_len < 0 ? 0 : match_matrix[j-base_len + prev_row_index]);
	  block_probs.push_back(rep_info->log_prob_pcr_artifact(block_option, artifact_size) + prob + pre_prob);
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

	for (int j = 1; j < seq_len; ++j, ++matrix_index){
	  // Compute all match-related deletion probabilities (including normal read extension, where k = 1)
	  std::vector<double> match_probs;	
	  int del_index = matrix_index - 1 - seq_len;
	  if (stutter_L == -1){
	    for (int k = 1; k <= std::min(haplotype_index, MAX_SEQ_DEL); k++){
	      match_probs.push_back(match_matrix[del_index]+LOG_DEL_N[homopolymer_len][k]);
	      del_index -= seq_len;
	    }
	    // Add deletion transitions to L flank state
	    for (int k = haplotype_index+1; k <= MAX_SEQ_DEL; k++)
	      match_probs.push_back(R_log_probs[j-1]+LOG_DEL_N[homopolymer_len][k]);
	  }
	  else {
	    // Add deletion transitions up until stutter_L+1
	    for (int k = 1; k  <= std::min(haplotype_index-stutter_L-1, MAX_SEQ_DEL); k++){
	      match_probs.push_back(match_matrix[del_index]+LOG_DEL_N[homopolymer_len][k]);
	      del_index -= seq_len;
	    }
	  }

	  match_probs.push_back(insert_matrix[matrix_index-1]+LOG_MATCH_TO_INS[homopolymer_len]); // Add insertion-related probability
	  double match_emit           = (seq_n[-j] == hap_char ? base_log_correct[-j] : base_log_wrong[-j]);
	  match_matrix[matrix_index]  = match_emit          + log_sum_exp(match_probs); 
	  insert_matrix[matrix_index] = base_log_correct[-j] + log_sum_exp(match_matrix[matrix_index-seq_len-1]+LOG_INS_TO_MATCH, 
									   insert_matrix[matrix_index-1]+LOG_INS_TO_INS);
	}
      }
    }
  }
}


double HaplotypeAligner::align(const char* base_seq, 
			       int base_seq_len, 
			       int seed_base,
			       const double* base_log_wrong, 
			       const double* base_log_correct,
			       double mapping_quality,
			       double* l_match_matrix, 
			       double* l_insert_matrix,
			       double* r_match_matrix, 
			       double* r_insert_matrix){
  return 0.0;


  int offset = base_seq_len-1;  
  double l_prob, r_prob;
  int l_flank_len = seed_base;
  int r_flank_len = offset-seed_base;
  
  // Align each flank
  align_left_flank(base_seq, l_flank_len, base_log_wrong, base_log_correct, l_match_matrix, l_insert_matrix, l_prob);

  
  return 0.0;


  align_right_flank(base_seq+offset, r_flank_len, base_log_wrong+offset, base_log_correct+offset, r_match_matrix, r_insert_matrix, r_prob);

  std::vector<double> log_probs;
  double SEED_LOG_MATCH_PRIOR, SEED_LOG_INSERT_PRIOR;

  // Left flank entirely outside of haplotype window, seed aligned with 0
  log_probs.push_back(SEED_LOG_MATCH_PRIOR + (base_seq[seed_base] == haplotype_->get_first_char() ? base_log_correct[seed_base]: base_log_wrong[seed_base])
		      + l_prob
		      + LOG_DEL_N[haplotype_->homopolymer_length(0, 1)][1] + r_match_matrix[r_flank_len*(haplotype_->cur_size()-1)-1]); 

  // Right flank entirely outside of haplotype window, seed aligned with n-1
  int last_block = haplotype_->num_blocks()-1;
  int char_index = haplotype_->get_seq(last_block).size()-2;
  log_probs.push_back(SEED_LOG_MATCH_PRIOR + (base_seq[seed_base] == haplotype_->get_last_char() ? base_log_correct[seed_base]: base_log_wrong[seed_base]) 
		      + r_prob 
		      + LOG_DEL_N[haplotype_->homopolymer_length(last_block, char_index)][1] + l_match_matrix[l_flank_len*(haplotype_->cur_size()-1)-1]);

  // Seed base outside of haplotype window
  if (mapping_quality != 1.0)
    log_probs.push_back(log(1.0-mapping_quality)); 

  // Seed base aligned with each haplotype base
  double* l_match_ptr  = l_match_matrix  + (l_flank_len - 1);
  double* r_match_ptr  = r_match_matrix  + (r_flank_len*(haplotype_->cur_size()-2) - 1);
  for (int block_index = 0; block_index < haplotype_->num_blocks(); block_index++){
    const std::string& block_seq = haplotype_->get_seq(block_index);
    bool stutter_block           = haplotype_->get_block(block_index)->get_repeat_info() != NULL;
    if (stutter_block){

    }
    else {
      char prev_char  = ' ';
      int homopolymer_len;

      int coord_index = (block_index == 0 ? 1 : 0);
      for (; coord_index < block_seq.size(); ++coord_index){
	if (block_seq[coord_index] != prev_char){
          prev_char       = block_seq[coord_index];
          homopolymer_len = std::min(MAX_HOMOP_LEN, haplotype_->homopolymer_length(block_index, coord_index));
        }

	log_probs.push_back(SEED_LOG_MATCH_PRIOR + (base_seq[seed_base] == block_seq[coord_index] ? base_log_correct[seed_base] : base_log_wrong[seed_base])
			    + LOG_DEL_N[homopolymer_len][1] + *l_match_ptr
			    + LOG_DEL_N[homopolymer_len][1] + *r_match_ptr);
	l_match_ptr += l_flank_len;
	r_match_ptr -= r_flank_len;
      }
    }
  }

  printErrorAndDie("Align not implemented");
  return log_sum_exp(log_probs);
}

/* 
 * Identify the base with the largest minimum distance from an insertion, a deletion and a stutter block
 * as defined by its alignment to the reference genome
 */
int HaplotypeAligner::calc_seed_base(Alignment& alignment){
  int32_t pos = alignment.get_start();
  const std::vector<CigarElement>& cigar_list = alignment.get_cigar_list();

  std::vector<int32_t> repeat_starts, repeat_stops;
  int best_seed = -1, cur_base = 0, max_dist = MIN_SEED_DIST;
  for (auto cigar_iter = cigar_list.begin(); cigar_iter != cigar_list.end(); cigar_iter++){
    switch(cigar_iter->get_type()){
    case '=': {
      int32_t min_region = pos, max_region = pos + cigar_iter->get_num() - 1; 
      for (int repeat_index = 0; repeat_index < repeat_starts.size(); repeat_index++){
	if (min_region >= repeat_starts[repeat_index])
	  min_region = std::max(min_region, repeat_stops[repeat_index]);    // rep_stop is not inclusive
	if (max_region < repeat_stops[repeat_index])
	  max_region = std::min(max_region, repeat_starts[repeat_index]-1); // rep_start is inclusive
      }
      if (min_region <= max_region){
	int dist = 1 + (max_region-min_region)/2;
	if (dist > max_dist){
	  max_dist  = dist;
	  best_seed = cur_base + ((max_region-min_region)/2 - (pos-min_region));
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
  if (best_seed < -1 || best_seed == 0 || best_seed >= ((int)alignment.get_sequence().size())-1){
    std::cerr << "'" << best_seed << "' " << alignment.get_sequence().size() << " " << alignment.get_sequence().size()-1 << std::endl;
    printErrorAndDie("Invalid alignment seed " + std::to_string(best_seed));
  }
  return best_seed;
}

void HaplotypeAligner::align(std::vector<Alignment>& alignments, 
			     BaseQuality& base_quality){
  // Determine the set of unique samples and 0-index them
  std::map<std::string, int> sample_index_mapping;
  std::vector<int> sample_indexes;
  int sample_index = 0;
  for (auto iter = alignments.begin(); iter != alignments.end(); iter++){
    auto match = sample_index_mapping.find(iter->get_sample());
    if (match == sample_index_mapping.end()){
      sample_indexes.push_back(sample_index);
      sample_index_mapping.insert(std::pair<std::string,int>(iter->get_sample(), sample_index));
      sample_index++;
    }
    else
      sample_indexes.push_back(match->second);
  }

#ifdef DEBUG
  std::cerr << "Aligning " << alignments.size() << " reads for " << sample_index << " unique samples" << std::endl;
  std::cerr << "Haplotype has " << haplotype_->num_combs() << " valid combinations" << std::endl;
#endif

  // Compute probability for each haplotype/alignment combination
  int ncombs = haplotype_->num_combs();
  int index  = 0;

  double* align_probs = new double[alignments.size()*ncombs];
  std::vector<bool> valid_read;
  
  for (int i = 0; i < alignments.size(); i++){
    // Compute base quality statistics for current read
    if(alignments[i].get_sequence().size() != alignments[i].get_base_qualities().size())
      printErrorAndDie("Lengths of sequence and base quality strings don't match");
    double base_log_wrong[alignments[i].get_sequence().size()];   // log10(Prob(error))
    double base_log_correct[alignments[i].get_sequence().size()]; // log10(Prob(correct))
    const std::string& qual_string = alignments[i].get_base_qualities();
    for (int j = 0; j < qual_string.size(); j++){
      base_log_wrong[j]   = base_quality.log_prob_error(qual_string[j]);
      base_log_correct[j] = base_quality.log_prob_correct(qual_string[j]);
    }
        
    // Determine seed base
    int seed_base = calc_seed_base(alignments[i]);
    if (seed_base == -1){
      // No valid seed, so ignore read
      valid_read.push_back(false);
      index += ncombs;
      continue;
    }
    else
      valid_read.push_back(true);

#ifdef DEBUG
    const std::string& seq = alignments[i].get_sequence();
    std::cerr << "Read split: " << seq.substr(0, seed_base) << " " << seq[seed_base] << " " << seq.substr(seed_base+1) << std::endl;
#endif
    
    const char* base_seq = alignments[i].get_sequence().c_str();
    int base_seq_len     = (int)alignments[i].get_sequence().size();

    // Allocate alignment matrices based on largest haplotype
    int max_hap_size        = haplotype_->max_size();
    double* l_match_matrix  = new double [(seed_base-1)*max_hap_size];
    double* l_insert_matrix = new double [(seed_base-1)*max_hap_size];
    double* r_match_matrix  = new double [(base_seq_len-seed_base)*max_hap_size];
    double* r_insert_matrix = new double [(base_seq_len-seed_base)*max_hap_size];

    do {  
      // Align to the current haplotype and store the resulting probability
      double prob = align(base_seq, base_seq_len, seed_base, base_log_wrong, base_log_correct, alignments[i].get_mapping_quality(),
			  l_match_matrix, l_insert_matrix,   r_match_matrix, r_insert_matrix);
      align_probs[index++] = prob;
    } while(haplotype_->next());
    haplotype_->reset();
    
    // Deallocate alignment matrices
    delete [] l_match_matrix;
    delete [] l_insert_matrix;
    delete [] r_match_matrix;
    delete [] r_insert_matrix;
  }
  
  // Compute combined likelihood for each haplotype pair and each sample using all reads
  int npairs          = (ncombs*(ncombs+1))/2;
  double* hap_probs   = new double[npairs*sample_index]();
  int*    read_counts = new int[sample_index]();
  int prob_index      = 0;
  for (int i = 0; i < alignments.size(); i++){
    if (!valid_read[i])
      continue;
        
    read_counts[sample_indexes[i]]++;
    int comb_index = sample_indexes[i]*npairs;
    for (int j = 0; j < ncombs; j++){
      double p1 = align_probs[prob_index+j];
      for (int k = j; k < ncombs; k++){
	double p2 = align_probs[prob_index+k];
	if (j == k)
	  hap_probs[comb_index++] += p1;
	else
	  hap_probs[comb_index++] += (LOG_ONE_HALF + log_sum_exp(p1, p2));

	// TO DO: Modify above expression when allelic dropout doesn't
	// occur at the same rate for both alleles
      }
    }
    prob_index += ncombs;
  }

  delete [] read_counts;
  delete [] hap_probs;
  delete [] align_probs;
}
