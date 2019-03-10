#include <time.h>

#include <algorithm>
#include <climits>
#include <sstream>

#include "AlignmentModel.h"
#include "AlignmentState.h"
#include "HapBlock.h"
#include "../mathops.h"
#include "RepeatBlock.h"
#include "StutterAligner.h"

// Initialize static class members
const double AlignmentState::IMPOSSIBLE               = -1000000000;
const double AlignmentState::MIN_SNP_LOG_PROB_CORRECT = -0.0043648054;

// Threshold below which two alignments are treated as have identical likelihoods
const double TRACE_LL_TOL = 0.001;

// Various functions we use to select from alignments with ~ identical likelihoods but different left alignment properties                 
inline int triple_min_index(double v1, double v2, double v3){
  if (v1 > v2+TRACE_LL_TOL)
    return (v1 > v3+TRACE_LL_TOL ? 0 : 2);
  else
    return (v2 > v3+TRACE_LL_TOL ? 1 : 2);
}
inline int rv_triple_min_index(double v1, double v2, double v3){
  if (v3 > v2+TRACE_LL_TOL)
    return (v3 > v1+TRACE_LL_TOL ? 2 : 0);
  else
    return (v2 > v1+TRACE_LL_TOL ? 1 : 0);
}
inline int    pair_min_index(double v1, double v2){ return (v1 > v2+TRACE_LL_TOL ? 0 : 1); }
inline int rv_pair_min_index(double v1, double v2){ return (v2 > v1+TRACE_LL_TOL ? 1 : 0); }

std::string AlignmentState::retrace_helper(AlignmentTrace& trace){
  if (seq_index_ == -1)
    return "";

  // Extract information about the alignment's current end position
  int seq_index    = seq_index_;
  int matrix_type  = matrix_type_;
  int matrix_index = matrix_index_; 
  int block_index, base_index;
  assert(hap_index_ >= 0 && hap_index_ < hap_->cur_size());
  hap_->get_coordinates(hap_index_, block_index, base_index);

  // Set function to select from alignments with ~ identical likelihoods
  int (*pair_index_fn)(double, double);
  int (*triple_index_fn)(double, double, double);
  if (!hap_->reversed()){
    pair_index_fn   = &pair_min_index;
    triple_index_fn = &triple_min_index;
  }
  else {
    pair_index_fn   = &rv_pair_min_index;
    triple_index_fn = &rv_triple_min_index;
  }

  std::stringstream aln_ss;
  while (block_index >= 0){
    const bool stutter_block     = hap_->get_block(block_index)->get_repeat_info() != NULL;
    const std::string& block_seq = hap_->get_seq(block_index);
    const int block_len          = block_seq.size();

    if (stutter_block){
      assert(matrix_type == MATCH && base_index+1 == block_len);
      int* artifact_size_ptr = best_artifact_size_ + seq_len_*block_index;
      int* artifact_pos_ptr  = best_artifact_pos_  + seq_len_*block_index;
      const int stutter_pos  = artifact_pos_ptr[seq_index];
      const int stutter_size = artifact_size_ptr[seq_index];

      std::stringstream str_ss;
      int i = 0;
      for (; i < std::min(seq_index+1, stutter_pos); ++i){
	aln_ss << "M";
	str_ss << seq_[seq_index-i];
      }
      if (stutter_size < 0)
	aln_ss << std::string(-stutter_size, 'D');
      else
	for (; i < std::min(seq_index+1, stutter_pos+stutter_size); ++i){
	  aln_ss << "I";
	  str_ss << seq_[seq_index-i];
	}
      for (; i < std::min(block_len + stutter_size, seq_index+1); ++i){
	aln_ss << "M";
	str_ss << seq_[seq_index-i];
      }

      // Add STR data to trace instance
      std::string str_seq = str_ss.str();
      if (hap_->reversed()){
	// Alignment for sequence to right of seed. Block indexes are reversed, but the alignment is correct
	trace.add_str_data(hap_->num_blocks()-1-block_index, stutter_size, str_seq);
      }
      else {
	// Alignment for sequence to left of seed. Block indexes are correct, but the alignment is reversed
	std::reverse(str_seq.begin(), str_seq.end());
	trace.add_str_data(block_index, stutter_size, str_seq);
      }

      int num_str_bases = (int)str_seq.size();
      if (num_str_bases == seq_index+1)
	return aln_ss.str(); // Sequence doesn't span stutter block
      else if (num_str_bases < seq_index+1) {
	matrix_index -= (num_str_bases + seq_len_*std::max(1, block_len));
	matrix_type   = MATCH;
	seq_index    -= (num_str_bases);
      }
      else
	assert(false);
    }
    else {
      int prev_matrix_type    = NONE;
      int32_t pos             = hap_->get_block(block_index)->start() + (hap_->reversed() ? -base_index : base_index);
      const int32_t increment = (hap_->reversed() ? 1 : -1);
      int32_t indel_seq_index = -1, indel_position = -1;

      // Retrace flanks while tracking any indels that occur
      // Indels are ultimately reported as (position, size) tuples, where position is the left-most
      // start coordinate and sizes are +/- for insertions and deletions, respectively
      std::stringstream flank_ss;
      while (base_index >= 0 && seq_index >= 0){
	// Update the homopolymer tract length
	int homopolymer_len = std::min(MAX_HOMOP_LEN, std::max(hap_->homopolymer_length(block_index, base_index),
							       hap_->homopolymer_length(block_index, std::max(0, base_index-1))));

	if (matrix_type != prev_matrix_type){
	  // Record any processed indels
	  if (prev_matrix_type == DEL){
	    if (hap_->reversed())
	      trace.add_flank_indel(std::pair<int32_t,int32_t>(indel_position, indel_position - pos));
	    else
	      trace.add_flank_indel(std::pair<int32_t,int32_t>(pos+1, pos - indel_position));
	  }
	  else if (prev_matrix_type == INS)
	    trace.add_flank_indel(std::pair<int32_t,int32_t>(indel_position + (hap_->reversed() ? 0 : 1), indel_seq_index - seq_index));

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
	  if (block_seq[base_index] != seq_[seq_index] && log_right_[seq_index] > MIN_SNP_LOG_PROB_CORRECT)
	    trace.add_flank_snp(pos, seq_[seq_index]);
	  flank_ss << seq_[seq_index];
	  aln_ss << "M";
	  seq_index--;
	  base_index--;
	  pos += increment;
	  break;
	case DEL:
	  trace.inc_flank_del();
	  aln_ss << "D";
	  base_index--;
	  pos += increment;
	  break;
	case INS:
	  trace.inc_flank_ins();
	  flank_ss << seq_[seq_index];
	  aln_ss << "I";
	  seq_index--;
	  break;
	default:
	  printErrorAndDie("Invalid matrix type when retracing alignments");
	  break;
	}

	if (seq_index == -1 || (base_index == -1 && block_index == 0)){
	  while(seq_index != -1){
	    aln_ss << "S";
	    seq_index--;
	  }

	  std::string flank_seq = flank_ss.str();
	  if (hap_->reversed())
	    trace.add_flank_data(hap_->num_blocks()-1-block_index, flank_seq);
	  else {
	    std::reverse(flank_seq.begin(), flank_seq.end());
	    trace.add_flank_data(block_index, flank_seq);
	  }
	  return aln_ss.str();
	}

	int best_opt;
	switch (matrix_type){
	case MATCH:
	  assert(matrix_index-seq_len_-1 >= 0);
	  
	  best_opt = triple_index_fn(ins_matrix_[matrix_index-1]            + LOG_MATCH_TO_INS[homopolymer_len],
				     del_matrix_[matrix_index-seq_len_-1]   + LOG_MATCH_TO_DEL[homopolymer_len],
				     match_matrix_[matrix_index-seq_len_-1] + LOG_MATCH_TO_MATCH[homopolymer_len]);
	  if (best_opt == 0){
	    matrix_type   = INS;
	    matrix_index -= 1;
	  }
	  else if (best_opt == 1){
	    matrix_type   = DEL;
	    matrix_index -= (seq_len_ + 1);
	  }
	  else {
	    matrix_type   = MATCH;
	    matrix_index -= (seq_len_ + 1);
	  }
	  break;
	case DEL:
	  best_opt = pair_index_fn(del_matrix_[matrix_index-seq_len_]   + LOG_DEL_TO_DEL,
				   match_matrix_[matrix_index-seq_len_] + LOG_DEL_TO_MATCH);
	  if (best_opt == 0){
	    matrix_type   = DEL;
	    matrix_index -= seq_len_;
	  }
	  else {
	    matrix_type   = MATCH;
	    matrix_index -= seq_len_;
	  }
	  break;
	case INS:
	  best_opt = pair_index_fn(ins_matrix_[matrix_index-1]            + LOG_INS_TO_INS,
				   match_matrix_[matrix_index-seq_len_-1] + LOG_INS_TO_MATCH);
	  if (best_opt == 0){
	    matrix_type   = INS;
	    matrix_index -= 1;
	  }
	  else {
	    matrix_type   = MATCH;
	    matrix_index -= (seq_len_ + 1);
	  }
	  break;
	default:
	  printErrorAndDie("Invalid matrix type when retracing alignments");
	  break;
	}
      }

      std::string flank_seq = flank_ss.str();
      if (hap_->reversed())
	trace.add_flank_data(hap_->num_blocks()-1-block_index, flank_seq);
      else {
	std::reverse(flank_seq.begin(), flank_seq.end());
	trace.add_flank_data(block_index, flank_seq);
      }
    }
    base_index = hap_->get_seq(--block_index).size()-1;
  }
  return aln_ss.str();
}

std::string AlignmentState::retrace(AlignmentTrace& trace){
  std::string res = retrace_helper(trace);
  if (!hap_->reversed()) // Alignment is backwards for Fw haplotypes
    std::reverse(res.begin(), res.end());
  assert(res.size() == (seq_index_ + 1 + std::count(res.begin(), res.end(), 'D')));
  return res;
}

std::string AlignmentState::stitch(const std::string& hap_aln, const std::string& read_aln, int h_index, int r_index, int increment){
  std::stringstream stitched_aln;
  while (r_index >= 0 && r_index < read_aln.size()){
    if (read_aln[r_index] == 'S'){
      stitched_aln << 'S';
      r_index += increment;
      continue;
    }

    assert(h_index >= 0 && h_index < hap_aln.size());
    if (hap_aln[h_index] == 'D'){
      if (read_aln[r_index] == 'I'){
	stitched_aln << 'M';
	r_index += increment;
	h_index += increment;
      }
      else {
	stitched_aln << 'D';
	h_index += increment;
      }
    }    
    else if (read_aln[r_index] == 'I'){
      stitched_aln << 'I';
      r_index += increment;
    }
    else if (read_aln[r_index] == 'D'){
      if (hap_aln[h_index] == 'M')
	stitched_aln << 'D';
      else if (hap_aln[h_index] == 'I')
	stitched_aln << "";
      else
	printErrorAndDie("Logical error in stitch()");
      r_index += increment;
      h_index += increment;
    }
    else if (read_aln[r_index] == 'M'){
      if (hap_aln[h_index] != 'M' && hap_aln[h_index] != 'I')
	printErrorAndDie("Logical error in stitch()");
      stitched_aln << hap_aln[h_index];
      r_index += increment;
      h_index += increment;
    }
    else
      printErrorAndDie("Logical error in stitch()");
   }
   return stitched_aln.str();
 }

void stitch(AlignmentState& fw_state, AlignmentState& rv_state, const Alignment& orig_aln,
	    AlignmentTrace& trace, bool debug){
   assert(fw_state.hap_index_ != -1 || rv_state.hap_index_ != -1);

   // Compute the full alignment of the read relative to the haplotype
   std::string left_aln        = fw_state.retrace(trace);
   std::string right_aln       = rv_state.retrace(trace);
   std::string read_aln_to_hap = left_aln + right_aln;

   if (debug)
     std::cerr << left_aln << " " << right_aln << " " << fw_state.seed_ << std::endl;

   // Determine the seed for the full alignment and the index of its matched haplotype base
   // NOTE: This is NOT the seed used during alignment, but rather a conveniently selected base that 
   // is matched to a haplotype base
   // Calculate the seed's index in the read vs. haplotype alignment string
   int seed_base, hap_index, read_aln_index;
   bool ok = false;
   if (fw_state.hap_index_ != -1){
     seed_base      = fw_state.seq_index_;
     hap_index      = fw_state.hap_index_;
     read_aln_index = (int)left_aln.size() - 1;

     // Find first matching base and use it as the seed
     for (auto iter = left_aln.rbegin(); iter != left_aln.rend(); ++iter, --read_aln_index){
       if (*iter == 'M'){
	 ok = true;
	 break;
       }
       
       if (*iter != 'D')
	 --seed_base;
       else
	 --hap_index;
     }
   }
   if (!ok && rv_state.hap_index_ != -1){
     seed_base      = rv_state.seq_len_         - rv_state.seq_index_ - 1;
     hap_index      = rv_state.hap_->cur_size() - rv_state.hap_index_ - 1;
     read_aln_index = left_aln.size();

     // Find first matching base and use it as the seed
     for (auto iter = right_aln.begin(); iter != right_aln.end(); ++iter, ++read_aln_index){
       if (*iter == 'M'){
	 ok = true;
	 break;
       }

       if (*iter != 'D')
	 ++seed_base;
       else
	 ++hap_index;
     }
   }
   assert(ok);
   assert(read_aln_index >= 0 && read_aln_index < read_aln_to_hap.size());

   // Retrieve pre-computed alignment of haplotype to reference genome
   std::string hap_aln_to_ref = fw_state.hap_->get_aln_info();

   // Determine the reference genome position of the seed
   // Calculate the seed's index in the haplotype vs. ref genome alignment string
   int32_t seed_pos  = fw_state.hap_->get_block(0)->start();
   int hap_aln_index = 0;
   while (hap_index > 0 && hap_aln_index < hap_aln_to_ref.size()){
     if (hap_aln_to_ref[hap_aln_index] == 'M' || hap_aln_to_ref[hap_aln_index] == 'I')
       hap_index--;
     if (hap_aln_to_ref[hap_aln_index] == 'M' || hap_aln_to_ref[hap_aln_index] == 'D')
       seed_pos++;
     hap_aln_index++;
   }
   while (hap_aln_index < hap_aln_to_ref.size() && hap_aln_to_ref[hap_aln_index] == 'D'){
     hap_aln_index++;
     // TO DO: Is there an error here? Shouldn't we be increasing seed_pos?
   }
   assert(hap_aln_index != hap_aln_to_ref.size());

   // Stitch the two sets of alignments together to generate read vs. ref genome alignments
   // Use the seed to separately handle the left and right sequences
   // Merge the two stitches together after omitting the double counted seed
   assert(read_aln_to_hap[read_aln_index] == 'M');
   std::string left_to_ref  = fw_state.stitch(hap_aln_to_ref, read_aln_to_hap, hap_aln_index, read_aln_index, -1);
   std::reverse(left_to_ref.begin(), left_to_ref.end());
   std::string right_to_ref = rv_state.stitch(hap_aln_to_ref, read_aln_to_hap, hap_aln_index, read_aln_index,  1);
   std::string full_to_ref  = left_to_ref + right_to_ref.substr(1);
   
   // Convert leading insertions to soft clips
   for (int i = 0; i < full_to_ref.size(); ++i){
     if (full_to_ref[i] == 'I')
       full_to_ref[i] = 'S';
     else
       break;
   }
   
   // Determine alignment start and end coordinates
   int32_t start = seed_pos+1; // +1 as left  stitch contains seed base
   int32_t stop  = seed_pos-1; // -1 as right stitch contains seed base
   for (auto iter = left_to_ref.begin(); iter != left_to_ref.end(); ++iter)
     if (*iter == 'D' || *iter == 'M')
       start--;
   for (auto iter = right_to_ref.begin(); iter != right_to_ref.end(); ++iter)
     if (*iter == 'D' || *iter == 'M')
       stop++;
   
   // Construct the CIGAR string elements
   std::vector<CigarElement> cigar_list;
   char cigar_char = full_to_ref[0];
   int  num        = 1;
   for(unsigned int i = 1; i < full_to_ref.size(); ++i){
     if (full_to_ref[i] != cigar_char){
       cigar_list.push_back(CigarElement(cigar_char, num));
       num = 1;
       cigar_char = full_to_ref[i];
     }
     else
       num += 1;
   }
   cigar_list.push_back(CigarElement(cigar_char, num));
   
   // Construct the actual alignment string from the string describing the alignment operations
   int read_index = 0;
   std::stringstream aln_ss;
   const std::string& bases = orig_aln.get_sequence();
   for (unsigned int i = 0; i < full_to_ref.size(); ++i){
     switch (full_to_ref[i]){
     case 'S':
       read_index++;
       break;
     case 'M':
     case 'I':
       aln_ss << bases[read_index];
       read_index++;
       break;
     case 'D':
       aln_ss << "-";
       break;
     default:
       printErrorAndDie("Invalid character encountered in stitch_alignment_trace()");
       break;
     }
   }

   // Update the fields in the trace instance
   trace.set_hap_aln(read_aln_to_hap);   
   trace.traced_aln() = Alignment(start, stop, false, "TRACE", orig_aln.get_base_qualities(), orig_aln.get_sequence(), aln_ss.str());
   trace.traced_aln().set_cigar_list(cigar_list);
}

void AlignmentState::align_seq_to_stutter_block(const int block_index,
						double* match_matrix, int* best_artifact_size, int* best_artifact_pos,
						const double* prev_match_matrix, const bool nonspanning_only){
  assert(hap_->get_block(block_index)->get_repeat_info() != NULL);
  assert((nonspanning_only && (prev_match_matrix == NULL)) || (!nonspanning_only && (prev_match_matrix != NULL)));
  double time = clock();

  RepeatStutterInfo* rep_info     = hap_->get_block(block_index)->get_repeat_info();
  const int period                = rep_info->get_period();
  const int block_option          = hap_->cur_index(block_index);
  const std::string& block_seq    = hap_->get_seq(block_index);
  const int block_len             = block_seq.size();
  const int num_stutter_artifacts = (rep_info->max_insertion()-rep_info->max_deletion())/period + 1;
  StutterAligner* stutter_aligner = hap_->get_block(block_index)->get_stutter_aligner(block_option);
  stutter_aligner->load_read(read_id_, seq_len_, seq_+seq_len_-1, log_wrong_+seq_len_-1, log_right_+seq_len_-1);

  // Maximum base in read to consider before we know there can't be any valid alignments
  const int max_base_index = (!nonspanning_only ? seq_len_ : std::min(seq_len_, block_len+rep_info->max_insertion()-1));

  int offset = seq_len_-1;
  int j      = 0;
  for (; j < max_base_index; ++j, --offset){
    // Consider valid range of insertions and deletions (including no stutter artifact)
    int art_idx           = 0;
    double best_LL        = IMPOSSIBLE;
    best_artifact_size[j] = -10000;
    double block_prob;
    for (int artifact_size = rep_info->max_deletion(); artifact_size <= rep_info->max_insertion(); artifact_size += period){
      int art_pos  = -1;
      int base_len = std::min(block_len+artifact_size, j+1);
      if (base_len >= 0){
	double prior_prob = rep_info->log_prob_pcr_artifact(block_option, artifact_size);
	if (j-base_len < 0)
	  block_prob = prior_prob +
	    stutter_aligner->align_stutter_region_reverse(base_len, seq_+j, offset, log_wrong_+j, log_right_+j, artifact_size, art_pos);
	else if (!nonspanning_only)
	  block_prob = prior_prob + prev_match_matrix[j-base_len] +
	    stutter_aligner->align_stutter_region_reverse(base_len, seq_+j, offset, log_wrong_+j, log_right_+j, artifact_size, art_pos);
	else
	  block_prob = IMPOSSIBLE;
      }
      else
	block_prob = IMPOSSIBLE;
      if (block_prob > best_LL){
	best_artifact_size[j] = artifact_size;
	best_artifact_pos[j]  = art_pos;
	best_LL               = block_prob;
      }
      art_idx++;
    }

    match_matrix[j] = best_LL;
  }

  // Populate values for bases where no valid alignments were attainable due to the restrictions
  for (; j < seq_len_; ++j){
    best_artifact_size[j] = -10000;
    match_matrix[j]       = IMPOSSIBLE;
  }
  total_stutter_time_ += (clock() - time)/CLOCKS_PER_SEC;
}

void AlignmentState::align_seq_to_nonstutter_block(const int block_index,
						   double* match_matrix, double* ins_matrix, double* del_matrix,
						   const double* prev_match_matrix){
  assert(hap_->get_block(block_index)->get_repeat_info() == NULL);
  double time = clock();

  const std::string& block_seq = hap_->get_seq(block_index);
  const int block_len          = block_seq.size();

  int coord_index = 0, matrix_index = 0;
  if (block_index == 0){
    // Initialize first row of matrix (each base position matched with the leftmost haplotype base)
    assert(prev_match_matrix == NULL);
    char first_hap_base = hap_->get_first_char();
    double left_prob    = 0.0;
    for (int j = 0; j < seq_len_; ++j){
      if (j == seed_){
	match_matrix[j] = IMPOSSIBLE; // Avoid situation where seed base is aligned with leftmost hap base
	ins_matrix[j]   = IMPOSSIBLE;
      }
      else {
	match_matrix[j] = (seq_[j] == first_hap_base ? log_right_[j] : log_wrong_[j]) + left_prob;
	ins_matrix[j]   = log_right_[j] + left_prob;
      }
      del_matrix[j]   = IMPOSSIBLE;
      left_prob       = ins_matrix[j];
    }
    coord_index  = 1;
    matrix_index = seq_len_;
  }
  else if (hap_->get_block(block_index-1)->get_repeat_info() != NULL){
    // Stutter block must be followed by a match
    assert(prev_match_matrix != NULL);
    char hap_char   = block_seq[0];
    match_matrix[0] = (seq_[0] == hap_char ? log_right_[0] : log_wrong_[0]);
    ins_matrix[0]   = IMPOSSIBLE;
    del_matrix[0]   = IMPOSSIBLE;
    for (int j = 1; j < seq_len_; ++j){
      double match_emit = (seq_[j] == hap_char ? log_right_[j] : log_wrong_[j]);
      match_matrix[j]   = match_emit + prev_match_matrix[j-1];
      ins_matrix[j]     = IMPOSSIBLE;
      del_matrix[j]     = IMPOSSIBLE;
    }

    coord_index  = 1;
    matrix_index = seq_len_;
  }
  else
    assert(false);
  
  // Handle normal n -> n-1 transitions
  for (; coord_index < block_len; ++coord_index){
    const char hap_char = block_seq[coord_index];

    // Update the homopolymer tract length
    int homopolymer_len = std::min(MAX_HOMOP_LEN, std::max(hap_->homopolymer_length(block_index, coord_index),
							   hap_->homopolymer_length(block_index, std::max(0, coord_index-1))));

    // Boundary conditions for leftmost base in read
    match_matrix[matrix_index] = (seq_[0] == hap_char ? log_right_[0] : log_wrong_[0]);
    ins_matrix[matrix_index]   = log_right_[0];
    del_matrix[matrix_index]   = std::max(  del_matrix[matrix_index-seq_len_]+LOG_DEL_TO_DEL, 
					  match_matrix[matrix_index-seq_len_]+LOG_DEL_TO_MATCH);

    // TO DO: Explore how removing this clause (and the clause below) affects results
    if (seed_ == 0)
      ins_matrix[matrix_index] = IMPOSSIBLE;
    matrix_index++;

    std::vector<double> match_probs; match_probs.reserve(3); // Reuse for each iteration to avoid reallocation penalty
    for (int j = 1; j < seq_len_; ++j, ++matrix_index){
      // Compute all match-related deletion probabilities (including normal read extension)
      match_probs.push_back(ins_matrix[matrix_index-1]            + LOG_MATCH_TO_INS[homopolymer_len]);
      match_probs.push_back(match_matrix[matrix_index-seq_len_-1] + LOG_MATCH_TO_MATCH[homopolymer_len]);
      match_probs.push_back(del_matrix[matrix_index-seq_len_-1]   + LOG_MATCH_TO_DEL[homopolymer_len]);

      double match_emit          = (seq_[j] == hap_char ? log_right_[j] : log_wrong_[j]);
      match_matrix[matrix_index] = match_emit    + std::max(match_probs[0], std::max(match_probs[1], match_probs[2]));
      ins_matrix[matrix_index]   = log_right_[j] + std::max(match_matrix[matrix_index-seq_len_-1]  + LOG_INS_TO_MATCH,
							    ins_matrix[matrix_index-1]             + LOG_INS_TO_INS);
      del_matrix[matrix_index]  = std::max(match_matrix[matrix_index-seq_len_] + LOG_DEL_TO_MATCH,
					   del_matrix[matrix_index-seq_len_]   + LOG_DEL_TO_DEL);

      // TO DO: Explore how removing this clause (and the clause above for the leftmost base) affects results
      // For now, clamp the seed base down on the haplotype by disallowing it to be inserted
      if (j == seed_)
	ins_matrix[matrix_index] = IMPOSSIBLE;

      match_probs.clear();
    }
  }
  total_nonstutter_time_ += (clock() - time)/CLOCKS_PER_SEC;
}

void AlignmentState::align_seq_to_haplotype(AlignmentMatrixCache* matrix_cache){
  double time = clock();

  // NOTE: Matrix structure: Row = Haplotype position, Column = Read index, fill in row-by-row
  // First row is alignment to haplotype position -1, which serves as a padding row
  // that prevents out-of-bounds indexes during traceback

  // Fill in padding row
  std::fill(match_matrix_, match_matrix_+seq_len_, IMPOSSIBLE);
  std::fill(  ins_matrix_,   ins_matrix_+seq_len_, IMPOSSIBLE);
  std::fill(  del_matrix_,   del_matrix_+seq_len_, IMPOSSIBLE);

  int haplotype_index = 0, matrix_index = seq_len_; // seq_len due to padding row
  const int num_blocks = hap_->num_blocks();

  bool reuse_alns = true;

  // Fill in matrix row by row by iterating through each haplotype block
  for (int block_index = 0; block_index < num_blocks; ++block_index){
    const std::string& block_seq = hap_->get_seq(block_index);
    const bool stutter_block     = (hap_->get_block(block_index)->get_repeat_info()) != NULL;
    const int block_len          = block_seq.size();
    const int block_option       = hap_->cur_index(block_index);

    reuse_alns &= (block_indexes_[block_index] == block_option);
    block_indexes_[block_index] = block_option;

    // Skip any blocks to the left of the last changed block (as we can reuse the alignments) and
    // skip any irrelevant non-stutter blocks for the Rv haplotype that will be exclusively handled by the Fw haplotype and vice-versa
    if (reuse_alns ||
	(!hap_->reversed() && !stutter_block && (block_index+1 == num_blocks)) ||
	( hap_->reversed() && !stutter_block && (block_index != 0))){
      haplotype_index += block_len;
      matrix_index    += std::max(1, block_len)*seq_len_;
      continue;
    }

    // For appropriate haplotype blocks, attempt to load the alignment data for the current block from the cache
    // If it does not exist, generate it and store it for later reuse and populate the main haplotype-level matrix as well
    // This is applicable whenever a haplotype block only involves non-spannning alignments, making it independent of all other blocks
    if (hap_->reversed() || block_index == 0){
      if (stutter_block){
	double *sub_match_matrix;
	int *sub_artifact_size, *sub_artifact_pos;
	int sub_matrix_size = seq_len_;
	if (matrix_cache->has(block_index, block_option))
	  matrix_cache->get(block_index, block_option, sub_match_matrix, sub_artifact_size, sub_artifact_pos);
	else {
	  sub_match_matrix  = new double[sub_matrix_size];
	  sub_artifact_size = new int[sub_matrix_size];
	  sub_artifact_pos  = new int[sub_matrix_size];
	  align_seq_to_stutter_block(block_index, sub_match_matrix, sub_artifact_size, sub_artifact_pos, NULL, true);
	  matrix_cache->add(block_index, block_option, sub_match_matrix, sub_artifact_size, sub_artifact_pos);
	}

	// Copy the results to the main arrays
	double* match_ptr      = match_matrix_ + matrix_index + seq_len_*std::max(0, block_len-1);
	int* artifact_size_ptr = best_artifact_size_ + seq_len_*block_index;
	int* artifact_pos_ptr  = best_artifact_pos_  + seq_len_*block_index;
	std::copy(sub_match_matrix,   sub_match_matrix+sub_matrix_size, match_ptr);
	std::copy(sub_artifact_size, sub_artifact_size+sub_matrix_size, artifact_size_ptr);
	std::copy(sub_artifact_pos,   sub_artifact_pos+sub_matrix_size, artifact_pos_ptr);

	// Prevent illegal alignments
	double* ins_ptr = ins_matrix_ + matrix_index + seq_len_*std::max(0, block_len-1);
	double* del_ptr = del_matrix_ + matrix_index + seq_len_*std::max(0, block_len-1);
	for (int i = 0; i < seq_len_; ++i)
	  ins_ptr[i]= del_ptr[i] = IMPOSSIBLE;
      }
      else
	align_seq_to_nonstutter_block(block_index, match_matrix_+matrix_index, ins_matrix_+matrix_index, del_matrix_+matrix_index, NULL);
    }
    else {
      if (stutter_block){
	int* artifact_size_ptr       = best_artifact_size_ + seq_len_*block_index;
	int* artifact_pos_ptr        = best_artifact_pos_  + seq_len_*block_index;
	const double* prev_match_ptr = match_matrix_       + matrix_index - seq_len_;
	double* match_ptr            = match_matrix_ + matrix_index + seq_len_*std::max(0, block_len-1);
	align_seq_to_stutter_block(block_index, match_ptr, artifact_size_ptr, artifact_pos_ptr, prev_match_ptr, false);

	// Prevent illegal alignments
	double* ins_ptr = ins_matrix_ + matrix_index + seq_len_*std::max(0, block_len-1);
	double* del_ptr = del_matrix_ + matrix_index + seq_len_*std::max(0, block_len-1);
	for (int i = 0; i < seq_len_; ++i)
	  ins_ptr[i]= del_ptr[i] = IMPOSSIBLE;
      }
      else
	align_seq_to_nonstutter_block(block_index, match_matrix_+matrix_index, ins_matrix_+matrix_index, del_matrix_+matrix_index,
				      match_matrix_+(matrix_index-seq_len_));
    }

    // Adjust indices
    haplotype_index += block_len;
    matrix_index    += seq_len_*std::max(1, block_len);
  }
  assert(haplotype_index == hap_->cur_size());
  total_hap_aln_time_ += (clock() - time)/CLOCKS_PER_SEC;
}

double calc_maximum_likelihood_alignment(AlignmentState& fw_state, AlignmentState& rv_state){
  assert(!fw_state.hap_->reversed() && rv_state.hap_->reversed());
  assert( fw_state.hap_->num_blocks() > 1);

  // Our key assumption/requirement is that the seed base must lie within the haplotype boundaries
  // In earlier HipSTR versions, we required that the seed base match a flanking base and be preceded and followed by matches
  // Here, we allow the seed base to lie anywhere within the haplotype (not just the flanks) 
  // and allow for it to be a match or insertion (vs. just a match)
  // There are 2 possibile scenarios:
  // Case I)  The read fully aligns within a single haplotype block
  //       If this is the case, the last base from the read much either be a match or an insertion
  //       This can occur at each position within the block. We use the Fw alignments to assess each of these configurations, 
  //        except for the last block where the Rv alignment is used 
  //       Note: One subcase is that the read is fully within a stutter block, which is not currently handled
  // Case II) The read is aligned across 2 or more contiguous haplotype blocks
  //       If this is the case, one base in the read must match at a boundary element between 
  //        the blocks (as we require matches at block transitions)
  //       We therefore examine all configurations in which any of a read's bases is matched to the boundary element
  //       This even allows the seed base to lie within a stutter block (which was previously disallowed) 
  //        if the whole read does not lie in the stutter block

  const int seq_len    = fw_state.seq_len_;
  const int num_blocks = fw_state.hap_->num_blocks(); 
  const int hap_size   = fw_state.hap_->cur_size();

  // Initialize values for optimal traceback configuration
  fw_state.reset_traceback(); rv_state.reset_traceback();

  // Initialize values for best solution. MAX_DIR = 0 if Fw, 1 if Rv, 2 if both Fw and Rv used
  double max_LL = AlignmentState::IMPOSSIBLE;
  int max_dir   = -1;

  // Step 1: Enumerate all possible configurations for Case I using subsets of the Fw and Rv alignment matrices
  for (int i = 0; i < 2; ++i){
    AlignmentState& state      = (i == 0 ? fw_state : rv_state);
    const Haplotype* haplotype = state.hap_;
    double* match_matrix       = state.match_matrix_;
    double* ins_matrix         = state.ins_matrix_;

    int hap_index    = 0;
    int matrix_index = seq_len-1 + seq_len; // +seq_len due to padding row
    for (int block_index = 0; block_index < num_blocks; ++block_index){
      const std::string& block_seq = haplotype->get_seq(block_index);
      const int block_len          = block_seq.size();
      const bool stutter_block     = haplotype->get_block(block_index)->get_repeat_info() != NULL;
      if (stutter_block){
	hap_index    += std::max(0, block_len-1);
	matrix_index += std::max(0, block_len-1)*seq_len;
	
	double v = match_matrix[matrix_index] + LOG_MATCH_TO_MATCH[1];
	if (v > max_LL){
	  max_LL  = v;
	  max_dir = i;
	  state.hap_index_    = hap_index;
	  state.seq_index_    = seq_len-1;
	  state.matrix_index_ = matrix_index;
	  state.matrix_type_  = AlignmentState::MATCH;
	}

	hap_index    += (block_len > 0 ? 1 : 0);
	matrix_index += seq_len;
	continue;
      }
      else {
	// Skip alignments in these region as they're handled elsewhere or by the other orientation
	// For the Fw haplotype, we consider all but the last haplotype block as that's covered by the Rv haplotype
	// For the Rv haplotype, we only consider the last block (which is the 0th block in its orientation)
	bool skip_block =  (haplotype->reversed() && block_index != 0);
	skip_block     |= (!haplotype->reversed() && block_index+1 == num_blocks);
	if (skip_block){
	  hap_index    += block_len;
	  matrix_index += seq_len*block_len;
	  continue;
	}

	for (int j = 0; j < block_len; ++j, ++hap_index){
	  // Although technically incorrect, we use these transition penalties (with a homopolymer length of 1) as a prior
	  // to penalize the states that end with an insertion. This penalty is not captured in the scoring matrices themselves
	  // b/c of the reverse nature of the alignment HMM (only transitions from the insertion state are captured)
	  double v1 = match_matrix[matrix_index] + LOG_MATCH_TO_MATCH[1];
	  double v2 = ins_matrix[matrix_index]   + LOG_MATCH_TO_INS[1];

	  if (v1 > max_LL){
	    max_LL  = v1;
	    max_dir = i;
	    state.hap_index_    = hap_index;
	    state.seq_index_    = seq_len-1;
	    state.matrix_index_ = matrix_index;
	    state.matrix_type_  = AlignmentState::MATCH;
	  }
	  if (j != 0 && v2 > max_LL){
	    max_LL  = v2;
	    max_dir = i;
	    state.hap_index_    = hap_index-1; // -1 b/c the insertion matrix captures insertions directly preceding the block base j
	    state.seq_index_    = seq_len-1;
	    state.matrix_index_ = matrix_index;
	    state.matrix_type_  = AlignmentState::INS;
	  }
	  
	  matrix_index += seq_len;
	}
      }
    }
  }

  // Reset the alignment state for the unused direction
  if (max_dir == 0)
    rv_state.reset_traceback();
  else if (max_dir == 1)
    fw_state.reset_traceback();
  else
    assert(false);

  // Step 2: Enumerate all possible configurations for Case II using both the Fw and Rv alignments
  // Consider all boundary elements and how the Fw and Rv alignments could merge at those elements

  // Use loop below instead of haplotype size to compute rv_matrix_index b/c some stutter blocks may have sequences of 0 size
  int fw_matrix_index = seq_len, rv_matrix_index = seq_len-1; // seq_len due to padding row
  for (int block_index = 0; block_index < num_blocks; ++block_index)
    rv_matrix_index += seq_len*std::max(1, (int)(rv_state.hap_->get_seq(block_index).size()));
  rv_matrix_index -= seq_len; // Index for match with next haplotype character

  int hap_index = 0;
  for (int block_index = 0; block_index < num_blocks-1; ++block_index){
    const std::string& block_seq = fw_state.hap_->get_seq(block_index);
    const int block_len          = block_seq.size();

    // Merge Fw and Rev alignments
    // Base j is aligned with the last character in BLOCK_INDEX, while base j+1 is aligned with the first character in BLOCK_INDEX+1
    // NOTE: Special cases of j=-1 (Rv only) and j=seq_len_-1 (Fw only) are implicitly handled in Case I
    hap_index       += std::max(0, block_len-1);
    fw_matrix_index += std::max(0, block_len-1)*seq_len;
    rv_matrix_index -= std::max(0, block_len-1)*seq_len;
    for (int j = 0; j < seq_len-1; ++j){
      double LL = LOG_MATCH_TO_MATCH[1] + fw_state.match_matrix_[fw_matrix_index + j] + rv_state.match_matrix_[rv_matrix_index - (j+1)];
      if (LL > max_LL){
	max_LL  = LL;
	max_dir = 2;
	fw_state.hap_index_    = hap_index;
	fw_state.seq_index_    = j;
	fw_state.matrix_index_ = fw_matrix_index + j;
	fw_state.matrix_type_  = AlignmentState::MATCH;

	// Save here as we know Case II is optimal and the last best alignment will set these appropriately
	// Other Rv state characteristics will be set later to avoid extra assignments
	rv_state.matrix_index_ = rv_matrix_index - (j+1);
      }
    }
    fw_matrix_index += seq_len;
    rv_matrix_index -= seq_len;
    hap_index       += (block_len > 0 ? 1 : 0);
  }

  // Optimal alignment is part Fw and part Rv from Case II. Set the remaining reverse alignment state based on the fw alignment parameters
  if (max_dir == 2){ 
    rv_state.hap_index_   = hap_size - fw_state.hap_index_ - 2;
    rv_state.seq_index_   = seq_len  - fw_state.seq_index_ - 2;
    rv_state.matrix_type_ = AlignmentState::MATCH;
    // NOTE: Matrix index was already saved in the for-loop above 
  }

  // Ensure we've found an alignment whose raw likelihood is <= 1
  assert(max_LL < TOLERANCE);
  assert(max_LL > AlignmentState::IMPOSSIBLE + TOLERANCE);

  return max_LL;
}
