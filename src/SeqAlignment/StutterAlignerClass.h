#ifndef STUTTER_ALIGNER_CLASS_H_
#define STUTTER_ALIGNER_CLASS_H_

#include <string>
#include <vector>

#include "RepeatStutterInfo.h"

class StutterAlignerClass {
 private:
  std::vector<double> log_probs_;
  char* block_seq_;
  const int  block_len_;
  const int  period_;
  const bool left_align_;
  std::vector<int*> upstream_match_lengths_;

  int num_artifacts_;
  double* ins_probs_;
  double* del_probs_;
  double* match_probs_;
 
  double align_no_artifact_reverse(const int base_seq_len,       const char*   base_seq, const int offset,
				   const double* base_log_wrong, const double* base_log_correct);
  
  double align_pcr_insertion_reverse(const int base_seq_len,       const char*   base_seq, const int offset,
				     const double* base_log_wrong, const double* base_log_correct, const int D,
				     int& best_ins_pos);
  
  double align_pcr_deletion_reverse(const int base_seq_len,       const char*   base_seq, const int offset,
				    const double* base_log_wrong, const double* base_log_correct, const int D,
				    int& best_del_pos);

  int* num_upstream_matches(const std::string& seq, int period) const {
    int* match_lengths = new int[seq.size()];
    for (unsigned int i = 0; i < std::min(period, (int)seq.size()); i++)
      match_lengths[i] = 0;
    for (unsigned int i = period; i < seq.size(); i++)
      match_lengths[i] = (seq[i-period] != seq[i] ? 0 : 1 + match_lengths[i-1]);
    return match_lengths;
  }
  
 public:
 StutterAlignerClass(const std::string& block_seq, int period, bool left_align, const RepeatStutterInfo* stutter_info)
   : block_len_(block_seq.size()), period_(period), left_align_(left_align){
    log_probs_.reserve(block_len_+1);

    // Create and fill the array and make the pointer refer to the last character
    block_seq_ = new char[block_len_];
    for (int i = 0; i < block_len_; i++)
      block_seq_[i] = block_seq[i];
    block_seq_ += (block_len_-1);
    
    // Calculate periodicity info for relevant offsets
    for (int i = -period; i >= stutter_info->max_deletion(); i -= period)
      upstream_match_lengths_.push_back(num_upstream_matches(block_seq, -i));

    assert(stutter_info->max_insertion() == -1*stutter_info->max_deletion());
    assert(stutter_info->max_insertion()%period_ == 0 && block_len_+stutter_info->max_deletion() >= 0);
    num_artifacts_ = stutter_info->max_insertion()/period_;
    ins_probs_     = NULL;
    del_probs_     = NULL;
    match_probs_   = NULL;
  }

  ~StutterAlignerClass(){
    // Point back to front of array and delete allocated space
    block_seq_ -= (block_len_-1);
    delete [] block_seq_;

    for (unsigned int i = 0; i < upstream_match_lengths_.size(); i++)
      delete [] upstream_match_lengths_[i];
    upstream_match_lengths_.clear();

    delete [] ins_probs_;
    delete [] del_probs_;
    delete [] match_probs_;
  }
  
  int block_len()         const { return block_len_;  }
  int period()            const { return period_;     }
  bool left_align()       const { return left_align_; }
  const char* block_seq() const { return block_seq_;  }


  void load_read(const int base_seq_len,       const char* base_seq,
		 const double* base_log_wrong, const double* base_log_correct,
		 int max_deletion,             int max_insertion);
  
  /* Returns the total log-likelihood of the base sequence given the block sequence and the associated quality scores.
   * Assumes that an artifact of size D occurs with equal probability throughout the block sequence.
   * Progresses BACKWARDS along the sequences using the provided pointers, so the pointers refer to the rightmost
   * bases in both the base sequence and block sequence
   */
  double align_stutter_region_reverse(const int base_seq_len,             const char*   base_seq,   const int offset,
				      const double* base_log_wrong, const double* base_log_correct, const int D,
				      int& best_pos);
};

#endif
