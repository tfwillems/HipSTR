#ifndef STUTTER_ALIGNER_CLASS_H_
#define STUTTER_ALIGNER_CLASS_H_

#include <string>
#include <vector>

#include "RepeatStutterInfo.h"


class StutterAlignerClass {
 private:
  std::vector<double> log_probs_;
  const int   block_len_;
  char* block_seq_;
  const int   period_;
  const bool  left_align_;
  std::vector<int*> upstream_match_lengths_;
 
  double align_no_artifact_reverse(const int base_seq_len,       const char*   base_seq,
				   const double* base_log_wrong, const double* base_log_correct);
  
  double align_pcr_insertion_reverse(const int base_seq_len,       const char*   base_seq,
				     const double* base_log_wrong, const double* base_log_correct, const int D,
				     int& best_ins_pos);
  
  double align_pcr_deletion_reverse(const int base_seq_len,       const char*   base_seq,
				    const double* base_log_wrong, const double* base_log_correct, const int D,
				    int& best_del_pos);

  int* num_upstream_matches(std::string& seq, int period){
    int* match_lengths = new int[seq.size()];
    for (unsigned int i = 0; i < period; i++)
      match_lengths[i] = 0;
    for (unsigned int i = period; i < seq.size(); i++)
      match_lengths[i] = (seq[i-period] != seq[i] ? 0 : 1 + match_lengths[i-1]);
    return match_lengths;
  }
  
 public:
 StutterAlignerClass(std::string& block_seq, int period, bool left_align, RepeatStutterInfo* stutter_info)
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
  }

  ~StutterAlignerClass(){
    // Point back to front of array and delete allocated space
    block_seq_ -= (block_len_-1);
    delete [] block_seq_;

    for (unsigned int i = 0; i < upstream_match_lengths_.size(); i++)
      delete [] upstream_match_lengths_[i];
    upstream_match_lengths_.clear();
  }
  
  int block_len()         const { return block_len_;  }
  int period()            const { return period_;     }
  bool left_align()       const { return left_align_; }
  const char* block_seq() const { return block_seq_;  }
  
  /* Returns the total log-likelihood of the base sequence given the block sequence and the associated quality scores.
   * Assumes that an artifact of size D occurs with equal probability throughout the block sequence.
   * Progresses BACKWARDS along the sequences using the provided pointers, so the pointers refer to the rightmost
   * bases in both the base sequence and block sequence
   */
  double align_stutter_region_reverse(const int base_seq_len,             const char*   base_seq,
				      const double* base_log_wrong, const double* base_log_correct, const int D,
				      int& best_pos);
};

#endif
