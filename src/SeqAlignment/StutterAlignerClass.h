#ifndef STUTTER_ALIGNER_CLASS_H_
#define STUTTER_ALIGNER_CLASS_H_

#include <string>
#include <vector>

#include "RepeatStutterInfo.h"

class StutterAlignerClass {
 private:
  char* block_seq_;
  int*  mods_;
  const int  block_len_;
  const int  period_;
  const bool left_align_;
  std::vector<int*> upstream_match_lengths_;

  int num_insertions_, num_deletions_;
  int max_insertion_,  max_deletion_;

  double* ins_probs_;
  double* del_probs_;
  double* match_probs_;
 
  double align_no_artifact_reverse(const int offset);
  
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

  // Private unimplemented copy constructor to prevent operation
  StutterAlignerClass(const StutterAlignerClass& other);
  StutterAlignerClass& operator=(const StutterAlignerClass& other);
  
 public:
 StutterAlignerClass(const std::string& block_seq, int period, bool left_align, const RepeatStutterInfo* stutter_info)
   : block_len_(block_seq.size()), period_(period), left_align_(left_align){
    assert(stutter_info->max_insertion()%period_ == 0 && stutter_info->max_deletion()% period == 0);

    // Create and fill the array and make the pointer refer to the last character
    if (block_len_ == 0)
      block_seq_ = NULL;
    else {
      block_seq_ = new char[block_len_];
      for (int i = 0; i < block_len_; i++)
	block_seq_[i] = block_seq[i];
      block_seq_ += (block_len_-1);
    }

    num_insertions_ = stutter_info->max_insertion()/period_;
    num_deletions_  = -1*(stutter_info->max_deletion()/period_);
    while (num_deletions_*period > block_len_)
      num_deletions_--;
    max_insertion_ = period_*num_insertions_;
    max_deletion_  = -period_*num_deletions_;

    // Calculate periodicity info for relevant offsets
    for (int i = -period; i >= max_deletion_; i -= period)
      upstream_match_lengths_.push_back(num_upstream_matches(block_seq, -i));
    if (max_deletion_ == 0) // We require this for insertion calculations
      upstream_match_lengths_.push_back(block_seq.empty() ? 0 : num_upstream_matches(block_seq, period));

    ins_probs_   = NULL;
    del_probs_   = NULL;
    match_probs_ = NULL;

    // Precompute various modulus values to avoid expensive downstream operations
    mods_ = new int[500];
    for (int i = 0; i < 500; ++i)
      mods_[i] = (i % period_);
  }

  ~StutterAlignerClass(){
    // Point back to front of array and delete allocated space
    if (block_seq_ != NULL)
      block_seq_ -= (block_len_-1);
    delete [] block_seq_;

    for (unsigned int i = 0; i < upstream_match_lengths_.size(); i++)
      delete [] upstream_match_lengths_[i];
    upstream_match_lengths_.clear();

    delete [] ins_probs_;
    delete [] del_probs_;
    delete [] match_probs_;
    delete [] mods_;
  }

  void load_read(const int base_seq_len,       const char* base_seq,
		 const double* base_log_wrong, const double* base_log_correct);
  
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
