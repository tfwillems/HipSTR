#include <assert.h>
#include <math.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include "../error.h"
#include "../mathops.h"
#include "StutterAlignerClass.h"

void StutterAlignerClass::load_read(const int base_seq_len,       const char* base_seq,
				    const double* base_log_wrong, const double* base_log_correct){
  delete [] ins_probs_;
  delete [] del_probs_;
  delete [] match_probs_;
  ins_probs_   = new double[base_seq_len*num_insertions_];
  match_probs_ = new double[base_seq_len];
  if (num_deletions_ != 0)
    del_probs_ = new double[base_seq_len*num_deletions_];
  else
    del_probs_ = NULL;

  int ins_index = 0, del_index = 0, match_index = 0;
  for (int i = 0; i < base_seq_len; i++){
    int j;
    double log_prob = 0.0;
    for (j = 0; j < std::min(base_seq_len-i, -max_deletion_); j++){
      log_prob += (base_seq[-i-j] == block_seq_[-j] ? base_log_correct[-i-j] : base_log_wrong[-i-j]);
      if ((j+1) % period_ == 0)
	del_probs_[del_index++] = log_prob;
    }
    for (; j < -max_deletion_; j++)
      if ((j+1) % period_ == 0)
	del_index++;
    for (; j < std::min(base_seq_len-i, block_len_); j++)
      log_prob += (base_seq[-i-j] == block_seq_[-j] ? base_log_correct[-i-j] : base_log_wrong[-i-j]);
    match_probs_[match_index++] = log_prob;

    double log_ins_prob = 0.0;
    for (j = 0; j < std::min(max_insertion_, base_seq_len-i); j++){
      if (j % period_ < block_len_)
	log_ins_prob += (base_seq[-i-j] == block_seq_[-(j%period_)] ? base_log_correct[-i-j] : base_log_wrong[-i-j]);
      else
	log_ins_prob += base_log_correct[-i-j]; // No base to match to, so assume observed without error
                                                // We could potentially match to upstream flank?
      if ((j+1) % period_ == 0)
	ins_probs_[ins_index++] = log_ins_prob;
    }
    for (; j < max_insertion_; j++)
      if ((j+1) % period_ == 0)
	ins_probs_[ins_index++] = log_ins_prob;
  }
}

double StutterAlignerClass::align_no_artifact_reverse(const int offset){
  return match_probs_[offset];
}

double StutterAlignerClass::align_pcr_insertion_reverse(const int base_seq_len,       const char*   base_seq, const int offset,
							const double* base_log_wrong, const double* base_log_correct, const int D,
							int& best_ins_pos){
  assert(D > 0 && base_seq_len <= block_len_+D && D%period_ == 0);
  log_probs_.clear();
  double log_prior      = -int_log(block_len_+1);
  int* upstream_matches = upstream_match_lengths_[0] + block_len_ - 1;

  // Compute probability for i = 0
  double log_prob = log_prior + ins_probs_[num_insertions_*offset + D/period_ - 1] + (base_seq_len > D ? match_probs_[offset+D] : 0);
  best_ins_pos    = 0;
  double best_LL  = log_prob;
  log_probs_.push_back(log_prob);

  // Compute for all other i's, reusing previous result to accelerate computation
  int i = 0;
  for (; i > -std::min(std::max(0, base_seq_len-D), block_len_); i--){
    if (-i+period_ < block_len_) {
      if (upstream_matches[i] == 0){
        for (int index = i-period_; index >= i-D; index -= period_){
          log_prob -= (base_seq[index] == block_seq_[i]         ? base_log_correct[index] : base_log_wrong[index]);
          log_prob += (base_seq[index] == block_seq_[i-period_] ? base_log_correct[index] : base_log_wrong[index]);
        }
	log_probs_.push_back(log_prob);
      }
      else {
	log_probs_.push_back(int_log(upstream_matches[i])+log_prob);
	i -= (upstream_matches[i]-1);
      }
    }
    else
      log_probs_.push_back(log_prob);
    
    if (log_prob > best_LL || (left_align_ && (log_prob == best_LL))){
      best_ins_pos = 1-i;
      best_LL      = log_prob;
    }
  }  

  // Remaining configurations all have same likelihood so count all of them
  if (i > -block_len_)
    log_probs_.push_back(int_log(block_len_+i)+log_prob);

  // Convert to raw probabilities, add, take the log while avoiding underflow
  return fast_log_sum_exp(log_probs_, std::max(best_LL, log_probs_.back()));
}

double StutterAlignerClass::align_pcr_deletion_reverse(const int base_seq_len,       const char*   base_seq, const int offset,
						       const double* base_log_wrong, const double* base_log_correct, const int D,
						       int& best_del_pos){
  assert(D < 0 && block_len_+D >= 0 && base_seq_len <= block_len_+D);
  log_probs_.clear();
  int* upstream_matches = upstream_match_lengths_[-D/period_ - 1] + block_len_ - 1;
  double log_prior = -int_log(block_len_+D+1);
  double log_prob  = log_prior;

  // Compute probability for i = 0
  if (offset+D >= 0)
    log_prob += match_probs_[offset+D] - del_probs_[(offset+D)*num_deletions_ - D/period_ - 1];
  else
    for (int j = 0; j > -base_seq_len; j--)
      log_prob += (block_seq_[j+D] == base_seq[j] ? base_log_correct[j] : base_log_wrong[j]);
  best_del_pos   = 0;
  double best_LL = log_prob;
  log_probs_.push_back(log_prob);

  // Compute for all other i's, reusing previous result to accelerate computation
  int i;
  for (i = 0; i > -base_seq_len; i--){
    if (upstream_matches[i] == 0){
      log_prob    -= (block_seq_[i+D] == base_seq[i] ? base_log_correct[i] : base_log_wrong[i]);
      log_prob    += (block_seq_[i]   == base_seq[i] ? base_log_correct[i] : base_log_wrong[i]);
      log_probs_.push_back(log_prob);
    }
    else {
      log_probs_.push_back(int_log(upstream_matches[i])+log_prob);
      i -= (upstream_matches[i]-1);
    }
    
    if (log_prob > best_LL || (left_align_ && (log_prob == best_LL))){
      best_del_pos = 1-i;
      best_LL      = log_prob;
    }
  }

  // Remaining configurations all have same likelihood so count all of them
  if(-i < block_len_+D)
    log_probs_.push_back(int_log(block_len_+D+i)+log_prob);

  // Convert to raw probabilities, add, take the log while avoiding underflow
  return fast_log_sum_exp(log_probs_, std::max(best_LL, log_probs_.back()));
}

double StutterAlignerClass::align_stutter_region_reverse(const int base_seq_len,       const char*   base_seq, const int offset,
							 const double* base_log_wrong, const double* base_log_correct, const int D,
							 int& best_pos){
  best_pos = -1;
  if (D == 0)
    return align_no_artifact_reverse(offset);
  else if (D > 0)
    return align_pcr_insertion_reverse(base_seq_len, base_seq, offset, base_log_wrong, base_log_correct, D, best_pos);
  else
    return align_pcr_deletion_reverse(base_seq_len, base_seq,  offset, base_log_wrong, base_log_correct, D, best_pos);
}
