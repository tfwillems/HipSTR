#include <assert.h>
#include <math.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include "../error.h"
#include "../mathops.h"
#include "StutterAlignerClass.h"


double StutterAlignerClass::align_no_artifact_reverse(const int base_seq_len,       const char*   base_seq,
						      const double* base_log_wrong, const double* base_log_correct){
  double log_prob = 0.0;
  for (int i = 0; i < base_seq_len; i++)
    log_prob += (block_seq_[-i] == base_seq[-i] ? base_log_correct[-i] : base_log_wrong[-i]);
  return log_prob;
}


double StutterAlignerClass::align_pcr_insertion_reverse(const int base_seq_len,       const char*   base_seq,
							const double* base_log_wrong, const double* base_log_correct, const int D,
							int& best_ins_pos){
  assert(D > 0 && base_seq_len <= block_len_+D && D%period_ == 0);
  log_probs_.clear();
  double log_prior = -int_log(block_len_+1);
  int* upstream_matches = upstream_match_lengths_[0] + block_len_ - 1;

  // Compute probability for i = 0
  double log_prob = log_prior;
  // Bases matched with insertion. Calculate emission probability according to agreement with proximal haplotype sequence
  for (int j = 0; j < std::min(D, base_seq_len); j++)
    log_prob += (base_seq[-j] == block_seq_[-(j%period_)] ? base_log_correct[-j] : base_log_wrong[-j]);
  for (int j = D; j < base_seq_len; j++)
    log_prob += (block_seq_[-j+D] == base_seq[-j] ? base_log_correct[-j] : base_log_wrong[-j]); // Bases matched with block characters
  log_probs_.push_back(log_prob);
  best_ins_pos = 0;
  double best_LL = log_prob;

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
  return fast_log_sum_exp(log_probs_);
}


double StutterAlignerClass::align_pcr_deletion_reverse(const int base_seq_len,       const char*   base_seq,
						       const double* base_log_wrong, const double* base_log_correct, const int D,
						       int& best_del_pos){
  assert(D < 0 && block_len_+D >= 0 && base_seq_len <= block_len_+D);
  log_probs_.clear();
  double log_prior = -int_log(block_len_+D+1);
  int* upstream_matches = upstream_match_lengths_[-D/period_ - 1] + block_len_ - 1;
  
  // Compute probability for i = 0
  double log_prob = log_prior;
  for (int j = 0; j > -base_seq_len; j--)
    log_prob += (block_seq_[j+D] == base_seq[j] ? base_log_correct[j] : base_log_wrong[j]);
  log_probs_.push_back(log_prob);
  best_del_pos = 0;
  double best_LL = log_prob;

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
  //if (base_seq_len < block_len_+D){
  if(-i < block_len_+D)
    log_probs_.push_back(int_log(block_len_+D+i)+log_prob);

  // Convert to raw probabilities, add, take the log while avoiding underflow
  return fast_log_sum_exp(log_probs_);
}

double StutterAlignerClass::align_stutter_region_reverse(const int base_seq_len,       const char*   base_seq,
							 const double* base_log_wrong, const double* base_log_correct, const int D,
							 int& best_pos){
  best_pos = -1;
  if (D == 0)
    return align_no_artifact_reverse(base_seq_len, base_seq, base_log_wrong, base_log_correct);
  else if (D > 0)
    return align_pcr_insertion_reverse(base_seq_len, base_seq, base_log_wrong, base_log_correct, D, best_pos);
  else
    return align_pcr_deletion_reverse(base_seq_len, base_seq, base_log_wrong, base_log_correct, D, best_pos);
}
