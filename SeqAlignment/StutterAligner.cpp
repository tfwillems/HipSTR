#include <assert.h>
#include <math.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include "../error.h"
#include "../mathops.h"
#include "StutterAligner.h"

double align_no_artifact_reverse(int block_len,                const char*   block_seq,
				 int base_seq_len,             const char*   base_seq,
				 const double* base_log_wrong, const double* base_log_correct){
  double log_prob = 0.0;
  for (int i = 0; i < base_seq_len; i++)
    log_prob += (block_seq[-i] == base_seq[-i] ? base_log_correct[-i] : base_log_wrong[-i]);
  return log_prob;
}

double align_pcr_insertion_reverse(int block_len,                const char*   block_seq,
				   int base_seq_len,             const char*   base_seq,
				   const double* base_log_wrong, const double* base_log_correct,
				   bool left_align, int D, int period,
				   int& best_ins_pos){
  assert(D > 0 && base_seq_len <= block_len+D && D%period == 0);
  std::vector<double> log_probs; log_probs.reserve(block_len+1);
  double log_prior = -int_log(block_len+1);

  // Compute probability for i = 0
  double log_prob = log_prior;
  // Bases matched with insertion. Calculate emission probability according to agreement with proximal haplotype sequence
  for (int j = 0; j < std::min(D, base_seq_len); j++)
    log_prob += (base_seq[-j] == block_seq[-(j%period)] ? base_log_correct[-j] : base_log_wrong[-j]);
  for (int j = D; j < base_seq_len; j++)
    log_prob += (block_seq[-j+D] == base_seq[-j] ? base_log_correct[-j] : base_log_wrong[-j]); // Bases matched with block characters
  log_probs.push_back(log_prob);
  best_ins_pos = 0;

  // Compute for all other i's, reusing previous result to accelerate computation
  int i = 0;
  for (; i > -std::min(std::max(0, base_seq_len-D), block_len); i--){
    if (-i+period < block_len) {
      if (block_seq[i] != block_seq[i-period]){
        for (int index = i-period; index >= i-D; index -= period){
          log_prob -= (base_seq[index] == block_seq[i]        ? base_log_correct[index] : base_log_wrong[index]);
          log_prob += (base_seq[index] == block_seq[i-period] ? base_log_correct[index] : base_log_wrong[index]);
        }
      }
    }
    log_probs.push_back(log_prob); // Store LL for configuration
    if (log_prob > log_probs[best_ins_pos] || (left_align && (log_prob == log_probs[best_ins_pos])))
      best_ins_pos = log_probs.size()-1;
  }

  // Remaining configurations all have same likelihood so count all of them
  if (base_seq_len < block_len+D)
    log_probs.push_back(int_log(block_len+i)+log_prob);

  // Convert to raw probabilities, add, take the log while avoiding underflow
  return fast_log_sum_exp(log_probs);
}

double align_pcr_deletion_reverse(int block_len,                const char*   block_seq,
				  int base_seq_len,             const char*   base_seq,
				  const double* base_log_wrong, const double* base_log_correct, bool left_align, int D,
				  int& best_del_pos){
  assert(D < 0 && block_len+D >= 0 && base_seq_len <= block_len+D);
  std::vector<double> log_probs; log_probs.reserve(block_len+D+1);
  double log_prior = -int_log(block_len+D+1);

  // Compute probability for i = 0
  double log_prob = log_prior;
  for (int j = 0; j > -base_seq_len; j--)
    log_prob += (block_seq[j+D] == base_seq[j] ? base_log_correct[j] : base_log_wrong[j]);
  log_probs.push_back(log_prob);
  best_del_pos = 0;

  // Compute for all other i's, reusing previous result to accelerate computation
  for (int i = 0; i > -base_seq_len; i--){
    log_prob    -= (block_seq[i+D] == base_seq[i] ? base_log_correct[i] : base_log_wrong[i]);
    log_prob    += (block_seq[i]   == base_seq[i] ? base_log_correct[i] : base_log_wrong[i]);
    log_probs.push_back(log_prob);
    if (log_prob > log_probs[best_del_pos] || (left_align && (log_prob == log_probs[best_del_pos])))
      best_del_pos = log_probs.size()-1;
  }

  // Remaining configurations all have same likelihood so count all of them
  if (base_seq_len < block_len+D)
    log_probs.push_back(int_log(block_len+D-base_seq_len)+log_prob);

  // Convert to raw probabilities, add, take the log while avoiding underflow
  return fast_log_sum_exp(log_probs);
}

double align_stutter_region_reverse(int block_len,                const char*   block_seq,
				    int base_seq_len,             const char*   base_seq,
				    const double* base_log_wrong, const double* base_log_correct,
				    bool left_align, int D, int period,
				    int& best_pos){
  best_pos = -1;
  if (D == 0)
    return align_no_artifact_reverse(block_len, block_seq, base_seq_len, base_seq, base_log_wrong, base_log_correct);
  else if (D > 0)
    return align_pcr_insertion_reverse(block_len, block_seq, base_seq_len, base_seq, base_log_wrong, base_log_correct, left_align, D, period, best_pos);
  else
    return align_pcr_deletion_reverse(block_len, block_seq, base_seq_len, base_seq, base_log_wrong, base_log_correct, left_align, D, best_pos);
}
