#ifndef STUTTER_ALIGNER_H_
#define STUTTER_ALIGNER_H_

#include <string>

double align_no_artifact_reverse(int block_len,                const char*   block_seq,
                                 int base_seq_len,             const char*   base_seq,
                                 const double* base_log_wrong, const double* base_log_correct);

double align_pcr_insertion_reverse(int block_len,                const char*   block_seq,
                                   int base_seq_len,             const char*   base_seq,
                                   const double* base_log_wrong, const double* base_log_correct,
                                   int D, int period,
				   int& best_ins_pos);

double align_pcr_deletion_reverse(int block_len,                const char*   block_seq,
                                  int base_seq_len,             const char*   base_seq,
                                  const double* base_log_wrong, const double* base_log_correct,
                                  int D,
				  int& best_del_pos);

/* Returns the total log-likelihood of the base sequence given the block sequence and the associated quality scores.
 * Assumes that an artifact of size D occurs with equal probability throughout the block sequence.
 * Progresses BACKWARDS along the sequences using the provided pointers, so the pointers refer to the rightmost
 * bases in both the base sequence and block sequence
 */
double align_stutter_region_reverse(int block_len,                const char*   block_seq,
				    int base_seq_len,             const char*   base_seq,
				    const double* base_log_wrong, const double* base_log_correct,
				    bool left_align, int D, int period,
				    int& best_pos);

#endif

