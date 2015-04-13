#ifndef STUTTER_ALIGNER_H_
#define STUTTER_ALIGNER_H_

#include <string>


/* 
   Returns the total log-likelihood of the sequence of bases given the block sequence
   and that no indels occur within the block. 
   Requires that BASE_SEQ_LEN <= BLOCK_LEN
*/
double align_no_artifact(int block_len,                const std::string& block_seq,
			 int base_seq_len,             const char*   base_seq,
			 const double* base_log_wrong, const double* base_log_correct);

/*
  Returns the total log-likelihood of the sequence of bases given the block sequence
  and that a single insertion of size D(> 0) occurred within the block.
  Assumes that the insertion is equally likely to occur anywhere in the block (including along ends).
  Requires that BASE_SEQ_LEN <= BLOCK_LEN + D
*/
double align_pcr_insertion(int block_len,                const std::string& block_seq,
			   int base_seq_len,             const char*   base_seq,
			   const double* base_log_wrong, const double* base_log_correct,
			   int D);

/*
  Returns the total log-likelihood of the sequence of bases given the block sequence
  and that a single deletion of size D(< 0) occurred within the block.
  Assumes that the deletion is equally likely to occur anywhere within the block sequence.
  Requires that BASE_SEQ_LEN <= BLOCK_LEN + D
*/
double align_pcr_deletion(int block_len,                const std::string& block_seq,
			  int base_seq_len,             const char*   base_seq,
			  const double* base_log_wrong, const double* base_log_correct,
			  int D);

double align_stutter_region(int block_len,                const std::string& block_seq,
			    int base_seq_len,             const char* base_seq,
			    const double* base_log_wrong, const double* base_log_correct,
			    int D);
#endif

