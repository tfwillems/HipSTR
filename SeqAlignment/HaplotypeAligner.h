#ifndef HAPLOTYPEALIGNER_H_
#define HAPLOTYPEALIGNER_H_

#include <string>
#include <vector>

#include "AlignmentData.h"
#include "BaseQuality.h"
#include "Haplotype.h"

class HaplotypeAligner {
 private:
  Haplotype* haplotype_;

 public:
  HaplotypeAligner(Haplotype* haplotype){
    haplotype_ = haplotype;
  }

  double align(const char* base_seq, 
	       int base_seq_len, 
	       int seed_base,
	       const double* base_log_wrong, 
	       const double* base_log_correct,
	       double mapping_quality,
	       double* l_match_matrix, 
	       double* l_insert_matrix,
	       double* r_match_matrix, 
	       double* r_insert_matrix);


  void align(std::vector<Alignment>& alignments, 
	     BaseQuality& base_quality);

  /** 
   * Returns the 0-based index into the sequence string that should be 
   * used as the seed for alignment. Returns -1 iff no valid seed exists 
   **/
  int calc_seed_base(Alignment& alignment);

  /** 
   *  Align the sequence contained in SEQ_0 -> SEQ_N using the recursion
   *  0 -> 1 -> 2 ... N
  **/
  void align_left_flank(const char* seq_0, 
			int seq_len, 
			const double* base_log_wrong, 
			const double* base_log_correct, 
			double* match_matrix, 
			double* insert_matrix,
			double& left_prob);
  
  /** 
   * Align the sequence contained in SEQ_N -> SEQ_END using the recursion
   * END -> END-1 -> END-2 ... N
   **/
  void align_right_flank(const char* seq_n, 
			 int seq_len, 
			 const double* base_log_wrong, 
			 const double* base_log_correct,
			 double* match_matrix, 
			 double* insert_matrix,
			 double& right_prob);
};

#endif
