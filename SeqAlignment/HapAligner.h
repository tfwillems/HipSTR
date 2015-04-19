#ifndef HAP_ALIGNER_H_
#define HAP_ALIGNER_H_

#include <string>
#include <vector>

#include "AlignmentData.h"
#include "../base_quality.h"
#include "../region.h"
#include "Haplotype.h"

class HapAligner {
 private:
  Haplotype* haplotype_;
  int total_reads_;
  Region* region_;
  int32_t max_ref_flank_len_; // Maximum length of reference sequences flanking the repetitive region
  BaseQuality* base_quality_;
  double* log_align_probs_;   // Iterates over haplotype and then read

  /**
   * Align the sequence contained in SEQ_0 -> SEQ_N using the recursion
   * 0 -> 1 -> 2 ... N
   **/
  void align_left_flank(const char* seq_0, int seq_len,
			const double* base_log_wrong, const double* base_log_correct,
			double* match_matrix, double* insert_matrix, double& left_prob);

  /**                                                                                                                                                                         
   * Align the sequence contained in SEQ_N -> SEQ_END using the recursion
   * END -> END-1 -> END-2 ... N
   **/
  void align_right_flank(const char* seq_n, int seq_len,
			 const double* base_log_wrong, const double* base_log_correct,
                         double* match_matrix, double* insert_matrix, double& right_prob);

  /**
   * Compute the log-probability of the alignment given the 
   * alignment matrices for the left and right segments
   **/
  double compute_aln_logprob(int base_seq_len, int seed_base,
			     char seed_char, double log_seed_wrong, double log_seed_correct,
			     double* l_match_matrix, double* l_insert_matrix, double l_prob,
			     double* r_match_matrix, double* r_insert_matrix, double r_prob);

 public:
  HapAligner(Haplotype* haplotype, Region region, int32_t max_ref_flank_len, BaseQuality* base_quality, int total_reads){
    haplotype_         = haplotype;
    total_reads_       = total_reads;
    region_            = new Region(region.chrom(), region.start(), region.stop(), region.period());
    max_ref_flank_len_ = max_ref_flank_len;
    base_quality_      = base_quality;
    log_align_probs_   = new double[haplotype_->num_combs()*total_reads_];
  }

  ~HapAligner(){
    delete region_;
    delete [] log_align_probs_;
  }

  /** 
   * Returns the 0-based index into the sequence string that should be 
   * used as the seed for alignment. Returns -1 iff no valid seed exists 
   **/
  int calc_seed_base(Alignment& alignment);


  void process_reads(std::vector<Alignment>& alignments, int init_read_index, double* aln_probs);
};

#endif
