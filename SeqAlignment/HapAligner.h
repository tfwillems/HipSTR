#ifndef HAP_ALIGNER_H_
#define HAP_ALIGNER_H_

#include <assert.h>
#include <string>
#include <vector>

#include "AlignmentData.h"
#include "AlignmentTraceback.h"
#include "../base_quality.h"
#include "Haplotype.h"

class HapAligner {
 private:
  Haplotype* fw_haplotype_;
  Haplotype* rev_haplotype_;
  std::vector<bool> realign_to_hap_;
  std::vector<HapBlock*> rev_blocks_;

  std::vector<int32_t> repeat_starts_;
  std::vector<int32_t> repeat_ends_;

  /**
   * Align the sequence contained in SEQ_0 -> SEQ_N using the recursion
   * 0 -> 1 -> 2 ... N
   **/
  void align_seq_to_hap(Haplotype* haplotype, bool reuse_alns,
			const char* seq_0, int seq_len,
			const double* base_log_wrong, const double* base_log_correct,
			double* match_matrix, double* insert_matrix, double* deletion_matrix,
			int* best_artifact_size, int* best_artifact_pos, double& left_prob);

  /**
   * Compute the log-probability of the alignment given the alignment matrices for the left and right segments.
   * Stores the index of the haplotype position with which the seed base is aligned in the maximum likelihood alignment
   **/
  double compute_aln_logprob(int base_seq_len, int seed_base,
			     char seed_char, double log_seed_wrong, double log_seed_correct,
			     double* l_match_matrix, double* l_insert_matrix, double* l_deletion_matrix, double l_prob,
			     double* r_match_matrix, double* r_insert_matrix, double* r_deletion_matrix, double r_prob,
			     int& max_index);

  std::string retrace(Haplotype* haplotype, const char* read_seq, const double* base_log_correct,
		      int seq_len, int block_index, int base_index, int matrix_index, double* l_match_matrix,
		      double* l_insert_matrix, double* l_deletion_matrix, int* best_artifact_size, int* best_artifact_pos,
		      AlignmentTrace& trace);

  void calc_best_seed_position(int32_t region_start, int32_t region_end,
			       int32_t& best_dist, int32_t& best_pos);

 public:
  HapAligner(Haplotype* haplotype, std::vector<bool>& realign_to_haplotype){
    assert(realign_to_haplotype.size() == haplotype->num_combs());
    fw_haplotype_   = haplotype;
    rev_haplotype_  = haplotype->reverse(rev_blocks_);
    realign_to_hap_ = realign_to_haplotype;


    for (int i = 0; i < fw_haplotype_->num_blocks(); i++){
      HapBlock* block = fw_haplotype_->get_block(i);
      if (block->get_repeat_info() != NULL){
	repeat_starts_.push_back(block->start());
	repeat_ends_.push_back(block->end());
      }
    }
  }

  ~HapAligner(){
    for (unsigned int i = 0; i < rev_blocks_.size(); i++)
      delete rev_blocks_[i];
    rev_blocks_.clear();
    delete rev_haplotype_;
  }

  /** 
   * Returns the 0-based index into the sequence string that should be used as the seed for alignment or -1 if no valid seed exists
   **/
  int calc_seed_base(Alignment& alignment);

  void process_read(Alignment& aln, int seed_base, BaseQuality* base_quality, bool retrace_aln,
		    double* prob_ptr, AlignmentTrace& traced_aln);

  void process_reads(std::vector<Alignment>& alignments, int init_read_index, BaseQuality* base_quality, std::vector<bool>& realign_read,
		     double* aln_probs, int* seed_positions);

  /*
    Retraces the Alignment's optimal alignment to the provided haplotype.
    Returns the result as a new Alignment relative to the reference haplotype
   */
  AlignmentTrace* trace_optimal_aln(Alignment& orig_aln, int seed_base, int best_haplotype, BaseQuality* base_quality);
};

#endif
