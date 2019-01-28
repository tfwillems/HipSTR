#ifndef HAP_ALIGNER_H_
#define HAP_ALIGNER_H_

#include <assert.h>
#include <vector>

#include "AlignmentData.h"
#include "AlignmentMatrixCache.h"
#include "AlignmentTraceback.h"
#include "../base_quality.h"
#include "../error.h"
#include "HapBlock.h"
#include "Haplotype.h"

class HapAligner {
 private:
  Haplotype *fw_haplotype_, *rv_haplotype_;
  std::vector<bool> realign_to_hap_;
  std::vector<HapBlock*> rv_blocks_;
  std::vector<int32_t> repeat_starts_, repeat_ends_;

  void calc_best_seed_position(int32_t region_start, int32_t region_end,
			       int32_t& best_dist, int32_t& best_pos);

  // Returns the 0-based index into the sequence string that should be used as the seed for alignment or -1 if no valid seed exists
  int calc_seed_base(const Alignment& alignment);

  void process_read(const Alignment& aln, int seed_base, const BaseQuality* base_quality, bool retrace_aln,
		    double* prob_ptr, AlignmentTrace& traced_aln, bool debug);

  // Private unimplemented copy constructor and assignment operator to prevent operations
  HapAligner(const HapAligner& other);
  HapAligner& operator=(const HapAligner& other);

 public:
  HapAligner(Haplotype* haplotype, std::vector<bool>& realign_to_haplotype){
    assert(realign_to_haplotype.size() == haplotype->num_combs());
    fw_haplotype_   = haplotype;
    rv_haplotype_   = haplotype->reverse(rv_blocks_);
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
    for (unsigned int i = 0; i < rv_blocks_.size(); i++)
      delete rv_blocks_[i];
    rv_blocks_.clear();
    delete rv_haplotype_;
  }

  void process_reads(const std::vector<Alignment>& alignments,
		     const BaseQuality* base_quality, const std::vector<bool>& realign_read,
		     double* aln_probs, int* seed_positions);

  /*
    Retraces the Alignment's optimal alignment to the provided haplotype.
    Returns the result as a new Alignment relative to the reference haplotype
   */
  AlignmentTrace* trace_optimal_aln(const Alignment& orig_aln, int seed_base, int best_haplotype, const BaseQuality* base_quality,
				    bool debug);
};

#endif
