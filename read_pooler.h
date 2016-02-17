#ifndef READ_POOLER_H_
#define READ_POOLER_H_

#include <assert.h>
#include <map>
#include <string>
#include <vector>

#include "bamtools/include/api/BamAlignment.h"

#include "base_quality.h"
#include "error.h"
#include "SeqAlignment/AlignmentData.h"

class ReadPooler {
 private:
  std::vector<Alignment> pooled_alns_;
  std::vector< std::vector<const std::string*> > qualities_by_pool_;
  std::map<std::string, int32_t> seq_to_pool_;
  bool pooled_;         // True iff pool() function has been invoked
  int32_t pool_index_;
  
 public:
  ReadPooler(){
    pool_index_ = 0;
    pooled_     = false;
  }

  ~ReadPooler(){
    for (unsigned int i = 0; i < qualities_by_pool_.size(); i++)
      for (unsigned int j = 0; j < qualities_by_pool_[i].size(); j++)
	delete qualities_by_pool_[i][j];
    qualities_by_pool_.clear();
  }

  int32_t num_pools(){ return pool_index_; }

  int32_t add_alignment(Alignment& aln);

  void pool(BaseQuality& base_quality){
    // For each pooled set of reads, set the base quality at each position to be the median across the set
    assert(pooled_alns_.size() == qualities_by_pool_.size());
    for (unsigned int i = 0; i < pooled_alns_.size(); i++)
      pooled_alns_[i].set_base_qualities(base_quality.median_base_qualities(qualities_by_pool_[i]));
    pooled_ = true;
  }

  std::vector<Alignment>& get_alignments(){
    return pooled_alns_;
  }
};

#endif
