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

  ~ReadPooler(){ }


  int32_t num_pools(){ return pool_index_; }

  int32_t add_alignment(Alignment& aln){
    if (pooled_)
      printErrorAndDie("Cannot call add_alignment function once pool() function has been invoked.");

    auto pool_iter = seq_to_pool_.find(aln.get_sequence());
    if (pool_iter == seq_to_pool_.end()){
      seq_to_pool_[aln.get_sequence()] = pool_index_;
      pooled_alns_.push_back(Alignment(aln.get_start(), aln.get_stop(), "POOL", "", aln.get_sequence(), aln.get_alignment()));
      pooled_alns_.back().set_cigar_list(aln.get_cigar_list());
      qualities_by_pool_.push_back(std::vector<const std::string*>());
      qualities_by_pool_.back().push_back(&(aln.get_base_qualities()));
      return pool_index_++;
    }
    else{
      qualities_by_pool_[pool_iter->second].push_back(&(aln.get_base_qualities()));
      return pool_iter->second;
    }  
  }

  void pool(BaseQuality& base_quality){
    // For each pooled set of reads, set the base quality at each position to be the median across the set
    assert(pooled_alns_.size() == qualities_by_pool_.size());
    for (unsigned int i = 0; i < pooled_alns_.size(); i++)
      pooled_alns_[i].set_base_qualities(base_quality.median_base_qualities(qualities_by_pool_[i]));
    qualities_by_pool_.clear();
    pooled_ = true;
  }

  std::vector<Alignment>& get_alignments(){
    return pooled_alns_;
  }

  Alignment& get_alignment(int32_t pool_index){
    assert(pool_index < pooled_alns_.size() && pool_index >= 0);
    if(!pooled_)
      printErrorAndDie("Cannot retrive pooled alignments before invoking pool() function");
    return pooled_alns_[pool_index];
  }
};

#endif
