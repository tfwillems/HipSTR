#include "read_pooler.h"

int32_t ReadPooler::add_alignment(Alignment& aln){
  if (pooled_)
    printErrorAndDie("Cannot call add_alignment function once pool() function has been invoked");
  
  auto pool_iter = seq_to_pool_.find(aln.get_sequence());
  if (pool_iter == seq_to_pool_.end()){
    seq_to_pool_[aln.get_sequence()] = pool_index_;
    pooled_alns_.push_back(Alignment(aln.get_start(), aln.get_stop(), "READPOOL", "", aln.get_sequence(), aln.get_alignment()));
    pooled_alns_.back().set_cigar_list(aln.get_cigar_list());
    qualities_by_pool_.push_back(std::vector<const std::string*>());
    qualities_by_pool_.back().push_back(new std::string(aln.get_base_qualities()));
    return pool_index_++;
  }
  else{
    qualities_by_pool_[pool_iter->second].push_back(new std::string(aln.get_base_qualities()));
    return pool_iter->second;
  }  
}
