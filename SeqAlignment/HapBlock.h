#ifndef HAPBLOCK_H_
#define HAPBLOCK_H_

#include <assert.h>
#include <stdint.h>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <string>
#include <vector>

#include "../error.h"
#include "../stringops.h"
#include "RepeatStutterInfo.h"

class HapBlock {
 protected:
  std::string ref_seq_;
  std::vector<std::string> alt_seqs_;
  int32_t start_;  // Start of region (inclusive)
  int32_t end_;    // End   of region (not inclusive)
  int min_size_;
  int max_size_;

  std::vector<int*> l_homopolymer_lens_;
  std::vector<int*> r_homopolymer_lens_;
  std::vector<int>  suffix_matches_;

  void calc_homopolymer_lengths();

 public:
  HapBlock(int32_t start, int32_t end, std::string ref_seq) {
    start_    = start;
    end_      = end;
    ref_seq_  = ref_seq;
    min_size_ = (int)ref_seq.size();
    max_size_ = (int)ref_seq.size();
    alt_seqs_ = std::vector<std::string>();
    suffix_matches_ = std::vector<int>();
    suffix_matches_.push_back(0);
  }

  virtual ~HapBlock() {
    for (unsigned int i = 0; i < l_homopolymer_lens_.size(); i++)
      delete [] l_homopolymer_lens_[i];
    for (unsigned int i = 0; i < r_homopolymer_lens_.size(); i++)
      delete [] r_homopolymer_lens_[i];
    l_homopolymer_lens_.clear();
    r_homopolymer_lens_.clear();
  }

  virtual RepeatStutterInfo* get_repeat_info() {
    return NULL;
  }

  int32_t start()   const { return start_; }
  int32_t end()     const { return end_;   }
  int num_options() const { return 1 + alt_seqs_.size(); }

  virtual void add_alternate(std::string& alt) {
    alt_seqs_.push_back(alt);
    min_size_ = std::min(min_size_, (int)alt.size());
    max_size_ = std::max(max_size_, (int)alt.size());
    if (alt_seqs_.size() == 1)
      suffix_matches_.push_back(length_suffix_match(ref_seq_, alt));
    else
      suffix_matches_.push_back(length_suffix_match(alt_seqs_[alt_seqs_.size()-2], alt));
  }

  int min_size() const { return min_size_; }
  int max_size() const { return max_size_; }

  void print(std::ostream& out);

  int size(int index) const {
    if (index == 0)
      return ref_seq_.size();
    else if (index-1 < alt_seqs_.size())
      return alt_seqs_[index-1].size();
    else
      throw std::out_of_range("Index out of bounds in HapBlock::size()");
  }

  const std::string& get_seq(unsigned int index) {
    if (index == 0)
      return ref_seq_;
    else if (index-1 < alt_seqs_.size())
      return alt_seqs_[index-1];
    else
      throw std::out_of_range("Index out of bounds in HapBlock::get_seq()");
  }

  inline int suffix_match_len(unsigned int seq_index){
    assert(seq_index < suffix_matches_.size());
    return suffix_matches_[seq_index];
  }

  inline unsigned int left_homopolymer_len(unsigned int seq_index, int base_index) {
    assert(seq_index < l_homopolymer_lens_.size());
    return l_homopolymer_lens_[seq_index][base_index];
  }

  inline unsigned int right_homopolymer_len(unsigned int seq_index, int base_index) {
    assert(seq_index < r_homopolymer_lens_.size());
    return r_homopolymer_lens_[seq_index][base_index];
  }
  
  void initialize();

  virtual HapBlock* reverse(){
    std::string rev_ref_seq = ref_seq_;
    std::reverse(rev_ref_seq.begin(), rev_ref_seq.end());
    HapBlock* rev_block = new HapBlock(end_, start_, rev_ref_seq);
    for (unsigned int i = 0; i < alt_seqs_.size(); i++) {
      std::string alt = alt_seqs_[i];
      std::reverse(alt.begin(), alt.end());
      rev_block->add_alternate(alt);
    }
    rev_block->initialize();
    return rev_block;
  }
};

#endif
