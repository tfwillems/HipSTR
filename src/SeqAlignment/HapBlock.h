#ifndef HAPBLOCK_H_
#define HAPBLOCK_H_

#include <assert.h>
#include <stdint.h>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "../stringops.h"
#include "RepeatStutterInfo.h"
#include "StutterAlignerClass.h"

class HapBlock {
 protected:
  const std::string ref_seq_;
  std::vector<std::string> alt_seqs_;
  std::set<std::string> seq_set_;
  int32_t start_; // Start of region (inclusive)
  int32_t end_;   // End of region (not inclusive)
  int min_size_, max_size_;

  std::vector<int*> l_homopolymer_lens_;
  std::vector<int*> r_homopolymer_lens_;
  std::vector<int>  suffix_matches_;

  // Compute the homopolymer lengths and store them in the resulting vectors
  void calc_homopolymer_lengths(const std::string& seq, std::vector<int*>& llen_vec, std::vector<int*>& rlen_vec);

  // Private unimplemented copy constructor and assignment operator to prevent operations
  HapBlock(const HapBlock& other);
  HapBlock& operator=(const HapBlock& other);

 public:
  HapBlock(int32_t start, int32_t end, const std::string& ref_seq)
    : ref_seq_(ref_seq){
    start_    = start;
    end_      = end;
    min_size_ = (int)ref_seq.size();
    max_size_ = (int)ref_seq.size();
    suffix_matches_.push_back(0);
    seq_set_.insert(ref_seq);
    calc_homopolymer_lengths(ref_seq, l_homopolymer_lens_, r_homopolymer_lens_);
  }

  virtual ~HapBlock() {
    for (unsigned int i = 0; i < l_homopolymer_lens_.size(); i++)
      delete [] l_homopolymer_lens_[i];
    for (unsigned int i = 0; i < r_homopolymer_lens_.size(); i++)
      delete [] r_homopolymer_lens_[i];
    l_homopolymer_lens_.clear();
    r_homopolymer_lens_.clear();
  }

  virtual RepeatStutterInfo* get_repeat_info()                    { return NULL; }
  virtual StutterAlignerClass* get_stutter_aligner(int seq_index) { return NULL; }

  int32_t start()   const { return start_; }
  int32_t end()     const { return end_;   }
  int num_options() const { return 1 + alt_seqs_.size(); }
  int min_size()    const { return min_size_; }
  int max_size()    const { return max_size_; }
  bool contains(const std::string& seq) const { return seq_set_.find(seq) != seq_set_.end(); }

  virtual void add_alternate(const std::string& alt) {
    alt_seqs_.push_back(alt);
    min_size_ = std::min(min_size_, (int)alt.size());
    max_size_ = std::max(max_size_, (int)alt.size());
    if (alt_seqs_.size() == 1)
      suffix_matches_.push_back(length_suffix_match(ref_seq_, alt));
    else
      suffix_matches_.push_back(length_suffix_match(alt_seqs_[alt_seqs_.size()-2], alt));
    seq_set_.insert(alt);
    calc_homopolymer_lengths(alt, l_homopolymer_lens_, r_homopolymer_lens_);
  }

  void print(std::ostream& out) const;

  int size(int index) const {
    if (index == 0)
      return ref_seq_.size();
    else if (index-1 < alt_seqs_.size())
      return alt_seqs_[index-1].size();
    else
      throw std::out_of_range("Index out of bounds in HapBlock::size()");
  }

  const std::string& get_seq(unsigned int index) const {
    if (index == 0)
      return ref_seq_;
    else if (index-1 < alt_seqs_.size())
      return alt_seqs_[index-1];
    else
      throw std::out_of_range("Index out of bounds in HapBlock::get_seq()");
  }

  inline int suffix_match_len(unsigned int seq_index) const {
    assert(seq_index < suffix_matches_.size());
    return suffix_matches_[seq_index];
  }

  inline unsigned int left_homopolymer_len(unsigned int seq_index, int base_index) const {
    assert(seq_index < l_homopolymer_lens_.size());
    return (l_homopolymer_lens_[seq_index] == NULL ? 0 : l_homopolymer_lens_[seq_index][base_index]);
  }

  inline unsigned int right_homopolymer_len(unsigned int seq_index, int base_index) const {
    assert(seq_index < r_homopolymer_lens_.size());
    return (r_homopolymer_lens_[seq_index] == NULL ? 0 : r_homopolymer_lens_[seq_index][base_index]);
  }
  
  virtual HapBlock* reverse(){
    std::string rev_ref_seq = ref_seq_;
    std::reverse(rev_ref_seq.begin(), rev_ref_seq.end());
    HapBlock* rev_block = new HapBlock(end_-1, start_-1, rev_ref_seq);
    for (unsigned int i = 0; i < alt_seqs_.size(); i++) {
      std::string alt = alt_seqs_[i];
      std::reverse(alt.begin(), alt.end());
      rev_block->add_alternate(alt);
    }
    return rev_block;
  }

  int index_of(const std::string& seq) const {
    if (seq.compare(ref_seq_) == 0)
      return 0;
    for (unsigned int i = 0; i < alt_seqs_.size(); i++)
      if (seq.compare(alt_seqs_[i]) == 0)
	return i+1;
    return -1;
  }

  virtual HapBlock* remove_alleles(const std::vector<int>& allele_indices){
    std::set<int> bad_alleles(allele_indices.begin(), allele_indices.end());
    assert(bad_alleles.find(0) == bad_alleles.end());

    HapBlock* new_block = new HapBlock(start_, end_, ref_seq_);
    for (unsigned int i = 0; i < alt_seqs_.size(); i++)
      if (bad_alleles.find(i+1) == bad_alleles.end())
        new_block->add_alternate(alt_seqs_[i]);
    return new_block;
  }
};

#endif
