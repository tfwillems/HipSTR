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
#include "StutterAligner.h"

class HapBlock {
 protected:
  std::vector<std::string> seqs_;
  std::set<std::string> seq_set_;
  int32_t start_; // Start of region (inclusive)
  int32_t end_;   // End of region (not inclusive)
  int min_size_, max_size_;

  std::vector<int*> l_homopolymer_lens_;
  std::vector<int*> r_homopolymer_lens_;

  // Compute the homopolymer lengths and store them in the resulting vectors
  void calc_homopolymer_lengths(const std::string& seq, std::vector<int*>& llen_vec, std::vector<int*>& rlen_vec);

  // Private unimplemented copy constructor and assignment operator to prevent operations
  HapBlock(const HapBlock& other);
  HapBlock& operator=(const HapBlock& other);

 public:
  HapBlock(int32_t start, int32_t end, const std::string& ref_seq){
    start_    = start;
    end_      = end;
    min_size_ = (int)ref_seq.size();
    max_size_ = (int)ref_seq.size();
    seqs_.push_back(ref_seq);
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

  virtual RepeatStutterInfo* get_repeat_info()               { return NULL; }
  virtual StutterAligner* get_stutter_aligner(int seq_index) { return NULL; }

  int32_t start()   const { return start_;       }
  int32_t end()     const { return end_;         }
  int num_options() const { return seqs_.size(); }
  int min_size()    const { return min_size_;    }
  int max_size()    const { return max_size_;    }
  bool contains(const std::string& seq) const { return seq_set_.find(seq) != seq_set_.end(); }

  virtual void add_alternate(const std::string& alt) {
    seqs_.push_back(alt);
    min_size_ = std::min(min_size_, (int)alt.size());
    max_size_ = std::max(max_size_, (int)alt.size());
    seq_set_.insert(alt);
    calc_homopolymer_lengths(alt, l_homopolymer_lens_, r_homopolymer_lens_);
  }

  void print(std::ostream& out) const;

  int size(int index) const {
    if (index < seqs_.size() && index >= 0)
      return seqs_[index].size();
    else
      throw std::out_of_range("Index out of bounds in HapBlock::size()");
  }

  const std::string& get_seq(unsigned int index) const {
    if (index < seqs_.size())
      return seqs_[index];
    else
      throw std::out_of_range("Index out of bounds in HapBlock::get_seq()");
  }

  inline unsigned int left_homopolymer_len(unsigned int seq_index, int base_index) const {
    return (l_homopolymer_lens_[seq_index] == NULL ? 0 : l_homopolymer_lens_[seq_index][base_index]);
  }

  inline unsigned int right_homopolymer_len(unsigned int seq_index, int base_index) const {
    return (r_homopolymer_lens_[seq_index] == NULL ? 0 : r_homopolymer_lens_[seq_index][base_index]);
  }
  
  virtual HapBlock* reverse(){
    std::string rev_ref_seq = seqs_[0];
    std::reverse(rev_ref_seq.begin(), rev_ref_seq.end());
    HapBlock* rev_block = new HapBlock(end_-1, start_-1, rev_ref_seq);
    for (unsigned int i = 1; i < seqs_.size(); ++i){
      std::string alt = seqs_[i];
      std::reverse(alt.begin(), alt.end());
      rev_block->add_alternate(alt);
    }
    return rev_block;
  }

  int index_of(const std::string& seq) const {
    for (unsigned int i = 0; i < seqs_.size(); ++i)
      if (seq.compare(seqs_[i]) == 0)
	return i;
    return -1;
  }

  virtual HapBlock* remove_alleles(const std::vector<int>& allele_indices){
    std::set<int> bad_alleles(allele_indices.begin(), allele_indices.end());
    assert(bad_alleles.find(0) == bad_alleles.end());

    HapBlock* new_block = new HapBlock(start_, end_, seqs_[0]);
    for (unsigned int i = 1; i < seqs_.size(); ++i)
      if (bad_alleles.find(i) == bad_alleles.end())
        new_block->add_alternate(seqs_[i]);
    return new_block;
  }
};

#endif
