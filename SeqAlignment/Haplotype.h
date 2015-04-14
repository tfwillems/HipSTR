#ifndef HAPLOTYPE_H_
#define HAPLOTYPE_H_

#include <iostream>
#include <string>
#include <vector>

#include "HapBlock.h"

class Haplotype {
 private:
  std::vector<HapBlock*> blocks_;
  std::vector<int> nopts_;
  int ncombs_, cur_size_;

  // Variables for graycode-based iterator
  int counter_, last_changed_, max_size_;
  std::vector<int> dirs_, factors_, counts_, nchanges_;

  // Information about repetitive regions [start, stop)
  std::vector<int> repeat_starts_, repeat_stops_;
 
  void init();
  
  int left_homopolymer_len(char c, int block_index);
  int right_homopolymer_len(char c, int block_index);

 public:
  Haplotype(std::vector<HapBlock*>& blocks) {
    int32_t ref_coord = blocks[0]->start();
    max_size_ = 0;
    for (int i = 0; i < blocks.size(); i++) {
      if (blocks[i]->get_repeat_info() != NULL) {
	repeat_starts_.push_back(ref_coord);
	repeat_stops_.push_back(ref_coord + blocks[i]->get_seq(0).size());
      }
      blocks_.push_back(blocks[i]);
      nopts_.push_back(blocks[i]->num_options());
      max_size_ += blocks[i]->max_size();
      ref_coord += blocks[i]->get_seq(0).size();
    }

    dirs_.resize(blocks_.size());
    factors_.resize(blocks_.size());
    counts_.resize(blocks_.size());
    nchanges_.resize(blocks_.size());
    init();
  }

  void print_nchanges(std::ostream& out) {
    for (int i = 0; i < blocks_.size(); i++)
      out << nchanges_[i] << " ";
    out << std::endl;
  }
  
  void print_counter_state(std::ostream& out) {
    for (int i = 0; i < blocks_.size(); i++)
      out << counts_[i] << " ";
    out << std::endl;
  }
  
  inline const std::string& get_seq(int block_index) {
    return blocks_[block_index]->get_seq(counts_[block_index]);
  }

  inline const std::vector<int32_t>& get_repeat_starts() {
    return repeat_starts_;
  }

  inline const std::vector<int32_t>& get_repeat_stops() {
    return repeat_stops_;
  }

  inline HapBlock* get_block(int block_index) {
    return blocks_[block_index];
  }

  inline char get_first_char() {
    return blocks_[0]->get_seq(counts_[0])[0];
  }

  inline char get_last_char() {
    return blocks_.back()->get_seq(counts_[blocks_.size()-1]).back();
  }

  inline int num_blocks()   const { return blocks_.size(); }
  inline int num_combs()    const { return ncombs_; }
  inline int last_changed() const { return last_changed_; }
  inline int max_size()     const { return max_size_; }
  inline int cur_size()     const { return cur_size_; }
  inline int cur_index(int block_index) const { return counts_[block_index]; }

  int num_options(int block_index) const {
    return blocks_[block_index]->num_options();
  }

  void print(std::ostream& out) {
    for (int i = 0; i < num_blocks(); i++)
      out << get_seq(i);
    out << std::endl;
  }
  
  void print_block_structure(int max_ref_len,
			     int max_other_len, 
			     std::ostream& out);

  void print_padded(std::ostream& out) {
    for (int i = 0; i < num_blocks(); i++) {
      const std::string& seq = get_seq(i);
      out << seq;
      for (int j = 0; j < blocks_[i]->max_size()+1-seq.size(); j++)
	out << " ";
    }
    out << std::endl;
  }
    
  void reset();

  bool next();

  int homopolymer_length(int block_index, int base_index);
};

#endif


