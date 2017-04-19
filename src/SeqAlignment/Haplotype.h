#ifndef HAPLOTYPE_H_
#define HAPLOTYPE_H_

#include <assert.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "HapBlock.h"

class Haplotype {
 private:
  std::vector<HapBlock*> blocks_;
  std::vector<int> nopts_;
  int ncombs_, cur_size_;
  bool fixed_;

  // Variables for graycode-based iterator
  int counter_, last_changed_, max_size_;
  std::vector<int> dirs_, factors_, counts_, nchanges_;
  bool inc_rev_; // Iff true, increment from back to front

  void init();
  
  unsigned int left_homopolymer_len(char c, int block_index) const;
  unsigned int right_homopolymer_len(char c, int block_index) const;

  std::vector<std::string> hap_aln_info_;
  void aln_haps_to_ref();
  void adjust_indels(std::string& ref_hap_al, std::string& alt_hap_al);

 public:
  explicit Haplotype(std::vector<HapBlock*>& blocks){
    max_size_ = 0;
    for (unsigned int i = 0; i < blocks.size(); i++) {
      blocks_.push_back(blocks[i]);
      nopts_.push_back(blocks[i]->num_options());
      max_size_ += blocks[i]->max_size();
    }
    fixed_ = false;

    dirs_.resize(blocks_.size());
    factors_.resize(blocks_.size());
    counts_.resize(blocks_.size());
    nchanges_.resize(blocks_.size());
    inc_rev_ = false;
    init();
    aln_haps_to_ref();
  }

  void print_nchanges(std::ostream& out) const {
    for (unsigned int i = 0; i < blocks_.size(); i++)
      out << nchanges_[i] << " ";
    out << std::endl;
  }
  
  void print_counter_state(std::ostream& out) const {
    for (unsigned int i = 0; i < blocks_.size(); i++)
      out << counts_[i] << " ";
    out << std::endl;
  }
  
  inline const std::string& get_seq(int block_index) const { return blocks_[block_index]->get_seq(counts_[block_index]); }
  inline const std::string& get_aln_info()           const { return hap_aln_info_[counter_]; }
  inline char get_first_char()                       const { return blocks_[0]->get_seq(counts_[0])[0]; }
  inline char get_last_char()                        const {
    const std::string& seq = blocks_.back()->get_seq(counts_[blocks_.size()-1]);
    return seq[seq.size()-1];
  }
  inline HapBlock* get_block(int block_index)        const { return blocks_[block_index]; }
  inline HapBlock* get_first_block()                 const { return blocks_.front();      }
  inline HapBlock* get_last_block()                  const { return blocks_.back();       }
  inline int num_blocks()                            const { return blocks_.size();       }
  inline int num_combs()                             const { return ncombs_;              }
  inline int last_changed()                          const { return last_changed_;        }
  inline int max_size()                              const { return max_size_;            }
  inline int cur_size()                              const { return cur_size_;            }
  inline int cur_index()                             const { return counter_;             }
  inline int cur_index(int block_index)              const { return counts_[block_index]; }
  inline bool reversed()                             const { return inc_rev_;             }
  int num_options(int block_index)                   const { return blocks_[block_index]->num_options(); }

  void get_coordinates(int hap_pos, int& block, int& block_pos) const {
    assert(hap_pos >= 0 && hap_pos < cur_size_);
    for (int i = 0; i < blocks_.size(); i++){
      if (hap_pos < blocks_[i]->size(counts_[i])){
	block     = i;
	block_pos = hap_pos;
	return;
      }
      hap_pos -= blocks_[i]->size(counts_[i]);
    }
    assert(false);
  }

  bool position_to_haplotype_index(int32_t pos, int& haplotype_index) const;

  std::string get_seq() const {
    std::stringstream ss;
    for (int i = 0; i < num_blocks(); i++)
      ss << get_seq(i);
    return ss.str();
  }

  void print(std::ostream& out) const {
    for (int i = 0; i < num_blocks(); i++)
      out << get_seq(i);
    out << std::endl;
  }
  
  void print_block_structure(int max_ref_len, int max_other_len, bool indent,
			     std::ostream& out) const;

  // Prevent haplotype from being changed using next()
  void fix()  { fixed_ = true;  }
  void unfix(){ fixed_ = false; }

  void reset();

  bool next();

  void go_to(int hap_index);

  unsigned int homopolymer_length(int block_index, int base_index) const;

  Haplotype* reverse(std::vector<HapBlock*>& rev_blocks);

  void check_indel_clobbering(const std::string& marker, std::vector<bool>& clobbered);
};

#endif
