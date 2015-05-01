#include "Haplotype.h"

void Haplotype::init(){
  ncombs_   = 1;
  cur_size_ = 0;
  for (int i = 0; i < blocks_.size(); i++){
    factors_[i]  = ncombs_;
    ncombs_     *= nopts_[i];
    dirs_[i]     = 1;
    counts_[i]   = 0;
    nchanges_[i] = 0;
    cur_size_   += blocks_[i]->size(0);
  }
  counter_      = 0;
  last_changed_ = -1;
}


void Haplotype::reset(){
  init();
  counter_      = 0;
  last_changed_ = -1;
}

bool Haplotype::next(){
  if (counter_ == ncombs_-1)
    return false;
  
  int t = counter_+1;
  for (int j = blocks_.size()-1; j >= 0; j--){
    t %= factors_[j];
    if (t == 0){
      last_changed_ = j;
      nchanges_[j] += 1;
      cur_size_    -= blocks_[j]->size(counts_[j]);
      counts_[j]   += dirs_[j];
      cur_size_    += blocks_[j]->size(counts_[j]);
      if (counts_[j] == 0 || counts_[j] == nopts_[j]-1)
	dirs_[j] *= -1;
      break;
    }
  }
  counter_++;
  return true;
}

void Haplotype::print_block_structure(int max_ref_len, 
				      int max_other_len, 
				      std::ostream& out){
  int max_rows = 0;
  for (int i = 0; i < num_blocks(); i++)
    max_rows = std::max(max_rows, blocks_[i]->num_options());
  
  for (int n = 0; n < max_rows; n++){
    for (int i = 0; i < num_blocks(); i++){
      int char_limit = blocks_[i]->num_options() == 1 ? max_ref_len : max_other_len;
      int num_chars  = 0;
      if (n < blocks_[i]->num_options()){
	num_chars = blocks_[i]->get_seq(n).size();
	if (num_chars > char_limit){
	  int v1 = char_limit/2;
	  int v2 = char_limit - v1 - 3;
	  out << blocks_[i]->get_seq(n).substr(0, v1) << "..." << blocks_[i]->get_seq(n).substr(blocks_[i]->get_seq(n).size()-v2, v2);
	  num_chars = char_limit;
	}
	else
	  out << blocks_[i]->get_seq(n);
      }
      while (num_chars <= std::min(blocks_[i]->max_size(), char_limit)){
	out << " ";
	num_chars++;
      }
    }
    out << std::endl;
  }
}

unsigned int Haplotype::left_homopolymer_len(char c, int block_index){
  unsigned int total = 0;
  while (block_index >= 0){
    const std::string& seq = get_seq(block_index);
    if (*seq.rbegin() == c){
      unsigned int llen = blocks_[block_index]->left_homopolymer_len(counts_[block_index], seq.size()-1);
      total += (1 + llen);
      if (llen != seq.size())
	break;
    }
    else
      break;
    block_index--;
  }
  return total;
}
 
unsigned int Haplotype::right_homopolymer_len(char c, int block_index){
  unsigned int total = 0;
  while (block_index < blocks_.size()){
    const std::string& seq = get_seq(block_index);
    if (seq[0] == c){
      unsigned int rlen = blocks_[block_index]->right_homopolymer_len(counts_[block_index], 0);
      total   += (1 + rlen); 
      if (rlen != seq.size())
	break;
    }
    else
      break;
    block_index++;
  }
  return total;
}


unsigned int Haplotype::homopolymer_length(int block_index, int base_index){
  HapBlock* block        = blocks_[block_index];
  const std::string& seq = block->get_seq(counts_[block_index]);
  unsigned int llen = block->left_homopolymer_len(counts_[block_index],  base_index);
  unsigned int rlen = block->right_homopolymer_len(counts_[block_index], base_index);
  if (base_index-llen == 0)
    llen += left_homopolymer_len(seq[base_index], block_index-1);
  if (base_index+rlen == seq.size()-1)
    rlen += right_homopolymer_len(seq[base_index], block_index+1);
  return llen + rlen + 1;
}
