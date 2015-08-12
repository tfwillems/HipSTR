#include <assert.h>

#include "../error.h"
#include "Haplotype.h"
#include "NWNoRefEndPenalty.h"


void Haplotype::adjust_indels(std::string& ref_hap_al, std::string& alt_hap_al){
  assert(blocks_.size() == 3);
  int32_t ref_pos = blocks_[0]->start(), str_pos = blocks_[1]->start();
  int aln_index   = 0;
  while (aln_index < alt_hap_al.size()){
    if (alt_hap_al[aln_index] == '-' && ref_pos < str_pos){
      // If necessary and possible, move deletion to the right until it lies fully within the repeat block
      int index = aln_index;
      while (index < alt_hap_al.size() && alt_hap_al[index] == '-')
	index++;
      int pos       = ref_pos;
      int del_index = aln_index;
      int del_size  = index - aln_index;
      while (index < alt_hap_al.size() && pos < str_pos && (ref_hap_al[del_index] == ref_hap_al[index])){
	alt_hap_al[del_index] = alt_hap_al[index];
	alt_hap_al[index]     = '-';
	index++;
	del_index++;
	pos++;
      }

      aln_index = index;
      ref_pos   = pos+del_size;
    }
    else if (ref_hap_al[aln_index] == '-' && ref_pos < str_pos){
      // If necessary and possible, move insertion to the right until it lies directly in front of the repeat block
      int index = aln_index;
      while (index < ref_hap_al.size() && ref_hap_al[index] == '-')
	index++;
      int pos       = ref_pos;
      int ins_index = aln_index;
      while (index < ref_hap_al.size() && pos < str_pos && (alt_hap_al[ins_index] == alt_hap_al[index])){
	ref_hap_al[ins_index] = ref_hap_al[index];
	ref_hap_al[index]     = '-';
	index++;
	ins_index++;
	pos++;
      }

      aln_index = index;
      ref_pos   = pos;
    }
    else {
      if (ref_hap_al[aln_index] != '-')
	ref_pos++;
      aln_index++;
    }
  }
}


void Haplotype::aln_haps_to_ref(){
  std::string ref_hap_seq = get_seq(), alt_hap_seq;
  std::string ref_hap_al, alt_hap_al;
  float score;
  std::vector<BamTools::CigarOp> cigar_list;

  do {
    alt_hap_seq = get_seq();
    if (!NWNoRefEndPenalty::LeftAlign(ref_hap_seq, alt_hap_seq, ref_hap_al, alt_hap_al, &score, cigar_list))
      printErrorAndDie("Failed to left-align haplotype sequence to reference allele");
    cigar_list.clear();

    // Attempt to merge indels inside of repeat block
    adjust_indels(ref_hap_al, alt_hap_al);

    std::string aln_info = ""; aln_info.reserve(alt_hap_al.size());
    for (unsigned int i = 0; i < alt_hap_al.size(); i++){
      if (ref_hap_al[i] == '-')
	aln_info += 'I';
      else if (alt_hap_al[i] == '-')
	aln_info += 'D';
      else
	aln_info += 'M';
    }
    hap_aln_info_.push_back(aln_info);
  }
  while (next());
  reset();
}

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
  if (fixed_)
    return false;
  if (counter_ == ncombs_-1)
    return false;

  int index = -1;
  int t     = counter_+1;
  for (int j = blocks_.size()-1; j >= 0; j--){
    t %= factors_[j];
    if (t == 0){
      index = j;
      break;
    }
  }
  if (inc_rev_)
    index = blocks_.size()-1-index;

  assert(index != -1);
  last_changed_     = index;
  nchanges_[index] += 1;
  cur_size_        -= blocks_[index]->size(counts_[index]);
  counts_[index]   += dirs_[index];
  cur_size_        += blocks_[index]->size(counts_[index]);
  if (counts_[index] == 0 || counts_[index] == nopts_[index]-1)
    dirs_[index] *= -1;

  counter_++;
  return true;
}

void Haplotype::go_to(int hap_index){
  if (hap_index < 0 || hap_index >= ncombs_)
    printErrorAndDie("Invalid haplotype index in go_to()");
  if (hap_index < counter_)
    reset();
  while (counter_ < hap_index)
    next();
  last_changed_ = -1; // Don't want to reuse haplotype info as we're no longer just incrementing
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

Haplotype* Haplotype::reverse( std::vector<HapBlock*>& rev_blocks){
  assert(rev_blocks.size() == 0);
  for (unsigned int i = 0; i < blocks_.size(); i++)
    rev_blocks.push_back(blocks_[i]->reverse());
  std::reverse(rev_blocks.begin(), rev_blocks.end());
  Haplotype* rev_hap = new Haplotype(rev_blocks);
  rev_hap->inc_rev_  = true;

  // Set the haplotype alignments to the reverse of those in the current haplotype
  // Can't use the values obtained from the constructor as the reverse
  // haplotype would have right align indels instead of left aligning them
  rev_hap->hap_aln_info_ = hap_aln_info_;
  for (auto iter = rev_hap->hap_aln_info_.begin(); iter != rev_hap->hap_aln_info_.end(); iter++)
    std::reverse(iter->begin(), iter->end());
  return rev_hap;
}
