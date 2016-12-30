#include <algorithm>
#include <iostream>

#include "HapBlock.h"
#include "../stringops.h"

void HapBlock::calc_homopolymer_lengths(std::string& seq, std::vector<int*>& llen_vec, std::vector<int*>& rlen_vec){
  int* llens = new int[seq.size()];
  int* rlens = new int[seq.size()];
  llens[0]   = 0;
  int count  = 0;
  for (unsigned int j = 1; j < seq.size(); j++){
    count    = (seq[j-1] == seq[j] ? count + 1 : 0);
    llens[j] = count;
  }

  rlens[seq.size()-1] = 0;
  for (int j = seq.size()-2; j >= 0; j--){
    count    = (seq[j+1] == seq[j] ? count + 1: 0);
    rlens[j] = count;
  }
  llen_vec.push_back(llens);
  rlen_vec.push_back(rlens);
}

bool compareStringLength(const std::string& s1, const std::string& s2){
  return s1.size() < s2.size();
}

void HapBlock::print(std::ostream& out){
  out << "Haplotype Block: {"          << std::endl 
      << start_ << " -> " << end_      << std::endl
      << "\t"   << "Reference seq:   " << std::endl
      << "\t\t" << ref_seq_            << std::endl
      << "\t"   << "Alternate seq(s):" << std::endl;
  for (unsigned int i = 0; i < alt_seqs_.size(); i++)
    out << "\t\t" << alt_seqs_[i] << std::endl;
  out << "}" << std::endl;
}
