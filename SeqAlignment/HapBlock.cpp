#include <algorithm>
#include <iostream>

#include "HapBlock.h"

void HapBlock::calc_homopolymer_lengths(){
  for (int i = 0; i <= alt_seqs_.size(); ++i){
    std::string& seq = (i == 0? ref_seq_ : alt_seqs_[i-1]);
    int* llens = new int[seq.size()];
    int* rlens = new int[seq.size()]; 

    llens[0] = 0;
    int count = 0;
    for (int j = 1; j < seq.size(); j++){
      count    = (seq[j-1] == seq[j] ? count + 1 : 0);
      llens[j] = count;
    }
    l_homopolymer_lens.push_back(llens);

    rlens[seq.size()-1] = 0;
    for (int j = seq.size()-2; j >= 0; j--){
      count    = (seq[j+1] == seq[j] ? count + 1: 0);
      rlens[j] = count;
    }
    r_homopolymer_lens.push_back(rlens);
  }
}


bool compareStringLength(const std::string& s1, const std::string& s2){
  return s1.size() < s2.size();
}

void HapBlock::order_alternates_by_length(){
  std::sort(alt_seqs_.begin(), alt_seqs_.end(), compareStringLength);
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
