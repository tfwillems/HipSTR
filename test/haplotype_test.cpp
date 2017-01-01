#include <algorithm>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <string>

#include "../src/SeqAlignment/HapBlock.h"
#include "../src/SeqAlignment/Haplotype.h"
#include "../src/SeqAlignment/RepeatBlock.h"
#include "../src/stutter_model.h"

int main(){
  std::string l1   = "ACGGTATC",   l2 = "ACGGTCTC";
  std::string r1   = "GAATCCC" ,   r2 = "GAATTCC";
  std::string rep1 = "ATAT",     rep2 = "ATATAT", rep3 = "ATATATAT";

  HapBlock left_flank(0, 8, l1); 
  left_flank.add_alternate(l2);

  StutterModel stutter_model(0.9,  0.01,  0.02, 0.7, 0.001, 0.001, 2);
  RepeatBlock rep_block(8, 12, rep1, 2, &stutter_model); 
  rep_block.add_alternate(rep2); 
  rep_block.add_alternate(rep3);
 
  HapBlock right_flank(12, 19, r1);
  right_flank.add_alternate(r2);

  std::vector<HapBlock*> hap_blocks;
  hap_blocks.push_back(&left_flank);
  hap_blocks.push_back(&rep_block);
  hap_blocks.push_back(&right_flank);

  Haplotype haplotype(hap_blocks);
  std::vector<HapBlock*> rev_blocks;
  Haplotype* rev_haplotype = haplotype.reverse(rev_blocks);
  haplotype.print_block_structure(100, 100, true, std::cout);
  rev_haplotype->print_block_structure(100, 100, true, std::cout);
  do {
    std::string s1 = haplotype.get_seq();
    std::string s2 = rev_haplotype->get_seq();
    std::cout << s1 << std::endl
	      << s2 << std::endl;
    std::reverse(s2.begin(), s2.end());
    assert(s2.compare(s1) == 0);
  }
  while(haplotype.next() && rev_haplotype->next());
  assert(rev_haplotype->get_block(1)->get_repeat_info() != NULL);

  for (unsigned int i = 0; i < 100; i++){
    int hap_index = rand() % haplotype.num_combs();
    haplotype.go_to(hap_index);
    rev_haplotype->go_to(hap_index);
    std::string s1 = haplotype.get_seq();
    std::string s2 = rev_haplotype->get_seq();
    std::reverse(s2.begin(), s2.end());
    assert(s2.compare(s1) == 0);
  }

  // Delete datastructures for reverse haplotype
  for (unsigned int i = 0; i < rev_blocks.size(); i++)
    delete rev_blocks[i];
  delete rev_haplotype;
}
