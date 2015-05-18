#include "../SeqAlignment/STRAlleleExpansion.h"

#include <iostream>
#include <set>
#include <string>
#include <vector>

int main(){
  std::set<std::string> new_str_seqs;
  std::vector<std::string> str_seqs, read_seqs;

  str_seqs.push_back( "CGCGCGTATATGGG");
  read_seqs.push_back("CGCGCGTATATATGGG");
  std::cerr << str_seqs.back() << std::endl
	    << "vs" << std::endl
	    << read_seqs.back() << std::endl;
  //get_candidates(str_seqs, read_seqs, 2, new_str_seqs);
  std::cerr << std::endl;
  
  new_str_seqs.clear();
  read_seqs.clear();
  read_seqs.push_back("ATATGGGACCGCGCGTATGGG");
  std::cerr << str_seqs.back() << std::endl
	    << "vs" << std::endl
	    << read_seqs.back() << std::endl;
  //get_candidates(str_seqs, read_seqs, 2, new_str_seqs);
  std::cerr << std::endl;  

  new_str_seqs.clear();
  read_seqs.clear();
  read_seqs.push_back("CGGATTCGCGCGTATATATATATGGGTTA");
  std::cerr << str_seqs.back() << std::endl
	    << "vs" << std::endl
	    << read_seqs.back() << std::endl;
  //get_candidates(str_seqs, read_seqs, 2, new_str_seqs);
  std::cerr << std::endl;  
}
