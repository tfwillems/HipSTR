#ifndef STR_ALLELE_EXPANSION_H_
#define STR_ALLELE_EXPANSION_H_

#include <map>
#include <set>
#include <string>
#include <vector>

void get_insertion_matches(const std::string& haplotype_allele, const std::string& read_seq, int ins_size, int period,
			   std::vector<int>& prefix_matches, std::vector<int>& suffix_matches, std::set<std::string>& existing,
			   std::map<std::string, int>& new_seq_counts);

void get_deletion_matches(const std::string& haplotype_allele, const std::string& read_seq, int del_size, int period,
			  std::vector<int>& prefix_matches, std::vector<int>& suffix_matches,
			  std::set<std::string>& new_seqs);

void get_candidates(std::vector<std::string>& str_seqs, std::vector< std::vector<std::string> >& read_seqs, int period,
		    std::set<std::string>& new_str_seqs);

#endif
