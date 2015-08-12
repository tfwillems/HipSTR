#include "STRAlleleExpansion.h"

#include <assert.h>
#include <iostream>
#include <map>
#include "../zalgorithm.h"

bool matches_pattern(std::string& hap_seq, std::string& read_seq, int period){
  std::vector<int> prefix_matches;
  ZAlgorithm::GetPrefixMatchCounts(read_seq, hap_seq, prefix_matches);
  for (unsigned int i = 0; i < prefix_matches.size(); i += period)
    if (prefix_matches[i] == read_seq.size())
      return true;
  return false;
}

void get_insertion_matches(const std::string& haplotype_allele, const std::string& read_seq, int ins_size, int period,
			   std::vector<int>& prefix_matches, std::vector<int>& suffix_matches, std::set<std::string>& existing,
			   std::map<std::string, int>& new_seq_counts){
  assert(ins_size > 0);
  assert(read_seq.size() == prefix_matches.size() && read_seq.size() == suffix_matches.size());
  int hap_len = haplotype_allele.size();
  if (read_seq.size() < hap_len + ins_size)
    return;

  // Iterate through all potential matching sequences
  for (int start_index = 0; start_index < read_seq.size(); start_index++){
    int end_index = start_index + hap_len + ins_size - 1;
    if (end_index >= read_seq.size())
      break;

    // Check if insertion comprises prefix of candidate sequence and matches haplotype prefix
    if (prefix_matches[start_index+ins_size] >= hap_len){
      if (prefix_matches[start_index] >= ins_size){
	std::string match = read_seq.substr(start_index, end_index-start_index+1);
	if (existing.find(match) == existing.end())
	  new_seq_counts[match] += 1;
	continue;
      }
    }

    // Check if insertion comprises suffix of candidate sequence and matches haplotype suffix
    if (prefix_matches[start_index] >= hap_len){
      if (suffix_matches[end_index] >= ins_size){
	std::string match = read_seq.substr(start_index, end_index-start_index+1);
	if (existing.find(match) == existing.end())
	  new_seq_counts[match] += 1;
	continue;
      }
    }

    // Check if insertion in interior of sequence
    int min_start = std::max(start_index+1, end_index - suffix_matches[end_index] - ins_size + 1);
    int max_start = std::min(end_index-1,  start_index + prefix_matches[start_index]);
    for (int ins_start = min_start; ins_start <= max_start; ins_start++){
      std::string ins_seq      = read_seq.substr(ins_start, ins_size);
      int min_surr_hap         = std::max(0,         ins_start-start_index-ins_size);
      int max_surr_hap         = std::min(hap_len-1, ins_start-start_index+ins_size-1);
      std::string surr_hap_seq = haplotype_allele.substr(min_surr_hap, max_surr_hap-min_surr_hap+1);
      if (matches_pattern(surr_hap_seq, ins_seq, period)){
	std::string match = read_seq.substr(start_index, end_index-start_index+1);
	if (existing.find(match) == existing.end())
	  new_seq_counts[match] += 1;
	break;
      }
    }
  }
}

void get_deletion_matches(const std::string& haplotype_allele, const std::string& read_seq, int del_size,
			  std::vector<int>& prefix_matches, std::vector<int>& suffix_matches, 
			  std::set<std::string>& new_seqs){
  assert(del_size < 0 && (del_size + haplotype_allele.size() > 0));
  assert(read_seq.size() == prefix_matches.size() && read_seq.size() == suffix_matches.size());
  int hap_len = haplotype_allele.size();
  if (read_seq.size() <= hap_len + del_size)
    return;

  // Iterate through all potential matching sequences
  for (int start_index = 0; start_index < read_seq.size(); start_index++){
    int end_index = start_index + hap_len + del_size - 1;
    if (end_index >= read_seq.size())
      break;

    // NOTE: Don't allow haplotype prefixes or suffixes as deletions as this could lead to way too many hits
    int min_pos = std::max(start_index+1, end_index - suffix_matches[end_index]);
    int max_pos = std::min(end_index-1, start_index + prefix_matches[start_index] - 1);
    for (int del_pos = min_pos; del_pos <= max_pos; del_pos++){
      // TO DO: Check if deleted sequence matches surrounding haplotype sequence
      if (true){
	new_seqs.insert(read_seq.substr(start_index, end_index-start_index+1));
	break;
      }
    }
  }
}


void get_candidates(std::vector<std::string>& str_seqs, std::vector< std::vector<std::string> >& read_seqs, int period,
		    std::set<std::string>& new_str_seqs, std::ostream& logger){
  assert(new_str_seqs.size() == 0);
  std::set<std::string> existing(str_seqs.begin(), str_seqs.end());
  logger << "Beginning additional allele identification..." << std::endl;
  for (unsigned int i = 0; i < read_seqs.size(); i++){
    std::map<std::string, int> match_counts;
    for (unsigned int j = 0; j < read_seqs[i].size(); j++){
      for (unsigned int k = 0; k < str_seqs.size(); k++){
	const std::string& hap_seq = str_seqs[k];
	std::vector<int> prefix_matches, suffix_matches;
	ZAlgorithm::GetPrefixMatchCounts(hap_seq, read_seqs[i][j], prefix_matches);
	ZAlgorithm::GetSuffixMatchCounts(hap_seq, read_seqs[i][j], suffix_matches);
	for (int ins_size = period; ins_size <= 3*period; ins_size += period)
	  get_insertion_matches(hap_seq, read_seqs[i][j], ins_size, period, prefix_matches, suffix_matches, existing, match_counts);
      }
    }

    // Check sample frequency counts
    for (auto iter = match_counts.begin(); iter != match_counts.end(); iter++)
      if (iter->second >= 2 && (1.0*iter->second/read_seqs[i].size() >= 0.2))
	new_str_seqs.insert(iter->first);
  }

  logger << "Existing alleles: " << std::endl;
  for (auto iter = str_seqs.begin(); iter != str_seqs.end(); iter++)
    logger << *iter << std::endl;

  logger << "Identified " << new_str_seqs.size() << " additional candidate STR alleles:" << std::endl;
  for (auto iter = new_str_seqs.begin(); iter != new_str_seqs.end(); iter++)
    logger << *iter << std::endl;
}
