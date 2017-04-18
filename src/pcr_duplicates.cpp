#include "pcr_duplicates.h"

#include <algorithm>
#include <assert.h>
#include <iostream>
#include <string>

std::string get_library(const BamAlignment& aln, const std::map<std::string, std::string>& rg_to_library){
  std::string rg;
  if (!aln.GetStringTag("RG", rg))
    printErrorAndDie("Failed to retrieve BAM alignment's RG tag");
  auto iter = rg_to_library.find(aln.Filename() + rg);
  if (iter == rg_to_library.end())
    printErrorAndDie("No library found for read group " + rg + " in BAM file headers");
  return iter->second;
}

void remove_pcr_duplicates(const BaseQuality& base_quality, bool use_bam_rgs,
			   const std::map<std::string, std::string>& rg_to_library,
			   std::vector< std::vector<BamAlignment> >& paired_strs_by_rg,
			   std::vector< std::vector<BamAlignment> >& mate_pairs_by_rg,
			   std::vector< std::vector<BamAlignment> >& unpaired_strs_by_rg, std::ostream& logger){
  int32_t dup_count = 0;
  assert(paired_strs_by_rg.size() == mate_pairs_by_rg.size() && paired_strs_by_rg.size() == unpaired_strs_by_rg.size());
  for (size_t i = 0; i < paired_strs_by_rg.size(); i++){
    assert(paired_strs_by_rg[i].size() == mate_pairs_by_rg[i].size());

    std::vector<ReadPair> read_pairs;
    for (size_t j = 0; j < paired_strs_by_rg[i].size(); j++){
      std::string library = use_bam_rgs ? get_library(paired_strs_by_rg[i][j], rg_to_library): rg_to_library.find(paired_strs_by_rg[i][j].Filename())->second;
      read_pairs.push_back(ReadPair(paired_strs_by_rg[i][j], mate_pairs_by_rg[i][j], library));
    }
    for (size_t j = 0; j < unpaired_strs_by_rg[i].size(); j++){
      std::string library = use_bam_rgs ? get_library(unpaired_strs_by_rg[i][j], rg_to_library): rg_to_library.find(unpaired_strs_by_rg[i][j].Filename())->second;
      read_pairs.push_back(ReadPair(unpaired_strs_by_rg[i][j], library));
    }
    std::sort(read_pairs.begin(), read_pairs.end());

    paired_strs_by_rg[i].clear();
    mate_pairs_by_rg[i].clear();
    unpaired_strs_by_rg[i].clear();
    if (read_pairs.size() == 0)
      continue;

    // When both mates in a pair overlap the STR, they generate pseudo PCR duplicates because the read pair is included twice in the input (but reversed).
    // To use both reads for genotyping, we don't want to remove these duplicates for downstream analysis. Instead, we use this flag to track
    // if this issue has occurred and undo the duplicate removal when saving the alignments
    bool include_rev  = false;
    size_t best_index = 0;
    for (size_t j = 1; j < read_pairs.size(); j++){
      if (read_pairs[j].duplicate(read_pairs[best_index])){
	dup_count++;
	// Update index if new pair's STR read has a higher total base quality
	if (base_quality.sum_log_prob_correct(read_pairs[j].aln_one().Qualities()) >
	    base_quality.sum_log_prob_correct(read_pairs[best_index].aln_one().Qualities())){
	  best_index  = j;
	  include_rev = (read_pairs[best_index].name().compare(read_pairs[j-1].name()) == 0);
	}
	else if (j == best_index+1)
	  include_rev |= (read_pairs[best_index].name().compare(read_pairs[j].name()) == 0);
      }
      else {
	// Keep best pair from prior set of duplicates
	if (read_pairs[best_index].single_ended())
	  unpaired_strs_by_rg[i].push_back(read_pairs[best_index].aln_one());
	else {
	  paired_strs_by_rg[i].push_back(read_pairs[best_index].aln_one());
	  mate_pairs_by_rg[i].push_back(read_pairs[best_index].aln_two());
	  if (include_rev){
	    dup_count--;
	    paired_strs_by_rg[i].push_back(read_pairs[best_index].aln_two());
	    mate_pairs_by_rg[i].push_back(read_pairs[best_index].aln_one());
	  }
	}
	best_index  = j; // Update index for new set of duplicates
	include_rev = false;
      }
    }

    // Keep best pair for last set of duplicates
    if (read_pairs[best_index].single_ended())
      unpaired_strs_by_rg[i].push_back(read_pairs[best_index].aln_one());
    else {
      paired_strs_by_rg[i].push_back(read_pairs[best_index].aln_one());
      mate_pairs_by_rg[i].push_back(read_pairs[best_index].aln_two());
      if (include_rev){
	dup_count--;
	paired_strs_by_rg[i].push_back(read_pairs[best_index].aln_two());
	mate_pairs_by_rg[i].push_back(read_pairs[best_index].aln_one());
      }
    }
  }
  logger << "Removed " << dup_count << " sets of PCR duplicate reads" << std::endl;
}
