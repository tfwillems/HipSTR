#include "pcr_duplicates.h"

#include <algorithm>
#include <assert.h>
#include <iostream>
#include <string>

class ReadPair {
private:
  int32_t min_read_start_;
  int32_t max_read_start_;
  BamTools::BamAlignment aln_1_;
  BamTools::BamAlignment aln_2_;
  std::string library_;

public:
  ReadPair(BamTools::BamAlignment& aln_1){
    aln_1_          = aln_1;
    min_read_start_ = -1;
    max_read_start_ = aln_1.Position;
    library_        = ""; // TO DO: Get from read group
  }

  ReadPair(BamTools::BamAlignment& aln_1, BamTools::BamAlignment& aln_2){
    aln_1_          = aln_1;
    aln_2_          = aln_2;
    min_read_start_ = std::min(aln_1.Position, aln_2.Position);
    max_read_start_ = std::max(aln_1.Position, aln_2.Position);
    library_        = ""; // TO DO: Get from read group
  }
  
  BamTools::BamAlignment& aln_one(){ return aln_1_; }
  BamTools::BamAlignment& aln_two(){ return aln_2_; }

  bool single_ended(){
    return min_read_start_ == -1;
  }

  bool duplicate (const ReadPair& pair) const {
    return (library_.compare(pair.library_) == 0)
      && (min_read_start_ == pair.min_read_start_)
      && (max_read_start_ == pair.max_read_start_);
  }

  bool operator < (const ReadPair& pair) const {
    int lib_comp = library_.compare(pair.library_);
    if (lib_comp != 0)
      return lib_comp < 0;
    if (min_read_start_ != pair.min_read_start_)
      return min_read_start_ < pair.min_read_start_;
    return max_read_start_ < pair.max_read_start_;
  }
};


void remove_pcr_duplicates(BaseQuality& base_quality,
			   std::vector< std::vector<BamTools::BamAlignment> >& paired_strs_by_rg,
			   std::vector< std::vector<BamTools::BamAlignment> >& mate_pairs_by_rg,
			   std::vector< std::vector<BamTools::BamAlignment> >& unpaired_strs_by_rg){
  int32_t dup_count = 0;
  assert(paired_strs_by_rg.size() == mate_pairs_by_rg.size() && paired_strs_by_rg.size() == unpaired_strs_by_rg.size());
  for (unsigned int i = 0; i < paired_strs_by_rg.size(); i++){
    assert(paired_strs_by_rg[i].size() == mate_pairs_by_rg[i].size());

    std::vector<ReadPair> read_pairs;
    for (unsigned int j = 0; j < paired_strs_by_rg[i].size(); j++)
      read_pairs.push_back(ReadPair(paired_strs_by_rg[i][j], mate_pairs_by_rg[i][j]));
    for (unsigned int j = 0; j < unpaired_strs_by_rg[i].size(); j++)
      read_pairs.push_back(ReadPair(unpaired_strs_by_rg[i][j]));
    std::sort(read_pairs.begin(), read_pairs.end());

    paired_strs_by_rg[i].clear();
    mate_pairs_by_rg[i].clear();
    unpaired_strs_by_rg[i].clear();
    if (read_pairs.size() == 0)
      continue;
    int best_index = 0;
    for (unsigned int j = 1; j < read_pairs.size(); j++){
      if (read_pairs[j].duplicate(read_pairs[best_index])){
	dup_count++;
	// Update index if new pair's STR read has a higher total base quality
	if (base_quality.sum_log_prob_correct(read_pairs[j].aln_one().Qualities) > 
	    base_quality.sum_log_prob_correct(read_pairs[best_index].aln_one().Qualities))
	  best_index = j;
      }
      else {
	// Keep best pair from prior set of duplicates
	if (read_pairs[best_index].single_ended())
	  unpaired_strs_by_rg[i].push_back(read_pairs[best_index].aln_one());
	else {
	  paired_strs_by_rg[i].push_back(read_pairs[best_index].aln_one());
	  mate_pairs_by_rg[i].push_back(read_pairs[best_index].aln_two());
	}
	best_index = j; // Update index for new set of duplicates
      }
    }

    // Keep best pair for last set of duplicates
    if (read_pairs[best_index].single_ended())
      unpaired_strs_by_rg[i].push_back(read_pairs[best_index].aln_one());
    else {
      paired_strs_by_rg[i].push_back(read_pairs[best_index].aln_one());
      mate_pairs_by_rg[i].push_back(read_pairs[best_index].aln_two());
    }
  }
  std::cerr << "Removed " << dup_count << " sets of PCR duplicate reads" << std::endl;
}
