#ifndef PCR_DUPLICATES_H_
#define PCR_DUPLICATES_H_

#include <iostream>
#include <map>
#include <vector>

#include "bam_io.h"
#include "base_quality.h"

class ReadPair {
 private:
  int32_t min_read_start_;
  int32_t max_read_start_;
  BamAlignment aln_1_;
  BamAlignment aln_2_;
  std::string library_;
  std::string name_;

 public:
  ReadPair(BamAlignment& aln_1, std::string& library)
    : aln_1_(aln_1), library_(library), name_(aln_1.Name()){
    min_read_start_ = -1;
    max_read_start_ = aln_1.Position();
  }

  ReadPair(BamAlignment& aln_1, BamAlignment& aln_2, std::string& library)
    : aln_1_(aln_1), aln_2_(aln_2), library_(library), name_(aln_1.Name()){
    assert(aln_1.Name().compare(aln_2.Name()) == 0);
    min_read_start_ = std::min(aln_1.Position(), aln_2.Position());
    max_read_start_ = std::max(aln_1.Position(), aln_2.Position());
  }

  BamAlignment& aln_one(){ return aln_1_; }
  BamAlignment& aln_two(){ return aln_2_; }
  const std::string& name() const { return name_;  }
  bool single_ended()       const { return min_read_start_ == -1; }

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
    if (max_read_start_ != pair.max_read_start_)
      return max_read_start_ < pair.max_read_start_;
    return (name_.compare(pair.name_) < 0);
  }
};

void remove_pcr_duplicates(const BaseQuality& base_quality, bool use_bam_rgs,
			   const std::map<std::string, std::string>& rg_to_library,
			   std::vector< std::vector<BamAlignment> >& paired_strs_by_rg,
			   std::vector< std::vector<BamAlignment> >& mate_pairs_by_rg,
			   std::vector< std::vector<BamAlignment> >& unpaired_strs_by_rg, std::ostream& logger);

#endif
