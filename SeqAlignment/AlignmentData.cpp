#include <algorithm>
#include <vector>

#include "AlignmentData.h"
#include "error.h"

bool compareAl(const Alignment& alignment_1, const Alignment& alignment_2){
  int sample_comp = alignment_1.get_sample().compare(alignment_2.get_sample());
  if (sample_comp != 0)
    return sample_comp < 0;
  if (alignment_1.get_start() != alignment_2.get_start())
    return alignment_1.get_start() < alignment_2.get_start();
  return alignment_1.get_stop() < alignment_2.get_stop();
}

void sortAlignments(std::vector<Alignment>& alignments){
  std::sort(alignments.begin(), alignments.end(), compareAl);
}

void Alignment::get_deletion_boundaries(std::vector<int32_t>& starts, std::vector<int32_t>& stops) const {
  int32_t pos = start_;
  for (auto iter = cigar_list_.begin(); iter != cigar_list_.end(); iter++){
    switch(iter->get_type()){
    case 'M': case 'X': case '=':
      pos += iter->get_num();
      break;
    case 'I': case 'S':
      break;
    case 'D':
      starts.push_back(pos);
      stops.push_back(pos+iter->get_num()-1);
      pos += iter->get_num();
      break;
    default:
      printErrorAndDie("Invalid CIGAR char detected in get_deletion_boundaries for alignment with CIGAR " + getCigarString() + "and alignment " + alignment_);
    }
  }
}

void Alignment::get_insertion_positions(std::vector<int32_t>& positions, std::vector<int32_t>& sizes) const {
  int32_t pos = start_;
  for (auto iter = cigar_list_.begin(); iter != cigar_list_.end(); iter++){
    switch(iter->get_type()){
    case 'M': case 'X': case '=':
      pos += iter->get_num();
      break;
    case 'S':
      break;
    case 'I':
      positions.push_back(pos);
      sizes.push_back(iter->get_num());
      break;
    case 'D':
      pos += iter->get_num();
      break;
    default:
      printErrorAndDie("Invalid CIGAR char detected in get_insertion_positions");
    }
  }
}
