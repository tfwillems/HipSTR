#include <algorithm>
#include <vector>

#include "AlignmentData.h"

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
