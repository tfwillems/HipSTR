#include <algorithm>
#include <assert.h>
#include <map>
#include <sstream>

#include "base_quality.h"
#include "error.h"
#include "mathops.h"
#include "stringops.h"

std::string BaseQuality::median_base_qualities(const std::vector<const std::string*>& qualities){
  assert(qualities.size() > 0);

  // Check that all base quality strings are of the same length
  for (unsigned int i = 0; i < qualities.size(); i++)
    if (qualities[i]->size() != qualities[0]->size())
      printErrorAndDie("All base quality strings must be of the same length when averaging probabilities");

  std::string median_qualities(qualities[0]->size(), 'N');
  for (unsigned int i = 0; i < qualities[0]->size(); i++){
    std::vector<char> quals;
    for (unsigned int j = 0; j < qualities.size(); j++)
      quals.push_back(qualities[j]->at(i));
    std::sort(quals.begin(), quals.end());
    median_qualities[i] = quals[quals.size()/2];
  }
  return median_qualities;
}

void printBaseCounts(int* counts, std::ostream& out){
  for (unsigned int i = 0; i < 256; i++)
    if (counts[i] != 0)
      out << (char)i << " " << counts[i] << std::endl;
}
