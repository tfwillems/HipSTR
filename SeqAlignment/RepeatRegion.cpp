#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

#include "RepeatRegion.h"


bool compareRepeatRegions(const RepeatRegion& r1, const RepeatRegion& r2){
  if (r1.get_start() != r2.get_start())
    return r1.get_start() < r2.get_start();
  else
    return r1.get_stop()  > r2.get_stop();
}

void sortRepeatRegions(std::vector<RepeatRegion>& repeat_regions){
  std::sort(repeat_regions.begin(), repeat_regions.end(), compareRepeatRegions);
}
