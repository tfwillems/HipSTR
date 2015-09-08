#ifndef INSERT_SIZE_COUNTER_H_
#define INSERT_SIZE_COUNTER_H_

#include <algorithm>
#include <iostream>
#include <vector>

#include "bamtools/include/api/BamAlignment.h"

constexpr int64_t bin_boundaries[21] = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1250, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 5000, 10000};

class InsertSizeCounter {
 private:
  //int64_t bin_boundaries[21];
  int64_t total_comps_, valid_comps_;
  double sum_diffs_, sum_square_diffs_;
  std::vector<int64_t> bins_, counts_;
  int64_t max_count_, diff_chrom_count_, unmapped_count_, max_diff_;

  bool add_diff(int abs_diff){
    int count_index = std::lower_bound(bins_.begin(), bins_.end(), abs_diff) - bins_.begin();
    if (count_index == 0)
      return false;
    counts_[count_index-1]++;
    return true;
  }

  /*
  void init_bins(){
    int index = 0;
    // 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000
    for (unsigned int i = 0; i < 11; i++)
      bin_boundaries[index++] = i*100;
    // 1250, 1500, 1750, 2000
    for (unsigned int i = 0; i < 4; i++)
      bin_boundaries[index++] = 1250 + i*250;
    // 2500, 3000, 3500, 4000
    for (unsigned int i = 0; i < 4; i++)
      bin_boundaries[index++] = 2500 + i*500;
    bin_boundaries[index++] = 5000;
    bin_boundaries[index++] = 10000;
  }
  */

 public:
  InsertSizeCounter(int64_t max_diff){
    total_comps_      = 0;
    valid_comps_      = 0;
    sum_diffs_        = 0;
    sum_square_diffs_ = 0;
    max_count_        = 0;
    diff_chrom_count_ = 0;
    unmapped_count_   = 0;
    max_diff_         = max_diff;
    counts_           = std::vector<int64_t>(20, 0);
    bins_             = std::vector<int64_t>(bin_boundaries, bin_boundaries+21);
  }

  void process_alignment(BamTools::BamAlignment& aln);

  void output_summary_statistics(std::ostream& out);
};

#endif
