#include "RepeatBlock.h"

#include "../stringops.h"

void RepeatBlock::add_alternate(const std::string& alt){
  HapBlock::add_alternate(alt);
  repeat_info_->add_alternate_allele(alt);
  stutter_aligners_.push_back(new StutterAligner(alt, repeat_info_->get_period(), !reversed_, repeat_info_));

  const int num_alleles = num_options();
  suffix_match_lengths_.push_back(std::vector<int>(num_alleles));
  for (int i = 0; i < num_alleles-1; ++i){
    int suffix_len = length_suffix_match(get_seq(i), alt);
    suffix_match_lengths_[i].push_back(suffix_len);
    suffix_match_lengths_[num_alleles-1][i] = suffix_len;
  }

  suffix_match_lengths_[num_alleles-1][num_alleles-1] = 0;
}
