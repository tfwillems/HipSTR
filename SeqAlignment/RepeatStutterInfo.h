#ifndef REPEAT_STUTTER_INFO_H_
#define REPEAT_STUTTER_INFO_H_

#include <algorithm>
#include <string>

#include "../error.h"
#include "../stutter_model.h"

// Maximum number of repeats that PCR stutter can add or remove from block's sequence
const int    MAX_STUTTER_REPEAT_INS = 3;
const int    MAX_STUTTER_REPEAT_DEL = -3;
const double LARGE_NEGATIVE = -10e6;

class RepeatStutterInfo {
 private:
  int period_;
  int max_ins_, max_del_;
  StutterModel* stutter_model_;
  std::vector<int> allele_sizes_;

 public:
  RepeatStutterInfo(int period, std::string& ref_allele, StutterModel* stutter_model){
    period_        = period;
    max_ins_       = MAX_STUTTER_REPEAT_INS*period_;
    max_del_       = MAX_STUTTER_REPEAT_DEL*period_;
    stutter_model_ = stutter_model;
    allele_sizes_.push_back(ref_allele.size());
  }
  
  inline const int get_period()         const  { return period_;  }
  inline const int max_insertion()      const  { return max_ins_; }
  inline const int max_deletion()       const  { return max_del_; }

  void add_alternate_allele(std::string& alt_allele){
    allele_sizes_.push_back((int)alt_allele.size());
  }

  inline double log_prob_pcr_artifact(int seq_index, int artifact_size) const {
    int read_size = allele_sizes_[seq_index]+artifact_size;
    if (artifact_size == 0)
      return stutter_model_->log_stutter_pmf(allele_sizes_[seq_index], read_size);
    else if (artifact_size > 0)
      return (artifact_size > max_ins_ ? LARGE_NEGATIVE : stutter_model_->log_stutter_pmf(allele_sizes_[seq_index], read_size));
    else
      return (artifact_size < max_del_ || read_size < 0 ? LARGE_NEGATIVE : stutter_model_->log_stutter_pmf(allele_sizes_[seq_index], read_size));
  }
};
#endif
