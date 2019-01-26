#ifndef REPEAT_STUTTER_INFO_H_
#define REPEAT_STUTTER_INFO_H_

#include <algorithm>
#include <string>

#include "../stutter_model.h"

// Maximum number of repeats that PCR stutter can add or remove from block's sequence
const int    MAX_STUTTER_REPEAT_INS = 6;
const int    MAX_STUTTER_REPEAT_DEL = -6;
const double LARGE_NEGATIVE         = -10e6;

class RepeatStutterInfo {
 private:
  const int period_, max_ins_, max_del_;
  StutterModel* stutter_model_;
  std::vector<int> allele_sizes_;
  std::vector<double*> stutter_probs_;

  // Private unimplemented copy constructor and assignment operator to prevent operations
  RepeatStutterInfo(const RepeatStutterInfo& other);
  RepeatStutterInfo& operator=(const RepeatStutterInfo& other);

  void add_stutter_parameters(int allele_size){
    double* probs = new double[max_ins_-max_del_+1];
    int idx       = 0;
    for (int artifact_size = max_del_; artifact_size <= max_ins_; ++artifact_size)
      probs[idx++] = (allele_size + artifact_size < 0 ? LARGE_NEGATIVE :
		      stutter_model_->log_stutter_pmf(allele_size, allele_size + artifact_size));
    stutter_probs_.push_back(probs);
  }

 public:
 RepeatStutterInfo(int period, const std::string& ref_allele, const StutterModel* stutter_model) :
  period_(period), max_ins_(MAX_STUTTER_REPEAT_INS*period_), max_del_(MAX_STUTTER_REPEAT_DEL*period_){
    assert(stutter_model != NULL && period > 0);
    stutter_model_ = stutter_model->copy();
    allele_sizes_.push_back(ref_allele.size());
    add_stutter_parameters(allele_sizes_.back());
  }

  ~RepeatStutterInfo(){
    delete stutter_model_;
    for (unsigned int i = 0; i < stutter_probs_.size(); ++i)
      delete [] stutter_probs_[i];
    stutter_probs_.clear();
  }

  void set_stutter_model(const StutterModel* model){
    assert(model != NULL);
    delete stutter_model_;
    stutter_model_ = model->copy();

    for (unsigned int i = 0; i < stutter_probs_.size(); ++i)
      delete [] stutter_probs_[i];
    stutter_probs_.clear();

    for (unsigned int i = 0; i < allele_sizes_.size(); ++i)
      add_stutter_parameters(allele_sizes_[i]);
  }

  inline StutterModel* get_stutter_model() const { return stutter_model_; }
  inline int get_period()                  const { return period_;        }
  inline int max_insertion()               const { return max_ins_;       }
  inline int max_deletion()                const { return max_del_;       }

  void add_alternate_allele(const std::string& alt_allele){
    allele_sizes_.push_back((int)alt_allele.size());
    add_stutter_parameters(allele_sizes_.back());
  }

  inline double log_prob_pcr_artifact(int seq_index, int artifact_size) const {
    if (artifact_size < max_del_ || artifact_size > max_ins_)
      return LARGE_NEGATIVE;
    return stutter_probs_[seq_index][artifact_size-max_del_];
  }
};

#endif
