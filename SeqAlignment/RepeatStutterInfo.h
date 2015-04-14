#ifndef REPEAT_STUTTER_INFO_H_
#define REPEAT_STUTTER_INFO_H_

#include <algorithm>
#include <string>

#include "../error.h"

// Maximum number of repeats that PCR stutter can add or remove from block's sequence
const int MAX_STUTTER_REPEAT_INS = 3;
const int MAX_STUTTER_REPEAT_DEL = -3;

class RepeatStutterInfo {
 private:
  int period_;
  int max_ins_, max_del_;

 public:
  RepeatStutterInfo(int period, std::string& ref_allele ){
    period_  = period;
    max_ins_ = MAX_STUTTER_REPEAT_INS*period_;
    max_del_ = std::max(MAX_STUTTER_REPEAT_DEL*period_, -(int)ref_allele.size());
  }
  
  inline const int get_period()         const  { return period_;  }
  inline const int max_insertion()      const  { return max_ins_; }
  inline const int max_deletion()       const  { return max_del_; }

  void add_alternate_allele(std::string& alt_allele){
    max_del_ = std::max(max_del_, -(int)alt_allele.size());
  }

  /*
   * Returns the probability that PCR stutter deletes NUM_REPEATS from observed sequence
   * for sequence SEQ_INDEX (0 = reference allele)
   */
  inline double log_prob_stutter_deletion (int num_repeats, int seq_index) const {
    printErrorAndDie("log_prob_stutter_deletion() not implemented");
    return 0.0;
  }

  /*
   * Returns the probability that PCR stutter inserts NUM_REPEATS from observed sequence
   * for sequence SEQ_INDEX (0 = reference allele)
   */
  inline double log_prob_stutter_insertion (int num_repeats, int seq_index) const {
    printErrorAndDie("log_prob_stutter_insertion() not implemented");
    return 0.0;
  }

  /*
   * Returns the probability that PCR stutter does not modify the observed sequence
   * for sequence SEQ_INDEX (0 = reference allele)
   */
  inline double log_prob_no_stutter (int seq_index) const {
    printErrorAndDie("log_prob_no_stutter() not implemented");
    return 0.0;
  }

  inline double log_prob_pcr_artifact(int seq_index, int artifact_size) const {
    if (artifact_size == 0)
      return log_prob_no_stutter(seq_index);
    else if (artifact_size < 0)
      return log_prob_stutter_deletion(seq_index, artifact_size);
    else
      return log_prob_stutter_insertion(seq_index, artifact_size);
  }
};
#endif
