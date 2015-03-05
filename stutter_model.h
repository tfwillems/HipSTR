#ifndef STUTTER_MODEL_H_
#define STUTTER_MODEL_H_

#include <assert.h>
#include <math.h>
#include <iostream>

class StutterModel {
 private:
  double p_geom_, p_up_, p_down_;
  double log_nostep_, log_step_;
  double log_up_, log_down_, log_equal_;
  
 public:
  StutterModel(double p_geom, double p_up, double p_down){
    assert(p_geom > 0 && p_geom < 1.0);
    assert(p_up > 0.0 && p_down > 0.0 && p_up+p_down < 1.0);
    p_geom_     = p_geom; 
    log_step_   = log(1-p_geom); 
    log_nostep_ = log(p_geom);
    p_up_       = p_up; 
    log_up_     = log(p_up);
    p_down_     = p_down;
    log_down_   = log(p_down);
    log_equal_  = log(1.0-p_up-p_down);
 }

  friend std::ostream& operator<< (std::ostream &out, StutterModel& model);

  inline double calc_log_stutter(int sample_repeats, int read_repeats){
    if (read_repeats == sample_repeats)
      return log_equal_;
    if (read_repeats < sample_repeats)
      return log_down_ + log_nostep_ + log_step_*(sample_repeats-read_repeats);
    else
      return log_up_   + log_nostep_ + log_step_*(read_repeats-sample_repeats);
  }
};

#endif
