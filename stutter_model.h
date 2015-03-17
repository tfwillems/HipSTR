#ifndef STUTTER_MODEL_H_
#define STUTTER_MODEL_H_

#include <assert.h>
#include <math.h>
#include <iostream>

class StutterModel {
 private:
  // Parameters for in-frame mutations. Geometric parameter in units of repeats
  double in_geom_, in_up_, in_down_;
  double in_log_nostep_, in_log_step_;
  double in_log_up_, in_log_down_, in_log_equal_;

  // Parameters for out-of-frame mutations. Geometric parameter in units of base pairs
  double out_geom_, out_up_, out_down_;
  double out_log_nostep_, out_log_step_;
  double out_log_up_, out_log_down_;
  
  int motif_len_;
  
 public:
  StutterModel(double inframe_geom,  double inframe_up,  double inframe_down, 
	       double outframe_geom, double outframe_up, double outframe_down, int motif_len){
    assert(inframe_geom > 0 && inframe_geom < 1.0);
    assert(inframe_up > 0.0 && inframe_down > 0.0 && inframe_up+inframe_down < 1.0);
    assert(outframe_geom > 0 && outframe_geom < 1.0);
    assert(outframe_up > 0.0 && outframe_down > 0.0 && outframe_up+outframe_down < 1.0);

    in_geom_        = inframe_geom;
    in_log_step_    = log(1-inframe_geom);
    in_log_nostep_  = log(inframe_geom);
    in_up_          = inframe_up;
    in_log_up_      = log(inframe_up);
    in_down_        = inframe_down;
    in_log_down_    = log(inframe_down);
    in_log_equal_   = log(1-inframe_up-inframe_down);

    out_geom_       = outframe_geom;
    out_log_step_   = log(1-outframe_geom);
    out_log_nostep_ = log(outframe_geom);
    out_up_         = outframe_up;
    out_log_up_     = log(outframe_up);
    out_down_       = outframe_down;
    out_log_down_   = log(outframe_down);

    motif_len_      = motif_len;
 }

  friend std::ostream& operator<< (std::ostream &out, StutterModel& model);

  double calc_log_stutter(int sample_bps, int read_bps){
    int bp_diff = read_bps - sample_bps;
    if (bp_diff % motif_len_ != 0){
      int eff_diff = bp_diff - (bp_diff/motif_len_);
      if (eff_diff < 0)
	return out_log_down_ + out_log_nostep_ + out_log_step_*(-eff_diff-1);
      else
	return out_log_up_   + out_log_nostep_ + out_log_step_*(eff_diff-1);
    }
    else {
      int rep_diff = bp_diff/motif_len_;
      if (rep_diff == 0)
	return in_log_equal_;
      if (rep_diff < 0)
	return in_log_down_ + in_log_nostep_ + in_log_step_*(-rep_diff-1);
      else
	return in_log_up_   + in_log_nostep_ + in_log_step_*(rep_diff-1);
    }
  }
};

#endif
