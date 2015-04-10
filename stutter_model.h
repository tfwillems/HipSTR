#ifndef STUTTER_MODEL_H_
#define STUTTER_MODEL_H_

#include <algorithm>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <vector>

#include "error.h"

class StutterModel {
 private:
  // Parameters for in-frame mutations. Geometric parameter in units of repeats
  double in_geom_, in_up_, in_down_;
  double in_log_nostep_, in_log_step_;
  double in_log_up_, in_log_down_;

  // Parameter for no stutter change
  double log_equal_;

  // Parameters for out-of-frame mutations. Geometric parameter in units of base pairs
  double out_geom_, out_up_, out_down_;
  double out_log_nostep_, out_log_step_;
  double out_log_up_, out_log_down_;
  
  int motif_len_;

  inline double log_sum_exp(std::vector<double>& log_vals){
    double max_val = *std::max_element(log_vals.begin(), log_vals.end());
    double total   = 0;
    for (auto iter = log_vals.begin(); iter != log_vals.end(); iter++)
      total += exp(*iter - max_val);
    return max_val + log(total);
  }

  inline double log_geom_geq(double p, int n){
    return (n-1)*log(1-p);
  }
  
  inline double log_geom_leq(double p, int n){
    return log(1-pow(1-p,n));
  }

 public:
  StutterModel(double inframe_geom,  double inframe_up,  double inframe_down, 
	       double outframe_geom, double outframe_up, double outframe_down, int motif_len){
    assert(inframe_geom  > 0.0 && inframe_geom  < 1.0);
    assert(outframe_geom > 0.0 && outframe_geom < 1.0);
    assert(inframe_up    > 0.0 && inframe_down  > 0.0);
    assert(outframe_up   > 0.0 && outframe_down > 0.0);
    assert(inframe_up + inframe_down + outframe_up + outframe_down < 1.0);

    in_geom_        = inframe_geom;
    in_log_step_    = log(1-inframe_geom);
    in_log_nostep_  = log(inframe_geom);
    in_up_          = inframe_up;
    in_log_up_      = log(inframe_up);
    in_down_        = inframe_down;
    in_log_down_    = log(inframe_down);

    out_geom_       = outframe_geom;
    out_log_step_   = log(1-outframe_geom);
    out_log_nostep_ = log(outframe_geom);
    out_up_         = outframe_up;
    out_log_up_     = log(outframe_up);
    out_down_       = outframe_down;
    out_log_down_   = log(outframe_down);

    log_equal_      = log(1-inframe_up-inframe_down-outframe_up-outframe_down);
    motif_len_      = motif_len;
 }

  friend std::ostream& operator<< (std::ostream &out, StutterModel& model);

  double get_parameter(bool in_frame, char parameter){
    switch(parameter){
    case 'U':
      return in_frame ? in_up_   : out_up_;
    case 'D':
      return in_frame ? in_down_ : out_down_;
    case 'P':
      return in_frame ? in_geom_ : out_geom_;
    default:
      printErrorAndDie("Invalid such stutter model parameter requested");
      return -1.0;
    }
  }


  /*
    Returns the read's log-likeliood given that it contains exactly the provided number of base pairs
  */
  double log_stutter_pmf(int sample_bps, int read_bps){
    double log_pmf;
    int bp_diff = read_bps - sample_bps;
    if (bp_diff % motif_len_ != 0){
      int eff_diff = bp_diff - (bp_diff/motif_len_);
      if (eff_diff < 0)
	log_pmf = out_log_down_ + out_log_nostep_ + out_log_step_*(-eff_diff-1);
      else
	log_pmf = out_log_up_   + out_log_nostep_ + out_log_step_*(eff_diff-1);
    }
    else {
      int rep_diff = bp_diff/motif_len_;
      if (rep_diff == 0)
	log_pmf = log_equal_;
      else {
	if (rep_diff < 0)
	  log_pmf = in_log_down_ + in_log_nostep_ + in_log_step_*(-rep_diff-1);
	else
	  log_pmf = in_log_up_   + in_log_nostep_ + in_log_step_*(rep_diff-1);
      }
    }
    assert(log_pmf <= 0);
    return log_pmf;
  }


  /*
     Returns the read's log-likelihood given that it contains at least the provided number of base pairs
   */
  double log_stutter_geq(int sample_bps, int min_read_bps){
    std::vector<double> log_probs;
    int min_bp_diff = min_read_bps - sample_bps;

    // Incorporate all potential in-frame stutters
    int next_rep_diff = (min_bp_diff < 0 || min_bp_diff % motif_len_ == 0 ? min_bp_diff/motif_len_:  1 + min_bp_diff/motif_len_);
    if (next_rep_diff < 0){
      log_probs.push_back(in_log_down_ + log_geom_leq(in_geom_, -next_rep_diff));
      log_probs.push_back(log_equal_);
      log_probs.push_back(in_log_up_);    
    }
    else if (next_rep_diff == 0){
      log_probs.push_back(log_equal_);
      log_probs.push_back(in_log_up_);
    }
    else
      log_probs.push_back(in_log_up_ + log_geom_geq(in_geom_, next_rep_diff));
    
    // Incorporate all potential out-of-frame stutters
    int next_outframe_diff = min_bp_diff + (min_bp_diff % motif_len_ == 0);
    int eff_diff = next_outframe_diff - (next_outframe_diff/motif_len_);
    if (eff_diff < 0){
      log_probs.push_back(out_log_down_ + log_geom_leq(out_geom_, -eff_diff));
      log_probs.push_back(out_log_up_);
    }
    else
      log_probs.push_back(out_log_up_ + log_geom_geq(out_geom_, eff_diff));

    return log_sum_exp(log_probs);
  }
};

#endif
