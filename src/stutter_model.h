#ifndef STUTTER_MODEL_H_
#define STUTTER_MODEL_H_

#include <algorithm>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <map>
#include <vector>

#include "error.h"
#include "mathops.h"
#include "region.h"

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

 public:
  StutterModel(double inframe_geom,  double inframe_up,  double inframe_down, 
	       double outframe_geom, double outframe_up, double outframe_down, int motif_len){
    assert(inframe_geom  > 0.0 && inframe_geom  < 1.0);
    assert(outframe_geom > 0.0 && outframe_geom < 1.0);
    assert(inframe_up    > 0.0 && inframe_down  > 0.0);
    assert(outframe_up   > 0.0 && outframe_down > 0.0);
    assert(inframe_up + inframe_down + outframe_up + outframe_down < 1.0);
    assert(motif_len > 0 && motif_len < 10);

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

  StutterModel* copy() const { return new StutterModel(in_geom_, in_up_, in_down_, out_geom_, out_up_, out_down_, motif_len_); }

  double get_parameter(bool in_frame, char parameter) const;

  double log_stutter_pmf(int sample_bps, int read_bps) const;

  int period() const { return motif_len_; }
  void set_period(int period){ motif_len_ = period; }

  void write(std::ostream& output) const;

  void write_model(const std::string& chrom, int32_t start, int32_t end, std::ostream& output) const;

  static StutterModel* read(std::istream& input);

  static void read_models(std::istream& input, std::map<Region, StutterModel*>& models);
};

#endif
