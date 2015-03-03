#include "stutter_model.h"

std::ostream& operator<< (std::ostream &out, StutterModel& model){
  out << "(Stutter Model: P_GEOM=" << model.p_geom_ << ", P_DOWN=" << model.p_down_ << ", P_UP=" << model.p_up_ << std::endl; 
  return out;
}

