#include "stutter_model.h"

std::ostream& operator<< (std::ostream &out, StutterModel& model){
  out << "(Stutter Model:"  << "\n"
      << "IN_FRAME [P_GEOM(rep) =" << model.in_geom_  << ", P_DOWN =" << model.in_down_  << ", P_UP =" << model.in_up_  << "]" << "\n"
      << "OUT_FRAME[P_GEOM(bp)  =" << model.out_geom_ << ", P_DOWN =" << model.out_down_ << ", P_UP =" << model.out_up_ << "]" << std::endl;
  return out;
}

