#include <sstream>

#include "stutter_model.h"

std::ostream& operator<< (std::ostream &out, StutterModel& model){
  out << "\n"
      << "\t" << "IN_FRAME [P_GEOM(rep)=" << model.in_geom_  << ", P_DOWN=" << model.in_down_  << ", P_UP=" << model.in_up_  << "]" << "\n"
      << "\t" << "OUT_FRAME[P_GEOM(bp) =" << model.out_geom_ << ", P_DOWN=" << model.out_down_ << ", P_UP=" << model.out_up_ << "]" << std::endl;
  return out;
}

double StutterModel::get_parameter(bool in_frame, char parameter) const {
  switch(parameter){
  case 'U':
    return in_frame ? in_up_   : out_up_;
  case 'D':
    return in_frame ? in_down_ : out_down_;
  case 'P':
      return in_frame ? in_geom_ : out_geom_;
  default:
    printErrorAndDie("Invalid such stutter model parameter requested");
  }
}


/*
  Returns the read's log-likeliood given that it contains exactly the provided number of base pairs
*/
double StutterModel::log_stutter_pmf(int sample_bps, int read_bps) const {
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

void StutterModel::write(std::ostream& output) const {
  output << in_geom_  << "\t" << in_down_  << "\t" << in_up_  << "\t" 
	 << out_geom_ << "\t" << out_down_ << "\t" << out_up_ << "\t" << motif_len_ << std::endl;
}

void StutterModel::write_model(const std::string& chrom, int32_t start, int32_t end, std::ostream& output) const {
  output << chrom << "\t" << start+1 << "\t" << end << "\t";
  write(output);
}

StutterModel* StutterModel::read(std::istream& input){
  double inframe_geom,  inframe_up,  inframe_down;
  double outframe_geom, outframe_up, outframe_down;
  int motif_len;
  std::string line;
  std::getline(input, line);
  std::istringstream ss(line);
  if (!(ss >> inframe_geom >> inframe_down >> inframe_up >> outframe_geom >> outframe_down >> outframe_up >> motif_len))
    printErrorAndDie("Improperly formatted stutter model file");
  if (motif_len < 1)
    printErrorAndDie("Improperly formatted stutter model file. One or more entries have a motif length < 1");
  if (motif_len > 9)
    printErrorAndDie("Improperly formatted stutter model file. One or more entries have a motif length > 9");
  return new StutterModel(inframe_geom, inframe_up, inframe_down, outframe_geom, outframe_up, outframe_down, motif_len);
}

void StutterModel::read_models(std::istream& input, std::map<Region, StutterModel*>& models){
  assert(models.size() == 0);
  while (input){
    std::string chrom;
    int32_t start, end;
    if (input >> chrom >> start >> end){
      StutterModel* model = read(input);
      models.insert(std::pair<Region,StutterModel*>(Region(chrom, start-1, end, model->period()), model));
    }
    else
      break;
  }
}

