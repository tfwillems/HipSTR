#include "RepeatBlock.h"

void RepeatBlock::add_alternate(const std::string& alt){
  HapBlock::add_alternate(alt);
  repeat_info_->add_alternate_allele(alt);
  stutter_aligners_.push_back(new StutterAligner(alt, repeat_info_->get_period(), !reversed_, repeat_info_));
}
