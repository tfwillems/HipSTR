#ifndef REPEAT_REGION_H_
#define REPEAT_REGION_H_

#include <iostream>
#include <string>
#include <vector>

class RepeatRegion {
 private:
  int32_t start_; // Inclusive
  int32_t stop_;  // Not inclusive
  int period_;
  std::string sequence_;

 public:
  RepeatRegion(int32_t start, int32_t stop, int period, std::string sequence){
    start_    = start;
    stop_     = stop;
    period_   = period;
    sequence_ = sequence;
  }

  inline int32_t get_start()               const{ return start_;    }
  inline int32_t get_stop()                const{ return stop_;     }
  inline int     get_period()              const{ return period_;   }
  inline const std::string& get_sequence() const{ return sequence_; }

  void print(std::ostream& out){
    out << start_ << "\t" << stop_ << "\t" << period_ << "\t" << sequence_ << std::endl; 
  }
};

bool compareRepeatRegions(const RepeatRegion& r1, const RepeatRegion& r2);

void sortRepeatRegions(std::vector<RepeatRegion>& repeat_regions);

#endif
