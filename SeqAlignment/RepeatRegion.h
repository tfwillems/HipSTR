#ifndef REPEAT_REGION_H_
#define REPEAT_REGION_H_

#include <iostream>
#include <string>
#include <vector>

class RepeatRegion {
 private:
  int32_t start_; // Inclusive
  int32_t stop_;  // Not inclusive
  std::string motif_;
  std::string sequence_;

 public:
  RepeatRegion(int32_t start, int32_t stop, std::string motif, std::string sequence){
    start_    = start;
    stop_     = stop;
    motif_    = motif;
    sequence_ = sequence;
  }

  inline int32_t get_start()               const{ return start_;    }
  inline int32_t get_stop()                const{ return stop_;     }
  inline const std::string& get_motif()    const{ return motif_;    }
  inline const std::string& get_sequence() const{ return sequence_; }

  void print(std::ostream& out){
    out << start_ << "\t" << stop_ << "\t" << motif_ << "\t" << sequence_ << std::endl; 
  }
};


void get_repeat_regions(std::string& seq, int32_t seq_start, std::vector<RepeatRegion>& regions);

bool compareRepeatRegions(const RepeatRegion& r1, const RepeatRegion& r2);

void sortRepeatRegions(std::vector<RepeatRegion>& repeat_regions);

#endif
