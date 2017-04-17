#ifndef PROCESS_TIMER_H_
#define PROCESS_TIMER_H_

#include <map>
#include <string>

class ProcessTimer {
 private:
  std::map<std::string, double> total_times_;

 public:
  ProcessTimer(){}
  
  void add_time(std::string key, double time){
    if (total_times_.find(key) == total_times_.end())
      total_times_[key] = 0.0;
    total_times_[key] += time;
  }

  double get_total_time(std::string key) const {
    auto iter = total_times_.find(key);
    if (iter == total_times_.end())
      return 0.0;
    return iter->second;
  }
};

#endif
