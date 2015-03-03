#ifndef BASE_QUALITY_H_
#define BASE_QUALITY_H_

#include <iostream>
#include <string>

#include "error.h"
#include "math.h"

class BaseQuality {
 private:
  // Based on the Illumina 1.8 Phred+33 system
  const static char MIN_BASE_QUALITY='!';
  const static char MAX_BASE_QUALITY='J';

  // Log likelihoods for quality scores
  double log_correct_[256];
  double log_error_[256];

 public:
  BaseQuality(){
    // Precalculate log likelihoods
    log_correct_[MIN_BASE_QUALITY] = -100000;
    log_error_[MIN_BASE_QUALITY]   = 0;

    int count = 1;
    for (int i = MIN_BASE_QUALITY+1; i <= MAX_BASE_QUALITY; ++i, ++count){
      log_correct_[i] = log(1.0 - pow(10.0, count/(-10.0)));
      log_error_[i]   = log(pow(10.0, count/(-10.0))/3.0);
    }
  }

  /*
   * Returns the log likelihood that the base with the
   * provided quality score should match a different base
   */
  double log_prob_error(char quality) {
    if (quality < MIN_BASE_QUALITY){
      std::cerr << "WARNING: Base quality " << quality << " outside of expected range. Proceeding using minimum expected score " << std::endl; 
      return log_error_[MIN_BASE_QUALITY];
    }
    else if (quality > MAX_BASE_QUALITY) {
      std::cerr << "WARNING: Base quality " << quality << " outside of expected range. Proceeding using maximum expected score " << std::endl;
      return log_error_[MAX_BASE_QUALITY];
    }
    else
      return log_error_[quality];
  }


  /*
   * Returns the log likelihood that the base with the
   * provided quality score was observed without error
   */  
  double log_prob_correct(char quality){
    if (quality < MIN_BASE_QUALITY){
      std::cerr << "WARNING: Base quality " << quality << " outside of expected range. Proceeding using minimum expected score " << std::endl; 
      return log_correct_[MIN_BASE_QUALITY];
    }
    else if (quality > MAX_BASE_QUALITY) {
      std::cerr << "WARNING: Base quality " << quality << " outside of expected range. Proceeding using maximum expected score " << std::endl;
      return log_correct_[MAX_BASE_QUALITY];
    }
    else
      return log_correct_[quality];
  }  
};

#endif
