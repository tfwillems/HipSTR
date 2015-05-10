#ifndef MATH_OPS_H_
#define MATH_OPS_H_

#define cast_uint32_t static_cast<uint32_t>

#include <algorithm>
#include <vector>

extern const double LOG_ONE_HALF;
extern const double TOLERANCE;


double sum(double* begin, double* end);

double sum(std::vector<double>& vals);

int sum (std::vector<bool>& vals);

double log_sum_exp(double* begin, double* end);

double log_sum_exp(double log_v1, double log_v2);

double log_sum_exp(double log_v1, double log_v2, double log_v3);

double log_sum_exp(std::vector<double>& log_vals);

double expected_value(double* log_probs, std::vector<int>& vals); 

// To accelerate logsumexp, ignore values if they're 1/1000th or less than the maximum value
const double LOG_THRESH   = log(0.001);

static inline float fasterlog (float x){
  union { float f; uint32_t i; } vx = { x };
  float y = vx.i;
  y *= 8.2629582881927490e-8f;
  return y - 87.989971088f;
}

static inline float fasterpow2 (float p){
  float clipp = (p < -126) ? -126.0f : p;
  union { uint32_t i; float f; } v = { cast_uint32_t ( (1 << 23) * (clipp + 126.94269504f) ) };
  return v.f;
}

static inline float fasterexp (float p) {
  return fasterpow2 (1.442695040f * p);
}

static inline double fast_log_sum_exp(double log_v1, double log_v2){
  if (log_v1 > log_v2){
    double diff = log_v2-log_v1;
    return diff < LOG_THRESH ? log_v1 : log_v1 + fasterlog(1 + fasterexp(diff));
  }
  else {
    double diff = log_v1-log_v2;
    return diff < LOG_THRESH ? log_v2 : log_v2 + fasterlog(1 + fasterexp(diff));
  }       
}

static inline double fast_log_sum_exp(std::vector<double>& log_vals){
  double max_val = *std::max_element(log_vals.begin(), log_vals.end());
  double total   = 0;
  for (auto iter = log_vals.begin(); iter != log_vals.end(); iter++){
    double diff = *iter - max_val;
    if (diff > LOG_THRESH)
      total += fasterexp(diff);
  }
  return max_val + fasterlog(total);
}

#endif
