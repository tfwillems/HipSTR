#ifndef MATH_OPS_H_
#define MATH_OPS_H_

#define cast_uint32_t static_cast<uint32_t>

#include <vector>

extern const double LOG_ONE_HALF;

double sum(double* begin, double* end);

double sum(std::vector<double>& vals);

double log_sum_exp(double* begin, double* end);

double log_sum_exp(double log_v1, double log_v2);

double log_sum_exp(double log_v1, double log_v2, double log_v3);

double log_sum_exp(std::vector<double>& log_vals);


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


#endif
