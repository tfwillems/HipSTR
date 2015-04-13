#ifndef MATH_OPS_H_
#define MATH_OPS_H_

#include <vector>

extern const double LOG_ONE_HALF;

double sum(double* begin, double* end);

double sum(std::vector<double>& vals);

double log_sum_exp(double* begin, double* end);

double log_sum_exp(double log_v1, double log_v2);

double log_sum_exp(double log_v1, double log_v2, double log_v3);

double log_sum_exp(std::vector<double>& log_vals);

#endif
