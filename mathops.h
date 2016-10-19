#ifndef MATH_OPS_H_
#define MATH_OPS_H_

#include <math.h>
#include <cstdint>
#include <algorithm>
#include <vector>

extern const double LOG_ONE_HALF;
extern const double TOLERANCE;
extern const double LOG_E_BASE_10;

void precompute_integer_logs();

double int_log(int val);

double sum(double* begin, double* end);

double sum(std::vector<double>& vals);

int sum (std::vector<bool>& vals);

double log_sum_exp(double* begin, double* end);

double log_sum_exp(double log_v1, double log_v2);

double log_sum_exp(double log_v1, double log_v2, double log_v3);

double log_sum_exp(std::vector<double>& log_vals);

double expected_value(double* log_probs, std::vector<int>& vals); 

double expected_value(std::vector<double>& log_likelihoods, std::vector<int>& vals);

void update_streaming_log_sum_exp(double log_val, double& max_val, double& total);

double finish_streaming_log_sum_exp(double max_val, double total);

void fast_update_streaming_log_sum_exp(double log_val, double& max_val, double& total);

double fast_finish_streaming_log_sum_exp(double max_val, double total);

// To accelerate logsumexp, ignore values if they're 1/1000th or less than the maximum value
const double LOG_THRESH   = log(0.001);

double fast_log_sum_exp(double log_v1, double log_v2);
double fast_log_sum_exp(std::vector<double>& log_vals);

#endif
