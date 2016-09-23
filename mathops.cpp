#include <algorithm>
#include <assert.h>
#include <math.h>

#include "mathops.h"

const double LOG_ONE_HALF  = log(0.5);
const double TOLERANCE     = 1e-10;
const double LOG_E_BASE_10 = 0.4342944819;

double INT_LOGS[10000];

void precompute_integer_logs(){
  INT_LOGS[0] = -1000;
  for (unsigned int i = 1; i < 10000; i++)
    INT_LOGS[i] = log(i);
}

double int_log(int val){ return INT_LOGS[val]; }

double sum(double* begin, double* end){
  double total = 0.0;
  for (double* iter = begin; iter != end; iter++)
    total += *iter;
  return total;
}

double sum(std::vector<double>& vals){
  double total = 0.0;
  for (auto iter = vals.begin(); iter != vals.end(); iter++)
    total += *iter;
  return total;
}

int sum(std::vector<bool>& vals){
  int total = 0;
  for (auto iter = vals.begin(); iter != vals.end(); iter++)
    total += *iter;
  return total;
}

double log_sum_exp(double* begin, double* end){
  double max_val = *std::max_element(begin, end);
  double total   = 0.0;
  for (double* iter = begin; iter != end; iter++)
    total += exp(*iter - max_val);
  return max_val + log(total);
}

double log_sum_exp(double log_v1, double log_v2){
  if (log_v1 > log_v2)
    return log_v1 + log(1 + exp(log_v2-log_v1));
  else
    return log_v2 + log(1 + exp(log_v1-log_v2));
}

double log_sum_exp(double log_v1, double log_v2, double log_v3){
  double max_val = std::max(std::max(log_v1, log_v2), log_v3);
  return max_val + log(exp(log_v1-max_val) + exp(log_v2-max_val) + exp(log_v3-max_val));
}

double log_sum_exp(std::vector<double>& log_vals){
  double max_val = *std::max_element(log_vals.begin(), log_vals.end());
  double total   = 0;
  for (auto iter = log_vals.begin(); iter != log_vals.end(); iter++)
    total += exp(*iter - max_val);
  return max_val + log(total);
}

double expected_value(double* log_likelihoods, std::vector<int>& vals){
  double total_LL = log_sum_exp(log_likelihoods, log_likelihoods+vals.size());
  double total    = 0.0;
  for (unsigned int i = 0; i < vals.size(); i++)
    total += exp(log_likelihoods[i]-total_LL)*vals[i];
  return total;
}

double expected_value(std::vector<double>& log_likelihoods, std::vector<int>& vals){
  assert(log_likelihoods.size() == vals.size());
  double total_LL = log_sum_exp(log_likelihoods);
  double total    = 0.0;
  for (unsigned int i = 0; i < vals.size(); i++)
    total += exp(log_likelihoods[i]-total_LL)*vals[i];
  return total;
}


void fast_update_streaming_log_sum_exp(double log_val, double& max_val, double& total){
  if (log_val <= max_val)
    total += fasterexp(log_val - max_val);
  else {
    total  *= fasterexp(max_val-log_val);
    total  += 1.0;
    max_val = log_val;
  }
}

double fast_finish_streaming_log_sum_exp(double max_val, double total){
  return max_val + fasterlog(total);
}

void update_streaming_log_sum_exp(double log_val, double& max_val, double& total){
  if (log_val <= max_val)
    total += exp(log_val - max_val);
  else {
    total  *= exp(max_val-log_val);
    total  += 1.0;
    max_val = log_val;
  }
}

double finish_streaming_log_sum_exp(double max_val, double total){
  return max_val + log(total);
}
