#include <assert.h>

#include "base_quality.h"
#include "error.h"
#include "mathops.h"

double BaseQuality::average_base_qualities(std::vector<std::string*> qualities, double* log_perror, double* log_pcorrect){
  assert(qualities.size() > 0);

  // Check that all base quality strings are of the same length
  for (unsigned int i = 0; i < qualities.size(); i++){
    if (qualities[i]->size() != qualities[0]->size())
      printErrorAndDie("All base quality strings must be of the same length when averaging probabilities");
  }

  // Average raw probabilities for each base
  std::vector<double> log_probs(qualities.size());
  for (unsigned int i = 0; i < qualities[0]->size(); i++){
    for (unsigned int j = 0; j < qualities.size(); j++)
      log_probs[j] = log_prob_correct(qualities[j]->at(i));
    double log_mean_prob = log_sum_exp(log_probs) - log(qualities.size());
    log_pcorrect[i] = log_mean_prob;
    log_perror[i]   = log(1.0 - exp(log_mean_prob));
  }
}
