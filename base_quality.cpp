#include <assert.h>

#include "base_quality.h"
#include "error.h"
#include "mathops.h"

std::string BaseQuality::average_base_qualities(std::vector<const std::string*> qualities){
  assert(qualities.size() > 0);

  // Check that all base quality strings are of the same length
  for (unsigned int i = 0; i < qualities.size(); i++){
    if (qualities[i]->size() != qualities[0]->size())
      printErrorAndDie("All base quality strings must be of the same length when averaging probabilities");
  }

  // Average raw error probabilities for each base and convert
  // to the closest quality score
  std::string avg_qualities('N', qualities[0]->size());
  std::vector<double> log_probs(qualities.size());
  for (unsigned int i = 0; i < qualities[0]->size(); i++){
    for (unsigned int j = 0; j < qualities.size(); j++)
      log_probs[j] = log_prob_error(qualities[j]->at(i));
    double log_mean_prob = log_sum_exp(log_probs) - log(qualities.size());
    avg_qualities[i]     = closest_char(log_mean_prob);
  }
  return avg_qualities;
}
