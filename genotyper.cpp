#include <assert.h>
#include <cfloat>
#include <cstring>
#include <time.h>

#include <algorithm>

#include "genotyper.h"
#include "mathops.h"

// Each genotype has an equal total prior, but heterozygotes have two possible phasings. Therefore,
// i)   Phased heterozygotes have a prior of 1/(n(n+1))
// ii)  Homozygotes have a prior of 2/(n(n+1))
// iii) Total prior is n*2/(n(n+1)) + n(n-1)*1/(n(n+1)) = 2/(n+1) + (n-1)/(n+1) = 1

double Genotyper::log_homozygous_prior(){
  if (haploid_)
    return -int_log(num_alleles_);
  else
    return int_log(2) - int_log(num_alleles_) - int_log(num_alleles_+1);
}

double Genotyper::log_heterozygous_prior(){
  if (haploid_)
    return -DBL_MAX/2;
  else
    return -int_log(num_alleles_) - int_log(num_alleles_+1);
}

void Genotyper::init_log_sample_priors(double* log_sample_ptr){
  if (log_allele_priors_ != NULL)
    std::memcpy(log_sample_ptr, log_allele_priors_, num_alleles_*num_alleles_*num_samples_*sizeof(double));
  else {
    // Set all elements to the het prior
    std::fill(log_sample_ptr, log_sample_ptr+(num_alleles_*num_alleles_*num_samples_), log_heterozygous_prior());

    // Fix homozygotes
    double log_homoz_prior = log_homozygous_prior();
    for (unsigned int i = 0; i < num_alleles_; i++){
      double* LL_ptr = log_sample_ptr + i*num_alleles_*num_samples_ + i*num_samples_;
      std::fill(LL_ptr, LL_ptr+num_samples_, log_homoz_prior);
    }
  }
}

double Genotyper::calc_log_sample_posteriors(std::vector<int>& read_weights){
  double posterior_time = clock();
  assert(read_weights.size() == num_reads_);
  std::vector<double> sample_max_LLs(num_samples_, -DBL_MAX);
  double* sample_LL_ptr = log_sample_posteriors_;
  init_log_sample_priors(log_sample_posteriors_);

  for (int index_1 = 0; index_1 < num_alleles_; ++index_1){
    for (int index_2 = 0; index_2 < num_alleles_; ++index_2){
      double* read_LL_ptr = log_aln_probs_;
      for (int read_index = 0; read_index < num_reads_; ++read_index){
        sample_LL_ptr[sample_label_[read_index]] += read_weights[read_index]*logsumexp_agg(LOG_ONE_HALF + log_p1_[read_index] + read_LL_ptr[index_1], 
											   LOG_ONE_HALF + log_p2_[read_index] + read_LL_ptr[index_2]);
        assert(sample_LL_ptr[sample_label_[read_index]] <= TOLERANCE);
        read_LL_ptr += num_alleles_;
      }
      // Update the per-sample maximum LLs
      for (int sample_index = 0; sample_index < num_samples_; ++sample_index)
        sample_max_LLs[sample_index] = std::max(sample_max_LLs[sample_index], sample_LL_ptr[sample_index]);

      sample_LL_ptr += num_samples_;
    }
  }

  // Compute the normalizing factor for each sample using logsumexp trick
  std::fill(sample_total_LLs_, sample_total_LLs_ + num_samples_, 0.0);
  sample_LL_ptr = log_sample_posteriors_;
  for (int index_1 = 0; index_1 < num_alleles_; ++index_1)
    for (int index_2 = 0; index_2 < num_alleles_; ++index_2)
      for (int sample_index = 0; sample_index < num_samples_; ++sample_index, ++sample_LL_ptr)
        sample_total_LLs_[sample_index] += exp(*sample_LL_ptr - sample_max_LLs[sample_index]);
  for (int sample_index = 0; sample_index < num_samples_; ++sample_index){
    sample_total_LLs_[sample_index] = sample_max_LLs[sample_index] + log(sample_total_LLs_[sample_index]);
    assert(sample_total_LLs_[sample_index] <= TOLERANCE);
  }
  // Compute the total log-likelihood given the current parameters
  double total_LL = sum(sample_total_LLs_, sample_total_LLs_ + num_samples_);

  // Normalize each genotype LL to generate valid log posteriors
  sample_LL_ptr = log_sample_posteriors_;
  for (int index_1 = 0; index_1 < num_alleles_; ++index_1)
    for(int index_2 = 0; index_2 < num_alleles_; ++index_2)
      for (int sample_index = 0; sample_index < num_samples_; ++sample_index, ++sample_LL_ptr)
        *sample_LL_ptr -= sample_total_LLs_[sample_index];

  posterior_time         = (clock() - posterior_time)/CLOCKS_PER_SEC;
  total_posterior_time_ += posterior_time;
  return total_LL;
}

void Genotyper::get_optimal_genotypes(double* log_posterior_ptr, std::vector< std::pair<int, int> >& gts){
  assert(gts.size() == 0);
  gts = std::vector< std::pair<int,int> > (num_samples_, std::pair<int,int>(-1,-1));
  std::vector<double> log_phased_posteriors(num_samples_, -DBL_MAX);
  for (int index_1 = 0; index_1 < num_alleles_; ++index_1){
    for (int index_2 = 0; index_2 < num_alleles_; ++index_2){
      for (unsigned int sample_index = 0; sample_index < num_samples_; ++sample_index, ++log_posterior_ptr){
        if (*log_posterior_ptr > log_phased_posteriors[sample_index]){
          log_phased_posteriors[sample_index] = *log_posterior_ptr;
          gts[sample_index] = std::pair<int,int>(index_1, index_2);
        }
      }
    }
  }
}
