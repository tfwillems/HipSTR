#include <math.h>
#include <string>
#include <vector>

#include "denovo_allele_priors.h"

void PopulationGenotypePrior::compute_allele_freqs(VCF::Variant& variant, std::vector<NuclearFamily>& families){
  allele_freqs_ = std::vector<double>(num_alleles_, 1.0); // Use a one sample pseudocount

  // Iterate over all founders in the families to compute allele counts
  double total_count = num_alleles_;
  int gt_a, gt_b;
  for (auto family_iter = families.begin(); family_iter != families.end(); family_iter++){
    for (int i = 0; i < 2; i++){
      std::string sample = (i == 0 ? family_iter->get_mother() : family_iter->get_father());
      if (variant.sample_call_missing(sample))
	continue;
      variant.get_genotype(sample, gt_a, gt_b);
      allele_freqs_[gt_a]++;
      allele_freqs_[gt_b]++;
      total_count += 2;
    }
  }

  // Normalize the allele counts to obtain frequencies
  for (int i = 0; i < allele_freqs_.size(); i++)
    allele_freqs_[i] /= total_count;

  // Precompute the logs of the allele frequencies
  log_allele_freqs_.clear();
  for (int i = 0; i < allele_freqs_.size(); i++)
    log_allele_freqs_.push_back(log10(allele_freqs_[i]));
}

void UniformGenotypePrior::compute_allele_freqs(VCF::Variant& variant, std::vector<NuclearFamily>& families){
  allele_freqs_     = std::vector<double>(num_alleles_, 1.0/num_alleles_);
  log_allele_freqs_ = std::vector<double>(num_alleles_, -log10(num_alleles_));
}

