#ifndef DENOVO_ALLELE_PRIORS_H_
#define DENOVO_ALLELE_PRIORS_H_

#include <assert.h>
#include <math.h>

#include <vector>
#include <string>

#include "error.h"
#include "pedigree.h"
#include "vcf_reader.h"

class DiploidGenotypePrior {
 protected:
  int num_alleles_;
  std::vector<double> allele_freqs_, log_allele_freqs_;
  double LOG_2;

  DiploidGenotypePrior(VCF::Variant& str_variant, std::vector<NuclearFamily>& families){
    num_alleles_ = str_variant.num_alleles();
    assert(num_alleles_ > 0);
    LOG_2 = log10(2);
  }

 public:
  /* Returns the log10 prior for the given unphased genotype, assuming Hardy-Weinberg equilibrium */
  double log_unphased_genotype_prior(int gt_a, int gt_b, const std::string& sample){
    if (gt_a < 0 || gt_a >= num_alleles_)
      printErrorAndDie("Invalid genotype index for log genotype prior");
    if (gt_b < 0 || gt_b >= num_alleles_)
      printErrorAndDie("Invalid genotype index for log genotype prior");
    if (gt_a == gt_b)
      return log_allele_freqs_[gt_a] + log_allele_freqs_[gt_b];
    else
      return log_allele_freqs_[gt_a] + log_allele_freqs_[gt_b] + LOG_2;
  }

  /* Returns the log10 prior for the given phased genotype, assuming Hardy-Weinberg equilibrium */
  double log_phased_genotype_prior(int gt_a, int gt_b, const std::string& sample){
    if (gt_a < 0 || gt_a >= num_alleles_)
      printErrorAndDie("Invalid genotype index for log genotype prior");
    if (gt_b < 0 || gt_b >= num_alleles_)
      printErrorAndDie("Invalid genotype index for log genotype prior");
    return log_allele_freqs_[gt_a] + log_allele_freqs_[gt_b];
  }
};


/*
 * Simple class the computes an STR's allele frequencies from a VCF and uses the
 * resulting values to produce genotype priors, assuming Hardy-Weinberg equilibrium
 */
class PopulationGenotypePrior : public DiploidGenotypePrior {
 protected:
  void compute_allele_freqs(VCF::Variant& variant, std::vector<NuclearFamily>& families);

 public:
 PopulationGenotypePrior(VCF::Variant& str_variant, std::vector<NuclearFamily>& families)
   : DiploidGenotypePrior(str_variant, families){
    compute_allele_freqs(str_variant, families);
  }
};

/*
 * Simple class that assigns all alleles a uniform frequency and uses the resulting
 * values to produce genotype priors, assuming Hardy-Weinberg equilibrium.
 * Results in unphased homozygous genotypes having half the prior likelihood of unphased heterozygous genotypes
 */
class UniformGenotypePrior : public DiploidGenotypePrior {
 protected:
  void compute_allele_freqs(VCF::Variant& variant, std::vector<NuclearFamily>& families);

 public:
  UniformGenotypePrior(VCF::Variant& str_variant, std::vector<NuclearFamily>& families)
    : DiploidGenotypePrior(str_variant, families){
    compute_allele_freqs(str_variant, families);
  }
};

#endif
