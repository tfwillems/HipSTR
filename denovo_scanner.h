#ifndef DENOVO_SCANNER_H_
#define DENOVO_SCANNER_H_

#include <assert.h>

#include <vector>
#include <set>
#include <string>

#include "pedigree.h"
#include "vcflib/src/Variant.h"

class DiploidGenotypePrior {
 private:
  int num_alleles_;
  std::vector<double> allele_freqs_, log_allele_freqs_;

  void compute_allele_freqs(vcflib::Variant& variant, std::vector<NuclearFamily>& families);
 public:
  DiploidGenotypePrior(vcflib::Variant& str_variant, std::vector<NuclearFamily>& families){
    num_alleles_ = str_variant.alleles.size();
    if (str_variant.alleles.back().compare(".") == 0)
      num_alleles_--;
    assert(num_alleles_ > 0);
    compute_allele_freqs(str_variant, families);
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

class DenovoScanner {
  const static int MIN_SECOND_BEST_SCORE = 50;
  const static int MAX_BEST_SCORE        = 10;

 private:
  int32_t window_size_;
  std::vector<NuclearFamily> families_;

 public:
  DenovoScanner(std::vector<NuclearFamily>& families){
    families_    = families;
    window_size_ = 500000;
  }
  
  void scan(vcflib::VariantCallFile& snp_vcf, vcflib::VariantCallFile& str_vcf, std::set<std::string>& sites_to_skip,
	    std::ostream& logger);
};

#endif
