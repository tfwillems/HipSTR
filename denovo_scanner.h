#ifndef DENOVO_SCANNER_H_
#define DENOVO_SCANNER_H_

#include <assert.h>

#include <iostream>
#include <vector>
#include <set>
#include <sstream>
#include <string>

#include "bgzf_streams.h"
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
  const static int MIN_SECOND_BEST_SCORE = 100;
  const static int MAX_BEST_SCORE        = 10;

 private:
  int32_t window_size_;
  std::vector<NuclearFamily> families_;
  bgzfostream denovo_vcf_;

  void write_vcf_header(std::string& full_command);
  void initialize_vcf_record(vcflib::Variant& str_variant);
  void add_family_to_record(NuclearFamily& family, double total_ll_no_denovo, std::vector<double>& total_lls_one_denovo, std::vector<double>& total_lls_one_other);


 public:
  DenovoScanner(std::vector<NuclearFamily>& families, std::string& output_file, std::string& full_command){
    families_    = families;
    window_size_ = 500000;
    denovo_vcf_.open(output_file.c_str());
    denovo_vcf_.precision(3);
    denovo_vcf_.setf(std::ios::fixed, std::ios::floatfield);
    write_vcf_header(full_command);
  }

  void scan(std::string& snp_vcf_file, vcflib::VariantCallFile& str_vcf, std::set<std::string>& sites_to_skip,
	    std::ostream& logger);

  void finish(){ denovo_vcf_.close(); }
};

#endif
