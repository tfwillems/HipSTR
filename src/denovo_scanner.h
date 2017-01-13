#ifndef DENOVO_SCANNER_H_
#define DENOVO_SCANNER_H_

#include <assert.h>
#include <math.h>

#include <iostream>
#include <vector>
#include <set>
#include <string>

#include "bgzf_streams.h"
#include "pedigree.h"
#include "vcf_reader.h"

class DenovoScanner {
 private:
  const static int MIN_SECOND_BEST_SCORE = 100;
  const static int MAX_BEST_SCORE        = 10;
  static std::string BPDIFFS_KEY, START_KEY, END_KEY, PERIOD_KEY;
  bool use_pop_priors_;

  int32_t window_size_;
  std::vector<NuclearFamily> families_;
  bgzfostream denovo_vcf_;

  void write_vcf_header(std::string& full_command);
  void initialize_vcf_record(VCF::Variant& str_variant);
  void add_family_to_record(NuclearFamily& family, double total_ll_no_denovo, std::vector<double>& total_lls_one_denovo, std::vector<double>& total_lls_one_other);

 public:
  DenovoScanner(std::vector<NuclearFamily>& families, std::string& output_file, std::string& full_command, bool use_pop_priors){
    families_       = families;
    use_pop_priors_ = use_pop_priors;
    window_size_    = 500000;
    denovo_vcf_.open(output_file.c_str());
    denovo_vcf_.precision(3);
    denovo_vcf_.setf(std::ios::fixed, std::ios::floatfield);
    write_vcf_header(full_command);
  }

  void scan(std::string& snp_vcf_file, VCF::VCFReader& str_vcf, std::set<std::string>& sites_to_skip,
	    std::ostream& logger);

  void finish(){ denovo_vcf_.close(); }
};

#endif
