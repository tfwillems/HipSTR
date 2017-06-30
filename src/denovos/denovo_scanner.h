#ifndef DENOVO_SCANNER_H_
#define DENOVO_SCANNER_H_

#include <assert.h>
#include <math.h>

#include <iostream>
#include <vector>
#include <set>
#include <string>

#include "../bgzf_streams.h"
#include "../pedigree.h"
#include "../vcf_reader.h"

class DenovoScanner {
 public:
  const static int MIN_SECOND_BEST_SCORE = 100;
  const static int MAX_BEST_SCORE        = 10;

 private:
  static std::string BPDIFFS_KEY, START_KEY, END_KEY, PERIOD_KEY;
  bool use_pop_priors_;

  int32_t window_size_;
  std::vector<NuclearFamily> families_;
  bgzfostream denovo_vcf_;

  void write_vcf_header(const std::string& full_command);
  void initialize_vcf_record(const VCF::Variant& str_variant);
  void add_family_to_record(const NuclearFamily& family, double total_ll_no_denovo,
			    const std::vector<double>& total_lls_one_denovo, const std::vector<double>& total_lls_one_other);

  // Private unimplemented copy constructor and assignment operator to prevent operations
  DenovoScanner(const DenovoScanner& other);
  DenovoScanner& operator=(const DenovoScanner& other);

 public:
 DenovoScanner(const std::vector<NuclearFamily>& families, const std::string& output_file, const std::string& full_command, bool use_pop_priors)
   : families_(families){
    use_pop_priors_ = use_pop_priors;
    window_size_    = 500000;
    denovo_vcf_.open(output_file.c_str());
    denovo_vcf_.precision(3);
    denovo_vcf_.setf(std::ios::fixed, std::ios::floatfield);
    write_vcf_header(full_command);
  }

  void scan(const std::string& snp_vcf_file, VCF::VCFReader& str_vcf, const std::set<std::string>& sites_to_skip,
	    std::ostream& logger);

  void finish(){ denovo_vcf_.close(); }
};

#endif
