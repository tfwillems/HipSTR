#ifndef TRIO_DENOVO_SCANNER_H_
#define TRIO_DENOVO_SCANNER_H_

#include <assert.h>

#include <iostream>
#include <vector>
#include <string>

#include "bgzf_streams.h"
#include "pedigree.h"
#include "vcf_reader.h"

class TrioDenovoScanner {
 private:
  static std::string BPDIFFS_KEY, START_KEY, END_KEY, PERIOD_KEY;

  std::vector<NuclearFamily> families_;
  bgzfostream denovo_vcf_;
  bool use_pop_priors_;

  void write_vcf_header(const std::string& full_command);
  void initialize_vcf_record(const VCF::Variant& str_variant);
  void add_child_to_record(double total_ll_no_denovo, double total_ll_one_denovo, double total_ll_one_other);

 public:
  TrioDenovoScanner(const std::vector<NuclearFamily>& families, const std::string& output_file, const std::string& full_command, bool use_pop_priors)
    : families_(families){
    use_pop_priors_ = use_pop_priors;
    denovo_vcf_.open(output_file.c_str());
    denovo_vcf_.precision(3);
    denovo_vcf_.setf(std::ios::fixed, std::ios::floatfield);
    write_vcf_header(full_command);
  }

  void scan(VCF::VCFReader& str_vcf, std::ostream& logger);

  void finish(){ denovo_vcf_.close(); }
};

#endif
