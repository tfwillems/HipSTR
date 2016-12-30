#ifndef VCF_INPUT_H_
#define VCF_INPUT_H_

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "error.h"
#include "region.h"
#include "vcf_reader.h"

extern const std::string GENOTYPE_KEY;
extern const std::string PHASED_GL_KEY;
extern std::string START_INFO_TAG;
extern std::string STOP_INFO_TAG;
extern const int32_t pad;

bool read_vcf_alleles(VCF::VCFReader* ref_vcf, const Region& region, std::vector<std::string>& alleles, int32_t& pos);

class PhasedGL{
 private:
  int num_alleles_;
  int num_samples_;
  std::map<std::string, int> sample_indices_;
  std::vector< std::vector<float> > phased_gls_;
  std::vector< std::vector<float> > max_gls_one_;
  std::vector< std::vector<float> > max_gls_two_;

  bool build(VCF::Variant& variant);

 public:
  PhasedGL(){
    num_alleles_ = 0;
    num_samples_ = 0;
  }

  PhasedGL(VCF::Variant& variant){
    if (!build(variant))
      printErrorAndDie("Failed to construct PhasedGL instance from VCF record");
  }

  bool has_sample(const std::string& sample){
    return sample_indices_.find(sample) != sample_indices_.end();
  }

  int get_sample_index(const std::string& sample){
    auto sample_iter = sample_indices_.find(sample);
    return (sample_iter == sample_indices_.end() ? -1 : sample_iter->second);
  }

  float get_gl(const std::string& sample, int gt_a, int gt_b){
    auto sample_iter = sample_indices_.find(sample);
    if (sample_iter == sample_indices_.end())
      printErrorAndDie("No data available for sample " + sample + " in PhasedGL instance");
    if (gt_a >= num_alleles_)
      printErrorAndDie("Genotype index exceeds the number of alleles present in PhasedGL instance");
    if (gt_b >= num_alleles_)
      printErrorAndDie("Genotype index exceeds the number of alleles present in PhasedGL instance");

    int gt_index = gt_a*num_alleles_ + gt_b;
    return phased_gls_[sample_iter->second][gt_index];
  }

  float get_gl(int sample_index, int gt_a, int gt_b){
    return phased_gls_[sample_index][gt_a*num_alleles_ + gt_b];
  }

  /*
   * For the relevant sample, returns the maximum phased GL of all genotypes
   * of the form GT_A | X, where X is any valid allele
   */
  float get_max_gl_allele_one_fixed(int sample_index, int gt_a){
    return max_gls_one_[sample_index][gt_a];
  }

  /*
   * For the relevant sample, returns the maximum phased GL of all genotypes
   * of the form X | GT_A, where X is any valid allele
   */
  float get_max_gl_allele_two_fixed(int sample_index, int gt_b){
    return max_gls_two_[sample_index][gt_b];
  }

};

#endif
