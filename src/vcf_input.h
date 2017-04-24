#ifndef VCF_INPUT_H_
#define VCF_INPUT_H_

#include <assert.h>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "error.h"
#include "region.h"
#include "vcf_reader.h"

extern const std::string GENOTYPE_KEY;
extern const std::string UNPHASED_GL_KEY;
extern const std::string PHASED_GL_KEY;
extern std::string START_INFO_TAG;
extern std::string STOP_INFO_TAG;
extern const int32_t pad;

bool read_vcf_alleles(VCF::VCFReader* ref_vcf, const Region& region, std::vector<std::string>& alleles, int32_t& pos);

class GL {
 protected:
  int num_alleles_;
  int num_samples_;
  std::map<std::string, int> sample_indices_;

 public:
  GL(){
    num_alleles_ = 0;
    num_samples_ = 0;
  }

  bool has_sample(const std::string& sample) const {
    return sample_indices_.find(sample) != sample_indices_.end();
  }

  int get_sample_index(const std::string& sample) const {
    auto sample_iter = sample_indices_.find(sample);
    return (sample_iter == sample_indices_.end() ? -1 : sample_iter->second);
  }

  virtual float get_gl(int sample_index, int gt_a, int gt_b) const = 0;
};

class UnphasedGL : public GL {
 private:
  std::vector< std::vector<float> > unphased_gls_;
  std::vector< std::vector<float> > max_gls_;

  bool build(const VCF::Variant& variant);

 public:
  explicit UnphasedGL(const VCF::Variant& variant){
    if (!variant.has_format_field(UNPHASED_GL_KEY))
      printErrorAndDie("Required FORMAT field " + UNPHASED_GL_KEY + " not present in VCF");
    if (!build(variant))
      printErrorAndDie("Failed to construct UnphasedGL instance from VCF record");
  }

  float get_gl(int sample_index, int min_gt, int max_gt) const {
    assert(min_gt <= max_gt);
    return unphased_gls_[sample_index][max_gt*(max_gt+1)/2 + min_gt];
  }

  /*
   * For the relevant sample, returns the maximum unphased GL of all genotypes
   * that contain GT_A as an allele
   */
  float get_max_gl_allele_fixed(int sample_index, int gt_a) const {
    return max_gls_[sample_index][gt_a];
  }
};


class PhasedGL : public GL {
 private:
  std::vector< std::vector<float> > phased_gls_;
  std::vector< std::vector<float> > max_gls_one_;
  std::vector< std::vector<float> > max_gls_two_;

  bool build(const VCF::Variant& variant);

 public:
  explicit PhasedGL(const VCF::Variant& variant){
   if (!variant.has_format_field(PHASED_GL_KEY))
     printErrorAndDie("Required FORMAT field " + PHASED_GL_KEY + " not present in VCF");
   if (!build(variant))
     printErrorAndDie("Failed to construct PhasedGL instance from VCF record");
  }

  float get_gl(int sample_index, int gt_a, int gt_b) const {
    return phased_gls_[sample_index][gt_a*num_alleles_ + gt_b];
  }

  /*
   * For the relevant sample, returns the maximum phased GL of all genotypes
   * of the form GT_A | X, where X is any valid allele
   */
  float get_max_gl_allele_one_fixed(int sample_index, int gt_a) const {
    return max_gls_one_[sample_index][gt_a];
  }

  /*
   * For the relevant sample, returns the maximum phased GL of all genotypes
   * of the form X | GT_A, where X is any valid allele
   */
  float get_max_gl_allele_two_fixed(int sample_index, int gt_b) const {
    return max_gls_two_[sample_index][gt_b];
  }
};

#endif
