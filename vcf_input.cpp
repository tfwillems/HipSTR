#include <algorithm>
#include <assert.h>
#include <cfloat>
#include <math.h>

#include "error.h"
#include "region.h"
#include "vcf_input.h"

const std::string GT_KEY        = "GT";
const std::string PHASED_GL_KEY = "PHASEDGL";
std::string PGP_KEY             = "PGP";
std::string START_INFO_TAG      = "START";
std::string STOP_INFO_TAG       = "END";
const float MIN_ALLELE_PRIOR    = 0.0001;

// Because HipSTR extends putative STR regions if there are nearby indels, the STR coordinates in the VCF may
// not exactly match the original reference region coordinates. As a result, when looking for a particular STR region,
// we look for entries a window around the locus. The size of this window is controlled by this parameter
const int32_t pad = 50;

bool read_vcf_alleles(VCF::VCFReader* ref_vcf, const Region& region, std::vector<std::string>& alleles, int32_t& pos){
    assert(alleles.size() == 0 && ref_vcf != NULL);
    int32_t pad_start = (region.start() < pad ? 0 : region.start()-pad);
    if (!ref_vcf->set_region(region.chrom(), pad_start, region.stop()+pad)){
      // Retry setting region if chr is in chromosome name
      if (region.chrom().size() <= 3 || region.chrom().substr(0, 3).compare("chr") != 0
	  || !ref_vcf->set_region(region.chrom().substr(3), pad_start, region.stop()+pad)){
	pos     = -1;
	return false;
      }
    }
   
    // Extract STR and ensure the coordinates match
    VCF::Variant variant;
    while (ref_vcf->get_next_variant(variant)){
      // Skip variants without the appropriate INFO fields (as they're not STRs)
      if (!variant.has_info_field(START_INFO_TAG) || !variant.has_info_field(STOP_INFO_TAG))
	continue;

      int32_t str_start, str_stop;
      variant.get_INFO_value_single_int(START_INFO_TAG, str_start);
      variant.get_INFO_value_single_int(STOP_INFO_TAG, str_stop);
      if (str_start == region.start()+1 && str_stop == region.stop()){
	pos     = variant.get_position()-1;
	alleles.insert(alleles.end(), variant.get_alleles().begin(), variant.get_alleles().end());
	return true;
      }
      if (variant.get_position() > region.start()+pad)
	break;
    }

    pos = -1;
    return false;
}

bool PhasedGL::build(VCF::Variant& variant){
  if (!variant.has_format_field(PHASED_GL_KEY))
    return false;

  std::vector< std::vector<float> > values;
  variant.get_FORMAT_value_multiple_floats(PHASED_GL_KEY, values);
  num_samples_         = 0;
  num_alleles_         = variant.num_alleles();
  int vcf_sample_index = 0;

  const std::vector<std::string>& samples = variant.get_samples();
  for (auto sample_iter = samples.begin(); sample_iter != samples.end(); ++sample_iter, ++vcf_sample_index){
    if (variant.sample_call_missing(vcf_sample_index))
      continue;
    phased_gls_.push_back(values[vcf_sample_index]);
    sample_indices_[*sample_iter] = num_samples_++;

    std::vector<float> max_gl_one(num_alleles_, -DBL_MAX/2), max_gl_two(num_alleles_, -DBL_MAX/2);
    int gl_index = 0;
    for (int i = 0; i < num_alleles_; ++i){
      for (int j = 0; j < num_alleles_; ++j, ++gl_index){
	max_gl_one[i] = std::max(max_gl_one[i], phased_gls_.back()[gl_index]);
	max_gl_two[j] = std::max(max_gl_two[j], phased_gls_.back()[gl_index]);
      }
    }

    max_gls_one_.push_back(max_gl_one);
    max_gls_two_.push_back(max_gl_two);
  }

  return true;
}
