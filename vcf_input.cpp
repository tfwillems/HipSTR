#include <algorithm>
#include <assert.h>
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

void read_vcf_alleles(VCF::VCFReader* ref_vcf, Region* region, std::vector<std::string>& alleles, int32_t& pos, bool& success){
    assert(alleles.size() == 0 && ref_vcf != NULL);
    int32_t pad_start = (region->start() < pad ? 0 : region->start()-pad);
    if (!ref_vcf->set_region(region->chrom(), pad_start, region->stop()+pad)){
      // Retry setting region if chr is in chromosome name
      if (region->chrom().size() <= 3 || region->chrom().substr(0, 3).compare("chr") != 0 
	  || !ref_vcf->set_region(region->chrom().substr(3), pad_start, region->stop()+pad)){
	success = false;
	pos     = -1;
	return;
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
      if (str_start == region->start()+1 && str_stop == region->stop()){
	success = true;
	pos     = variant.get_position()-1;
	alleles.insert(alleles.end(), variant.get_alleles().begin(), variant.get_alleles().end());
	return;
      }
      if (variant.get_position() > region->start()+pad)
	break;
    }

    success = false;
    pos     = -1;
}


/*
 * Searchs for an entry in the provided VCF that matches the region. If found, stores the alleles in the provided vector
 * and returns an array of size NUM_ALLELES*NUM_ALLELES*NUM_SAMPLES containing the log prior for each sample's diploid genotype.
 * The array iterates over allele_1, allele_2 and then each sample.
 * Method exits with an error if no VCF entry is found, if the VCF doesn't containg PGP allele priors in the FORMAT field or if it only contains a subset of the samples.
 * The user is responsible for freeing the returned array when it is no longer needed.
 */

double* extract_vcf_alleles_and_log_priors(VCF::VCFReader* ref_vcf, Region* region, std::map<std::string, int>& sample_indices,
					   std::vector<std::string>& alleles, std::vector<bool>& got_priors, int32_t& pos, bool& success, std::ostream& logger){
  assert(alleles.size() == 0 && got_priors.size() == 0);
  got_priors.resize(sample_indices.size(), false);
  int32_t pad_start = (region->start() < pad ? 0 : region->start()-pad);

  if (!ref_vcf->set_region(region->chrom(), pad_start, region->stop()+pad)){
    // Retry setting region if chr is in chromosome name
    if (region->chrom().size() <= 3 || region->chrom().substr(0, 3).compare("chr") != 0 
	|| !ref_vcf->set_region(region->chrom().substr(3), pad_start, region->stop()+pad)){
      success = false;
      pos     = -1;
      return NULL;
    }
  }
  VCF::Variant variant;
  bool matches_region = false;
  while(ref_vcf->get_next_variant(variant)){
    // Skip variants without the appropriate INFO fields (as they're not STRs)
    if (!variant.has_info_field(START_INFO_TAG) || !variant.has_info_field(STOP_INFO_TAG) || !variant.has_format_field(PGP_KEY))
      continue;

    int32_t str_start, str_stop;
    variant.get_INFO_value_single_int(START_INFO_TAG, str_start);
    variant.get_INFO_value_single_int(STOP_INFO_TAG, str_stop);
    if (str_start == region->start()+1 && str_stop == region->stop()){
      matches_region = true;
      break;
    }
  }
  if (!matches_region){
    success = false;
    pos     = -1;
    return NULL;
  }

  // Extract and store the number of alleles and each of their sequences
  success = true;
  pos     = variant.get_position()-1;
  alleles.insert(alleles.end(), variant.get_alleles().begin(), variant.get_alleles().end());
  
  // Allocate allele prior storage
  size_t num_samples  = sample_indices.size();
  size_t num_alleles  = variant.get_alleles().size();
  double* log_allele_priors = new double[num_alleles*num_alleles*num_samples];

  // Initialize array with what is equivalent to log of uniform prior
  // Results in valid prior for samples without VCF priors
  std::fill(log_allele_priors, log_allele_priors+(num_alleles*num_alleles*num_samples), -2*log(num_alleles));

  // Extract priors for each sample
  size_t sample_count = 0;
  std::vector< std::vector<float> > gp_probs;
  variant.get_FORMAT_value_multiple_floats(PGP_KEY, gp_probs);
  int vcf_sample_index = 0;
  int exp_num_vals     = num_alleles*num_alleles;
  for (auto sample_iter = variant.get_samples().begin(); sample_iter != variant.get_samples().end(); ++sample_iter, ++vcf_sample_index){
    if (sample_indices.find(*sample_iter) == sample_indices.end())
      continue;
    int sample_index = sample_indices[*sample_iter];
    float total      = 0.0;
    got_priors[sample_index] = true;

    if (gp_probs[vcf_sample_index].size() != exp_num_vals)
      printErrorAndDie("Number of items in PGP FORMAT field does not match the expected value");

    std::vector<float>& probs = gp_probs[vcf_sample_index];
    for (int i = 0; i < probs.size(); i++){
      probs[i] = std::max(probs[i], MIN_ALLELE_PRIOR);
      total += probs[i];
    }

    // Normalize and log-transform priors and store at appropriate index
    size_t gp_index = 0;
    double* log_prior_ptr = log_allele_priors + sample_index;
    for (size_t i = 0; i < num_alleles; ++i){
      for (size_t j = 0; j < num_alleles; ++j, ++gp_index){
	*log_prior_ptr = log(probs[gp_index]/total);
	log_prior_ptr += num_samples;
      }
    }

    sample_count++;
  }
    
  // Warn if the VCF did not contain priors for all samples
  if (sample_count != num_samples)
    logger << "WARNING: VCF only contained allele priors for " << sample_count << " out of " << num_samples << " samples";
  return log_allele_priors;
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
  }

  return true;
}
