#include <assert.h>

#include "error.h"
#include "region.h"
#include "vcf_input.h"

std::string PGP_KEY           = "PGP";
std::string START_INFO_TAG    = "START";
std::string STOP_INFO_TAG     = "END";
const double MIN_ALLELE_PRIOR = 0.001;

// Because HipSTR extends putative STR regions if there are nearby indels, the STR coordinates in the VCF may
// not exactly match the original reference region coordinates. As a result, when looking for a particular STR region,
// we look for entries a window around the locus. The size of this window is controlled by this parameter
const int32_t pad = 50;

void read_vcf_alleles(vcf::VariantCallFile* ref_vcf, Region* region, std::vector<std::string>& alleles, int32_t& pos){
    assert(alleles.size() == 0 && ref_vcf != NULL);
    if (!ref_vcf->setRegion(region->chrom(), region->start()-pad, region->stop()+pad)){
      // Retry setting region if chr is in chromosome name
      if (region->chrom().size() <= 3 || region->chrom().substr(0, 3).compare("chr") != 0 
	  || !ref_vcf->setRegion(region->chrom().substr(3), region->start()-pad, region->stop()+pad))
	printErrorAndDie("Failed to set VCF region when reading candidate STR alleles for region " + region->str());
    }
   
    // Extract STR and ensure the coordinates match
    vcf::Variant variant(*ref_vcf);
    while (ref_vcf->getNextVariant(variant)){
      // Skip variants without the appropriate INFO fields (as they're not STRs)
      if (ref_vcf->infoTypes.find(START_INFO_TAG) == ref_vcf->infoTypes.end())
	continue;
      if (ref_vcf->infoTypes.find(STOP_INFO_TAG) == ref_vcf->infoTypes.end())
	continue;

      int32_t str_start = (int32_t)variant.getInfoValueFloat(START_INFO_TAG);
      int32_t str_stop  = (int32_t)variant.getInfoValueFloat(STOP_INFO_TAG);
      if (str_start == region->start() && str_stop == region->stop()){
	pos = variant.position;
	alleles.insert(alleles.end(), variant.alleles.begin(), variant.alleles.end());
	return;
      }
      if (variant.position > region->start()+pad)
	break;
    }
    printErrorAndDie("Failed to extract matching VCF entry when reading candidate STR alleles for region " + region->str());
}


/*
 * Searchs for an entry in the provided VCF that matches the region. If found, stores the alleles in the provided vector
 * and returns an array of size NUM_ALLELES*NUM_ALLELES*NUM_SAMPLES containing the log prior for each sample's diploid genotype.
 * The array iterates over allele_1, allele_2 and then each sample.
 * Method exits with an error if no VCF entry is found, if the VCF doesn't containg PGP allele priors in the FORMAT field or if it only contains a subset of the samples.
 * The user is responsible for freeing the returned array when it is no longer needed.
 */
double* extract_vcf_alleles_and_log_priors(vcf::VariantCallFile* ref_vcf, Region* region, std::map<std::string, int>& sample_indices,
					   std::vector<std::string>& alleles, int32_t& pos){
  assert(alleles.size() == 0);
  if (ref_vcf->formatTypes.find(PGP_KEY) == ref_vcf->formatTypes.end())
    printErrorAndDie("VCF doesn't contain the PGP format field required for setting allele priors");
  if (!ref_vcf->setRegion(region->chrom(), region->start()-pad, region->stop()+pad)){
    // Retry setting region if chr is in chromosome name
    if (region->chrom().size() <= 3 || region->chrom().substr(0, 3).compare("chr") != 0 
	|| !ref_vcf->setRegion(region->chrom().substr(3), region->start()-pad, region->stop()+pad))
      printErrorAndDie("Failed to set VCF region when obtaining allele priors for region " + region->str());
  }
  vcf::Variant variant(*ref_vcf);
  bool matches_region = false;
  while(ref_vcf->getNextVariant(variant)){
    // Skip variants without the appropriate INFO fields (as they're not STRs)
    if (ref_vcf->infoTypes.find(START_INFO_TAG) == ref_vcf->infoTypes.end())
      continue;
    if (ref_vcf->infoTypes.find(STOP_INFO_TAG) == ref_vcf->infoTypes.end())
      continue;
    
    int32_t str_start = (int32_t)variant.getInfoValueFloat(START_INFO_TAG);
    int32_t str_stop  = (int32_t)variant.getInfoValueFloat(STOP_INFO_TAG);
    if (str_start == region->start() && str_stop == region->stop()){
      matches_region = true;
      break;
    }
  }
  if (!matches_region)
    printErrorAndDie("Failed to extract VCF entry when obtaining allele priors for region " + region->str());

  // Extract and store the number of alleles and each of their sequences
  pos = variant.position;
  alleles.insert(alleles.end(), variant.alleles.begin(), variant.alleles.end());
  
  // Allocate allele prior storage
  int num_samples = sample_indices.size();
  int num_alleles = variant.alleles.size();
  double* log_allele_priors = new double[num_alleles*num_alleles*num_samples];

  // Extract priors for each sample
  int sample_count = 0;
  std::vector<double> gp_probs; gp_probs.reserve(num_alleles*num_alleles);
  for (auto sample_iter = variant.sampleNames.begin(); sample_iter != variant.sampleNames.end(); ++sample_iter){
    if (sample_indices.find(*sample_iter) == sample_indices.end())
      continue;
    int sample_index = sample_indices.find(*sample_iter)->second;
    int gp_index     = 0;
    double total     = 0.0;
    for (unsigned int i = 0; i < num_alleles; ++i){
      for (unsigned int j = 0; j < num_alleles; ++j, ++gp_index){
	  double prob = variant.getSampleValueFloat(PGP_KEY, *sample_iter, gp_index);
	  gp_probs.push_back(std::max(prob, MIN_ALLELE_PRIOR));
	  total += gp_probs.back();
      }
    }

    // Normalize and log-transform priors and store at appropriate index
    gp_index = 0;
    double* log_prior_ptr = log_allele_priors + sample_index;
    for (unsigned int i = 0; i < num_alleles; ++i){
      for (unsigned int j = 0; j < num_alleles; ++j, ++gp_index){
	*log_prior_ptr = log(gp_probs[gp_index]/total);
	log_prior_ptr += num_samples;
      }
    }

    gp_probs.clear();
    sample_count++;
  }
    
  // Ensure that the VCF contained priors for all samples
  if (sample_count != num_samples)
    printErrorAndDie("VCF only contained allele priors for a subset of samples");

  return log_allele_priors;
}
