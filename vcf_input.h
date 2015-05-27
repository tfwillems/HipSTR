#ifndef VCF_INPUT_H_
#define VCF_INPUT_H_

#include <map>
#include <string>
#include <vector>

#include "vcflib/src/Variant.h"

#include "region.h"

extern std::string PGP_KEY;
extern std::string START_INFO_TAG;
extern std::string STOP_INFO_TAG;
extern const double MIN_ALLELE_PRIOR;
extern const int32_t pad;

void read_vcf_alleles(vcf::VariantCallFile* ref_vcf, Region* region, std::vector<std::string>& alleles, int32_t& pos, bool& success);


double* extract_vcf_alleles_and_log_priors(vcf::VariantCallFile* ref_vcf, Region* region, std::map<std::string, int>& sample_indices,
					   std::vector<std::string>& alleles, std::vector<bool>& got_priors, int32_t& pos, bool& success);
 
#endif
