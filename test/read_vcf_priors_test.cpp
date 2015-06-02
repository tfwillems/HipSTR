#include <iostream>
#include <map>
#include <vector>

#include "../vcflib/src/Variant.h"

#include "../error.h"
#include "../region.h"
#include "../vcf_input.h"

int main(int argc, char* argv[]){
  if (argc != 3)
    printErrorAndDie("Script requires exactly 2 arguments");
  std::string region_file = std::string(argv[1]);
  std::string vcf_file    = std::string(argv[2]);

  // Read list of regions
  std::vector<Region> regions;  
  readRegions(region_file, regions, 1000);

  vcflib::VariantCallFile ref_vcf;
  if(!ref_vcf.open(vcf_file))
    printErrorAndDie("Failed to open VCF");

  // Populate map with samples in VCF header
  std::map<std::string, int> sample_indices;
  for (unsigned int i = 0; i < ref_vcf.sampleNames.size(); i++)
    sample_indices[ref_vcf.sampleNames[i]] = i;

  std::vector<std::string> alleles;
  std::vector<bool> got_priors;
  int32_t pos;
  for (unsigned int i = 0; i < regions.size(); i++){
    bool success;
    double* priors = extract_vcf_alleles_and_log_priors(&ref_vcf, &(regions[i]), sample_indices, alleles, got_priors, pos, success);

    if (success){
      std::cerr << "Position=" << pos << std::endl;
      std::cerr << "Alleles:" << std::endl;
      for (unsigned int j = 0; j < alleles.size(); j++)
	std::cerr << alleles[j] << std::endl;
    }
    else {
      std::cerr << "Failed to read alleles and priors for region " << regions[i].str() << std::endl;
    }
      
    alleles.clear();
    got_priors.clear();
    delete [] priors;
  }

}
