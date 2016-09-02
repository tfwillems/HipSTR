#include <iostream>
#include <map>
#include <vector>

#include "../error.h"
#include "../region.h"
#include "../vcf_input.h"
#include "../vcf_reader.h"

int main(int argc, char* argv[]){
  if (argc != 3)
    printErrorAndDie("Script requires exactly 2 arguments");
  std::string region_file = std::string(argv[1]);
  std::string vcf_file    = std::string(argv[2]);

  // Read list of regions
  std::vector<Region> regions;  
  readRegions(region_file, regions, 1000, "", std::cerr);

  VCF::VCFReader ref_vcf(vcf_file);

  // Populate map with samples in VCF header
  const std::vector<std::string>& samples = ref_vcf.get_samples();
  std::map<std::string, int> sample_indices;
  for (unsigned int i = 0; i < samples.size(); i++)
    sample_indices[samples[i]] = i;

  std::vector<std::string> alleles;
  std::vector<bool> got_priors;
  int32_t pos;
  for (unsigned int i = 0; i < regions.size(); i++){
    bool success;
    double* priors = extract_vcf_alleles_and_log_priors(&ref_vcf, &(regions[i]), sample_indices, alleles, got_priors, pos, success, std::cerr);

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
