#include <string>
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

  std::vector<std::string> alleles;
  int32_t pos;
  for (unsigned int i = 0; i < regions.size(); i++){
    bool success;
    read_vcf_alleles(&ref_vcf, &regions[i], alleles, pos, success);
    if (success){
      std::cerr << "Position=" << pos << std::endl;
      std::cerr << "Alleles:" << std::endl;
      for (unsigned int j = 0; j < alleles.size(); j++)
	std::cerr << alleles[j] << std::endl;
    }
    else {
      std::cerr << "Failed to read alleles for region " << regions[i].str() << std::endl; 
    }
    alleles.clear();
  }
}
