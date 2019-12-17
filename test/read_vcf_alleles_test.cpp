#include <iostream>
#include <string>
#include <vector>

#include "../src/error.h"
#include "../src/region.h"
#include "../src/vcf_input.h"
#include "../src/vcf_reader.h"

int main(int argc, char* argv[]){
  if (argc != 3)
    printErrorAndDie("Script requires exactly 2 arguments");
  std::string region_file = std::string(argv[1]);
  std::string vcf_file    = std::string(argv[2]);

  // Read list of regions
  std::vector<Region> regions;  
  readRegions(region_file, 1000, "", regions, std::cerr);

  VCF::VCFReader ref_vcf(vcf_file);
  std::vector<std::string> alleles, lflanks, rflanks;
  int32_t pos;
  for (unsigned int i = 0; i < regions.size(); i++){
    bool success = read_vcf_alleles(&ref_vcf, regions[i], alleles, lflanks, rflanks, pos);
    if (success){
      std::cerr << "Position=" << pos << std::endl;
      std::cerr << "Alleles:" << std::endl;
      for (unsigned int j = 0; j < alleles.size(); j++)
	std::cerr << alleles[j] << std::endl;

      if (!lflanks.empty()){
	std::cerr << "Left Flanks:" << std::endl;
	for (unsigned int j = 0; j < lflanks.size(); j++)
	  std::cerr << lflanks[j] << std::endl;
      }

      if (!rflanks.empty()){
	std::cerr << "Right Flanks:" << std::endl;
	for (unsigned int j = 0; j < rflanks.size(); j++)
	  std::cerr << rflanks[j] << std::endl;
      }
    }
    else {
      std::cerr << "Failed to read alleles for region " << regions[i].str() << std::endl; 
    }
    alleles.clear();
    std::cerr << std::endl << std::endl;
  }
}
