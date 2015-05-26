#include "../vcflib/src/Variant.h"

#include "../error.h"
#include "../region.h"
#include "../vcf_input.h"

int main(){
  std::string region_file = "input/chr1_regions.bed", vcf_file = "input/1kg.chr1.imputed.vcf.gz";

  // Read list of regions
  std::vector<Region> regions;  
  readRegions(region_file, regions, 1000);

  vcf::VariantCallFile ref_vcf;
  if(!ref_vcf.open(vcf_file))
    printErrorAndDie("Failed to open VCF");

  std::vector<std::string> alleles;
  int32_t pos;
  for (unsigned int i = 0; i < regions.size(); i++){
    read_vcf_alleles(&ref_vcf, &regions[i], alleles, pos);
    std::cerr << "Position=" << pos << std::endl;
    std::cerr << "Alleles:" << std::endl;
    for (unsigned int j = 0; j < alleles.size(); j++)
      std::cerr << alleles[j] << std::endl;
    alleles.clear();
  }
}
