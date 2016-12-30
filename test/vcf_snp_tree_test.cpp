#include <stdlib.h>
#include <iostream>

#include "../src/region.h"
#include "../src/vcf_reader.h"
#include "../src/snp_tree.h"

int main(int argc, char** argv) {
  std::string filename = argv[1];
  VCF::VCFReader vcf_reader(filename);

  std::string chrom = "22";
  uint32_t start    = 10000000; 
  uint32_t end      = 20000000;
  std::vector<SNPTree*> snp_trees;
  std::map<std::string, unsigned int> sample_indices;
  std::vector<Region> skip_regions;
  int32_t skip_pad = 0;
  create_snp_trees(chrom, start, end, skip_regions, skip_pad, &vcf_reader, NULL, sample_indices, snp_trees, std::cerr);
  destroy_snp_trees(snp_trees);
  return 0;
}
