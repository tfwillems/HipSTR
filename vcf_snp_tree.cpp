#include <stdlib.h>

#include "vcflib/src/Variant.h"

#include "snp_tree.h"

int main(int argc, char** argv) {
  //std::cin.sync_with_stdio(false);
  vcf::VariantCallFile variant_file;
  if (argc > 1){
    std::string filename = argv[1];
    variant_file.open(filename);
  }
  else
    variant_file.open(std::cin);

  std::string chrom = "X";
  uint32_t start    = 50000000; 
  uint32_t end      = 75000000;
  std::vector<SNPTree*> snp_trees;
  std::map<std::string, unsigned int> sample_indices;
  create_snp_trees(chrom, start, end, variant_file, sample_indices, snp_trees);
  destroy_snp_trees(snp_trees);
  return 0;
}
