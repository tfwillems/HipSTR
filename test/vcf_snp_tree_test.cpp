#include <stdlib.h>
#include <iostream>

#include "../vcflib/src/Variant.h"

#include "../snp_tree.h"

int main(int argc, char** argv) {
  vcflib::VariantCallFile variant_file;
  if (argc > 1){
    std::string filename = argv[1];
    variant_file.open(filename);
  }
  else
    variant_file.open(std::cin);

  std::string chrom = "22";
  uint32_t start    = 10000000; 
  uint32_t end      = 20000000;
  std::vector<SNPTree*> snp_trees;
  std::map<std::string, unsigned int> sample_indices;
  create_snp_trees(chrom, start, end, 1, 1, variant_file, NULL, sample_indices, snp_trees, std::cerr);
  destroy_snp_trees(snp_trees);
  return 0;
}
