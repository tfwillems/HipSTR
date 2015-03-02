#include <stdlib.h>

#include "vcflib/src/Variant.h"
#include "SNPTree.h"

bool is_biallelic_snp(vcf::Variant& variant){
  if (variant.alleles.size() != 2)
    return false;

  std::vector<std::string> alleles = variant.alleles;
  for (auto iter = variant.alleles.begin(); iter != variant.alleles.end(); iter++)
    if (iter->size() != 1)
      return false;
  return true;
}


int main(int argc, char** argv) {
  //std::cin.sync_with_stdio(false);
  vcf::VariantCallFile variant_file;
  if (argc > 1){
    std::string filename = argv[1];
    variant_file.open(filename);
  }
  else
    variant_file.open(std::cin);

  assert(variant_file.is_open());
  
  uint32_t start = 0;
  uint32_t end   = 10000000;
  assert(variant_file.setRegion("", start, end));

  // Index samples
  std::map<std::string, unsigned int> sample_indices;
  unsigned int sample_count = 0;
  for (auto sample_iter = variant_file.sampleNames.begin(); sample_iter != variant_file.sampleNames.end(); ++sample_iter)
    sample_indices[*sample_iter] = sample_count++;
  
  // Iterate through all VCF entries
  std::vector< std::vector<SNP> > snps_by_sample(variant_file.sampleNames.size());
  vcf::Variant variant(variant_file);
  uint32_t locus_count = 0;
  while(variant_file.getNextVariant(variant)){
    ++locus_count;
    if (locus_count % 1000 == 0)   std::cout << "\rProcessing locus #" << locus_count << std::flush;
    if(!is_biallelic_snp(variant)) continue;

    for (auto sample_iter = variant.sampleNames.begin(); sample_iter != variant.sampleNames.end(); ++sample_iter){
      std::string gts = variant.getGenotype(*sample_iter);
      assert(gts.size() == 3);
      if (gts[1] == '|'){
	int gt_1 = gts[0]-'0';
	int gt_2 = gts[2]-'0';

	// Only Heterozygous SNPs are informative
	if (gt_1 != gt_2){
	  char a1 = variant.alleles[gt_1][0];
	  char a2 = variant.alleles[gt_2][0];
	  snps_by_sample[sample_indices[*sample_iter]].push_back(SNP(variant.position, a1, a2));
	}
      }
    }
  }

  // Create SNP trees
  std::vector< SNPTree* > snp_trees;
  for(unsigned int i = 0; i < snps_by_sample.size(); i++){
    std::cout << "Building interval tree for " << variant_file.sampleNames[i] << " containing " << snps_by_sample[i].size() << " heterozygous SNPs" << std::endl;
    snp_trees.push_back(new SNPTree(snps_by_sample[i]));
  }
  
  // Discard SNPs
  snps_by_sample.clear();

  // Delete SNP trees
  for (int i = 0; i < snp_trees.size(); i++)
    delete snp_trees[i];

  std::cout << std::endl;

  return 0;
}
