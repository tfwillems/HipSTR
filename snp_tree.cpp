#include <assert.h>

#include "error.h"
#include "snp_tree.h"

std::ostream& operator<<(std::ostream& out, SNP& snp) {
  out << "SNP: (" << snp.pos() << ", " << snp.base_one() << "|" << snp.base_two() << ")";
  return out;
}

bool is_biallelic_snp(vcflib::Variant& variant){
  if (variant.alleles.size() != 2)
    return false;
  for (auto iter = variant.alleles.begin(); iter != variant.alleles.end(); ++iter)
    if (iter->size() != 1)
      return false;
  return true;
}

bool create_snp_trees(const std::string& chrom, uint32_t start, uint32_t end, vcflib::VariantCallFile& variant_file,
                      std::map<std::string, unsigned int>& sample_indices, std::vector<SNPTree*>& snp_trees, std::ostream& logger){
  logger << "Building SNP tree for region " << chrom << ":" << start << "-" << end << std::endl;
  assert(sample_indices.size() == 0 && snp_trees.size() == 0);
  assert(variant_file.is_open());

  if (!variant_file.setRegion(chrom, start, end)){
    // Retry setting region if chr is in chromosome name
    if (chrom.size() <= 3 || chrom.substr(0, 3).compare("chr") != 0 || !variant_file.setRegion(chrom.substr(3), start, end))
      return false;
  }

  // Index samples
  unsigned int sample_count = 0;
  for (auto sample_iter = variant_file.sampleNames.begin(); sample_iter != variant_file.sampleNames.end(); ++sample_iter)
    sample_indices[*sample_iter] = sample_count++;

  // Iterate through all VCF entries
  std::vector< std::vector<SNP> > snps_by_sample(variant_file.sampleNames.size());
  vcflib::Variant variant(variant_file);
  uint32_t locus_count = 0, skip_count = 0;
  while(variant_file.getNextVariant(variant)){
    //if (locus_count % 1000 == 0)   std::cout << "\rProcessing locus #" << locus_count << " (skipped " << skip_count << ") at position " << variant.position << std::flush;
    if(!is_biallelic_snp(variant)){
      skip_count++;
      continue;
    }
    ++locus_count;
    for (auto sample_iter = variant.sampleNames.begin(); sample_iter != variant.sampleNames.end(); ++sample_iter){
      std::string gts = variant.getGenotype(*sample_iter);
      if (gts.size() == 0)
	continue;
      assert(gts.size() == 3);
      if (gts[1] == '|'){
        int gt_1 = gts[0]-'0';
        int gt_2 = gts[2]-'0';

        // Only Heterozygous SNPs are informative
        if (gt_1 != gt_2){
          char a1 = variant.alleles[gt_1][0];
          char a2 = variant.alleles[gt_2][0];
	  
	  // IMPORTANT NOTE: VCFs are 1-based, but BAMs are 0-based. Decrease VCF coordinate by 1 for consistency
          snps_by_sample[sample_indices[*sample_iter]].push_back(SNP(variant.position-1, a1, a2)); 
        }
      }
      else {
	//printErrorAndDie("SNP panel VCF must contain phased genotypes and therefore utilize the | genotype separator");
      }
    }
  }
  logger << "Region contained a total of " << locus_count << " valid SNPs" << std::endl;
  
  // Create SNP trees
  for(unsigned int i = 0; i < snps_by_sample.size(); i++){
    //logger << "Building interval tree for " << variant_file.sampleNames[i] << " containing " << snps_by_sample[i].size() << " heterozygous SNPs" << std::endl;
    snp_trees.push_back(new SNPTree(snps_by_sample[i]));
  }

  // Discard SNPs
  snps_by_sample.clear();
  
  return true;
}

void destroy_snp_trees(std::vector<SNPTree*>& snp_trees){
  for (unsigned int i = 0; i < snp_trees.size(); i++)
    delete snp_trees[i];
  snp_trees.clear();
}
