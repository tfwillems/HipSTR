#include <assert.h>

#include <set>

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

bool in_region(vcflib::Variant& variant, uint32_t region_start, uint32_t region_end){
  return variant.position >= region_start && variant.position <= region_end;
}

void filter_snps(std::vector<SNP>& snps, std::set<int32_t>& bad_sites){
  int insert_index = 0;
  for (int i = 0; i < snps.size(); i++)
    if (bad_sites.find(snps[i].pos()+1) == bad_sites.end()) // +1 required b/c bad sites are 1-based, while SNPs are 0-based
      snps[insert_index++] = snps[i];
  snps.resize(insert_index);
}

bool create_snp_trees(const std::string& chrom, uint32_t start, uint32_t end, uint32_t skip_start, uint32_t skip_stop, vcflib::VariantCallFile& variant_file, HaplotypeTracker* tracker,
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

  std::vector< std::set<int32_t> > bad_sites_by_family(tracker != NULL ? tracker->families().size() : 0);

  // Iterate through all VCF entries
  std::vector< std::vector<SNP> > snps_by_sample(variant_file.sampleNames.size());
  vcflib::Variant variant(variant_file);
  uint32_t locus_count = 0, skip_count = 0;
  while (variant_file.getNextVariant(variant)){
    //if (locus_count % 1000 == 0)   std::cout << "\rProcessing locus #" << locus_count << " (skipped " << skip_count << ") at position " << variant.position << std::flush;
    if (!is_biallelic_snp(variant) || in_region(variant, skip_start, skip_stop)){
      skip_count++;
      continue;
    }

    // When performing pedigree-based filtering, we need to identify sites with any Mendelian
    // inconsistencies or missing genotypes as these won't be detected by the haplotype tracker
    if (tracker != NULL){
      const std::vector<NuclearFamily>& families = tracker->families();
      int family_index = 0;
      for (auto family_iter = families.begin(); family_iter != families.end(); ++family_iter, ++family_index)
	if (family_iter->is_missing_genotype(variant) || !family_iter->is_mendelian(variant))
	  bad_sites_by_family[family_index].insert(variant.position);
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

  // Filter out SNPs on a per-sample basis using any available pedigree information
  int MAX_BEST_SCORE = 10;
  int MIN_SECOND_BEST_SCORE = 50;
  if (tracker != NULL){
    int32_t filt_count = 0, unfilt_count = 0;
    const std::vector<NuclearFamily>& families = tracker->families();
    int family_index = 0;
    for (auto family_iter = families.begin(); family_iter != families.end(); ++family_iter, ++family_index){
      std::vector<int> maternal_indices, paternal_indices;
      bool good_haplotypes = tracker->infer_haplotype_inheritance(*family_iter, MAX_BEST_SCORE, MIN_SECOND_BEST_SCORE, maternal_indices, paternal_indices, bad_sites_by_family[family_index]);

      // If the family haplotypes aren't good enough, clear all of the sample's SNPs. Otherwise, remove only the bad sites from each sample's list
      for (auto sample_iter = family_iter->get_samples().begin(); sample_iter != family_iter->get_samples().end(); sample_iter++){
	auto sample_index = sample_indices.find(*sample_iter);
	if (sample_index != sample_indices.end()){
	  filt_count += snps_by_sample[sample_index->second].size();
	  if (!good_haplotypes)
	    snps_by_sample[sample_index->second].clear();
	  else
	    filter_snps(snps_by_sample[sample_index->second], bad_sites_by_family[family_index]);
	  filt_count   -= snps_by_sample[sample_index->second].size();
	  unfilt_count += snps_by_sample[sample_index->second].size();
	}
      }
    }
    logger << "Removed " << filt_count << " out of " << filt_count+unfilt_count << " individual heterozygous SNP calls due to pedigree uncertainties or inconsistencies" << std::endl;
  }
  

  // Create SNP trees
  for (unsigned int i = 0; i < snps_by_sample.size(); i++){
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
