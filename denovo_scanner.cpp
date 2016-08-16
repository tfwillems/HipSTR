#include "denovo_scanner.h"
#include "error.h"
#include "haplotype_tracker.h"
#include "vcf_input.h"

#include <vector>

void DenovoScanner::scan(vcflib::VariantCallFile& snp_vcf, vcflib::VariantCallFile& str_vcf, std::set<std::string>& sites_to_skip,
			 std::ostream& logger){
  HaplotypeTracker haplotype_tracker(families_);
  vcflib::Variant snp_variant(snp_vcf), str_variant(str_vcf);

  std::string chrom = "";
  int32_t num_strs = 0;
  while (str_vcf.getNextVariant(str_variant)){
    num_strs++;
    //PhasedGL phased_gls(str_vcf, str_variant);

    if (str_variant.sequenceName.compare(chrom) != 0){
      chrom = str_variant.sequenceName;
      haplotype_tracker.reset();
      if(!snp_vcf.setRegion(chrom, 1))
	printErrorAndDie("Failed to set the region to chromosome " + chrom + " in the SNP VCF. Please check the SNP VCF and rerun the analysis");
    }

    int32_t start_of_window = str_variant.position - window_size_;
    int32_t end_of_window   = str_variant.position + window_size_;
    if (start_of_window < 0)
      start_of_window = 0;

    // Incorporate new SNPs within the window
    while (haplotype_tracker.last_snp_position() < end_of_window && snp_vcf.getNextVariant(snp_variant)){
      std::string key = snp_variant.sequenceName + ":" + std::to_string(snp_variant.position);
      if (sites_to_skip.find(key) != sites_to_skip.end())
	continue;
      haplotype_tracker.add_snp(snp_variant);
    }

    // Remove SNPs to left of window
    while (haplotype_tracker.next_snp_position() < start_of_window && haplotype_tracker.next_snp_position() != -1)
      haplotype_tracker.remove_next_snp();

    std::cerr << str_variant.position << " " << haplotype_tracker.num_stored_snps()  << std::endl;  // TO DO: Add this to VCF INFO field?

    // Analyze edit distances between the phased SNP haplotypes of each child and its parents
    int d11, d12, d21, d22;
    for (auto family_iter = families_.begin(); family_iter != families_.end(); family_iter++){
      for (auto child_iter = family_iter->get_children().begin(); child_iter != family_iter->get_children().end(); child_iter++){
	DiploidEditDistance maternal_distance = haplotype_tracker.edit_distances(*child_iter, family_iter->get_mother());
	std::cout << *child_iter << maternal_distance;
	//std::cout << *child_iter << " " << str_variant.position << " " << d11 << " " << d12 << " " << d21 << " " << d22 << "\t";

	DiploidEditDistance paternal_distance = haplotype_tracker.edit_distances(*child_iter, family_iter->get_father());
	std::cout << paternal_distance << std::endl;
	//std::cout << d11 << " " << d12 << " " << d21 << " " << d22 << std::endl;
      }
    }
  }
}
