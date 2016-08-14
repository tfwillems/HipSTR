#include "denovo_scanner.h"
#include "haplotype_tracker.h"
#include "vcf_input.h"

void DenovoScanner::scan(vcflib::VariantCallFile& snp_vcf, vcflib::VariantCallFile& str_vcf, std::set<std::string>& sites_to_skip,
			 std::ostream& logger){
  HaplotypeTracker haplotype_tracker(families_);
  vcflib::Variant variant(snp_vcf);



  PhasedGL phased_gls;
  // Iterate through the SNP VCF to determine haplotype sharing at each position

  int32_t count = 0;
  while (snp_vcf.getNextVariant(variant)){
    std::string key = variant.sequenceName + ":" + std::to_string(variant.position);
    if (sites_to_skip.find(key) != sites_to_skip.end())
      continue;

    haplotype_tracker.add_snp(variant);

    if (++count % 1000 == 0){
      std::cerr << variant.position << std::endl;
      int32_t position = variant.position;
      while (haplotype_tracker.next_snp_position() < position-window_size_ && haplotype_tracker.next_snp_position() != -1)
        haplotype_tracker.remove_next_snp();
      std::cerr << haplotype_tracker.num_stored_snps()  << std::endl;

      int d11, d12, d21, d22;
      for (auto family_iter = families_.begin(); family_iter != families_.end(); family_iter++){
        for (auto child_iter = family_iter->get_children().begin(); child_iter != family_iter->get_children().end(); child_iter++){
          haplotype_tracker.edit_distances(*child_iter, family_iter->get_mother(), d11, d12, d21, d22);
	  std::cout << *child_iter << " " << position << " " << d11 << " " << d12 << " " << d21 << " " << d22 << "\t";

          haplotype_tracker.edit_distances(*child_iter, family_iter->get_father(), d11, d12, d21, d22);
	  std::cout << d11 << " " << d12 << " " << d21 << " " << d22 << std::endl;
        }
      }
    }
  }


}
