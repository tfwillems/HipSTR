#ifndef DENOVO_SCANNER_H_
#define DENOVO_SCANNER_H_

#include <vector>
#include <set>
#include <string>

#include "pedigree.h"
#include "vcflib/src/Variant.h"

class DenovoScanner {
 private:
  int32_t window_size_;
  std::vector<NuclearFamily> families_;

 public:
  DenovoScanner(std::vector<NuclearFamily>& families){
    families_    = families;
    window_size_ = 500000;
  }
  
  void scan(vcflib::VariantCallFile& snp_vcf, vcflib::VariantCallFile& str_vcf, std::set<std::string>& sites_to_skip,
	    std::ostream& logger);



};

#endif
