#ifndef SNP_BAM_PROCESSOR_H_
#define SNP_BAM_PROCESSOR_H_

#include <iostream>
#include <string>
#include <vector>

#include "bamtools/include/api/BamAlignment.h"
#include "vcflib/src/Variant.h"

#include "bam_processor.h"
#include "base_quality.h"
#include "error.h"
#include "region.h"

class SNPBamProcessor : public BamProcessor {
private:
  bool have_snp_vcf;
  vcflib::VariantCallFile phased_snp_vcf;
  int32_t match_count_, mismatch_count_;

  // Timing statistics (in seconds)
  double total_snp_phase_info_time_;
  double locus_snp_phase_info_time_;

public:
 SNPBamProcessor(bool use_bam_rgs, bool check_mate_chroms, bool remove_pcr_dups):BamProcessor(use_bam_rgs, check_mate_chroms, remove_pcr_dups){
    have_snp_vcf     = false;
    match_count_     = 0;
    mismatch_count_  = 0;
    total_snp_phase_info_time_  = 0;;
    locus_snp_phase_info_time_  = -1;
  }

  double total_snp_phase_info_time() { return total_snp_phase_info_time_; }
  double locus_snp_phase_info_time() { return locus_snp_phase_info_time_; }

  void process_reads(std::vector< std::vector<BamTools::BamAlignment> >& paired_strs_by_rg,
		     std::vector< std::vector<BamTools::BamAlignment> >& mate_pairs_by_rg,
		     std::vector< std::vector<BamTools::BamAlignment> >& unpaired_strs_by_rg,
		     std::vector<std::string>& rg_names, Region& region, std::string& ref_allele, std::string& chrom_seq,
		     std::ostream& out);


  virtual void analyze_reads_and_phasing(std::vector< std::vector<BamTools::BamAlignment> >& alignments,
					 std::vector< std::vector<double> >& log_p1s, 
					 std::vector< std::vector<double> >& log_p2s,
					 std::vector<std::string>& rg_names, Region& region, std::string& ref_allele, std::string& chrom_seq){
    std::cerr << "Ignoring read phasing probabilties" << std::endl;
  }

  void set_input_snp_vcf(std::string& vcf_file){
    if(!phased_snp_vcf.open(vcf_file))
      printErrorAndDie("Failed to open input SNP VCF file");
    have_snp_vcf = true;
  }

  void finish(){
    std::cerr << "SNP matching statistics: "   << match_count_     << "\t" << mismatch_count_ << std::endl;
  }
};


#endif
