#ifndef LENGTH_EM_BAM_PROCESSOR_H_
#define LENGTH_EM_BAM_PROCESSOR_H_

#include <fstream>
#include <string>
#include <vector>

#include "bamtools/include/api/BamAlignment.h"
#include "snp_bam_processor.h"

class LengthEMBamProcessor : public SNPBamProcessor {
private:  
  // Counters for EM convergence
  int num_em_converge_, num_em_fail_;

  // Output file for STR genotypes
  bool output_str_gts_;
  std::ofstream str_vcf_;
  std::vector<std::string> samples_to_genotype_;

public:
 LengthEMBamProcessor(bool use_lobstr_rg, bool check_mate_chroms):SNPBamProcessor(use_lobstr_rg, check_mate_chroms){
    output_str_gts_ = false;
  }

  void set_output_str_vcf(std::string& vcf_file, std::set<std::string>& samples_to_output){
    output_str_gts_ = true;
    str_vcf_.open(vcf_file, std::ofstream::out);
    if (!str_vcf_.is_open())
      printErrorAndDie("Failed to open VCF file for STR genotypes");
    
    // Print floats with exactly 3 decimal places
    str_vcf_.precision(3);
    str_vcf_.setf(std::ios::fixed, std::ios::floatfield);
    
    // Assemble a list of sample names for genotype output
    std::copy(samples_to_output.begin(), samples_to_output.end(), std::back_inserter(samples_to_genotype_));
    std::sort(samples_to_genotype_.begin(), samples_to_genotype_.end());
    
    // Write VCF header
    EMStutterGenotyper::write_vcf_header(samples_to_genotype_, str_vcf_);
  }

  void analyze_reads_and_phasing(std::vector< std::vector<BamTools::BamAlignment> >& alignments,
				 std::vector< std::vector<double> >& log_p1s,
				 std::vector< std::vector<double> >& log_p2s,
				 std::vector<std::string>& rg_names, Region& region, std::string& ref_allele, std::string& chrom_seq);
  void finish(){
    SNPBamProcessor::finish();
    if (output_str_gts_)
      str_vcf_.close();  
  }

  int MAX_EM_ITER         = 100;
  double ABS_LL_CONVERGE  = 0.01;  // For EM convergence, new_LL - prev_LL < ABS_LL_CONVERGE
  double FRAC_LL_CONVERGE = 0.001; // For EM convergence, -(new_LL-prev_LL)/prev_LL < FRAC_LL_CONVERGE
};



#endif
