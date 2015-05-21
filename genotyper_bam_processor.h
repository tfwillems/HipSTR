#ifndef GENOTYPER_BAM_PROCESSOR_H_
#define GENOTYPER_BAM_PROCESSOR_H_

#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "bamtools/include/api/BamAlignment.h"
#include "vcflib/src/Variant.h"

#include "em_stutter_genotyper.h"
#include "region.h"
#include "seq_stutter_genotyper.h"
#include "snp_bam_processor.h"
#include "stutter_model.h"
#include "SeqAlignment/HTMLCreator.h"

class GenotyperBamProcessor : public SNPBamProcessor {
private:  
  // Counters for EM convergence
  int num_em_converge_, num_em_fail_;

  // Parameters for stutter models read from file
  bool read_stutter_models_;
  std::map<Region, StutterModel*> stutter_models_;

  // Output file for stutter models
  bool output_stutter_models_;
  std::ofstream stutter_model_out_;

  // Outupt file for STR alleles (w/o genotypes)
  bool output_alleles_;
  std::ofstream allele_vcf_;

  // Output file for STR genotypes
  bool output_str_gts_;
  std::ofstream str_vcf_;
  std::vector<std::string> samples_to_genotype_;

  // Flag for type of genotyper to use
  bool use_seq_aligner_;

  // Counters for genotyping success;
  int num_genotype_success_, num_genotype_fail_;

  // VCF containg SNP and STR genotypes for a reference panel
  bool have_ref_vcf_;
  vcf::VariantCallFile ref_vcf_;

  bool output_viz_;
  std::ofstream viz_out_;

public:
 GenotyperBamProcessor(bool use_lobstr_rg, bool check_mate_chroms, bool use_seq_aligner):SNPBamProcessor(use_lobstr_rg, check_mate_chroms){
    output_stutter_models_ = false;
    output_alleles_        = false;
    output_str_gts_        = false;
    read_stutter_models_   = false;
    have_ref_vcf_          = false;
    use_seq_aligner_       = use_seq_aligner;
    num_em_converge_       = 0;
    num_em_fail_           = 0;
    num_genotype_success_  = 0;
    num_genotype_fail_     = 0;
    MAX_EM_ITER            = 100;
    ABS_LL_CONVERGE        = 0.01;
    FRAC_LL_CONVERGE       = 0.001;
    MIN_TOTAL_READS        = 100;
    
    // Exploratory
    output_viz_ = true;
    viz_out_.open("boo.html", std::ofstream::out);
    writeHeader(viz_out_);
  }

  ~GenotyperBamProcessor(){
    for (auto iter = stutter_models_.begin(); iter != stutter_models_.end(); iter++)
      delete iter->second;
    stutter_models_.clear();
  }

  void use_seq_aligner(){
    use_seq_aligner_ = true;
  }

  void set_ref_vcf(std::string& ref_vcf_file){
    if(!ref_vcf_.open(ref_vcf_file))
      printErrorAndDie("Failed to open VCF file for reference panel");
    have_ref_vcf_ = true;
  }

  void set_input_stutter(std::string& model_file){
    std::ifstream input;
    input.open(model_file, std::ifstream::in);
    if (!input.is_open())
      printErrorAndDie("Failed to open input file for stutter models");
    StutterModel::read_models(input, stutter_models_);
    std::cerr << stutter_models_.size() << std::endl;
    read_stutter_models_ = true;
    input.close();
  }
  
  void set_output_stutter(std::string& model_file){
    output_stutter_models_ = true;
    stutter_model_out_.open(model_file, std::ofstream::out);
    if (!stutter_model_out_.is_open())
      printErrorAndDie("Failed to open output file for stutter models");
  }

  void set_output_allele_vcf(std::string& vcf_file){
    output_alleles_ = true;
    allele_vcf_.open(vcf_file, std::ofstream::out);
    if (!allele_vcf_.is_open())
      printErrorAndDie("Failed to open VCF file for STR alleles");
    
    std::vector<std::string> no_samples;
    if (use_seq_aligner_)
      SeqStutterGenotyper::write_vcf_header(no_samples, allele_vcf_);
    else
      printErrorAndDie("Cannot output STR allele VCF if --seq-genotyper option has not been specified");
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
    if (use_seq_aligner_)
      SeqStutterGenotyper::write_vcf_header(samples_to_genotype_, str_vcf_);
    else
      EMStutterGenotyper::write_vcf_header(samples_to_genotype_, str_vcf_);
  }

  void analyze_reads_and_phasing(std::vector< std::vector<BamTools::BamAlignment> >& alignments,
				 std::vector< std::vector<double> >& log_p1s,
				 std::vector< std::vector<double> >& log_p2s,
				 std::vector<std::string>& rg_names, Region& region, std::string& ref_allele, std::string& chrom_seq);
  void finish(){
    SNPBamProcessor::finish();
    if (output_alleles_)
      allele_vcf_.close();
    if (output_str_gts_)
      str_vcf_.close();
    if (output_stutter_models_)
      stutter_model_out_.close();
    if (output_viz_)
      viz_out_.close();

    std::cerr << "Stutter model training succeeded for " << num_em_converge_ << " out of " << num_em_converge_+num_em_fail_ << " loci" << std::endl;
    std::cerr << "Genotyping succeeded for " << num_genotype_success_ << " out of " << num_genotype_success_+num_genotype_fail_ << " loci" << std::endl;
  }

  // EM parameters for length-based stutter learning
  int MAX_EM_ITER;
  double ABS_LL_CONVERGE;  // For EM convergence, new_LL - prev_LL < ABS_LL_CONVERGE
  double FRAC_LL_CONVERGE; // For EM convergence, -(new_LL-prev_LL)/prev_LL < FRAC_LL_CONVERGE
  int MIN_TOTAL_READS;     // Minimum total reads required to genotype locus
};



#endif
