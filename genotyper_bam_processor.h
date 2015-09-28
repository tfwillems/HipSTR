#ifndef GENOTYPER_BAM_PROCESSOR_H_
#define GENOTYPER_BAM_PROCESSOR_H_

#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "bamtools/include/api/BamAlignment.h"
#include "vcflib/src/Variant.h"
#include "bgzf_streams.h"

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

  // Output file for STR alleles (w/o genotypes)
  bool output_alleles_;
  bgzfostream allele_vcf_;

  // Output file for STR genotypes
  bool output_str_gts_;
  bgzfostream str_vcf_;
  std::vector<std::string> samples_to_genotype_;

  // Flag for type of genotyper to use
  bool use_seq_aligner_;

  // Counters for genotyping success;
  int num_genotype_success_, num_genotype_fail_;

  // VCF containg SNP and STR genotypes for a reference panel
  bool have_ref_vcf_;
  vcflib::VariantCallFile ref_vcf_;

  bool output_viz_;
  bgzfostream viz_out_;

  bool output_bstrap_quals_;    // Output the BQ FORMAT field to the VCF
  bool output_gls_;             // Output the GL FORMAT field to the VCF
  bool output_pls_;             // Output the PL FORMAT field to the VCF
  bool output_all_reads_;       // Output the ALLREADS  FORMAT field to the VCF
  bool output_pall_reads_;      // Output the PALLREADS FORMAT field to the VCF
  bool output_mall_reads_;      // Output the MALLREADS FORMAT field to the VCF

  std::set<std::string> haploid_chroms_;

  // Timing statistics (in seconds)
  double total_stutter_time_,  locus_stutter_time_;
  double total_genotype_time_, locus_genotype_time_;

  // True iff we should recalculate the stutter model after performing haplotype alignments
  // The idea is that the haplotype-based alignments should be far more accurate, and reperforming
  // the stutter analysis will result in a better stutter model
  bool recalc_stutter_model_;

  // If this flag is set, HTML alignments are written for both the haplotype alignments and Needleman-Wunsch left alignments
  bool viz_left_alns_;

  // If true, the seqeunce-based genotyper will collapse reads with identical sequences
  // and merge their base quality scores. Results in a large reduction in computation time
  bool pool_seqs_;

public:
 GenotyperBamProcessor(bool use_bam_rgs, bool check_mate_chroms, bool remove_pcr_dups,
		       bool use_seq_aligner):SNPBamProcessor(use_bam_rgs, check_mate_chroms, remove_pcr_dups){
    output_stutter_models_ = false;
    output_alleles_        = false;
    output_str_gts_        = false;
    output_viz_            = false;
    read_stutter_models_   = false;
    have_ref_vcf_          = false;
    viz_left_alns_         = false;
    pool_seqs_             = false;
    use_seq_aligner_       = use_seq_aligner;
    haploid_chroms_        = std::set<std::string>();
    num_em_converge_       = 0;
    num_em_fail_           = 0;
    num_genotype_success_  = 0;
    num_genotype_fail_     = 0;
    MAX_EM_ITER            = 100;
    ABS_LL_CONVERGE        = 0.01;
    FRAC_LL_CONVERGE       = 0.001;
    MIN_TOTAL_READS        = 100;
    MAX_TOTAL_READS        = 25000;
    output_bstrap_quals_   = true;
    output_gls_            = false;
    output_pls_            = false;
    output_all_reads_      = true;
    output_pall_reads_     = true;
    output_mall_reads_     = true;
    total_stutter_time_    = 0;
    locus_stutter_time_    = -1;
    total_genotype_time_   = 0;
    locus_genotype_time_   = -1;
    recalc_stutter_model_  = false;
  }

  ~GenotyperBamProcessor(){
    for (auto iter = stutter_models_.begin(); iter != stutter_models_.end(); iter++)
      delete iter->second;
    stutter_models_.clear();
  }

  double total_stutter_time()  { return total_stutter_time_;  }
  double locus_stutter_time()  { return locus_stutter_time_;  }
  double total_genotype_time() { return total_genotype_time_; }
  double locus_genotype_time() { return locus_genotype_time_; }

  void output_gls(){
    output_gls_ = true;
  }

  void output_pls(){
    output_pls_ = true;
  }

  void hide_all_reads(){
    output_all_reads_ = false;
  }

  void hide_pall_reads(){
    output_pall_reads_ = false;
  }

  void hide_mall_reads(){
    output_mall_reads_ = false;
  }

  void add_haploid_chrom(std::string chrom){
    haploid_chroms_.insert(chrom);
  }

  void use_seq_aligner(){
    use_seq_aligner_ = true;
  }

  void use_len_model(){
    use_seq_aligner_ = false;
  }

  void visualize_left_alns(){
    viz_left_alns_ = true;
  }

  void pool_sequences(){
    pool_seqs_ = true;
  }

  void set_output_viz(std::string& viz_file){
    output_viz_ = true;
    viz_out_.open(viz_file.c_str());
    writeHeader(viz_out_);
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
      printErrorAndDie("Failed to open input file for stutter models. Filename = " + model_file);
    StutterModel::read_models(input, stutter_models_);
    log("Read stutter models for " + std::to_string(stutter_models_.size()) + " loci");
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
    allele_vcf_.open(vcf_file.c_str());
    std::vector<std::string> no_samples;
    if (use_seq_aligner_)
      SeqStutterGenotyper::write_vcf_header(no_samples, output_gls_, output_pls_, allele_vcf_);
    else
      printErrorAndDie("Cannot output STR allele VCF if --seq-genotyper option has not been specified");
  }

  void set_output_str_vcf(std::string& vcf_file, std::set<std::string>& samples_to_output){
    output_str_gts_ = true;
    str_vcf_.open(vcf_file.c_str());

    // Print floats with exactly 3 decimal places
    str_vcf_.precision(3);
    str_vcf_.setf(std::ios::fixed, std::ios::floatfield);
    
    // Assemble a list of sample names for genotype output
    std::copy(samples_to_output.begin(), samples_to_output.end(), std::back_inserter(samples_to_genotype_));
    std::sort(samples_to_genotype_.begin(), samples_to_genotype_.end());
    
    // Write VCF header
    if (use_seq_aligner_)
      SeqStutterGenotyper::write_vcf_header(samples_to_genotype_, output_gls_, output_pls_, str_vcf_);
    else
      EMStutterGenotyper::write_vcf_header(samples_to_genotype_, output_gls_, output_pls_, str_vcf_);
  }

  void analyze_reads_and_phasing(std::vector< std::vector<BamTools::BamAlignment> >& alignments,
				 std::vector< std::vector<double> >& log_p1s,
				 std::vector< std::vector<double> >& log_p2s,
				 std::vector<std::string>& rg_names, Region& region, std::string& ref_allele, std::string& chrom_seq, int iter);
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

    log("Stutter model training succeeded for " + std::to_string(num_em_converge_) + " out of " + std::to_string(num_em_converge_+num_em_fail_) + " loci");
    log("Genotyping succeeded for " + std::to_string(num_genotype_success_) + " out of " + std::to_string(num_genotype_success_+num_genotype_fail_) + " loci");
  }

  // EM parameters for length-based stutter learning
  int MAX_EM_ITER;
  double ABS_LL_CONVERGE;  // For EM convergence, new_LL - prev_LL < ABS_LL_CONVERGE
  double FRAC_LL_CONVERGE; // For EM convergence, -(new_LL-prev_LL)/prev_LL < FRAC_LL_CONVERGE
  int32_t MIN_TOTAL_READS; // Minimum total reads required to genotype locus
  int32_t MAX_TOTAL_READS; // Don't genotype loci where the number of STR reads exceeds this limit
};



#endif
