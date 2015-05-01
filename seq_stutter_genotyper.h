#ifndef SEQ_STUTTER_GENOTYPER_H_
#define SEQ_STUTTER_GENOTYPER_H_

#include <algorithm>
#include <iostream>
#include <map>
#include <math.h>
#include <set>
#include <string>
#include <vector>

#include "bamtools/include/api/BamAlignment.h"
#include "vcflib/src/Variant.h"

#include "base_quality.h"
#include "region.h"
#include "stutter_model.h"

#include "SeqAlignment/AlignmentData.h"
#include "SeqAlignment/Haplotype.h"
#include "SeqAlignment/HapBlock.h"

class SeqStutterGenotyper{
 private:
  std::string END_KEY;

  // Locus information
  Region* region_;

  unsigned int num_reads_;   // Total number of reads across all samples
  int num_samples_; // Total number of samples
  int motif_len_;   // # bp in STR motif
  int num_alleles_; // Number of valid alleles

  int MAX_REF_FLANK_LEN;

  double* log_p1_;                 // Log of SNP phasing likelihoods for each read
  double* log_p2_;
  int* sample_label_;              // Sample index for each read
  StutterModel* stutter_model_;
  BaseQuality base_quality_;

  std::vector<int> bp_diffs_;                  // Base pair difference of each read from reference
  std::vector< std::vector<Alignment> > alns_; // Vector of left-aligned alignments  
  std::vector<std::string> sample_names_;      // List of sample names
  std::map<std::string, int> sample_indices_;  // Mapping from sample name to index
  std::vector<int> reads_per_sample_;          // Number of reads for each sample
  std::vector<HapBlock*> hap_blocks_;          // Haplotype blocks
  Haplotype* haplotype_;                       // Potential STR haplotypes

  bool alleles_from_bams_; // Flag that determines if we examine BAMs for candidate alleles

  std::vector<std::string> alleles_; // Vector of indexed alleles
  int32_t pos_;                      // Position of reported alleles in VCF     

  // Iterates through reads and then alleles by their indices
  double* log_aln_probs_;

  // Iterates through allele_1, allele_2 and then samples by their indices
  double* log_sample_posteriors_; 
  
  // Iterates through allele_1, allele_2 and then samples by their indices
  // Only used if per-allele priors have been specified for each sample
  double* log_allele_priors_;

  // VCF containing STR and SNP genotypes for a reference panel
  vcf::VariantCallFile* ref_vcf_;

  /* Compute the alignment probabilites between each read and each haplotype */
  double calc_align_probs();

  /* Compute the posteriors for each sample, given the haplotype probabilites and stutter model */
  double calc_log_sample_posteriors();  

  // Set up the relevant data structures. Invoked by the constructor 
  void init(std::vector< std::vector<BamTools::BamAlignment> >& alignments,
	    std::vector< std::vector<double> >& log_p1,
	    std::vector< std::vector<double> >& log_p2,
	    std::vector<std::string>& sample_names, std::string& chrom_seq);

  // Extract the sequences for each allele and the VCF start position
  void get_alleles(std::string& chrom_seq, std::vector<std::string>& alleles);

  void debug_sample(int sample_index);
  
  void read_ref_vcf_alleles(std::vector<std::string>& alleles);
  

 public:
  SeqStutterGenotyper(Region& region,
		      std::vector< std::vector<BamTools::BamAlignment> >& alignments,
		      std::vector< std::vector<double> >& log_p1, 
		      std::vector< std::vector<double> >& log_p2, 
		      std::vector<std::string>& sample_names, std::string& chrom_seq, 
		      StutterModel& stutter_model, vcf::VariantCallFile* ref_vcf){
    assert(alignments.size() == log_p1.size() && alignments.size() == log_p2.size() && alignments.size() == sample_names.size());
    log_p1_                = NULL;
    log_p2_                = NULL;
    log_aln_probs_         = NULL;
    log_sample_posteriors_ = NULL;
    log_allele_priors_     = NULL;
    sample_label_          = NULL;
    haplotype_             = NULL;
    MAX_REF_FLANK_LEN      = 30;
    END_KEY                = "END";
    pos_                   = -1;
    
    region_       = region.copy();
    num_samples_  = alignments.size();
    sample_names_ = sample_names;
    for (unsigned int i = 0; i < sample_names.size(); i++){
      sample_indices_.insert(std::pair<std::string,int>(sample_names[i], i));
      reads_per_sample_.push_back(alignments[i].size());
    }
    stutter_model_      = stutter_model.copy();
    ref_vcf_            = ref_vcf;
    alleles_from_bams_  = true;
    init(alignments, log_p1, log_p2, sample_names, chrom_seq);
  }

  /*
  SeqStutterGenotyper(Region& region,
		      std::vector< std::vector<BamTools::BamAlignment> >& alignments,
		      std::vector< std::vector<double> >& log_p1, 
		      std::vector< std::vector<double> >& log_p2, 
		      std::vector<std::string>& sample_names, std::string& chrom_seq, 
		      StutterModel& stutter_model, vcf::VariantCallFile& beagle_imp_vcf){

  }
  */



  ~SeqStutterGenotyper(){
    delete region_;
    delete stutter_model_;
    delete [] log_p1_;
    delete [] log_p2_;
    delete [] sample_label_;
    delete [] log_aln_probs_;
    delete [] log_sample_posteriors_;
    delete [] log_allele_priors_;
    for (unsigned int i = 0; i < hap_blocks_.size(); i++)
      delete hap_blocks_[i];
    hap_blocks_.clear();
    delete haplotype_;
  }
  
  static void write_vcf_header(std::vector<std::string>& sample_names, std::ostream& out);
  
  void write_vcf_record(std::vector<std::string>& sample_names, std::ostream& out);
  

  bool genotype();  
  
};

#endif



