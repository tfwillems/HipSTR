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
  // Locus information
  Region* region_;

  unsigned int num_reads_; // Total number of reads across all samples
  int num_samples_;        // Total number of samples
  int motif_len_;          // # bp in STR motif
  int num_alleles_;        // Number of valid alleles
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
  std::vector<HapBlock*> hap_blocks_;          // Haplotype blocks
  Haplotype* haplotype_;                       // Potential STR haplotypes
  std::vector<bool> call_sample_;              // True iff we should try to genotype the sample with the associated index
                                               // Based on the deletion boundaries in the sample's reads
  std::vector<bool> got_priors_;               // True iff we obtained the priors for the sample with the associated index
                                               // from the VCF. If false, the sample will not be genotyped. This data structure
                                               // isn't used if priors aren't read from a VCF

  bool alleles_from_bams_; // Flag that determines if we examine BAMs for candidate alleles

  std::vector<std::string> alleles_; // Vector of indexed alleles
  int32_t pos_;                      // Position of reported alleles in VCF     

  // 0-based seed index for each read
  int* seed_positions_;

  // Iterates through reads and then alleles by their indices
  double* log_aln_probs_;

  // Iterates through allele_1, allele_2 and then samples by their indices
  double* log_sample_posteriors_; 
  
  // Iterates through allele_1, allele_2 and then samples by their indices
  // Only used if per-allele priors have been specified for each sample
  double* log_allele_priors_;

  // VCF containing STR and SNP genotypes for a reference panel
  vcflib::VariantCallFile* ref_vcf_;

  // If this flag is set, reads with identical sequences are pooled and their base emission error
  // probabilities averaged. Each unique sequence is then only aligned once using these
  // probabilities. Should result in significant speedup but may introduce genotyping errors
  bool pool_identical_seqs_;

  /* Combine reads with identical base sequences into a single representative alignment */
  void combine_reads(std::vector<Alignment>& alignments, Alignment& pooled_aln);

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
  
  // Attempt to identify additional haplotypes given all alignments and the current
  // haplotype structure. Modifies the underlying haplotype and haplotype blocks accordingly
  void expand_haplotype();


  // Identify a list of alleles that aren't the MAP genotype for any sample
  // Doesn't include the reference allele (index = 0)
  void get_uncalled_alleles(std::vector<int>& allele_indices);

  // Modify the internal structures to remove the alleles at the associated indices
  // Designed to remove alleles who aren't the MAP genotype of any samples
  // However, it does not modify any of the haplotype-related data structures
  void remove_alleles(std::vector<int>& allele_indices);

  std::set<std::string> expanded_alleles_;

  // True iff we only report genotypes for samples with >= 1 read
  // In an imputation-only setting, this should be set to false
  bool require_one_read_;

 public:
  SeqStutterGenotyper(Region& region,
		      std::vector< std::vector<BamTools::BamAlignment> >& alignments,
		      std::vector< std::vector<double> >& log_p1, 
		      std::vector< std::vector<double> >& log_p2, 
		      std::vector<std::string>& sample_names, std::string& chrom_seq, 
		      StutterModel& stutter_model, vcflib::VariantCallFile* ref_vcf){
    assert(alignments.size() == log_p1.size() && alignments.size() == log_p2.size() && alignments.size() == sample_names.size());
    log_p1_                = NULL;
    log_p2_                = NULL;
    seed_positions_        = NULL;
    log_aln_probs_         = NULL;
    log_sample_posteriors_ = NULL;
    log_allele_priors_     = NULL;
    sample_label_          = NULL;
    haplotype_             = NULL;
    MAX_REF_FLANK_LEN      = 30;
    pos_                   = -1;
    pool_identical_seqs_   = false;

    // True iff no allele priors are available (for imputation)
    require_one_read_ = (ref_vcf == NULL);
    
    region_       = region.copy();
    num_samples_  = alignments.size();
    sample_names_ = sample_names;
    for (unsigned int i = 0; i < sample_names.size(); i++)
      sample_indices_.insert(std::pair<std::string,int>(sample_names[i], i));
    stutter_model_      = stutter_model.copy();
    ref_vcf_            = ref_vcf;
    alleles_from_bams_  = true;
    init(alignments, log_p1, log_p2, sample_names, chrom_seq);
  }

  ~SeqStutterGenotyper(){
    delete region_;
    delete stutter_model_;
    delete [] log_p1_;
    delete [] log_p2_;
    delete [] sample_label_;
    delete [] seed_positions_;
    delete [] log_aln_probs_;
    delete [] log_sample_posteriors_;
    delete [] log_allele_priors_;
    for (unsigned int i = 0; i < hap_blocks_.size(); i++)
      delete hap_blocks_[i];
    hap_blocks_.clear();
    delete haplotype_;
  }
  
  static void write_vcf_header(std::vector<std::string>& sample_names, bool output_gls, bool output_pls, std::ostream& out);

  /*
   *  When aligning to each haplotype, align each unique sequence instead of each read.
   *  As quality scores, the genotype utilizes the average of the base quality scores (raw probabilities) for
   *  reads with identical sequences. Should result in significant speedup if many reads have the same sequence.
   *  By default, each read is aligned using its own quality scores.
   */
  void pool_identical_sequences(){
    pool_identical_seqs_ = true;
  }

  void write_vcf_record(std::vector<std::string>& sample_names, bool print_info, std::string& chrom_seq, bool output_gls, bool output_pls,
			bool output_viz, std::ostream& html_output, std::ostream& out);
  
  bool genotype();
};

#endif



