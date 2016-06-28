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
#include "read_pooler.h"
#include "region.h"
#include "stutter_model.h"
#include "vcf_input.h"

#include "SeqAlignment/AlignmentData.h"
#include "SeqAlignment/AlignmentTraceback.h"
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
  int* pool_index_;                // Pool index for each read
  BaseQuality base_quality_;
  ReadPooler pooler_;

  std::vector<int> bp_diffs_;                     // Base pair difference of each read from reference
  std::vector< std::vector<Alignment> > alns_;    // Vector of left-aligned alignments
  std::vector< std::vector<bool> > use_for_haps_; // True iff we should use the alignment for identifying candidate haplotypes
  std::vector<std::string> sample_names_;         // List of sample names
  std::map<std::string, int> sample_indices_;     // Mapping from sample name to index
  std::vector<HapBlock*> hap_blocks_;             // Haplotype blocks
  Haplotype* haplotype_;                          // Potential STR haplotypes
  std::vector<bool> call_sample_;                 // True iff we should try to genotype the sample with the associated index
                                                  // Based on the deletion boundaries in the sample's reads

  bool alleles_from_bams_; // Flag that determines if we examine BAMs for candidate alleles

  std::vector<std::string> alleles_; // Vector of indexed alleles
  int32_t pos_;                      // Position of reported alleles in VCF     

  // 0-based seed index for each read
  // -1 denotes that no seed position was determined for the read
  int* seed_positions_;

  // Iterates through reads and then alleles by their indices
  double* log_aln_probs_;

  // Iterates through allele_1, allele_2 and then samples by their indices
  double* log_sample_posteriors_; 

  // Total log-likelihoods for each sample
  double* sample_total_LLs_;
  
  // Iterates through allele_1, allele_2 and then samples by their indices
  // Only used if per-allele priors have been specified for each sample
  double* log_allele_priors_;

  // VCF containing STR and SNP genotypes for a reference panel
  vcflib::VariantCallFile* ref_vcf_;

  // If this flag is set, reads with identical sequences are pooled and their base emission error
  // probabilities averaged. Each unique sequence is then only aligned once using these
  // probabilities. Should result in significant speedup but may introduce genotyping errors
  bool pool_identical_seqs_;

  // True iff we only report genotypes for samples with >= 1 read
  // In an imputation-only setting, this should be set to false
  bool require_one_read_;

  // True iff the underlying marker is haploid
  bool haploid_;

  // Timing statistics (in seconds)
  double total_hap_build_time_;
  double total_left_aln_time_;
  double total_hap_aln_time_;
  double total_posterior_time_;
  double total_aln_trace_time_;
  double total_bootstrap_time_;

  // Cache of traced back alignments
  std::map<std::pair<int,int>, AlignmentTrace*> trace_cache_;

  // True iff both the indexed read and its mate overlap the STR and the current read's index is greater
  bool* second_mate_;

  // Read weights used to calculate posteriors (See calc_log_sample_posteriors function)
  // Used to account for special cases in which both reads in a pair overlap the STR by setting
  // the weight for the second read to zero. Elsewhere, the alignments probabilities for the two reads are summed
  std::vector<int> read_weights_;

  /* Compute the alignment probabilites between each read and each haplotype */
  double calc_align_probs();

  /* Initialize the priors required for computed sample genotype posteriors */
  void init_log_sample_priors(double* log_sample_ptr);

  /* Compute the posteriors for each sample using the haplotype probabilites, stutter model and read weights */
  double calc_log_sample_posteriors();  
  double calc_log_sample_posteriors(std::vector<int>& read_weights);

  // Set up the relevant data structures. Invoked by the constructor 
  void init(std::vector< std::vector<BamTools::BamAlignment> >& alignments,
	    std::vector< std::vector<double> >& log_p1,
	    std::vector< std::vector<double> >& log_p2,
	    std::vector<std::string>& sample_names, StutterModel& stutter_model, std::string& chrom_seq, std::ostream& logger);

  // Extract the sequences for each allele and the VCF start position
  void get_alleles(std::string& chrom_seq, std::vector<std::string>& alleles);

  void debug_sample(int sample_index, std::ostream& logger);
  
  // Identify a list of alleles that aren't the MAP genotype for any sample
  // Doesn't include the reference allele (index = 0)
  void get_uncalled_alleles(std::vector<int>& allele_indices);

  // Modify the internal structures to remove the alleles at the associated indices
  // Designed to remove alleles who aren't the MAP genotype of any samples
  // However, it does not modify any of the haplotype-related data structures
  void remove_alleles(std::vector<int>& allele_indices);

  // Determine the genotype associated with each sample based on the current genotype posteriors
  void get_optimal_genotypes(double* log_posterior_ptr, std::vector< std::pair<int, int> >& gts);

  // Compute bootstrapped quality scores by resampling reads and determining how frequently
  // the genotypes match the ML genotype
  void compute_bootstrap_qualities(int num_iter, std::vector<double>& bootstrap_qualities);

  // Convert a list of integers into a string with key|count pairs separated by semicolons
  // eg. -1,0,-1,2,2,1 will be converted into -1|2;0|1;1|1;2|2
  std::string condense_read_counts(std::vector<int>& read_diffs);

  // Retrace the alignment for each read and store the associated pointers in the provided vector
  // Reads which were unaligned will have a NULL pointer
  void retrace_alignments(std::ostream& logger, std::vector<AlignmentTrace*>& traced_alns);

  // Filter reads based on their retraced ML alignments
  void filter_alignments(std::ostream& logger, std::vector<int>& masked_reads);

  // Identify additional candidate STR alleles using the sequences observed
  // in reads with stutter artifacts
  void get_stutter_candidate_alleles(std::ostream& logger, std::vector<std::string>& candidate_seqs);

  // Align each read to each of the candidate alleles, and store the results in the provided arrays
  void calc_hap_aln_probs(Haplotype* haplotype, double* log_aln_probs, int* seed_positions);

  // Identify alleles present in stutter artifacts
  // Align each read to these alleles and incorporate these alignment probabilities and
  // alleles into the relevant data structures
  bool id_and_align_to_stutter_alleles(std::string& chrom_seq, std::ostream& logger);

  // Exploratory function related to identifying indels in the flanking sequences
  void analyze_flank_indels(std::ostream& logger);

  // Exploratory function related to identifying alleles without any spanning reads
  // These alleles should likely be removed
  void get_unspanned_alleles(std::vector<int>& allele_indices, std::ostream& logger);


 public:
  
  // In the VCF format fields for ALLREADS and MALLREADS, condense the fields into size|count
  // instead of a long comma-separated list of sizes e.g. -2,-2,0,-2,0 will be converted to -2|3;0|2
  static bool condense_read_count_fields;

  SeqStutterGenotyper(Region& region, bool haploid,
		      std::vector< std::vector<BamTools::BamAlignment> >& alignments,
		      std::vector< std::vector<double> >& log_p1, 
		      std::vector< std::vector<double> >& log_p2, 
		      std::vector<std::string>& sample_names, std::string& chrom_seq,
		      bool pool_identical_seqs,
		      StutterModel& stutter_model, vcflib::VariantCallFile* ref_vcf, std::ostream& logger){
    assert(alignments.size() == log_p1.size() && alignments.size() == log_p2.size() && alignments.size() == sample_names.size());
    log_p1_                = NULL;
    log_p2_                = NULL;
    seed_positions_        = NULL;
    log_aln_probs_         = NULL;
    log_sample_posteriors_ = NULL;
    sample_total_LLs_      = NULL;
    log_allele_priors_     = NULL;
    sample_label_          = NULL;
    pool_index_            = NULL;
    haplotype_             = NULL;
    second_mate_           = NULL;
    MAX_REF_FLANK_LEN      = 30;
    pos_                   = -1;
    pool_identical_seqs_   = pool_identical_seqs;
    haploid_               = haploid;
    total_hap_build_time_  = 0;
    total_left_aln_time_   = 0;
    total_hap_aln_time_    = 0;
    total_aln_trace_time_  = 0;
    total_posterior_time_  = 0;
    total_bootstrap_time_  = 0;

    // True iff no allele priors are available (for imputation)
    if (ref_vcf == NULL)
      require_one_read_ = true;
    else
      require_one_read_ = (ref_vcf->formatTypes.find(PGP_KEY) == ref_vcf->formatTypes.end());
    
    region_       = region.copy();
    num_samples_  = alignments.size();
    sample_names_ = sample_names;
    for (unsigned int i = 0; i < sample_names.size(); i++)
      sample_indices_.insert(std::pair<std::string,int>(sample_names[i], i));
    ref_vcf_            = ref_vcf;
    alleles_from_bams_  = true;
    init(alignments, log_p1, log_p2, sample_names, stutter_model, chrom_seq, logger);
  }

  ~SeqStutterGenotyper(){
    delete region_;
    delete [] log_p1_;
    delete [] log_p2_;
    delete [] sample_label_;
    delete [] seed_positions_;
    delete [] log_aln_probs_;
    delete [] log_sample_posteriors_;
    delete [] sample_total_LLs_;
    delete [] log_allele_priors_;
    delete [] pool_index_;
    delete [] second_mate_;
    for (auto trace_iter = trace_cache_.begin(); trace_iter != trace_cache_.end(); trace_iter++)
      delete trace_iter->second;
    trace_cache_.clear();
    for (unsigned int i = 0; i < hap_blocks_.size(); i++)
      delete hap_blocks_[i];
    hap_blocks_.clear();
    delete haplotype_;
  }
  
  static void write_vcf_header(std::string& full_command, std::vector<std::string>& sample_names, bool output_gls, bool output_pls, std::ostream& out);

  /*
   *  Returns true iff the read with the associated retraced maximum log-likelihood alignment should be used in genotyping
   *  Considers factors such as indels in the regions flanking the STR block and the total number of matched bases
   */
  bool use_read(AlignmentTrace* trace);

  void write_vcf_record(std::vector<std::string>& sample_names, bool print_info, std::string& chrom_seq,
			bool output_bootstrap_qualities, bool output_gls, bool output_pls,
			bool output_allreads, bool output_pallreads, bool output_mallreads, bool output_viz, float max_flank_indel_frac,
			bool visualize_left_alns,
			std::ostream& html_output, std::ostream& out, std::ostream& logger);


  double hap_build_time() { return total_hap_build_time_;  }
  double left_aln_time()  { return total_left_aln_time_;   }
  double hap_aln_time()   { return total_hap_aln_time_;    }
  double posterior_time() { return total_posterior_time_;  }
  double aln_trace_time() { return total_aln_trace_time_;  }
  double bootstrap_time() { return total_bootstrap_time_;  }

  bool genotype(std::string& chrom_seq, std::ostream& logger);

  /*
   * Recompute the stutter model(s) using the PCR artifacts obtained from the ML alignments
   * and regenotype the samples using this new model
  */
  bool recompute_stutter_models(std::string& chrom_seq, std::ostream& logger, int max_em_iter, double abs_ll_converge, double frac_ll_converge);
};

#endif
