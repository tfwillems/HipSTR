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

#include "base_quality.h"
#include "genotyper.h"
#include "read_pooler.h"
#include "region.h"
#include "stutter_model.h"
#include "vcf_input.h"
#include "vcf_reader.h"

#include "SeqAlignment/AlignmentData.h"
#include "SeqAlignment/AlignmentTraceback.h"
#include "SeqAlignment/Haplotype.h"
#include "SeqAlignment/HapBlock.h"

class SeqStutterGenotyper : public Genotyper {
 private:
  int MAX_REF_FLANK_LEN;
  BaseQuality base_quality_;
  ReadPooler pooler_;
  int* pool_index_;                               // Pool index for each read
  std::vector<int> bp_diffs_;                     // Base pair difference of each read from reference

  typedef std::vector<Alignment> AlnList;
  AlnList alns_;                                  // Vector of left-aligned alignments
  std::vector<bool> use_for_haps_;                // True iff we should use the alignment for identifying candidate haplotypes
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

  // VCF containing STR and SNP genotypes for a reference panel
  VCF::VCFReader* ref_vcf_;

  // If this flag is set, reads with identical sequences are pooled and their base emission error
  // probabilities averaged. Each unique sequence is then only aligned once using these
  // probabilities. Should result in significant speedup but may introduce genotyping errors
  bool pool_identical_seqs_;

  // True iff we only report genotypes for samples with >= 1 read
  // In an imputation-only setting, this should be set to false
  bool require_one_read_;

  // Timing statistics (in seconds)
  double total_hap_build_time_;
  double total_hap_aln_time_;
  double total_aln_trace_time_;
  double total_bootstrap_time_;

  // Cache of traced back alignments
  std::map<std::pair<int,int>, AlignmentTrace*> trace_cache_;

  // True iff both the indexed read and its mate overlap the STR and the current read's index is greater
  bool* second_mate_;

  /* Compute the alignment probabilites between each read and each haplotype */
  double calc_align_probs();

  // Set up the relevant data structures. Invoked by the constructor 
  void init(StutterModel& stutter_model, std::string& chrom_seq, std::ostream& logger);

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

  // Compute bootstrapped quality scores by resampling reads and determining how frequently
  // the genotypes match the ML genotype
  void compute_bootstrap_qualities(int num_iter, std::vector<double>& bootstrap_qualities);

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
  SeqStutterGenotyper(Region& region, bool haploid,
		      std::vector<Alignment>& alignments, std::vector<bool>& use_to_generate_haps, std::vector<int>& bp_diffs,
		      std::vector< std::vector<double> >& log_p1, std::vector< std::vector<double> >& log_p2,
		      std::vector<std::string>& sample_names, std::string& chrom_seq,
		      bool pool_identical_seqs,
		      StutterModel& stutter_model, VCF::VCFReader* ref_vcf, std::ostream& logger): Genotyper(region, haploid, false, sample_names, log_p1, log_p2){
    alns_                  = alignments;
    bp_diffs_              = bp_diffs;
    use_for_haps_          = use_to_generate_haps;
    seed_positions_        = NULL;
    pool_index_            = NULL;
    haplotype_             = NULL;
    second_mate_           = NULL;
    ref_vcf_               = ref_vcf;
    MAX_REF_FLANK_LEN      = 30;
    pos_                   = -1;
    pool_identical_seqs_   = pool_identical_seqs;
    total_hap_build_time_  = total_hap_aln_time_    = 0;
    total_aln_trace_time_  = total_bootstrap_time_  = 0;
    ref_vcf_               = ref_vcf;
    alleles_from_bams_     = true;

    require_one_read_      = true;
    /* TO DO: Properly set this flag based on whether the VCF has the required FORMAT fields
    // True iff no allele priors are available (for imputation)
    if (ref_vcf == NULL)
      require_one_read_ = true;
    else
      require_one_read_ = (ref_vcf->formatTypes.find(PGP_KEY) == ref_vcf->formatTypes.end());
    */
    assert(num_reads_ == alns_.size() && num_reads_ == bp_diffs_.size() && num_reads_ == use_for_haps_.size());
    init(stutter_model, chrom_seq, logger);
  }

  ~SeqStutterGenotyper(){
    delete [] seed_positions_;
    delete [] pool_index_;
    delete [] second_mate_;
    if (ref_vcf_ != NULL)
      delete ref_vcf_;
    for (auto trace_iter = trace_cache_.begin(); trace_iter != trace_cache_.end(); trace_iter++)
      delete trace_iter->second;
    for (unsigned int i = 0; i < hap_blocks_.size(); i++)
      delete hap_blocks_[i];
    delete haplotype_;
  }
  
  static void write_vcf_header(std::string& full_command, std::vector<std::string>& sample_names, bool output_gls, bool output_pls, bool output_phased_gls, std::ostream& out);

  /*
   *  Returns true iff the read with the associated retraced maximum log-likelihood alignment should be used in genotyping
   *  Considers factors such as indels in the regions flanking the STR block and the total number of matched bases
   */
  bool use_read(AlignmentTrace* trace);

  void write_vcf_record(std::vector<std::string>& sample_names, bool print_info, std::string& chrom_seq,
			bool output_bootstrap_qualities, bool output_gls, bool output_pls, bool output_phased_gls,
			bool output_allreads, bool output_pallreads, bool output_mallreads, bool output_viz, float max_flank_indel_frac,
			bool visualize_left_alns,
			std::ostream& html_output, std::ostream& out, std::ostream& logger);


  double hap_build_time() { return total_hap_build_time_;  }
  double hap_aln_time()   { return total_hap_aln_time_;    }
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
