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
  bool initialized_;       // True iff initialization succeeded and genotyping can proceed

  // 0-based seed index for each read
  // -1 denotes that no seed position was determined for the read
  int* seed_positions_;

  // VCF containing STR and SNP genotypes for a reference panel
  VCF::VCFReader* ref_vcf_;

  // If this flag is set, reads with identical sequences are pooled and their base emission error
  // probabilities averaged. Each unique sequence is then only aligned once using these
  // probabilities. Should result in significant speedup but may introduce genotyping errors
  bool pool_identical_seqs_;

  // Timing statistics (in seconds)
  double total_hap_build_time_;
  double total_hap_aln_time_;
  double total_aln_trace_time_;
  double total_bootstrap_time_;

  // Cache of traced back alignments
  std::map<std::pair<int,int>, AlignmentTrace*> trace_cache_;

  // True iff both the indexed read and its mate overlap the STR and the current read's index is greater
  bool* second_mate_;

  // Set up the relevant data structures. Invoked by the constructor 
  void init(StutterModel& stutter_model, std::string& chrom_seq, std::ostream& logger);

  void reorder_alleles(std::vector<std::string>& alleles,
		       std::vector<int>& old_to_new, std::vector<int>& new_to_old);

  // Extract the sequences for each allele and the VCF start position
  void get_alleles(Region& region, int block_index, std::string& chrom_seq,
		   int32_t& pos, std::vector<std::string>& alleles);

  void debug_sample(int sample_index, std::ostream& logger);
  
  // Modify the internal structures to remove the alleles at the associated indices
  // Designed to remove alleles who aren't the MAP genotype of any samples
  // However, it does not modify any of the haplotype-related data structures
  void remove_alleles(std::vector< std::vector<int> >& allele_indices);

  // Compute bootstrapped quality scores by resampling reads and determining how frequently
  // the genotypes match the ML genotype
  void compute_bootstrap_qualities(int num_iter, std::vector< std::vector<double> >& bootstrap_qualities);

  // Retrace the alignment for each read and store the associated pointers in the provided vector
  // Reads which were unaligned will have a NULL pointer
  void retrace_alignments(std::vector<AlignmentTrace*>& traced_alns);

  // Filter reads based on their retraced ML alignments
  void filter_alignments(std::ostream& logger, std::vector<int>& masked_reads);

  // Identify additional candidate STR alleles using the sequences observed
  // in reads with stutter artifacts
  void get_stutter_candidate_alleles(int block_index, std::ostream& logger, std::vector<std::string>& candidate_seqs);

  // Align each read to each of the candidate alleles, and store the results in the provided arrays
  void calc_hap_aln_probs(std::vector<bool>& realign_to_haplotype);

  // Identify alleles present in stutter artifacts
  // Align each read to these alleles and incorporate these alignment probabilities and
  // alleles into the relevant data structures
  bool id_and_align_to_stutter_alleles(std::string& chrom_seq, std::ostream& logger);

  // Exploratory function related to identifying indels in the flanking sequences
  void analyze_flank_indels(std::ostream& logger);

  // Determines the allele index in the given haplotype block that is associated with each haplotype configuration
  // Stores the results in the provided vector
  void haps_to_alleles(int hap_block_index, std::vector<int>& allele_indices);


  // If CHECK_CALLED,  identifies any alleles that are not in any samples' ML genotypes
  // If CHECK_SPANNED, identifies any alleles that are not spanned by any reads without stutter artifacts
  // Adds any identifed alleles to the provided vector, separated by haplotype block index
  void get_unused_alleles(bool check_spanned, bool check_called,
			  std::vector< std::vector<int> >& allele_indices, int& num_aff_blocks, int& num_aff_alleles);

  // Add and/or remove the provided alleles to the underlying haplotype structures. Realigns each read to any novel alleles
  // and updates the alignment probabilities and genotype posteriors accordingly.
  void add_and_remove_alleles(std::vector< std::vector<int> >& alleles_to_remove,
			      std::vector< std::vector<std::string> >& alleles_to_add);

  void write_vcf_record(std::vector<std::string>& sample_names, int hap_block_index, Region& region, std::string& chrom_seq, std::vector<double>& bootstrap_qualities,
			bool output_bootstrap_qualities, bool output_gls, bool output_pls, bool output_phased_gls,
			bool output_allreads, bool output_pallreads, bool output_mallreads, bool output_viz, float max_flank_indel_frac, bool viz_left_alns,
			std::ostream& html_output, std::ostream& out, std::ostream& logger);

  Region* region_;

 public:
  SeqStutterGenotyper(Region& region, bool haploid,
		      std::vector<Alignment>& alignments, std::vector<bool>& use_to_generate_haps, std::vector<int>& bp_diffs,
		      std::vector< std::vector<double> >& log_p1, std::vector< std::vector<double> >& log_p2,
		      std::vector<std::string>& sample_names, std::string& chrom_seq,
		      bool pool_identical_seqs,
		      StutterModel& stutter_model, VCF::VCFReader* ref_vcf, std::ostream& logger): Genotyper(haploid, false, sample_names, log_p1, log_p2){
    region_                = region.copy();
    alns_                  = alignments;
    bp_diffs_              = bp_diffs;
    use_for_haps_          = use_to_generate_haps;
    seed_positions_        = NULL;
    pool_index_            = NULL;
    haplotype_             = NULL;
    second_mate_           = NULL;
    ref_vcf_               = ref_vcf;
    MAX_REF_FLANK_LEN      = 30;
    initialized_           = false;
    pool_identical_seqs_   = pool_identical_seqs;
    total_hap_build_time_  = total_hap_aln_time_    = 0;
    total_aln_trace_time_  = total_bootstrap_time_  = 0;
    ref_vcf_               = ref_vcf;
    alleles_from_bams_     = true;
    assert(num_reads_ == alns_.size() && num_reads_ == bp_diffs_.size() && num_reads_ == use_for_haps_.size());
    init(stutter_model, chrom_seq, logger);
  }

  ~SeqStutterGenotyper(){
    delete region_;
    delete [] seed_positions_;
    delete [] pool_index_;
    delete [] second_mate_;
    for (auto trace_iter = trace_cache_.begin(); trace_iter != trace_cache_.end(); trace_iter++)
      delete trace_iter->second;
    for (unsigned int i = 0; i < hap_blocks_.size(); i++)
      delete hap_blocks_[i];
    delete haplotype_;
  }
  
  /*
   *  Returns true iff the read with the associated retraced maximum log-likelihood alignment should be used in genotyping
   *  Considers factors such as indels in the regions flanking the STR block and the total number of matched bases
   */
  bool use_read(AlignmentTrace* trace);

  void write_vcf_record(std::vector<std::string>& sample_names, std::string& chrom_seq,
			bool output_bootstrap_qualities, bool output_gls, bool output_pls, bool output_phased_gls,
			bool output_allreads, bool output_pallreads, bool output_mallreads, bool output_viz, float max_flank_indel_frac, bool viz_left_alns,
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
