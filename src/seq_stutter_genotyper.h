#ifndef SEQ_STUTTER_GENOTYPER_H_
#define SEQ_STUTTER_GENOTYPER_H_

#include <algorithm>
#include <iostream>
#include <map>
#include <math.h>
#include <set>
#include <string>
#include <vector>

#include "bam_io.h"
#include "base_quality.h"
#include "genotyper.h"
#include "read_pooler.h"
#include "region.h"
#include "stutter_model.h"
#include "vcf_input.h"
#include "vcf_reader.h"
#include "vcf_writer.h"

#include "SeqAlignment/AlignmentData.h"
#include "SeqAlignment/AlignmentMatrixCache.h"
#include "SeqAlignment/AlignmentTraceback.h"
#include "SeqAlignment/Haplotype.h"
#include "SeqAlignment/HapBlock.h"

class SeqStutterGenotyper : public Genotyper {
 private:
  int REF_FLANK_LEN;
  double STRAND_TOLERANCE;

  BaseQuality base_quality_;
  ReadPooler pooler_;
  int* pool_index_;                                 // Pool index for each read

  typedef std::vector<Alignment> AlnList;
  AlnList alns_;                                        // Vector of left-aligned alignments
  std::vector<AlignmentMatrixCache*> fw_matrix_caches_; // Matrix caches for each pool
  std::vector<AlignmentMatrixCache*> rv_matrix_caches_; // Matrix caches for each pool
  std::vector<HapBlock*> hap_blocks_;                   // Haplotype blocks
  Haplotype* haplotype_;                                // Potential STR haplotypes
  std::vector<std::string> call_sample_;                // True iff we should try to genotype the sample with the associated index
                                                        // Based on the deletion boundaries in the sample's reads

  bool initialized_; // True iff initialization succeeded and genotyping can proceed

  // 0-based seed index for each read
  // -1 denotes that no seed position was determined for the read
  int* seed_positions_;

  // VCF containing STR and SNP genotypes for a reference panel
  VCF::VCFReader* ref_vcf_;

  // If this flag is set, the genotyper will reassemble the flanking sequences after an initial round of genotyping
  bool reassemble_flanks_;

  // If true, clear matrix caches between each round of genotyping. Otherwise, allow them to persist
  // to reduce computation time but at the expense of increased memory usage
  bool clear_caches_;

  // Timing statistics (in seconds)
  double total_hap_build_time_;
  double total_hap_aln_time_;
  double total_aln_trace_time_;
  double total_assembly_time_;

  // Used to identify candidate haplotypes during flank reassembly
  int MIN_PATH_WEIGHT, MIN_KMER, MAX_KMER;

  // Cache of traced back alignments
  std::map<std::pair<int,int>, AlignmentTrace*> trace_cache_;

  // True iff both the indexed read and its mate overlap the STR and the current read's index is greater
  bool* second_mate_;

  // Set up the relevant data structures. Invoked by the constructor 
  bool build_haplotype(const std::string& chrom_seq, std::vector<StutterModel*>& stutter_models, std::ostream& logger);
  void init(std::vector<StutterModel *>& stutter_models, const std::string& chrom_seq, std::ostream& logger);

  void reorder_alleles(std::vector<std::string>& alleles,
		       std::vector<int>& old_to_new, std::vector<int>& new_to_old);

  // Extract the sequences for each allele and the VCF start position
  std::pair<int,int> get_alleles(const Region& region, int block_index, const std::string& chrom_seq,
				 int32_t& pos, std::vector<std::string>& alleles);

  void debug_sample(int sample_index, std::ostream& logger);
  
  // Modify the internal structures to remove the alleles at the associated indices
  // Designed to remove alleles who aren't the MAP genotype of any samples
  void remove_alleles(std::vector< std::vector<int> >& allele_indices);

  // Retrace the alignment for each read and store the associated pointers in the provided vector
  // Reads which were unaligned will have a NULL pointer
  void retrace_alignments(std::vector<AlignmentTrace*>& traced_alns);

  // Identify additional candidate STR alleles using the sequences observed in reads with stutter artifacts
  void get_stutter_candidate_alleles(int block_index, std::ostream& logger, std::vector<std::string>& candidate_seqs);

  // Aligns each read to each of the candidate haplotypes and stores the results in internal arrays
  void calc_hap_aln_probs(std::vector<bool>& realign_to_haplotype);
  void calc_hap_aln_probs(std::vector<bool>& realign_to_haplotype, std::vector<bool>& realign_pool, std::vector<bool>& copy_read);

  // Identify alleles present in stutter artifacts. Align each read to the new haplotypes
  // containing these alleles and incorporate these alignment probabilities into the relevant data structures
  bool id_and_align_to_stutter_alleles(int max_total_haplotypes, std::ostream& logger);

  // Exploratory function related to identifying indels in the flanking sequences
  void analyze_flank_indels(std::ostream& logger);

  // Exploratory function related to identifying SNPs in the flanking sequences
  void analyze_flank_snps(std::ostream& logger);

  // Use local assembly to identify variants in the flanking sequences, generate new haplotypes and realign reads
  bool assemble_flanks(int max_total_haplotypes, int max_flank_haplotypes, double min_flank_freq, std::ostream& logger);

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
  void add_and_remove_alleles(std::vector< std::vector<int> >& alleles_to_remove,
			      std::vector< std::vector<std::string> >& alleles_to_add,
			      std::vector<bool>& realign_pool, std::vector<bool>& copy_read);

  bool genotype_samples(bool first_round, int max_total_haplotypes, int max_flank_haplotypes, double min_flank_freq,
			std::ostream& logger);

  double compute_allele_bias(int hap_a_read_count, int hap_b_read_count);

  void write_vcf_record(const std::vector<std::string>& sample_names, int hap_block_index, const Region& region,
			const std::string& chrom_seq, bool output_viz, bool viz_left_alns,
			std::ostream& html_output, VCFWriter* vcf_writer, std::ostream& logger);

  RegionGroup* region_group_;


  // Private unimplemented copy constructor and assignment operator to prevent operations
  SeqStutterGenotyper(const SeqStutterGenotyper& other);
  SeqStutterGenotyper& operator=(const SeqStutterGenotyper& other);

 public:
 SeqStutterGenotyper(const RegionGroup& region_group, bool haploid, bool reassemble_flanks, int ref_flank_len,
		     std::vector<Alignment>& alignments, std::vector< std::vector<double> >& log_p1,
		     std::vector< std::vector<double> >& log_p2,
		     const std::vector<std::string>& sample_names, const std::string& chrom_seq,
		     std::vector<StutterModel*>& stutter_models,
		     VCF::VCFReader* ref_vcf, std::ostream& logger): Genotyper(haploid, sample_names, log_p1, log_p2){
    region_group_          = region_group.copy();
    alns_                  = alignments;
    seed_positions_        = NULL;
    pool_index_            = NULL;
    haplotype_             = NULL;
    second_mate_           = NULL;
    MIN_PATH_WEIGHT        = 2;
    MIN_KMER               = 10;
    MAX_KMER               = 15;
    STRAND_TOLERANCE       = 0.1;
    REF_FLANK_LEN          = ref_flank_len;
    initialized_           = false;
    reassemble_flanks_     = reassemble_flanks;
    clear_caches_          = false;
    total_hap_build_time_  = total_hap_aln_time_  = 0;
    total_aln_trace_time_  = total_assembly_time_ = 0;
    ref_vcf_               = ref_vcf;
    assert(num_reads_ == alns_.size());
    init(stutter_models, chrom_seq, logger);
  }

  ~SeqStutterGenotyper(){
    delete region_group_;
    delete [] seed_positions_;
    delete [] pool_index_;
    delete [] second_mate_;

    for (int pool_index = 0; pool_index < pooler_.num_pools(); ++pool_index){
      // Clear the alignment matrix cache as it stores allocated matrices that won't otherwise be freed
      fw_matrix_caches_[pool_index]->clear();
      rv_matrix_caches_[pool_index]->clear();
      delete fw_matrix_caches_[pool_index];
      delete rv_matrix_caches_[pool_index];
    }

    for (auto trace_iter = trace_cache_.begin(); trace_iter != trace_cache_.end(); trace_iter++)
      delete trace_iter->second;
    for (unsigned int i = 0; i < hap_blocks_.size(); i++)
      delete hap_blocks_[i];
    delete haplotype_;
  }
  
  void write_vcf_record(const std::vector<std::string>& sample_names, const std::string& chrom_seq,
			bool output_viz, bool viz_left_alns,
			std::ostream& html_output, VCFWriter* vcf_writer, std::ostream& logger);

  double hap_build_time() { return total_hap_build_time_;  }
  double hap_aln_time()   { return total_hap_aln_time_;    }
  double aln_trace_time() { return total_aln_trace_time_;  }
  double assembly_time()  { return total_assembly_time_;   }

  bool genotype(int max_total_haplotypes, int max_flank_haplotypes, double min_flank_freq, std::ostream& logger);

  /*
   * Recompute the stutter model(s) using the PCR artifacts obtained from the ML alignments
   * and regenotype the samples using this new model
  */
  bool recompute_stutter_models(std::ostream& logger, int max_total_haplotypes, int max_flank_haplotypes, double min_flank_freq,
				int max_em_iter, double abs_ll_converge, double frac_ll_converge);
};

#endif
