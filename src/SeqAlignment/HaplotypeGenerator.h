#ifndef HAPLOTYPE_GENERATOR_H_
#define HAPLOTYPE_GENERATOR_H_

#include <iostream>
#include <string>
#include <vector>

#include "AlignmentData.h"
#include "../error.h"
#include "../region.h"
#include "../stutter_model.h"
#include "Haplotype.h"
#include "HapBlock.h"

class HaplotypeGenerator {
 private:
  // Criteria used to determine whether a candidate sequence should be identified as an allele
  double MIN_FRAC_READS;
  double MIN_FRAC_SAMPLES;
  double MIN_FRAC_STRONG_SAMPLE;
  double MIN_READS_STRONG_SAMPLE;
  double MIN_STRONG_SAMPLES;

  // When extracting alleles in regions, we pad by these amounts to improve the capture of proximal indels
  int32_t LEFT_PAD, RIGHT_PAD;

  int32_t MIN_BLOCK_SPACING; // Minimum distance (bp) between variant haplotype blocks
  int32_t REF_FLANK_LEN;     // Maximum length of reference sequences flanking the variant haplotype blocks

  bool finished_; // True iff the underlying haplotype blocks are ready for downstream use
  std::string failure_msg_;
  int32_t min_aln_start_, max_aln_stop_;
  std::vector<HapBlock*> hap_blocks_;

  void trim(int ideal_min_length,
	    int32_t& region_start, int32_t& region_end, std::vector<std::string>& sequences) const;

  bool extract_sequence(const Alignment& aln, int32_t start, int32_t end, std::string& seq) const;

  void gen_candidate_seqs(const std::string& ref_seq, int ideal_min_length,
			  const std::vector< std::vector<Alignment> >& alignments, const std::vector<std::string>& vcf_alleles,
			  int32_t& region_start, int32_t& region_end, std::vector<std::string>& sequences) const;

  void get_aln_bounds(const std::vector< std::vector<Alignment> >& alignments,
		      int32_t& min_aln_start, int32_t& max_aln_stop) const;

 public:
  HaplotypeGenerator(int32_t min_aln_start, int32_t max_aln_stop){
    finished_                = false;
    MIN_FRAC_READS           = 0.05;
    MIN_FRAC_SAMPLES         = 0.05;
    MIN_FRAC_STRONG_SAMPLE   = 0.2;
    MIN_READS_STRONG_SAMPLE  = 2;
    MIN_STRONG_SAMPLES       = 1;
    LEFT_PAD                 = 5;
    RIGHT_PAD                = 5;
    MIN_BLOCK_SPACING        = 10;
    REF_FLANK_LEN            = 35;
    min_aln_start_           = min_aln_start;
    max_aln_stop_            = max_aln_stop;
  }

  bool add_vcf_haplotype_block(int32_t pos, const std::string& chrom_seq,
			       const std::vector<std::string>& vcf_alleles, const StutterModel* stutter_model);

  bool add_haplotype_block(const Region& region, const std::string& chrom_seq, const std::vector< std::vector<Alignment> >& alignments,
			   const std::vector<std::string>& vcf_alleles, const StutterModel* stutter_model);

  bool fuse_haplotype_blocks(const std::string& chrom_seq);

  const std::string& failure_msg(){ return failure_msg_; }

  const std::vector<HapBlock*> get_haplotype_blocks() const {
    if (!finished_)
      printErrorAndDie("Haplotype blocks are not ready for downstream use");
    return hap_blocks_;
  }
};

#endif
