#ifndef SNP_PHASING_QUALITY_H_
#define SNP_PHASING_QUALITY_H_

#include <vector>

#include "bam_reader.h"
#include "base_quality.h"
#include "snp_tree.h"

void extract_bases_and_qualities(BamAlignment& aln, std::vector<SNP>& snps,
				 std::vector<char>& bases,  std::vector<char>& quals);

void add_log_phasing_probs(BamAlignment& aln, SNPTree* tree, BaseQuality& base_qualities,
			   double& log_p1, double& log_p2, int32_t& p1_match_count, int32_t& p2_match_count, int32_t& mismatch_count);

void calc_het_snp_factors(std::vector<BamAlignment>& str_reads, std::vector<BamAlignment>& mate_reads,
			  BaseQuality& base_qualities, SNPTree* snp_tree,
			  std::vector<double>& log_p1s, std::vector<double>& log_p2s, int32_t& match_count, int32_t& mismatch_count);

void calc_het_snp_factors(std::vector<BamAlignment>& str_reads, BaseQuality& base_qualities, SNPTree* snp_tree,
			  std::vector<double>& log_p1s, std::vector<double>& log_p2s, int32_t& match_count, int32_t& mismatch_count);

#endif
