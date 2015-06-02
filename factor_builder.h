#ifndef FACTOR_BUILDER_H_
#define FACTOR_BUILDER_H_

#include <map>
#include <vector>

#include <dai/factorgraph.h>
#include "bamtools/include/api/BamAlignment.h"
#include "vcflib/src/Variant.h"

#include "base_quality.h"
#include "snp_tree.h"
#include "stutter_model.h"

inline double log_sum_exp(double v1, double v2);

void add_GT_factors(vcflib::Variant& variant, std::map<std::string, int>& sample_indices,
		    int min_bp_length, int max_bp_length, bool ignore_phase,
		    std::vector<dai::Var>& node_1s, std::vector<dai::Var>& node_2s,
		    std::vector<dai::Factor>& factors, std::vector< std::pair<dai::Var,int> >& clamps);

void add_GL_factors(vcflib::Variant& variant, std::map<std::string, int>& sample_indices,
		    int min_allele, int max_allele, std::vector<dai::Var>& node_1s, std::vector<dai::Var>& node_2s,
		    std::vector<dai::Factor>& factors);

void add_read_factors(std::vector<int>& nrepeats, std::vector<double>& snp_log_p1s, std::vector<double>& snp_log_p2s,
		      int min_allele, int max_allele, StutterModel& stutter_model,
		      dai::Var& node_id1, dai::Var& node_id2,
		      std::vector<dai::Factor>& factors);

void add_read_factors(std::vector<int>& nrepeats, int min_allele, int max_allele, StutterModel& stutter_model,
		      dai::Var& node_id1, dai::Var& node_id2,
		      std::vector<dai::Factor>& factors);

typedef std::vector<BamTools::BamAlignment> AlnVector;

void add_read_factors(AlnVector& paired_str_reads, AlnVector& mate_reads, AlnVector& unpaired_str_reads,
		      BaseQuality& base_qualities, SNPTree* snp_tree, int min_allele, int max_allele, StutterModel& model,
		      dai::Var node_id1, dai::Var node_id2,
		      std::vector<dai::Factor>& factors);

void add_sample_factors(std::vector<std::string>& samples, int num_alleles,
                        std::vector<dai::Var>& node_1s, std::vector<dai::Var>& node_2s, std::map<std::string, int>& sample_indices);

#endif
