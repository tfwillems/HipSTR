#ifndef HAPLOTYPE_GENERATOR_H_
#define HAPLOTYPE_GENERATOR_H_

#include <string>
#include <vector>

#include "AlignmentData.h"
#include "../region.h"
#include "../stutter_model.h"
#include "Haplotype.h"
#include "HapBlock.h"

Haplotype* generate_haplotype(Region& str_region, int32_t max_ref_flank_len, std::string& chrom_seq,
                              std::vector< std::vector<Alignment> >& alignments, std::vector<std::string>& vcf_alleles,
			      StutterModel* stutter_model, bool search_bams_for_alleles,
			      std::vector<HapBlock*>& blocks);


#endif
