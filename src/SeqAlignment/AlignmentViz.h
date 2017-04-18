#ifndef ALIGNMENT_VIZ_H_
#define ALIGNMENT_VIZ_H_

#include <climits>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <vector>

#include "AlignmentData.h"
#include "HapBlock.h"

void visualizeAlignments(const std::vector< std::vector<Alignment> >& alns, const std::vector<std::string>& sample_names,
			 const std::map<std::string, std::string>& sample_info, const std::vector<HapBlock*>& hap_blocks,
			 const std::string& chrom_seq, const std::string& locus_id, bool draw_locus_id,
			 std::ostream& output);

#endif
