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

void visualizeAlignments(std::vector< std::vector<Alignment> >& alns, std::vector<std::string>& sample_names,
			 std::map<std::string, std::string>& sample_info, std::vector<HapBlock*>& hap_blocks,
			 std::string& chrom_seq, std::string locus_id, bool draw_locus_id,
			 std::ostream& output);

#endif
