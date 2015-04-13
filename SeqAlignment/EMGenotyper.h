#ifndef EMGENOTYPER_H
#define EMGENOTYPER_H

#include <climits>
#include <map>
#include <string>
#include <sstream>
#include <vector>

#include "bamtools/include/api/BamAlignment.h"
#include "bamtools/include/api/BamAux.h"

#include "BaseQuality.h"


void getMaxInsertionSizes(std::vector<BamTools::BamAlignment>& alignments, 
			  std::map<int32_t,int>& max_insertions);

void overlayAlignments(std::vector<BamTools::BamAlignment>& alignments,
		       std::map<int32_t,int>& max_insertions, 
		       std::vector<std::string>& results, 
		       int32_t& min_start, 
		       int32_t& max_stop);
						
void arrangeAlignments(std::vector<BamTools::BamAlignment>& alignments, 
		       std::string& ref_sequence, 
		       std::string locus_id, 
		       std::ostream& output,
		       bool draw_locus_id,
		       vcf::VariantCallFile* vcf_data,
		       std::map<std::string, std::string>& sample_info,
		       BaseQuality& base_quality);

#endif
