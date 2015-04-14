#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "RepeatRegion.h"
#include "../bamtools/include/api/BamAlignment.h"
#include "../vcflib/src/Variant.h"
#include "../base_quality.h"


void arrangeAlignments(RepeatRegion& rep_region, 
		       std::string& ref_sequence,
		       std::vector<BamTools::BamAlignment>& alignments, 
		       std::string locus_id, bool draw_locus_id,
		       std::ostream& output, 
		       std::map<std::string, std::string>& sample_info,
		       BaseQuality& base_quality,
		       vcf::VariantCallFile* vcf_data);
