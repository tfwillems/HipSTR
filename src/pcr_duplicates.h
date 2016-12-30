#ifndef PCR_DUPLICATES_H_
#define PCR_DUPLICATES_H_

#include <iostream>
#include <map>
#include <vector>

#include "bamtools/include/api/BamAlignment.h"
#include "base_quality.h"

void remove_pcr_duplicates(BaseQuality& base_quality, bool use_bam_rgs,
			   std::map<std::string, std::string>& rg_to_library,
			   std::vector< std::vector<BamTools::BamAlignment> >& paired_strs_by_rg,
			   std::vector< std::vector<BamTools::BamAlignment> >& mate_pairs_by_rg,
			   std::vector< std::vector<BamTools::BamAlignment> >& unpaired_strs_by_rg, std::ostream& logger);

#endif
