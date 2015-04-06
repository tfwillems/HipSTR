#ifndef EXTRACT_INDELS_H_
#define EXTRACT_INDELS_H_

#include <vector>

#include "bamtools/include/api/BamAlignment.h"

bool ExtractCigar(std::vector<BamTools::CigarOp>& cigar_data,
		  const int& cigar_start, const int& region_start, const int& region_end,
		  int& bp_diff_from_ref);

bool ExtractCigarSingleFlank(std::vector<BamTools::CigarOp>& cigar_data,
			     const int& cigar_start, const int& region_start, const int& region_end,
			     int& min_bp_length);
#endif
