#ifndef EXTRACT_INDELS_H_
#define EXTRACT_INDELS_H_

#include <vector>

#include "bam_reader.h"
#include "SeqAlignment/AlignmentData.h"

bool ExtractCigar(const std::vector<CigarElement>& cigar_data,
		  const int& cigar_start, const int& region_start, const int& region_end,
		  int& bp_diff_from_ref);

bool ExtractCigar(const std::vector<CigarOp>& cigar_data,
		  const int& cigar_start, const int& region_start, const int& region_end,
		  int& bp_diff_from_ref);

#endif
