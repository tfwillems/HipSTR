#ifndef ALIGNMENT_TRACEBACK_H_
#define ALIGNMENT_TRACEBACK_H_

#include <string>

#include "AlignmentData.h"

std::string stitch(const std::string& hap_aln, const std::string& read_aln, int h_index, int r_index, int increment);

void stitch_alignment_trace(int32_t hap_start, const std::string& hap_aln_to_ref, 
			    const std::string& read_aln_to_hap, int hap_index, int seed_base, Alignment& orig_aln,
			    Alignment& new_aln);

#endif
