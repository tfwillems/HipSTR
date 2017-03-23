#ifndef ALIGNMENT_OPS_H_
#define ALIGNMENT_OPS_H_

#include <string>
#include <vector>

#include "../bam_io.h"
#include "AlignmentData.h"

extern const int ALIGN_WINDOW_WIDTH;

bool realign(BamAlignment& alignment, std::string& ref_sequence, Alignment& new_alignment);

void convertAlignment(BamAlignment& alignment, std::string& ref_sequence, Alignment& new_alignment);

#endif
