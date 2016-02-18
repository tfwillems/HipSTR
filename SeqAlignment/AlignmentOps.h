#ifndef ALIGNMENT_OPS_H_
#define ALIGNMENT_OPS_H_

#include <string>
#include <vector>

#include "../bamtools/include/api/BamAlignment.h"
#include "AlignmentData.h"


extern const int ALIGN_WINDOW_WIDTH;

bool GetIntBamTag(const BamTools::BamAlignment& aln, const std::string& tag_name, int* destination);

bool realign(BamTools::BamAlignment& alignment, std::string& ref_sequence, Alignment& new_alignment);

#endif
