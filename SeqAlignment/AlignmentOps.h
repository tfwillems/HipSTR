#ifndef ALIGNMENT_OPS_H_
#define ALIGNMENT_OPS_H_

#include <string>
#include <vector>

#include "../bamtools/include/api/BamAlignment.h"
#include "AlignmentData.h"

extern const int ALIGN_WINDOW_WIDTH;

bool GetIntBamTag(const BamTools::BamAlignment& aln, const std::string& tag_name, int* destination);

bool realign(BamTools::BamAlignment& alignment, std::string& ref_sequence, Alignment& new_alignment);

void convertAlignment(BamTools::BamAlignment& alignment, std::string& ref_sequence, Alignment& new_alignment);

bool startsWithSoftClip(const BamTools::BamAlignment& aln);
bool endsWithSoftClip(const BamTools::BamAlignment& aln);
bool startsWithHardClip(const BamTools::BamAlignment& aln);
bool endsWithHardClip(const BamTools::BamAlignment& aln);
bool matchesReference(const BamTools::BamAlignment& aln);


void trimAlignment(BamTools::BamAlignment& aln, int32_t min_read_start, int32_t max_read_stop, char min_base_qual='~');

void trimLowQualityEnds(BamTools::BamAlignment& aln, char min_base_qual);

#endif
