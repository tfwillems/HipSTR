#ifndef ALIGNMENT_OPS_H_
#define ALIGNMENT_OPS_H_

#include <string>
#include <vector>

#include "../bamtools/include/api/BamAlignment.h"
#include "AlignmentData.h"


extern const std::string START_TAG;
extern const std::string STOP_TAG;
extern const std::string RG_TAG;
extern const std::string SAMPLE_TAG;
extern const std::string MOTIF_TAG;
extern const int ALIGN_WINDOW_WIDTH;

bool compareAlignments(const BamTools::BamAlignment& alignment_1, const BamTools::BamAlignment& alignment_2);


void sortAlignments(std::vector<BamTools::BamAlignment>& alignments);


std::string getCigarString(BamTools::BamAlignment& alignment);


std::string getCigarString(std::vector<BamTools::CigarOp>& cigar_list);

std::string getSampleName(const BamTools::BamAlignment& alignment);

bool GetIntBamTag(const BamTools::BamAlignment& aln, const std::string& tag_name, int* destination);

bool GetStringBamTag(const BamTools::BamAlignment& alignment, const std::string& tag_name, std::string& value);

bool realign(BamTools::BamAlignment& alignment, std::string& ref_sequence, Alignment& new_alignment);

#endif
