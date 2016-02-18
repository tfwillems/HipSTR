#include "AlignmentOps.h"
#include "../error.h"
#include "NWNoRefEndPenalty.h"
#include "../stringops.h"

extern const int ALIGN_WINDOW_WIDTH = 75;

bool GetIntBamTag(const BamTools::BamAlignment& aln, const std::string& tag_name, int* destination) {
  char tag_type;
  if (!aln.GetTagType(tag_name, tag_type)) {return false;}
  switch (tag_type) {
  case (BamTools::Constants::BAM_TAG_TYPE_INT32):
    return aln.GetTag(tag_name, *destination);
  case (BamTools::Constants::BAM_TAG_TYPE_INT8):
    int8_t d8;
    if (!aln.GetTag(tag_name, d8))
      return false;
    *destination = static_cast<int>(d8);
    return true;
  case (BamTools::Constants::BAM_TAG_TYPE_UINT8):
    uint8_t ud8;
    if (!aln.GetTag(tag_name, ud8))
      return false;
    *destination = static_cast<int>(ud8);
    return true;
  case (BamTools::Constants::BAM_TAG_TYPE_INT16):
    int16_t d16;
    if (!aln.GetTag(tag_name, d16))
      return false;
    *destination = static_cast<int>(d16);
    return true;
  case (BamTools::Constants::BAM_TAG_TYPE_UINT16):
    uint16_t ud16;
    if (!aln.GetTag(tag_name, ud16))
      return false;
    *destination = static_cast<int>(ud16);
    return true;
  case (BamTools::Constants::BAM_TAG_TYPE_UINT32):
    uint32_t ud32;
    if (!aln.GetTag(tag_name, ud32))
      return false;
    *destination = static_cast<int>(ud32);
    return true;
  default:
    printErrorAndDie("Invalid BAM integer flag type");
    break;
  }
  return false;
}

/*
 Realign read to reference region using left alignment variant. Store the new alignment information using
 the provided Alignment reference. Converts the bases to their upper case variants
 */
bool realign(BamTools::BamAlignment& alignment, std::string& ref_sequence, Alignment& new_alignment){
    int32_t start        = std::max(alignment.Position-ALIGN_WINDOW_WIDTH-1, 0);
    int32_t stop         = std::min(alignment.GetEndPosition()+ALIGN_WINDOW_WIDTH-1, (int32_t)(ref_sequence.size()-1));
    int32_t length       = stop-start+1;
    std::string ref_seq  = ref_sequence.substr(start, length);
    std::string read_seq = alignment.QueryBases;
    
    // Realign read using left-alignment method
    std::string ref_al, read_al;
    float score;
    std::vector<BamTools::CigarOp> cigar_list;
    bool left_aligned = NWNoRefEndPenalty::LeftAlign(ref_seq, read_seq, ref_al, read_al, &score, cigar_list);
    
    // Calculate number of leading spaces in read's alignment and start position
    unsigned int num_lead = 0;
    while (num_lead < read_al.size() && read_al[num_lead] == '-')
        num_lead++;
    int32_t start_position = start + num_lead;
    
    // Calculate number of trailing spaces in read's alignment
    int trail_index = read_al.size()-1;
    while (trail_index >= 0 && read_al[trail_index] == '-')
        trail_index--;
    int num_trail = read_al.size()-1-trail_index;
    
    // Calculate alignment end position
    int32_t end_position = start_position;
    bool halt = false;
    for (std::vector<BamTools::CigarOp>::iterator cigar_iter = cigar_list.begin(); cigar_iter != cigar_list.end() && !halt; cigar_iter++){
      switch(cigar_iter->Type){
      case 'X': case '=': case 'D':
	end_position += cigar_iter->Length;
	break;
      case 'I':
	break;
      case 'S':
	halt = true;
	break;
      default:
	std::stringstream msg;
	msg << "Invalid CIGAR character " << cigar_iter->Type << " in adjustAlignment()";
	printErrorAndDie(msg.str());
	break;
      }
    }
    end_position--;
    
    // Determine if any leading and trailing soft-clips should be applied in new alignment
    int num_head_sclips = 0, num_back_sclips = ref_al.size()-1;
    while (num_head_sclips < ref_al.size() && ref_al[num_head_sclips] == '-')
        num_head_sclips++;
    while (num_back_sclips > 0 && ref_al[num_back_sclips] == '-')
        num_back_sclips--;
    num_back_sclips = ref_al.size()-1-num_back_sclips;

    // Store new alignment information
    if (alignment.QueryBases.size() != alignment.Qualities.size())
        printErrorAndDie("Lengths of sequence and quality strings don't match");
    std::string base_qualities = alignment.Qualities.substr(num_head_sclips, read_seq.size()-num_head_sclips-num_back_sclips);
    std::string sequence       = uppercase(read_seq.substr(num_head_sclips, read_seq.size()-num_head_sclips-num_back_sclips));
    std::string alignment_seq  = uppercase(read_al.substr(num_head_sclips+num_lead, read_al.size()-num_head_sclips-num_lead-num_trail-num_back_sclips));
    new_alignment = Alignment(start_position, end_position, base_qualities, sequence, alignment_seq);

    // Add CIGAR data while approriately trimming for clipped bases
    int head = num_head_sclips, tail = num_back_sclips;
    std::vector<BamTools::CigarOp>::iterator end_iter = cigar_list.end()-1;
    while (tail > end_iter->Length && end_iter != cigar_list.begin()){
      tail -= end_iter->Length;
      end_iter--;
    }
    for (std::vector<BamTools::CigarOp>::iterator cigar_iter = cigar_list.begin(); cigar_iter != end_iter; cigar_iter++){
      if (head >= cigar_iter->Length)
	head -= cigar_iter->Length;
      else if (head > 0){
	new_alignment.add_cigar_element(CigarElement(cigar_iter->Type, cigar_iter->Length-head));
	head = 0;
      }
      else
	new_alignment.add_cigar_element(CigarElement(cigar_iter->Type, cigar_iter->Length));
    }
    if (head+tail > end_iter->Length)
      printErrorAndDie("Can't trim CIGAR character as the trim amount exceeds the CIGAR's length");
    if (head+tail < end_iter->Length)
      new_alignment.add_cigar_element(CigarElement(end_iter->Type, end_iter->Length-head-tail));

    return left_aligned;
}
