#include "AlignmentOps.h"
#include "../error.h"
#include "NWNoRefEndPenalty.h"
#include "../stringops.h"

extern const std::string START_TAG  = "XS";
extern const std::string STOP_TAG   = "XE";
extern const std::string RG_TAG     = "RG";
extern const std::string SAMPLE_TAG = "SN";
extern const std::string MOTIF_TAG  = "XR";
extern const int ALIGN_WINDOW_WIDTH = 75;


bool compareAlignments(const BamTools::BamAlignment& alignment_1, const BamTools::BamAlignment& alignment_2){
  std::string samp_1, samp_2;
  bool got_tag = alignment_1.GetTag(SAMPLE_TAG, samp_1);
  if (!got_tag) printErrorAndDie("Failed to retrieve sample tag for alignment");
  got_tag = alignment_2.GetTag(SAMPLE_TAG, samp_2);
  if (!got_tag) printErrorAndDie("Failed to retrieve sample tag for alignment");

  int sample_comp = samp_1.compare(samp_2);
  if (sample_comp != 0)
    return sample_comp < 0;

  if (alignment_1.Position != alignment_2.Position)
    return alignment_1.Position < alignment_2.Position;

  return alignment_1.GetEndPosition() < alignment_2.GetEndPosition();
}


void sortAlignments(std::vector<BamTools::BamAlignment>& alignments){
  std::sort(alignments.begin(), alignments.end(), compareAlignments); 
}

std::string getSampleName(const BamTools::BamAlignment& alignment){
  if (!alignment.HasTag("RG"))
    printErrorAndDie("Alignment does not have read group tag (RG).");

  std::string rg;
  alignment.GetTag("RG", rg);

  // Split RG by semicolon. Should result in three substrings
  std::vector<std::string> tokens;
  split_by_delim(rg, ';', tokens);
  if (tokens.size() != 3)
    printErrorAndDie("RG tag split by ';' delimiter should result in 3 tokens.");
  return tokens[1];
}

std::string getCigarString(BamTools::BamAlignment& alignment){
  std::stringstream cigar_str;
  for (std::vector<BamTools::CigarOp>::iterator iter = alignment.CigarData.begin(); iter != alignment.CigarData.end(); iter++)
    cigar_str << iter->Length << iter->Type;
  return cigar_str.str();
}

std::string getCigarString(std::vector<BamTools::CigarOp>& cigar_list){
  std::stringstream cigar_str;
  for (std::vector<BamTools::CigarOp>::iterator iter = cigar_list.begin(); iter != cigar_list.end(); iter++)
    cigar_str << iter->Length << iter->Type;
  return cigar_str.str();
}


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


bool GetStringBamTag(const BamTools::BamAlignment& alignment, const std::string& tag_name, std::string& value){
  char tag_type = 'Z';
  if (!alignment.GetTagType(tag_name, tag_type)) return false;
  if (tag_type == BamTools::Constants::BAM_TAG_TYPE_STRING)
    return alignment.GetTag(tag_name, value);
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
    std::string sample;
    alignment.GetTag(SAMPLE_TAG, sample);
    
    if (alignment.QueryBases.size() != alignment.Qualities.size())
        printErrorAndDie("Lengths of sequence and quality strings don't match");
    std::string base_qualities = alignment.Qualities.substr(num_head_sclips, read_seq.size()-num_head_sclips-num_back_sclips);
    std::string sequence       = uppercase(read_seq.substr(num_head_sclips, read_seq.size()-num_head_sclips-num_back_sclips));
    std::string alignment_seq  = uppercase(read_al.substr(num_head_sclips+num_lead, read_al.size()-num_head_sclips-num_lead-num_trail-num_back_sclips));
    new_alignment = Alignment(start_position, end_position, sample, base_qualities, sequence, alignment_seq);

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
    

    if (alignment.Name.compare("HS2000-1266_147:5:2102:3686:41029") == 0){
      std::cerr << "DEBUG INFO:" << std::endl
		<< alignment.QueryBases << std::endl
		<< num_head_sclips << " " << num_back_sclips << " " << num_lead << " " << num_trail << std::endl
		<< read_seq << std::endl
		<< read_al << std::endl
		<< ref_al  << std::endl
		<< start_position << " " << end_position << std::endl
		<< left_aligned << std::endl;
      for (std::vector<BamTools::CigarOp>::iterator cigar_iter = cigar_list.begin(); cigar_iter != cigar_list.end(); cigar_iter++)
	std::cerr <<  cigar_iter->Length <<cigar_iter->Type; 
      std::cerr << std::endl;
    }
    return left_aligned;
}
