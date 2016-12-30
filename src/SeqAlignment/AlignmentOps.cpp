#include <assert.h>

#include "AlignmentOps.h"
#include "../error.h"
#include "NeedlemanWunsch.h"
#include "../stringops.h"

extern const int ALIGN_WINDOW_WIDTH = 75;

bool GetFloatBamTag(const BamTools::BamAlignment& aln, const std::string& tag_name, float* destination){
  char tag_type;
  if (!aln.GetTagType(tag_name, tag_type))
    return false;
  if (tag_type != BamTools::Constants::BAM_TAG_TYPE_FLOAT)
    return false;
  float val;
  if (!aln.GetTag(tag_name, val))
    return false;
  *destination = val;
  return true;
}

bool GetIntBamTag(const BamTools::BamAlignment& aln, const std::string& tag_name, int* destination){
  char tag_type;
  if (!aln.GetTagType(tag_name, tag_type)) return false;
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
    bool aligned = NeedlemanWunsch::Align(ref_seq, read_seq, ref_al, read_al, &score, cigar_list);
    
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
	msg << "Invalid CIGAR character " << cigar_iter->Type << " in realign() for alignment " << alignment.Name;
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
    new_alignment = Alignment(start_position, end_position, alignment.Name, base_qualities, sequence, alignment_seq);

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

    return aligned;
}


bool startsWithSoftClip(const BamTools::BamAlignment& aln){
  if (aln.CigarData.size() == 0)
    return false;
  return aln.CigarData.front().Type == 'S';
}

bool endsWithSoftClip(const BamTools::BamAlignment& aln){
  if (aln.CigarData.size() == 0)
    return false;
  return aln.CigarData.back().Type == 'S';
}

bool startsWithHardClip(const BamTools::BamAlignment& aln){
  if (aln.CigarData.size() == 0)
    return false;
  return aln.CigarData.front().Type == 'H';
}

bool endsWithHardClip(const BamTools::BamAlignment& aln){
  if (aln.CigarData.size() == 0)
    return false;
  return aln.CigarData.back().Type == 'H';
}

/*
 *  Trim an alignment that extends too far upstream or downstream of the provided region or has low base qualities on the ends
 *  First trims MIN_LTRIM and MIN_RTRIM bases from the left and right ends of the alignment.
 *  It then trims until either i) the base quality exceeds the provided threshold or ii) the alignment is fully within the provided region bounds
 *  Modifies the following members of the BamAlignment: Position, CigarData, Qualities, AlignedBases, QueryBases, Length
 */
void trimAlignment(BamTools::BamAlignment& aln, int32_t min_read_start, int32_t max_read_stop, char min_base_qual, int32_t min_ltrim, int32_t min_rtrim){
  assert(min_ltrim + min_rtrim <= aln.QueryBases.size());
  int ltrim = 0, aligned_bases_ltrim = 0;
  int32_t start_pos = aln.Position;
  while ((ltrim < min_ltrim || start_pos < min_read_start) && aln.CigarData.size() > 0){
    // Check if we should stop trimming b/c the quality score is above the threshold
    bool qual_above_thresh = false;
    switch(aln.CigarData.front().Type){
    case 'M': case '=': case 'X': case 'I': case 'S':
      qual_above_thresh = (aln.Qualities[ltrim] > min_base_qual);
      break;
    case 'D': case 'H':
      break;
    default:
      printErrorAndDie("Invalid CIGAR option encountered in trimAlignment");
      break;
    }
    if (qual_above_thresh && ltrim >= min_ltrim)
      break;

    switch(aln.CigarData.front().Type){
    case 'M': case '=': case 'X':
      ltrim++;
      aligned_bases_ltrim++;
      start_pos++;
      break;
    case 'D':
      aligned_bases_ltrim++;
      start_pos++;
      break;
    case 'I':
      ltrim++;
      aligned_bases_ltrim++;
      break;
    case 'S':
      ltrim++;
      break;
    case 'H':
      break;
    default:
      printErrorAndDie("Invalid CIGAR option encountered in trimAlignment");
      break;
    }
    if (aln.CigarData.front().Length == 1)
      aln.CigarData.erase(aln.CigarData.begin(), aln.CigarData.begin()+1);
    else
      aln.CigarData.front().Length--;
  }
  aln.Position = start_pos;

  int rtrim = 0, aligned_bases_rtrim = 0, qual_string_len = aln.Qualities.size()-1;
  int32_t end_pos = aln.GetEndPosition();
  while ((rtrim < min_rtrim || end_pos > max_read_stop) && aln.CigarData.size() > 0){
    // Check if we should stop trimming b/c the quality score is above the threshold
    bool qual_above_thresh = false;
    switch(aln.CigarData.back().Type){
    case 'M': case '=': case 'X': case 'I': case 'S':
      qual_above_thresh = (aln.Qualities[qual_string_len-rtrim] > min_base_qual);
      break;
    case 'D': case 'H':
      break;
    default:
      printErrorAndDie("Invalid CIGAR option encountered in trimAlignment");
      break;
    }
    if (qual_above_thresh && rtrim >= min_rtrim)
      break;

    switch(aln.CigarData.back().Type){
    case 'M': case '=': case 'X':
      rtrim++;
      aligned_bases_rtrim++;
      end_pos--;
      break;
    case 'D':
      aligned_bases_rtrim++;
      end_pos--;
      break;
    case 'I':
      rtrim++;
      aligned_bases_rtrim++;
      break;
    case 'S':
      rtrim++;
      break;
    case 'H':
      break;
    default:
      printErrorAndDie("Invalid CIGAR option encountered in trimAlignment");
      break;
    }
    if(aln.CigarData.back().Length == 1)
      aln.CigarData.pop_back();
    else
      aln.CigarData.back().Length--;
  }

  assert(ltrim+rtrim <= aln.QueryBases.size());
  assert(aligned_bases_ltrim + aligned_bases_rtrim <= aln.AlignedBases.size());
  aln.QueryBases   = aln.QueryBases.substr(ltrim, aln.QueryBases.size()-ltrim-rtrim);
  aln.AlignedBases = aln.AlignedBases.substr(aligned_bases_ltrim, aln.AlignedBases.size()-aligned_bases_ltrim-aligned_bases_rtrim);
  aln.Qualities    = aln.Qualities.substr(ltrim, aln.Qualities.size()-ltrim-rtrim);
  aln.Length      -= (ltrim + rtrim);
}

void trimLowQualityEnds(BamTools::BamAlignment& aln, char min_base_qual){
  int32_t min_read_start  = aln.GetEndPosition()+1;
  int32_t max_read_stop   = aln.Position-1;
  return trimAlignment(aln, min_read_start, max_read_stop, min_base_qual);
}

bool matchesReference(const BamTools::BamAlignment& aln){
  for (auto cigar_iter = aln.CigarData.begin(); cigar_iter != aln.CigarData.end(); cigar_iter++)
    if (cigar_iter->Type != 'M' && cigar_iter->Type != '=')
      return false;
  return true;
}

void convertAlignment(BamTools::BamAlignment& alignment, std::string& ref_sequence, Alignment& new_alignment){
  std::string read_sequence = uppercase(alignment.QueryBases);
  int32_t seq_index = 0, ref_index = alignment.Position;
  std::stringstream aln_ss;
  new_alignment = Alignment(alignment.Position, alignment.GetEndPosition()-1, alignment.Name, alignment.Qualities, read_sequence, "");
  for (auto cigar_iter = alignment.CigarData.begin(); cigar_iter != alignment.CigarData.end(); cigar_iter++){
    int32_t cigar_index    = 0;
    char prev_cigar_type   = '=';
    int32_t prev_cigar_num = 0;

    switch (cigar_iter->Type){
    case 'H':
      break;
    case 'S':
      new_alignment.add_cigar_element(CigarElement(cigar_iter->Type, cigar_iter->Length));
      seq_index += cigar_iter->Length;
      break;
    case 'I':
      new_alignment.add_cigar_element(CigarElement(cigar_iter->Type, cigar_iter->Length));
      aln_ss << read_sequence.substr(seq_index, cigar_iter->Length);
      seq_index += cigar_iter->Length;
      break;
    case 'D':
      new_alignment.add_cigar_element(CigarElement(cigar_iter->Type, cigar_iter->Length));
      aln_ss << std::string(cigar_iter->Length, '-');
      ref_index += cigar_iter->Length;
      break;
    case 'M': case '=': case 'X':
      while (cigar_index < cigar_iter->Length){
	if (read_sequence[seq_index] == static_cast<char>(toupper(ref_sequence[ref_index]))){
	  if (prev_cigar_type == '=')
	    prev_cigar_num++;
	  else {
	    if (prev_cigar_num != 0)
	      new_alignment.add_cigar_element(CigarElement(prev_cigar_type, prev_cigar_num));
	    prev_cigar_type = '=';
	    prev_cigar_num  = 1;
	  }
	}
	else {
	  if (prev_cigar_type == 'X')
	    prev_cigar_num++;
	  else {
	    if (prev_cigar_num != 0)
	      new_alignment.add_cigar_element(CigarElement(prev_cigar_type, prev_cigar_num));
	    prev_cigar_type = 'X';
	    prev_cigar_num  = 1;
	  }
	}
	aln_ss << read_sequence[seq_index];
	cigar_index++; ref_index++; seq_index++;
      }

      if (prev_cigar_num != 0)
	new_alignment.add_cigar_element(CigarElement(prev_cigar_type, prev_cigar_num));
      break;
    default:
      printErrorAndDie("Invalid CIGAR option encountered in convertAlignment");
      break;
    }
  }

  new_alignment.set_alignment(aln_ss.str());
  assert(seq_index == read_sequence.size());
  assert(ref_index == alignment.GetEndPosition());
}
