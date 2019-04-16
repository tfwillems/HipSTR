#include <assert.h>

#include "AlignmentOps.h"
#include "../error.h"
#include "NeedlemanWunsch.h"
#include "../stringops.h"

extern const int ALIGN_WINDOW_WIDTH = 75;

/*
 Realign read to reference region using left alignment variant. Store the new alignment information using
 the provided Alignment reference. Converts the bases to their upper case variants
 */
bool realign(BamAlignment& alignment, const std::string& ref_sequence, Alignment& new_alignment){
    int32_t start        = std::max(alignment.Position()-ALIGN_WINDOW_WIDTH-1, 0);
    int32_t stop         = std::min(alignment.GetEndPosition()+ALIGN_WINDOW_WIDTH-1, (int32_t)(ref_sequence.size()-1));
    int32_t length       = stop-start+1;
    std::string ref_seq  = ref_sequence.substr(start, length);
    std::string read_seq = alignment.QueryBases();
    
    // Realign read using left-alignment method
    std::string ref_al, read_al;
    float score;
    std::vector<CigarOp> cigar_list;
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
    for (auto cigar_iter = cigar_list.begin(); cigar_iter != cigar_list.end() && !halt; cigar_iter++){
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
	msg << "Invalid CIGAR character " << cigar_iter->Type << " in realign() for alignment " << alignment.Name();
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
    if (alignment.QueryBases().size() != alignment.Qualities().size())
        printErrorAndDie("Lengths of sequence and quality strings don't match");
    std::string base_qualities = alignment.Qualities().substr(num_head_sclips, read_seq.size()-num_head_sclips-num_back_sclips);
    std::string sequence       = uppercase(read_seq.substr(num_head_sclips, read_seq.size()-num_head_sclips-num_back_sclips));
    std::string alignment_seq  = uppercase(read_al.substr(num_head_sclips+num_lead, read_al.size()-num_head_sclips-num_lead-num_trail-num_back_sclips));
    new_alignment = Alignment(start_position, end_position, alignment.IsReverseStrand(), alignment.Name(), base_qualities, sequence, alignment_seq);

    // Add CIGAR data while approriately trimming for clipped bases
    int head = num_head_sclips, tail = num_back_sclips;
    auto end_iter = cigar_list.end()-1;
    while (tail > end_iter->Length && end_iter != cigar_list.begin()){
      tail -= end_iter->Length;
      end_iter--;
    }
    for (auto cigar_iter = cigar_list.begin(); cigar_iter != end_iter; cigar_iter++){
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

void convertAlignment(BamAlignment& alignment, const std::string& ref_sequence, Alignment& new_alignment){
  std::string read_sequence = uppercase(alignment.QueryBases());
  int32_t seq_index = 0, ref_index = alignment.Position();
  std::stringstream aln_ss;
  new_alignment = Alignment(alignment.Position(), alignment.GetEndPosition()-1, alignment.IsReverseStrand(), alignment.Name(), alignment.Qualities(), read_sequence, "");
  for (auto cigar_iter = alignment.CigarData().begin(); cigar_iter != alignment.CigarData().end(); cigar_iter++){
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
