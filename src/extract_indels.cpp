#include "error.h"
#include "extract_indels.h"

bool ExtractCigar(const std::vector<CigarElement>& cigar_data, const int& cigar_start,
		  const int& region_start, const int& region_end,
		  int& bp_diff_from_ref) {
  std::vector<CigarOp> cigar_ops;
  for (auto cigar_iter = cigar_data.begin(); cigar_iter != cigar_data.end(); cigar_iter++)
    cigar_ops.push_back(CigarOp(cigar_iter->get_type(), cigar_iter->get_num()));
  return ExtractCigar(cigar_ops, cigar_start, region_start, region_end, bp_diff_from_ref);
}

/* Returns true iff the function successfully extracted the base pair difference
   of the read from the reference genome within the provided region.
   If successful, stores the difference using the provided int reference.
   Adapted from same function in lobSTR
 */
bool ExtractCigar(const std::vector<CigarOp>& cigar_data, const int& cigar_start,
		  const int& region_start, const int& region_end,
		  int& bp_diff_from_ref) {
  assert(cigar_start >= 0 && region_end >= region_start);

  int pos = cigar_start;        // beginning position of current cigar item
  size_t cigar_start_index = 0; // position in cigar string where region starts
  size_t last_match_index = 0;  // position of last 'M' cigar
  int bp = 0;
  char cigar_type; // type of the cigar item

  // Determine cigar end
  int cigar_region_length = 0;
  for (size_t i = 0; i < cigar_data.size(); i++) {
    bp = cigar_data[i].Length;
    cigar_type = cigar_data[i].Type;
    if (cigar_type == 'M' || cigar_type == '=' || cigar_type == 'X' || cigar_type == 'D')
      cigar_region_length += bp;
  }

  // Checks for boundaries
  if (region_start < cigar_start) return false;
  if (region_end   >= cigar_start+cigar_region_length) return false;

  // Set start index
  while (pos < region_start && cigar_start_index < cigar_data.size()) {
    bp = cigar_data[cigar_start_index].Length;
    cigar_type = cigar_data[cigar_start_index].Type;
    // If match or del, increment position. All other CIGAR strings (I, P, H, S, N) don't increment
    if (cigar_type == 'M' || cigar_type == '=' || cigar_type == 'X' || cigar_type == 'D')
      pos += bp;
    if (cigar_type == 'M' || cigar_type == '=' || cigar_type == 'X')
      last_match_index = cigar_start_index;
    cigar_start_index++;
  }
  cigar_start_index = last_match_index;
  if (cigar_start_index == 0){
    char type = cigar_data[cigar_start_index].Type;
    if (!(type == 'M' || type == '=' || type == 'X'))
      return false;
  }
    
  // *** Set end index *** //
  size_t cigar_end_index = cigar_data.size() - 1;
  last_match_index = cigar_data.size() - 1;
  pos = cigar_start + cigar_region_length;
  while (pos > region_end) {
    bp = cigar_data[cigar_end_index].Length;
    cigar_type = cigar_data[cigar_end_index].Type;
    // If match or del, decrement position. All other CIGAR strings (I, P, H, S, N) don't decrement
    if (cigar_type == 'M' || cigar_type == '=' || cigar_type == 'X' || cigar_type == 'D')
      pos -= bp;
    if (cigar_type == 'M' || cigar_type == '=' || cigar_type == 'X')
      last_match_index = cigar_end_index;
    if (cigar_end_index == 0) break;
    cigar_end_index -= 1;
  }
  cigar_end_index = last_match_index;
  if (cigar_end_index == cigar_data.size()-1){
    char type = cigar_data[cigar_end_index].Type;
    if (!(type == 'M' || type == '=' || type == 'X'))
      return false;
  }

  bp_diff_from_ref = 0;
  for (size_t i = cigar_start_index; i <= cigar_end_index; i++){
    if (cigar_data[i].Type == 'D')
      bp_diff_from_ref -= cigar_data[i].Length;
    else if (cigar_data[i].Type == 'I')
      bp_diff_from_ref += cigar_data[i].Length;
  }
  return true;
}
