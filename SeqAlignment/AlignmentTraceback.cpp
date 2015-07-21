#include <algorithm>
#include <assert.h>

#include "AlignmentTraceback.h"
#include "AlignmentData.h"

std::string stitch(const std::string& hap_aln, const std::string& read_aln, int h_index, int r_index, int increment){
  std::stringstream stitched_aln;
  while (r_index >= 0 && r_index < read_aln.size()){
    if (read_aln[r_index] == 'S'){
      stitched_aln << 'S';
      r_index += increment;
      continue;
    }
    
    assert(h_index >= 0 && h_index < hap_aln.size());
    if (hap_aln[h_index] == 'D'){
      if (read_aln[r_index] == 'I'){
	stitched_aln << 'M';
	r_index += increment;
	h_index += increment;
      }
      else {
	stitched_aln << 'D';
	h_index += increment;
      }
    }    
    else if (read_aln[r_index] == 'I'){
      stitched_aln << 'I';
      r_index += increment;
    }
    else if (read_aln[r_index] == 'D'){
      if (hap_aln[h_index] == 'M')
	stitched_aln << 'D';
      else if (hap_aln[h_index] == 'I')
	stitched_aln << "";//'M';
      else
	printErrorAndDie("Logical error in stitch_alignment_trace()");
      r_index += increment;
      h_index += increment;
    }
    else if (read_aln[r_index] == 'M'){
      if (hap_aln[h_index] != 'M' && hap_aln[h_index] != 'I')
	printErrorAndDie("Logical error in stitch_alignment_trace()");
      stitched_aln << hap_aln[h_index];
      r_index += increment;
      h_index += increment;
    }
    else
      printErrorAndDie("Logical error in stitch_alignment_trace()");
   }
   return stitched_aln.str();
}

void stitch_alignment_trace(int32_t hap_start, const std::string& hap_aln_to_ref, const std::string& read_aln_to_hap, 
			    int hap_index, int seed_base, Alignment& orig_aln,
			    Alignment& new_aln){
  int hap_aln_index = 0;
  int32_t seed_pos  = hap_start;

  while (hap_index > 0 && hap_aln_index < hap_aln_to_ref.size()){
    if (hap_aln_to_ref[hap_aln_index] == 'M' || hap_aln_to_ref[hap_aln_index] == 'I')
      hap_index--;
    if (hap_aln_to_ref[hap_aln_index] == 'M' || hap_aln_to_ref[hap_aln_index] == 'D')
      seed_pos++;
    hap_aln_index++;
  }
  assert(hap_aln_index != hap_aln_to_ref.size());

  int read_aln_index = 0;
  while (seed_base > 0 && read_aln_index < read_aln_to_hap.size()){
    if (read_aln_to_hap[read_aln_index] == 'M' || read_aln_to_hap[read_aln_index] == 'I' || read_aln_to_hap[read_aln_index] == 'S')
      seed_base--;
    read_aln_index++;
  }
  assert(read_aln_index != read_aln_to_hap.size());
    
  std::string left_aln  = stitch(hap_aln_to_ref, read_aln_to_hap, hap_aln_index-1, read_aln_index-1, -1);
  std::reverse(left_aln.begin(), left_aln.end());
  std::string right_aln = stitch(hap_aln_to_ref, read_aln_to_hap, hap_aln_index+1, read_aln_index+1,  1);
  std::string full_aln  = left_aln + "M" + right_aln;

  for (int i = 0; i < full_aln.size(); i++){
    if (full_aln[i] == 'I')
      full_aln[i] = 'S';
    else
      break;
  }

  // Determine alignment start and end coordinates
  int32_t start = seed_pos;
  int32_t stop  = seed_pos;
  for (auto iter = left_aln.begin(); iter != left_aln.end(); iter++)
    if (*iter == 'D' || *iter == 'M')
      start--;
  for (auto iter = right_aln.begin(); iter != right_aln.end(); iter++)
    if (*iter == 'D' || *iter == 'M')
      stop++;

  // Construct the CIGAR string elements
  std::vector<CigarElement> cigar_list;
  char cigar_char = full_aln[0];
  int  num        = 1;
  char new_cigar_char;
  for(unsigned int i = 1; i < full_aln.size(); i++){
    new_cigar_char = full_aln[i];
    if (new_cigar_char != cigar_char){
      cigar_list.push_back(CigarElement(cigar_char, num));
      num = 1;
      cigar_char = new_cigar_char;
    }
    else
      num += 1;
  }
  cigar_list.push_back(CigarElement(cigar_char, num));

  // Construct the actual alignment string (from the string describing the alignment operations)
  int read_index = 0;
  std::stringstream aln_ss;
  const std::string& bases = orig_aln.get_sequence();
  for (unsigned int i = 0; i < full_aln.size(); i++){
    switch (full_aln[i]){
    case 'S':
      read_index++;
      break;
    case 'M':
    case 'I':
      aln_ss << bases[read_index];
      read_index++;
      break;
    case 'D':
      aln_ss << "-";
      break;
    default:
      printErrorAndDie("Invalid character encountered in stitch_alignment_trace()");
      break;
    }
  }

  new_aln = Alignment(start, stop, orig_aln.get_sample(), orig_aln.get_base_qualities(), orig_aln.get_sequence(), aln_ss.str());
  new_aln.set_cigar_list(cigar_list);
}
