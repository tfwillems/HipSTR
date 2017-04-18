#include <algorithm>
#include <climits>
#include <map>
#include <string>
#include <sstream>
#include <vector>

#include "../error.h"
#include "../stringops.h"

#include "AlignmentData.h"
#include "AlignmentViz.h"
#include "HapBlock.h"
#include "HTMLCreator.h"

void getMaxInsertionSizes(const std::vector<Alignment>& alignments,
			  std::map<int32_t,int>& max_insertions){
  max_insertions.clear();
  for (auto align_iter = alignments.begin(); align_iter != alignments.end(); ++align_iter){
    int32_t position = align_iter->get_start();
    const std::vector<CigarElement>& cigar_ops = align_iter->get_cigar_list();
    std::map<int32_t, int>::iterator insertion_iter;
    for (auto cigar_iter = cigar_ops.begin(); cigar_iter != cigar_ops.end(); ++cigar_iter){
      char type  = cigar_iter->get_type();
      int length = cigar_iter->get_num();
      switch(type){
      case 'M': case '=': case 'X': case 'D':
	position += length;
	break;
      case 'I':
	insertion_iter = max_insertions.find(position);
	if(insertion_iter == max_insertions.end())
	  max_insertions.insert(std::pair<long,int>(position, length));
	else
	  insertion_iter->second = std::max(length, insertion_iter->second);
	break;
      case 'S': case 'H':
	break;
      default:
	printErrorAndDie("Invalid cigar option encountered");
	break;	
      }
    }
  }
}

void overlayAlignments(const std::vector<Alignment>& alignments, std::map<int32_t,int>& max_insertions,
		       std::vector<std::string>& results, int32_t& min_start, int32_t& max_stop){
  min_start = INT_MAX;
  max_stop  = INT_MIN;
  for (auto align_iter = alignments.begin(); align_iter != alignments.end(); ++align_iter)
    min_start = std::min(min_start, align_iter->get_start());

  results.clear();
  max_insertions.clear();
  getMaxInsertionSizes(alignments, max_insertions);
  max_insertions.insert(std::pair<int32_t,int>(INT_MAX, -1)); // Dummy insertion so that there's always an insertion after the end of the reads
  
  for (auto align_iter = alignments.begin(); align_iter != alignments.end(); ++align_iter){
    std::stringstream result;
    const std::string& nucleotides = align_iter->get_sequence();
    const std::vector<CigarElement>& cigar_ops = align_iter->get_cigar_list();
    std::vector<CigarElement>::const_iterator cigar_iter = cigar_ops.begin();
    int32_t position = align_iter->get_start();
    int nuc_index    = 0;
    char type        = cigar_iter->get_type();
    int  length      = cigar_iter->get_num();
    
    // Left pad with space characters
    std::map<int32_t,int>::iterator insertion_iter = max_insertions.begin();
    for (int32_t i = min_start; i <= position; i++){
      if (i == insertion_iter->first){
	for (int j = 0; j < insertion_iter->second; j++)
	  result << SPACE_CHAR;
	++insertion_iter;
      }
      if (i != position)
	result << SPACE_CHAR;
    }
    
    while (true){
      if (length == 0){
	++cigar_iter;
	if (cigar_iter == cigar_ops.end())
	  break;
	type   = cigar_iter->get_type();
	length = cigar_iter->get_num();
      }
      
      int num_rem_ins = 0;
      if (position == insertion_iter->first){
	int num_ins = insertion_iter->second;
	if (type == 'I'){
	  num_ins    -= length;
	  num_rem_ins = num_ins;
	}
	else
	  for (int i = 0; i < num_ins; i++)
	    result << NOT_APP_CHAR;
	++insertion_iter;
      }
      
      if (type == 'M' || type == '=' || type == 'X'){
	int num_match = std::min(length, insertion_iter->first - position);
	for (int i = 0; i < num_match; i++){
	  result << nucleotides[nuc_index++];
	  position++;
	  length--;
	}
      }
      else if (type == 'I'){
	for (int i = 0; i < length; i++)
	  result << nucleotides[nuc_index++];
	length = 0;
	if (cigar_iter+1 != cigar_ops.end())
	  for (int i = 0; i < num_rem_ins; i++)
	    result << NOT_APP_CHAR;
	num_rem_ins = 0;
      }
      else if (type == 'D'){
	int num_del = std::min(length, insertion_iter->first - position);
	for (int i = 0; i < num_del; i++){
	  result << DELETION_CHAR;
	  position++;
	  length--;
	}
      }
      else if (type == 'H')
	length = 0;
      else if (type == 'S'){
	nuc_index += length;
	length = 0;
      }
      else
	printErrorAndDie("Invalid cigar option encountered.");
    }
    max_stop = std::max(max_stop, position-1);
    results.push_back(result.str());
  }
}

std::string arrangeReferenceString(const std::string& chrom_seq, std::map<int32_t,int>& max_insertions,
				   const std::string& locus_id, int32_t str_start, int32_t str_stop,
				   int32_t min_start, int32_t max_stop, bool draw_locus_id,
				   std::ostream& output){
  // Hack to deal with lobSTR VCF off by 1
  int offset = 1;

  std::stringstream ref_result;
  std::vector<bool> within_locus;
  std::map<int32_t,int>::iterator insertion_iter = max_insertions.begin();
  for (int32_t i = min_start; i <= max_stop; i++){
    if (i == insertion_iter->first){
      for (int j = 0; j < insertion_iter->second; j++){
	ref_result << NOT_APP_CHAR;
	within_locus.push_back(false);
      }
      ++insertion_iter;
    }
    ref_result << chrom_seq[i];
    within_locus.push_back((i>=str_start-offset && i <= str_stop-offset));
  } 
  
  std::string ref_alignment = ref_result.str();
  writeReferenceString(ref_alignment, output, locus_id, within_locus, draw_locus_id);
  return ref_alignment;
}

void get_alignment_bounds(const std::vector<Alignment>& alignments,
			  int32_t& min_start, int32_t& max_stop){
  min_start = INT_MAX; max_stop = INT_MIN;  
  for (auto iter = alignments.begin(); iter != alignments.end(); ++iter){
    min_start = std::min(min_start, iter->get_start());
    max_stop  = std::max(max_stop,  iter->get_stop());
  }
}

void visualizeAlignments(const std::vector< std::vector<Alignment> >& alns, const std::vector<std::string>& sample_names,
			 const std::map<std::string, std::string>& sample_info, const std::vector<HapBlock*>& hap_blocks,
			 const std::string& chrom_seq, const std::string& locus_id, bool draw_locus_id,
			 std::ostream& output) {
  assert(hap_blocks.size() == 3 && alns.size() == sample_names.size());


  // Sort samples by name
  std::vector<std::pair<std::string, int> > sample_ordering;
  for (unsigned int i = 0; i < sample_names.size(); i++)
    sample_ordering.push_back(std::pair<std::string,int>(sample_names[i], i));
  std::sort(sample_ordering.begin(), sample_ordering.end());

  std::vector<Alignment> alignments;
  std::vector<std::string> alignment_samples;
  for (unsigned int i = 0; i < sample_ordering.size(); i++){
    std::string& sample = sample_ordering[i].first;
    const std::vector<Alignment>& sample_alns = alns[sample_ordering[i].second];
    alignments.insert(alignments.end(), sample_alns.begin(), sample_alns.end());
    for (unsigned int j = 0; j < sample_alns.size(); j++)
      alignment_samples.push_back(sample);
  }

  // Minimum and maximum coordinates of alignments
  int32_t min_coord, max_coord;
  get_alignment_bounds(alignments, min_coord, max_coord);

  // Arrange alignments
  std::vector<std::string> align_results;
  std::map<int32_t,int> max_insertions;
  int32_t min_start, max_stop;
  overlayAlignments(alignments, max_insertions, align_results, min_start, max_stop);

  // Arrange reference sequence
  std::string ref_alignment = arrangeReferenceString(chrom_seq, max_insertions, locus_id, 
						     hap_blocks[1]->start(), hap_blocks[1]->end(), min_start, max_stop,
						     draw_locus_id, output);

  // Write to HTML
  writeAlignmentStrings(ref_alignment, output, locus_id, align_results, alignment_samples, sample_info, true);
  output << locus_id << "\t" << "ALL" << "\t"
	 << "\t</table> "
	 << "<br> "
	 << "</div>" << std::endl;
}
