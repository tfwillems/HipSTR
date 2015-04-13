#include <algorithm>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <vector>

#include "AlignmentData.h"
#include "HaplotypeCreator.h"
#include "error.h"
#include "RepeatRegion.h"
#include "StringOps.h"

// TO DO: New region identification strategy
// Identify regions on a per-sample basis:
//     1. Min fraction
//     2. Min count (>= 2)
//     3. Distance from end of alignment...
// Merge per-sample regions into global regions
// Extract sequences
// Refilter sequences based on global stats:
//     1. Min fraction (5%)
//     2. Min count
// Consider merging close mismatch windows (< 5 bp)


const int MIN_COUNTS_PER_SAMPLE = 2;



/* 
   Determines the disjoint set of regions within which indels, mismatches or deletions occur. Regions that are overlapping or adjacent
   are combined into a single larger region. Each region is represented as [start, end) and stored in the provided
   vector in order of increasing start coordinate
*/
void getHaplotypeRegions(std::vector<Alignment>& alignments, 
			 std::vector< std::pair<int32_t, int32_t> >& regions){
  if (alignments.size() == 0)
    return;
  
  std::map<int32_t,int32_t> max_regions;
  std::map<int32_t, int>::iterator riter;
  std::map<int32_t, int> mismatch_counts; 
  std::string current_sample = alignments[0].get_sample();

  for (auto align_iter = alignments.begin(); ; align_iter++){
    // Process set of mismatches for previous sample
    if (align_iter == alignments.end() || align_iter->get_sample().compare(current_sample) != 0){
      for (auto mm_iter = mismatch_counts.begin(); mm_iter != mismatch_counts.end(); mm_iter++){
	if (mm_iter->second >= MIN_COUNTS_PER_SAMPLE){
	  if (max_regions.find(mm_iter->first) == max_regions.end())
	    max_regions[mm_iter->first] = mm_iter->first+1; 
	}
      }
      
      if (align_iter == alignments.end())
	break;
      else
	current_sample = align_iter->get_sample();
      mismatch_counts.clear();
    }
    
    int32_t pos = align_iter->get_start();
    for (auto cigar_iter = align_iter->get_cigar_list().begin(); cigar_iter != align_iter->get_cigar_list().end(); cigar_iter++){
      switch(cigar_iter->get_type()){
      case '=':
	pos += cigar_iter->get_num();
	break;
      case 'I':
	riter = max_regions.find(pos);
	if (riter == max_regions.end())
	  max_regions[pos] = pos;
	break;
      case 'X':
	for (int i = 0; i < cigar_iter->get_num(); i++)
	  mismatch_counts[i+pos]++;
	pos += cigar_iter->get_num();
	break;
      case 'D':
	riter = max_regions.find(pos);
	if (riter != max_regions.end())
	  riter->second = std::max(riter->second, pos+cigar_iter->get_num());
	else
	  max_regions[pos] = pos+cigar_iter->get_num();
	pos += cigar_iter->get_num();
	break;
      default:
	printErrorAndDie("Invalid CIGAR char in getHaplotypeRegions()");
	break;
      }
    }
  }

 
  // Reduce set of regions by merging ones that overlap or are adjacent to one another
  regions.clear();
  if (max_regions.size() != 0){
    std::map<int32_t,int32_t>::iterator riter = max_regions.begin();
    int32_t start = riter->first;
    int32_t stop  = riter->second;
    riter++;
    while (riter != max_regions.end()){
      if (riter->first <= stop)
	stop = std::max(stop, riter->second);
      else {
	regions.push_back(std::pair<int,int>(start, stop));
	start = riter->first;
	stop  = riter->second;
      }
      riter++;
    }
    regions.push_back(std::pair<int32_t,int32_t>(start, stop));
  }
}

/* 
   Extracts the set of sequences observed for each of the regions in the provided alignments.
   An alignment's sequence for a region is only considered if the alignment 
   fully spans the region. By default, the reference sequence for each region is also added
   to each region's set.
 */
void extractRegionSequences(std::vector<Alignment>& alignments, 
			    std::vector< std::pair<int32_t,int32_t> >& regions, 
			    std::string& ref_sequence, 
			    std::vector< std::vector<std::string> >& region_seqs){ 
  std::vector< std::map<std::string, int> > regseq_counts;
  int region_counts[regions.size()];
  for (int i = 0; i < regions.size(); i++){
    region_counts[i] = 0; 
    regseq_counts.push_back(std::map<std::string, int>());
  }
 
  int alignment_count = 0;
  for (std::vector<Alignment>::iterator align_iter = alignments.begin(); align_iter != alignments.end(); align_iter++){
    alignment_count++;

    int32_t pos = align_iter->get_start();
    std::vector<CigarElement>::const_iterator cigar_iter = align_iter->get_cigar_list().begin();
    std::vector< std::pair<int32_t,int32_t> >::iterator region_iter = regions.begin();
    
    int align_index = 0; // Index into alignment string
    int char_index  = 0; // Index of current base in current CIGAR element
    int reg_index   = 0; // Index of current region

    // Skip any leading regions not fully spanned by alignment
    while (region_iter != regions.end() && (region_iter->first < pos)){
      region_iter++;
      reg_index++;
    }

    // Extract each region sequence, if fully spanned by alignment
    std::stringstream reg_seq;
    while (region_iter != regions.end() && cigar_iter != align_iter->get_cigar_list().end()){
      if (char_index == cigar_iter->get_num()){
	cigar_iter++;
	char_index = 0;
      }
      else if (pos > region_iter->second){
	if (reg_seq.str() == "")
	    reg_seq << "X";
	region_counts[reg_index]++;

	std::string seq = uppercase(reg_seq.str());
	std::map<std::string, int>::iterator iter = regseq_counts[reg_index].find(seq);
	if (iter != regseq_counts[reg_index].end())
	  iter->second++;
	else
	  regseq_counts[reg_index].insert(std::pair<std::string, int>(seq, 1));

	reg_index++;
	region_iter++;
	reg_seq.str("");
      }
      else if (pos == region_iter->second){
	if (cigar_iter->get_type() == 'I'){
	  reg_seq << align_iter->get_alignment().substr(align_index, cigar_iter->get_num());
	  align_index += cigar_iter->get_num();
	  char_index = 0;
	  cigar_iter++;
	}
	else {
	  if (reg_seq.str() == "")
	    reg_seq << "X";
	  region_counts[reg_index]++;

	  std::string seq = reg_seq.str();
	  std::map<std::string, int>::iterator iter = regseq_counts[reg_index].find(seq);
	  if (iter != regseq_counts[reg_index].end())
	    iter->second++;
	  else
	    regseq_counts[reg_index].insert(std::pair<std::string, int>(seq, 1));

	  reg_index++;
	  region_iter++;
	  reg_seq.str("");
	}
      }
      else if (pos >= region_iter->first){
	int32_t num_bases = std::min(region_iter->second-pos, cigar_iter->get_num()-char_index);
	switch(cigar_iter->get_type()){	 
	case 'I':
	  // Insertion within region
	  num_bases = cigar_iter->get_num();
	  reg_seq << align_iter->get_alignment().substr(align_index, num_bases);
	  break;
	case '=': case 'X':
	  reg_seq << align_iter->get_alignment().substr(align_index, num_bases);
	  pos += num_bases;
	  break;
	case 'D':
	  pos += num_bases;
	  break;
	default:
	  printErrorAndDie("Invalid CIGAR char in extractRegionSequences()");
	  break;
	}
	align_index += num_bases;
	char_index  += num_bases;
      }
      else {
	int32_t num_bases;
	if (cigar_iter->get_type() == 'I')
	  num_bases = cigar_iter->get_num()-char_index;
	else {
	  num_bases = std::min(region_iter->first-pos, cigar_iter->get_num()-char_index);
	  pos      += num_bases;
	}
	align_index  += num_bases;
	char_index   += num_bases;
      }
    }
  }

  region_seqs.clear();
  for (unsigned int i = 0; i < regions.size(); i++)
    region_seqs.push_back(std::vector<std::string>());

  // Add reference sequence for each region to the respective set
  for (unsigned int i = 0; i < regions.size(); i++){
    if (regions[i].first == regions[i].second)
      region_seqs[i].push_back("X");
    else
      region_seqs[i].push_back(uppercase(ref_sequence.substr(regions[i].first, regions[i].second-regions[i].first)));
  }
   
  // Filter potential elements based on the fraction of reads in which they're observed and the number of times they're observed
  // TO DO: Add more refined filters on a per-sample basis
  double min_frac = 0.01;
  int min_count   = 2;
  for (int i = 0; i < regions.size(); i++){
    for (std::map<std::string, int>::iterator iter = regseq_counts[i].begin(); iter != regseq_counts[i].end(); iter++){
      if (iter->first.compare(region_seqs[i][0]) != 0 && iter->second >= min_frac*region_counts[i] && iter->second >= min_count)
	region_seqs[i].push_back(iter->first);
    }
  }
}

