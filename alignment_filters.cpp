/*
Copyright (C) 2014 Thomas Willems <twillems@mit.edu>
*/

#include <vector>
#include <string>
#include <sstream>

#include "alignment_filters.h"
#include "error.h"
#include "zalgorithm.h"

using namespace std;

namespace AlignmentFilters {
  template<typename CigarIterator> int GetDistToIndel(CigarIterator iter, CigarIterator end){
    // Process leading clipping ops
    if (iter != end && iter->Type == 'H')
      iter++;
    if (iter != end && iter->Type == 'S')
      iter++;
    
    int dist = 0;
    while (iter != end){
      char type = iter->Type;
      if (type == 'M')
	dist += iter->Length;
      else if (type == 'I' || type == 'D')
	return dist;
      else if (type == 'S' || type == 'H')
	return -1;
      else {
	string msg = "Invalid CIGAR char";
	msg += type;
	printErrorAndDie(msg);
      }
      iter++;
    }
    return -1;
  }

  std::string GetCigarString(std::vector<BamTools::CigarOp>& cigar_ops){
    std::stringstream ss;
    for (auto iter = cigar_ops.begin(); iter != cigar_ops.end(); iter++)
      ss << iter->Length << iter->Type;
    return ss.str();
  }
  
  pair<int,int> GetEndDistToIndel(BamTools::BamAlignment& aln){
    vector<int> vals;
    int head_dist = GetDistToIndel(aln.CigarData.begin(),  aln.CigarData.end());
    int tail_dist = GetDistToIndel(aln.CigarData.rbegin(), aln.CigarData.rend());
    return pair<int,int>(head_dist, tail_dist);
  }
  
  pair<int,int> GetNumEndMatches(BamTools::BamAlignment& aln, const string& ref_seq, int ref_seq_start){
    if (aln.Position < ref_seq_start)
      return pair<int,int>(-1,-1);
    
    unsigned int read_index = 0;
    unsigned int ref_index  = aln.Position-ref_seq_start;
    vector<BamTools::CigarOp>::iterator cigar_iter = aln.CigarData.begin();
    bool beginning = true;
    int match_run  = 0;
    int head_match = 0;

    // Process leading clip CIGAR types
    if (cigar_iter != aln.CigarData.end() && cigar_iter->Type == 'H')
      cigar_iter++;
    if (cigar_iter != aln.CigarData.end() && cigar_iter->Type == 'S'){
      read_index += cigar_iter->Length;
      cigar_iter++;
    }
    
    // Process CIGAR items as long as read region lies within reference sequence bounds
    while (cigar_iter != aln.CigarData.end() && ref_index < ref_seq.size() && read_index < aln.QueryBases.size()){
      if (cigar_iter->Type == 'M'){
	if (ref_index + cigar_iter->Length > ref_seq.size()) 
	  return pair<int,int>(-1, -1);
	if (read_index + cigar_iter->Length > aln.QueryBases.size())
	  printErrorAndDie("Nucleotides for aligned read don't correspond to the CIGAR string");
	for (unsigned int len = cigar_iter->Length; len > 0; len--){
	  if ((char)tolower(ref_seq[ref_index]) == (char)tolower(aln.QueryBases[read_index]))
	    match_run++;
	  else {
	    if (beginning) head_match = match_run;
	    beginning = false;
	    match_run = 0;
	  }
	  read_index++;
	  ref_index++;
	}
      }
      else if (cigar_iter->Type == 'I'){
	if (beginning) head_match = match_run;
	beginning   = false;
	match_run   = 0;
	read_index += cigar_iter->Length;
      }
      else if (cigar_iter->Type == 'D'){
	if (beginning) head_match = match_run;
	beginning  = false;
	match_run  = 0;
	ref_index += cigar_iter->Length;
      }
      else if (cigar_iter->Type == 'S' || cigar_iter->Type == 'H')
	break;
      else {
	string msg = "Invalid CIGAR char";
	msg += cigar_iter->Type;
	printErrorAndDie(msg);
      }
      cigar_iter++;
    }

    // Process trailing clip CIGAR types
    if (cigar_iter != aln.CigarData.end() && cigar_iter->Type == 'S'){
      read_index += cigar_iter->Length;
      cigar_iter++;
    }
    if (cigar_iter != aln.CigarData.end() && cigar_iter->Type == 'H')
      cigar_iter++;
    
    // Ensure that we processed all CIGAR options
    if (cigar_iter != aln.CigarData.end()){
      if (ref_index >= ref_seq.size())
	return pair<int,int>(-1,-1);
      else
	printErrorAndDie("Improperly formatted CIGAR string");
    }
    
    // Ensure that CIGAR string corresponded to aligned bases
    if (read_index != aln.QueryBases.size()){
      if (ref_index >= ref_seq.size())
	return pair<int,int>(-1,-1);
      else
	printErrorAndDie("CIGAR string does not correspond to alignment bases");
    }
    
    if (beginning)
      return pair<int,int>(match_run, match_run);
    else
      return pair<int,int>(head_match, match_run);
  } 


  /* 
     Stores the sequence, start and end position of the read after removing clipped bases
     using the provided references
   */
  void GetUnclippedInfo(BamTools::BamAlignment& aln, string& bases, int32_t& unclipped_start, int32_t& unclipped_end){
    unclipped_start = aln.Position;
    unclipped_end   = aln.Position-1;
    bool begin      = true;
    int start_index = 0, num_bases = 0;
    for(vector<BamTools::CigarOp>::iterator cigar_iter = aln.CigarData.begin(); cigar_iter != aln.CigarData.end(); cigar_iter++){
      switch(cigar_iter->Type) {
      case 'D':
	unclipped_end += cigar_iter->Length;
	begin          = false;
	break;
      case 'H':
	break;
      case 'S':
	if (begin) start_index += cigar_iter->Length;
	break;
      case 'M':
	unclipped_end += cigar_iter->Length;
	num_bases     += cigar_iter->Length;
	begin          = false;
	break;
      case 'I':
	num_bases += cigar_iter->Length;
	begin      = false;
	break;
      default:
	string msg = "Invalid CIGAR char ";
	msg += cigar_iter->Type;
	printErrorAndDie(msg);
	break;
      }
    }
    bases = aln.QueryBases.substr(start_index, num_bases);
  }

 
  bool HasLargestEndMatches(BamTools::BamAlignment& aln, const string& ref_seq, int ref_seq_start, int max_external, int max_internal){
    // Extract sequence, start and end coordinates of read after clipping
    string bases;
    int start, end;
    GetUnclippedInfo(aln, bases, start, end);

    // Check that the prefix match is the longest
    if (start >= ref_seq_start && start < ref_seq_start + static_cast<int>(ref_seq.size())){
      int start_index = start - ref_seq_start;
      int start       = max(0, start_index - max_external);
      int stop        = min(static_cast<int>((ref_seq.size()-1)), start_index + max_internal);
      vector<int> match_counts;
      ZAlgorithm::GetPrefixMatchCounts(bases, ref_seq, start, stop, match_counts);

      int align_index = start_index - start;
      int num_matches = match_counts[align_index];
      for (int i = 0; i < static_cast<int>(match_counts.size()); i++){
	if (i == align_index)
	  continue;
	if (match_counts[i] >= num_matches)
	  return false;
      }
    }

    // Check that the suffix match is the longest
    if (end >= ref_seq_start && end < ref_seq_start + static_cast<int>(ref_seq.size())){
      int end_index = end - ref_seq_start;
      int start     = max(0, end_index - max_internal);
      int stop      = min(static_cast<int>(ref_seq.size()-1), end_index + max_external);
      vector<int> match_counts;
      ZAlgorithm::GetSuffixMatchCounts(bases, ref_seq, start, stop, match_counts);
      
      int align_index = end_index - start;
      int num_matches = match_counts[align_index];
      for (int i = 0; i < static_cast<int>(match_counts.size()); i++){
	if (i == align_index)
	  continue;
	if (match_counts[i] >= num_matches)
	  return false;
      }
    }       
    return true;
  }

  void GetNumClippedBases(BamTools::BamAlignment& aln, int& num_hard_clips, int& num_soft_clips){
    num_hard_clips = 0;
    num_soft_clips = 0;
    for (std::vector<BamTools::CigarOp>::iterator cigar_iter = aln.CigarData.begin(); cigar_iter != aln.CigarData.end(); cigar_iter++){
      switch(cigar_iter->Type){
      case 'H':
	num_hard_clips += cigar_iter->Length;
	break;
      case 'S':
	num_soft_clips += cigar_iter->Length;
	break;
      default:
	break;
      }
    }
  }
}
