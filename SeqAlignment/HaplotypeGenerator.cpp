#include <algorithm>
#include <assert.h>
#include <climits>
#include <iostream>
#include <string>
#include <vector>

#include "HaplotypeGenerator.h"
#include "RepeatBlock.h"
#include "../stringops.h"

const double MIN_FRAC_READS   = 0.01;
const double MIN_FRAC_SAMPLES = 0.01;

void trim_from_left(int32_t& rep_region_start, int32_t& rep_region_end, std::vector<std::string>& sequences){
  int min_len = INT_MAX;
  for (unsigned int i = 0; i < sequences.size(); i++)
    min_len = std::min(min_len, (int)sequences[i].size());
 
  int left_trim = 0;
  while (left_trim < min_len){
    unsigned int j = 1;
    while (j < sequences.size()){
      if (sequences[j][left_trim] != sequences[j-1][left_trim])
	break;
      j++;
    }
    if (j != sequences.size()) 
      break;
    left_trim++;
  }

  // Adjust sequences and position
  for (unsigned int i = 0; i < sequences.size(); i++)
    sequences[i] = sequences[i].substr(left_trim);
  rep_region_start += left_trim;
}


void trim_from_right(int32_t& rep_region_start, int32_t& rep_region_end, std::vector<std::string>& sequences){
  int min_len = INT_MAX;
  for (unsigned int i = 0; i < sequences.size(); i++)
    min_len = std::min(min_len, (int)sequences[i].size());

  int right_trim = 0;
  while (right_trim < min_len){
    char c = sequences[0][sequences[0].size()-1-right_trim];
    unsigned int j = 1;
    while (j < sequences.size()){
      if (sequences[j][sequences[j].size()-1-right_trim] != c)
	break;
      j++;
    }
    if (j != sequences.size())
      break;
    right_trim++;
  }
  for (unsigned int i = 0; i < sequences.size(); i++)
    sequences[i] = sequences[i].substr(0, sequences[i].size()-right_trim);
  rep_region_end -= right_trim;
}


bool extract_sequence(Alignment& aln, int32_t start, int32_t end, std::string& seq){    
  if (aln.get_start() >= start) return false;
  if (aln.get_stop()  <= end)   return false;

  int align_index = 0; // Index into alignment string
  int char_index  = 0; // Index of current base in current CIGAR element
  int32_t pos     = aln.get_start();
  auto cigar_iter = aln.get_cigar_list().begin();
  
  // Extract region sequence if fully spanned by alignment
  std::stringstream reg_seq;
  while (cigar_iter != aln.get_cigar_list().end()){
    if (char_index == cigar_iter->get_num()){
      cigar_iter++;
      char_index = 0;
    }
    else if (pos > end){
      if (reg_seq.str() == "")
	reg_seq << "X";
      seq = uppercase(reg_seq.str());
      return true;
    }
    else if (pos == end){
      if (cigar_iter->get_type() == 'I'){
	reg_seq << aln.get_alignment().substr(align_index, cigar_iter->get_num());
	align_index += cigar_iter->get_num();
	char_index = 0;
	cigar_iter++;
      }
      else {
	if (reg_seq.str() == "")
	  reg_seq << "X";
	seq = uppercase(reg_seq.str());
	return true;
      }
    }
    else if (pos >= start){
      int32_t num_bases = std::min(end-pos, cigar_iter->get_num()-char_index);
      switch(cigar_iter->get_type()){	 
      case 'I':
	// Insertion within region
	num_bases = cigar_iter->get_num();
	reg_seq << aln.get_alignment().substr(align_index, num_bases);
	break;
      case '=': case 'X':
	reg_seq << aln.get_alignment().substr(align_index, num_bases);
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
	num_bases = std::min(start-pos, cigar_iter->get_num()-char_index);
	pos      += num_bases;
      }
      align_index  += num_bases;
      char_index   += num_bases;
    }
  }
}


bool stringLengthLT(const std::string& s1, const std::string& s2){
  return s1.size() < s2.size();
}

void generate_candidate_str_seqs(std::string& ref_seq,
				 std::vector< std::vector<Alignment> >& paired_strs_by_rg,
				 std::vector< std::vector<Alignment> >& unpaired_strs_by_rg,
				 int32_t& rep_region_start, int32_t& rep_region_end, std::vector<std::string>& sequences){
  assert(paired_strs_by_rg.size() == unpaired_strs_by_rg.size());
  std::map<std::string, double> sample_counts;
  std::map<std::string, int>    read_counts;

  // Determine the number of reads and number of samples supporting each allele
  int num_reads = 0;
  for (unsigned int i = 0; i < paired_strs_by_rg.size(); i++){
    int num_reads = paired_strs_by_rg[i].size() + unpaired_strs_by_rg[i].size();
    for (unsigned int j = 0; j < 2; j++){
      std::vector<Alignment>& reads = (j == 0 ? paired_strs_by_rg[i] : unpaired_strs_by_rg[i]);
      for (unsigned int k = 0; k < reads.size(); k++){
	std::string subseq;
	if (extract_sequence(reads[k], rep_region_start, rep_region_end, subseq)){
	  sample_counts[subseq] += 1.0/num_reads;
	  read_counts[subseq]   += 1;
	  num_reads++;
	}
      } 
    }   
  }

  // Identify reads satisfying thresholds
  int num_samples = paired_strs_by_rg.size();
  int ref_index   = -1;
  sequences.clear();
  for (auto iter = sample_counts.begin(); iter != sample_counts.end(); iter++){
    if (iter->second > MIN_FRAC_SAMPLES*num_samples && read_counts[iter->first] > MIN_FRAC_READS*num_reads){
      sequences.push_back(iter->first);
      if (ref_index == -1 && (iter->first.compare(ref_seq) == 0))
	ref_index = sequences.size()-1;
    }
  }

  // Arrange reference sequence as first element
  if (ref_index == -1)
    sequences.insert(sequences.begin(), ref_seq);
  else {
    sequences[ref_index] = sequences[0];
    sequences[0]         = ref_seq;
  }

  // Sort regions by length (apart from reference sequence)
  std::sort(sequences.begin()+1, sequences.end(), stringLengthLT);

  // Trim to remove identical sequences
  trim_from_left(rep_region_start, rep_region_end, sequences);
  trim_from_right(rep_region_start, rep_region_end, sequences);
}

Haplotype* generate_haplotype(Region& str_region, std::string& chrom_seq,
			      std::vector< std::vector<Alignment> >& paired_strs_by_rg,
			      std::vector< std::vector<Alignment> >& unpaired_strs_by_rg,
			      StutterModel* stutter_model,
			      std::vector<HapBlock*>& blocks){
  // Determine the minimum and maximum alignment boundaries
  int32_t min_start = INT_MAX, max_stop = INT_MIN;
  for (auto vec_iter = paired_strs_by_rg.begin(); vec_iter != paired_strs_by_rg.end(); vec_iter++){
    for (auto aln_iter = vec_iter->begin(); aln_iter != vec_iter->end(); aln_iter++){
      min_start = std::min(min_start, aln_iter->get_start());
      max_stop  = std::max(max_stop,  aln_iter->get_stop());
    }
  }
  for (auto vec_iter = unpaired_strs_by_rg.begin(); vec_iter != unpaired_strs_by_rg.end(); vec_iter++){
    for (auto aln_iter = vec_iter->begin(); aln_iter != vec_iter->end(); aln_iter++){
      min_start = std::min(min_start, aln_iter->get_start());
      max_stop  = std::max(max_stop,  aln_iter->get_stop());
    }
  }
  
  // Extract candidate STR sequences  (use some padding to ensure indels near STR ends are included)
  std::vector<std::string> str_seqs;
  int32_t rep_region_start = str_region.start() < 5 ? 0 : str_region.start()-5;
  int32_t rep_region_end   = str_region.stop() + 5;
  std::string ref_seq      = uppercase(chrom_seq.substr(rep_region_start, rep_region_end-rep_region_start));
  generate_candidate_str_seqs(ref_seq, paired_strs_by_rg, unpaired_strs_by_rg, rep_region_start, rep_region_end, str_seqs);
  
  // Create a set of haplotype regions, consisting of STR sequence block flanked by two reference sequence stretches
  assert(rep_region_start > min_start && rep_region_end < max_stop);
  assert(str_seqs[0].compare(uppercase(chrom_seq.substr(rep_region_start, rep_region_end-rep_region_start))) == 0);
  blocks.clear();
  blocks.push_back(new HapBlock(min_start,        rep_region_start, uppercase(chrom_seq.substr(min_start, rep_region_start-min_start))));    // Ref sequence preceding STRS
  blocks.push_back(new RepeatBlock(rep_region_start, rep_region_end, 
				   uppercase(chrom_seq.substr(rep_region_start, rep_region_end-rep_region_start)), str_region.period(), stutter_model));
  blocks.push_back(new HapBlock(rep_region_end,   max_stop,         uppercase(chrom_seq.substr(rep_region_end, max_stop-rep_region_end))));  // Ref sequence following STRs
  for (unsigned int j = 1; j < str_seqs.size(); j++)
    blocks[1]->add_alternate(str_seqs[j]);

  // Initialize each block's data structures, namely the homopolymer length information
  for (unsigned int i = 0; i < blocks.size(); i++)
    blocks[i]->initialize();
 
  std::cerr << "Constructing haplotype" << std::endl;
  Haplotype* haplotype = new Haplotype(blocks);
  haplotype->print_block_structure(30, 100, std::cerr);  
  return haplotype;
}

