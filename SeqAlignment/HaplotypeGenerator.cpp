#include <algorithm>
#include <assert.h>
#include <climits>
#include <iostream>
#include <string>
#include <vector>

#include "HaplotypeGenerator.h"
#include "RepeatBlock.h"
#include "../stringops.h"

const double MIN_FRAC_READS   = 0.05;
const double MIN_FRAC_SAMPLES = 0.05;

// Minimum fraction and minimum number of a sample's reads an allele 
// must be present in for the allele to be strongly supported by the sample
const double MIN_FRAC_STRONG_SAMPLE  = 0.2; 
const double MIN_READS_STRONG_SAMPLE = 2;
const double MIN_STRONG_SAMPLES      = 1;

// Maximum fraction of samples for which a deletion can extend into the STR block
// If more samples have one such deletion, we'll retry to rebuild the STR region
// using a large padding window
const double MAX_FRAC_SAMPLE_DEL_FAIL = 0.01; 

void trim(int32_t left_padding, int32_t right_padding, int ideal_min_length, int32_t& rep_region_start, int32_t& rep_region_end, std::vector<std::string>& sequences){
  int min_len = INT_MAX;
  for (unsigned int i = 0; i < sequences.size(); i++)
    min_len = std::min(min_len, (int)sequences[i].size());
  if (min_len <= ideal_min_length)
    return;

  int max_left_trim = 0, max_right_trim = 0;
  while (max_left_trim < min_len-ideal_min_length){
    unsigned int j = 1;
    while (j < sequences.size()){
      if (sequences[j][max_left_trim] != sequences[j-1][max_left_trim])
	break;
      j++;
    }
    if (j != sequences.size()) 
      break;
    max_left_trim++;
  }
  while (max_right_trim < min_len-ideal_min_length){
    char c = sequences[0][sequences[0].size()-1-max_right_trim];
    unsigned int j = 1;
    while (j < sequences.size()){
      if (sequences[j][sequences[j].size()-1-max_right_trim] != c)
	break;
      j++;
    }
    if (j != sequences.size())
      break;
    max_right_trim++;
  }

  // Don't trim past the flanks
  max_left_trim  = std::min(left_padding,  max_left_trim);
  max_right_trim = std::min(right_padding, max_right_trim); 

  // Don't trim past the padding flanks
  max_left_trim  = std::max(0, std::min(min_len-right_padding, max_left_trim));
  max_right_trim = std::max(0, std::min(min_len-left_padding,  max_right_trim));

  // Determine the left and right trims that clip as much as possible
  // but are as equal in size as possible
  int left_trim, right_trim;
  if (min_len - 2*std::min(max_left_trim, max_right_trim) <= ideal_min_length){
    left_trim = right_trim = std::min(max_left_trim, max_right_trim);
    while (min_len - left_trim - right_trim < ideal_min_length){
      if (left_trim > right_trim)
	left_trim--;
      else
	right_trim--;
    }
  }
  else {
    if (max_left_trim > max_right_trim){
      right_trim = max_right_trim;
      left_trim  = std::min(max_left_trim, min_len-ideal_min_length-max_right_trim);
    }
    else {
      left_trim  = max_left_trim;
      right_trim = std::min(max_right_trim, min_len-ideal_min_length-max_left_trim);
    }
  }

  // Adjust sequences and position
  for (unsigned int i = 0; i < sequences.size(); i++)
    sequences[i] = sequences[i].substr(left_trim, sequences[i].size()-left_trim-right_trim);
  rep_region_start += left_trim;
  rep_region_end   -= right_trim;
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
  printErrorAndDie("Logical error in extract_sequence");
  return false;
}

bool stringLengthLT(const std::string& s1, const std::string& s2){
  return s1.size() < s2.size();
}

void generate_candidate_str_seqs(std::string& ref_seq, std::string& chrom_seq, int32_t left_padding, int32_t right_padding, int ideal_min_length,
				 std::vector< std::vector<Alignment> >& alignments, std::vector<std::string>& vcf_alleles, bool search_bams_for_alleles,
				 int32_t& rep_region_start, int32_t& rep_region_end, std::vector<std::string>& sequences){
  assert(sequences.size() == 0);

  std::map<std::string, double> sample_counts;
  std::map<std::string, int>    read_counts;
  std::map<std::string, int>    must_inc;
  int tot_reads   = 0;
  int tot_samples = 0;
  if (search_bams_for_alleles){
    // Determine the number of reads and number of samples supporting each allele
    for (unsigned int i = 0; i < alignments.size(); i++){
      int samp_reads = 0;
      std::map<std::string, int> counts;
      for (unsigned int j = 0; j < alignments[i].size(); j++){
	std::string subseq;
	if (extract_sequence(alignments[i][j], rep_region_start, rep_region_end, subseq)){
	  read_counts[subseq]   += 1;
	  counts[subseq]        += 1;
	  tot_reads++;
	  samp_reads++;
	}
      } 
      
      // Identify alleles strongly supported by sample
      for (auto iter = counts.begin(); iter != counts.end(); iter++){
	if (iter->second > MIN_READS_STRONG_SAMPLE && iter->second > MIN_FRAC_STRONG_SAMPLE*samp_reads)
	  must_inc[iter->first] += 1;
	sample_counts[iter->first] += iter->second*1.0/samp_reads;
      }
      
      if (samp_reads > 0)
	tot_samples++;
    } 
  }

  // Add VCF alleles to list (apart from reference sequence) and remove from other sets
  int ref_index = -1;
  for (unsigned int i = 0; i < vcf_alleles.size(); i++){
    sequences.push_back(vcf_alleles[i]);
    auto iter_1 = sample_counts.find(vcf_alleles[i]);
    if (iter_1 != sample_counts.end()){
      sample_counts.erase(iter_1);
      read_counts.erase(vcf_alleles[i]);
    }
    auto iter_2 = must_inc.find(vcf_alleles[i]);
    if (iter_2 != must_inc.end())
      must_inc.erase(iter_2);
    if (vcf_alleles[i].compare(ref_seq) == 0)
      ref_index = i;
  }
  
  if (search_bams_for_alleles) {
    // Add alleles with strong support from a subset of samples
    for (auto iter = must_inc.begin(); iter != must_inc.end(); iter++){
      if (iter->second >= MIN_STRONG_SAMPLES){
	auto iter_1 = sample_counts.find(iter->first);
	auto iter_2 = read_counts.find(iter->first);
	std::cerr << "Strong   stats: " << iter->first << " " << iter->second << " " << iter_1->second/alignments.size()
		  << " " << iter_2->second << " " << iter_2->second*1.0/tot_reads << std::endl;
	sample_counts.erase(iter_1);
	read_counts.erase(iter_2);
	sequences.push_back(iter->first);
	if (iter->first.compare(ref_seq) == 0)
	  ref_index = sequences.size()-1;
      }
  }
  
    // Identify additional alleles satisfying thresholds
    for (auto iter = sample_counts.begin(); iter != sample_counts.end(); iter++){
      if (iter->second > MIN_FRAC_SAMPLES*tot_samples || read_counts[iter->first] > MIN_FRAC_READS*tot_reads){
	std::cerr << "Sequence stats: " << iter->first << " " << iter->first.size() << " " 
		  << iter->second*1.0/tot_samples << " " << read_counts[iter->first] << " " << read_counts[iter->first]*1.0/tot_reads << std::endl;
	sequences.push_back(iter->first);
	if (ref_index == -1 && (iter->first.compare(ref_seq) == 0))
	  ref_index = sequences.size()-1;
      }
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

  // Clip identical regions
  trim(left_padding, right_padding, ideal_min_length, rep_region_start, rep_region_end, sequences); 
}


int check_deletion_bounds(std::vector< std::vector<Alignment> >& alignments, int32_t start, int32_t end,
			  int32_t& min_del_start, int32_t& max_del_stop, std::vector<bool>& call_sample){
  assert(call_sample.size() == 0);
  int sample_fail_count = 0; // Number of samples with 1 or more non-enclosed deletions
  min_del_start = (start+end)/2;
  max_del_stop  = (start+end)/2;

  // Extract deletion coordinates
  std::vector<int32_t> del_starts, del_ends;
  for (auto vec_iter = alignments.begin(); vec_iter != alignments.end(); vec_iter++){
    bool sample_fail = false;
    for (auto aln_iter = vec_iter->begin(); aln_iter != vec_iter->end(); aln_iter++)
      aln_iter->get_deletion_boundaries(del_starts, del_ends);

    // Check that they all are contained within the region or outside the region
    assert(del_starts.size() == del_ends.size());
    for (unsigned i = 0; i < del_starts.size(); i++){
      bool fail = false;
      if (del_starts[i] <= start && del_ends[i] >= start){
	fail = true;
	min_del_start = std::min(min_del_start, del_starts[i]);
      }
      if (del_starts[i] <= end && del_ends[i] >= end){
	fail = true;
	max_del_stop = std::max(max_del_stop, del_ends[i]);
      }      
      sample_fail |= fail;
    }
    del_starts.clear(); del_ends.clear();	
    if (sample_fail){
      sample_fail_count++;
      call_sample.push_back(false);
    }
    else
      call_sample.push_back(true);
  }
  return sample_fail_count;
}

Haplotype* generate_haplotype(Region& str_region, int32_t max_ref_flank_len, std::string& chrom_seq,
			      std::vector< std::vector<Alignment> >& alignments, std::vector<std::string>& vcf_alleles,
			      StutterModel* stutter_model, bool search_bams_for_alleles,
			      std::vector<HapBlock*>& blocks, std::vector<bool>& call_sample){
  assert(blocks.size() == 0);
  assert(call_sample.size() == 0);

  std::map<std::string, int> seqs;

  // Determine the minimum and maximum alignment boundaries
  int32_t min_start = INT_MAX, max_stop = INT_MIN;
  for (auto vec_iter = alignments.begin(); vec_iter != alignments.end(); vec_iter++){
    for (auto aln_iter = vec_iter->begin(); aln_iter != vec_iter->end(); aln_iter++){
      min_start = std::min(min_start, aln_iter->get_start());
      max_stop  = std::max(max_stop,  aln_iter->get_stop());
      seqs[aln_iter->get_sequence()]++;
    }
  }
  
  std::cerr << "# Seqs = " << seqs.size() << std::endl;
			 
  // Extract candidate STR sequences (use some padding to ensure indels near STR ends are included)
  std::vector<std::string> str_seqs;
  int32_t left_padding = 5, right_padding = 5;
  int32_t rep_region_start = str_region.start() < left_padding ? 0 : str_region.start()-left_padding;
  int32_t rep_region_end   = str_region.stop() + right_padding;
  std::string ref_seq      = uppercase(chrom_seq.substr(rep_region_start, rep_region_end-rep_region_start));
  int ideal_min_length     = 3*str_region.period(); // Would ideally have at least 3 repeat units in each allele after trimming

  // TO DO: Use frequency of deletions not contained within window to 
  //   i) Identify problematic regions
  //  ii) Retry with increased window padding?
  int32_t min_del_start, max_del_stop;
  int sample_fail_count = check_deletion_bounds(alignments, rep_region_start, rep_region_end, min_del_start, max_del_stop, call_sample);
  
  // Recheck deletions after expanding the STR window
  if(1.0*sample_fail_count/alignments.size() <= MAX_FRAC_SAMPLE_DEL_FAIL){
    std::cerr << "PASS SAMPLE DEL STATS: " << rep_region_start << "\t" << rep_region_end << "\t" << sample_fail_count 
	      << "\t" << min_del_start << "\t" << max_del_stop << "\t" << std::endl;
  }
  else {
    std::cerr << "FAIL SAMPLE DEL STATS: " << rep_region_start << "\t" << rep_region_end << "\t" << sample_fail_count 
	      << "\t" << min_del_start << "\t" << max_del_stop << "\t" << std::endl;

    //if ((rep_region_start-min_del_start <= 10) && (max_del_stop-rep_region_end <= 10)){
      call_sample.clear();
      rep_region_start  = std::min(rep_region_start, min_del_start-5);
      rep_region_end    = std::max(rep_region_end,   max_del_stop+5);
      left_padding      = str_region.start() - rep_region_start;
      right_padding     = rep_region_end     - str_region.stop();
      ref_seq           = uppercase(chrom_seq.substr(rep_region_start, rep_region_end-rep_region_start));
      sample_fail_count = check_deletion_bounds(alignments, rep_region_start, rep_region_end, min_del_start, max_del_stop, call_sample);
      if(1.0*sample_fail_count/alignments.size() <= MAX_FRAC_SAMPLE_DEL_FAIL){
	std::cerr << "REDONE_PASS SAMPLE DEL STATS: " << rep_region_start << "\t" << rep_region_end << "\t" << sample_fail_count
		  << "\t" << min_del_start << "\t" << max_del_stop << "\t" << std::endl;
      }
      else {
	std::cerr << "REDONE_FAIL SAMPLE DEL STATS: " << rep_region_start << "\t" << rep_region_end << "\t" << sample_fail_count
		  << "\t" << min_del_start << "\t" << max_del_stop << "\t" << std::endl;
      }
      //}
  }
  
  // Trim boundaries so that the reference flanks aren't too long
  if (rep_region_start > max_ref_flank_len)
    min_start = std::max(rep_region_start-max_ref_flank_len, min_start);
  max_stop = std::min(rep_region_end+max_ref_flank_len, max_stop);

  // Extend each VCF allele by padding size
  std::vector<std::string> ext_vcf_alleles;
  if (vcf_alleles.size() != 0){
    std::string lflank = uppercase(chrom_seq.substr(rep_region_start,  str_region.start()-rep_region_start));
    std::string rflank = uppercase(chrom_seq.substr(str_region.stop(), rep_region_end-str_region.stop())); 
    for (unsigned int i = 0; i < vcf_alleles.size(); i++){
      ext_vcf_alleles.push_back(lflank + vcf_alleles[i] + rflank);
      std::cerr << ext_vcf_alleles.back() << std::endl;
    }
    std::cerr << ext_vcf_alleles[0] << std::endl << ref_seq << std::endl;
    assert(ext_vcf_alleles[0].compare(ref_seq) == 0);
  }
  generate_candidate_str_seqs(ref_seq, chrom_seq, left_padding, right_padding, ideal_min_length, alignments, 
			      ext_vcf_alleles, search_bams_for_alleles, rep_region_start, rep_region_end, str_seqs);
  
  // Create a set of haplotype regions, consisting of STR sequence block flanked by two reference sequence stretches
  assert(rep_region_start > min_start && rep_region_end < max_stop);

  if(str_seqs[0].compare(uppercase(chrom_seq.substr(rep_region_start, rep_region_end-rep_region_start))) != 0){
    std::cerr << str_seqs[0] << std::endl
	      << uppercase(chrom_seq.substr(rep_region_start, rep_region_end-rep_region_start)) << std::endl;
    assert(str_seqs[0].compare(uppercase(chrom_seq.substr(rep_region_start, rep_region_end-rep_region_start))) == 0);
  }
  blocks.clear();
  blocks.push_back(new HapBlock(min_start, rep_region_start, uppercase(chrom_seq.substr(min_start, rep_region_start-min_start))));    // Ref sequence preceding STRS
  blocks.push_back(new RepeatBlock(rep_region_start, rep_region_end, 
				   uppercase(chrom_seq.substr(rep_region_start, rep_region_end-rep_region_start)), str_region.period(), stutter_model));
  blocks.push_back(new HapBlock(rep_region_end, max_stop, uppercase(chrom_seq.substr(rep_region_end, max_stop-rep_region_end))));  // Ref sequence following STRs
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


Haplotype* generate_haplotype(std::vector<std::string>& vcf_alleles, StutterModel* stutter_model, std::vector<HapBlock*>& blocks){
  // TO DO: Implement function designed for cases in which we're reading alleles and their priors from a VCF
  printErrorAndDie("Variant of generate_haplotype not implemented");
  return NULL;
}
