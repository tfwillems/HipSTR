#include <algorithm>
#include <assert.h>
#include <climits>
#include <iostream>
#include <string>
#include <vector>

#include "HaplotypeGenerator.h"
#include "RepeatBlock.h"
#include "../stringops.h"

void HaplotypeGenerator::trim(int ideal_min_length, int32_t& region_start, int32_t& region_end, std::vector<std::string>& sequences) const {
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
  max_left_trim  = std::min(LEFT_PAD,  max_left_trim);
  max_right_trim = std::min(RIGHT_PAD, max_right_trim);

  // Don't trim past the padding flanks
  max_left_trim  = std::max(0, std::min(min_len-RIGHT_PAD, max_left_trim));
  max_right_trim = std::max(0, std::min(min_len-LEFT_PAD,  max_right_trim));

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
  region_start += left_trim;
  region_end   -= right_trim;
}

bool HaplotypeGenerator::extract_sequence(const Alignment& aln, int32_t region_start, int32_t region_end, std::string& seq) const {
  if (aln.get_start() >= region_start) return false;
  if (aln.get_stop()  <= region_end)   return false;

  int align_index = 0; // Index into alignment string
  int char_index  = 0; // Index of current operation in current CIGAR element
  int32_t pos     = aln.get_start();
  auto cigar_iter = aln.get_cigar_list().begin();
  
  // Extract region sequence if fully spanned by alignment
  std::stringstream reg_seq;
  while (cigar_iter != aln.get_cigar_list().end()){
    if (char_index == cigar_iter->get_num()){
      cigar_iter++;
      char_index = 0;
    }
    else if (pos == region_end){
      if (cigar_iter->get_type() == 'I'){
	reg_seq << aln.get_alignment().substr(align_index, cigar_iter->get_num());
	align_index += cigar_iter->get_num();
	char_index = 0;
	cigar_iter++;
      }
      else {
	seq = uppercase(reg_seq.str());
	return true;
      }
    }
    else if (pos >= region_start){
      int32_t num_bases = std::min(region_end-pos, cigar_iter->get_num()-char_index);
      switch(cigar_iter->get_type()){	 
      case 'I':
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
    else if (pos > region_end)
      assert(false);
    else {
      int32_t num_bases;
      if (cigar_iter->get_type() == 'I')
	num_bases = cigar_iter->get_num()-char_index;
      else {
	num_bases = std::min(region_start-pos, cigar_iter->get_num()-char_index);
	pos      += num_bases;
      }
      align_index  += num_bases;
      char_index   += num_bases;
    }
  }
  printErrorAndDie("Logical error in extract_sequence");
  return false;
}

void HaplotypeGenerator::select_candidate_alleles(const std::string& ref_seq, int ideal_min_length,
						  const std::vector< std::vector<Alignment> >& alignments, const std::vector<std::string>& vcf_alleles,
						  int32_t& region_start, int32_t& region_end,
						  std::vector<std::string>& alleles, std::vector<bool>& from_vcf) const {
  assert(alleles.empty() && from_vcf.empty());
  std::map<std::string, double> total_frac_counts;
  std::map<std::string, int> total_read_counts, must_inc;
  int total_reads = 0, total_samples = 0;

  // Initialize sequences with reference sequence and VCF alleles
  alleles.push_back(ref_seq);
  from_vcf.push_back(true);
  bool have_ref = false;
  for (unsigned int i = 0; i < vcf_alleles.size(); ++i)
    if (vcf_alleles[i].compare(ref_seq) != 0){
      alleles.push_back(vcf_alleles[i]);
      from_vcf.push_back(true);
    }
    else
      have_ref = true;
  assert(vcf_alleles.empty() || have_ref); // VCF should always contain the reference allele

  // Set of currently selected alleles
  std::set<std::string> allele_set(alleles.begin(), alleles.end());

  // Determine the number of reads and samples supporting each allele
  for (unsigned int i = 0; i < alignments.size(); ++i){
    int sample_reads = 0;
    std::map<std::string, int> counts;
    for (unsigned int j = 0; j < alignments[i].size(); ++j){
      std::string subseq;
      if (extract_sequence(alignments[i][j], region_start, region_end, subseq)){
	total_read_counts[subseq] += 1;
	counts[subseq]      += 1;
	total_reads++;
	sample_reads++;
      }
    }

    // Identify alleles strongly supported by the current sample
    for (auto iter = counts.begin(); iter != counts.end(); ++iter){
      if ((iter->second >= MIN_READS_STRONG_SAMPLE) && (iter->second >= MIN_FRAC_STRONG_SAMPLE*sample_reads))
	must_inc[iter->first] += 1;
      total_frac_counts[iter->first] += iter->second*1.0/sample_reads;
    }

    if (sample_reads > 0)
      total_samples++;
  }
  
  // Add alleles with strong support from a subset of samples
  for (auto iter = must_inc.begin(); iter != must_inc.end(); ++iter){
    if ((iter->second >= MIN_STRONG_SAMPLES) && (allele_set.find(iter->first) == allele_set.end())){
      alleles.push_back(iter->first);
      allele_set.insert(iter->first);
      from_vcf.push_back(false);
    }
  }

  // Add alleles supported by other criteria
  for (auto iter = total_frac_counts.begin(); iter != total_frac_counts.end(); ++iter){
    if ((iter->second > MIN_FRAC_SAMPLES*total_samples) || (total_read_counts[iter->first] > MIN_FRAC_READS*total_reads)){
      if (allele_set.find(iter->first) == allele_set.end()){
	alleles.push_back(iter->first);
	allele_set.insert(iter->first);
	from_vcf.push_back(false);
      }
    }
  }
  
  // Sort regions by length and then by sequence (apart from reference sequence)
  std::sort(alleles.begin()+1, alleles.end(), orderByLengthAndSequence);

  // Clip identical regions
  trim(ideal_min_length, region_start, region_end, alleles);
}

void HaplotypeGenerator::get_aln_bounds(const std::vector< std::vector<Alignment> >& alignments,
					int32_t& min_aln_start, int32_t& max_aln_stop) const {
  // Determine the minimum and maximum alignment boundaries
  min_aln_start = INT_MAX;
  max_aln_stop  = INT_MIN;
  for (auto vec_iter = alignments.begin(); vec_iter != alignments.end(); ++vec_iter){
    for (auto aln_iter = vec_iter->begin(); aln_iter != vec_iter->end(); ++aln_iter){
      min_aln_start = std::min(min_aln_start, aln_iter->get_start());
      max_aln_stop  = std::max(max_aln_stop,  aln_iter->get_stop());
    }
  }
}

bool HaplotypeGenerator::add_vcf_nonrepeat_haplotype_block(int32_t pos, const std::string& chrom_seq,
							   const std::vector<std::string>& vcf_alleles){
  add_vcf_haplotype_block(pos, chrom_seq, vcf_alleles, NULL);
}


bool HaplotypeGenerator::add_vcf_haplotype_block(int32_t pos, const std::string& chrom_seq,
						 const std::vector<std::string>& vcf_alleles, const StutterModel* stutter_model){
  assert(!vcf_alleles.empty());
  if (!failure_msg_.empty())
    printErrorAndDie("Unable to add a VCF haplotype block as a previous addition failed");
  if (finished_)
    printErrorAndDie("Unable to add a VCF haplotype block as the haplotype generation process is already finished");
  
  // Verify that the first VCF allele matches the reference genome
  std::string ref = vcf_alleles.front();
  int32_t ref_len = (int32_t) ref.size();
  assert(ref.compare(uppercase(chrom_seq.substr(pos, ref_len))) == 0);

  // Ensure that the new haplotype block won't overlap with previous blocks
  if (!hap_blocks_.empty() && (pos < hap_blocks_.back()->end())){
    failure_msg_ = "Haplotype blocks are too near to one another";
    return false;
  }

  // Create the new haplotype block, whose type depends on whether the stutter model is NULL
  if (stutter_model == NULL)
    hap_blocks_.push_back(new    HapBlock(pos, pos+ref_len, ref, true));
  else
    hap_blocks_.push_back(new RepeatBlock(pos, pos+ref_len, ref, true, stutter_model->period(), stutter_model));

  // Add each alternative allele to the haplotype block
  for (unsigned int i = 1; i < vcf_alleles.size(); ++i){
    std::string seq = uppercase(vcf_alleles[i]);
    hap_blocks_.back()->add_alternate(seq, false); // VCF-based alleles are designated as not removable
  }

  return true;
}

bool HaplotypeGenerator::add_repeat_haplotype_block(const Region& region, const std::string& chrom_seq, const std::vector< std::vector<Alignment> >& alignments,
						    const StutterModel* stutter_model){
  if (!failure_msg_.empty())
    printErrorAndDie("Unable to add a haplotype block, as a previous addition failed");
  if (finished_)
    printErrorAndDie("Unable to add a VCF haplotype block as the haplotype generation process is already finished");

  // Ensure that we don't exceed the chromosome bounds
  if (region.start() < REF_FLANK_LEN + LEFT_PAD || region.stop() + REF_FLANK_LEN + RIGHT_PAD > chrom_seq.size()){
    failure_msg_ = "Haplotype blocks are too near to the chromosome ends";
    return false;
  }

  // Extract the alignment boundaries
  int32_t min_aln_start, max_aln_stop;
  get_aln_bounds(alignments, min_aln_start, max_aln_stop);

  // Ensure that we have some spanning alignments to work with
  int32_t region_start = region.start() - LEFT_PAD;
  int32_t region_end   = region.stop()  + RIGHT_PAD;
  std::string ref_seq  = uppercase(chrom_seq.substr(region_start, region_end-region_start));
  if (min_aln_start + 5 >= region_start || max_aln_stop - 5 <= region_end){
    failure_msg_ = "No spanning alignments";
    return false;
  }

  // Extract candidate STR sequences (using some padding to ensure indels near STR ends are included)
  std::vector<std::string> vcf_alleles; // Use an empty list as no VCF input here
  std::vector<std::string> alleles;
  std::vector<bool> from_vcf;
  int ideal_min_length = 3*region.period(); // Would ideally have at least 3 repeat units in each allele after trimming
  select_candidate_alleles(ref_seq, ideal_min_length, alignments, vcf_alleles, region_start, region_end, alleles, from_vcf);

  // Ensure that the new haplotype block won't overlap with previous blocks
  if (!hap_blocks_.empty() && (region_start < hap_blocks_.back()->end())){
    failure_msg_ = "Haplotype blocks are too near to one another";
    return false;
  }

  // Add the new haplotype block
  hap_blocks_.push_back(new RepeatBlock(region_start, region_end, alleles.front(), !vcf_alleles.empty(), stutter_model->period(), stutter_model));
  for (unsigned int i = 1; i < alleles.size(); i++)
    hap_blocks_.back()->add_alternate(alleles[i], from_vcf[i]);

  return true;
}

bool HaplotypeGenerator::fuse_haplotype_blocks(const std::string& chrom_seq){
  if (!failure_msg_.empty())
    printErrorAndDie("Unable to fuse haplotype blocks as previous additions failed");
  if (hap_blocks_.empty())
    printErrorAndDie("Unable to fuse haplotype blocks as none have been added");
  if (finished_)
    return true;

  assert(REF_FLANK_LEN > 10);
  assert(hap_blocks_.front()->start() >= REF_FLANK_LEN);
  assert(hap_blocks_.back()->end() + REF_FLANK_LEN <= chrom_seq.size());

  // Determine the boundaries for the fixed length reference-only flanking blocks
  int32_t min_start = std::min(hap_blocks_.front()->start()-10, std::max(hap_blocks_.front()->start() - REF_FLANK_LEN, min_aln_start_));
  int32_t max_stop  = std::max(hap_blocks_.back()->end()+10,    std::min(hap_blocks_.back()->end()    + REF_FLANK_LEN, max_aln_stop_));
  //int32_t min_start = hap_blocks_.front()->start() - REF_FLANK_LEN;
  //int32_t max_stop  = hap_blocks_.back()->end()    + REF_FLANK_LEN;

  // Interleave the existing variant blocks with new reference-only haplotype blocks
  std::vector<HapBlock*> fused_blocks;
  int32_t start = min_start;
  for (int i = 0; i < hap_blocks_.size(); ++i){
    int32_t end = hap_blocks_[i]->start();
    if (start < end)
      fused_blocks.push_back(new HapBlock(start, end, uppercase(chrom_seq.substr(start, end-start)), false));
    fused_blocks.push_back(hap_blocks_[i]);
    start = hap_blocks_[i]->end();
  }
  if (start < max_stop)
    fused_blocks.push_back(new HapBlock(start, max_stop, uppercase(chrom_seq.substr(start, max_stop-start)), false));

  hap_blocks_ = fused_blocks;
  finished_   = true;
  return true;
}
