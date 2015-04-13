#include <algorithm>
#include <climits>
#include <map>
#include <string>
#include <sstream>
#include <vector>

#include "constants.h"
#include "error.h"

#include "AlignmentData.h"
#include "AlignmentOps.h"
#include "BaseQuality.h"
#include "HapBlock.h"
#include "Haplotype.h"
#include "HaplotypeAligner.h"
#include "HaplotypeCreator.h"
#include "HTMLCreator.h"
#include "NWNoRefEndPenalty.h"
#include "RepeatRegion.h"
#include "StringOps.h"




void getMaxInsertionSizes(std::vector<Alignment>& alignments, 
			  std::map<int32_t,int>& max_insertions){
  max_insertions.clear();
  for (std::vector<Alignment>::iterator align_iter = alignments.begin(); align_iter != alignments.end(); align_iter++){
    int32_t position = align_iter->get_start();
    const std::vector<CigarElement>& cigar_ops = align_iter->get_cigar_list();
    std::map<int32_t, int>::iterator insertion_iter;
    for(std::vector<CigarElement>::const_iterator cigar_iter = cigar_ops.begin(); cigar_iter != cigar_ops.end(); cigar_iter++){
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


bool strcomp(std::string s1, std::string s2){
  return s1.size() > s2.size();
}

void overlayHaplotypeRegions(std::vector<HapBlock*> blocks,
			     std::map<int32_t,int>& max_insertions,
			     int32_t min_start,
			     std::vector<std::string>& results){
  int max_num_seqs = 0;
  for (unsigned int i = 0; i < blocks.size(); i++)
    max_num_seqs = std::max(max_num_seqs, blocks[i]->num_options());

  std::vector<int> max_reg_sizes;
  std::map<int32_t,int>::iterator insertion_iter = max_insertions.begin();
  for (unsigned int i = 0; i < blocks.size(); i++){
    // Skip reference-only haplotype blocks
    if (blocks[i]->num_options() == 1){
      max_reg_sizes.push_back(-1);
      continue;
    }

    while (insertion_iter->first < blocks[i]->start())
      insertion_iter++;
    int size = blocks[i]->end()-blocks[i]->start();

    while(insertion_iter->first <= blocks[i]->end()){
      size += insertion_iter->second;
      insertion_iter++;
    }
    max_reg_sizes.push_back(size);
  }

  // Create N rows, where N is the maximum number of sequences across all regions
  for (int n = 0; n < max_num_seqs; n++){
    std::stringstream result;
    int32_t position = min_start;
    int block_index  = 0;
    std::map<int32_t,int>::iterator insertion_iter = max_insertions.begin();
    for (std::vector<HapBlock*>::iterator block_iter = blocks.begin(); block_iter != blocks.end(); block_iter++){
      // Skip haplotype blocks that only have 1 potential sequence (stretch of reference sequence)
      if ((*block_iter)->num_options() == 1){
	block_index++;
	continue;
      }

      // Spaces before/between regions
      for(int32_t i = position; i < (*block_iter)->start(); i++){
	if (insertion_iter != max_insertions.end() && (i == insertion_iter->first)){
	  for (int j = 0; j < insertion_iter->second; j++)
	    result << SPACE_CHAR;
	  insertion_iter++;
	}
	result << SPACE_CHAR;
      }
      
      // Write region sequence
      if (n < blocks[block_index]->num_options())
	result << blocks[block_index]->get_seq(n);
      
      // Right pad region with spaces for non-maximal length sequences
      if (n < blocks[block_index]->num_options())
	for (int i = 0; i < max_reg_sizes[block_index]-blocks[block_index]->get_seq(n).size();i++)
	  result << NOT_APP_CHAR;
      else
	for(int i = 0; i < max_reg_sizes[block_index]; i++)
	  result << SPACE_CHAR;
      
      // Discard insertions within region
      if ((*block_iter)->start() == (*block_iter)->end())
	while(insertion_iter != max_insertions.end() && insertion_iter->first <= (*block_iter)->end())
	  insertion_iter++;
      else
	while (insertion_iter != max_insertions.end() && insertion_iter->first < (*block_iter)->end())
	  insertion_iter++;
      
      // Update variables for next iteration
      position = (*block_iter)->end();
      block_index++;
    }
    results.push_back(result.str());
  }
}

void overlayAlignments(std::vector<Alignment>& alignments, 
		       std::map<int32_t,int>& max_insertions, 
		       std::vector<std::string>& results, 
		       int32_t& min_start, int32_t& max_stop){						
  min_start = INT_MAX;
  max_stop  = INT_MIN;
  for (std::vector<Alignment>::iterator align_iter = alignments.begin(); align_iter != alignments.end(); align_iter++)
    min_start = std::min(min_start, align_iter->get_start());

  results.clear();
  max_insertions.clear();
  getMaxInsertionSizes(alignments, max_insertions);
  max_insertions.insert(std::pair<int32_t,int>(INT_MAX, -1)); // Dummy insertion so that there's always an insertion after the end of the reads
  
  for (std::vector<Alignment>::iterator align_iter = alignments.begin(); align_iter != alignments.end(); align_iter++){
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
    for(int32_t i = min_start; i < position; i++){
      if (i == insertion_iter->first){
	for (int j = 0; j < insertion_iter->second; j++)
	  result << SPACE_CHAR;
	insertion_iter++;
      }
      result <<  SPACE_CHAR;
    }
    
    while (true){
      if (length == 0){
	cigar_iter++;
	if (cigar_iter == cigar_ops.end())
	  break;
	type   = cigar_iter->get_type();
	length = cigar_iter->get_num();
      }
      
      if (position == insertion_iter->first){
	int num_ins = insertion_iter->second;
	if (type == 'I')
	  num_ins -= length;
	for(int i = 0; i < num_ins; i++)
	  result << NOT_APP_CHAR;
	insertion_iter++;
      }
      
      if (type == 'M' || type == '=' || type == 'X'){
	int num_match = std::min(length, insertion_iter->first - position);
	for(int i = 0; i < num_match; i++){
	  result << nucleotides[nuc_index++];
	  position++;
	  length--;
	}
      }
      else if (type == 'I'){
	for(int i = 0; i < length; i++)
	  result << nucleotides[nuc_index++];
	length = 0;
      }
      else if (type == 'D'){
	int num_del = std::min(length, insertion_iter->first - position);
	for(int i = 0; i < num_del; i++){
	  result << DELETION_CHAR;
	  position++;
	  length--;
	}
      }
      else if (type == 'H'){
	length = 0;
	continue;
	printErrorAndDie("Cigar string H option not implemeted.");
      }
      else if (type == 'S'){
	nuc_index += length;
	length = 0;
	continue;
	printErrorAndDie("Cigar string S option not implemeted.");
      }
      else
	printErrorAndDie("Invalid cigar option encountered.");
    }
    max_stop = std::max(max_stop, position-1);
    results.push_back(result.str());
  }
}


std::string arrangeReferenceString(std::string& ref_sequence, 
				   std::map<int32_t,int>& max_insertions,
				   std::string& locus_id,
				   int32_t str_start,
				   int32_t str_stop,
				   int32_t min_start,
				   int32_t max_stop,
				   bool draw_locus_id,
				   std::ostream& output){
  // Hack to deal with lobSTR VCF off by 1
  int offset = 1;

  std::stringstream ref_result;
  std::vector<bool> within_locus;
  std::map<int32_t,int>::iterator insertion_iter = max_insertions.begin();
  for(int32_t i = min_start; i <= max_stop; i++){
    if (i == insertion_iter->first){
      for (int j = 0; j < insertion_iter->second; j++){
	ref_result << NOT_APP_CHAR;
	within_locus.push_back(false);
      }
      insertion_iter++;
    }
    ref_result << ref_sequence[i];
    within_locus.push_back((i>=str_start-offset && i <= str_stop-offset));
  } 
  
  std::string ref_alignment = ref_result.str();
  writeReferenceString(ref_alignment, output, locus_id, within_locus, draw_locus_id);
  return ref_alignment;
}

void get_alignment_bounds(std::vector<Alignment>& alignments, 
			  int32_t& min_start, 
			  int32_t& max_stop){
  min_start = INT_MAX; max_stop = INT_MIN;  
  for (std::vector<Alignment>::const_iterator iter = alignments.begin(); iter != alignments.end(); iter++){
    min_start = std::min(min_start, iter->get_start());
    max_stop  = std::max(max_stop,  iter->get_stop());
  }
}

void arrangeAlignments(std::vector<BamTools::BamAlignment>& alignments, 
		       std::string& ref_sequence, 
		       std::string locus_id, 
		       std::ostream& output, 
		       bool draw_locus_id,
		       vcf::VariantCallFile* vcf_data, 
		       std::map<std::string, std::string>& sample_info,
		       BaseQuality& base_quality) {
  // Generate realigned reads
  std::vector<Alignment> new_alignments;
  for (std::vector<BamTools::BamAlignment>::iterator iter = alignments.begin(); iter != alignments.end(); iter++){
    Alignment new_alignment;
    realign(*iter, ref_sequence, new_alignment);
    new_alignments.push_back(new_alignment);
  }

  // Sort alignments by (sample name, start position)
  sortAlignments(new_alignments);

  // Determine STR coordinates
  int32_t str_start, str_stop;
  bool got_tag = GetIntBamTag(alignments[0], START_TAG, &str_start);
  if (!got_tag) printErrorAndDie("Failed to retrieve start tag for STR locus");
  got_tag = GetIntBamTag(alignments[0], STOP_TAG, &str_stop);
  if (!got_tag) printErrorAndDie("Failed to retrieve stop tag for STR locus");

  // Minimum and maximum coordinates of alignments
  int32_t min_coord, max_coord;
  get_alignment_bounds(new_alignments, min_coord, max_coord);

  // Calculate repetitive regions in reference sequence
  std::string region_seq = ref_sequence.substr(min_coord, max_coord-min_coord+1);
  std::vector<RepeatRegion> repeat_regions;
  get_repeat_regions(region_seq, min_coord, repeat_regions);
#ifdef DEBUG  
  std::cerr << "Repeat Regions:" << std::endl;
  for (unsigned int i = 0; i < repeat_regions.size(); i++)
    repeat_regions[i].print(std::cerr);
#endif

  std::vector<std::string> hap_results;
  std::vector< std::pair<int32_t, int32_t> > regions;
  std::vector< std::vector<std::string> > region_seqs;
  getHaplotypeRegions(new_alignments, regions);
  extractRegionSequences(new_alignments, regions, ref_sequence, region_seqs);

  // Discard regions with only 1 candidate sequence (i.e. the reference sequence)
  int cur_index = 0, ins_index = 0;
  while (cur_index < regions.size()){
    if (region_seqs[cur_index].size() > 1){
      regions[ins_index]     = regions[cur_index];
      region_seqs[ins_index] = region_seqs[cur_index];
      ins_index++;
    }  
    cur_index++;
  }
  regions.resize(ins_index);
  region_seqs.resize(ins_index);

  // TO DO: Identify stutter regions by allele size pattern and merge them with repetitive regions
  /*
  int reg_index = 0, rep_index = 0;
  while (reg_index < regions.size() && rep_index < repeat_regions.size()){
    
  }
  */

  if (regions.size() != 0 && (regions.front().first < min_coord || regions.back().second > max_coord))
    printErrorAndDie("Regions not fully contained within aligment window bounds");
  regions.push_back(std::pair<int32_t,int32_t>(max_coord, max_coord+1)); // Dummy region

  // Create an integrated set of haplotype regions, consisting of alternating blocks of 
  // reference sequence and alternate sequences
  std::vector<HapBlock*> hap_blocks;
  int region_index = 0;
  int32_t coord    = min_coord;
  while (coord < max_coord){
    if (coord < regions[region_index].first){
      hap_blocks.push_back(new HapBlock(coord, regions[region_index].first, uppercase(ref_sequence.substr(coord, regions[region_index].first-coord))));
      coord = regions[region_index].first;
    }
    else {
      HapBlock* block = new HapBlock(regions[region_index].first, regions[region_index].second, region_seqs[region_index][0]);
      for (int j=1; j < region_seqs[region_index].size(); j++)
	block->add_alternate(region_seqs[region_index][j]);
      hap_blocks.push_back(block);
      coord = regions[region_index].second;
      region_index++;
    }
  }

  // Order alternate sequences by length
  for (unsigned int i = 0; i < hap_blocks.size(); i++)
    hap_blocks[i]->order_alternates_by_length();

  // Output HTML rendering of the reference sequence, haplotype blocks, and alignments 
  if (true){
    // Arrange alignments
    std::vector<std::string> align_results;
    std::map<int32_t,int> max_insertions;
    int32_t min_start, max_stop;
    overlayAlignments(new_alignments, max_insertions, align_results, min_start, max_stop);
    std::vector<std::string> sample_names;
    for (std::vector<Alignment>::iterator iter = new_alignments.begin(); iter != new_alignments.end(); iter++)
      sample_names.push_back(iter->get_sample());

    // Arrange reference sequence
    std::string ref_alignment = arrangeReferenceString(ref_sequence, max_insertions, locus_id, 
						       str_start, str_stop, min_start, max_stop,
						       draw_locus_id, output);
    // Arrange haplotype blocks
    std::vector<std::string> hap_results;
    overlayHaplotypeRegions(hap_blocks, max_insertions, min_start, hap_results);
    std::vector<std::string> hap_samples;
    for (int i = 0; i < hap_results.size(); i++)
      hap_samples.push_back("Blocks");
    
    // Write to HTML
    output << "<div class='alignments'>"      << "\n"
	   << "\t<table class=\"readtable\">" << "\n";
    writeAlignmentStrings(ref_alignment, output, locus_id, hap_results,   hap_samples,  NULL,     sample_info, false);
    writeAlignmentStrings(ref_alignment, output, locus_id, align_results, sample_names, vcf_data, sample_info, true);
    output << "\t</table>" << "\n"
	   << "<br>"       << "\n"
	   << "</div>"     << std::endl;
  }

  // Create Haplotype from constituent blocks
  Haplotype hap(hap_blocks);
#ifdef DEBUG
  std::cerr << "Haplotype block structure:" << std::endl;
  hap.print_block_structure(20, 45, std::cerr);
#endif

  // Align all reads across all samples
  // HaplotypeAligner hap_aligner(&hap);
  // hap_aligner.align(new_alignments, base_quality);

  // Delete haplotype block instances
  for (int i = 0; i < hap_blocks.size(); i++)
    delete hap_blocks[i];
  hap_blocks.clear();
}
