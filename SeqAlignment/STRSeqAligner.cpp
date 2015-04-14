#include <algorithm>
#include <climits>
#include <map>
#include <string>
#include <sstream>
#include <vector>

#include "constants.h"
#include "../error.h"

#include "AlignmentData.h"
#include "AlignmentOps.h"
#include "../base_quality.h"
#include "HapBlock.h"
#include "Haplotype.h"
#include "HaplotypeAligner.h"
#include "HaplotypeCreator.h"
#include "HTMLCreator.h"
#include "NWNoRefEndPenalty.h"
#include "RepeatRegion.h"
#include "../stringops.h"

#include "HaplotypeGenerator.h"


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

void get_alignment_bounds(std::vector<Alignment>& alignments, 
			  int32_t& min_start, 
			  int32_t& max_stop){
  min_start = INT_MAX; max_stop = INT_MIN;  
  for (std::vector<Alignment>::const_iterator iter = alignments.begin(); iter != alignments.end(); iter++){
    min_start = std::min(min_start, iter->get_start());
    max_stop  = std::max(max_stop,  iter->get_stop());
  }
}

void arrangeAlignments(RepeatRegion& rep_region, 
		       std::string& ref_sequence,
		       std::vector<BamTools::BamAlignment>& alignments, 
		       std::string locus_id, bool draw_locus_id,
		       std::ostream& output, 
		       std::map<std::string, std::string>& sample_info,
		       BaseQuality& base_quality,
		       vcf::VariantCallFile* vcf_data) {
  // Generate realigned reads
  std::cerr << "Realigning reads" << std::endl;
  std::vector<Alignment> new_alignments;
  for (std::vector<BamTools::BamAlignment>::iterator iter = alignments.begin(); iter != alignments.end(); iter++){
    Alignment new_alignment;
    realign(*iter, ref_sequence, new_alignment);
    new_alignments.push_back(new_alignment);
  }

  // Sort alignments by (sample name, start position)
  std::cerr << "Sorting alignments" << std::endl;
  sortAlignments(new_alignments);
  
  // Minimum and maximum coordinates of alignments
  int32_t min_coord, max_coord;
  get_alignment_bounds(new_alignments, min_coord, max_coord);

  std::string region_seq = ref_sequence.substr(min_coord, max_coord-min_coord+1);

  // Generate candidate haplotype window regions
  std::vector<std::string> hap_results;
  std::vector< std::pair<int32_t, int32_t> > regions;
  std::vector< std::vector<std::string> > region_seqs;
  std::cerr << "Extracting haplotype regions" << std::endl;
  getHaplotypeRegions(new_alignments, regions);

  // Merge regions that overlap the repetitive tract
  mergeRegions(rep_region.get_start(), rep_region.get_stop(), regions);


  // Extract sequences for each remaining region
  std::cerr << "Extracting region sequences" << std::endl;
  extractRegionSequences(new_alignments, regions, ref_sequence, region_seqs);

  // Discard regions with only 1 candidate sequence (i.e. the reference sequence)
  // or regions in which all sequences have the same length
  int cur_index = 0, ins_index = 0;
  while (cur_index < regions.size()){
    bool keep_region = region_seqs[cur_index].size() > 1;
    if (keep_region){
      int i = 1;
      for (; i < region_seqs[cur_index].size(); i++){
	if (region_seqs[cur_index][i].size() != region_seqs[cur_index][i-1].size())
	  break;
      }
      keep_region &= (i != region_seqs[cur_index].size());
    }

    if (keep_region){
      regions[ins_index]     = regions[cur_index];
      region_seqs[ins_index] = region_seqs[cur_index];
      ins_index++;
    }  
    cur_index++;
  }
  regions.resize(ins_index);
  region_seqs.resize(ins_index);

  if (regions.size() != 0 && (regions.front().first < min_coord || regions.back().second > max_coord))
    printErrorAndDie("Regions not fully contained within aligment window bounds");
  regions.push_back(std::pair<int32_t,int32_t>(max_coord, max_coord+1)); // Dummy region

  // Create an integrated set of haplotype regions, consisting of alternating blocks of 
  // reference sequence and alternate sequences
  std::cerr << "Constructing haplotype blocks" << std::endl;
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

  std::cerr << "Constructing haplotype" << std::endl;
  Haplotype hap(hap_blocks);
  hap.print_block_structure(30, 100, std::cerr);


  // Align all reads across all samples
  // HaplotypeAligner hap_aligner(&hap);
  // hap_aligner.align(new_alignments, base_quality);



  // Delete haplotype block instances
  for (int i = 0; i < hap_blocks.size(); i++)
    delete hap_blocks[i];
  hap_blocks.clear();

  
  return;
}
