#include <deque>
#include <vector>

#include "error.h"
#include "filter_bams.h"
#include "interval_tree.h"

#include "bamtools/include/api/BamWriter.h"


std::string trim_alignment_name(BamTools::BamAlignment& aln){
  std::string aln_name = aln.Name;
  if (aln_name.size() > 2){
    if (aln_name[aln_name.size()-2] == '/')
      aln_name.resize(aln_name.size()-2);
  }
  return aln_name;
}


/*
  Filter BAMs in which paired end reads are adjacent to one another in the BAM file.
  Utilizes an interval tree to determine if each pair overlaps a region of interest.
 */
void filter_bam_paired_mode(BamTools::BamReader& reader,
			    std::vector< std::vector<Region> >& regions,
			    std::map<std::string, int>& chrom_order,
			    std::string& output_filename,
			    InsertSizeCounter& counter,
			    bool analyze_insert_size){
  BamTools::RefVector ref_vector = reader.GetReferenceData();

  // Construct interval trees for each chromosome's regions indexed by the reference
  // id in the BAM file + 1
  std::vector< IntervalTree<bool> > interval_trees;
  interval_trees.push_back(IntervalTree<bool>()); // Empty interval tree for RefID=-1 (unmapped)
  std::cerr << "Constructing interval trees" << std::endl;
  for (int i = 0; i < ref_vector.size(); i++){
    std::vector< Interval<bool> > region_intervals;

    if (chrom_order.find(ref_vector[i].RefName) == chrom_order.end())
      interval_trees.push_back(IntervalTree<bool>(region_intervals));
    else {
      std::vector<Region>& chrom_regions = regions[chrom_order.find(ref_vector[i].RefName)->second];
      for (int j = 0; j < chrom_regions.size(); j++)
	region_intervals.push_back(Interval<bool>(chrom_regions[j].start(), chrom_regions[j].stop(), true));
      interval_trees.push_back(IntervalTree<bool>(region_intervals));
    }
  }

  BamTools::BamWriter writer;
  bool file_open = writer.Open(output_filename, reader.GetHeaderText(), ref_vector);
  if (!file_open) printErrorAndDie("Failed to open output BAM file");
  BamTools::BamAlignment alignments[2];
   
  int64_t read_count = 0;
  int aln_count      = 0;
  while (reader.GetNextAlignment(alignments[aln_count])){
    read_count++;
    if (analyze_insert_size)
      counter.process_alignment(alignments[aln_count]);

    if (aln_count == 0)
      aln_count++;
    else {
      std::string key_a = trim_alignment_name(alignments[0]);
      std::string key_b = trim_alignment_name(alignments[1]);

      if (key_a.compare(key_b) == 0){
	// Process set of paired reads
	bool overlaps = interval_trees[alignments[0].RefID+1].overlaps(alignments[0].Position, alignments[0].GetEndPosition()); 
	overlaps     |= interval_trees[alignments[1].RefID+1].overlaps(alignments[1].Position, alignments[1].GetEndPosition());
	if (overlaps){
	  if (!writer.SaveAlignment(alignments[0]))
	    printErrorAndDie("Failed to save alignment");   
	  if (!writer.SaveAlignment(alignments[1]))
	    printErrorAndDie("Failed to save alignment");
	}
	aln_count = 0;
      }
      else {
	// Process single read
	bool overlaps = interval_trees[alignments[0].RefID+1].overlaps(alignments[0].Position, alignments[0].GetEndPosition());
	if (overlaps)
	  if(!writer.SaveAlignment(alignments[0]))
	    printErrorAndDie("Failed to save alignment");
	    
	// Move second read to first position
	alignments[0] = alignments[1];
      }
    }

    if (read_count % 1000000 == 0){
      std::cerr << "\tProcessing read # " << read_count << std::endl;
      if (analyze_insert_size)
	counter.output_summary_statistics(std::cerr);
    }
  }
  writer.Close();
}


void filter_bam(BamTools::BamReader& reader,
                std::vector< std::vector<Region> >&regions,
                std::map<std::string, int>& chrom_order,
                std::string& output_filename,
		InsertSizeCounter& counter,
		bool analyze_insert_size){

  BamTools::RefVector ref_vector = reader.GetReferenceData();
  BamTools::BamWriter writer;
  bool file_open = writer.Open(output_filename, reader.GetHeaderText(), ref_vector);
  if (!file_open) printErrorAndDie("Failed to open output BAM file");
  BamTools::BamAlignment alignment;

  int chrom_id = -2, chrom_idx = -1, region_idx = -1;  // Use -2 for chrom_id b/c *, the reference for unmapped reads, has a RefID of -1 
  int64_t read_count = 0;
  int32_t prev_pos   = -1, max_regions = 0;
  std::set<int> proc_chroms;
  std::map<std::string, int32_t> aligned_mate_pairs;
  std::map<std::string, BamTools::BamAlignment> unaligned_mate_pairs;

  while (reader.GetNextAlignment(alignment)){
    read_count++;
    if (analyze_insert_size)
      counter.process_alignment(alignment);

    if (read_count % 1000000 == 0){
      std::cerr << "\tProcessing read # " << read_count << std::endl;
      if (analyze_insert_size)
	counter.output_summary_statistics(std::cerr);
      
           
      // Clear unaligned reads without identified mate pairs and likely don't have any mapped mates
      auto unaln_iter = unaligned_mate_pairs.begin();
      while (unaln_iter != unaligned_mate_pairs.end()){
	if (prev_pos - unaln_iter->second.Position > 100000)
	  unaligned_mate_pairs.erase(unaln_iter++);
	else
	  unaln_iter++;
      }
    }

    if (alignment.RefID != chrom_id){
      proc_chroms.insert(chrom_id);
      if (proc_chroms.find(alignment.RefID) != proc_chroms.end())
	printErrorAndDie("Chromosomes in BAM file must be in sorted order. Out of order at read " + alignment.Name);
      chrom_id   = alignment.RefID;
      if (alignment.RefID == -1 || chrom_order.find(ref_vector[alignment.RefID].RefName) == chrom_order.end()) {
	chrom_idx   = -1;
	max_regions = 0;
      }
      else {
	chrom_idx  = chrom_order[ref_vector[alignment.RefID].RefName];
	max_regions = regions[chrom_idx].size();
      }

      region_idx = 0;
      prev_pos   = alignment.Position;
      aligned_mate_pairs.clear();
      unaligned_mate_pairs.clear();
    }

    if (alignment.Position < prev_pos)
      printErrorAndDie("BAM files must be sorted. Out of order at read " + alignment.Name);
    else
      prev_pos = alignment.Position;
    
    // Check if it overlaps a region
    while (region_idx < max_regions && regions[chrom_idx][region_idx].stop() < alignment.Position)
      region_idx++;

    std::string aln_key = trim_alignment_name(alignment);

    if (region_idx < max_regions && regions[chrom_idx][region_idx].start() <= alignment.GetEndPosition()){
      // Process region overlap
      if (!writer.SaveAlignment(alignment))  printErrorAndDie("Failed to save alignment for STR-containing read");

      auto iter = unaligned_mate_pairs.find(aln_key);
      if (iter != unaligned_mate_pairs.end()){
	if (!writer.SaveAlignment(iter->second))
	  printErrorAndDie("Failed to save alignment for mate pair read");
	unaligned_mate_pairs.erase(iter);
      }
      else {
	auto iter = aligned_mate_pairs.find(aln_key);
	if (iter != aligned_mate_pairs.end())
	  aligned_mate_pairs.erase(iter);
	else
	  aligned_mate_pairs.insert(std::pair<std::string, int32_t>(aln_key, alignment.Position));
      }
    }
    else {
      // No region overlap but examine any relevant mate pair information
      auto iter = aligned_mate_pairs.find(aln_key);
      if (iter != aligned_mate_pairs.end()){
	aligned_mate_pairs.erase(iter);
	if (!writer.SaveAlignment(alignment))
          printErrorAndDie("Failed to save alignment for mate pair read");
      }
      else {
	auto iter = unaligned_mate_pairs.find(aln_key);
	if (iter != unaligned_mate_pairs.end())
	  unaligned_mate_pairs.erase(iter);
	else 
	  unaligned_mate_pairs.insert(std::pair<std::string, BamTools::BamAlignment>(aln_key, alignment));
      }
    } 
  }
  writer.Close();
}
