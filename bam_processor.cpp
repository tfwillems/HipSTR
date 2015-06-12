#include <fstream>
#include <iostream>
#include <locale>
#include <sstream>
#include <stdlib.h>

#include "bam_processor.h"
#include "alignment_filters.h"
#include "error.h"
#include "seqio.h"
#include "stringops.h"

std::string BamProcessor::parse_lobstr_rg(BamTools::BamAlignment& aln){
  std::string rg;
  std::string rg_tag = "RG";
  char tag_type = 'Z';
  if (!aln.GetTagType(rg_tag, tag_type))
    printErrorAndDie("Failed to retrieve BAM alignment's RG tag");
  aln.GetTag("RG", rg);
  std::vector<std::string> rg_tokens;
  split_by_delim(rg, ';', rg_tokens);
  if (rg_tokens.size() != 3 || rg_tokens[0].compare("lobSTR") != 0){
    printErrorAndDie("Improperly formatted lobSTR RG tag: " + rg);
  }
  return rg_tokens[1];
}

std::string BamProcessor::trim_alignment_name(BamTools::BamAlignment& aln){
  std::string aln_name = aln.Name;
  if (aln_name.size() > 2){
    if (aln_name[aln_name.size()-2] == '/')
      aln_name.resize(aln_name.size()-2);
  }
  return aln_name;
}

std::string get_str_ref_allele(uint32_t start, uint32_t end, std::string& chrom_seq){
  std::locale loc;
  std::string seq = chrom_seq.substr(start-1, end-start+1);
  return uppercase(seq);
}

void BamProcessor::read_and_filter_reads(BamTools::BamMultiReader& reader, std::string& chrom_seq, 
					 std::vector<Region>::iterator region_iter, std::map<std::string, std::string>& file_read_groups,
					 std::vector<std::string>& rg_names, 
					 std::vector< std::vector<BamTools::BamAlignment> >& paired_strs_by_rg,
					 std::vector< std::vector<BamTools::BamAlignment> >& mate_pairs_by_rg,
					 std::vector< std::vector<BamTools::BamAlignment> >& unpaired_strs_by_rg,
					 BamTools::BamWriter& bam_writer){
  std::vector<BamTools::BamAlignment> region_alignments;
  int read_count = 0;
  int diff_chrom_mate = 0, unmapped_mate = 0, not_spanning = 0; // Counts for filters that are always applied
  int insert_size = 0, multimapped = 0, flank_len = 0, bp_before_indel = 0, end_match_window = 0, num_end_matches = 0; // Counts for filters that are user-controlled
  BamTools::BamAlignment alignment;

  std::vector<BamTools::BamAlignment> paired_str_alns, mate_alns, unpaired_str_alns;
  std::map<std::string, BamTools::BamAlignment> potential_strs, potential_mates;

  while (reader.GetNextAlignment(alignment)){
    read_count++;
    
    if (check_mate_info_){
      // Ignore read if its mate pair chromosome doesn't match
      if (alignment.RefID != alignment.MateRefID){
	diff_chrom_mate++;
	continue;
      }
      // Ignore read if its mate pair is unmapped
      if (alignment.InsertSize == 0){
	unmapped_mate++;
	continue;
      }
    }

    // Ignore read if its mate pair distance exceeds the threshold
    if (abs(alignment.InsertSize) > MAX_MATE_DIST){
      insert_size++;
      continue;
    }
    // Ignore read if multimapper and filter specified
    if (REMOVE_MULTIMAPPERS && alignment.HasTag("XA")){
      multimapped++;
      continue;
    }

    bool pass = true;

    // Simple test to exclude mate pairs
    if (alignment.Position > region_iter->stop() || alignment.GetEndPosition() < region_iter->start()){
      pass = false;
    }
    // Ignore read if it does not span the STR
    if (pass && REQUIRE_SPANNING && (alignment.Position > region_iter->start() || alignment.GetEndPosition() < region_iter->stop())){
      not_spanning++;
      pass = false;
    }
    // Ignore read if it has insufficient flanking bases on either side of the STR
    if (pass && (alignment.Position > (region_iter->start()-MIN_FLANK) || alignment.GetEndPosition() < (region_iter->stop()+MIN_FLANK))){
      flank_len++;
      pass = false;
    }
    // Ignore read if there is an indel within the first MIN_BP_BEFORE_INDEL bps from each end
    if (pass && MIN_BP_BEFORE_INDEL > 0){
      std::pair<int, int> num_bps = AlignmentFilters::GetEndDistToIndel(alignment);
      if ((num_bps.first != -1 && num_bps.first < MIN_BP_BEFORE_INDEL) || (num_bps.second != -1 && num_bps.second < MIN_BP_BEFORE_INDEL)){
	bp_before_indel++;
	pass = false;
      }
    }
    // Ignore read if there is another location within MAXIMAL_END_MATCH_WINDOW bp for which it has a longer end match
    if (pass && MAXIMAL_END_MATCH_WINDOW > 0){
      bool maximum_end_matches = AlignmentFilters::HasLargestEndMatches(alignment, chrom_seq, 0, MAXIMAL_END_MATCH_WINDOW, MAXIMAL_END_MATCH_WINDOW);
      if (!maximum_end_matches){
	end_match_window++;
	pass = false;
      }
    }
    // Ignore read if it doesn't match perfectly for at least MIN_READ_END_MATCH bases on each end
    if (pass && MIN_READ_END_MATCH > 0){
      std::pair<int,int> match_lens = AlignmentFilters::GetNumEndMatches(alignment, chrom_seq, 0);
      if (match_lens.first < MIN_READ_END_MATCH || match_lens.second < MIN_READ_END_MATCH){ 
	num_end_matches++;
	pass = false;
      }
    }

    if (pass)
      region_alignments.push_back(alignment);

    std::string aln_key = trim_alignment_name(alignment);
    if (pass){
      auto aln_iter = potential_mates.find(aln_key);
      if (aln_iter != potential_mates.end()){
	paired_str_alns.push_back(alignment);
	mate_alns.push_back(aln_iter->second);
	potential_mates.erase(aln_iter);
      }
      else 
	potential_strs.insert(std::pair<std::string, BamTools::BamAlignment>(aln_key, alignment));
    }
    else {
      auto aln_iter = potential_strs.find(aln_key);
      if (aln_iter != potential_strs.end()){
	paired_str_alns.push_back(aln_iter->second);
	mate_alns.push_back(alignment);
	potential_strs.erase(aln_iter);
      }
      else 
	potential_mates.insert(std::pair<std::string, BamTools::BamAlignment>(aln_key, alignment));
    }
  }

  for (auto aln_iter = potential_strs.begin(); aln_iter != potential_strs.end(); ++aln_iter)
    unpaired_str_alns.push_back(aln_iter->second);
  potential_strs.clear();
  potential_mates.clear();
  std::cerr << "Found " << paired_str_alns.size() << " fully paired reads and " << unpaired_str_alns.size() << " unpaired reads" << std::endl;
  
  std::cerr << read_count << " reads overlapped region, of which " 
	    << "\n\t" << diff_chrom_mate  << " had mates on a different chromosome"
	    << "\n\t" << unmapped_mate    << " had unmapped mates"
	    << "\n\t" << not_spanning     << " did not span the STR"
	    << "\n\t" << insert_size      << " failed the insert size filter"      
	    << "\n\t" << multimapped      << " were removed due to multimapping"
	    << "\n\t" << flank_len        << " had too bps in one or more flanks"
	    << "\n\t" << bp_before_indel  << " had too few bp before the first indel"
	    << "\n\t" << end_match_window << " did not have the maximal number of end matches within the specified window"
	    << "\n\t" << num_end_matches  << " had too few bp matches along the ends"
  	    << "\n"   << region_alignments.size() << " PASSED ALL FILTERS" << "\n" << std::endl; 
    
  // Output the spanning reads to a BAM file, if requested
  if (bam_writer.IsOpen()){
    for (auto read_iter = region_alignments.begin(); read_iter != region_alignments.end(); read_iter++){
      // Add RG to BAM record based on file
      if (!use_lobstr_rg_){
	std::string rg_tag = "lobSTR;" + file_read_groups[read_iter->Filename] + ";" + file_read_groups[read_iter->Filename];
	read_iter->AddTag("RG", "Z", rg_tag);
      }
      
      // Add STR start and stop tags
      if (read_iter->HasTag("XS"))
	read_iter->RemoveTag("XS");
      if(!read_iter->AddTag("XS",  "i", region_iter->start()))
	printErrorAndDie("Failed to modify XS tag");
      if (read_iter->HasTag("XE"))
	read_iter->RemoveTag("XE");
      if(!read_iter->AddTag("XE", "i", region_iter->stop()))
	printErrorAndDie("Failed to modify XE tag");
      if (!bam_writer.SaveAlignment(*read_iter))  
	printErrorAndDie("Failed to save alignment for STR-spanning read");
    }
  }

  // Separate the reads based on their associated read groups
  std::map<std::string, int> rg_indices;
  for (unsigned int i = 0; i < paired_str_alns.size(); ++i){
    std::string rg = use_lobstr_rg_ ? parse_lobstr_rg(paired_str_alns[i]) : file_read_groups[paired_str_alns[i].Filename];
    int rg_index;
    auto index_iter = rg_indices.find(rg);
    if (index_iter == rg_indices.end()){
      rg_index = rg_indices.size();
      rg_indices[rg] = rg_index;
      rg_names.push_back(rg);
      paired_strs_by_rg.push_back(std::vector<BamTools::BamAlignment>());
      unpaired_strs_by_rg.push_back(std::vector<BamTools::BamAlignment>());
      mate_pairs_by_rg.push_back(std::vector<BamTools::BamAlignment>());
    }
    else
      rg_index = index_iter->second;

    // Record STR read and its mate pair
    paired_strs_by_rg[rg_index].push_back(paired_str_alns[i]);
    mate_pairs_by_rg[rg_index].push_back(mate_alns[i]);
  }
  for (unsigned int i = 0; i < unpaired_str_alns.size(); ++i){
    std::string rg = use_lobstr_rg_ ? parse_lobstr_rg(unpaired_str_alns[i]) : file_read_groups[unpaired_str_alns[i].Filename];
    int rg_index;
    auto index_iter = rg_indices.find(rg);
    if (index_iter == rg_indices.end()){
      rg_index = rg_indices.size();
      rg_indices[rg] = rg_index;
      rg_names.push_back(rg);
      paired_strs_by_rg.push_back(std::vector<BamTools::BamAlignment>());
      unpaired_strs_by_rg.push_back(std::vector<BamTools::BamAlignment>());
      mate_pairs_by_rg.push_back(std::vector<BamTools::BamAlignment>());
    }
    else
      rg_index = index_iter->second;

    // Record unpaired STR read
    unpaired_strs_by_rg[rg_index].push_back(unpaired_str_alns[i]);
  }
}

void BamProcessor::process_regions(BamTools::BamMultiReader& reader, 
				   std::string& region_file, std::string& fasta_dir,
				   std::map<std::string, std::string>& file_read_groups,
				   BamTools::BamWriter& bam_writer, std::ostream& out, int32_t max_regions, std::string chrom){
  std::vector<Region> regions;
  readRegions(region_file, regions, max_regions, chrom);
  orderRegions(regions);
  
  std::string ref_seq;
  BamTools::RefVector ref_vector = reader.GetReferenceData();
  int cur_chrom_id = -1; std::string chrom_seq;
  for (auto region_iter = regions.begin(); region_iter != regions.end(); region_iter++){
    std::cerr << "Processing region " << region_iter->chrom() << " " << region_iter->start() << " " << region_iter->stop() << std::endl;
    int chrom_id = reader.GetReferenceID(region_iter->chrom());
    
    // Read FASTA sequence for chromosome 
    if (cur_chrom_id != chrom_id){
      cur_chrom_id      = chrom_id;
      std::string chrom = ref_vector[cur_chrom_id].RefName;
      std::cerr << "Reading fasta file for " << chrom << std::endl;
      readFasta(chrom+".fa", fasta_dir, chrom_seq);
    }
    
    if(!reader.SetRegion(chrom_id, (region_iter->start() < MAX_MATE_DIST ? 0: region_iter->start()-MAX_MATE_DIST), 
			 chrom_id, region_iter->stop() + MAX_MATE_DIST)){
      printErrorAndDie("One or more BAM files failed to set the region properly");
    }

    std::vector<std::string> rg_names;
    std::vector< std::vector<BamTools::BamAlignment> > paired_strs_by_rg, mate_pairs_by_rg, unpaired_strs_by_rg;
    read_and_filter_reads(reader, chrom_seq, region_iter, file_read_groups, rg_names, paired_strs_by_rg, mate_pairs_by_rg, unpaired_strs_by_rg, bam_writer);
    std::string ref_allele = get_str_ref_allele(region_iter->start(), region_iter->stop(), chrom_seq);
    process_reads(paired_strs_by_rg, mate_pairs_by_rg, unpaired_strs_by_rg, rg_names, *region_iter, ref_allele, chrom_seq, out);
  }
}
