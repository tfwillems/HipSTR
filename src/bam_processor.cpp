#include <fstream>
#include <iostream>
#include <locale>
#include <sstream>
#include <stdlib.h>
#include <time.h>

#include "bam_processor.h"
#include "adapter_trimmer.h"
#include "alignment_filters.h"
#include "error.h"
#include "fasta_reader.h"
#include "pcr_duplicates.h"
#include "stringops.h"
#include "SeqAlignment/AlignmentOps.h"

const std::string ALT_MAP_TAG           = "XA";
const std::string PRIMARY_ALN_SCORE_TAG = "AS";
const std::string SUBOPT_ALN_SCORE_TAG  = "XS";

void BamProcessor::add_passes_filters_tag(BamAlignment& aln, const std::string& passes){
  if (aln.HasTag("PF"))
    if (!aln.RemoveTag("PF"))
      printErrorAndDie("Failed to remove existing passes filters tag from BAM alignment");
  if (!aln.AddStringTag("PF", passes))
    printErrorAndDie("Failed to add passes filters tag to BAM alignment");
}

void BamProcessor::passes_filters(BamAlignment& aln, std::vector<bool>& region_passes){
  assert(region_passes.empty());
  std::string passes;
  if (!aln.GetStringTag("PF", passes))
    printErrorAndDie("Failed to extract passes filters tag from BAM alignment");
  for (int i = 0; i < passes.size(); i++)
    region_passes.push_back(passes[i] == '1' ? true : false);
}

void BamProcessor::write_passing_alignment(BamAlignment& aln, BamWriter* writer){
  if (writer == NULL)
    return;
  if (!writer->SaveAlignment(aln))
    printErrorAndDie("Failed to save alignment");
}

void BamProcessor::write_filtered_alignment(BamAlignment& aln, std::string filter, BamWriter* writer){
  if (writer == NULL)
    return;

  if (aln.HasTag("FT"))
    if (!aln.RemoveTag("FT"))
      printErrorAndDie("Failed to remove alignment's FT tag");

  if (!aln.AddStringTag("FT", filter))
    printErrorAndDie("Failed to add filter tag to alignment");

  if (!writer->SaveAlignment(aln))
    printErrorAndDie("Failed to save alignment");
}

void BamProcessor::extract_mappings(BamAlignment& aln,
				    std::vector< std::pair<std::string, int32_t> >& chrom_pos_pairs) const {
  assert(chrom_pos_pairs.size() == 0);
  if (aln.Ref().compare("*") == 0 || aln.CigarData().size() == 0)
    return;
  chrom_pos_pairs.push_back(std::pair<std::string, int32_t>(aln.Ref(), aln.Position()));
  std::string aln_cigar_string = "";

  for (unsigned int i = 0; i < 2; i++){
    std::string tag = (i == 0 ? "XA" : "SA");
    if (!aln.HasTag(tag.c_str()))
      continue;
    std::string alt_info;
    if (!aln.GetStringTag(tag.c_str(), alt_info))
      printErrorAndDie("Failed to extract XA or SA tag from BAM alignment");
    std::vector<std::string> alts;
    split_by_delim(alt_info, ';', alts);
    for (unsigned int j = 0; j < alts.size(); j++){
      std::vector<std::string> tokens;
      split_by_delim(alts[j], ',', tokens);
      int32_t pos = std::abs(std::stol(tokens[1]));
      if (tokens[0].compare(chrom_pos_pairs[0].first) != 0 || std::abs(pos - chrom_pos_pairs[0].second) > 200){
	// Appropriately handle alt contigs in GRCh38
	// Don't count alternate mappings if they i) are from an alt contig that matches the alignment's chromosome and ii) have the same CIGAR string
	if (i == 0 && string_ends_with(tokens[0], "_alt") && string_starts_with(tokens[0], chrom_pos_pairs[0].first + "_")){
	  if (aln_cigar_string.empty())
	    aln_cigar_string = BuildCigarString(aln.CigarData());
	  if (tokens[2].compare(aln_cigar_string) == 0)
	    continue;
	}

	chrom_pos_pairs.push_back(std::pair<std::string, int32_t>(tokens[0], pos));
      }
    }
  }
}

void BamProcessor::get_valid_pairings(BamAlignment& aln_1, BamAlignment& aln_2,
				      std::vector< std::pair<std::string, int32_t> >& p1, std::vector< std::pair<std::string, int32_t> >& p2) const {
  assert(p1.size() == 0 && p2.size() == 0);
  if (aln_1.Ref().compare("*") == 0 || aln_2.Ref().compare("*") == 0)
    return;

  // BWA-MEM sometimes doesn't report the alternate alignment tag if there are too many alternate mappings
  // To avoid including mismapped reads, we need to be careful about drawing any conclusions
  // regarding alternate mappings. We do this by examining the alignment score of the read or its mate pair relative
  // to the suboptimal alignment score.
  //  i) If the mate pair has no XA tag, we require that it has a decent alignment score relative to its best alternate score.
  //     Otherwise, the mate mapping information is not informative and we discard the pair of reads.
  // ii) If the mate pair has an XA tag but the read has no XA tag, we require that the read has a decent alignment score relative to its best alternate score.
  //     Otherwise, even with an informative mate, there may be other equally valid alignments for the read
  if (!aln_2.HasTag(ALT_MAP_TAG.c_str())){
    if (aln_2.HasTag(PRIMARY_ALN_SCORE_TAG.c_str()) && aln_2.HasTag(SUBOPT_ALN_SCORE_TAG.c_str())){
      int64_t mate_primary_score, mate_subopt_score;
      if(!aln_2.GetIntTag(PRIMARY_ALN_SCORE_TAG.c_str(), mate_primary_score)) printErrorAndDie("Failed to extract the primary alignment score from the BAM record");
      if(!aln_2.GetIntTag(SUBOPT_ALN_SCORE_TAG.c_str(),  mate_subopt_score))  printErrorAndDie("Failed to extract the suboptimal alignment score from the BAM record");
      if (mate_primary_score - mate_subopt_score < 10)
	return;
    }
  }
  else if (!aln_1.HasTag(ALT_MAP_TAG.c_str())){
    if (aln_1.HasTag(PRIMARY_ALN_SCORE_TAG.c_str()) && aln_1.HasTag(SUBOPT_ALN_SCORE_TAG.c_str())){
      int64_t primary_score, subopt_score;
      if(!aln_1.GetIntTag(PRIMARY_ALN_SCORE_TAG.c_str(), primary_score))      printErrorAndDie("Failed to extract the primary alignment score from the BAM record");
      if(!aln_1.GetIntTag(SUBOPT_ALN_SCORE_TAG.c_str(),  subopt_score))       printErrorAndDie("Failed to extract the suboptimal alignment score from the BAM record");
      if (primary_score - subopt_score < 10)
	return;
    }
  }

  std::vector< std::pair<std::string, int32_t> > pairs_1, pairs_2;
  extract_mappings(aln_1, pairs_1);
  extract_mappings(aln_2, pairs_2);
  std::sort(pairs_1.begin(), pairs_1.end());
  std::sort(pairs_2.begin(), pairs_2.end());

  unsigned int min_j = 0;
  for (unsigned int i = 0; i < pairs_1.size(); i++){
    for (unsigned int j = min_j; j < pairs_2.size(); j++){
      int chrom_comp = pairs_1[i].first.compare(pairs_2[j].first);
      if (chrom_comp < 0)
	break;
      else if (chrom_comp > 0)
	min_j = j+1;
      else {
	if (abs(pairs_1[i].second - pairs_2[j].second) < MAX_MATE_DIST){
	  p1.push_back(pairs_1[i]);
	  p2.push_back(pairs_2[j]);
	}
      }
    }
  }
}

std::string BamProcessor::get_read_group(const BamAlignment& aln, const std::map<std::string, std::string>& read_group_mapping) const {
  std::string rg;
  if (!aln.GetStringTag("RG", rg))
    printErrorAndDie("Failed to retrieve BAM alignment's RG tag");
  auto iter = read_group_mapping.find(aln.Filename() + rg);
  if (iter == read_group_mapping.end())
    printErrorAndDie("No sample found for read group " + rg + " in BAM file headers");
  return iter->second;
}

std::string BamProcessor::trim_alignment_name(const BamAlignment& aln) const {
  std::string aln_name = aln.Name();
  if (aln_name.size() > 2){
    if (aln_name[aln_name.size()-2] == '/')
      aln_name.resize(aln_name.size()-2);
  }
  return aln_name;
}

void BamProcessor::read_and_filter_reads(BamCramMultiReader& reader, const std::string& chrom_seq, const RegionGroup& region_group,
					 const std::map<std::string, std::string>& rg_to_sample, std::vector<std::string>& rg_names,
					 std::vector<BamAlnList>& paired_strs_by_rg, std::vector<BamAlnList>& mate_pairs_by_rg, std::vector<BamAlnList>& unpaired_strs_by_rg,
					 BamWriter* pass_writer, BamWriter* filt_writer){
  locus_read_filter_time_ = clock();
  assert(reader.get_merge_type() == BamCramMultiReader::ORDER_ALNS_BY_FILE);

  int32_t read_count = 0, not_spanning = 0, unique_mapping = 0, read_has_N = 0, hard_clip = 0, low_qual_score = 0, num_filt_unpaired_reads = 0;
  BamAlignment alignment;
  BamAlnList paired_str_alns, mate_alns, unpaired_str_alns;
  std::map<std::string, BamAlignment> potential_strs, potential_mates;
  TOO_MANY_READS = false;

  const std::vector<Region>& regions = region_group.regions();
  std::string prev_file  = "";
  int32_t file_index     = 0;
  std::string file_label = "0_";

  while (reader.GetNextAlignment(alignment)){
    // Discard reads where the 1st/2nd mate info isn't clear
    if (alignment.IsPaired() && (!alignment.IsFirstMate() && !alignment.IsSecondMate()))
      continue;

    // Discard reads that don't overlap the STR region and whose mate pair has no chance of overlapping the region
    if (alignment.Position() > region_group.stop() || alignment.GetEndPosition() < region_group.start()){
      if (!alignment.IsPaired() || alignment.MatePosition() == alignment.Position())
	continue;
      if (alignment.MatePosition() > region_group.stop())
	continue;
      if (alignment.MatePosition()+alignment.Length()+100 < region_group.start())
	continue;
    }

    // Stop parsing reads if we've already exceeded the maximum number for downstream analyses
    if (paired_str_alns.size() > MAX_TOTAL_READS){
      TOO_MANY_READS = true;
      break;
    }

    if (!alignment.IsMapped() || alignment.Position() == 0 || alignment.CigarData().size() == 0 || alignment.Length() == 0)
	continue;
    assert(alignment.CigarData().size() > 0 && alignment.Ref().compare("*") != 0);

    // If requested, trim any reads that potentially overlap the STR regions
    if (alignment.Position() < region_group.stop() && alignment.GetEndPosition() >= region_group.start()){
      if (BASE_QUAL_TRIM > ' '){
	// Read trimming doesn't work if hard clipping has been applied
	if (alignment.StartsWithHardClip() || alignment.EndsWithHardClip()){
	  read_count++;
	  hard_clip++;
	  write_filtered_alignment(alignment, "HARD_CLIPPED", filt_writer);
	  continue;
	}

	int32_t length = alignment.Length();
	alignment.TrimLowQualityEnds(BASE_QUAL_TRIM);
	if (alignment.Position() < region_group.stop() && alignment.GetEndPosition() >= region_group.start())
	  if ((alignment.Length() == 0) || (alignment.Length() < length/2))
	    continue;
      }

      // Apply adapter trimming
      adapter_trimmer_.trim_adapters(alignment);

      if (alignment.CigarData().size() == 0 || alignment.Length() == 0)
	continue;
    }

    // Clear out mate alignment cache if we've switched to a new file to reduce memory usage
    // and update the file label for later use
    if (prev_file.compare(alignment.Filename()) != 0){
      prev_file = alignment.Filename();
      potential_mates.clear();

      std::stringstream ss;
      ss << ++file_index << "_";
      file_label = ss.str();
    }

    // Only apply tests to putative STR reads that overlap the STR region
    if (alignment.Position() < region_group.stop() && alignment.GetEndPosition() >= region_group.start()){
      bool pass_one = false; // Denotes if read passed first set of simpler filters
      std::string pass_two(regions.size(), '0'); // Denotes if read passed second set of additional filters for each region
                                                 // Meant to signify if reads that pass first set should be used to generate haplotypes
      std::string filter = "";
      read_count++;

      // Ignore reads with N bases
      if (alignment.QueryBases().find('N') != std::string::npos){
	read_has_N++;
	filter.append("HAS_N_BASES");
      }
      // Ignore reads with a very low overall base quality score
      // Want to avoid situations in which it's more advantageous to have misalignments b/c the scores are so low
      else if (base_quality_.sum_log_prob_correct(alignment.Qualities()) < MIN_SUM_QUAL_LOG_PROB){
	low_qual_score++;
	filter.append("LOW_BASE_QUALS");
      }
      else
	pass_one = true;

      if (pass_one){
	// Determine whether we can use the read for haplotype generation for each region in the group
	int region_index = 0;
	for (auto region_iter = regions.begin(); region_iter != regions.end(); ++region_iter, ++region_index){
	  if ((MIN_FLANK > 0) && (alignment.Position() > (region_iter->start()-MIN_FLANK) || alignment.GetEndPosition() < (region_iter->stop()+MIN_FLANK)))
	    continue;

	  // Ignore read if there is another location within MAXIMAL_END_MATCH_WINDOW bp for which it has a longer end match
	  if (MAXIMAL_END_MATCH_WINDOW > 0){
	    bool maximum_end_matches = AlignmentFilters::HasLargestEndMatches(alignment, chrom_seq, 0, MAXIMAL_END_MATCH_WINDOW, MAXIMAL_END_MATCH_WINDOW);
	    if (!maximum_end_matches){
	      pass_two = std::string(regions.size(), '0');
	      break;
	    }
	  }
	  // Ignore read if it doesn't match perfectly for at least MIN_READ_END_MATCH bases on each end
	  if (MIN_READ_END_MATCH > 0){
	    std::pair<int,int> match_lens = AlignmentFilters::GetNumEndMatches(alignment, chrom_seq, 0);
	    if (match_lens.first < MIN_READ_END_MATCH || match_lens.second < MIN_READ_END_MATCH){
	      pass_two = std::string(regions.size(), '0');
	      break;
	    }
	  }
	  // Ignore read if there is an indel within the first MIN_BP_BEFORE_INDEL bps from each end
	  if (MIN_BP_BEFORE_INDEL > 0){
	    std::pair<int, int> num_bps = AlignmentFilters::GetEndDistToIndel(alignment);
	    if ((num_bps.first != -1 && num_bps.first < MIN_BP_BEFORE_INDEL) || (num_bps.second != -1 && num_bps.second < MIN_BP_BEFORE_INDEL)){
	      pass_two = std::string(regions.size(), '0');
	      break;
	    }
	  }
	  pass_two[region_index] = '1';
	}
      }

      std::string aln_key = file_label + trim_alignment_name(alignment);
      if (pass_one){
	add_passes_filters_tag(alignment, pass_two);
	auto aln_iter = potential_mates.find(aln_key);
	if (aln_iter != potential_mates.end()){
	  if (alignment.IsFirstMate() == aln_iter->second.IsFirstMate()){
	    potential_mates.erase(aln_iter);
	    potential_strs.insert(std::pair<std::string, BamAlignment>(aln_key, alignment));
	    continue;
	  }

	  std::vector< std::pair<std::string, int32_t> > p_1, p_2;
	  get_valid_pairings(alignment, aln_iter->second, p_1, p_2);
	  if (p_1.size() == 1 && p_1[0].second == alignment.Position()){
	    paired_str_alns.push_back(alignment);
	    mate_alns.push_back(aln_iter->second);
	    write_passing_alignment(alignment, pass_writer);
	    write_passing_alignment(aln_iter->second, pass_writer);
	  }
	  else {
	    unique_mapping++;
	    filter.append("NO_UNIQUE_MAPPING");
	    write_filtered_alignment(alignment, filter, filt_writer);
	  }
	  potential_mates.erase(aln_iter);
	}
	else {
	  // Check if read's mate pair also overlaps the STR
	  auto str_iter = potential_strs.find(aln_key);
	  if (str_iter != potential_strs.end()){
	    if (alignment.IsFirstMate() == str_iter->second.IsFirstMate()){
	      read_count--;
	      continue;
	    }

	    std::vector< std::pair<std::string, int32_t> > p_1, p_2;
	    get_valid_pairings(alignment, str_iter->second, p_1, p_2);
	    if (p_1.size() == 1 && p_1[0].second == alignment.Position()){
	      paired_str_alns.push_back(alignment);
	      mate_alns.push_back(str_iter->second);
	      write_passing_alignment(alignment, pass_writer);

	      paired_str_alns.push_back(str_iter->second);
	      mate_alns.push_back(alignment);
	      write_passing_alignment(str_iter->second, pass_writer);
	    }
	    else {
	      unique_mapping += 2;
	      std::string filter = "NO_UNIQUE_MAPPING";
	      write_filtered_alignment(alignment, filter, filt_writer);
	      write_filtered_alignment(str_iter->second, filter, filt_writer);
	    }
	    potential_strs.erase(str_iter);
	  }
	  else
	    potential_strs.insert(std::pair<std::string, BamAlignment>(aln_key, alignment));
	}
      }
      else {
	assert(!filter.empty());
	write_filtered_alignment(alignment, filter, filt_writer);
	potential_mates.insert(std::pair<std::string, BamAlignment>(aln_key, alignment));
      }
    }
    else {
      std::string aln_key = file_label + trim_alignment_name(alignment);
      auto aln_iter = potential_strs.find(aln_key);
      if (aln_iter != potential_strs.end()){
	if (alignment.IsFirstMate() == aln_iter->second.IsFirstMate())
	  continue;

	std::vector< std::pair<std::string, int32_t> > p_1, p_2;
	get_valid_pairings(aln_iter->second, alignment, p_1, p_2);
	if (p_1.size() == 1 && p_1[0].second == aln_iter->second.Position()){
	  paired_str_alns.push_back(aln_iter->second);
	  mate_alns.push_back(alignment);
	  write_passing_alignment(aln_iter->second, pass_writer);
	  write_passing_alignment(alignment, pass_writer);
	}
	else {
	  unique_mapping++;
	  std::string filter = "NO_UNIQUE_MAPPING";
	  write_filtered_alignment(aln_iter->second, filter, filt_writer);
	}
	potential_strs.erase(aln_iter);
      }
      else {
	auto other_iter = potential_mates.find(aln_key);
	if (other_iter != potential_mates.end()){
	  if (alignment.IsFirstMate() == other_iter->second.IsFirstMate())
	    continue;
	  potential_mates.erase(other_iter);
	}
	else
	  potential_mates.insert(std::pair<std::string, BamAlignment>(aln_key, alignment));
      }
    }
  }

  for (auto aln_iter = potential_strs.begin(); aln_iter != potential_strs.end(); ++aln_iter){
    std::string filter = "";
    if (aln_iter->second.HasTag(ALT_MAP_TAG.c_str())){
      unique_mapping++;
      filter = "NO_UNIQUE_MAPPING";
    }
    else if (REQUIRE_PAIRED_READS){
      num_filt_unpaired_reads++;
      filter = "NO_MATE_PAIR";
    }

    if (filter.empty()){
      unpaired_str_alns.push_back(aln_iter->second);
      write_passing_alignment(aln_iter->second, pass_writer);
    }
    else
      write_filtered_alignment(aln_iter->second, filter, filt_writer);
  }
  potential_strs.clear(); potential_mates.clear();

  selective_logger() << adapter_trimmer_.get_trimming_stats_msg() << "\n"
		     << read_count << " reads overlapped region, of which "
		     << "\n\t" << hard_clip      << " were hard clipped"
		     << "\n\t" << read_has_N     << " had an 'N' base call"
		     << "\n\t" << low_qual_score << " had low base quality scores"
		     << "\n\t" << unique_mapping << " did not have a unique mapping";
  if (REQUIRE_PAIRED_READS)
    selective_logger() << "\n\t" << num_filt_unpaired_reads << " did not have a mate pair";
  selective_logger() << "\n\t" << (paired_str_alns.size()+unpaired_str_alns.size()) << " PASSED ALL FILTERS" << "\n"
		     << "Found " << paired_str_alns.size() << " fully paired reads and " << unpaired_str_alns.size() << " unpaired reads for downstream analyses" << std::endl;
    
  // Separate the reads based on their associated read groups
  std::map<std::string, int> rg_indices;
  for (unsigned int type = 0; type < 2; ++type){
    BamAlnList& aln_src  = (type == 0 ? paired_str_alns : unpaired_str_alns);
    while (!aln_src.empty()){
      const BamAlignment& aln = aln_src.back();
      std::string rg = use_bam_rgs_ ? get_read_group(aln, rg_to_sample): rg_to_sample.find(aln.Filename())->second;
      int rg_index;
      auto index_iter = rg_indices.find(rg);
      if (index_iter == rg_indices.end()){
	rg_index = rg_indices.size();
	rg_indices[rg] = rg_index;
	rg_names.push_back(rg);
	paired_strs_by_rg.push_back(BamAlnList());
	unpaired_strs_by_rg.push_back(BamAlnList());
	mate_pairs_by_rg.push_back(BamAlnList());
      }
      else
	rg_index = index_iter->second;

      // Record STR read and its mate pair
      if (type == 0){
	paired_strs_by_rg[rg_index].push_back(aln);
	mate_pairs_by_rg[rg_index].push_back(mate_alns.back());
	mate_alns.pop_back();
      }
      // Record unpaired STR read
      else
	unpaired_strs_by_rg[rg_index].push_back(aln);
      aln_src.pop_back();
    }
  }

  locus_read_filter_time_  = (clock() - locus_read_filter_time_)/CLOCKS_PER_SEC;
  total_read_filter_time_ += locus_read_filter_time_;
}

// Ensure that all of the chromosomes are present in i) the FASTA file, ii) the BAM files and iii) the SNP VCF file, if provided
void BamProcessor::verify_chromosomes(const std::vector<std::string>& chroms, const BamHeader* bam_header, FastaReader& fasta_reader){
  for (auto chrom_iter = chroms.begin(); chrom_iter != chroms.end(); chrom_iter++){
    std::string chrom = (*chrom_iter);

    std::vector<std::string> alt_names(1, "chr" + chrom);
    if (chrom.size() > 3 && chrom.substr(0, 3).compare("chr") == 0)
      alt_names.push_back(chrom.substr(3));

    // 1. Check FASTA
    if (fasta_reader.get_sequence_length(chrom) == -1){
      std::stringstream err_msg;
      err_msg << "No sequence for chromosome " << chrom << " found in the FASTA file" << "\n"
	      << "\t" << "Please ensure that the chromosome names in your region BED file match those in you FASTA file";
      full_logger() << "\n" << "ERROR: " << err_msg.str() << std::endl;

      // Prompt if simple changes to the chromosome name would solve the issue
      for (auto alt_iter = alt_names.begin(); alt_iter != alt_names.end(); alt_iter++)
	if (fasta_reader.get_sequence_length(*alt_iter) != -1)
	  full_logger() << "\t" << "NOTE: Found chromosome " << (*alt_iter) << " in the FASTA, but not chromosome " << chrom << std::endl;

      // Abort execution
      printErrorAndDie("Terminating HipSTR as chromosomes in the region file are missing from the FASTA file. Please see the log for details");
    }

    // 2. Check BAMs
    if (bam_header->ref_id(chrom) == -1){
      std::stringstream err_msg;
      err_msg << "No entries for chromosome " << chrom << " found in the BAM/CRAM(s)" << "\n"
	      << "\t" << "Please ensure that the chromosome names in your region BED file match those in your BAM/CRAM(s)";
      full_logger() << "\n" << "ERROR: " << err_msg.str() << std::endl;

      // Prompt if simple changes to the chromosome name would solve the issue
      for (auto alt_iter = alt_names.begin(); alt_iter != alt_names.end(); alt_iter++)
	if (bam_header->ref_id(*alt_iter) != -1)
	  full_logger() << "\t" << "NOTE: Found chromosome " << (*alt_iter) << " in the BAM/CRAM(s), but not chromosome " << chrom << std::endl;

      // Abort execution
      printErrorAndDie("Terminating HipSTR as chromosomes in the region file are missing from the BAM/CRAM(s). Please see the log for details");
    }
  }

  // 3. Check SNP VCF file
  verify_vcf_chromosomes(chroms);
}


void BamProcessor::process_regions(BamCramMultiReader& reader, const std::string& region_file, const std::string& fasta_file,
				   const std::map<std::string, std::string>& rg_to_sample, const std::map<std::string, std::string>& rg_to_library, const std::string& full_command,
				   BamWriter* pass_writer, BamWriter* filt_writer, int32_t max_regions, const std::string& chrom){
  std::vector<Region> regions;
  readRegions(region_file, max_regions, chrom, regions, full_logger());
  orderRegions(regions);

  FastaReader fasta_reader(fasta_file);
  const BamHeader* bam_header = reader.bam_header();

  // Collect a list of all chromosomes present in the region file
  std::vector<std::string> chroms;
  std::string prev_chrom = "";
  for (auto region_iter = regions.begin(); region_iter != regions.end(); region_iter++){
    if (region_iter->chrom().compare(prev_chrom) != 0){
      prev_chrom = region_iter->chrom();
      chroms.push_back(prev_chrom);
    }
  }

  // Ensure consistent chromosome naming between the relevant input files
  verify_chromosomes(chroms, bam_header, fasta_reader);

  // Add the chromosome information to the VCF
  init_output_vcf(fasta_file, chroms, full_command);

  std::string cur_chrom = "", chrom_seq = "";
  for (auto region_iter = regions.begin(); region_iter != regions.end(); region_iter++){
    full_logger() << "" << "Processing region " << region_iter->chrom() << " " << region_iter->start() << " " << region_iter->stop() << std::endl;

    if (region_iter->stop() - region_iter->start() > MAX_STR_LENGTH){
      num_too_long_++;
      full_logger() << "Skipping region as the reference allele length exceeds the threshold (" 
		    << region_iter->stop()-region_iter->start() << " vs " << MAX_STR_LENGTH << ")" << "\n"
		    << "You can increase this threshold using the --max-str-len option" << std::endl;
      continue;
    }
    
    // Read FASTA sequence for chromosome 
    if (region_iter->chrom().compare(cur_chrom) != 0){
      cur_chrom = region_iter->chrom();
      fasta_reader.get_sequence(cur_chrom, chrom_seq);
      assert(chrom_seq.size() != 0);
    }

    if (region_iter->start() < 50 || region_iter->stop()+50 >= chrom_seq.size()){
      full_logger() << "Skipping region within 50bp of the end of the contig" << std::endl;
      continue;
    }

    locus_bam_seek_time_ = clock();
    if (!reader.SetRegion(cur_chrom, (region_iter->start() < MAX_MATE_DIST ? 0: region_iter->start()-MAX_MATE_DIST),
			  region_iter->stop() + MAX_MATE_DIST))
      printErrorAndDie("One or more BAM files failed to set the region properly");

    locus_bam_seek_time_  =  (clock() - locus_bam_seek_time_)/CLOCKS_PER_SEC;
    total_bam_seek_time_ += locus_bam_seek_time_;

    std::vector<std::string> rg_names;
    std::vector<BamAlnList> paired_strs_by_rg, mate_pairs_by_rg, unpaired_strs_by_rg;
    RegionGroup region_group(*region_iter); // TO DO: Extend region groups to have multiple regions
    read_and_filter_reads(reader, chrom_seq, region_group, rg_to_sample, rg_names,
			  paired_strs_by_rg, mate_pairs_by_rg, unpaired_strs_by_rg, pass_writer, filt_writer);

    // The user specified a list of samples to which we need to restrict the analyses
    // Discard reads for any samples not in this set
    if (!sample_set_.empty()){
      selective_logger() << "Restricting reads to the " << sample_set_.size() << " samples in the specified sample list" << std::endl;
      unsigned int ins_index = 0;
      for (unsigned int i = 0; i < rg_names.size(); i++){
	if (sample_set_.find(rg_names[i]) != sample_set_.end()){
	  if (i != ins_index){
	    rg_names[ins_index]            = rg_names[i];
	    paired_strs_by_rg[ins_index]   = paired_strs_by_rg[i];
	    mate_pairs_by_rg[ins_index]    = mate_pairs_by_rg[i];
	    unpaired_strs_by_rg[ins_index] = unpaired_strs_by_rg[i];
	  }
	  ins_index++;
	}
      }
      if (ins_index != rg_names.size()){
	rg_names.resize(ins_index);
	paired_strs_by_rg.resize(ins_index);
	mate_pairs_by_rg.resize(ins_index);
	unpaired_strs_by_rg.resize(ins_index);
      }
    }

    if (REMOVE_PCR_DUPS == 1)
      remove_pcr_duplicates(base_quality_, use_bam_rgs_, rg_to_library, paired_strs_by_rg, mate_pairs_by_rg, unpaired_strs_by_rg, selective_logger());

    process_reads(paired_strs_by_rg, mate_pairs_by_rg, unpaired_strs_by_rg, rg_names, region_group, chrom_seq);

    adapter_trimmer_.mark_new_locus(); // Inform the trimmer that future alignments will be for a new STR
  }
}
