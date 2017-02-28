#include <fstream>
#include <iostream>
#include <locale>
#include <sstream>
#include <stdlib.h>
#include <time.h>

#include "bam_processor.h"
#include "alignment_filters.h"
#include "error.h"
#include "fasta_reader.h"
#include "pcr_duplicates.h"
#include "stringops.h"
#include "SeqAlignment/AlignmentOps.h"

const std::string ALT_MAP_TAG           = "XA";
const std::string PRIMARY_ALN_SCORE_TAG = "AS";
const std::string SUBOPT_ALN_SCORE_TAG  = "XS";

void BamProcessor::add_passes_filters_tag(BamAlignment& aln, std::string& passes){
  if (aln.HasTag("PF"))
    aln.RemoveTag("PF");
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

void BamProcessor::add_filtered_alignment(BamAlignment& alignment, std::string filter, BamAlnList& filtered_alignments){
  filtered_alignments.push_back(alignment);
  if (!filtered_alignments.back().AddStringTag("FT", filter))
    printErrorAndDie("Failed to add filter tag to alignment");
}

void BamProcessor::extract_mappings(BamAlignment& aln, const BamHeader* bam_header,
				    std::vector< std::pair<std::string, int32_t> >& chrom_pos_pairs){
  assert(chrom_pos_pairs.size() == 0);
  if (aln.RefID() == -1 || aln.CigarData().size() == 0)
    return;
  assert(aln.RefID() < bam_header->num_seqs());
  chrom_pos_pairs.push_back(std::pair<std::string, int32_t>(bam_header->ref_name(aln.RefID()), aln.Position()));

  if (aln.HasTag(ALT_MAP_TAG.c_str())){
    std::string alt_info;
    if (!aln.GetStringTag(ALT_MAP_TAG.c_str(), alt_info))
      printErrorAndDie("Failed to extract XA tag from BAM alignment");
    std::vector<std::string> alts;
    split_by_delim(alt_info, ';', alts);
    for (unsigned int i = 0; i < alts.size(); i++){
      std::vector<std::string> tokens;
      split_by_delim(alts[i], ',', tokens);
      chrom_pos_pairs.push_back(std::pair<std::string, int32_t>(tokens[0], std::abs(std::stol(tokens[1]))));
    }
  }
}

void BamProcessor::get_valid_pairings(BamAlignment& aln_1, BamAlignment& aln_2, const BamHeader* bam_header,
				      std::vector< std::pair<std::string, int32_t> >& p1, std::vector< std::pair<std::string, int32_t> >& p2){
  assert(p1.size() == 0 && p2.size() == 0);
  if (aln_1.RefID() == -1 || aln_2.RefID() == -1)
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
  extract_mappings(aln_1, bam_header, pairs_1);
  extract_mappings(aln_2, bam_header, pairs_2);
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

std::string BamProcessor::get_read_group(BamAlignment& aln, std::map<std::string, std::string>& read_group_mapping){
  std::string rg;
  if (!aln.GetStringTag("RG", rg))
    printErrorAndDie("Failed to retrieve BAM alignment's RG tag");
  auto iter = read_group_mapping.find(aln.Filename() + rg);
  if (iter == read_group_mapping.end())
    printErrorAndDie("No sample found for read group " + rg + " in BAM file headers");
  return iter->second;
}

std::string BamProcessor::trim_alignment_name(BamAlignment& aln){
  std::string aln_name = aln.Name();
  if (aln_name.size() > 2){
    if (aln_name[aln_name.size()-2] == '/')
      aln_name.resize(aln_name.size()-2);
  }
  return aln_name;
}

void BamProcessor::modify_and_write_alns(BamAlnList& alignments, std::map<std::string, std::string>& rg_to_sample,
					 BamWriter* writer){
  for (auto read_iter = alignments.begin(); read_iter != alignments.end(); read_iter++){
    // Add RG to BAM record based on file
    if (!use_bam_rgs_){
      std::string rg_tag = "HipSTR;" + rg_to_sample[read_iter->Filename()] + ";" + rg_to_sample[read_iter->Filename()];
      read_iter->AddStringTag("RG", rg_tag);
    }
    if (!writer->SaveAlignment(*read_iter))
      printErrorAndDie("Failed to save alignment for STR-spanning read");
  }
}

// Returns true if the alignment spans at least one region in the group.
// Reads that don't actually span the region but start/end with a soft clip are also counted as spanning,
// as the clipped sequence frequently extends past the region
bool BamProcessor::spans_a_region(const std::vector<Region>& regions, BamAlignment& alignment){
  for (auto region_iter = regions.begin(); region_iter != regions.end(); region_iter++){
    if (alignment.Position() > region_iter->stop() || alignment.GetEndPosition() < region_iter->start())
      continue;
    if (alignment.Position() > region_iter->start() && !alignment.StartsWithSoftClip())
      continue;
    if (alignment.GetEndPosition() < region_iter->stop() && !alignment.EndsWithSoftClip())
      continue;
    return true;
  }
  return false;
}

void BamProcessor::read_and_filter_reads(BamCramMultiReader& reader, std::string& chrom_seq, RegionGroup& region_group,
					 std::map<std::string, std::string>& rg_to_sample, std::map<std::string, std::string>& rg_to_library, std::vector<std::string>& rg_names,
					 std::vector<BamAlnList>& paired_strs_by_rg, std::vector<BamAlnList>& mate_pairs_by_rg, std::vector<BamAlnList>& unpaired_strs_by_rg,
					 BamWriter* pass_writer, BamWriter* filt_writer){
  locus_read_filter_time_ = clock();
  assert(reader.get_merge_type() == BamCramMultiReader::ORDER_ALNS_BY_FILE);

  bool pass_to_bam     = (pass_writer != NULL);
  bool filtered_to_bam = (filt_writer != NULL);
  BamAlnList region_alignments, filtered_alignments;
  int32_t read_count = 0, not_spanning = 0, unique_mapping = 0, read_has_N = 0, hard_clip = 0, split_alignment = 0, low_qual_score = 0, num_filt_unpaired_reads = 0;
  BamAlignment alignment;
  const BamHeader* bam_header = reader.bam_header();
  BamAlnList paired_str_alns, mate_alns, unpaired_str_alns;
  std::map<std::string, BamAlignment> potential_strs, potential_mates;
  TOO_MANY_READS = false;

  const std::vector<Region>& regions = region_group.regions();
  std::string prev_file = "";
  while (reader.GetNextAlignment(alignment)){
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
    assert(alignment.CigarData().size() > 0 && alignment.RefID() != -1);

    // If requested, trim any reads that potentially overlap the STR regions
    if (alignment.Position() < region_group.stop() && alignment.GetEndPosition() >= region_group.start()){
      if (BASE_QUAL_TRIM > ' '){
	// Read trimming doesn't work if hard clipping has been applied
	if (alignment.StartsWithHardClip() || alignment.EndsWithHardClip()){
	  read_count++;
	  hard_clip++;
	  if (filtered_to_bam) add_filtered_alignment(alignment, "HARD_CLIPPED", filtered_alignments);
	  continue;
	}

	int32_t length = alignment.Length();
	alignment.TrimLowQualityEnds(BASE_QUAL_TRIM);
	if (alignment.Position() < region_group.stop() && alignment.GetEndPosition() >= region_group.start())
	  if ((alignment.Length() == 0) || (alignment.Length() < length/2))
	    continue;
      }
    }

    // Clear out mate alignment cache if we've switched to a new file to reduce memory
    if (prev_file.compare(alignment.Filename()) != 0){
      prev_file = alignment.Filename();
      potential_mates.clear();
    }

    // Only apply tests to putative STR reads that overlap the STR region
    if (alignment.Position() < region_group.stop() && alignment.GetEndPosition() >= region_group.start()){
      bool pass_one = false; // Denotes if read passed first set of simpler filters
      std::string pass_two(regions.size(), '0'); // Denotes if read passed sceond set of additional filters for each region
                                                 // Meant to signify if reads that pass first set should be used to generate haplotypes
      std::string filter = "";
      read_count++;

      // Ignore chimeric alignments
      if (alignment.HasTag("SA")){
	split_alignment++;
	filter.append("HAS_SA_TAG");
      }
      // Ignore reads with N bases
      else if (alignment.QueryBases().find('N') != std::string::npos){
	read_has_N++;
	filter.append("HAS_N_BASES");
      }
      // Ignore reads with a very low overall base quality score
      // Want to avoid situations in which it's more advantageous to have misalignments b/c the scores are so low
      else if (base_quality_.sum_log_prob_correct(alignment.Qualities()) < MIN_SUM_QUAL_LOG_PROB){
	low_qual_score++;
	filter.append("LOW_BASE_QUALS");
      }
      else if (REQUIRE_SPANNING && !spans_a_region(regions, alignment)){
	not_spanning++;
	filter.append("NOT_SPANNING");
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

      std::string aln_key = trim_alignment_name(alignment);
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
	  get_valid_pairings(alignment, aln_iter->second, bam_header, p_1, p_2);
	  if (p_1.size() == 1 && p_1[0].second == alignment.Position()){
	    paired_str_alns.push_back(alignment);
	    mate_alns.push_back(aln_iter->second);
	    if (pass_to_bam){
	      region_alignments.push_back(alignment);
	      region_alignments.push_back(aln_iter->second);
	    }
	  }
	  else {
	    unique_mapping++;
	    filter.append("NO_UNIQUE_MAPPING");
	    if (filtered_to_bam)
	      add_filtered_alignment(alignment, filter, filtered_alignments);
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
	    get_valid_pairings(alignment, str_iter->second, bam_header, p_1, p_2);
	    if (p_1.size() == 1 && p_1[0].second == alignment.Position()){
	      paired_str_alns.push_back(alignment);
	      mate_alns.push_back(str_iter->second);
	      if (pass_to_bam) region_alignments.push_back(alignment);

	      paired_str_alns.push_back(str_iter->second);
	      mate_alns.push_back(alignment);
	      if (pass_to_bam) region_alignments.push_back(str_iter->second);
	    }
	    else {
	      unique_mapping += 2;
	      std::string filter = "NO_UNIQUE_MAPPING";
	      if (filtered_to_bam){
		add_filtered_alignment(alignment, filter, filtered_alignments);
		add_filtered_alignment(str_iter->second, filter, filtered_alignments);
	      }
	    }
	    potential_strs.erase(str_iter);
	  }
	  else
	    potential_strs.insert(std::pair<std::string, BamAlignment>(aln_key, alignment));
	}
      }
      else {
	assert(!filter.empty());
	if (filtered_to_bam)
	  add_filtered_alignment(alignment, filter, filtered_alignments);
	potential_mates.insert(std::pair<std::string, BamAlignment>(aln_key, alignment));
      }
    }
    else {
      std::string aln_key = trim_alignment_name(alignment);
      auto aln_iter = potential_strs.find(aln_key);
      if (aln_iter != potential_strs.end()){
	if (alignment.IsFirstMate() == aln_iter->second.IsFirstMate())
	  continue;

	std::vector< std::pair<std::string, int32_t> > p_1, p_2;
	get_valid_pairings(aln_iter->second, alignment, bam_header, p_1, p_2);
	if (p_1.size() == 1 && p_1[0].second == aln_iter->second.Position()){
	  paired_str_alns.push_back(aln_iter->second);
	  mate_alns.push_back(alignment);
	  if (pass_to_bam){
	    region_alignments.push_back(aln_iter->second);
	    region_alignments.push_back(alignment);
	  }
	}
	else {
	  unique_mapping++;
	  std::string filter = "NO_UNIQUE_MAPPING";
	  if (filtered_to_bam)
	    add_filtered_alignment(aln_iter->second, filter, filtered_alignments);
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
      if (pass_to_bam) region_alignments.push_back(aln_iter->second);
    }
    else {
      if (filtered_to_bam)
	add_filtered_alignment(aln_iter->second, filter, filtered_alignments);
    }
  }
  potential_strs.clear(); potential_mates.clear();
  
  logger() << read_count << " reads overlapped region, of which "
	   << "\n\t" << hard_clip        << " were hard clipped"
	   << "\n\t" << split_alignment  << " had an SA (split alignment) BAM tag"
	   << "\n\t" << read_has_N       << " had an 'N' base call"
	   << "\n\t" << low_qual_score   << " had low base quality scores"
	   << "\n\t" << not_spanning     << " did not span the STR"
	   << "\n\t" << unique_mapping   << " did not have a unique mapping";
  if (REQUIRE_PAIRED_READS)
    logger() << "\n\t" << num_filt_unpaired_reads << " did not have a mate pair";
  logger() << "\n\t" << (paired_str_alns.size()+unpaired_str_alns.size()) << " PASSED ALL FILTERS" << "\n"
	   << "Found " << paired_str_alns.size() << " fully paired reads and " << unpaired_str_alns.size() << " unpaired reads for downstream analyses" << std::endl;
    
  // Output the reads passing all filters to a BAM file (if requested)
  if (pass_writer != NULL)
    modify_and_write_alns(region_alignments, rg_to_sample, pass_writer);

  // Output reads that overlapped the STR but were filtered to a BAM file (if requested)
  if (filt_writer != NULL)
    modify_and_write_alns(filtered_alignments, rg_to_sample, filt_writer);

  // Separate the reads based on their associated read groups
  std::map<std::string, int> rg_indices;
  for (unsigned int type = 0; type < 2; ++type){
    BamAlnList& aln_src  = (type == 0 ? paired_str_alns : unpaired_str_alns);
    for (unsigned int i = 0; i < aln_src.size(); ++i){
      std::string rg = use_bam_rgs_ ? get_read_group(aln_src[i], rg_to_sample): rg_to_sample[aln_src[i].Filename()];
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
	paired_strs_by_rg[rg_index].push_back(aln_src[i]);
	mate_pairs_by_rg[rg_index].push_back(mate_alns[i]);
      }
      // Record unpaired STR read
      else
	unpaired_strs_by_rg[rg_index].push_back(aln_src[i]);
    }
  }

  locus_read_filter_time_  = (clock() - locus_read_filter_time_)/CLOCKS_PER_SEC;
  total_read_filter_time_ += locus_read_filter_time_;
}

void BamProcessor::process_regions(BamCramMultiReader& reader, std::string& region_file, std::string& fasta_dir,
				   std::map<std::string, std::string>& rg_to_sample, std::map<std::string, std::string>& rg_to_library,
				   BamWriter* pass_writer, BamWriter* filt_writer,
				   std::ostream& out, int32_t max_regions, std::string chrom){
  std::vector<Region> regions;
  readRegions(region_file, regions, max_regions, chrom, logger());
  orderRegions(regions);

  FastaReader fasta_reader(fasta_dir);
  std::string ref_seq;
  const BamHeader* bam_header = reader.bam_header();
  int cur_chrom_id = -1; std::string chrom_seq;
  for (auto region_iter = regions.begin(); region_iter != regions.end(); region_iter++){
    logger() << "\n\n" << "Processing region " << region_iter->chrom() << " " << region_iter->start() << " " << region_iter->stop() << std::endl;
    int chrom_id = bam_header->ref_id(region_iter->chrom());
    if (chrom_id == -1 && region_iter->chrom().size() > 3 && region_iter->chrom().substr(0, 3).compare("chr") == 0)
      chrom_id = bam_header->ref_id(region_iter->chrom().substr(3));

    if (chrom_id == -1){
      logger() << "\n" << "WARNING: No reference sequence for chromosome " << region_iter->chrom() << " found in BAMs"  << "\n"
	       << "\t" << "Please ensure that the names of reference sequences in your BED file match those in you BAMs" << "\n"
	       << "\t" << "Skipping region " << region_iter->chrom() << " " << region_iter->start() << " " << region_iter->stop() << "\n" << std::endl;
      continue;
    }

    if (region_iter->stop() - region_iter->start() > MAX_STR_LENGTH){
      logger() << "Skipping region as the reference allele length exceeds the threshold (" << region_iter->stop()-region_iter->start() << " vs " << MAX_STR_LENGTH << ")" << "\n"
	       << "You can increase this threshold using the --max-str-length option" << std::endl;
      continue;
    }
    
    // Read FASTA sequence for chromosome 
    if (cur_chrom_id != chrom_id){
      cur_chrom_id      = chrom_id;
      std::string chrom = region_iter->chrom();
      fasta_reader.get_sequence(chrom, chrom_seq);
      assert(chrom_seq.size() != 0);
    }

    if (region_iter->start() < 50 || region_iter->stop()+50 >= chrom_seq.size()){
      logger() << "Skipping region within 50bp of the end of the contig" << std::endl;
      continue;
    }

    locus_bam_seek_time_ = clock();
    if (!reader.SetRegion(bam_header->ref_name(chrom_id), (region_iter->start() < MAX_MATE_DIST ? 0: region_iter->start()-MAX_MATE_DIST),
			  region_iter->stop() + MAX_MATE_DIST))
      printErrorAndDie("One or more BAM files failed to set the region properly");

    locus_bam_seek_time_  =  (clock() - locus_bam_seek_time_)/CLOCKS_PER_SEC;
    total_bam_seek_time_ += locus_bam_seek_time_;

    std::vector<std::string> rg_names;
    std::vector<BamAlnList> paired_strs_by_rg, mate_pairs_by_rg, unpaired_strs_by_rg;
    RegionGroup region_group(*region_iter); // TO DO: Extend region groups to have multiple regions
    read_and_filter_reads(reader, chrom_seq, region_group, rg_to_sample, rg_to_library, rg_names,
			  paired_strs_by_rg, mate_pairs_by_rg, unpaired_strs_by_rg, pass_writer, filt_writer);

    // The user specified a list of samples to which we need to restrict the analyses
    // Discard reads for any samples not in this set
    if (!sample_set_.empty()){
      logger() << "Restricting reads to the " << sample_set_.size() << " samples in the specified sample list" << std::endl;
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

    if (rem_pcr_dups_)
      remove_pcr_duplicates(base_quality_, use_bam_rgs_, rg_to_library, paired_strs_by_rg, mate_pairs_by_rg, unpaired_strs_by_rg, logger());

    process_reads(paired_strs_by_rg, mate_pairs_by_rg, unpaired_strs_by_rg, rg_names, region_group, chrom_seq, out);
  }
}
