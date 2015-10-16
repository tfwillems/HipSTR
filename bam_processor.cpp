#include <fstream>
#include <iostream>
#include <locale>
#include <sstream>
#include <stdlib.h>
#include <time.h>

#include "fastahack/Fasta.h"

#include "bam_processor.h"
#include "alignment_filters.h"
#include "error.h"
#include "pcr_duplicates.h"
#include "seqio.h"
#include "stringops.h"

const std::string ALT_MAP_TAG = "XA";

void BamProcessor::extract_mappings(BamTools::BamAlignment& aln, const BamTools::RefVector& ref_vector,
				    std::vector< std::pair<std::string, int32_t> >& chrom_pos_pairs){
  assert(chrom_pos_pairs.size() == 0);
  if (aln.RefID == -1 || aln.CigarData.size() == 0)
    return;
  assert(aln.RefID < ref_vector.size());
  chrom_pos_pairs.push_back(std::pair<std::string, int32_t>(ref_vector[aln.RefID].RefName, aln.Position));

  if (aln.HasTag(ALT_MAP_TAG)){
    std::string alt_info;
    if (!aln.GetTag(ALT_MAP_TAG, alt_info))
      printErrorAndDie("Failed to extract XA tag from BAM alignment");
    std::vector<std::string> alts;
    split_by_delim(alt_info, ';', alts);
    for (unsigned int i = 0; i < alts.size(); i++){
      std::vector<std::string> tokens;
      split_by_delim(alts[i], ',', tokens);
      chrom_pos_pairs.push_back(std::pair<std::string, int32_t>(tokens[0], abs(std::stol(tokens[1]))));
    }
  }
}

void BamProcessor::get_valid_pairings(BamTools::BamAlignment& aln_1, BamTools::BamAlignment& aln_2, const BamTools::RefVector& ref_vector,
				      std::vector< std::pair<std::string, int32_t> >& p1, std::vector< std::pair<std::string, int32_t> >& p2){
  assert(p1.size() == 0 && p2.size() == 0);
  if (aln_1.RefID == -1 || aln_2.RefID == -1)
    return;

  std::vector< std::pair<std::string, int32_t> > pairs_1, pairs_2;
  extract_mappings(aln_1, ref_vector, pairs_1);
  extract_mappings(aln_2, ref_vector, pairs_2);
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

std::string BamProcessor::get_read_group(BamTools::BamAlignment& aln, std::map<std::string, std::string>& read_group_mapping){
  std::string rg;
  std::string rg_tag = "RG";
  char tag_type = 'Z';
  if (!aln.GetTagType(rg_tag, tag_type))
    printErrorAndDie("Failed to retrieve BAM alignment's RG tag");
  aln.GetTag("RG", rg);
  auto iter = read_group_mapping.find(rg);
  if (iter == read_group_mapping.end())
    printErrorAndDie("No sample found for read group " + rg + " in BAM file headers");
  return iter->second;
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


void BamProcessor::modify_and_write_alns(std::vector<BamTools::BamAlignment>& alignments,
					 std::map<std::string, std::string>& rg_to_sample,
					 Region& region,
					 BamTools::BamWriter& writer){
  for (auto read_iter = alignments.begin(); read_iter != alignments.end(); read_iter++){
    // Add RG to BAM record based on file
    if (!use_bam_rgs_){
      std::string rg_tag = "lobSTR;" + rg_to_sample[read_iter->Filename] + ";" + rg_to_sample[read_iter->Filename];
      read_iter->AddTag("RG", "Z", rg_tag);
    }

    // Add STR start and stop tags
    if (read_iter->HasTag("XS"))
      read_iter->RemoveTag("XS");
    if(!read_iter->AddTag("XS",  "i", region.start()))
      printErrorAndDie("Failed to modify XS tag");
    if (read_iter->HasTag("XE"))
      read_iter->RemoveTag("XE");
    if(!read_iter->AddTag("XE", "i", region.stop()))
      printErrorAndDie("Failed to modify XE tag");
    if (!writer.SaveAlignment(*read_iter))
      printErrorAndDie("Failed to save alignment for STR-spanning read");
  }

}

// TO DO: Track reads from observed from both strands, whether or not they're filtered
void BamProcessor::read_and_filter_reads(BamTools::BamMultiReader& reader, std::string& chrom_seq, 
					 std::vector<Region>::iterator region_iter,
					 std::map<std::string, std::string>& rg_to_sample, std::map<std::string, std::string>& rg_to_library,
					 std::vector<std::string>& rg_names,
					 std::vector< std::vector<BamTools::BamAlignment> >& paired_strs_by_rg,
					 std::vector< std::vector<BamTools::BamAlignment> >& mate_pairs_by_rg,
					 std::vector< std::vector<BamTools::BamAlignment> >& unpaired_strs_by_rg,
					 BamTools::BamWriter& pass_writer, BamTools::BamWriter& filt_writer){
  locus_read_filter_time_ = clock();

  std::vector<BamTools::BamAlignment> region_alignments, filtered_alignments;
  int32_t read_count = 0;
  int32_t diff_chrom_mate = 0, unmapped_mate = 0, not_spanning = 0; // Counts for filters that are always applied
  int32_t insert_size = 0, multimapped = 0, mapping_quality = 0, flank_len = 0; // Counts for filters that are user-controlled
  int32_t bp_before_indel = 0, end_match_window = 0, num_end_matches = 0, read_has_N = 0, hard_clip = 0, soft_clip = 0, split_alignment = 0, low_qual_score = 0;
  int32_t unique_mapping = 0;
  BamTools::BamAlignment alignment;

  const BamTools::RefVector& ref_vector = reader.GetReferenceData();
  std::vector<BamTools::BamAlignment> paired_str_alns, mate_alns, unpaired_str_alns;
  std::map<std::string, BamTools::BamAlignment> potential_strs, potential_mates;
  const std::string FILTER_TAG_NAME = "FT";
  const std::string FILTER_TAG_TYPE = "Z";

  while (reader.GetNextAlignment(alignment)){
    // Stop parsing reads if we've already exceeded the maximum number for downstream analyses
    if (paired_str_alns.size() > MAX_TOTAL_READS)
      break;

    if (!alignment.IsMapped() || alignment.Position == 0 || alignment.CigarData.size() == 0)
	continue;
    assert(alignment.CigarData.size() > 0 && alignment.RefID != -1);

    // Only apply tests to putative STR reads that overlap the STR region
    if (alignment.Position < region_iter->stop() && alignment.GetEndPosition() >= region_iter->start()){
      bool pass = true;
      std::string filter = "";
      read_count++;

      if (check_mate_info_){
	// Ignore read if its mate pair chromosome doesn't match
	if (pass && alignment.RefID != alignment.MateRefID){
	  diff_chrom_mate++;
	  pass = false;
	  filter.append("MATE_DIFF_CHROM");
	}
	// Ignore read if its mate pair is unmapped
	// TO DO: Check that insert size != 0 ?
	if (pass && (!alignment.IsMateMapped() || alignment.MatePosition == 0)){
	  unmapped_mate++;
	  pass = false;
	  filter.append("MATE_UNMAPPED");
	}
      }
      // Ignore read if its mate pair distance exceeds the threshold
      if (pass && abs(alignment.InsertSize) > MAX_MATE_DIST){
	insert_size++;
	pass = false;
	filter.append("INSERT_SIZE");
      }
      // Ignore read if multimapper and filter specified
      if (pass && REMOVE_MULTIMAPPERS && alignment.HasTag("XA")){
	multimapped++;
	pass = false;
	filter.append("MULTIMAPPED");
      }
      // Ignore chimeric alignments
      if (pass && alignment.HasTag("SA")){
	split_alignment++;
	pass = false;
	filter.append("HAS_SA_TAG");
      }
      // Ignore reads with N bases
      if (pass && REMOVE_READS_WITH_N && (alignment.QueryBases.find('N') != std::string::npos)){
	read_has_N++;
	pass = false;
	filter.append("HAS_N_BASES");
      }
      // Ignore read if its mapping quality is too low
      if (pass && MIN_MAPPING_QUALITY > alignment.MapQuality){
	mapping_quality++;
	pass = false;
	filter.append("LOW_MAPQ");
      }
      // Ignore reads with too many clipped bases
      int num_hard_clips, num_soft_clips;
      AlignmentFilters::GetNumClippedBases(alignment, num_hard_clips, num_soft_clips);
      if (pass && num_hard_clips > MAX_HARD_CLIPS){
	hard_clip++;
	pass = false;
	filter.append("NUM_HCLIPS");
      }
      if (pass && num_soft_clips > MAX_SOFT_CLIPS){
	soft_clip++;
	pass = false;
	filter.append("NUM_SCLIPS");
      }
      // Ignore read if it does not span the STR
      if (pass && REQUIRE_SPANNING && (alignment.Position > region_iter->start() || alignment.GetEndPosition() < region_iter->stop())){
	not_spanning++;
	pass = false;
	filter.append("NOT_SPANNING");
      }
      // Ignore read if it has insufficient flanking bases on either side of the STR
      if (pass && (alignment.Position > (region_iter->start()-MIN_FLANK) || alignment.GetEndPosition() < (region_iter->stop()+MIN_FLANK))){
	flank_len++;
	pass = false;
	filter.append("FLANK_LEN");
      }
      // Ignore read if there is an indel within the first MIN_BP_BEFORE_INDEL bps from each end
      if (pass && MIN_BP_BEFORE_INDEL > 0){
	std::pair<int, int> num_bps = AlignmentFilters::GetEndDistToIndel(alignment);
	if ((num_bps.first != -1 && num_bps.first < MIN_BP_BEFORE_INDEL) || (num_bps.second != -1 && num_bps.second < MIN_BP_BEFORE_INDEL)){
	  bp_before_indel++;
	  pass = false;
	  filter.append("BP_BEFORE_INDEL");
	}
      }
      // Ignore read if there is another location within MAXIMAL_END_MATCH_WINDOW bp for which it has a longer end match
      if (pass && MAXIMAL_END_MATCH_WINDOW > 0){
	bool maximum_end_matches = AlignmentFilters::HasLargestEndMatches(alignment, chrom_seq, 0, MAXIMAL_END_MATCH_WINDOW, MAXIMAL_END_MATCH_WINDOW);
	if (!maximum_end_matches){
	  end_match_window++;
	  pass = false;
	  filter.append("END_MATCHES_NOT_MAXIMAL");
	}
      }
      // Ignore read if it doesn't match perfectly for at least MIN_READ_END_MATCH bases on each end
      if (pass && MIN_READ_END_MATCH > 0){
	std::pair<int,int> match_lens = AlignmentFilters::GetNumEndMatches(alignment, chrom_seq, 0);
	if (match_lens.first < MIN_READ_END_MATCH || match_lens.second < MIN_READ_END_MATCH){
	  num_end_matches++;
	  pass = false;
	  filter.append("NUM_END_MATCHES");
	}
      }
      // Ignore reads with a very low overall base quality score
      // Want to avoid situations in which it's more advantageous to have misalignments b/c the scores are so low
      if (pass && base_quality_.sum_log_prob_correct(alignment.Qualities) < MIN_SUM_QUAL_LOG_PROB){
	low_qual_score++;
	pass = false;
	filter.append("LOW_BASE_QUALS");
      }

      if (pass){
	std::string aln_key = trim_alignment_name(alignment);
	auto aln_iter = potential_mates.find(aln_key);
	if (aln_iter != potential_mates.end()){
	  bool add = true;
	  if (check_unique_mapping_){
	    std::vector< std::pair<std::string, int32_t> > p_1, p_2;
	    get_valid_pairings(alignment, aln_iter->second, ref_vector, p_1, p_2);
	    if (p_1.size() != 1 || p_1[0].second != alignment.Position){
	      unique_mapping++;
	      add = false;
	      filter.append("NO_UNIQUE_MAPPING");
	    }
	  }

	  if (add){
	    paired_str_alns.push_back(alignment);
	    mate_alns.push_back(aln_iter->second);
	    region_alignments.push_back(alignment);
	  }
	  else {
	    assert(!filter.empty());
	    filtered_alignments.push_back(alignment);
	    if(!filtered_alignments.back().AddTag(FILTER_TAG_NAME, FILTER_TAG_TYPE, filter))
	      printErrorAndDie("Failed to add filter tag to alignment");
	  }
	  potential_mates.erase(aln_iter);
	}
	else
	  potential_strs.insert(std::pair<std::string, BamTools::BamAlignment>(aln_key, alignment));
      }
      else {
	assert(!filter.empty());
	filtered_alignments.push_back(alignment);
	if(!filtered_alignments.back().AddTag(FILTER_TAG_NAME, FILTER_TAG_TYPE, filter))
	  printErrorAndDie("Failed to add filter tag to alignment");
      }
    }
    else {
      std::string aln_key = trim_alignment_name(alignment);
      auto aln_iter = potential_strs.find(aln_key);
      if (aln_iter != potential_strs.end()){
	bool add = true;
	if (check_unique_mapping_){
	  std::vector< std::pair<std::string, int32_t> > p_1, p_2;
          get_valid_pairings(alignment, aln_iter->second, ref_vector, p_1, p_2);
          if (p_2.size() != 1 || p_2[0].second != aln_iter->second.Position){
            add = false;
            unique_mapping++;
          }
	}

	if (add){
	  paired_str_alns.push_back(aln_iter->second);
	  mate_alns.push_back(alignment);
	  region_alignments.push_back(aln_iter->second);
	}
	else {
	  std::string filter = "NO_UNIQUE_MAPPING";
	  filtered_alignments.push_back(aln_iter->second);
	  if (!filtered_alignments.back().AddTag(FILTER_TAG_NAME, FILTER_TAG_TYPE, filter))
	    printErrorAndDie("Failed to add filter tag to alignment");
	}
	potential_strs.erase(aln_iter);
      }
      else {
	auto other_iter = potential_mates.find(aln_key);
	if (other_iter != potential_mates.end())
	  potential_mates.erase(other_iter);
	else
	  potential_mates.insert(std::pair<std::string, BamTools::BamAlignment>(aln_key, alignment));
      }
    }
  }

  int32_t num_filt_unpaired_reads = 0;
  for (auto aln_iter = potential_strs.begin(); aln_iter != potential_strs.end(); ++aln_iter){
    if (check_unique_mapping_){
      if (!aln_iter->second.HasTag(ALT_MAP_TAG)){
	unpaired_str_alns.push_back(aln_iter->second);
	region_alignments.push_back(aln_iter->second);
      }
      else {
	unique_mapping++;
	std::string filter = "NO_UNIQUE_MAPPING";
	filtered_alignments.push_back(aln_iter->second);
	if(!filtered_alignments.back().AddTag(FILTER_TAG_NAME, FILTER_TAG_TYPE, filter))
	  printErrorAndDie("Failed to add filter tag to alignment");
      }
    }
    else {
      if (require_paired_reads_)
	num_filt_unpaired_reads++;
      else {
	unpaired_str_alns.push_back(aln_iter->second);
	region_alignments.push_back(aln_iter->second);
      }
    }
  }
  potential_strs.clear();
  potential_mates.clear();
  logger() << "Found " << paired_str_alns.size() << " fully paired reads and " << unpaired_str_alns.size() << " unpaired reads" << std::endl;
  
  logger() << read_count << " reads overlapped region, of which "
	   << "\n\t" << diff_chrom_mate  << " had mates on a different chromosome"
	   << "\n\t" << unmapped_mate    << " had unmapped mates"
	   << "\n\t" << insert_size      << " failed the insert size filter"
	   << "\n\t" << multimapped      << " were removed due to multimapping"
	   << "\n\t" << split_alignment  << " had an SA (split alignment) BAM tag"
	   << "\n\t" << read_has_N       << " had an 'N' base call"
	   << "\n\t" << mapping_quality  << " had too low of a mapping quality"
	   << "\n\t" << hard_clip        << " had too many hard clipped bases"
	   << "\n\t" << soft_clip        << " had too many soft clipped bases"
	   << "\n\t" << not_spanning     << " did not span the STR"
	   << "\n\t" << flank_len        << " had too bps in one or more flanks"
	   << "\n\t" << bp_before_indel  << " had too few bp before the first indel"
	   << "\n\t" << end_match_window << " did not have the maximal number of end matches within the specified window"
	   << "\n\t" << num_end_matches  << " had too few bp matches along the ends"
	   << "\n\t" << low_qual_score   << " had low base quality scores";
  if (check_unique_mapping_)
    logger() << "\n\t" << unique_mapping << " did not have a unique mapping";
  if (require_paired_reads_)
    logger() << "\n\t" << num_filt_unpaired_reads << " did not have a mate pair";
  logger() << "\n" << region_alignments.size() << " PASSED ALL FILTERS" << "\n" << std::endl;
    
  // Output the reads passing all filters to a BAM file (if requested)
  if (pass_writer.IsOpen())
    modify_and_write_alns(region_alignments, rg_to_sample, *region_iter, pass_writer);

  // Output reads that overlapped the STR but were filtered to a BAM file (if requested)
  if (filt_writer.IsOpen())
    modify_and_write_alns(filtered_alignments, rg_to_sample, *region_iter, filt_writer);

  // Separate the reads based on their associated read groups
  std::map<std::string, int> rg_indices;
  for (unsigned int i = 0; i < paired_str_alns.size(); ++i){
    std::string rg = use_bam_rgs_ ? get_read_group(paired_str_alns[i], rg_to_sample): rg_to_sample[paired_str_alns[i].Filename];

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
    std::string rg = use_bam_rgs_ ? get_read_group(unpaired_str_alns[i], rg_to_sample): rg_to_sample[unpaired_str_alns[i].Filename];

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

  locus_read_filter_time_  = (clock() - locus_read_filter_time_)/CLOCKS_PER_SEC;
  total_read_filter_time_ += locus_read_filter_time_;
}

void BamProcessor::process_regions(BamTools::BamMultiReader& reader, 
				   std::string& region_file, std::string& fasta_dir,
				   std::map<std::string, std::string>& rg_to_sample, std::map<std::string, std::string>& rg_to_library,
				   BamTools::BamWriter& pass_writer, BamTools::BamWriter& filt_writer,
				   std::ostream& out, int32_t max_regions, std::string chrom){
  std::vector<Region> regions;
  readRegions(region_file, regions, max_regions, chrom, logger());
  orderRegions(regions);

  FastaReference* fasta_ref = NULL;
  if (is_file(fasta_dir)){
    fasta_ref = new FastaReference();
    log("Fasta file exists... " + fasta_dir);
    fasta_ref->open(fasta_dir);
  }

  std::string ref_seq;
  BamTools::RefVector ref_vector = reader.GetReferenceData();
  int cur_chrom_id = -1; std::string chrom_seq;
  for (auto region_iter = regions.begin(); region_iter != regions.end(); region_iter++){
    logger() << "Processing region " << region_iter->chrom() << " " << region_iter->start() << " " << region_iter->stop() << std::endl;
    int chrom_id = reader.GetReferenceID(region_iter->chrom());
    if (chrom_id == -1 && region_iter->chrom().size() > 3 && region_iter->chrom().substr(0, 3).compare("chr") == 0)
      chrom_id = reader.GetReferenceID(region_iter->chrom().substr(3));

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
      if (fasta_ref != NULL)
	chrom_seq = fasta_ref->getSequence(chrom);
      else
	readFastaFromDir(chrom+".fa", fasta_dir, chrom_seq);
      assert(chrom_seq.size() != 0);
    }

    locus_bam_seek_time_ = clock();
    if(!reader.SetRegion(chrom_id, (region_iter->start() < MAX_MATE_DIST ? 0: region_iter->start()-MAX_MATE_DIST), 
			 chrom_id, region_iter->stop() + MAX_MATE_DIST)){
      printErrorAndDie("One or more BAM files failed to set the region properly");
    }
    locus_bam_seek_time_  =  (clock() - locus_bam_seek_time_)/CLOCKS_PER_SEC;
    total_bam_seek_time_ += locus_bam_seek_time_;

    std::vector<std::string> rg_names;
    std::vector< std::vector<BamTools::BamAlignment> > paired_strs_by_rg, mate_pairs_by_rg, unpaired_strs_by_rg;
    read_and_filter_reads(reader, chrom_seq, region_iter, rg_to_sample, rg_to_library, rg_names,
			  paired_strs_by_rg, mate_pairs_by_rg, unpaired_strs_by_rg, pass_writer, filt_writer);

    if (rem_pcr_dups_)
      remove_pcr_duplicates(base_quality_, use_bam_rgs_, rg_to_library, paired_strs_by_rg, mate_pairs_by_rg, unpaired_strs_by_rg, logger());

    std::string ref_allele = get_str_ref_allele(region_iter->start(), region_iter->stop(), chrom_seq);
    process_reads(paired_strs_by_rg, mate_pairs_by_rg, unpaired_strs_by_rg, rg_names, *region_iter, ref_allele, chrom_seq, out);
  }

  if (fasta_ref != NULL)
    delete fasta_ref;
}

