#include "../bamtools/include/api/BamAlignment.h"
#include "../bamtools/include/api/BamMultiReader.h"
#include "../bamtools/include/api/BamWriter.h"

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>

#include "../error.h"
#include "exon_info.h"
#include "../stringops.h"
#include "../SeqAlignment/AlignmentOps.h"

const std::string SCORE_TAG     = "AS";
const std::string ALT_MAP_TAG   = "XA";
const std::string SPLIT_ALN_TAG = "SA";

std::string getCigarString(std::vector<BamTools::CigarOp>& cigar_list){
  std::stringstream cigar_str;
  for (std::vector<BamTools::CigarOp>::iterator iter = cigar_list.begin(); iter != cigar_list.end(); iter++)
    cigar_str << iter->Length << iter->Type;
  return cigar_str.str();
}

std::string trim_alignment_name(BamTools::BamAlignment& aln){
  std::string aln_name = aln.Name;
  if (aln_name.size() > 2){
    if (aln_name[aln_name.size()-2] == '/')
      aln_name.resize(aln_name.size()-2);
  }
  return aln_name;
}

int32_t get_end_coord(int32_t start_coord, std::string cigar_string){
  std::stringstream ss; ss << cigar_string;
  int len; char type;
  while (ss >> len >> type){
    switch(type){
    case 'M': case 'X': case 'D':
      start_coord += len;
      break;
    case 'S': case 'H': case 'I':
      break;
    default:
      printErrorAndDie("Invalid CIGAR character encountered");
      break;
    }
  }
  return start_coord;
}


void extract_transcript_aln_info(BamTools::BamAlignment& aln, ExonInfo& exon_info, BamTools::RefVector& ref_vector,
				 std::vector<std::string>& transcripts, std::vector<int32_t>& starts, std::vector<std::string>& cigar_strings){
  assert(transcripts.size() == 0 && starts.size() == 0 && cigar_strings.size() == 0);
  assert(aln.RefID != -1);

  if (exon_info.is_transcript(ref_vector[aln.RefID].RefName)){
    transcripts.push_back(ref_vector[aln.RefID].RefName);
    starts.push_back(aln.Position);
    cigar_strings.push_back(getCigarString(aln.CigarData));
  }

  if (aln.HasTag(ALT_MAP_TAG)){
    std::string alt_info;
    if (!aln.GetTag(ALT_MAP_TAG, alt_info))
      printErrorAndDie("Failed to extract XA tag from BAM alignment");
    std::vector<std::string> alts;
    split_by_delim(alt_info, ';', alts);
    for (unsigned int i = 0; i < alts.size(); i++){
      std::vector<std::string> tokens;
      split_by_delim(alts[i], ',', tokens);
      if (exon_info.is_transcript(tokens[0])){
	transcripts.push_back(tokens[0]);
	starts.push_back(abs(std::stol(tokens[1])));
	cigar_strings.push_back(tokens[2]);
      }
    }
  }
}


void populate_cigar_data(std::string& cigar_string, std::vector<BamTools::CigarOp>& cigar_data){
  assert(cigar_data.size() == 0);
  std::istringstream iss(cigar_string);
  char Type;
  uint32_t Length;
  while (iss >> Length >> Type){
    cigar_data.push_back(BamTools::CigarOp(Type, Length));
  }
}

bool is_hard_clipped(BamTools::BamAlignment& aln){
  for (auto iter = aln.CigarData.begin(); iter != aln.CigarData.end(); iter++)
    if (iter->Type == 'H')
      return true;
  return false;
}

bool is_soft_clipped(BamTools::BamAlignment& aln){
  for (auto iter = aln.CigarData.begin(); iter != aln.CigarData.end(); iter++)
    if (iter->Type == 'S')
      return true;
  return false;
}

int calc_num_bases(std::vector<BamTools::CigarOp>& cigar_ops){
  int num_bases = 0;
  for (auto iter = cigar_ops.begin(); iter != cigar_ops.end(); iter++){
    switch (iter->Type){
    case 'M': case '=': case 'X': case 'I': case'S':
      num_bases += iter->Length;
      break;
    case 'H': case 'D':
      break;
    default:
      printErrorAndDie("Invalid CIGAR char encountered");
      break;
    }
  }
  return num_bases;
}

void build_transcript_alns(BamTools::BamAlignment& aln_1, BamTools::BamAlignment& aln_2, ExonInfo& exon_info, 
			   BamTools::RefVector& ref_vector, std::map<std::string, int32_t>& ref_indexes,
			   std::vector<BamTools::BamAlignment>& alns_1_out, std::vector<BamTools::BamAlignment>& alns_2_out){
  assert(alns_1_out.size() == 0 && alns_2_out.size() == 0);

  std::vector<std::string> transcripts_1, transcripts_2;
  std::vector<int32_t> starts_1, starts_2;
  std::vector<std::string> cigar_strings_1, cigar_strings_2;
  extract_transcript_aln_info(aln_1, exon_info, ref_vector, transcripts_1, starts_1, cigar_strings_1);
  extract_transcript_aln_info(aln_2, exon_info, ref_vector, transcripts_2, starts_2, cigar_strings_2);

  std::set<std::string> bad_transcripts;
  std::map<std::string, std::pair<int, int32_t> > trans_to_pos;
  for (unsigned int i = 0; i < transcripts_1.size(); i++){
    if (bad_transcripts.find(transcripts_1[i]) != bad_transcripts.end())
      continue;
    if (trans_to_pos.find(transcripts_1[i]) != trans_to_pos.end()){
      trans_to_pos.erase(transcripts_1[i]);
      bad_transcripts.insert(transcripts_1[i]);
      continue;
    }
    trans_to_pos[transcripts_1[i]] = std::pair<int, int32_t>(i, starts_1[i]);
  }
  
  for (unsigned int i = 0; i < transcripts_2.size(); i++){
    auto t1_iter = trans_to_pos.find(transcripts_2[i]);
    if (t1_iter != trans_to_pos.end() && abs(t1_iter->second.second - starts_2[i]) < 5000){
      int t1_index = t1_iter->second.first;
      int t2_index = i;

      int32_t template_start = std::min(starts_1[t1_index], starts_2[t2_index]);
      int32_t template_end   = std::max(get_end_coord(starts_1[t1_index], cigar_strings_1[t1_index]), 
					get_end_coord(starts_2[t2_index], cigar_strings_2[t2_index])); 

      // Process valid pair of reads for transcript
      BamTools::BamAlignment out_1, out_2;
      out_1.Name         = aln_1.Name;
      out_1.Length       = aln_1.Length;
      out_1.QueryBases   = aln_1.QueryBases;
      out_1.Qualities    = aln_1.Qualities;
      out_1.RefID        = ref_indexes[transcripts_1[t1_index]];
      out_1.Position     = starts_1[t1_index];
      out_1.MapQuality   = 255;
      out_1.MateRefID    = ref_indexes[transcripts_2[t2_index]];
      out_1.MatePosition = starts_2[t2_index];
      out_1.InsertSize   = (starts_1[t1_index] <= starts_2[t2_index] ? template_end-template_start : template_start-template_end);
      populate_cigar_data(cigar_strings_1[t1_index], out_1.CigarData);
      
      out_2.Name         = aln_2.Name;
      out_2.Length       = aln_2.Length;
      out_2.QueryBases   = aln_2.QueryBases;
      out_2.Qualities    = aln_2.Qualities;
      out_2.RefID        = ref_indexes[transcripts_2[t2_index]];
      out_2.Position     = starts_2[t2_index];
      out_2.MapQuality   = 255;
      out_2.MateRefID    = ref_indexes[transcripts_1[t1_index]];
      out_2.MatePosition = starts_1[t1_index];
      out_2.InsertSize   = (starts_2[t2_index] < starts_1[t1_index] ? template_end-template_start: template_start-template_end);
      populate_cigar_data(cigar_strings_2[t2_index], out_2.CigarData);

      // Make sure that the number of bases in the CIGAR string matches the number of the bases in the reads
      if (out_1.QueryBases.size() != calc_num_bases(out_1.CigarData)){
	std::cerr << out_1.QueryBases.size() << " " << calc_num_bases(out_1.CigarData) << std::endl
		  << getCigarString(aln_1.CigarData) << std::endl
		  << aln_1.QueryBases << std::endl;

	if (aln_1.HasTag(ALT_MAP_TAG)){
	  std::string alt_info;
	  aln_1.GetTag(ALT_MAP_TAG, alt_info);
	  std::cerr << alt_info << std::endl;	
	}
      }
      assert(out_1.QueryBases.size() == calc_num_bases(out_1.CigarData));
      assert(out_2.QueryBases.size() == calc_num_bases(out_2.CigarData));
      
      alns_1_out.push_back(out_1);
      alns_2_out.push_back(out_2);
    }
  }
}


void build_transcript_alns(BamTools::BamAlignment& aln, ExonInfo& exon_info, BamTools::RefVector& ref_vector, std::map<std::string, int32_t>& ref_indexes,
			   std::vector<BamTools::BamAlignment>& alns_out){
  assert(alns_out.size() == 0);  
  std::vector<std::string> transcripts;
  std::vector<int32_t> starts;
  std::vector<std::string> cigar_strings;
  extract_transcript_aln_info(aln, exon_info, ref_vector, transcripts, starts, cigar_strings);
  for (unsigned int i = 0; i < transcripts.size(); i++){
    // Write valid read for transcript
    BamTools::BamAlignment out;
    out.Name         = aln.Name;
    out.Length       = aln.Length;
    out.QueryBases   = aln.QueryBases;
    out.Qualities    = aln.Qualities;
    out.RefID        = ref_indexes[transcripts[i]];
    out.Position     = starts[i];
    out.MapQuality   = 255;
    out.MateRefID    = -1;
    out.MatePosition = 0;
    out.InsertSize   = 0; // Since it's single ended
    populate_cigar_data(cigar_strings[i], out.CigarData);
    // Filename = aln.Filename;
    //, AlignedBases(other.AlignedBases)
    // , TagData(other.TagData)
    // , Bin(other.Bin)
    //  , AlignmentFlag(other.AlignmentFlag)
    //, SupportData(other.SupportData)
    
    assert(out.QueryBases.size() == calc_num_bases(out.CigarData));
    alns_out.push_back(out);
  }
}


void extract_gene_clusters(BamTools::BamAlignment& aln, ExonInfo& exon_info, 
			   BamTools::RefVector& ref_vector, std::map<std::string, int>& neg_chrom_indexes, 
			   std::set<int>& gene_clusters){
  assert(gene_clusters.size() == 0);
  
  // Unmapped read
  if (aln.RefID == -1)
    return;
  if (aln.CigarData.size() == 0)
    return;

  bool has_trans = false;

  std::vector<std::string> ref_names;
  ref_names.push_back(ref_vector[aln.RefID].RefName);
  if (exon_info.is_transcript(ref_names.back())){
    gene_clusters.insert(exon_info.gene_cluster_for_transcript(ref_names.back()));
    has_trans = true;
  }
  else {
    // Check to see if region intersects gene body/exons using interval tree
    std::vector<int> cluster_ids;
    exon_info.gene_clusters_from_region(ref_names.back(), aln.Position, aln.GetEndPosition(), cluster_ids);
    if (cluster_ids.size() != 0)
      for (auto iter = cluster_ids.begin(); iter != cluster_ids.end(); iter++)
	gene_clusters.insert(*iter);
    else
      gene_clusters.insert(neg_chrom_indexes[ref_names.back()]);
  }

  if (aln.HasTag(ALT_MAP_TAG)){
    std::string alt_info;
    if (!aln.GetTag(ALT_MAP_TAG, alt_info))
      printErrorAndDie("Failed to extract XA tag from BAM alignment");
    std::vector<std::string> alts;
    split_by_delim(alt_info, ';', alts);
    for (unsigned int i = 0; i < alts.size(); i++){
      std::vector<std::string> tokens;
      split_by_delim(alts[i], ',', tokens);
      assert(tokens.size() == 4);      
      ref_names.push_back(tokens[0]);
      if (exon_info.is_transcript(ref_names.back())){
	gene_clusters.insert(exon_info.gene_cluster_for_transcript(ref_names.back()));
	has_trans = true;
      }
      else {
	int32_t start = std::stol(tokens[1]);
	if (start < 0) start = -start;
	int32_t end   = get_end_coord(start, tokens[2]);

	// Check to see if region intersects gene body/exons using interval tree
	std::vector<int> cluster_ids;
	exon_info.gene_clusters_from_region(ref_names.back(), start, end, cluster_ids);
	if (cluster_ids.size() != 0)
	  for (auto iter = cluster_ids.begin(); iter != cluster_ids.end(); iter++)
	    gene_clusters.insert(*iter);
	else
	  gene_clusters.insert(neg_chrom_indexes[ref_names.back()]);
      }
    }
  }
  
  // A transcript wasn't listed in the primary alignment or any of the alternates
  // If one alignment overlapped with the gene body, we'll still assign it to the gene cluster
  // We should instead remove all alignments
  if (!has_trans){
    bool gt_zero = false;
    for (auto iter = gene_clusters.begin(); iter != gene_clusters.end(); iter++)
      gt_zero |= (*iter > 0);
    if (gt_zero)
      gene_clusters.clear();
  }
} 


/*
 * Filter BAMs in which paired end reads are adjacent to one another in the BAM file.
 */
void filter_bam_paired_mode(BamTools::BamReader& reader,
			    ExonInfo& exon_info,
			    std::string& keep_trans_filename, std::string& keep_dna_filename, std::string& filt_filename){
  BamTools::RefVector ref_vector = reader.GetReferenceData();

  // Index the chromosome names starting at -2
  std::map<std::string, int> neg_chrom_indexes;
  neg_chrom_indexes["*"] = -1;
  int chrom_index = -2;
  for (unsigned int i = 0; i < ref_vector.size(); ++i, --chrom_index)
    neg_chrom_indexes[ref_vector[i].RefName] = chrom_index;

  // Index sequences in BAM by id
  std::map<std::string, int32_t> ref_indexes;
  for (int32_t i = 0; i < ref_vector.size(); i++)
    ref_indexes[ref_vector[i].RefName] = i;

  BamTools::BamWriter keep_dna_writer, keep_trans_writer, filt_writer;
  if (!keep_dna_writer.Open(keep_dna_filename, reader.GetHeaderText(), ref_vector)) 
    printErrorAndDie("Failed to open output BAM file for DNA reads to keep");
  if (!keep_trans_writer.Open(keep_trans_filename, reader.GetHeaderText(), ref_vector)) 
    printErrorAndDie("Failed to open output BAM file for transcript reads to keep");
  if (!filt_writer.Open(filt_filename, reader.GetHeaderText(), ref_vector)) 
    printErrorAndDie("Failed to open output BAM file for reads to filter");
  BamTools::BamAlignment alignments[2];
   
  int32_t read_count      = 0;
  int aln_count           = 0;
  int32_t pairs_saved     = 0;
  int32_t pairs_removed   = 0;
  int32_t singles_saved   = 0;
  int32_t singles_removed = 0;
  int32_t pair_transcript = 0;
  int32_t pair_dna        = 0;
 
  while (reader.GetNextAlignment(alignments[aln_count])){
    read_count++;
    if (aln_count == 0)
      aln_count++;
    else {
      std::string key_a = trim_alignment_name(alignments[0]);
      std::string key_b = trim_alignment_name(alignments[1]);

      if (key_a.compare(key_b) == 0){
	// Process paired end read
	bool save_alns = true;

	// Check that the optimal alignments for both reads map to the same gene
	std::set<int> aln_a_ids, aln_b_ids;
	std::vector<int> lst_intersection_ids;
	extract_gene_clusters(alignments[0], exon_info, ref_vector, neg_chrom_indexes, aln_a_ids);
	extract_gene_clusters(alignments[1], exon_info, ref_vector, neg_chrom_indexes, aln_b_ids);
	std::set_intersection(aln_a_ids.begin(), aln_a_ids.end(), aln_b_ids.begin(), aln_b_ids.end(), std::back_inserter(lst_intersection_ids));
	std::set<int> intersection_ids(lst_intersection_ids.begin(), lst_intersection_ids.end());
	
	if (intersection_ids.size() == 1 && *(intersection_ids.begin()) != -1){
	  if (*intersection_ids.begin() > 0)
	    save_alns = true;
	  else {
	    // Check that there's only 1 DNA location and that the mate pairs are within a certain distance
	    save_alns =  (!alignments[0].HasTag(ALT_MAP_TAG) && !alignments[0].HasTag(ALT_MAP_TAG) && 
			  std::abs(alignments[0].Position - alignments[1].Position) < 5000);
	  }
	}
	else
	  save_alns = false;

	// Alignments with this tag lead to big problems, as the primary alignment might have a hard clip
	// where the secondary ones have soft clips. But when we're trying to build the secondary alignments,
	// the bases are gone
	if (alignments[0].HasTag(SPLIT_ALN_TAG) || alignments[1].HasTag(SPLIT_ALN_TAG))
	  save_alns = false;


	/*
	std::cerr << "A IDs" << std::endl;
	for (auto i = aln_a_ids.begin(); i != aln_a_ids.end(); i++)
	  exon_info.print_cluster_info(*i, std::cerr);
	std::cerr << std::endl;
	std::cerr << "B IDs" << std::endl;
	for (auto i = aln_b_ids.begin(); i != aln_b_ids.end(); i++)
	  exon_info.print_cluster_info(*i, std::cerr);
	std::cerr << std::endl;
	*/

	if (save_alns){
	  if (*(intersection_ids.begin()) < 0){
	    pair_dna++;
	    if (!keep_dna_writer.SaveAlignment(alignments[0])) printErrorAndDie("Failed to save alignment");   
	    if (!keep_dna_writer.SaveAlignment(alignments[1])) printErrorAndDie("Failed to save alignment");
	  }
	  else {
	    std::vector<BamTools::BamAlignment> alns_1_out, alns_2_out;
	    build_transcript_alns(alignments[0], alignments[1], exon_info, ref_vector, ref_indexes, alns_1_out, alns_2_out);
	    pair_transcript++;
	    for (unsigned int i = 0; i < alns_1_out.size(); i++){
	      if (!keep_trans_writer.SaveAlignment(alns_1_out[i])) printErrorAndDie("Failed to save alignment");   
	      if (!keep_trans_writer.SaveAlignment(alns_2_out[i])) printErrorAndDie("Failed to save alignment");   
	    }
	  }
	  pairs_saved++;
	}
	else {
	  pairs_removed++;
	  if (!filt_writer.SaveAlignment(alignments[0])) printErrorAndDie("Failed to save alignment");   
	  if (!filt_writer.SaveAlignment(alignments[1])) printErrorAndDie("Failed to save alignment");
	}
	aln_count = 0;
      }
      else {
	// Process single read
	bool save_aln = true;

	// Check that the alignment for the read is unique
	std::set<int> aln_a_ids;
	extract_gene_clusters(alignments[0], exon_info, ref_vector, neg_chrom_indexes, aln_a_ids);
	save_aln = false;

	/*
	for (auto i = aln_a_ids.begin(); i != aln_a_ids.end(); i++)
	  exon_info.print_cluster_info(*i, std::cerr);
	std::cerr << *(aln_a_ids.begin()) << " " << (*aln_a_ids.begin() > 0) << " " << *aln_a_ids.begin() << std::endl;
	printErrorAndDie("boo");
	*/

	if (aln_a_ids.size() == 1 && *(aln_a_ids.begin()) != -1){
	  if (*aln_a_ids.begin() > 0)
	    save_aln = true;
	  else
	    save_aln = (!alignments[0].HasTag(ALT_MAP_TAG)); // Check that there's only 1 DNA location and that the mate pairs are within a certain distance
	}
	else
	  save_aln = false;
	
	if (alignments[0].HasTag(SPLIT_ALN_TAG))
	  save_aln = false;

	if (save_aln){
	  singles_saved++;
	  if ((*aln_a_ids.begin() > 0)){
	    std::vector<BamTools::BamAlignment> alns_out;
	    build_transcript_alns(alignments[0], exon_info, ref_vector, ref_indexes, alns_out);	   
	    for (unsigned int i = 0; i < alns_out.size(); i++)
	      if (!keep_trans_writer.SaveAlignment(alns_out[i])) printErrorAndDie("Failed to save alignment");   
	  }
	  else if(!keep_dna_writer.SaveAlignment(alignments[0])) 
	    printErrorAndDie("Failed to save alignment");
	}
	else {
	  singles_removed++;
	  if(!filt_writer.SaveAlignment(alignments[0])) printErrorAndDie("Failed to save alignment");
	}

	// Move second read to first position
	alignments[0] = alignments[1];
      }
    }

    if (read_count % 1000000 == 0)
      std::cerr << "\tProcessing read # "    << read_count      << "\t" 
		<< "PAIRED(#filt, #kept) = " << pairs_removed   << " " << pairs_saved   << "\t"
		<< "SINGLE(#filt, #kept) = " << singles_removed << " " << singles_saved << std::endl;
  }
  keep_dna_writer.Close();
  keep_trans_writer.Close();
  filt_writer.Close();
  std::cerr << "\tProcessing read # "    << read_count      << "\t" 
	    << "PAIRED(#filt, #kept) = " << pairs_removed   << " " << pairs_saved      << "\n"
	    << "            #DNA = " << pair_dna << " #TRANSCRIPT = " << pair_transcript << " #FILTER = " << pairs_removed << "\n"
	    << "SINGLE(#filt, #kept) = " << singles_removed << " " << singles_saved << std::endl;
}


int main(int argc, char** argv){
  ExonInfo exon_info(argv[1]);
  std::string input_file          = argv[2];
  std::string keep_dna_filename   = argv[3];
  std::string keep_trans_filename = argv[4];
  std::string filt_filename       = argv[5];
  BamTools::BamReader reader;
  if (!reader.Open(input_file)) printErrorAndDie("Failed to open BAM file");
  filter_bam_paired_mode(reader, exon_info, keep_trans_filename, keep_dna_filename, filt_filename);
}

