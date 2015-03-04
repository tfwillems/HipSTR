#include <climits>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>

#include "bamtools/include/api/BamAlignment.h"
#include "bamtools/include/api/BamMultiReader.h"
#include "bamtools/include/api/BamWriter.h"

#include "alignment_filters.h"
#include "error.h"
#include "extract_indels.h"
#include "filter_bams.h"
#include "region.h"
#include "seqio.h"
#include "stringops.h"

int  MAX_MATE_DIST            = 1000;
int  MIN_BP_BEFORE_INDEL      = 7;
int  MIN_FLANK                = 5;
int  MIN_READ_END_MATCH       = 10;
int  MAXIMAL_END_MATCH_WINDOW = 15;
int  REMOVE_MULTIMAPPERS      = 0;

void printCounts(std::map<int,int>& counts, std::ostream& out){
  auto iter = counts.begin();
  while (iter != counts.end()){
    if (std::next(iter) == counts.end())
      out << iter->first << ":" << iter->second;
    else
      out << iter->first << ":" << iter->second << ",";
    iter++;
  }
}

void analyzeStutter(std::vector< std::vector<BamTools::BamAlignment> >& alignment_lists,
		    std::vector<std::string>& rg_names, Region& region,
		    std::ofstream& stutter_output){
  stutter_output << region.chrom() << "\t" << region.start() << "\t" << region.stop() << "\t" << region.period() << "\t";

  // Process each set of reads separately
  for (int i = 0; i < alignment_lists.size(); i++){
    // Determine the indel size associated with each read
    std::map<int, int> indel_counts;
    int fail_count = 0;
    for (int j = 0; j < alignment_lists[i].size(); j++){
      BamTools::BamAlignment& alignment = alignment_lists[i][j];
      int bp_diff;
      bool got_size = ExtractCigar(alignment.CigarData, alignment.Position, region.start(), region.stop(), bp_diff);
      if (got_size) 
	indel_counts[bp_diff]++;
      else
	fail_count++;
    }
    std::cerr << "RG " << rg_names[i] << "\t"
	      << "Fail count = " << fail_count << "\t"
	      << "Bp counts" << "\t";
    printCounts(indel_counts, std::cerr);
    std::cerr << "\n";
    
    stutter_output << rg_names[i] << "\t";
    printCounts(indel_counts, stutter_output);
    stutter_output << "\t";
  }
  stutter_output << "\n";
  std::cerr << std::endl;
}


void processRegions(BamTools::BamMultiReader& reader, 
		    std::string& region_file, std::string& fasta_dir,
		    std::map<std::string, std::string>& file_read_groups,
		    std::ofstream& stutter_output, BamTools::BamWriter& bam_writer){
  std::vector<Region> regions;
  readRegions(region_file, regions);
  
  std::string ref_seq;
  BamTools::BamAlignment alignment;
  BamTools::RefVector ref_vector = reader.GetReferenceData();
  int32_t str_start, str_stop;
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

    if(!reader.SetRegion(chrom_id, region_iter->start(), chrom_id, region_iter->stop()))
      printErrorAndDie("One or more BAM files failed to set the region properly");
   
    std::vector<BamTools::BamAlignment> region_alignments;
    int read_count = 0;
    int diff_chrom_mate = 0, unmapped_mate = 0, not_spanning = 0; // Counts for filters that are always applied
    int insert_size = 0, multimapped = 0, flank_len = 0, bp_before_indel = 0, end_match_window = 0, num_end_matches = 0; // Counts for filters that are user-controlled
    while (reader.GetNextAlignment(alignment)){
      read_count++;

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
      // Ignore read if it does not span the STR
      if (alignment.Position > region_iter->start() || alignment.GetEndPosition() < region_iter->stop()){
	not_spanning++;
	continue;
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
      // Ignore read if it has insufficient flanking bases on either side of the STR
      if (alignment.Position > (region_iter->start()-MIN_FLANK) || alignment.GetEndPosition() < (region_iter->stop()+MIN_FLANK)){
	flank_len++;
	continue;
      }
      // Ignore read if there is an indel within the first MIN_BP_BEFORE_INDEL bps from each end
      if (MIN_BP_BEFORE_INDEL > 0){
	std::pair<int, int> num_bps = AlignmentFilters::GetEndDistToIndel(alignment);
	if ((num_bps.first != -1 && num_bps.first < MIN_BP_BEFORE_INDEL) || (num_bps.second != -1 && num_bps.second < MIN_BP_BEFORE_INDEL)){
	  bp_before_indel++;
	  continue;
	}
      }
      // Ignore read if there is another location within MAXIMAL_END_MATCH_WINDOW bp for which it has a longer end match
      if (MAXIMAL_END_MATCH_WINDOW > 0){
	bool maximum_end_matches = AlignmentFilters::HasLargestEndMatches(alignment, chrom_seq, 0, MAXIMAL_END_MATCH_WINDOW, MAXIMAL_END_MATCH_WINDOW);
	if (!maximum_end_matches){
	  end_match_window++;
	  continue;
	}
      }
      // Ignore read if it doesn't match perfectly for at least MIN_READ_END_MATCH bases on each end
      if (MIN_READ_END_MATCH > 0){
	std::pair<int,int> match_lens = AlignmentFilters::GetNumEndMatches(alignment, chrom_seq, 0);
	if (match_lens.first < MIN_READ_END_MATCH || match_lens.second < MIN_READ_END_MATCH){ 
	  num_end_matches++;
	  continue;
	}
      }
      region_alignments.push_back(alignment);
    }
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
	std::string rg_tag = "lobSTR;" + file_read_groups[read_iter->Filename] + ";" + file_read_groups[read_iter->Filename];
	read_iter->AddTag("RG", "Z", rg_tag);

	// Add STR start and stop tags
	bool success;
	if (read_iter->HasTag("XS"))
	  read_iter->RemoveTag("XS");
	if(!read_iter->AddTag("XS",  "I", region_iter->start()))
	  printErrorAndDie("Failed to modify XS tag");
	if (read_iter->HasTag("XE"))
	  read_iter->RemoveTag("XE");
	if(!read_iter->EditTag("XE", "I", region_iter->stop()))
	  printErrorAndDie("Failed to modify XE tag");

	if (!bam_writer.SaveAlignment(*read_iter))  
	  printErrorAndDie("Failed to save alignment for STR-spanning read");
      }
    }

    // Separate the reads based on their associated read groups
    std::map<std::string, int> rg_indices;
    std::vector< std::vector<BamTools::BamAlignment> > alignments_by_rg;
    std::vector<std::string> rg_names;
    for (auto align_iter = region_alignments.begin(); align_iter != region_alignments.end(); align_iter++){
      std::string rg = file_read_groups[align_iter->Filename];
      int rg_index;
      auto index_iter = rg_indices.find(rg);
      if (index_iter == rg_indices.end()){
	rg_index = rg_indices.size(); 
	rg_indices[rg] = rg_index;
	rg_names.push_back(rg);
	alignments_by_rg.push_back(std::vector<BamTools::BamAlignment>());
      }
      else
	rg_index = index_iter->second;
      alignments_by_rg[rg_index].push_back(*align_iter);
    }

    // Analze the PCR stutter levels
    analyzeStutter(alignments_by_rg, rg_names, *region_iter, stutter_output);
  } 
}

void parse_command_line_args(int argc, char** argv, 
			     std::string& bamfile_string, std::string& bamindex_string, std::string& rg_string,
			     std::string& fasta_dir, std::string& region_file,  std::string& stutter_out_file,
			     std::string& bam_out_file){
   if (argc == 1){
     std::cerr << "Usage: StutterTrainer --bams <list_of_bams> --indexes <list_of_bam_indexes> --rgs <list_of_read_groups>" << "\n"
	       << "                      --stutter-out <stutter.txt> --regions <region_file.bed> --fasta <dir> [--bam-out <spanning_reads.bam>] [--rem-multimaps]"   << "\n\n"
	       << "\t" << "--bams          <list_of_bams>        "  << "\t" << "Comma separated list of .bam files"                                                  << "\n"
	       << "\t" << "--indexes       <list_of_bam_indexes> "  << "\t" << "Comma separated list of .bai files in same order as .bam files"                      << "\n"
	       << "\t" << "--rgs           <list_of_read_groups> "  << "\t" << "Comma separated list of read groups in same order as .bam files"                     << "\n" 
	       << "\t" << "--fasta         <dir>                 "  << "\t" << "Directory in which FASTA files for each chromosome are located"                      << "\n"
	       << "\t" << "--stutter-out   <stutter.txt>         "  << "\t" << "File to which the stutter models will be written"                                    << "\n"
	       << "\t" << "--regions       <region_file.bed>     "  << "\t" << "BED file containing coordinates for regions to render "                              << "\n"
	       << "\t" << "--bam-out       <spanning_reads.bam   "  << "\t" << "Output a BAM file containing the reads spanning each region to the provided file"    << "\n"
	       << "\t" << "--rem-multimaps                       "  << "\t" << "Remove reads that map to multiple locations (Default = False)"                       << "\n"
	       << "\t" << "--max-mate-dist <max_bp>              "  << "\t" << "Remove reads whose mate pair distance is > MAX_BP (Default = " << MAX_MATE_DIST << ")" << "\n"
	       << "\n";
     exit(0);
  }
 
  static struct option long_options[] = {
    {"bams",            required_argument, 0, 'b'},
    {"max-mate-dist",   required_argument, 0, 'd'},
    {"fasta",           required_argument, 0, 'f'},
    {"rgs",             required_argument, 0, 'g'},
    {"indexes",         required_argument, 0, 'i'},
    {"stutter-out",     required_argument, 0, 'o'},
    {"regions",         required_argument, 0, 'r'},
    {"bam-out",         required_argument, 0, 'w'},
    {"rem-multimaps",   no_argument, &REMOVE_MULTIMAPPERS, 1},
    {0, 0, 0, 0}
  };

  int c;
  while (true){
    int option_index = 0;
    c = getopt_long(argc, argv, "nb:f:i:m:o:r:v:", long_options, &option_index);
    if (c == -1)
      break;

    switch(c){
    case 0:
      break;
    case 'b':
      bamfile_string = std::string(optarg);
      break;
    case 'd':
      MAX_MATE_DIST = atoi(optarg);
      break;
    case 'f':
      fasta_dir = std::string(optarg);
      break;
    case 'g':
      rg_string = std::string(optarg);
      break;
    case 'i':
      bamindex_string = std::string(optarg);
      break;
    case 'o':
      stutter_out_file = std::string(optarg);
      break;
    case 'r':
      region_file = std::string(optarg);
      break;
    case 'w':
      bam_out_file = std::string(optarg);
      break;
    case '?':
      printErrorAndDie("Unrecognized command line option");
      break;
    default:
      abort();
      break;
    }
  }
}

int main(int argc, char** argv){
  std::string bamfile_string= "", bamindex_string="", rg_string="", region_file="", stutter_out_file="", fasta_dir="", bam_out_file="";
  parse_command_line_args(argc, argv, bamfile_string, bamindex_string, rg_string, fasta_dir, 
			  region_file, stutter_out_file, bam_out_file);
  int num_flank = 0;
  if (bamfile_string.empty())
    printErrorAndDie("--bams option required");
  else if (bamindex_string.empty())
    printErrorAndDie("--indexes option required");
  else if (rg_string.empty())
    printErrorAndDie("--rgs option required");
  else if (stutter_out_file.empty())
    printErrorAndDie("--out option required");
  else if (region_file.empty())
    printErrorAndDie("--region option required");
  else if (fasta_dir.empty())
    printErrorAndDie("--fasta option required");

  if (fasta_dir.back() != '/')
    fasta_dir += "/";
  std::cerr << "--bams         " << bamfile_string   << "\n"
	    << "--indexes      " << bamindex_string  << "\n"
	    << "--rgs          " << rg_string        << "\n"
	    << "--fasta        " << fasta_dir        << "\n"
	    << "--stutter-out  " << stutter_out_file << "\n"
	    << "--regions      " << region_file      << "\n" << std::endl;

  std::vector<std::string> bam_files;
  split_by_delim(bamfile_string, ',', bam_files);
  std::vector<std::string> bam_indexes;
  split_by_delim(bamindex_string, ',', bam_indexes);
  std::vector<std::string> read_groups;
  split_by_delim(rg_string, ',', read_groups);
  std::cerr << "Detected " << bam_files.size() << " BAM files" << std::endl;
  if (bam_files.size() != bam_indexes.size())
    printErrorAndDie("Number of .bam and .bai files must match");

  // Open all .bam and .bai files
  BamTools::BamMultiReader reader;
  if (!reader.Open(bam_files)) printErrorAndDie("Failed to open one or more BAM files");
  if (!reader.OpenIndexes(bam_indexes)) printErrorAndDie("Failed to open one or more BAM index files");

  // Open output
  std::ofstream stutter_output(stutter_out_file.c_str());
  if (!stutter_output.is_open()) printErrorAndDie("Failed to open output file");

  // Construct filename->read group map
  std::map<std::string, std::string> file_read_groups;
  for (int i = 0; i < bam_files.size(); i++)
    file_read_groups[bam_files[i]] = read_groups[i];

  BamTools::BamWriter bam_writer;
  if (!bam_out_file.empty()){
    BamTools::RefVector ref_vector = reader.GetReferenceData();
    bool file_open = bam_writer.Open(bam_out_file, reader.GetHeaderText(), ref_vector);
    if (!file_open) printErrorAndDie("Failed to open output BAM file");
  }

  // Run analysis
  processRegions(reader, region_file, fasta_dir, file_read_groups, stutter_output, bam_writer);

  if (!bam_out_file.empty()) bam_writer.Close();
  stutter_output.close();
  reader.Close();
  return 0;  
}
