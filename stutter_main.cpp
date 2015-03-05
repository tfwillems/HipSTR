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

#include "bam_processor.h"
#include "error.h"
#include "extract_indels.h"
#include "region.h"
#include "stringops.h"

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

class StutterBamProcessor : public BamProcessor {
public:
  void process_reads(std::vector< std::vector<BamTools::BamAlignment> >& alignments_by_rg, std::vector<std::string>& rg_names, Region& region, std::ostream& out){
    out << region.chrom() << "\t" << region.start() << "\t" << region.stop() << "\t" << region.period() << "\t";
    
    // Process each set of reads separately
    for (unsigned int i = 0; i < alignments_by_rg.size(); i++){
      // Determine the indel size associated with each read
      std::map<int, int> indel_counts;
      int fail_count = 0;
      for (unsigned int j = 0; j < alignments_by_rg[i].size(); j++){
	BamTools::BamAlignment& alignment = alignments_by_rg[i][j];
	int bp_diff;
	bool got_size = ExtractCigar(alignment.CigarData, alignment.Position, region.start(), region.stop(), bp_diff);
	if (got_size) 
	  indel_counts[bp_diff]++;
	else
	  fail_count++;
      }
      std::cerr << "RG " << rg_names[i] << "\t" << "Fail count = " << fail_count << "\t" << "Bp counts" << "\t";
      printCounts(indel_counts, std::cerr);
      std::cerr << "\n";
      
      out << rg_names[i] << "\t";
      printCounts(indel_counts, out);
      out << "\t";
    }
    out << "\n";
    std::cerr << std::endl; 
  }
};

void parse_command_line_args(int argc, char** argv, 
			     std::string& bamfile_string, std::string& bamindex_string, std::string& rg_string,
			     std::string& fasta_dir, std::string& region_file,  std::string& stutter_out_file,
			     std::string& bam_out_file, BamProcessor& bam_processor){
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
	       << "\t" << "--max-mate-dist <max_bp>              "  << "\t" << "Remove reads whose mate pair distance is > MAX_BP (Default = " << bam_processor.MAX_MATE_DIST << ")" << "\n"
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
    {"rem-multimaps",   no_argument, &(bam_processor.REMOVE_MULTIMAPPERS), 1},
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
      bam_processor.MAX_MATE_DIST = atoi(optarg);
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
  StutterBamProcessor bam_processor;
  std::string bamfile_string= "", bamindex_string="", rg_string="", region_file="", stutter_out_file="", fasta_dir="", bam_out_file="";
  parse_command_line_args(argc, argv, bamfile_string, bamindex_string, rg_string, fasta_dir, 
			  region_file, stutter_out_file, bam_out_file, bam_processor);
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
  bam_processor.process_regions(reader, region_file, fasta_dir, file_read_groups, bam_writer, stutter_output);

  if (!bam_out_file.empty()) bam_writer.Close();
  stutter_output.close();
  reader.Close();
  return 0;  
}
