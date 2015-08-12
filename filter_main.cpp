#include <climits>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>

#include "bamtools/include/api/BamReader.h"

#include "error.h"
#include "filter_bams.h"
#include "insert_size.h"
#include "region.h"


void parse_command_line_args(int argc, char** argv,    std::string& input_file, std::string& output_file, 
			     std::string& region_file, std::string& insert_stats_file, int& max_diff, int& paired_mode, int& region_pad){
   if (argc == 1){
    std::cerr << "Usage: BamSieve --in <in.bam> --out <out.bam> --regions <region_file.bed> [--insert-stats <stat_file.txt>]]"                            << "\n"
	      << "\t" << "--in            <in.bam>         " << "\t"  << "Input BAM file to filter"                                                       << "\n"
	      << "\t" << "--out           <out.bam>        " << "\t"  << "Output BAM file containing filtered reads and their mate pairs"                 << "\n"
	      << "\t" << "--regions       <region_file.bed>" << "\t"  << "BED file containing coordinates for regions to filter "                         << "\n"
	      << "\t" << "--insert-stats  <stat_file.txt>  " << "\t"  << "Compute insert size stats from all reads and write them to the provided file"   << "\n"
	      << "\t" << "--max-diff      <max_size>       " << "\t"  << "Only compute insert size stats for pairs with distance < MAX_SIZE. Default = "  << max_diff << "\n"
	      << "\t" << "--paired                         " << "\t"  << "Paired end reads are adjacent in the BAM file. By default, it is assumed that"  << "\n"
	      << "\t" << "                                 " << "\t"  << "this is not the case and that the BAM is sorted by position"                    << "\n"
	      << "\t" << "--pad           <bp_pad>         " << "\t"  << "Extend each region by BP_PAD base pairs. By default, each region is unmodified" << "\n"
	      << "\n";
    exit(0);
  }
 
  static struct option long_options[] = {
    {"max-diff",     required_argument, 0, 'd'},
    {"in",           required_argument, 0, 'i'},
    {"out",          required_argument, 0, 'o'},
    {"pad",          required_argument, 0, 'p'},
    {"regions",      required_argument, 0, 'r'},
    {"insert-stats", required_argument, 0, 's'},
    {"paired",       no_argument,  &paired_mode, 1},
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
    case 'd':
      max_diff = atoi(optarg);
      break;
    case 'i':
      input_file = std::string(optarg);
      break;
    case 'o':
      output_file = std::string(optarg);
      break;
    case 'p':
      region_pad = atoi(optarg);
      break;
    case 'r':
      region_file = std::string(optarg);
      break;
    case 's':
      insert_stats_file = std::string(optarg);
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
  std::string input_file="", output_file="", region_file="", insert_stats_file="";
  int max_insert_size = 2000, paired_mode = 0, region_pad = 0;
  parse_command_line_args(argc, argv, input_file, output_file, region_file, insert_stats_file, max_insert_size, paired_mode, region_pad);
 
  if (input_file.empty())
    printErrorAndDie("--in option required");
  else if (output_file.empty())
    printErrorAndDie("--out option required");
  else if (region_file.empty())
    printErrorAndDie("--region option required");

  std::cerr << "--in      " << input_file  << "\n"
	    << "--out     " << output_file << "\n"
	    << "--regions " << region_file << "\n";
  if (region_pad != 0)
    std::cerr << "--pad     " << region_pad << "\n";
  std::cerr << std::endl;
  
  // Open the BAM file
  BamTools::BamReader reader;
  if (!reader.Open(input_file)) printErrorAndDie("Failed to open BAM file");

  // Read and arrange regions
  std::vector<Region> regions;
  readRegions(region_file, regions, -1, "", std::cerr);

  // Extend regions by padding
  if (region_pad != 0){
    for (unsigned int i = 0; i < regions.size(); i++){
      regions[i].set_start(regions[i].start() - region_pad);
      regions[i].set_stop(regions[i].stop()   + region_pad);
    }
  }

  // Sort and arrange regions
  std::vector< std::vector<Region> > ordered_regions;
  std::map<std::string, int> chrom_order;
  orderRegions(regions, ordered_regions, chrom_order);

  // Filter BAM
  std::cerr << "Filtering BAM" << std::endl;
  InsertSizeCounter insert_counter(max_insert_size);

  if (paired_mode == 0)
    filter_bam(reader, ordered_regions, chrom_order, output_file, insert_counter, !(insert_stats_file.empty()));
  else
    filter_bam_paired_mode(reader, ordered_regions, chrom_order, output_file, insert_counter, (!insert_stats_file.empty()));

  if (!insert_stats_file.empty()){
    std::ofstream stats_output (insert_stats_file, std::ofstream::out);
    insert_counter.output_summary_statistics(stats_output);
    stats_output.close();
  }
 
  reader.Close();
  return 0;  
}
