#include <getopt.h>
#include <stdlib.h>
#include <time.h>

#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <sstream>
#include <vector>

#include "bgzf_streams.h"
#include "denovos/denovo_scanner.h"
#include "error.h"
#include "haplotype_tracker.h"
#include "pedigree.h"
#include "region.h"
#include "stringops.h"
#include "version.h"
#include "vcf_reader.h"

bool file_exists(const std::string& path){
  return (access(path.c_str(), F_OK) != -1);
}

void print_usage(){
  std::cerr << "Usage: PhasingChecker --fam <fam_file> --snp-vcf <phased_snps.vcf.gz> --regions <region_file.bed> --out <edit_distances.txt>" << "\n" << "\n"
    
	    << "Required parameters:" << "\n"
	    << "\t" << "--fam        <fam_file>            "  << "\t" << "FAM file containing pedigree information for samples of interest"                     << "\n"
	    << "\t" << "--snp-vcf    <phased_snps.vcf.gz>  "  << "\t" << "Bgzipped input VCF file containing phased SNP genotypes for the samples."             << "\n"
	    << "\t" << "                                   "  << "\t" << " File should be identical to --snp-vcf argument provided to HipSTR during genotyping" << "\n"
	    << "\t" << "--regions    <region_file.bed>     "  << "\t" << "BED file containing coordinates for each STR region"                                  << "\n"
	    << "\t" << "--out        <edit_distances.gz>   "  << "\t" << "Path to which edit distance file will be written"                                     << "\n" << "\n"

	    << "Other optional parameters:" << "\n"
	    << "\t" << "--help                             "  << "\t" << "Print this help message and exit"                                                     << "\n"
	    << "\n";
}
  
void parse_command_line_args(int argc, char** argv, std::string& fam_file, std::string& snp_vcf_file, std::string& region_file, std::string& output_file){
  if (argc == 1 || (argc == 2 && std::string("-h").compare(std::string(argv[1])) == 0)){
    print_usage();
    exit(0);
  }

  int print_help = 0;
  static struct option long_options[] = {
    {"fam",             required_argument, 0, 'f'},
    {"out",             required_argument, 0, 'o'},
    {"snp-vcf",         required_argument, 0, 'v'},
    {"regions",         required_argument, 0, 'r'},
    {"h",               no_argument, &print_help, 1},
    {"help",            no_argument, &print_help, 1},
    {0, 0, 0, 0}
  };

  std::string filename;
  while (true){
    int option_index = 0;
    int c = getopt_long(argc, argv, "f:o:r:v:", long_options, &option_index);
    if (c == -1)
      break;

    switch(c){
    case 0:
      break;
    case 'f':
      fam_file = std::string(optarg);
      break;
    case 'o':
      output_file = std::string(optarg);
      break;
    case 'r':
      region_file = std::string(optarg);
      break;
    case 'v':
      snp_vcf_file = std::string(optarg);
      break;
    case '?':
      printErrorAndDie("Unrecognized command line option");
      break;
    default:
      abort();
      break;
    }
  }

  if (optind < argc) {
    std::stringstream msg;
    msg << "Did not recognize the following command line arguments:" << "\n";
    while (optind < argc)
      msg << "\t" << argv[optind++] << "\n";
    msg << "Please check your command line syntax or type ./PhasingChecker --help for additional information" << "\n";
    printErrorAndDie(msg.str());
  }

  if (print_help){
    print_usage();
    exit(0);
  }
}


int main(int argc, char** argv){
  double total_time = clock();

  std::string fam_file = "", snp_vcf_file = "", region_file = "", output_file = "", chrom="";
  parse_command_line_args(argc, argv, fam_file, snp_vcf_file, region_file, output_file);

  if (fam_file.empty())
    printErrorAndDie("--fam option required");
  else if (snp_vcf_file.empty())
    printErrorAndDie("--snp-vcf option required");
  else if (region_file.empty())
    printErrorAndDie("--regions option required");
  else if (output_file.empty())
    printErrorAndDie("--out option required");

  // Check that the FAM file exists
  if (!file_exists(fam_file))
    printErrorAndDie("FAM file " + fam_file + " does not exist. Please ensure that the path provided to --fam is valid");

  // Check that the SNP VCF file exists, has a tabix index and then open it
  if (!file_exists(snp_vcf_file))
    printErrorAndDie("SNP VCF file " + snp_vcf_file + " does not exist. Please ensure that the path provided to --snp-vcf is valid");
  if (!file_exists(snp_vcf_file + ".tbi"))
    printErrorAndDie("No .tbi index found for the SNP VCF file. Please index using tabix and rerun the analysis");
  VCF::VCFReader snp_vcf(snp_vcf_file);
  
  std::ostream& logger = std::cerr;

  // Open the output file
  if (!string_ends_with(output_file, ".gz"))
    printErrorAndDie("Output file must end in .gz");
  bgzfostream output;
  output.open(output_file.c_str());

  std::set<std::string> sites_to_skip;
  
  // Determine which samples have SNP data
  std::set<std::string> snp_samples(snp_vcf.get_samples().begin(), snp_vcf.get_samples().end());
  
  // Extract the nuclear families with SNP data
  std::vector<NuclearFamily> families;
  extract_pedigree_nuclear_families(fam_file, snp_samples, families, logger);

  // Read the STR regions
  std::vector<Region> regions;
  int32_t max_regions = 1000000000;
  readRegions(region_file, max_regions, chrom, regions, logger);
  orderRegions(regions);
  
  // Build the haplotype tracker
  int32_t window_size = 500000;
  HaplotypeTracker haplotype_tracker(families, snp_vcf_file, window_size);

  // Output the header line
  output << "#CHROM\tPOS";
  for (auto family_iter = families.begin(); family_iter != families.end(); ++family_iter)
      for (auto child_iter = family_iter->get_children().begin(); child_iter != family_iter->get_children().end(); ++child_iter)
	output << "\t" << *child_iter;
  output << "\n";

  std::string prev_chrom    = "";
  int min_second_best_score = DenovoScanner::MIN_SECOND_BEST_SCORE;
  int max_best_score        = DenovoScanner::MAX_BEST_SCORE;

  // Output the edit distances for each region
  for (auto region_iter = regions.begin(); region_iter != regions.end(); region_iter++){
    if (region_iter->chrom().compare(prev_chrom) != 0){
      logger << "Processing chromosome " << region_iter->chrom() << std::endl;
      prev_chrom = region_iter->chrom();
    }

    output << region_iter->chrom() << "\t" << region_iter->start();
    haplotype_tracker.advance(region_iter->chrom(), region_iter->start(), sites_to_skip);
    for (auto family_iter = families.begin(); family_iter != families.end(); ++family_iter){
      std::string mother = family_iter->get_mother();
      std::string father = family_iter->get_father();
      bool all_pass      = true;

      for (auto child_iter = family_iter->get_children().begin(); child_iter != family_iter->get_children().end(); ++child_iter){	
	DiploidEditDistance maternal_distance = haplotype_tracker.edit_distances(*child_iter, mother);
	DiploidEditDistance paternal_distance = haplotype_tracker.edit_distances(*child_iter, father);

	int min_mat_dist, min_mat_index, second_mat_dist, second_mat_index;
	maternal_distance.min_distance(min_mat_dist, min_mat_index);
	maternal_distance.second_min_distance(second_mat_dist, second_mat_index);
	if (min_mat_dist > max_best_score || second_mat_dist < min_second_best_score)
	  all_pass = false;

	int min_pat_dist, min_pat_index, second_pat_dist, second_pat_index;
	paternal_distance.min_distance(min_pat_dist, min_pat_index);
	paternal_distance.second_min_distance(second_pat_dist, second_pat_index);
	if (min_pat_dist > max_best_score || second_pat_dist < min_second_best_score)
	  all_pass = false;

	// Ensure that one of the parental best matches involves haplotype #1 and the other involves haplotype #2
	if (min_mat_index == 0 || min_mat_index == 1){
	  if (min_pat_index != 2 && min_pat_index != 3)
	    all_pass = false;
	}
	else if (min_pat_index != 0 && min_pat_index != 1)
	  all_pass = false;
      }

      for (auto child_iter = family_iter->get_children().begin(); child_iter != family_iter->get_children().end(); ++child_iter){	
	DiploidEditDistance maternal_distance = haplotype_tracker.edit_distances(*child_iter, mother);
	DiploidEditDistance paternal_distance = haplotype_tracker.edit_distances(*child_iter, father);
	output << "\t" << (all_pass ? "PASS" : "FAIL")
	       << ":"  << maternal_distance.distance(0,0) << "," << maternal_distance.distance(0,1)
	       << ","  << maternal_distance.distance(1,0) << "," << maternal_distance.distance(1,1)
	       << ":"  << paternal_distance.distance(0,0) << "," << paternal_distance.distance(0,1)
	       << ","  << paternal_distance.distance(1,0) << "," << paternal_distance.distance(1,1);
      }
    }
    output << "\n";
    
  }
  output.close();
  total_time = (clock() - total_time)/CLOCKS_PER_SEC;
  logger << "Execution finished: Total runtime = " << total_time << " sec" << std::endl;
  
  return 0;  
}
