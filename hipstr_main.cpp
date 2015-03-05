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
#include "vcflib/src/Variant.h"

#include "bam_processor.h"
#include "error.h"
#include "extract_indels.h"
#include "region.h"
#include "seqio.h"
#include "snp_tree.h"
#include "stringops.h"

class SNPBamProcessor : public BamProcessor {
private:
  bool have_vcf;
  vcf::VariantCallFile phased_vcf;

public:
  SNPBamProcessor(){
    have_vcf = false;
  }

  SNPBamProcessor(std::string& vcf_file){
    set_vcf(vcf_file);
  }

  void set_vcf(std::string& vcf_file){
    if(!phased_vcf.open(vcf_file))
      printErrorAndDie("Failed to open VCF file");
    have_vcf = true;
  }

  void process_reads(std::vector< std::vector<BamTools::BamAlignment> >& alignments_by_rg, std::vector<std::string>& rg_names, Region& region, std::ostream& out){
    if (have_vcf){
      std::vector<SNPTree*> snp_trees;
      std::map<std::string, unsigned int> sample_indices;
      
      create_snp_trees(region.chrom(), (region.start() > MAX_MATE_DIST ? region.start()-MAX_MATE_DIST : 1), 
		       region.stop()+MAX_MATE_DIST, phased_vcf, sample_indices, snp_trees);
      
      destroy_snp_trees(snp_trees);
    }
    else {

    }
  }
};


void parse_command_line_args(int argc, char** argv, 
			     std::string& bamfile_string, std::string& bamindex_string, std::string& rg_string,
			     std::string& fasta_dir, std::string& region_file,  std::string& vcf_file, std::string& chrom, 
			     std::string& bam_out_file, BamProcessor& bam_processor){
   if (argc == 1){
     std::cerr << "Usage: StutterTrainer --bams <list_of_bams> --indexes <list_of_bam_indexes> --rgs <list_of_read_groups>" << "\n"
	       << "                      --stutter-out <stutter.txt> --regions <region_file.bed> --fasta <dir> [--bam-out <spanning_reads.bam>] [--rem-multimaps]"   << "\n\n"
	       << "\t" << "--bams          <list_of_bams>        "  << "\t" << "Comma separated list of .bam files"                                                  << "\n"
	       << "\t" << "--indexes       <list_of_bam_indexes> "  << "\t" << "Comma separated list of .bai files in same order as .bam files"                      << "\n"
	       << "\t" << "--rgs           <list_of_read_groups> "  << "\t" << "Comma separated list of read groups in same order as .bam files"                     << "\n" 
	       << "\t" << "--fasta         <dir>                 "  << "\t" << "Directory in which FASTA files for each chromosome are located"                      << "\n"
	       << "\t" << "--chrom         <chrom>               "  << "\t" << "Only consider STRs on the provided chromosome"                                       << "\n"
	       << "\t" << "--vcf           <phased_gts.vcf>      "  << "\t" << "VCF file containing phased genotypes"                                                << "\n"
	       << "\t" << "--regions       <region_file.bed>     "  << "\t" << "BED file containing coordinates for each STR region"                                 << "\n"
	       << "\t" << "--bam-out       <spanning_reads.bam   "  << "\t" << "Output a BAM file containing the reads spanning each region to the provided file"    << "\n"
	       << "\t" << "--rem-multimaps                       "  << "\t" << "Remove reads that map to multiple locations (Default = False)"                       << "\n"
	       << "\t" << "--max-mate-dist <max_bp>              "  << "\t" << "Remove reads whose mate pair distance is > MAX_BP (Default = " << bam_processor.MAX_MATE_DIST << ")" << "\n"
	       << "\n";
     exit(0);
  }
 
  static struct option long_options[] = {
    {"bams",            required_argument, 0, 'b'},
    {"chrom",           required_argument, 0, 'c'},
    {"max-mate-dist",   required_argument, 0, 'd'},
    {"fasta",           required_argument, 0, 'f'},
    {"rgs",             required_argument, 0, 'g'},
    {"indexes",         required_argument, 0, 'i'},
    {"regions",         required_argument, 0, 'r'},
    {"vcf",             required_argument, 0, 'v'},
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
    case 'c':
      chrom = std::string(optarg);
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
    case 'r':
      region_file = std::string(optarg);
      break;
    case 'v':
      vcf_file = std::string(optarg);
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
  SNPBamProcessor bam_processor;
  std::string bamfile_string= "", bamindex_string="", rg_string="", region_file="", fasta_dir="", chrom="", vcf_file="", bam_out_file="";
  parse_command_line_args(argc, argv, bamfile_string, bamindex_string, rg_string, fasta_dir, region_file, vcf_file, chrom, bam_out_file, bam_processor);
  int num_flank = 0;
  if (bamfile_string.empty())
    printErrorAndDie("--bams option required");
  else if (bamindex_string.empty())
    printErrorAndDie("--indexes option required");
  else if (rg_string.empty())
    printErrorAndDie("--rgs option required");
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
  if (bam_files.size() != read_groups.size())
    printErrorAndDie("Number of .bam and RGs must match");

  // Open all .bam and .bai files
  BamTools::BamMultiReader reader;
  if (!reader.Open(bam_files)) printErrorAndDie("Failed to open one or more BAM files");
  if (!reader.OpenIndexes(bam_indexes)) printErrorAndDie("Failed to open one or more BAM index files");

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

  if (!vcf_file.empty())
    bam_processor.set_vcf(vcf_file);

  // Run analysis
  bam_processor.process_regions(reader, region_file, fasta_dir, file_read_groups, bam_writer, std::cout);

  if (!bam_out_file.empty()) bam_writer.Close();
  reader.Close();
  return 0;  
}
