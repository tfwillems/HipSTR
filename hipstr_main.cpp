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

#include "error.h"
#include "genotyper_bam_processor.h"
#include "stringops.h"

void parse_command_line_args(int argc, char** argv, 
			     std::string& bamfile_string, std::string& bamindex_string, std::string& rg_string,
			     std::string& fasta_dir, std::string& region_file,  std::string& snp_vcf_file, std::string& chrom, 
			     std::string& bam_out_file, std::string& str_vcf_out_file,  std::string& allele_vcf_out_file, std::string& viz_out_file,
			     std::string& stutter_in_file, std::string& stutter_out_file, int& use_hap_aligner, int& remove_all_filters,
			     std::string& ref_vcf_file,
			     BamProcessor& bam_processor){
  int def_mdist = bam_processor.MAX_MATE_DIST;
  if (argc == 1){
    std::cerr << "Usage: HipSTR --bams  <list_of_bams>  --indexes <list_of_bam_indexes>" << "\n"
	      << "              --fasta <dir>           --regions <region_file.bed>"     << "\n" << "\n"
	      << "Required parameters:" << "\n"
	      << "\t" << "--bams          <list_of_bams>        "  << "\t" << "Comma separated list of .bam files"                                                  << "\n"
	      << "\t" << "--indexes       <list_of_bam_indexes> "  << "\t" << "Comma separated list of .bai files in same order as .bam files"                      << "\n"
	      << "\t" << "--fasta         <dir>                 "  << "\t" << "Directory in which FASTA files for each chromosome are located"                      << "\n"
	      << "\t" << "--regions       <region_file.bed>     "  << "\t" << "BED file containing coordinates for each STR region"                                 << "\n" << "\n"
	    
	      << "Optional input parameters:" << "\n"
	      << "\t" << "--ref-vcf    <str_snp_ref_gts.vcf.gz> "  << "\t" << "Bgzipped input VCF file containing STR and SNP genotypes for a reference panel"      << "\n" 
	      << "\t" << "--snp-vcf    <phased_snp_gts.vcf.gz>  "  << "\t" << "Bgzipped input VCF file containing phased SNP genotypes for the samples"             << "\n" 
	      << "\t" << "                                      "  << "\t" << " that are going to be genotyped"                                                     << "\n"
	      << "\t" << "--stutter-in <stutter_models.txt>     "  << "\t" << "Input file containing stutter models for each locus. By default, an EM algorithm "   << "\n"
      	      << "\t" << "                                      "  << "\t" << "  will be used to learn locus-specific models"                                       << "\n" << "\n"

      	      << "Optional output parameters:" << "\n"
	      << "\t" << "--bam-out       <spanning_reads.bam   "  << "\t" << "Output a BAM file containing the reads spanning each region to the provided file"    << "\n"
	      << "\t" << "--str-vcf       <str_gts.vcf.gz>      "  << "\t" << "Output a gzipped VCF file containing phased STR genotypes"                           << "\n"
	      << "\t" << "--allele-vcf    <str_alleles.vcf>     "  << "\t" << "Output a VCF file containing alleles with strong evidence in the BAMs"               << "\n"
	      << "\t" << "--stutter-out   <stutter_models.txt>  "  << "\t" << "Output stutter models learned by the EM algorithm to the provided file"              << "\n"
	      << "\t" << "--viz-out       <aln_viz.html>        "  << "\t" << "Output an HTML file containing Needleman-Wunsch alignments for each genotyped locus" << "\n"
	      << "\t" << "                                      "  << "\t" << " Option only available when the --seq-genotyper flag has been specified"             << "\n"
	      
	      << "Other optional parameters:" << "\n"
      	      << "\t" << "--chrom         <chrom>               "  << "\t" << "Only consider STRs on the provided chromosome"                                       << "\n"
	      << "\t" << "--no-filters                          "  << "\t" << "Don't filter any putative STR reads"                                                 << "\n" 
	      << "\t" << "--max-mate-dist <max_bp>              "  << "\t" << "Remove reads whose mate pair distance is > MAX_BP (Default = " << def_mdist << ")"   << "\n"
      	      << "\t" << "--rem-multimaps                       "  << "\t" << "Remove reads that map to multiple locations (Default = False)"                       << "\n"
	      << "\t" << "--rgs           <list_of_read_groups> "  << "\t" << "Comma separated list of read groups in same order as .bam files. "                   << "\n"
	      << "\t" << "                                      "  << "\t" << "  Assign each read the RG tag corresponding to its file. By default, "               << "\n"
	      << "\t" << "                                      "  << "\t" << "  each read must have an RG flag from lobSTR and this is used instead"               << "\n"
      	      << "\t" << "--seq-genotyper                       "  << "\t" << "Use a haplotype-based aligment model to genotype each STR. This option is much "     << "\n"
	      << "\t" << "                                      "  << "\t" << "  more accurate than the default length-based model but also requires substantially" << "\n"
	      << "\t" << "                                      "  << "\t" << "  more computation time"                                                             << "\n"

	      << "\n";
    exit(0);
  }
 
  static struct option long_options[] = {
    {"allele-vcf",      required_argument, 0, 'a'},
    {"bams",            required_argument, 0, 'b'},
    {"chrom",           required_argument, 0, 'c'},
    {"max-mate-dist",   required_argument, 0, 'd'},
    {"fasta",           required_argument, 0, 'f'},
    {"rgs",             required_argument, 0, 'g'},
    {"ref-vcf",         required_argument, 0, 'h'},
    {"indexes",         required_argument, 0, 'i'},
    {"no-filters",      no_argument, &remove_all_filters, 1},
    {"str-vcf",         required_argument, 0, 'o'},
    {"regions",         required_argument, 0, 'r'},
    {"snp-vcf",         required_argument, 0, 'v'},
    {"seq-genotyper",   no_argument, &use_hap_aligner, 1},
    {"stutter-in",      required_argument, 0, 'm'},
    {"stutter-out",     required_argument, 0, 's'},
    {"bam-out",         required_argument, 0, 'w'},
    {"viz-out",         required_argument, 0, 'z'},
    {"rem-multimaps",   no_argument, &(bam_processor.REMOVE_MULTIMAPPERS), 1},
    {0, 0, 0, 0}
  };

  int c;
  while (true){
    int option_index = 0;
    c = getopt_long(argc, argv, "a:b:c:d:f:g:h:i:o:r:v:m:s:w:", long_options, &option_index);
    if (c == -1)
      break;

    switch(c){
    case 0:
      break;
    case 'a':
      allele_vcf_out_file = std::string(optarg);
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
    case 'h':
      ref_vcf_file = std::string(optarg);
      break;
    case 'i':
      bamindex_string = std::string(optarg);
      break;
    case 'm':
      stutter_in_file = std::string(optarg);
      break;
    case 'o':
      str_vcf_out_file = std::string(optarg);
      break;
    case 'r':
      region_file = std::string(optarg);
      break;
    case 's':
      stutter_out_file = std::string(optarg);
      break;
    case 'v':
      snp_vcf_file = std::string(optarg);
      break;
    case 'w':
      bam_out_file = std::string(optarg);
      break;
    case 'z':
      viz_out_file = std::string(optarg);
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
  bool check_mate_chroms = false;
  GenotyperBamProcessor bam_processor(false, check_mate_chroms, false);
  
  int use_hap_aligner = 0, remove_all_filters = 0;
  std::string bamfile_string= "", bamindex_string="", rg_string="", region_file="", fasta_dir="", chrom="", snp_vcf_file="";
  std::string bam_out_file="", str_vcf_out_file="", allele_vcf_out_file="", stutter_in_file="", stutter_out_file="", viz_out_file="";
  std::string ref_vcf_file="";
  parse_command_line_args(argc, argv, bamfile_string, bamindex_string, rg_string, fasta_dir, region_file, snp_vcf_file, chrom, 
			  bam_out_file, str_vcf_out_file, allele_vcf_out_file, viz_out_file, stutter_in_file, stutter_out_file, use_hap_aligner, remove_all_filters, 
			  ref_vcf_file, bam_processor);
  if (use_hap_aligner)
    bam_processor.use_seq_aligner();

  if (bamfile_string.empty())
    printErrorAndDie("--bams option required");
  else if (bamindex_string.empty())
    printErrorAndDie("--indexes option required");
  else if (region_file.empty())
    printErrorAndDie("--region option required");
  else if (fasta_dir.empty())
    printErrorAndDie("--fasta option required");

  if (fasta_dir.back() != '/')
    fasta_dir += "/";

  std::vector<std::string> bam_files;
  split_by_delim(bamfile_string, ',', bam_files);
  std::vector<std::string> bam_indexes;
  split_by_delim(bamindex_string, ',', bam_indexes);
  std::vector<std::string> read_groups;
  if (!rg_string.empty())
    split_by_delim(rg_string, ',', read_groups);
  std::cerr << "Detected " << bam_files.size() << " BAM files and " << bam_indexes.size() << " BAI files" << std::endl;
  if (bam_files.size() != bam_indexes.size())
    printErrorAndDie("Number of .bam and .bai files must match");

  // Open all .bam and .bai files
  BamTools::BamMultiReader reader;
  if (!reader.Open(bam_files)) {
    std::cerr << reader.GetErrorString() << std::endl;
    printErrorAndDie("Failed to open one or more BAM files");
  }
  if (!reader.OpenIndexes(bam_indexes)) {
    std::cerr << reader.GetErrorString() << std::endl;
    printErrorAndDie("Failed to open one or more BAM index files");
  }

  // Construct filename->read group map (if one has been specified) 
  // and determine the list of samples of interest based on either
  // the specified names or the RG tags in the BAM headers
  std::set<std::string> rg_samples;
  std::map<std::string, std::string> file_read_groups;
  if (!rg_string.empty()){
    if(bam_files.size() != read_groups.size())
      printErrorAndDie("Number of .bam and RGs must match");
    for (unsigned int i = 0; i < bam_files.size(); i++){
      file_read_groups[bam_files[i]] = read_groups[i];
      rg_samples.insert(read_groups[i]);
    }
    std::cerr << "User-specified read groups for " << rg_samples.size() << " unique samples" << std::endl;
  }
  else {
    bam_processor.set_lobstr_rg_usage(true);
    if (!reader.GetHeader().HasReadGroups())
      printErrorAndDie("Provided BAM files don't contain read groups in the header and the --rgs flag was not specified");

    BamTools::SamReadGroupDictionary rg_dict = reader.GetHeader().ReadGroups;
    for (auto rg_iter = rg_dict.Begin(); rg_iter != rg_dict.End(); rg_iter++){
      if (!rg_iter->HasID() || !rg_iter->HasSample())
	printErrorAndDie("RG in BAM header is lacking the ID or SM tag");
      rg_samples.insert(rg_iter->Sample);
    }
    std::cerr << "BAMs contain read groups for " << rg_samples.size() << " unique samples" << std::endl;
  }
  
  BamTools::BamWriter bam_writer;
  if (!bam_out_file.empty()){
    BamTools::RefVector ref_vector = reader.GetReferenceData();
    bool file_open = bam_writer.Open(bam_out_file, reader.GetHeaderText(), ref_vector);
    if (!file_open) printErrorAndDie("Failed to open output BAM file");
  }

  if (!ref_vcf_file.empty()){
    if (!string_ends_with(ref_vcf_file, ".gz"))
      printErrorAndDie("Ref VCF file must be bgzipped (and end in .gz)");
    bam_processor.set_ref_vcf(ref_vcf_file);
  }
  if (!snp_vcf_file.empty()){
    if (!string_ends_with(snp_vcf_file, ".gz"))
      printErrorAndDie("SNP VCF file must be bgzipped (and end in .gz)");
    bam_processor.set_input_snp_vcf(snp_vcf_file);
  }

  if(!allele_vcf_out_file.empty())
    bam_processor.set_output_allele_vcf(allele_vcf_out_file);
  if(!str_vcf_out_file.empty()){
    if (!string_ends_with(str_vcf_out_file, ".gz"))
      printErrorAndDie("Path for STR VCF output file must end in .gz as it will be gzipped");
    bam_processor.set_output_str_vcf(str_vcf_out_file, rg_samples);
  }
  if (!stutter_in_file.empty())
    bam_processor.set_input_stutter(stutter_in_file);
  if (!stutter_out_file.empty())
    bam_processor.set_output_stutter(stutter_out_file);
  if(!viz_out_file.empty())
    bam_processor.set_output_viz(viz_out_file);
  
  if (remove_all_filters)
    bam_processor.remove_all_filters();

  // Run analysis
  bam_processor.process_regions(reader, region_file, fasta_dir, file_read_groups, bam_writer, std::cout, 1000000);

  bam_processor.finish();

  if (!bam_out_file.empty()) bam_writer.Close();
  reader.Close();
  return 0;  
}
