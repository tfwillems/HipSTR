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

bool file_exists(std::string path){
  return (access(path.c_str(), F_OK) != -1);
}

void print_usage(int def_mdist){
  std::cerr << "Usage: HipSTR --bams  <list_of_bams> --fasta <dir> --regions <region_file.bed>" << "\n" << "\n"
    
	    << "Required parameters:" << "\n"
	    << "\t" << "--bams          <list_of_bams>        "  << "\t" << "Comma separated list of .bam files"                                                  << "\n"
	    << "\t" << "--fasta         <dir>                 "  << "\t" << "Directory in which FASTA files for each chromosome are located"                      << "\n"
	    << "\t" << "--regions       <region_file.bed>     "  << "\t" << "BED file containing coordinates for each STR region"                         << "\n" << "\n"
    
	    << "Optional input parameters:" << "\n"
	    << "\t" << "--ref-vcf    <str_ref_panel.vcf.gz>   "  << "\t" << "Bgzipped input VCF file containing STR (and possibly SNP) genotypes for a ref panel" << "\n" 
	    << "\t" << "                                      "  << "\t" << " This option is not available when the --len-genotyper option has been specified"    << "\n"
	    << "\t" << "--snp-vcf    <phased_snps.vcf.gz>     "  << "\t" << "Bgzipped input VCF file containing phased SNP genotypes for the samples"             << "\n" 
	    << "\t" << "                                      "  << "\t" << " that are going to be genotyped. These SNPs will be used to physically phase any "   << "\n"
	    << "\t" << "                                      "  << "\t" << " STRs in which a read or its mate pair overlaps a heterozygous site"                 << "\n"
	    << "\t" << "--stutter-in <stutter_models.txt>     "  << "\t" << "Input file containing stutter models for each locus. By default, an EM algorithm "   << "\n"
	    << "\t" << "                                      "  << "\t" << "  will be used to learn locus-specific models"                               << "\n" << "\n"
    
	    << "Optional output parameters:" << "\n"
	    << "\t" << "--bam-out       <used_reads.bam>      "  << "\t" << "Output a BAM file containing the reads used to genotype each region"                 << "\n"
	    << "\t" << "--str-vcf       <str_gts.vcf.gz>      "  << "\t" << "Output a bgzipped VCF file containing phased STR genotypes"                          << "\n"
	    << "\t" << "--allele-vcf    <str_alleles.vcf>     "  << "\t" << "Output a bgzipped VCF file containing alleles with strong evidence in the BAMs"      << "\n"
	    << "\t" << "--stutter-out   <stutter_models.txt>  "  << "\t" << "Output stutter models learned by the EM algorithm to the provided file"              << "\n"
	    << "\t" << "--viz-out       <aln_viz.html.gz>     "  << "\t" << "Output a bgzipped file containing Needleman-Wunsch alignments for each locus"        << "\n"
	    << "\t" << "                                      "  << "\t" << " The resulting file can be readily visualized with VizAln"                           << "\n"
	    << "\t" << "                                      "  << "\t" << " Option only available when the --len-genotyper flag has not been specified"         << "\n"
	    << "\t" << "--output-gls                          "  << "\t" << "Write genotype likelihoods to VCF (default = False)"                                 << "\n"
	    << "\t" << "--output-pls                          "  << "\t" << "Write phred-scaled genotype likelihoods to VCF (default = False)"                    << "\n" << "\n"
    
	    << "Other optional parameters:" << "\n"
	    << "\t" << "--help                                "  << "\t" << "Print this help message and exit"
	    << "\t" << "--chrom         <chrom>               "  << "\t" << "Only consider STRs on the provided chromosome"                                       << "\n"
	    << "\t" << "--haploid-chrs  <list_of_chroms>      "  << "\t" << "Comma separated list of chromosomes to treat as haploid"                             << "\n"
	    << "\t" << "                                      "  << "\t" << " By default, all chromosomes are treated as diploid"                                 << "\n"
	    << "\t" << "--no-filters                          "  << "\t" << "Don't filter any putative STR reads"                                                 << "\n"
	    << "\t" << "--no-rmdup                            "  << "\t" << "Don't remove PCR duplicates. By default, they'll be removed"                         << "\n"
	    << "\t" << "--max-mate-dist <max_bp>              "  << "\t" << "Remove reads whose mate pair distance is > MAX_BP (Default = " << def_mdist << ")"   << "\n"
	    << "\t" << "--rem-multimaps                       "  << "\t" << "Remove reads that map to multiple locations (Default = False)"                       << "\n"
	    << "\t" << "--rgs           <list_of_read_groups> "  << "\t" << "Comma separated list of read groups in same order as BAM files. "                    << "\n"
	    << "\t" << "                                      "  << "\t" << "  Assign each read the RG tag corresponding to its file. By default, "               << "\n"
	    << "\t" << "                                      "  << "\t" << "  each read must have an RG flag and this is used instead"                           << "\n"
	    << "\t" << "--lbs           <list_of_read_groups> "  << "\t" << "Comma separated list of libraries in same order as BAM files. "                      << "\n"
	    << "\t" << "                                      "  << "\t" << "  Assign each read the library (LB tag) corresponding to its file. By default, "     << "\n"
	    << "\t" << "                                      "  << "\t" << "  each read must have an RG flag and the associated library is used instead"         << "\n"
	    << "\t" << "                                      "  << "\t" << "  NOTE: This option is required when --rgs has been specified"                       << "\n"
	    << "\t" << "--len-genotyper                       "  << "\t" << "Use a length-based model to genotype each STR. This option is much"                  << "\n"
	    << "\t" << "                                      "  << "\t" << "  faster than the default sequence-based model but does not model the underlying"    << "\n"
	    << "\t" << "                                      "  << "\t" << "  STR sequence. As a result, it cannot detect homoplasy and all STR alleles output"  << "\n"
	    << "\t" << "                                      "  << "\t" << "  in the VCF assume that indels are perfect copies of the repeat motif"              << "\n"
	    << "\t" << "                                      "  << "\t" << "  The length-based approach is very similar to lobSTR's allelotype module"           << "\n"
	    << "\n";
}
  
void parse_command_line_args(int argc, char** argv, 
			     std::string& bamfile_string,  std::string& rg_string,        std::string& haploid_chr_string,
			     std::string& fasta_dir,       std::string& region_file,      std::string& snp_vcf_file,        std::string& chrom,
			     std::string& bam_out_file,    std::string& str_vcf_out_file, std::string& allele_vcf_out_file, std::string& viz_out_file,
			     std::string& stutter_in_file, std::string& stutter_out_file, int& use_hap_aligner, int& remove_all_filters, int& remove_pcr_dups,
			     int& output_gls, int& output_pls, std::string& ref_vcf_file,
			     BamProcessor& bam_processor){
  int def_mdist = bam_processor.MAX_MATE_DIST;
  if (argc == 1){
    print_usage(def_mdist);
    exit(0);
  }

  int print_help = 0;
  
  static struct option long_options[] = {
    {"allele-vcf",      required_argument, 0, 'a'},
    {"bams",            required_argument, 0, 'b'},
    {"chrom",           required_argument, 0, 'c'},
    {"max-mate-dist",   required_argument, 0, 'd'},
    {"fasta",           required_argument, 0, 'f'},
    {"rgs",             required_argument, 0, 'g'},
    {"h",               no_argument, &print_help, 1},
    {"help",            no_argument, &print_help, 1},
    {"len-genotyper",   no_argument, &use_hap_aligner,    0},
    {"no-filters",      no_argument, &remove_all_filters, 1},
    {"no-rmdup",        no_argument, &remove_pcr_dups,    0},
    {"output-gls",      no_argument, &output_gls, 1},
    {"output-pls",      no_argument, &output_pls, 1},
    {"str-vcf",         required_argument, 0, 'o'},
    {"ref-vcf",         required_argument, 0, 'p'},
    {"regions",         required_argument, 0, 'r'},
    {"snp-vcf",         required_argument, 0, 'v'},
    {"stutter-in",      required_argument, 0, 'm'},
    {"stutter-out",     required_argument, 0, 's'},
    {"haploid-chrs",    required_argument, 0, 't'},
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
    case 'm':
      stutter_in_file = std::string(optarg);
      break;
    case 'o':
      str_vcf_out_file = std::string(optarg);
      break;
    case 'p':
      ref_vcf_file = std::string(optarg);
      break;
    case 'r':
      region_file = std::string(optarg);
      break;
    case 's':
      stutter_out_file = std::string(optarg);
      break;
    case 't':
      haploid_chr_string = std::string(optarg);
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

  if (print_help){
    print_usage(def_mdist);
    exit(0);
  }
}

int main(int argc, char** argv){
  bool check_mate_chroms = false;
  GenotyperBamProcessor bam_processor(true, check_mate_chroms, true, true);
  
  int use_hap_aligner = 1, remove_all_filters = 0, remove_pcr_dups = 1;
  std::string bamfile_string= "", rg_string="", lb_string="", hap_chr_string="", region_file="", fasta_dir="", chrom="", snp_vcf_file="";
  std::string bam_out_file="", str_vcf_out_file="", allele_vcf_out_file="", stutter_in_file="", stutter_out_file="", viz_out_file="";
  int output_gls = 0, output_pls = 0;
  std::string ref_vcf_file="";
  parse_command_line_args(argc, argv, bamfile_string, rg_string, hap_chr_string, fasta_dir, region_file, snp_vcf_file, chrom,
			  bam_out_file, str_vcf_out_file, allele_vcf_out_file, viz_out_file, stutter_in_file, stutter_out_file, use_hap_aligner, remove_all_filters, 
			  remove_pcr_dups, output_gls, output_pls, ref_vcf_file, bam_processor);
  if (output_gls) bam_processor.output_gls();
  if (output_pls) bam_processor.output_pls();
  if (remove_pcr_dups == 0)
    bam_processor.allow_pcr_dups();

  if (!use_hap_aligner) {
    bam_processor.use_len_model();
    if (!ref_vcf_file.empty())
      printErrorAndDie("--ref-vcf option is not compatible with the --len-genotyper option");
    if (!viz_out_file.empty())
      printErrorAndDie("--viz-out option is not compatible with the --len-genotyper option");
  }
    
  
  if (bamfile_string.empty())
    printErrorAndDie("--bams option required");
  else if (region_file.empty())
    printErrorAndDie("--region option required");
  else if (fasta_dir.empty())
    printErrorAndDie("--fasta option required");

  if (fasta_dir.back() != '/')
    fasta_dir += "/";

  std::vector<std::string> bam_files;
  split_by_delim(bamfile_string, ',', bam_files);
  std::cerr << "Detected " << bam_files.size() << " BAM files" << std::endl;

  // Open all BAM files
  BamTools::BamMultiReader reader;
  if (!reader.Open(bam_files)) {
    std::cerr << reader.GetErrorString() << std::endl;
    printErrorAndDie("Failed to open one or more BAM files");
  }

  // Open BAM index files, assuming they're the same path with a .bai suffix
  std::vector<std::string> bam_indexes;
  for (unsigned int i = 0; i < bam_files.size(); i++){
    std::string bai_file = bam_files[i] + ".bai";
    if (!file_exists(bai_file))
      printErrorAndDie("BAM index file " + bai_file + " does not exist. Please ensure that each BAM has been sorted by position and indexed using samtools");
    bam_indexes.push_back(bai_file);
  }
  if (!reader.OpenIndexes(bam_indexes)) {
    std::cerr << reader.GetErrorString() << std::endl;
    printErrorAndDie("Failed to open one or more BAM index files");
  }

  // Construct filename->read group map (if one has been specified) and determine the list
  // of samples of interest based on either the specified names or the RG tags in the BAM headers
  std::set<std::string> rg_samples, rg_libs;
  std::map<std::string, std::string> rg_ids_to_sample, rg_ids_to_library;
  if (!rg_string.empty()){
    if (lb_string.empty())
      printErrorAndDie("--lbs option required when --rgs option specified");

    std::vector<std::string> read_groups, libraries;
    split_by_delim(rg_string, ',', read_groups);
    split_by_delim(lb_string, ',', libraries);
    if (bam_files.size() != read_groups.size())
      printErrorAndDie("Number of BAM files in --bams and RGs in --rgs must match");
    if (bam_files.size() != libraries.size())
      printErrorAndDie("Number of BAM files in --bams and LBs in --lbs must match");

    for (unsigned int i = 0; i < bam_files.size(); i++){
      rg_ids_to_sample[bam_files[i]]  = read_groups[i];
      rg_ids_to_library[bam_files[i]] = read_groups[i];
      rg_samples.insert(read_groups[i]);
    }
    bam_processor.use_custom_read_groups();
    std::cerr << "User-specified read groups for " << rg_samples.size() << " unique samples" << std::endl;
  }
  else {
    if (!reader.GetHeader().HasReadGroups())
      printErrorAndDie("Provided BAM files don't contain read groups in the header and the --rgs flag was not specified");

    BamTools::SamReadGroupDictionary rg_dict = reader.GetHeader().ReadGroups;
    for (auto rg_iter = rg_dict.Begin(); rg_iter != rg_dict.End(); rg_iter++){
      if (!rg_iter->HasID() || !rg_iter->HasSample() || !rg_iter->HasLibrary())
	printErrorAndDie("RG in BAM header is lacking the ID, SM or LB tag");

      // Ensure that there aren't identical read group ids that map to different samples or libraries
      if (rg_ids_to_sample.find(rg_iter->ID) != rg_ids_to_sample.end())
	if (rg_ids_to_sample[rg_iter->ID].compare(rg_iter->Sample) != 0)
	  printErrorAndDie("Read group id " + rg_iter->ID + " maps to more than one sample");
      if (rg_ids_to_library.find(rg_iter->ID) != rg_ids_to_library.end())
	if (rg_ids_to_library[rg_iter->ID].compare(rg_iter->Library) != 0)
	  printErrorAndDie("Read group id " + rg_iter->ID + " maps to more than one library");

      rg_ids_to_sample[rg_iter->ID]  = rg_iter->Sample;
      rg_ids_to_library[rg_iter->ID] = rg_iter->Library; 
      rg_samples.insert(rg_iter->Sample);
      rg_libs.insert(rg_iter->Library);
    }
    std::cerr << "BAMs contain " << rg_ids_to_sample.size() << " unique read group IDs for "
	      << rg_libs.size()    << " unique libraries and "
	      << rg_samples.size() << " unique samples" << std::endl;
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

    // Check that the VCF exists
    if (!file_exists(ref_vcf_file)) 
      printErrorAndDie("Ref VCF file " + ref_vcf_file + " does not exist. Please ensure that the path provided to --ref-vcf is valid");

    // Check that tabix index exists
    if (!file_exists(ref_vcf_file + ".tbi"))
	printErrorAndDie("No .tbi index found for the ref VCF file. Please index using tabix and rerun HipSTR");

    bam_processor.set_ref_vcf(ref_vcf_file);
  }
  if (!snp_vcf_file.empty()){
    if (!string_ends_with(snp_vcf_file, ".gz"))
      printErrorAndDie("SNP VCF file must be bgzipped (and end in .gz)");
    
    // Check that the VCF exists
    if (!file_exists(snp_vcf_file))
      printErrorAndDie("SNP VCF file " + snp_vcf_file + " does not exist. Please ensure that the path provided to --snp-vcf is valid");

    // Check that tabix index exists
    if (!file_exists(snp_vcf_file + ".tbi"))
	printErrorAndDie("No .tbi index found for the SNP VCF file. Please index using tabix and rerun HipSTR");

    bam_processor.set_input_snp_vcf(snp_vcf_file);
  }

  if(!allele_vcf_out_file.empty())
    bam_processor.set_output_allele_vcf(allele_vcf_out_file);
  if(!str_vcf_out_file.empty()){
    if (!string_ends_with(str_vcf_out_file, ".gz"))
      printErrorAndDie("Path for STR VCF output file must end in .gz as it will be bgzipped");
    bam_processor.set_output_str_vcf(str_vcf_out_file, rg_samples);
  }
  if (!stutter_in_file.empty())
    bam_processor.set_input_stutter(stutter_in_file);
  if (!stutter_out_file.empty())
    bam_processor.set_output_stutter(stutter_out_file);
  if(!viz_out_file.empty()){
    if (!string_ends_with(viz_out_file, ".gz"))
      printErrorAndDie("Path for alignment visualization file must end in .gz as it will be bgzipped");
    bam_processor.set_output_viz(viz_out_file);
  }
  
  if (remove_all_filters)
    bam_processor.remove_all_filters();

  if (!hap_chr_string.empty()){
    std::vector<std::string> haploid_chroms;
    split_by_delim(hap_chr_string, ',', haploid_chroms);
    for (auto chrom_iter = haploid_chroms.begin(); chrom_iter != haploid_chroms.end(); chrom_iter++)
      bam_processor.add_haploid_chrom(*chrom_iter);
  }

  // Run analysis
  bam_processor.process_regions(reader, region_file, fasta_dir, rg_ids_to_sample, rg_ids_to_library, bam_writer, std::cout, 1000000, chrom);
  bam_processor.finish();

  if (!bam_out_file.empty()) bam_writer.Close();
  reader.Close();
  return 0;  
}
