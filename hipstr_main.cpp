#include <climits>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <time.h>

#include "bamtools/include/api/BamAlignment.h"
#include "vcflib/src/Variant.h"

#include "error.h"
#include "genotyper_bam_processor.h"
#include "seqio.h"
#include "stringops.h"
#include "version.h"

bool file_exists(std::string path){
  return (access(path.c_str(), F_OK) != -1);
}

void print_usage(int def_mdist, int def_min_reads, int def_max_reads, int def_max_str_len){
  std::cerr << "Usage: HipSTR --bams <list_of_bams> --fasta <dir> --regions <region_file.bed> [OPTIONS]" << "\n" << "\n"
    
	    << "Required parameters:" << "\n"
	    << "\t" << "--bams          <list_of_bams>        "  << "\t" << "Comma separated list of BAM files. Either --bams or --bam-files must be specified"   << "\n"
	    << "\t" << "--fasta         <dir>                 "  << "\t" << "Directory in which FASTA files for each chromosome are located or the path to a"     << "\n"
	    << "\t" << "                                      "  << "\t" << " single FASTA file that contains all of the relevant sequences"                      << "\n"
	    << "\t" << "--regions       <region_file.bed>     "  << "\t" << "BED file containing coordinates for each STR region"                         << "\n" << "\n"
    
	    << "Optional input parameters:" << "\n"
	    << "\t" << "--bam-files  <bam_files.txt>         "  << "\t" << "File containing BAM files to analyze, one per line."                                 << "\n"
	    << "\t" << "--ref-vcf    <str_ref_panel.vcf.gz>   "  << "\t" << "Bgzipped input VCF file containing a reference panel of STR genotypes"               << "\n"
	    << "\t" << "                                      "  << "\t" << " For each locus, alleles in the VCF will be used as candidate STR variants instead"  << "\n"
	    << "\t" << "                                      "  << "\t" << " of trying to identify candidates from the BAMs (default)"                           << "\n"
	    << "\t" << "                                      "  << "\t" << " This option is not available when the --len-genotyper option has been specified"    << "\n"
	    << "\t" << "--snp-vcf    <phased_snps.vcf.gz>     "  << "\t" << "Bgzipped input VCF file containing phased SNP genotypes for the samples"             << "\n" 
	    << "\t" << "                                      "  << "\t" << " that are going to be genotyped. These SNPs will be used to physically phase any "   << "\n"
	    << "\t" << "                                      "  << "\t" << " STRs when a read or its mate pair overlaps a heterozygous site"                     << "\n"
	    << "\t" << "--stutter-in <stutter_models.txt>     "  << "\t" << "Input file containing stutter models for each locus. By default, an EM algorithm "   << "\n"
	    << "\t" << "                                      "  << "\t" << "  will be used to learn locus-specific models"                               << "\n" << "\n"
    
	    << "Optional output parameters:" << "\n"
	    << "\t" << "--str-vcf       <str_gts.vcf.gz>      "  << "\t" << "Output a bgzipped VCF file containing STR genotypes"                                 << "\n"
	    << "\t" << "                                      "  << "\t" << " NOTE: If you don't specify this option, no genotyping will be performed"            << "\n"
	    << "\t" << "--stutter-out   <stutter_models.txt>  "  << "\t" << "Output stutter models learned by the EM algorithm to the provided file"              << "\n"
	    << "\t" << "--log <log.txt>                       "  << "\t" << "Output the log information to the provided file. By default, the log will be "       << "\n"
	    << "\t" << "                                      "  << "\t" << " written to standard err"                                                            << "\n"
	    << "\t" << "--viz-out       <aln_viz.html.gz>     "  << "\t" << "Output a bgzipped file containing haplotype alignments for each locus"               << "\n"
	    << "\t" << "                                      "  << "\t" << " The resulting file can be readily visualized with VizAln"                           << "\n"
	    << "\t" << "                                      "  << "\t" << " Option only available when the --len-genotyper flag has not been specified"         << "\n"
	    << "\t" << "--viz-left-alns                       "  << "\t" << "Output the original left aligned reads to the HTML output in addition to the "       << "\n"
	    << "\t" << "                                      "  << "\t" << " haplotype alignments. By default, only the latter is output"                        << "\n"
	    << "\t" << "--pass-bam      <used_reads.bam>      "  << "\t" << "Output a BAM file containing the reads used to genotype each region"                 << "\n"
	    << "\t" << "--filt-bam      <filt_reads.bam>      "  << "\t" << "Output a BAM file containing the reads filtered in each region. Each BAM entry"      << "\n"
	    << "\t" << "                                      "  << "\t" << " has an FT tag specifying the reason for filtering"                                  << "\n"
	    << "\t" << "--max-flank-indel <max_flank_frac>    "  << "\t" << "Don't output genotypes for a sample if the fraction of reads containing an indel"    << "\n"
	    << "\t" << "                                      "  << "\t" << " in the sequence flanking the STR is greater than MAX_FLANK_FRAC (Default = 0.15)"   << "\n"
	    << "\t" << "--hide-allreads                       "  << "\t" << "Don't output the ALLREADS  FORMAT field to the VCF. By default, it will be output"   << "\n"
	    << "\t" << "--hide-mallreads                      "  << "\t" << "Don't output the MALLREADS FORMAT field to the VCF. By default, it will be output"   << "\n"
	    << "\t" << "--output-pallreads                    "  << "\t" << "Output the PALLREADS FORMAT field to the VCF. By default, it will not be output"     << "\n"
	    << "\t" << "--output-gls                          "  << "\t" << "Write genotype likelihoods to VCF (Default = False)"                                 << "\n"
	    << "\t" << "--output-pls                          "  << "\t" << "Write phred-scaled genotype likelihoods to VCF (Default = False)"                    << "\n" << "\n"

	    << "Optional read filtering parameters:" << "\n"
	    << "\t" << "--no-rmdup                            "  << "\t" << "Don't remove PCR duplicates. By default, they'll be removed"                         << "\n"
	    << "\t" << "--max-mate-dist <max_bp>              "  << "\t" << "Remove reads whose mate pair distance is > MAX_BP (Default = " << def_mdist << ")"   << "\n"
	    << "\t" << "--use-unpaired                        "  << "\t" << "Use unpaired reads when genotyping. (Default = False)"                               << "\n"
	    << "\t" << "--min-mapq      <min_qual>            "  << "\t" << "Remove reads with a mapping quality below this threshold (Default = 0)"              << "\n" << "\n"

	    << "Other optional parameters:" << "\n"
	    << "\t" << "--help                                "  << "\t" << "Print this help message and exit"                                                    << "\n"
	    << "\t" << "--chrom         <chrom>               "  << "\t" << "Only consider STRs on the provided chromosome"                                       << "\n"
	    << "\t" << "--haploid-chrs  <list_of_chroms>      "  << "\t" << "Comma separated list of chromosomes to treat as haploid"                             << "\n"
	    << "\t" << "                                      "  << "\t" << " By default, all chromosomes are treated as diploid"                                 << "\n"
	    << "\t" << "--min-reads     <num_reads>           "  << "\t" << "Minimum total reads required to genotype a locus (Default = " << def_min_reads << ")" << "\n"
	    << "\t" << "--max-reads     <num_reads>           "  << "\t" << "Skip a locus if it has more than NUM_READS reads (Default = " << def_max_reads << ")" << "\n"
	    << "\t" << "--max-str-len   <max_bp>              "  << "\t" << "Only genotype STRs in the provided BED file with length < MAX_BP (Default = " << def_max_str_len << ")" << "\n"
	    << "\t" << "--bam-samps     <list_of_samples>     "  << "\t" << "Comma separated list of read groups in same order as BAM files. "                    << "\n"
	    << "\t" << "                                      "  << "\t" << "  Assign each read the read group corresponding to its file. By default, "           << "\n"
	    << "\t" << "                                      "  << "\t" << "  each read must have an RG tag and the sample is determined from the SM field"      << "\n"
	    << "\t" << "--bam-libs      <list_of_libraries>   "  << "\t" << "Comma separated list of libraries in same order as BAM files. "                      << "\n"
	    << "\t" << "                                      "  << "\t" << "  Assign each read the library corresponding to its file. By default, "              << "\n"
	    << "\t" << "                                      "  << "\t" << "  each read must have an RG tag and the library is determined from the LB field"     << "\n"
	    << "\t" << "                                      "  << "\t" << "  NOTE: This option is required when --bam-samps has been specified"                 << "\n"
	    << "\t" << "--lib-from-samp                       "  << "\t" << " Assign each read the library corresponding to its sample name. By default,  "       << "\n"
	    << "\t" << "                                      "  << "\t" << "  each read must have an RG tag and the library is determined from the LB field"     << "\n"
	    << "\t" << "--10x-bams                            "  << "\t" << "BAM files were generated by 10X Genomics. HipSTR will utilize haplotype tags in the" << "\n"
	    << "\t" << "                                      "  << "\t" << "  BAMs to phase and more accurately genotype STRs (Experimental)"                    << "\n"
	    << "\t" << "--len-genotyper                       "  << "\t" << "Use a length-based model to genotype each STR. This option is much"                  << "\n"
	    << "\t" << "                                      "  << "\t" << "  faster than the default sequence-based model but does not model the underlying"    << "\n"
	    << "\t" << "                                      "  << "\t" << "  STR sequence. As a result, it cannot detect homoplasy and all STR alleles output"  << "\n"
	    << "\t" << "                                      "  << "\t" << "  in the VCF assume that indels are perfect copies of the repeat motif"              << "\n"
	    << "\t" << "                                      "  << "\t" << "  The length-based approach is very similar to lobSTR's allelotype module"           << "\n"
	    << "\t" << "--no-pool-seqs                        "  << "\t" << "Do not merge reads with identical sequences and combine their base quality scores."  << "\n"
	    << "\t" << "                                      "  << "\t" << "  By default, pooled reads will be aligned using the haplotype aligner instead"      << "\n"
	    << "\t" << "                                      "  << "\t" << "  of the reads themselves, resulting in a large speedup."                            << "\n"
	    << "\t" << "--def-stutter-model                   "  << "\t" << "For each locus, use a stutter model with PGEOM=0.9, UP=0.05, DOWN=0.05 for "         << "\n"
	    << "\t" << "                                      "  << "\t" << " in-frame artifacts and PGEOM=0.9, UP=0.01, DOWN=0.01 for out-of-frame artifacts"    << "\n"
	    << "\t" << "--use-all-reads                       "  << "\t" << "Use all reads overlapping the region to genotype the STR, even those that are "      << "\n"
	    << "\t" << "                                      "  << "\t" << " unlikely to be informative. By default, HipSTR only utilizes the reads it thinks"   << "\n"
	    << "\t" << "                                      "  << "\t" << " will be informative. Enabling this option can slightly increase accuracy but at"    << "\n"
	    << "\t" << "                                      "  << "\t" << " the expense of significantly longer runtimes (~4-5 times longer)."                  << "\n"
	    << "\t" << "--read-qual-trim <min_qual>           "  << "\t" << "Trim both ends of the read until reaching a base with quality > MIN_QUAL"            << "\n"
	    << "\t" << "                                      "  << "\t" << " By default, no quality trimmming is performed"                                      << "\n"
	    << "\t" << "--version                             "  << "\t" << "Print HipSTR version and exit"                                                       << "\n"  
	    << "\n";
}
  
void parse_command_line_args(int argc, char** argv, 
			     std::string& bamfile_string,     std::string& bamlist_string,    std::string& rg_sample_string,      std::string& rg_lib_string,
			     std::string& haploid_chr_string, std::string& fasta_dir,         std::string& region_file,           std::string& snp_vcf_file,
			     std::string& chrom,              std::string& bam_pass_out_file, std::string& bam_filt_out_file,
			     std::string& str_vcf_out_file,   std::string& log_file,      int& use_all_reads,
			     int& use_hap_aligner, int& remove_pcr_dups,  int& bams_from_10x,    int& bam_lib_from_samp,     int& def_stutter_model,
			     int& output_gls,      int& output_pls,       int& output_all_reads, int& output_pall_reads,     int& output_mall_reads, std::string& ref_vcf_file,
			     GenotyperBamProcessor& bam_processor){
  int def_mdist       = bam_processor.MAX_MATE_DIST;
  int def_min_reads   = bam_processor.MIN_TOTAL_READS;
  int def_max_reads   = bam_processor.MAX_TOTAL_READS;
  int def_max_str_len = bam_processor.MAX_STR_LENGTH;
  if (argc == 1 || (argc == 2 && std::string("-h").compare(std::string(argv[1])) == 0)){
    print_usage(def_mdist, def_min_reads, def_max_reads, def_max_str_len);
    exit(0);
  }

  int print_help           = 0;
  int condense_read_fields = 1;
  int pool_seqs            = 1;
  int viz_left_alns        = 0;
  int print_version        = 0;

  static struct option long_options[] = {
    {"10x-bams",        no_argument, &bams_from_10x, 1},
    {"bams",            required_argument, 0, 'b'},
    {"bam-files",       required_argument, 0, 'B'},
    {"chrom",           required_argument, 0, 'c'},
    {"max-mate-dist",   required_argument, 0, 'd'},
    {"fasta",           required_argument, 0, 'f'},
    {"bam-samps",       required_argument, 0, 'g'},
    {"bam-libs",        required_argument, 0, 'q'},
    {"lib-from-samp",    no_argument, &bam_lib_from_samp,    1},
    {"full-read-fields", no_argument, &condense_read_fields, 0},
    {"min-mapq",        required_argument, 0, 'e'},
    {"min-reads",       required_argument, 0, 'i'},
    {"read-qual-trim",  required_argument, 0, 'j'},
    {"log",             required_argument, 0, 'l'},
    {"max-reads",       required_argument, 0, 'n'},
    {"h",               no_argument, &print_help, 1},
    {"help",            no_argument, &print_help, 1},
    {"hide-allreads",   no_argument, &output_all_reads,   0},
    {"hide-mallreads",  no_argument, &output_mall_reads,  0},
    {"output-pallreads",no_argument, &output_pall_reads,  1},
    {"len-genotyper",   no_argument, &use_hap_aligner,    0},
    {"no-rmdup",        no_argument, &remove_pcr_dups,    0},
    {"output-gls",      no_argument, &output_gls, 1},
    {"output-pls",      no_argument, &output_pls, 1},
    {"no-pool-seqs",    no_argument, &pool_seqs,  0},
    {"version",         no_argument, &print_version, 1},
    {"max-flank-indel", required_argument, 0, 'F'},
    {"str-vcf",         required_argument, 0, 'o'},
    {"ref-vcf",         required_argument, 0, 'p'},
    {"regions",         required_argument, 0, 'r'},
    {"use-unpaired",    no_argument, &(bam_processor.REQUIRE_PAIRED_READS), 0},
    {"use-all-reads",   no_argument, &use_all_reads, 1},
    {"def-stutter-model", no_argument, &def_stutter_model, 1},
    {"snp-vcf",         required_argument, 0, 'v'},
    {"stutter-in",      required_argument, 0, 'm'},
    {"stutter-out",     required_argument, 0, 's'},
    {"haploid-chrs",    required_argument, 0, 't'},
    {"pass-bam",        required_argument, 0, 'w'},
    {"max-str-len",     required_argument, 0, 'x'},
    {"filt-bam",        required_argument, 0, 'y'},
    {"viz-left-alns",   no_argument, &viz_left_alns, 1},
    {"viz-out",         required_argument, 0, 'z'},
    {0, 0, 0, 0}
  };

  std::string filename;
  int c;
  while (true){
    int option_index = 0;
    c = getopt_long(argc, argv, "b:B:c:d:e:f:g:i:j:k:l:m:n:o:p:q:r:s:t:v:w:x:y:z:", long_options, &option_index);
    if (c == -1)
      break;

    switch(c){
    case 0:
      break;
    case 'b':
      bamlist_string = std::string(optarg);
      break;
    case 'B':
      bamfile_string = std::string(optarg);
      break;
    case 'c':
      chrom = std::string(optarg);
      break;
    case 'd':
      bam_processor.MAX_MATE_DIST = atoi(optarg);
      break;
    case 'e':
      bam_processor.MIN_MAPPING_QUALITY = atoi(optarg);
      break;
    case 'f':
      fasta_dir = std::string(optarg);
      break;
    case 'g':
      rg_sample_string = std::string(optarg);
      break;
    case 'i':
      bam_processor.MIN_TOTAL_READS = atoi(optarg);
      if (bam_processor.MIN_TOTAL_READS < 1)
	printErrorAndDie("--min-total-reads must be greater than 0");
      break;
    case 'j':
      if (std::string(optarg).size() != 1)
	printErrorAndDie("--read-qual-trim requires a single character argument");
      bam_processor.BASE_QUAL_TRIM = std::string(optarg)[0];
      break;
    case 'l':
      log_file = std::string(optarg);
      break;
    case 'm':
      filename = std::string(optarg);
      bam_processor.set_input_stutter(filename);
      break;
    case 'n':
      bam_processor.MAX_TOTAL_READS = atoi(optarg);
      break;
    case 'o':
      str_vcf_out_file = std::string(optarg);
      break;
    case 'p':
      ref_vcf_file = std::string(optarg);
      break;
    case 'q':
      rg_lib_string = std::string(optarg);
      break;
    case 'r':
      region_file = std::string(optarg);
      break;
    case 's':
      filename = std::string(optarg);
      bam_processor.set_output_stutter(filename);
      break;
    case 't':
      haploid_chr_string = std::string(optarg);
      break;
    case 'v':
      snp_vcf_file = std::string(optarg);
      break;
    case 'w':
      bam_pass_out_file = std::string(optarg);
      break;
    case 'x':
      bam_processor.MAX_STR_LENGTH = atoi(optarg);
      break;
    case 'y':
      bam_filt_out_file = std::string(optarg);
      break;
    case 'z':
      filename = std::string(optarg);
      if (!string_ends_with(filename, ".gz"))
	printErrorAndDie("Path for alignment visualization file must end in .gz as it will be bgzipped");
      bam_processor.set_output_viz(filename);
      break;
    case 'F':
      bam_processor.set_max_flank_indel_frac(atof(optarg));
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
    msg << "Did not recognizee the following command line arguments:" << "\n";
    while (optind < argc)
      msg << "\t" << argv[optind++] << "\n";
    msg << "Please check your command line syntax or type ./HipSTR --help for additional information" << "\n";
    printErrorAndDie(msg.str());
  }

  if (print_version == 1){
    std::cerr << "HipSTR version " << VERSION << std::endl;
    exit(0);
  }
  if (pool_seqs == 1)
    bam_processor.pool_sequences();
  if (print_help){
    print_usage(def_mdist, def_min_reads, def_max_reads,  def_max_str_len);
    exit(0);
  }
  if (viz_left_alns)
    bam_processor.visualize_left_alns();
  SeqStutterGenotyper::condense_read_count_fields = (condense_read_fields == 0 ? false : true);
}

int main(int argc, char** argv){
  double total_time = clock();
  precompute_integer_logs(); // Calculate and cache log of integers from 1 -> 999

  std::stringstream full_command_ss;
  full_command_ss << "HipSTR-" << VERSION;
  for (int i = 1; i < argc; i++)
    full_command_ss << " " << argv[i];
  std::string full_command = full_command_ss.str();

  GenotyperBamProcessor bam_processor(true, true, true);
  bam_processor.set_max_flank_indel_frac(0.15);

  int use_all_reads = 0, use_hap_aligner = 1, remove_pcr_dups = 1, bams_from_10x = 0, bam_lib_from_samp = 0, def_stutter_model = 0;
  std::string bamfile_string= "", bamlist_string = "", rg_sample_string="", rg_lib_string="", hap_chr_string="", region_file="", fasta_dir="", chrom="", snp_vcf_file="";
  std::string bam_pass_out_file="", bam_filt_out_file="", str_vcf_out_file="", log_file = "";
  int output_gls = 0, output_pls = 0, output_all_reads = 1, output_pall_reads = 0, output_mall_reads = 1;
  std::string ref_vcf_file="";
  parse_command_line_args(argc, argv, bamfile_string, bamlist_string, rg_sample_string, rg_lib_string, hap_chr_string, fasta_dir, region_file, snp_vcf_file, chrom,
			  bam_pass_out_file, bam_filt_out_file, str_vcf_out_file, log_file, use_all_reads, use_hap_aligner, remove_pcr_dups, bams_from_10x,
			  bam_lib_from_samp, def_stutter_model, output_gls, output_pls, output_all_reads, output_pall_reads, output_mall_reads,
			  ref_vcf_file, bam_processor);

  if (!log_file.empty())
    bam_processor.set_log(log_file);
  if (bams_from_10x){
    bam_processor.use_10x_bam_tags();
    bam_processor.logger() << "Using 10X BAM tags to genotype and phase STRs (WARNING: Any arguments provided to --snp-vcf will be ignored)" << std::endl;
  }
  if (use_all_reads == 1)     bam_processor.REQUIRE_SPANNING = false;
  if (output_gls)             bam_processor.output_gls();
  if (output_pls)             bam_processor.output_pls();
  if (output_all_reads == 0)  bam_processor.hide_all_reads();
  if (output_pall_reads == 0) bam_processor.hide_pall_reads();
  if (output_mall_reads == 0) bam_processor.hide_mall_reads();
  if (remove_pcr_dups == 0)   bam_processor.allow_pcr_dups();
  if (def_stutter_model == 1) bam_processor.set_default_stutter_model(0.95, 0.05, 0.05, 0.95, 0.01, 0.01);

  if (!use_hap_aligner) {
    bam_processor.use_len_model();
    if (!ref_vcf_file.empty())
      printErrorAndDie("--ref-vcf option is not compatible with the --len-genotyper option");
    //if (!viz_out_file.empty())
    //printErrorAndDie("--viz-out option is not compatible with the --len-genotyper option");
  }

  if (bamfile_string.empty() && bamlist_string.empty())
    printErrorAndDie("You must specify either the --bams or --bam-files option");
  else if (region_file.empty())
    printErrorAndDie("--region option required");
  else if (fasta_dir.empty())
    printErrorAndDie("--fasta option required");

  if (fasta_dir.back() != '/' && !is_file(fasta_dir))
    fasta_dir += "/";

  std::vector<std::string> bam_files;
  if (!bamlist_string.empty())
    split_by_delim(bamlist_string, ',', bam_files);

  // Read file containing BAM paths line-by-line
  if (!bamfile_string.empty()){
    if (!file_exists(bamfile_string))
      printErrorAndDie("File containing BAM file names does not exist: " + bamfile_string);
    std::ifstream input(bamfile_string.c_str());
    if (!input.is_open())
      printErrorAndDie("Failed to open file containing BAM file names: " + bamfile_string);
    std::string line;
    while (std::getline(input, line))
      if (!line.empty())
	bam_files.push_back(line);
    input.close();
  }
  bam_processor.logger() << "Detected " << bam_files.size() << " BAM files" << std::endl;

  // Open all BAM files
  BamTools::BamMultiReader reader;
  if (!reader.Open(bam_files)) {
    std::cerr << reader.GetErrorString() << std::endl;
    printErrorAndDie("Failed to open one or more BAM files");
  }
  //reader.SetExplicitMergeOrder(BamTools::BamMultiReader::MergeOrder::RoundRobinMerge);


  // Construct filename->read group map (if one has been specified) and determine the list
  // of samples of interest based on either the specified names or the RG tags in the BAM headers
  std::set<std::string> rg_samples, rg_libs;
  std::map<std::string, std::string> rg_ids_to_sample, rg_ids_to_library;
  if (!rg_sample_string.empty()){
    if (rg_lib_string.empty())
      printErrorAndDie("--bam-libs option required when --bam-samps option specified");

    std::vector<std::string> read_groups, libraries;
    split_by_delim(rg_sample_string, ',', read_groups);
    split_by_delim(rg_lib_string, ',', libraries);
    if (bam_files.size() != read_groups.size())
      printErrorAndDie("Number of BAM files in --bams and samples in --bam-samps must match");
    if (bam_files.size() != libraries.size())
      printErrorAndDie("Number of BAM files in --bams and libraries in --bam-libs must match");

    for (unsigned int i = 0; i < bam_files.size(); i++){
      rg_ids_to_sample[bam_files[i]]  = read_groups[i];
      rg_ids_to_library[bam_files[i]] = libraries[i];
      rg_samples.insert(read_groups[i]);
    }
    bam_processor.use_custom_read_groups();
    bam_processor.logger() << "User-specified read groups for " << rg_samples.size() << " unique samples" << std::endl;
  }
  else {
    if (!reader.GetHeader().HasReadGroups())
      printErrorAndDie("Provided BAM files don't contain read groups in the header and the --bam-samps flag was not specified");

    for (unsigned int i = 0; i < bam_files.size(); i++){
      // Although it would be ideal to use the dictionary provided by the BamMultiReader, it doesn't appropriately handle
      // read group ID collisions across BAMs (as it just overwrites the previous entry).
      // Instead, we can open an individual reader for each BAM and allow for conficting IDs, as long as they lie in separate files

      BamTools::BamReader single_file_reader;
      if (!single_file_reader.Open(bam_files[i]))
	printErrorAndDie("Failed to open one or more BAM files");

      BamTools::SamReadGroupDictionary rg_dict = single_file_reader.GetHeader().ReadGroups;
      for (auto rg_iter = rg_dict.Begin(); rg_iter != rg_dict.End(); rg_iter++){
	if (!rg_iter->HasID() || !rg_iter->HasSample() || ((bam_lib_from_samp == 0) && !rg_iter->HasLibrary()))
	  printErrorAndDie("RG in BAM header is lacking the ID, SM or LB tag");

	std::string rg_library = (bam_lib_from_samp == 0 ? rg_iter->Library : rg_iter->Sample);

	// Ensure that there aren't identical read group ids that map to different samples or libraries
	if (rg_ids_to_sample.find(rg_iter->ID) != rg_ids_to_sample.end())
	  if (rg_ids_to_sample[rg_iter->ID].compare(rg_iter->Sample) != 0)
	    printErrorAndDie("Read group id " + rg_iter->ID + " maps to more than one sample");
	if (rg_ids_to_library.find(rg_iter->ID) != rg_ids_to_library.end())
	  if (rg_ids_to_library[rg_iter->ID].compare(rg_library) != 0)
	    printErrorAndDie("Read group id " + rg_iter->ID + " maps to more than one library");

	rg_ids_to_sample[bam_files[i] + rg_iter->ID]  = rg_iter->Sample;
	rg_ids_to_library[bam_files[i] + rg_iter->ID] = rg_library;
	rg_samples.insert(rg_iter->Sample);
	rg_libs.insert(rg_library);
      }
      single_file_reader.Close();
    }
    bam_processor.logger() << "BAMs contain unique read group IDs for "
			   << rg_libs.size()    << " unique libraries and "
			   << rg_samples.size() << " unique samples" << std::endl;
  }

  // Open BAM index files, assuming they're either the same path with a .bai suffix or a path where .bai replaces .bam
  std::vector<std::string> bam_indexes;
  for (unsigned int i = 0; i < bam_files.size(); i++){
    bool have_index      = false;
    std::string bai_file = bam_files[i] + ".bai";
    if (!file_exists(bai_file)){
      int filename_len = bam_files[i].size();
      if (filename_len > 4 && bam_files[i].substr(filename_len-4).compare(".bam") == 0){
	bai_file = bam_files[i].substr(0, filename_len-4) + ".bai";
	if (file_exists(bai_file)){
	  have_index = true;
	  bam_indexes.push_back(bai_file);
	}
      }
    }
    else {
      have_index = true;
      bam_indexes.push_back(bai_file);
    }

    if(!have_index){
      std::stringstream error_msg;
      error_msg << "Unable to find a BAM index file for " << bam_files[i] << "\n"
		<< "Please ensure that each BAM has been sorted by position and indexed using samtools.";
      printErrorAndDie(error_msg.str());
    }
  }
  if (!reader.OpenIndexes(bam_indexes))
    printErrorAndDie("Failed to open one or more BAM index files");

  BamTools::BamWriter bam_pass_writer;
  if (!bam_pass_out_file.empty()){
    BamTools::RefVector ref_vector = reader.GetReferenceData();
    bool file_open = bam_pass_writer.Open(bam_pass_out_file, reader.GetHeaderText(), ref_vector);
    if (!file_open) printErrorAndDie("Failed to open output BAM file for reads used to genotype region");
  }
  BamTools::BamWriter bam_filt_writer;
  if (!bam_filt_out_file.empty()){
    BamTools::RefVector ref_vector = reader.GetReferenceData();
    bool file_open = bam_filt_writer.Open(bam_filt_out_file, reader.GetHeaderText(), ref_vector);
    if (!file_open) printErrorAndDie("Failed to open output BAM file for reads filtered for each region");
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

  if(!str_vcf_out_file.empty()){
    if (!string_ends_with(str_vcf_out_file, ".gz"))
      printErrorAndDie("Path for STR VCF output file must end in .gz as it will be bgzipped");
    bam_processor.set_output_str_vcf(str_vcf_out_file, full_command, rg_samples);
  }

  if (!hap_chr_string.empty()){
    std::vector<std::string> haploid_chroms;
    split_by_delim(hap_chr_string, ',', haploid_chroms);
    for (auto chrom_iter = haploid_chroms.begin(); chrom_iter != haploid_chroms.end(); chrom_iter++)
      bam_processor.add_haploid_chrom(*chrom_iter);
  }

  // Run analysis
  bam_processor.process_regions(reader, region_file, fasta_dir, rg_ids_to_sample, rg_ids_to_library, bam_pass_writer, bam_filt_writer, std::cout, 1000000, chrom);
  bam_processor.finish();

  if (!bam_pass_out_file.empty()) bam_pass_writer.Close();
  if (!bam_filt_out_file.empty()) bam_filt_writer.Close();
  reader.Close();

  total_time = (clock() - total_time)/CLOCKS_PER_SEC;
  bam_processor.logger() << "HipSTR execution finished: Total runtime = " << total_time << " sec" << std::endl;
  return 0;  
}
