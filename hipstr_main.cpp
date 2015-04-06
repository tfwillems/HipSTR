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
#include "base_quality.h"
#include "em_stutter_genotyper.h"
#include "error.h"
#include "extract_indels.h"
#include "region.h"
#include "seqio.h"
#include "snp_phasing_quality.h"
#include "snp_tree.h"
#include "stringops.h"

int MAX_EM_ITER         = 100;
double ABS_LL_CONVERGE  = 0.01;  // For EM convergence, new_LL - prev_LL < ABS_LL_CONVERGE
double FRAC_LL_CONVERGE = 0.001; // For EM convergence, -(new_LL-prev_LL)/prev_LL < FRAC_LL_CONVERGE

class SNPBamProcessor : public BamProcessor {
private:
  bool have_vcf;
  vcf::VariantCallFile phased_vcf;
  BaseQuality base_qualities;
  int32_t match_count_, mismatch_count_;

  // Settings controlling EM algorithm convergence
  int max_em_iter_;
  double LL_abs_change_, LL_frac_change_;
  
  // Counters for EM convergence
  int num_em_converge_, num_em_fail_;

public:
  SNPBamProcessor(bool use_lobstr_rg, int max_iter, double LL_abs_change, double LL_frac_change):BamProcessor(use_lobstr_rg){
    have_vcf         = false;
    match_count_     = 0;
    mismatch_count_  = 0;
    max_em_iter_     = max_iter;
    LL_abs_change_   = LL_abs_change;
    LL_frac_change_  = LL_frac_change;
    num_em_converge_ = 0;
    num_em_fail_     = 0;
  }

  SNPBamProcessor(bool use_lobstr_rg, int max_iter, double LL_abs_change, double LL_frac_change, std::string& vcf_file):BamProcessor(use_lobstr_rg){
    set_vcf(vcf_file);
    match_count_    = 0;
    mismatch_count_ = 0;
    max_em_iter_    = max_iter;
    LL_abs_change_  = LL_abs_change;
    LL_frac_change_ = LL_frac_change;
    num_em_converge_ = 0;
    num_em_fail_     = 0;
  }

  void set_vcf(std::string& vcf_file){
    if(!phased_vcf.open(vcf_file))
      printErrorAndDie("Failed to open VCF file");
    have_vcf = true;
  }

  void process_reads(std::vector< std::vector<BamTools::BamAlignment> >& paired_strs_by_rg,
		     std::vector< std::vector<BamTools::BamAlignment> >& mate_pairs_by_rg,
		     std::vector< std::vector<BamTools::BamAlignment> >& unpaired_strs_by_rg,
		     std::vector<std::string>& rg_names, Region& region,
		     std::ostream& out){
    assert(paired_strs_by_rg.size() == mate_pairs_by_rg.size() && paired_strs_by_rg.size() == unpaired_strs_by_rg.size());
    if(paired_strs_by_rg.size() == 0 && unpaired_strs_by_rg.size() == 0)
      return;

    std::vector< std::vector<double> > log_p1s, log_p2s;
    if (have_vcf){
      std::vector<SNPTree*> snp_trees;
      std::map<std::string, unsigned int> sample_indices;      
      if(create_snp_trees(region.chrom(), (region.start() > MAX_MATE_DIST ? region.start()-MAX_MATE_DIST : 1), 
			  region.stop()+MAX_MATE_DIST, phased_vcf, sample_indices, snp_trees)){
	for (unsigned int i = 0; i < paired_strs_by_rg.size(); ++i){
	  assert(sample_indices.find(rg_names[i]) != sample_indices.end());
	  std::vector<double> log_p1, log_p2;
	  SNPTree* snp_tree = snp_trees[sample_indices[rg_names[i]]];
	  calc_het_snp_factors(paired_strs_by_rg[i], mate_pairs_by_rg[i], base_qualities, snp_tree, log_p1, log_p2, match_count_, mismatch_count_);
	  calc_het_snp_factors(unpaired_strs_by_rg[i], base_qualities, snp_tree, log_p1, log_p2, match_count_, mismatch_count_);
	  log_p1s.push_back(log_p1);
	  log_p2s.push_back(log_p2);
	}
      }
      destroy_snp_trees(snp_trees);      
    }

    // Extract STR sizes for each read (if possible) and their associated phasing likelihoods
    std::vector< std::vector<int> > str_bp_lengths(paired_strs_by_rg.size());
    std::vector< std::vector<double> > str_log_p1s(paired_strs_by_rg.size()), str_log_p2s(paired_strs_by_rg.size());
    for (unsigned int i = 0; i < paired_strs_by_rg.size(); i++){
      for (int read_type = 0; read_type < 2; read_type++){
	std::vector<BamTools::BamAlignment>& reads = (read_type == 0 ? paired_strs_by_rg[i] : unpaired_strs_by_rg[i]);
	unsigned int read_index = 0;
	for (unsigned int j = 0; j < reads.size(); ++j, ++read_index){
	  int bp_diff;
	  bool got_size = ExtractCigar(reads[read_index].CigarData, reads[read_index].Position, region.start(), region.stop(), bp_diff);
	  if (got_size){
	    str_bp_lengths[i].push_back(bp_diff);
	    if (log_p1s.size() == 0){
	      str_log_p1s[i].push_back(-1); str_log_p2s[i].push_back(-1); // Assign equal phasing LLs as no SNP info is available
	    }
	    else {
	      str_log_p1s[i].push_back(log_p1s[i][read_index]); str_log_p2s[i].push_back(log_p2s[i][read_index]);
	    }
	  }
	}
      }
    }
	
    // Train stutter model and genotype each sample
    std::cerr << "Building EM stutter genotyper" << std::endl;
    EMStutterGenotyper stutter_genotyper(str_bp_lengths, str_log_p1s, str_log_p2s, rg_names, region.period());
    std::cerr << "Training EM sutter genotyper" << std::endl;
    bool trained = stutter_genotyper.train(max_em_iter_, LL_abs_change_, LL_frac_change_);
    if (trained){
      num_em_converge_++;
      stutter_genotyper.genotype();
      std::cerr << "Learned stutter model: " << (*stutter_genotyper.get_stutter_model()) << std::endl;
    }
    else {
      num_em_fail_++;
      std::cerr << "Stutter model failed to converge" << std::endl;
    }
  }

  void finish(){
    std::cerr << "SNP matching statistics: "   << match_count_     << "\t" << mismatch_count_ << "\n"
	      << "EM convergence statistics: " << num_em_converge_ << "\t" << num_em_fail_ << std::endl;
  }
};


void parse_command_line_args(int argc, char** argv, 
			     std::string& bamfile_string, std::string& bamindex_string, std::string& rg_string,
			     std::string& fasta_dir, std::string& region_file,  std::string& vcf_file, std::string& chrom, 
			     std::string& bam_out_file, BamProcessor& bam_processor){
   if (argc == 1){
     std::cerr << "Usage: HipSTR --bams  <list_of_bams>  --indexes <list_of_bam_indexes> --rgs <list_of_read_groups>" << "\n"
	       << "              --fasta <dir>           --regions <region_file.bed>" << "\n"
	       << "              [--bam-out <spanning_reads.bam>] [--rem-multimaps] [--chrom <chrom>] [--vcf <phased_snp_gts.vcf>]" << "\n\n"
	       << "Required parameters:" << "\n"
	       << "\t" << "--bams          <list_of_bams>        "  << "\t" << "Comma separated list of .bam files"                                                  << "\n"
	       << "\t" << "--indexes       <list_of_bam_indexes> "  << "\t" << "Comma separated list of .bai files in same order as .bam files"                      << "\n"
	       << "\t" << "--fasta         <dir>                 "  << "\t" << "Directory in which FASTA files for each chromosome are located"                      << "\n"
	       << "\t" << "--regions       <region_file.bed>     "  << "\t" << "BED file containing coordinates for each STR region"                                 << "\n" << "\n"
	       << "Optional parameters:" << "\n"
	       << "\t" << "--bam-out       <spanning_reads.bam   "  << "\t" << "Output a BAM file containing the reads spanning each region to the provided file"    << "\n"
	       << "\t" << "--rem-multimaps                       "  << "\t" << "Remove reads that map to multiple locations (Default = False)"                       << "\n"
	       << "\t" << "--max-mate-dist <max_bp>              "  << "\t" << "Remove reads whose mate pair distance is > MAX_BP (Default = " << bam_processor.MAX_MATE_DIST << ")" << "\n"
       	       << "\t" << "--chrom         <chrom>               "  << "\t" << "Only consider STRs on the provided chromosome"                                       << "\n"
	       << "\t" << "--vcf           <phased_snp_gts.vcf>  "  << "\t" << "VCF file containing phased SNP genotypes"                                            << "\n"
	       << "\t" << "--rgs           <list_of_read_groups> "  << "\t" << "Comma separated list of read groups in same order as .bam files. "                   << "\n"
	       << "\t" << "                                      "  << "\t" << "Assign each read the RG tag corresponding to its file. By default, "                 << "\n"
	       << "\t" << "                                      "  << "\t" << "each read must have an RG flag from lobSTR and this is used instead"                 << "\n"
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
  SNPBamProcessor bam_processor(MAX_EM_ITER, ABS_LL_CONVERGE, FRAC_LL_CONVERGE, false);
  std::string bamfile_string= "", bamindex_string="", rg_string="", region_file="", fasta_dir="", chrom="", vcf_file="", bam_out_file="";
  parse_command_line_args(argc, argv, bamfile_string, bamindex_string, rg_string, fasta_dir, region_file, vcf_file, chrom, bam_out_file, bam_processor);
  int num_flank = 0;
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
  if (!rg_string.empty())
    split_by_delim(rg_string, ',', read_groups);
  std::cerr << "Detected " << bam_files.size() << " BAM files" << std::endl;
  if (bam_files.size() != bam_indexes.size())
    printErrorAndDie("Number of .bam and .bai files must match");

  // Open all .bam and .bai files
  BamTools::BamMultiReader reader;
  if (!reader.Open(bam_files)) printErrorAndDie("Failed to open one or more BAM files");
  if (!reader.OpenIndexes(bam_indexes)) printErrorAndDie("Failed to open one or more BAM index files");

  // Construct filename->read group map (if one has been specified)
  std::map<std::string, std::string> file_read_groups;
  if (!rg_string.empty()){
    if(bam_files.size() != read_groups.size())
      printErrorAndDie("Number of .bam and RGs must match");
    for (int i = 0; i < bam_files.size(); i++)
      file_read_groups[bam_files[i]] = read_groups[i];
  }
  else
    bam_processor.set_lobstr_rg_usage(true);
    
  BamTools::BamWriter bam_writer;
  if (!bam_out_file.empty()){
    BamTools::RefVector ref_vector = reader.GetReferenceData();
    bool file_open = bam_writer.Open(bam_out_file, reader.GetHeaderText(), ref_vector);
    if (!file_open) printErrorAndDie("Failed to open output BAM file");
  }

  if (!vcf_file.empty())
    bam_processor.set_vcf(vcf_file);

  // Run analysis
  bam_processor.process_regions(reader, region_file, fasta_dir, file_read_groups, bam_writer, std::cout, 1000);

  bam_processor.finish();

  if (!bam_out_file.empty()) bam_writer.Close();
  reader.Close();
  return 0;  
}
