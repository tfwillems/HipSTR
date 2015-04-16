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

// Exploratory
#include "SeqAlignment/AlignmentOps.h"
#include "SeqAlignment/HaplotypeGenerator.h"
#include "SeqAlignment/Haplotype.h"
#include "SeqAlignment/HapAligner.h"

int MAX_EM_ITER         = 100;
double ABS_LL_CONVERGE  = 0.01;  // For EM convergence, new_LL - prev_LL < ABS_LL_CONVERGE
double FRAC_LL_CONVERGE = 0.001; // For EM convergence, -(new_LL-prev_LL)/prev_LL < FRAC_LL_CONVERGE

class SNPBamProcessor : public BamProcessor {
private:
  bool have_snp_vcf;
  vcf::VariantCallFile phased_snp_vcf;
  BaseQuality base_qualities;
  int32_t match_count_, mismatch_count_;

  // Settings controlling EM algorithm convergence
  int max_em_iter_;
  double LL_abs_change_, LL_frac_change_;
  
  // Counters for EM convergence
  int num_em_converge_, num_em_fail_;

  // Output file for STR genotypes
  bool output_str_gts_;
  std::ofstream str_vcf_;
  std::vector<std::string> samples_to_genotype_;

  // Output file to visualize alignments
  bool output_viz_;
  std::ofstream viz_pdf_;

public:
  SNPBamProcessor(bool use_lobstr_rg, bool check_mate_chroms, 
		  int max_iter, double LL_abs_change, double LL_frac_change):BamProcessor(use_lobstr_rg, check_mate_chroms){
    have_snp_vcf     = false;
    match_count_     = 0;
    mismatch_count_  = 0;
    max_em_iter_     = max_iter;
    LL_abs_change_   = LL_abs_change;
    LL_frac_change_  = LL_frac_change;
    num_em_converge_ = 0;
    num_em_fail_     = 0;
    output_str_gts_  = false;
    output_viz_      = false;
  }

  void set_viz_output(std::string& pdf_file){
    output_viz_ = true;
    viz_pdf_.open(pdf_file, std::ofstream::out);
    if (!viz_pdf_.is_open())
      printErrorAndDie("Failed to open PDF file for alignment visualization");
  }

  void set_output_str_vcf(std::string& vcf_file, std::set<std::string>& samples_to_output){
    output_str_gts_ = true;
    str_vcf_.open(vcf_file, std::ofstream::out);
    if (!str_vcf_.is_open())
      printErrorAndDie("Failed to open VCF file for STR genotypes");

    // Print floats with exactly 3 decimal places
    str_vcf_.precision(3);
    str_vcf_.setf(std::ios::fixed, std::ios::floatfield);

    // Assemble a list of sample names for genotype output
    std::copy(samples_to_output.begin(), samples_to_output.end(), std::back_inserter(samples_to_genotype_));
    std::sort(samples_to_genotype_.begin(), samples_to_genotype_.end());

    // Write VCF header
    EMStutterGenotyper::write_vcf_header(samples_to_genotype_, str_vcf_);
  }

  void set_input_snp_vcf(std::string& vcf_file){
    if(!phased_snp_vcf.open(vcf_file))
      printErrorAndDie("Failed to open input SNP VCF file");
    have_snp_vcf = true;
  }

  void process_reads(std::vector< std::vector<BamTools::BamAlignment> >& paired_strs_by_rg,
		     std::vector< std::vector<BamTools::BamAlignment> >& mate_pairs_by_rg,
		     std::vector< std::vector<BamTools::BamAlignment> >& unpaired_strs_by_rg,
		     std::vector<std::string>& rg_names, Region& region, 
		     std::string& ref_allele, std::string& chrom_seq, std::ostream& out){
    assert(paired_strs_by_rg.size() == mate_pairs_by_rg.size() && paired_strs_by_rg.size() == unpaired_strs_by_rg.size());
    if(paired_strs_by_rg.size() == 0 && unpaired_strs_by_rg.size() == 0)
      return;

    // Exploratory code related to sequence-based genotyping
    if (output_viz_){
      std::map<std::string, std::string> empty_sample_info;
      std::stringstream ss; ss << region.chrom() << ":" << region.start() << "-" << region.stop();

      std::cerr << "Realigning reads" << std::endl;
      std::vector< std::vector<Alignment> > paired, unpaired;
      int total_reads = 0;
      for (unsigned int i = 0; i < paired_strs_by_rg.size(); i++){
	paired.push_back(std::vector<Alignment>());
	for (unsigned int j = 0; j < paired_strs_by_rg[i].size(); j++){
	  Alignment new_alignment;
	  realign(paired_strs_by_rg[i][j], chrom_seq, new_alignment);
	  paired.back().push_back(new_alignment);
	  total_reads++;
	}
	unpaired.push_back(std::vector<Alignment>());
	for (unsigned int j = 0; j < unpaired_strs_by_rg[i].size(); j++){
          Alignment new_alignment;
          realign(unpaired_strs_by_rg[i][j], chrom_seq, new_alignment);
          unpaired.back().push_back(new_alignment);
	  total_reads++;
        }
      }

      // TO DO: Replace parameters with values learned from length-based EM algorithm
      StutterModel stutter_model(0.9, 0.05, 0.05, 0.7, 0.005, 0.005, region.period());

      std::vector<HapBlock*> blocks;
      int max_ref_flank_len = 20;
      Haplotype* haplotype = generate_haplotype(region, max_ref_flank_len, chrom_seq, paired, unpaired, &stutter_model, blocks);
      std::cerr << "Max block sizes: ";
      for (unsigned int i = 0; i < haplotype->num_blocks(); i++)
	std::cerr << haplotype->get_block(i)->max_size() << " ";
      std::cerr << std::endl;
      
      do {
	haplotype->print(std::cerr);

	for (unsigned int i = 0; i < haplotype->num_blocks(); i++){
	  const std::string& block_seq = haplotype->get_seq(i);
	  std::cerr << i << " " << " " << block_seq << " IS_REPEAT? " << (haplotype->get_block(i)->get_repeat_info() != NULL) << std::endl;
	  /*
	  for (unsigned int j = 0; j < block_seq.size(); j++){
	    std::cerr << haplotype->homopolymer_length(i,j) << " ";
	  }
	  std::cerr << std::endl;	  
	  */
	}
      } while (haplotype->next());

      
      std::cerr << "Beginning sequence-specific functions" << std::endl;
      HapAligner hap_aligner(haplotype, region, max_ref_flank_len, &base_qualities, total_reads);
      int read_index = 0;
      for (unsigned int i = 0; i < paired_strs_by_rg.size(); i++){
	hap_aligner.process_reads(paired[i],   read_index); read_index += paired[i].size();
	hap_aligner.process_reads(unpaired[i], read_index); read_index += unpaired[i].size();
      }
      assert(read_index == total_reads);
      

      // Clean up created data structures
      for (int i = 0; i < blocks.size(); i++)
	delete blocks[i];
      blocks.clear();
      delete haplotype;

      return; // Temporary
    }
    // End of exploratory section


    std::vector< std::vector<double> > log_p1s, log_p2s;
    if (have_snp_vcf){
      std::vector<SNPTree*> snp_trees;
      std::map<std::string, unsigned int> sample_indices;      
      if(create_snp_trees(region.chrom(), (region.start() > MAX_MATE_DIST ? region.start()-MAX_MATE_DIST : 1), 
			  region.stop()+MAX_MATE_DIST, phased_snp_vcf, sample_indices, snp_trees)){
	std::set<std::string> bad_samples, good_samples;
	for (unsigned int i = 0; i < paired_strs_by_rg.size(); ++i){
	  if (sample_indices.find(rg_names[i]) != sample_indices.end()){
	    good_samples.insert(rg_names[i]);
	    std::vector<double> log_p1, log_p2;
	    SNPTree* snp_tree = snp_trees[sample_indices[rg_names[i]]];
	    calc_het_snp_factors(paired_strs_by_rg[i], mate_pairs_by_rg[i], base_qualities, snp_tree, log_p1, log_p2, match_count_, mismatch_count_);
	    calc_het_snp_factors(unpaired_strs_by_rg[i], base_qualities, snp_tree, log_p1, log_p2, match_count_, mismatch_count_);
	    log_p1s.push_back(log_p1); log_p2s.push_back(log_p2);
	  }
	  else {
	    std::vector<double> log_p1, log_p2;
	    for (unsigned int j = 0; j < paired_strs_by_rg[i].size()+unpaired_strs_by_rg[i].size(); ++j){
	      log_p1.push_back(0); log_p2.push_back(0); // Assign equal phasing LLs as no SNP info is available
	    }
	    log_p1s.push_back(log_p1); log_p2s.push_back(log_p2);
	    bad_samples.insert(rg_names[i]);
	  }
	}
	std::cerr << "Found VCF info for " << good_samples.size() << " out of " << good_samples.size()+bad_samples.size() << " samples with STR reads" << std::endl;
      }
      else 
	std::cerr << "Warning: Failed to construct SNP trees for " << region.chrom() << ":" << region.start() << "-" << region.stop() << std::endl; 
      destroy_snp_trees(snp_trees);      
    }

    // Extract STR sizes for each read (if possible) and their associated phasing likelihoods
    int snp_power = 0, no_snp_power = 0, phased_samples = 0;

    std::vector< std::vector<int> > str_bp_lengths(paired_strs_by_rg.size());
    std::vector< std::vector<double> > str_log_p1s(paired_strs_by_rg.size()), str_log_p2s(paired_strs_by_rg.size());
    std::vector<std::string> names;
    for (unsigned int i = 0; i < paired_strs_by_rg.size(); i++){
      bool sample_phased = false;
      unsigned int read_index = 0;
      for (int read_type = 0; read_type < 2; read_type++){
	std::vector<BamTools::BamAlignment>& reads = (read_type == 0 ? paired_strs_by_rg[i] : unpaired_strs_by_rg[i]);
	for (unsigned int j = 0; j < reads.size(); ++j, ++read_index){
	  int bp_diff;
	  //bool got_size = ExtractCigar(reads[j].CigarData, reads[j].Position, region.start(), region.stop(), bp_diff);
	  bool got_size = ExtractCigar(reads[j].CigarData, reads[j].Position, region.start()-region.period(), region.stop()+region.period(), bp_diff);
	  if (got_size){
	    if (bp_diff < -(int)(region.stop()-region.start()+1)) {
	      std::cerr << "WARNING: Excluding read with bp difference greater than reference allele: " << reads[j].Name << std::endl;
	      continue;
	    }
	    str_bp_lengths[i].push_back(bp_diff);
	    if (log_p1s.size() == 0){
	      str_log_p1s[i].push_back(0); str_log_p2s[i].push_back(0); // Assign equal phasing LLs as no SNP info is available
	    }
	    else {
	      str_log_p1s[i].push_back(log_p1s[i][read_index]); str_log_p2s[i].push_back(log_p2s[i][read_index]);
	    }
	    if (str_log_p1s[i].back() != str_log_p2s[i].back()){
	      snp_power++;
	      sample_phased = true;
	    }
	    else
	      no_snp_power++;
	  }
	}
      }
      if (sample_phased) phased_samples++;
    }
    std::cout << "Phased SNPs add info for " << snp_power << " out of " << snp_power+no_snp_power << " reads" 
	      << " and " << phased_samples << " samples" << std::endl;
	
    // Train stutter model and genotype each sample
    std::cerr << "Building EM stutter genotyper" << std::endl;
    EMStutterGenotyper stutter_genotyper(region.chrom(), region.start(), region.stop(), str_bp_lengths, str_log_p1s, str_log_p2s, rg_names, region.period(), 0);

    std::cerr << "Training EM sutter genotyper" << std::endl;
    bool trained = stutter_genotyper.train(max_em_iter_, LL_abs_change_, LL_frac_change_);
    if (trained){
      num_em_converge_++;
      std::cerr << "Learned stutter model: " << *(stutter_genotyper.get_stutter_model()) << std::endl;
      bool use_pop_freqs = false;
      stutter_genotyper.genotype(use_pop_freqs);

      if (output_str_gts_)
	stutter_genotyper.write_vcf_record(ref_allele, samples_to_genotype_, str_vcf_);
    }
    else {
      num_em_fail_++;
      std::cerr << "Stutter model training failed for locus " << region.chrom() << ":" << region.start() << "-" << region.stop() 
		<< " with " << snp_power+no_snp_power << " informative reads" << std::endl;
    }
  }

  void finish(){
    if (output_str_gts_)
      str_vcf_.close();

    if (output_viz_)
      viz_pdf_.close();

    std::cerr << "SNP matching statistics: "   << match_count_     << "\t" << mismatch_count_ << "\n"
	      << "EM convergence statistics: " << num_em_converge_ << "\t" << num_em_fail_ << std::endl;
  }
};


void parse_command_line_args(int argc, char** argv, 
			     std::string& bamfile_string, std::string& bamindex_string, std::string& rg_string,
			     std::string& fasta_dir, std::string& region_file,  std::string& vcf_file, std::string& chrom, 
			     std::string& bam_out_file, std::string& str_vcf_out_file, std::string& viz_out_file, BamProcessor& bam_processor){
  int def_mdist = bam_processor.MAX_MATE_DIST;
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
	      << "\t" << "--str-vcf       <str_gts.vcf>         "  << "\t" << "Output a VCF file containing phased STR genotypes"                                   << "\n"
	      << "\t" << "--viz-out       <viz.pdf>             "  << "\t" << "Output a PDF file visualizing the alignments for each locus"                         << "\n"
	      << "\t" << "--rem-multimaps                       "  << "\t" << "Remove reads that map to multiple locations (Default = False)"                       << "\n"
	      << "\t" << "--max-mate-dist <max_bp>              "  << "\t" << "Remove reads whose mate pair distance is > MAX_BP (Default = " << def_mdist << ")"   << "\n"
	      << "\t" << "--chrom         <chrom>               "  << "\t" << "Only consider STRs on the provided chromosome"                                       << "\n"
	      << "\t" << "--snp-vcf       <phased_snp_gts.vcf>  "  << "\t" << "Input VCF file containing phased SNP genotypes"                                      << "\n"
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
    {"str-vcf",         required_argument, 0, 'o'},
    {"viz-out",         required_argument, 0, 'p'},
    {"regions",         required_argument, 0, 'r'},
    {"snp-vcf",         required_argument, 0, 'v'},
    {"bam-out",         required_argument, 0, 'w'},
    {"rem-multimaps",   no_argument, &(bam_processor.REMOVE_MULTIMAPPERS), 1},
    {0, 0, 0, 0}
  };

  int c;
  while (true){
    int option_index = 0;
    c = getopt_long(argc, argv, "b:c:d:f:g:i:o:r:v:w:", long_options, &option_index);
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
    case 'o':
      str_vcf_out_file = std::string(optarg);
      break;
    case 'p':
      viz_out_file = std::string(optarg);
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
  bool check_mate_chroms = false;
  SNPBamProcessor bam_processor(false, check_mate_chroms, MAX_EM_ITER, ABS_LL_CONVERGE, FRAC_LL_CONVERGE);
  std::string bamfile_string= "", bamindex_string="", rg_string="", region_file="", fasta_dir="", chrom="", vcf_file="";
  std::string bam_out_file="", str_vcf_out_file="", viz_out_file="";
  parse_command_line_args(argc, argv, bamfile_string, bamindex_string, rg_string, fasta_dir, region_file, vcf_file, chrom, 
			  bam_out_file, str_vcf_out_file, viz_out_file, bam_processor);
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
  // and determine the list of samples of interest based on either
  // the specified names or the RG tags in the BAM headers
  std::set<std::string> rg_samples;
  std::map<std::string, std::string> file_read_groups;
  if (!rg_string.empty()){
    if(bam_files.size() != read_groups.size())
      printErrorAndDie("Number of .bam and RGs must match");
    for (int i = 0; i < bam_files.size(); i++){
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

  if (!vcf_file.empty())
    bam_processor.set_input_snp_vcf(vcf_file);
  if(!str_vcf_out_file.empty())
    bam_processor.set_output_str_vcf(str_vcf_out_file, rg_samples);
  if(!viz_out_file.empty())
    bam_processor.set_viz_output(viz_out_file);

  // Run analysis
  bam_processor.process_regions(reader, region_file, fasta_dir, file_read_groups, bam_writer, std::cout, 1000);

  bam_processor.finish();

  if (!bam_out_file.empty()) bam_writer.Close();
  reader.Close();
  return 0;  
}
