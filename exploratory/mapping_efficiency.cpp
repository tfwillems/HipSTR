#include <assert.h>
#include <getopt.h>
#include <iostream>
#include <random>
#include <string>

#include "../vcflib/src/Variant.h"
#include "../fastahack/Fasta.h"

#include "../bgzf_streams.h"
#include "../error.h"
#include "../seqio.h"
#include "../stringops.h"

bool file_exists(std::string path){
  return (access(path.c_str(), F_OK) != -1);
}

void parse_command_line_args(int argc, char** argv, std::string& vcf, std::string& fasta, std::string& out, int& read_length){
   if (argc == 1){
     std::cerr << "Usage: Mapper --vcf <str_gts.vcf.gz> --fasta <fasta_dir> --out <paired_reads.fq.gz> --read-length <num_bps>"                     << "\n\n"
	       << "\t" << "--fasta       <fasta_dir>           " << "\t" << "Directory in which FASTA files for each chromosome are located"        << "\n"
	       << "\t" << "--out         <paired_reads.fq.gz>  " << "\t" << "Output path where bgzipped FASTQ for paired-end reads will be written" << "\n"
	       << "\t" << "--read-length <num_bps>             " << "\t" << "Base pair length of simulated reads"                                   << "\n"
	       << "\t" << "--vcf         <str_gts.vcf.gz>      " << "\t" << "Bgzipped input VCF file containing STR alleles"                        << "\n"
	       << "\t" << "                                    " << "\t" << " for which STR reads will be simulated"                                << "\n\n";
     exit(0);
  }

  static struct option long_options[] = {
    {"fasta",        required_argument, 0, 'f'},
    {"vcf",          required_argument, 0, 'i'},
    {"read-length",  required_argument, 0, 'l'},
    {"out",          required_argument, 0, 'o'},
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
    case 'f':
      fasta = std::string(optarg);
      break;
    case 'i':
      vcf = std::string(optarg);
      break;
    case 'l':
      read_length = atoi(optarg);
      break;
    case 'o':
      out = std::string(optarg);
      break;
    case '?':
      printErrorAndDie("Unrecognized command line option");
      break;
    default:
      abort();
      break;
    }
  }

  if (vcf.empty())   printErrorAndDie("Argument --vcf is required");
  if (fasta.empty()) printErrorAndDie("Argument --fasta is required");
  if (out.empty())   printErrorAndDie("Argument --out is required");
  if (!string_ends_with(vcf, ".gz")) printErrorAndDie("Input VCF file must be bgzipped (and end in .gz)");
  if (!string_ends_with(out, ".gz")) printErrorAndDie("Output file path must end in .gz");
  if (!file_exists(vcf + ".tbi"))    printErrorAndDie("No .tbi index found for the VCF file. Please index using tabix and rerun");
  if (read_length <= 0) printErrorAndDie("Argument --read-length is required and must be a positive integer");
}

std::string get_qual_string(int read_length){
  return std::string(read_length, 'a');
}

int main(int argc, char** argv){
  std::string vcf_file = "", fasta_dir = "", fastq_out = "";
  int read_length = -1;
  parse_command_line_args(argc, argv, vcf_file, fasta_dir, fastq_out, read_length);

  bgzfostream fastq_writer;
  fastq_writer.open(fastq_out.c_str());
  vcflib::VariantCallFile vcf;
  if(!vcf.open(vcf_file))
    printErrorAndDie("Failed to open input VCF file: " + vcf_file);

  FastaReference* fasta_ref = NULL;
  if (is_file(fasta_dir)){
    fasta_ref = new FastaReference();
    fasta_ref->open(fasta_dir);
  }

  std::string chrom_seq  = "";
  std::string prev_chrom = "";
  vcflib::Variant variant(vcf);
  int64_t read_count = 0;
  while (vcf.getNextVariant(variant)){
    std::string var_chrom = variant.sequenceName;

    // Load new FASTA sequence if necessary
    if (var_chrom.compare(prev_chrom) != 0){
     if (fasta_ref != NULL)
	chrom_seq = fasta_ref->getSequence(var_chrom);
      else
	readFastaFromDir(var_chrom+".fa", fasta_dir, chrom_seq);
      prev_chrom = var_chrom;
      assert(chrom_seq.size() != 0);
    }

    int32_t str_start       = variant.position-1;
    int32_t str_stop        = variant.position-1 + variant.ref.size()-1; // Inclusive stop coordinate
    int32_t lstart          = std::max(0, str_start-read_length);
    int32_t rend            = std::min((int32_t)(chrom_seq.size()-1), str_stop+read_length);
    std::string left_flank  = uppercase(chrom_seq.substr(lstart, str_start-lstart));
    std::string right_flank = uppercase(chrom_seq.substr(str_stop+1, rend-str_stop));

    std::default_random_engine generator;
    std::normal_distribution<double> mate_dist_distribution(500, 150);
    std::binomial_distribution<int> mate_dir_distribution;

    int gt = 0;
    std::string qual_string = get_qual_string(read_length);
    for (auto iter = variant.alleles.begin(); iter != variant.alleles.end(); iter++){
      if (iter->compare(".") == 0)
	continue;

      std::string haplotype = left_flank + *iter + right_flank;
      int32_t lflank_len    = read_length - iter->size();
      int32_t pos           = (int)left_flank.size() - lflank_len;
      while (lflank_len > 0){
	// Write out the STR read information
	fastq_writer << "@" << var_chrom << "_" << variant.position << "_" << gt << "_" << read_count << "\n"
		     << haplotype.substr(pos, read_length) << "\n"
		     << "+\n" 
		     << qual_string << "\n";


	// Determine a suitable mate pair position
	int32_t mate_pos = -1;
	while (mate_pos < 0 || mate_pos+read_length >= chrom_seq.size()){
	  int32_t mate_dir = (mate_dir_distribution(generator) == 0 ? -1 : 1);
	  mate_pos = lstart + pos + mate_dir*((int32_t)mate_dist_distribution(generator));
	}

	// Write out the mate pair information
	fastq_writer << "@" << var_chrom << "_" << variant.position << "_" << gt << "_" << read_count << "\n"
		     << uppercase(chrom_seq.substr(mate_pos, read_length)) << "\n"
		     << "+\n"
		     << qual_string << "\n";

	lflank_len--;
	pos++;
	read_count++;
      }
      gt++;
    }
  }

  if (fasta_ref != NULL)
    delete fasta_ref;
}

