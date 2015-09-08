#include <assert.h>
#include <iostream>
#include <string>

#include "../vcflib/src/Variant.h"
#include "../fastahack/Fasta.h"

#include "../bgzf_streams.h"
#include "../error.h"
#include "../seqio.h"
#include "../stringops.h"

std::string get_qual_string(int read_length){
  return std::string(read_length, 'a');
}

int main(int argc, char** argv){
  // Input arguments
  std::string vcf_file  = "/data/twillems/STR_datasets/illumina_200x/init_candidates.vcf.gz";
  int read_length       = 100;
  std::string fasta_dir = "/data/dbase/human/hg19/chromFa/";
  std::string fastq_out = "reads_for_candidates.gz";

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

    //std::cerr << left_flank << " " << variant.ref << " " << right_flank << std::endl;
    //std::cerr << uppercase(chrom_seq.substr(lstart, 250)) << std::endl << std::endl;

    int gt = 0;
    std::string qual_string = get_qual_string(read_length);
    for (auto iter = variant.alleles.begin(); iter != variant.alleles.end(); iter++){
      std::string haplotype = left_flank + *iter + right_flank;
      int32_t lflank_len    = read_length - iter->size();
      int32_t pos           = (int)left_flank.size() - lflank_len;
      while (lflank_len > 0){
	fastq_writer << "@" << var_chrom << "_" << variant.position << "_" << gt << "_" << read_count << "\n"
		     << haplotype.substr(pos, read_length) << "\n"
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

