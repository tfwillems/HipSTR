#include <sys/stat.h>

#include "vcf_reader.h"

namespace VCF {

  const std::vector<std::string>& Variant::get_samples(){
    return vcf_reader_->get_samples();
  }

  void Variant::get_genotype(std::string& sample, int& gt_a, int& gt_b){
    int sample_index = vcf_reader_->get_sample_index(sample);
    if (sample_index == -1)
      gt_a = gt_b = -1;
    else {
      gt_a = gt_1_[sample_index];
      gt_b = gt_2_[sample_index];
    }
  }

  bool Variant::sample_call_missing(const std::string& sample){
    int sample_index = vcf_reader_->get_sample_index(sample);
    return (sample_index == -1 ? true : missing_[sample_index]);
  }

  void Variant::extract_alleles(){
    for (int i = 0; i < vcf_record_->n_allele; i++)
      alleles_.push_back(vcf_record_->d.allele[i]);
  }

  void Variant::extract_genotypes(){
    int   mem = 0;
    int* gts_ = NULL;
    std::string GT_KEY = "GT";
    if (bcf_get_format_int32(vcf_header_, vcf_record_, GT_KEY.c_str(), &gts_, &mem) <= 0)
      printErrorAndDie("Failed to extract the genotypes from the VCF record");
    
    missing_.reserve(num_samples_);
    gt_1_.reserve(num_samples_);
    gt_2_.reserve(num_samples_);
    
    int gt_index = 0;
    for (int i = 0; i < num_samples_; i++){
      if (bcf_gt_is_missing(gts_[gt_index]) || bcf_gt_is_missing(gts_[gt_index+1])){
	missing_.push_back(true);
	phased_.push_back(false);
	gt_1_.push_back(-1);
	gt_2_.push_back(-1);
      }
      else {
	missing_.push_back(false);
	phased_.push_back(bcf_gt_is_phased(gts_[gt_index+1]));
	gt_1_.push_back(bcf_gt_allele(gts_[gt_index]));
	gt_2_.push_back(bcf_gt_allele(gts_[gt_index+1]));
      }
      gt_index += 2;
    }
    free(gts_);
  }

void VCFReader::open(std::string& filename){
  const char* cfilename = filename.c_str();
    
  if (bgzf_is_bgzf(cfilename) != 1)
    printErrorAndDie("VCF file is not in a valid bgzipped file. Please ensure that bgzip was used to compress it");
  
  char *fnidx = (char*) calloc(strlen(cfilename) + 5, 1);
  strcat(strcpy(fnidx, cfilename), ".tbi");
  struct stat stat_tbi, stat_vcf;
  stat(fnidx, &stat_tbi);
  stat(cfilename, &stat_vcf);
  if (stat_vcf.st_mtime > stat_tbi.st_mtime)
    printErrorAndDie("The tabix index for the VCF file is older than the VCF itself. Please reindex the VCF with tabix");
  free(fnidx);
  
  if ((vcf_input_ = hts_open(cfilename, "r")) == NULL)
    printErrorAndDie("Failed to open the VCF file");
  if ((tbx_input_ = tbx_index_load(cfilename)) == NULL)
    printErrorAndDie("Failed to open the VCF file's tabix index");
  
  int nseq;
  const char** seq = tbx_seqnames(tbx_input_, &nseq);
  for (int i = 0; i < nseq; i++)
    chroms_.push_back(seq[i]);
  free(seq);
  
  if (chroms_.size() == 0)
    printErrorAndDie("VCF does not contain any chromosomes");
  
  vcf_header_  = bcf_hdr_read(vcf_input_);
  tbx_iter_    = tbx_itr_querys(tbx_input_, chroms_.front().c_str());
  chrom_index_ = 0;
  
  for (int i = 0; i < bcf_hdr_nsamples(vcf_header_); i++){
    samples_.push_back(vcf_header_->samples[i]);
    sample_indices_[vcf_header_->samples[i]] = i;
  }
}

bool VCFReader::get_next_variant(Variant& variant){
  if ((tbx_iter_ != NULL) && tbx_itr_next(vcf_input_, tbx_input_, tbx_iter_, &vcf_line_) >= 0){
    if (vcf_parse(&vcf_line_, vcf_header_, vcf_record_) < 0)
      printErrorAndDie("Failed to parse VCF record");
    variant = Variant(vcf_header_, vcf_record_, this);
    return true;
  }
  
  if (jumped_)
    return false;
  
  while (chrom_index_+1 < chroms_.size()){
    chrom_index_++;
    tbx_itr_destroy(tbx_iter_);
    tbx_iter_ = tbx_itr_querys(tbx_input_, chroms_[chrom_index_].c_str());
    
    if ((tbx_iter_ != NULL) && tbx_itr_next(vcf_input_, tbx_input_, tbx_iter_, &vcf_line_) >= 0){
      if (vcf_parse(&vcf_line_, vcf_header_, vcf_record_) < 0)
	printErrorAndDie("Failed to parse VCF record");
      variant = Variant(vcf_header_, vcf_record_, this);
      return true;
    }
  }
  return false;
}

};
