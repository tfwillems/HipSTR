#ifndef VCF_READER_H_
#define VCF_READER_H_

#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "htslib/htslib/bgzf.h" 
#include "htslib/htslib/tbx.h" 
#include "htslib/htslib/vcf.h" 

#include "error.h"

namespace VCF {

class Variant {
private:
  bcf_hdr_t* vcf_header_;
  bcf1_t*    vcf_record_;
  
  std::vector<std::string> alleles_;

  int num_samples_;
  std::vector<bool> missing_;
  std::vector<bool> phased_;
  std::vector<int> gt_1_, gt_2_;
  
  void extract_alleles();
  void extract_genotypes();

public:
  Variant(){
    vcf_record_ = NULL;
    vcf_header_ = NULL;
  }

  Variant(bcf_hdr_t* vcf_header, bcf1_t* vcf_record){
    vcf_header_  = vcf_header;
    vcf_record_  = vcf_record;
    num_samples_ = bcf_hdr_nsamples(vcf_header_);
    bcf_unpack(vcf_record_, BCF_UN_ALL);
    extract_genotypes();
  }
  
  ~Variant(){ }

  const std::string& get_allele(int allele)    { return alleles_[allele]; }

  bool is_biallelic_snp(){
    if (vcf_record_ != NULL)
      return (vcf_record_->n_allele == 2) && bcf_is_snp(vcf_record_);
    return false;
  }

  int32_t get_position(){
    if (vcf_record_ != NULL)
      return vcf_record_->pos;
    else
      return -1;
  }

  std::string get_chromosome(){
    if (vcf_record_ != NULL)
      return bcf_seqname(vcf_header_, vcf_record_);
    else
      return "";
  }

  int get_format_field_index(const std::string &fieldname) const{
    bcf_fmt_t* fmt = bcf_get_fmt(vcf_header_, vcf_record_, fieldname.c_str());
    return (fmt == NULL ? -1 : fmt->id);
  }

  int get_info_field_index(const std::string &fieldname) const{
    bcf_info_t* info = bcf_get_info(vcf_header_, vcf_record_, fieldname.c_str());
    return (info == NULL ? -1 : info->key);
  }

  bool has_format_field(const std::string& fieldname) const{
    return (bcf_get_fmt(vcf_header_, vcf_record_, fieldname.c_str()) != NULL);
  }

  bool has_info_field(const std::string& fieldname) const{
    return (bcf_get_info(vcf_header_, vcf_record_, fieldname.c_str()) != NULL);
  }

  bool sample_call_phased(int sample_index) const{
    return phased_[sample_index];
  }

  bool sample_call_missing(int sample_index) const{
    return missing_[sample_index];
  }

  void get_genotype(int sample_index, int& gt_a, int& gt_b){
    gt_a = gt_1_[sample_index];
    gt_b = gt_2_[sample_index];
  }
};

class VCFReader {
private:
  htsFile*    vcf_input_;
  bcf_hdr_t*  vcf_header_;
  kstring_t   vcf_line_;
  bcf1_t*     vcf_record_;
  tbx_t*      tbx_input_;
  hts_itr_t*  tbx_iter_;
  bool        jumped_;
  int         chrom_index_;
  std::vector<std::string> samples_;
  std::vector<std::string> chroms_;
  std::map<std::string, int> sample_indices_;

  void open(std::string& filename);

public:
  VCFReader(std::string& filename){
    vcf_input_  = NULL;
    tbx_input_  = NULL;
    tbx_iter_   = NULL;
    jumped_     = false;
    vcf_line_.l = 0;
    vcf_line_.m = 0;
    vcf_line_.s = NULL;
    vcf_record_ = bcf_init();
    open(filename);
  }

  ~VCFReader(){
    if (vcf_input_  != NULL)   ;
    if (vcf_header_ != NULL)   bcf_hdr_destroy(vcf_header_);
    if (tbx_iter_   != NULL)   tbx_itr_destroy(tbx_iter_);
    if (tbx_input_  != NULL)   tbx_destroy(tbx_input_);
    if (vcf_line_.s != NULL)   free(vcf_line_.s);
    bcf_destroy(vcf_record_);
  }

  bool has_sample(std::string& sample){
    return sample_indices_.find(sample) != sample_indices_.end();
  }

  int get_sample_index(std::string& sample){
    auto sample_iter = sample_indices_.find(sample);
    if (sample_iter == sample_indices_.end())
      return -1;
    else
      return sample_iter->second;
  }
  
  bool set_region(const std::string& region){
    tbx_itr_destroy(tbx_iter_);
    tbx_iter_ = tbx_itr_querys(tbx_input_, region.c_str());
    jumped_  = true;
    return tbx_iter_ != NULL;
  }

  bool set_region(const std::string& chrom, int32_t start, int32_t end = 0){
    std::stringstream region_str;
    if (end) region_str << chrom << ":" << start << "-" << end;
    else     region_str << chrom << ":" << start;
    return set_region(region_str.str());
  }

  const std::vector<std::string>& get_samples(){ return samples_; }

  bool get_next_variant(Variant& variant);
};

};
#endif
