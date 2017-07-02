#ifndef VCF_WRITER_H_
#define VCF_WRITER_H_

#include <algorithm>
#include <string>

#include "bgzf_streams.h"
#include "error.h"

class RecordTuple {
 private:
  int32_t pos_;
  std::string text_;

 public:  
 RecordTuple(int32_t pos, const std::string& text) : pos_(pos), text_(text) {}
  
  int32_t pos()            { return pos_;  }
  const std::string& text(){ return text_; }
};

bool tuple_comparator(RecordTuple* r1, RecordTuple* r2);

class VCFWriter {
 private:
  bgzfostream str_vcf_;
  bool open_;

  std::string chrom_;
  std::vector<RecordTuple*> record_heap_;

  // We assume regions are processed in sorted order, but that regions
  // can end up with start positions minus this amount of padding at the most
  int32_t MAX_RECORD_PAD;

  void write_all_records(){
    while (!record_heap_.empty()){
      std::pop_heap(record_heap_.begin(), record_heap_.end(), tuple_comparator);
      RecordTuple* best = record_heap_.back(); record_heap_.pop_back();
      str_vcf_ << best->text() << std::endl;
      delete best;
    }
  }

  // Private unimplemented copy constructor and assignment operator to prevent operations
  VCFWriter(const VCFWriter& other);
  VCFWriter& operator=(const VCFWriter& other);

 public:
  VCFWriter(){
    open_          = false;
    MAX_RECORD_PAD = 50;
    chrom_         = "";
  }

  ~VCFWriter(){
    write_all_records();
  }

  bool is_open() const { return open_; }

  void open(const std::string& vcf_file){
    if (open_)
      printErrorAndDie("Cannot reopen an open VCFWriter");
    open_ = true;
    str_vcf_.open(vcf_file.c_str());
  }

  void write_header(const std::string& header_text){
    if (!open_)
      printErrorAndDie("Cannot invoke write_header() on a non-open VCFWriter");
    str_vcf_ << header_text;
  }

  void add_vcf_record(const std::string& chrom, int32_t record_pos, const std::string& record_text);

  void close(){
    write_all_records();
    open_ = false;
    str_vcf_.close();
  }
};

#endif
