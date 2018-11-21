#ifndef ADAPTER_TRIMMER_H_
#define ADAPTER_TRIMMER_H_

#include <string>
#include <vector>

#include "bam_io.h"

class AdapterTrimmer {
 private:
  std::vector<std::string> r1_fw_adapters_, r2_fw_adapters_;
  std::vector<std::string> r1_rc_adapters_, r2_rc_adapters_;
  bool trim_; // True iff trimming should even be performed

  // Counters to track trimming statistics
  int64_t r1_trimmed_bases_, r2_trimmed_bases_;
  int64_t r1_trimmed_reads_, r2_trimmed_reads_;
  int64_t r1_total_reads_,   r2_total_reads_;

  // Variables to track timing
  double locus_trimming_time_;
  double total_trimming_time_;

  std::string reverse_complement(const std::string& pattern);

  void init();

 public:
  // Minimum overlap between the adapter sequence and the read sequence for trimming to be considered
  const static int MIN_OVERLAP;

  // Commonly used Illumina adapter pairs (see https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html)
  // i)  Applicable to most WGS experiments
  const static std::string TRUSEQ_R1_ADAPTER;
  const static std::string TRUSEQ_R2_ADAPTER;
  // ii) Applicable to most exome/capture experiments
  const static std::string NEXTERA_R1_ADAPTER;
  const static std::string NEXTERA_R2_ADAPTER;

  AdapterTrimmer(){
    // Use both Truseq and Nextera adapters for trimming by default to be comprehensive
    r1_fw_adapters_.push_back(TRUSEQ_R1_ADAPTER);
    r2_fw_adapters_.push_back(TRUSEQ_R2_ADAPTER);
    r1_fw_adapters_.push_back(NEXTERA_R1_ADAPTER);
    r2_fw_adapters_.push_back(NEXTERA_R2_ADAPTER);
    trim_ = true;
    init();
  }

  AdapterTrimmer(const std::string& read1_adapter, const std::string& read2_adapter, bool trim){
    r1_fw_adapters_.push_back(read1_adapter);
    r2_fw_adapters_.push_back(read2_adapter);
    trim_ = trim;
    init();
  }
  
  void mark_new_locus(){
    total_trimming_time_ += locus_trimming_time_;
    locus_trimming_time_  = 0;
  }

  double locus_trimming_time(){ return locus_trimming_time_; }
  double total_trimming_time(){ return total_trimming_time_; }

  void trim_adapters(BamAlignment& aln);
};

#endif
