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
  int64_t locus_r1_trimmed_bases_, locus_r2_trimmed_bases_;
  int64_t locus_r1_trimmed_reads_, locus_r2_trimmed_reads_;
  int64_t locus_r1_total_reads_,   locus_r2_total_reads_;

  int64_t total_r1_trimmed_bases_, total_r2_trimmed_bases_;
  int64_t total_r1_trimmed_reads_, total_r2_trimmed_reads_;
  int64_t total_r1_total_reads_,   total_r2_total_reads_;

  // Variables to track trimming timing
  double locus_trimming_time_;
  double total_trimming_time_;

  std::string reverse_complement(const std::string& pattern);
  void init();

  /* Identifies the rightmost index in the alignment's sequence that matches any adapter sequence with at most 1 mismatch
     Bases in and to the left of the index are then discarded
     Configurations are only considered if both
     i)  the adapter is fully within the alignment sequence or has a left-side overhang
     ii) the adapter and alignment sequence overlap by at least MIN_OVERLAP bases
     Left-side overhangs are not penalized as mismatches
     Returns the number of trimmed bases or 0 if no bases were trimmed */
  int64_t trim_five_prime(BamAlignment& aln, const std::vector<std::string>& adapters);

  
  /* Identifies the leftmost index in the alignment's sequence that matches any adapter sequence with at most 1 mismatch
     Bases in and to the right of the index are then discarded
     Configurations are only considered if both
     i)  the adapter is fully within the alignment sequence or has a right-side overhang
     ii) the adapter and alignment sequence overlap by at least MIN_OVERLAP bases
     Right-side overhangs are not penalized as mismatches
     Returns the number of trimmed bases or 0 if no bases were trimmed */
  int64_t trim_three_prime(BamAlignment& aln, const std::vector<std::string>& adapters);

  // Private unimplemented copy constructor and assignment operator to prevent operations
  AdapterTrimmer(const AdapterTrimmer& other);
  AdapterTrimmer& operator=(const AdapterTrimmer& other);
  
 public:
  // Minimum overlap between the adapter sequence and the read sequence for trimming to be considered
  const static int    MIN_OVERLAP;
  const static double MAX_ERROR_RATE;

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
  
  /* Save the current locus' statistics and make room for the next locus */
  void mark_new_locus(){
    total_r1_trimmed_bases_ += locus_r1_trimmed_bases_;
    total_r2_trimmed_bases_ += locus_r2_trimmed_bases_;
    total_r1_trimmed_reads_ += locus_r1_trimmed_reads_;
    total_r2_trimmed_reads_ += locus_r2_trimmed_reads_;
    total_r1_total_reads_   += locus_r1_total_reads_;
    total_r2_total_reads_   += locus_r2_total_reads_;
    locus_r1_trimmed_bases_  = locus_r2_trimmed_bases_ = 0;
    locus_r1_trimmed_reads_  = locus_r2_trimmed_reads_ = 0;
    locus_r1_total_reads_    = locus_r2_total_reads_   = 0;

    total_trimming_time_    += locus_trimming_time_;
    locus_trimming_time_     = 0;
  }

  std::string get_trimming_stats_msg();

  double locus_trimming_time(){ return locus_trimming_time_; }
  double total_trimming_time(){ return total_trimming_time_; }

  /* Identify any exact or 1 mismatch adapter sequences present in the alignment's sequence
     Uses information about the alignment's orientation to check for 5' or 3' adapters (as reverse alignments have already been reverse complemented)
     Removes the adapter sequence from the alignment's bases and modifies the alignment properties accordingly
  */
  void trim_adapters(BamAlignment& aln);
};

#endif
