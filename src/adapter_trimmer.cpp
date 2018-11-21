#include <assert.h>
#include <time.h>
#include <sstream>

#include "adapter_trimmer.h"
#include "error.h"

void AdapterTrimmer::init(){
  for (auto iter = r1_fw_adapters_.begin(); iter != r1_fw_adapters_.end(); ++iter)
    r1_rc_adapters_.push_back(reverse_complement(*iter));
  for (auto iter = r2_fw_adapters_.begin(); iter != r2_fw_adapters_.end(); ++iter)
    r2_rc_adapters_.push_back(reverse_complement(*iter));

  r1_trimmed_bases_    = r2_trimmed_bases_ = 0;
  r1_trimmed_reads_    = r2_trimmed_reads_ = 0;
  r1_total_reads_      = r2_total_reads_   = 0;
  locus_trimming_time_ = 0;
  total_trimming_time_ = 0;
}

std::string AdapterTrimmer::reverse_complement(const std::string& pattern){
  std::stringstream res;
  for (auto iter = pattern.rbegin(); iter != pattern.rend(); ++iter){
    switch(*iter){
    case 'A': case 'a':
      res << 'T';
      break;
    case 'C': case 'c':
      res << 'G';
      break;
    case 'G': case 'g':
      res << 'C';
      break;
    case 'T': case 't':
      res << 'A';
      break;
    default:
      printErrorAndDie("Invalid character in pattern argument to reverse_complement(): " + *iter);
      break;
    }
  }
  return res.str();
}

void AdapterTrimmer::trim_adapters(BamAlignment& aln){
  // In case trimming isn't actually desired, do nothing
  if (!trim_)
    return;

  // Start the clock
  double start_time = clock();

  if (aln.IsFirstMate()){
    if (aln.IsReverseStrand()){
      // TO DO: Implement
      // Look for adapter near 5' end, scanning from right-to-left
      r1_rc_adapters_;
    }
    else {
      // TO DO: Implement
      // Look for adapter near 3' end, scanning from left-to-right
      r1_fw_adapters_;
    }

    // TO DO: Update these
    r1_trimmed_bases_;
    r1_trimmed_reads_;
    r1_total_reads_++;
  }
  else if (aln.IsSecondMate()){
    if (aln.IsReverseStrand()){
      // TO DO: Implement
      // Look for adapter near 5' end, scanning from right-to-left
      r2_rc_adapters_;
    }
    else {
      // TO DO: Implement
      // Look for adapter near 3' end, scanning from left-to-right
      r2_fw_adapters_;
    }

    // TO DO: Update these
    r2_trimmed_bases_;
    r2_trimmed_reads_;
    r2_total_reads_++;
  }
  else
    assert(false);
  

  // TO DO : Use this to extract the bases for alignment
  //QueryBases();

  // TO DO: Use this to execute the actual trimming
  //aln.TrimNumBases(int left_trim, int right_trim){

  // Turn off the clock
  locus_trimming_time_ += (clock() - start_time)/CLOCKS_PER_SEC;
}


// Minimum overlap between the alignment sequence and adapter sequence
const int AdapterTrimmer::MIN_OVERLAP = 5;

// Initialize various commonly used adapter sequences
const std::string AdapterTrimmer::TRUSEQ_R1_ADAPTER  = "AGATCGGAAGAGCACACGTCTGAAC";
const std::string AdapterTrimmer::TRUSEQ_R2_ADAPTER  = "AGATCGGAAGAGCGTCGTGTAGGGA";
const std::string AdapterTrimmer::NEXTERA_R1_ADAPTER = "CTGTCTCTTATACACATCT";
const std::string AdapterTrimmer::NEXTERA_R2_ADAPTER = "CTGTCTCTTATACACATCT";

