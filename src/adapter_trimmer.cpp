#include <assert.h>
#include <time.h>
#include <sstream>

#include "adapter_trimmer.h"
#include "error.h"
#include "zalgorithm.h"

void AdapterTrimmer::init(){
  assert(MIN_OVERLAP > 0);
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

int64_t AdapterTrimmer::trim_five_prime(BamAlignment& aln, const std::vector<std::string>& adapters){
  const std::string& bases = aln.QueryBases();
  const int read_length    = bases.size();

  // Identify the rightmost trimming index across all adapters
  int trim_index = -1;
  std::vector<int> prefix_matches, suffix_matches;
  for (auto adapter_iter = adapters.begin(); adapter_iter != adapters.end(); ++adapter_iter){
    const int adapter_length = adapter_iter->size();
    std::string query        = std::string(adapter_length, '*') + bases;
    ZAlgorithm::GetPrefixMatchCounts(*adapter_iter, query, prefix_matches);
    ZAlgorithm::GetSuffixMatchCounts(*adapter_iter, query, suffix_matches);
    for (int index = read_length - 1 + adapter_length; index >= MIN_OVERLAP - 1 + adapter_length; --index){
      // Check if a full match or 1 mismatch configuration exists      
      bool valid = (suffix_matches[index] == adapter_length);
      valid     |= ((suffix_matches[index] + 1 + prefix_matches[index-adapter_length+1]) == adapter_length);
      if (valid){
	if (index-adapter_length > trim_index)
	  trim_index = index-adapter_length;
	break;
      }
    }
  }

  // Trim the bases
  if (trim_index >= 0)
    aln.TrimNumBases(trim_index+1, 0);

  return trim_index+1;
}

int64_t AdapterTrimmer::trim_three_prime(BamAlignment& aln, const std::vector<std::string>& adapters){
  const std::string& bases = aln.QueryBases();
  const int read_length    = bases.size();

  // Identify the leftmost trimming index across all adapters
  int trim_index = read_length;
  std::vector<int> prefix_matches, suffix_matches;
  for (auto adapter_iter = adapters.begin(); adapter_iter != adapters.end(); ++adapter_iter){
    const int adapter_length = adapter_iter->size();
    std::string query        = bases + std::string(adapter_length, '*');
    ZAlgorithm::GetPrefixMatchCounts(*adapter_iter, query, prefix_matches);
    ZAlgorithm::GetSuffixMatchCounts(*adapter_iter, query, suffix_matches);
    for (int index = 0; index <= read_length-MIN_OVERLAP; ++index){
      // Check if a full match or 1 mismatch configuration exists
      bool valid = (prefix_matches[index] == adapter_length);
      valid     |= ((prefix_matches[index] + 1 + suffix_matches[index+adapter_length-1]) == adapter_length);
      if (valid){
	if (index < trim_index)
	  trim_index = index;
	break;
      }
    }
  }
  
  // Trim the bases
  if (trim_index < read_length)
    aln.TrimNumBases(0, read_length-trim_index);

  return read_length-trim_index;
}

void AdapterTrimmer::trim_adapters(BamAlignment& aln){
  // In case trimming isn't actually desired, do nothing
  if (!trim_) return;

  double start_time = clock(); // Start the clock
  int64_t num_trim;
  if (aln.IsFirstMate()){
    if (aln.IsReverseStrand())
      num_trim = trim_five_prime(aln, r1_rc_adapters_);
    else
      num_trim = trim_three_prime(aln, r1_fw_adapters_);

    r1_trimmed_bases_ += num_trim;
    r1_trimmed_reads_ += (num_trim > 0 ? 1 : 0);
    r1_total_reads_++;
  }
  else if (aln.IsSecondMate()){
    if (aln.IsReverseStrand())
      num_trim = trim_five_prime(aln, r2_rc_adapters_);
    else
      num_trim = trim_three_prime(aln, r2_fw_adapters_);

    r2_trimmed_bases_ += num_trim;
    r2_trimmed_reads_ += (num_trim > 0 ? 1 : 0);
    r2_total_reads_++;
  }
  else
    assert(false);
  locus_trimming_time_ += (clock() - start_time)/CLOCKS_PER_SEC; // Turn off the clock
}

// Minimum overlap between the alignment sequence and adapter sequence
const int AdapterTrimmer::MIN_OVERLAP = 5;

// Initialize various commonly used adapter sequences. 
// Don't use the full sequences as we only allow for 1 mismatch
// Checking for full matches might decrease sensitivity but should have little effect on specificity
const std::string AdapterTrimmer::TRUSEQ_R1_ADAPTER  = "AGATCGGAAGAGCAC";
const std::string AdapterTrimmer::TRUSEQ_R2_ADAPTER  = "AGATCGGAAGAGCGT";
const std::string AdapterTrimmer::NEXTERA_R1_ADAPTER = "CTGTCTCTTATACAC";
const std::string AdapterTrimmer::NEXTERA_R2_ADAPTER = "CTGTCTCTTATACAC";
