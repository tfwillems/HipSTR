#include <algorithm>
#include <assert.h>
#include <map>
#include <sstream>

#include "bamtools/include/api/BamAlignment.h"
#include "base_quality.h"
#include "error.h"
#include "mathops.h"
#include "stringops.h"

std::string BaseQuality::average_base_qualities(const std::vector<const std::string*>& qualities){
  assert(qualities.size() > 0);

  // Check that all base quality strings are of the same length
  for (unsigned int i = 0; i < qualities.size(); i++)
    if (qualities[i]->size() != qualities[0]->size())
      printErrorAndDie("All base quality strings must be of the same length when averaging probabilities");

  // Average raw error probabilities for each base and convert
  // to the closest quality score
  std::string avg_qualities(qualities[0]->size(), 'N');
  std::vector<double> log_probs(qualities.size());
  for (unsigned int i = 0; i < qualities[0]->size(); i++){
    for (unsigned int j = 0; j < qualities.size(); j++)
      log_probs[j] = log_prob_error(qualities[j]->at(i));
    double log_mean_prob = log_sum_exp(log_probs) - int_log(qualities.size());
    avg_qualities[i]     = closest_char(log_mean_prob);
  }
  return avg_qualities;
}

std::string BaseQuality::median_base_qualities(const std::vector<const std::string*>& qualities){
  assert(qualities.size() > 0);

  // Check that all base quality strings are of the same length
  for (unsigned int i = 0; i < qualities.size(); i++)
    if (qualities[i]->size() != qualities[0]->size())
      printErrorAndDie("All base quality strings must be of the same length when averaging probabilities");

  std::string median_qualities(qualities[0]->size(), 'N');
  for (unsigned int i = 0; i < qualities[0]->size(); i++){
    std::vector<char> quals;
    for (unsigned int j = 0; j < qualities.size(); j++)
      quals.push_back(qualities[j]->at(i));
    std::sort(quals.begin(), quals.end());
    median_qualities[i] = quals[quals.size()/2];
  }
  return median_qualities;
}

void printBaseCounts(int* counts, std::ostream& out){
  for (unsigned int i = 0; i < 256; i++)
    if (counts[i] != 0)
      out << (char)i << " " << counts[i] << std::endl;
}

void BaseQuality::deduce_quality_encodings(BamTools::BamMultiReader& reader){
  int32_t max_read_count = 1000000; // Number of reads used to determine quality encoding
  int32_t read_count     = 0;
  BamTools::BamAlignment alignment;
  std::map<std::string, int*> qual_counts;
  while (read_count < max_read_count && reader.GetNextAlignment(alignment)){
    read_count++;
  
    size_t colon_pos  = alignment.Name.find_first_of(":");
    std::string id    = (colon_pos != std::string::npos ? alignment.Name.substr(0, colon_pos) : alignment.Name);
    size_t period_pos = alignment.Name.find_first_of(".");
    std::string id_2  = (period_pos != std::string::npos ? id.substr(0, period_pos) : id);

    auto count_iter = qual_counts.find(id_2);
    if (count_iter == qual_counts.end()){
      count_iter = qual_counts.insert(std::pair<std::string, int*>(id_2, new int [256])).first;
      memset(count_iter->second, 0, 256*sizeof(int));
    }

    // Count base qualities
    for (auto qual_iter = alignment.Qualities.begin(); qual_iter != alignment.Qualities.end(); qual_iter++)
      count_iter->second[(int)(*qual_iter)]++;
  }

  int Illumina_18_count = 0, Illumina_15_count = 0, Illumina_13_count = 0, invalid_count = 0;
  for (auto count_iter = qual_counts.begin(); count_iter != qual_counts.end(); count_iter++){
    int* count_ptr = count_iter->second;
    char min_qual  = '}', max_qual = '!';

    // Determine minimum and maximum observed quality scores
    for (int i = 0; i < 256; i++){
      if (count_ptr[i] != 0){
	min_qual = (char)i;
	break;
      }
    }
    for (int i = 255; i >= 0; i--){
      if (count_ptr[i] != 0){
	max_qual = (char)i;
	break;
      }
    }

    // Assess compatibility with different encoding schemes
    if (min_qual >= '!' && max_qual <= 'J'){
      // i) Illumina 1.8+ (! to J)
      Illumina_18_count++;
    }
    else if (min_qual >= 'B' && max_qual <= 'h'){
      // ii) Illumina 1.5 (C to h where B encodes a bad base)
      Illumina_15_count++;
    }
    else if (min_qual >= '@' && max_qual <= 'h'){
      // iii) Illumina 1.3 (@ to h)
      Illumina_13_count++;
    }
    else {
      // Invalid quality score scale (for now)
      invalid_count++;
      std::stringstream msg;
      msg << "WARNING: Invalid quality score for group "  << count_iter->first
	  << " with a minimum quality score of " << min_qual
	  << " and a maximum quality score of "  << max_qual << std::endl
	  << "Only Illumina 1.8, 1.5 and 1.3 score are currently supported" << std::endl
	  << "All reads from this group will be ignored" << std::endl;
      std::cerr << msg.str() << std::endl;
    }
    delete [] (count_iter->second);
  }
  std::cerr << Illumina_18_count << " " << Illumina_15_count << " " << Illumina_13_count << " " << invalid_count << std::endl;
  qual_counts.clear();
}
