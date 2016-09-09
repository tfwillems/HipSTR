#ifndef GENOTYPER_H_
#define GENOTYPER_H_

#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "region.h"

class Genotyper{
 protected:
  Region* region_;            // Locus information
  unsigned int num_reads_;    // Total number of reads across all samples
  int num_samples_;           // Total number of samples
  int num_alleles_;           // Number of valid alleles
  double* log_p1_, *log_p2_;  // Log of SNP phasing likelihoods for each read
  int* sample_label_;         // Sample index for each read
  bool haploid_;              // True iff the underlying marker is haploid

  std::vector<std::string> sample_names_;      // List of sample names
  std::map<std::string, int> sample_indices_;  // Mapping from sample name to index

  // Iterates through allele_1, allele_2 and then samples by their indices
  // Only used if per-allele priors have been specified for each sample
  double* log_allele_priors_;

  // Iterates through allele_1, allele_2 and then samples by their indices
  double* log_sample_posteriors_; 

  // Total log-likelihoods for each sample
  double* sample_total_LLs_;

  // Convert a list of integers into a string with key|count pairs separated by semicolons
  // e.g. -1,0,-1,2,2,1 will be converted into -1|2;0|1;1|1;2|2
  std::string condense_read_counts(std::vector<int>& read_diffs){
    if (read_diffs.size() == 0)
      return ".";
    std::map<int, int> diff_counts;
    for (unsigned int i = 0; i < read_diffs.size(); i++)
      diff_counts[read_diffs[i]]++;
    std::stringstream res;
    for (auto iter = diff_counts.begin(); iter != diff_counts.end(); iter++){
      if (iter != diff_counts.begin())
	res << ";";
      res << iter->first << "|" << iter->second;
    }
    return res.str();
  }

 public:
  Genotyper(Region& region, bool haploid, std::vector<std::string>& sample_names,
	    std::vector< std::vector<double> >& log_p1, std::vector< std::vector<double> >& log_p2){
    assert(log_p1.size() == log_p2.size() && log_p1.size() == sample_names.size());
    num_reads_ = 0;
    for (unsigned int i = 0; i < log_p1.size(); i++)
      num_reads_ += log_p1[i].size();

    num_alleles_  = -1;
    region_       = region.copy();
    haploid_      = haploid;
    num_samples_  = log_p1.size();
    sample_names_ = sample_names;
    for (unsigned int i = 0; i < sample_names.size(); i++)
      sample_indices_.insert(std::pair<std::string,int>(sample_names[i], i));

    log_p1_                = new double[num_reads_];
    log_p2_                = new double[num_reads_];
    sample_label_          = new int[num_reads_];
    sample_total_LLs_      = new double[num_samples_];

    unsigned int read_index = 0;
    for (unsigned int i = 0; i < log_p1.size(); ++i){
      for (unsigned int j = 0; j < log_p1[i].size(); ++j, ++read_index){
	assert(log_p1[i][j] <= 0.0 && log_p2[i][j] <= 0.0);
	log_p1_[read_index]       = log_p1[i][j];
	log_p2_[read_index]       = log_p2[i][j];
	sample_label_[read_index] = i;
      }
    }

    // These data structures need to be allocated once the number of alleles is known
    // within the derived classes
    log_sample_posteriors_ = NULL;
    log_allele_priors_     = NULL;
  }

  ~Genotyper(){
    delete region_;
    delete [] log_p1_;
    delete [] log_p2_;
    delete [] sample_label_;
    delete [] sample_total_LLs_;
    
    if (log_sample_posteriors_ != NULL)
      delete [] log_sample_posteriors_;
    if (log_allele_priors_ != NULL)
      delete [] log_allele_priors_;
  }

  virtual bool genotype(std::string& chrom_seq, std::ostream& logger) = 0;
};

#endif
