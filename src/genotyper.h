#ifndef GENOTYPER_H_
#define GENOTYPER_H_

#include <assert.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "mathops.h"

class Genotyper {
 private:
  // Private unimplemented copy constructor and assignment operator to prevent operations
  Genotyper(const Genotyper& other);
  Genotyper& operator=(const Genotyper& other);

 protected:
  unsigned int num_reads_;    // Total number of reads across all samples
  int num_samples_;           // Total number of samples
  int num_alleles_;           // Number of valid alleles
  double* log_p1_, *log_p2_;  // Log of SNP phasing likelihoods for each read
  int* sample_label_;         // Sample index for each read
  bool haploid_;              // True iff the underlying marker is haploid

  std::vector<std::string> sample_names_;      // List of sample names
  std::map<std::string, int> sample_indices_;  // Mapping from sample name to index

  // Iterates through samples and then through allele_1 and allele_2
  double* log_sample_posteriors_; 

  // Iterates through reads and then alleles by their indices
  double* log_aln_probs_;

  // Total log-likelihoods for each sample
  double* sample_total_LLs_;

  // Total time spent computing posteriors (seconds)
  double total_posterior_time_;

  // Read weights used to calculate posteriors (See calc_log_sample_posteriors function)
  // Used to account for special cases in which both reads in a pair overlap the STR by setting
  // the weight for the second read to zero. Elsewhere, the alignments probabilities for the two reads are summed
  std::vector<int> read_weights_;

  // Convert a list of integers into a string with key|count pairs separated by semicolons
  // e.g. -1,0,-1,2,2,1 will be converted into -1|2;0|1;1|1;2|2
  std::string condense_read_counts(const std::vector<int>& read_diffs) const {
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

  double log_homozygous_prior() const;

  double log_heterozygous_prior() const;

  virtual void init_log_sample_priors(double* log_sample_ptr);

  /* Compute the posteriors for each sample using the haplotype probabilites, stutter model and read weights */
  double calc_log_sample_posteriors(std::vector<int>& read_weights);

  double calc_log_sample_posteriors(){
    return calc_log_sample_posteriors(read_weights_);
  }

  // Determine the genotype associated with each sample based on the current genotype posteriors
  void get_optimal_haplotypes(std::vector< std::pair<int, int> >& gts) const;

 public:
  Genotyper(bool haploid,
	    const std::vector<std::string>& sample_names,
	    const std::vector< std::vector<double> >& log_p1,
	    const std::vector< std::vector<double> >& log_p2){
    assert(log_p1.size() == log_p2.size() && log_p1.size() == sample_names.size());
    num_reads_ = 0;
    for (unsigned int i = 0; i < log_p1.size(); i++)
      num_reads_ += log_p1[i].size();

    num_alleles_      = -1;
    haploid_          = haploid;
    num_samples_      = log_p1.size();
    sample_names_     = sample_names;
    for (unsigned int i = 0; i < sample_names.size(); i++)
      sample_indices_.insert(std::pair<std::string,int>(sample_names[i], i));

    total_posterior_time_  = 0;
    log_p1_                = new double[num_reads_];
    log_p2_                = new double[num_reads_];
    sample_label_          = new int[num_reads_];
    sample_total_LLs_      = new double[num_samples_];
    read_weights_          = std::vector<int>(num_reads_, 1);
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
    log_aln_probs_         = NULL;
  }

  virtual ~Genotyper(){
    delete [] log_p1_;
    delete [] log_p2_;
    delete [] sample_label_;
    delete [] sample_total_LLs_;
    
    if (log_sample_posteriors_ != NULL)
      delete [] log_sample_posteriors_;
    if (log_aln_probs_ != NULL)
      delete [] log_aln_probs_;
  }

  double posterior_time() const { return total_posterior_time_;  }

  static std::string get_vcf_header(const std::string& fasta_path, const std::string& full_command, const std::vector<std::string>& chroms, const std::vector<std::string>& sample_names,
				    bool output_gls, bool output_pls, bool output_phased_gls);

  void calc_PLs(const std::vector<double>& gls, std::vector<int>& pls) const;

  double calc_gl_diff(const std::vector<double>& gls, int gt_a, int gt_b) const;

  void extract_genotypes_and_likelihoods(int num_variants, std::vector<int>& hap_to_allele,
					 std::vector< std::pair<int,int>  >& best_haplotypes,
					 std::vector< std::pair<int,int>  >& best_gts,
					 std::vector<double>& log_phased_posteriors, std::vector<double>& log_unphased_posteriors,
					 bool calc_gls,        std::vector< std::vector<double> >& gls, std::vector<double>& gl_diffs,
					 bool calc_pls,        std::vector< std::vector<int> >& pls,
					 bool calc_phased_gls, std::vector< std::vector<double> >& phased_gls);
};

#endif
