#ifndef EM_STUTTER_GENOTYPER_H_
#define EM_STUTTER_GENOTYPER_H_

#include <algorithm>
#include <iostream>
#include <map>
#include <math.h>
#include <set>
#include <string>
#include <vector>

#include "error.h"
#include "genotyper.h"
#include "stutter_model.h"

class EMStutterGenotyper: public Genotyper {
 private:
  int motif_len_;                     // Number of base pairs in STR motif
  int* allele_index_;                 // Index of each read's STR size
  StutterModel* stutter_model_;
  std::vector<int> bps_per_allele_;   // Size of each STR allele in bps
  std::vector<int> reads_per_sample_; // Number of reads for each sample
  double* log_gt_priors_;

  bool use_pop_freqs_;

  // Iterates through reads and then allele_1, allele_2, and phase 1 or 2 by their indices
  double* log_read_phase_posteriors_; 

  void calc_hap_aln_probs(double* log_aln_probs);

  void init_log_sample_priors(double* log_sample_ptr);
  
  // Initialization functions for the EM algorithm
  void init_log_gt_priors();
  void init_stutter_model();
  
  // Functions for the M step of the EM algorithm
  void recalc_log_gt_priors();
  void recalc_stutter_model();
  
  // Functions for the E step of the EM algorithm
  void recalc_log_read_phase_posteriors();

 public:
 EMStutterGenotyper(bool haploid, int motif_length,
		    const std::vector< std::vector<int> >& num_bps,
		    const std::vector< std::vector<double> >& log_p1,
		    const std::vector< std::vector<double> >& log_p2,
		    const std::vector<std::string>& sample_names, int ref_allele): Genotyper(haploid, sample_names, log_p1, log_p2){
    assert(num_bps.size() == log_p1.size() && num_bps.size() == log_p2.size() && num_bps.size() == sample_names.size());
    motif_len_     = motif_length;
    use_pop_freqs_ = false;

    // Compute the set of allele sizes
    std::set<int> allele_sizes;
    for (unsigned int i = 0; i < num_bps.size(); i++){
      assert(num_bps[i].size() == log_p1[i].size() && num_bps[i].size() == log_p2[i].size());
      for (unsigned int j = 0; j < num_bps[i].size(); j++)
	allele_sizes.insert(num_bps[i][j]);
    }

    // Remove the reference allele if it's in the list of alleles
    if (allele_sizes.find(ref_allele) != allele_sizes.end())
      allele_sizes.erase(ref_allele);

    // Construct a list of allele sizes (including the reference allele)
    // The reference allele is stored as the first element but the remaining elements are sorted
    bps_per_allele_ = std::vector<int>(allele_sizes.begin(), allele_sizes.end());
    std::sort(bps_per_allele_.begin(), bps_per_allele_.end());
    bps_per_allele_.insert(bps_per_allele_.begin(), ref_allele);

    // Construct a mapping from allele size to allele index
    num_alleles_ = bps_per_allele_.size();
    std::map<int, int> allele_indices;
    for (unsigned int i = 0; i < bps_per_allele_.size(); i++)
      allele_indices[bps_per_allele_[i]] = i;

    // Allocate the relevant data structures
    allele_index_              = new int[num_reads_];
    log_gt_priors_             = new double[num_alleles_]; 
    log_sample_posteriors_     = new double[num_samples_*num_alleles_*num_alleles_];
    log_read_phase_posteriors_ = new double[num_reads_*num_alleles_*num_alleles_*2];
    log_aln_probs_             = new double[num_reads_*num_alleles_];

    // Iterate through all reads and store the relevant information
    unsigned int read_index = 0;
    for (unsigned int i = 0; i < num_bps.size(); i++){
      reads_per_sample_.push_back(num_bps[i].size());
      for (unsigned int j = 0; j < num_bps[i].size(); ++j, ++read_index){
	assert(log_p1[i][j] <= 0.0 && log_p2[i][j] <= 0.0);
        allele_index_[read_index] = allele_indices[num_bps[i][j]];
      }
    }
    assert(read_index == num_reads_);
    stutter_model_     = NULL;
  }

  ~EMStutterGenotyper(){
    delete [] allele_index_;
    delete [] log_gt_priors_;
    delete [] log_read_phase_posteriors_;
    delete stutter_model_;
  }  
  
  bool train(int max_iter, double min_LL_abs_change, double min_LL_frac_change, bool disp_stats, std::ostream& logger);

  StutterModel* get_stutter_model() const {
    if (stutter_model_ == NULL)
      printErrorAndDie("No stutter model has been specified or learned");
    return stutter_model_;
  }
};

#endif
