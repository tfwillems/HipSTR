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
#include "stutter_model.h"


class EMStutterGenotyper{
 private:
  int num_reads_;   // Total number of reads across all samples
  int num_samples_; // Total number of samples
  int num_alleles_; // Total number of valid alleles
  int motif_len_;   // # bp in STR motif
  
  int* allele_index_ = NULL;                  // Index of each read's STR size
  double* log_p1_    = NULL, *log_p2_ = NULL; // Log of SNP phasing likelihoods for each read
  int* sample_label_ = NULL;                  // Sample index for each read

  StutterModel* stutter_model_ = NULL;

  std::vector<std::string> sample_names_;     // List of sample names
  std::map<std::string, int> sample_indices_; // Mapping from sample name to index

  std::vector<int> bps_per_allele_;   // Size of each STR allele in bps
  std::vector<int> reads_per_sample_;
  double* log_gt_priors_ = NULL;

  // Iterates through allele_1, allele_2 and then samples by their indices
  double* log_sample_posteriors_ = NULL; 
  
  // Iterates through allele_1, allele_2, and then reads and phases 1 or 2 by their indices
  double* log_read_phase_posteriors_ = NULL; 
  
  // Various functions to compute the sum and logsumexp of values and arrays of values
  inline double sum(double* begin, double* end);
  inline double sum(std::vector<double>& vals);
  inline double log_sum_exp(double* begin, double* end);
  inline double log_sum_exp(double log_v1, double log_v2);
  inline double log_sum_exp(double log_v1, double log_v2, double log_v3);
  inline double log_sum_exp(std::vector<double>& log_vals);

  // Initialization functions for the EM algorithm
  void init_log_gt_priors();
  void init_stutter_model();
  
  // Functions for the M step of the EM algorithm
  void recalc_log_gt_priors();
  void recalc_stutter_model();
  
  // Functions for the E step of the EM algorithm
  double recalc_log_sample_posteriors();  
  void recalc_log_read_phase_posteriors();

 public:
  EMStutterGenotyper(std::vector< std::vector<int> >& num_bps, 
		     std::vector< std::vector<double> >& log_p1, 
		     std::vector< std::vector<double> >& log_p2, 
		     std::vector<std::string>& sample_names, int motif_len){
    assert(num_bps.size() == log_p1.size() && num_bps.size() == log_p2.size() && num_bps.size() == sample_names.size());
    num_samples_  = num_bps.size();
    motif_len_    = motif_len;

    sample_names_ = sample_names;
    for (unsigned int i = 0; i < sample_names.size(); i++)
      sample_indices_.insert(std::pair<std::string,int>(sample_names[i], i));

    // Compute total number of reads and set of allele sizes
    std::set<int> allele_sizes;
    num_reads_ = 0;
    for (unsigned int i = 0; i < num_bps.size(); i++){
      assert(num_bps[i].size() == log_p1[i].size() && num_bps[i].size() == log_p2[i].size());
      for (unsigned int j = 0; j < num_bps[i].size(); j++)
	allele_sizes.insert(num_bps[i][j]);
      num_reads_ += num_bps[i].size();
    }

    // Construct a sorted list of allele sizes and a mapping from allele size to allele index
    bps_per_allele_ = std::vector<int>(allele_sizes.begin(), allele_sizes.end());
    std::sort(bps_per_allele_.begin(), bps_per_allele_.end());
    num_alleles_ = bps_per_allele_.size();
    std::map<int, int> allele_indices;
    for (unsigned int i = 0; i < bps_per_allele_.size(); i++)
      allele_indices[bps_per_allele_[i]] = i;

    // Allocate the relevant data structures
    allele_index_ = new int[num_reads_];
    log_p1_       = new double[num_reads_];
    log_p2_       = new double[num_reads_];
    sample_label_ = new int[num_reads_];
    log_gt_priors_             = new double[num_alleles_]; 
    log_sample_posteriors_     = new double[num_alleles_*num_alleles_*num_samples_]; 
    log_read_phase_posteriors_ = new double[num_alleles_*num_alleles_*num_reads_*2]; 

    // Iterate through all reads and store the relevant information
    int read_index = 0;
    for (unsigned int i = 0; i < num_bps.size(); i++){
      reads_per_sample_.push_back(num_bps[i].size());
      for (unsigned int j = 0; j < num_bps[i].size(); ++j, ++read_index){
        allele_index_[read_index] = allele_indices[num_bps[i][j]];
	log_p1_[read_index]       = log_p1[i][j];
	log_p2_[read_index]       = log_p2[i][j];
	sample_label_[read_index] = i;
      }
    }
  }

  ~EMStutterGenotyper(){
    delete [] allele_index_;
    delete [] log_p1_;
    delete [] log_p2_;
    delete [] sample_label_;
    delete [] log_gt_priors_;
    delete [] log_sample_posteriors_;
    delete [] log_read_phase_posteriors_;
    delete stutter_model_;
  }  


  void write_vcf_record(std::string chrom, uint32_t pos, std::vector<std::string>& sample_names, std::ostream& out);
  
  void set_stutter_model(double inframe_geom,  double inframe_up,  double inframe_down,
			 double outframe_geom, double outframe_up, double outframe_down){
    delete stutter_model_;
    stutter_model_ = new StutterModel(inframe_geom,  inframe_up,  inframe_down, outframe_geom, outframe_up, outframe_down, motif_len_);
  }

  void genotype();
  
  bool train(int max_iter, double min_LL_abs_change, double min_LL_frac_change);

  StutterModel* get_stutter_model(){
    if (stutter_model_ == NULL)
      printErrorAndDie("No stutter model has been specified or learned");
    return stutter_model_;
  }
};

#endif
