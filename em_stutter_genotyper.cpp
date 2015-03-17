#include <algorithm>
#include <cfloat>

#include "em_stutter_genotyper.h"

inline double EMStutterGenotyper::sum(double* begin, double* end){
  double total = 0.0;
  for (double* iter = begin; iter != end; iter++)
    total += *iter;
  return total;
}

inline double EMStutterGenotyper::sum(std::vector<double>& vals){
  double total = 0.0;
  for (auto iter = vals.begin(); iter != vals.end(); iter++)
    total += *iter;
  return total;
}

inline double EMStutterGenotyper::log_sum_exp(double* begin, double* end){
  double max_val = *std::max_element(begin, end);
  double total   = 0.0;
  for (double* iter = begin; iter != end; iter++)
    total += exp(*iter - max_val);
  return max_val + log(total);
}

inline double EMStutterGenotyper::log_sum_exp(double log_v1, double log_v2){
  if (log_v1 > log_v2)
    return log_v1 + log(1 + exp(log_v2-log_v1));
  else
    return log_v2 + log(1 + exp(log_v1-log_v2));
}

inline double EMStutterGenotyper::log_sum_exp(double log_v1, double log_v2, double log_v3){
  double max_val = std::max(std::max(log_v1, log_v2), log_v3);
  return max_val + log(exp(log_v1-max_val) + exp(log_v2-max_val) + exp(log_v3-max_val));
}

inline double EMStutterGenotyper::log_sum_exp(std::vector<double>& log_vals){
  double max_val = *std::max_element(log_vals.begin(), log_vals.end());
  double total   = 0;
  for (auto iter = log_vals.begin(); iter != log_vals.end(); iter++)
    total += exp(*iter - max_val);
  return max_val + log(total);
}

void EMStutterGenotyper::init_log_gt_priors(){
  std::fill(log_gt_priors_, log_gt_priors_+num_alleles_, 1); // Use 1 sample pseudocount                                                                                  
  for (int i = 0; i < num_reads_; i++)
    log_gt_priors_[allele_index_[i]] += 1.0/reads_per_sample_[sample_label_[i]];
  double log_total = log(sum(log_gt_priors_, log_gt_priors_+num_alleles_));
  for (int i = 0; i < num_alleles_; i++)
    log_gt_priors_[i] = log(log_gt_priors_[i]) - log_total;
}

void EMStutterGenotyper::recalc_log_gt_priors(){
  std::fill(log_gt_priors_, log_gt_priors_+num_alleles_, 0);
  for (int i = 0; i < num_alleles_; i++)
    log_gt_priors_[i] = log_sum_exp(log_sample_posteriors_[i*num_alleles_], log_sample_posteriors_[(i+1)*num_alleles_]);
  double log_total = log_sum_exp(log_gt_priors_, log_gt_priors_+num_alleles_);
  for (int i = 0; i < num_alleles_; i++)
    log_gt_priors_[i] -= log_total;
}

void EMStutterGenotyper::init_stutter_model(){
  delete stutter_model_;
  stutter_model_ = new StutterModel(0.9, 0.1, 0.1, 0.8, 0.01, 0.01, motif_len_);
}
  
void EMStutterGenotyper::recalc_stutter_model(){
  std::vector<double> in_log_up,  in_log_down,  in_log_eq, in_log_diffs; // In-frame values
  std::vector<double> out_log_up, out_log_down, out_log_diffs;           // Out-of-frame values
  
  // Add various pseudocounts such that p_geom < 1 for both in-frame and out-of-frame stutter models
  in_log_up.push_back(0.0);  in_log_down.push_back(0.0);  in_log_diffs.push_back(0.0);  in_log_diffs.push_back(log(2));
  out_log_up.push_back(0.0); out_log_down.push_back(0.0); out_log_diffs.push_back(0.0); out_log_diffs.push_back(log(2));

  double* log_posterior_ptr = log_sample_posteriors_;
  double* log_phase_ptr     = log_read_phase_posteriors_;
  for (int index_1 = 0; index_1 < num_alleles_; ++index_1){
    for (int index_2 = 0; index_2 < num_alleles_; ++index_2){
      for (int read_index = 0; read_index < num_reads_; ++read_index){
	double log_gt_posterior = log_posterior_ptr[sample_label_[read_index]];
	for (int phase = 0; phase < 2; ++phase, ++log_phase_ptr){
	  int gt_index = (phase == 0 ? index_1 : index_2);
	  int bp_diff  = bps_per_allele_[allele_index_[read_index]] - bps_per_allele_[gt_index];
	  
	  if (allele_index_[read_index] == gt_index)
	    in_log_eq.push_back(log_gt_posterior + *log_phase_ptr);
	  else {
	    if (bp_diff % motif_len_ != 0){
	      int eff_diff = bp_diff - bp_diff/motif_len_; // Effective stutter bp difference (excludes unit changes)
	      out_log_diffs.push_back(log_gt_posterior + *log_phase_ptr + log(abs(eff_diff)));
	      if (bp_diff > 0)
		out_log_up.push_back(log_gt_posterior + *log_phase_ptr);
	      else
		out_log_down.push_back(log_gt_posterior + *log_phase_ptr);
 	    }
	    else {
	      int eff_diff = bp_diff/motif_len_; // Effective stutter repeat difference
	      in_log_diffs.push_back(log_gt_posterior + *log_phase_ptr + log(abs(eff_diff)));
	      if (bp_diff > 0)
		in_log_up.push_back(log_gt_posterior + *log_phase_ptr);
	      else
		in_log_down.push_back(log_gt_posterior + *log_phase_ptr);
	    }
	  }
	}
      }
      log_posterior_ptr += num_samples_;
    }
  }

  // Compute new in-frame parameter estimates
  double in_log_total_up = log_sum_exp(in_log_up), in_log_total_down  = log_sum_exp(in_log_down);
  double in_log_total_eq = log_sum_exp(in_log_eq), in_log_total_diffs = log_sum_exp(in_log_diffs);
  double in_log_total    = log_sum_exp(in_log_total_up, in_log_total_down, in_log_total_eq);
  double in_pgeom_hat    = exp(log_sum_exp(in_log_total_up, in_log_total_down) - in_log_total_diffs);
  double in_pup_hat      = exp(in_log_total_up   - in_log_total);
  double in_pdown_hat    = exp(in_log_total_down - in_log_total);

  // Compute new out-of-frame parameter estimates
  double out_log_total_up = log_sum_exp(out_log_up), out_log_total_down = log_sum_exp(out_log_down), out_log_total_diffs = log_sum_exp(out_log_diffs);
  double out_log_total    = log_sum_exp(out_log_total_up, out_log_total_down);
  double out_pgeom_hat    = exp(out_log_total      - out_log_total_diffs);
  double out_pup_hat      = exp(out_log_total_up   - out_log_total);
  double out_pdown_hat    = exp(out_log_total_down - out_log_total);

  // Update stutter model
  delete stutter_model_;
  stutter_model_ = new StutterModel(in_pgeom_hat, in_pup_hat, in_pdown_hat, out_pgeom_hat, out_pup_hat, out_pdown_hat, motif_len_);
}

/*
  Returns the total log-likelihood given the current stutter model
 */
double EMStutterGenotyper::recalc_log_sample_posteriors(){
  std::vector<double> sample_max_LLs(num_samples_, DBL_MIN);
  double* LL_ptr = log_sample_posteriors_;
  for (int index_1 = 0; index_1 < num_alleles_; ++index_1){
    int len_1 = bps_per_allele_[index_1];
    for (int index_2 = 0; index_2 < num_alleles_; ++index_2){
      int len_2 = bps_per_allele_[index_2];
      std::fill(LL_ptr, LL_ptr+num_samples_, log_gt_priors_[index_1]+log_gt_priors_[index_2]); // Initialize LL's with log genotype priors  
      for (int read_index = 0; read_index < num_reads_; ++read_index)
	LL_ptr[sample_label_[read_index]] += log_sum_exp(log_p1_[read_index]+stutter_model_->calc_log_stutter(len_1, bps_per_allele_[allele_index_[read_index]]), 
							 log_p2_[read_index]+stutter_model_->calc_log_stutter(len_2, bps_per_allele_[allele_index_[read_index]]));
      // Update the per-sample maximum LLs
      for (int sample_index = 0; sample_index < num_samples_; ++sample_index)
	sample_max_LLs[sample_index] = std::max(sample_max_LLs[sample_index], LL_ptr[sample_index]);

      LL_ptr += num_samples_;
    }
  }

  // Compute the normalizing factor for each sample using logsumexp trick
  std::vector<double> sample_total_LLs;
  LL_ptr = log_sample_posteriors_;
  for (int index_1 = 0; index_1 < num_alleles_; ++index_1)
    for (int index_2 = 0; index_2 < num_alleles_; ++index_2)
      for (int sample_index = 0; sample_index < num_samples_; ++sample_index, ++LL_ptr)
	sample_total_LLs[sample_index] += exp(*LL_ptr - sample_max_LLs[sample_index]);
  for (int sample_index = 0; sample_index < num_samples_; ++sample_index)
    sample_total_LLs[sample_index] = sample_max_LLs[sample_index] + log(sample_total_LLs[sample_index]);

  // Compute the total log-likelihood given the current parameters
  double total_LL = sum(sample_total_LLs);

  // Normalize each genotype LL to generate valid log posteriors
  for (int index_1 = 0; index_1 < num_alleles_;++index_1)
    for(int index_2 = 0; index_2 < num_alleles_; ++index_2)
      for (int sample_index = 0; sample_index < num_samples_; ++sample_index, ++LL_ptr)
	*LL_ptr -= sample_total_LLs[sample_index];  

  return total_LL;
}

 
void EMStutterGenotyper::recalc_log_read_phase_posteriors(){
  double* log_phase_ptr = log_read_phase_posteriors_;
  for (int index_1 = 0; index_1 < num_alleles_; ++index_1){
    int len_1 = bps_per_allele_[index_1];
    for (int index_2 = 0; index_2 < num_alleles_; ++index_2){
      int len_2 = bps_per_allele_[index_2];
      for (int read_index = 0; read_index < num_reads_; ++read_index){
	double log_phase_one   = log_p1_[read_index] + stutter_model_->calc_log_stutter(len_1, bps_per_allele_[allele_index_[read_index]]);
	double log_phase_two   = log_p2_[read_index] + stutter_model_->calc_log_stutter(len_2, bps_per_allele_[allele_index_[read_index]]);
	double log_phase_total = log_sum_exp(log_phase_one, log_phase_two);
	log_phase_ptr[0] = log_phase_one-log_phase_total;
	log_phase_ptr[1] = log_phase_two-log_phase_total;
	log_phase_ptr   += 2;
      }
    }
  }
}

bool EMStutterGenotyper::run_EM(int max_iter, double min_LL_abs_change, double min_LL_frac_change){
  // Initialization
  init_log_gt_priors();
  init_stutter_model();

  int num_iter   = 1;
  bool converged = false;
  double LL      = DBL_MIN;
  
  while (num_iter <= max_iter && !converged){
    // E-step
    double new_LL = recalc_log_sample_posteriors();
    recalc_log_read_phase_posteriors();
    std::cerr << "Iteration " << num_iter << ": LL = " << new_LL << "\n" << stutter_model_;

    // M-step                                                                                                                                                               
    recalc_log_gt_priors();
    recalc_stutter_model();
    

    double abs_change  = new_LL - LL;
    double frac_change = -(new_LL - LL)/LL;
    if (abs_change < min_LL_abs_change && frac_change < min_LL_frac_change){
      converged = true;
      return true;
    }

    LL = new_LL;
    num_iter++;
  }
  return false;
}
