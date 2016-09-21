#include <algorithm>
#include <cfloat>
#include <cstring>
#include <sstream>

#include "em_stutter_genotyper.h"
#include "error.h"
#include "mathops.h"

void EMStutterGenotyper::init_log_gt_priors(){
  std::fill(log_gt_priors_, log_gt_priors_+num_alleles_, 1); // Use 1 sample pseudocount                                                                                  
  for (int i = 0; i < num_reads_; i++)
    log_gt_priors_[allele_index_[i]] += 1.0/reads_per_sample_[sample_label_[i]];
  double log_total = log(sum(log_gt_priors_, log_gt_priors_+num_alleles_));
  for (int i = 0; i < num_alleles_; i++){
    log_gt_priors_[i] = log(log_gt_priors_[i]) - log_total;
    assert(log_gt_priors_[i] <= TOLERANCE);
  }
}

void EMStutterGenotyper::recalc_log_gt_priors(){
  // Compute log diploid counts
  double* LL_ptr = log_sample_posteriors_;
  std::vector<double> log_dip_counts;
  for (int index_1 = 0; index_1 < num_alleles_; ++index_1){
    for (int index_2 = 0; index_2 < num_alleles_; ++index_2){
      double log_count = log_sum_exp(LL_ptr, LL_ptr+num_samples_);
      log_dip_counts.push_back(log_count);      
      LL_ptr += num_samples_;
    }
  }

  // Compute log haploid counts from log diploid counts
  std::fill(log_gt_priors_, log_gt_priors_+num_alleles_, -DBL_MAX);
  auto count_iter = log_dip_counts.begin();
  for (int index_1 = 0; index_1 < num_alleles_; ++index_1){
    for (int index_2 = 0; index_2 < num_alleles_; ++index_2, ++count_iter){
      log_gt_priors_[index_1] = log_sum_exp(log_gt_priors_[index_1], *count_iter);
      log_gt_priors_[index_2] = log_sum_exp(log_gt_priors_[index_2], *count_iter);
    }
  }
	
  // Normalize log counts to log probabilities
  double log_total = log_sum_exp(log_gt_priors_, log_gt_priors_+num_alleles_);
  for (int i = 0; i < num_alleles_; i++){
    log_gt_priors_[i] -= log_total;
    assert(log_gt_priors_[i] <= TOLERANCE);
  }
}

void EMStutterGenotyper::init_stutter_model(){
  delete stutter_model_;
  stutter_model_ = new StutterModel(0.9, 0.1, 0.1, 0.8, 0.01, 0.01, motif_len_);
}
  
void EMStutterGenotyper::recalc_stutter_model(){
  std::vector<double> in_log_up,  in_log_down,  in_log_eq, in_log_diffs; // In-frame values
  std::vector<double> out_log_up, out_log_down, out_log_diffs;           // Out-of-frame values
  
  // Add various pseudocounts such that p_geom < 1 for both in-frame and out-of-frame stutter models
  in_log_up.push_back(0.0);  in_log_down.push_back(0.0);  in_log_diffs.push_back(0.0);  in_log_diffs.push_back(log(1.1));
  out_log_up.push_back(0.0); out_log_down.push_back(0.0); out_log_diffs.push_back(0.0); out_log_diffs.push_back(log(1.1));
  in_log_eq.push_back(0.0);

  double* log_posterior_ptr = log_sample_posteriors_;
  double* log_phase_ptr     = log_read_phase_posteriors_;
  for (int index_1 = 0; index_1 < num_alleles_; ++index_1){
    for (int index_2 = 0; index_2 < num_alleles_; ++index_2){
      for (int read_index = 0; read_index < num_reads_; ++read_index){
	double log_gt_posterior = log_posterior_ptr[sample_label_[read_index]];
	for (int phase = 0; phase < 2; ++phase, ++log_phase_ptr){
	  int gt_index = (phase == 0 ? index_1 : index_2);
	  int bp_diff  = bps_per_allele_[allele_index_[read_index]] - bps_per_allele_[gt_index];
	  
	  if (bp_diff == 0)
	    in_log_eq.push_back(log_gt_posterior + *log_phase_ptr);
	  else {
	    if (bp_diff % motif_len_ != 0){
	      int eff_diff = bp_diff - bp_diff/motif_len_; // Effective stutter bp difference (excludes unit changes)
	      out_log_diffs.push_back(log_gt_posterior + *log_phase_ptr + int_log(abs(eff_diff)));
	      if (bp_diff > 0)
		out_log_up.push_back(log_gt_posterior + *log_phase_ptr);
	      else
		out_log_down.push_back(log_gt_posterior + *log_phase_ptr);
 	    }
	    else {
	      int eff_diff = bp_diff/motif_len_; // Effective stutter repeat difference
	      in_log_diffs.push_back(log_gt_posterior + *log_phase_ptr + int_log(abs(eff_diff)));
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

  // Compute new parameter estimates
  double in_log_total_up     = fast_log_sum_exp(in_log_up);
  double in_log_total_down   = fast_log_sum_exp(in_log_down);
  double in_log_total_eq     = fast_log_sum_exp(in_log_eq);
  double in_log_total_diffs  = fast_log_sum_exp(in_log_diffs);
  double out_log_total_up    = fast_log_sum_exp(out_log_up);
  double out_log_total_down  = fast_log_sum_exp(out_log_down);
  double out_log_total_diffs = fast_log_sum_exp(out_log_diffs);
  double out_log_total       = fast_log_sum_exp(out_log_total_up, out_log_total_down);
  double in_pgeom_hat        = std::min(0.999, exp(log_sum_exp(in_log_total_up, in_log_total_down) - in_log_total_diffs));
  double out_pgeom_hat       = std::min(0.999, exp(out_log_total - out_log_total_diffs));
  double log_total           = log_sum_exp(log_sum_exp(in_log_total_up, in_log_total_down, in_log_total_eq), out_log_total);
  double in_pup_hat          = exp(in_log_total_up    - log_total);
  double in_pdown_hat        = exp(in_log_total_down  - log_total);
  double out_pup_hat         = exp(out_log_total_up   - log_total);
  double out_pdown_hat       = exp(out_log_total_down - log_total);

  // Update stutter model
  delete stutter_model_;
  stutter_model_ = new StutterModel(in_pgeom_hat, in_pup_hat, in_pdown_hat, out_pgeom_hat, out_pup_hat, out_pdown_hat, motif_len_);
}

void EMStutterGenotyper::init_log_sample_priors(double* log_sample_ptr){
  // Initialize using the standard Genotyper approach
  Genotyper::init_log_sample_priors(log_sample_ptr);

  // If we're using priors based on allele frequencies, we need to reinitialize the values
  if (use_pop_freqs_ && log_allele_priors_ == NULL){
    double* LL_ptr = log_sample_ptr;
    for (int index_1 = 0; index_1 < num_alleles_; ++index_1){
      for (int index_2 = 0; index_2 < num_alleles_; ++index_2){
	if (!haploid_)
	  std::fill(LL_ptr, LL_ptr+num_samples_, log_gt_priors_[index_1]+log_gt_priors_[index_2]); // Initialize LL's with log genotype priors
	else
	  std::fill(LL_ptr, LL_ptr+num_samples_, (index_1 == index_2 ? log_gt_priors_[index_1] : -DBL_MAX/2)); // Homoz prior is the allele frequency while hetz genotypes are disallowed
	LL_ptr += num_samples_;
      }
    }
  }
}

void EMStutterGenotyper::calc_hap_aln_probs(double* log_aln_probs){
  for (int read_index = 0; read_index < num_reads_; ++read_index)
    for (int allele_id = 0; allele_id < num_alleles_; ++allele_id, ++log_aln_probs)
      *log_aln_probs = stutter_model_->log_stutter_pmf(bps_per_allele_[allele_id], bps_per_allele_[allele_index_[read_index]]);
}

void EMStutterGenotyper::recalc_log_read_phase_posteriors(){
  double* log_phase_ptr = log_read_phase_posteriors_;
  for (int index_1 = 0; index_1 < num_alleles_; ++index_1){
    int len_1 = bps_per_allele_[index_1];
    for (int index_2 = 0; index_2 < num_alleles_; ++index_2){
      int len_2 = bps_per_allele_[index_2];
      for (int read_index = 0; read_index < num_reads_; ++read_index){
	double log_phase_one   = LOG_ONE_HALF + log_p1_[read_index] + stutter_model_->log_stutter_pmf(len_1, bps_per_allele_[allele_index_[read_index]]);
	double log_phase_two   = LOG_ONE_HALF + log_p2_[read_index] + stutter_model_->log_stutter_pmf(len_2, bps_per_allele_[allele_index_[read_index]]);
	double log_phase_total = log_sum_exp(log_phase_one, log_phase_two);
	log_phase_ptr[0] = log_phase_one-log_phase_total;
	log_phase_ptr[1] = log_phase_two-log_phase_total;
	log_phase_ptr   += 2;
      }
    }
  }
}

bool EMStutterGenotyper::train(int max_iter, double min_LL_abs_change, double min_LL_frac_change, bool disp_stats, std::ostream& logger){
  // Initialization
  if (log_allele_priors_ == NULL)
    init_log_gt_priors();
  init_stutter_model();

  int num_iter   = 1;
  bool converged = false;
  double LL      = -DBL_MAX;
  use_pop_freqs_ = true;

  while (num_iter <= max_iter && !converged){
    // E-step
    calc_hap_aln_probs(log_aln_probs_);
    double new_LL = calc_log_sample_posteriors();
    recalc_log_read_phase_posteriors();
    if (disp_stats){
      logger << "Iteration " << num_iter << ": LL = " << new_LL << "\n" << *stutter_model_;
      logger << "Pop freqs: ";
      for (unsigned int i = 0; i < num_alleles_; i++)
	logger << exp(log_gt_priors_[i]) << " ";
      logger << std::endl;
    }

    assert(new_LL <= TOLERANCE);
    if (new_LL < LL+TOLERANCE){
      // Occasionally the LL isn't monotonically increasing b/c of the pseudocounts in
      // recalc_stutter_model(). Let's return true anyways
      return true;
    }

    // M-step
    if (log_allele_priors_ == NULL)
      recalc_log_gt_priors();
    recalc_stutter_model();
    
    double abs_change  = new_LL - LL;
    double frac_change = -(new_LL - LL)/LL;
    if (disp_stats)
      logger << abs_change << " " << min_LL_abs_change << " " << frac_change << " " << min_LL_frac_change << std::endl;
    if (abs_change < min_LL_abs_change && frac_change < min_LL_frac_change){
      converged = true;
      return true;
    }

    LL = new_LL;
    num_iter++;
  }
  return false;
}

bool EMStutterGenotyper::genotype(std::string& chrom_seq, std::ostream& logger){
  use_pop_freqs_ = false;
  if (stutter_model_ == NULL)
    printErrorAndDie("Must specify stutter model before running genotype()");
  calc_hap_aln_probs(log_aln_probs_);
  calc_log_sample_posteriors();
  recalc_log_read_phase_posteriors();
  return true;
}
