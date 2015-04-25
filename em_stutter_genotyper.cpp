#include <algorithm>
#include <cfloat>
#include <sstream>

#include "em_stutter_genotyper.h"
#include "error.h"
#include "mathops.h"

void EMStutterGenotyper::write_vcf_header(std::vector<std::string>& sample_names, std::ostream& out){
  out << "##fileformat=VCFv4.1" << "\n";

  // Info field descriptors
  out << "##INFO=<ID=" << "INFRAME_PGEOM,"  << "Number=1,Type=Float,Description=\""   << "Parameter for in-frame geometric step size distribution"                   << "\">\n"
      << "##INFO=<ID=" << "INFRAME_UP,"     << "Number=1,Type=Float,Description=\""   << "Probability that stutter causes an in-frame increase in obs. STR size"     << "\">\n"
      << "##INFO=<ID=" << "INFRAME_DOWN,"   << "Number=1,Type=Float,Description=\""   << "Probability that stutter causes an in-frame decrease in obs. STR size"     << "\">\n"
      << "##INFO=<ID=" << "OUTFRAME_PGEOM," << "Number=1,Type=Float,Description=\""   << "Parameter for out-of-frame geometric step size distribution"               << "\">\n"
      << "##INFO=<ID=" << "OUTFRAME_UP,"    << "Number=1,Type=Float,Description=\""   << "Probability that stutter causes an out-of-frame increase in obs. STR size" << "\">\n"
      << "##INFO=<ID=" << "OUTFRAME_DOWN,"  << "Number=1,Type=Float,Description=\""   << "Probability that stutter causes an out-of-frame decrease in obs. STR size" << "\">\n"
      << "##INFO=<ID=" << "BPDIFFS,"        << "Number=A,Type=Integer,Description=\"" << "Base pair difference of each alternate allele from the reference allele"   << "\">\n"
      << "##INFO=<ID=" << "END,"            << "Number=1,Type=Integer,Description=\"" << "Inclusive end coordinate for STR's reference allele"                       << "\">\n";

  // Format field descriptors
  out << "##FORMAT=<ID=" << "GT"          << ",Number=1,Type=String,Description=\""   << "Genotype" << "\">" << "\n"
      << "##FORMAT=<ID=" << "GB"          << ",Number=1,Type=String,Description=\""   << "Base pair differences of genotype from reference" << "\">" << "\n"
      << "##FORMAT=<ID=" << "POSTERIOR"   << ",Number=1,Type=Float,Description=\""    << "Posterior probability of phased genotype"                      << "\">" << "\n"
      << "##FORMAT=<ID=" << "DP"          << ",Number=1,Type=Integer,Description=\""  << "Total observed reads for sample"                               << "\">" << "\n"
      << "##FORMAT=<ID=" << "DSNP"        << ",Number=1,Type=Integer,Description=\""  << "Total observed reads for sample with SNP phasing information"  << "\">" << "\n"
      << "##FORMAT=<ID=" << "PDP"         << ",Number=1,Type=String,Description=\""   << "Fractional reads supporting each haploid genotype"             << "\">" << "\n"
      << "##FORMAT=<ID=" << "ALLREADS"    << ",Number=.,Type=Integer,Description=\""  << "Base pair difference observed in each read"                    << "\">" << "\n";

  // Sample names
  out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  for (unsigned int i = 0; i < sample_names.size(); i++)
    out << "\t" << sample_names[i];
  out << "\n";
}

void EMStutterGenotyper::set_allele_priors(vcf::VariantCallFile& variant_file){
  delete [] log_allele_priors_;
  log_allele_priors_ = new double[num_alleles_*num_alleles_*num_samples_];
  
  std::string GP_KEY = "GP";
  

  if (!variant_file.setRegion(chrom_, start_, end_)){
    // Retry setting region if chr is in chromosome name
    if (chrom_.size() <= 3 || chrom_.substr(0, 3).compare("chr") != 0 || !variant_file.setRegion(chrom_.substr(3), start_, end_))
      printErrorAndDie("Failed to set VCF region when obtaining allele priors");
  }
  vcf::Variant variant(variant_file);
  if (!variant_file.getNextVariant(variant))
    printErrorAndDie("Failed to extract VCF entry when obtaining allele priors");

  int sample_count = 0;
  for (auto sample_iter = variant.sampleNames.begin(); sample_iter != variant.sampleNames.end(); ++sample_iter){
    if (sample_indices_.find(*sample_iter) == sample_indices_.end())
      continue;
    int sample_index = sample_indices_.find(*sample_iter)->second;
    sample_count++;

    int gp_index = 0;
    for (unsigned int i = 0; i < num_alleles_; i++){
      for (unsigned int j = 0; j <= i; ++j, ++gp_index){
	double prob = variant.getSampleValueFloat(GP_KEY, *sample_iter, gp_index);
      }
    }
    
    // TO DO: Parse BEAGLE format fields to set priors
    printErrorAndDie("set_allele_priors() function not fully implemented");
  }

  // Ensure that the VCF contained priors for all samples
  if (sample_count != num_samples_)
    printErrorAndDie("BEAGLE VCF only contained allele priors for a subset of samples");
}

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
  in_log_up.push_back(0.0);  in_log_down.push_back(0.0);  in_log_diffs.push_back(0.0);  in_log_diffs.push_back(log(2));
  out_log_up.push_back(0.0); out_log_down.push_back(0.0); out_log_diffs.push_back(0.0); out_log_diffs.push_back(log(2));
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

  // Compute new parameter estimates
  double in_log_total_up     = log_sum_exp(in_log_up);
  double in_log_total_down   = log_sum_exp(in_log_down);
  double in_log_total_eq     = log_sum_exp(in_log_eq);
  double in_log_total_diffs  = log_sum_exp(in_log_diffs);
  double out_log_total_up    = log_sum_exp(out_log_up);
  double out_log_total_down  = log_sum_exp(out_log_down);
  double out_log_total_diffs = log_sum_exp(out_log_diffs);
  double out_log_total       = log_sum_exp(out_log_total_up, out_log_total_down);
  double in_pgeom_hat        = exp(log_sum_exp(in_log_total_up, in_log_total_down) - in_log_total_diffs);
  double out_pgeom_hat       = exp(out_log_total - out_log_total_diffs);
  double log_total           = log_sum_exp(log_sum_exp(in_log_total_up, in_log_total_down, in_log_total_eq), out_log_total);
  double in_pup_hat          = exp(in_log_total_up    - log_total);
  double in_pdown_hat        = exp(in_log_total_down  - log_total);
  double out_pup_hat         = exp(out_log_total_up   - log_total);
  double out_pdown_hat       = exp(out_log_total_down - log_total);

  // Update stutter model
  delete stutter_model_;
  stutter_model_ = new StutterModel(in_pgeom_hat, in_pup_hat, in_pdown_hat, out_pgeom_hat, out_pup_hat, out_pdown_hat, motif_len_);
}

/*
  Returns the total log-likelihood given the current stutter model
 */
double EMStutterGenotyper::recalc_log_sample_posteriors(bool use_pop_freqs){
  std::vector<double> sample_max_LLs(num_samples_, -DBL_MAX);
  double* LL_ptr = log_sample_posteriors_;

  if (use_pop_freqs){
    // If per-allele priors have been set for each sample, use them
    if (log_allele_priors_ != NULL)
      memcpy(log_sample_posteriors_, log_allele_priors_, num_alleles_*num_alleles_*num_samples_*sizeof(double));

    // Otherwise we'll set them in the for loop below on a per-allele basis
  }
  else
    std::fill(log_sample_posteriors_, log_sample_posteriors_+(num_alleles_*num_alleles_*num_samples_), -2*log(num_alleles_));

  for (int index_1 = 0; index_1 < num_alleles_; ++index_1){
    int len_1 = bps_per_allele_[index_1];
    for (int index_2 = 0; index_2 < num_alleles_; ++index_2){
      int len_2 = bps_per_allele_[index_2];
      if (use_pop_freqs && log_allele_priors_ == NULL)
	  std::fill(LL_ptr, LL_ptr+num_samples_, log_gt_priors_[index_1]+log_gt_priors_[index_2]); // Initialize LL's with log genotype priors  

      for (int read_index = 0; read_index < num_reads_; ++read_index){
	LL_ptr[sample_label_[read_index]] += log_sum_exp(LOG_ONE_HALF + log_p1_[read_index]+stutter_model_->log_stutter_pmf(len_1, bps_per_allele_[allele_index_[read_index]]), 
							 LOG_ONE_HALF + log_p2_[read_index]+stutter_model_->log_stutter_pmf(len_2, bps_per_allele_[allele_index_[read_index]]));
	assert(LL_ptr[sample_label_[read_index]] <= TOLERANCE);
      }
      // Update the per-sample maximum LLs
      for (int sample_index = 0; sample_index < num_samples_; ++sample_index)
	sample_max_LLs[sample_index] = std::max(sample_max_LLs[sample_index], LL_ptr[sample_index]);

      LL_ptr += num_samples_;
    }
  }

  // Compute the normalizing factor for each sample using logsumexp trick
  std::vector<double> sample_total_LLs(num_samples_, 0.0);
  LL_ptr = log_sample_posteriors_;
  for (int index_1 = 0; index_1 < num_alleles_; ++index_1)
    for (int index_2 = 0; index_2 < num_alleles_; ++index_2)
      for (int sample_index = 0; sample_index < num_samples_; ++sample_index, ++LL_ptr)
	sample_total_LLs[sample_index] += exp(*LL_ptr - sample_max_LLs[sample_index]);
  for (int sample_index = 0; sample_index < num_samples_; ++sample_index){
    sample_total_LLs[sample_index] = sample_max_LLs[sample_index] + log(sample_total_LLs[sample_index]);    
    assert(sample_total_LLs[sample_index] <= TOLERANCE);
  }
  // Compute the total log-likelihood given the current parameters
  double total_LL = sum(sample_total_LLs);

  // Normalize each genotype LL to generate valid log posteriors
  LL_ptr = log_sample_posteriors_;
  for (int index_1 = 0; index_1 < num_alleles_; ++index_1)
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

bool EMStutterGenotyper::train(int max_iter, double min_LL_abs_change, double min_LL_frac_change){
  // Initialization
  init_log_gt_priors();
  init_stutter_model();

  int num_iter   = 1;
  bool converged = false;
  double LL      = -DBL_MAX;

  while (num_iter <= max_iter && !converged){
    // E-step
    double new_LL = recalc_log_sample_posteriors(true);
    recalc_log_read_phase_posteriors();
    std::cerr << "Iteration " << num_iter << ": LL = " << new_LL << "\n" << *stutter_model_;
    std::cerr << "Pop freqs: ";
    for (unsigned int i = 0; i < num_alleles_; i++)
      std::cerr << exp(log_gt_priors_[i]) << " ";
    std::cerr << std::endl;
    
    assert(new_LL <= TOLERANCE);
    if (new_LL < LL+TOLERANCE)
      return false;

    // M-step
    if (log_allele_priors_ == NULL) 
      recalc_log_gt_priors();
    recalc_stutter_model();
    
    double abs_change  = new_LL - LL;
    double frac_change = -(new_LL - LL)/LL;
    std::cerr << abs_change << " " << min_LL_abs_change << " " << frac_change << " " << min_LL_frac_change << std::endl;
    if (abs_change < min_LL_abs_change && frac_change < min_LL_frac_change){
      converged = true;
      return true;
    }

    LL = new_LL;
    num_iter++;
  }
  return false;
}

bool EMStutterGenotyper::genotype(bool use_pop_freqs){
  if (stutter_model_ == NULL)
    printErrorAndDie("Must specify stutter model before running genotype()");
  recalc_log_sample_posteriors(use_pop_freqs);
  recalc_log_read_phase_posteriors();
  return true;
}

std::string EMStutterGenotyper::get_allele(std::string& ref_allele, int bp_diff){
  if (bp_diff < -((int)ref_allele.size())){
    std::cerr << "Invalid bp diff size: " << ref_allele << " " << ref_allele.size() << " " << bp_diff << std::endl;
    assert(bp_diff >= -((int)ref_allele.size()));
  }
    
  if (bp_diff == 0)
    return ref_allele;
  else if (bp_diff < 0) {
    if (bp_diff + ref_allele.size() == 0)
      return "*";
    else
      return ref_allele.substr(0, ref_allele.size()+bp_diff);
  }
  else {
    std::stringstream ss;
    std::string motif = ref_allele.substr(ref_allele.size()-motif_len_);
    ss << ref_allele;
    while (bp_diff >= motif_len_){
      ss << motif;
      bp_diff -= motif_len_;
    }
    if (bp_diff > 0)
      ss << motif.substr(0, bp_diff);
    return ss.str();
  }
}

void EMStutterGenotyper::write_vcf_record(std::string& ref_allele, std::vector<std::string>& sample_names, std::ostream& out){
  std::vector< std::pair<int,int> > gts(num_samples_, std::pair<int,int>(-1,-1));
  std::vector<double> log_phased_posteriors(num_samples_, -DBL_MAX);
  std::vector< std::vector<double> > log_read_phases(num_samples_);
  std::vector<double> log_unphased_posteriors, phase_probs;

  // TO DO: Consider selecting GT based on genotype with maximum UNPHASED posterior instead of maximum PHASED posterior
  // Are we then double-counting het GTs vs hom GTs?

  // Extract each sample's MAP genotype and the associated posterior
  double* log_post_ptr = log_sample_posteriors_;
  for (int index_1 = 0; index_1 < num_alleles_; ++index_1)
    for (int index_2 = 0; index_2 < num_alleles_; ++index_2)
      for (unsigned int sample_index = 0; sample_index < num_samples_; ++sample_index, ++log_post_ptr){
	if (*log_post_ptr > log_phased_posteriors[sample_index]){
	  log_phased_posteriors[sample_index] = *log_post_ptr;
	  gts[sample_index] = std::pair<int,int>(index_1, index_2);
	}
      }

  // Extract the phasing probability conditioned on the determined sample genotypes
  for (unsigned int sample_index = 0; sample_index < num_samples_; sample_index++){
    int gt_a = gts[sample_index].first, gt_b = gts[sample_index].second;
    if (gt_a == gt_b){
      log_unphased_posteriors.push_back(log_phased_posteriors[sample_index]);
      phase_probs.push_back(1.0);
    }
    else {
      double log_p1  = log_phased_posteriors[sample_index];
      double log_p2  = log_sample_posteriors_[gt_b*num_alleles_*num_samples_ + gt_a*num_samples_ + sample_index];
      double log_tot = log_sum_exp(log_p1, log_p2);
      log_unphased_posteriors.push_back(log_tot);
      phase_probs.push_back(exp(log_p1-log_tot));
    }
  }

  // Determine the number of reads with SNP information for each sample
  std::vector<int> num_reads_with_snps(num_samples_, 0);
  for (unsigned int read_index = 0; read_index < num_reads_; read_index++)
    num_reads_with_snps[sample_label_[read_index]] += (log_p1_[read_index] != log_p2_[read_index]);
  
  // Extract each read's phase posterior conditioned on the determined sample genotypes
  for (unsigned int read_index = 0; read_index < num_reads_; read_index++){
    int gt_a = gts[sample_label_[read_index]].first;
    int gt_b = gts[sample_label_[read_index]].second;
    log_read_phases[sample_label_[read_index]].push_back(log_read_phase_posteriors_[2*num_reads_*(gt_a*num_alleles_ + gt_b) + 2*read_index]);
  }

  // Determine the bp differences observed in reads for each sample
  std::vector< std::vector<int> > bps_per_sample(num_samples_);
  for (unsigned int read_index = 0; read_index < num_reads_; read_index++){
    bps_per_sample[sample_label_[read_index]].push_back(bps_per_allele_[allele_index_[read_index]]);
  }

  //VCF line format = CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE_1 SAMPLE_2 ... SAMPLE_N
  out << chrom_ << "\t" << start_ << "\t" << ".";

  // Add reference allele and alternate alleles
  out << "\t" << get_allele(ref_allele, bps_per_allele_[0]) << "\t";
  if (num_alleles_ == 1)
    out << ".";
  else {
    for (int i = 1; i < num_alleles_-1; i++)
      out << get_allele(ref_allele, bps_per_allele_[i]) << ",";
    out << get_allele(ref_allele, bps_per_allele_[num_alleles_-1]);
  }

  // Add QUAL and FILTER fields
  out << "\t" << "." << "\t" << ".";

  // Add INFO field items
  out << "\tINFRAME_PGEOM=" << stutter_model_->get_parameter(true,  'P') << ";" 
      << "INFRAME_UP="      << stutter_model_->get_parameter(true,  'U') << ";" 
      << "INFRAME_DOWN="    << stutter_model_->get_parameter(true,  'D') << ";" 
      << "OUTFRAME_PGEOM="  << stutter_model_->get_parameter(false, 'P') << ";" 
      << "OUTFRAME_UP="     << stutter_model_->get_parameter(false, 'U') << ";" 
      << "OUTFRAME_DOWN="   << stutter_model_->get_parameter(false, 'D') << ";"
      << "END="             << end_ << ";";
  if (num_alleles_ > 1){
   out << "BPDIFFS=" << bps_per_allele_[1];
    for (unsigned int i = 2; i < num_alleles_; i++)
      out << "," << bps_per_allele_[i];
    out << ";";
  }

  // Add FORMAT field
  out << "\tGT:GB:POSTERIOR:DP:DSNP:PDP:ALLREADS";

  for (unsigned int i = 0; i < sample_names.size(); i++){
    out << "\t";
    auto sample_iter = sample_indices_.find(sample_names[i]);
    if (sample_iter == sample_indices_.end() || reads_per_sample_[sample_iter->second] == 0){
      out << ".";
      continue;
    }

    int sample_index    = sample_iter->second;
    int total_reads     = reads_per_sample_[sample_index];
    double phase1_reads = exp(log_sum_exp(log_read_phases[sample_index]));
    double phase2_reads = total_reads - phase1_reads;
    std::sort(bps_per_sample[sample_index].begin(), bps_per_sample[sample_index].end());

    out << gts[sample_index].first << "|" << gts[sample_index].second     // Genotype
	<< ":" << bps_per_allele_[gts[sample_index].first] << "|" << bps_per_allele_[gts[sample_index].second] // Bp diffs from reference
	<< ":" << exp(log_phased_posteriors[sample_index])                // Posterior
	<< ":" << total_reads                                             // Total reads
	<< ":" << num_reads_with_snps[sample_index]                       // Total reads with SNP information
	<< ":" << phase1_reads << "|" << phase2_reads                     // Reads per allele
	<< ":" << bps_per_sample[sample_index][0];
    for (unsigned int j = 1; j < bps_per_sample[sample_index].size(); j++)
      out << "," << bps_per_sample[sample_index][j]; 
  }

  out << "\n";
}

