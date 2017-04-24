#include <assert.h>
#include <cfloat>
#include <cstring>
#include <time.h>

#include <algorithm>
#include <cfloat>
#include <cmath>

#include "genotyper.h"
#include "mathops.h"

// Each genotype has an equal total prior, but heterozygotes have two possible phasings. Therefore,
// i)   Phased heterozygotes have a prior of 1/(n(n+1))
// ii)  Homozygotes have a prior of 2/(n(n+1))
// iii) Total prior is n*2/(n(n+1)) + n(n-1)*1/(n(n+1)) = 2/(n+1) + (n-1)/(n+1) = 1

double Genotyper::log_homozygous_prior() const {
  if (haploid_)
    return -int_log(num_alleles_);
  else
    return int_log(2) - int_log(num_alleles_) - int_log(num_alleles_+1);
}

double Genotyper::log_heterozygous_prior() const {
  if (haploid_)
    return -DBL_MAX/2;
  else
    return -int_log(num_alleles_) - int_log(num_alleles_+1);
}

void Genotyper::init_log_sample_priors(double* log_sample_ptr){
  const double log_homoz_prior = log_homozygous_prior();
  const double log_hetz_prior  = log_heterozygous_prior();
  double* LL_ptr = log_sample_ptr;
  for (unsigned int i = 0; i < num_samples_; ++i)
    for (unsigned int j = 0; j < num_alleles_; ++j)
      for (unsigned int k = 0; k < num_alleles_; ++k, ++LL_ptr)
	  *LL_ptr = (j == k ? log_homoz_prior : log_hetz_prior);
}

double Genotyper::calc_log_sample_posteriors(std::vector<int>& read_weights){
  double posterior_time = clock();
  assert(read_weights.size() == num_reads_);
  init_log_sample_priors(log_sample_posteriors_);

  const int num_diplotypes = num_alleles_*num_alleles_;
  double* read_LL_ptr      = log_aln_probs_;
  for (int read_index = 0; read_index < num_reads_; ++read_index){
    double* sample_LL_ptr = log_sample_posteriors_ + num_diplotypes*sample_label_[read_index];
    for (int index_1 = 0; index_1 < num_alleles_; ++index_1){
      for (int index_2 = 0; index_2 < num_alleles_; ++index_2, ++sample_LL_ptr){
        *sample_LL_ptr += read_weights[read_index]*fast_log_sum_exp(LOG_ONE_HALF + log_p1_[read_index] + read_LL_ptr[index_1],
								    LOG_ONE_HALF + log_p2_[read_index] + read_LL_ptr[index_2]);
        assert(*sample_LL_ptr <= TOLERANCE);
      }
    }
    read_LL_ptr += num_alleles_;
  }

  // Compute each sample's total LL and normalize each genotype LL to generate valid log posteriors
  double* sample_LL_ptr = log_sample_posteriors_;
  for (int sample_index = 0; sample_index < num_samples_; ++sample_index){
    const double sample_total_LL = log_sum_exp(sample_LL_ptr, sample_LL_ptr+num_diplotypes);
    sample_total_LLs_[sample_index] = sample_total_LL;
    assert(sample_total_LL <= TOLERANCE);
    for (int index_1 = 0; index_1 < num_alleles_; ++index_1)
      for (int index_2 = 0; index_2 < num_alleles_; ++index_2, ++sample_LL_ptr)
	*sample_LL_ptr -= sample_total_LL;
  }

  // Compute the total log-likelihood given the current parameters
  double total_LL = sum(sample_total_LLs_, sample_total_LLs_ + num_samples_);

  posterior_time         = (clock() - posterior_time)/CLOCKS_PER_SEC;
  total_posterior_time_ += posterior_time;
  return total_LL;
}

void Genotyper::get_optimal_haplotypes(std::vector< std::pair<int, int> >& gts) const {
  assert(gts.size() == 0);
  gts = std::vector< std::pair<int,int> > (num_samples_, std::pair<int,int>(-1,-1));
  double* log_posterior_ptr = log_sample_posteriors_;
  std::vector<double> log_phased_posteriors(num_samples_, -DBL_MAX);
  for (unsigned int sample_index = 0; sample_index < num_samples_; ++sample_index){
    for (int index_1 = 0; index_1 < num_alleles_; ++index_1){
      for (int index_2 = 0; index_2 < num_alleles_; ++index_2, ++log_posterior_ptr){
        if (*log_posterior_ptr > log_phased_posteriors[sample_index]){
          log_phased_posteriors[sample_index] = *log_posterior_ptr;
          gts[sample_index] = std::pair<int,int>(index_1, index_2);
        }
      }
    }
  }
}

void Genotyper::calc_PLs(const std::vector<double>& gls, std::vector<int>& pls) const {
  assert(pls.empty());
  double max_gl = *(std::max_element(gls.begin(), gls.end()));
  for (unsigned int i = 0; i < gls.size(); i++)
    pls.push_back(std::min(999, (int)(-10*(gls[i]-max_gl))));
}

double Genotyper::calc_gl_diff(const std::vector<double>& gls, int gt_a, int gt_b) const {
  if (num_alleles_ == 1)
    return -1000;

  double max_gl    = *(std::max_element(gls.begin(), gls.end()));
  double second_gl = -DBL_MAX;
  for (unsigned int i = 0; i < gls.size(); i++)
    if (gls[i] < max_gl)
      second_gl = std::max(second_gl, gls[i]);
  if (second_gl == -DBL_MAX)
    second_gl = max_gl;

  int gl_index;
  if (haploid_)
    gl_index = gt_a;
  else {
    int min_gt = std::min(gt_a, gt_b);
    int max_gt = std::max(gt_a, gt_b);
    gl_index   = max_gt*(max_gt+1)/2 + min_gt;
  }
  return ((std::abs(max_gl-gls[gl_index]) < TOLERANCE) ? (max_gl-second_gl) : gls[gl_index]-max_gl);
}

void Genotyper::extract_genotypes_and_likelihoods(int num_variants, std::vector<int>& hap_to_allele,
						  std::vector< std::pair<int,int>  >& best_haplotypes,
						  std::vector< std::pair<int,int>  >& best_gts,
						  std::vector<double>& log_phased_posteriors, std::vector<double>& log_unphased_posteriors,
						  bool calc_gls,        std::vector< std::vector<double> >& gls, std::vector<double>& gl_diffs,
						  bool calc_pls,        std::vector< std::vector<int> >& pls,
						  bool calc_phased_gls, std::vector< std::vector<double> >& phased_gls){
  assert(log_phased_posteriors.empty() && log_unphased_posteriors.empty() && gl_diffs.empty());
  assert(best_haplotypes.empty() && best_gts.empty() && gls.empty() && pls.empty() && phased_gls.empty());

  // Use the standard genotyper approach to find the ML combination of haplotypes
  get_optimal_haplotypes(best_haplotypes);

  // Extract the ML alleles for the variant
  for (int sample_index = 0; sample_index < num_samples_; sample_index++)
    best_gts.push_back(std::pair<int,int>(hap_to_allele[best_haplotypes[sample_index].first], hap_to_allele[best_haplotypes[sample_index].second]));

  // Marginalize over all haplotypes to compute the genotype posteriors and use streaming log-sum-exp to aggregate values
  std::vector< std::vector<double>  > max_log_phased_posteriors   (num_samples_, std::vector<double>(num_variants*num_variants, -DBL_MAX/2));
  std::vector< std::vector<double>  > total_log_phased_posteriors (num_samples_, std::vector<double>(num_variants*num_variants, 0.0));
  double* log_posterior_ptr = log_sample_posteriors_;
  for (unsigned int sample_index = 0; sample_index < num_samples_; ++sample_index)
    for (int index_1 = 0; index_1 < num_alleles_; ++index_1){
      for (int index_2 = 0; index_2 < num_alleles_; ++index_2, ++log_posterior_ptr){
	int gt_index = num_variants*hap_to_allele[index_1] + hap_to_allele[index_2];
	update_streaming_log_sum_exp(*log_posterior_ptr, max_log_phased_posteriors[sample_index][gt_index], total_log_phased_posteriors[sample_index][gt_index]);
    }
  }
  int gt_index = 0;
  for (int index_1 = 0; index_1 < num_variants; ++index_1){
    for (int index_2 = 0; index_2 < num_variants; ++index_2, ++gt_index){
      for (unsigned int sample_index = 0; sample_index < num_samples_; ++sample_index){
	double total = finish_streaming_log_sum_exp(max_log_phased_posteriors[sample_index][gt_index], total_log_phased_posteriors[sample_index][gt_index]);
	total_log_phased_posteriors[sample_index][gt_index] = total;
      }
    }
  }

  // Store the aggregated posterior values in the provided vectors
  for (int sample_index = 0; sample_index < num_samples_; sample_index++){
    int gt_a = best_gts[sample_index].first, gt_b = best_gts[sample_index].second;
    double log_phased_prob = total_log_phased_posteriors[sample_index][num_variants*gt_a + gt_b];
    log_phased_posteriors.push_back(log_phased_prob);
    if (gt_a == gt_b)
      log_unphased_posteriors.push_back(log_phased_prob);
    else {
      double alt_log_phased_prob =  total_log_phased_posteriors[sample_index][num_variants*gt_b + gt_a];
      log_unphased_posteriors.push_back(log_sum_exp(log_phased_prob, alt_log_phased_prob));
    }
  }
  assert(best_haplotypes.size() == num_samples_ && best_gts.size() == num_samples_);
  assert(log_phased_posteriors.size() == num_samples_ && log_unphased_posteriors.size() == num_samples_);

  // Compute GLs and phased GLs if necessary
  if (calc_gls || calc_phased_gls || calc_pls){
    // The genotype likelihoods should not contain the priors we used during the posterior calculation
    // To obtain the true likelihoods, we subtract out the priors from the posteriors using these values.
    gls = std::vector< std::vector<double> >(num_samples_);
    if (calc_phased_gls)
      phased_gls = std::vector< std::vector<double> >(num_samples_);
    double hom_ll_correction  = log_homozygous_prior();
    double het_ll_correction  = (haploid_ ? 0 : log_heterozygous_prior());                    // If haploid, don't correct hetz genotypes as they're impossible
    double nconfig_correction = int_log(2) + 2*(int_log(num_alleles_)-int_log(num_variants)); // Need to correct for the number of haplotypes whose GL's we're averaging
    int gt_index = 0;
    for (int index_1 = 0; index_1 < num_variants; ++index_1){
      for (int index_2 = 0; index_2 < num_variants; ++index_2, ++gt_index){
	int alt_gt_index              = index_2*num_variants + index_1;
	double gl_ll_correction       = (index_1 == index_2 ? hom_ll_correction : het_ll_correction) + nconfig_correction;
	double phasedgl_ll_correction = (index_1 == index_2 ? hom_ll_correction : het_ll_correction);
	for (int sample_index = 0; sample_index < num_samples_; sample_index++){
	  if (index_2 <= index_1){
	    if (!haploid_ || (index_1 == index_2)){
	      double gl_base_e = sample_total_LLs_[sample_index] - gl_ll_correction + fast_log_sum_exp(total_log_phased_posteriors[sample_index][gt_index],
												       total_log_phased_posteriors[sample_index][alt_gt_index]);
	      gls[sample_index].push_back(gl_base_e*LOG_E_BASE_10); // Convert from ln to log10
	    }
	  }
	  if (calc_phased_gls)
	    if (!haploid_ || (index_1 == index_2))
	      phased_gls[sample_index].push_back((total_log_phased_posteriors[sample_index][gt_index] + sample_total_LLs_[sample_index] - phasedgl_ll_correction)*LOG_E_BASE_10);
	}
      }
    }

    // GLDIFFs
    for (int sample_index = 0; sample_index < num_samples_; sample_index++)
      gl_diffs.push_back(calc_gl_diff(gls[sample_index], best_gts[sample_index].first, best_gts[sample_index].second));

    // PLs
    if (calc_pls){
      pls = std::vector< std::vector<int> >(num_samples_);
      for (int sample_index = 0; sample_index < num_samples_; sample_index++)
	calc_PLs(gls[sample_index], pls[sample_index]);
    }

    if (!calc_gls)
      gls.clear();
  }
}

void Genotyper::write_vcf_header(const std::string& full_command, const std::vector<std::string>& sample_names,
				 bool output_gls, bool output_pls, bool output_phased_gls, std::ostream& out){
  out << "##fileformat=VCFv4.1" << "\n"
      << "##command=" << full_command << "\n";

  // Info field descriptors
  out << "##INFO=<ID=" << "INFRAME_PGEOM"  << ",Number=1,Type=Float,Description=\""   << "Parameter for in-frame geometric step size distribution"                      << "\">\n"
      << "##INFO=<ID=" << "INFRAME_UP"     << ",Number=1,Type=Float,Description=\""   << "Probability that stutter causes an in-frame increase in obs. STR size"        << "\">\n"
      << "##INFO=<ID=" << "INFRAME_DOWN"   << ",Number=1,Type=Float,Description=\""   << "Probability that stutter causes an in-frame decrease in obs. STR size"        << "\">\n"
      << "##INFO=<ID=" << "OUTFRAME_PGEOM" << ",Number=1,Type=Float,Description=\""   << "Parameter for out-of-frame geometric step size distribution"                  << "\">\n"
      << "##INFO=<ID=" << "OUTFRAME_UP"    << ",Number=1,Type=Float,Description=\""   << "Probability that stutter causes an out-of-frame increase in read's STR size"  << "\">\n"
      << "##INFO=<ID=" << "OUTFRAME_DOWN"  << ",Number=1,Type=Float,Description=\""   << "Probability that stutter causes an out-of-frame decrease in read's STR size"  << "\">\n"
      << "##INFO=<ID=" << "BPDIFFS"        << ",Number=A,Type=Integer,Description=\"" << "Base pair difference of each alternate allele from the reference allele"      << "\">\n"
      << "##INFO=<ID=" << "START"          << ",Number=1,Type=Integer,Description=\"" << "Inclusive start coodinate for the repetitive portion of the reference allele" << "\">\n"
      << "##INFO=<ID=" << "END"            << ",Number=1,Type=Integer,Description=\"" << "Inclusive end coordinate for the repetitive portion of the reference allele"  << "\">\n"
      << "##INFO=<ID=" << "PERIOD"         << ",Number=1,Type=Integer,Description=\"" << "Length of STR motif"                                                          << "\">\n"
      << "##INFO=<ID=" << "AN"             << ",Number=1,Type=Integer,Description=\"" << "Total number of alleles in called genotypes"                                  << "\">\n"
      << "##INFO=<ID=" << "REFAC"          << ",Number=1,Type=Integer,Description=\"" << "Reference allele count"                                                       << "\">\n"
      << "##INFO=<ID=" << "AC"             << ",Number=A,Type=Integer,Description=\"" << "Alternate allele counts"                                                      << "\">\n"
      << "##INFO=<ID=" << "NSKIP"          << ",Number=1,Type=Integer,Description=\"" << "Number of samples not genotyped due to various issues"                        << "\">\n"
      << "##INFO=<ID=" << "NFILT"          << ",Number=1,Type=Integer,Description=\"" << "Number of samples whose genotypes were filtered due to various issues"        << "\">\n"
      << "##INFO=<ID=" << "DP"             << ",Number=1,Type=Integer,Description=\"" << "Total number of valid reads used to genotype all samples"                     << "\">\n"
      << "##INFO=<ID=" << "DSNP"           << ",Number=1,Type=Integer,Description=\"" << "Total number of reads with SNP phasing information"                           << "\">\n"
      << "##INFO=<ID=" << "DSTUTTER"       << ",Number=1,Type=Integer,Description=\"" << "Total number of reads with a stutter indel in the STR region"                 << "\">\n"
      << "##INFO=<ID=" << "DFLANKINDEL"    << ",Number=1,Type=Integer,Description=\"" << "Total number of reads with an indel in the regions flanking the STR"          << "\">\n";

  // Format field descriptors
  out << "##FORMAT=<ID=" << "GT"          << ",Number=1,Type=String,Description=\""  << "Genotype" << "\">" << "\n"
      << "##FORMAT=<ID=" << "GB"          << ",Number=1,Type=String,Description=\""  << "Base pair differences of genotype from reference"              << "\">" << "\n"
      << "##FORMAT=<ID=" << "Q"           << ",Number=1,Type=Float,Description=\""   << "Posterior probability of unphased genotype"                    << "\">" << "\n"
      << "##FORMAT=<ID=" << "PQ"          << ",Number=1,Type=Float,Description=\""   << "Posterior probability of phased genotype"                      << "\">" << "\n"
      << "##FORMAT=<ID=" << "DP"          << ",Number=1,Type=Integer,Description=\"" << "Number of valid reads used for sample's genotype"              << "\">" << "\n"
      << "##FORMAT=<ID=" << "DSNP"        << ",Number=1,Type=Integer,Description=\"" << "Number of reads with SNP phasing information"                  << "\">" << "\n"
      << "##FORMAT=<ID=" << "PSNP"        << ",Number=1,Type=String,Description=\""  << "Number of reads with SNPs supporting each haploid genotype"    << "\">" << "\n"
      << "##FORMAT=<ID=" << "PDP"         << ",Number=1,Type=String,Description=\""  << "Fractional reads supporting each haploid genotype"             << "\">" << "\n"
      << "##FORMAT=<ID=" << "GLDIFF"      << ",Number=1,Type=Float,Description=\""   << "Difference in likelihood between the reported and next best genotypes"  << "\">" << "\n"
      << "##FORMAT=<ID=" << "DSTUTTER"    << ",Number=1,Type=Integer,Description=\"" << "Number of reads with a stutter indel in the STR region"        << "\">" << "\n"
      << "##FORMAT=<ID=" << "DFLANKINDEL" << ",Number=1,Type=Integer,Description=\"" << "Number of reads with an indel in the regions flanking the STR" << "\">" << "\n"
      << "##FORMAT=<ID=" << "AB"          << ",Number=1,Type=Float,Description=\""
      << "log10 of the allele bias pvalue, where 0 is no bias and more negative values are increasingly biased. Not applicable to homozygous genotypes" << "\">" << "\n"
      << "##FORMAT=<ID=" << "DAB"         << ",Number=1,Type=Integer,Description=\"" << "Number of reads used in the allele bias calculation"           << "\">" << "\n"
      << "##FORMAT=<ID=" << "ALLREADS"    << ",Number=1,Type=String,Description=\""  << "Base pair difference observed in each read's Needleman-Wunsch alignment" << "\">" << "\n"
      << "##FORMAT=<ID=" << "MALLREADS"   << ",Number=1,Type=String,Description=\""
      << "Maximum likelihood bp diff in each read based on haplotype alignments for reads that span the repeat region by at least 5 base pairs"            << "\">" << "\n";

  if (output_gls)
    out << "##FORMAT=<ID=" << "GL"       << ",Number=G,Type=Float,Description=\""   << "log-10 genotype likelihoods" << "\">" << "\n";
  if (output_pls)
    out << "##FORMAT=<ID=" << "PL"       << ",Number=G,Type=Integer,Description=\"" << "Phred-scaled genotype likelihoods" << "\">" << "\n";
  if (output_phased_gls)
    out << "##FORMAT=<ID=" << "PHASEDGL" << ",Number=.,Type=Float,Description=\""
	<< "log-10 genotype likelihood for each phased genotype. Value for phased genotype X|Y is stored at a 0-based index of X*A + Y, where A is the number of alleles. Not applicable to haploid genotypes"
	<< "\">" << "\n";

  // Sample names
  out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  for (unsigned int i = 0; i < sample_names.size(); i++)
    out << "\t" << sample_names[i];
  out << "\n";
}
