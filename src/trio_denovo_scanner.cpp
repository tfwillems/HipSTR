#include <stdlib.h>

#include <cfloat>
#include <vector>

#include "trio_denovo_scanner.h"
#include "denovo_allele_priors.h"
#include "error.h"
#include "mathops.h"
#include "mutation_model.h"
#include "vcf_input.h"

std::string TrioDenovoScanner::BPDIFFS_KEY = "BPDIFFS";
std::string TrioDenovoScanner::START_KEY   = "START";
std::string TrioDenovoScanner::END_KEY     = "END";
std::string TrioDenovoScanner::PERIOD_KEY  = "PERIOD";

void TrioDenovoScanner::write_vcf_header(std::string& full_command){
  denovo_vcf_ << "##fileformat=VCFv4.1" << "\n"
	      << "##command=" << full_command << "\n";

  // Info field descriptors
  denovo_vcf_ << "##INFO=<ID=" << "BPDIFFS" << ",Number=A,Type=Integer,Description=\"" << "Base pair difference of each alternate allele from the reference allele"      << "\">\n"
	      << "##INFO=<ID=" << "START"   << ",Number=1,Type=Integer,Description=\"" << "Inclusive start coodinate for the repetitive portion of the reference allele" << "\">\n"
	      << "##INFO=<ID=" << "END"     << ",Number=1,Type=Integer,Description=\"" << "Inclusive end coordinate for the repetitive portion of the reference allele"  << "\">\n"
	      << "##INFO=<ID=" << "PERIOD"  << ",Number=1,Type=Integer,Description=\"" << "Length of STR motif"                                                          << "\">\n";

  // Format field descriptors
  denovo_vcf_ << "##FORMAT=<ID=" << "NOMUT"  << ",Number=1,Type=Float,Description=\"" << "Log10-likelihood that no mutations occurred in the child"              << "\">\n"
	      << "##FORMAT=<ID=" << "DENOVO" << ",Number=1,Type=Float,Description=\"" << "Log10-likelihood that a single de novo mutation occurred in the child" << "\">\n"
	      << "##FORMAT=<ID=" << "OTHER"  << ",Number=1,Type=Float,Description=\"" << "Log10-likelihood that a single other mutation occurred in the child"   << "\">\n";

  denovo_vcf_ << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  for (auto family_iter = families_.begin(); family_iter != families_.end(); family_iter++)
    for (auto child_iter = family_iter->get_children().begin(); child_iter != family_iter->get_children().end(); ++child_iter)
      denovo_vcf_ << "\t" << *child_iter;
  denovo_vcf_ << "\n";
}

void TrioDenovoScanner::initialize_vcf_record(VCF::Variant& str_variant){
  // VCF line format = CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE_1 SAMPLE_2 ... SAMPLE_N
  denovo_vcf_ << str_variant.get_chromosome() << "\t" << str_variant.get_position() << "\t" << str_variant.get_id() << "\t" << str_variant.get_allele(0) << "\t";
  if (str_variant.num_alleles() > 1){
    denovo_vcf_ << str_variant.get_allele(1);
    for (int i = 2; i < str_variant.num_alleles(); i++)
      denovo_vcf_ << "," << str_variant.get_allele(i);
  }
  else
    denovo_vcf_ << ".";
  denovo_vcf_ << "\t" << "." << "\t" << "." << "\t";

  // INFO field
  int32_t start;  str_variant.get_INFO_value_single_int(START_KEY, start);
  int32_t end;    str_variant.get_INFO_value_single_int(END_KEY, end);
  int32_t period; str_variant.get_INFO_value_single_int(PERIOD_KEY, period);
  std::vector<int32_t> bp_diffs; str_variant.get_INFO_value_multiple_ints(BPDIFFS_KEY, bp_diffs);
  assert(bp_diffs.size()+1 == str_variant.num_alleles());

  denovo_vcf_ << "BPDIFFS=" << bp_diffs[0];
  for (int i = 2; i < str_variant.num_alleles(); i++)
    denovo_vcf_ << "," <<  bp_diffs[i-1];
  denovo_vcf_ << ";START="  << start
	      << ";END="    << end
	      << ";PERIOD=" << period;

  // FORMAT field
  denovo_vcf_ << "\t" << "NOMUT:DENOVO:OTHER";
}

void TrioDenovoScanner::add_child_to_record(double total_ll_no_mutation, double total_ll_one_denovo, double total_ll_one_other){
  denovo_vcf_ << "\t" << total_ll_no_mutation << ":" << total_ll_one_denovo << ":" << total_ll_one_other;
}

void TrioDenovoScanner::scan(VCF::VCFReader& str_vcf, std::ostream& logger){
  VCF::Variant str_variant;
  while (str_vcf.get_next_variant(str_variant)){
    int num_alleles = str_variant.num_alleles();
    if (num_alleles <= 1)
      continue;

    int32_t start;  str_variant.get_INFO_value_single_int(START_KEY, start);
    int32_t end;    str_variant.get_INFO_value_single_int(END_KEY, end);
    logger << "Processing STR region " << str_variant.get_chromosome() << ":" << start << "-" << end << " with " << num_alleles << " alleles" << "\n";

    UnphasedGL unphased_gls(str_variant);
    MutationModel mut_model(str_variant);
    DiploidGenotypePrior* dip_gt_priors;
    if (use_pop_priors_)
      dip_gt_priors = new PopulationGenotypePrior(str_variant, families_);
    else
      dip_gt_priors = new UniformGenotypePrior(str_variant, families_);
    initialize_vcf_record(str_variant);
    const double LOG_ONE_FOURTH = -log10(4);
    const double LOG_TWO        = log10(2);

    logger << "\t" << "Computing log-likelihoods for mutation scenarios" << "\n";
    for (auto family_iter = families_.begin(); family_iter != families_.end(); family_iter++){
      bool scan_for_denovo = unphased_gls.has_sample(family_iter->get_mother()) && unphased_gls.has_sample(family_iter->get_father());
      for (auto child_iter = family_iter->get_children().begin(); child_iter != family_iter->get_children().end(); ++child_iter){
	if (!scan_for_denovo || !unphased_gls.has_sample(*child_iter)){
	  denovo_vcf_ << "\t" << ".";
	  continue;
	}

	// To accelerate computations, we will ignore configurations that make a neglible contribution (< 0.01%) to the total LL
	// For mutational scenarios, we aggregate 1/4*A^2*(A+1)^2*4*2*A values. Therefore, to ignore a configuration with LL=X:
	// X*A^3*(A+1)^2*2 < TOTAL/10000;
	// logX < log(TOTAL) - log(10000*A^3*(A+1)^2*2) = log(TOTAL) - [log(10000) + 3log(A) + 2log(A+1) + log(2)];
	double MIN_CONTRIBUTION   = 4 + 3*log10(num_alleles) + 2*log(num_alleles+1) + LOG_TWO;
	double ll_no_mutation_max = -DBL_MAX/2, ll_no_mutation_total = 0.0;
	double ll_one_denovo_max  = -DBL_MAX/2, ll_one_denovo_total  = 0.0;
	double ll_one_other_max   = -DBL_MAX/2, ll_one_other_total   = 0.0;
	int mother_gl_index       = unphased_gls.get_sample_index(family_iter->get_mother());
	int father_gl_index       = unphased_gls.get_sample_index(family_iter->get_father());
	int child_gl_index        = unphased_gls.get_sample_index(*child_iter);

	// Iterate over all maternal genotypes
	for (int mat_i = 0; mat_i < num_alleles; mat_i++){
	  for (int mat_j = 0; mat_j <= mat_i; mat_j++){
	    double mat_ll = dip_gt_priors->log_unphased_genotype_prior(mat_j, mat_i, family_iter->get_mother()) + unphased_gls.get_gl(mother_gl_index, mat_j, mat_i);

	    // Iterate over all paternal genotypes
	    for (int pat_i = 0; pat_i < num_alleles; pat_i++){
	      for (int pat_j = 0; pat_j <= pat_i; pat_j++){
		double pat_ll    = dip_gt_priors->log_unphased_genotype_prior(pat_j, pat_i, family_iter->get_father()) + unphased_gls.get_gl(father_gl_index, pat_j, pat_i);
		double config_ll = mat_ll + pat_ll + LOG_ONE_FOURTH;

		// Iterate over all 4 possible inheritance patterns for the child
		for (int mat_index = 0; mat_index < 2; ++mat_index){
		  int mat_allele = (mat_index == 0 ? mat_i : mat_j);
		  for (int pat_index = 0; pat_index < 2; ++pat_index){
		    int pat_allele = (pat_index == 0 ? pat_i : pat_j);

		    double no_mutation_config_ll = config_ll + unphased_gls.get_gl(child_gl_index, std::min(mat_allele, pat_allele), std::max(mat_allele, pat_allele));
		    update_streaming_log_sum_exp(no_mutation_config_ll, ll_no_mutation_max, ll_no_mutation_total);

		    // All putative mutations to the maternal allele
		    double max_ll_mat_mut = config_ll + unphased_gls.get_max_gl_allele_fixed(child_gl_index, pat_allele) + mut_model.max_log_prior_mutation(mat_allele);
		    if (max_ll_mat_mut > std::min(ll_one_denovo_max, ll_one_other_max)-MIN_CONTRIBUTION){
		      for (int mut_allele = 0; mut_allele < num_alleles; mut_allele++){
			if (mut_allele == mat_allele)
			  continue;
			double prob = config_ll + unphased_gls.get_gl(child_gl_index, std::min(mut_allele, pat_allele), std::max(mut_allele, pat_allele))
			  + mut_model.log_prior_mutation(mat_allele, mut_allele);
			if (mut_allele != mat_i && mut_allele != mat_j && mut_allele != pat_i && mut_allele != pat_j)
			  update_streaming_log_sum_exp(prob, ll_one_denovo_max, ll_one_denovo_total);
			else
			  update_streaming_log_sum_exp(prob, ll_one_other_max, ll_one_other_total);
		      }
		    }

		    // All putative mutations to the paternal allele
		    double max_ll_pat_mut = config_ll + unphased_gls.get_max_gl_allele_fixed(child_gl_index, mat_allele) + mut_model.max_log_prior_mutation(pat_allele);
		    if (max_ll_pat_mut > std::min(ll_one_denovo_max, ll_one_other_max)-MIN_CONTRIBUTION){
		      for (int mut_allele = 0; mut_allele < num_alleles; mut_allele++){
			if (mut_allele == pat_allele)
			  continue;
			double prob = config_ll + unphased_gls.get_gl(child_gl_index, std::min(mat_allele, mut_allele), std::max(mat_allele, mut_allele))
			  + mut_model.log_prior_mutation(pat_allele, mut_allele);
			if (mut_allele != mat_i && mut_allele != mat_j && mut_allele != pat_i && mut_allele != pat_j)
			  update_streaming_log_sum_exp(prob, ll_one_denovo_max, ll_one_denovo_total);
			else
			  update_streaming_log_sum_exp(prob, ll_one_other_max, ll_one_other_total);
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}

	// Compute total LL for each scenario and add it to the VCF
	double total_ll_no_mutation = finish_streaming_log_sum_exp(ll_no_mutation_max, ll_no_mutation_total);
	double total_ll_one_denovo  = finish_streaming_log_sum_exp(ll_one_denovo_max,  ll_one_denovo_total);
	double total_ll_one_other   = finish_streaming_log_sum_exp(ll_one_other_max,   ll_one_other_total);
	add_child_to_record(total_ll_no_mutation, total_ll_one_denovo, total_ll_one_other);
      }
    }

    // End of VCF record line
    denovo_vcf_ << "\n";
    delete dip_gt_priors;
  }
}
