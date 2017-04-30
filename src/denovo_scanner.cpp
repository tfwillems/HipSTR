#include <stdlib.h>

#include <cfloat>
#include <vector>

#include "denovo_scanner.h"
#include "denovo_allele_priors.h"
#include "error.h"
#include "haplotype_tracker.h"
#include "mathops.h"
#include "mutation_model.h"
#include "vcf_input.h"

std::string DenovoScanner::BPDIFFS_KEY = "BPDIFFS";
std::string DenovoScanner::START_KEY   = "START";
std::string DenovoScanner::END_KEY     = "END";
std::string DenovoScanner::PERIOD_KEY  = "PERIOD";

void DenovoScanner::write_vcf_header(const std::string& full_command){
  denovo_vcf_ << "##fileformat=VCFv4.1" << "\n"
	      << "##command=" << full_command << "\n";

  // Info field descriptors
  denovo_vcf_ << "##INFO=<ID=" << "BPDIFFS"        << ",Number=A,Type=Integer,Description=\"" << "Base pair difference of each alternate allele from the reference allele"      << "\">\n"
	      << "##INFO=<ID=" << "START"          << ",Number=1,Type=Integer,Description=\"" << "Inclusive start coodinate for the repetitive portion of the reference allele" << "\">\n"
	      << "##INFO=<ID=" << "END"            << ",Number=1,Type=Integer,Description=\"" << "Inclusive end coordinate for the repetitive portion of the reference allele"  << "\">\n"
	      << "##INFO=<ID=" << "PERIOD"         << ",Number=1,Type=Integer,Description=\"" << "Length of STR motif"                                                          << "\">\n";

  // Format field descriptors
  denovo_vcf_ << "##FORMAT=<ID=" << "CHILDREN" << ",Number=.,Type=String,Description=\""  << "Ordered list of children in family that were tested for mutations. Specifies order of values for AFF, DENOVO and OTHER FORMAT fields" << "\">\n"
	      << "##FORMAT=<ID=" << "NOMUT"    << ",Number=1,Type=Float,Description=\""   << "Log10-likelihood that no mutations occurred in any of the family members" << "\">\n"
	      << "##FORMAT=<ID=" << "ANYMUT"   << ",Number=1,Type=Float,Description=\""   << "Log10-likelihood that a mutation occurred in any of the family members"   << "\">\n"
	      << "##FORMAT=<ID=" << "DENOVO"   << ",Number=.,Type=Float,Description=\""
	      << "Log10-likelihood that a single de novo mutation occurred in the family, and it occurred in the current child" << "\">\n"
	      << "##FORMAT=<ID=" << "OTHER"    << ",Number=.,Type=Float,Description=\""
	      << "Log10-likelihood that a single mutation occurred in the family, and it occurred in the current child. In contrast to DENOVO, the mutated allele is also present in a parental genotype" << "\">\n";

  denovo_vcf_ << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  for (auto family_iter = families_.begin(); family_iter != families_.end(); family_iter++)
    denovo_vcf_ << "\t" << family_iter->get_family_id();
  denovo_vcf_ << "\n";
}


void DenovoScanner::initialize_vcf_record(const VCF::Variant& str_variant){
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

  std::vector<int32_t> bp_diffs;
  if (str_variant.num_alleles() > 2)
    str_variant.get_INFO_value_multiple_ints(BPDIFFS_KEY, bp_diffs);
  else {
    int32_t diff;
    str_variant.get_INFO_value_single_int(BPDIFFS_KEY, diff);
    bp_diffs.push_back(diff);
  }
  assert(bp_diffs.size()+1 == str_variant.num_alleles());

  denovo_vcf_ << "BPDIFFS=" << bp_diffs[0];
  for (int i = 2; i < str_variant.num_alleles(); i++)
    denovo_vcf_ << "," <<  bp_diffs[i-1];
  denovo_vcf_ << ";START="  << start
	      << ";END="    << end
	      << ";PERIOD=" << period;

  // FORMAT field
  denovo_vcf_ << "\t" << "CHILDREN:NOMUT:ANYMUT:DENOVO:OTHER";
}

void DenovoScanner::add_family_to_record(const NuclearFamily& family, double total_ll_no_mutation,
					 const std::vector<double>& total_lls_one_denovo, const std::vector<double>& total_lls_one_other){
  assert(total_lls_one_denovo.size() == total_lls_one_other.size() && total_lls_one_denovo.size() == family.get_children().size());
  const std::vector<std::string>& children = family.get_children();

  // Names of children
  denovo_vcf_ << "\t" << children.at(0);
  for (int i = 1; i < children.size(); i++)
    denovo_vcf_ << "," << children[i];

  // LL no mutation
  denovo_vcf_ << ":" << total_ll_no_mutation;

  // LL a mutation
  denovo_vcf_ << ":" << fast_log_sum_exp(fast_log_sum_exp(total_lls_one_denovo), fast_log_sum_exp(total_lls_one_other));

  // LL one denovo, for each child
  denovo_vcf_ << ":" << total_lls_one_denovo.at(0);
  for (int i = 1; i < total_lls_one_denovo.size(); i++)
    denovo_vcf_ << "," << total_lls_one_denovo[i];

  // LL one other mutation, for each child
  denovo_vcf_ << ":" << total_lls_one_other.at(0);
  for (int i = 1; i < total_lls_one_other.size(); i++)
    denovo_vcf_ << "," << total_lls_one_other[i];
}

void DenovoScanner::scan(const std::string& snp_vcf_file, VCF::VCFReader& str_vcf, const std::set<std::string>& sites_to_skip,
			 std::ostream& logger){
  HaplotypeTracker haplotype_tracker(families_, snp_vcf_file, window_size_);
  VCF::Variant str_variant;
  while (str_vcf.get_next_variant(str_variant)){
    int num_alleles = str_variant.num_alleles();
    if (num_alleles <= 1)
      continue;
    if (str_variant.num_samples() == str_variant.num_missing())
      continue;

    int32_t start;  str_variant.get_INFO_value_single_int(START_KEY, start);
    int32_t end;    str_variant.get_INFO_value_single_int(END_KEY, end);
    logger << "Processing STR region " << str_variant.get_chromosome() << ":" << start << "-" << end << " with " << num_alleles << " alleles" << "\n";

    PhasedGL phased_gls(str_variant);
    logger << "\t";
    haplotype_tracker.advance(str_variant.get_chromosome(), str_variant.get_position(), sites_to_skip);

    MutationModel mut_model(str_variant);
    DiploidGenotypePrior* dip_gt_priors;
    if (use_pop_priors_)
      dip_gt_priors = new PopulationGenotypePrior(str_variant, families_);
    else
      dip_gt_priors = new UniformGenotypePrior(str_variant, families_);
    initialize_vcf_record(str_variant);

    logger << "\t" << "Computing log-likelihoods for mutation scenarios" << "\n";
    for (auto family_iter = families_.begin(); family_iter != families_.end(); family_iter++){
      // Determine if all samples have well-phased SNP haplotypes and infer the inheritance pattern
      std::vector<int> maternal_indices, paternal_indices;
      std::set<int32_t> bad_sites;
      bool scan_for_denovo = haplotype_tracker.infer_haplotype_inheritance(*family_iter, MAX_BEST_SCORE, MIN_SECOND_BEST_SCORE,
									   maternal_indices, paternal_indices, bad_sites);

      // Don't look for de novos if any of the family members are missing genotype likelihoods
      scan_for_denovo &= phased_gls.has_sample(family_iter->get_mother());
      scan_for_denovo &= phased_gls.has_sample(family_iter->get_father());
      if (scan_for_denovo)
	for (auto child_iter = family_iter->get_children().begin(); child_iter != family_iter->get_children().end(); ++child_iter)
	  scan_for_denovo &= phased_gls.has_sample(*child_iter);

      if (!scan_for_denovo)
	denovo_vcf_ << "\t" << ".";
      else {
	// To accelerate computations, we will ignore configurations that make a neglible contribution (< 0.01%) to the total LL
	// For mutational scenarios, we aggregate A^5*2*NUM_CHILDREN values. Therefore, to ignore a configuration with LL=X:
	// X*A^5*2*NUM_CHILDREN < TOTAL/10000;
	// logX < log(TOTAL) - log(10000*A^5*2*NUM_CHILDREN) = log(TOTAL) - [log(10000) + 5log(A) + log(2) + log(NUM_CHILDREN)];
	float MIN_CONTRIBUTION = 4 + 5*log10(num_alleles) + log10(2) + log10(family_iter->get_children().size());

	assert(family_iter->get_children().size() == maternal_indices.size() && maternal_indices.size() == paternal_indices.size());
	double ll_no_mutation_max = -DBL_MAX/2, ll_no_mutation_total = 0.0;
	std::vector<double> ll_one_denovo_max(family_iter->get_children().size(), -DBL_MAX/2), ll_one_denovo_total(family_iter->get_children().size(), 0.0);
	std::vector<double>  ll_one_other_max(family_iter->get_children().size(), -DBL_MAX/2),  ll_one_other_total(family_iter->get_children().size(), 0.0);

	int mother_gl_index = phased_gls.get_sample_index(family_iter->get_mother());
	int father_gl_index = phased_gls.get_sample_index(family_iter->get_father());
	std::vector<int> children_gl_index;
	for (auto child_iter = family_iter->get_children().begin(); child_iter != family_iter->get_children().end(); ++child_iter)
	  children_gl_index.push_back(phased_gls.get_sample_index(*child_iter));

	// Iterate over all maternal genotypes
	for (int mat_i = 0; mat_i < num_alleles; mat_i++){
	  for (int mat_j = 0; mat_j < num_alleles; mat_j++){
	    double mat_ll = dip_gt_priors->log_phased_genotype_prior(mat_i, mat_j, family_iter->get_mother()) + phased_gls.get_gl(mother_gl_index, mat_i, mat_j);

	    // Iterate over all paternal genotypes
	    for (int pat_i = 0; pat_i < num_alleles; pat_i++){
	      for (int pat_j = 0; pat_j < num_alleles; pat_j++){
		double pat_ll = dip_gt_priors->log_phased_genotype_prior(pat_i, pat_j, family_iter->get_father()) + phased_gls.get_gl(father_gl_index, pat_i, pat_j);

		double no_mutation_config_ll = mat_ll + pat_ll;

		// Iterate over all the children to compute the likelihood for no denovos
		int child_index = 0;
		for (auto child_iter = family_iter->get_children().begin(); child_iter != family_iter->get_children().end(); ++child_iter, ++child_index){
		  int child_i = -1, child_j = -1;

		  if (maternal_indices[child_index] == 0)       child_i = mat_i;
		  else if (maternal_indices[child_index] == 1)  child_i = mat_j;
		  else if (maternal_indices[child_index] == 2)  child_j = mat_i;
		  else                                          child_j = mat_j;

		  if (paternal_indices[child_index] == 0)       child_i = pat_i;
		  else if (paternal_indices[child_index] == 1)  child_i = pat_j;
		  else if (paternal_indices[child_index] == 2)  child_j = pat_i;
		  else                                          child_j = pat_j;

		  assert(child_i != -1 && child_j != -1);
		  no_mutation_config_ll += phased_gls.get_gl(children_gl_index[child_index], child_i, child_j);
		}
		update_streaming_log_sum_exp(no_mutation_config_ll, ll_no_mutation_max, ll_no_mutation_total);

		// Iterate over all the children to compute the likelihood that a single mutation occurs, and it occurs in the current child
		child_index = 0;
                for (auto child_iter = family_iter->get_children().begin(); child_iter != family_iter->get_children().end(); ++child_iter, ++child_index){
		  int child_i, child_j;

		  if (maternal_indices[child_index] == 0)       child_i = mat_i;
		  else if (maternal_indices[child_index] == 1)  child_i = mat_j;
		  else if (maternal_indices[child_index] == 2)  child_j = mat_i;
		  else                                          child_j = mat_j;

		  if (paternal_indices[child_index] == 0)       child_i = pat_i;
		  else if (paternal_indices[child_index] == 1)  child_i = pat_j;
		  else if (paternal_indices[child_index] == 2)  child_j = pat_i;
		  else                                          child_j = pat_j;

		  double config_ll = no_mutation_config_ll - phased_gls.get_gl(children_gl_index[child_index], child_i, child_j);

		  // All putative mutations on haplotype #1 (if they can contribute to the total LL)
		  double max_ll_hap_one = config_ll + phased_gls.get_max_gl_allele_two_fixed(children_gl_index[child_index], child_j) + mut_model.max_log_prior_mutation(child_i);
		  if (max_ll_hap_one > std::min(ll_one_denovo_max[child_index], ll_one_other_max[child_index])-MIN_CONTRIBUTION){
		    for (int mut_allele = 0; mut_allele < num_alleles; mut_allele++){
		      if (mut_allele == child_i)
			continue;
		      double prob = config_ll + phased_gls.get_gl(children_gl_index[child_index], mut_allele, child_j) + mut_model.log_prior_mutation(child_i, mut_allele);
		      if (mut_allele != mat_i && mut_allele != mat_j && mut_allele != pat_i && mut_allele != pat_j)
			update_streaming_log_sum_exp(prob, ll_one_denovo_max[child_index], ll_one_denovo_total[child_index]);
		      else
			update_streaming_log_sum_exp(prob, ll_one_other_max[child_index], ll_one_other_total[child_index]);
		    }
		  }

		  // All putative mutations on haplotype #2 (if they can contribute to the total LL)
		  double max_ll_hap_two = config_ll + phased_gls.get_max_gl_allele_one_fixed(children_gl_index[child_index], child_i) + mut_model.max_log_prior_mutation(child_j);
		  if (max_ll_hap_two > std::min(ll_one_denovo_max[child_index], ll_one_other_max[child_index])-MIN_CONTRIBUTION){
		    for (int mut_allele = 0; mut_allele < num_alleles; mut_allele++){
		      if (mut_allele == child_j)
			continue;
		      double prob = config_ll + phased_gls.get_gl(children_gl_index[child_index], child_i, mut_allele) + mut_model.log_prior_mutation(child_j, mut_allele);
		      if (mut_allele != mat_i && mut_allele != mat_j && mut_allele != pat_i && mut_allele != pat_j)
			update_streaming_log_sum_exp(prob, ll_one_denovo_max[child_index], ll_one_denovo_total[child_index]);
		      else
			update_streaming_log_sum_exp(prob, ll_one_other_max[child_index], ll_one_other_total[child_index]);
		    }
		  }
		}
	      }
	    }
	  }
	}

	// Compute total LL for each scenario
	double total_ll_no_mutation = finish_streaming_log_sum_exp(ll_no_mutation_max, ll_no_mutation_total);

	std::vector<double> total_lls_one_denovo, total_lls_one_other;
	int child_index = 0;
	for (auto child_iter = family_iter->get_children().begin(); child_iter != family_iter->get_children().end(); ++child_iter, ++child_index){
	  total_lls_one_denovo.push_back(finish_streaming_log_sum_exp(ll_one_denovo_max[child_index], ll_one_denovo_total[child_index]));
	  total_lls_one_other.push_back(finish_streaming_log_sum_exp(ll_one_other_max[child_index], ll_one_other_total[child_index]));
	}

	// Add family's mutation likelihoods to the VCF record
	add_family_to_record(*family_iter, total_ll_no_mutation, total_lls_one_denovo, total_lls_one_other);
      }
    }

    // End of VCF record line
    denovo_vcf_ << "\n";
    delete dip_gt_priors;
  }
}
