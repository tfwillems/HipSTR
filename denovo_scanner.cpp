#include <stdlib.h>

#include "denovo_scanner.h"
#include "error.h"
#include "haplotype_tracker.h"
#include "mathops.h"
#include "mutation_model.h"
#include "vcf_input.h"

#include <vector>

void DiploidGenotypePrior::compute_allele_freqs(vcflib::Variant& variant, std::vector<NuclearFamily>& families){
  allele_freqs_ = std::vector<double>(num_alleles_, 1.0); // Use a one sample pseudocount

  // Iterate over all founders in the families to compute allele counts
  double total_count = num_alleles_;
  for (auto family_iter = families.begin(); family_iter != families.end(); family_iter++){
    for (int i = 0; i < 2; i++){
      std::string sample = (i == 0 ? family_iter->get_mother() : family_iter->get_father());
      std::string gts = variant.getGenotype(sample);
      if (gts.size() == 0)
	continue;

      size_t separator_index = gts.find("|");
      if (separator_index == std::string::npos)
	separator_index = gts.find("/");
      if (separator_index == std::string::npos || separator_index+1 == gts.size())
	printErrorAndDie("Failed to find valid separator in genotype: " + gts);
      int gt_a = std::atoi(gts.substr(0, separator_index).c_str());
      int gt_b = std::atoi(gts.substr(separator_index+1).c_str());
      allele_freqs_[gt_a]++;
      allele_freqs_[gt_b]++;
      total_count += 2;
    }
  }

  // Normalize the allele counts to obtain frequencies
  for (int i = 0; i < allele_freqs_.size(); i++)
    allele_freqs_[i] /= total_count;

  // Precompute the logs of the allele frequencies
  log_allele_freqs_.clear();
  for (int i = 0; i < allele_freqs_.size(); i++)
    log_allele_freqs_.push_back(log10(allele_freqs_[i]));
}


void DenovoScanner::scan(std::string& snp_vcf_file, vcflib::VariantCallFile& str_vcf, std::set<std::string>& sites_to_skip,
			 std::ostream& logger){
  HaplotypeTracker haplotype_tracker(families_, snp_vcf_file, window_size_);
  vcflib::Variant str_variant(str_vcf);
  int32_t num_strs  = 0;
  while (str_vcf.getNextVariant(str_variant)){
    num_strs++;
    PhasedGL phased_gls(str_vcf, str_variant);
    haplotype_tracker.advance(str_variant.sequenceName, str_variant.position, sites_to_skip, logger);

    int num_alleles = str_variant.alleles.size();
    if (num_alleles <= 1)
      continue;

    MutationModel mut_model(str_variant);
    DiploidGenotypePrior dip_gt_priors(str_variant, families_);

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
	1;
      else {
	std::vector<double> lls_no_denovo;
	std::vector< std::vector<double> > lls_one_denovo_mut(family_iter->get_children().size());
	std::vector< std::vector<double> > lls_one_other_mut(family_iter->get_children().size());

	// Iterate over all maternal genotypes
	for (int mat_i = 0; mat_i < num_alleles; mat_i++){
	  for (int mat_j = 0; mat_j < num_alleles; mat_j++){
	    double mat_ll = dip_gt_priors.log_phased_genotype_prior(mat_i, mat_j, family_iter->get_mother()) + phased_gls.get_gl(family_iter->get_mother(), mat_i, mat_j);

	    // Iterate over all paternal genotypes
	    for (int pat_i = 0; pat_i < num_alleles; pat_i++){
	      for (int pat_j = 0; pat_j < num_alleles; pat_j++){
		double pat_ll = dip_gt_priors.log_phased_genotype_prior(pat_i, pat_j, family_iter->get_father()) + phased_gls.get_gl(family_iter->get_father(), pat_i, pat_j);

		double no_denovo_config_ll = mat_ll + pat_ll;

		// Iterate over all the children to compute the likelihood for no denovos
		int child_index = 0;
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

		  no_denovo_config_ll += phased_gls.get_gl(*child_iter, child_i, child_j);
		}
		lls_no_denovo.push_back(no_denovo_config_ll);

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

		  double config_ll = no_denovo_config_ll - phased_gls.get_gl(*child_iter, child_i, child_j);

		  // All putative mutations on haplotype #1
		  for (int mut_allele = 0; mut_allele < num_alleles; mut_allele++){
		    if (mut_allele == child_i)
		      continue;
		    double prob = config_ll + phased_gls.get_gl(*child_iter, mut_allele, child_j) + mut_model.log_prior_mutation(child_i, mut_allele);
		    if (mut_allele != mat_i && mut_allele != mat_j && mut_allele != pat_i && mut_allele != pat_j)
		      lls_one_denovo_mut[child_index].push_back(prob);
		    else
		      lls_one_other_mut[child_index].push_back(prob);
		  }

		  // All putative mutations on haplotype #2
		  for (int mut_allele = 0; mut_allele < num_alleles; mut_allele++){
		    if (mut_allele == child_j)
		      continue;
		    double prob = config_ll + phased_gls.get_gl(*child_iter, child_i, mut_allele) + mut_model.log_prior_mutation(child_j, mut_allele);
		    if (mut_allele != mat_i && mut_allele != mat_j && mut_allele != pat_i && mut_allele != pat_j)
		      lls_one_denovo_mut[child_index].push_back(prob);
		    else
		      lls_one_other_mut[child_index].push_back(prob);
		  }
		}
	      }
	    }
	  }
	}

	// Compute total LL for each scenario
	double total_ll_no_denovo = log_sum_exp(lls_no_denovo);
	std::vector<double> total_lls_one_denovo, total_lls_one_other;
	int child_index = 0;
	for (auto child_iter = family_iter->get_children().begin(); child_iter != family_iter->get_children().end(); ++child_iter, ++child_index){
	  total_lls_one_denovo.push_back(log_sum_exp(lls_one_denovo_mut[child_index]));
	  total_lls_one_other.push_back(log_sum_exp(lls_one_other_mut[child_index]));
	}
      }
      std::cout << str_variant.position << " " << (scan_for_denovo ? 1 : 0) << "\n";
    }
  }
}
