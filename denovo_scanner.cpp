#include "denovo_scanner.h"
#include "error.h"
#include "haplotype_tracker.h"
#include "mathops.h"
#include "mutation_model.h"
#include "vcf_input.h"

#include <vector>

double log_genotype_prior(vcflib::Variant& str_variant, int gt_a, int gt_b, const std::string& sample_name){
  printErrorAndDie("log_genotype_prior() function not implemented");
  return 0.0;
}

void DenovoScanner::scan(vcflib::VariantCallFile& snp_vcf, vcflib::VariantCallFile& str_vcf, std::set<std::string>& sites_to_skip,
			 std::ostream& logger){
  HaplotypeTracker haplotype_tracker(families_);
  vcflib::Variant snp_variant(snp_vcf), str_variant(str_vcf);

  std::string chrom = "";
  int32_t num_strs = 0;
  while (str_vcf.getNextVariant(str_variant)){
    num_strs++;
    PhasedGL phased_gls(str_vcf, str_variant);

    if (str_variant.sequenceName.compare(chrom) != 0){
      chrom = str_variant.sequenceName;
      haplotype_tracker.reset();
      if(!snp_vcf.setRegion(chrom, 1))
	printErrorAndDie("Failed to set the region to chromosome " + chrom + " in the SNP VCF. Please check the SNP VCF and rerun the analysis");
    }

    int32_t start_of_window = str_variant.position - window_size_;
    int32_t end_of_window   = str_variant.position + window_size_;
    if (start_of_window < 0)
      start_of_window = 0;

    // Incorporate new SNPs within the window
    while (haplotype_tracker.last_snp_position() < end_of_window && snp_vcf.getNextVariant(snp_variant)){
      std::string key = snp_variant.sequenceName + ":" + std::to_string(snp_variant.position);
      if (sites_to_skip.find(key) != sites_to_skip.end())
	continue;
      haplotype_tracker.add_snp(snp_variant);
    }

    // Remove SNPs to left of window
    while (haplotype_tracker.next_snp_position() < start_of_window && haplotype_tracker.next_snp_position() != -1)
      haplotype_tracker.remove_next_snp();

    int num_alleles = str_variant.alleles.size();
    MutationModel mut_model(str_variant);

    // Analyze edit distances between the phased SNP haplotypes of each child and its parents
    int d11, d12, d21, d22;
    for (auto family_iter = families_.begin(); family_iter != families_.end(); family_iter++){

      // Determine if all samples have well-phased SNP haplotypes and infer the inheritance pattern
      bool scan_for_denovo = true;
      std::vector<int> maternal_indices, paternal_indices;
      for (auto child_iter = family_iter->get_children().begin(); child_iter != family_iter->get_children().end(); child_iter++){
	DiploidEditDistance maternal_distance = haplotype_tracker.edit_distances(*child_iter, family_iter->get_mother());
	int min_mat_dist, min_mat_index, second_mat_dist, second_mat_index;
	maternal_distance.min_distance(min_mat_dist, min_mat_index);
	maternal_distance.second_min_distance(second_mat_dist, second_mat_index);
	if (min_mat_dist > MAX_BEST_SCORE || second_mat_dist < MIN_SECOND_BEST_SCORE){
	  scan_for_denovo = false;
	  break;
	}

	DiploidEditDistance paternal_distance = haplotype_tracker.edit_distances(*child_iter, family_iter->get_father());
	int min_pat_dist, min_pat_index, second_pat_dist, second_pat_index;
	paternal_distance.min_distance(min_pat_dist, min_pat_index);
	paternal_distance.second_min_distance(second_pat_dist, second_pat_index);

	if (min_pat_dist > MAX_BEST_SCORE || second_pat_dist < MIN_SECOND_BEST_SCORE){
	  scan_for_denovo = false;
	  break;
	}

	// Ensure that one of the best matches involves haplotype #1 and the other involves haplotype #2
	if (min_mat_index == 0 || min_mat_index == 1){
	  if (min_pat_index != 2 && min_pat_index != 3){
	    scan_for_denovo = false;
	    break;
	  }
	}
	else {
	  if (min_pat_index != 0 && min_pat_index != 1){
	    scan_for_denovo = false;
	    break;
	  }
	}
	maternal_indices.push_back(min_mat_index);
	paternal_indices.push_back(min_pat_index);
	assert(maternal_indices.back() >= 0 && maternal_indices.back() < 4);
	assert(paternal_indices.back() >= 0 && paternal_indices.back() < 4);
      }

      if (!scan_for_denovo || num_alleles <= 1)
	1;
      else {
	std::vector<double> lls_no_denovo;
	std::vector< std::vector<double> > lls_one_denovo_mut(family_iter->get_children().size());
	std::vector< std::vector<double> > lls_one_other_mut(family_iter->get_children().size());

	// Iterate over all maternal genotypes
	for (int mat_i = 0; mat_i < num_alleles; mat_i++){
	  for (int mat_j = 0; mat_j < num_alleles; mat_j++){
	    double mat_ll = log_genotype_prior(str_variant, mat_i, mat_j, family_iter->get_mother()) + phased_gls.get_gl(family_iter->get_mother(), mat_i, mat_j);

	    // Iterate over all paternal genotypes
	    for (int pat_i = 0; pat_i < num_alleles; pat_i++){
	      for (int pat_j = 0; pat_j < num_alleles; pat_j++){
		double pat_ll = log_genotype_prior(str_variant, pat_i, pat_j, family_iter->get_father()) + phased_gls.get_gl(family_iter->get_father(), pat_i, pat_j);

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
