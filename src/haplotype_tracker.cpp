#include "haplotype_tracker.h"

void DiploidHaplotype::add_snp(int gt_a, int gt_b){
  assert(snps_1_.size() > 0 && snps_2_.size() > 0);
  if (gt_a == 1) snps_1_.back() |= set_mask_;
  if (gt_b == 1) snps_2_.back() |= set_mask_;
  
  set_mask_ <<= 1;
  if (set_mask_ == MAX_DIGIT){
    snps_1_.push_back(0);
    snps_2_.push_back(0);
    set_mask_ = 1;
  }
}

void DiploidHaplotype::extract_set_bits(int64_t value, int offset, std::set<int>& mismatch_indices) const {
  while (value){
    if (value & 1)
      mismatch_indices.insert(offset);
    value >>= 1;
    offset++;
  }
}

void DiploidHaplotype::add_mismatched_sites(int hap_index, DiploidHaplotype& other_hap, int other_index,
					    std::set<int>& mismatch_indices) const {
  assert((hap_index == 0 || hap_index == 1) && (other_index == 0 || other_index == 1));
  const std::deque<int64_t>& hap_a = (hap_index == 0 ? snps_1_ : snps_2_);
  const std::deque<int64_t>& hap_b = (other_index == 0 ? other_hap.snps_1_ : other_hap.snps_2_);
  assert(hap_a.size() == hap_b.size());
  auto iter_a = hap_a.begin();
  auto iter_b = hap_b.begin();

  int offset = -num_removed_;
  while (iter_a != hap_a.end()){
    int64_t set_bits = (*iter_a)^(*iter_b);
    if (set_bits)
      extract_set_bits(set_bits, offset, mismatch_indices);
    offset += 63;
    iter_a++; iter_b++;
  }
}

void DiploidHaplotype::remove_next_snp(){
  assert(snps_1_.size() > 0 && snps_2_.size() > 0);
  snps_1_.front() &= erase_mask_;
  snps_2_.front() &= erase_mask_;
  erase_mask_ <<= 1;
  num_removed_++;
  if (erase_mask_ == 0){
    snps_1_.pop_front();
    snps_2_.pop_front();
    erase_mask_  = -2;
    num_removed_ = 0;
  } 
}

void HaplotypeTracker::add_snp(const VCF::Variant& variant){
  num_snps_++;
  positions_.push_back(variant.get_position());
  
  int sample_index = 0;
  for (unsigned int i = 0; i < families_.size(); i++){
    NuclearFamily& family = families_[i];
    bool use_gts          = true;
    
    if (family.is_missing_genotype(variant))
      use_gts = false; // Ignore a SNP if any samples in the family are missing a genotype
    else if (!family.is_mendelian(variant))
      use_gts = false; // Ignore a SNP if any samples in the family have a Mendelian inconsistency

    int gt_a, gt_b;
    for (int j = 0; j < family.size(); j++){
      if (use_gts){
	variant.get_genotype(vcf_indices_[sample_index], gt_a, gt_b);
	snp_haplotypes_[sample_index].add_snp(gt_a, gt_b);
      }
      else
	snp_haplotypes_[sample_index].add_snp(0, 0);
      sample_index++;
    }      
  }
}

void HaplotypeTracker::advance(const std::string& chrom, int32_t position, const std::set<std::string>& sites_to_skip){
  int32_t start_of_window = (position >= window_size_ ? position - window_size_ : 0);
  int32_t end_of_window   = position + window_size_;
  if (chrom.compare(chrom_) != 0){
    chrom_ = chrom;
    reset();
    if (!snp_vcf_.set_region(chrom, start_of_window))
      printErrorAndDie("Failed to set the region to chromosome " + chrom + " in the SNP VCF. Please check the SNP VCF and rerun the analysis");
  }
  else {
    if (start_of_window < prev_window_start_)
      printErrorAndDie("Haplotype tracker's advance() function can only consecutively process regions sorted by position");
    if (start_of_window > prev_window_end_){
      reset();
      if (!snp_vcf_.set_region(chrom, start_of_window))
	printErrorAndDie("Failed to set the region in the SNP VCF. Please check the SNP VCF and rerun the analysis");
    }
  }
  prev_window_start_ = start_of_window;
  prev_window_end_   = end_of_window;

  // Incorporate new SNPs within the window
  VCF::Variant snp_variant;
  while (last_snp_position() < end_of_window && snp_vcf_.get_next_variant(snp_variant)){
    std::stringstream ss;
    ss << snp_variant.get_chromosome() << ":" << snp_variant.get_position();
    std::string key = ss.str();
    if (sites_to_skip.find(key) != sites_to_skip.end())
      continue;
    add_snp(snp_variant);
  }

  // Remove SNPs to left of window
  while (next_snp_position() < start_of_window && next_snp_position() != -1)
    remove_next_snp();
  //logger << " done" << std::endl;
}

/* Analyze edit distances between the phased SNP haplotypes of each child and its parents. Returns true iff all the children in the family
 * have a valid match, which is controlled by the parameters MAX_BEST_SCORE and MIN_SECOND_BEST_SCORE.
 * For a given child, a valid match occurs if a child's haplotype matches a parental haplotype with distance <= MAX_BEST_SCORE, no other child-parent haplotype
 * pairings have a distance < MIN_SECOND_BEST_SCORE, and the optimal maternal and paternal matches involve both of the child's haplotypes.
 *
 * If all children satisfy these conditions, the method stores the inheritance pattern for each child in the provided vectors,
 * where 0, 1, 2 and 3 denote that the child-parent haplotype pairing is 1+1, 1+2, 2+1 or 2+2, respectively.
 * Stores the positions of any SNPs that are inconsistent with the inhertiance pattern in the set.
 * NOTE: this set can ONLY contain SNPs that are Mendelian and have no missing genotypes in the family.
 */
bool HaplotypeTracker::infer_haplotype_inheritance(const NuclearFamily& family, int max_best_score, int min_second_best_score,
						   std::vector<int>& maternal_indices, std::vector<int>& paternal_indices, std::set<int32_t>& bad_sites){
  assert(maternal_indices.size() == 0 && paternal_indices.size() == 0);
  DiploidHaplotype& mat_haplotypes = snp_haplotypes_[sample_indices_[family.get_mother()]];
  DiploidHaplotype& pat_haplotypes = snp_haplotypes_[sample_indices_[family.get_father()]];
  std::set<int> mismatch_indices;

  for (auto child_iter = family.get_children().begin(); child_iter != family.get_children().end(); child_iter++){
    DiploidEditDistance maternal_distance = edit_distances(*child_iter, family.get_mother());
    int min_mat_dist, min_mat_index, second_mat_dist, second_mat_index;
    maternal_distance.min_distance(min_mat_dist, min_mat_index);
    maternal_distance.second_min_distance(second_mat_dist, second_mat_index);
    if (min_mat_dist > max_best_score || second_mat_dist < min_second_best_score)
      return false;

    DiploidEditDistance paternal_distance = edit_distances(*child_iter, family.get_father());
    int min_pat_dist, min_pat_index, second_pat_dist, second_pat_index;
    paternal_distance.min_distance(min_pat_dist, min_pat_index);
    paternal_distance.second_min_distance(second_pat_dist, second_pat_index);
    if (min_pat_dist > max_best_score || second_pat_dist < min_second_best_score)
      return false;

    // Ensure that one of the parental best matches involves haplotype #1 and the other involves haplotype #2
    if (min_mat_index == 0 || min_mat_index == 1){
      if (min_pat_index != 2 && min_pat_index != 3)
	return false;
    }
    else if (min_pat_index != 0 && min_pat_index != 1)
      return false;
    assert(min_mat_index >= 0 && min_mat_index < 4 && min_pat_index >= 0 && min_pat_index < 4);

    // Identify the indices of sites that are inconsistent with the inheritance structure
    // Only identifies sites that are Mendelian and have no missing genotypes
    DiploidHaplotype& child_haplotypes = snp_haplotypes_[sample_indices_[*child_iter]];
    int idx_a = (min_mat_index == 0 || min_mat_index == 1 ? 0 : 1);
    int idx_b = (min_mat_index == 0 || min_mat_index == 2 ? 0 : 1);
    child_haplotypes.add_mismatched_sites(idx_a, mat_haplotypes, idx_b, mismatch_indices);
    idx_a = (min_pat_index == 0 || min_pat_index == 1 ? 0 : 1);
    idx_b = (min_pat_index == 0 || min_pat_index == 2 ? 0 : 1);
    child_haplotypes.add_mismatched_sites(idx_a, pat_haplotypes, idx_b, mismatch_indices);

    // Store the best indices
    maternal_indices.push_back(min_mat_index);
    paternal_indices.push_back(min_pat_index);
  }

  // Convert from internal SNP indices to SNP positions
  for (auto snp_index_iter = mismatch_indices.begin(); snp_index_iter != mismatch_indices.end(); snp_index_iter++)
    bad_sites.insert(positions_.at(*snp_index_iter));
  return true;
}
