#include "haplotype_tracker.h"

std::ostream& operator<< (std::ostream &out, DiploidEditDistance& edit_distance){
  out << "\t" << edit_distance.distances_[0] << " " << edit_distance.distances_[1]
      << " "  << edit_distance.distances_[2] << " " << edit_distance.distances_[3];
  return out;
}

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

void DiploidHaplotype::remove_next_snp(){
  assert(snps_1_.size() > 0 && snps_2_.size() > 0);
  snps_1_.front() &= erase_mask_;
  snps_2_.front() &= erase_mask_;
  erase_mask_ <<= 1;
  if (erase_mask_ == 0){
    snps_1_.pop_front();
    snps_2_.pop_front();
    erase_mask_ = -2;
  } 
}

void HaplotypeTracker::add_snp(vcflib::Variant& variant){
  num_snps_++;
  positions_.push(variant.position);
  
  int sample_index = 0;
  for (unsigned int i = 0; i < families_.size(); i++){
    NuclearFamily& family = families_[i];
    bool use_gts          = true;
    
    if (family.is_missing_genotype(variant))
      use_gts = false; // Ignore a SNP if any samples in the family are missing a genotype
    else if (!family.is_mendelian(variant))
      use_gts = false; // Ignore a SNP if any samples in the family have a Mendelian inconsistency


    for (int j = 0; j < family.size(); j++){
      if (use_gts){
	std::string gt = variant.getGenotype(samples_[sample_index]);
	assert(gt.size() == 3);
	assert(gt[1] == '|');
	int gt_a = gt[0]-'0', gt_b = gt[2]-'0';
	snp_haplotypes_[sample_index].add_snp(gt_a, gt_b);
      }
      else
	snp_haplotypes_[sample_index].add_snp(0, 0);
      sample_index++;
    }      
  }
}
