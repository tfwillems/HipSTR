#ifndef HAPLOTYPE_TRACKER_H_
#define HAPLOTYPE_TRACKER_H_

#include <deque>
#include <string>
#include <vector>

#include "vcflib/src/Variant.h"
#include "unused/Pedigree.h"

class DiploidHaplotype {
 private:
  std::deque<int64_t> snps_1_, snps_2_;
  int64_t erase_mask_, set_mask_;

  int num_set_bits(int64_t value){
    int count = 0;
    while (value != 0){
      value &= (value-1);
      count++;
    }
    return count;
  }

  int edit_distance(std::deque<int64_t>& snp_hap_1, std::deque<int64_t>& snp_hap_2){
    int distance = 0;
    auto iter_a = snp_hap_1.begin();
    auto iter_b = snp_hap_2.begin();
    while (iter_a != snp_hap_1.end()){
      distance += num_set_bits((*iter_a) ^ (*iter_b));
      iter_a++; iter_b++;
    }
    return distance;
  }

 public:  
  DiploidHaplotype(){
    snps_1_        = std::deque<int64_t>();
    snps_2_        = std::deque<int64_t>();
    snps_1_.push_back(0);
    snps_2_.push_back(0);
    erase_mask_    = -2;
    set_mask_      = 1;
  }

  void edit_distances(DiploidHaplotype& other_hap, int& d11, int& d12, int& d21, int& d22){
    d11 = edit_distance(snps_1_, other_hap.snps_1_);
    d12 = edit_distance(snps_1_, other_hap.snps_2_);
    d21 = edit_distance(snps_2_, other_hap.snps_1_);
    d22 = edit_distance(snps_2_, other_hap.snps_2_);
  }

  void add_snp(int gt_a, int gt_b);

  void remove_next_snp();
};







class HaplotypeTracker {
 private:
  std::vector<NuclearFamily> families_;
  std::vector<std::string> samples_;
  std::map<std::string, int> sample_indices_;
  std::vector<DiploidHaplotype> snp_haplotypes_;

  int32_t num_snps_;
  std::queue<int32_t> positions_;

 public:  
  HaplotypeTracker(std::vector<NuclearFamily>& families){
    families_ = families;
    samples_  = std::vector<std::string>();
    for (auto family_iter = families.begin(); family_iter != families.end(); family_iter++){
      samples_.push_back(family_iter->get_mother());
      samples_.push_back(family_iter->get_father());
      for (auto child_iter = family_iter->get_children().begin(); child_iter != family_iter->get_children().end(); child_iter++)
	samples_.push_back(*child_iter);
    }    
    for (unsigned int i = 0; i < samples_.size(); i++)
      sample_indices_[samples_[i]] = i;
    
    snp_haplotypes_ = std::vector<DiploidHaplotype>(samples_.size(), DiploidHaplotype());
    num_snps_       = 0;
    positions_      = std::queue<int32_t>();
  }

  int32_t next_snp_position(){
    if (num_snps_ == 0)
      return -1;
    return positions_.front();
  }

  void remove_next_snp(){
    if (num_snps_ == 0)
      return;
    num_snps_--;
    for (unsigned int i = 0; i < snp_haplotypes_.size(); i++)
      snp_haplotypes_[i].remove_next_snp();
    positions_.pop();
  }

  void add_snp(vcflib::Variant& variant);

  int32_t num_stored_snps() { return num_snps_; }

  void edit_distances(const std::string& sample_1, const std::string& sample_2,
		      int& d11, int& d12, int& d21, int& d22){
    int index_1 = sample_indices_[sample_1];
    int index_2 = sample_indices_[sample_2];
    snp_haplotypes_[index_1].edit_distances(snp_haplotypes_[index_2], d11, d12, d21, d22);
  }
};

#endif
