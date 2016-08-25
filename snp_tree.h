#ifndef SNP_TREE_H_
#define SNP_TREE_H_

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "haplotype_tracker.h"

#include "vcflib/src/Variant.h"

class SNP {
 private:
  uint32_t pos_;
  char base_1_;
  char base_2_;
 public:
  SNP(){
    pos_    = 0;
    base_1_ = 'N';
    base_2_ = 'N';
  }

  SNP(uint32_t pos, char base_1, char base_2){
    pos_    = pos;
    base_1_ = base_1;
    base_2_ = base_2;
  }
  
  friend std::ostream& operator<< (std::ostream &out, SNP& snp);

  inline uint32_t pos()  const { return pos_;    }
  inline char base_one() const { return base_1_; }
  inline char base_two() const { return base_2_; }
};

class SNPSorter {
 public:
  bool operator() (const SNP& snp_a, const SNP& snp_b){
    return snp_a.pos() < snp_b.pos();
  }
};

class SNPTree {
    std::vector<SNP> snps_;
    SNPTree* left_;
    SNPTree* right_;
    uint32_t center_;

 public:
    SNPTree(){
      left_   = right_ = NULL;
      center_ = 0;
    }

    SNPTree(const SNPTree& other) {
        center_ = other.center_;
        snps_   = other.snps_;
	left_   = (other.left_  ? new SNPTree(*other.left_)  : NULL);
	right_  = (other.right_ ? new SNPTree(*other.right_) : NULL);
    }
    
    SNPTree& operator=(const SNPTree& other) {
      center_ = other.center_;
      snps_   = other.snps_;
      if (other.left_)
	left_ = new SNPTree(*other.left_);
      else {
	if (left_) delete left_;
	left_ = NULL;
      }

      if (other.right_)
	right_ = new SNPTree(*other.right_);
      else {
	if (right_) delete right_;
	right_ = NULL;
      }
      return *this;
    }

 SNPTree(std::vector<SNP>& snp_vals,
	 unsigned int depth = 16,
	 unsigned int minbucket = 64,
	 int leftextent  = 0,
	 int rightextent = 0,
	 unsigned int maxbucket = 512) {
   left_ = right_ = NULL;
   
   --depth;
   SNPSorter snp_sorter;
   if (depth == 0 || (snp_vals.size() < minbucket && snp_vals.size() < maxbucket)) {
     std::sort(snp_vals.begin(), snp_vals.end(), snp_sorter);
     snps_ = snp_vals;
   } else {
     if (leftextent == 0 && rightextent == 0)
       std::sort(snp_vals.begin(), snp_vals.end(), snp_sorter); // sort SNPs by position
     
     int leftp = 0, rightp = 0, centerp = 0;
     if (leftextent || rightextent) {
       leftp  = leftextent;
       rightp = rightextent;
     } else {
       leftp  = snp_vals.front().pos();
       rightp = snp_vals.back().pos();
     }
     
     centerp = snp_vals.at(snp_vals.size() / 2).pos();
     center_ = centerp;
     
     std::vector<SNP> lefts, rights;
     for (auto snp_iter = snp_vals.begin(); snp_iter != snp_vals.end(); ++snp_iter){
       SNP snp = *snp_iter;
       if (snp.pos() < center_)
	 lefts.push_back(snp);
       else if (snp.pos() > center_)
	 rights.push_back(snp);
       else
	 snps_.push_back(snp);
     }
     
     if (!lefts.empty())
       left_ = new SNPTree(lefts, depth, minbucket, leftp, centerp);
     if (!rights.empty())
       right_ = new SNPTree(rights, depth, minbucket, centerp, rightp);
   }
 }

 void findContained(uint32_t start, uint32_t stop, std::vector<SNP>& overlapping) const {
   if (left_ && start <= center_)
     left_->findContained(start, stop, overlapping);
   
   if (!snps_.empty() && ! (stop < snps_.front().pos())) {
     for (auto snp_iter = snps_.begin(); snp_iter != snps_.end(); ++snp_iter)
       if (snp_iter->pos() >= start && snp_iter->pos() <= stop)
	 overlapping.push_back(*snp_iter);
   }
   
   if (right_ && stop >= center_)
     right_->findContained(start, stop, overlapping);
 }
 
 ~SNPTree(void) {
   // traverse the left and right
   // delete them all the way down
   if (left_)
     delete left_;
   if (right_)
     delete right_;
 }
 
};


bool create_snp_trees(const std::string& chrom, uint32_t start, uint32_t end, uint32_t skip_start, uint32_t skip_stop, vcflib::VariantCallFile& variant_file, HaplotypeTracker* tracker,
                      std::map<std::string, unsigned int>& sample_indices, std::vector<SNPTree*>& snp_trees, std::ostream& logger);

void destroy_snp_trees(std::vector<SNPTree*>& snp_trees);

#endif
