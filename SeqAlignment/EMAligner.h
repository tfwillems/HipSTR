#ifndef EM_ALIGNER_H_
#define EM_ALIGNER_H_

#include <vector>

#include "error.h"
#include "AlignmentData.h"
#include "Haplotype.h"

class PathTracker {
 private:
  const static int MAX_NUM_BITS = 32; // Based on size of long
  std::vector<int> bit_indices;
  
  inline int num_bits(int val){
    int n = 0;
    while (val >>= 1)
      n++;
    return n;
  }

  void calc_bit_indices(Haplotype* haplotype){
    int nbits = 0;
    for (unsigned int block_index = 0; block_index < haplotype->num_blocks(); block_index++){
      int block_size = haplotype->num_options(block_index);
      bit_indices.push_back(nbits);
      nbits += num_bits(block_size); 
    }
    if (nbits > MAX_NUM_BITS)
      printErrorAndDie("Too many bits required to track haplotype path");
  }

 public:
  PathTracker(Haplotype* haplotype){
    calc_bit_indices(haplotype);
  }

  inline long add_track_to_path(int block_index, int block_option, long path){ 
    return (path || (block_option << bit_indices[block_index]));
  }
};


class EMAligner {
 private:
  Haplotype* haplotype_;
 public:
  EMAligner(Haplotype* haplotype){
    haplotype_ = haplotype;
  }

  void run_EM(std::vector<Alignment>& all_alignments);
  
  void compute_seed_bases(std::vector<Alignment>::iterator begin_iter,
			  std::vector<Alignment>::iterator end_iter,
			  std::vector<int>& seed_bases);

  void compute_forward_probs(std::vector<Alignment>::iterator begin_iter, 
			     std::vector<Alignment>::iterator end_iter,
			     std::vector<int>& seed_bases);

  void compute_backward_probs(std::vector<Alignment>::iterator begin_iter,
			      std::vector<Alignment>::iterator end_iter,
			      std::vector<int>& seed_bases);

  void compute_posterior_probs();
};


#endif
