#include "AlignmentOps.h"
#include "EMAligner.h"

// Mode 1: Stutter model fixed
// Want to estimate haplotype frequencies and underlying genotypes
// 1. Calculate P(Read/Hap) for all samples and all reads using Viterbi/Sum-Product algorithm?
// 2. Run EM until convergence,  without recalculating P(Read/Hap)
// 3. Marginalize across haplotypes to determine total STR probabilities 

// Mode 2: Stutter model unknown
// Want to estimate stutter probabilities
// 1. For each sample:
//      i.   Calculate forward probs for all reads
//      ii.  One read at a time, calculate backward probs
//      iii. Use the backwards probs to determine the posteriors and
//           relevant statistics for stutter estimation
// 2. Pool statistics across all samples to reestimate stutter probabilities
//  and haplotype frequencies
// 3. Repeat 1-2 until convergence


void EMAligner::run_EM(std::vector<Alignment>& all_alignments){
  // Sort alignments by sample
  sortAlignments(all_alignments);

  bool converged = false;
  std::vector<int> seed_bases;
  while (!converged){
    auto begin_iter = all_alignments.begin();
    auto end_iter   = all_alignments.begin();
    while (true){
      // Compute per-sample posterior probabilities using forward-backward algorithm
      if (end_iter == all_alignments.end() || (begin_iter->get_sample().compare(end_iter->get_sample()) != 0)){
	seed_bases.clear();
	compute_seed_bases(begin_iter, end_iter, seed_bases);
	compute_backward_probs(begin_iter, end_iter, seed_bases);
	compute_forward_probs(begin_iter, end_iter, seed_bases);
	compute_posterior_probs();

	if (end_iter == all_alignments.end())
	  break;
      }
      end_iter++;
    }

    // TO DO: Use the posterior probabilities across all samples to reestimate
    // 1. Haplotype   frequencies
    // 2. PCR Stutter probabilities


    // TO DO: Test for convergence
  }
}


void EMAligner::compute_seed_bases(std::vector<Alignment>::iterator begin_iter,
				   std::vector<Alignment>::iterator end_iter,
				   std::vector<int>& seed_bases){
  while (begin_iter != end_iter){
    


    begin_iter++;
  }
}


void EMAligner::compute_forward_probs(std::vector<Alignment>::iterator begin_iter, 
				      std::vector<Alignment>::iterator end_iter,
				      std::vector<int>& seed_bases){
  int cur_align_index=0, num_skip=0, num_proc=0;
  while (begin_iter != end_iter){
    if (seed_bases[cur_align_index] == -1)
      num_skip++;
    else{
      num_proc++;
      


      // Process region to left of seed
      if (seed_bases[cur_align_index] != 0){

      }

      // Process region to right of seed
      if (seed_bases[cur_align_index] != begin_iter->get_base_qualities().size()-1){

      }
    }
    begin_iter++;
    cur_align_index++;
  }
}

void EMAligner::compute_backward_probs(std::vector<Alignment>::iterator begin_iter,
				       std::vector<Alignment>::iterator end_iter,
				       std::vector<int>& seed_bases){
  while (begin_iter != end_iter){

    for (unsigned int hap_index = haplotype_->num_blocks()-1; hap_index >= 0; hap_index--){
      
    }
    begin_iter++;
  }
}

void EMAligner::compute_posterior_probs(){

}
