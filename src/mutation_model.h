#ifndef MUTATION_MODEL_H_
#define MUTATION_MODEL_H_

#include <math.h>

#include "vcf_reader.h"

class MutationModel {
  double log_mut_prior_;

 public:
  MutationModel(VCF::Variant& str_variant){
    assert(str_variant.num_alleles() > 1);
    
    // The allele on each haplotype can mutate to N-1 alleles, so assuming a 
    // uniform prior each mutation has a prior of 1/(2*(N-1))
    log_mut_prior_ = -log10(2) - log10(str_variant.num_alleles()-1);
  }

  /*
   * Log10-likelihood of mutating from the parental to the child allele,
   * given that a mutation occurred
   */
  double log_prior_mutation(int parental_allele, int child_allele){
    return log_mut_prior_;
  }
};

#endif
