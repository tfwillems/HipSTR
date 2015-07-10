#include "math.h"

#include <iomanip>
#include <iostream>

#include "AlignmentModel.h"

// Obtained from Dindel source code v 1.01
double dindel_probs[10] = {2.9e-5, 2.9e-5, 2.9e-5, 2.9e-5, 4.3e-5, 1.1e-4, 2.4e-4, 5.7e-4, 1.0e-3, 1.4e-3};

double LOG_MATCH_TO_MATCH[MAX_HOMOP_LEN+1];
double LOG_MATCH_TO_INS[MAX_HOMOP_LEN+1];
double LOG_MATCH_TO_DEL[MAX_HOMOP_LEN+1];

/*
  Initialize each transition array
  LOG_MATCH_TO_INS and LOG_MATCH_TO_DEL values are obtained by logging the values utilized in Dindel
  LOG_MATCH_TO_MATCH is equal to log(1-P(M->I)-P(M->D))
 */
void init_alignment_model(){
  // Values won't be used anyways as homopolymer length >= 1
  LOG_MATCH_TO_INS[0]   = 0;
  LOG_MATCH_TO_DEL[0]   = 0;
  LOG_MATCH_TO_MATCH[0] = 0;

  // Set values for each homopolymer length
  for (unsigned int i = 1; i <= MAX_HOMOP_LEN; i++){
    LOG_MATCH_TO_INS[i]   = (i <= 10 ? log(dindel_probs[i-1]) : log(dindel_probs[9]+(4.3e-4)*(i-10)));
    LOG_MATCH_TO_DEL[i]   = (i <= 10 ? log(dindel_probs[i-1]) : log(dindel_probs[9]+(4.3e-4)*(i-10)));
    LOG_MATCH_TO_MATCH[i] = log(1.0 - exp(LOG_MATCH_TO_INS[i]) - exp(LOG_MATCH_TO_DEL[i]));
  }
}


void print_alignment_model(){
  std::cerr << "Match->Insertion transition probabilities:\n"
	    << " homopolymer_len: log_P \n";
  for (unsigned int i = 0; i <= MAX_HOMOP_LEN; i++)
    std::cerr <<  " " << i << ":" << LOG_MATCH_TO_INS[i] << "\n";
  std::cerr << "\n";

  std::cerr << "Match->Deletion transition probabilities:\n"
	    << " homopolymer_len: log_P \n";
  for (unsigned int i = 0; i <= MAX_HOMOP_LEN; i++)
    std::cerr << " " << i << ":" << LOG_MATCH_TO_DEL[i] << "\n";
  std::cerr << std::endl;
}
