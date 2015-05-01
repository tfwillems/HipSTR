#include "math.h"

#include <iomanip>
#include <iostream>

#include "AlignmentModel.h"

// Obtained from Dindel source code v 1.01
double dindel_probs[10] = {2.9e-5, 2.9e-5, 2.9e-5, 2.9e-5, 4.3e-5, 1.1e-4, 2.4e-4, 5.7e-4, 1.0e-3, 1.4e-3};

double LOG_MATCH_TO_INS[MAX_HOMOP_LEN+1];
double LOG_DEL_N[MAX_HOMOP_LEN+1][MAX_SEQ_DEL+1];

/*
  Initialize LOG_MATCH_TO_INS and LOG_DEL_N arrays. 
  LOG_MATCH_TO_INS values are obtained by logging the values utilized in Dindel.
  LOG_DEL_N values are obtained by scaling P(d=1) to be equal to 1-2P(M->I). 
  The remaining deletions then have total probability P(M->I) where P(d) is proportional to e^(-(d-1)). 
  This results in a symmetric treatment of insertions and deletions w.r.t. transtion probabilities
 */
void init_alignment_model(){
  LOG_MATCH_TO_INS[0] = 0;
  for (int j = 0; j <= MAX_SEQ_DEL; j++)
    LOG_DEL_N[0][j] = 0.0;

  for (unsigned int i = 1; i <= MAX_HOMOP_LEN; i++){
    LOG_MATCH_TO_INS[i] = (i <= 10 ? log(dindel_probs[i-1]) : log(dindel_probs[9]+(4.3e-4)*(i-10)));
    double ins_prob = exp(LOG_MATCH_TO_INS[i]);
    LOG_DEL_N[i][0] = 0.0;
    LOG_DEL_N[i][1] = log(1-2*ins_prob);
    double del_prob = ins_prob;
    double tot = 0.0;
    for (int j = 2; j<= MAX_SEQ_DEL; j++)
      tot += exp(-(j-1));
    double log_norm_factor = log(tot/del_prob); 
    for (int j = 2; j <= MAX_SEQ_DEL; j++)
      LOG_DEL_N[i][j] = (-(j-1) - log_norm_factor);
  }
}


void print_alignment_model(){
  std::cout << "Match->Insertion transition probabilities:\n" 
	    << " homopolymer_len: log_P \n";
  for (unsigned int i = 0; i <= MAX_HOMOP_LEN; i++)
    std::cout <<  " " << i << ":" << LOG_MATCH_TO_INS[i] << "\n";
  std::cout << "\n";

  std::cout << "Deletion transition probabilities:\n"
	    << " homopolymer_len deletion_size: log_P\n";
  for (unsigned int i = 0; i <= MAX_HOMOP_LEN; i++){
    for (int j = 0; j <= MAX_SEQ_DEL; j++){
      std::cout << " " << i << " " << j << ": " << std::setw(5) << LOG_DEL_N[i][j] << "\n";
    }
    std::cout << std::endl;
  }
}
