#ifndef ALIGNMENT_MODEL_H_
#define ALIGNMENT_MODEL_H_
  
const int MAX_HOMOP_LEN       = 15;
const int MAX_SEQ_DEL         = 5;
const double LOG_INS_TO_INS   = -1.0; // log(e^-1)
const double LOG_INS_TO_MATCH = -0.4586751453870818910216436; // log(1-e^-1)

// Obtained by taking the log_e of the values in the insertion probabilities in the Dindel source code
extern double LOG_MATCH_TO_INS[MAX_HOMOP_LEN+1];

// Incorporates both MATCH->MATCH transition probability and P(d)
// P(no deletion) = P(d=1) = 1-2*P(M->I) such that total probability of deletions is equal to P(M->I)
// Scales remaining deletion sizes in same way as insertions (e^-(d-1)) so that the model is entirely symmetrical
// w.r.t. insertion and deletion transition probabilities
extern double LOG_DEL_N[MAX_HOMOP_LEN+1][MAX_SEQ_DEL+1];


void init_probabilities();
void print_parameters();

#endif
