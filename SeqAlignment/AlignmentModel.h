#ifndef ALIGNMENT_MODEL_H_
#define ALIGNMENT_MODEL_H_
  
const unsigned int MAX_HOMOP_LEN = 15;
const double LOG_INS_TO_INS      = -1.0; // log(e^-1)
const double LOG_INS_TO_MATCH    = -0.4586751453870818910216436; // log(1-e^-1)
const double LOG_DEL_TO_DEL      = -1.0; // log(e^-1)
const double LOG_DEL_TO_MATCH    = -0.4586751453870818910216436; // log(1-e^-1)

extern double LOG_MATCH_TO_MATCH[MAX_HOMOP_LEN+1];
extern double LOG_MATCH_TO_INS[MAX_HOMOP_LEN+1];
extern double LOG_MATCH_TO_DEL[MAX_HOMOP_LEN+1];

void init_alignment_model();
void print_alignment_model();

#endif
