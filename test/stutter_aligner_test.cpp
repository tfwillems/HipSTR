#include <iostream>
#include <math.h>

#include "../SeqAlignment/StutterAligner.h"

void fill_qual_scores(double* arr, int num_elements, double value){
  for (int i = 0; i < num_elements; i++)
    arr[i] = value;
}

int main(){
  std::string block_seq = "ATATATATAT"; 
  std::string base_seq  = "ATATATATAT";
  
  double base_log_wrong[block_seq.size()]; 
  double base_log_right[base_seq.size()];
  fill_qual_scores(base_log_wrong, base_seq.size(), log(1e-5));
  fill_qual_scores(base_log_right, base_seq.size(), log(1-1e-5));

  double prob_a, prob_b, prob_c;
  prob_a = align_no_artifact(  block_seq.size(),   block_seq, base_seq.size(),   base_seq.c_str(), base_log_wrong, base_log_right);
  prob_b = align_pcr_insertion(block_seq.size()-4, block_seq, base_seq.size(),   base_seq.c_str(), base_log_wrong, base_log_right, 4);
  prob_c = align_pcr_deletion( block_seq.size(),   block_seq, base_seq.size()-4, base_seq.c_str(), base_log_wrong, base_log_right, -4);
  std::cerr << prob_a << " " << prob_b << " " << prob_c << std::endl;

  prob_a = align_no_artifact(  block_seq.size(),   block_seq, base_seq.size()-3, base_seq.c_str(), base_log_wrong, base_log_right);
  prob_b = align_pcr_insertion(block_seq.size()-4, block_seq, base_seq.size()-3, base_seq.c_str(), base_log_wrong, base_log_right, 4);
  prob_c = align_pcr_deletion( block_seq.size(),   block_seq, base_seq.size()-7, base_seq.c_str(), base_log_wrong, base_log_right, -4);
  std::cerr << prob_a << " " << prob_b << " " << prob_c << std::endl;

  block_seq = "ATATATATAC";
  prob_a = align_no_artifact(  block_seq.size(),   block_seq, base_seq.size(),   base_seq.c_str(), base_log_wrong, base_log_right);
  prob_b = align_pcr_insertion(block_seq.size()-4, block_seq, base_seq.size(),   base_seq.c_str(), base_log_wrong, base_log_right, 4);
  prob_c = align_pcr_deletion( block_seq.size(),   block_seq, base_seq.size()-4, base_seq.c_str(), base_log_wrong, base_log_right, -4);
  std::cerr << prob_a << " " << prob_b << " " << prob_c << std::endl;

  prob_a = align_no_artifact(  block_seq.size(),   block_seq, base_seq.size()-3, base_seq.c_str(), base_log_wrong, base_log_right);
  prob_b = align_pcr_insertion(block_seq.size()-4, block_seq, base_seq.size()-3, base_seq.c_str(), base_log_wrong, base_log_right, 4);
  prob_c = align_pcr_deletion( block_seq.size(),   block_seq, base_seq.size()-7, base_seq.c_str(), base_log_wrong, base_log_right, -4);
  std::cerr << prob_a << " " << prob_b << " " << prob_c << std::endl;  
}
