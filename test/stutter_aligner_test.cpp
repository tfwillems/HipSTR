#include <iostream>
#include <math.h>
#include <vector>

#include "../SeqAlignment/StutterAligner.h"

void fill_qual_scores(std::vector<double>& vals, int num_elements, double value){
  vals.clear();
  for (int i = 0; i < num_elements; i++)
    vals.push_back(value);
}

int main(){
  std::string block_seq = "ATATATATAT"; 
  std::string base_seq  = "ATATATATAT";
  std::vector<double> base_log_wrong, base_log_right;
  fill_qual_scores(base_log_wrong, base_seq.size(), log(1e-5));
  fill_qual_scores(base_log_right, base_seq.size(), log(1-1e-5));
  //double prob_a, prob_b, prob_c, prob_d, prob_e;
  /*
  prob_a = align_no_artifact_forward(  block_seq.size(),   block_seq.c_str(), base_seq.size(),   base_seq.c_str(), &base_log_wrong[0], &base_log_right[0]);
  prob_b = align_pcr_insertion_forward(block_seq.size()-4, block_seq.c_str(), base_seq.size(),   base_seq.c_str(), &base_log_wrong[0], &base_log_right[0], 4, 2);
  prob_c = align_pcr_deletion_forward( block_seq.size(),   block_seq.c_str(), base_seq.size()-4, base_seq.c_str(), &base_log_wrong[0], &base_log_right[0], -4);
  prob_d = align_pcr_insertion_reverse(block_seq.size()-4, block_seq.c_str()+block_seq.size()-5, 
				       base_seq.size(), base_seq.c_str()+base_seq.size()-1, &base_log_wrong[base_seq.size()], &base_log_right[base_seq.size()], 4, 2);
  prob_e = align_pcr_deletion_reverse( block_seq.size(),   block_seq.c_str()+block_seq.size()-1, 
				       base_seq.size()-4, base_seq.c_str()+base_seq.size()-5, &base_log_wrong[base_seq.size()-4], &base_log_right[base_seq.size()-4], -4);
  std::cerr << prob_a << " " << prob_b << " " << prob_c << " " << prob_d << " " << prob_e << std::endl;

  prob_a = align_no_artifact_forward(  block_seq.size(),   block_seq.c_str(), base_seq.size()-3, base_seq.c_str(), &base_log_wrong[0], &base_log_right[0]);
  prob_b = align_pcr_insertion_forward(block_seq.size()-4, block_seq.c_str(), base_seq.size()-3, base_seq.c_str(), &base_log_wrong[0], &base_log_right[0], 4, 2);
  prob_c = align_pcr_deletion_forward( block_seq.size(),   block_seq.c_str(), base_seq.size()-7, base_seq.c_str(), &base_log_wrong[0], &base_log_right[0], -4);
  std::cerr << prob_a << " " << prob_b << " " << prob_c << std::endl;

  block_seq = "ATATATATAC";
  prob_a = align_no_artifact_forward(  block_seq.size(),   block_seq.c_str(), base_seq.size(),   base_seq.c_str(), &base_log_wrong[0], &base_log_right[0]);
  prob_b = align_pcr_insertion_forward(block_seq.size()-4, block_seq.c_str(), base_seq.size(),   base_seq.c_str(), &base_log_wrong[0], &base_log_right[0], 4, 2);
  prob_c = align_pcr_deletion_forward( block_seq.size(),   block_seq.c_str(), base_seq.size()-4, base_seq.c_str(), &base_log_wrong[0], &base_log_right[0], -4);
  std::cerr << prob_a << " " << prob_b << " " << prob_c << std::endl;

  prob_a = align_no_artifact_forward(  block_seq.size(),   block_seq.c_str(), base_seq.size()-3, base_seq.c_str(), &base_log_wrong[0], &base_log_right[0]);
  prob_b = align_pcr_insertion_forward(block_seq.size()-4, block_seq.c_str(), base_seq.size()-3, base_seq.c_str(), &base_log_wrong[0], &base_log_right[0], 4, 2);
  prob_c = align_pcr_deletion_forward( block_seq.size(),   block_seq.c_str(), base_seq.size()-7, base_seq.c_str(), &base_log_wrong[0], &base_log_right[0], -4);
  std::cerr << prob_a << " " << prob_b << " " << prob_c << std::endl;

  base_seq  = "ATATACATATAT";
  block_seq = "ATATATATAT";
  prob_a = align_no_artifact_forward(  block_seq.size(),   block_seq.c_str(), base_seq.size()-2, base_seq.c_str(), &base_log_wrong[0], &base_log_right[0]);
  prob_b = align_pcr_insertion_forward(block_seq.size(),   block_seq.c_str(), base_seq.size(),   base_seq.c_str(), &base_log_wrong[0], &base_log_right[0], 2, 2);
  prob_c = align_pcr_deletion_forward( block_seq.size(),   block_seq.c_str(), base_seq.size()-4, base_seq.c_str(), &base_log_wrong[0], &base_log_right[0], -2);
  std::cerr << prob_a << " " << prob_b << " " << prob_c << std::endl;
  

  base_seq  = "ATATACACACATAT";
  block_seq = "ATATACATAT";
  fill_qual_scores(base_log_wrong, base_seq.size(), log(1e-5));
  fill_qual_scores(base_log_right, base_seq.size(), log(1-1e-5));
  prob_b = align_pcr_insertion_forward(block_seq.size(), block_seq.c_str(), base_seq.size(), base_seq.c_str(), &base_log_wrong[0], &base_log_right[0], 4, 2);
  prob_d = align_pcr_insertion_reverse(block_seq.size(), block_seq.c_str()+block_seq.size()-1, base_seq.size(), base_seq.c_str()+base_seq.size()-1, &base_log_wrong[base_seq.size()-1], &base_log_right[base_seq.size()-1], 4, 2);
  std::cerr << prob_b << " " << prob_d << std::endl;*/
}
