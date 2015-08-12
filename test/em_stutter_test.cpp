#include <algorithm>
#include <assert.h>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "../em_stutter_genotyper.h"
#include "../stringops.h"
#include "../stutter_model.h"

void read_bp_info(std::string input_file, int motif_len,
		  std::vector<std::string>& sample_names, std::vector< std::vector<int> >& num_bps){
  assert(sample_names.size() == 0 && num_bps.size() == 0);
  std::ifstream input(input_file.c_str());
  if (!input.is_open())
    printErrorAndDie("Failed to open stutter data file " + input_file);

  // Each line is of the form SAMPLE_NUM REPEAT_DIFF_STRING COUNT_STRING
  std::string line;
  while (std::getline(input, line)){
    std::istringstream iss(line);
    std::string sample, repeat_diffs, counts;
    std::vector<std::string> rep_diff_lst, count_lst;
    if (!(iss >> sample >> repeat_diffs >> counts))
      printErrorAndDie("Improperly formatted file for EM stutter test");

    sample_names.push_back(sample);
    num_bps.push_back(std::vector<int>());
    split_by_delim(repeat_diffs, ',', rep_diff_lst);
    split_by_delim(counts, ',', count_lst);
    assert(rep_diff_lst.size() == count_lst.size());
    for (unsigned int i = 0; i < rep_diff_lst.size(); i++){
      int bp_diff = atoi(rep_diff_lst[i].c_str())*motif_len;
      for (unsigned int j = 0; j < atoi(count_lst[i].c_str()); j++){
	num_bps.back().push_back(bp_diff);
      }
    }
  }
  input.close();
}

int main(int argc, char* argv[]){
  if (argc != 5)
    printErrorAndDie("EM stutter test requires exactly four arguments: INPUT_FILE HAPLOID MOTIF_LEN FREQ_PHASE_INFO");

  std::string filename = std::string(argv[1]);
  bool haploid         = (std::string(argv[2]).compare("TRUE") == 0) || (std::string(argv[2]).compare("True") == 0);
  int motif_len        = atoi(argv[3]);
  double phasing_freq  = atof(argv[4]);
  std::string chrom    = "chr1";
  std::vector<std::string> sample_names;
  std::vector< std::vector<int> > num_bps;
  std::vector< std::vector<double> > log_p1s, log_p2s; 

  // Read input file data
  read_bp_info(filename, motif_len, sample_names, num_bps);

  if (haploid){
    std::cerr << "Using HAPLOID model" << std::endl;
    // Fill in phasing info with equal phasing values (for compatibility)
    for (unsigned int i = 0; i < num_bps.size(); i++){
      log_p1s.push_back(std::vector<double>());
      log_p2s.push_back(std::vector<double>());
      for (unsigned int j = 0; j < num_bps[i].size(); j++){
	log_p1s.back().push_back(0);
	log_p2s.back().push_back(0);
      }
    }
  }
  else {
    std::cerr << "Using DIPLOID model" << std::endl;
    // Randomly pair haploid individuals to create diploid individuals
    std::vector<int> indices;
    for (unsigned int i = 0; i < (sample_names.size()/2)*2; i++)
      indices.push_back(i);
    std::random_shuffle(indices.begin(), indices.end());

    std::vector<std::string> merged_sample_names;
    std::vector< std::vector<int> > merged_num_bps;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    for (unsigned int i = 0; i < indices.size(); i += 2){
      merged_sample_names.push_back(sample_names[indices[i]] + "_" + sample_names[indices[i+1]]);
      merged_num_bps.push_back(std::vector<int>());
      log_p1s.push_back(std::vector<double>());
      log_p2s.push_back(std::vector<double>());

      for (unsigned int j = 0; j < num_bps[indices[i]].size(); j++){
	merged_num_bps.back().push_back(num_bps[indices[i]][j]);

	// If we have SNP phasing info, assign read to strand #1
	bool have_phasing_info = distribution(generator) < phasing_freq;
	log_p1s.back().push_back(have_phasing_info ?   0.0 : 0.0);
	log_p2s.back().push_back(have_phasing_info ? -10.0 : 0.0);

      }
      for (unsigned int j = 0; j < num_bps[indices[i+1]].size(); j++){
	merged_num_bps.back().push_back(num_bps[indices[i+1]][j]);

	// If we have SNP phasing info, assign read to strand #2
	bool have_phasing_info = distribution(generator) < phasing_freq;
	log_p1s.back().push_back(have_phasing_info ? -10.0 : 0.0);
	log_p2s.back().push_back(have_phasing_info ?   0.0 : 0.0);
      }
    }
    sample_names = merged_sample_names;
    num_bps      = merged_num_bps;
  }

  // Would ideally obtain this from genotyper_bam_processor.h, but there's some annoying namespace collisions between bamtools and vcflib's tabix
  int MAX_EM_ITER         = 100;
  double ABS_LL_CONVERGE  = 0.01;
  double FRAC_LL_CONVERGE = 0.001;
  EMStutterGenotyper genotyper(chrom, 0, 100, haploid, num_bps, log_p1s, log_p2s, sample_names, motif_len, 0);

  if (!genotyper.train(MAX_EM_ITER, ABS_LL_CONVERGE, FRAC_LL_CONVERGE, false, std::cerr)){
    std::cout << "EM_FAILED_TO_CONVERGE" << std::endl;
  }
  else {
    StutterModel* stutter_model = genotyper.get_stutter_model()->copy();
    stutter_model->write_model(chrom, 0, 100, std::cout);
  }
}
