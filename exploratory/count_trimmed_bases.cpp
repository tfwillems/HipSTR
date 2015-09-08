#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>

#include "../bamtools/include/api/BamAlignment.h"
#include "../bamtools/include/api/BamMultiReader.h"

#include "../bgzf_streams.h"
#include "../error.h"
#include "../zalgorithm.h"

struct compare {
  bool operator()(const std::string& first, const std::string& second) {
    return first.size() > second.size();
  }
};

std::string reverse_complement(std::string& seq){
  std::stringstream res;
  for (auto iter = seq.rbegin(); iter != seq.rend(); iter++){
    switch(*iter){
    case 'A': case 'a':
      res << 'T';
      break;
    case 'C': case 'c':
      res << 'G';
      break;
    case 'G': case 'g':
	res << 'C';
      break;
    case 'T': case 't':
      res << 'A';
      break;
    case 'N': case 'n':
      res << 'N';
      break;
    default:
      printErrorAndDie("Invalid base character");
      break;
    }
  }
  return res.str();
}


void reduce_clip_seqs(std::map<std::string, int>& clip_counts, std::vector<std::string>& final_seqs, bool five_prime){
  std::vector<std::string> seqs;
  for (auto iter = clip_counts.begin(); iter != clip_counts.end(); iter++)
    if (iter->second > 1000)
      seqs.push_back(iter->first);
  std::sort(seqs.begin(), seqs.end(), compare());

  for (unsigned int i = 0; i < seqs.size(); i++){
    bool found = false;
    for (unsigned int j = 0; j < final_seqs.size(); j++){
      if ((five_prime && final_seqs[j].substr(final_seqs[j].size()-seqs[i].size()).compare(seqs[i]) == 0) ||
	  (!five_prime && final_seqs[j].substr(0, seqs[i].size()).compare(seqs[i]) == 0)){
	found = true;
	break;
      }
    }
    if (!found)
      final_seqs.push_back(seqs[i]);
  }
}




int trim_five_prime_end(std::vector<std::string>& adapter_seqs, std::string& bases, std::string& qualities){
  int32_t max_coord = -1;
  for (unsigned int i = 0; i < adapter_seqs.size(); i++){
    std::vector<int> num_matches;
    ZAlgorithm::GetSuffixMatchCounts(adapter_seqs[i], bases, num_matches);
    for (unsigned int j = 0; j < num_matches.size(); j++){
      if (num_matches[j] == j+1)
	max_coord = std::max(max_coord, (int32_t)j);
    }
  }
  if (max_coord < 6)
    max_coord = -1;
  int trim     = max_coord + 1;
  bases        = bases.substr(max_coord+1);
  qualities    = qualities.substr(max_coord+1);
  return trim;
}

int trim_three_prime_end(std::vector<std::string>& adapter_seqs, std::string& bases, std::string& qualities){
  int32_t min_coord = bases.size();
  for (unsigned int i = 0; i < adapter_seqs.size(); i++){
    std::vector<int> num_matches;
    ZAlgorithm::GetPrefixMatchCounts(adapter_seqs[i], bases, num_matches);
    for (unsigned int j = 0; j < num_matches.size(); j++){
      if (j+num_matches[j] == bases.size())
	min_coord = std::min(min_coord, (int32_t)j);
    }
  }
  if (min_coord + 7 > bases.size())
    min_coord = bases.size();
  int trim = (bases.size() - min_coord);
  bases        = bases.substr(0, min_coord);
  qualities    = qualities.substr(0, min_coord);
  return trim;
}




int main(int argc, char** argv){
  std::string input_bam_file    = argv[1];
  std::string input_fastq_file  = argv[2];
  std::string output_fastq_file = argv[3];
  BamTools::BamReader reader;
  if (!reader.Open(input_bam_file)) printErrorAndDie("Failed to open BAM file");

  BamTools::BamAlignment alignment;
  int32_t read_count = 0;
  std::map<std::string, int32_t> three_clip_counts, five_clip_counts;

  int32_t fw_five = 0, fw_three = 0, bw_five = 0, bw_three = 0;

  while (reader.GetNextAlignment(alignment)){
    read_count++;
    if (read_count % 1000000 == 0)
      std::cerr << "Processing read #" << read_count << std::endl;

    std::vector<BamTools::CigarOp> cigar_data = alignment.CigarData;
    if (cigar_data.size() <= 1)
      continue;
    else {
      std::string clipped_seq = "";
      if (cigar_data[0].Type == 'H' && cigar_data[1].Type == 'S')
	clipped_seq = alignment.QueryBases.substr(0, cigar_data[1].Length);
      else if (cigar_data[0].Type == 'S')
	clipped_seq = alignment.QueryBases.substr(0, cigar_data[0].Length);
      if (clipped_seq.size() != 0){
	if (alignment.IsReverseStrand()){
	  if (clipped_seq.size() >= 10)
	    three_clip_counts[reverse_complement(clipped_seq)]++;
	  bw_five++;
	}
	else {
	  if (clipped_seq.size() >= 10)
	    five_clip_counts[clipped_seq]++;
	  fw_five++;
	}
      }

      clipped_seq = "";
      if (cigar_data.back().Type == 'H' && cigar_data[cigar_data.size()-2].Type == 'S')
	clipped_seq = alignment.QueryBases.substr(alignment.QueryBases.size()-cigar_data[cigar_data.size()-2].Length);
      else if (cigar_data.back().Type == 'S')
	clipped_seq = alignment.QueryBases.substr(alignment.QueryBases.size()-cigar_data.back().Length);
      if (clipped_seq.size() != 0){
	if (alignment.IsReverseStrand()){
	  bw_three++;
	  if (clipped_seq.size() >= 10)
	    five_clip_counts[reverse_complement(clipped_seq)]++;
	}
	else {
	  if (clipped_seq.size() >= 10)
	    three_clip_counts[clipped_seq]++;
	  fw_three++;
	}
      }

    }
  }
  reader.Close();

  std::vector<std::string> five_final_seqs, three_final_seqs;
  reduce_clip_seqs(five_clip_counts, five_final_seqs, true);
  reduce_clip_seqs(three_clip_counts, three_final_seqs, false);

  std::cout << "5'-seqs:" << " ";
  for (unsigned int i = 0; i < five_final_seqs.size(); i++)
    std::cout << five_final_seqs[i] << " ";
  std::cout << std::endl;

  std::cout << "3'-seqs:" << " ";
  for (unsigned int i = 0; i < three_final_seqs.size(); i++)
    std::cout << three_final_seqs[i] << " ";
  std::cout << std::endl;
  
  // Process FASTQ
  bgzfistream input_fastq(input_fastq_file.c_str());
  bgzfostream output_fastq(output_fastq_file.c_str());

  int32_t clip_count = 0;
  read_count = 0;

  std::string id, bases, sep, qualities;  
  while (std::getline(input_fastq, id)){
    std::getline(input_fastq, bases);
    std::getline(input_fastq, sep);
    std::getline(input_fastq, qualities);

    int clip_five  = trim_five_prime_end(five_final_seqs,  bases, qualities);
    int clip_three = trim_three_prime_end(three_final_seqs, bases, qualities);
    if (clip_five || clip_three)
      clip_count++;
    
    output_fastq << id        << "\n" 
		 << bases     << "\n" 
		 << sep       << "\n" 
		 << qualities << "\n";
  }
  input_fastq.close();
  output_fastq.close();
  std::cout << clip_count << std::endl;
}
