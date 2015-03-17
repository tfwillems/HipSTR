#include "snp_phasing_quality.h"
#include "error.h"

void printCigarString(BamTools::BamAlignment& aln){
  for (unsigned int i = 0; i < aln.CigarData.size(); ++i)
    std::cerr << aln.CigarData[i].Length << aln.CigarData[i].Type;
}

void extract_bases_and_qualities(BamTools::BamAlignment& aln, std::vector<SNP>& snps, 
				 std::vector<char>& bases,  std::vector<char>& quals){
  assert(bases.size() == 0 && quals.size() == 0);
  assert(aln.CigarData.size() > 0);
  assert(std::is_sorted(snps.begin(), snps.end(), SNPSorter()));

  int32_t pos = aln.Position;
  unsigned int snp_index = 0, cigar_index = 0, base_index = 0;
  while (snp_index < snps.size() && cigar_index < aln.CigarData.size()){
    switch(aln.CigarData[cigar_index].Type){
    case 'M': case '=': case 'X': 
      if (snps[snp_index].pos() < pos + aln.CigarData[cigar_index].Length){
	bases.push_back(aln.QueryBases.at(snps[snp_index].pos() - pos + base_index));
	quals.push_back(aln.Qualities.at(snps[snp_index].pos() - pos + base_index));
	snp_index++;
      }
      else {
	pos += aln.CigarData[cigar_index].Length;
	base_index += aln.CigarData[cigar_index].Length;
	cigar_index++;
      }
      break;
    case 'D':
      if (snps[snp_index].pos() < pos + aln.CigarData[cigar_index].Length){
	bases.push_back('-');
	quals.push_back('-');
	snp_index++;
      }
      else {
	pos += aln.CigarData[cigar_index].Length;
	cigar_index++;
      }
      break;
    case 'I':
      base_index += aln.CigarData[cigar_index].Length;
      cigar_index++;
      break;
    case 'S': 
      if (snps[snp_index].pos() < pos){
	bases.push_back('-'); // Ignore bases in soft clips by marking them as dummy deletions
	quals.push_back('-');
	snp_index++;
      }
      else {
	base_index += aln.CigarData[cigar_index].Length;
	cigar_index++;
      }
      break;
    case 'H':
      cigar_index++;
      break;
    default:
      printErrorAndDie("Invalid CIGAR option encountered");
      break;
    }
  }
  assert(bases.size() == snps.size() && snp_index == snps.size());
}

void add_log_phasing_probs(BamTools::BamAlignment& aln, SNPTree* tree, BaseQuality& base_qualities, 
			   double& log_p1, double& log_p2, int& match_count, int& mismatch_count){
  std::vector<SNP> snps;  
  // NOTE: GetEndPosition() returns a non-inclusive position. Use -1 to only find SNPs overlapped by read
  tree->findContained(aln.Position, aln.GetEndPosition()-1, snps);
  if (snps.size() != 0){
    std::vector<char> bases, quals;
  
    extract_bases_and_qualities(aln, snps, bases, quals);
    assert(snps.size() == bases.size());
    for (unsigned int i = 0; i < snps.size(); ++i){
      if (bases[i] != '-'){
	  log_p1 += (bases[i] == snps[i].base_one() ? base_qualities.log_prob_correct(quals[i]) : base_qualities.log_prob_error(quals[i]));
	  log_p2 += (bases[i] == snps[i].base_two() ? base_qualities.log_prob_correct(quals[i]) : base_qualities.log_prob_error(quals[i]));
	  
	  if (bases[i] != snps[i].base_one() && bases[i] != snps[i].base_two())
	    mismatch_count++;
	  else
	    match_count++;
      }
    }
  }
}


void calc_het_snp_factors(std::vector<BamTools::BamAlignment>& str_reads, std::vector<BamTools::BamAlignment>& mate_reads, 
			  BaseQuality& base_qualities, SNPTree* snp_tree,
			  std::vector<double>& log_p1s, std::vector<double>& log_p2s, int& match_count, int& mismatch_count) {
  assert(str_reads.size() == mate_reads.size());
  log_p1s.resize(str_reads.size(), 0.0);
  log_p2s.resize(str_reads.size(), 0.0);
  for (unsigned int i = 0; i < str_reads.size(); i++){
    double log_p1 = 0.0, log_p2 = 0.0;
    add_log_phasing_probs(str_reads[i],  snp_tree, base_qualities, log_p1, log_p2, match_count, mismatch_count);
    add_log_phasing_probs(mate_reads[i], snp_tree, base_qualities, log_p1, log_p2, match_count, mismatch_count);
    log_p1s.push_back(log_p1);
    log_p2s.push_back(log_p2);
  }
}

void calc_het_snp_factors(std::vector<BamTools::BamAlignment>& str_reads, BaseQuality& base_qualities, SNPTree* snp_tree, 
			  std::vector<double>& log_p1s, std::vector<double>& log_p2s, int& match_count, int& mismatch_count){
  for (unsigned int i = 0; i < str_reads.size(); ++i){
    double log_p1 = 0.0, log_p2 = 0.0;
    add_log_phasing_probs(str_reads[i], snp_tree, base_qualities, log_p1, log_p2, match_count, mismatch_count);
    log_p1s.push_back(log_p1);
    log_p2s.push_back(log_p2);
  }
}
