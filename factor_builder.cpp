#include <algorithm>

#include "factor_builder.h"
#include "snp_phasing_quality.h"

inline double log_sum_exp(double v1, double v2){
  return (v1 > v2 ? v1 + log(1 + exp(v2-v1)) : v2 + log(1 + exp(v1-v2)));
}

void add_GT_factors(vcf::Variant& variant, std::map<std::string, int>& sample_indices,
		    int min_bp_length, int max_bp_length, bool ignore_phase,
		    std::vector<dai::Var>& node_1s, std::vector<dai::Var>& node_2s,
		    std::vector<dai::Factor>& factors, std::vector< std::pair<dai::Var,int> >& clamps){
  std::string motif_key = "MOTIF";
  int motif_len         = variant.getInfoValueString(motif_key).size();
  int num_alleles       = (max_bp_length - min_bp_length)/motif_len + 1;
  int num_out_frame     = 0;

  for (auto sample_iter = variant.sampleNames.begin(); sample_iter != variant.sampleNames.end(); ++sample_iter){
    assert(sample_indices.find(*sample_iter) != sample_indices.end());
    int sample_index = sample_indices[*sample_iter];
    std::string gts  = variant.getGenotype(*sample_iter);
    size_t pos;
    if ((pos=gts.find("|") != std::string::npos) || (pos=gts.find("/") != std::string::npos)){
      std::string& allele_1 = variant.alleles.at(atoi(gts.substr(0, pos).c_str()));
      std::string& allele_2 = variant.alleles.at(atoi(gts.substr(pos+1).c_str()));
      assert(allele_1.size() >= min_bp_length && allele_1.size() <= max_bp_length);
      assert(allele_2.size() >= min_bp_length && allele_2.size() <= max_bp_length);
      
      bool out_of_frame = ((allele_1.size()-min_bp_length) % motif_len != 0);
      out_of_frame     |= ((allele_2.size()-min_bp_length) % motif_len != 0);
      if (out_of_frame){
	num_out_frame++;
	continue;
      }
      
      int nrepeats_one = (allele_1.size()-min_bp_length)/motif_len;
      int nrepeats_two = (allele_2.size()-min_bp_length)/motif_len;
      if (ignore_phase || (gts.find("|") == std::string::npos)){
	factors.push_back(dai::Factor(dai::VarSet(node_1s[sample_index], node_2s[sample_index])));
	int index_one = nrepeats_one*num_alleles + nrepeats_two;
	int index_two = nrepeats_two*num_alleles + nrepeats_one;
	if (index_one != index_two){
	  factors.back().set(index_one, 0.5);
	  factors.back().set(index_two, 0.5);
	}
	else
	  factors.back().set(index_one, 1.0);
      }
      else {
	clamps.push_back(std::pair<dai::Var,int>(node_1s.at(sample_index), nrepeats_one));
	clamps.push_back(std::pair<dai::Var,int>(node_2s.at(sample_index), nrepeats_two));
      }
    }
  }
  
  if (num_out_frame != 0)
    std::cerr << "WARNING: Skipped " << num_out_frame << " samples with out-of-frame STR alleles when adding GT factors" << std::endl;
}

void add_GL_factors(vcf::Variant& variant, std::map<std::string, int>& sample_indices,
		    int min_allele, int max_allele, std::vector<dai::Var>& node_1s, std::vector<dai::Var>& node_2s,
		    std::vector<dai::Factor>& factors){
  /*
  std::vector<bool> in_frame;
  std::vector<int>  allele_indices;
  for (unsigned int i = 0; i < variant.alleles.size(); ++i){
    in_frame.append((variant.alleles[i].size()-min_bp_length)%motif_len == 0);
    if (in_frame.back())
      allele_indices.append((variant.alleles[i].size()-min_bp_length)/motif_len);
    else
      allele_indices.append(-1);
  }

  for (auto sample_iter = variant.sampleNames.begin(); sample_iter != variant.sampleNames.end(); ++sample_iter){
    assert(sample_indices.find(*sample_iter) != sample_indices.end());
    int sample_index = sample_indices[*sample_iter];
    std::string gl_string = variant.getSampleValueString("GL", *sample_iter);

    std::vector<double> potentials;
    std::vector<int> potential_indices;
    int gl_index = 0;
    for (int allele_1 = 0; allele_1 < variant.alleles.size(); ++allele_1){
      if (!in_frame[allele_1]){
	gl_index += allele_1 + 1;
	continue;
      }
	  
      for (int allele_2 = 0; allele_2 <= allele_1; ++allele_2){
	if(in_frame[allele_2]){
	  potentials.append();
	  potential_indices.append();
	}
	gl_index++;
      }
    }
  }
  */
}

void add_read_factors(std::vector<int>& nrepeats, 
		      std::vector<double>& snp_log_p1s, std::vector<double>& snp_log_p2s,
		      int min_allele, int max_allele, StutterModel& model,
		      dai::Var& node_id1, dai::Var& node_id2,
		      std::vector<dai::Factor>& factors){
  assert(nrepeats.size() == snp_log_p1s.size() && nrepeats.size() == snp_log_p2s.size());
      
  // TO DO: Incorporate dropout probabilites by replacing LOG_ONE_HALF with unequal values
  std::vector<double> log_vals;
  for (int gt_1 = min_allele; gt_1 <= max_allele; ++gt_1){
    for (int gt_2 = min_allele; gt_2 <= max_allele; ++gt_2){
      double log_val = 0.0;
      for (unsigned int k = 0; k < nrepeats.size(); ++k)
	log_val += log_sum_exp(LOG_ONE_HALF + snp_log_p1s[k] + model.log_stutter_pmf(gt_1, nrepeats[k]), 
			       LOG_ONE_HALF + snp_log_p2s[k] + model.log_stutter_pmf(gt_2, nrepeats[k]));
      log_vals.push_back(log_val);
    }
  }
  
  // Construct new factor and record probability table
  // Value for node_2 changes fastest in table
  factors.push_back(dai::Factor(dai::VarSet(node_id1, node_id2)));
  double max_log_val = *std::max_element(log_vals.begin(), log_vals.end());
  for (unsigned int count = 0; count < log_vals.size(); ++count)
    factors.back().set(count, exp(log_vals[count]-max_log_val));
}


void add_read_factors(AlnVector& paired_str_reads, AlnVector& mate_reads, AlnVector& unpaired_str_reads,
		      BaseQuality& base_qualities, SNPTree* snp_tree, 
		      int min_allele, int max_allele, int motif_len, StutterModel& model,
		      dai::Var node_id1, dai::Var node_id2,
		      std::vector<dai::Factor>& factors){
  std::vector<int> nrepeats;
  std::vector<bool> success;
  std::vector<double> snp_log_p1s, snp_log_p2s;
  // Extract the number of repeats associated with each STR-spanning read
  /*
  int fail_count = 0, out_frame_count = 0;
  for (unsigned int i = 0; i < paired_str_reads.size(); ++i){
    int bp_diff;
    bool got_size = ExtractCigar(paired_str_reads[i].CigarData, paired_str_reads[i].Position, region.start(), region.stop(), bp_diff);
    if (got_size){
      int bp_size;
      assert(bp_size >= min_allele && bp_size <= max_allele);
      if ((bp_size-min_allele)%motif_len == 0){
	success.push_back(true);
	nrepeats.push_back((bp_size-min_allele)/motif_len);
      }
      else {
	success.push_back(false);
	nrepeats.push_back(-1);
	out_frame_count++;
      }
    }
    else{
      success.push_back(false);
      nrepeats.push_back(-1);
      fail_count++;
    }
  }
  */

  // Extract the phasing qualities determined from overlapping heterozygous SNPs
  //calc_het_snp_factors(paired_str_reads, mate_reads, base_qualities, snp_tree, snp_log_p1s, snp_log_p2s);
  //calc_het_snp_factors(unpaired_str_reads, base_qualities, snp_tree, snp_log_p1s, snp_log_p2s);

  // TO DO: Filter vectors to only include successful reads

  // Add appropriate factors
  add_read_factors(nrepeats, snp_log_p1s, snp_log_p2s, min_allele, max_allele, model, node_id1, node_id2, factors);
}






void add_sample_factors(std::vector<std::string>& samples, int num_alleles, 
			std::vector<dai::Var>& node_1s, std::vector<dai::Var>& node_2s, std::map<std::string, int>& sample_indices){
  assert(node_1s.size() == 0 && node_2s.size() == 0 && sample_indices.size() == 0);
  int node_index = 0;
  for (unsigned int i = 0; i < samples.size(); ++i){
    sample_indices[samples[i]] = i;
    node_1s.push_back(dai::Var(node_index++, num_alleles));
    node_2s.push_back(dai::Var(node_index++, num_alleles));
  }
}




/*
// Outline
// Create variables for samples
// i)   Add GTs if provided
// ii)  Add GLs if provided
// iii) Add reads if provided
// Enusre that only 1 of i,ii,iii was added for each sample


// Input
vcf::Variant variant;
int min_bp_length, int max_bp_length, bool ignore_phase;
int min_allele, int max_allele, motif_len;
std::vector<AlnVector> paired_str_reads, mate_reads, unpaired_str_reads;
std::vector<SNPTree*> snp_trees;
StutterModel model;



// Output
std::vector<dai::Var> node_1s, node_2s;
std::map<std::string, int>& sample_indices;
std::vector<dai::Factor>& factors;
std::vector< std::pair<dai::Var,int> >& clamps;


int num_alleles = max_allele - min_allele + 1;

add_sample_factors(samples, num_alleles, node_1s, node_2s, sample_indices);

if (variant != NULL)
  add_GT_factors(*variant, sample_indices, min_bp_length, max_bp_length, ignore_phase, node_1s, node_2s, factors, clamps);

if (variant != NULL)
  add_GL_factors(*variant, sample_indices, min_allele, max_allele, node_1s, node_2s, factors);

for (unsigned int i = 0; i < sample.size(); ++i){
  int sample_index = sample_indices[samples[i]];
  add_read_factors(paired_str_reads[i], mate_reads[i], unpaired_str_reads, base_qualities, snp_trees[i],
		   min_allele, max_allele, motif_len, model, node_1s[sample_index], node_2s[sample_index],
		   factors);
 }
*/
