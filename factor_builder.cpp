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

  /*
  for (auto sample_iter = variant.sampleNames.begin(); sample_iter != variant.sampleNames.end(); ++sample_iter){
    assert(sample_indices.find(*sample_iter) != sample_indices.end());
    int sample_index = sample_indices[*sample_iter];
    std::string gts  = variant.getGenotype(*sample_iter);
    size_t pos;
    if ((pos=gts.find("|") != std::string::npos) || (pos=gts.find("/") != std::stringnpos)){
      std::string allele_1 = variant.alleles.at(atoi(gts.substr(0, pos).c_str()));
      std::string allele_2 = variant.alleles.at(atoi(gts.substr(pos+1).c_str()));

      assert(allele_1.size() >= min_bp_length && allele_1.size() <= max_bp_length);
      assert(allele_2.size() >= min_bp_length && allele_2.size() <= max_bp_length);
      
      bool out_of_frame = ((allele_1.size()-min_bp_length) % motif_len != 0);
      out_of_frame     |= ((allele_2.size()-min_bp_length) % motif_len != 0);
      if (out_out_frame){
	num_out_frame++;
	continue;
      }
      
      int nrepeats_one = (allele_1.size() - min_bp_length)/motif_len;
      int nrepeats_two = (allele_2.size() - min_bp_length)/motif_len;
      if (ignore_phase || (gts.find("|") == std::string::npos)){
	factors.push_back(dai::Factor(dai::VarSet(node_1s[i], node_2s[i])));
	int index_one = nrepeats_one*num_alleles + nrepeats_two;
	int index_two = nrepeats_two*num_alleles + nrepeats_one;
	if (index_one != index_two){
	  factors.back().set(index_one, 0.5);
	  factors.back().set(index_two, 0.5);
	}
	else
	  factors.back.set(index_one, 1.0);
      }
      else {
	clamps.push_back(std::pair<dai::Var,int>(node_1s.at(sample_index), nrepeats_one));
	clamps.push_back(std::pair<dai::Var,int>(node_2s.at(sample_index), nrepeats_two));
      }
    }
  }
  */
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

void add_read_factors(std::vector<int>& nrepeats, std::vector<double>& snp_log_p1s, std::vector<double>& snp_log_p2s,
		      int min_allele, int max_allele, StutterModel& model,
		      dai::Var& node_id1, dai::Var& node_id2,
		      std::vector<dai::Factor>& factors){
  assert(nrepeats.size() == snp_log_p1s.size() && nrepeats.size() == snp_log_p2s.size());
  factors.push_back(dai::Factor(dai::VarSet(node_id1, node_id2)));

  const double LOG_ONE_HALF = -0.69314718056;

    
  // Value for node_2 changes fastest in table
  int count = 0;
  for (int gt_1 = min_allele; gt_1 <= max_allele; ++gt_1){
    for (int gt_2 = min_allele; gt_2 <= max_allele; ++gt_2){
      double val = 0.0;
      
      // Iterate over all reads for the sample
      for (unsigned int k = 0; k < nrepeats.size(); ++k){
	// TO DO: Incorporate dropout probabilites by replacing LOG_ONE_HALF with unequal values
	val += log_sum_exp(LOG_ONE_HALF + snp_log_p1s[k] + model.calc_log_stutter(gt_1, nrepeats[k]), 
			   LOG_ONE_HALF + snp_log_p2s[k] + model.calc_log_stutter(gt_2, nrepeats[k]));
      }

      // Record factor probability for gt_2, gt_1
      factors.back().set(count, val);
      count++;
    }
  }
}


void add_read_factors(AlnVector& paired_str_reads, AlnVector& mate_reads, AlnVector& unpaired_str_reads,
		      BaseQuality& base_qualities, SNPTree* snp_tree, int min_allele, int max_allele, StutterModel& model,
		      dai::Var node_id1, dai::Var node_id2,
		      std::vector<dai::Factor>& factors){
  std::vector<int> nrepeats;
  std::vector<double> snp_log_p1s, snp_log_p2s;

  // Extract the number of repeats associated with each STR-spanning read

  // Extract the phasing qualities determined from overlapping heterozygous SNPs
  calc_het_snp_factors(paired_str_reads, mate_reads, base_qualities, snp_tree, snp_log_p1s, snp_log_p2s);
  calc_het_snp_factors(unpaired_str_reads, base_qualities, snp_tree, snp_log_p1s, snp_log_p2s);

  // Add appropriate factor
  add_read_factors(nrepeats, snp_log_p1s, snp_log_p2s, min_allele, max_allele, model, node_id1, node_id2, factors);
}
