#include "snp_bam_processor.h"
#include "snp_phasing_quality.h"
#include "snp_tree.h"

void SNPBamProcessor::process_reads(std::vector< std::vector<BamTools::BamAlignment> >& paired_strs_by_rg,
				    std::vector< std::vector<BamTools::BamAlignment> >& mate_pairs_by_rg,
				    std::vector< std::vector<BamTools::BamAlignment> >& unpaired_strs_by_rg,
				    std::vector<std::string>& rg_names, Region& region, 
				    std::string& ref_allele, std::string& chrom_seq, std::ostream& out){
  assert(paired_strs_by_rg.size() == mate_pairs_by_rg.size() && paired_strs_by_rg.size() == unpaired_strs_by_rg.size());
  if (paired_strs_by_rg.size() == 0 && unpaired_strs_by_rg.size() == 0)
    return;
  
  std::vector<  std::vector<BamTools::BamAlignment> > alignments(paired_strs_by_rg.size());
  std::vector< std::vector<double> > log_p1s, log_p2s;
  bool got_snp_info = false;
  if (have_snp_vcf){
    std::vector<SNPTree*> snp_trees;
    std::map<std::string, unsigned int> sample_indices;      
    if(create_snp_trees(region.chrom(), (region.start() > MAX_MATE_DIST ? region.start()-MAX_MATE_DIST : 1), 
			region.stop()+MAX_MATE_DIST, phased_snp_vcf, sample_indices, snp_trees)){
      got_snp_info = true;
      std::set<std::string> bad_samples, good_samples;
      for (unsigned int i = 0; i < paired_strs_by_rg.size(); ++i){
	if (sample_indices.find(rg_names[i]) != sample_indices.end()){
	  good_samples.insert(rg_names[i]);
	  std::vector<double> log_p1, log_p2;
	  SNPTree* snp_tree = snp_trees[sample_indices[rg_names[i]]];
	  calc_het_snp_factors(paired_strs_by_rg[i], mate_pairs_by_rg[i], base_quality_, snp_tree, log_p1, log_p2, match_count_, mismatch_count_);
	  calc_het_snp_factors(unpaired_strs_by_rg[i], base_quality_, snp_tree, log_p1, log_p2, match_count_, mismatch_count_);
	  log_p1s.push_back(log_p1); log_p2s.push_back(log_p2);
	}
	else {
	  std::vector<double> log_p1, log_p2;
	  for (unsigned int j = 0; j < paired_strs_by_rg[i].size()+unpaired_strs_by_rg[i].size(); ++j){
	    log_p1.push_back(0); log_p2.push_back(0); // Assign equal phasing LLs as no SNP info is available
	  }
	  log_p1s.push_back(log_p1); log_p2s.push_back(log_p2);
	  bad_samples.insert(rg_names[i]);
	}
	
	// Copy alignments
	alignments[i].insert(alignments[i].end(), paired_strs_by_rg[i].begin(),   paired_strs_by_rg[i].end());
	alignments[i].insert(alignments[i].end(), unpaired_strs_by_rg[i].begin(), unpaired_strs_by_rg[i].end());
      }
      std::cerr << "Found VCF info for " << good_samples.size() << " out of " << good_samples.size()+bad_samples.size() << " samples with STR reads" << std::endl;
    }
    else 
      std::cerr << "Warning: Failed to construct SNP trees for " << region.chrom() << ":" << region.start() << "-" << region.stop() << std::endl; 
    destroy_snp_trees(snp_trees);      
  }
  if (!got_snp_info){
    for (unsigned int i = 0; i < paired_strs_by_rg.size(); i++){
      // Copy alignments                                                                                                                                             
      alignments[i].insert(alignments[i].end(), paired_strs_by_rg[i].begin(),   paired_strs_by_rg[i].end());
      alignments[i].insert(alignments[i].end(), unpaired_strs_by_rg[i].begin(), unpaired_strs_by_rg[i].end());
      
      // Assign equal phasing LLs as no SNP info is available
      log_p1s.push_back(std::vector<double>(paired_strs_by_rg[i].size()+unpaired_strs_by_rg[i].size(), 0.0));
      log_p2s.push_back(std::vector<double>(paired_strs_by_rg[i].size()+unpaired_strs_by_rg[i].size(), 0.0));
    }
  }
  
  int phased_samples = 0, phased_reads = 0, total_reads = 0;
  for (unsigned int i = 0; i < alignments.size(); i++){
    bool sample_phased = false;
    for (unsigned int j = 0; j < alignments[i].size(); j++){
      sample_phased |=  (log_p1s[i][j] != log_p2s[i][j]);
      phased_reads  += (log_p1s[i][j] != log_p2s[i][j]);
    }
    total_reads    += alignments[i].size();
    phased_samples += sample_phased;
  }

  std::cout << "Phased SNPs add info for " << phased_reads << " out of " << total_reads << " reads" 
	    << " and " << phased_samples << " out of " << rg_names.size() <<  " samples" << std::endl;

  // Run any additional analyses using phasing probabilities
  analyze_reads_and_phasing(alignments, log_p1s, log_p2s, rg_names, region, ref_allele, chrom_seq);
}




