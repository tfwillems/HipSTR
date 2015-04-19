#include <iostream>

#include "em_stutter_genotyper.h"
#include "extract_indels.h"
#include "length_em_bam_processor.h"

#include "seq_stutter_genotyper.h"

void LengthEMBamProcessor::analyze_reads_and_phasing(std::vector< std::vector<BamTools::BamAlignment> >& alignments,
						     std::vector< std::vector<double> >& log_p1s,
						     std::vector< std::vector<double> >& log_p2s,
						     std::vector<std::string>& rg_names, Region& region, std::string& ref_allele, std::string& chrom_seq){
  
  // Exploratory
  
  StutterModel stutter_model(0.9, 0.05, 0.05, 0.7, 0.005, 0.005, region.period());
  SeqStutterGenotyper seq_genotyper(region, alignments, log_p1s, log_p2s, rg_names, chrom_seq, stutter_model);
  seq_genotyper.genotype();
  seq_genotyper.write_vcf_record(samples_to_genotype_, str_vcf_);
  return;
  
  // End of exploratory section
  

  std::vector< std::vector<int> > str_bp_lengths(alignments.size());
  std::vector< std::vector<double> > str_log_p1s(alignments.size()), str_log_p2s(alignments.size());
  int inf_reads = 0;
  for (unsigned int i = 0; i < alignments.size(); ++i){
    for (unsigned int j = 0; j < alignments[i].size(); ++j){
      int bp_diff;
      bool got_size = ExtractCigar(alignments[i][j].CigarData, alignments[i][j].Position, region.start()-region.period(), region.stop()+region.period(), bp_diff);
      if (got_size){
	if (bp_diff < -(int)(region.stop()-region.start()+1)) {
	  std::cerr << "WARNING: Excluding read with bp difference greater than reference allele: " << alignments[i][j].Name << std::endl;
	  continue;
	}
	inf_reads++;
	str_bp_lengths[i].push_back(bp_diff);
	if (log_p1s.size() == 0){
	  str_log_p1s[i].push_back(0); str_log_p2s[i].push_back(0); // Assign equal phasing LLs as no SNP info is available
	}
	else {
	  str_log_p1s[i].push_back(log_p1s[i][j]); str_log_p2s[i].push_back(log_p2s[i][j]);
	}
      }
    }
  }
  
  // Train stutter model and genotype each sample
  std::cerr << "Building EM stutter genotyper" << std::endl;
  EMStutterGenotyper stutter_genotyper(region.chrom(), region.start(), region.stop(), str_bp_lengths, str_log_p1s, str_log_p2s, rg_names, region.period(), 0);
  
  std::cerr << "Training EM stutter genotyper" << std::endl;
  bool trained = stutter_genotyper.train(MAX_EM_ITER, ABS_LL_CONVERGE, FRAC_LL_CONVERGE);
  if (trained){
    num_em_converge_++;
    std::cerr << "Learned stutter model: " << *(stutter_genotyper.get_stutter_model()) << std::endl;
    bool use_pop_freqs = false;
    stutter_genotyper.genotype(use_pop_freqs);
    
    if (output_str_gts_)
      stutter_genotyper.write_vcf_record(ref_allele, samples_to_genotype_, str_vcf_);
  }
  else {
    num_em_fail_++;
    std::cerr << "Stutter model training failed for locus " << region.chrom() << ":" << region.start() << "-" << region.stop() 
	      << " with " << inf_reads << " informative reads" << std::endl;
  }
}

