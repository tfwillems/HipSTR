#include <iomanip>
#include <iostream>
#include <time.h>

//#include "sys/sysinfo.h"
//#include "sys/types.h"

#include "extract_indels.h"
#include "genotyper_bam_processor.h"

int parseLine(char* line){
  int i = strlen(line);
  while (*line < '0' || *line > '9') line++;
  line[i-3] = '\0';
  i = atoi(line);
  return i;
}

int getUsedPhysicalMemoryKB(){
  FILE* file = fopen("/proc/self/status", "r");
  int result = -1;
  char line[128];

  while (fgets(line, 128, file) != NULL){
    if (strncmp(line, "VmRSS:", 6) == 0){
      result = parseLine(line);
      break;
    }
  }
  fclose(file);
  return result;
}

void GenotyperBamProcessor::analyze_reads_and_phasing(std::vector< std::vector<BamTools::BamAlignment> >& alignments,
						      std::vector< std::vector<double> >& log_p1s,
						      std::vector< std::vector<double> >& log_p2s,
						      std::vector<std::string>& rg_names, Region& region, std::string& ref_allele, std::string& chrom_seq, int iter){
  int32_t total_reads = 0;
  for (unsigned int i = 0; i < alignments.size(); i++)
    total_reads += alignments[i].size();
  if (total_reads < MIN_TOTAL_READS){
    logger() << "Skipping locus with too few reads: TOTAL=" << total_reads << ", MIN=" << MIN_TOTAL_READS << std::endl;
    return;
  }
  if (total_reads > MAX_TOTAL_READS){
    logger() << "Skipping locus with too many reads: TOTAL=" << total_reads << ", MAX=" << MAX_TOTAL_READS << std::endl;
    return;
  }

  assert(alignments.size() == log_p1s.size() && alignments.size() == log_p2s.size() && alignments.size() == rg_names.size());
  std::vector< std::vector<int> > str_bp_lengths(alignments.size());
  std::vector< std::vector<double> > str_log_p1s(alignments.size()), str_log_p2s(alignments.size());
  int inf_reads = 0;

  // Extract bp differences and phasing probabilities for each read if we 
  // need to utilize the length-based EM genotyper for stutter model training or genotyping
  int skip_count = 0;
  if (!read_stutter_models_ || !use_seq_aligner_){
    for (unsigned int i = 0; i < alignments.size(); ++i){
      for (unsigned int j = 0; j < alignments[i].size(); ++j){
	int bp_diff;
	bool got_size = ExtractCigar(alignments[i][j].CigarData, alignments[i][j].Position, region.start()-region.period(), region.stop()+region.period(), bp_diff);
	if (got_size){
	  if (bp_diff < -(int)(region.stop()-region.start()+1)) {
	    log("WARNING: Excluding read with bp difference greater than reference allele: " +alignments[i][j].Name);
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
	else
	  skip_count++;
      }
    }
  }

  if (total_reads-skip_count < MIN_TOTAL_READS){
    logger() << "Skipping locus with too few reads: TOTAL=" << total_reads-skip_count << ", MIN=" << MIN_TOTAL_READS << std::endl;
    return;
  }

  bool haploid = (haploid_chroms_.find(region.chrom()) != haploid_chroms_.end());
  bool trained = false;
  StutterModel* stutter_model          = NULL;
  EMStutterGenotyper* length_genotyper = NULL;
  locus_stutter_time_ = clock();
  if (read_stutter_models_){
    // Attempt to extact model from dictionary
    auto model_iter = stutter_models_.find(region);
    if (model_iter != stutter_models_.end())
      stutter_model = model_iter->second->copy();
    else
      logger() << "WARNING: No stutter model found for " << region.chrom() << ":" << region.start() << "-" << region.stop() << std::endl;
  }
  else {
    // Learn stutter model using length-based EM algorithm
    log("Building EM stutter genotyper");
    length_genotyper = new EMStutterGenotyper(region.chrom(), region.start(), region.stop(), haploid, str_bp_lengths,
					      str_log_p1s, str_log_p2s, rg_names, region.period(), 0);
    log("Training EM stutter genotyper");
    trained = length_genotyper->train(MAX_EM_ITER, ABS_LL_CONVERGE, FRAC_LL_CONVERGE, false, logger());
    if (trained){
      if (output_stutter_models_)
	length_genotyper->get_stutter_model()->write_model(region.chrom(), region.start(), region.stop(), stutter_model_out_);
      num_em_converge_++;
      stutter_model = length_genotyper->get_stutter_model()->copy();
      logger() << "Learned stutter model: " << *stutter_model << std::endl;
    }
    else {
      num_em_fail_++;
      logger() << "Stutter model training failed for locus " << region.chrom() << ":" << region.start() << "-" << region.stop()
	       << " with " << inf_reads << " informative reads" << std::endl;
    }
  }
  locus_stutter_time_  = (clock() - locus_stutter_time_)/CLOCKS_PER_SEC;;
  total_stutter_time_ += locus_stutter_time_;

  SeqStutterGenotyper* seq_genotyper = NULL;
  locus_genotype_time_ = clock();
  if (output_str_gts_){
    if (stutter_model != NULL) {
      if (use_seq_aligner_){
	// Use sequence-based genotyper
	vcflib::VariantCallFile* reference_panel_vcf = NULL;
	if (have_ref_vcf_)
	  reference_panel_vcf = &ref_vcf_;

	seq_genotyper = new SeqStutterGenotyper(region, haploid, alignments, log_p1s, log_p2s, rg_names, chrom_seq, pool_seqs_,
						*stutter_model, reference_panel_vcf, logger());
	if (output_alleles_){
	  std::vector<std::string> no_samples;
	  std::vector<int> read_str_sizes;
	  seq_genotyper->write_vcf_record(no_samples, false, chrom_seq, false, false, false,
					  false, false, false, false, false, read_str_sizes, viz_out_, allele_vcf_, logger());
	}

	if (output_str_gts_){
	  if (seq_genotyper->genotype(logger())) {
	    bool pass = true;

	    // If appropriate, recalculate the stutter model using the haplotype ML alignments,
	    // realign the reads and regenotype the samples
	    if (recalc_stutter_model_)
	      pass = seq_genotyper->recompute_stutter_model(chrom_seq, logger(), MAX_EM_ITER, ABS_LL_CONVERGE, FRAC_LL_CONVERGE);

	    if (pass){
	      num_genotype_success_++;
	      std::vector<int> read_str_sizes;
	      seq_genotyper->write_vcf_record(samples_to_genotype_, true, chrom_seq, output_bstrap_quals_, output_gls_, output_pls_,
					      output_all_reads_, output_pall_reads_, output_mall_reads_, output_viz_, viz_left_alns_,
					      read_str_sizes, viz_out_, str_vcf_, logger());
	    }
	    else
	      num_genotype_fail_++;
	  }
	  else
	    num_genotype_fail_++;
	}
      }
      else {
	// Use length-based genotyper
	if (length_genotyper == NULL){
	  length_genotyper = new EMStutterGenotyper(region.chrom(), region.start(), region.stop(), haploid,
						    str_bp_lengths, str_log_p1s, str_log_p2s, rg_names, region.period(), 0);
	  length_genotyper->set_stutter_model(*stutter_model);
	}

	if (output_str_gts_){
	  bool use_pop_freqs = false;
	  if (length_genotyper->genotype(use_pop_freqs)){
	    num_genotype_success_++;
	    if (output_str_gts_)
	      length_genotyper->write_vcf_record(ref_allele, samples_to_genotype_, output_gls_, output_pls_, output_all_reads_, str_vcf_);
	  }
	  else
	    num_genotype_fail_++;
	}
      }
    }
  }
  locus_genotype_time_  = (clock() - locus_genotype_time_)/CLOCKS_PER_SEC;
  total_genotype_time_ += locus_genotype_time_;

  logger() << "Locus timing:"                                          << "\n"
	   << " BAM seek time       = " << locus_bam_seek_time()       << " seconds\n"
	   << " Read filtering      = " << locus_read_filter_time()    << " seconds\n"
	   << " SNP info extraction = " << locus_snp_phase_info_time() << " seconds\n"
	   << " Stutter estimation  = " << locus_stutter_time()        << " seconds\n";
  if (stutter_model != NULL){
    logger() << " Genotyping          = " << locus_genotype_time()       << " seconds\n";
    if (use_seq_aligner_ && output_str_gts_){
      assert(seq_genotyper != NULL);
      logger() << "\t" << " Left alignment       = "  << seq_genotyper->left_aln_time()   << " seconds\n"
	       << "\t" << " Haplotype generation = "  << seq_genotyper->hap_build_time()  << " seconds\n"
	       << "\t" << " Haplotype alignment  = "  << seq_genotyper->hap_aln_time()    << " seconds\n"
	       << "\t" << " Alignment filtering  = "  << seq_genotyper->aln_filter_time() << " seconds\n"
	       << "\t" << " Alignment traceback  = "  << seq_genotyper->aln_trace_time()  << " seconds\n";

      process_timer_.add_time("Left alignment",       seq_genotyper->left_aln_time());
      process_timer_.add_time("Haplotype generation", seq_genotyper->hap_build_time());
      process_timer_.add_time("Haplotype alignment",  seq_genotyper->hap_aln_time());
      process_timer_.add_time("Alignment filtering",  seq_genotyper->aln_filter_time());
      process_timer_.add_time("Alignment traceback",  seq_genotyper->aln_trace_time());
    }
  }

  /*
  logger() << "Total memory in use = " << getUsedPhysicalMemoryKB() << " KB"
	   << std::endl;
  */

  delete seq_genotyper;
  delete stutter_model;
  delete length_genotyper;
}
 

