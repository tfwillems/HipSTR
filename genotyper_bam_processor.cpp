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

/*
  Left align BamAlignments in the provided vector and store those that successfully realign in the provided vector.
  Also extracts other information for successfully realigned reads into provided vectors.
 */
void GenotyperBamProcessor::left_align_reads(Region& region, std::string& chrom_seq, std::vector< std::vector<BamTools::BamAlignment> >& alignments,
					     std::vector< std::vector<double> >& log_p1,       std::vector< std::vector<double> >& log_p2,
					     std::vector< std::vector<double> >& filt_log_p1,  std::vector< std::vector<double> >& filt_log_p2,
					     std::vector< Alignment>& left_alns, std::vector<int>& bp_diffs, std::vector<bool>& use_for_hap_generation,
					     std::ostream& logger){
  locus_left_aln_time_ = clock();
  logger << "Left aligning reads..." << std::endl;
  std::map<std::string, int> seq_to_alns;
  int32_t align_fail_count = 0, total_reads = 0;
  int bp_diff;
  left_alns.clear(); filt_log_p1.clear(); filt_log_p2.clear();
  bp_diffs.clear(); use_for_hap_generation.clear();

  for (unsigned int i = 0; i < alignments.size(); ++i){
    filt_log_p1.push_back(std::vector<double>());
    filt_log_p2.push_back(std::vector<double>());

    for (unsigned int j = 0; j < alignments[i].size(); ++j, ++total_reads){
      // Trim alignment if it extends very far upstream or downstream of the STR. For tractability, we limit it to 40bp
      trimAlignment(alignments[i][j], (region.start() > 40 ? region.start()-40 : 1), region.stop()+40);
      if (alignments[i][j].Length == 0)
        continue;

      auto iter      = seq_to_alns.find(alignments[i][j].QueryBases);
      bool have_prev = (iter != seq_to_alns.end());
      if (have_prev)
        have_prev &= left_alns[iter->second].get_sequence().size() == alignments[i][j].QueryBases.size();

      if (!have_prev){
        left_alns.push_back(Alignment(alignments[i][j].Name));
        if (matchesReference(alignments[i][j]))
          convertAlignment(alignments[i][j], chrom_seq, left_alns.back());
        else if (!realign(alignments[i][j], chrom_seq, left_alns.back())){
	  // Failed to realign read
          align_fail_count++;
          left_alns.pop_back();
          continue;
	}
	seq_to_alns[alignments[i][j].QueryBases] = left_alns.size()-1;
      }
      else {
        // Reuse alignments if the sequence has already been observed and didn't lead to a soft-clipped alignment
        // Soft-clipping is problematic because it complicates base quality extration (but not really that much)
        Alignment& prev_aln = left_alns[iter->second];
        assert(prev_aln.get_sequence().size() == alignments[i][j].QueryBases.size());
	std::string bases = uppercase(alignments[i][j].QueryBases);
        Alignment new_aln(prev_aln.get_start(), prev_aln.get_stop(), alignments[i][j].Name, alignments[i][j].Qualities, bases, prev_aln.get_alignment());
        new_aln.set_cigar_list(prev_aln.get_cigar_list());
        left_alns.push_back(new_aln);
      }

      left_alns.back().check_CIGAR_string(alignments[i][j].Name); // Ensure alignment is properly formatted
      filt_log_p1[i].push_back(log_p1[i][j]);
      filt_log_p2[i].push_back(log_p2[i][j]);
      bool got_size = ExtractCigar(alignments[i][j].CigarData, alignments[i][j].Position, region.start()-region.period(), region.stop()+region.period(), bp_diff);
      bp_diffs.push_back(got_size ? bp_diff : -999);
      use_for_hap_generation.push_back(BamProcessor::passes_filters(alignments[i][j]));
    }
  }

  locus_left_aln_time_  = (clock() - locus_left_aln_time_)/CLOCKS_PER_SEC;
  total_left_aln_time_ += locus_left_aln_time_;
  if (align_fail_count != 0)
    logger << "Failed to left align " << align_fail_count << " out of " << total_reads << " reads" << std::endl;
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
  // Can't simply check the total number of reads because the bam processor may have stopped reading at the threshold and then removed PCR duplicates
  // Instead, we check this flag which it sets when too many reads are encountered during filtering
  if (TOO_MANY_READS){
    logger() << "Skipping locus with too many reads: MAX=" << MAX_TOTAL_READS << std::endl;
    return;
  }

  assert(alignments.size() == log_p1s.size() && alignments.size() == log_p2s.size() && alignments.size() == rg_names.size());
  std::vector< std::vector<int> > str_bp_lengths(alignments.size());
  std::vector< std::vector<double> > str_log_p1s(alignments.size()), str_log_p2s(alignments.size());
  int inf_reads = 0;

  // Extract bp differences and phasing probabilities for each read if we 
  // need to utilize the length-based EM genotyper for stutter model training or genotyping
  int skip_count = 0;
  if ((def_stutter_model_ == NULL) && !read_stutter_models_){
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
  if (def_stutter_model_ != NULL){
    log("Using default stutter model");
    stutter_model = def_stutter_model_->copy();
    stutter_model->set_period(region.period());
  }
  else if (read_stutter_models_){
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
    length_genotyper = new EMStutterGenotyper(region, haploid, str_bp_lengths, str_log_p1s, str_log_p2s, rg_names, 0);
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
      VCF::VCFReader* reference_panel_vcf = NULL;
      if (ref_vcf_ != NULL)
	reference_panel_vcf = ref_vcf_;

      std::vector<Alignment> left_alignments;
      std::vector< std::vector<double> > filt_log_p1s, filt_log_p2s;
      std::vector<bool> use_to_generate_haps;
      std::vector<int> bp_diffs;
      left_align_reads(region, chrom_seq, alignments, log_p1s, log_p2s, filt_log_p1s,
		       filt_log_p2s, left_alignments, bp_diffs, use_to_generate_haps, logger());

      seq_genotyper = new SeqStutterGenotyper(region, haploid, left_alignments, use_to_generate_haps, bp_diffs, filt_log_p1s, filt_log_p2s, rg_names, chrom_seq, pool_seqs_,
					      *stutter_model, reference_panel_vcf, logger());

      if (output_str_gts_){
	if (seq_genotyper->genotype(chrom_seq, logger())) {
	  bool pass = true;

	  // If appropriate, recalculate the stutter model using the haplotype ML alignments,
	  // realign the reads and regenotype the samples
	  if (recalc_stutter_model_)
	    pass = seq_genotyper->recompute_stutter_models(chrom_seq, logger(), MAX_EM_ITER, ABS_LL_CONVERGE, FRAC_LL_CONVERGE);

	  if (pass){
	    num_genotype_success_++;
	    seq_genotyper->write_vcf_record(samples_to_genotype_, true, chrom_seq, output_bstrap_quals_, output_gls_, output_pls_, output_phased_gls_,
					    output_all_reads_, output_pall_reads_, output_mall_reads_, output_viz_, max_flank_indel_frac_,
					    viz_left_alns_, viz_out_, str_vcf_, logger());
	  }
	  else
	    num_genotype_fail_++;
	}
	else
	  num_genotype_fail_++;
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
    if (output_str_gts_){
      assert(seq_genotyper != NULL);
      logger() << "\t" << " Left alignment        = "  << locus_left_aln_time_             << " seconds\n"
	       << "\t" << " Haplotype generation  = "  << seq_genotyper->hap_build_time()  << " seconds\n"
	       << "\t" << " Haplotype alignment   = "  << seq_genotyper->hap_aln_time()    << " seconds\n"
	       << "\t" << " Posterior computation = "  << seq_genotyper->posterior_time()  << " seconds\n"
	       << "\t" << " Alignment traceback   = "  << seq_genotyper->aln_trace_time()  << " seconds\n"
	       << "\t" << " Bootstrap computation = "  << seq_genotyper->bootstrap_time()  << " seconds\n";

      process_timer_.add_time("Left alignment",        locus_left_aln_time_);
      process_timer_.add_time("Haplotype generation",  seq_genotyper->hap_build_time());
      process_timer_.add_time("Haplotype alignment",   seq_genotyper->hap_aln_time());
      process_timer_.add_time("Posterior computation", seq_genotyper->posterior_time());
      process_timer_.add_time("Alignment traceback",   seq_genotyper->aln_trace_time());
      process_timer_.add_time("Bootstrap computation", seq_genotyper->bootstrap_time());
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
 

