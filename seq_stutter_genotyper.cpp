#include <algorithm>
#include <cfloat>
#include <random>
#include <sstream>
#include <time.h>

#include "seq_stutter_genotyper.h"
#include "em_stutter_genotyper.h"
#include "error.h"
#include "extract_indels.h"
#include "mathops.h"

#include "read_pooler.h"

#include "stringops.h"
#include "vcf_input.h"

#include "SeqAlignment/AlignmentOps.h"
#include "SeqAlignment/AlignmentData.h"
#include "SeqAlignment/AlignmentModel.h"
#include "SeqAlignment/AlignmentViz.h"
#include "SeqAlignment/HaplotypeGenerator.h"
#include "SeqAlignment/HapAligner.h"
#include "SeqAlignment/RepeatStutterInfo.h"
#include "SeqAlignment/RepeatBlock.h"
#include "SeqAlignment/STRAlleleExpansion.h"

bool SeqStutterGenotyper::condense_read_count_fields = true;

int max_index(double* vals, unsigned int num_vals){
  int best_index = 0;
  for (unsigned int i = 1; i < num_vals; i++)
    if (vals[i] > vals[best_index])
      best_index = i;
  return best_index;
}

void SeqStutterGenotyper::get_uncalled_alleles(std::vector<int>& allele_indices){
  assert(allele_indices.size() == 0);
 
  // Determine which samples have >= 1 aligned read
  std::vector<bool> aligned_read(num_samples_, false);
  for (unsigned int read_index = 0; read_index < num_reads_; read_index++){
    if (seed_positions_[read_index] >= 0)
      aligned_read[sample_label_[read_index]] = true;
  }

  // Extract each sample's MAP genotype
  std::vector< std::pair<int,int> > gts;
  get_optimal_genotypes(gts);

  // Mark all alleles with a call by a valid sample
  std::vector<bool> called(num_alleles_, false);
  for (unsigned int i = 0; i < gts.size(); i++){
    if ((!require_one_read_ || aligned_read[i]) && call_sample_[i]){
      called[gts[i].first]  = true;
      called[gts[i].second] = true;
    }
  }

  // Unmarked alleles are uncalled (apart from reference allele which must always be kept)
  for (unsigned int i = 1; i < called.size(); i++)
    if (!called[i])
      allele_indices.push_back(i);
}

void SeqStutterGenotyper::remove_alleles(std::vector<int>& allele_indices){
  assert(log_allele_priors_ == NULL);           // Can't use this option if priors have been set
  assert(allele_indices.size() < num_alleles_); // Make sure we'll have at least 1 allele

  std::vector<bool> keep_allele(num_alleles_, true);
  for (auto iter = allele_indices.begin(); iter != allele_indices.end(); iter++){
    assert(*iter < keep_allele.size() && *iter >= 0);
    assert(keep_allele[*iter] == true);
    keep_allele[*iter] = false;
  }

  int fixed_num_alleles = num_alleles_ - allele_indices.size();
  std::vector<std::string> fixed_alleles;
  for (unsigned int i = 0; i < alleles_.size(); i++)
    if (keep_allele[i])
      fixed_alleles.push_back(alleles_[i]);
  
  // Fix read alignment probability array
  double* fixed_log_aln_probs = new double[fixed_num_alleles*num_reads_];
  double* old_log_aln_ptr     = log_aln_probs_;
  double* new_log_aln_ptr     = fixed_log_aln_probs;;
  for (unsigned int i = 0; i < num_reads_; ++i){
    for (unsigned int j = 0; j < num_alleles_; ++j, ++old_log_aln_ptr){
      if (keep_allele[j]){
	*new_log_aln_ptr = *old_log_aln_ptr;
	new_log_aln_ptr++;
      }
    }
  }
  delete [] log_aln_probs_;
  log_aln_probs_ = fixed_log_aln_probs;

  // Replace other variables
  num_alleles_ = fixed_num_alleles;
  alleles_     = fixed_alleles;

  // Rebuild the haplotype
  assert(haplotype_->num_blocks() == 3);
  assert(haplotype_->get_block(1)->get_repeat_info() != NULL);
  RepeatBlock* new_str_block = ((RepeatBlock*)(haplotype_->get_block(1)))->remove_alleles(allele_indices);
  delete hap_blocks_[1];
  delete haplotype_;
  hap_blocks_[1] = new_str_block;
  haplotype_     = new Haplotype(hap_blocks_);


  // TO DO: Consider reshrinking allele bounds as in the get_alleles function

  // Resize and recalculate genotype posterior array
  delete [] log_sample_posteriors_;
  log_sample_posteriors_ = new double[fixed_num_alleles*fixed_num_alleles*num_samples_];
  calc_log_sample_posteriors();
}

void SeqStutterGenotyper::combine_reads(std::vector<Alignment>& alignments, Alignment& pooled_aln){
  assert(alignments.size() > 0);
  pooled_aln.set_start(alignments[0].get_start());
  pooled_aln.set_stop(alignments[0].get_stop());
  pooled_aln.set_sample("");
  pooled_aln.set_sequence(alignments[0].get_sequence());
  pooled_aln.set_alignment(alignments[0].get_alignment());
  pooled_aln.set_cigar_list(alignments[0].get_cigar_list());

  // Utilize mean base quality scores for pooled alignment
  std::vector<const std::string*> qual_ptrs;
  for (unsigned int i = 0; i < alignments.size(); i++)
    qual_ptrs.push_back(&(alignments[i].get_base_qualities()));
  std::string mean_base_quals = base_quality_.average_base_qualities(qual_ptrs);
  assert(mean_base_quals.size() == alignments[0].get_sequence().size());
  pooled_aln.set_base_qualities(mean_base_quals);
}


bool SeqStutterGenotyper::expand_haplotype(std::ostream& logger){
  assert(haplotype_->num_blocks() == 3);
  std::vector< std::vector<std::string> > read_seqs(alns_.size());
  for (unsigned int i = 0; i < alns_.size(); ++i){
    read_seqs.push_back(std::vector<std::string>());
    for (unsigned int j = 0; j < alns_[i].size(); ++j)
      read_seqs[i].push_back(alns_[i][j].get_sequence());
  }

  // Check that the flanking sequences are sufficiently large to support the operations
  if (std::min(haplotype_->get_block(0)->get_seq(0).size(),  haplotype_->get_block(2)->get_seq(0).size()) <= 5)
    return false;

  // Identify additional candidate STR alleles
  int flank_match = 5;
  std::vector<std::string> str_seqs;
  std::string lflank = haplotype_->get_block(0)->get_seq(0).substr(haplotype_->get_block(0)->get_seq(0).size()-flank_match, flank_match);
  std::string rflank = haplotype_->get_block(2)->get_seq(0).substr(0, flank_match);
  for (unsigned int i = 0; i < haplotype_->get_block(1)->num_options(); i++)
    str_seqs.push_back(lflank + haplotype_->get_block(1)->get_seq(i) + rflank);
  std::set<std::string> new_str_seqs;
  get_candidates(str_seqs, read_seqs, region_->period(), new_str_seqs, logger);

  // Extract actual STR sequences from existing haplotype and new alleles
  std::vector<std::string> final_str_seqs;
  for (unsigned int i = 0; i < haplotype_->get_block(1)->num_options(); i++)
    final_str_seqs.push_back(haplotype_->get_block(1)->get_seq(i));
  for (auto iter = new_str_seqs.begin(); iter != new_str_seqs.end(); iter++){
    if (iter->substr(0, flank_match).compare(lflank) == 0 && iter->substr(iter->size()-flank_match, flank_match).compare(rflank) == 0){
      std::string seq = iter->substr(flank_match, iter->size()-2*flank_match);
      final_str_seqs.push_back(seq);
      expanded_alleles_.insert(seq);
    }
  }
  std::sort(final_str_seqs.begin()+1, final_str_seqs.end(), stringLengthLT);
  
  // Construct the new STR block
  HapBlock* str_block = new RepeatBlock(hap_blocks_[1]->start(), hap_blocks_[1]->end(), hap_blocks_[1]->get_seq(0), region_->period(), stutter_model_);
  for (unsigned int i = 1; i < final_str_seqs.size(); i++)
    str_block->add_alternate(final_str_seqs[i]);
  str_block->initialize();
  
  // Reconstruct the haplotype
  delete haplotype_;
  delete hap_blocks_[1];
  hap_blocks_[1] = str_block;
  haplotype_     = new Haplotype(hap_blocks_);
  num_alleles_   = haplotype_->num_combs();
  logger << "Haplotype after additional allele identification:" << std::endl;
  haplotype_->print_block_structure(30, 100, logger);
  return true;
}

void SeqStutterGenotyper::init(std::vector< std::vector<BamTools::BamAlignment> >& alignments, 
			       std::vector< std::vector<double> >& log_p1, 
			       std::vector< std::vector<double> >& log_p2,
			       std::vector<std::string>& sample_names,
			       std::string& chrom_seq, std::ostream& logger){
  // Compute the total number of reads
  num_reads_ = 0;
  for (unsigned int i = 0; i < alignments.size(); ++i)
    num_reads_ += alignments[i].size();

  // Allocate some data structures
  log_p1_           = new double[num_reads_];
  log_p2_           = new double[num_reads_];
  sample_label_     = new int[num_reads_];
  sample_total_LLs_ = new double[num_samples_];

  double locus_left_aln_time = clock();
  logger << "Left aligning reads..." << std::endl;
  std::map<std::string, std::pair<int,int> > seq_to_alns;
  int read_index = 0, align_fail_count = 0, qual_filt_count = 0;
  int bp_diff;
  for (unsigned int i = 0; i < alignments.size(); ++i){
    alns_.push_back(std::vector<Alignment>());
    for (unsigned int j = 0; j < alignments[i].size(); ++j, ++read_index){
      // Ignore/remove reads with a very low overall base quality score
      // Want to avoid situations in which it's more advantageous to have misalignments
      // b/c the scores are so low
      if (base_quality_.sum_log_prob_correct(alignments[i][j].Qualities) < MIN_SUM_QUAL_LOG_PROB){
	qual_filt_count++;
	num_reads_--;
	read_index--;
	continue;
      }

      auto iter      = seq_to_alns.find(alignments[i][j].QueryBases);
      bool have_prev = (iter != seq_to_alns.end());
      if (have_prev)
	have_prev &= alns_[iter->second.first][iter->second.second].get_sequence().size() == alignments[i][j].QueryBases.size();

      if (!have_prev){
	alns_.back().push_back(Alignment());
	if (realign(alignments[i][j], chrom_seq, alns_.back().back())){
	  seq_to_alns[alignments[i][j].QueryBases] = std::pair<int,int>(i, alns_[i].size()-1);
	  alns_.back().back().check_CIGAR_string(alignments[i][j].Name);
	}
	else {
	  // Failed to realign read
	  align_fail_count++;
	  alns_.back().pop_back();
	  num_reads_--;
	  read_index--;
	  continue;
	}
      }
      else {
	// Reuse alignments if the sequence has already been observed and didn't lead to a soft-clipped alignment
	// Soft-clipping is problematic because it complicates base quality extration (but not really that much)
	Alignment& prev_aln = alns_[iter->second.first][iter->second.second];
	assert(prev_aln.get_sequence().size() == alignments[i][j].QueryBases.size());
	std::string sample; alignments[i][j].GetTag(SAMPLE_TAG, sample);
	std::string bases = uppercase(alignments[i][j].QueryBases);
	Alignment new_aln(prev_aln.get_start(), prev_aln.get_stop(), sample, alignments[i][j].Qualities, bases, prev_aln.get_alignment());
	new_aln.set_cigar_list(alns_[iter->second.first][iter->second.second].get_cigar_list());
	new_aln.check_CIGAR_string(alignments[i][j].Name);
	alns_.back().push_back(new_aln);
      }
      bool got_size = ExtractCigar(alignments[i][j].CigarData, alignments[i][j].Position,
				   region_->start()-region_->period(), region_->stop()+region_->period(), bp_diff);
      bp_diffs_.push_back(got_size ? bp_diff : -999);
      log_p1_[read_index]       = log_p1[i][j];
      log_p2_[read_index]       = log_p2[i][j];
      sample_label_[read_index] = i; 
    }
  }
  locus_left_aln_time  = (clock() - locus_left_aln_time)/CLOCKS_PER_SEC;
  total_left_aln_time_ += locus_left_aln_time;

  if (align_fail_count != 0)
    logger << "Failed to left align " << align_fail_count << " out of " << num_reads_ << " reads" << std::endl;
  if (qual_filt_count != 0)
    logger << "Filtered " << qual_filt_count << " reads due to low overall base qualities." << std::endl
	   << "If this value is high (>1% of reads), there may be an issue to with the base quality score encoding" << std::endl;

  // Replace base quality scores for 'N' bases with minimum quality score
  for (unsigned int i = 0; i < alns_.size(); ++i)
    for (unsigned int j = 0; j < alns_[i].size(); ++j)
      alns_[i][j].fix_N_base_qualities(base_quality_);

  double locus_hap_build_time = clock();
  std::vector<std::string> vcf_alleles;
  if (ref_vcf_ != NULL){
    bool success = false;
    if (ref_vcf_->formatTypes.find(PGP_KEY) == ref_vcf_->formatTypes.end()){
      // Read alleles from VCF
      logger << "Reading STR alleles from VCF" << std::endl;
      read_vcf_alleles(ref_vcf_, region_, alleles_, pos_, success);
      assert(log_allele_priors_ == NULL);
    }
    else {
      // Read alleles and priors for each sample's genotypes from VCF generated by PhasedBEAGLE
      logger << "Reading STR alleles and priors from VCF" << std::endl;
      log_allele_priors_ = extract_vcf_alleles_and_log_priors(ref_vcf_, region_, sample_indices_, alleles_, got_priors_, pos_, success, logger);
      assert(got_priors_.size() == num_samples_);
    }

    num_alleles_ = alleles_.size();
    if (success){
      assert(num_alleles_ > 0);

      // Construct the haplotype using the set of VCF alleles
      haplotype_ = generate_haplotype(pos_, *region_, MAX_REF_FLANK_LEN, chrom_seq, alleles_, stutter_model_, hap_blocks_, logger);
      
      // TO DO: Set vector based on proximity of indels to haplotype
      call_sample_ = std::vector<bool>(num_samples_, true);

      // If priors in the VCF, don't call samples without allele priors
      if (log_allele_priors_ != NULL){
	for (unsigned int i = 0; i < num_samples_; i++)
	  call_sample_[i] = (call_sample_[i] && got_priors_[i]);
      }
    }
    else
      pos_ = -1;
  }
  else {
    // Generate putative haplotypes and determine the number of alleles
    logger << "Generating putative haplotypes..." << std::endl;
    haplotype_   = generate_haplotype(*region_, MAX_REF_FLANK_LEN, chrom_seq, alns_, vcf_alleles, stutter_model_, alleles_from_bams_, hap_blocks_, call_sample_, logger);
    num_alleles_ = haplotype_->num_combs();
    assert(call_sample_.size() == num_samples_);

    // Reconstruct haplotype after looking for additional alleles
    if(expand_haplotype(logger))
      get_alleles(chrom_seq, alleles_); // Extract full STR sequence for each allele using annotated repeat region and the haplotype above
    else
      pos_ = -1;
  }
  locus_hap_build_time  = (clock() - locus_hap_build_time)/CLOCKS_PER_SEC;
  total_hap_build_time_ += locus_hap_build_time;

  if (pos_ != -1){
    // Print information about the haplotype and the stutter model
    logger << "Max block sizes: ";
    for (unsigned int i = 0; i < haplotype_->num_blocks(); i++)
      logger << haplotype_->get_block(i)->max_size() << " ";
    logger << std::endl << "Stutter model information" << std::endl;
    RepeatStutterInfo* stutter_info = hap_blocks_[1]->get_repeat_info();
    for (int i = stutter_info->max_deletion(); i <= stutter_info->max_insertion(); i += stutter_info->get_period())
      logger << i << " " << stutter_info->log_prob_pcr_artifact(0, i) << std::endl;
    logger << std::endl;
    
    // Allocate the remaining data structures
    log_sample_posteriors_ = new double[num_alleles_*num_alleles_*num_samples_];
    log_aln_probs_         = new double[num_reads_*num_alleles_];
    seed_positions_        = new int[num_reads_];
  }
  else
    logger << "WARNING: Unsuccessful initialization. " << std::endl;
}

bool SeqStutterGenotyper::genotype(std::ostream& logger){
  // Unsuccessful initialization. May be due to
  // 1) Failing to find the corresponding allele priors in the VCF (if one has been provided)
  // 2) Large deletion extending past STR
  if (pos_ == -1)
    return false;

  // If the smallest stutter block sequence is smaller than the maximum deletion size, the stutter aligner will fail
  // We could extend the stutter block to prevent this, but if this happens, the locus is probably not very high quality
  // As a result, for now, just abort the genotyping for this locus
  RepeatBlock* rep_block = (RepeatBlock *)haplotype_->get_block(1);
  if (rep_block->min_size() < std::abs(rep_block->get_repeat_info()->max_deletion()))
    return false;

  init_alignment_model();
  double locus_hap_aln_time = clock();
  HapAligner hap_aligner(haplotype_);

  if (pool_identical_seqs_){
    // Merge reads with identical sequences into pools
    logger << "Pooling reads with identical sequences..." << std::endl;
    std::vector< std::vector<int> > pool_indices(alns_.size());
    ReadPooler pooler;
    for (unsigned int i = 0; i < alns_.size(); i++){
      pool_indices[i].reserve(alns_[i].size());
      for (unsigned int j = 0; j < alns_[i].size(); j++)
	pool_indices[i].push_back(pooler.add_alignment(alns_[i][j]));
    }
    pooler.pool(base_quality_);

    // Align each pooled read to each haplotype
    logger << "Aligning pooled reads to each candidate haplotype..." << std::endl;
    std::vector<Alignment>& pooled_alns = pooler.get_alignments();
    double* log_pool_aln_probs = new double[pooled_alns.size()*num_alleles_];
    int* pool_seed_positions   = new int[pooled_alns.size()];
    hap_aligner.process_reads(pooled_alns, 0, &base_quality_, log_pool_aln_probs, pool_seed_positions);

    // Copy each pool reads's alignment probabilities to the entries for its constituent reads
    int* seed_ptr       = seed_positions_;
    double* log_aln_ptr = log_aln_probs_;
    for (unsigned int i = 0; i < alns_.size(); ++i){
      for (unsigned int j = 0; j < alns_[i].size(); ++j, ++seed_ptr){
	int pool_index = pool_indices[i][j];
	*seed_ptr = pool_seed_positions[pool_index];
	memcpy(log_aln_ptr, log_pool_aln_probs + num_alleles_*pool_index, num_alleles_*sizeof(double));
	log_aln_ptr += num_alleles_;
      }
    }
    delete [] log_pool_aln_probs;
    delete [] pool_seed_positions;
  }
  else {
    // Align each read against each candidate haplotype
    logger << "Aligning reads to each candidate haplotype..." << std::endl;
    int read_index = 0;
    for (unsigned int i = 0; i < alns_.size(); i++){
      hap_aligner.process_reads(alns_[i], read_index, &base_quality_, log_aln_probs_, seed_positions_); 
      read_index += alns_[i].size();
    }  
  }
  locus_hap_aln_time   = (clock() - locus_hap_aln_time)/CLOCKS_PER_SEC;
  total_hap_aln_time_ += locus_hap_aln_time;

  logger << "Genotyping samples..." << std::endl;
  if (stutter_model_ == NULL)
    printErrorAndDie("Must specify stutter model before running genotype()");
  calc_log_sample_posteriors();

  //debug_sample(sample_indices_["NA12878"]);

  // Remove alleles with no MAP genotype calls and recompute the posteriors
  if (log_allele_priors_ == NULL){
    std::vector<int> uncalled_indices;
    get_uncalled_alleles(uncalled_indices);
    if (uncalled_indices.size() != 0){
      logger << "Recomputing sample posteriors after removing " << uncalled_indices.size() << " uncalled alleles" << std::endl;
      remove_alleles(uncalled_indices);
    }
  }
  
  return true;
}

void SeqStutterGenotyper::write_vcf_header(std::vector<std::string>& sample_names, bool output_gls, bool output_pls, std::ostream& out){
  out << "##fileformat=VCFv4.1" << "\n";

  // Info field descriptors
  out << "##INFO=<ID=" << "INFRAME_PGEOM,"  << "Number=1,Type=Float,Description=\""   << "Parameter for in-frame geometric step size distribution"                      << "\">\n"
      << "##INFO=<ID=" << "INFRAME_UP,"     << "Number=1,Type=Float,Description=\""   << "Probability that stutter causes an in-frame increase in obs. STR size"        << "\">\n"
      << "##INFO=<ID=" << "INFRAME_DOWN,"   << "Number=1,Type=Float,Description=\""   << "Probability that stutter causes an in-frame decrease in obs. STR size"        << "\">\n"
      << "##INFO=<ID=" << "OUTFRAME_PGEOM," << "Number=1,Type=Float,Description=\""   << "Parameter for out-of-frame geometric step size distribution"                  << "\">\n"
      << "##INFO=<ID=" << "OUTFRAME_UP,"    << "Number=1,Type=Float,Description=\""   << "Probability that stutter causes an out-of-frame increase in obs. STR size"    << "\">\n"
      << "##INFO=<ID=" << "OUTFRAME_DOWN,"  << "Number=1,Type=Float,Description=\""   << "Probability that stutter causes an out-of-frame decrease in obs. STR size"    << "\">\n"
      << "##INFO=<ID=" << "BPDIFFS,"        << "Number=A,Type=Integer,Description=\"" << "Base pair difference of each alternate allele from the reference allele"      << "\">\n"
      << "##INFO=<ID=" << "START,"          << "Number=1,Type=Integer,Description=\"" << "Inclusive start coodinate for the repetitive portion of the reference allele" << "\">\n"
      << "##INFO=<ID=" << "END,"            << "Number=1,Type=Integer,Description=\"" << "Inclusive end coordinate for the repetitive portion of the reference allele"  << "\">\n"
      << "##INFO=<ID=" << "PERIOD,"         << "Number=1,Type=Integer,Description=\"" << "Length of STR motif"                                                          << "\">\n"
      << "##INFO=<ID=" << "REFAC,"          << "Number=1,Type=Integer,Description=\"" << "Reference allele count"                                                       << "\">\n"
      << "##INFO=<ID=" << "AC,"             << "Number=A,Type=Integer,Description=\"" << "Alternate allele counts"                                                      << "\">\n"
      << "##INFO=<ID=" << "NSKIP,"          << "Number=1,Type=Integer,Description=\"" << "Number of samples not genotyped due to various issues"                        << "\">\n"
      << "##INFO=<ID=" << "NFILT,"          << "Number=1,Type=Integer,Description=\"" << "Number of samples whose genotypes were filtered due to various issues"        << "\">\n";

  // Format field descriptors
  out << "##FORMAT=<ID=" << "GT"          << ",Number=1,Type=String,Description=\""  << "Genotype" << "\">" << "\n"
      << "##FORMAT=<ID=" << "GB"          << ",Number=1,Type=String,Description=\""  << "Base pair differences of genotype from reference"              << "\">" << "\n"
      << "##FORMAT=<ID=" << "Q"           << ",Number=1,Type=Float,Description=\""   << "Posterior probability of unphased genotype"                    << "\">" << "\n"
      << "##FORMAT=<ID=" << "PQ"          << ",Number=1,Type=Float,Description=\""   << "Posterior probability of phased genotype"                      << "\">" << "\n"
      << "##FORMAT=<ID=" << "DP"          << ",Number=1,Type=Integer,Description=\"" << "Number of valid reads used for sample's genotype"              << "\">" << "\n"
      << "##FORMAT=<ID=" << "DSNP"        << ",Number=1,Type=Integer,Description=\"" << "Number of reads with SNP phasing information"                  << "\">" << "\n"
      << "##FORMAT=<ID=" << "PSNP"        << ",Number=1,Type=String,Description=\""  << "Number of reads with SNPs supporting each haploid genotype"    << "\">" << "\n"
      << "##FORMAT=<ID=" << "PDP"         << ",Number=1,Type=String,Description=\""  << "Fractional reads supporting each haploid genotype"             << "\">" << "\n"
      << "##FORMAT=<ID=" << "BQ"          << ",Number=1,Type=Float,Description=\""   << "Bootstrapped quality score"                                    << "\">" << "\n"
      << "##FORMAT=<ID=" << "GLDIFF"      << ",Number=1,Type=Float,Description=\""   << "Difference in likelihood between the reported and next best genotypes" << "\">" << "\n"
      << "##FORMAT=<ID=" << "DFILT"       << ",Number=1,Type=Integer,Description=\"" << "Number of reads filtered due to various issues"                << "\">" << "\n"
      << "##FORMAT=<ID=" << "DSTUTTER"    << ",Number=1,Type=Integer,Description=\"" << "Number of reads with a stutter indel in the STR region"        << "\">" << "\n"
      << "##FORMAT=<ID=" << "DFLANKINDEL" << ",Number=1,Type=Integer,Description=\"" << "Number of reads with an indel in the regions flanking the STR" << "\">" << "\n"
      << "##FORMAT=<ID=" << "BPDOSE"      << ",Number=1,Type=Float,Description=\""   << "Posterior mean base pair difference from reference"            << "\">" << "\n";

  if (condense_read_count_fields)
    out << "##FORMAT=<ID=" << "ALLREADS"  << ",Number=1,Type=String,Description=\"" << "Base pair difference observed in each read's Needleman-Wunsch alignment" << "\">" << "\n"
	<< "##FORMAT=<ID=" << "MALLREADS" << ",Number=1,Type=String,Description=\"" << "Maximum likelihood bp diff in each read based on haplotype alignments"   << "\">" << "\n";
  else
    out << "##FORMAT=<ID=" << "ALLREADS"  << ",Number=.,Type=Integer,Description=\"" << "Base pair difference observed in each read's Needleman-Wunsch alignment" << "\">" << "\n"
	<< "##FORMAT=<ID=" << "MALLREADS" << ",Number=.,Type=Integer,Description=\"" << "Maximum likelihood bp diff in each read based on haplotype alignments"   << "\">" << "\n";
  out << "##FORMAT=<ID=" << "PALLREADS"   << ",Number=.,Type=Float,Description=\""   << "Expected bp diff in each read based on haplotype alignment probs"        << "\">" << "\n";

  if (output_gls)
    out << "##FORMAT=<ID=" << "GL" << ",Number=G,Type=Float,Description=\""   << "log-10 genotype likelihoods" << "\">" << "\n";
  if (output_pls)
    out << "##FORMAT=<ID=" << "PL" << ",Number=G,Type=Integer,Description=\"" << "Phred-scaled genotype likelihoods" << "\">" << "\n";

  // Sample names
  out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  for (unsigned int i = 0; i < sample_names.size(); i++)
    out << "\t" << sample_names[i];
  out << "\n";
}


void SeqStutterGenotyper::get_alleles(std::string& chrom_seq, std::vector<std::string>& alleles){
  assert(hap_blocks_.size() == 3 && alleles.size() == 0);
  HapBlock* block = hap_blocks_[1];
  int32_t start   = block->start(), end = block->end();

  std::string left_flank  = (start >= region_->start() ? uppercase(chrom_seq.substr(region_->start(), start-region_->start())) : "");
  std::string right_flank = (end <= region_->stop()    ? uppercase(chrom_seq.substr(end, region_->stop()-end)) : "");
  pos_ = std::min((int32_t)region_->start(), start);
  
  // If necessary, add 1bp on the left so that all the alleles match the reference sequence
  if (left_flank.empty()){
    bool pad_left = false;
    const std::string& ref_seq = block->get_seq(0);
    for (unsigned int i = 1; i < block->num_options(); i++){
      if (block->get_seq(i)[0] != ref_seq[0]){
	pad_left = true;
	break;
      }
    }

    if (pad_left){
      pos_ -= 1;
      left_flank = uppercase(chrom_seq.substr(pos_, 1));
    }
  }

  for (unsigned int i = 0; i < block->num_options(); i++){
    std::stringstream ss; 
    ss << left_flank << block->get_seq(i) << right_flank;
    alleles.push_back(ss.str());
  }

  pos_ += 1; // Fix off-by-1 VCF error
}
 
void SeqStutterGenotyper::debug_sample(int sample_index){
  std::cerr << "DEBUGGING SAMPLE..." << std::endl;
  std::cerr << "READ LL's:" << std::endl;
  double* read_LL_ptr = log_aln_probs_;
  int read_index = 0;
  for (unsigned int i = 0; i < num_reads_; ++i){
    if(sample_label_[i] == sample_index){
      std::cerr << "\t" << "READ #" << read_index << ", SEED BASE=" << seed_positions_[i] 
		<< ", TOTAL QUAL CORRECT= " << alns_[sample_index][read_index].sum_log_prob_correct(base_quality_) << ", " 
		<< bp_diffs_[i] << " " << max_index(read_LL_ptr, num_alleles_) << ", "
		<< log_p1_[read_index] << " " << log_p2_[read_index] <<  ", "
		<< alns_[sample_index][read_index].get_sequence().substr(0, seed_positions_[i]) 
		<< " " << alns_[sample_index][read_index].get_sequence().substr(seed_positions_[i]+1) << std::endl;
      for (unsigned int j = 0; j < num_alleles_; ++j, ++read_LL_ptr)
	std::cerr << "\t\t" << j << " " << *read_LL_ptr << std::endl;
      read_index++;
    }
    else
      read_LL_ptr += num_alleles_;
  }

  std::cerr << std::endl << "SAMPLE LL's:" << std::endl;
  double* sample_LL_ptr = log_sample_posteriors_ + sample_index;
  for (int index_1 = 0; index_1 < num_alleles_; ++index_1)
    for (int index_2 = 0; index_2 < num_alleles_; ++index_2){
      std::cerr << index_1 << " " << index_2 << " " << *sample_LL_ptr << "(" << exp(*sample_LL_ptr) << ")" << std::endl;
      sample_LL_ptr += num_samples_;
    }
  
  std::cerr << "END OF SAMPLE DEBUGGING..." << std::endl;
}


double SeqStutterGenotyper::calc_log_sample_posteriors(std::vector<int>& read_weights){
  assert(read_weights.size() == num_reads_);
  std::vector<double> sample_max_LLs(num_samples_, -DBL_MAX);
  double* sample_LL_ptr = log_sample_posteriors_;

  if (log_allele_priors_ != NULL)
    memcpy(log_sample_posteriors_, log_allele_priors_, num_alleles_*num_alleles_*num_samples_*sizeof(double));
  else {
    if (!haploid_){
      // Each genotype has an equal total prior, but heterozygotes have two possible phasings. Therefore,
      // i)   Phased heterozygotes have a prior of 1/(n(n+1))
      // ii)  Homozygotes have a prior of 2/(n(n+1))
      // iii) Total prior is n*2/(n(n+1)) + n(n-1)*1/(n(n+1)) = 2/(n+1) + (n-1)/(n+1) = 1

      // Set all elements to het prior
      double log_hetz_prior = -log(num_alleles_) - log(num_alleles_+1);
      std::fill(log_sample_posteriors_, log_sample_posteriors_+(num_alleles_*num_alleles_*num_samples_), log_hetz_prior);

      // Fix homozygotes
      double log_homoz_prior = log(2) - log(num_alleles_) - log(num_alleles_+1);
      for (unsigned int i = 0; i < num_alleles_; i++){
	double* LL_ptr = log_sample_posteriors_ + i*num_alleles_*num_samples_ + i*num_samples_;
	std::fill(LL_ptr, LL_ptr+num_samples_, log_homoz_prior);
      }
    }
    else {
      // Set all elements to impossible
      std::fill(log_sample_posteriors_, log_sample_posteriors_+(num_alleles_*num_alleles_*num_samples_), -DBL_MAX/2);

      // Fix homozygotes using a uniform prior
      double log_homoz_prior = -log(num_alleles_);
      for (unsigned int i = 0; i < num_alleles_; i++){
	double* LL_ptr = log_sample_posteriors_ + i*num_alleles_*num_samples_ + i*num_samples_;
	std::fill(LL_ptr, LL_ptr+num_samples_, log_homoz_prior);
      }
    }
  }

  for (int index_1 = 0; index_1 < num_alleles_; ++index_1){
    for (int index_2 = 0; index_2 < num_alleles_; ++index_2){
      double* read_LL_ptr = log_aln_probs_;
      for (int read_index = 0; read_index < num_reads_; ++read_index){
	sample_LL_ptr[sample_label_[read_index]] += read_weights[read_index]*log_sum_exp(LOG_ONE_HALF + log_p1_[read_index] + read_LL_ptr[index_1],
											 LOG_ONE_HALF + log_p2_[read_index] + read_LL_ptr[index_2]);
	assert(sample_LL_ptr[sample_label_[read_index]] <= TOLERANCE);
	read_LL_ptr += num_alleles_;
      }
      // Update the per-sample maximum LLs
      for (int sample_index = 0; sample_index < num_samples_; ++sample_index)
	sample_max_LLs[sample_index] = std::max(sample_max_LLs[sample_index], sample_LL_ptr[sample_index]);

      sample_LL_ptr += num_samples_;
    }
  }

  // Compute the normalizing factor for each sample using logsumexp trick
  std::fill(sample_total_LLs_, sample_total_LLs_ + num_samples_, 0.0);
  sample_LL_ptr = log_sample_posteriors_;
  for (int index_1 = 0; index_1 < num_alleles_; ++index_1)
    for (int index_2 = 0; index_2 < num_alleles_; ++index_2)
      for (int sample_index = 0; sample_index < num_samples_; ++sample_index, ++sample_LL_ptr)
	sample_total_LLs_[sample_index] += exp(*sample_LL_ptr - sample_max_LLs[sample_index]);
  for (int sample_index = 0; sample_index < num_samples_; ++sample_index){
    sample_total_LLs_[sample_index] = sample_max_LLs[sample_index] + log(sample_total_LLs_[sample_index]);
    assert(sample_total_LLs_[sample_index] <= TOLERANCE);
  }
  // Compute the total log-likelihood given the current parameters
  double total_LL = sum(sample_total_LLs_, sample_total_LLs_ + num_samples_);

  // Normalize each genotype LL to generate valid log posteriors
  sample_LL_ptr = log_sample_posteriors_;
  for (int index_1 = 0; index_1 < num_alleles_; ++index_1)
    for(int index_2 = 0; index_2 < num_alleles_; ++index_2)
      for (int sample_index = 0; sample_index < num_samples_; ++sample_index, ++sample_LL_ptr)
	*sample_LL_ptr -= sample_total_LLs_[sample_index];

  return total_LL;
}


double SeqStutterGenotyper::calc_log_sample_posteriors(){
  std::vector<int> weights(num_reads_, 1);
  return calc_log_sample_posteriors(weights);
}

bool SeqStutterGenotyper::use_read(Alignment& max_LL_aln, int num_flank_ins, int num_flank_del){
  return true;
  int32_t left_flank_bases  = haplotype_->get_block(1)->start() - max_LL_aln.get_start();
  int32_t right_flank_bases = max_LL_aln.get_stop() - haplotype_->get_block(2)->start() + 1;
  if (std::min(left_flank_bases, right_flank_bases) < 5)
    return false;
  else
    return true;

  /*
  if (max_LL_aln.num_matched_bases() < 25)
    return false;
  return true;
  */
}

void SeqStutterGenotyper::get_optimal_genotypes(std::vector< std::pair<int, int> >& gts){
  assert(gts.size() == 0);
  gts = std::vector< std::pair<int,int> > (num_samples_, std::pair<int,int>(-1,-1));
  double* log_post_ptr = log_sample_posteriors_;
  std::vector<double> log_phased_posteriors(num_samples_, -DBL_MAX);
  for (int index_1 = 0; index_1 < num_alleles_; ++index_1)
    for (int index_2 = 0; index_2 < num_alleles_; ++index_2)
      for (unsigned int sample_index = 0; sample_index < num_samples_; ++sample_index, ++log_post_ptr)
        if (*log_post_ptr > log_phased_posteriors[sample_index]){
          log_phased_posteriors[sample_index] = *log_post_ptr;
          gts[sample_index] = std::pair<int,int>(index_1, index_2);
        }
}

std::string SeqStutterGenotyper::condense_read_counts(std::vector<int>& read_diffs){
  if (read_diffs.size() == 0)
    return ".";
  std::map<int, int> diff_counts;
  for (unsigned int i = 0; i < read_diffs.size(); i++)
    diff_counts[read_diffs[i]]++;
  std::stringstream res;
  for (auto iter = diff_counts.begin(); iter != diff_counts.end(); iter++){
    if (iter != diff_counts.begin())
      res << ";";
    res << iter->first << "|" << iter->second;
  }
  return res.str();
}

void SeqStutterGenotyper::filter_alignments(std::ostream& logger, std::vector<int>& masked_reads){
  assert(masked_reads.size() == 0);
  masked_reads = std::vector<int>(num_samples_, 0);

  double filter_start = clock();
  std::vector< std::pair<int, int> > gts;
  get_optimal_genotypes(gts);

  int32_t filt_count = 0, keep_count = 0;
  std::vector<int> num_proc_alns(num_samples_, 0);
  HapAligner hap_aligner(haplotype_);
  double* read_LL_ptr = log_aln_probs_;
  for (unsigned int read_index = 0; read_index < num_reads_; read_index++){
    if (seed_positions_[read_index] < 0){
      masked_reads[sample_label_[read_index]]++;
      read_LL_ptr += num_alleles_;
      num_proc_alns[sample_label_[read_index]]++;
      continue;
    }

    int gt_a    = gts[sample_label_[read_index]].first;
    int gt_b    = gts[sample_label_[read_index]].second;
    int best_gt = ((LOG_ONE_HALF+log_p1_[read_index]+read_LL_ptr[gt_a] >  LOG_ONE_HALF+log_p2_[read_index]+read_LL_ptr[gt_b]) ? gt_a : gt_b);

    // Retrace alignment and ensure that it's of sufficient quality
    Alignment traced_aln;
    int idx_1 = sample_label_[read_index], idx_2 = num_proc_alns[sample_label_[read_index]], num_flank_ins, num_flank_del, stutter_size;
    hap_aligner.trace_optimal_aln(alns_[idx_1][idx_2], seed_positions_[read_index], best_gt, &base_quality_, traced_aln, num_flank_ins, num_flank_del, stutter_size);
    num_proc_alns[idx_1]++;

    // Zero out alignment probabilities for filtered reads
    if (!use_read(traced_aln, num_flank_ins, num_flank_del)){
      seed_positions_[read_index] = -2;
      for (unsigned int i = 0; i < haplotype_->num_combs(); ++i)
	read_LL_ptr[i] = 0;
      filt_count++;
    }
    else
      keep_count++;
    read_LL_ptr += num_alleles_;
  }
  calc_log_sample_posteriors();
  logger << "Filtered " << filt_count << " out of " << filt_count+keep_count << " reads based on their ML alignment tracebacks" << "\n";
  total_aln_filter_time_ = (clock() - filter_start)/CLOCKS_PER_SEC;
}

void SeqStutterGenotyper::write_vcf_record(std::vector<std::string>& sample_names, bool print_info, std::string& chrom_seq,
					   bool output_bootstrap_qualities, bool output_gls, bool output_pls,
					   bool output_allreads, bool output_pallreads, bool output_mallreads, bool output_viz,
					   bool visualize_left_alns, std::vector<int>& read_str_sizes,
					   std::ostream& html_output, std::ostream& out, std::ostream& logger){
  assert(haplotype_->num_blocks() == 3);
  assert(read_str_sizes.size() == 0);

  if(log_allele_priors_ != NULL)
    assert(!output_gls && !output_pls); // These fields only make sense in the context of MLE estimation, not MAP estimation

  // Compute the base pair differences from the reference
  std::vector<int> allele_bp_diffs;
  for (unsigned int i = 0; i < alleles_.size(); i++)
    allele_bp_diffs.push_back((int)alleles_[i].size() - (int)alleles_[0].size());
  
  // Filter reads with questionable alignments
  std::vector<int> masked_reads;
  if (true)
    filter_alignments(logger, masked_reads);

  // Extract each sample's posterior base pair dosage, MAP genotype, the associated phased genotype posterior
  // and the genotype likelihoods
  std::vector< std::pair<int,int> > gts(num_samples_, std::pair<int,int>(-1,-1));
  std::vector<double> log_phased_posteriors(num_samples_, -DBL_MAX), bp_dosages;
  std::vector<int> dip_bpdiffs;
  std::vector< std::vector<double> > log_post_probs(num_samples_), gls(num_samples_);
  std::vector< std::vector<int> > pls(num_samples_);
  double* log_post_ptr = log_sample_posteriors_;
  for (int index_1 = 0; index_1 < num_alleles_; ++index_1){
    for (int index_2 = 0; index_2 < num_alleles_; ++index_2){
      dip_bpdiffs.push_back(allele_bp_diffs[index_1]+allele_bp_diffs[index_2]);
      for (unsigned int sample_index = 0; sample_index < num_samples_; ++sample_index, ++log_post_ptr){
	if (*log_post_ptr > log_phased_posteriors[sample_index]){
	  log_phased_posteriors[sample_index] = *log_post_ptr;
	  gts[sample_index] = std::pair<int,int>(index_1, index_2);
	}
	log_post_probs[sample_index].push_back(*log_post_ptr);
	if (index_2 <= index_1){
	  double gl_base_e =  sample_total_LLs_[sample_index] + LOG_ONE_HALF
	    + log_sum_exp(*log_post_ptr, log_sample_posteriors_[index_2*num_alleles_*num_samples_ + index_1*num_samples_ + sample_index]);

	  if (!haploid_ || (index_1 == index_2))
	    gls[sample_index].push_back(gl_base_e*LOG_E_BASE_10); // Convert from ln to log10
	}
      }
    }
  }

  std::vector<double> gl_diffs;
  for (unsigned int sample_index = 0; sample_index < num_samples_; sample_index++){
    bp_dosages.push_back((!haploid_ ? 1.0 : 0.5)*expected_value(log_post_probs[sample_index], dip_bpdiffs));
    double max_gl    = *(std::max_element(gls[sample_index].begin(), gls[sample_index].end()));
    double second_gl = -DBL_MAX;
    for (unsigned int j = 0; j < gls[sample_index].size(); j++){
      pls[sample_index].push_back((int)(gls[sample_index][j]-max_gl));
      if (gls[sample_index][j] < max_gl)
	second_gl = std::max(second_gl, gls[sample_index][j]);
    }

    if (num_alleles_ == 1)
      gl_diffs.push_back(-1000);
    else {
      int gl_index;
      if (haploid_)
	gl_index = gts[sample_index].first;
      else {
	int min_gt = std::min(gts[sample_index].first, gts[sample_index].second);
	int max_gt = std::max(gts[sample_index].first, gts[sample_index].second);
	gl_index   = max_gt*(max_gt+1)/2 + min_gt;
      }
      gl_diffs.push_back((abs(max_gl-gls[sample_index][gl_index]) < TOLERANCE) ? (max_gl-second_gl) : gls[sample_index][gl_index]-max_gl);
    }
  }

  // Extract the genotype phasing probability conditioned on the determined sample genotypes
  std::vector<double> log_unphased_posteriors, phase_probs;
  for (unsigned int sample_index = 0; sample_index < num_samples_; sample_index++){
    int gt_a = gts[sample_index].first, gt_b = gts[sample_index].second;
    if (gt_a == gt_b){
      log_unphased_posteriors.push_back(log_phased_posteriors[sample_index]);
      phase_probs.push_back(1.0);
    }
    else {
      double log_p1  = log_phased_posteriors[sample_index];
      double log_p2  = log_sample_posteriors_[gt_b*num_alleles_*num_samples_ + gt_a*num_samples_ + sample_index];
      double log_tot = log_sum_exp(log_p1, log_p2);
      log_unphased_posteriors.push_back(log_tot);
      phase_probs.push_back(exp(log_p1-log_tot));
    }
  }

  // Extract information about each read and group by sample
  assert(bp_diffs_.size() == num_reads_);
  std::vector<int> num_aligned_reads(num_samples_, 0), num_reads_with_snps(num_samples_, 0), num_proc_alns(num_samples_, 0);
  std::vector<int> num_reads_with_stutter(num_samples_, 0), num_reads_with_flank_indels(num_samples_, 0);
  std::vector<int> num_reads_strand_one(num_samples_, 0), num_reads_strand_two(num_samples_, 0);
  std::vector< std::vector<int> > bps_per_sample(num_samples_), ml_bps_per_sample(num_samples_);
  std::vector< std::vector<double> > log_read_phases(num_samples_), posterior_bps_per_sample(num_samples_);

  assert(max_LL_alns_.size() == 0);
  max_LL_alns_ = std::vector< std::vector<Alignment> >(num_samples_);
  HapAligner hap_aligner(haplotype_);
  double* read_LL_ptr   = log_aln_probs_;
  for (unsigned int read_index = 0; read_index < num_reads_; read_index++){
    if (seed_positions_[read_index] < 0){
      read_LL_ptr += num_alleles_;
      num_proc_alns[sample_label_[read_index]]++;
      read_str_sizes.push_back(-999);
      continue;
    }

    // Extract read's phase posterior conditioned on the determined sample genotype
    int gt_a = gts[sample_label_[read_index]].first;
    int gt_b = gts[sample_label_[read_index]].second;
    double total_read_LL = log_sum_exp(LOG_ONE_HALF+log_p1_[read_index]+read_LL_ptr[gt_a], LOG_ONE_HALF+log_p2_[read_index]+read_LL_ptr[gt_b]);
    double log_phase_one = LOG_ONE_HALF + log_p1_[read_index] + read_LL_ptr[gt_a] - total_read_LL; 
    log_read_phases[sample_label_[read_index]].push_back(log_phase_one);

    // Retrace alignment and ensure that it's of sufficient quality
    double trace_start = clock();
    int best_gt = (log_phase_one > LOG_ONE_HALF ? gt_a : gt_b);
    Alignment traced_aln;
    int idx_1 = sample_label_[read_index], idx_2 = num_proc_alns[sample_label_[read_index]], num_flank_ins, num_flank_del, stutter_size;
    hap_aligner.trace_optimal_aln(alns_[idx_1][idx_2], seed_positions_[read_index], best_gt, &base_quality_, traced_aln, num_flank_ins, num_flank_del, stutter_size);
    num_proc_alns[idx_1]++;

    if (stutter_size != 0)
      num_reads_with_stutter[sample_label_[read_index]]++;
    if (num_flank_ins != 0 || num_flank_del != 0)
      num_reads_with_flank_indels[sample_label_[read_index]]++;
    read_str_sizes.push_back(allele_bp_diffs[best_gt]+stutter_size);

    if (visualize_left_alns)
      max_LL_alns_[idx_1].push_back(alns_[idx_1][idx_2]);
    max_LL_alns_[idx_1].push_back(traced_aln);
    total_aln_trace_time_ += (clock() - trace_start)/CLOCKS_PER_SEC;

    // Adjust number of aligned reads per sample
    num_aligned_reads[sample_label_[read_index]]++;

    // Adjust number of reads with SNP information for each sample
    if (abs(log_p1_[read_index] - log_p2_[read_index]) > TOLERANCE){
      num_reads_with_snps[sample_label_[read_index]]++;
      if (log_p1_[read_index] > log_p2_[read_index])
	num_reads_strand_one[sample_label_[read_index]]++;
      else
	num_reads_strand_two[sample_label_[read_index]]++;
    }

    // Extract the bp difference observed in read from left-alignment
    bps_per_sample[sample_label_[read_index]].push_back(bp_diffs_[read_index]);

    // Extract the posterior bp differences observed in read from haplotype alignment
    posterior_bps_per_sample[sample_label_[read_index]].push_back(expected_value(read_LL_ptr, allele_bp_diffs));

    // Extract the ML bp difference observed in read based on the ML genotype
    ml_bps_per_sample[sample_label_[read_index]].push_back(read_str_sizes.back());

    read_LL_ptr += num_alleles_;
  }

  // Compute bootstrap qualities if flag set
  std::vector<double> bootstrap_qualities;
  int bootstrap_iter = 100;
  if (output_bootstrap_qualities)
    compute_bootstrap_qualities(bootstrap_iter, bootstrap_qualities);
 
  // Compute allele counts for samples of interest
  std::set<std::string> samples_of_interest(sample_names.begin(), sample_names.end());
  std::vector<int> allele_counts(num_alleles_);
  int sample_index = 0, skip_count = 0, filt_count = 0;
  for (auto gt_iter = gts.begin(); gt_iter != gts.end(); ++gt_iter, ++sample_index){
    if (samples_of_interest.find(sample_names_[sample_index]) == samples_of_interest.end())
      continue;
    if (require_one_read_ && num_aligned_reads[sample_index] == 0)
      continue;
    if (call_sample_[sample_index]) {
      if (haploid_){
	assert(gt_iter->first == gt_iter->second);
	allele_counts[gt_iter->first]++;
      }
      else {
	allele_counts[gt_iter->first]++;
	allele_counts[gt_iter->second]++;
      }
    }
    else
      skip_count++;
  }

  if (print_info){
    logger << "Allele counts" << std::endl;
    for (unsigned int i = 0; i < alleles_.size(); i++){
      bool expanded = (expanded_alleles_.find(haplotype_->get_block(1)->get_seq(i)) != expanded_alleles_.end());
      logger << alleles_[i] << " " << allele_counts[i] << " " << expanded << std::endl;
    }
    logger << std::endl;
  }

  //VCF line format = CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE_1 SAMPLE_2 ... SAMPLE_N
  out << region_->chrom() << "\t" << pos_ << "\t" << (region_->name().empty() ? "." : region_->name());

  // Add reference allele and alternate alleles
  out << "\t" << alleles_[0] << "\t";
  if (num_alleles_ == 1)
    out << ".";
  else {
    for (int i = 1; i < num_alleles_-1; i++)
      out << alleles_[i] << ",";
    out << alleles_[num_alleles_-1];
  }

  // Add QUAL and FILTER fields
  out << "\t" << "." << "\t" << ".";

  // Add INFO field items
  out << "\tINFRAME_PGEOM=" << stutter_model_->get_parameter(true,  'P') << ";" 
      << "INFRAME_UP="      << stutter_model_->get_parameter(true,  'U') << ";" 
      << "INFRAME_DOWN="    << stutter_model_->get_parameter(true,  'D') << ";" 
      << "OUTFRAME_PGEOM="  << stutter_model_->get_parameter(false, 'P') << ";" 
      << "OUTFRAME_UP="     << stutter_model_->get_parameter(false, 'U') << ";" 
      << "OUTFRAME_DOWN="   << stutter_model_->get_parameter(false, 'D') << ";"
      << "START="           << region_->start()  << ";"
      << "END="             << region_->stop()   << ";"
      << "PERIOD="          << region_->period() << ";"
      << "NSKIP="           << skip_count        << ";"
      << "NFILT="           << filt_count        << ";";
  if (num_alleles_ > 1){
    out << "BPDIFFS=" << allele_bp_diffs[1];
    for (unsigned int i = 2; i < num_alleles_; i++)
      out << "," << allele_bp_diffs[i];
    out << ";";
  }
  // Add allele counts
  out << "REFAC=" << allele_counts[0] << ";";
  if (allele_counts.size() > 1){
    out << "AC=";
    for (unsigned int i = 1; i < allele_counts.size()-1; i++)
      out << allele_counts[i] << ",";
    out << allele_counts.back() << ";";
  }

  // Add FORMAT field
  out << (!haploid_ ? "\tGT:GB:Q:PQ:DP:DSNP:DFILT:DSTUTTER:DFLANKINDEL:PDP:PSNP:BPDOSE:GLDIFF" : "\tGT:GB:Q:DP:DFILT:DSTUTTER:DFLANKINDEL:BPDOSE:GLDIFF");
  if (output_bootstrap_qualities) out << ":BQ";
  if (output_allreads)            out << ":ALLREADS";
  if (output_pallreads)           out << ":PALLREADS";
  if (output_mallreads)           out << ":MALLREADS";
  if (output_gls)                 out << ":GL";
  if (output_pls)                 out << ":PL";

  std::map<std::string, std::string> sample_results;
  for (unsigned int i = 0; i < sample_names.size(); i++){
    out << "\t";
    auto sample_iter = sample_indices_.find(sample_names[i]);
    if (sample_iter == sample_indices_.end()){
      out << ".";
      continue;
    }
    
    // Don't report information for a sample if none of its reads were successfully realigned
    // and we require at least one read
    if (require_one_read_ && num_aligned_reads[sample_iter->second] == 0){
      out << ".";
      continue;
    }

    // Don't report information for a sample if flag has been set to false
    if (!call_sample_[sample_iter->second]){
      out << ".";
      continue;
    }
    
    int sample_index    = sample_iter->second;
    double phase1_reads = (num_aligned_reads[sample_index] == 0 ? 0 : exp(log_sum_exp(log_read_phases[sample_index])));
    double phase2_reads = num_aligned_reads[sample_index] - phase1_reads;

    std::stringstream samp_info;
    samp_info << allele_bp_diffs[gts[sample_index].first] << "|" << allele_bp_diffs[gts[sample_index].second];
    sample_results[sample_names[i]] = samp_info.str();

    if (!haploid_){
      out << gts[sample_index].first << "|" << gts[sample_index].second                             // Genotype
	  << ":" << allele_bp_diffs[gts[sample_index].first]
	  << "|" << allele_bp_diffs[gts[sample_index].second]                                       // Base pair differences from reference
	  << ":" << exp(log_unphased_posteriors[sample_index])                                      // Unphased posterior
	  << ":" << exp(log_phased_posteriors[sample_index])                                        // Phased posterior
	  << ":" << num_aligned_reads[sample_index]                                                 // Total reads used to genotype (after filtering)
	  << ":" << num_reads_with_snps[sample_index]                                               // Total reads with SNP information
	  << ":" << masked_reads[sample_index]                                                      // Total masked reads
	  << ":" << num_reads_with_stutter[sample_index]                                            // Total reads with a non-zero stutter artifact in ML alignment
	  << ":" << num_reads_with_flank_indels[sample_index]                                       // Total reads with an indel in flank in ML alignment
	  << ":" << phase1_reads << "|" << phase2_reads                                             // Reads per allele
	  << ":" << num_reads_strand_one[sample_index] << "|" << num_reads_strand_two[sample_index] // Reads with SNPs supporting each haploid genotype
	  << ":" << bp_dosages[sample_index];                                                       // Posterior STR dosage (in base pairs)

      // Difference in GL between the current and next best genotype
      if (num_alleles_ == 1)
	out << ":" << ".";
      else
	out << ":" << gl_diffs[sample_index];
    }
    else {
      out << gts[sample_index].first                                                                // Genotype
	  << ":" << allele_bp_diffs[gts[sample_index].first]                                        // Base pair differences from reference
	  << ":" << exp(log_unphased_posteriors[sample_index])                                      // Unphased posterior
	  << ":" << num_aligned_reads[sample_index]                                                 // Total reads used to genotype (after filtering)
	  << ":" << masked_reads[sample_index]                                                      // Total masked reads
	  << ":" << num_reads_with_stutter[sample_index]                                            // Total reads with a non-zero stutter artifact in ML alignment
	  << ":" << num_reads_with_flank_indels[sample_index]                                       // Total reads with an indel in flank in ML alignment
	  << ":" << bp_dosages[sample_index];                                                       // Posterior STR dosage (in base pairs)

      // Difference in GL between the current and next best genotype
      if (num_alleles_ == 1)
	out << ":" << ".";
      else
	out << ":" << gl_diffs[sample_index];
    }

    if (output_bootstrap_qualities)
      out << ":" << bootstrap_qualities[sample_index];

    // Add bp diffs from regular left-alignment
    if (output_allreads){
      if (condense_read_count_fields)
	out << ":" << condense_read_counts(bps_per_sample[sample_index]);
      else {
	if (bps_per_sample[sample_index].size() != 0){
	  out << ":" << bps_per_sample[sample_index][0];
	  for (unsigned int j = 1; j < bps_per_sample[sample_index].size(); j++)
	    out << "," << bps_per_sample[sample_index][j];
	}
	else
	  out << ":" << ".";
      }
    }

    // Expected base pair differences from alignment probabilities
    if (output_pallreads){
      if (posterior_bps_per_sample[sample_index].size() != 0){
	out << ":" << posterior_bps_per_sample[sample_index][0];
	for (unsigned int j = 1; j < posterior_bps_per_sample[sample_index].size(); j++)
	  out << "," << posterior_bps_per_sample[sample_index][j];
      }
      else
	out << ":" << ".";
    }

    // Maximum likelihood base pair differences in each read from alignment probabilites
    if (output_mallreads){
      if (condense_read_count_fields)
	out << ":" << condense_read_counts(ml_bps_per_sample[sample_index]);
      else {
	if (ml_bps_per_sample[sample_index].size() != 0){
	  out << ":" << ml_bps_per_sample[sample_index][0];
	  for (unsigned int j = 1; j < ml_bps_per_sample[sample_index].size(); j++)
	    out << "," << ml_bps_per_sample[sample_index][j];
	}
	else
	  out << ":" << ".";
      }
    }

    // Genotype and phred-scaled likelihoods
    if (output_gls){
      out << ":" << gls[sample_index][0];
      for (unsigned int j = 1; j < gls[sample_index].size(); j++)
	out << "," << gls[sample_index][j];
    }
    if (output_pls){
      out << ":" << pls[sample_index][0];
      for (unsigned int j = 1; j < pls[sample_index].size(); j++)
	out << "," << pls[sample_index][j];
    }
  }
  out << "\n";

  // Render HTML of Smith-Waterman alignments (or haplotype alignments)
  if (output_viz){
    std::stringstream locus_info;
    locus_info << region_->chrom() << "\t" << region_->start() << "\t" << region_->stop();
    visualizeAlignments(max_LL_alns_, sample_names_, sample_results, hap_blocks_, chrom_seq, locus_info.str(), true, html_output);
  }
}


bool SeqStutterGenotyper::recompute_stutter_model(std::string& chrom_seq, std::ostream& logger,
						  int max_em_iter, double abs_ll_converge, double frac_ll_converge){
  logger << "Retraining EM stutter genotyper using maximum likelihood alignments" << std::endl;

  // Get the artifact sizes observed in each read
  std::vector<std::string> empty_sample_names;
  std::vector<int> read_str_sizes;
  write_vcf_record(empty_sample_names, false, chrom_seq, false, false, false, false, false, false, false, false, read_str_sizes, std::cerr, std::cerr, logger);
  max_LL_alns_.clear(); // Need to clear this data structure for a future call to write_vcf_record to work
  assert(read_str_sizes.size() == num_reads_);

  std::vector< std::vector<int> > str_num_bps(num_samples_);
  std::vector< std::vector<double> > str_log_p1s(num_samples_), str_log_p2s(num_samples_);
  for (unsigned int read_index = 0; read_index < num_reads_; read_index++){
    if (read_str_sizes[read_index] != -999){
      str_num_bps[sample_label_[read_index]].push_back(read_str_sizes[read_index]);
      str_log_p1s[sample_label_[read_index]].push_back(log_p1_[read_index]);
      str_log_p2s[sample_label_[read_index]].push_back(log_p2_[read_index]);
    }
  }

  EMStutterGenotyper length_genotyper(region_->chrom(), region_->start(), region_->stop(), haploid_, str_num_bps, str_log_p1s, str_log_p2s, sample_names_, region_->period(), 0);
  bool trained = length_genotyper.train(max_em_iter, abs_ll_converge, frac_ll_converge, false, logger);
  if (trained){
    delete stutter_model_;
    stutter_model_ = length_genotyper.get_stutter_model()->copy();
    logger << "Learned stutter model: " << *stutter_model_ << std::endl;

    // Replace the stutter model in the repeat block
    assert(haplotype_->num_blocks() == 3);
    assert(haplotype_->get_block(1)->get_repeat_info() != NULL);
    ((RepeatBlock*)(haplotype_->get_block(1)))->get_repeat_info()->set_stutter_model(stutter_model_);
    return genotype(logger);
  }
  else {
    logger << "Retraining stutter model training failed for locus " << region_->chrom() << ":" << region_->start() << "-" << region_->stop() << std::endl;
    return false;
  }
}


void SeqStutterGenotyper::compute_bootstrap_qualities(int num_iter, std::vector<double>& bootstrap_qualities){
  assert(bootstrap_qualities.size() == 0);
  double bootstrap_start = clock();

  // Extract the original ML genotypes
  std::vector< std::pair<int, int> > ML_gts;
  get_optimal_genotypes(ML_gts);

  // Partition the aligned reads by sample
  std::vector< std::vector<int> > reads_by_sample(num_samples_);
  for (unsigned int i = 0; i < num_reads_; i++)
    if (seed_positions_[i] >= 0)
      reads_by_sample.at(sample_label_[i]).push_back(i);

  std::vector<int> ML_gt_counts(num_samples_, 0);
  std::uniform_int_distribution<int> unif_dist;
  std::default_random_engine gen;
  for (unsigned int i = 0; i < num_iter; i++){
    std::vector<int> bootstrap_weights(num_reads_, 0);

    // Bootstrap reads for each sample
    for (unsigned int j = 0; j < num_samples_; j++){
      int mod = reads_by_sample[j].size();
      for (unsigned int k = 0; k < mod; k++)
	bootstrap_weights[reads_by_sample[j][unif_dist(gen)% mod]]++;
    }

    // Recompute the posteriors using bootstrapped read weights
    calc_log_sample_posteriors(bootstrap_weights);

    // Increment count if bootstrapped ML genotype (unordered) matches the ML genotype
    std::vector< std::pair<int, int> > bootstrap_gts;
    get_optimal_genotypes(bootstrap_gts);
    for (unsigned int i = 0; i < num_samples_; i++){
      if (bootstrap_gts[i].first == ML_gts[i].first && bootstrap_gts[i].second == ML_gts[i].second)
	ML_gt_counts[i]++;
      else if (bootstrap_gts[i].first == ML_gts[i].second && bootstrap_gts[i].second == ML_gts[i].first)
	ML_gt_counts[i]++;
    }
  }

  // Recalculate the posteriors using the original read indices to restore everything to normalcy
  calc_log_sample_posteriors();

  // Compute the boostrapped qualities as the fraction of iterations in which the genotype matched
  for (unsigned int i = 0; i < num_samples_; i++)
    bootstrap_qualities.push_back(1.0*ML_gt_counts[i]/num_iter);

  double bootstrap_time = (clock() - bootstrap_start)/CLOCKS_PER_SEC;
  //std::cerr << "Bootstrapping time = " << bootstrap_time << std::endl;
}
