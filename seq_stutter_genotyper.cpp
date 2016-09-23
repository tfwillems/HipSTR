#include <algorithm>
#include <cfloat>
#include <cstring>
#include <random>
#include <string>
#include <sstream>
#include <time.h>

#include "seq_stutter_genotyper.h"
#include "bam_processor.h"
#include "em_stutter_genotyper.h"
#include "error.h"
#include "extract_indels.h"
#include "mathops.h"
#include "stringops.h"
#include "vcf_input.h"

#include "SeqAlignment/AlignmentData.h"
#include "SeqAlignment/AlignmentModel.h"
#include "SeqAlignment/AlignmentViz.h"
#include "SeqAlignment/HaplotypeGenerator.h"
#include "SeqAlignment/HapAligner.h"
#include "SeqAlignment/RepeatStutterInfo.h"
#include "SeqAlignment/RepeatBlock.h"

#include "cephes/cephes.h"

int max_index(double* vals, unsigned int num_vals){
  int best_index = 0;
  for (unsigned int i = 1; i < num_vals; i++)
    if (vals[i] > vals[best_index])
      best_index = i;
  return best_index;
}

void SeqStutterGenotyper::haps_to_alleles(int hap_block_index, std::vector<int>& allele_indices){
  assert(allele_indices.empty());
  allele_indices.reserve(haplotype_->num_combs());
  haplotype_->reset();
  do {
    allele_indices.push_back(haplotype_->cur_index(hap_block_index));
  } while (haplotype_->next());
  haplotype_->reset();
}

void SeqStutterGenotyper::get_unspanned_alleles(std::vector<int>& allele_indices, std::ostream& logger){
  assert(allele_indices.size() == 0);

  // Extract each sample's MAP genotype
  std::vector< std::pair<int,int> > haps;
  get_optimal_haplotypes(log_sample_posteriors_, haps);

  std::vector<AlignmentTrace*> traced_alns;
  retrace_alignments(logger, traced_alns);

  assert(haplotype_->num_blocks() == 3);
  int str_block_index = 1;
  HapBlock* str_block = haplotype_->get_block(str_block_index);

  std::vector<bool> spanned(num_alleles_, false);
  spanned[0] = true;
  double* read_LL_ptr = log_aln_probs_;
  for (unsigned int read_index = 0; read_index < num_reads_; read_index++){
    if (seed_positions_[read_index] < 0){
      read_LL_ptr += num_alleles_;
      continue;
    }
    AlignmentTrace* trace = traced_alns[read_index];
    if (trace->traced_aln().get_start() < str_block->start())
      if (trace->traced_aln().get_stop() > str_block->end())
	if (trace->stutter_size(str_block_index) == 0){
	  int hap_a    = haps[sample_label_[read_index]].first;
	  int hap_b    = haps[sample_label_[read_index]].second;
	  int best_hap = hap_a;
	  if (!haploid_ && (hap_a != hap_b)){
	    double v1 = log_p1_[read_index]+read_LL_ptr[hap_a], v2 = log_p2_[read_index]+read_LL_ptr[hap_b];
	    if (abs(v1-v2) > TOLERANCE)
	      best_hap = (v1 > v2 ? hap_a : hap_b);
	  }
	  spanned[best_hap] = true;
	}
    read_LL_ptr += num_alleles_;
  }

  int count = 0;
  for (unsigned int i = 0; i < num_alleles_; ++i)
    if (!spanned[i]){
      allele_indices.push_back(i);
      count++;
    }
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
  get_optimal_haplotypes(log_sample_posteriors_, gts);

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
  std::vector<int> allele_mapping;
  int keep_count = 0;
  for (unsigned int i = 0; i < alleles_.size(); i++){
    if (keep_allele[i]){
      fixed_alleles.push_back(alleles_[i]);
      allele_mapping.push_back(keep_count++);
    }
    else
      allele_mapping.push_back(-1);
  }
  
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

  // Fix alignment traceback cache (as allele indices have changed)
  std::map<std::pair<int,int>, AlignmentTrace*> new_trace_cache;
  for (auto cache_iter = trace_cache_.begin(); cache_iter != trace_cache_.end(); cache_iter++){
    int new_allele_index = allele_mapping[cache_iter->first.second];
    if (new_allele_index != -1)
      new_trace_cache[std::pair<int,int>(cache_iter->first.first, new_allele_index)] = cache_iter->second;
  }
  trace_cache_ = new_trace_cache;

  // Resize and recalculate genotype posterior array
  delete [] log_sample_posteriors_;
  log_sample_posteriors_ = new double[fixed_num_alleles*fixed_num_alleles*num_samples_];
  calc_log_sample_posteriors();
}

void SeqStutterGenotyper::init(StutterModel& stutter_model, std::string& chrom_seq, std::ostream& logger){
  // Allocate and initiate additional data structures
  read_weights_.clear();
  pool_index_   = new int[num_reads_];
  second_mate_  = new bool[num_reads_];
  int32_t min_start = INT_MAX, max_stop = INT_MIN;
  std::string prev_aln_name = "";
  for (unsigned int read_index = 0; read_index < num_reads_; read_index++){
    if (use_for_haps_[read_index]){
      min_start = std::min(min_start, alns_[read_index].get_start());
      max_stop  = std::max(max_stop,  alns_[read_index].get_stop());
    }

    pool_index_[read_index]   = (pool_identical_seqs_ ? pooler_.add_alignment(alns_[read_index]) : read_index);
    second_mate_[read_index]  = (alns_[read_index].get_name().compare(prev_aln_name) == 0);
    read_weights_.push_back(second_mate_[read_index] ? 0 : 1);
    prev_aln_name = alns_[read_index].get_name();
  }

  double locus_hap_build_time = clock();
  std::vector<std::string> vcf_alleles;
  if (min_start >= region_->start()-5 || max_stop < region_->stop()+5){
    // No reads extend 5bp upstream and downstream of the STR
    logger << "Skipping region as no reads extend +- 5bp from the STR boundary" << std::endl;
    pos_ = -1;
  }
  else if (ref_vcf_ != NULL){
    bool success = false;
    std::vector<bool> got_priors;

    if (true){
      // Read alleles from VCF
      logger << "Reading STR alleles from VCF" << std::endl;
      read_vcf_alleles(ref_vcf_, region_, alleles_, pos_, success);
      assert(log_allele_priors_ == NULL);
    }
    else {
      // Read alleles and priors for each sample's genotypes from VCF generated by PhasedBEAGLE
      logger << "Reading STR alleles and priors from VCF" << std::endl;
      log_allele_priors_ = extract_vcf_alleles_and_log_priors(ref_vcf_, region_, sample_indices_, alleles_, got_priors, pos_, success, logger);
      assert(got_priors.size() == num_samples_);
    }

    num_alleles_ = alleles_.size();
    if (success){
      assert(num_alleles_ > 0);

      // Construct the haplotype using the set of VCF alleles
      haplotype_ = generate_haplotype(pos_, *region_, MAX_REF_FLANK_LEN, chrom_seq, alleles_, &stutter_model, hap_blocks_, logger);
      
      // If priors in the VCF, don't call samples without allele priors
      call_sample_ = std::vector<bool>(num_samples_, true);
      if (log_allele_priors_ != NULL)
	for (unsigned int i = 0; i < num_samples_; i++)
	  call_sample_[i] = (call_sample_[i] && got_priors[i]);
    }
    else
      pos_ = -1;
  }
  else {
    // Generate putative haplotypes and determine the number of alleles
    logger << "Generating putative haplotypes..." << std::endl;

    // Select only those alignments marked as good for haplotype generation
    std::vector<AlnList> gen_hap_alns(num_samples_);
    for (unsigned int read_index = 0; read_index < num_reads_; read_index++)
      if (use_for_haps_[read_index])
	gen_hap_alns[sample_label_[read_index]].push_back(alns_[read_index]);

    haplotype_   = generate_haplotype(*region_, MAX_REF_FLANK_LEN, chrom_seq, gen_hap_alns, vcf_alleles, &stutter_model,
				      alleles_from_bams_, hap_blocks_, call_sample_, logger);
    call_sample_ =  std::vector<bool>(num_samples_, true);
    num_alleles_ = haplotype_->num_combs();
    assert(call_sample_.size() == num_samples_);

    // Extract full STR sequence for each allele using annotated repeat region and the haplotype above
    get_alleles(chrom_seq, alleles_);
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

void SeqStutterGenotyper::calc_hap_aln_probs(Haplotype* haplotype, double* log_aln_probs, int* seed_positions){
  double locus_hap_aln_time = clock();
  HapAligner hap_aligner(haplotype);
  int num_alleles = haplotype->num_combs();

  if (pool_identical_seqs_){
    // Align each pooled read to each haplotype
    AlnList& pooled_alns = pooler_.get_alignments();
    double* log_pool_aln_probs = new double[pooled_alns.size()*num_alleles];
    int* pool_seed_positions   = new int[pooled_alns.size()];
    hap_aligner.process_reads(pooled_alns, 0, &base_quality_, log_pool_aln_probs, pool_seed_positions);

    // Copy each pool's alignment probabilities to the entries for its constituent reads
    double* log_aln_ptr = log_aln_probs;
    for (unsigned int i = 0; i < num_reads_; i++){
      seed_positions[i] = pool_seed_positions[pool_index_[i]];
      std::memcpy(log_aln_ptr, log_pool_aln_probs + num_alleles*pool_index_[i], num_alleles*sizeof(double));
      log_aln_ptr += num_alleles;
    }

    delete [] log_pool_aln_probs;
    delete [] pool_seed_positions;
  }
  else {
    // Align each read against each candidate haplotype
    int read_index = 0;
    hap_aligner.process_reads(alns_, read_index, &base_quality_, log_aln_probs, seed_positions);
  }

  // If both mate pairs overlap the STR region, they share the same phasing probabilities
  // We therefore need to avoid treating them as independent reads
  // To do so, we combine the alignment probabilities here and set the read weight
  // for the second in the pair to zero during posterior calculation
  for (unsigned int i = 0; i < num_reads_; ++i){
    if (!second_mate_[i])
      continue;
    double* mate_one_ptr = log_aln_probs + (i-1)*num_alleles;
    double* mate_two_ptr = log_aln_probs + i*num_alleles;
    for (unsigned int j = 0; j < num_alleles; ++j, ++mate_one_ptr, ++mate_two_ptr){
      double total  = *mate_one_ptr + *mate_two_ptr;
      *mate_one_ptr = total;
      *mate_two_ptr = total;
    }
  }

  locus_hap_aln_time   = (clock() - locus_hap_aln_time)/CLOCKS_PER_SEC;
  total_hap_aln_time_ += locus_hap_aln_time;
}

bool SeqStutterGenotyper::id_and_align_to_stutter_alleles(std::string& chrom_seq, std::ostream& logger){
  assert(haplotype_->num_blocks() == 3);
  assert(hap_blocks_[1]->get_repeat_info() != NULL);

  // Look for candidate alleles present in stutter artifacts
  std::vector<std::string> stutter_seqs;
  get_stutter_candidate_alleles(logger, stutter_seqs);
  while (stutter_seqs.size() != 0){
    std::sort(stutter_seqs.begin(), stutter_seqs.end(), stringLengthLT);
    RepeatBlock* rep_block = (RepeatBlock *)haplotype_->get_block(1);
    if (stutter_seqs[0].size() < std::abs(rep_block->get_repeat_info()->max_deletion()))
      return false;

    // Construct a new haplotype containing only stutter alleles and align each read to it
    std::vector<HapBlock*> blocks;
    blocks.push_back(hap_blocks_[0]);
    blocks.push_back(new RepeatBlock(hap_blocks_[1]->start(), hap_blocks_[1]->end(), stutter_seqs[0], region_->period(),
				     hap_blocks_[1]->get_repeat_info()->get_stutter_model()));
    blocks.push_back(hap_blocks_[2]);
    for (unsigned int i = 1; i < stutter_seqs.size(); i++)
      blocks[1]->add_alternate(stutter_seqs[i]);
    Haplotype* haplotype      = new Haplotype(blocks);
    double* new_log_aln_probs = new double[num_reads_*stutter_seqs.size()];
    calc_hap_aln_probs(haplotype, new_log_aln_probs, seed_positions_);
    delete blocks[1];
    delete haplotype;

    // Create a new sorted list of alleles and an STR haplotype block with all alleles
    std::vector<std::string> str_seqs;
    for (unsigned int i = 0; i < haplotype_->get_block(1)->num_options(); i++)
      str_seqs.push_back(haplotype_->get_block(1)->get_seq(i));
    for (unsigned int i = 0; i < stutter_seqs.size(); i++)
      str_seqs.push_back(stutter_seqs[i]);
    std::sort(str_seqs.begin()+1, str_seqs.end(), stringLengthLT);
    HapBlock* str_block = new RepeatBlock(hap_blocks_[1]->start(), hap_blocks_[1]->end(), hap_blocks_[1]->get_seq(0), region_->period(),
					  hap_blocks_[1]->get_repeat_info()->get_stutter_model());
    for (unsigned int i = 1; i < str_seqs.size(); i++)
      str_block->add_alternate(str_seqs[i]);

    // Determine the mapping from each allele to its new index
    std::vector<int> original_indices, stutter_indices;
    for (unsigned int i = 0; i < num_alleles_; i++)
      original_indices.push_back(str_block->index_of(hap_blocks_[1]->get_seq(i)));
    for (unsigned int i = 0; i < stutter_seqs.size(); i++)
      stutter_indices.push_back(str_block->index_of(stutter_seqs[i]));

    // Combine alignment probabilities by copying them to their new indices
    int total_alleles            = num_alleles_ + stutter_seqs.size();
    double* fixed_log_aln_probs  = new double[total_alleles*num_reads_];
    double* log_aln_ptr_original = log_aln_probs_;
    double* log_aln_ptr_stutter  = new_log_aln_probs;
    double* log_aln_ptr_all      = fixed_log_aln_probs;
    for (unsigned int i = 0; i < num_reads_; ++i){
      for (unsigned int j = 0; j < num_alleles_; ++j, ++log_aln_ptr_original)
	log_aln_ptr_all[original_indices[j]] = *log_aln_ptr_original;
      for (unsigned int j = 0; j < stutter_seqs.size(); ++j, ++log_aln_ptr_stutter)
	log_aln_ptr_all[stutter_indices[j]] = *log_aln_ptr_stutter;
      log_aln_ptr_all += total_alleles;
    }
    delete [] log_aln_probs_;
    delete [] new_log_aln_probs;
    log_aln_probs_ = fixed_log_aln_probs;

    // Fix the trace cache indexing
    std::map<std::pair<int,int>, AlignmentTrace*> new_trace_cache;
    for (auto cache_iter = trace_cache_.begin(); cache_iter != trace_cache_.end(); cache_iter++){
      int new_allele_index = original_indices[cache_iter->first.second];
      new_trace_cache[std::pair<int,int>(cache_iter->first.first, new_allele_index)] = cache_iter->second;
    }
    trace_cache_ = new_trace_cache;

    // Construct a haplotype that includes all the alleles
    delete haplotype_;
    delete hap_blocks_[1];
    num_alleles_   = total_alleles;
    hap_blocks_[1] = str_block;
    haplotype_     = new Haplotype(hap_blocks_);

    // Reextract the allele info
    alleles_.clear();
    get_alleles(chrom_seq, alleles_);

    // Reallocate and recompute genotype posteriors
    delete [] log_sample_posteriors_;
    log_sample_posteriors_ = new double[num_alleles_*num_alleles_*num_samples_];
    calc_log_sample_posteriors();

    stutter_seqs.clear();
    get_stutter_candidate_alleles(logger, stutter_seqs);
  }
  return true;
}

bool SeqStutterGenotyper::genotype(std::string& chrom_seq, std::ostream& logger){
  // Unsuccessful initialization. May be due to
  // 1) Failing to find the corresponding allele priors in the VCF (if one has been provided)
  // 2) Large deletion extending past STR
  if (pos_ == -1)
    return false;

  // If the smallest stutter block sequence is smaller than the maximum deletion size, the stutter aligner will fail
  // We could extend the stutter block to prevent this, but if this happens, the locus is probably not very high quality
  // As a result, for now, just abort the genotyping for this locus
  for (int i = 0; i < haplotype_->num_blocks(); ++i){
    HapBlock* hap_block = haplotype_->get_block(i);
    if (hap_block->get_repeat_info() == NULL)
      continue;
    if (hap_block->min_size() < std::abs(hap_block->get_repeat_info()->max_deletion()))
      return false;
  }

  init_alignment_model();
  if (pool_identical_seqs_){
    logger << "Pooling reads with identical sequences..." << std::endl;
    pooler_.pool(base_quality_);
  }

  // Align each read to each candidate haplotype and store them in the provided arrays
  logger << "Aligning reads to each candidate haplotype..." << std::endl;
  calc_hap_aln_probs(haplotype_, log_aln_probs_, seed_positions_);
  calc_log_sample_posteriors();

  // Look for additional alleles in stutter artifacts and align to them (if necessary)
  if (ref_vcf_ == NULL){
    if(!id_and_align_to_stutter_alleles(chrom_seq, logger))
      return false;
  }

  // Remove alleles with no MAP genotype calls and recompute the posteriors
  if (ref_vcf_ == NULL && log_allele_priors_ == NULL){
    std::vector<int> uncalled_indices;
    get_uncalled_alleles(uncalled_indices);
    if (uncalled_indices.size() != 0){
      logger << "Recomputing sample posteriors after removing " << uncalled_indices.size() << " uncalled alleles" << std::endl;
      remove_alleles(uncalled_indices);
    }

    uncalled_indices.clear();
    get_unspanned_alleles(uncalled_indices, logger);
    if (uncalled_indices.size() != 0){
      logger << "Recomputing sample posteriors after removing " << uncalled_indices.size() << " alleles with no spanning reads" << std::endl;
      remove_alleles(uncalled_indices);
    }
  }
  
  if (ref_vcf_ != NULL)
    pos_ += 1;
  return true;
}

void SeqStutterGenotyper::get_alleles(std::string& chrom_seq, std::vector<std::string>& alleles){
  assert(alleles.size() == 0);

  // Extract all the alleles
  std::vector<std::string> allele_seqs;
  do {
    allele_seqs.push_back(uppercase(haplotype_->get_seq()));
  } while(haplotype_->next());
  haplotype_->reset();

  // Trim from the left until the region boundary or a mismatched character
  int32_t left_trim = 0;
  int32_t start     = hap_blocks_.front()->start();
  while (start + left_trim < region_->start()){
    bool trim = true;
    for (unsigned int i = 0; i < allele_seqs.size(); ++i){
      if ((left_trim+1 >= allele_seqs[i].size()) || (allele_seqs[i][left_trim] != allele_seqs[0][left_trim])){
	trim = false;
	break;
      }
    }
    if (!trim) break;
    left_trim++;
  }
  start += left_trim;
  for (unsigned int i = 0; i < allele_seqs.size(); ++i)
    allele_seqs[i] = allele_seqs[i].substr(left_trim);

  // Trim from the right until the region boundary or a mismatched character
  int32_t right_trim = 0;
  int32_t end        = hap_blocks_.back()->end();
  while (end - right_trim > region_->stop()){
    bool trim    = true;
    int ref_size = allele_seqs[0].size();
    for (unsigned int i = 0; i < allele_seqs.size(); ++i){
      int alt_size = allele_seqs[i].size();
      if ((right_trim+1 >= allele_seqs[i].size()) || (allele_seqs[i][alt_size-right_trim-1] != allele_seqs[0][ref_size-right_trim-1])){
	trim = false;
	break;
      }
    }
    if (!trim) break;
    right_trim++;
  }
  end -= right_trim;
  for (unsigned int i = 0; i < allele_seqs.size(); ++i)
    allele_seqs[i] = allele_seqs[i].substr(0, allele_seqs[i].size()-right_trim);

  std::string left_flank  = (start >= region_->start() ? uppercase(chrom_seq.substr(region_->start(), start-region_->start())) : "");
  std::string right_flank = (end <= region_->stop()    ? uppercase(chrom_seq.substr(end, region_->stop()-end)) : "");
  pos_ = std::min((int32_t)region_->start(), start);

  // If necessary, add 1bp on the left so that all the alleles match the reference sequence
  if (left_flank.empty()){
    bool pad_left = false;
    for (unsigned int i = 1; i < allele_seqs.size(); ++i){
      if (allele_seqs[i][0] != allele_seqs[0][0]){
	pad_left = true;
	break;
      }
    }
    if (pad_left){
      pos_ -= 1;
      left_flank = uppercase(chrom_seq.substr(pos_, 1));
    }
  }

  for (unsigned int i = 0; i < allele_seqs.size(); ++i){
    std::stringstream ss;
    ss << left_flank << allele_seqs[i] << right_flank;
    alleles.push_back(ss.str());
  }

  pos_ += 1; // Fix off-by-1 VCF error
}
 
void SeqStutterGenotyper::debug_sample(int sample_index, std::ostream& logger){
  std::vector<AlignmentTrace*> traced_alns;
  retrace_alignments(logger, traced_alns);

  std::cerr << "DEBUGGING SAMPLE..." << std::endl;
  std::cerr << "READ LL's:" << std::endl;
  double* read_LL_ptr = log_aln_probs_;
  int read_index = 0;
  for (unsigned int i = 0; i < num_reads_; ++i){
    if(sample_label_[i] == sample_index){
      std::cerr << "\t" << "READ #" << read_index << ", SEED BASE=" << seed_positions_[i]
		<< ", TOTAL QUAL CORRECT= " << alns_[i].sum_log_prob_correct(base_quality_) << ", "
		<< bp_diffs_[i] << " " << max_index(read_LL_ptr, num_alleles_) << ", "
		<< log_p1_[read_index] << " " << log_p2_[read_index] <<  ", "
		<< alns_[i].get_sequence().substr(0, seed_positions_[i])
		<< " " << alns_[i].get_sequence().substr(seed_positions_[i]+1) << std::endl
		<< traced_alns[i]->hap_aln() << std::endl
		<< traced_alns[i]->traced_aln().get_alignment()  << std::endl
		<< traced_alns[i]->traced_aln().getCigarString() << std::endl;

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

bool SeqStutterGenotyper::use_read(AlignmentTrace* trace){
  return true;
}

void SeqStutterGenotyper::filter_alignments(std::ostream& logger, std::vector<int>& masked_reads){
  assert(masked_reads.size() == 0);
  masked_reads = std::vector<int>(num_samples_, 0);
  std::vector<AlignmentTrace*> traced_alns;
  retrace_alignments(logger, traced_alns);
  assert(traced_alns.size() == num_reads_);

  int32_t filt_count = 0, keep_count = 0;
  double* read_LL_ptr = log_aln_probs_;
  for (unsigned int read_index = 0; read_index < num_reads_; read_index++){
    if (seed_positions_[read_index] < 0){
      masked_reads[sample_label_[read_index]]++;
      read_LL_ptr += num_alleles_;
      continue;
    }
    assert(traced_alns[read_index] != NULL);

    // Zero out alignment probabilities for filtered reads
    if (!use_read(traced_alns[read_index])){
      seed_positions_[read_index] = -2;
      for (unsigned int i = 0; i < num_alleles_; ++i)
	read_LL_ptr[i] = 0;
      filt_count++;
      masked_reads[sample_label_[read_index]]++;
    }
    else
      keep_count++;
    read_LL_ptr += num_alleles_;
  }

  calc_log_sample_posteriors();
  if (filt_count != 0)
    logger << "Filtered " << filt_count << " out of " << filt_count+keep_count << " reads based on their ML alignment tracebacks" << "\n";
}

void SeqStutterGenotyper::retrace_alignments(std::ostream& logger, std::vector<AlignmentTrace*>& traced_alns){
  double trace_start = clock();
  assert(traced_alns.size() == 0);
  traced_alns.reserve(num_reads_);
  std::vector< std::pair<int, int> > haps;
  get_optimal_haplotypes(log_sample_posteriors_, haps);

  HapAligner hap_aligner(haplotype_);
  double* read_LL_ptr = log_aln_probs_;
  for (unsigned int read_index = 0; read_index < num_reads_; read_index++){
    if (seed_positions_[read_index] < 0){
      read_LL_ptr += num_alleles_;
      traced_alns.push_back(NULL);
      continue;
    }

    int hap_a    = haps[sample_label_[read_index]].first;
    int hap_b    = haps[sample_label_[read_index]].second;
    int best_hap = ((LOG_ONE_HALF+log_p1_[read_index]+read_LL_ptr[hap_a] >  LOG_ONE_HALF+log_p2_[read_index]+read_LL_ptr[hap_b]) ? hap_a : hap_b);

    AlignmentTrace* trace = NULL;
    std::pair<int,int> trace_key(pool_index_[read_index], best_hap);
    auto trace_iter = trace_cache_.find(trace_key);
    if (trace_iter == trace_cache_.end()){
      trace  = hap_aligner.trace_optimal_aln(alns_[read_index], seed_positions_[read_index], best_hap, &base_quality_);
      trace_cache_[trace_key] = trace;
    }
    else
      trace = trace_iter->second;

    traced_alns.push_back(trace);
    read_LL_ptr += num_alleles_;
  }
  total_aln_trace_time_ += (clock() - trace_start)/CLOCKS_PER_SEC;
}

void SeqStutterGenotyper::get_stutter_candidate_alleles(std::ostream& logger, std::vector<std::string>& candidate_seqs){
  assert(candidate_seqs.size() == 0);

  // Get the index for the STR haplotype block and check that only one block exists
  // TO DO: We need to modify this function to accomodate multiple STR blocks, but for now, let's just die
  int str_block_index = -1;
  for (int i = 0; i < haplotype_->num_blocks(); ++i)
    if (haplotype_->get_block(i)->get_repeat_info() != NULL){
      assert(str_block_index == -1);
      str_block_index = i;
    }
  assert(str_block_index != -1);


  std::vector<AlignmentTrace*> traced_alns;
  retrace_alignments(logger, traced_alns);

  std::vector<int> sample_counts(num_samples_, 0);
  std::vector< std::map<std::string, int> > sample_stutter_counts(num_samples_);

  HapBlock* str_block = haplotype_->get_block(str_block_index);
  for (unsigned int read_index = 0; read_index < num_reads_; read_index++){
    if (traced_alns[read_index] == NULL)
      continue;
    AlignmentTrace* trace = traced_alns[read_index];
    if (trace->traced_aln().get_start() < str_block->start()){
      if (trace->traced_aln().get_stop() > str_block->end()){
	if (trace->stutter_size(str_block_index) != 0)
	  sample_stutter_counts[sample_label_[read_index]][trace->str_seq(str_block_index)]++;
	sample_counts[sample_label_[read_index]]++;
      }
    }
  }

  std::set<std::string> candidate_set;
  for (unsigned int i = 0; i < num_samples_; i++)
    for (auto seq_iter = sample_stutter_counts[i].begin(); seq_iter != sample_stutter_counts[i].end(); seq_iter++)
      if (!str_block->contains(seq_iter->first))
	if (seq_iter->second >= 2 && 1.0*seq_iter->second/sample_counts[i] >= 0.15)
	  candidate_set.insert(seq_iter->first);
  candidate_seqs = std::vector<std::string>(candidate_set.begin(), candidate_set.end());

  logger << "Identified " << candidate_seqs.size() << " additional candidate alleles from stutter artifacts" << "\n";
  for (unsigned int i = 0; i < candidate_seqs.size(); i++)
    logger << "\t" << candidate_seqs[i] << "\n";
}

void SeqStutterGenotyper::analyze_flank_indels(std::ostream& logger){
  std::vector<AlignmentTrace*> traced_alns;
  retrace_alignments(logger, traced_alns);
  std::vector<int> sample_counts(num_samples_, 0);
  std::vector< std::map<std::pair<int,int>, int> > sample_flank_indel_counts(num_samples_);

  for (unsigned int read_index = 0; read_index < num_reads_; read_index++){
    if (traced_alns[read_index] == NULL)
      continue;
    AlignmentTrace* trace = traced_alns[read_index];

    if (!trace->has_stutter()){
      bool use = false;
      use     |= (trace->flank_ins_size() == 0 && trace->flank_del_size() != 0);
      use     |= (trace->flank_ins_size() != 0 && trace->flank_del_size() == 0);
      use     &= (trace->flank_indel_data().size() == 1);
      if (use)
	sample_flank_indel_counts[sample_label_[read_index]][trace->flank_indel_data()[0]]++;
    }
    sample_counts[sample_label_[read_index]]++;
  }

  std::map< std::pair<int,int>, int>  candidate_set;
  for (unsigned int i = 0; i < num_samples_; i++)
    for (auto indel_iter = sample_flank_indel_counts[i].begin(); indel_iter != sample_flank_indel_counts[i].end(); indel_iter++)
      if (indel_iter->second >= 2 && 1.0*indel_iter->second/sample_counts[i] >= 0.15){
	  candidate_set[indel_iter->first]++;
	  //std::cerr << sample_names_[i] << " " << indel_iter->first.first << " " << indel_iter->first.second << std::endl;
      }

  if (candidate_set.size() != 0){
    for (auto candidate_iter = candidate_set.begin(); candidate_iter != candidate_set.end(); candidate_iter++){
      std::cerr << region_->chrom() << "\t" << pos_ << "\t" << (region_->name().empty() ? "." : region_->name()) << " ";
      std::cerr << candidate_iter->first.first << " " << candidate_iter->first.second << " " << candidate_iter->second << std::endl;
    }
  }
}


void SeqStutterGenotyper::write_vcf_record(std::vector<std::string>& sample_names, bool print_info, std::string& chrom_seq,
					   bool output_bootstrap_qualities, bool output_gls, bool output_pls, bool output_phased_gls,
					   bool output_allreads, bool output_pallreads, bool output_mallreads, bool output_viz, float max_flank_indel_frac,
					   bool visualize_left_alns,
					   std::ostream& html_output, std::ostream& out, std::ostream& logger){
  assert(haplotype_->num_blocks() == 3);

  //analyze_flank_indels(logger);
  //debug_sample(sample_indices_["SSC11604"], logger);

  if(log_allele_priors_ != NULL)
    assert(!output_gls && !output_pls && !output_phased_gls); // These fields only make sense in the context of MLE estimation, not MAP estimation

  // Compute the base pair differences from the reference
  std::vector<int> allele_bp_diffs;
  for (unsigned int i = 0; i < alleles_.size(); i++)
    allele_bp_diffs.push_back((int)alleles_[i].size() - (int)alleles_[0].size());
  
  // Filter reads with questionable alignments
  std::vector<int> masked_reads;
  filter_alignments(logger, masked_reads);

  // Extract the optimal genotyps and their associated likelihoods
  std::vector< std::pair<int,int> > haplotypes, gts;
  std::vector<double> log_phased_posteriors, log_unphased_posteriors, gl_diffs;
  std::vector< std::vector<double> > gls, phased_gls;
  std::vector< std::vector<int> > pls;
  int hap_block_index = 1; // TO DO: Determine this dynamically
  std::vector<int> hap_to_allele;
  haps_to_alleles(hap_block_index, hap_to_allele);
  int num_variants = haplotype_->get_block(hap_block_index)->num_options();
  extract_genotypes_and_likelihoods(num_variants, hap_to_allele, log_sample_posteriors_, haplotypes, gts, log_phased_posteriors, log_unphased_posteriors,
				    true, gls, gl_diffs, output_pls, pls, output_phased_gls, phased_gls);

  // Extract information about each read and group by sample
  assert(bp_diffs_.size() == num_reads_);
  std::vector<int> num_aligned_reads(num_samples_, 0), num_reads_with_snps(num_samples_, 0);
  std::vector<int> num_reads_with_stutter(num_samples_, 0), num_reads_with_flank_indels(num_samples_, 0);
  std::vector<int> num_reads_strand_one(num_samples_, 0), num_reads_strand_two(num_samples_, 0);
  std::vector< std::vector<int> > bps_per_sample(num_samples_), ml_bps_per_sample(num_samples_);
  std::vector< std::vector<double> > log_read_phases(num_samples_), posterior_bps_per_sample(num_samples_);
  std::vector<AlnList> max_LL_alns_strand_one(num_samples_), left_alns_strand_one(num_samples_), orig_alns_strand_one(num_samples_);
  std::vector<AlnList> max_LL_alns_strand_two(num_samples_), left_alns_strand_two(num_samples_), orig_alns_strand_two(num_samples_);

  HapAligner hap_aligner(haplotype_);
  double* read_LL_ptr = log_aln_probs_;
  for (unsigned int read_index = 0; read_index < num_reads_; read_index++){
    if (seed_positions_[read_index] < 0){
      read_LL_ptr += num_alleles_;
      continue;
    }

    // Extract read's phase posterior conditioned on the determined sample genotype
    int hap_a = gts[sample_label_[read_index]].first;
    int hap_b = gts[sample_label_[read_index]].second;
    double total_read_LL = log_sum_exp(LOG_ONE_HALF+log_p1_[read_index]+read_LL_ptr[hap_a], LOG_ONE_HALF+log_p2_[read_index]+read_LL_ptr[hap_b]);
    double log_phase_one = LOG_ONE_HALF + log_p1_[read_index] + read_LL_ptr[hap_a] - total_read_LL; 
    log_read_phases[sample_label_[read_index]].push_back(log_phase_one);

    // Determine which of the two genotypes each read is associated with
    int read_strand = 0;
    if (!haploid_ && ((hap_a != hap_b) || (abs(log_p1_[read_index]- log_p2_[read_index]) > TOLERANCE))){
      double v1 = log_p1_[read_index]+read_LL_ptr[hap_a], v2 = log_p2_[read_index]+read_LL_ptr[hap_b];
      if (abs(v1-v2) > TOLERANCE)
	read_strand = (v1 > v2 ? 0 : 1);
    }

    // Retrace alignment and ensure that it's of sufficient quality
    double trace_start = clock();
    int best_hap = (read_strand == 0 ? hap_a : hap_b);
    int best_gt  = (read_strand == 0 ? gts[sample_label_[read_index]].first: gts[sample_label_[read_index]].second);
    AlignmentTrace* trace = NULL;
    std::pair<int,int> trace_key(pool_index_[read_index], best_hap);
    auto trace_iter = trace_cache_.find(trace_key);
    if (trace_iter == trace_cache_.end()){
      trace  = hap_aligner.trace_optimal_aln(alns_[read_index], seed_positions_[read_index], best_hap, &base_quality_);
      trace_cache_[trace_key] = trace;
    }
    else
      trace = trace_iter->second;

    if (trace->has_stutter())
      num_reads_with_stutter[sample_label_[read_index]]++;
    if (trace->flank_ins_size() != 0 || trace->flank_del_size() != 0)
      num_reads_with_flank_indels[sample_label_[read_index]]++;

    if (visualize_left_alns)
      (read_strand == 0 ? left_alns_strand_one : left_alns_strand_two)[sample_label_[read_index]].push_back(alns_[read_index]);
    (read_strand == 0 ? max_LL_alns_strand_one : max_LL_alns_strand_two)[sample_label_[read_index]].push_back(trace->traced_aln());
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

    // Extract the ML bp difference observed in read based on the ML genotype,
    // but only for reads that span the original repeat region by 5 bp
    if (trace->traced_aln().get_start() < (region_->start() > 4 ? region_->start()-4 : 0))
      if (trace->traced_aln().get_stop() > region_->stop() + 4)
	ml_bps_per_sample[sample_label_[read_index]].push_back(allele_bp_diffs[best_gt]+trace->total_stutter_size());

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
  int sample_index = 0, skip_count = 0, filt_count = 0, allele_number = 0;
  for (auto gt_iter = gts.begin(); gt_iter != gts.end(); ++gt_iter, ++sample_index){
    if (samples_of_interest.find(sample_names_[sample_index]) == samples_of_interest.end())
      continue;
    if (require_one_read_ && num_aligned_reads[sample_index] == 0)
      continue;
    if (num_aligned_reads[sample_index] > 0 &&
        (num_reads_with_flank_indels[sample_index] > max_flank_indel_frac*num_aligned_reads[sample_index])){
      filt_count++;
      continue;
    }
    if (call_sample_[sample_index]) {
      if (haploid_){
	assert(gt_iter->first == gt_iter->second);
	allele_counts[gt_iter->first]++;
	allele_number++;
      }
      else {
	allele_counts[gt_iter->first]++;
	allele_counts[gt_iter->second]++;
	allele_number += 2;
      }
    }
    else
      skip_count++;
  }

  if (print_info){
    logger << "Allele counts" << std::endl;
    for (unsigned int i = 0; i < alleles_.size(); i++)
      logger << alleles_[i] << " " << allele_counts[i] <<  std::endl;
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

  // Obtain relevant stutter model. For now, get it from first repeat block
  // TO DO: Generalize this
  assert(haplotype_->get_block(1)->get_repeat_info() != NULL);
  StutterModel* stutter_model = haplotype_->get_block(1)->get_repeat_info()->get_stutter_model();

  // Add INFO field items
  out << "\tINFRAME_PGEOM=" << stutter_model->get_parameter(true,  'P') << ";"
      << "INFRAME_UP="      << stutter_model->get_parameter(true,  'U') << ";"
      << "INFRAME_DOWN="    << stutter_model->get_parameter(true,  'D') << ";"
      << "OUTFRAME_PGEOM="  << stutter_model->get_parameter(false, 'P') << ";"
      << "OUTFRAME_UP="     << stutter_model->get_parameter(false, 'U') << ";"
      << "OUTFRAME_DOWN="   << stutter_model->get_parameter(false, 'D') << ";"
      << "START="           << region_->start()+1 << ";"
      << "END="             << region_->stop()    << ";"
      << "PERIOD="          << region_->period()  << ";"
      << "NSKIP="           << skip_count         << ";"
      << "NFILT="           << filt_count         << ";";
  if (num_alleles_ > 1){
    out << "BPDIFFS=" << allele_bp_diffs[1];
    for (unsigned int i = 2; i < num_alleles_; i++)
      out << "," << allele_bp_diffs[i];
    out << ";";
  }

  // Compute INFO field values for DP, DFILT, DSTUTTER and DFLANKINDEL and add them to the VCF
  int32_t tot_dp = 0, tot_dsnp = 0, tot_dfilt = 0, tot_dstutter = 0, tot_dflankindel = 0;
  for (unsigned int i = 0; i < sample_names.size(); i++){
    auto sample_iter = sample_indices_.find(sample_names[i]);
    if (sample_iter == sample_indices_.end())
      continue;
    if (!call_sample_[sample_iter->second])
      continue;
    if (num_aligned_reads[sample_iter->second] > 0 &&
	(num_reads_with_flank_indels[sample_iter->second] > num_aligned_reads[sample_iter->second]*max_flank_indel_frac))
      continue;

    int sample_index = sample_iter->second;
    tot_dp          += num_aligned_reads[sample_index];
    tot_dsnp        += num_reads_with_snps[sample_index];
    tot_dfilt       += masked_reads[sample_index];
    tot_dstutter    += num_reads_with_stutter[sample_index];
    tot_dflankindel += num_reads_with_flank_indels[sample_index];
  }
  out << "DP="          << tot_dp          << ";"
      << "DSNP="        << tot_dsnp        << ";"
      << "DFILT="       << tot_dfilt       << ";"
      << "DSTUTTER="    << tot_dstutter    << ";"
      << "DFLANKINDEL=" << tot_dflankindel << ";";

  // Add allele counts
  out << "AN=" << allele_number << ";" << "REFAC=" << allele_counts[0];
  if (allele_counts.size() > 1){
    out << ";AC=";
    for (unsigned int i = 1; i < allele_counts.size()-1; i++)
      out << allele_counts[i] << ",";
    out << allele_counts.back();
  }

  // Add FORMAT field
  out << (!haploid_ ? "\tGT:GB:Q:PQ:DP:DSNP:DFILT:DSTUTTER:DFLANKINDEL:PDP:PSNP:BPDOSE:GLDIFF" : "\tGT:GB:Q:DP:DFILT:DSTUTTER:DFLANKINDEL:BPDOSE:GLDIFF");
  if (output_bootstrap_qualities) out << ":BQ";
  if (output_allreads)            out << ":ALLREADS";
  if (output_pallreads)           out << ":PALLREADS";
  if (output_mallreads)           out << ":MALLREADS";
  if (output_gls)                 out << ":GL";
  if (output_pls)                 out << ":PL";
  if (output_phased_gls)          out << ":PHASEDGL";

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

    // Don't report genotype for a sample if it exceeds the flank indel fraction
    if (num_aligned_reads[sample_iter->second] > 0 &&
	(num_reads_with_flank_indels[sample_iter->second] > num_aligned_reads[sample_iter->second]*max_flank_indel_frac)){
      out << ".";
      continue;
    }
    
    int sample_index    = sample_iter->second;
    double phase1_reads = (num_aligned_reads[sample_index] == 0 ? 0 : exp(log_sum_exp(log_read_phases[sample_index])));
    double phase2_reads = num_aligned_reads[sample_index] - phase1_reads;

    std::stringstream samp_info;
    samp_info << allele_bp_diffs[gts[sample_index].first] << "|" << allele_bp_diffs[gts[sample_index].second];
    sample_results[sample_names[i]] = samp_info.str();

    // TO DO: Compute p-value for allele read depth bias
    // i)  Spanning reads
    // ii) All reads if --use-all-reads is specified
    // We will use the  bdtr(k, N, p) function from the cephes directory, which computes the CDF for a binomial distribution
    // e.g.: double val = bdtr (24, 50, 0.5);

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
	  << ":" << num_reads_strand_one[sample_index] << "|" << num_reads_strand_two[sample_index]; // Reads with SNPs supporting each haploid genotype

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
	  << ":" << num_reads_with_flank_indels[sample_index];                                      // Total reads with an indel in flank in ML alignment

      // Difference in GL between the current and next best genotype
      if (num_alleles_ == 1)
	out << ":" << ".";
      else
	out << ":" << gl_diffs[sample_index];
    }

    if (output_bootstrap_qualities)
      out << ":" << bootstrap_qualities[sample_index];

    // Add bp diffs from regular left-alignment
    if (output_allreads)
	out << ":" << condense_read_counts(bps_per_sample[sample_index]);

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
    if (output_mallreads)
	out << ":" << condense_read_counts(ml_bps_per_sample[sample_index]);

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
    if (output_phased_gls){
      out << ":" << phased_gls[sample_index][0];
      for (unsigned int j = 1; j < phased_gls[sample_index].size(); j++)
	out << "," << phased_gls[sample_index][j];
    }
  }
  out << "\n";

  // Render HTML of Smith-Waterman alignments (or haplotype alignments)
  if (output_viz){
    // Combine alignments from both strands after ordering them by position independently
    std::vector<AlnList> max_LL_alns(num_samples_);
    for (unsigned int i = 0; i < num_samples_; i++){
      for (unsigned int j = 0; j < 3; j++){
	AlnList& aln_ref_one = (j == 0 ? orig_alns_strand_one[i] : (j == 1 ? left_alns_strand_one[i] : max_LL_alns_strand_one[i]));
	AlnList& aln_ref_two = (j == 0 ? orig_alns_strand_two[i] : (j == 1 ? left_alns_strand_two[i] : max_LL_alns_strand_two[i]));
	std::sort(aln_ref_one.begin(), aln_ref_one.end());
	std::sort(aln_ref_two.begin(), aln_ref_two.end());
	max_LL_alns[i].insert(max_LL_alns[i].end(), aln_ref_one.begin(), aln_ref_one.end());
	max_LL_alns[i].insert(max_LL_alns[i].end(), aln_ref_two.begin(), aln_ref_two.end());
	aln_ref_one.clear(); aln_ref_two.clear();
      }
    }

    std::stringstream locus_info;
    locus_info << region_->chrom() << "\t" << region_->start()+1 << "\t" << region_->stop();
    visualizeAlignments(max_LL_alns, sample_names_, sample_results, hap_blocks_, chrom_seq, locus_info.str(), true, html_output);
  }
}

bool SeqStutterGenotyper::recompute_stutter_models(std::string& chrom_seq, std::ostream& logger,
						  int max_em_iter, double abs_ll_converge, double frac_ll_converge){
  logger << "Retraining EM stutter genotyper using maximum likelihood alignments" << std::endl;
  std::vector<AlignmentTrace*> traced_alns;
  retrace_alignments(logger, traced_alns);
  assert(traced_alns.size() == num_reads_);

  for (int block_index = 0; block_index < haplotype_->num_blocks(); ++block_index){
    HapBlock* block = haplotype_->get_block(block_index);
    if (block->get_repeat_info() == NULL)
      continue;

    std::vector< std::vector<int> > str_num_bps(num_samples_);
    std::vector< std::vector<double> > str_log_p1s(num_samples_), str_log_p2s(num_samples_);
    for (unsigned int read_index = 0; read_index < num_reads_; ++read_index){
      AlignmentTrace* trace = traced_alns[read_index];
      if (trace != NULL){
	if (trace->traced_aln().get_start() < block->start()){
	  if (trace->traced_aln().get_stop() > block->end()){
	    str_num_bps[sample_label_[read_index]].push_back(((int)trace->str_seq(block_index).size())+trace->stutter_size(block_index));
	    str_log_p1s[sample_label_[read_index]].push_back(log_p1_[read_index]);
	    str_log_p2s[sample_label_[read_index]].push_back(log_p2_[read_index]);
	  }
	}
      }
    }

    int period = block->get_repeat_info()->get_period();
    EMStutterGenotyper length_genotyper(*region_, haploid_, str_num_bps, str_log_p1s, str_log_p2s, sample_names_, 0);
    bool trained = length_genotyper.train(max_em_iter, abs_ll_converge, frac_ll_converge, false, logger);
    if (!trained){
      logger << "Retraining stutter model training failed for locus " << region_->chrom() << ":" << region_->start() << "-" << region_->stop() << std::endl;
      return false;
    }

    logger << "Learned stutter model: " << (*length_genotyper.get_stutter_model()) << std::endl;
    block->get_repeat_info()->set_stutter_model(length_genotyper.get_stutter_model());
  }
  trace_cache_.clear();
  return genotype(chrom_seq, logger);
}

void SeqStutterGenotyper::compute_bootstrap_qualities(int num_iter, std::vector<double>& bootstrap_qualities){
  assert(bootstrap_qualities.size() == 0);
  double bootstrap_start = clock();

  // Extract the original ML genotypes
  std::vector< std::pair<int, int> > ML_gts;
  get_optimal_haplotypes(log_sample_posteriors_, ML_gts);

  // Partition the aligned reads by sample
  std::vector< std::vector<int> > reads_by_sample(num_samples_);
  for (unsigned int i = 0; i < num_reads_; i++)
    if (seed_positions_[i] >= 0)
      reads_by_sample.at(sample_label_[i]).push_back(i);

  std::vector<int> ML_gt_counts(num_samples_, 0);
  std::uniform_int_distribution<int> unif_dist;
  std::default_random_engine gen;
  double log_homoz_prior = log_homozygous_prior(), log_hetz_prior = log_heterozygous_prior();
  for (unsigned int i = 0; i < num_samples_; i++){
    int num_sample_reads = reads_by_sample[i].size();

    // Precompute all read LLs for each of the sample's diploid genotypes
    double* read_gt_LLs = new double[num_alleles_*num_alleles_*num_sample_reads];
    double* ptr         = read_gt_LLs;
    for (int index_1 = 0; index_1 < num_alleles_; ++index_1){
      for (int index_2 = 0; index_2 < num_alleles_; ++index_2){
	for (auto read_index = reads_by_sample[i].begin(); read_index != reads_by_sample[i].end(); ++read_index, ++ptr){
	  double* read_LL_ptr = log_aln_probs_ + *read_index*num_alleles_;
	  *ptr = log_sum_exp(LOG_ONE_HALF + log_p1_[*read_index] + read_LL_ptr[index_1],
			     LOG_ONE_HALF + log_p2_[*read_index] + read_LL_ptr[index_2]);
	}
      }
    }

    for (int j = 0; j < num_iter; j++){
      // Bootstrap reads for the sample
      std::vector<int> bootstrap_weights(num_sample_reads, 0);
      for (unsigned int k = 0; k < num_sample_reads; k++)
	bootstrap_weights[unif_dist(gen) % num_sample_reads]++;

      // Recompute the genotype posteriors using bootstrapped read weights
      double* read_LL_ptr = read_gt_LLs;
      std::pair<int, int> bootstrap_gt(0, 0);
      double bootstrap_max_LL = -10e6;
      for (int index_1 = 0; index_1 < num_alleles_; ++index_1){
	for (int index_2 = 0; index_2 < num_alleles_; ++index_2){
	  if (haploid_ && (index_1 != index_2))
	    continue;
	  double gt_LL = (index_1 == index_2 ? log_homoz_prior: log_hetz_prior);
	  for (int read_index = 0; read_index < num_sample_reads; ++read_index, ++read_LL_ptr)
	    gt_LL += bootstrap_weights[read_index]*(*read_LL_ptr);
	  if (gt_LL > bootstrap_max_LL){
	    bootstrap_gt     = std::pair<int,int>(index_1, index_2);
	    bootstrap_max_LL = gt_LL;
	  }
	}
      }

      // Increment count if bootstrapped ML genotype (unordered) matches the ML genotype
      if (bootstrap_gt.first == ML_gts[i].first && bootstrap_gt.second == ML_gts[i].second)
	ML_gt_counts[i]++;
      else if (bootstrap_gt.first == ML_gts[i].second && bootstrap_gt.second == ML_gts[i].first)
	ML_gt_counts[i]++;
    }
    delete [] read_gt_LLs;
  }

  // Compute the boostrapped qualities as the fraction of iterations in which the genotype matched
  for (unsigned int i = 0; i < num_samples_; i++)
    bootstrap_qualities.push_back(1.0*ML_gt_counts[i]/num_iter);

  double bootstrap_time  = (clock() - bootstrap_start)/CLOCKS_PER_SEC;
  total_bootstrap_time_ += bootstrap_time;
}
