#include <algorithm>
#include <cfloat>
#include <cstring>
#include <math.h>
#include <random>
#include <string>
#include <sstream>
#include <time.h>

#include "seq_stutter_genotyper.h"
#include "bam_processor.h"
#include "debruijn_graph.h"
#include "em_stutter_genotyper.h"
#include "error.h"
#include "extract_indels.h"
#include "mathops.h"
#include "stringops.h"
#include "vcf_input.h"
#include "zalgorithm.h"

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

bool SeqStutterGenotyper::assemble_flanks(std::ostream& logger){
  std::vector<AlignmentTrace*> traced_alns;
  retrace_alignments(traced_alns);

  double locus_assembly_time = clock();
  logger << "Reassembling flanking sequences" << std::endl;
  std::vector< std::vector<std::string> > alleles_to_add (haplotype_->num_blocks());
  std::vector<bool> realign_sample(num_samples_, false);

  for (int flank = 0; flank < 2; flank++){
    int block_index     = (flank == 0 ? 0 : haplotype_->num_blocks()-1);
    std::string ref_seq = hap_blocks_[block_index]->get_seq(0);
    int max_k           = std::min(MAX_KMER, ref_seq.size() == 0 ? -1 : (int)ref_seq.size()-1);
    int kmer_length;
    if (!DebruijnGraph::calc_kmer_length(ref_seq, MIN_KMER, max_k, kmer_length))
      return false;

    std::map<std::string, int> haplotype_counts;
    std::vector< std::pair<std::string,int> > assembly_data;
    int min_read_index = 0, read_index;
    for (int sample_index = 0; sample_index < num_samples_; sample_index++){
      assembly_data.clear();
      bool acyclic = false;
      for (int k = kmer_length; k <= max_k; k++){
	DebruijnGraph assembler(k, ref_seq);
	for (read_index = min_read_index; read_index < num_reads_; read_index++){
	  if (sample_label_[read_index] != sample_index)
	    break;
	  if (traced_alns[read_index] == NULL)
	    continue;

	  std::string seq = traced_alns[read_index]->flank_seq(block_index);
	  if (!seq.empty())
	    assembler.add_string(seq);
	}

	assembler.prune_edges(0.02, 2);
	if (!assembler.has_cycles() && assembler.is_source_ok() && assembler.is_sink_ok()){
	  acyclic = true;
	  assembler.enumerate_paths(MIN_PATH_WEIGHT, 10, assembly_data);
	  break;
	}
      }
      min_read_index = read_index;

      if (acyclic){
	if (call_sample_[sample_index] && assembly_data.size() > 1){
	  int total_depth = 0;
	  for (unsigned int i = 0; i < assembly_data.size(); i++)
	    total_depth += assembly_data[i].second;

	  for (unsigned int i = 0; i < assembly_data.size(); i++){
	    if (assembly_data[i].first.compare(ref_seq) != 0){
	      if (assembly_data[i].second*1.0/total_depth > 0.25){
		if (ref_seq.size() != assembly_data[i].first.size()){
		  // (i) Mask a sample if it supports a flank that contains an indel and (ii) ignore the candidate flanking sequence
		  // NOTE: We may want to incorporate alternate flanks with indels at a later stage,
		  // but right now there's an issue with indels in the flanks clobbering indels in the STR region
		  // We need to be more intelligent about how we resolve this before we can include these flanks
		  // For now, we'll avoid genotyping any samples with support for these indels
		  call_sample_[sample_index]   = false;
		  realign_sample[sample_index] = false;
		}
		else {
		  //logger << sample_names_[sample_index] << " " << total_depth << " " << assembly_data[i].second << std::endl;
		  realign_sample[sample_index] = true;
		  haplotype_counts[assembly_data[i].first]++;
		}
	      }
	    }
	  }
	}
      }
      else
	call_sample_[sample_index] = false;
    }

    if (!haplotype_counts.empty()){
      logger << "Identified " << haplotype_counts.size() << " new flank haplotype(s)" << "\n";
      for (auto hap_iter = haplotype_counts.begin(); hap_iter != haplotype_counts.end(); hap_iter++){
	logger << "\t" << hap_iter->first << "\t" << hap_iter->second << "\n";
	alleles_to_add[block_index].push_back(hap_iter->first);
      }
      logger << "\t" << ref_seq << "\t" << "REF_SEQ" << "\n" << std::endl;
    }
  }
  locus_assembly_time   = (clock() - locus_assembly_time)/CLOCKS_PER_SEC;
  total_assembly_time_ += locus_assembly_time;

  // Determine which read pools we need to realign and which read's probabilities we should update
  // We need to realign a pool if any of its associated reads has a sample with a new candidate flank
  // We only want to update a read's probabilities if all of the other read's for the sample were also realigned
  std::vector<bool> realign_pools(pooler_.num_pools(), false);
  std::vector<bool> copy_reads(num_reads_, false);
  for (unsigned int read_index = 0; read_index < num_reads_; read_index++){
    bool sample_flag = realign_sample[sample_label_[read_index]];
    realign_pools[pool_index_[read_index]] = (realign_pools[pool_index_[read_index]] || sample_flag);
    copy_reads[read_index] = sample_flag;
  }

  // Realign reads for samples with new candidate flanking sequences
  int realign_count = 0;
  for (unsigned int i = 0; i < realign_pools.size(); i++)
    if (realign_pools[i])
      realign_count++;
  if (realign_count > 0){
    logger << "Realigning " << realign_count << " out of " << realign_pools.size() << " read pools to polish flanking sequences" << std::endl;
    std::vector< std::vector<int> > alleles_to_remove(haplotype_->num_blocks());
    add_and_remove_alleles(alleles_to_remove, alleles_to_add, realign_pools, copy_reads);


    // Remove alleles with no MAP genotype calls and recompute the posteriors
    if (ref_vcf_ == NULL){
      std::vector< std::vector<int> > unused_indices;
      int num_aff_blocks = 0, num_aff_alleles = 0;
      get_unused_alleles(false, true, unused_indices, num_aff_blocks, num_aff_alleles);
      if (num_aff_alleles != 0){
	logger << "Recomputing sample posteriors after removing " << num_aff_alleles
	       << " uncalled alleles across " << num_aff_blocks << " blocks" << std::endl;
	remove_alleles(unused_indices);
      }
    }
  }

  return true;
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

void SeqStutterGenotyper::get_unused_alleles(bool check_spanned, bool check_called,
					     std::vector< std::vector<int> >& allele_indices, int& num_aff_blocks, int& num_aff_alleles){
  assert(allele_indices.size() == 0);
  num_aff_blocks  = 0;
  num_aff_alleles = 0;

  // Extract each sample's optimal haplotypes
  std::vector< std::pair<int,int> > haps;
  get_optimal_haplotypes(haps);

  // Retrace the alignments
  std::vector<AlignmentTrace*> traced_alns;
  retrace_alignments(traced_alns);

  // Determine which samples have >= 1 aligned read
  std::vector<bool> aligned_read(num_samples_, false);
  for (unsigned int read_index = 0; read_index < num_reads_; read_index++)
    if (seed_positions_[read_index] >= 0)
      aligned_read[sample_label_[read_index]] = true;

  // Iterate over each block in the haplotype
  for (int block_index = 0; block_index < haplotype_->num_blocks(); block_index++){
    allele_indices.push_back(std::vector<int>());

    // Skip blocks in which only the reference allele is present
    HapBlock* block = haplotype_->get_block(block_index);
    if (block->num_options() == 1)
      continue;

    // Determine the mapping from haplotype->allele for the current block
    std::vector<int> hap_to_allele;
    haps_to_alleles(block_index, hap_to_allele);

    std::vector<bool> spanned(block->num_options(), false);
    std::vector<bool>  called(block->num_options(), false);

    // Mark all alleles spanned by at least one ML alignment
    if (check_spanned){
      double* read_LL_ptr = log_aln_probs_;
      for (unsigned int read_index = 0; read_index < num_reads_; read_index++){
	if (seed_positions_[read_index] < 0){
	  read_LL_ptr += num_alleles_;
	  continue;
	}
	AlignmentTrace* trace = traced_alns[read_index];
	if (trace->traced_aln().get_start() < block->start()){
	  if (trace->traced_aln().get_stop() > block->end()){
	    if (trace->stutter_size(block_index) == 0){
	      int hap_a    = haps[sample_label_[read_index]].first;
	      int hap_b    = haps[sample_label_[read_index]].second;
	      int best_hap = hap_a;
	      if (!haploid_ && (hap_a != hap_b)){
		double v1 = log_p1_[read_index]+read_LL_ptr[hap_a], v2 = log_p2_[read_index]+read_LL_ptr[hap_b];
		if (abs(v1-v2) > TOLERANCE)
		  best_hap = (v1 > v2 ? hap_a : hap_b);
	      }
	      spanned[hap_to_allele[best_hap]] = true;
	    }
	  }
	}
	read_LL_ptr += num_alleles_;
      }
    }

    // Mark all alleles where at least one sample contains the allele in its optimal haplotypes
    if (check_called){
      for (unsigned int sample_index = 0; sample_index < haps.size(); sample_index++){
	if (aligned_read[sample_index] && call_sample_[sample_index]){
	  called[hap_to_allele[haps[sample_index].first]]  = true;
	  called[hap_to_allele[haps[sample_index].second]] = true;
	}
      }
    }

    // Add non-reference alleles that fail the requirements
    bool block_affected = false;
    for (int allele_index = 1; allele_index < block->num_options(); allele_index++){
      if ((check_spanned && !spanned[allele_index]) || (check_called && !called[allele_index])){
	allele_indices.back().push_back(allele_index);
	block_affected = true;
	num_aff_alleles++;
      }
    }
    if (block_affected)
      num_aff_blocks++;
  }
}

void SeqStutterGenotyper::add_and_remove_alleles(std::vector< std::vector<int> >& alleles_to_remove,
                                                 std::vector< std::vector<std::string> >& alleles_to_add){
  std::vector<bool> realign_pool = std::vector<bool>(pooler_.num_pools(), true);
  std::vector<bool> copy_read    = std::vector<bool>(num_reads_, true);
  add_and_remove_alleles(alleles_to_remove, alleles_to_add, realign_pool, copy_read);
}

void SeqStutterGenotyper::add_and_remove_alleles(std::vector< std::vector<int> >& alleles_to_remove,
						 std::vector< std::vector<std::string> >& alleles_to_add,
						 std::vector<bool>& realign_pool, std::vector<bool>& copy_read){
  assert(alleles_to_remove.size() == hap_blocks_.size() && alleles_to_add.size() == hap_blocks_.size());

  // Store the initial set of haplotype sequences
  haplotype_->reset();
  std::map<std::string, int> hap_indices;
  std::vector<std::string> hap_seqs;
  do {
    hap_seqs.push_back(haplotype_->get_seq());
    hap_indices[hap_seqs.back()] = hap_seqs.size()-1;
  } while(haplotype_->next());
  haplotype_->reset();

  // Construct new haplotype blocks by removing the unwanted alleles and adding the new alleles
  std::vector<HapBlock*> updated_blocks;
  for (int i = 0; i < hap_blocks_.size(); i++)
    updated_blocks.push_back(hap_blocks_[i]->remove_alleles(alleles_to_remove[i]));

  bool added_seq = false;
  for (int i = 0; i < updated_blocks.size(); i++){
    for (int j = 0; j < alleles_to_add[i].size(); j++){
      updated_blocks[i]->add_alternate(alleles_to_add[i][j]);
      added_seq = true;
    }
  }

  // Construct the new haplotype and record its set of haplotype sequences
  // Determine the mapping from old sequences to new sequences, if they're still present
  Haplotype* updated_haplotype = new Haplotype(updated_blocks);
  std::vector<std::string> updated_hap_seqs;
  std::vector<int> allele_mapping(num_alleles_, -1);
  std::vector<bool> realign_to_haplotype;
  do {
    updated_hap_seqs.push_back(updated_haplotype->get_seq());
    auto match = hap_indices.find(updated_hap_seqs.back());
    if (match == hap_indices.end())
      realign_to_haplotype.push_back(true);
    else {
      realign_to_haplotype.push_back(false);
      allele_mapping[match->second] = updated_hap_seqs.size()-1;
    }
  } while(updated_haplotype->next());
  updated_haplotype->reset();
  assert(updated_hap_seqs.front().compare(hap_seqs.front()) == 0);

  // Copy over the alignment probabilities for old sequences present in the new haplotype
  int new_num_alleles         = updated_haplotype->num_combs();
  double* fixed_log_aln_probs = new double[num_reads_*new_num_alleles];
  std::fill_n(fixed_log_aln_probs, num_reads_*new_num_alleles, -100000);
  double* old_log_aln_ptr     = log_aln_probs_;
  double* new_log_aln_ptr     = fixed_log_aln_probs;
  for (unsigned int i = 0; i < num_reads_; ++i){
    for (unsigned int j = 0; j < num_alleles_; ++j, ++old_log_aln_ptr)
      if (allele_mapping[j] != -1){
	assert(!realign_to_haplotype[allele_mapping[j]]);
	new_log_aln_ptr[allele_mapping[j]] = *old_log_aln_ptr;
      }
    new_log_aln_ptr += new_num_alleles;
  }
  delete [] log_aln_probs_;
  log_aln_probs_ = fixed_log_aln_probs;

  // Delete the old haplotype data structures and replace them with the updated ones
  delete haplotype_;
  for (int i = 0; i < hap_blocks_.size(); i++)
    delete hap_blocks_[i];
  haplotype_   = updated_haplotype;
  hap_blocks_  = updated_blocks;
  num_alleles_ = new_num_alleles;

  // Compute the alignment probabilites for new haplotype sequences if there are any
  if (added_seq)
    calc_hap_aln_probs(realign_to_haplotype, realign_pool, copy_read);

  // Fix alignment traceback cache (as allele indices have changed)
  std::map<std::pair<int,int>, AlignmentTrace*> new_trace_cache;
  for (auto cache_iter = trace_cache_.begin(); cache_iter != trace_cache_.end(); cache_iter++){
    int new_allele_index = allele_mapping[cache_iter->first.second];
    if (new_allele_index != -1)
      new_trace_cache[std::pair<int,int>(cache_iter->first.first, new_allele_index)] = cache_iter->second;
    else
      delete cache_iter->second;
  }
  trace_cache_ = new_trace_cache;

  // Resize and recalculate the genotype posterior array
  delete [] log_sample_posteriors_;
  log_sample_posteriors_ = new double[num_samples_*num_alleles_*num_alleles_];
  calc_log_sample_posteriors();
}

void SeqStutterGenotyper::remove_alleles(std::vector< std::vector<int> >& allele_indices){
  std::vector< std::vector<std::string> > alleles_to_add(hap_blocks_.size());
  add_and_remove_alleles(allele_indices, alleles_to_add);
}

bool SeqStutterGenotyper::build_haplotype(std::string& chrom_seq, std::vector<StutterModel*>& stutter_models, std::ostream& logger){
  double locus_hap_build_time = clock();
  assert(hap_blocks_.empty() && haplotype_ == NULL);
  logger << "Generating candidate haplotypes" << std::endl;

  // Determine the minimum and maximum alignment boundaries
  int32_t min_aln_start = INT_MAX, max_aln_stop = INT_MIN;
  for (unsigned int read_index = 0; read_index < num_reads_; read_index++){
    min_aln_start = std::min(min_aln_start, alns_[read_index].get_start());
    max_aln_stop  = std::max(max_aln_stop,  alns_[read_index].get_stop());
  }

  HaplotypeGenerator hap_generator(min_aln_start, max_aln_stop);
  const std::vector<Region>& regions = region_group_->regions();
  bool success = true;
  for (int region_index = 0; region_index < regions.size(); region_index++){
    // Select only those alignments marked as good for haplotype generation
    std::vector<AlnList> gen_hap_alns(num_samples_);
    for (unsigned int read_index = 0; read_index < num_reads_; read_index++)
      if (alns_[read_index].use_for_hap_generation(region_index))
	gen_hap_alns[sample_label_[read_index]].push_back(alns_[read_index]);

    std::vector<std::string> vcf_alleles;
    if (ref_vcf_ != NULL){
      int32_t pos;
      if (!read_vcf_alleles(ref_vcf_, regions[region_index], vcf_alleles, pos)){
	logger << "Haplotype construction failed: The alleles could not be extracted from the reference VCF" << std::endl;
	success = false;
	break;
      }

      // Add the haplotype block based on the extracted VCF alleles
      if (!hap_generator.add_vcf_haplotype_block(pos, chrom_seq, vcf_alleles, stutter_models[region_index])){
	logger << "Haplotype construction failed: " << hap_generator.failure_msg() << std::endl;
	success = false;
	break;
      }
    }
    else {
      // Add the haplotype block in which alleles are derived from the alignments
      if (!hap_generator.add_haplotype_block(regions[region_index], chrom_seq, gen_hap_alns, vcf_alleles, stutter_models[region_index])){
	logger << "Haplotype construction failed: " << hap_generator.failure_msg() << std::endl;
	success = false;
	break;
      }
    }
  }

  if (success){
    if (hap_generator.fuse_haplotype_blocks(chrom_seq)){
      // Copy over the constructed haplotype blocks and build the haplotype
      hap_blocks_  = hap_generator.get_haplotype_blocks();
      haplotype_   = new Haplotype(hap_blocks_);
      num_alleles_ = haplotype_->num_combs();
      call_sample_ = std::vector<bool>(num_samples_, true);
      haplotype_->print_block_structure(30, 100, logger);
    }
    else {
      logger << "Haplotype construction failed: " << hap_generator.failure_msg() << std::endl;
      success = false;
    }
  }

  locus_hap_build_time   = (clock() - locus_hap_build_time)/CLOCKS_PER_SEC;
  total_hap_build_time_ += locus_hap_build_time;
  return success;
}

void SeqStutterGenotyper::init(std::vector<StutterModel*>& stutter_models, std::string& chrom_seq, std::ostream& logger){
  // Allocate and initiate additional data structures
  read_weights_.clear();
  pool_index_   = new int[num_reads_];
  second_mate_  = new bool[num_reads_];
  std::string prev_aln_name = "";

  for (unsigned int read_index = 0; read_index < num_reads_; read_index++){
    pool_index_[read_index]   = pooler_.add_alignment(alns_[read_index]);
    second_mate_[read_index]  = (alns_[read_index].get_name().compare(prev_aln_name) == 0);
    read_weights_.push_back(second_mate_[read_index] ? 0 : 1);
    prev_aln_name = alns_[read_index].get_name();
  }

  initialized_ = build_haplotype(chrom_seq, stutter_models, logger);
  if (initialized_){
    // Allocate the remaining data structures
    log_sample_posteriors_ = new double[num_samples_*num_alleles_*num_alleles_];
    log_aln_probs_         = new double[num_reads_*num_alleles_];
    seed_positions_        = new int[num_reads_];
  }
}

void SeqStutterGenotyper::calc_hap_aln_probs(std::vector<bool>& realign_to_haplotype){
  std::vector<bool> realign_pool = std::vector<bool>(pooler_.num_pools(), true);
  std::vector<bool> copy_read    = std::vector<bool>(num_reads_, true);
  calc_hap_aln_probs(realign_to_haplotype, realign_pool, copy_read);
}

void SeqStutterGenotyper::calc_hap_aln_probs(std::vector<bool>& realign_to_haplotype, std::vector<bool>& realign_pool, std::vector<bool>& copy_read){
  double locus_hap_aln_time = clock();
  assert(haplotype_->num_combs() == realign_to_haplotype.size() && haplotype_->num_combs() == num_alleles_);
  HapAligner hap_aligner(haplotype_, realign_to_haplotype);

  // Align each pooled read to each haplotype
  AlnList& pooled_alns       = pooler_.get_alignments();
  double* log_pool_aln_probs = new double[pooled_alns.size()*num_alleles_];
  int* pool_seed_positions   = new int[pooled_alns.size()];
  hap_aligner.process_reads(pooled_alns, 0, &base_quality_, realign_pool, log_pool_aln_probs, pool_seed_positions);

  // Copy each pool's alignment probabilities to the entries for its constituent reads, but only for realigned haplotypes
  double* log_aln_ptr = log_aln_probs_;
  for (unsigned int i = 0; i < num_reads_; i++){
    if (!copy_read[i]){
      log_aln_ptr += num_alleles_;
      continue;
    }

    seed_positions_[i] = pool_seed_positions[pool_index_[i]];
    double* src_ptr = log_pool_aln_probs + num_alleles_*pool_index_[i];
    for (unsigned int j = 0; j < num_alleles_; ++j, ++log_aln_ptr, ++src_ptr)
      if (realign_to_haplotype[j])
	*log_aln_ptr = *src_ptr;
  }
  delete [] log_pool_aln_probs;
  delete [] pool_seed_positions;

  // If both mate pairs overlap the STR region, they share the same phasing probabilities and we need to avoid treating them as independent
  // To do so, we combine the alignment probabilities here and set the read weight for the second in the pair to zero during the posterior calculation
  // NOTE: It's very important that we don't recombine the values for haplotypes that have already been aligned,
  // or we'll effectively keep doubling those values with each iteration
  for (unsigned int i = 0; i < num_reads_; ++i){
    if (!second_mate_[i] || !copy_read[i])
      continue;

    double* mate_one_ptr = log_aln_probs_ + (i-1)*num_alleles_;
    double* mate_two_ptr = log_aln_probs_ + i*num_alleles_;
    for (unsigned int j = 0; j < num_alleles_; ++j, ++mate_one_ptr, ++mate_two_ptr){
      if (realign_to_haplotype[j]){
	double total  = *mate_one_ptr + *mate_two_ptr;
	*mate_one_ptr = total;
	*mate_two_ptr = total;
      }
    }
  }

  locus_hap_aln_time   = (clock() - locus_hap_aln_time)/CLOCKS_PER_SEC;
  total_hap_aln_time_ += locus_hap_aln_time;
}

bool SeqStutterGenotyper::id_and_align_to_stutter_alleles(std::string& chrom_seq, std::ostream& logger){
  std::vector< std::vector<int> > alleles_to_remove(haplotype_->num_blocks());
  while (true){
    // Look for candidate alleles present in stutter artifacts
    bool added_alleles = false;
    std::vector< std::vector<std::string> > stutter_seqs(haplotype_->num_blocks());
    for (int i = 0; i < haplotype_->num_blocks(); i++){
      HapBlock* block = haplotype_->get_block(i);
      if (block->get_repeat_info() != NULL){
	get_stutter_candidate_alleles(i, logger, stutter_seqs[i]);
	for (auto allele_iter = stutter_seqs[i].begin(); allele_iter != stutter_seqs[i].end(); allele_iter++)
	  if (allele_iter->size() < std::abs(block->get_repeat_info()->max_deletion()))
	    return false;
	added_alleles |= !stutter_seqs[i].empty();
	std::sort(stutter_seqs[i].begin(), stutter_seqs[i].end(), orderByLengthAndSequence);
      }
    }
    // Terminate if no new alleles identified in any of the blocks
    if (!added_alleles) break;

    // Otherwise, add the new alleles to the haplotype and recompute the relevant values
    add_and_remove_alleles(alleles_to_remove, stutter_seqs);
  }
  return true;
}

bool SeqStutterGenotyper::genotype(std::string& chrom_seq, std::ostream& logger){
  // Unsuccessful initialization. May be due to
  // 1) Failing to find the corresponding alleles in the VCF (if one has been provided)
  // 2) Large deletion extending past STR
  if (!initialized_)
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

  // Check if we can assemble the sequences flanking the STR
  // If not, it's likely that the flanks are too repetitive and will introduce genotyping errors
  for (int flank = 0; flank < 2; flank++){
    int block_index     = (flank == 0 ? 0 : haplotype_->num_blocks()-1);
    std::string ref_seq = hap_blocks_[block_index]->get_seq(0);
    int max_k           = std::min(MAX_KMER, ref_seq.size() == 0 ? -1 : (int)ref_seq.size()-1);
    int kmer_length;
    if (!DebruijnGraph::calc_kmer_length(ref_seq, MIN_KMER, max_k, kmer_length)){
      logger << "Aborting genotyping of the locus as the sequence " << (flank == 0 ? "upstream" : "downstream")
	     << " of the repeat is too repetitive for accurate genotyping" << "\n";
      logger << "\tFlanking sequence = " << ref_seq << std::endl;
      return false;
    }
  }

  init_alignment_model();
  logger << "Pooling reads with identical sequences" << std::endl;
  pooler_.pool(base_quality_);

  // Align each read to each candidate haplotype and store them in the provided arrays
  logger << "Aligning reads to each candidate haplotype" << std::endl;
  std::vector<bool> realign_to_haplotype(num_alleles_, true);
  assert(realign_to_haplotype.size() == haplotype_->num_combs());
  calc_hap_aln_probs(realign_to_haplotype);
  calc_log_sample_posteriors();

  if (ref_vcf_ == NULL){
    // Look for additional alleles in stutter artifacts and align to them (if necessary)
    if (!id_and_align_to_stutter_alleles(chrom_seq, logger))
      return false;

    // Remove alleles with no MAP genotype calls and recompute the posteriors
    std::vector< std::vector<int> > unused_indices;
    int num_aff_blocks = 0, num_aff_alleles = 0;
    get_unused_alleles(false, true, unused_indices, num_aff_blocks, num_aff_alleles);
    if (num_aff_alleles != 0){
      logger << "Recomputing sample posteriors after removing " << num_aff_alleles
	     << " uncalled alleles across " << num_aff_blocks << " blocks" << std::endl;
      remove_alleles(unused_indices);
    }

    // Remove alleles with no spanning reads and recompute the posteriors
    unused_indices.clear();
    get_unused_alleles(true, false, unused_indices, num_aff_blocks, num_aff_alleles);
    if (num_aff_alleles != 0){
      logger << "Recomputing sample posteriors after removing " << num_aff_alleles
	     << " alleles with no spanning reads across " << num_aff_blocks << " blocks" << std::endl;
      remove_alleles(unused_indices);
    }
  }

  if (reassemble_flanks_)
    if (!assemble_flanks(logger))
      return false;

  return true;
}

void SeqStutterGenotyper::reorder_alleles(std::vector<std::string>& alleles,
					  std::vector<int>& old_to_new, std::vector<int>& new_to_old){
  assert(old_to_new.empty() && new_to_old.empty());
  std::vector< std::pair<std::string, int> > tuples;
  for (int i = 0; i < alleles.size(); i++)
    tuples.push_back(std::pair<std::string, int>(alleles[i], i));
  std::sort(tuples.begin()+1, tuples.end(), [](const std::pair<std::string, int>& left, const std::pair<std::string, int>& right) { 
      return orderByLengthAndSequence(left.first, right.first);});

  old_to_new = std::vector<int>(alleles.size(), -1);
  for (int i = 0; i < tuples.size(); i++){
    new_to_old.push_back(tuples[i].second);
    old_to_new[tuples[i].second] = i;
  }
}

void SeqStutterGenotyper::get_alleles(const Region& region, int block_index, std::string& chrom_seq,
				      int32_t& pos, std::vector<std::string>& alleles){
  assert(alleles.size() == 0);

  // Extract all the alleles
  HapBlock* block = haplotype_->get_block(block_index);
  for (int i = 0; i < block->num_options(); i++)
    alleles.push_back(block->get_seq(i));

  // Trim from the left until the region boundary or a mismatched character
  int32_t left_trim = 0;
  int32_t start     = block->start();
  while (start + left_trim < region.start()){
    bool trim = true;
    for (unsigned int i = 0; i < alleles.size(); ++i){
      if ((left_trim+1 >= alleles[i].size()) || (alleles[i][left_trim] != alleles[0][left_trim])){
	trim = false;
	break;
      }
    }
    if (!trim) break;
    left_trim++;
  }
  start += left_trim;
  for (unsigned int i = 0; i < alleles.size(); ++i)
    alleles[i] = alleles[i].substr(left_trim);

  // Trim from the right until the region boundary or a mismatched character
  int32_t right_trim = 0;
  int32_t end        = block->end();
  while (end - right_trim > region.stop()){
    bool trim    = true;
    int ref_size = alleles[0].size();
    for (unsigned int i = 0; i < alleles.size(); ++i){
      int alt_size = alleles[i].size();
      if ((right_trim+1 >= alleles[i].size()) || (alleles[i][alt_size-right_trim-1] != alleles[0][ref_size-right_trim-1])){
	trim = false;
	break;
      }
    }
    if (!trim) break;
    right_trim++;
  }
  end -= right_trim;
  for (unsigned int i = 0; i < alleles.size(); ++i)
    alleles[i] = alleles[i].substr(0, alleles[i].size()-right_trim);

  std::string left_flank  = (start >= region.start() ? uppercase(chrom_seq.substr(region.start(), start-region.start())) : "");
  std::string right_flank = (end <= region.stop()    ? uppercase(chrom_seq.substr(end, region.stop()-end)) : "");
  pos = std::min((int32_t)region.start(), start);

  // If necessary, add 1bp on the left so that all the alleles match the reference sequence
  if (left_flank.empty()){
    bool pad_left = false;
    for (unsigned int i = 1; i < alleles.size(); ++i){
      if (alleles[i][0] != alleles[0][0]){
	pad_left = true;
	break;
      }
    }
    if (pad_left){
      pos -= 1;
      left_flank = uppercase(chrom_seq.substr(pos, 1));
    }
  }

  for (unsigned int i = 0; i < alleles.size(); ++i){
    std::stringstream ss;
    ss << left_flank << alleles[i] << right_flank;
    alleles[i] = ss.str();
  }

  pos += 1; // Fix off-by-1 VCF error
}
 
void SeqStutterGenotyper::debug_sample(int sample_index, std::ostream& logger){
  std::vector<AlignmentTrace*> traced_alns;
  retrace_alignments(traced_alns);

  std::cerr << "DEBUGGING SAMPLE..." << std::endl;
  std::cerr << "READ LL's:" << std::endl;
  double* read_LL_ptr = log_aln_probs_;
  int read_index = 0;
  for (unsigned int i = 0; i < num_reads_; ++i){
    if(sample_label_[i] == sample_index && seed_positions_[i] >= 0){
      std::cerr << "\t" << "READ #" << i << ", SEED BASE=" << seed_positions_[i] << ", POOL INDEX=" << pool_index_[i] << ", IS_SECOND_MATE=" << second_mate_[i]
		<< ", TOTAL QUAL CORRECT= " << alns_[i].sum_log_prob_correct(base_quality_) << ", "
		<< max_index(read_LL_ptr, num_alleles_) << ", "
		<< log_p1_[i] << " " << log_p2_[i] <<  ", "
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
  double* sample_LL_ptr = log_sample_posteriors_ + num_alleles_*num_alleles_*sample_index;
  for (int index_1 = 0; index_1 < num_alleles_; ++index_1)
    for (int index_2 = 0; index_2 < num_alleles_; ++index_2, ++sample_LL_ptr)
      std::cerr << index_1 << " " << index_2 << " " << *sample_LL_ptr << "(" << exp(*sample_LL_ptr) << ")" << std::endl;
   std::cerr << "END OF SAMPLE DEBUGGING" << std::endl;
}

void SeqStutterGenotyper::retrace_alignments(std::vector<AlignmentTrace*>& traced_alns){
  assert(traced_alns.size() == 0);
  double trace_start = clock();
  traced_alns.reserve(num_reads_);
  std::vector< std::pair<int, int> > haps;
  get_optimal_haplotypes(haps);

  std::vector<bool> realign_to_haplotype(num_alleles_, true);
  HapAligner hap_aligner(haplotype_, realign_to_haplotype);
  double* read_LL_ptr = log_aln_probs_;
  for (unsigned int read_index = 0; read_index < num_reads_; read_index++){
    if (seed_positions_[read_index] < 0){
      read_LL_ptr += num_alleles_;
      traced_alns.push_back(NULL);
      continue;
    }

    int hap_a    = haps[sample_label_[read_index]].first;
    int hap_b    = haps[sample_label_[read_index]].second;
    int best_hap = ((LOG_ONE_HALF+log_p1_[read_index]+read_LL_ptr[hap_a] > LOG_ONE_HALF+log_p2_[read_index]+read_LL_ptr[hap_b]) ? hap_a : hap_b);

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

void SeqStutterGenotyper::get_stutter_candidate_alleles(int str_block_index, std::ostream& logger, std::vector<std::string>& candidate_seqs){
  assert(candidate_seqs.size() == 0 && haplotype_->get_block(str_block_index)->get_repeat_info() != NULL);
  HapBlock* str_block = haplotype_->get_block(str_block_index);
  std::vector<AlignmentTrace*> traced_alns;
  retrace_alignments(traced_alns);

  // Count the number of reads spanning the STR and the number of reads with stutter
  std::vector<int> sample_counts(num_samples_, 0);
  std::vector< std::map<std::string, int> > sample_stutter_counts(num_samples_);
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

  // Add frequently observed stutter artifacts as candidate sequences
  std::set<std::string> candidate_set;
  for (unsigned int i = 0; i < num_samples_; i++)
    for (auto seq_iter = sample_stutter_counts[i].begin(); seq_iter != sample_stutter_counts[i].end(); seq_iter++)
      if (seq_iter->second >= 2 && 1.0*seq_iter->second/sample_counts[i] >= 0.15)
	if (!str_block->contains(seq_iter->first))
	  candidate_set.insert(seq_iter->first);
  candidate_seqs = std::vector<std::string>(candidate_set.begin(), candidate_set.end());

  if (candidate_seqs.size() != 0){
    logger << "Identified " << candidate_seqs.size() << " additional candidate alleles from stutter artifacts" << "\n";
    for (unsigned int i = 0; i < candidate_seqs.size(); i++)
      logger << "\t" << candidate_seqs[i] << "\n";
  }
}

/*
void SeqStutterGenotyper::analyze_flank_snps(std::ostream& logger){
  std::vector<AlignmentTrace*> traced_alns;
  std::vector< std::pair<int, int> > haps;
  retrace_alignments(traced_alns);
  get_optimal_haplotypes(haps);

  double* read_LL_ptr = log_aln_probs_;
  int min_read_index = 0, read_index;
  for (int sample_index = 0; sample_index < num_samples_; sample_index++){
    int hap_a = haps[sample_index].first;
    int hap_b = haps[sample_index].second;
    std::map<std::pair<int32_t, char>, int> hap_snp_counts;
    std::vector<std::pair<int32_t, int32_t> > hap_read_coords;

    for (read_index = min_read_index; read_index < num_reads_; read_index++){
      if (traced_alns[read_index] == NULL){
	read_LL_ptr += num_alleles_;
	continue;
      }
      if (sample_label_[read_index] != sample_index)
	break;

      AlignmentTrace* trace = traced_alns[read_index];
      hap_read_coords.push_back(std::pair<int32_t,int32_t>(trace->traced_aln().get_start(), trace->traced_aln().get_stop()));
      std::vector<std::pair<int32_t, char> >& snp_data = trace->flank_snp_data();
      for (auto snp_iter = snp_data.begin(); snp_iter != snp_data.end(); snp_iter++)
	hap_snp_counts[*snp_iter]++;
      read_LL_ptr += num_alleles_;
    }
    min_read_index = read_index;

    for (auto snp_iter = hap_snp_counts.begin(); snp_iter != hap_snp_counts.end(); snp_iter++){
      if (snp_iter->second >= 2){
	int32_t pos = snp_iter->first.first;
	int num_span = 0;
	for (int i = 0; i < hap_read_coords.size(); i++)
	  if (hap_read_coords[i].first < pos && hap_read_coords[i].second > pos)
	    num_span++;
	const Region& region = region_group_->regions()[0];
	std::cerr << "SNP INFO " << region.chrom() << " " << region.start()+1 << " " << sample_names_[sample_index] << " " << pos << " " << snp_iter->second << " " << num_span << std::endl;
      }
    }
  }
}


void SeqStutterGenotyper::analyze_flank_indels(std::ostream& logger){
  std::vector<AlignmentTrace*> traced_alns;
  retrace_alignments(traced_alns);
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
    for (auto candidate_iter = candidate_set.begin(); candidate_iter != candidate_set.end(); candidate_iter++)
      std::cerr << "FLANK INDEL " << region_group_->chrom() << " " << candidate_iter->first.first << " " << candidate_iter->first.second << " " << candidate_iter->second << std::endl;
  }
}
*/

double SeqStutterGenotyper::compute_allele_bias(int hap_a_read_count, int hap_b_read_count){
  // Compute p-value for allele read depth bias according to the counts of reads uniquely assigned to each haplotype
  // We use the bdtr(k, N, p) function from the cephes directory, which computes the CDF for a binomial distribution
  // e.g.: double val = bdtr (24, 50, 0.5);
  int total = hap_a_read_count + hap_b_read_count;

  // Not applicable
  if (total == 0)
    return 1;

  // p-value is 1
  if (hap_a_read_count == hap_b_read_count)
    return 0.0;

  int min_count = std::min(hap_a_read_count, hap_b_read_count);
  double pvalue = 2*bdtr(min_count, total, 0.5); // Two-sided pvalue
  return log10(std::min(1.0, pvalue));
}

void SeqStutterGenotyper::write_vcf_record(std::vector<std::string>& sample_names, std::string& chrom_seq,
					   bool output_gls, bool output_pls, bool output_phased_gls, bool output_allreads,
					   bool output_mallreads, bool output_viz, float max_flank_indel_frac, bool viz_left_alns,
                                           std::ostream& html_output, std::ostream& out, std::ostream& logger){
  int region_index = 0;
  for (int block_index = 0; block_index < haplotype_->num_blocks(); block_index++)
    if (haplotype_->get_block(block_index)->get_repeat_info() != NULL)
      write_vcf_record(sample_names, block_index, region_group_->regions()[region_index++], chrom_seq, output_gls, output_pls, output_phased_gls,
		       output_allreads, output_mallreads, output_viz, max_flank_indel_frac, viz_left_alns, html_output, out, logger);
  assert(region_index == region_group_->num_regions());
}

void SeqStutterGenotyper::write_vcf_record(std::vector<std::string>& sample_names, int hap_block_index, const Region& region, std::string& chrom_seq, bool output_gls,
					   bool output_pls, bool output_phased_gls, bool output_allreads, bool output_mallreads,
					   bool output_viz, float max_flank_indel_frac, bool viz_left_alns,
					   std::ostream& html_output, std::ostream& out, std::ostream& logger){
  // Extract the alleles and position for the current haplotype block
  int32_t pos;
  std::vector<std::string> alleles;
  get_alleles(region, hap_block_index, chrom_seq, pos, alleles);

  //debug_sample(sample_indices_["LP6005442-DNA_D08"], logger);
  //std::vector<bool> clobbered;
  //haplotype_->check_indel_clobbering(region.name(), clobbered);

  // Compute the base pair differences from the reference
  std::vector<int> allele_bp_diffs;
  for (unsigned int i = 0; i < alleles.size(); i++)
    allele_bp_diffs.push_back((int)alleles[i].size() - (int)alleles[0].size());

  // Extract the optimal genotypes and their associated likelihoods
  std::vector< std::pair<int,int> > haplotypes, gts;
  std::vector<double> log_phased_posteriors, log_unphased_posteriors, gl_diffs;
  std::vector< std::vector<double> > gls, phased_gls;
  std::vector< std::vector<int> > pls;
  std::vector<int> hap_to_allele;
  haps_to_alleles(hap_block_index, hap_to_allele);
  int num_variants = haplotype_->get_block(hap_block_index)->num_options();
  extract_genotypes_and_likelihoods(num_variants, hap_to_allele, haplotypes, gts, log_phased_posteriors, log_unphased_posteriors,
				    true, gls, gl_diffs, output_pls, pls, output_phased_gls, phased_gls);

  // Extract information about each read and group by sample
  std::vector<int> num_aligned_reads(num_samples_, 0), num_reads_with_snps(num_samples_, 0);
  std::vector<int> num_reads_with_stutter(num_samples_, 0), num_reads_with_flank_indels(num_samples_, 0);
  std::vector<int> num_reads_strand_one(num_samples_, 0), num_reads_strand_two(num_samples_, 0);
  std::vector<int> unique_reads_hap_one(num_samples_, 0), unique_reads_hap_two(num_samples_, 0);
  std::vector< std::vector<int> > bps_per_sample(num_samples_), ml_bps_per_sample(num_samples_);
  std::vector< std::vector<double> > log_read_phases(num_samples_);
  std::vector<AlnList> max_LL_alns_strand_one(num_samples_), left_alns_strand_one(num_samples_), orig_alns_strand_one(num_samples_);
  std::vector<AlnList> max_LL_alns_strand_two(num_samples_), left_alns_strand_two(num_samples_), orig_alns_strand_two(num_samples_);
  std::vector<bool> realign_to_haplotype(num_alleles_, true);
  HapAligner hap_aligner(haplotype_, realign_to_haplotype);
  double* read_LL_ptr = log_aln_probs_;
  int bp_diff; bool got_size;
  for (unsigned int read_index = 0; read_index < num_reads_; read_index++){
    if (seed_positions_[read_index] < 0){
      read_LL_ptr += num_alleles_;
      continue;
    }

    // Extract read's phase posterior conditioned on the determined sample genotype
    int hap_a            = haplotypes[sample_label_[read_index]].first;
    int hap_b            = haplotypes[sample_label_[read_index]].second;
    double total_read_LL = log_sum_exp(LOG_ONE_HALF+log_p1_[read_index]+read_LL_ptr[hap_a], LOG_ONE_HALF+log_p2_[read_index]+read_LL_ptr[hap_b]);
    double log_phase_one = LOG_ONE_HALF + log_p1_[read_index] + read_LL_ptr[hap_a] - total_read_LL; 
    log_read_phases[sample_label_[read_index]].push_back(log_phase_one);

    // Determine which of the two genotypes each read is associated with
    int read_strand = 0;
    if (!haploid_ && ((hap_a != hap_b) || (abs(log_p1_[read_index]-log_p2_[read_index]) > TOLERANCE))){
      double v1 = log_p1_[read_index]+read_LL_ptr[hap_a], v2 = log_p2_[read_index]+read_LL_ptr[hap_b];
      if (abs(v1-v2) > TOLERANCE){
	read_strand = (v1 > v2 ? 0 : 1);
	if (read_strand == 0)
	  unique_reads_hap_one[sample_label_[read_index]]++;
	else
	  unique_reads_hap_two[sample_label_[read_index]]++;
      }
    }

    // Retrace alignment and ensure that it's of sufficient quality
    double trace_start = clock();
    int best_hap = (read_strand == 0 ? hap_a : hap_b);
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

    if (viz_left_alns)
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
    got_size = ExtractCigar(alns_[read_index].get_cigar_list(), alns_[read_index].get_start(), region.start()-region.period(), region.stop()+region.period(), bp_diff);
    if (got_size) bps_per_sample[sample_label_[read_index]].push_back(bp_diff);

    // Extract the ML bp difference observed in read based on the ML genotype,
    // but only for reads that span the original repeat region by 5 bp
    if (trace->traced_aln().get_start() < (region.start() > 4 ? region.start()-4 : 0))
      if (trace->traced_aln().get_stop() > region.stop() + 4)
	ml_bps_per_sample[sample_label_[read_index]].push_back(allele_bp_diffs[hap_to_allele[best_hap]]+trace->total_stutter_size());

    read_LL_ptr += num_alleles_;
  }
 
  // Compute allele counts for samples of interest
  std::set<std::string> samples_of_interest(sample_names.begin(), sample_names.end());
  std::vector<int> allele_counts(alleles.size());
  int sample_index = 0, skip_count = 0, filt_count = 0, allele_number = 0;
  for (auto gt_iter = gts.begin(); gt_iter != gts.end(); ++gt_iter, ++sample_index){
    if (samples_of_interest.find(sample_names_[sample_index]) == samples_of_interest.end())
      continue;
    if (num_aligned_reads[sample_index] == 0)
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

  // Determine an ordering for the alleles, in which they're sorted by allele length
  std::vector<int> old_to_new, new_to_old;
  reorder_alleles(alleles, old_to_new, new_to_old);

  // Print the allele count information
  logger << "Allele counts" << std::endl;
  for (unsigned int i = 0; i < alleles.size(); i++)
    logger << alleles[new_to_old[i]] << " " << allele_counts[new_to_old[i]] <<  std::endl;
  logger << std::endl;

  //VCF line format = CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE_1 SAMPLE_2 ... SAMPLE_N
  out << region.chrom() << "\t" << pos << "\t" << (region.name().empty() ? "." : region.name());

  // Add reference allele and alternate alleles
  out << "\t" << alleles[new_to_old[0]] << "\t";
  if (alleles.size() == 1)
    out << ".";
  else {
    for (int i = 1; i < alleles.size()-1; i++)
      out << alleles[new_to_old[i]] << ",";
    out << alleles[new_to_old.back()];
  }

  // Add QUAL and FILTER fields
  out << "\t" << "." << "\t" << ".";

  // Obtain relevant stutter model
  assert(haplotype_->get_block(hap_block_index)->get_repeat_info() != NULL);
  StutterModel* stutter_model = haplotype_->get_block(hap_block_index)->get_repeat_info()->get_stutter_model();

  // Add INFO field items
  out << "\tINFRAME_PGEOM=" << stutter_model->get_parameter(true,  'P') << ";"
      << "INFRAME_UP="      << stutter_model->get_parameter(true,  'U') << ";"
      << "INFRAME_DOWN="    << stutter_model->get_parameter(true,  'D') << ";"
      << "OUTFRAME_PGEOM="  << stutter_model->get_parameter(false, 'P') << ";"
      << "OUTFRAME_UP="     << stutter_model->get_parameter(false, 'U') << ";"
      << "OUTFRAME_DOWN="   << stutter_model->get_parameter(false, 'D') << ";"
      << "START="           << region.start()+1 << ";"
      << "END="             << region.stop()    << ";"
      << "PERIOD="          << region.period()  << ";"
      << "NSKIP="           << skip_count       << ";"
      << "NFILT="           << filt_count       << ";";
  if (alleles.size() > 1){
    out << "BPDIFFS=" << allele_bp_diffs[new_to_old[1]];
    for (unsigned int i = 2; i < alleles.size(); i++)
      out << "," << allele_bp_diffs[new_to_old[i]];
    out << ";";
  }

  // Compute INFO field values for DP, DSTUTTER and DFLANKINDEL and add them to the VCF
  int32_t tot_dp = 0, tot_dsnp = 0, tot_dstutter = 0, tot_dflankindel = 0;
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
    tot_dstutter    += num_reads_with_stutter[sample_index];
    tot_dflankindel += num_reads_with_flank_indels[sample_index];
  }
  out << "DP="          << tot_dp          << ";"
      << "DSNP="        << tot_dsnp        << ";"
      << "DSTUTTER="    << tot_dstutter    << ";"
      << "DFLANKINDEL=" << tot_dflankindel << ";";

  // Add allele counts
  out << "AN=" << allele_number << ";" << "REFAC=" << allele_counts[0];
  if (allele_counts.size() > 1){
    out << ";AC=";
    for (unsigned int i = 1; i < allele_counts.size()-1; i++)
      out << allele_counts[new_to_old[i]] << ",";
    out << allele_counts[new_to_old.back()];
  }

  // If we used all reads during genotyping and performed assembly, we'll output the allele bias
  bool output_allele_bias = (!haploid_ && reassemble_flanks_);

  // Add FORMAT field
  out << (!haploid_ ? "\tGT:GB:Q:PQ:DP:DSNP:DSTUTTER:DFLANKINDEL:PDP:PSNP:GLDIFF" : "\tGT:GB:Q:DP:DSTUTTER:DFLANKINDEL:GLDIFF");
  if (output_allele_bias)         out << ":AB:DAB";
  if (output_allreads)            out << ":ALLREADS";
  if (output_mallreads)           out << ":MALLREADS";
  if (output_gls)                 out << ":GL";
  if (output_pls)                 out << ":PL";
  if (!haploid_ && output_phased_gls) out << ":PHASEDGL";

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
    if (num_aligned_reads[sample_iter->second] == 0){
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

    double allele_bias = 1;
    if (!haploid_ && (gts[sample_index].first != gts[sample_index].second))
      allele_bias = compute_allele_bias(unique_reads_hap_one[sample_index], unique_reads_hap_two[sample_index]);

    if (!haploid_){
      out << old_to_new[gts[sample_index].first] << "|" << old_to_new[gts[sample_index].second]     // Genotype
	  << ":" << allele_bp_diffs[gts[sample_index].first]
	  << "|" << allele_bp_diffs[gts[sample_index].second]                                       // Base pair differences from reference
	  << ":" << exp(log_unphased_posteriors[sample_index])                                      // Unphased posterior
	  << ":" << exp(log_phased_posteriors[sample_index])                                        // Phased posterior
	  << ":" << num_aligned_reads[sample_index]                                                 // Total reads used to genotype (after filtering)
	  << ":" << num_reads_with_snps[sample_index]                                               // Total reads with SNP information
	  << ":" << num_reads_with_stutter[sample_index]                                            // Total reads with a non-zero stutter artifact in ML alignment
	  << ":" << num_reads_with_flank_indels[sample_index]                                       // Total reads with an indel in flank in ML alignment
	  << ":" << phase1_reads << "|" << phase2_reads                                             // Reads per allele
	  << ":" << num_reads_strand_one[sample_index] << "|" << num_reads_strand_two[sample_index]; // Reads with SNPs supporting each haploid genotype

      // Difference in GL between the current and next best genotype
      if (alleles.size() == 1)
        out << ":" << ".";
      else
        out << ":" << gl_diffs[sample_index];
    }
    else {
      out << old_to_new[gts[sample_index].first]                                                    // Genotype
	  << ":" << allele_bp_diffs[gts[sample_index].first]                                        // Base pair differences from reference
	  << ":" << exp(log_unphased_posteriors[sample_index])                                      // Unphased posterior
	  << ":" << num_aligned_reads[sample_index]                                                 // Total reads used to genotype (after filtering)
	  << ":" << num_reads_with_stutter[sample_index]                                            // Total reads with a non-zero stutter artifact in ML alignment
	  << ":" << num_reads_with_flank_indels[sample_index];                                      // Total reads with an indel in flank in ML alignment

      // Difference in GL between the current and next best genotype
      if (alleles.size() == 1)
	out << ":" << ".";
      else
	out << ":" << gl_diffs[sample_index];
    }

    // Output the log-10 value of the allele bias p-value
    if (output_allele_bias){
      if (std::abs(allele_bias-1) < TOLERANCE)
	out << ":.:.";
      else
	out << ":" << allele_bias << ":" << (unique_reads_hap_one[sample_index] + unique_reads_hap_two[sample_index]);
    }

    // Add bp diffs from regular left-alignment
    if (output_allreads)
	out << ":" << condense_read_counts(bps_per_sample[sample_index]);

    // Maximum likelihood base pair differences in each read from alignment probabilites
    if (output_mallreads)
	out << ":" << condense_read_counts(ml_bps_per_sample[sample_index]);

    // Genotype and phred-scaled likelihoods, taking into account new allele ordering
    if (haploid_){
      if (output_gls){
	out << ":" << gls[sample_index][0];
	for (int i = 1; i < new_to_old.size(); i++)
	  out << "," << gls[sample_index][new_to_old[i]];
      }

      if (output_pls){
	out << ":" << pls[sample_index][0];
	for (int i = 1; i < new_to_old.size(); i++)
	  out << "," << pls[sample_index][new_to_old[i]];
      }
    }
    else {
      if (output_gls){
	out << ":" << gls[sample_index][0];
	for (int i = 1; i < new_to_old.size(); i++){
	  for (int j = 0; j <= i; j++){
	    int index_a = std::min(new_to_old[i], new_to_old[j]);
	    int index_b = std::max(new_to_old[i], new_to_old[j]);
	    out << "," << gls[sample_index][index_b*(index_b+1)/2 + index_a];
	  }
	}
      }

      if (output_pls){
	out << ":" << pls[sample_index][0];
	for (int i = 1; i < new_to_old.size(); i++){
	  for (int j = 0; j <= i; j++){
	    int index_a = std::min(new_to_old[i], new_to_old[j]);
	    int index_b = std::max(new_to_old[i], new_to_old[j]);
	    out << "," << pls[sample_index][index_b*(index_b+1)/2 + index_a];
	  }
	}
      }

      if (output_phased_gls){
	out << ":" << phased_gls[sample_index][0];
	for (int i = 0; i < new_to_old.size(); i++){
	  for (int j = 0; j < new_to_old.size(); j++){
	    if (i == 0 && j == 0)
	      continue;
	    out << "," << phased_gls[sample_index][new_to_old[i]*new_to_old.size() + new_to_old[j]];
	  }
	}
      }
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
    locus_info << region.chrom() << "\t" << region.start()+1 << "\t" << region.stop();
    visualizeAlignments(max_LL_alns, sample_names_, sample_results, hap_blocks_, chrom_seq, locus_info.str(), true, html_output);
  }
}

bool SeqStutterGenotyper::recompute_stutter_models(std::string& chrom_seq, std::ostream& logger,
						  int max_em_iter, double abs_ll_converge, double frac_ll_converge){
  logger << "Retraining EM stutter genotyper using maximum likelihood alignments" << std::endl;
  std::vector<AlignmentTrace*> traced_alns;
  retrace_alignments(traced_alns);

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
    EMStutterGenotyper length_genotyper(haploid_, period, str_num_bps, str_log_p1s, str_log_p2s, sample_names_, 0);
    bool trained = length_genotyper.train(max_em_iter, abs_ll_converge, frac_ll_converge, false, logger);
    if (!trained){
      logger << "Retraining stutter model training failed" << std::endl;
      return false;
    }

    logger << "Learned stutter model for block #" << block_index << ":" << (*length_genotyper.get_stutter_model()) << std::endl;
    block->get_repeat_info()->set_stutter_model(length_genotyper.get_stutter_model());
  }
  trace_cache_.clear();
  return genotype(chrom_seq, logger);
}
