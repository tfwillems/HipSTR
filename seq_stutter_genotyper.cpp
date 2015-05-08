#include <algorithm>
#include <cfloat>
#include <sstream>

#include "seq_stutter_genotyper.h"
#include "error.h"
#include "extract_indels.h"
#include "mathops.h"
#include "stringops.h"

#include "SeqAlignment/AlignmentOps.h"
#include "SeqAlignment/AlignmentData.h"
#include "SeqAlignment/AlignmentModel.h"
#include "SeqAlignment/HaplotypeGenerator.h"
#include "SeqAlignment/HapAligner.h"
#include "SeqAlignment/RepeatStutterInfo.h"

void SeqStutterGenotyper::read_ref_vcf_alleles(std::vector<std::string>& alleles){
    assert(alleles.size() == 0);
    assert(ref_vcf_ != NULL);
    if (!ref_vcf_->setRegion(region_->chrom(), region_->start(), region_->stop())){
      // Retry setting region if chr is in chromosome name
      if (region_->chrom().size() <= 3 || region_->chrom().substr(0, 3).compare("chr") != 0 
	  || !ref_vcf_->setRegion(region_->chrom().substr(3), region_->start(), region_->stop()))
	printErrorAndDie("Failed to set VCF region when reading candidate STR alleles");
    }                                                                                                                                                                                  
    // Extract STR and ensure the coordinates match
    vcf::Variant variant(*ref_vcf_);
    while (ref_vcf_->getNextVariant(variant)){
      if (variant.position < region_->start())
	continue;
      else if (variant.position == region_->start()){
	int32_t end = (int32_t)variant.getInfoValueFloat(END_KEY);
	if (end == region_->stop()) {
	  for (auto iter = variant.alleles.begin(); iter != variant.alleles.end(); iter++)
	    alleles.push_back(*iter);
	  return;
	}
      }
      else 
	break;
    }
    printErrorAndDie("Failed to extract matching VCF entry when reading candidate STR alleles");
}

void SeqStutterGenotyper::init(std::vector< std::vector<BamTools::BamAlignment> >& alignments, 
			       std::vector< std::vector<double> >& log_p1, 
			       std::vector< std::vector<double> >& log_p2,
			       std::vector<std::string>& sample_names,
			       std::string& chrom_seq){
  // Compute the total number of reads
  num_reads_ = 0;
  for (unsigned int i = 0; i < alignments.size(); ++i)
    num_reads_ += alignments[i].size();

  // Allocate some data structures
  log_p1_       = new double[num_reads_];
  log_p2_       = new double[num_reads_];
  sample_label_ = new int[num_reads_];
  
  // Extract base pair differences (for debugging purposes)
  int bp_diff;
  for (unsigned int i = 0; i < alignments.size(); ++i){
    for (unsigned int j = 0; j < alignments[i].size(); ++j){
      bool got_size = ExtractCigar(alignments[i][j].CigarData, alignments[i][j].Position, region_->start()-region_->period(), region_->stop()+region_->period(), bp_diff);
      if (got_size)
	bp_diffs_.push_back(bp_diff);
      else
	bp_diffs_.push_back(-999);
    }
  }

  std::cerr << "Left aligning reads..." << std::endl;
  std::map<std::string, Alignment*> seq_to_alns;
  int read_index = 0;
  // TO DO: Test memoized vs. non-memoized speed
  for (unsigned int i = 0; i < alignments.size(); ++i){
    alns_.push_back(std::vector<Alignment>());
    for (unsigned int j = 0; j < alignments[i].size(); ++j, ++read_index){
      auto iter = seq_to_alns.find(alignments[i][j].QueryBases);
      if (iter == seq_to_alns.end()){
	alns_.back().push_back(Alignment());
	realign(alignments[i][j], chrom_seq, alns_.back().back());
      }
      else {
	alns_.back().push_back(*(iter->second));
	seq_to_alns[alignments[i][j].QueryBases] = &(alns_.back().back());
      }
      log_p1_[read_index]       = log_p1[i][j];
      log_p2_[read_index]       = log_p2[i][j];
      sample_label_[read_index] = i; 
    }
  }

  std::vector<std::string> vcf_alleles;
  if (ref_vcf_ != NULL){
    std::cerr << "Reading STR alleles from reference VCF" << std::endl;
    read_ref_vcf_alleles(vcf_alleles);
  }

  // Generate putative haplotypes and determine the number of alleles
  std::cerr << "Generating putative haplotypes..." << std::endl;
  haplotype_   = generate_haplotype(*region_, MAX_REF_FLANK_LEN, chrom_seq, alns_, vcf_alleles, stutter_model_, alleles_from_bams_, hap_blocks_, call_sample_);
  num_alleles_ = haplotype_->num_combs();
  assert(call_sample_.size() == num_samples_);

  // Print information about the haplotype and the stutter model 
  std::cerr << "Max block sizes: ";
  for (unsigned int i = 0; i < haplotype_->num_blocks(); i++)
    std::cerr << haplotype_->get_block(i)->max_size() << " ";
  std::cerr << std::endl << "Stutter model information" << std::endl;
  RepeatStutterInfo* stutter_info = hap_blocks_[1]->get_repeat_info();
  for (int i = stutter_info->max_deletion(); i <= stutter_info->max_insertion(); i += stutter_info->get_period())
    std::cerr << i << " " << stutter_info->log_prob_pcr_artifact(1, i) << std::endl;
  std::cerr << std::endl;

  // Allocate the remaining data structures
  log_sample_posteriors_ = new double[num_alleles_*num_alleles_*num_samples_];
  log_aln_probs_         = new double[num_reads_*num_alleles_];
  seed_positions_        = new int[num_reads_];

  // Extract full STR sequence for each allele using annotated repeat region
  // and the haplotype above
  get_alleles(chrom_seq, alleles_);
}

bool SeqStutterGenotyper::genotype(){
  // Failed to extract alleles (likely due to large deletion extending past STR)
  if (pos_ == -1)
    return false;

  // Align each read against each candidate haplotype
  std::cerr << "Aligning reads to each candidate haplotype..." << std::endl;
  init_alignment_model();
  HapAligner hap_aligner(haplotype_, MAX_REF_FLANK_LEN, &base_quality_, num_reads_);
  int read_index = 0;
  for (unsigned int i = 0; i < alns_.size(); i++){
    // TO DO: Add check to see if sequence already encountered
    // If so, reuse alignment probs(even though qual scores may differ)
    // Should result in 5-7x speedup

    hap_aligner.process_reads(alns_[i], read_index, log_aln_probs_, seed_positions_); 
    read_index += alns_[i].size();
  }  

  std::cerr << "Genotyping samples..." << std::endl;
  if (stutter_model_ == NULL)
    printErrorAndDie("Must specify stutter model before running genotype()");
  calc_log_sample_posteriors();


  debug_sample(sample_indices_["LP6005441-DNA_E08"]);

  return true;
}

void SeqStutterGenotyper::write_vcf_header(std::vector<std::string>& sample_names, std::ostream& out){
  out << "##fileformat=VCFv4.1" << "\n";

  // Info field descriptors
  out << "##INFO=<ID=" << "INFRAME_PGEOM,"  << "Number=1,Type=Float,Description=\""   << "Parameter for in-frame geometric step size distribution"                      << "\">\n"
      << "##INFO=<ID=" << "INFRAME_UP,"     << "Number=1,Type=Float,Description=\""   << "Probability that stutter causes an in-frame increase in obs. STR size"        << "\">\n"
      << "##INFO=<ID=" << "INFRAME_DOWN,"   << "Number=1,Type=Float,Description=\""   << "Probability that stutter causes an in-frame decrease in obs. STR size"        << "\">\n"
      << "##INFO=<ID=" << "OUTFRAME_PGEOM," << "Number=1,Type=Float,Description=\""   << "Parameter for out-of-frame geometric step size distribution"                  << "\">\n"
      << "##INFO=<ID=" << "OUTFRAME_UP,"    << "Number=1,Type=Float,Description=\""   << "Probability that stutter causes an out-of-frame increase in obs. STR size"    << "\">\n"
      << "##INFO=<ID=" << "OUTFRAME_DOWN,"  << "Number=1,Type=Float,Description=\""   << "Probability that stutter causes an out-of-frame decrease in obs. STR size"    << "\">\n"
      << "##INFO=<ID=" << "BPDIFFS,"        << "Number=A,Type=Integer,Description=\"" << "Base pair difference of each alternate allele from the reference allele"      << "\">\n"
      << "##INFO=<ID=" << "START,"          << "Number=1,Type=Integer,Description=\"" << "Inclusive start coodinate for the repetitive potrion of the reference allele" << "\">\n"
      << "##INFO=<ID=" << "END,"            << "Number=1,Type=Integer,Description=\"" << "Inclusive end coordinate for the repetitive portion of the reference allele"  << "\">\n"
      << "##INFO=<ID=" << "AC,"             << "Number=A,Type=Integer,Description=\"" << "Alternate allele counts"                                                      << "\">\n"
      << "##INFO=<ID=" << "DELSKIP,"        << "Number=1,Type=Integer,Description=\"" << "Number of samples not genotyped due to problematic deletion boundaries"       << "\">\n";

  // Format field descriptors
  out << "##FORMAT=<ID=" << "GT"          << ",Number=1,Type=String,Description=\""  << "Genotype" << "\">" << "\n"
      << "##FORMAT=<ID=" << "GB"          << ",Number=1,Type=String,Description=\""  << "Base pair differences of genotype from reference"              << "\">" << "\n"
      << "##FORMAT=<ID=" << "Q"           << ",Number=1,Type=Float,Description=\""   << "Posterior probability of unphased genotype"                    << "\">" << "\n"
      << "##FORMAT=<ID=" << "PQ"          << ",Number=1,Type=Float,Description=\""   << "Posterior probability of phased genotype"                      << "\">" << "\n"
      << "##FORMAT=<ID=" << "DP"          << ",Number=1,Type=Integer,Description=\"" << "Read depth"                                                    << "\">" << "\n"
      << "##FORMAT=<ID=" << "DSNP"        << ",Number=1,Type=Integer,Description=\"" << "Number of reads with SNP phasing information"                  << "\">" << "\n"
      << "##FORMAT=<ID=" << "PDP"         << ",Number=1,Type=String,Description=\""  << "Fractional reads supporting each haploid genotype"             << "\">" << "\n"
      << "##FORMAT=<ID=" << "ALLREADS"    << ",Number=.,Type=Integer,Description=\"" << "Base pair difference observed in each read"                    << "\">" << "\n"
      << "##FORMAT=<ID=" << "PALLREADS"   << ",Number=.,Type=Float,Description=\""   << "Expected bp diff in each read based on haplotype alignment probabilities"<< "\">" << "\n";

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

  pos_ += 1; // First off-by-1 VCF error
}



/*
void SeqStutterGenotyper::set_allele_priors(vcf::VariantCallFile& variant_file){
  delete [] log_allele_priors_;
  log_allele_priors_  = new double[num_alleles_*num_alleles_*num_samples_];
  std::string PGP_KEY = "PGP";
  
  if (!variant_file.setRegion(chrom_, start_, end_)){
    // Retry setting region if chr is in chromosome name
    if (chrom_.size() <= 3 || chrom_.substr(0, 3).compare("chr") != 0 || !variant_file.setRegion(chrom_.substr(3), start_, end_))
      printErrorAndDie("Failed to set VCF region when obtaining allele priors");
  }
  vcf::Variant variant(variant_file);
  if (!variant_file.getNextVariant(variant))
    printErrorAndDie("Failed to extract VCF entry when obtaining allele priors");
  if (variant_file.formatTypes.find(GP_KEY) == variant_file.formatTypes.end())
    printErrorAndDie("VCF doesn't contain the PGP format field required for setting allele priors");

  int sample_count = 0;
  for (auto sample_iter = variant.sampleNames.begin(); sample_iter != variant.sampleNames.end(); ++sample_iter){
    if (sample_indices_.find(*sample_iter) == sample_indices_.end())
      continue;
    int sample_index = sample_indices_.find(*sample_iter)->second;
    sample_count++;

    int gp_index = 0;
    for (unsigned int i = 0; i < num_alleles_; ++i){
      for (unsigned int j = 0; j < num_alleles_; ++j, ++gp_index){
	  double prob = variant.getSampleValueFloat(PGP_KEY, *sample_iter, gp_index);
	  
	  
	}
      }
    }
    
    // TO DO: Parse BEAGLE format fields to set priors
    printErrorAndDie("set_allele_priors() function not fully implemented");
  }

  // Ensure that the VCF contained priors for all samples
  if (sample_count != num_samples_)
    printErrorAndDie("VCF only contained allele priors for a subset of samples");
}
*/
 
void SeqStutterGenotyper::debug_sample(int sample_index){
  std::cerr << "DEBUGGING SAMPLE..." << std::endl;
  std::cerr << "READ LL's:" << std::endl;
  double* read_LL_ptr = log_aln_probs_;
  int read_index = 0;
  for (unsigned int i = 0; i < num_reads_; ++i){
    if(sample_label_[i] == sample_index){
      std::cerr << "\t" << "READ #" << read_index << ", SEED BASE=" << seed_positions_[i] << " " 
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
      std::cerr << index_1 << " " << index_2 << " " << *sample_LL_ptr << std::endl;
      sample_LL_ptr += num_samples_;
    }
  
  std::cerr << "END OF SAMPLE DEBUGGING..." << std::endl;
}


double SeqStutterGenotyper::calc_log_sample_posteriors(){
  std::vector<double> sample_max_LLs(num_samples_, -DBL_MAX);
  double* sample_LL_ptr = log_sample_posteriors_;

  if (log_allele_priors_ != NULL)
    memcpy(log_sample_posteriors_, log_allele_priors_, num_alleles_*num_alleles_*num_samples_*sizeof(double));
  else
    std::fill(log_sample_posteriors_, log_sample_posteriors_+(num_alleles_*num_alleles_*num_samples_), -2*log(num_alleles_));

  for (int index_1 = 0; index_1 < num_alleles_; ++index_1){
    for (int index_2 = 0; index_2 < num_alleles_; ++index_2){
      double* read_LL_ptr = log_aln_probs_;
      for (int read_index = 0; read_index < num_reads_; ++read_index){
	sample_LL_ptr[sample_label_[read_index]] += log_sum_exp(LOG_ONE_HALF + log_p1_[read_index] + read_LL_ptr[index_1], 
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
  std::vector<double> sample_total_LLs(num_samples_, 0.0);
  sample_LL_ptr = log_sample_posteriors_;
  for (int index_1 = 0; index_1 < num_alleles_; ++index_1)
    for (int index_2 = 0; index_2 < num_alleles_; ++index_2)
      for (int sample_index = 0; sample_index < num_samples_; ++sample_index, ++sample_LL_ptr)
	sample_total_LLs[sample_index] += exp(*sample_LL_ptr - sample_max_LLs[sample_index]);
  for (int sample_index = 0; sample_index < num_samples_; ++sample_index){
    sample_total_LLs[sample_index] = sample_max_LLs[sample_index] + log(sample_total_LLs[sample_index]);
    assert(sample_total_LLs[sample_index] <= TOLERANCE);
  }
  // Compute the total log-likelihood given the current parameters
  double total_LL = sum(sample_total_LLs);

  // Normalize each genotype LL to generate valid log posteriors
  sample_LL_ptr = log_sample_posteriors_;
  for (int index_1 = 0; index_1 < num_alleles_; ++index_1)
    for(int index_2 = 0; index_2 < num_alleles_; ++index_2)
      for (int sample_index = 0; sample_index < num_samples_; ++sample_index, ++sample_LL_ptr)
	*sample_LL_ptr -= sample_total_LLs[sample_index];  

  return total_LL;
}


void SeqStutterGenotyper::write_vcf_record(std::vector<std::string>& sample_names, std::ostream& out){
  std::vector< std::pair<int,int> > gts(num_samples_, std::pair<int,int>(-1,-1));
  std::vector<double> log_phased_posteriors(num_samples_, -DBL_MAX);
  std::vector< std::vector<double> > log_read_phases(num_samples_);
  std::vector<double> log_unphased_posteriors, phase_probs;

  // TO DO: Consider selecting GT based on genotype with maximum UNPHASED posterior instead of maximum PHASED posterior
  // Are we then double-counting het GTs vs hom GTs?
  // TO DO: How do we pool STRs with identical lengths?
  
  // Extract each sample's MAP genotype and the associated posterior
  double* log_post_ptr = log_sample_posteriors_;
  for (int index_1 = 0; index_1 < num_alleles_; ++index_1)
    for (int index_2 = 0; index_2 < num_alleles_; ++index_2)
      for (unsigned int sample_index = 0; sample_index < num_samples_; ++sample_index, ++log_post_ptr){
	if (*log_post_ptr > log_phased_posteriors[sample_index]){
	  log_phased_posteriors[sample_index] = *log_post_ptr;
	  gts[sample_index] = std::pair<int,int>(index_1, index_2);
	}
      }
  
  // Compute allele counts for samples of interest
  std::set<std::string> samples_of_interest(sample_names.begin(), sample_names.end());
  std::vector<int> allele_counts(num_alleles_);
  int sample_index   = 0;
  int del_skip_count = 0;
  for (auto gt_iter = gts.begin(); gt_iter != gts.end(); ++gt_iter, ++sample_index){
    if (samples_of_interest.find(sample_names_[sample_index]) == samples_of_interest.end())
      continue;
    if (call_sample_[sample_index]) {
      allele_counts[gt_iter->first]++;
      allele_counts[gt_iter->second]++;
    }
    else
      del_skip_count++;
  }

  std::cerr << "Allele counts" << std::endl;
  for (unsigned int i = 0; i < alleles_.size(); i++)
    std::cerr << alleles_[i] << " " << allele_counts[i] << std::endl;
  std::cerr << std::endl;

  // Extract the phasing probability conditioned on the determined sample genotypes
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

  // Determine the number of reads with SNP information for each sample
  std::vector<int> num_reads_with_snps(num_samples_, 0);
  for (unsigned int read_index = 0; read_index < num_reads_; read_index++)
    num_reads_with_snps[sample_label_[read_index]] += (log_p1_[read_index] != log_p2_[read_index]);
  
  // Extract each read's phase posterior conditioned on the determined sample genotypes
  double* read_LL_ptr = log_aln_probs_;
  for (unsigned int read_index = 0; read_index < num_reads_; read_index++){
    int gt_a = gts[sample_label_[read_index]].first;
    int gt_b = gts[sample_label_[read_index]].second;
    double total_read_LL = log_sum_exp(LOG_ONE_HALF+log_p1_[read_index]+read_LL_ptr[gt_a], LOG_ONE_HALF+log_p2_[read_index]+read_LL_ptr[gt_b]);
    double log_phase_one = LOG_ONE_HALF + log_p1_[read_index] + read_LL_ptr[gt_a] - total_read_LL; 
    log_read_phases[sample_label_[read_index]].push_back(log_phase_one);
    read_LL_ptr += num_alleles_;
  }

  // Determine the bp differences observed in reads for each sample
  assert(bp_diffs_.size() == num_reads_);
  std::vector< std::vector<int> > bps_per_sample(num_samples_);
  for (unsigned int read_index = 0; read_index < num_reads_; read_index++){
    bps_per_sample[sample_label_[read_index]].push_back(bp_diffs_[read_index]);
  }

  assert(haplotype_->num_blocks() == 3);

  //VCF line format = CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE_1 SAMPLE_2 ... SAMPLE_N
  out << region_->chrom() << "\t" << pos_ << "\t";
  
  if (region_->name().empty())
    out << ".";
  else
    out << region_->name();

  // Compute the base pair differences from the reference
  std::vector<int> bp_diffs;
  for (unsigned int i = 0; i < alleles_.size(); i++)
    bp_diffs.push_back((int)alleles_[i].size() - (int)alleles_[0].size());

  // Determine the posterior bp differences observed in reads for each sample
  read_LL_ptr = log_aln_probs_;
  std::vector< std::vector<double> > posterior_bps_per_sample(num_samples_);
  for (unsigned int read_index = 0;  read_index < num_reads_; read_index++){
    posterior_bps_per_sample[sample_label_[read_index]].push_back(expected_value(read_LL_ptr, bp_diffs));
    read_LL_ptr += num_alleles_;
  }

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
      << "START="           << region_->start() << ";"
      << "END="             << region_->stop()  << ";"
      << "DELSKIP="         << del_skip_count   << ";";
  if (num_alleles_ > 1){
    out << "BPDIFFS=" << bp_diffs[1];
    for (unsigned int i = 2; i < num_alleles_; i++)
      out << "," << bp_diffs[i];
    out << ";";
  }
  // Add allele counts
  if (allele_counts.size() > 1){
    out << "AC=";
    for (unsigned int i = 1; i < allele_counts.size()-1; i++)
      out << allele_counts[i] << ",";
    out << allele_counts.back() << ";";
  }

  // Add FORMAT field
  out << "\tGT:GB:Q:PQ:DP:DSNP:PDP:ALLREADS:PALLREADS";

  for (unsigned int i = 0; i < sample_names.size(); i++){
    out << "\t";
    auto sample_iter = sample_indices_.find(sample_names[i]);
    if (sample_iter == sample_indices_.end() || reads_per_sample_[sample_iter->second] == 0){
      out << ".";
      continue;
    }
    
    // Don't report information for a sample if flag has been set to false
    // Likely due to problematic deletion boundaries
    if (!call_sample_[sample_iter->second]){
      out << ".";
      continue;
    }

    int sample_index    = sample_iter->second;
    int total_reads     = reads_per_sample_[sample_index];
    double phase1_reads = exp(log_sum_exp(log_read_phases[sample_index]));
    double phase2_reads = total_reads - phase1_reads;

    out << gts[sample_index].first << "|" << gts[sample_index].second                            // Genotype
	<< ":" << bp_diffs[gts[sample_index].first] << "|" << bp_diffs[gts[sample_index].second] // Base pair differences from reference
	<< ":" << exp(log_unphased_posteriors[sample_index])                                     // Unphased posterior
	<< ":" << exp(log_phased_posteriors[sample_index])                                       // Phased posterior
	<< ":" << total_reads                                                                    // Total reads
	<< ":" << num_reads_with_snps[sample_index]                                              // Total reads with SNP information
	<< ":" << phase1_reads << "|" << phase2_reads;                                           // Reads per allele

    // Add bp diffs back in for debugging purposes
    if (bps_per_sample[sample_index].size() != 0){
      std::sort(bps_per_sample[sample_index].begin(), bps_per_sample[sample_index].end());
      out << ":" << bps_per_sample[sample_index][0];
      for (unsigned int j = 1; j < bps_per_sample[sample_index].size(); j++)
	out << "," << bps_per_sample[sample_index][j];
    }
    else
      out << ":" << ".";

    // Expected base pair differences from alignment probabilities
    if (posterior_bps_per_sample[sample_index].size() != 0){
      std::sort(posterior_bps_per_sample[sample_index].begin(), posterior_bps_per_sample[sample_index].end());
      out << ":" << posterior_bps_per_sample[sample_index][0];
      for (unsigned int j = 1; j < posterior_bps_per_sample[sample_index].size(); j++)
        out << "," << posterior_bps_per_sample[sample_index][j];
    }
    else
      out << ":" << ".";
  }
  out << "\n";
}

