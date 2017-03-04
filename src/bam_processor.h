#ifndef BAM_PROCESSOR_H_
#define BAM_PROCESSOR_H_

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "bam_io.h"
#include "base_quality.h"
#include "error.h"
#include "region.h"
#include "stringops.h"

class BamProcessor {
 protected:
  typedef std::vector<BamAlignment> BamAlnList;

 private:
  bool use_bam_rgs_;
  bool rem_pcr_dups_;

  // Timing statistics (in seconds)
  double total_bam_seek_time_;
  double locus_bam_seek_time_;
  double total_read_filter_time_;
  double locus_read_filter_time_;


  void  write_passing_alignment(BamAlignment& aln, std::map<std::string, std::string>& rg_to_sample, BamWriter* writer);
  void write_filtered_alignment(BamAlignment& aln, std::string filter, std::map<std::string, std::string>& rg_to_sample, BamWriter* writer);

  void extract_mappings(BamAlignment& aln, const BamHeader* bam_header,
			std::vector< std::pair<std::string, int32_t> >& chrom_pos_pairs);

  void get_valid_pairings(BamAlignment& aln_1, BamAlignment& aln_2, const BamHeader* bam_header,
			  std::vector< std::pair<std::string, int32_t> >& p1, std::vector< std::pair<std::string, int32_t> >& p2);

  void read_and_filter_reads(BamCramMultiReader& reader, std::string& chrom_seq, RegionGroup& region,
			     std::map<std::string, std::string>& rg_to_sample, std::map<std::string, std::string>& rg_to_library, std::vector<std::string>& rg_names,
			     std::vector<BamAlnList>& paired_strs_by_rg, std::vector<BamAlnList>& mate_pairs_by_rg, std::vector<BamAlnList>& unpaired_strs_by_rg,
			     BamWriter* pass_writer, BamWriter* filt_writer);

 std::string get_read_group(BamAlignment& aln, std::map<std::string, std::string>& read_group_mapping);

 std::string trim_alignment_name(BamAlignment& aln);

 bool spans_a_region(const std::vector<Region>& regions, BamAlignment& alignment);

 protected:
 BaseQuality base_quality_;

 bool bams_from_10x_; // True iff BAMs were generated from 10X GEMCODE platform

 bool log_to_file_;
 std::ofstream log_;

 std::set<std::string> sample_set_;

  public:
 BamProcessor(bool use_bam_rgs, bool remove_pcr_dups){
   use_bam_rgs_             = use_bam_rgs;
   rem_pcr_dups_            = remove_pcr_dups;
   MAX_MATE_DIST            = 1000;
   MIN_BP_BEFORE_INDEL      = 7;
   MIN_FLANK                = 5;
   MIN_READ_END_MATCH       = 10;
   MAXIMAL_END_MATCH_WINDOW = 15;
   REQUIRE_SPANNING         = true;
   REQUIRE_PAIRED_READS     = 1;
   total_bam_seek_time_     = 0;
   locus_bam_seek_time_     = -1;
   total_read_filter_time_  = 0;
   locus_read_filter_time_  = -1;
   MAX_STR_LENGTH           = 100;
   MIN_SUM_QUAL_LOG_PROB    = -10;
   log_to_file_             = false;
   MAX_TOTAL_READS          = 1000000;
   BASE_QUAL_TRIM           = '5';
   bams_from_10x_           = false;
 }

 ~BamProcessor(){
   if (log_to_file_)
     log_.close();
 }

 double total_bam_seek_time()    { return total_bam_seek_time_;    }
 double locus_bam_seek_time()    { return locus_bam_seek_time_;    }
 double total_read_filter_time() { return total_read_filter_time_; }
 double locus_read_filter_time() { return locus_read_filter_time_; }
 void use_custom_read_groups()   { use_bam_rgs_ = false;           }
 void allow_pcr_dups()           { rem_pcr_dups_ = false;          }

 void process_regions(BamCramMultiReader& reader,
		      std::string& region_file, std::string& fasta_file,
		      std::map<std::string, std::string>& rg_to_sample, std::map<std::string, std::string>& rg_to_library,
		      BamWriter* pass_writer, BamWriter* filt_writer,
		      std::ostream& out, int32_t max_regions, std::string chrom);
  
 virtual void process_reads(std::vector<BamAlnList>& paired_strs_by_rg,
			    std::vector<BamAlnList>& mate_pairs_by_rg,
			    std::vector<BamAlnList>& unpaired_strs_by_rg,
			    std::vector<std::string>& rg_names, RegionGroup& region_group, std::string& chrom_seq,
			    std::ostream& out){
   log("Doing nothing with reads");
 }


 void set_log(std::string log_file){
   if (log_to_file_)
     printErrorAndDie("Cannot reset the log file multiple times");
   log_to_file_ = true;
   log_.open(log_file, std::ofstream::out);
   if (!log_.is_open())
     printErrorAndDie("Failed to open the log file: " + log_file);
 }

 inline void log(std::string msg){
   if (log_to_file_)
     log_ << msg << std::endl;
   else
     std::cerr << msg << std::endl;
 }

 inline std::ostream& logger(){
   return (log_to_file_ ? log_ : std::cerr);
 }

 void set_sample_set(std::string sample_names){
   std::vector<std::string> sample_list;
   split_by_delim(sample_names, ',', sample_list);
   sample_set_ = std::set<std::string>(sample_list.begin(), sample_list.end());
 }

 void use_10x_bam_tags(){
   bams_from_10x_ = true;
 }

 static void add_passes_filters_tag(BamAlignment& aln, std::string& passes);

 static void passes_filters(BamAlignment& aln, std::vector<bool>& region_passes);

 int32_t MAX_MATE_DIST;
 int32_t MIN_BP_BEFORE_INDEL;
 int32_t MIN_FLANK;
 int32_t MIN_READ_END_MATCH;
 int32_t MAXIMAL_END_MATCH_WINDOW;
 int32_t MAX_STR_LENGTH;
 bool    REQUIRE_SPANNING;
 int     REQUIRE_PAIRED_READS;  // Only utilize paired STR reads to genotype individuals
 double  MIN_SUM_QUAL_LOG_PROB;
 int32_t MAX_TOTAL_READS;       // Skip loci where the number of STR reads passing all filters exceeds this limit
 char    BASE_QUAL_TRIM;        // Trim boths ends of the read until encountering a base with quality greater than this threshold
 bool    TOO_MANY_READS;        // Flag set if the current locus being processed as too many reads
};

#endif
