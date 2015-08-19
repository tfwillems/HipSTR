#ifndef BAM_PROCESSOR_H_
#define BAM_PROCESSOR_H_

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "bamtools/include/api/BamAlignment.h"
#include "bamtools/include/api/BamMultiReader.h"
#include "bamtools/include/api/BamWriter.h"

#include "base_quality.h"
#include "error.h"
#include "region.h"

class BamProcessor {
 private:
  bool use_bam_rgs_;
  bool check_mate_info_;
  bool rem_pcr_dups_;

  // Timing statistics (in seconds)
  double total_bam_seek_time_;
  double locus_bam_seek_time_;
  double total_read_filter_time_;
  double locus_read_filter_time_;

 void read_and_filter_reads(BamTools::BamMultiReader& reader, std::string& chrom_seq,
			    std::vector<Region>::iterator region_iter,
			    std::map<std::string, std::string>& rg_to_sample, std::map<std::string, std::string>& rg_to_library,
			    std::vector<std::string>& rg_names,
			    std::vector< std::vector<BamTools::BamAlignment> >& paired_strs_by_rg,                                                             
			    std::vector< std::vector<BamTools::BamAlignment> >& mate_pairs_by_rg,                                                              
			    std::vector< std::vector<BamTools::BamAlignment> >& unpaired_strs_by_rg,                                                           
			    BamTools::BamWriter& bam_writer);

 std::string get_read_group(BamTools::BamAlignment& aln, std::map<std::string, std::string>& read_group_mapping);

 std::string trim_alignment_name(BamTools::BamAlignment& aln);

 protected:
 BaseQuality base_quality_;

 bool log_to_file_;
 std::ofstream log_;

  public:
 BamProcessor(bool use_bam_rgs, bool filter_by_mate, bool remove_pcr_dups){
   use_bam_rgs_             = use_bam_rgs;
   check_mate_info_         = filter_by_mate;
   rem_pcr_dups_            = remove_pcr_dups;
   MAX_MATE_DIST            = 1000;
   MIN_BP_BEFORE_INDEL      = 7;
   MIN_FLANK                = 5;
   MIN_READ_END_MATCH       = 10;
   MAXIMAL_END_MATCH_WINDOW = 15;
   REMOVE_MULTIMAPPERS      = 0;
   REQUIRE_SPANNING         = true;
   MIN_MAPPING_QUALITY      = 0;
   REMOVE_READS_WITH_N      = 1;
   total_bam_seek_time_     = 0;
   locus_bam_seek_time_     = -1;
   total_read_filter_time_  = 0;
   locus_read_filter_time_  = -1;
   MAX_SOFT_CLIPS           = 100000;
   MAX_HARD_CLIPS           = 100000;
   log_to_file_             = false;
 }

 ~BamProcessor(){
   if (log_to_file_)
     log_.close();
 }

 void remove_all_filters(){
   MIN_BP_BEFORE_INDEL      = 0;
   MIN_FLANK                = 0;
   MIN_READ_END_MATCH       = 0;
   MAXIMAL_END_MATCH_WINDOW = 0;
   REMOVE_MULTIMAPPERS      = 0;
   REQUIRE_SPANNING         = false;
   MIN_MAPPING_QUALITY      = 0;
   REMOVE_READS_WITH_N      = 0;
   MAX_SOFT_CLIPS           = 100000;
   MAX_HARD_CLIPS           = 100000;
 }

 double total_bam_seek_time()    { return total_bam_seek_time_;    }
 double locus_bam_seek_time()    { return locus_bam_seek_time_;    }
 double total_read_filter_time() { return total_read_filter_time_; }
 double locus_read_filter_time() { return locus_read_filter_time_; }
 void use_custom_read_groups()   { use_bam_rgs_ = false;           }
 void allow_pcr_dups()           { rem_pcr_dups_ = false;          }

 void set_min_mapping_quality(int quality) { MIN_MAPPING_QUALITY = quality; }

 void process_regions(BamTools::BamMultiReader& reader,
		      std::string& region_file, std::string& fasta_dir,
		      std::map<std::string, std::string>& rg_to_sample, std::map<std::string, std::string>& rg_to_library,
		      BamTools::BamWriter& bam_writer, std::ostream& out, int32_t max_regions, std::string chrom);
  
 virtual void process_reads(std::vector< std::vector<BamTools::BamAlignment> >& paired_strs_by_rg,
			    std::vector< std::vector<BamTools::BamAlignment> >& mate_pairs_by_rg,
			    std::vector< std::vector<BamTools::BamAlignment> >& unpaired_strs_by_rg,
			    std::vector<std::string>& rg_names, Region& region, std::string& ref_allele, std::string& chrom_seq,
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

 int32_t MAX_MATE_DIST;
 int32_t MIN_BP_BEFORE_INDEL;
 int32_t MIN_FLANK;
 int32_t MIN_READ_END_MATCH;
 int32_t MAXIMAL_END_MATCH_WINDOW;
 int32_t MIN_MAPPING_QUALITY;
 int32_t MAX_SOFT_CLIPS;
 int32_t MAX_HARD_CLIPS;
 int REMOVE_MULTIMAPPERS;
 int REMOVE_READS_WITH_N;
 bool REQUIRE_SPANNING;
};


#endif
