#ifndef BAM_PROCESSOR_H_
#define BAM_PROCESSOR_H_

#include <iostream>
#include <string>
#include <vector>

#include "bamtools/include/api/BamAlignment.h"
#include "bamtools/include/api/BamMultiReader.h"
#include "bamtools/include/api/BamWriter.h"

#include "error.h"
#include "region.h"

class BamProcessor {
 private:
  bool use_lobstr_rg_;
  bool check_mate_info_;

 void read_and_filter_reads(BamTools::BamMultiReader& reader, std::string& chrom_seq,
			    std::vector<Region>::iterator region_iter, std::map<std::string, std::string>& file_read_groups,
			    std::vector<std::string>& rg_names,
			    std::vector< std::vector<BamTools::BamAlignment> >& paired_strs_by_rg,                                                             
			    std::vector< std::vector<BamTools::BamAlignment> >& mate_pairs_by_rg,                                                              
			    std::vector< std::vector<BamTools::BamAlignment> >& unpaired_strs_by_rg,                                                           
			    BamTools::BamWriter& bam_writer);

 std::string parse_lobstr_rg(BamTools::BamAlignment& aln);

 std::string trim_alignment_name(BamTools::BamAlignment& aln);

  public:
 BamProcessor(bool use_lobstr_rg, bool filter_by_mate){
   use_lobstr_rg_    = use_lobstr_rg;
   check_mate_info_  = filter_by_mate;
   MAX_MATE_DIST            = 1000;
   MIN_BP_BEFORE_INDEL      = 7;
   MIN_FLANK                = 5;
   MIN_READ_END_MATCH       = 10;
   MAXIMAL_END_MATCH_WINDOW = 15;
   REMOVE_MULTIMAPPERS      = 0;
 }

  void process_regions(BamTools::BamMultiReader& reader, 
		       std::string& region_file, std::string& fasta_dir,
		       std::map<std::string, std::string>& file_read_groups,
		       BamTools::BamWriter& bam_writer, std::ostream& out, int32_t max_regions);

  
  virtual void process_reads(std::vector< std::vector<BamTools::BamAlignment> >& paired_strs_by_rg,
			     std::vector< std::vector<BamTools::BamAlignment> >& mate_pairs_by_rg,
			     std::vector< std::vector<BamTools::BamAlignment> >& unpaired_strs_by_rg,
			     std::vector<std::string>& rg_names, Region& region, std::string& ref_allele, std::string& chrom_seq,
			     std::ostream& out){
    std::cerr << "Doing nothing with reads" << std::endl;
  }

  void set_lobstr_rg_usage(bool use_lobstr_rg){
    use_lobstr_rg_ = use_lobstr_rg;
  }


  int MAX_MATE_DIST;
  int MIN_BP_BEFORE_INDEL;
  int MIN_FLANK;
  int MIN_READ_END_MATCH;
  int MAXIMAL_END_MATCH_WINDOW;
  int REMOVE_MULTIMAPPERS;
};


#endif
