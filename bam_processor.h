#ifndef BAM_PROCESSOR_H_
#define BAM_PROCESSOR_H_

#include <iostream>
#include <string>
#include <vector>

#include "bamtools/include/api/BamAlignment.h"
#include "bamtools/include/api/BamMultiReader.h"
#include "bamtools/include/api/BamWriter.h"

#include "region.h"

class BamProcessor {
 private:
  void read_and_filter_reads(BamTools::BamMultiReader& reader, std::string& chrom_seq, 
			     std::vector<Region>::iterator region_iter, std::map<std::string, std::string>& file_read_groups, 
			     std::vector<std::string>& rg_names, std::vector< std::vector<BamTools::BamAlignment> >& alignments_by_rg, 
			     BamTools::BamWriter& bam_writer);
  public:
  void process_regions(BamTools::BamMultiReader& reader, 
		       std::string& region_file, std::string& fasta_dir,
		       std::map<std::string, std::string>& file_read_groups,
		       BamTools::BamWriter& bam_writer, std::ostream& out);
  
  void process_reads(std::vector< std::vector<BamTools::BamAlignment> >& alignments_by_rg, std::vector<std::string>& rg_names, Region& region, std::ostream& out);
   
  int MAX_MATE_DIST            = 1000;
  int MIN_BP_BEFORE_INDEL      = 7;
  int MIN_FLANK                = 5;
  int MIN_READ_END_MATCH       = 10;
  int MAXIMAL_END_MATCH_WINDOW = 15;
  int REMOVE_MULTIMAPPERS      = 0;
};


#endif
