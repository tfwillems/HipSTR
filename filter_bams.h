#ifndef FILTER_BAMS_H_
#define FILTER_BAMS_H_

#include <map>
#include <string>
#include <vector>

#include "bamtools/include/api/BamAlignment.h"
#include "bamtools/include/api/BamReader.h"

#include "insert_size.h"
#include "region.h"


void filter_bam_paired_mode(BamTools::BamReader& reader,
                            std::vector< std::vector<Region> >& regions,
                            std::map<std::string, int>& chrom_order,
                            std::string& output_filename,
                            InsertSizeCounter& counter,
			    bool analyze_insert_size);

void filter_bam(BamTools::BamReader& reader,
                std::vector< std::vector<Region> >&regions,
                std::map<std::string, int>& chrom_order,
                std::string& output_filename,
		InsertSizeCounter& counter,
		bool analyze_insert_size);
#endif
