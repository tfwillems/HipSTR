#ifndef FILTER_BAMS_H_
#define FILTER_BAMS_H_

#include <map>
#include <string>
#include <vector>

#include "bam_io.h"
#include "region.h"

void filter_bam_paired_mode(BamCramReader& reader,
                            std::vector< std::vector<Region> >& regions,
                            std::map<std::string, int>& chrom_order,
                            std::string& output_filename);

void filter_bam(BamCramReader& reader,
                std::vector< std::vector<Region> >&regions,
                std::map<std::string, int>& chrom_order,
                std::string& output_filename);
#endif
