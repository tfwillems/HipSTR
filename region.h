#ifndef REGION_H_
#define REGION_H_

#include <map>
#include <set>
#include <string>
#include <vector>
#include <stdlib.h>

#include "error.h"

class Region{
private:
  std::string chrom_;
  uint32_t start_, stop_;
  std::set<std::string> sample_set_;
  std::map<std::string, std::string> sample_annots_;
  int period_;
public:
  Region(std::string chrom, uint32_t start, uint32_t stop, int period){
    chrom_ = chrom; start_ = start; stop_ = stop; period_ = period;
  }
  const std::string& chrom() const { return chrom_;}
  uint32_t start() const { return start_;}
  uint32_t  stop() const { return stop_; }
  int     period() const { return period_; }     
};

void readRegions(std::string& input_file, std::vector<Region>& regions);

void orderRegions(std::vector<Region>& input_regions, 
		  std::vector< std::vector<Region> >& output_regions, 
		  std::map<std::string, int>& chrom_order);

#endif
