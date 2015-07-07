#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>

#include "error.h"
#include "region.h"

void readRegions(std::string& input_file, std::vector<Region>& regions, uint32_t max_regions, std::string chrom_limit){
  std::cerr << "Reading region file " << input_file << std::endl;
  std::ifstream input(input_file.c_str());
  if (!input.is_open()) 
    printErrorAndDie("Failed to open region file");

  regions.clear();
  std::string line;
  while (std::getline(input, line) && regions.size() < max_regions){
    std::istringstream iss(line);
    std::string chrom, name;
    int32_t start, stop;
    int period;
    double ref_copy;
    if (!(iss >> chrom >> start >> stop >> period >> ref_copy))
      printErrorAndDie("Improperly formatted region file. Required format is tab-delimited columns CHROM START STOP PERIOD NCOPIES");
    if (!chrom_limit.empty() && chrom.compare(chrom_limit) != 0)
      continue;
    if (iss >> name)
      regions.push_back(Region(chrom, start, stop, period, name));
    else
      regions.push_back(Region(chrom, start, stop, period));
  }
  input.close();
  std::cerr << "Region file contains " << regions.size() << " regions" << std::endl;
}

bool region_lt(const Region& r1, const Region& r2){
  int chrom_comp = (r1.chrom().compare(r2.chrom()));
  if (chrom_comp != 0)
    return chrom_comp < 0;
  else
    return r1.start() < r2.start();
}

void orderRegions(std::vector<Region>& regions){
  std::sort(regions.begin(), regions.end(), region_lt);
}

void orderRegions(std::vector<Region>& input_regions, std::vector< std::vector<Region> >& output_regions, std::map<std::string, int>& chrom_order){
  output_regions.clear();
  chrom_order.clear();
  int chrom_count = 0;
  for (auto iter = input_regions.begin(); iter != input_regions.end(); iter++){
    if (chrom_order.find(iter->chrom()) == chrom_order.end()){
      chrom_order[iter->chrom()] = chrom_count++;
      output_regions.push_back(std::vector<Region>());
    }
    output_regions[chrom_order[iter->chrom()]].push_back(*iter);
  }
  for (unsigned int i = 0; i < output_regions.size(); i++)
    std::sort(output_regions[i].begin(), output_regions[i].end(), region_lt);
} 



