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

void readRegions(std::string& input_file, std::vector<Region>& regions, uint32_t max_regions, std::string chrom_limit, std::ostream& logger){
  logger << "Reading region file " << input_file << std::endl;
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
      printErrorAndDie("Improperly formatted region file. \nRequired format is tab-delimited columns CHROM START STOP PERIOD NCOPIES\n Bad line: " + line);
    if (start < 1)      printErrorAndDie("Improperly formatted region file. \n Region has a START < 1, but START must be >= 1\n Bad line: " + line);
    if (stop <= start)  printErrorAndDie("Improperly formatted region file. \n Region has a STOP <= START. Bad line: " + line);
    if (period < 1)     printErrorAndDie("Improperly formatted region file. \n Region has a PERIOD < 1. Bad line: " + line);
    if (period > 9)     printErrorAndDie("Improperly formatted region file. \n Region has a PERIOD > 9. Bad line: " + line);

    if (!chrom_limit.empty() && chrom.compare(chrom_limit) != 0)
      continue;
    if (iss >> name)
      regions.push_back(Region(chrom, start-1, stop, period, name));
    else
      regions.push_back(Region(chrom, start-1, stop, period));
  }
  input.close();
  logger << "Region file contains " << regions.size() << " regions" << std::endl;
}

void orderRegions(std::vector<Region>& regions){
  std::sort(regions.begin(), regions.end());
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
    std::sort(output_regions[i].begin(), output_regions[i].end());
} 



