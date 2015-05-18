#include <iostream>
#include <string>
#include <vector>

#include "../base_quality.h"

int main(){
  std::vector<const std::string*> base_qualities;
  base_qualities.push_back(new std::string("A!3=E!EF"));
  base_qualities.push_back(new std::string("AEAAE!ED"));
  base_qualities.push_back(new std::string("A!AAE!ED"));
  BaseQuality base_qual;  
  std::string res = base_qual.average_base_qualities(base_qualities);
  std::cerr << *base_qualities[0] << std::endl
	    << *base_qualities[1] << std::endl
	    << *base_qualities[2] << std::endl
	    << res << std::endl;
  for (unsigned int i = 0; i < base_qualities.size(); i++)
    delete base_qualities[i];
}
