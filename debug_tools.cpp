#include "debug_tools.h"

#include <iostream>

void print_vector(std::vector<double>& vals){
  for (unsigned int i = 0; i < vals.size(); i++)
    std::cerr << vals[i] << " ";
  std::cerr << std::endl;
}

