#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "error.h"
#include "seqio.h"

void readFasta(std::string input_file, std::string fasta_dir, std::string& res){
  std::stringstream res_stream;
  std::string path = fasta_dir + input_file;
  std::ifstream data(path.c_str());
  if (!data.is_open())
    printErrorAndDie("Unable to open FASTA file " + path);
  std::string line;
  getline(data, line);
  if (line.size() == 0) printErrorAndDie("FASTA file " + input_file + " is empty");
  if (line[0] != '>') printErrorAndDie("First line of FASTA file " + input_file + " must begin with a '>'.");
  
  while(getline(data, line))
    res_stream << line;
  res = res_stream.str();
  data.close();
}

