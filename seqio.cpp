#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>

#include "error.h"
#include "seqio.h"

bool is_file(const std::string& name){
  struct stat st_buf;
  int status = stat(name.c_str(), &st_buf);
  return (S_ISREG (st_buf.st_mode));
}

void readFastaFromDir(std::string input_file, std::string fasta_dir, std::string& res){
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

