#ifndef SEQ_IO_H_
#define SEQ_IO_H_

#include <string>

bool is_file(const std::string& name);

void readFastaFromDir(std::string input_file, std::string fasta_dir,
		      std::string& res);

#endif

