#include "fasta_reader.h"

#include <assert.h>
#include <dirent.h>
#include <sstream>

#include "error.h"
#include "stringops.h"

void FastaReader::add_index(std::string& path){
  // Check if path for file exists
  if (!file_exists(path))
    printErrorAndDie("FASTA file " + path + " does not exist");

  // Check if path for index exists
  if (!file_exists(path + ".fai")){
    std::stringstream error_msg;
    error_msg << "No FASTA index file exists for " <<  path << "\n"
	      << "Please rerun the analysis after generating the index using the command:\n" 
	      << "\tsamtools faidx " << path << "\n";
    printErrorAndDie(error_msg.str());
  }

  // Use htslib to load the index
  faidx_t* index = fai_load(path.c_str());
  if (index == NULL)
    printErrorAndDie("Failed to load FASTA index file for " + path);

  // Add all of the chromosomes in the index to the map
  int num_seqs = faidx_nseq(index);
  for (int i = 0; i < num_seqs; i++){
    std::string seq_name = faidx_iseq(index, i);
    if (chrom_to_index_.find(seq_name) != chrom_to_index_.end())
      printErrorAndDie("Multiple entries for chromosome " + seq_name + " exist in FASTA files");
    chrom_to_index_[seq_name] = index;
  }
  fasta_indices_.push_back(index);
}

void FastaReader::init(std::string& path){
  assert(chrom_to_index_.empty() && fasta_indices_.empty());

  if (is_file(path))
    add_index(path);
  else {
    DIR* dir = opendir(path.c_str());
    if (dir == NULL)
      printErrorAndDie("Failed to access directory " + path);
    
    struct dirent* file_info;
    while ((file_info = readdir(dir)) != NULL){
      std::string filename(file_info->d_name);
      if (string_ends_with(filename, ".fa")){
	std::string full_path = path + "/" + filename;
	add_index(full_path);
      }
    }
    closedir(dir);
    
    if (fasta_indices_.empty())
      printErrorAndDie("Failed to locate any FASTA files in the provided directory: \n\t" + path);
  }
}

