#ifndef FASTA_READER_H_
#define FASTA_READER_H_

#include <assert.h>
#include <map>
#include <stdlib.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#include "error.h"
#include "stringops.h"

extern "C" {
#include "htslib/htslib/faidx.h"
}

/*
 * Simple class that facilitates extracting FASTA sequences from one or more indexed FASTA files.
 * Internally, it uses htslib functions to provide this functionality
 * The constructor accepts either a path to a directory containing indexed FASTA files or a path
 * to a single indexed FASTA file that contains one or more chromosomes.
 */
class FastaReader {
 private:
  std::map<std::string, faidx_t*> chrom_to_index_;
  std::vector<faidx_t*> fasta_indices_;

  bool file_exists(std::string path) const {
    return (access(path.c_str(), F_OK) != -1);
  }

  bool is_file(const std::string& name) const {
    struct stat st_buf;
    stat(name.c_str(), &st_buf);
    return (S_ISREG (st_buf.st_mode));
  }

  void add_index(const std::string& path);
  void init(const std::string& path);

  // Private unimplemented copy constructor and assignment operator to prevent operations
  FastaReader(const FastaReader& other);
  FastaReader& operator=(const FastaReader& other);

 public:
  /*
   * PATH is either (i)  a single indexed FASTA file, containing one or more chromosomes
   *             or (ii) a directory containing one or more indexed FASTA files
   * Sequences from all chromosomes in the relevant FASTA file(s) will be available for queries
   */
  explicit FastaReader(const std::string& path){
    init(path);
  }

  ~FastaReader(){
    for (unsigned int i = 0; i < fasta_indices_.size(); i++)
      fai_destroy(fasta_indices_[i]);
  }

  /*
   * Retrieves the sequence with name CHROM from the relevant FASTA file and stores it in SEQ
   */
  void get_sequence(const std::string& chrom, std::string& seq){
    std::string chrom_key = chrom;
    auto index_iter = chrom_to_index_.find(chrom);
    if (index_iter == chrom_to_index_.end()){
      if (chrom.size() > 3 && string_starts_with(chrom, "chr")){
	chrom_key  = chrom.substr(3);
	index_iter = chrom_to_index_.find(chrom_key);
      }
      if (index_iter == chrom_to_index_.end())
	printErrorAndDie("No entry for chromosome " + chrom + " found in FASTA files"); 
    }

    int length;
    char* result = fai_fetch(index_iter->second, chrom_key.c_str(), &length);
    assert(result != NULL);
    seq.assign(result, length);
    free((void *)result);
  }

  /*
   * Retrieves the sequence with name CHROM from the relevant FASTA file
   * and stores the 0-index based substring from START -> END (inclusive) in SEQ
   */
  void get_sequence(const std::string& chrom, int32_t start, int32_t end, std::string& seq){
    std::string chrom_key = chrom;
    auto index_iter = chrom_to_index_.find(chrom);
    if (index_iter == chrom_to_index_.end()){
      if (chrom.size() > 3 && string_starts_with(chrom, "chr")){
	chrom_key  = chrom.substr(3);
	index_iter = chrom_to_index_.find(chrom_key);
      }
      if (index_iter == chrom_to_index_.end())
      printErrorAndDie("No entry for chromosome " + chrom + " found in FASTA files");
    }
      
    int length;
    char* result = faidx_fetch_seq(index_iter->second, chrom_key.c_str(), start, end, &length);
    assert(result != NULL);
    seq.assign(result, length);
    free((void *)result);
  }
};

#endif
