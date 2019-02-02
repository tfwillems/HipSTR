#ifndef ALIGNMENT_MATRIX_CACHE_H_
#define ALIGNMENT_MATRIX_CACHE_H_

#include <map>
#include <tuple>
#include <vector>

#include "../error.h"

class AlignmentMatrixCache {
 private:
  // Key = (block_index, block_option), Value = (vector index, vector type)
  std::map< std::pair<int, int>, std::pair<int, int> > entries_;

  // Private unimplemented copy constructor and assignment operator to prevent operations
  AlignmentMatrixCache(const AlignmentMatrixCache& other);
  AlignmentMatrixCache& operator=(const AlignmentMatrixCache& other);
  
   // Cached matrices
  std::vector< std::tuple<double*, double*, double*> > nonstutter_matrices_; // Type 0
  std::vector< std::tuple<double*,    int*,    int*> > stutter_matrices_;    // Type 1

 public:
  AlignmentMatrixCache(){}

  bool has(int block_index, int block_option){
    std::pair<int, int> key(block_index, block_option);
    return entries_.find(key) != entries_.end();
  }

  void add(int block_index, int block_option, double* match_matrix, double* insertion_matrix, double* deletion_matrix);

  void add(int block_index, int block_option, double* match_matrix, int* artifact_sizes, int* artifact_positions);

  void get(int block_index, int block_option,
	   double*& match_matrix, double*& insertion_matrix, double*& deletion_matrix);

  void get(int block_index, int block_option,
	   double*& match_matrix, int*& artifact_sizes, int*& artifact_positions);

  void reindex(const std::map< std::pair<int, int>, int>& block_allele_mapping);

  void clear();
  
  ~AlignmentMatrixCache(){
    clear();
  }
};

#endif
