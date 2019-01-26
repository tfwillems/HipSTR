#include "AlignmentMatrixCache.h"

void AlignmentMatrixCache::add(int block_index, int block_option, double* match_matrix, double* insertion_matrix, double* deletion_matrix){
  if (has(block_index, block_option)) 
    printErrorAndDie("AlignmentMatrixCache instance already contains the requested entry");
  entries_[std::pair<int, int>(block_index, block_option)] = std::pair<int, int>(nonstutter_matrices_.size(), 0);
  nonstutter_matrices_.push_back(std::tuple<double*, double*, double*>(match_matrix, insertion_matrix, deletion_matrix));
}

void AlignmentMatrixCache::add(int block_index, int block_option, double* match_matrix, int* artifact_sizes, int* artifact_positions){
  if (has(block_index, block_option)) 
    printErrorAndDie("AlignmentMatrixCache instance already contains the requested entry");
  entries_[std::pair<int, int>(block_index, block_option)] = std::pair<int, int>(stutter_matrices_.size(), 1);
  stutter_matrices_.push_back(std::tuple<double*, int*, int*>(match_matrix, artifact_sizes, artifact_positions));
}

void AlignmentMatrixCache::get(int block_index, int block_option,
			       double*& match_matrix, double*& insertion_matrix, double*& deletion_matrix){
  std::pair<int, int>key(block_index, block_option);
  auto iter = entries_.find(key);
  if (iter == entries_.end())
    printErrorAndDie("AlignmentMatrixCache instance does not contain the requested entry");
  if (iter->second.second != 0)
    printErrorAndDie("Cached matrices don't have the requested types");
  
  int index        = iter->second.first;
  match_matrix     = std::get<0>(nonstutter_matrices_[index]);
  insertion_matrix = std::get<1>(nonstutter_matrices_[index]);
  deletion_matrix  = std::get<2>(nonstutter_matrices_[index]);
}

void AlignmentMatrixCache::get(int block_index, int block_option,
			       double*& match_matrix, int*& artifact_sizes, int*& artifact_positions){
  std::pair<int, int>key(block_index, block_option);
  auto iter = entries_.find(key);
  if (iter == entries_.end())
    printErrorAndDie("AlignmentMatrixCache instance does not contain the requested entry");
  if (iter->second.second != 1)
    printErrorAndDie("Cached matrices don't have the requested types");
  
  int index          = iter->second.first;
  match_matrix       = std::get<0>(stutter_matrices_[index]);
  artifact_sizes     = std::get<1>(stutter_matrices_[index]);
  artifact_positions = std::get<2>(stutter_matrices_[index]);
}

void AlignmentMatrixCache::clear(){
  entries_.clear();
  
  // Clear the non-stutter data structures
  for (auto tuple_iter = nonstutter_matrices_.begin(); tuple_iter != nonstutter_matrices_.end(); ++tuple_iter){
    delete [] std::get<0>(*tuple_iter);
    delete [] std::get<1>(*tuple_iter);
    delete [] std::get<2>(*tuple_iter);
  }
  nonstutter_matrices_.clear();
  
  // Clear the stutter data structures
  for (auto tuple_iter = stutter_matrices_.begin(); tuple_iter != stutter_matrices_.end(); ++tuple_iter){
    delete [] std::get<0>(*tuple_iter);
    delete [] std::get<1>(*tuple_iter);
    delete [] std::get<2>(*tuple_iter);
  }
  stutter_matrices_.clear();
}
