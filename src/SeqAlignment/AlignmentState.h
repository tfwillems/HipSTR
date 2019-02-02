#ifndef ALIGNMENT_STATE_H_
#define ALIGNMENT_STATE_H_

#include <string>
#include <vector>

#include "AlignmentMatrixCache.h"
#include "AlignmentTraceback.h"
#include "../error.h"
#include "Haplotype.h"

class AlignmentState {
 private:
  // Haplotype relative to which the alignment is made
  Haplotype* hap_;

  // Sequence and quality score information for the read
  const char* seq_;
  const int seq_len_;
  const int read_id_;
  const double* log_wrong_;
  const double* log_right_;
  const int seed_;

  // Counters that describe the current endpoint of the alignment if one has been set
  int hap_index_;
  int seq_index_;
  int matrix_index_;
  int matrix_type_;

  // Pointers to the start of various alignment matrices
  double* match_matrix_;
  double* ins_matrix_;
  double* del_matrix_;
  int* best_artifact_size_;
  int* best_artifact_pos_;

  // Static members
  const static double IMPOSSIBLE;               // Large negative value to prevent impossible or undesirable configurations
  const static double MIN_SNP_LOG_PROB_CORRECT; // Only consider a base as a SNP if the log prob is above this threshold
  
  // Types of alignment matrices
  const static int NONE = -1;
  const static int MATCH = 0;
  const static int DEL   = 1;
  const static int INS   = 2;

  void align_seq_to_stutter_block(const int block_index,
				  double* match_matrix, int* best_artifact_size, int* best_artifact_pos,
				  const double* prev_match_matrix, const bool nonspanning_only);

  void align_seq_to_nonstutter_block(const int block_index,
				     double* match_matrix, double* ins_matrix, double* del_matrix,
				     const double* prev_match_matrix);

  std::string retrace_helper(AlignmentTrace& trace);
  std::string retrace(AlignmentTrace& trace);

  std::string stitch(const std::string& hap_aln, const std::string& read_aln, int h_index, int r_index, int increment);

  // Private unimplemented copy constructor and assignment operator to prevent operations
  AlignmentState(const AlignmentState& other);
  AlignmentState& operator=(const AlignmentState& other);

 public:
 AlignmentState(Haplotype* haplotype, const char* seq, const int seq_len, const int read_id,
		const double* log_wrong, const double* log_right, int seed):
  hap_(haplotype), seq_(seq), seq_len_(seq_len), read_id_(read_id), log_wrong_(log_wrong), log_right_(log_right), seed_(seed){
    int max_hap_size    = hap_->max_size();
    int num_hap_blocks  = hap_->num_blocks();
    int max_matrix_size = seq_len_*(max_hap_size + num_hap_blocks + 1); // +1 for padding row
    match_matrix_       = new double [max_matrix_size];
    ins_matrix_         = new double [max_matrix_size];
    del_matrix_         = new double [max_matrix_size];
    best_artifact_size_ = new int    [seq_len_*num_hap_blocks];
    best_artifact_pos_  = new int    [seq_len_*num_hap_blocks];
    reset_traceback();
  }

  void reset_traceback(){
    hap_index_    = -1;
    seq_index_    = -1;
    matrix_index_ = -1;
    matrix_type_  = NONE;
  }

  void align_seq_to_haplotype(const bool reuse_alns, AlignmentMatrixCache* matrix_cache);

  friend double calc_maximum_likelihood_alignment(AlignmentState& fw_state, AlignmentState& rv_state);

  friend void stitch(AlignmentState& fw_state, AlignmentState& rv_state, const Alignment& orig_aln,
		     AlignmentTrace& trace, bool debug);

  ~AlignmentState(){
    // Deallocate full scoring matrices
    delete [] match_matrix_;
    delete [] ins_matrix_;
    delete [] del_matrix_;
    delete [] best_artifact_size_;
    delete [] best_artifact_pos_;
  }
};

#endif
