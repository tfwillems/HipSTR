#ifndef ALIGNMENT_TRACEBACK_H_
#define ALIGNMENT_TRACEBACK_H_

#include <string>
#include <vector>

#include "AlignmentData.h"
#include "Haplotype.h"

class AlignmentTrace {
 private:
  std::string hap_aln_;      // Alignment string for read against its genotype's haplotype
  Alignment trace_vs_ref_;   // Alignment trace relative to the reference allele
  int gt_;                   // Index for genotype against which the alignment was originally traced
  int flank_ins_size_;       // Number of inserted base pairs in sequences flanking the STR (positive)
  int flank_del_size_;       // Number of deleted base pairs in sequences flanking the STR (positive)
  int stutter_size_;         // Size of stutter artifact in STR block
  std::string str_seq_;      // Sequence in STR region
  std::string full_str_seq_; // Hypothetical sequence in STR region if the read fully spanned the stutter block
                             // Identical to str_seq if read spans STR region
  std::vector< std::pair<int,int> > flank_indel_data_;

 public:
  AlignmentTrace(int gt){
    gt_             = gt;
    hap_aln_        = "";
    flank_ins_size_ = -1;
    flank_del_size_ = -1;
    stutter_size_   = -999;
    str_seq_        = "";
  }

  int gt_index()              { return gt_;             }
  int flank_ins_size()        { return flank_ins_size_; }
  int flank_del_size()        { return flank_del_size_; }
  int stutter_size()          { return stutter_size_;   }
  std::string& hap_aln()      { return hap_aln_;        }
  Alignment& traced_aln()     { return trace_vs_ref_;   }
  std::string& str_seq()      { return str_seq_;        }
  std::string& full_str_seq() { return full_str_seq_;   }
  std::vector< std::pair<int,int> >& flank_indel_data() { return flank_indel_data_; }

  void set_flank_ins_size(int num_ins_bp) { flank_ins_size_ = num_ins_bp;   }
  void set_flank_del_size(int num_del_bp) { flank_del_size_ = num_del_bp;   }
  void set_stutter_size(int stutter_size) { stutter_size_   = stutter_size; }
  void set_hap_aln(std::string& aln)      { hap_aln_        = aln;          }
  void set_str_seq(std::string& seq)      { str_seq_        = seq;          }
  void set_full_str_seq(std::string& seq) { full_str_seq_   = seq;          }
  void set_flank_indel_data(std::vector< std::pair<int,int> >& data){ flank_indel_data_ = data; }
};


std::string stitch(const std::string& hap_aln, const std::string& read_aln, int h_index, int r_index, int increment);

void stitch_alignment_trace(int32_t hap_start, const std::string& hap_aln_to_ref, 
			    const std::string& read_aln_to_hap, int hap_index, int seed_base, Alignment& orig_aln,
			    Alignment& new_aln);

#endif
