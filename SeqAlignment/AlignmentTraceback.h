#ifndef ALIGNMENT_TRACEBACK_H_
#define ALIGNMENT_TRACEBACK_H_

#include <string>

#include "AlignmentData.h"

class AlignmentTrace {
 private:
  std::string hap_aln_;    // Alignment string for read against its genotype's haplotype
  Alignment trace_vs_ref_; // Alignment trace relative to the reference allele
  int gt_;                 // Index for genotype against which the alignment was originally traced
  int num_flank_ins_;      // Number of inserted base pairs in sequences flanking the STR
  int num_flank_del_;      // Number of deleted base pairs in sequences flanking the STR
  int stutter_size_;       // Size of stutter artifact in STR block

 public:
  AlignmentTrace(int gt){
    gt_            = gt;
    hap_aln_       = "";
    num_flank_ins_ = -1;
    num_flank_del_ = -1;
    stutter_size_  = -999;
  }

  int gt_index()              { return gt_;            }
  int num_flank_ins()         { return num_flank_ins_; }
  int num_flank_del()         { return num_flank_del_; }
  int stutter_size()          { return stutter_size_;  }
  std::string& hap_aln()      { return hap_aln_;       }
  Alignment& traced_aln()     { return trace_vs_ref_;  }


  void set_num_flank_ins(int num_ins)     { num_flank_ins_;               }
  void set_num_flank_del(int num_del)     { num_flank_del_;               }
  void set_stutter_size(int stutter_size) { stutter_size_ = stutter_size; }
  void set_hap_aln(std::string& aln)      { hap_aln_ = aln;               }
};


std::string stitch(const std::string& hap_aln, const std::string& read_aln, int h_index, int r_index, int increment);

void stitch_alignment_trace(int32_t hap_start, const std::string& hap_aln_to_ref, 
			    const std::string& read_aln_to_hap, int hap_index, int seed_base, Alignment& orig_aln,
			    Alignment& new_aln);

#endif
