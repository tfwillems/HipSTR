#ifndef ALIGNMENT_TRACEBACK_H_
#define ALIGNMENT_TRACEBACK_H_

#include <string>
#include <vector>

#include "AlignmentData.h"
#include "Haplotype.h"

class AlignmentTrace {
 private:
  // Simple class for traceback data related to STR blocks
  class STRTraceData {
  private:
    int stutter_size_;         // Size of stutter artifact in STR block
    std::string str_seq_;      // Sequence in STR region
  public:
    STRTraceData(int stutter_size, std::string& str_seq){
      stutter_size_ = stutter_size;
      str_seq_      = str_seq;
    }

    int stutter_size()          { return stutter_size_; }
    std::string& str_seq()      { return str_seq_;      }
  };

  std::string hap_aln_;      // Alignment string for read against its genotype's haplotype
  Alignment trace_vs_ref_;   // Alignment trace relative to the reference allele
  int flank_ins_size_;       // Number of inserted base pairs in sequences flanking the STR (positive)
  int flank_del_size_;       // Number of deleted base pairs in sequences flanking the STR (positive)
  std::vector<STRTraceData*> str_data_;
  std::vector<std::string> flank_seqs_;
  std::vector< std::pair<int32_t,int32_t> > flank_indel_data_;
  std::vector< std::pair<int32_t, char> > flank_snp_data_;

 public:
 AlignmentTrace(int num_haplotype_blocks): trace_vs_ref_("TRACE"){
    hap_aln_        = "";
    flank_ins_size_ = 0;
    flank_del_size_ = 0;
    str_data_       = std::vector<STRTraceData*>(num_haplotype_blocks, NULL);
    flank_seqs_     = std::vector<std::string>(num_haplotype_blocks, "");
  }

  ~AlignmentTrace(){
    for (auto str_trace_iter = str_data_.begin(); str_trace_iter != str_data_.end(); ++str_trace_iter)
      if (*str_trace_iter != NULL)
	delete *str_trace_iter;
    str_data_.clear();
  }

  int flank_ins_size()     { return flank_ins_size_; }
  int flank_del_size()     { return flank_del_size_; }
  std::string& hap_aln()   { return hap_aln_;        }
  Alignment& traced_aln()  { return trace_vs_ref_;   }

  void add_flank_indel(std::pair<int32_t, int32_t> indel){ flank_indel_data_.push_back(indel); }
  void add_flank_snp(int32_t pos, char base){
    flank_snp_data_.push_back(std::pair<int32_t, char>(pos, base));
  }
  void inc_flank_ins()                    { flank_ins_size_++;   }
  void inc_flank_del()                    { flank_del_size_++;   }
  void set_hap_aln(std::string& aln)      { hap_aln_ = aln;      }

  void add_flank_data(int block_index, std::string& flank_seq){
    flank_seqs_[block_index].append(flank_seq);
  }

  void add_str_data(int block_index, int stutter_size, std::string& str_seq){
    assert(str_data_[block_index] == NULL);
    str_data_[block_index] = new STRTraceData(stutter_size, str_seq);
  }

  std::vector< std::pair<int32_t,int32_t> >& flank_indel_data() { return flank_indel_data_; }
  std::vector< std::pair<int32_t, char> >& flank_snp_data()     { return flank_snp_data_;   }

  bool has_stutter(){
    for (unsigned int i = 0; i < str_data_.size(); ++i)
      if (str_data_[i] != NULL)
	if (str_data_[i]->stutter_size() != 0)
	  return true;
    return false;
  }

  int total_stutter_size(){
    int total_size = 0;
    for (unsigned int i = 0; i < str_data_.size(); ++i)
      if (str_data_[i] != NULL)
	total_size += str_data_[i]->stutter_size();
    return total_size;
  }

  int stutter_size(int block_index){
    assert(str_data_[block_index] != NULL);
    return str_data_[block_index]->stutter_size();
  }

  std::string& flank_seq(int block_index){
    return flank_seqs_[block_index];
  }

  std::string& str_seq(int block_index){
    assert(str_data_[block_index] != NULL);
    return str_data_[block_index]->str_seq();
  }
};


std::string stitch(const std::string& hap_aln, const std::string& read_aln, int h_index, int r_index, int increment);

void stitch_alignment_trace(int32_t hap_start, const std::string& hap_aln_to_ref, 
			    const std::string& read_aln_to_hap, int hap_index, int seed_base, Alignment& orig_aln,
			    Alignment& new_aln);

#endif
