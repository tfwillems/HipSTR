#ifndef ALIGNMENT_DATA_H_
#define ALIGNMENT_DATA_H_

#include <sstream>
#include <string>
#include <vector>

#include "../error.h"

class CigarElement {
 private:
  char type_;
  int  num_;

 public:
  CigarElement(char type, int num){
    type_ = type;
    num_  = num;
  }

  inline void set_type(char type){ type_ = type;}
  inline void set_num(int num)   {  num_ = num; }
  inline char get_type()   const { return type_; }
  inline int  get_num()    const { return num_;  }
};

class Alignment {
 private:
  int32_t start_;
  int32_t stop_;
  std::string sample_;
  std::vector<CigarElement> cigar_list_;
  std::string base_qualities_;
  std::string sequence_;
  std::string alignment_;
  double mapq_;

 public:
  Alignment(int32_t start, int32_t stop,
	    std::string sample,
	    std::string base_qualities,
	    std::string sequence,
	    std::string alignment,
	    double mapq){
    start_          = start; 
    stop_           = stop;
    sample_         = sample; 
    base_qualities_ = base_qualities;
    sequence_       = sequence;
    alignment_      = alignment;
    mapq_           = mapq;
  }

  Alignment(){
    start_          = 0; 
    stop_           = -1;
    sample_         = ""; 
    base_qualities_ = "";
    sequence_       = "";
    alignment_      = "";
    mapq_           = 0.0;
  }

  inline int32_t get_start()             const { return start_;  }
  inline int32_t get_stop()              const { return stop_;   }
  inline std::string get_sample()        const { return sample_; }
  inline double get_mapping_quality()    const { return mapq_;   }

  inline void set_start(int32_t start)       { start_ = start;   }
  inline void set_stop(int32_t stop)         { stop_  = stop;    }
  inline void set_sample(std::string sample) { sample_ = sample; }

  int num_indels() const{
    int num = 0;
    for (std::vector<CigarElement>::const_iterator iter = cigar_list_.begin(); iter != cigar_list_.end(); iter++)
      if (iter->get_type() == 'I' || iter->get_type() == 'D')
	num++;
    return num;
  }
  
  int num_mismatches() const{
    int num = 0;
    for (std::vector<CigarElement>::const_iterator iter = cigar_list_.begin(); iter != cigar_list_.end(); iter++)
      if (iter->get_type() == 'X')
	num++;
    return num;
  }

  inline void set_base_qualities(std::string& base_qualities){ base_qualities_.assign(base_qualities); }
  inline void set_sequence(std::string& sequence)            { sequence_.assign(sequence);             }
  inline void set_alignment(std::string& alignment)          { alignment_.assign(alignment);           }
  inline void add_cigar_element(CigarElement e)              { cigar_list_.push_back(e);               }

  inline const std::string& get_base_qualities()           const { return base_qualities_; }
  inline const std::string& get_sequence()                 const { return sequence_;       }
  inline const std::string& get_alignment()                const { return alignment_;      }
  inline const std::vector<CigarElement>& get_cigar_list() const { return cigar_list_;     }

  void get_deletion_boundaries(std::vector<int32_t>& starts, std::vector<int32_t>& stops) const {
    int32_t pos = start_;
    for (auto iter = cigar_list_.begin(); iter != cigar_list_.end(); iter++){
      switch(iter->get_type()){
      case 'M': case 'X': case '=':
	pos += iter->get_num();
	break;
      case 'I': case 'S':
	break;
      case 'D':
	starts.push_back(pos);
	stops.push_back(pos+iter->get_num()-1);
	pos += iter->get_num();
	break;
      default:
	printErrorAndDie("Invalid CIGAR char detected in get_deletion_boundaries");
      }
    }
  }

  void get_insertion_positions(std::vector<int32_t>& positions, std::vector<int32_t>& sizes) const {
    int32_t pos = start_;
    for (auto iter = cigar_list_.begin(); iter != cigar_list_.end(); iter++){
      switch(iter->get_type()){
      case 'M': case 'X': case '=':
	pos += iter->get_num();
	break;
      case 'S':
	break;
      case 'I':
	positions.push_back(pos);
	sizes.push_back(iter->get_num());
	break;
      case 'D':
	pos += iter->get_num();
	break;
      default:
	printErrorAndDie("Invalid CIGAR char detected in get_deletion_boundaries");
      }
    }
  }
  
  std::string getCigarString(){
    std::stringstream cigar_str;
    for (std::vector<CigarElement>::iterator iter = cigar_list_.begin(); iter != cigar_list_.end(); iter++)
      cigar_str << iter->get_num() << iter->get_type();
    return cigar_str.str();
  }
};


bool compareAl(const Alignment& alignment_1, const Alignment& alignment_2);

void sortAlignments(std::vector<Alignment>& alignments);

#endif
