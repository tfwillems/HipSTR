#ifndef SRC_NWNOREFENDPENALTY_H_
#define SRC_NWNOREFENDPENALTY_H_

#include <stdlib.h>
#include <stdio.h>

#include <algorithm>
#include <string>
#include <vector>

#include "bamtools/include/api/BamAux.h"

namespace NWNoRefEndPenalty { 
  void Align(const std::string& ref_seq, 
	     const std::string& read_seq,
	     std::string& ref_seq_al, 
	     std::string& read_seq_al,
	     float* score, 
	     std::vector<BamTools::CigarOp>& cigar_list);

  bool LeftAlign(const std::string& ref_seq, 
		 const std::string& read_seq,
		 std::string& ref_seq_al, 
		 std::string& read_seq_al,
		 float* score, 
		 std::vector<BamTools::CigarOp>& cigar_list);
}
#endif

