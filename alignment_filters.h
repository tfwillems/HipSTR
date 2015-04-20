/*
Copyright (C) 2014 Thomas Willems <twillems@mit.edu>
*/


#ifndef SRC_ALIGNMENTFILTERS_H_
#define SRC_ALIGNMENTFILTERS_H_

#include <string>
#include <vector>

#include "bamtools/include/api/BamAlignment.h"

namespace AlignmentFilters {
  /* Returns the CIGAR string corresponding to the vector of CigarOps. */
  std::string GetCigarString(std::vector<BamTools::CigarOp>& cigar_ops); 
  
  /* Length of perfect base matches at 5' and 3' end of read. */
  std::pair<int,int> GetNumEndMatches(BamTools::BamAlignment& aln, const std::string& ref_seq, int ref_seq_start);
  
  /* Minimum distances from 5' and 3' end of reads to first indel. If no such indel exists, returns (-1,-1). */
  std::pair<int,int> GetEndDistToIndel(BamTools::BamAlignment& aln);

  /* Returns true iff the alignment has: 
     1) a maximal matching prefix compared to alignments that start [-max_upstream, max_downstream] from the 5' alignment position of the read
     2) a maximal matching suffix compared to alignments that end   [-max_downstream, max_upstream] from the 3' alignment position of the read
     Ignores clipped bases when performing these comparions 
  */
  bool HasLargestEndMatches(BamTools::BamAlignment& aln, const std::string& ref_seq, int ref_seq_start, int max_upstream, int max_downstream);
}

#endif
