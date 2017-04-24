#ifndef SRC_ALIGNMENTFILTERS_H_
#define SRC_ALIGNMENTFILTERS_H_

#include <string>
#include <vector>

#include "bam_io.h"

namespace AlignmentFilters {
  /* Returns the CIGAR string corresponding to the vector of CigarOps. */
  std::string GetCigarString(const std::vector<CigarOp>& cigar_ops);
  
  /* Length of perfect base matches at 5' and 3' end of read. */
  std::pair<int,int> GetNumEndMatches(BamAlignment& aln, const std::string& ref_seq, int ref_seq_start);
  
  /* Minimum distances from 5' and 3' end of reads to first indel. If no such indel exists, returns (-1,-1). */
  std::pair<int,int> GetEndDistToIndel(BamAlignment& aln);

  /* Returns true iff the alignment has: 
     1) a maximal matching prefix compared to alignments that start [-max_upstream, max_downstream] from the 5' alignment position of the read
     2) a maximal matching suffix compared to alignments that end   [-max_downstream, max_upstream] from the 3' alignment position of the read
     Ignores clipped bases when performing these comparions 
  */
  bool HasLargestEndMatches(BamAlignment& aln, const std::string& ref_seq, int ref_seq_start, int max_upstream, int max_downstream);


  /* Determines the number of hard and soft clipped bases based on the CIGAR string and stores them in the provided integer references */
  void GetNumClippedBases(BamAlignment& aln, int& num_hard_clips, int& num_soft_clips);
}

#endif
