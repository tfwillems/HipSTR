#ifndef HAPLOTYPE_CREATOR_H_
#define HAPLOTYPE_CREATOR_H_

#include <string>
#include <vector>

#include "AlignmentData.h"
#include "RepeatRegion.h"

/* 
   Determines the disjoint set of regions within which indels, mismatches or deletions occur. Regions that are overlapping or adjacent
   are combined into a single larger region. Each region is represented as [start, end) and stored in the provided
   vector in order of increasing start coordinate
*/
void getHaplotypeRegions(std::vector<Alignment>& alignments, 
			 std::vector< std::pair<int32_t, int32_t> >& regions);
 
/* 
   Extracts the set of sequences observed for each of the regions in the provided alignments.
   An alignment's sequence for a region is only considered if the alignment 
   fully spans the region and the reference sequence for each region is also added
   to each set.
 */
void extractRegionSequences(std::vector<Alignment>& alignments, 
			    std::vector< std::pair<int,int> >& regions, 
			    std::string& ref_sequence, 
			    std::vector< std::vector<std::string> >& region_seqs);


/*
   Merge regions which overlap the provided repeat interval into a single region
 */
void mergeRegions(int32_t str_start, int32_t str_end, std::vector< std::pair<int32_t, int32_t> >& regions);


#endif
