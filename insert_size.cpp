#include <iomanip>  
#include <math.h>

#include "error.h"
#include "insert_size.h"

void InsertSizeCounter::process_alignment(BamTools::BamAlignment& aln){
  total_comps_++;
  if (aln.RefID == -1 || aln.MateRefID == -1 || aln.InsertSize == 0)
    unmapped_count_++;
  else if (aln.RefID != aln.MateRefID)
    diff_chrom_count_++;
  else if (abs(aln.InsertSize) > max_diff_)
    max_count_++;
  else {
    if (!add_diff(abs(aln.InsertSize)))
      printErrorAndDie("Failed to add insert size for alignment " + aln.Name);
    valid_comps_++;
    sum_diffs_        += abs(aln.InsertSize);
    sum_square_diffs_ += aln.InsertSize*aln.InsertSize;
  }
}

void InsertSizeCounter::output_summary_statistics(std::ostream& out){
  out << "UNMAPPED"    << "\t" << unmapped_count_   << "\t" << std::setprecision(3) << 100.0*unmapped_count_/total_comps_   << "\n"
      << "DIFF_CHROM"  << "\t" << diff_chrom_count_ << "\t" << std::setprecision(3) << 100.0*diff_chrom_count_/total_comps_ << "\n"
      << "GT_MAX_DIST" << "\t" << max_count_        << "\t" << std::setprecision(3) << 100.0*max_count_/total_comps_        << "\n"
      << "VALID_PAIRS" << "\t" << valid_comps_      << "\t" << std::setprecision(3) << 100.0*valid_comps_/total_comps_      << "\n"
      << "TOTAL_PAIRS" << "\t" << total_comps_      << "\t" << std::setprecision(3) << 100.0*total_comps_/total_comps_      << "\n";
  
  out << "STATISTICS FOR VALID PAIRS:" << "\n";
  
  // Mean
  double mean = 1.0*sum_diffs_/valid_comps_;
  out << "\t" << "MEAN" << "\t" << mean << "\n";
  
  // Variance
  double variance = 1.0*sum_square_diffs_/valid_comps_ - mean*mean; 
  out << "\t" << "STD" << "\t" << sqrt(variance) << "\n";
  
  // Histogram counts
  for (int i = 0; i < 20; i++)
    out << bins_[i+1] << "\t" << counts_[i] << "\t" << std::setprecision(3) << 100.0*counts_[i]/valid_comps_ << "\n";
  out << std::endl;
}
