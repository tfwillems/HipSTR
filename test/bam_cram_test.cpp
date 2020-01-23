#include <time.h>

#include <iostream>
#include <string>
#include <vector>

#include "../src/bam_io.h"
#include "../src/null_ostream.h"
#include "../src/region.h"

double test_readers(const std::vector<Region>& regions, std::string fasta_ref_path,
			std::string bam_cram_path_one, bool reopen_one,
			std::string bam_cram_path_two, bool reopen_two,
			bool verbose, int32_t& total_alns){
  double start_time = clock();
  total_alns = 0;
  BamAlignment aln_one, aln_two;

  BamCramReader* reader_one = new BamCramReader(bam_cram_path_one, fasta_ref_path);
  BamCramReader* reader_two = new BamCramReader(bam_cram_path_two, fasta_ref_path);  

  // Iterate over each region in the file
  for (auto region_iter = regions.begin(); region_iter != regions.end(); ++region_iter){
    if (verbose)
      std::cerr << "Processing region " << region_iter->chrom() << " " << region_iter->start() << " " << region_iter->stop() << "\t";

    if (reopen_one){
      delete reader_one;
      reader_one = new BamCramReader(bam_cram_path_one, fasta_ref_path);
    }

    if (reopen_two){
      delete reader_two;
      reader_two = new BamCramReader(bam_cram_path_two, fasta_ref_path);
    }

    if (!reader_one->SetRegion(region_iter->chrom(), region_iter->start(), region_iter->stop()))
      printErrorAndDie("1st reader failed to set the region");
    if (!reader_two->SetRegion(region_iter->chrom(), region_iter->start(), region_iter->stop()))
      printErrorAndDie("2nd reader failed to set the region");

    // Ensure the BAM and CRAM readers result in identical numbers of alignments
    int32_t bam_count = 0;
    while (reader_one->GetNextAlignment(aln_one)){
      bam_count++;

      if (!reader_two->GetNextAlignment(aln_two))
        printErrorAndDie("1st reader has more alignments than 2nd reader");

      // Compare each pair of alignments to ensure that their contents are identical
      if (aln_one.Name().compare(aln_two.Name()) != 0)
        printErrorAndDie("Alignments' names differ");
      if (aln_one.QueryBases().compare(aln_two.QueryBases()) != 0){
	std::cerr << aln_one.Name()       << "\n"
		  << aln_two.Name()       << "\n"
                  << aln_one.QueryBases() << "\n"
                  << aln_two.QueryBases() << std::endl;
        printErrorAndDie("Alignments' sequences differ");
      }
    }
    while (reader_two->GetNextAlignment(aln_two))
      printErrorAndDie("2nd reader has more alignments than 1st reader");
    if (verbose)
      std::cerr << bam_count << " " << std::endl;
    total_alns += bam_count;
  }
  
  delete reader_one;
  delete reader_two;

  // Test passed
  double seconds = (clock() - start_time)/CLOCKS_PER_SEC;
  return seconds;
}

double test_single_reader(const std::vector<Region>& regions, std::string fasta_ref_path, std::string bam_cram_path, bool verbose, int32_t& total_alns){
  double start_time = clock();
  total_alns = 0;
  BamAlignment aln;
  BamCramReader reader(bam_cram_path, fasta_ref_path);

  // Iterate over each region in the file
  for (auto region_iter = regions.begin(); region_iter != regions.end(); ++region_iter){
    if (verbose)
      std::cerr << "Processing region " << region_iter->chrom() << " " << region_iter->start() << " " << region_iter->stop() << "\t";

    if (!reader.SetRegion(region_iter->chrom(), region_iter->start(), region_iter->stop()))
      printErrorAndDie("Reader failed to set the region");

    int32_t count = 0;
    while (reader.GetNextAlignment(aln))
      count++;

    if (verbose)
      std::cerr << count << std::endl;
    total_alns += count;
  }

  // Test passed
  double seconds = (clock() - start_time)/CLOCKS_PER_SEC;
  return seconds;
}

void run_tests(const std::vector<Region>& regions, bool verbose, std::string& bam_path, std::string cram_path, std::string& cram_fasta_path, std::string& region_file){
  double time;
  int32_t aln_count;
  std::cerr << "Running tests using a list of " << regions.size() << " regions" << std::endl;

  // Test that a single BAM reader works
  time = test_single_reader(regions, "", bam_path, verbose, aln_count);
  std::cerr << "\t" << "BAM test completed in " << time << " seconds and read " << aln_count << " records" << std::endl;

  // Test that a single CRAM reader works
  time = test_single_reader(regions, cram_fasta_path, cram_path, verbose, aln_count);
  std::cerr << "\t" << "CRAM test completed in " << time << " seconds and read " << aln_count << " records" << std::endl;

  // Test that a BAM and CRAM reader give identical alignments
  time = test_readers(regions, cram_fasta_path, bam_path, false, cram_path, false, verbose, aln_count);
  std::cerr << "\t" << "BAM vs. CRAM test completed in " << time << " seconds and read " << aln_count << " records" << std::endl;

  // Test that reopening a BAM reader for each region is equivalent to keeping the BAM reader open
  time = test_readers(regions, "", bam_path, false, bam_path, true, verbose, aln_count);
  std::cerr << "\t" << "BAM vs. BAM reopening test completed in " << time << " seconds and read " << aln_count << " records" << std::endl;
  
  std::cerr << std::endl;
}

int main(){
  std::string bam_path("/cluster/moredata/willems/Polaris/Backup_HipSTR/NA12878.bam");
  std::string cram_path("/cluster/moredata/willems/Polaris/Backup_HipSTR/NA12878.final.cram");
  std::string cram_fasta_path("/cluster/home/willems/reference_genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa");
  std::string region_file("/cluster/moredata/willems/Benchmarking/reference/all_annot_markers.hg38.fixed.pad_1kb.bed");

  // Read file containing regions
  std::vector<Region> regions;
  std::string chrom_limit = "";
  NullOstream logger;  
  readRegions(region_file, 1000000, chrom_limit, regions, logger);

  // Run tests on initial regions
  bool verbose = false;
  run_tests(regions, verbose, bam_path, cram_path, cram_fasta_path, region_file);

  // Add more overlapping regions
  std::vector<Region> regions_2, regions_3;
  size_t num_regions = regions.size();
  for (size_t i = 0; i < num_regions; ++i){
    Region r = regions[i];

    regions_2.push_back(r);
    regions_2.push_back(Region(r.chrom(), r.start()+10, r.stop()+50, r.period()));

    regions_3.push_back(r);
    regions_3.push_back(Region(r.chrom(), r.start()+10, r.stop()+50, r.period()));
    regions_3.push_back(Region(r.chrom(), r.start()+20, r.stop()+70, r.period()));
  }
  orderRegions(regions_2);
  orderRegions(regions_3);

  // Re-run tests using overlapping region sets
  run_tests(regions_2, verbose, bam_path, cram_path, cram_fasta_path, region_file);
  run_tests(regions_3, verbose, bam_path, cram_path, cram_fasta_path, region_file);
}

