#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "../bamtools/include/api/BamAlignment.h"
#include "../bamtools/include/api/BamMultiReader.h"

#include "../base_quality.h"
#include "../error.h"
#include "../interval_tree.h"
#include "../snp_phasing_quality.h"
#include "../snp_tree.h"

const std::string BARCODE_TAG  = "BX";
const std::string MULTIMAP_TAG = "XA";
const uint32_t WINDOW_BP_SIZE  = 1000000;
const uint32_t READ_PAD        = 1000;
const int      MIN_MAPQ        = 30;

class BarcodePhasing {
public:
  std::string barcode;
  double log_p1, log_p2; 
  int32_t p1_match_count, p2_match_count,  mismatch_count;
  int32_t num_reads, num_reads_mapq;

  BarcodePhasing(std::string name){
    barcode        = name;
    log_p1         = 0;
    log_p2         = 0;
    p1_match_count = 0;
    p2_match_count = 0;
    mismatch_count = 0;
    num_reads      = 0;
    num_reads_mapq = 0;
  }

  void print(std::ostream& out){
    out << barcode         << "\t"   << log_p1         << "\t" << log_p2         << "\t" 
	<< p1_match_count  << "\t"   << p2_match_count << "\t" << mismatch_count << "\t" 
	<< num_reads       << "\t"   << num_reads_mapq << "\n";
  }
};


void build_PS_interval_trees(std::string& input_file, std::map<std::string, int>& chrom_indexes, std::vector< IntervalTree<int64_t> >& ps_interval_trees){
  assert(chrom_indexes.size() == 0 && ps_interval_trees.size() == 0);
  std::ifstream input(input_file.c_str());
  if (!input.is_open())
    printErrorAndDie("Failed to open PS region file");

  std::set<int64_t> ps_id_set;
  std::string line;
  std::vector< std::vector< Interval<int64_t> > > ps_intervals;
  int num_chroms = 0;
  while (std::getline(input, line)){
    std::istringstream iss(line);
    std::string chrom;
    int32_t start, stop;
    int64_t ps_id;
    if (!(iss >> chrom >> start >> stop >> ps_id)){
      std::cerr << line << std::endl;
      printErrorAndDie("Improperly formatted PS region file.");
    }
    if (ps_id_set.find(ps_id) != ps_id_set.end())
      printErrorAndDie("PS region file contains multiple regions with the same ID");
    else
      ps_id_set.insert(ps_id);    

    if (chrom_indexes.find(chrom) == chrom_indexes.end()){
      chrom_indexes[chrom] = num_chroms++;
      ps_intervals.push_back(std::vector< Interval<int64_t> >());
    }
    ps_intervals[chrom_indexes[chrom]].push_back(Interval<int64_t>(start, stop, ps_id));
  }
  input.close();

  for (unsigned int i = 0; i < ps_intervals.size(); i++)
    ps_interval_trees.push_back(IntervalTree<int64_t>(ps_intervals[i]));  
}


void print_barcodes(std::map<std::pair<int64_t, std::string>, BarcodePhasing>& barcode_info){
  for (auto iter = barcode_info.begin(); iter != barcode_info.end(); iter++){
    std::cout << iter->first.first << "_";
    iter->second.print(std::cout);
  }
  std::cout << std::endl;
}

int main(int argc, char** argv){
  bool filter_multimappers = true;

  std::string ps_interval_filename = "/data2/thomas/10x/filt_ps_regions.bed";
  std::vector< IntervalTree<int64_t> > ps_interval_trees;
  std::map<std::string, int> ps_chrom_indexes;
  build_PS_interval_trees(ps_interval_filename, ps_chrom_indexes, ps_interval_trees);

  std::string bam_filename = "/data2/thomas/10x/NA12878_phased_possorted_bam.bam";
  std::string vcf_filename = "/data2/thomas/10x/filt.filt_ps_regions.vcf.gz"; //"/data2/thomas/10x/filt.vcf.gz"; //"/data2/thomas/10x/NA12878_phased_variants.vcf.gz";
  BaseQuality base_qualities;

  std::vector<std::string> bam_files;
  bam_files.push_back(bam_filename);
  
  std::cerr << "Opening BAM" << std::endl;
  // Open all BAM files
  BamTools::BamMultiReader reader;
  if (!reader.Open(bam_files))
    printErrorAndDie("Failed to open one or more BAM files");

  std::cerr << "Opening VCF" << std::endl;
  // Open the VCF
  VCF::VCFReader vcf_reader(vcf_filename);

  std::cerr << "Iterating over regions" << std::endl;
  // Iterate over all reference sequences in the BAM in chunks
  const BamTools::RefVector ref_seqs = reader.GetReferenceData();
  std::vector<SNPTree*> snp_trees;
  std::map<std::pair<int64_t, std::string>, BarcodePhasing> barcode_info;
  int64_t read_count = 0;
  for (unsigned int i = 0; i < ref_seqs.size(); i++){
    uint32_t ref_len  = (uint32_t)ref_seqs[i].RefLength;
    std::string chrom = ref_seqs[i].RefName;
    
    bool sex_chrom = false;
    if (chrom.compare("chrX") == 0 || chrom.compare("X") == 0 || chrom.compare("x") == 0)
      sex_chrom = true;
    if (chrom.compare("chrY") == 0 || chrom.compare("y") == 0 || chrom.compare("y") == 0)
      sex_chrom = true;
    if (sex_chrom)
      continue;

    if (ps_chrom_indexes.find(chrom) == ps_chrom_indexes.end()){
      std::cerr << "Skipping chromosome " << chrom << " because no PS region information was available" << std::endl;
      continue;
    }
    int ps_index = ps_chrom_indexes[chrom];

    std::cerr << "Setting region" << std::endl;
    if(!reader.SetRegion(i, 0, i, ref_len)){
      printErrorAndDie("One or more BAM files failed to set the region properly");
    }
    std::cerr << "Done setting region" << std::endl;

    for (uint32_t start = 0; start < ref_len; start += WINDOW_BP_SIZE){
      uint32_t end = start + WINDOW_BP_SIZE;
      std::cerr << "Processing region " << chrom << " " << start << " " << end << ", # reads = " << read_count << std::endl;

      std::map<std::string, unsigned int> sample_indices;
      std::vector<Region> skip_regions;
      int32_t skip_pad = 0;
      if (create_snp_trees(chrom, start, end + READ_PAD, skip_regions, skip_pad, &vcf_reader, NULL, sample_indices, snp_trees, std::cerr)){
	std::cerr << "Built SNP tree" << std::endl;
	assert(snp_trees.size() > 0);  
	
	BamTools::BamAlignment alignment;
	while (reader.GetNextAlignment(alignment)){
	  read_count++;
	  if (alignment.Position < start)
	    continue;
	  if (alignment.Position > end)
	    break;

	  
	  std::vector< Interval<int64_t> > ps_ids;
	  ps_interval_trees[ps_index].findOverlapping(alignment.Position, alignment.GetEndPosition(), ps_ids);
	  if (ps_ids.size() != 1)
	    continue;
	  int64_t ps_id = ps_ids[0].value;

	  if (alignment.HasTag(BARCODE_TAG)){
	    std::string barcode;
	    alignment.GetTag(BARCODE_TAG, barcode);
	    std::pair<int64_t,  std::string> key(ps_id, barcode);
	    
	    auto barcode_iter = barcode_info.find(key);
	    if (barcode_iter == barcode_info.end()){
	      barcode_info.insert(std::pair< std::pair<int64_t,  std::string>, BarcodePhasing>(key, BarcodePhasing(barcode)));
	      barcode_iter          = barcode_info.find(key);
	    }

	    BarcodePhasing& bc = barcode_iter->second;
	    bc.num_reads++;

	    // Skip if filtering multimappers
	    if (filter_multimappers && alignment.HasTag(MULTIMAP_TAG))
	      continue;

	    // Use SNP tree to look for overlapping heterozygous SNPs
	    if (alignment.MapQuality > MIN_MAPQ){
	      bc.num_reads_mapq++;
	      add_log_phasing_probs(alignment, snp_trees[0], base_qualities,
				    bc.log_p1, bc.log_p2, bc.p1_match_count, bc.p2_match_count, bc.mismatch_count);
	    }
	  }
	}
      }
      else {
	std::cerr << "Failed to build SNP tree" << std::endl;
      }
      destroy_snp_trees(snp_trees);
    }
    //print_barcodes(barcode_info);
  }
  print_barcodes(barcode_info);
}
