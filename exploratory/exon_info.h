#ifndef EXON_INFO_H_
#define EXON_INFO_H_

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "../error.h"
#include "../interval_tree.h"

class EXON {
 private:
  std::string chrom_;
  int32_t start_, end_;
  std::string gene_, transcript_;

 public:
  EXON(std::string chrom, int32_t start, int32_t end, std::string gene, std::string transcript){
    chrom_      = chrom;
    start_      = start;
    end_        = end;
    gene_       = gene;
    transcript_ = transcript;
  }

  const std::string& get_gene()    { return gene_; }
  void  set_gene(std::string gene) { gene_ = gene; }
};

class ExonInfo {
private:
  std::vector<std::string> gene_names_;
  std::vector<int> gene_cluster_ids_;                              // Gene cluster associated with each gene
  std::map<std::string, int> gene_indexes_;                        // Map from gene to associated index
  std::map<std::string, std::string> transcript_to_gene_;          // Map from transcript to its associated gene. Must be 1 gene per transcript
  std::vector< std::vector<std::string> > transcripts_by_gene_;    // Transcripts associated with each gene
  std::map<std::string, int> chrom_indexes_;                       // Map from chromosome to associated index
  std::vector< IntervalTree<EXON> > exon_trees_;                   // Interval tree for exons on each chromosome
  ExonInfo(){}
  
  // Parse text file to populate each of the fields
  void read_exon_info(std::string input_file);

  // Cluster genes with overlapping exons
  void cluster_genes(std::vector< std::vector< Interval<EXON> > >& exon_intervals_by_chrom);
   
public:
  ExonInfo(std::string input_file)     { read_exon_info(input_file);                                          }
  int  num_genes()                     { return gene_indexes_.size();                                         }
  int  num_transcripts()               { return transcript_to_gene_.size();                                   }
  bool is_gene(std::string& name)      { return gene_indexes_.find(name) != gene_indexes_.end();              }
  bool is_transcript(std::string& name){ return transcript_to_gene_.find(name) != transcript_to_gene_.end();  }

  std::string gene_for_transcript(std::string& transcript){
    auto trans_iter = transcript_to_gene_.find(transcript);
    if (trans_iter == transcript_to_gene_.end())
      printErrorAndDie("No transcript found with name " + transcript);
    return trans_iter->second;
  }

  int gene_cluster_for_gene(std::string gene){
    auto gene_iter = gene_indexes_.find(gene);
    if (gene_iter == gene_indexes_.end())
      printErrorAndDie("No gene information for gene " + gene);
    return gene_cluster_ids_[gene_iter->second];
  }

  int gene_cluster_for_transcript(std::string transcript){
    auto trans_iter = transcript_to_gene_.find(transcript);
    if (trans_iter == transcript_to_gene_.end())
      printErrorAndDie("No transcript information for transcript " + transcript);
    return gene_cluster_for_gene(trans_iter->second);
  }

  void gene_clusters_from_region(std::string& chrom, int32_t start, int32_t end, std::vector<int>& cluster_ids){
    assert(cluster_ids.size() == 0);
    auto chrom_iter = chrom_indexes_.find(chrom);
    if (chrom_iter == chrom_indexes_.end())
      return;
    std::vector< Interval<EXON> > exon_overlaps;
    exon_trees_[chrom_iter->second].findOverlapping(start, end, exon_overlaps);
    for (unsigned int i = 0; i < exon_overlaps.size(); i++)
      cluster_ids.push_back(gene_cluster_ids_[gene_indexes_[exon_overlaps[i].value.get_gene()]]);
  }
  
  void print_cluster_info(int cluster_id, std::ostream& out){
    out << "CLUSTER #" << cluster_id << ":";
    for (unsigned int i = 0; i < gene_cluster_ids_.size(); i++)
      if (gene_cluster_ids_[i] == cluster_id)
	out  << " " << gene_names_[i];
    out << "\n";
  }
};


#endif
