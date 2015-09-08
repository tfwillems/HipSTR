#include <algorithm>
#include <assert.h>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "../error.h"
#include "exon_info.h"

void ExonInfo::cluster_genes(std::vector< std::vector< Interval<EXON> > >& exon_intervals_by_chrom){
  assert(gene_cluster_ids_.size() == 0);
  
  // Initialize each gene with its own cluster
  for (unsigned int i = 0; i < gene_indexes_.size(); i++)
    gene_cluster_ids_.push_back(i);
  
  // Combine genes that overlap into clusters
  for (unsigned int i = 0; i < exon_intervals_by_chrom.size(); i++){
    std::sort(exon_intervals_by_chrom[i].begin(), exon_intervals_by_chrom[i].end(), IntervalStartSorter<EXON>());
    int32_t max_stop = -1;
    std::vector<std::string> genes_in_cluster;
    for (unsigned int j = 0; j < exon_intervals_by_chrom[i].size(); j++){
      if (intervalStart(exon_intervals_by_chrom[i][j]) < max_stop){
	genes_in_cluster.push_back(exon_intervals_by_chrom[i][j].value.get_gene());
	max_stop = std::max(max_stop, intervalStop(exon_intervals_by_chrom[i][j]));
      }
      else {
	while (genes_in_cluster.size() > 1){
	  int id_1 = gene_cluster_ids_[gene_indexes_[genes_in_cluster[0]]];	  
	  while (gene_cluster_ids_[id_1] != id_1){
	    gene_cluster_ids_[id_1] = gene_cluster_ids_[gene_cluster_ids_[id_1]];
	    id_1 = gene_cluster_ids_[id_1];
	  }
	  int id_2 = gene_cluster_ids_[gene_indexes_[genes_in_cluster.back()]];
	  while (gene_cluster_ids_[id_2] != id_2){
	    gene_cluster_ids_[id_2] = gene_cluster_ids_[gene_cluster_ids_[id_2]];
	    id_2 = gene_cluster_ids_[id_2];
	  }
	  gene_cluster_ids_[id_2] = id_1;
	  genes_in_cluster.pop_back();
	}
	genes_in_cluster.clear();
	max_stop = intervalStop(exon_intervals_by_chrom[i][j]);
      }
    }
  }

  // Fully flatten cluster id structure
  std::set<int> final_ids;
  for (unsigned int i = 0; i < gene_cluster_ids_.size(); i++){
    int id = i;
    while (gene_cluster_ids_[id] != id){
      gene_cluster_ids_[id] = gene_cluster_ids_[gene_cluster_ids_[id]];
      id = gene_cluster_ids_[id];
    }
    gene_cluster_ids_[i] = id;
    final_ids.insert(id);
  }

  // Reindex the cluster ids to go from 0 to N-1
  std::vector<int> final_id_list(final_ids.begin(), final_ids.end());
  std::map<int,int> final_id_mapping;
  std::sort(final_id_list.begin(), final_id_list.end());
  for (unsigned int i = 0; i < final_id_list.size(); i++)
    final_id_mapping[final_id_list[i]] = i;
  for (unsigned int i = 0; i < gene_cluster_ids_.size(); i++)
    gene_cluster_ids_[i] = final_id_mapping[gene_cluster_ids_[i]];
  
  std::cerr << "After clustering overlapping genes, there are " << final_id_list.size() << " distinct clusters" << std::endl;
}

void ExonInfo::read_exon_info(std::string input_file){
  std::ifstream input(input_file.c_str());
  if (!input.is_open()) 
    printErrorAndDie("Failed to open exon annotation file " + input_file);

  int num_chroms = 0, num_genes = 0;
  std::string line;
  std::vector< std::vector< Interval<EXON> > > exon_intervals_by_chrom;
  std::vector<std::string> gene_names;
  while (std::getline(input, line)){
    std::istringstream iss(line);
    std::string chrom, gene, transcript;
    int32_t start, end;
    if (!(iss >> chrom >> start >> end >> gene >> transcript))
      printErrorAndDie("Improperly formatted exon annotation file. Required format is tab-delimited columns CHROM START END GENE TRANSCRIPT");
    if (chrom_indexes_.find(chrom) == chrom_indexes_.end()){
      chrom_indexes_[chrom] = num_chroms++;
      exon_intervals_by_chrom.push_back(std::vector< Interval<EXON> >());
    }
    if (gene_indexes_.find(gene) == gene_indexes_.end()){
      gene_names_.push_back(gene);
      gene_indexes_[gene] = num_genes++;
      transcripts_by_gene_.push_back(std::vector<std::string>());
    }
    if (transcript_to_gene_.find(transcript) == transcript_to_gene_.end()){
      transcript_to_gene_[transcript] = gene;
      transcripts_by_gene_[gene_indexes_[gene]].push_back(transcript);
    }
    else if (transcript_to_gene_[transcript].compare(gene) != 0)
      printErrorAndDie("Transcript " + transcript + " is already associated with a gene");

    int chrom_index = chrom_indexes_.find(chrom)->second;
    exon_intervals_by_chrom[chrom_index].push_back(Interval<EXON>(start, end, EXON(chrom, start, end, gene, transcript)));
  }
  input.close();
  std::cerr << "Exon annotation file contains " << num_genes << " genes across " << num_chroms << " chromosomes" << std::endl;

  // Cluster genes with overlapping exons
  std::vector<int> gene_cluster_ids;
  cluster_genes(exon_intervals_by_chrom);

  // Construct the interval trees
  assert(exon_trees_.size() == 0);
  for (unsigned int i = 0; i < exon_intervals_by_chrom.size(); i++)
    exon_trees_.push_back(IntervalTree<EXON>(exon_intervals_by_chrom[i]));
}
