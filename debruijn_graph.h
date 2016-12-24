#ifndef DEBRUIJN_GRAPH_
#define DEBRUIJN_GRAPH_

#include <assert.h>
#include <algorithm>
#include <string>
#include <vector>

#include "directed_graph.h"

class DebruijnPath;

class DebruijnGraph : public DirectedGraph {
 protected:
  int k_;
  std::string ref_seq_;
  std::string source_kmer_;
  std::string sink_kmer_;
  int32_t num_strings_;
  std::vector<bool> ref_edge_; // True iff the edge at the corresponding index is from the reference sequence

  void get_alt_kmer_nodes(std::string& kmer, bool source, bool sink, std::vector<Node*>& nodes);

  void prune_edges(std::vector<bool>& remove_edges);

 public:
  DebruijnGraph(int k, std::string& ref_seq){
    assert(ref_seq.size() > k);
    k_           = k;
    ref_seq_     = ref_seq;
    source_kmer_ = ref_seq.substr(0, k_);
    sink_kmer_   = ref_seq.substr(ref_seq.size()-k, k_);
    num_strings_ = 0;

    // Add the reference path with a weight of 2
    add_string(ref_seq, 2);
    ref_edge_.clear();
    ref_edge_.resize(edges_.size(), true);
  }

  void add_string(std::string& seq, int weight=1);

  void enumerate_paths(int min_weight, int max_paths, std::vector<std::pair<std::string, int> >& paths);

  static bool calc_kmer_length(std::string& ref_seq, int min_kmer, int max_kmer, int& kmer);

  bool is_source_ok();

  bool is_sink_ok();

  void prune_edges(double min_edge_freq, int min_weight);
};

class DebruijnPath {
 private:
  DebruijnPath* parent_;
  int node_id_;
  int min_weight_;

 public:
  DebruijnPath(int node_id){
    parent_     = NULL;
    min_weight_ = 1000000;
    node_id_    = node_id;
  }

  int get_min_weight() { return min_weight_; }
  int get_node_id()    { return node_id_;    } 

  DebruijnPath* add_edge(Edge* edge){
    DebruijnPath* new_path = new DebruijnPath(edge->get_destination());
    new_path->parent_      = this;
    new_path->min_weight_  = std::min(min_weight_, edge->get_weight());
    return new_path;
  }
  
  std::string get_sequence(DebruijnGraph* graph);
};



#endif
