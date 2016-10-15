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
  std::string source_kmer;
  std::string sink_kmer;

 public:
  DebruijnGraph(int k, std::string& ref_seq){
    assert(ref_seq.size() > k);
    k_          = k;
    source_kmer = ref_seq.substr(0, k_);
    sink_kmer   = ref_seq.substr(ref_seq.size()-k, k_);

    // Add it twice so that the reference path has a weight of at least 2
    add_string(ref_seq); add_string(ref_seq);
  }

  void add_string(std::string& seq){
    if (seq.size() <= k_)
      return;

    std::string prev_kmer = seq.substr(0, k_);
    for (int i = 1; i < seq.size()+1-k_; i++){
      std::string next_kmer = seq.substr(i, k_);
      increment_edge(prev_kmer, next_kmer);
      prev_kmer = next_kmer;
    }
  }

  void enumerate_paths(int min_weight, int max_paths, std::vector<std::pair<std::string, int> >& paths);
};

class DebruijnPath {
 private:
  DebruijnPath* parent_;
  int node_id_;
  int min_weight_;

 public:
  DebruijnPath(Node* node){
    parent_     = NULL;
    min_weight_ = 1000000;
    node_id_    = node->get_id();
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
