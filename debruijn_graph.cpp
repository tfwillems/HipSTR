#include "debruijn_graph.h"

#include <math.h>
#include <algorithm>
#include <set>
#include <sstream>

bool DebruijnGraph::is_source_ok(){
  Node* source = get_node(source_kmer_);
  return (source->num_departing_edges() > 0) && (source->num_incident_edges() == 0);
}

bool DebruijnGraph::is_sink_ok(){
  Node* sink = get_node(sink_kmer_);
  return (sink->num_incident_edges() > 0) && (sink->num_departing_edges() == 0);
}

bool DebruijnGraph::calc_kmer_length(std::string& ref_seq, int min_kmer, int max_kmer, int& kmer){
  for (kmer = min_kmer; kmer <= max_kmer; kmer++){
    DebruijnGraph graph(kmer, ref_seq);
    if (!graph.has_cycles())
      return true;
  }
  return false;
}

void DebruijnGraph::add_string(std::string& seq, int weight){
  if (seq.size() <= k_)
    return;

  num_strings_++;
  std::string prev_kmer = seq.substr(0, k_);
  for (int i = 1; i < seq.size()+1-k_; i++){
    std::string next_kmer = seq.substr(i, k_);
    increment_edge(prev_kmer, next_kmer, weight);
    prev_kmer = next_kmer;
  }
}

void DebruijnGraph::prune_edges(double min_edge_freq, int min_weight){
  min_weight    = std::max(min_weight, (int)ceil(min_edge_freq*num_strings_));
  int ins_index = 0;
  std::vector<bool> keep_node(nodes_.size(), false);
  keep_node[get_node(source_kmer_)->get_id()] = true;
  keep_node[get_node(sink_kmer_)->get_id()]   = true;

  // Filter all edges with a weight below the threshold
  for (unsigned int i = 0; i < edges_.size(); i++){
    if (edges_[i]->get_weight() >= min_weight){
      if (i != ins_index)
	edges_[ins_index] = edges_[i];
      keep_node[edges_[i]->get_source()] = true;
      keep_node[edges_[i]->get_destination()] = true;
      ins_index++;
    }
    else {
      if (edges_[i]->get_source() == edges_[i]->get_destination())
	nodes_[edges_[i]->get_source()]->remove_edge(edges_[i]);
      else {
	nodes_[edges_[i]->get_source()]->remove_edge(edges_[i]);
	nodes_[edges_[i]->get_destination()]->remove_edge(edges_[i]);
      }
      delete edges_[i];
    }
  }
  edges_.resize(ins_index);

  // Filter and reindex all of the nodes with at least one edge
  std::vector<int> node_indices(nodes_.size(), -1);
  num_nodes_ = 0;
  node_map_.clear();
  for (unsigned int i = 0; i < nodes_.size(); i++){
    if (keep_node[i]){
      nodes_[i]->set_id(num_nodes_);
      node_indices[i] = num_nodes_;
      if (num_nodes_ != i){
	nodes_[num_nodes_]       = nodes_[i];
	node_labels_[num_nodes_] = node_labels_[i];
      }
      node_map_[node_labels_[num_nodes_]] = num_nodes_;
      num_nodes_++;
    }
    else
      delete nodes_[i];
  }
  nodes_.resize(num_nodes_);
  node_labels_.resize(num_nodes_);

  // Fix the node indices in each edge
  for (unsigned int i = 0; i < edges_.size(); i++){
    edges_[i]->set_source(node_indices[edges_[i]->get_source()]);
    edges_[i]->set_destination(node_indices[edges_[i]->get_destination()]);
  }

  // Add the reference path with a weight of 2
  add_string(ref_seq_, 2);
}

bool path_comparator(DebruijnPath* p1, DebruijnPath* p2){
  return p1->get_min_weight() < p2->get_min_weight();
}

/*
 * Generate all kmers that differ from the input kmer by a 1 bp mismatch
 * Add them to the list of nodes if they're present in the graph
 * and they satisfy the source/sink requirements
 */
void DebruijnGraph::get_alt_kmer_nodes(std::string& kmer, bool source, bool sink, std::vector<Node*>& nodes){
  assert(nodes.empty());
  std::vector<char> bases(4, 'N'); bases[0]='A'; bases[1]='C'; bases[2]='G'; bases[3]='T';
  for (unsigned int i = 0; i < kmer.size(); ++i){
    char orig = kmer[i];
    for (unsigned int j = 0; j < 4; ++j){
      if (bases[j] != orig){
	kmer[i] = bases[j];
	if (has_node(kmer)){
	  Node* node = get_node(kmer);
	  if (source && node->num_incident_edges() > 0)
	    continue;
	  if (sink && node->num_departing_edges() > 0)
	    continue;
	  nodes.push_back(node);
	}
      }
    }
    kmer[i] = orig;
  }
}

void DebruijnGraph::enumerate_paths(int min_weight, int max_paths, std::vector<std::pair<std::string, int> >& paths){
  assert(paths.empty());

  // Create a heap containing the source node
  Node* source = get_node(source_kmer_);
  Node* sink   = get_node(sink_kmer_);
  int sink_id  = sink->get_id();
  std::vector<DebruijnPath*> all_paths(1, new DebruijnPath(source->get_id()));
  std::vector<DebruijnPath*> heap(1, all_paths.back());
  std::make_heap(heap.begin(), heap.end(), path_comparator);

  // Add all kmers that differ by a 1 bp mismatch from the source kmer to the heap
  std::vector<Node*> alt_source_nodes;
  get_alt_kmer_nodes(source_kmer_, true, false, alt_source_nodes);
  for (unsigned int i = 0; i < alt_source_nodes.size(); i++){
    all_paths.push_back(new DebruijnPath(alt_source_nodes[i]->get_id()));
    heap.push_back(all_paths.back());
    std::push_heap(heap.begin(), heap.end(), path_comparator);
  }

  // Construct a set of sink nodes based on the sink kmer and all of its 1bp mismatches
  std::vector<Node*> alt_sink_nodes;
  std::set<int> sink_ids; sink_ids.insert(sink_id);
  get_alt_kmer_nodes(sink_kmer_, false, true, alt_sink_nodes);
  for (unsigned int i = 0; i < alt_sink_nodes.size(); i++)
    sink_ids.insert(alt_sink_nodes[i]->get_id());

  while (!heap.empty()){
    if (paths.size() == max_paths)
      break;

    std::pop_heap(heap.begin(), heap.end(), path_comparator);
    DebruijnPath* best = heap.back(); heap.pop_back();

    // If we reached a sink, record the weight and sequence of the path
    if (sink_ids.find(best->get_node_id()) != sink_ids.end())
      paths.push_back(std::pair<std::string, int>(best->get_sequence(this), best->get_min_weight()));

    std::vector<Edge*>& edges = nodes_[best->get_node_id()]->get_departing_edges();
    for (unsigned int i = 0; i < edges.size(); i++){
      if (edges[i]->get_weight() < min_weight)
	continue;
      all_paths.push_back(best->add_edge(edges[i]));
      heap.push_back(all_paths.back());
      std::push_heap(heap.begin(), heap.end(), path_comparator);
    }
  }

  for (unsigned int i = 0; i < all_paths.size(); i++)
    delete all_paths[i];
}


std::string DebruijnPath::get_sequence(DebruijnGraph* graph){
  std::stringstream ss;
  const std::string& kmer = graph->get_node_label(node_id_);
  for (auto base_iter = kmer.rbegin(); base_iter != kmer.rend(); base_iter++)
    ss << *base_iter;

  DebruijnPath* parent = parent_;
  while (parent != NULL){
    ss << graph->get_node_label(parent->node_id_).front();
    parent = parent->parent_;
  }
  std::string result = ss.str();
  std::reverse(result.begin(), result.end());
  return result;
}
