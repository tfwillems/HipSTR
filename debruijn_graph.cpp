#include "debruijn_graph.h"

#include <set>
#include <sstream>

bool DebruijnGraph::is_source_ok(){
  Node* source = get_node(source_kmer);
  assert(source->num_departing_edges() > 0);
  return source->num_incident_edges() == 0;
}

bool DebruijnGraph::is_sink_ok(){
  Node* sink = get_node(sink_kmer);
  assert(sink->num_incident_edges() > 0);
  return sink->num_departing_edges() == 0;
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

  std::string prev_kmer = seq.substr(0, k_);
  for (int i = 1; i < seq.size()+1-k_; i++){
    std::string next_kmer = seq.substr(i, k_);
    increment_edge(prev_kmer, next_kmer, weight);
    prev_kmer = next_kmer;
  }
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
  Node* source = get_node(source_kmer);
  Node* sink   = get_node(sink_kmer);
  int sink_id  = sink->get_id();
  std::vector<DebruijnPath*> all_paths(1, new DebruijnPath(source->get_id()));
  std::vector<DebruijnPath*> heap(1, all_paths.back());
  std::make_heap(heap.begin(), heap.end(), path_comparator);

  // Add all kmers that differ by a 1 bp mismatch from the source kmer to the heap
  std::vector<Node*> alt_source_nodes;
  get_alt_kmer_nodes(source_kmer, true, false, alt_source_nodes);
  for (unsigned int i = 0; i < alt_source_nodes.size(); i++){
    all_paths.push_back(new DebruijnPath(alt_source_nodes[i]->get_id()));
    heap.push_back(all_paths.back());
    std::push_heap(heap.begin(), heap.end(), path_comparator);
  }

  // Construct a set of sink nodes based on the sink kmer and all of its 1bp mismatches
  std::vector<Node*> alt_sink_nodes;
  std::set<int> sink_ids; sink_ids.insert(sink_id);
  get_alt_kmer_nodes(sink_kmer, false, true, alt_sink_nodes);
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
