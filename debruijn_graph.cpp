#include "debruijn_graph.h"

#include <sstream>

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

bool path_comparator(DebruijnPath* p1, DebruijnPath* p2){
  return p1->get_min_weight() < p2->get_min_weight();
}

void DebruijnGraph::enumerate_paths(int min_weight, int max_paths, std::vector<std::pair<std::string, int> >& paths){
  assert(paths.empty());

  // Create a heap containing the source node
  Node* source = get_node(source_kmer);
  Node* sink   = get_node(sink_kmer);
  std::vector<DebruijnPath*> all_paths(1, new DebruijnPath(source));
  std::vector<DebruijnPath*> heap(1, all_paths.back());
  std::make_heap(heap.begin(), heap.end(), path_comparator);

  while (!heap.empty()){
    if (paths.size() == max_paths)
      break;

    std::pop_heap(heap.begin(), heap.end(), path_comparator);
    DebruijnPath* best = heap.back(); heap.pop_back();

    // If we reached the sink, record the weight and sequence of the path
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
