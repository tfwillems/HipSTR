#include "error.h"
#include "directed_graph.h"

void DirectedGraph::increment_edge(const std::string& val_1, const std::string& val_2, int delta){
  Node* source  = get_node(val_1);
  Node* dest    = get_node(val_2);
  int source_id = source->get_id();

  std::vector<Edge*>& edges = dest->get_incident_edges();
  for (unsigned int i = 0; i < edges.size(); i++){
    if (edges[i]->get_source() == source_id){
      edges[i]->inc_weight(delta);
      return;
    }
  }

  Edge* new_edge = new Edge(edges_.size(), source_id, dest->get_id(), delta);
  if (source_id == dest->get_id())
    source->add_edge(new_edge);
  else {
    source->add_edge(new_edge);
    dest->add_edge(new_edge);
  }
  edges_.push_back(new_edge);
}

bool DirectedGraph::can_sort_topologically() const {
  std::map<int, int> parent_counts;
  std::vector<int> sources;
  for (unsigned int i = 0; i < nodes_.size(); i++){
    int count = nodes_[i]->num_incident_edges();
    if (count == 0)
      sources.push_back(i);
    else
      parent_counts[i] = count;
  }

  std::vector<int> ordered_nodes;
  std::vector<int> children;
  while (sources.size() != 0){
    int source = sources.back();
    nodes_[source]->get_child_nodes(children);
    ordered_nodes.push_back(source);
    sources.pop_back();

    for (auto child_iter = children.begin(); child_iter != children.end(); child_iter++){
      auto count_iter = parent_counts.find(*child_iter);
      if (count_iter == parent_counts.end())
 	printErrorAndDie("Logical error in topological_sort()");
      else if (count_iter->second == 1){
	sources.push_back(*child_iter);
	parent_counts.erase(count_iter);
      }
      else
	count_iter->second -= 1;
    }
    children.clear();
  }
 
  if (parent_counts.size() == 0)
    return true;
  else
    return false; // Only a DAG if no unprocessed individuals are left
}

void DirectedGraph::print(std::ostream& out) const {
  assert(nodes_.size() == node_labels_.size());

  out << "NODES" << "\n";
  for (unsigned int i = 0; i < nodes_.size(); i++)
    out << "\t" << i << "\t" << node_labels_[i] << "\n";
  out << "\n";

  out << "EDGES " << edges_.size() <<  "\n";
  for (unsigned int i = 0; i < edges_.size(); i++)
    out << "\t" << i << "\t" << edges_[i]->get_source() << "\t" << edges_[i]->get_destination()
	<< "\t" << edges_[i]->get_weight() << "\n";
  out << "\n" << std::endl;
}
