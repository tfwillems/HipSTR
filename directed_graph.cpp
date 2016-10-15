#include "error.h"
#include "directed_graph.h"

bool DirectedGraph::can_sort_topologically(){
  std::map<Node*, int> parent_counts;
  std::vector<Node*> sources;
  for (unsigned int i = 0; i < nodes_.size(); i++){
    int count = nodes_[i]->num_incident_edges();
    if (count == 0)
      sources.push_back(nodes_[i]);
    else
      parent_counts[nodes_[i]] = count;
  }

  std::vector<Node*> ordered_nodes;
  std::vector<Node*> children;
  while (sources.size() != 0){
    Node* source = sources.back();
    source->get_child_nodes(children);
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
