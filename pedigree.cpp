#include <assert.h>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>

#include "pedigree.h"
#include "error.h"

void PedigreeGraph::init_no_ancestors() {
  no_ancestors_.clear();
  for (int i = 0; i < nodes_.size(); i++)
    if (!nodes_[i]->has_mother() && !nodes_[i]->has_father())
      no_ancestors_.push_back(nodes_[i]);
}

void PedigreeGraph::init_no_descendants() {
  no_descendants_.clear();
  for (int i = 0; i < nodes_.size(); i++)
    if (nodes_[i]->get_children().size() == 0)
      no_descendants_.push_back(nodes_[i]);
}

bool PedigreeGraph::topological_sort(std::vector<PedigreeNode*>& nodes){
  no_ancestors_.clear();
  no_descendants_.clear();
  nodes_.clear();
  
  std::map<PedigreeNode*, int> parent_counts;
  std::vector<PedigreeNode*>   sources;
  for (int i = 0; i < nodes.size(); i++){
    int count = nodes[i]->has_mother() + nodes[i]->has_father();
    if (count == 0)
      sources.push_back(nodes[i]);
    else
      parent_counts[nodes[i]] = count;
  }

  while (sources.size() != 0){
    PedigreeNode* source = sources.back();
    std::vector<PedigreeNode*>& children = source->get_children();
    nodes_.push_back(source);
    sources.pop_back();

    for (auto child_iter = children.begin(); child_iter != children.end(); child_iter++){
      auto count_iter = parent_counts.find(*child_iter);
      if (count_iter == parent_counts.end()){
	source->print(std::cerr);
	(*child_iter)->print(std::cerr);
	printErrorAndDie("Logical error in topological_sort() for parent " + source->get_name() + " and child " + (*child_iter)->get_name());
      }
      else if (count_iter->second == 1){
	sources.push_back(*child_iter);
	parent_counts.erase(count_iter);
      }
      else
	count_iter->second -= 1;
    }
  }
  return parent_counts.size() == 0; // Only a DAG if no unprocessed individuals are left
}

bool PedigreeGraph::build(std::string filename) {
  std::ifstream input(filename.c_str());
  if (!input.is_open())
    printErrorAndDie("Failed to open pedigree file " + filename);
  
  std::map<std::string, PedigreeNode*> samples;
  std::vector<PedigreeNode*> nodes;
  std::string line;
  //std::getline(input, line));  // TO DO: Fix weird expection that occurs if I try to skip header  
  while (std::getline(input, line)){
    std::istringstream iss(line);
    std::string family, child, father, mother;
    if(! (iss >> family >> child >> father >> mother))
      printErrorAndDie("Improperly formated .ped pedigree file " + filename);

    if (child.compare("0") == 0)
      printErrorAndDie("Invalid individual id " + child);

    // Create new nodes for any previously unseen samples that have an identifier other than 0
    if (samples.find(child) == samples.end()){
      PedigreeNode* new_node = new PedigreeNode(child);
      nodes.push_back(new_node);
      samples[child] = new_node;
    }
    if (mother.compare("0") != 0 && samples.find(mother) == samples.end()){
      PedigreeNode* new_node = new PedigreeNode(mother);
      samples[mother] = new_node;
      nodes.push_back(new_node);
    }
    if (father.compare("0") != 0 && samples.find(father) == samples.end()){
      PedigreeNode* new_node = new PedigreeNode(father);
      samples[father] = new_node;
      nodes.push_back(new_node);
    }
    
    // Store relationships in node instance
    PedigreeNode* child_node  = samples.find(child)->second;
    PedigreeNode* mother_node = (mother.compare("0") == 0 ? NULL : samples.find(mother)->second);
    PedigreeNode* father_node = (father.compare("0") == 0 ? NULL : samples.find(father)->second);
    child_node->set_mother(mother_node);
    child_node->set_father(father_node);
    if (mother_node != NULL) mother_node->add_child(child_node);
    if (father_node != NULL) father_node->add_child(child_node);
  }
  input.close();
  
  // Sort nodes in pedigree graph topologically
  return topological_sort(nodes);
}

void PedigreeGraph::prune(std::set<std::string>& sample_set){
  // Determine if each node has an upstream requested sample
  std::map<PedigreeNode*, bool> upstream_status;
  for (int i = 0; i < nodes_.size(); i++){
    bool has_upstream = sample_set.find(nodes_[i]->get_name()) != sample_set.end();
    has_upstream     |= (nodes_[i]->has_father() && upstream_status[nodes_[i]->get_father()]);
    has_upstream     |= (nodes_[i]->has_mother() && upstream_status[nodes_[i]->get_mother()]);
    upstream_status[nodes_[i]] = has_upstream;
  }
  
  // Determine if each node has a downstream requested sample
  std::map<PedigreeNode*, bool> downstream_status;
  for (int i = nodes_.size()-1; i >= 0; i--) {
    bool has_downstream = sample_set.find(nodes_[i]->get_name()) != sample_set.end();
    for(auto iter = nodes_[i]->get_children().begin(); iter != nodes_[i]->get_children().end(); iter++)
      has_downstream |= downstream_status[(*iter)];
    downstream_status[nodes_[i]] = has_downstream;
  }

  // Determine if nodes have a requested sample both above and below
  // If not, mark them for removal
  std::map<PedigreeNode*, bool> removal_status;
  for (int i = 0; i < nodes_.size(); i++)
    removal_status[nodes_[i]] = (!upstream_status[nodes_[i]] || !downstream_status[nodes_[i]]);
      
  // Remove and modify nodes member data accordingly
  int insert_index = 0;
  for (int i = 0; i < nodes_.size(); i++){
    if (removal_status[nodes_[i]])
      delete nodes_[i];
    else {
      if (nodes_[i]->has_father() && removal_status[nodes_[i]->get_father()])
	nodes_[i]->set_father(NULL);
      if (nodes_[i]->has_mother() && removal_status[nodes_[i]->get_mother()])
	nodes_[i]->set_mother(NULL);
      int child_ins_index = 0;
      std::vector<PedigreeNode*>& children = nodes_[i]->get_children();;
      for (int j = 0; j < children.size(); j++)
	if (!removal_status[children[j]])
	  children[child_ins_index++] = children[j];
      children.resize(child_ins_index);
      nodes_[insert_index++] = nodes_[i];
    }
  }
  nodes_.resize(insert_index);

  // Rebuild internal structures
  no_ancestors_.clear();
  no_descendants_.clear();
  init_no_ancestors();
  init_no_descendants();
}

bool PedigreeGraph::build_subgraph(std::vector<PedigreeNode*>& subgraph_nodes){
  std::map<std::string, PedigreeNode*> samples;
  std::vector<PedigreeNode*> nodes;
  for (auto node_iter = subgraph_nodes.begin(); node_iter != subgraph_nodes.end(); node_iter++){
    PedigreeNode* child_node;
    PedigreeNode* mother_node = NULL;
    PedigreeNode* father_node = NULL;
    std::string child = (*node_iter)->get_name();

    // Create new nodes for any previously unseen samples
    if (samples.find(child) == samples.end()){
      child_node = new PedigreeNode(child);
      samples[child] = child_node;
      nodes.push_back(child_node);
    }
    else
      child_node = samples[child];

    if ((*node_iter)->has_mother()){
      std::string mother = (*node_iter)->get_mother()->get_name();
      if (samples.find(mother) == samples.end()){
	mother_node = new PedigreeNode(mother);
	nodes.push_back(mother_node);
	samples[mother] = mother_node;
      }
      else
	mother_node = samples[mother];
    }

    if ((*node_iter)->has_father()){
      std::string father = (*node_iter)->get_father()->get_name();
      if (samples.find(father) == samples.end()){
	father_node = new PedigreeNode(father);
	nodes.push_back(father_node);
	samples[father] = father_node;
      }
      else
	father_node = samples[father];
    }

    // Store relationships in node instance
    child_node->set_mother(mother_node);
    child_node->set_father(father_node);
    if (mother_node != NULL) mother_node->add_child(child_node);
    if (father_node != NULL) father_node->add_child(child_node);
  }
  return topological_sort(nodes);
}

void PedigreeGraph::split_into_connected_components(std::vector<PedigreeGraph>& components){
  // Determine the component to which each node belongs
  std::set<PedigreeNode*> visited;
  std::vector< std::vector<PedigreeNode*> > component_nodes;
  for (auto node_iter = nodes_.begin(); node_iter != nodes_.end(); node_iter++){
    if (visited.find(*node_iter) != visited.end())
      continue;

    component_nodes.push_back(std::vector<PedigreeNode*>());
    std::vector<PedigreeNode*> to_process(1, *node_iter);
    while (to_process.size() != 0){
      PedigreeNode* seed = to_process.back();
      to_process.pop_back();
      if (visited.find(seed) != visited.end())
	continue;

      visited.insert(seed);
      component_nodes.back().push_back(seed);

      if (seed->has_mother() && (visited.find(seed->get_mother()) == visited.end()))
	to_process.push_back(seed->get_mother());
      if (seed->has_father() && (visited.find(seed->get_father()) == visited.end()))
	to_process.push_back(seed->get_father());

      std::vector<PedigreeNode*>& children = seed->get_children();
      for (auto child_iter = children.begin(); child_iter != children.end(); child_iter++)
	if (visited.find(*child_iter) == visited.end())
	  to_process.push_back(*child_iter);
    }
  }
  assert(visited.size() == nodes_.size());

  // Construct individual graphs for each component
  assert(components.size() == 0);
  for (int i = 0; i < component_nodes.size(); i++)
    components.push_back(PedigreeGraph(component_nodes[i]));
}

bool PedigreeGraph::is_nuclear_family(){
  if (no_ancestors_.size() != 2)
    return false;
  if (no_descendants_.size() == 0)
    return false;
  if (no_ancestors_.size() + no_descendants_.size() != nodes_.size())
    return false;

  std::string p1 = no_ancestors_[0]->get_name();
  std::string p2 = no_ancestors_[1]->get_name();
  for (auto node_iter = no_descendants_.begin(); node_iter != no_descendants_.end(); node_iter++){
    if ((!(*node_iter)->has_mother()) || (!(*node_iter)->has_father()))
      return false;

    std::string mother = (*node_iter)->get_mother()->get_name();
    std::string father = (*node_iter)->get_father()->get_name();
    if ((mother.compare(p1) != 0) || (father.compare(p2) != 0))
      if ((mother.compare(p2) != 0) || (father.compare(p1) != 0))
	return false;
  }
  return true;
}

NuclearFamily PedigreeGraph::convert_to_nuclear_family(){
  assert(is_nuclear_family());
  std::string mother = no_descendants_[0]->get_mother()->get_name();;
  std::string father = no_descendants_[0]->get_father()->get_name();;
  std::vector<std::string> children;
  for (auto node_iter = no_descendants_.begin(); node_iter != no_descendants_.end(); node_iter++)
    children.push_back((*node_iter)->get_name());
  return NuclearFamily(mother, father, children);
}

void read_sample_list(std::string input_file, std::set<std::string>& sample_set){
  sample_set.clear();
  std::ifstream input(input_file);
  if (!input.is_open())
    printErrorAndDie("Unable to open sample list file " + input_file);
  std::string line;
  while (std::getline(input, line))
    sample_set.insert(line);
}

void extract_pedigree_nuclear_families(std::string pedigree_fam_file, std::set<std::string>& samples_with_data,
				       std::vector<NuclearFamily>& nuclear_families, std::ostream& logger){
  assert(nuclear_families.size() == 0);

  // Read the original pedigree
  PedigreeGraph pedigree(pedigree_fam_file);

  // Remove irrelevant samples from pedigree
  pedigree.prune(samples_with_data);

  // Identify simple nuclear families in the pedigree
  std::vector<PedigreeGraph> pedigree_components;
  pedigree.split_into_connected_components(pedigree_components);
  int num_others = 0;
  for (unsigned int i = 0; i < pedigree_components.size(); i++){
    if (pedigree_components[i].is_nuclear_family())
      nuclear_families.push_back(pedigree_components[i].convert_to_nuclear_family());
    else
      num_others++;
  }
  logger << "Detected " << nuclear_families.size() << " nuclear families and " << num_others << " other family structures\n";
}
