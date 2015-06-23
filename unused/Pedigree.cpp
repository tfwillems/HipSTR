#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>

#include "Pedigree.h"

#include "error.h"

void PEDIGREE_GRAPH::init_no_ancestors() {
  no_ancestors_.clear();
  for (int i = 0; i < nodes_.size(); i++)
    if (!nodes_[i]->has_mother() && !nodes_[i]->has_father())
      no_ancestors_.push_back(nodes_[i]);
}

void PEDIGREE_GRAPH::init_no_descendants() {
  no_descendants_.clear();
  for (int i = 0; i < nodes_.size(); i++)
    if (nodes_[i]->get_children().size() == 0)
      no_descendants_.push_back(nodes_[i]);
}

bool PEDIGREE_GRAPH::topological_sort(std::vector<PEDIGREE_NODE*>& nodes){
  no_ancestors_.clear();
  no_descendants_.clear();
  nodes_.clear();
  
  std::map<PEDIGREE_NODE*, int> parent_counts;
  std::vector<PEDIGREE_NODE*>   sources;
  for (int i = 0; i < nodes.size(); i++){
    int count = nodes[i]->has_mother() + nodes[i]->has_father();
    if (count == 0)
      sources.push_back(nodes[i]);
    else
      parent_counts[nodes[i]] = count;
  }

  while (sources.size() != 0){
    PEDIGREE_NODE* source = sources.back();
    std::vector<PEDIGREE_NODE*>& children = source->get_children();
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

bool PEDIGREE_GRAPH::build(std::string filename) {
  std::ifstream input(filename.c_str());
  if (!input.is_open())
    printErrorAndDie("Failed to open pedigree file " + filename);
  
  std::map<std::string, PEDIGREE_NODE*> samples;
  std::vector<PEDIGREE_NODE*> nodes;
  std::string line;
  //std::getline(input, line));  // TO DO: Fix weird expection that occurs if I try to skip header  
  while (std::getline(input, line)){
    std::istringstream iss(line);
    std::string family, child, father, mother;
    if(! (iss >> family >> child >> father >> mother))
      printErrorAndDie("Improperly formated .ped pedigree file " + filename);

    if (child.compare("0") == 0)
      printErrorAndDie("Invalid individual id " + child);

    // Create new nodes for any previously unseen samples that have 
    // an identifier other than 0
    if (samples.find(child) == samples.end()){
      PEDIGREE_NODE* new_node = new PEDIGREE_NODE(child);
      nodes.push_back(new_node);
      samples[child] = new_node;
    }
    if (mother.compare("0") != 0 && samples.find(mother) == samples.end()){
      PEDIGREE_NODE* new_node = new PEDIGREE_NODE(mother);
      samples[father] = new_node;
      nodes.push_back(new_node);
    }
    if (father.compare("0") != 0 && samples.find(father) == samples.end()){
      PEDIGREE_NODE* new_node = new PEDIGREE_NODE(father);
      samples[mother] = new_node;
      nodes.push_back(new_node);
    }
    
    // Store relationships in node instance
    PEDIGREE_NODE* child_node  = samples.find(child)->second;
    PEDIGREE_NODE* mother_node = (mother.compare("0") == 0 ? NULL : samples.find(mother)->second);
    PEDIGREE_NODE* father_node = (father.compare("0") == 0 ? NULL : samples.find(father)->second); 
    child_node->set_mother(mother_node);
    child_node->set_father(father_node);
    if (mother_node != NULL) mother_node->add_child(child_node);
    if (father_node != NULL) father_node->add_child(child_node);
  }
  input.close();
  
  // Sort nodes in pedigree graph topologically
  return topological_sort(nodes);
}

void PEDIGREE_GRAPH::prune(std::set<std::string>& sample_set){
  // Determine if each node has an upstream requested sample
  std::map<PEDIGREE_NODE*, bool> upstream_status;
  for (int i = 0; i < nodes_.size(); i++){
    bool has_upstream = sample_set.find(nodes_[i]->get_name()) != sample_set.end();
    has_upstream     |= (nodes_[i]->has_father() && upstream_status[nodes_[i]->get_father()]);
    has_upstream     |= (nodes_[i]->has_mother() && upstream_status[nodes_[i]->get_mother()]);
    upstream_status[nodes_[i]] = has_upstream;
  }
  
  // Determine if each node has a downstream requested sample
  std::map<PEDIGREE_NODE*, bool> downstream_status;
  for (int i = nodes_.size()-1; i >= 0; i--) {
    bool has_downstream = sample_set.find(nodes_[i]->get_name()) != sample_set.end();
    for(auto iter = nodes_[i]->get_children().begin(); iter != nodes_[i]->get_children().end(); iter++)
      has_downstream |= downstream_status[(*iter)];
    downstream_status[nodes_[i]] = has_downstream;
  }

  // Determine if nodes have a requested sample both above and below
  // If not, mark them for removal
  std::map<PEDIGREE_NODE*, bool> removal_status;
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
      std::vector<PEDIGREE_NODE*>& children = nodes_[i]->get_children();;
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


void read_sample_list(std::string input_file, std::set<std::string>& sample_set){
  sample_set.clear();
  std::ifstream input(input_file);
  if (!input.is_open())
    printErrorAndDie("Unable to open sample list file " + input_file);
  std::string line;
  while (std::getline(input, line))
    sample_set.insert(line);
}
