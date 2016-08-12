#ifndef PEDIGREE_H_
#define PEDIGREE_H_

#include <algorithm>
#include <assert.h>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "../error.h"
#include "../vcflib/src/Variant.h"

class NuclearFamily {
 private:
  std::string mother_, father_;
  std::vector<std::string> children_;

 public:
  NuclearFamily(std::string mother, std::string father, std::vector<std::string> children){
    mother_   = mother;
    father_   = father;
    children_ = children;
  }

  const std::string& get_mother() { return mother_; }
  const std::string& get_father() { return father_; }
  const std::vector<std::string>& get_children() { return children_; }
  const int size()         {  return 2 + children_.size(); }
  const int num_children() { return children_.size();      }

  bool is_missing_genotype(vcflib::Variant& variant){
    if (variant.getGenotype(mother_).empty() || variant.getGenotype(father_).empty())
      return true;
    for (auto child_iter = children_.begin(); child_iter != children_.end(); child_iter++)
      if (variant.getGenotype(*child_iter).empty())
	return true;
    return false;
  }

  bool is_mendelian(vcflib::Variant& variant){
    std::string father_gt = variant.getGenotype(father_);
    std::string mother_gt = variant.getGenotype(mother_);
    if (father_gt.empty() || mother_gt.empty())
      return false;
    assert(father_gt.size() == 3 && mother_gt.size() == 3);
    int f_1 = father_gt[0]-'0', f_2 = father_gt[2]-'0';
    int m_1 = mother_gt[0]-'0', m_2 = mother_gt[2]-'0';
    for (auto child_iter = children_.begin(); child_iter != children_.end(); child_iter++){
      std::string child_gt = variant.getGenotype(*child_iter);
      if (child_gt.empty())
	return false;
      assert(child_gt.size() == 3);
      int c_1 = child_gt[0]-'0', c_2 = child_gt[2]-'0';
      if ((c_1 != m_1 && c_1 != m_2) || (c_2 != f_1 && c_2 != f_2))
	if ((c_1 != f_1 && c_1 != f_2) || (c_2 != m_1 && c_2 != m_2))
	  return false;
    }
    return true;
  }
};

class PedigreeNode {
 private:
  std::string name_;
  PedigreeNode* mother_;
  PedigreeNode* father_;
  std::vector<PedigreeNode*> children_;

 public:
  PedigreeNode(std::string name){
    name_     = name;
    mother_   = NULL;
    father_   = NULL;
    children_ = std::vector<PedigreeNode*>();
  }

  ~PedigreeNode(){
    children_.clear();
  }

  bool has_mother() const    { return mother_ != NULL; }
  bool has_father() const    { return father_ != NULL; }
  PedigreeNode* get_mother() const { return mother_; }
  PedigreeNode* get_father() const { return father_; }
  std::string   get_name()   const { return name_;   }
  std::vector<PedigreeNode*>& get_children() { return children_; }

  void set_mother(PedigreeNode* mother) { mother_ = mother;           }
  void set_father(PedigreeNode* father) { father_ = father;           }
  void add_child (PedigreeNode* child)  { children_.push_back(child); }
  void del_child (PedigreeNode* child)  {
    auto iter = std::find(children_.begin(), children_.end(), child);;
    if (iter == children_.end())
      printErrorAndDie("Can't delete chid from node as it is not among the existing children");
    children_.erase(iter);
  }
  
  void print(std::ostream& out){
    out << "NAME:"     << name_
	<< "\tFATHER:" << (father_ == NULL ? "NONE": father_->get_name()) 
	<< "\tMOTHER:" << (mother_ == NULL ? "NONE" : mother_->get_name()) << std::endl; 
  }
};

class PedigreeGraph {
 private:
  // Nodes that don't have any ancestors 
  std::vector<PedigreeNode*> no_ancestors_;
  
  // Nodes that don't have any descendants
  std::vector<PedigreeNode*> no_descendants_;

  // Nodes in graph sorted in topological order
  std::vector<PedigreeNode*> nodes_;

  bool topological_sort(std::vector<PedigreeNode*>& nodes);
  bool build(std::string input_file);  
  void init_no_ancestors();
  void init_no_descendants();
  bool build_subgraph(std::vector<PedigreeNode*>& sorted_nodes);
    
 public:
  PedigreeGraph(){
    no_ancestors_   = std::vector<PedigreeNode*>();
    no_descendants_ = std::vector<PedigreeNode*>();
    nodes_          = std::vector<PedigreeNode*>();
  }

  PedigreeGraph(std::string input_file){
    no_ancestors_   = std::vector<PedigreeNode*>();
    no_descendants_ = std::vector<PedigreeNode*>();
    nodes_          = std::vector<PedigreeNode*>();
    bool success    = build(input_file);
    if (!success)
      printErrorAndDie("Supplied pedigree file " + input_file + " contains cycles");
    init_no_ancestors();
    init_no_descendants();
  }

  PedigreeGraph(std::vector<PedigreeNode*>& subgraph_nodes){
    if (!build_subgraph(subgraph_nodes))
      printErrorAndDie("Subgraph in pedigree contains a cycle");
    init_no_ancestors();
    init_no_descendants();
  }

  PedigreeGraph(const PedigreeGraph& other){
    // Create new nodes with identical names
    std::map<std::string, int> indices;
    for (int i = 0; i < other.nodes_.size(); i++){
      indices[other.nodes_[i]->get_name()] = i;
      nodes_.push_back(new PedigreeNode(other.nodes_[i]->get_name()));
    }

    // Restore links between cloned nodes
    for (int i = 0; i < nodes_.size(); i++){
      if (other.nodes_[i]->has_mother()){
	std::string mother = other.nodes_[i]->get_mother()->get_name();
	nodes_[i]->set_mother(nodes_[indices[mother]]);
      }

      if (other.nodes_[i]->has_father()){
	std::string father = other.nodes_[i]->get_father()->get_name();
	nodes_[i]->set_father(nodes_[indices[father]]);
      }

      const std::vector<PedigreeNode*> children = other.nodes_[i]->get_children();
      for (auto child_iter = children.begin(); child_iter != children.end(); child_iter++){
	std::string child = (*child_iter)->get_name();
	nodes_[i]->add_child(nodes_[indices[child]]);
      }
    }

    init_no_ancestors();
    init_no_descendants();
  }

  ~PedigreeGraph(){
    for (int i = 0; i < nodes_.size(); i++)
      delete nodes_[i];
    nodes_.clear();
    no_ancestors_.clear();
    no_descendants_.clear();
  }

  int size(){ return nodes_.size(); }

  void print(std::ostream& out){
    out << "Pedigree graph contains " << nodes_.size() << " nodes" << std::endl;
  }

  void prune(std::set<std::string>& sample_set);

  void split_into_connected_components(std::vector<PedigreeGraph>& components);

  bool is_nuclear_family();

  NuclearFamily convert_to_nuclear_family();
};

void read_sample_list(std::string input_file, std::set<std::string>& sample_set);

#endif
