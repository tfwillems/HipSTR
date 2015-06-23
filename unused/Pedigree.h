#ifndef PEDIGREE_H_
#define PEDIGREE_H_

#include <algorithm>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "error.h"

class PEDIGREE_NODE {
 private:
  std::string name_;
  PEDIGREE_NODE* mother_;
  PEDIGREE_NODE* father_;
  std::vector<PEDIGREE_NODE*> children_;


 public:
  PEDIGREE_NODE(std::string name){
    name_     = name;
    mother_   = NULL;
    father_   = NULL;
    children_ = std::vector<PEDIGREE_NODE*>();
  }

  PEDIGREE_NODE* get_mother() { return mother_; }
  PEDIGREE_NODE* get_father() { return father_; }
  std::string    get_name()   { return name_;   }
  std::vector<PEDIGREE_NODE*>& get_children() { return children_; }

  bool has_mother()      { return mother_ != NULL; }
  bool has_father()      { return father_ != NULL; }

  void set_mother(PEDIGREE_NODE* mother) { mother_ = mother;           }
  void set_father(PEDIGREE_NODE* father) { father_ = father;           }
  void add_child (PEDIGREE_NODE* child)  { children_.push_back(child); }
  void del_child (PEDIGREE_NODE* child)  {
    auto iter = std::find(children_.begin(), children_.end(), child);;
    if (iter == children_.end())
      printErrorAndDie("Can't delete chid from node as it is not among the existing children");
    children_.erase(iter);
  }
  
  void print(std::ostream& out){
    out << "NAME: "    << name_ 
	<< "\tFATHER:" << (father_ == NULL ? "NONE": father_->get_name()) 
	<< "\tMOTHER:" << (mother_ == NULL ? "NONE" : mother_->get_name()) << std::endl; 
  }
};

class PEDIGREE_GRAPH {
 private:
  // Nodes that don't have any ancestors 
  std::vector<PEDIGREE_NODE*> no_ancestors_;
  
  // Nodes that don't have any descendants
  std::vector<PEDIGREE_NODE*> no_descendants_;

  // Nodes in graph sorted in topological order
  std::vector<PEDIGREE_NODE*> nodes_;

  bool topological_sort(std::vector<PEDIGREE_NODE*>& nodes);
  bool build(std::string input_file);  
  void init_no_ancestors();
  void init_no_descendants();
    
 public:
  PEDIGREE_GRAPH(std::string input_file){
    no_ancestors_   = std::vector<PEDIGREE_NODE*>();
    no_descendants_ = std::vector<PEDIGREE_NODE*>();
    nodes_          = std::vector<PEDIGREE_NODE*>();
    bool success    = build(input_file);
    if (!success)
      printErrorAndDie("Supplied pedigree file " + input_file + " contains cycles");
    init_no_ancestors();
    init_no_descendants();
  }

  ~PEDIGREE_GRAPH(){
    for (int i = 0; i < nodes_.size(); i++)
      delete nodes_[i];
    nodes_.clear();
  }

  void print(std::ostream& out){
    out << "Pedigree graph contains " << nodes_.size() << " nodes" << std::endl;
  }

  void prune(std::set<std::string>& sample_set);
};


void read_sample_list(std::string input_file, std::set<std::string>& sample_set);

#endif
