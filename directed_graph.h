#ifndef DIRECTED_GRAPH_H_
#define DIRECTED_GRAPH_H_

#include <assert.h>
#include <iostream>
#include <map>
#include <vector>

class Node;

class Edge {
  int source_;
  int destination_;
  int weight_;

 public:
  Edge(int source, int destination, int weight){
    source_      = source;
    destination_ = destination;
    weight_      = weight;
  }

  int get_source()           { return source_;      }
  int get_destination()      { return destination_; }
  void inc_weight(int delta) { weight_ += delta;    }
  void dec_weight(int delta) { weight_ -= delta;    }
  int  get_weight()          { return weight_;      }
};

class Node {
 protected:
  int id_;
  std::vector<Edge*> arriving_;
  std::vector<Edge*> departing_;

  void remove_edge(Edge* edge, std::vector<Edge*>& edges){
    unsigned int ins_index = 0;
    for (unsigned int i = 0; i < edges.size(); i++)
      if (edges[i] != edge)
	edges[ins_index++] = edges[i];
    assert(ins_index+1 == edges.size());
    edges.pop_back();
  }

 public:
  Node (int id){
    id_ = id;
  }

  int get_id(){ return id_; }
  std::vector<Edge*>& get_incident_edges() { return arriving_;         }
  std::vector<Edge*>& get_departing_edges(){ return departing_;        }
  int num_departing_edges()                { return departing_.size(); }
  int num_incident_edges()                 { return arriving_.size();  }

  void add_edge(Edge* edge){
    bool added = false;
    if (edge->get_source() == id_){
      departing_.push_back(edge);
      added = true;
    }
    if (edge->get_destination() == id_){
      arriving_.push_back(edge);
      added = true;
    }
    assert(added);
  }

  void remove_edge(Edge* edge){
    bool removed = false;
    if (edge->get_source() == id_){
      remove_edge(edge, departing_);
      removed = true;
    }
    if (edge->get_destination() == id_){
      remove_edge(edge, arriving_);
      removed = true;
    }
    assert(removed);
  }

  void get_child_nodes(std::vector<int>& children){
    for (unsigned int i = 0; i < departing_.size(); i++)
      children.push_back(departing_[i]->get_destination());
  }
};


class DirectedGraph {
protected:
  int num_nodes_;
  std::vector<Node*> nodes_;
  std::vector<Edge*> edges_;
  std::vector<std::string> node_labels_;
  std::map<std::string, int> node_map_;
  
public:
  DirectedGraph(){
    num_nodes_ = 0;
  }

  bool can_sort_topologically();

  bool has_cycles(){
    return !can_sort_topologically();
  }

  ~DirectedGraph(){
    for (unsigned int i = 0; i < edges_.size(); i++)
      delete edges_[i];
  }

  bool has_node(std::string& value){
    return node_map_.find(value) != node_map_.end();
  }

  Node* get_node(std::string& value){
    auto node_iter = node_map_.find(value);
    if (node_iter != node_map_.end())
      return nodes_[node_iter->second];

    node_labels_.push_back(value);
    node_map_[value] = num_nodes_;
    nodes_.push_back(new Node(num_nodes_++));
    return nodes_.back();
  }

  const std::string& get_node_label(int node_id){
    return node_labels_[node_id];
  }

  void increment_edge(std::string& val_1, std::string& val_2, int delta=1);
  void print(std::ostream& out);
};

#endif
