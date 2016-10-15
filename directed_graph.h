#ifndef DIRECTED_GRAPH_H_
#define DIRECTED_GRAPH_H_

#include <assert.h>
#include <map>
#include <vector>

class Node;

class Edge {
  Node* source_;
  Node* destination_;
  int weight_;

 public:
  Edge(Node* source, Node* destination, int weight){
    source_      = source;
    destination_ = destination;
    weight_      = weight;
  }

  Node* get_source()     { return source_; }
  Node* get_destination(){ return destination_; }

  void inc_weight(){ weight_++;      }
  void dec_weight(){ weight_--;      }
  int  get_weight(){ return weight_; }
  void set_weight(int weight){
    weight_ = weight;
  }
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
  std::vector<Edge*>& get_incident_edges() { return arriving_;              }
  std::vector<Edge*>& get_departing_edges(){ return departing_;             }
  void remove_incident_edge(Edge* edge)    { remove_edge(edge, arriving_);  }
  void remove_departing_edge(Edge* edge)   { remove_edge(edge, departing_); }
  void add_incident_edge(Edge* edge)       { arriving_.push_back(edge);     }
  void add_departing_edge(Edge* edge)      { departing_.push_back(edge);    }
  int num_incident_edges()                 { return arriving_.size();       }
  
  bool has_parent(Node* node){
    for (unsigned int i = 0; i < arriving_.size(); i++)
      if (arriving_[i]->get_source() == node)
	return true;
    return false;
  }

  bool has_child(Node* node){
    for (unsigned int i = 0; i < departing_.size(); i++)
      if (departing_[i]->get_destination() == node)
	return true;
    return false;
  }

  void get_child_nodes(std::vector<Node*>& children){
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
  DirectedGraph(){}

  bool can_sort_topologically();

  bool has_cycles(){
    return !can_sort_topologically();
  }

  ~DirectedGraph(){
    for (unsigned int i = 0; i < nodes_.size(); i++)
      delete nodes_[i];
    for (unsigned int i = 0; i < edges_.size(); i++)
      delete edges_[i];
  }

  void prune_edges(int min_weight){
    unsigned int ins_index = 0;
    for (unsigned int i = 0; i < edges_.size(); i++){
      if (edges_[i]->get_weight() < min_weight){
	edges_[i]->get_source()->remove_departing_edge(edges_[i]);
	edges_[i]->get_destination()->remove_incident_edge(edges_[i]);
	delete edges_[i];
      }
      else
	edges_[ins_index++] = edges_[i];
    }
    edges_.resize(ins_index);
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

  const std::string& get_node_label(Node* node){
    return node_labels_[node->get_id()];
  }

  const std::string& get_node_label(int node_id){
    return node_labels_[node_id];
  }

  void increment_edge(std::string& val_1, std::string& val_2){
    Node* source = get_node(val_1);
    Node* dest   = get_node(val_2);
    
    std::vector<Edge*> edges = dest->get_incident_edges();
    for (unsigned int i = 0; i < edges.size(); i++){
      if (edges[i]->get_source() == source){
	edges[i]->inc_weight();
	return;
      }
    }

    Edge* new_edge = new Edge(source, dest, 1);
    edges.push_back(new_edge);
    source->add_departing_edge(new_edge);
    dest->add_incident_edge(new_edge);
  }
};

#endif
