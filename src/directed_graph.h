#ifndef DIRECTED_GRAPH_H_
#define DIRECTED_GRAPH_H_

#include <assert.h>
#include <iostream>
#include <map>
#include <vector>

class Node;

class Edge {
 protected:
  int id_;
  int source_;
  int destination_;
  int weight_;

 public:
  Edge(int id, int source, int destination, int weight){
    id_          = id;
    source_      = source;
    destination_ = destination;
    weight_      = weight;
  }

  int  get_id()            const { return id_;            }
  int  get_source()        const { return source_;        }
  int  get_destination()   const { return destination_;   }
  int  get_weight()        const { return weight_;        }
  void set_id(int id)            { id_          = id;     }
  void set_source(int source)    { source_      = source; }
  void set_destination(int dest) { destination_ = dest;   }
  void inc_weight(int delta)     { weight_     += delta;  }
  void dec_weight(int delta)     { weight_     -= delta;  }
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
  explicit Node (int id){
    id_ = id;
  }

  int get_id() const { return id_; }
  void set_id(int id){ id_ = id;   }
  std::vector<Edge*>& get_incident_edges()  { return arriving_;   }
  std::vector<Edge*>& get_departing_edges() { return departing_;  }
  int num_departing_edges()           const { return departing_.size(); }
  int num_incident_edges()            const { return arriving_.size();  }

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

  void get_child_nodes(std::vector<int>& children) const {
    for (unsigned int i = 0; i < departing_.size(); i++)
      children.push_back(departing_[i]->get_destination());
  }
};


class DirectedGraph {
 private:
  // Private unimplemented copy constructor and assignment operator to prevent operations
  DirectedGraph(const DirectedGraph& other);
  DirectedGraph& operator=(const DirectedGraph& other);

protected:
  std::vector<Node*> nodes_;
  std::vector<Edge*> edges_;
  std::vector<std::string> node_labels_;
  std::map<std::string, int> node_map_;
  
public:
  bool can_sort_topologically() const;

  bool has_cycles() const {
    return !can_sort_topologically();
  }

  DirectedGraph(){}

  ~DirectedGraph(){
    for (unsigned int i = 0; i < nodes_.size(); i++)
      delete nodes_[i];
    for (unsigned int i = 0; i < edges_.size(); i++)
      delete edges_[i];
  }

  bool has_node(const std::string& value) const {
    return node_map_.find(value) != node_map_.end();
  }

  Node* get_node(const std::string& value){
    auto node_iter = node_map_.find(value);
    if (node_iter != node_map_.end())
      return nodes_[node_iter->second];

    node_labels_.push_back(value);
    node_map_[value] = nodes_.size();
    nodes_.push_back(new Node(nodes_.size()));
    return nodes_.back();
  }

  const std::string& get_node_label(int node_id) const {
    return node_labels_[node_id];
  }

  void increment_edge(const std::string& val_1, const std::string& val_2, int delta=1);
  void print(std::ostream& out) const;
};

#endif
