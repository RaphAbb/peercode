#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <iostream>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
 private:

  struct node_element;
  struct edge_element;

 public:

  using graph_type = Graph;
  using size_type  = unsigned;

  class Node;
  class Edge;

  using node_type = Node;
  using edge_type = Edge;

  Graph() 
    : node_map_(), edge_map_(), nodes_to_edge_(), n_nodes_(0), n_edges_(0) {
  }

  /** Default destructor */
  ~Graph() = default;

  //  ---------------------------------- Nodes ----------------------------- //

  // @class Graph::Node. Proxy object to get info about Graph's nodes.
  class Node {
   public:

    // Invalid constructor. Valid ctor is private.
    Node() {
    }

    // Return this node's position.
    const Point& position() const {
      return get().position;
    }

    // Return this node's index, a number in the range [0, graph_size). //
    size_type index() const {
      return this->nid_;
    }

    // Return this node's container, a pointer to a Graph.
    const Graph* graph() const {
      return this->graph_;
    }

    // Test whether this node and @a n are equal.
    // Equal nodes have both the same graph and the same index.
    bool operator==(const Node& n) const {
      if (this->graph()    == n.graph()  &&
          this->index()    == n.index()  &&
          this->position() == n.position()) 
      {
        return true; 
      } else {
        return false;
      }

    }

    // Test whether this node is less than @a n in a global order.
    // It has no geometric meaning.
    bool operator<(const Node& n) const {
      if (this->index() == n.index()) {
        return true; 
      } else {
        return false;
      }
    }

   private:

    friend class Graph; // Proxy design pattern necessitates Graph access.

    Graph* graph_;      // Pointer to parent Graph class (container)
    size_type nid_;     // Unique id. Key for graph_'s unordered_map.

    // The only valid constructor. Accessible only through graph_.
    Node(const Graph* graph, size_type nid)
      : graph_(const_cast<Graph*>(graph)), nid_(nid) {
    }
    
    // Grab node with Node ID = id from Graph.
    node_element& get() const {
      return graph_->node_map_[this->nid_];
    }

  };

  // -------------------------- Graph Node Methods ------------------------- //

  // Return the number of nodes in the graph.
  size_type size() const {
    return n_nodes_;
  }
  size_type num_nodes() const {
    return size();
  }

  // Add a new Node to Graph. 
  // Does not check if a Node with identical position already exists.
  Node add_node(const Point& position) {
   
    // Create the new Node, with id = n_nodes_.
    Node new_node(this, n_nodes_);

    // Create and populate the associated node_element data structure.
    node_element node_el;
    node_el.position = position;

    // Store node in Graph's map of nodes, update n_nodes_.
    node_map_[n_nodes_] = node_el;
    ++n_nodes_;

    return new_node;  // Return a valid node.
  }

  // Determine if a Node belongs to Graph.
  // Checks for identical Node ID and position.
  bool has_node(const Node& n) const {

    // Check if Node ID exists
    assert(n.index() < size() and n.index() >= 0);

    // Compare positions of @a n and node_map[@a n.nid_].
    auto search = node_map_.find(n.index());
    if (search->second.position == n.position()) {
      return true;
    } else {
      return false;
    }
  }

  // Return the Node with index @a i, for 0 =< i < num_nodes().
  Node node(size_type i) const {
    assert(i < size() and i >= 0);
    return Node(this, i);  // Valid Node.
  }

  //  ---------------------------------- Edges ----------------------------- //

  /** @class Graph::Edge
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge {

   public:

    // Invalid constructor. Valid ctor is private.
    Edge() {
    }

    // STUDENT ADDED THE FOLLOWING MEMBER
    // Return this Edge's index.
    // Edge indices are keys in edge_map_, a member of Graph which holds
    //   Edge data. It's used for lookups in the same way the Node ID is.
    // Adds moderate interpretability to inequality methods.
    size_type index() const {
      return eid_;
    }

    // Return a Node of this Edge
    Node node1() const {
      size_type nid = graph_->edge_map_[eid_].node1;
      return graph_->node(nid);
    }

    // Return the other Node of this Edge 
    Node node2() const {
      size_type nid = graph_->edge_map_[eid_].node2;
      return graph_->node(nid);
    }

    // Test whether this edge and @a e are equal (i.e. equivalent).
    // Equal edges represent the same undirected edge between two nodes.
    bool operator==(const Edge& e) const {
      if  ((this->node1() == e.node1() and this->node2() == e.node2())
        or (this->node1() == e.node2() and this->node2() == e.node1())) {
          return true;
        } else {
          return false;
        }
    }

    // Test whether this edge is less than @a e in a global order.
    // Compares Edge IDs.
    // Result has no interpretive meaning.
    bool operator<(const Edge& e) const {
      if (this->index() < e.index()) {
        return true;
      } else {
        return false;
      }
    }

   private:

    friend class Graph; // Proxy design pattern necessitates Graph access.

    Graph* graph_;      // Pointer to parent Graph class (container)

    // STUDENT ADDED THE FOLLOWING MEMBER
    size_type eid_;     // Unique id. Key for graph_'s unordered_map.

    // The only valid constructor. Accessible only through graph_.
    Edge(const Graph* graph, size_type eid)
      : graph_(const_cast<Graph*>(graph)), eid_(eid) {
    }
    
    // Grab Edge with Edge ID = eid_ from Graph.
    edge_element& get() const {
      return graph_->edge_map_[this->eid_];
    } 

  };

  // -------------------------- Graph Edge Methods ------------------------- //

  // Return the total number of edges in the graph.
  size_type num_edges() const {
    return n_edges_;
  }

  // Return the Edge with index i.
  Edge edge(size_type i) const {
    if (i >= num_edges() or i < 0) {
      return Edge();         // Invalid Edge
    } else {
      return Edge(this, i);  // Valid Edge
    }
  }

  // Test whether two nodes are connected by an edge.
  // Complexity: O(num_edges()).
  bool has_edge(const Node& a, const Node& b) const {

    // Get sorted Node indices.
    std::pair <node_type, node_type> sorted = node_ID_sort(a, b);
    
    // Search existing nodes for an Edge connecting them.
    std::unordered_map<std::size_t, size_type>::const_iterator search;
    search = node_connection_search(sorted.first, sorted.second);

    if (search != nodes_to_edge_.end()) {
      return true;
    } else {
      return false;
    }
  }

  // Add an edge to the graph, or return the current edge if it already exists.
  // Invalidates neither Edge indices nor outstanding Edge objects.
  // Complexity: No more than O(num_nodes() + num_edges()), hopefully less.
  Edge add_edge(const Node& a, const Node& b) {

    assert(has_node(a) and has_node(b));

    // Get sorted Node indices.
    std::pair <node_type, node_type> sorted = node_ID_sort(a, b);
    size_type id1 = sorted.first.index(), id2 = sorted.second.index();
    
    // Search existing nodes for an Edge connecting them.
    std::unordered_map<std::size_t, size_type>::const_iterator search;
    search = node_connection_search(sorted.first, sorted.second);

    if (search != nodes_to_edge_.end()) {

      // Existing Edge is found, return it.
      return Edge(this, search->second);

    } else {
       
      // Create the new Edge, with id = n_edges_.
      Edge new_edge(this, n_edges_);

      // Create and populate the associated edge_element data structure.
      edge_element edge_el;
      edge_el.node1 = id1;
      edge_el.node2 = id2;

      // Store Edge in Graph's map of Edges, update n_edges_.
      edge_map_[n_edges_] = edge_el;

      std::size_t hash = int_hash(id1, id2);
      nodes_to_edge_[hash] = n_edges_;

      ++n_edges_;

      return new_edge;  // Return a new, valid Edge.
    }

  }

  // Remove all nodes and edges from this graph.
  // Invalidates all outstanding Node and Edge objects.
  void clear() {
    // Reset Counters
    n_nodes_ = 0;
    n_edges_ = 0;

    // Reset Maps
    node_map_.clear();
    edge_map_.clear();
    nodes_to_edge_.clear();
  }

  // Sort two given Nodes by their ID.
  std::pair <node_type, node_type> 
        node_ID_sort(const Node& a, const Node& b) const {

    std::pair <node_type, node_type> sorted;

    if (a.index() < b.index()) {
      sorted = std::make_pair(a, b);
    } else if (b.index() < a.index()) { 
      sorted = std::make_pair(b, a);
    } else {
      std::string err = "Warning: sort requested between two nodes both of ID ";
      std::cout << err << a.index() << std::endl;
      sorted = std::make_pair(Node(), Node());
    }
    
    return sorted;
  }

  // Search for an Edge connecting the two given Nodes, a and b.
  // Both has_edge and add_edge need the value of this search operation.
  std::unordered_map<std::size_t, size_type>::const_iterator
        node_connection_search(const Node& a, const Node& b) const {

    assert(a.index() < b.index());

    // Search edges for one joining this pair of Node indices.
    std::size_t hash = int_hash(a.index(), b.index());

    std::unordered_map<std::size_t, size_type>::const_iterator search;
    search = nodes_to_edge_.find(hash);
    return search;
  }


 private:

  // data structures for nodes and edges
  struct node_element {
    Point position;
  };

  struct edge_element {
    size_type node1;  // ID of Node 1
    size_type node2;  // ID of Node 2
  };

  // maps for nodes and edges. {n,e}id mapped to {node,edge}_element.
  std::unordered_map<size_type, node_element> node_map_;
  std::unordered_map<size_type, edge_element> edge_map_;

  // Map Node IDs to Edge ID. 
  // Allows for quick check of whether an Edge exists (i.e. for has_edge).
  std::unordered_map<std::size_t,  size_type> nodes_to_edge_;

  // count of each element.
  // accessible through num_nodes() and num_edges() public methods.
  // also used for assigning new indices to preserve the invariant:
  // For all i with 0 ≤ i < graph.n_nodes(), graph.Node(i).nid() == i.
  size_type n_nodes_;
  size_type n_edges_;

  // Forbid copy or assignment of Graphs
  Graph(const Graph&) = delete;
  Graph& operator=(const Graph&) = delete;

  
  // ------------------------ Miscellaneous Functions ----------------------- //

  // Hashes two Node IDs into a unique key for nodes_to_edge_ map.
  std::size_t int_hash(const size_type& nid1, const size_type& nid2) const {
      std::size_t seed = 0;
      boost::hash_combine(seed, nid1);
      boost::hash_combine(seed, nid2);

      return seed;
  }

};

#endif // CME212_GRAPH_HPP
