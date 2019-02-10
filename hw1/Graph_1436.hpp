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
template <typename V>
class Graph {

 private:

  struct node_element;
  struct edge_element;

 public:

  using graph_type = Graph<V>;

  // Type of indexes and sizes.
  //    Return type of Graph::Node::index(), Graph::num_nodes(),
  //    Graph::num_edges(), and argument type of Graph::node(size_type) 
  using size_type     = unsigned;
  using node_map_type = std::unordered_map<size_type, node_element>;
  using edge_map_type = std::unordered_map<size_type, edge_element>;

  // Predeclaration of Node class and associated types/classes 
  class Node;
    class NodeIterator;
    using node_type       = Node;
    using node_value_type = V;
    using node_iterator   = NodeIterator;

  // Predeclaration of Edge class and associated types/classes 
  class Edge;
    class EdgeIterator;
    using edge_type = Edge;
    using edge_iterator = EdgeIterator;

  // Predeclaration of IncidentIterator 
  class IncidentIterator;
    using incident_iterator = IncidentIterator;


  // Default constructor
  Graph() 
    : node_map_(), edge_map_(), nodes_to_edge_(), n_nodes_(0), n_edges_(0) {
  }

  // Default destructor 
  ~Graph() = default;

  //  ---------------------------------- Nodes ----------------------------- //

  // @class Graph::Node. Proxy object to get info about Graph's nodes.
  class Node : private totally_ordered<Node> {
   public:

    // Invalid constructor. Valid ctor is private.
    Node() {
    }

    // Grab node with Node ID = id from Graph.
    node_element& get() const {
      return graph_->node_map_[this->nid_];
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

    // Return reference to edit this node's value.
    node_value_type& value() {
      return get().value;
    }

    // Return this node's value.
    const node_value_type& value() const {
      return get().value;
    }

    // Return this node's value.
    size_type degree() const {
      return get().connections.size();
    }

    // Return this node's connections (a vector).
    std::vector<size_type> connections() const {
      return get().connections;
    }

    incident_iterator edge_begin() const {
      return IncidentIterator(this->index(), get().connections.begin());
    }

    incident_iterator edge_end() const {
      return IncidentIterator(this->index(), get().connections.end());
    }

    // Test whether this node and @a n are equal.
    // Equal nodes have both the same graph and the same index.
    bool operator==(const Node& n) const {
      return (this->graph() == n.graph() && this->index() == n.index());
    }

    // Test whether this node is less than _n_ in a global order.
    // It has no geometric meaning.
    bool operator<(const Node& n) const {
      assert(graph_->has_node(n));
      return (this->index() < n.index());
    }

   private:

    friend class Graph; // Proxy design pattern necessitates Graph access.

    Graph* graph_;      // Pointer to parent Graph class (container)
    size_type nid_;     // Unique id. Key for graph_'s unordered_map.

    // The only valid constructor. Accessible only through graph_.
    Node(const Graph* graph, size_type nid)
      : graph_(const_cast<Graph*>(graph)), nid_(nid) {
    }
    
  };

  // -------------------------- Graph Node Methods ------------------------- //

  // Return the number of nodes in the graph.
  size_type size() const {
    return n_nodes_;
  }
  
  // Synonym of size().
  size_type num_nodes() const {
    return size();
  }

  // Add a new Node to Graph. 
  // Does not check if a Node with identical position already exists.
  Node add_node(const Point& position, 
                  const node_value_type& value = node_value_type()) {
   
    // Create and populate the associated node_element data structure.
    node_element node_el;
    node_el.node     = Node(this, n_nodes_);   // id = n_nodes.
    node_el.index    = n_nodes_;
    node_el.position = position;
    node_el.value    = value;
    // Leave node_el.connections unassigned.
    //    add_edge updates this attribute.

    // Store node in Graph's map of nodes, update n_nodes_.
    node_map_[n_nodes_] = node_el;
    ++n_nodes_;

    return node_el.node;
  }

  // Determine if a Node is a member of this Graph, comparing pointers.
  bool has_node(const Node& n) const {
    return this == n.graph();
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
  class Edge : private totally_ordered<Edge> { 

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

    // Swap node indices of edge. Allows node1() > node2() for IncidentIterator.
    // Grabbing an edge from edge_map_ should still be done with sorted keys.
    void swap_node_indices() {
      
      size_type nid1 = graph_->edge_map_[eid_].node1;
      size_type nid2 = graph_->edge_map_[eid_].node2;
      graph_->edge_map_[eid_].node1 = nid2;
      graph_->edge_map_[eid_].node2 = nid1;

      return;
    }

  };

  // -------------------------- Graph Edge Methods ------------------------- //

  // Return the total number of edges in the graph.
  size_type num_edges() const {
    return n_edges_;
  }

  // Return the Edge with index i.
  Edge edge(size_type i) const {
    if (i >= num_edges()) {
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
      edge_element    edge_el;
      edge_el.edge  = new_edge;
      edge_el.node1 = id1;
      edge_el.node2 = id2;

      // Store Edge in Graph's map of Edges.
      edge_map_[n_edges_]  = edge_el;
      std::size_t hash     = int_hash(id1, id2);
      nodes_to_edge_[hash] = n_edges_;

      // Add this edge to "connections" containers of the nodes it connects.
      node_map_[id1].connections.push_back(new_edge);
      node_map_[id2].connections.push_back(new_edge);

      // Increment n_edges_.
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

  // ----------------- Iterator Members, Associated Aliases ----------------- //

  // ------------------------------- NodeIterator --------------------------- //
  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. Wrapper around std::unordered_map iterator.*/
  class NodeIterator : private totally_ordered<NodeIterator> {
   public:

    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    // Return current node object.
    Node operator*() const {
      return it_->second.node;
    };

    // Advance one element in the map.
    NodeIterator& operator++() {
      ++it_;
      return *this;
    };

    // Assess equality of current node with another iterator's node.
    bool operator==(const NodeIterator& x) const {
      return this->it_ == x.it_;
    };

   private:

    friend class Graph;
    typename node_map_type::const_iterator it_;

    // Construct the iterator, instantiating a begin() pointer.
    NodeIterator(typename node_map_type::const_iterator it)
      : it_(it) {
    }

  };

  // Type alias for NodeIterator, returning the begin() pointer.
  node_iterator node_begin() const {
    return NodeIterator(node_map_.begin());
  };
  
  // Type alias for NodeIterator, returning the end() pointer.
  node_iterator node_end() const {
    return NodeIterator(node_map_.end());
  };

  // ----------------------------- IncidentIterator ------------------------- //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : totally_ordered<IncidentIterator> {
   public:

    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    // Return the current Edge the iterator is on. 
    // If necessary, underlying Edge indices will be adjusted so
    //    base node will be current_edge.node1().
    Edge operator*() const {

      edge_type current_edge = *this->it_;
      
      // Swap node indices if necessary to return current node first.
      if (base_node_ != current_edge.node1().index()) {
        current_edge.swap_node_indices();
        assert(base_node_ == current_edge.node1().index()); // check self, don't wreck self.
      }

      return current_edge;
    }

    // Move on to the next element.
    IncidentIterator& operator++() {
      ++it_;
      return *this;
    }

    // Compare Edges of two iterators
    bool operator==(const IncidentIterator& x) const {
      return this->it_ == x.it_;
    }

   private:

    friend class Graph;

    // Iterator always returns the index of its spawning node first. 
    // Store the corresponding node index here.
    size_type base_node_;

    typename std::vector<edge_type>::iterator it_;

    // Construct a valid IncidentIterator.
    IncidentIterator(size_type n, typename std::vector<edge_type>::iterator it)
      : base_node_(n), it_(it) {
    }

  };

  // ------------------------------- EdgeIterator --------------------------- //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    // Return current edge object.
    Edge operator*() const {
      return Edge(graph_, it_->second.edge.index());
    };

    // Advance one element in the map.
    EdgeIterator& operator++() {
      ++it_;
      return *this;
    };

    // Assess equality of current iterator with another.
    bool operator==(const EdgeIterator& x) const {
      return this->it_ == x.it_;
    };

   private:
    friend class Graph;
    const Graph* graph_;
    typename edge_map_type::const_iterator it_;

    // Construct the iterator, instantiating a begin() pointer.
    EdgeIterator(const Graph* g, typename edge_map_type::const_iterator it)
      : graph_(g), it_(it) {
    }
  };

  // Type alias for EdgeIterator, returning the begin() pointer.
  edge_iterator edge_begin() const {
    return EdgeIterator(this, edge_map_.begin());
  };
  
  // Type alias for EdgeIterator, returning the end() pointer.
  edge_iterator edge_end() const {
    return EdgeIterator(this, edge_map_.end());
  };

 // -------------------------- Private Graph Methods ------------------------ //

 private:

  // Data structures for nodes and edges
  struct node_element {
    Node                   node;
    size_type              index;
    std::vector<edge_type> connections;
    Point                  position;
    node_value_type        value;
  };

  struct edge_element {
    Edge      edge;   // Pointer back to the Edge.
    size_type node1;  // ID of Node 1
    size_type node2;  // ID of Node 2
  };

  // Maps for nodes and edges. {n,e}id mapped to {node,edge}_element.
  node_map_type node_map_;
  edge_map_type edge_map_;

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
