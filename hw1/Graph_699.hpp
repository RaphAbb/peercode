#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <set>
#include <string>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = double>
class Graph {
 private:

  //
  // PRIVATE PREDECLARATIONS
  //

  struct internal_node;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
  /** Type of node value */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  using node_iterator = NodeIterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  using incident_iterator = IncidentIterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
    : nodes_(), size_(0), idx_to_edges_(), edges_to_idx_(), edge_size_(0),
      node_connectivity_() {
  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node : private totally_ordered<Node> {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node()
      : graph_(nullptr), id_() {
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_[id_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return id_;
    } 
    /** Returns the value of the node */
    node_value_type& value() {
      return graph_->nodes_[id_].value;
    }

    /** Overloading const version */
    const node_value_type& value() const {
      return graph_->nodes_[id_].value;
    }

    // Number of nodes this node is connected to O(1)
    size_type degree() const {
      return graph_->node_connectivity_[id_].size();
    };

    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, id_, 0);
    };

    incident_iterator edge_end() const {
      return IncidentIterator(graph_, id_, degree());
    };


    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (n.graph_ == graph_ && n.id_ == id_)
        return true;
      return false;
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
      if (id_ < n.id_)
        return true;
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Pointer back to the Graph container
    Graph* graph_;
    // This node's index number
    size_type id_;

    // Private Constructor
    Node(const Graph* graph, size_type id)
      : graph_(const_cast<Graph*>(graph)), id_(id) {
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return size_;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    internal_node temp_node { position, value};
    nodes_.push_back(temp_node);
    ++size_;
    return Node(this, size_-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.id_ < size_ && n.graph_ == this)
      return true;
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < size_);
    return Node(this, i);
  }

  //
  // END NODES
  //

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge()
      : graph_(nullptr), id_(0) {
    }

    /** Return a node of this Edge */
    Node node1() const {
      std::vector<size_type> node_pair = graph_->idx_to_edges_[id_];
      return Node(graph_, node_pair[0]);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      std::vector<size_type> node_pair = graph_->idx_to_edges_[id_];
      return Node(graph_, node_pair[1]);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (id_ == e.id_)
        return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (id_ < e.id_)
        return true;
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    
    Edge(const Graph* graph, size_type id)
      : graph_(const_cast<Graph*>(graph)), id_(id) {
    }

    // Pointer back to graph container
    Graph* graph_;
    // ID of the edge
    size_type id_;

  };

  // 
  // END EDGES
  //

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_size_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert(a.id_ < num_nodes() && b.id_ < num_nodes());
    // constructs a string, node_pair_str,
    // min(a.id_, b.id_) + " " max(a.id_, b.id_)
    // used to lookup the existence of that edge, o(1) complexity
    size_type small_id = std::min(a.id_, b.id_);
    size_type large_id = std::max(a.id_, b.id_);
    std::string node_pair_str;
    node_pair_str = std::to_string(small_id) + "_" + std::to_string(large_id);

    auto search = edges_to_idx_.find(node_pair_str);
    if ( search != edges_to_idx_.end() )
      return true;
    return false;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {
    // Check if nodes are valid
    assert(a.id_ < num_nodes() && b.id_ < num_nodes());

    // Constructs a string, node_pair_str,
    // min(a.id_, b.id_) + " " max(a.id_, b.id_)
    // used to lookup the existence of that edge, O(1) complexity
    size_type small_id = std::min(a.id_, b.id_);
    size_type large_id = std::max(a.id_, b.id_);
    std::vector<size_type> node_pair {small_id, large_id};
    std::string node_pair_str;
    node_pair_str = std::to_string(small_id) + "_" + std::to_string(large_id);

    // Check if edge already exists
    if ( has_edge(a, b) ) {
      unsigned edge_id = edges_to_idx_.at(node_pair_str);
      return Edge(this, edge_id);
    } 

    // Create new edge
    edges_to_idx_[node_pair_str] = edge_size_;
    idx_to_edges_.push_back(node_pair);
    // Storing the connectivity
    node_connectivity_[a.id_].push_back(b.id_);
    node_connectivity_[b.id_].push_back(a.id_);

    // Update num_edges()
    ++edge_size_;
    return Edge(this, edge_size_-1);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    idx_to_edges_.clear();
    edges_to_idx_.clear();
    edge_size_ = 0;
    size_ = 0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator> {
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

    Node operator*() const {
      return Node(graph_, p_);
    }

    NodeIterator& operator++() {
      p_++;
      return *this;
    }

    bool operator==(const node_iterator& x) const {
      return p_ == x.p_;
    }

   private:
    friend class Graph;
    Graph<V>* graph_;
    size_type p_;
    NodeIterator(const Graph* graph, size_type p) 
      : graph_(const_cast<Graph*>(graph)), p_(p) {
    }
  };

  //
  // END NODE ITERATOR
  //

  node_iterator node_begin() const { return NodeIterator(this, 0); }
  node_iterator node_end() const { return NodeIterator(this, num_nodes()); }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator> {
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

    Edge operator*() const {
      // Gets the adjacent node id from the iterator
      size_type a_node_id = graph_->node_connectivity_[h_node_id_][a_node_iter_];

      // Constructs a string, node_pair_str,
      // min(a.id_, b.id_) + " " max(a.id_, b.id_)
      // used to lookup the existence of that edge, O(1) complexity
      size_type small_id = std::min(h_node_id_, a_node_id);
      size_type large_id = std::max(h_node_id_, a_node_id);
      std::string node_pair_str;
      node_pair_str = std::to_string(small_id) + "_" + std::to_string(large_id);
      unsigned edge_id = graph_->edges_to_idx_.at(node_pair_str);
      Edge e = Edge(graph_, edge_id);


      // If Node1 isn't the "home node", swap the ordering such that it is
      if (e.node1().index() != h_node_id_) {
        std::vector<size_type> node_pair = graph_->idx_to_edges_[edge_id];
        std::swap(node_pair[0], node_pair[1]);
        graph_->idx_to_edges_[edge_id] = node_pair;
      }

      return Edge(graph_, edge_id);
    }

    IncidentIterator& operator++() {
      a_node_iter_++;
      return *this;
    }

    bool operator==(const IncidentIterator& x) const {
      return a_node_iter_ == x.a_node_iter_;
    }


   private:
    friend class Graph;

    Graph* graph_;
    size_type h_node_id_; // ID of the central node
    size_type a_node_iter_; // Indexing into node_connectivity_[h_node_id_]

    // Private constructor to make IncidetIterator tied to a certain graph and home node
    IncidentIterator(const Graph* graph, size_type h_node_id, size_type a_node_iter) 
      : graph_(const_cast<Graph*>(graph)), h_node_id_(h_node_id), 
        a_node_iter_(a_node_iter) {
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator> {
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

    Edge operator*() const {
      return Edge(graph_, edge_id_);
    };

    EdgeIterator& operator++() {
      edge_id_++;
      return *this;
    };

    bool operator==(const EdgeIterator& x) const {
      return edge_id_ == x.edge_id_;
    };

   private:
    friend class Graph;
    Graph* graph_;
    size_type edge_id_;

    EdgeIterator(const Graph* graph, size_type edge_id)
      : graph_(const_cast<Graph*>(graph)), edge_id_(edge_id) {
      };
  };

  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  };

  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  };

 private:

  /* 
   * Node helper functions, data members, etc.
   */

  struct internal_node {
    Point position;
    node_value_type value;
  };

  std::vector<internal_node> nodes_;

  size_type size_;

  /* 
   * Edge helper functions, data members, etc.
   */

  // A vector of edges, where each edge is represented by a node pair
  std::vector<std::vector<size_type>> idx_to_edges_;

  // A map to lookup the existence of an edge between 2 nodes
  // Key: node_a_id + "_" + node_b_id, Value: edge_id
  std::unordered_map<std::string, size_type> edges_to_idx_;

  size_type edge_size_;

  // Index of node, and all the nodes its connected to
  // Equivalent of a sparse matrix
  std::unordered_map<size_type, std::vector<size_type>> node_connectivity_;


};


#endif // CME212_GRAPH_HPP

