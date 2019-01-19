#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

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
class Graph {
 private:


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

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
    : nodes_(), size_(0), idx_to_edges_(), edges_to_idx_(), edge_size_(0) {
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
  class Node {
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
      return graph_->nodes_[id_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return id_;
    }

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
  Node add_node(const Point& position) {
    nodes_.push_back(position);
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
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge {
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

 private:

  /* 
   * Node helper functions, data members, etc.
   */

  std::vector<Point> nodes_;

  size_type size_;

  /* 
   * Edge helper functions, data members, etc.
   */

  std::vector<std::vector<size_type>> idx_to_edges_;
  // Key: node_a_id + "_" + node_b_id, Value: edge_id
  std::unordered_map<std::string, size_type> edges_to_idx_;

  size_type edge_size_;

};


#endif // CME212_GRAPH_HPP

