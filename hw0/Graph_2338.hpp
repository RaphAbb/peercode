#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
#include <set>

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

  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  // For storage of Nodes
  std::vector<Point> points_;

  unsigned next_index_;
  unsigned next_edge_index_;

  // For storage of Edges
  class Internal_Edge;
  std::vector<std::vector<unsigned>> edges_;
  std::vector<std::vector<Internal_Edge>> adjacency_list;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  using size_type = unsigned;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph(): next_index_(0), next_edge_index_(0){
    // HW0: YOUR CODE HERE
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
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->points_[index_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {

      // With the current functionality the index will never be greater than graph_size
      // since we are not enabling the user to delete any nodes.
      unsigned index = index_;

      return index;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      bool equals = (graph_ == n.graph_) and (index_ == n.index_);
      return equals;
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
      // HW0: YOUR CODE HERE

      return index_ < n.index();
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    Graph* graph_;
    size_type index_;

    /** Private Constructor **/
    Node(const Graph* graph, size_type index)
              : graph_(const_cast<Graph*>(graph)), index_(index){
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {

    return points_.size();
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
    // Add point to graph data structure
    points_.push_back(const_cast<Point&>(position));
    adjacency_list.push_back(std::vector<Internal_Edge>());

    // Make proxy node instance
    Node n = Node(this, next_index_);

    // Update graph info
    next_index_ = next_index_ + 1;

    return n;        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    // Todo: check if size is compatible
    return (points_[n.index()] == n.position()) and (this == n.graph_);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i >= 0 and i <= num_nodes());
    Node n = Node(this, i);
    return n;
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
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      unsigned node1_id = node1_id_;
      Node node = graph_->node(node1_id);

      return node;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      unsigned node2_id = node2_id_;
      Node node = graph_->node(node2_id);

      return node;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {

      bool equal_same_order = (node1_id_ == e.node1_id_) and (node2_id_ == e.node2_id_);
      bool equal_opposite_order = (node1_id_ == e.node2_id_) and (node2_id_ == e.node1_id_);

      bool equal = equal_same_order or equal_opposite_order;

      return equal;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {

      bool less = index_ < e.index_;

      return less;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    Graph* graph_;
    unsigned index_;

    unsigned node1_id_;
    unsigned node2_id_;

    /** Constructor for Edge */
    Edge(const Graph* graph, unsigned index, unsigned node1_id, unsigned node2_id):
      graph_(const_cast<Graph*>(graph)), index_(index), node1_id_(node1_id),
      node2_id_(node2_id) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {

    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {

    assert((i >= 0) and (i <= num_edges()));
    std::vector<unsigned> pair = edges_[i];
    unsigned node_a_index = pair[0];
    unsigned node_b_index = pair[1];

    Edge e = Edge(this, i, node_a_index, node_b_index);

    return e;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {

    // Check if a and b are in this graph
    assert(a.graph_ == b.graph_);

    unsigned node_a_index = a.index();
    unsigned node_b_index = b.index();

    std::vector<Internal_Edge> edges_to_neighbors_a = adjacency_list[node_a_index];

    bool found = false;

    for (auto e : edges_to_neighbors_a) {

      // Order of node_1 and node_2 is arbitrary so we have to check both orders
      if (e.node1_id_ == node_a_index and e.node2_id_ == node_b_index) {
        found = true;
      }

      if (e.node1_id_ == node_b_index and e.node2_id_ == node_a_index) {
        found = true;
      }
    }

    return found;
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

    assert(not(a == b));

    unsigned node_a_id = a.index();
    unsigned node_b_id = b.index();

    // If edge exists, find the edge index and return the edge
    if (this->has_edge(a, b)){
      std::vector<Internal_Edge> neigbors_a = adjacency_list[node_a_id];

      Internal_Edge edge_internal;

      // Find the appropriate edge
      for (auto e: neigbors_a) {
        if (e.node1_id_ == node_a_id and e.node2_id_ == node_b_id) {
          edge_internal = e;
        }

        if (e.node1_id_ == node_b_id and e.node2_id_ == node_a_id) {
          edge_internal = e;
        }
      }

      Edge edge = Edge(edge_internal.graph_, edge_internal.node1_id_,
              edge_internal.node2_id_, edge_internal.index_);

      return edge;
    }

    // If the edge does not exist yet

    // Add edge to adjacency list
    Internal_Edge e_internal = Internal_Edge(this, next_edge_index_, node_a_id, node_b_id);
    adjacency_list[node_a_id].push_back(e_internal);
    adjacency_list[node_b_id].push_back(e_internal);

    // Add edge to internal edge_ vector
    std::vector<unsigned> pair;
    pair.push_back(node_a_id);
    pair.push_back(node_b_id);
    edges_.push_back(pair);

    // Make new edge to return
    Edge e = Edge(this, next_edge_index_ , node_a_id, node_b_id);

    // Increment edge counter
    next_edge_index_++;

    return e;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    points_ = std::vector<Point>();

    next_index_ = 0;
    next_edge_index_ = 0;

    // For storage of Edges
    edges_ = std::vector<std::vector<unsigned>>();
    adjacency_list = std::vector<std::vector<Internal_Edge>>();
  }

 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  class Internal_Edge {

      friend class Graph;

      const Graph *graph_;
      unsigned index_;

      unsigned node1_id_;
      unsigned node2_id_;

      // Invalid constructor
      Internal_Edge() {
      }

      Internal_Edge(const Graph* graph, unsigned index, unsigned node1_id, unsigned node2_id):
              graph_(const_cast<Graph*>(graph)), index_(index), node1_id_(node1_id),
              node2_id_(node2_id) {
      }
  };

};

#endif // CME212_GRAPH_HPP
