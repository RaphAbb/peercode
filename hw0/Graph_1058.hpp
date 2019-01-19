#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <stdexcept>
#include <vector>

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
  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

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
  Graph(){
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
      return graph_->points_.at(idx_);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return this->idx_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if(n.graph_ == this->graph_ and n.index() == this->index()){
        return true;
      }
      else{
        return false;
      }
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
      return index() < n.index();
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    const graph_type* graph_;
    size_type idx_;

    Node(graph_type* graph, size_type idx)
      : graph_(graph), idx_(idx) {

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
    points_.push_back(position);     // Add node's point to graph
    next_idx_++;                     // Increment index
    adj_list_.push_back(std::vector<InternalEdge>());
    return Node(this, next_idx_- 1); // Create & return node instance
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return this == n.graph_;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(const_cast<graph_type*>(this), i);
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
      return Node(graph_, node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, node2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if(e.node1_ == node1_ and e.node2_ == node2_){
        return true;
      }
      else if(e.node1_ == node2_ and e.node2_ == node1_){
        return true;
      }
      else{
        return false;
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return idx_ < e.idx_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    size_type node1_;
    size_type node2_;
    size_type idx_;
    graph_type* graph_;

    /** Construct valid edge object */
    Edge(size_type node1, size_type node2, size_type idx, graph_type* graph)
      : node1_(node1), node2_(node2), idx_(idx), graph_(graph) {

    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return next_edge_idx_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    if(i < 0 || i > next_edge_idx_-1){
      throw std::out_of_range("Edge index out of range when calling edge.");
    }
    return Edge(edges_.at(i).node1_, edges_.at(i).node2_, i, const_cast<graph_type*>(this));
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Check that nodes belong to this graph
    if((has_node(a) and has_node(b)) == false){
      throw std::runtime_error("Node(s) do not belong to graph in has_edge()");
    }
    // Loop through edges that indclude a and look for (a,b) edge
    for(auto edge : adj_list_.at(a.index())){
      if (edge.node1_ == a.index() and edge.node2_ == b.index()){
        return true;
      }
      else if (edge.node1_ == b.index() and edge.node2_ == a.index()){
        return true;
      }
    }
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
    // Check that a and b are distinct valid nodes of this graph
    if(a.graph_ != this || b.graph_ != this){
      throw std::runtime_error("One of the nodes do not belong to the graph.");
    }
    if(a == b){
      throw std::runtime_error("add_edge() received non-distinct nodes. Self-loops are frowned upon.");
    }

    // If the edge exists loop over edges for node a, find index, return edge
    if (has_edge(a, b)){
      size_type index;
      for(auto edge_ : adj_list_.at(a.index())){
        if(edge_.node1_ == a.index() and edge_.node2_ == b.index()){
          index = edge_.idx_;
          break;
        }
        else if (edge_.node1_ == b.index() and edge_.node2_ == a.index()){
          index = edge_.idx_;
          break;
        }
      }
      return edge(index);
    }
    // Else create new edge and store information in appropriate containers
    else{
      next_edge_idx_++;
      InternalEdge edge_ab(a.index(), b.index(), next_edge_idx_-1);
      adj_list_[a.index()].push_back(edge_ab);
      adj_list_[b.index()].push_back(edge_ab);
      edges_.push_back(NodePair(a.index(), b.index()));
      return Edge(a.index(), b.index(), next_edge_idx_-1, this);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 and num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // Reset counters
    next_idx_      = 0;
    next_edge_idx_ = 0;
    // Clear graph containers
    points_.clear();
    adj_list_.clear();
    edges_.clear();
  }

 private:
  // Structure to store a pair of node indices
  struct NodePair{
    NodePair(size_type a, size_type b)
      : node1_(a), node2_(b) {

    }
    size_type node1_;
    size_type node2_;
  };

  // Structure to store the node indices and index of an edge
  struct InternalEdge{
    InternalEdge(size_type a, size_type b, size_type idx)
      : node1_(a), node2_(b), idx_(idx){

    }
    size_type node1_;
    size_type node2_;
    size_type idx_;
  };

  std::vector<Point> points_;                       // Node idx -> Point
  std::vector<std::vector<InternalEdge>> adj_list_; // Map node idx -> Edges
  std::vector<NodePair> edges_;                     // Map node idxs -> edge idx
  size_type next_idx_ = 0;
  size_type next_edge_idx_ = 0;
};

#endif // CME212_GRAPH_HPP
