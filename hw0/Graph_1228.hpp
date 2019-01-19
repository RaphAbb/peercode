#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
// Added by me
#include <map>
#include <algorithm>

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
  Graph() {}

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
    Node() {} // this invalid constructor does not do anything

    /** Return this node's position. */
    const Point& position() const {
      // proxy design: let the graph object handles everything
      assert(id_ < graph_->size());
      return graph_->point_list_[id_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return this->id_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return ((n.id_ == id_) && (n.graph_ == graph_));
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
      return (id_ < n.id_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type id_;
    // private constructor accessible in Graph class
    Node(const Graph* graph, const size_type id)
        : graph_(const_cast<Graph*>(graph)), id_(id) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return point_list_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return point_list_.size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position) {
    size_type s = size();
    Node new_node = Node(this, s);
    point_list_.push_back(position);
    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return n.id_ < size();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < size());
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
    Edge() {}

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, node1_id_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, node2_id_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (e.id_ == id_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (id_ < e.id_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type node1_id_; // id of node 1
    size_type node2_id_; // id of node 2
    size_type id_;       // id of the edge
    Edge(const Graph* graph, const size_type node1_id, const size_type node2_id, 
        const size_type id)
        : graph_(const_cast<Graph*>(graph)), node1_id_(node1_id), node2_id_(node2_id), 
          id_(id) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_list_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());
    return edge_list_[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert(has_node(a) && has_node(b));
    size_type id_min = std::min(a.id_, b.id_);
    size_type id_max = std::max(a.id_, b.id_);
    // check if an edge might exist in the first place
    if (adj_list_.find(id_min) != adj_list_.end()){
      for (size_type i : adj_list_.at(id_min))
        if (edge_list_[i].node2().id_ == id_max) return true;
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
    assert(has_node(a) && has_node(b) && !(a == b));
    size_type id_min = std::min(a.id_, b.id_);
    size_type id_max = std::max(a.id_, b.id_);
    size_type i {find_edge(a, b)};

    if (i == num_edges()){
      // the edge does not exist, create it
      Edge new_edge = Edge(this, id_min, id_max, i);
      edge_list_.push_back(new_edge);
      adj_list_[id_min].push_back(i);
      return new_edge;
    }  
    // the edge exists, return it
    else return edge_list_[i];
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    point_list_.clear();
    edge_list_.clear();
    adj_list_.clear();
  }

 private:

  /** For now, the graph is represented with an node and edge vectors. */
  std::vector<Point> point_list_;
  std::vector<edge_type> edge_list_;
  /** Kind of adjacency list
   * It relate a node with a vector of edge ids' of the strictly greater
   * nodes connected to it. Don't know yet if it is a good data structure,
   * will see in future assignments if I have to modify it.
   */
  std::map<size_type, std::vector<size_type>> adj_list_;


  /** Find an edge and return its id
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return either the id of the edge if it exists,
   *         otherwise the number of edges in the graph.
   */
  size_type find_edge(const Node& a, const Node& b){
    assert(has_node(a) && has_node(b) && !(a == b));    
    size_type id_min = std::min(a.id_, b.id_);
    size_type id_max = std::max(a.id_, b.id_);
    if (adj_list_.find(id_min) != adj_list_.end()){
      for (size_type i : adj_list_.at(id_min))
        if (edge_list_[i].node2().id_ == id_max) return i;
    }
    return num_edges();
  }

};

#endif // CME212_GRAPH_HPP
