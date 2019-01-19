#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
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

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  struct NodeInternal;
  struct EdgeInternal;

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

  // Index of the true Node in the node_list_
  using NodeListIndex = size_type;

  // Index of true Edge in the edge_list_
  using EdgeListIndex = size_type;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
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
      // HW0: YOUR CODE HERE
      node_list_index_ = std::numeric_limits<EdgeListIndex>::max();
      graph_ptr_ = nullptr;
    }

    Node(NodeListIndex node_list_index, const Graph *graph_ptr)
      : node_list_index_(node_list_index), graph_ptr_(const_cast<Graph*>(graph_ptr)){
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      // check to make sure node is valid
      assert(node_list_index_ < (graph_ptr_->node_list_).size());
      return (graph_ptr_->node_list_)[node_list_index_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      assert(node_list_index_ < (graph_ptr_->node_list_).size());
      return node_list_index_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      assert(node_list_index_ < (graph_ptr_->node_list_).size());
      return (graph_ptr_ == n.get_graph_pointer()) &&
              (node_list_index_ == n.get_node_index());
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
      assert(node_list_index_ < (graph_ptr_->node_list_).size());
      return (graph_ptr_ == n.get_graph_pointer()) &&
              (node_list_index_ < n.get_node_index());
    }

    // Access the index of the node in node_list_.
    // NOTE: I define and enforce this to be the same as the node's index.
    NodeListIndex get_node_index() const {
      return node_list_index_;
    }

    // Access the pointer to the graph to which this node belongs.
    Graph *get_graph_pointer() const {
      return graph_ptr_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Index of node in the node_list_
    // NOTE: I enforce each node's index to be the same as its position in node_list_.
    NodeListIndex node_list_index_;

    // Pointer to the graph to which this node belongs.
    Graph *graph_ptr_;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return node_list_.size();
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
    // HW0: YOUR CODE HERE

    // obtain index for this new node
    const size_type new_index = node_list_.size();

    // add NodeInternal to graph using emplace
    node_list_.emplace_back(new_index, position);

    // create node proxy
    Node node_proxy = Node(new_index, this);

    return node_proxy;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return n.get_graph_pointer() == this;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < node_list_.size());
    return Node(i, this);
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
      // HW0: YOUR CODE HERE
      edge_list_index_ = std::numeric_limits<EdgeListIndex>::max();
      graph_ptr_ = nullptr;
    }

    Edge(EdgeListIndex edge_list_index, const Graph *graph_ptr)
      : edge_list_index_(edge_list_index), graph_ptr_(const_cast<Graph*>(graph_ptr)) {
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      assert(edge_list_index_ < (graph_ptr_->edge_list_).size());
      return (graph_ptr_->edge_list_)[edge_list_index_].node1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      assert(edge_list_index_ < (graph_ptr_->edge_list_).size());
      return (graph_ptr_->edge_list_)[edge_list_index_].node2;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      assert(edge_list_index_ < (graph_ptr_->edge_list_).size());
      // again, assuming here that edge index is same as edge list index
      return (graph_ptr_ == e.get_graph_pointer()) && (edge_list_index_ == e.get_edge_index());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      assert(edge_list_index_ < (graph_ptr_->edge_list_).size());
      return (graph_ptr_ == e.get_graph_pointer()) && (edge_list_index_ < e.get_edge_index());
    }

    // Access the index to the edge in edge_list_.
    EdgeListIndex get_edge_index() const {
      return edge_list_index_;
    }

    // Access the pointer to the graph to which this edge belongs.
    Graph *get_graph_pointer() const {
      return graph_ptr_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Pointer to the true edge.
    // NOTE: Edge's index is enforced to be same as edge's position in edge_list_.
    EdgeListIndex edge_list_index_;

    // Pointer to the graph to which this node belongs.
    Graph *graph_ptr_;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_list_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < edge_list_.size());
    return Edge(i, this);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    Edge edge = find_connecting_edge(a, b);
    return edge.get_edge_index() < edge_list_.size();
  }

  // Helper function to obtain the connecting edge.
  // Improves efficiency of add_edge() fxn. Otherwise would have to potentially
  // iterate through adjacent_edges twice.
  Edge find_connecting_edge(const Node &a, const Node &b) const {

    const NodeInternal &node = node_list_[a.get_node_index()];

    for (const size_type edge_index : node.adjacent_edges) {
      const EdgeInternal &edge = edge_list_[edge_index];
      if (edge.node1 == b || edge.node2 == b) {
        return Edge(edge_index, this);
      }
    }
    return Edge();
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
    // HW0: YOUR CODE HERE

    Edge edge_proxy = find_connecting_edge(a, b);

    if (edge_proxy.get_edge_index() < edge_list_.size()) {
      return edge_proxy;
    } else {
      // create index for this new edge
      const size_type new_index = edge_list_.size();

      // add edge to graph
      edge_list_.emplace_back(new_index, a, b);

      // NOTE: do not specify these as const b/c need to add edges
      NodeInternal &node_a = node_list_[a.get_node_index()];
      NodeInternal &node_b = node_list_[b.get_node_index()];

      // add edge to nodes a and b
      (node_a.adjacent_edges).insert(new_index);
      (node_b.adjacent_edges).insert(new_index);

      // create edge proxy
      Edge edge_proxy = Edge(new_index, this);

      return edge_proxy;
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    node_list_.clear();
    edge_list_.clear();
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // Vector containing true nodes.
  using NodeList = std::vector<NodeInternal>;

  // Vector containing true edges.
  using EdgeList = std::vector<EdgeInternal>;

  // Stores the indices of the edges adjacent to the node.
  using AdjacentEdges = std::set<size_type>;

  struct NodeInternal {
    size_type index; // NOTE: this index == position of node in node_list_
    Point position;
    AdjacentEdges adjacent_edges {};
    NodeInternal(const size_type index_, const Point &position_)
      : index(index_), position(position_) {
    }
  };

  struct EdgeInternal {
    size_type index; // NOTE: this index == position of edge in edge_list_
    Node node1; // NOTE: storing node proxies in true edge
    Node node2;
    EdgeInternal(const size_type index_, const Node &node1_,
      const Node &node2_) : index(index_), node1(node1_), node2(node2_) {
    }
  };

  // Vector containing the true nodes.
  NodeList node_list_ {};

  // Vector containing the true edges.
  EdgeList edge_list_ {};

};

#endif // CME212_GRAPH_HPP
