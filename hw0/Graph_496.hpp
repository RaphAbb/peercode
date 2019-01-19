#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <tuple>
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
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->points[this-> idx_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return this-> idx_;
      return size_type(-1);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if ((n.graph_ == this->graph_) and (n.idx_ == this->idx_)) {
      	return true;
      }

      (void) n;          // Quiet compiler warning
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
      // HW0: YOUR CODE HERE
      if (this-> idx_ < n.idx_) {
      	return true;
      }
      (void) n;           // Quiet compiler warning
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Attribute to keep track of graph of node
    graph_type* graph_;

    // Attribute to keep track of index of node
    size_type idx_;

    // Constructor for Node
    Node(graph_type* g, size_type id)
      : graph_(g), idx_(id) {
      }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return this->next_idx;
    return 0;
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

    // Add node to the Graph
    this->points.push_back(position);

    // Increment next_idx
    this->next_idx++;

    // Add empty vector to edges
    connections.push_back(std::vector<edge_type>());

    // Create and return node instance
    return Node(this, next_idx - 1);

    (void) position;      // Quiet compiler warning
    return Node();        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return this == n.graph_;
    (void) n;            // Quiet compiler warning
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE

    // Test for precondition
    // Note: no need to check 0 <= i, as 0 <= size_type always true
    if (i < this->num_nodes()) {
      return Node(const_cast<graph_type*>(this), i);
    }

    // If precondition fails, return invalid node
    (void) i;             // Quiet compiler warning
    return Node();        // Invalid node
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      
      // Index into edges data structure
      return std::get<0>(graph_->edges.at(this->idy_));
      return Node();      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE

      // Index into edges data structure
      return std::get<1>(graph_->edges.at(this->idy_));
      return Node();      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if ((e.graph_ == this->graph_) and (e.idy_ == this->idy_)) {
        return true;
      }

      (void) e;           // Quiet compiler warning
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (this-> idy_ < e.idy_) {
        return true;
      }
      (void) e;           // Quiet compiler warning
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Attribute to keep track of graph of edge
    graph_type* graph_;

    // Attribute to keep track of index of edge
    size_type idy_;

    // Attribute to keep track of nodes of edge
    node_type node1_;
    node_type node2_;

    // Constructor for Edge
    Edge(graph_type* g, size_type id, node_type n1, node_type n2)
      : graph_(g), idy_(id), node1_(n1), node2_(n2) {
      }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return next_idy;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE

    // Check for precondition
    // Note: no need to check 0 <= i, as 0 <= size_type always true
    if (i < num_edges()) {

      return Edge(const_cast<graph_type*>(this), i, 
        std::get<0>(edges.at(i)), std::get<1>(edges.at(i)));
    }
    
    // If precondition fails, return an invalid edge
    (void) i;             // Quiet compiler warning
    return Edge();        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully lessls
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE

    // Test for precondition
    if (has_node(a) and has_node(b)) {
      for (auto edge : connections.at(a.index())) {
        
        // Check forward version of edge
        if ((edge.node1() == a) and (edge.node2() ==b)) {
          return true;
        }
        // Check backward version of edge
        if ((edge.node1() == b) and (edge.node2() == a)) {
          return true;
        }
      }
    }
    // If precondition fails, throw an exception
    else {
      throw std::invalid_argument( "Invalid nodes received" );
    }
    
    // Otherwise, if edge not in graph, return false
    (void) a; (void) b;   // Quiet compiler warning
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
    // HW0: YOUR CODE HERE

    //Test for precondition 
    if (a == b) {
      throw std::runtime_error("Self loop illegal");
    }
    if (!has_node(a) or !has_node(b)) {
      throw std::runtime_error("Invalid nodes for edge addition");
    }

    // If edge present in graph
    if (has_edge(a,b)) {

      // If input ordered edge is present in graph, return edge object
      for (auto edge : connections.at(a.index())) {
        
        // Check forward version of edge
        if ((edge.node1() == a) and (edge.node2() ==b)) {
          return Edge(this, edge.idy_, a, b);
        }
        // Check backward version of edge
        if ((edge.node1() == b) and (edge.node2() == a)) {
          std::tuple<node_type, node_type> new_edge = std::make_tuple(b, a);
          edges.at(edge.idy_) = new_edge;
          return Edge(this, edge.idy_, b, a);
        }
      }
    }

    // Otherwise, add edge to graph, return edge object
    std::tuple<node_type, node_type> new_tup = std::make_tuple(a, b);
    this->edges.push_back(new_tup);

    edge_type new_edge =Edge(this, next_idy, a, b);

    // Add edge to our adjacency list
    (connections.at(a.index())).push_back(new_edge);
    (connections.at(b.index())).push_back(new_edge);

    this->next_idy++;

    return Edge(this, next_idy - 1, a, b);
    

    (void) a, (void) b;   // Quiet compiler warning
    return Edge();        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE

    // Delete all nodes stored in graph
    points.clear();

    // Reset next_idx (num_nodes) to 0
    next_idx = 0;

    // Delete all edges stored in adjacency list
    connections.clear();

    // Delete all edges stored in vector
    edges.clear();

    // Reset next_idy (num_edges) to 0
    next_idy = 0;
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // Vector of Points representing our Graphs nodes
 	std::vector<Point> points;

  // Index for points
  size_type next_idx = 0;

  // Vector of vectors representing an adjacency list
  std::vector<std::vector<edge_type>> connections;

  // Vector of tuples of node objects representing edges
  std::vector<std::tuple<node_type, node_type>> edges;

  // Index for edges
  size_type next_idy = 0;

};

#endif // CME212_GRAPH_HPP