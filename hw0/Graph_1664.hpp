#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

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
      // Return Point position using index and Graph pointer
      return ((this->graph_)->listOfNodes[this->index_].position_);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // Return this Node's index
      return this->index_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // Return whether Nodes have the same Graph and the same index.
      bool sameGraph = ( this->graph_ == n.graph_ ); // Check if Nodes have same Graph pointer
      bool sameIndex = ( this->index() == n.index() ); // Check if same index
      return (sameIndex && sameGraph);
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
      // Define order based on Node's index
      return ( this->index() < n.index() );
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Private Variables tracked by Node
    graph_type* graph_; // Pointer back to Graph containing Node
    size_type index_; // Index value associated with Node
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // Return number of Nodes
    return size_type(this->listOfNodes.size());
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
    // Create a Node with pointer to Graph and index
    Node n = Node(); // Create (invalid) Node inside nodeInfo struct
    n.graph_ = this; // Assign this Graph to Node
    n.index_ = size_type(this->listOfNodes.size()); // Give Index to Node
    // Create a nodeInfo struct with Node and associated information
    nodeInfo n_info; // Create Empty nodeInfo struct
    n_info.position_ = position; // Give position to nodeInfo struct
    n_info.index_ = size_type(this->listOfNodes.size()); // Give Index to nodeInfo struct
    n_info.node_ = n; // Put Node inside nodeInfo struct
    // Add Node (with correct index) to Graph
    this->listOfNodes.push_back(n_info); // Assign index to Node
    return n; // Return Node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // Make sure that Node is inside Graph at index indicated by Node
    if (n.index_ < this->num_nodes()) { // Make sure index of Node is inside Graph
      return (this->listOfNodes[n.index_].node_ == n); // Match Nodes
    } else {
      return false;
    }
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // Return Node with index i
    return this->listOfNodes[i].node_;
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
      // Return reference to dereferenced pointer
      return (*((this->graph_)->listOfEdges[this->index_].node1_));
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // Return reference to dereferenced pointer
      return (*((this->graph_)->listOfEdges[this->index_].node2_));
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // Equal Edges represent the same undirected edge between two Nodes
      bool nodesMatchExactly = ((e.node1() == this->node1()) && (e.node2() == this->node2()));
      bool nodesMatchSwitched = ((e.node1() == this->node2()) && (e.node2() == this->node1()));
      return (nodesMatchExactly || nodesMatchSwitched);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // Inequality based on Edge index
      return (this->index_ < e.index_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    graph_type* graph_; // Pointer back to Graph containing Edge
    size_type index_; // Index value associated with Node
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // Return number of Edges in Graph
    return size_type(this->listOfEdges.size());
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // Return Edge from Graph
    return (this->listOfEdges[i]).edge_;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Iterate through the vector of Edges
    for(size_type i=0; i<this->num_edges(); i++) {
      edge_type e = this->edge(i);
      if ( ((e.node1()==a) && (e.node2()==b)) || ((e.node1()==b) && (e.node2()==a)) ) {
        return true; // Return true if this Edge connects Nodes
      }
    }
    return false; // Return false if no Edge connects those two Nodes
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
    // Create Edge with two associated Nodes, and add to Graph
    if (this->has_edge(a,b)) { // If Graph has Edge that connects these Nodes
      for(size_type i=0; i<this->num_edges(); i++) {
        edge_type e = this->edge(i);
        if ( ((e.node1()==a) && (e.node2()==b)) || ((e.node1()==b) && (e.node2()==a)) ) {
          return e; // Return Edge that connects Nodes
        }
      }
      assert(false); // Execution should not reach this point! If it does something is wrong
      return Edge(); // Gets rid of warning
    } else {
      // Create Edge
      edge_type e = Edge(); // Create (invalid) Edge
      e.graph_ = this; // Assign this Graph to Edge
      e.index_ = size_type(this->listOfEdges.size()); // Record index of Edge
      // Add Edge (and associated Nodes) to edgeInfo
      edgeInfo e_info; // Create edgeInfo struct
      e_info.index_ = e.index_; // Record index of Edge into edgeInfo struct
      e_info.node1_ = &a; // First Node associated with Edge
      e_info.node2_ = &b; // Second Node associated with Edge
      e_info.edge_ = e; // Record Edge into edgeInfo struct
      this->listOfEdges.push_back(e_info); // Add Edge and edgeInfo to Graph
      return e; // Return edge
    }

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // Remove all Nodes
    this->listOfNodes.clear();
    // Remove all Edges
    this->listOfEdges.clear();
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // Structs for Holding Node and Edge information
  struct nodeInfo {
    Point position_; // position as a Point object
    size_type index_; // Index value associated with Node
    node_type node_; // Node
  };
  struct edgeInfo {
    size_type index_; // Index value associated with Node
    const node_type* node1_; // Pointer to first Node
    const node_type* node2_; // Pointer to second Node
    edge_type edge_; // Edge
  };

  // Node information within Graph
  std::vector<nodeInfo> listOfNodes; // Vector of nodeInfo structs

  // Edge information within graph
  std::vector<edgeInfo> listOfEdges; // Vector of edgeInfo structs

};

#endif // CME212_GRAPH_HPP
