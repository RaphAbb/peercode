#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <tuple>
#include <vector>
#include <unordered_map>
#include <stdexcept>
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

// Publicly declare Node class for its use in the edgevec_ below
public:
  class Node;

 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  // Store the nodes using a vector of the positions of the nodes in the graph,
  // ordered by nodeID.
  std::vector<Point> nodevec_;

  // Store the edges with an adjacency list. adjlist_ maps nodeIDs to maps of
  // nodeIDs to edgeIDs. Thus if the edge with edgeID 4 connects the nodes with
  // nodeIDs 1 and 2, adjlist_ will contain the entries 1 -> (2 -> 4) and 2 -> (1 -> 4).
  // This allows O(1) lookup for both the presence AND the ID of an edge given two nodes, the
  // first of which is useful in has_edge() and the second of which is useful in add_edge().
  std::unordered_map<unsigned int, std::unordered_map<unsigned int,unsigned int>> adjlist_;

  // Also keep a vector of the edges in the graph, ordered by edgeID. This doesn't explicitly
  // contain any information different from the adjacency list, but it does allow for O(1)
  // lookup of an edge given an edgeID. This is useful in the edge() method. We could choose
  // to not use this vector and instead iterate over the adjacency list in search of a given
  // edgeID, but then the edge() method would not support constant time lookup. There is a
  // tradeoff between time and space here, and as no guidelines for managing this tradeoff
  // were given, I have chosen to use more space for less time.
  std::vector<std::tuple<unsigned int , unsigned int>> edgevec_;

  // Track the number of nodes and edges currently in the Graph.
  unsigned int numnodes_;
  unsigned int numedges_;

public:

 //
 // PUBLIC TYPE DEFINITIONS
 //

 /** Type of this graph. */
 using graph_type = Graph;

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
  Graph():
  nodevec_(), adjlist_(), edgevec_(), numnodes_(0), numedges_(0)
    // HW0: YOUR CODE HERE
  {}

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
      return this->graph_->nodevec_[nodeid_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return this->nodeid_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (n.graph_ == this->graph_ && n.index() == this->index()){return true;}
      else{return false;}
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
      // Base the decision on node IDs.
      // As stated on Piazza, we can assume the nodes are in the same Graph
        if (this->nodeid_ >= n.index()){return false;}
        else {return true;}
      }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    Graph* graph_; // Pointer back to graph
    size_type nodeid_; // Unique node id
    // Private constructor
    Node(const Graph* graph, size_type id):
    graph_(const_cast<Graph*>(graph)), nodeid_(id) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return numnodes_;
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
    Node new_node(this, numnodes_); // Create new node
    nodevec_.push_back(position); // Add position of new node to the nodevec
    numnodes_ += 1; // Increment size of graph
    return new_node;        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (n.graph_ == this && numnodes_ > n.index()){return true;}
    else{return false;}
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i >= numnodes_){
      // Out of bounds; throw error
    throw std::invalid_argument("Index out of range");
    }
    else{
    Node nodei = Node(this, i);
    return nodei;}       // Invalid node
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
      return Node(this->graph_, this->node1id_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(this->graph_, this->node2id_);
    }

    /** Return this edge's index, a number in the range [0, numedges_). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return this->edgeid_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */

    bool operator==(const Edge& e) const {
      if (e.graph_ == this->graph_){ // Verify the Graphs are the same
        // Verify that the pairs of nodes are identical
        if ((e.node1() == this->node1()) && (e.node2() == this->node2())){
          return true;
        }
        else if ((e.node1() == this->node2()) && (e.node2() == this->node1())){
          return true;
        }
        // Edge is between different pairs of nodes
        else{return false;}
      }
      // Edges are in different graphs
      else{return false;}
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */

    bool operator<(const Edge& e) const {
      // Again, can assume the edges are in the same graph, according to Piazza
      // Thus, we can just base this on edgeID.
      if (this->edgeid_ >= e.index()){return false;}
      else {return true;}
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_; // Pointer back to graph
    size_type edgeid_; // Unique edge id, the position in the edgevec_ vector
    size_type node1id_; // ID of one of the endpoints
    size_type node2id_; // ID of the other endpoint
    // Private constructor
    Edge(const Graph* graph, size_type id, size_type node1id, size_type node2id):
    graph_(const_cast<Graph*>(graph)), edgeid_(id), node1id_(node1id), node2id_(node2id) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return this->numedges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i >= numedges_){
      // Out of bounds; throw error
    throw std::invalid_argument("Index out of range");
    }
    else{
      // Return proxy. edgevec_ allows for O(1) lookup of the nodeIDs for this edgeIDs.
      size_type node1id = std::get<0>(this->edgevec_[i]);
      size_type node2id = std::get<1>(this->edgevec_[i]);
    Edge edgei = Edge(this, i, node1id, node2id);
    return edgei;}       // Invalid node
  }


  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const{
    // HW0: YOUR CODE HERE
    unsigned int ind1 = a.index();
    unsigned int ind2 = b.index();
    // test for (a,b)
    if (adjlist_.count(ind1) > 0){
      // a has neighbors; check for b as one of them
      if (adjlist_.at(ind1).count(ind2) > 0){
        return true; // b is a neighbor of a
      }
      else{return false;} // b is not a neighbor of a
    }
      else{return false;} // a has no neighbors at all
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
    if (has_edge(a,b)){
    // Return a proxy to the existing edge. First we need its index
    unsigned int indx = adjlist_[a.index()][b.index()];
    return Edge(this, indx, a.index(), b.index());
    }
    else{
    // Now, the case where the edge was not already present in Graph
    Edge new_edge(this, numedges_, a.index(), b.index()); // Create new edge
    // Update edgevec_
    std::tuple<unsigned int, unsigned int> given_edge_fwd(a.index(),b.index());
    edgevec_.push_back(given_edge_fwd); // Add node pair to the end of edgevec_
    // Update adjlist_
    adjlist_[a.index()][b.index()] = numedges_; // numedges_ is the new edgeID
    adjlist_[b.index()][a.index()] = numedges_; // numedges_ is the new edgeID
    numedges_ += 1; // Increment size of edge set
    return new_edge;
  }
}

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    this->numnodes_ = 0;
    this->numedges_ = 0;
    this->nodevec_.clear();
    this->edgevec_.clear();
    this->adjlist_.clear();
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
