#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <set>
#include <map>

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

  struct pseudoedge;

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
  Graph() : numnodes_(0), numedges_(0){
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
      //std::cout << "Node::position() begin" << std::endl;
      return (this->pointer_to_graph->vecpoints_[unid_]);
    }

    /*Need to finish*/ // can I check if index is always the same as unid?
    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      //std::cout << "Node::index() begin" << std::endl;
      assert(this->unid_ < this->pointer_to_graph->num_nodes());
      //std::cout << "Node::index() end" << std::endl;
      return this->unid_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      //std::cout << "Node::operator==() begin" << std::endl;
      if ((this->unid_ == n.unid_) & (this->pointer_to_graph == n.pointer_to_graph)){
        return true;
      }
      //(void) n;          // Quiet compiler warning
      //std::cout << "Node::operator==() end" << std::endl;
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
      //std::cout << "Node::operator <() begin" << std::endl;
      if ((this->unid_ < n.unid_) & (this->pointer_to_graph == n.pointer_to_graph)){
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* pointer_to_graph;
    size_type unid_;

    // private Node Constructor
    Node(const Graph* p_t_g, size_type unid_in)
      : pointer_to_graph(const_cast<Graph*>(p_t_g)), unid_(unid_in){
      }

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
  Node add_node(const Point& position) { /*Need to finish*/ //Should I check if the node already exists here?
    // HW0: YOUR CODE HERE
    //std::cout << "Graph::add_node() begin" << std::endl;
    this->vecpoints_.push_back(position);
    ++numnodes_;

    // Push back in the edges adjacency list
    this->vecedges_.push_back(std::vector<size_type>());
    return Node(this, numnodes_-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const { /*Need to finish*/
    // HW0: YOUR CODE HERE
    //std::cout << "Graph::has_node() begin" << std::endl;
    if ((n.pointer_to_graph == this) && (n.unid_ < this->numnodes_)){
      return true;
    }
    //std::cout << "Graph::has_node() end" << std::endl;
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const { /*Need to finish*/ // Should I return Point here because that's what the user expects?
    // HW0: YOUR CODE HERE
    assert(0 <= i);
    assert(i < numnodes_);

    return Node(this, i);        // Invalid node
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
      return Node(p_t_g_edge_, node1_ID_);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(p_t_g_edge_, node2_ID_);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //std::cout << "Edge::operator==() begin" << std::endl;
      if (((e.node1_ID_ == this->node1_ID_) && (e.node2_ID_ == this->node2_ID_)) ||
          ((e.node1_ID_ == this->node2_ID_) && (e.node2_ID_ == this->node1_ID_))){
            return true;
          }
      //std::cout << "Edge::operator==() end: returning false" << std::endl;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (this->edge_ID_ < e.edge_ID_){
        //std::cout << "Edge::operator <() end: returning false" << std::endl;
        return true;
      }
      //std::cout << "Edge::operator <() end: returning false" << std::endl;
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    Graph* p_t_g_edge_;
    size_type edge_ID_;
    size_type node1_ID_;
    size_type node2_ID_;

    // private Edge constructor
    Edge(const Graph* p_t_g_in, size_type edge_ID_in, size_type nodeval1, size_type nodeval2)
      : p_t_g_edge_(const_cast<Graph*>(p_t_g_in)), edge_ID_(edge_ID_in), node1_ID_(nodeval1), node2_ID_(nodeval2){
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return numedges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {/*Need to finish*/
    // HW0: YOUR CODE HERE
    assert(0 <= i);
    assert(i < num_edges());
    // loop over the Key values of mapedges_ and find the match for i
    std::pair<size_type, size_type> temppair = this->edge_index_[i];

    return Edge(this, i, temppair.first, temppair.second);
    //for(auto p: v)
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    //std::cout << "Graph::has_edge() begin" << std::endl;
    // Check if @a a and @a b are valid nodes of this graph
    assert(this == a.pointer_to_graph);
    assert(this == b.pointer_to_graph);
    assert(a.unid_ < numnodes_);
    assert(b.unid_ < numnodes_);

    // Search the vector of connected nodes to the first node passed in to see if an edge exists there.
    std::vector<size_type> tempvec = this->vecedges_[a.unid_];
    for (size_type i = 0; i < tempvec.size(); i++){
      if (tempvec[i] == b.unid_){
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
  Edge add_edge(const Node& a, const Node& b) { /*Need to finish*/ // What is current edge??
    // HW0: YOUR CODE HERE
    //assert(a and b are valid nodes.);
    assert(this == a.pointer_to_graph);
    assert(this == b.pointer_to_graph);
    assert(a.index() < this->numnodes_);
    assert(b.index() < this->numnodes_); // /*Need to finish*/ can I improve this?
    //assert that a and b are distinct
    assert(a.index() != b.index());

    // Only add the edge if the edge does not already exist in the graph
    bool edgeexists = has_edge(a,b);
    if (!edgeexists){
      this->vecedges_[a.unid_].push_back(b.unid_);
      this->vecedges_[b.unid_].push_back(a.unid_);
      this->edge_index_.push_back(std::pair<size_type,size_type>(a.unid_,b.unid_));
      numedges_++;
      return Edge(this, numedges_-1, a.unid_, b.unid_);
    }
    else{
      // edge exists already --> Return the existing edge index.
      // need to search through this->edge_index_ for it.
      //this->edge_index;
      //std::make_pair(a.unid_,b.unid_);
      std::pair<size_type, size_type> temppair(a.unid_, b.unid_);
      std::pair<size_type, size_type> temppair2(b.unid_, a.unid_);
      for (size_type i = 0; i < this->edge_index_.size(); i++){
        if ((this->edge_index_[i] == temppair) || (this->edge_index_[i] == temppair2)){
          return Edge(this, i, a.unid_, b.unid_);
        }
      }
      std::cout << "add_edge bad option happened" << std::endl;
      return Edge(this, numedges_-1, a.unid_, b.unid_); // this numedges_-1?
    }
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    // remove all nodes from the graph
    vecpoints_.clear();
    assert(vecpoints_.size() == 0);
    numnodes_ = 0;

    // remove all edges from the graph
    vecedges_.clear();
    edge_index_.clear();
    numedges_ = 0;
    assert(vecedges_.size()==0);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // Nodes:
  std::vector<Point> vecpoints_;
  size_type numnodes_;

  // Edges:
  std::vector<std::vector<size_type>> vecedges_;
  std::vector<std::pair<size_type,size_type>> edge_index_;
  size_type numedges_;


// Search to find all remaining parts to finish: /*Need to finish*/


};

#endif // CME212_GRAPH_HPP
