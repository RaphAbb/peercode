#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <iostream>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
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
  //for an empty graphy, we want the nodes_ and adj_ vectors to be empty
  Graph() 
      : nodes_(std::vector<Point*>(0)),
	adj_(std::vector<std::vector<size_type>>(0)) {
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

    //invalid node
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      /** Dereferences the Point* in the nodes_ vector belong to Graph object
       * at the index associated with the Node object
       */
      return *((*graph_).nodes_[node_idx_]);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return node_idx_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      //return true if nodes belong to same graph and have same index
      if (node_idx_ == n.node_idx_ and graph_ == n.graph_) return true;
      (void) n;          // Quiet compiler warning
      //return false otherwise
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
      /** Return true if the index of this node is less than that of n
       * and they belong to the same graph
       */
      if (node_idx_ < n.node_idx_ and graph_ == n.graph_) return true;
      (void) n;           // Quiet compiler warning
      //return false otherwise
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    //pointer to the location of the Graph object of which Node is a node
    Graph* graph_;
    //current index of the node in the node_ vector 
    size_type node_idx_;

    //valid Node constructor
    Node(const Graph* g, size_type idx) 
	: graph_(const_cast<Graph*>(g)), node_idx_(idx) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    //return the size of the vector holding all nodes
    return nodes_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size(); //reuse the size() member function
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position) {
    //dynamically allocated a new Point on the heap
    Point *p = new Point;
    //dereference and set equal to position argument
    *p = position;
    nodes_.push_back(p); //add pointer to p to nodes vector
    //update length of adj list; node is not connected any any other nodes yet
    adj_.push_back(std::vector<size_type>(0));
    return Node(this, num_nodes()-1);
    (void) position;      // Quiet compiler warning
    return Node();
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    //if n belongs to the same graph and has index less than # of graph nodes
    if (n.node_idx_ < num_nodes() and n.graph_ == this) return true;
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
    //return a Node constructor for this graph object and specified point 
    return Node(this, i);
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
    //invalid Edge constructor
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      //return a Node object constructed from graph and node 1 index
      return Node(graph_, node1_idx_);
      return Node();	  //Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      //return a Node object constructed from graph and node 2 index
      return Node(graph_, node2_idx_);
      return Node();      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //if edges in the same graph and node indices match; order not considered 
      if (graph_ == e.graph_) {
        if (node1_idx_ == e.node1_idx_ and node2_idx_ == e.node2_idx_) return true;
	if (node1_idx_ == e.node2_idx_ and node2_idx_ == e.node1_idx_) return true;
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
      //we will define Edge e1 < Edge e2 if e1.node1_idx_ < e2.node1_idx_
      if (graph_ == e.graph_ and node1_idx_ < e.node1_idx_) return true;
      (void) e;           // Quiet compiler warning
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    //pointer to Graph object to which the edge belongs
    Graph* graph_;
    //index of the first node to which the edge is adjacent
    //in terms of the adjacency list, the edge will beloing to adj[node_1_idx_] 
    size_type node1_idx_;
    //index of the second node to which the edge is adjacent
    size_type node2_idx_;

    //valid Edge contructor
    Edge(const Graph* g, size_type n1_idx, size_type n2_idx)
	: graph_(const_cast<Graph*>(g)),
	  node1_idx_(n1_idx),
	  node2_idx_(n2_idx) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    //initialize number of edges variable to zero
    size_type num_e = 0;
    //loop through adjacency list and add the length of each vec in adj vec
    for (auto v : adj_) {
      num_e += v.size();
    }
    //adjacency list encodes each edge twice, so divide by 2
    return num_e/2;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    /** adj_ is set up as an adjacency list, meaning edge (u,v) appears as 
     * v in adj_[u] and u in adj[v]. So edges are repeated.
     *
     * First we obtain the index of the ith unique element if all vectors in
     * were to be concatenated. Call this index k.
     */
    //counter to keep track of unique elements
    size_type counter = i+1;
    //initialize index k to zero
    size_type k = 0;
    //loop through all vectors
    for (size_type idx1 = 0; idx1 < adj_.size(); ++idx1) {
      //first vector will contain all unique edges
      if (idx1 == 0) {
        if (i < adj_[idx1].size()) {k = i; break;}
        else {k += adj_[idx1].size(); counter -= adj_[idx1].size();}
      }
      else {
	for (size_type idx2 = 0; idx2 < adj_[idx1].size(); ++idx2) {
	  //edge is already counted if the node 2 index is less than node 1
	  if (adj_[idx1][idx2] < idx1) {k += 1;}
          else {
            counter -= 1;
	    //when the counter reaches zero we have the element we want
            if (counter == 0) goto out_of_nested_loops;
	    k += 1;
          }
	} 
      }
    }
    //use goto to break out of nested loops
    out_of_nested_loops:
    //now extract the kth element of "concatenated" vectors
    size_type elems = 0; //sum of edges in each vector in adj_
    size_type v_idx = -1; //idx of vector whose size was just added
    //loop through all vectors in adj_
    for (auto v : adj_) {
      elems += v.size(); //updates elems
      v_idx += 1; //update v_idx
      //if the edge # we wanted is less than elem sum, it was in most recent vec 
      if (k < elems) break;
    }
    //calculate index of ith element within the inner adj_ vector
    size_type j = k-(elems-adj_[v_idx].size());
    //create resulting edge and return
    return Edge(this,v_idx,adj_[v_idx][j]);
    (void) i;             // Quiet compiler warning
    return Edge();        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    //if nodes belong to different graphs return false
    if (a.graph_ != b.graph_) return false;
    //isolate vector corresponding to index of node a
    std::vector<size_type> v = adj_[a.node_idx_];
    //if index of node b appears in this vector, edge exists
    if (std::find(v.begin(),v.end(),b.node_idx_) != v.end()) return true; 
    (void) a; (void) b;   // Quiet compiler warning
    //otherwise edge does not exist
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
    size_type a_idx = a.node_idx_;
    size_type b_idx = b.node_idx_;
    //if the edge already exists, return current edge
    std::vector<size_type> v = adj_[a_idx];
    //std::vector<size_type> w = adj_[b_idx];
    if (std::find(v.begin(),v.end(),b_idx) != v.end()) 
      return Edge(this,a_idx,b_idx);
    //else if (std::find(w.begin(),w.end(),a_idx) != w.end()) 
      //return Edge(this,b_idx,a_idx);
    //otherwise add edge
    else {
      adj_[a_idx].push_back(b_idx);
      adj_[b_idx].push_back(a_idx);
    }
    //return new edge object
    return Edge(this,a_idx,b_idx);
    (void) a, (void) b;   // Quiet compiler warning
    return Edge();        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    //make use of vector clear() function
    nodes_.clear();
    adj_.clear();
  }

 private:
   /** Vector of the pointers to the Point object associated with all nodes */
   std::vector<Point*> nodes_;

   /** Vector of vectors which contain the edge connections for each node.
    * Outer vector index corresponds to node index, and inner vectors contain
    * size_type elements, each corresponding to the node to which it is connected
    */
   std::vector<std::vector<size_type>> adj_;
};

#endif // CME212_GRAPH_HPP
