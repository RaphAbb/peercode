#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
//#include <iostream>

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

  // Node Point data (index is node id)
  std::vector<Point> Points;
  // Edge Node id's data (index is edge id)
  std::vector<unsigned> EdgeA_ids;
  std::vector<unsigned> EdgeB_ids;
  // adjacency list: Node id -> connected Node id's
  std::vector<std::vector<unsigned>> Nadj_list;
  // adjacency list: Node id -> corresponding Edge id's to the Node id's above
  std::vector<std::vector<unsigned>> Eadj_list;
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
      return graph_->Points[id_];
    }

    /** Return this node's index, a number in the range [0, num_Nodes). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return id_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if(graph_ == n.graph_)
        if(id_ == n.index())
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
      // HW0: YOUR CODE HERE
      if(id_ < n.index())
        return true;
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // internal Node variables
    Graph* graph_;
    size_type id_;
    /** Private Constructor */
    Node(const Graph* g_in, size_type id_in) : graph_(const_cast<Graph*>(g_in)), id_(id_in) {}
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return num_Nodes;
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
    // updates node information
    Points.push_back(position);

    // increases size of adjacency lists
    std::vector<unsigned> newvec;
    Nadj_list.push_back(newvec);
    Eadj_list.push_back(newvec);
    
    num_Nodes++;
    return Node(this,num_Nodes-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    // also makes sure if Node id is not too large
    if(n.graph_ == this)
      if(n.id_ < size())
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
    // HW0: YOUR CODE HERE
    assert(i < size());
    return Node(this,i);
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

    /** Return this edge's index, a number in the range [0, num_Edges). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return id_;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      // retrieves Node id
      size_type A_i = graph_->EdgeA_ids[id_];
      return graph_->node(A_i);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      // retrieves Node id
      size_type B_i = graph_->EdgeB_ids[id_];
      return graph_->node(B_i);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // compares node id's in both directions
      if(graph_->EdgeA_ids[id_] == graph_->EdgeA_ids[e.index()])
        if(graph_->EdgeB_ids[id_] == graph_->EdgeB_ids[e.index()])
          return true;
      if(graph_->EdgeA_ids[id_] == graph_->EdgeB_ids[e.index()])
        if(graph_->EdgeB_ids[id_] == graph_->EdgeA_ids[e.index()])
          return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if(id_ < e.index())
        return true;
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // internal Edge variables
    Graph* graph_;
    size_type id_;
    /** Private Constructor */
    Edge(const Graph* g_in, size_type id_in) : graph_(const_cast<Graph*>(g_in)), id_(id_in) {}
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return num_Edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i<num_edges());
    return Edge(this,i);        // Valid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    // Uses the node adjacency list to determine if edge already exists
    for(size_type i = 0; i < Nadj_list[a.index()].size(); i++)
    {
      if(Nadj_list[a.index()][i] == b.index())
      {
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
    // HW0: YOUR CODE HERE
    // Uses the node adjacency list to determine if edge exists
    for(size_type i = 0; i < Nadj_list[a.index()].size(); i++)
    {
      if(Nadj_list[a.index()][i] == b.index())
      {
        // determines the edge id corresponding to that position on the node adjacency list
        size_type E_ind = Eadj_list[a.index()][i];
        // returns existing edge
        if(EdgeA_ids[E_ind] == a.index())
          if(EdgeB_ids[E_ind] == b.index())
          {
            return Edge(this,E_ind);
          }
        // switches a,b and returns existing edge
        if(EdgeB_ids[E_ind] == a.index())
          if(EdgeA_ids[E_ind] == b.index())
          {
            EdgeA_ids[E_ind] = a.index();
            EdgeB_ids[E_ind] = b.index();
            return Edge(this,E_ind);
          } 
      }
    }

    // else, creates and returns new edge
    // updates adjacency lists
    Nadj_list[a.index()].push_back(b.index());
    Nadj_list[b.index()].push_back(a.index());
    Eadj_list[a.index()].push_back(num_Edges);
    Eadj_list[b.index()].push_back(num_Edges);
    // adds edge id data
    EdgeA_ids.push_back(a.index());
    EdgeB_ids.push_back(b.index());
    num_Edges++;
    return Edge(this,num_Edges-1);  
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    Points.clear();
    EdgeA_ids.clear();
    EdgeB_ids.clear();
    Nadj_list.clear();
    Eadj_list.clear();
    num_Nodes = 0;
    num_Edges = 0;
  }

 private:

  // HW0: YOUR CODE HERE
  size_type num_Nodes = 0;
  size_type num_Edges = 0;
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
