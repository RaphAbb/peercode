#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <map>
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
 struct internal_node;
 struct internal_edge;
 
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
      return graph_->nodes_[index_].position; //access position via proxy vector of internal_node's
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return index_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return ((graph_ == n.graph_)&&(index_ == n.index_));
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
      return (index_ < n.index_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    // Pointer back to the Graph
    Graph* graph_;
    
    //index of node in graph.nodes
    size_type index_;

    /** Private Constructor */
    Node(const Graph* graph, size_type index):graph_(const_cast<Graph*>(graph)), index_(index) {};

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes_.size(); //accessing member of vector is O(1)
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
    size_type index = num_nodes();
    nodes_.push_back({position, index}); //add internal_node proxy; inserting is O(1) amortized
    adjacency_.push_back(std::map<size_type, size_type>()); //map in adjacency vector will store <neighbor, edge> pairs
    return Node(this, index);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.index() >= num_nodes())  
      return false; //invalid index for graph
    else {
      internal_node n2 = nodes_[n.index()]; //internal_node proxy
      return ((n2.position==n.position())&&(n2.index==n.index())); 
    }
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i >= 0);
    assert(i < num_nodes());
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
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //equal edges should be in the same graph and connect the same nodes (regardless)
      //of the order of the connection They will also have the same index since each
      //edge is added only once.
      return ((graph_ == e.graph_)&&(((node1_ == e.node1_)&&(node2_ == e.node2_)) ||
          ((node1_ == e.node2_)&&(node2_ == e.node1_))));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (index_ < e.index_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objectsi
    
    // Pointer back to the Graph
    Graph* graph_;

    // index of edge in graph.edges_
    size_type index_;

    // index of node1 in graph.nodes_
    size_type node1_;

    // index of node2 in graph.nodes_
    size_type node2_;

    /** Private Constructor */
    Edge(const Graph* graph, size_type index, size_type node1, size_type node2):graph_(const_cast<Graph*>(graph)), index_(index), node1_(node1), node2_(node2) {};
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_.size(); //O(1) to get size of vector

  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i >= 0);
    assert(i < num_edges());
    return Edge(this, i, (edges_[i]).node1, (edges_[i]).node2); //determine edge's nodes via proxy vector of internal_edge's; O(1) for vector operator []
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    std::map<size_type, size_type> a_edges = adjacency_[a.index()]; //maps neighbors of a to index of connecting edge in proxy edges_
    return (a_edges.find(b.index()) != a_edges.end()); //edge exists iff b in neighbors of a
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
    size_type i; //index of connecting edge in proxy edges_
    if (has_edge(a, b))
      i = (adjacency_[a.index()])[b.index()];  
    else {
      i = num_edges(); 
      (adjacency_[a.index()]).insert(std::pair<size_type, size_type>(b.index(), i)); //add that b is a neighbor of a, connected by edge i; O(num_edges()) to insert into map of size O(num_edges())
      (adjacency_[b.index()]).insert(std::pair<size_type, size_type>(a.index(), i)); //add that a is a neighbor of b, connected by edge i; O(num_edges()) as well
      edges_.push_back({a.index(), b.index(), i}); //add internal_edge i to edges_ as proxy for edge i; O(1) for vector push_back
    }
    return Edge(this, i, a.index(), b.index());
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    adjacency_.clear();
  }

 private:
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals
  struct internal_node {
    const Point position;
    const size_type index; //index in proxy vector nodes_
  };
  struct internal_edge {
    const size_type node1; //index in proxy vector nodes_
    const size_type node2; //index in proxy vector nodes_
    const size_type index; //index in proxy vector edges_
  };

  std::vector<internal_edge> edges_; //proxy for Edge's
  std::vector<std::map<size_type, size_type>> adjacency_; //adjacency vector, map in index a with entry <b, i> denotes edge i between nodes a and b
  std::vector<internal_node> nodes_; //proxy for Node's
};

#endif // CME212_GRAPH_HPP
