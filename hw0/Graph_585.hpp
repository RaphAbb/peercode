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

  /** Pre-declaration for internal type for node */
  struct internal_node;

  /** Pre-declaration for internal type for edge */
  struct internal_edge;


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

  /** Construct an empty graph, initializing node and edge vectors. */
  Graph()
    : nodes_(), edges_() {
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
      return graph_->nodes_[index_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->nodes_[index_].index;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_->has_node(n) and this->index() == n.index());
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
      assert(graph_->has_node(n));
      return (this->index() < n.index());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    /** Pointer back to Graph */
    Graph* graph_ = nullptr;

    /** Index of node */
    size_type index_;

    /** Private constructor */
    Node(const Graph* graph, size_type index)
      : graph_(const_cast<Graph*>(graph)), index_(index) {
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes_.size();
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
  Node add_node(const Point& pos) {
    // Initialize internal node struct
    size_type id = num_nodes();
    internal_node new_node;
    new_node.position = pos;
    new_node.index = id;
    new_node.node = Node(this, id);
    // Add to vector of nodes
    nodes_.push_back(new_node);
    // return proxy
    return new_node.node;
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
    assert(i >= 0 and i < nodes_.size());
    return nodes_[i].node;
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
      return graph_->edges_[index_].node1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->edges_[index_].node2;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      node_type node1 = graph_->edges_[index_].node1;
      node_type node2 = graph_->edges_[index_].node2;
      // Check if both nodes match
      return (node1 == e.node1() and node2 == e.node2())
          or (node2 == e.node1() and node1 == e.node2());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // Vector of graph's edges
      std::vector<internal_edge> edges = graph_->edges_;
      // Loop through vector, if this edge < e then this will be first
      for(unsigned int i = 0; i < edges.size(); ++i) {
        if(edges[i].edge == *this) return true;
        if(edges[i].edge == e) return false;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    /** Graph that this edge belongs to */
    Graph* graph_ = nullptr;

    /** Index */
    size_type index_;

    /** Private constructor for a valid Edge */
    Edge(Graph* graph, size_type index)
      : graph_(const_cast<Graph*>(graph)), index_(index) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i >= 0 and i < edges_.size());
    return edges_[i].edge;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Assert Nodes A and B are on this Graph
    assert(this->has_node(a) and this->has_node(b));
    // Loop through Node A's edges to see if it is connected with B
    std::vector<edge_type> edge_list = nodes_[a.index()].connections;
    for(const auto& edge : edge_list) {
      if (edge.node1() == b or edge.node2() == b) return true;
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
    assert(this->has_node(a) and this->has_node(b)); // valid nodes of the graph
    assert(not (a == b));                              // distinct

    // Edge does not exist yet
    if(not has_edge(a,b)) {
      size_type id = this->num_edges();
      Edge new_edge(this,id);                            // Construct a valid Edge proxy
      edges_.emplace_back(a,b,new_edge);                 // Create the "real" internal struct
      nodes_[a.index()].connections.push_back(new_edge); // Add edge proxy to appropriate lists
      nodes_[b.index()].connections.push_back(new_edge);
      return new_edge;
    }

    // Edge exists
    std::vector<edge_type> edge_list = nodes_[a.index()].connections;
    // Loop through the list of edges associated with Node A to find current edge
    for(const auto& edge : edge_list) {
      if (edge.node1() == b or edge.node2() == b) return edge;
    }

    return Edge();
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
  }

 private:

  /** Internal type for Node */
  struct internal_node {
    Point position;
    size_type index;
    std::vector<edge_type> connections; // List of edges associated with this node
    Node node;
  };

  /** Internal type for Edge */
  struct internal_edge {
    Node node1;
    Node node2;
    Edge edge;
    // Constructor for emplace_back
    internal_edge(Node a, Node b, Edge e) : node1(a), node2(b), edge(e) {}
  };

  /** Vector containing node structs (in order of assignment/indexing) */
  std::vector<internal_node> nodes_;

  /** Vector containing edge structs (in order of assignment/indexing) */
  std::vector<internal_edge> edges_;
};

#endif // CME212_GRAPH_HPP
