#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
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
 public:

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
    using size_type = unsigned;

 private:

  /* Container for storing node data */
  typedef struct internal_node {
    Point point;   // The position of the node in (x, y, z)
    size_type uid; // The unique identifcation for a node
  } internal_node;

  /* Map from node index to node */
  typedef std::unordered_map<size_type, internal_node> Nodes;
  /* Map from node UID to node index */
  typedef std::unordered_map<size_type, size_type> IDMap;

  /* Map from edge index to edge nodes */
  typedef std::unordered_map< size_type, 
    std::pair<size_type, size_type>> Edges;
  /* Map from edge nodes to edge index (i.e. an 'inverted'
  version of the Edges map) */
  typedef std::map<
    std::pair<size_type, size_type>,
    size_type> EdgesInv;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() 
    : nodes_(), edges_(), idmap_(), edgeinv_(),
      next_uid_(0), size_(0), num_edges_(0) {}

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
      this->G_ = nullptr;
    }

    /** Return this node's position. */
    const Point& position() const {
      // Find the position by UID
      return this->G_->nodes_[this->uid_].point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // Invalid node
      if (this->G_ == nullptr) {
        std::cerr << "ERROR: This node does not belong to a graph.";
        exit(1);
      }
      // Find the index by UID
      return this->G_->idmap_[this->uid_];
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (n.uid_ == this->uid_ and n.G_ == this->G_) {
        return true;
      }
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
      if (this->uid_ < n.uid_) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Point to the owner of this node
    Graph* G_; 
    // Uniquely identify a node
    size_type uid_;

    // Private constructor for Graph to construct a valid Node object.
    Node(const Graph* G, size_type uid)
        : G_(const_cast<Graph*>(G)), uid_(uid) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return this->size_;
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

    // Build internal node for storage
    internal_node in;
    in.point = position;
    in.uid = this->next_uid_;

    // Add internal node
    this->nodes_.insert(std::make_pair(this->size_, in));
    // Add to ID map
    this->idmap_.insert(std::make_pair(this->next_uid_, this->size_));
    // Increment next UID
    this->next_uid_ += 1;
    // Increment size of Graph
    this->size_ += 1;
    // Return valid node
    return Node(this, this->next_uid_- 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.G_ == this and idmap_.count(n.uid_)) {
      return true;
    }
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // Find correct UID
    size_type uid = this->nodes_.at(i).uid;
    // Return a node with the correct UID
    return Node(this, uid);
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
      this->G_ = nullptr;
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(this->G_, this->nodepair_.first);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(this->G_, this->nodepair_.second);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {

      if (this->G_ == e.G_
      and this->nodepair_.first == e.nodepair_.first
      and this->nodepair_.second == e.nodepair_.second) {
        return true;
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (this->index() < e.index()) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Private data members
    Graph* G_;
    std::pair<size_type, size_type> nodepair_;

    // Private helper function to get the index of an Edge
    size_type index() const {
      return this->G_->edgeinv_.at(this->nodepair_);
    }
    // Private constructor for Graph to construct a valid Edge object
    Edge(const Graph* G, std::pair<size_type, size_type> nodepair)
        : G_(const_cast<Graph*>(G)), nodepair_(nodepair) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return this->num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // Get node UIDs from index i
    return Edge(this, this->edges_.at(i));
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {

    // Check if nodes are the same
    if (a == b) {
      return false;
    }
    // Check if nodes are in the graph
    if (has_node(a) and has_node(b)) {

      if (this->edgeinv_.count(std::make_pair(
          std::min(a.uid_, b.uid_), 
          std::max(a.uid_, b.uid_)))) {
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

    // Check if nodes are in the graph
    if (! has_node(a) or ! has_node(b)) {
        std::cerr << "ERROR: Nodes a and b must belong to this graph.";
        exit(1);
    }

    // Create pair of node UIDs
    std::pair<size_type, size_type> nodepair = std::make_pair(
          std::min(a.uid_, b.uid_), 
          std::max(a.uid_, b.uid_));

    // Check if edge is in the graph (if so, return that edge)
    if (has_edge(a, b)) {

      // return Edge();
      return Edge(this, nodepair);
      
    }
    //// Add edge to graph

    // (a) Add to list of edges
    this->edges_.insert(std::make_pair(this->num_edges_, nodepair));
    // (b) Add to inverted edge map
    this->edgeinv_.insert(std::make_pair(nodepair, this->num_edges_)); 
    // Update number of edges
    this->num_edges_ += 1;
    // New edge
    return Edge(this, nodepair);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // Reset all internal data structures
    this->nodes_ = Nodes();
    this->edges_ = Edges();
    this->idmap_ = IDMap();
    this-> edgeinv_ = EdgesInv();
    this->next_uid_ = 0;
    this->size_ = 0;
    this->num_edges_ = 0;
  }

 private:

  Nodes nodes_;
  Edges edges_;
  IDMap idmap_;
  EdgesInv edgeinv_;

  // Next unique UID for a node
  size_type next_uid_;
  // Number of nodes in the graph
  size_type size_;
  // Number of edges in the graph
  size_type num_edges_;

};

#endif // CME212_GRAPH_HPP
