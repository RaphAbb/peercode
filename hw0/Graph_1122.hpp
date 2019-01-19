#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>

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
  Graph() :
    next_uid_(0),
    points_(),
    node_uid_to_index_(),
    node_index_to_uid_(),
    edges_(),
    node_neighbors_() {}

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
    Node() : graph_(), uid_() {}

    /** Return this node's position. */
    const Point& position() const {
      return graph()->position(uid_);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph()->node_index(uid_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return n.graph() == this->graph() && n.index() == this->index();
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
      return this->index() < n.index();
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Pointer to parent graph
    const Graph* graph_;
    // Persistent ID for this node
    const size_type uid_;

    Node(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)), uid_(uid) {}

    const Graph* graph() const {
      return graph_;
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return this->num_nodes();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return points_.size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position) {
    size_type next_index = this->num_nodes();
    points_.push_back(position);
    node_uid_to_index_[next_uid_] = next_index;
    node_index_to_uid_[next_index] = next_uid_;
    node_neighbors_.push_back({});
    return Node(this, next_uid_++);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return node_uid_to_index_.count(n.uid_) > 0;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, node_uid(i));
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
    Edge() : graph_(), node1_(), node2_() {}

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->node(node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(node2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (this->node1_ == e.node1_ && this->node2_ == e.node2_)
        || (this->node1_ == e.node2_ && this->node2_ == e.node1_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      std::pair<size_type, size_type> this_edge = std::minmax(this->node1_, this->node2_);
      std::pair<size_type, size_type> other_edge = std::minmax(e.node1_, e.node2_);
      return this_edge.first < other_edge.first ||
        (this_edge.first == other_edge.first && this_edge.second < other_edge.second);
    }

  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Pointer to parent graph
    const Graph* graph_;
    // Persistent IDs of nodes that define this edge
    const size_type node1_;
    const size_type node2_;

    Edge(const Graph* graph, size_type node1, size_type node2)
      : graph_(const_cast<Graph*>(graph)), node1_(node1), node2_(node2) {}
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
    return Edge(this, edges_[i][0], edges_[i][1]);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    const std::vector<size_type> a_neighbors = node_neighbors_[a.index()];
    return std::find(a_neighbors.begin(), a_neighbors.end(), b.index()) != a_neighbors.end();
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
    if (!has_edge(a, b)) {
      edges_.push_back({a.index(), b.index()}),
      node_neighbors_[a.index()].push_back(b.index());
      node_neighbors_[b.index()].push_back(a.index());
    }
    return Edge(this, a.index(), b.index());
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    points_.clear();
    edges_.clear();
    node_uid_to_index_.clear();
    node_index_to_uid_.clear();
    node_neighbors_.clear();
  }

private:
  size_type next_uid_;
  std::vector<Point> points_;
  std::unordered_map<size_type, size_type> node_uid_to_index_;
  std::unordered_map<size_type, size_type> node_index_to_uid_;
  std::vector<std::array<size_type, 2>> edges_;
  std::vector<std::vector<size_type>> node_neighbors_;

  size_type node_index(const size_type uid) const {
    return node_uid_to_index_.at(uid);
  }

  size_type node_uid(const size_type index) const {
    return node_index_to_uid_.at(index);
  }

  const Point& position(const size_type uid) const {
    return points_[node_index(uid)];
  }
};

#endif // CME212_GRAPH_HPP
