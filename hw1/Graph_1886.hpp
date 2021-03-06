#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <unordered_set>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph {
private:
  /** Type of unique IDs for nodes, only used internally. */
  using uid_type = unsigned;

public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V>;

  /** Type of values stored in nodes. */
  using node_value_type = V;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  using node_iterator = NodeIterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  using incident_iterator = IncidentIterator;

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
    nodes_(),
    uid2index_(),
    index2uid_(),
    edges_(),
    adjacencies_() {}

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
  class Node : private totally_ordered<Node> {
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

    /** Return this node's value. */
    const node_value_type& value() const {
      return graph_->nodes_.at(uid_).value_;
    }

    /** Return this node's value. */
    node_value_type& value() {
      return graph_->nodes_.at(uid_).value_;
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_.at(uid_).position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->uid2index_.at(uid_);
    }

    /** Return this node's degree. */
    size_type degree() const {
      return graph_->adjacencies_.at(uid_).size();
    }

    /** Returns an iterator to the first edge incident to this node.
     *
     * @pre This node must be a valid node belonging to the parent graph.
     *
     * If no edges are incident to this node, the returned iterator will be equal
     * to @see edge_end().
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, uid_, graph_->adjacencies_.at(uid_).begin());
    }

    /** Returns an iterator to the element past the last edge incident to this node.
     *
     * @pre This node must be a valid node belonging to the parent graph.
     *
     * This element acts as a placeholder; attempting to access it results in
     * undefined behavior.
     */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, uid_, graph_->adjacencies_.at(uid_).end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return n.graph_ == this->graph_ && n.index() == this->index();
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
    Graph* graph_;
    // Persistent ID for this node
    const size_type uid_;

    Node(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)), uid_(uid) {}
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
   * @param[in] value The new node's value; defaults to the default value
   *   for the type.
   *
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: Amortized O(1).
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    uid2index_[next_uid_] = index2uid_.size();
    index2uid_.push_back(next_uid_);
    nodes_[next_uid_] = { position, value };
    adjacencies_[next_uid_] = {};
    return Node(this, next_uid_++);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return uid2index_.count(n.uid_) > 0;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, index2uid_.at(i));
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
  class Edge : private totally_ordered<Edge> {
  public:
    /** Construct an invalid Edge. */
    Edge() : graph_(), node1_uid_(), node2_uid_() {}

    /** Return a node of this Edge. */
    Node node1() const {
      return Node(graph_, node1_uid_);
    }

    /** Return the other node of this Edge. */
    Node node2() const {
      return Node(graph_, node2_uid_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return std::minmax(node1_uid_, node2_uid_) == std::minmax(e->node1_uid_, e->node2_uid_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return std::minmax(node1_uid_, node2_uid_) < std::minmax(e->node1_uid_, e->node2_uid_);
    }

  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Pointer to parent graph
    Graph* graph_;
    // Persistent IDs of nodes that define this edge
    const uid_type node1_uid_;
    const uid_type node2_uid_;

    Edge(const Graph* graph, const uid_type node1_uid, const uid_type node2_uid)
      : graph_(const_cast<Graph*>(graph)), node1_uid_(node1_uid), node2_uid_(node2_uid) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1).
   */
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(1).
   */
  Edge edge(size_type i) const {
    return Edge(this, edges_.at(i).at(0), edges_.at(i).at(1));
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: Average O(1), worst case O(num_edges()).
   */
  bool has_edge(const Node& a, const Node& b) const {
    return adjacencies_.at(a.uid_).count(b.uid_) > 0;
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
   * Complexity: Amortized O(1).
   */
  Edge add_edge(const Node& a, const Node& b) {
    if (!has_edge(a, b)) {
      edges_.push_back({a.uid_, b.uid_}),
      adjacencies_[a.uid_].insert(b.uid_);
      adjacencies_[b.uid_].insert(a.uid_);
    }
    return Edge(this, a.uid_, b.uid_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    uid2index_.clear();
    index2uid_.clear();
    adjacencies_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() : graph_(), node_uid_iterator_() {}

    /** Return the next node from this iterator.
     * @pre this != end after every previous invocation of @see operator++().
     */
    Node operator*() const {
      return Node(graph_, *node_uid_iterator_);
    }

    /** Increment this iterator. */
    NodeIterator& operator++() {
      node_uid_iterator_++;
      return *this;
    }

    /** Test for equality between this iterator and another NodeIterator.
     *
     * Note that the two iterators must pertain to the same graph to test equal.
     */
    bool operator==(const NodeIterator& other) const {
      return this->graph_ == other.graph_
          && this->node_uid_iterator_ == other.node_uid_iterator_;
    }

    /** Boolean inverse of @see operator==. */
    bool operator!=(const NodeIterator& other) const {
      return !(this->operator==(other));
    }

   private:
    friend class Graph;

    // Pointer to parent graph
    Graph* graph_;
    // Underlying iterator over the node UIDs
    std::vector<uid_type>::const_iterator node_uid_iterator_;

    NodeIterator(const Graph* graph, std::vector<uid_type>::const_iterator node_uid_iterator)
      : graph_(const_cast<Graph*>(graph)), node_uid_iterator_(node_uid_iterator) {}
  };

  /** Returns an iterator to the first node belonging to this graph.
   *
   * If there are no nodes in this graph, the returned iterator will be equal to @see node_end().
   */
  node_iterator node_begin() const {
    return NodeIterator(this, index2uid_.begin());
  }

  /** Returns an iterator to the element following the last node belonging to this graph.
   *
   * This element acts as a placeholder; attempting to access it results in undefined behavior.
   */
  node_iterator node_end() const {
    return NodeIterator(this, index2uid_.end());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    using uid_set_iterator_type = std::unordered_set<uid_type>::const_iterator;

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {}

    /** Return the next edge from this iterator.
     * @pre this != end after every previous invocation of @see operator++().
     */
    Edge operator*() const {
      return Edge(graph_, origin_uid_, *neighbor_);
    }

    /** Increment this iterator. */
    IncidentIterator& operator++() {
      ++neighbor_;
      return *this;
    }

    /** Test for equality between this and another IncidentIterator.
     *
     * Note that the two iterators must pertain to the same graph to test equal.
     */
    bool operator==(const IncidentIterator& other) const {
      return this->graph_ == other.graph_
          && this->origin_uid_ == other.origin_uid_
          && this->neighbor_ == other.neighbor_;
    }

    /** Boolean inverse of @see operator==. */
    bool operator!=(const IncidentIterator& other) const {
      return !(this->operator==(other));
    }

   private:
    friend class Graph;

    Graph* graph_;
    const uid_type origin_uid_;
    uid_set_iterator_type neighbor_;

    IncidentIterator(const Graph* graph, const uid_type origin_uid, const uid_set_iterator_type neighbor)
      : graph_(const_cast<Graph*>(graph)), origin_uid_(origin_uid), neighbor_(neighbor) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    using edge_iter_type = std::vector<std::array<uid_type, 2>>::const_iterator;

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {}

    /** Return the next edge in this graph.
     * @pre this != end after every previous invocation of @see operator++().
     */
    Edge operator*() const {
      std::array<uid_type, 2> nodes = *edge_iterator_;
      return Edge(graph_, nodes[0], nodes[1]);
    }

    /** Increment this iterator. */
    EdgeIterator& operator++() {
      ++edge_iterator_;
      return *this;
    }

    /** Test for equality between this and another EdgeIterator.
     *
     * Note that the two iterators must pertain to the same graph to test equal.
     */
    bool operator==(const EdgeIterator& other) const {
      return this->graph_ == other.graph_
          && this->edge_iterator_ == other.edge_iterator_;
    }

    /** Boolean inverse of @see operator==. */
    bool operator!=(const EdgeIterator& other) const {
      return !(this->operator==(other));
    }

   private:
    friend class Graph;
    // Pointer to parent graph
    Graph* graph_;
    // Underlying iterator over internal edge representation
    edge_iter_type edge_iterator_;

    EdgeIterator(const Graph* graph, edge_iter_type edge_iterator)
      : graph_(const_cast<Graph*>(graph)), edge_iterator_(edge_iterator) {}
  };

  /** Returns an iterator to the first edge in the graph.
   *
   * If there are no edges in this graph, the returned iterator will be equal to @see edge_end().
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, edges_.begin());
  }

  /** Returns an iterator to the element past the last edge in the graph.
   *
   * This element is a placeholder; attempting to access it results in undefined behavior.
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, edges_.end());
  }

 private:
  // A handy container to keep track of a node's position and value in the same place.
  struct NodeInfo {
    Point position_;
    node_value_type value_;
  };

  // To be used in generating a unique ID for each node. Implicitly limits the total
  // number of nodes in this graph to UINT_MAX, which seems like a reasonable
  // constraint right now.
  uid_type next_uid_;
  // Allows us to locate the NodeInfo for any node by its UID, in O(1) time on average.
  std::unordered_map<uid_type, NodeInfo> nodes_;
  // The master mapping of each node's UID to its corresponding "index" in the graph.
  std::unordered_map<uid_type, size_type> uid2index_;
  // The reverse mapping, allowing reverse lookups of a node's UID given its index.
  std::vector<uid_type> index2uid_;
  // A vector of all edges, represented as pairs of node UIDs. Allows for O(1) lookups
  // of edges by their "index", and for straightforward iteration over edges.
  std::vector<std::array<uid_type, 2>> edges_;
  // The adjacency list representation of the edges. Maintaining both the adjacency
  // list and a flat vector of edges increases our memory requirements, but enables
  // very fast access to both a single edge or a set of neighbors of a given node.
  std::unordered_map<uid_type, std::unordered_set<uid_type>> adjacencies_;
};

#endif // CME212_GRAPH_HPP
