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
template<typename V>
class Graph {
private:
  struct internal_edge;
  struct internal_node;
  const static int INVALID_INDEX = -1;

public:
  using graph_type = Graph;
  using node_value_type = V;

  class Node;

  using node_type = Node;

  class Edge;

  using edge_type = Edge;

  class NodeIterator;

  using node_iterator = NodeIterator;

  class EdgeIterator;

  using edge_iterator = EdgeIterator;

  class IncidentIterator;

  using incident_iterator = IncidentIterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = signed;
  using map_iterator = std::unordered_map<size_type, size_type>::const_iterator;

  /** Construct an empty graph. */
  Graph() = default;

  ~Graph() = default;

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
    Node() {
      this->index_ = INVALID_INDEX;
    }

    /** Return this node's position. */
    const Point &position() const {
      return graph_->nodes[index_].point;
    }

    /** Return this node's index, a number in the range [0, graph_size).
     * If an invalid node is constructed, the index will be -1
     */
    size_type index() const {
      return this->index_;
    }

    node_value_type &value() {
      return graph_->nodes[index_].value;
    }

    const node_value_type &value() const {
      return graph_->nodes[index_].value;
    }

    size_type degree() const {
      return graph_->degree(index_);
    }

    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, index_);
    }

    incident_iterator edge_end() const {
      return IncidentIterator(graph_, index_, true);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node &n) const {
      return this->index_ == n.index() && graph_->has_node(n);
    }

    /** Test whether this node is less than @a n in a global order.
     * Ordering is based on index and only is defined on Nodes from the
     * same graph.
     *
     * Throws assertion error if the nodes are not from the same graph.
     */
    bool operator<(const Node &n) const {
      assert(this->graph_->has_node(n));
      return this->index_ < n.index();
    }

  private:
    friend class Graph;

    Graph *graph_;
    size_type index_;

    Node(const Graph *graph, size_type index)
        : graph_(const_cast<Graph *>(graph)), index_(index) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes.size();
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
  Node add_node(const Point &position, const node_value_type &val = node_value_type()) {
    size_type index = num_nodes();
    internal_node int_node{position, val};
    nodes.push_back(int_node);

    std::unordered_map<size_type, size_type> hash_map{};
    node_edge_map[index] = hash_map;

    return {this, index};
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node &n) const {
    return n.graph_ == this;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    size_type index = i;
    if (i < 0 || i >= num_nodes()) {
      index = INVALID_INDEX;
    }
    return {this, index};
  }

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Edge> {
  public:
    /** Construct an invalid Edge. */
    Edge() {
      this->index_ = INVALID_INDEX;
    }

    /** Return a node of this Edge
     * Throws assertion error if edge has index -1 or no graph
     */
    Node node1() const {
      assert(index_ >= 0 && graph_ != nullptr);
      if (relative_to_node_b_) {
        return graph_->edges.at(index_).b_node;
      }
      return graph_->edges.at(index_).a_node;
    }

    /** Return the other node of this Edge
     * Throws assertion error if edge has index -1 or no graph
     */
    Node node2() const {
      assert(index_ >= 0 && graph_ != nullptr);
      if (relative_to_node_b_) {
        return graph_->edges.at(index_).a_node;
      }
      return graph_->edges.at(index_).b_node;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge &e) const {
      return this->index_ == e.index_ && graph_->has_edge(e.node1(), e.node2());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge &e) const {
      assert(this->graph_->has_edge(e.node1(), e.node2()));
      return this->index_ < e.index_;
    }

  private:
    friend class Graph;

    Graph *graph_;
    size_type index_;
    bool relative_to_node_b_ = false;

    Edge(const Graph *graph, size_type index)
        : graph_(const_cast<Graph *>(graph)), index_(index) {
    }

    Edge(const Graph *graph, size_type index, bool relative_to_node_b)
        : graph_(const_cast<Graph *>(graph)), index_(index), relative_to_node_b_(relative_to_node_b) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(1)
   */
  Edge edge(size_type i) const {
    size_type index = i;
    if (i < 0 || i >= num_edges()) {
      index = INVALID_INDEX;
    }
    return Edge(this, index);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: O(1)
   */
  bool has_edge(const Node &a, const Node &b) const {
    return node_edge_map.count(a.index()) > 0 && node_edge_map.at(a.index()).count(b.index());
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
   * Complexity: O(1)
   */
  Edge add_edge(const Node &a, const Node &b) {
    Edge edge;
    if (has_edge(a, b)) {
      edge = Edge(this, node_edge_map.at(a.index()).at(b.index()));
    } else {
      size_type index = num_edges();
      internal_edge int_edge{a, b};
      edges.push_back(int_edge);
      edge = Edge(this, index);
      add_edge_to_lookup_by_node(a, b, edge);
    }
    return edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();
    edges.clear();
    node_edge_map.clear();
  }

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node *;                    // Pointers to elements
    using reference         = Node &;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    /** @pre valid NodeIterator must be constructed */
    Node operator*() const {
      return graph_->node(iterator_index_);
    }

    NodeIterator &operator++() {
      iterator_index_++;
      return *this;
    }

    /** Tests iterator equality
     * @param node_iterator  Another node iterator to compare to
     * @return bool indicating whether the graphs are equal and iterator is pointing to the same index
     *
     * @pre valid NodeIterator must be constructed
    */
    bool operator==(const NodeIterator &node_iterator) const {
      return this->graph_ == node_iterator.graph_ && this->iterator_index_ == node_iterator.iterator_index_;

    }

  private:
    friend class Graph;

    Graph<V> *graph_;
    size_type iterator_index_;

    NodeIterator(const Graph<V> *graph, size_type node_index)
        : graph_(const_cast<Graph *>(graph)), iterator_index_(node_index) {
    }
  };

  node_iterator node_begin() const {
    return {this, 0};
  }

  node_iterator node_end() const {
    return {this, num_nodes()};
  }

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge *;                    // Pointers to elements
    using reference         = Edge &;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    Edge operator*() const {
      size_type edge_index = (*map_iterator_).second;
      return graph_->getEdge(edge_index, node_index_);
    }

    IncidentIterator &operator++() {
      map_iterator_++;
      return *this;
    }

    /** Tests iterator equality
     * @param other_iterator  Another incident iterator to compare to
     * @return bool indicating whether the graphs are equal, iterators are pointing to the element
     *
     * @pre valid IncidentIterator must be constructed
    */
    bool operator==(const IncidentIterator &other_iterator) const {
      bool graphs_equal = graph_ == other_iterator.graph_;
      bool node_index_equal = node_index_ == other_iterator.node_index_;
      bool map_it_equal = map_iterator_ == other_iterator.map_iterator_;

      return graphs_equal && node_index_equal && map_it_equal;
    }

  private:
    friend class Graph;

    map_iterator map_iterator_;
    Graph<V> *graph_;
    size_type node_index_;

    IncidentIterator(const Graph<V> *graph, const size_type node_index, bool isEnd = false)
        : graph_(const_cast<Graph *>(graph)), node_index_(node_index) {
      if (isEnd) {
        map_iterator_ = graph->get_end_node_neighbor_iterator(node_index);
      } else {
        map_iterator_ = graph->get_node_neighbor_iterator(node_index);
      }
    }

  };

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge *;                    // Pointers to elements
    using reference         = Edge &;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    Edge operator*() const {
      return graph_->edge(iterator_index_);
    }

    EdgeIterator &operator++() {
      iterator_index_++;
      return *this;
    }

    /** Tests iterator equality
     * @param edge_iterator  Another edge iterator to compare to
     * @return bool indicating whether the graphs are equal and iterator is pointing to the same index
     *
     * @pre valid EdgeIterator must be constructed
    */
    bool operator==(const EdgeIterator &edge_iterator) const {
      return this->graph_ == edge_iterator.graph_ && this->iterator_index_ == edge_iterator.iterator_index_;
    }

  private:
    friend class Graph;

    Graph<V> *graph_;
    size_type iterator_index_;

    EdgeIterator(const Graph<V> *graph, size_type index)
        : graph_(const_cast<Graph *>(graph)), iterator_index_(index) {
    }
  };

  edge_iterator edge_begin() const {
    return {this, 0};
  }

  edge_iterator edge_end() const {
    return {this, num_edges()};

  }

private:
  /** Provides easy lookup of edges by either node as the top level key*/
  std::unordered_map<size_type, std::unordered_map<size_type, size_type>> node_edge_map;
  std::vector<internal_node> nodes;
  /** Use for insertion only, reording or deletion requires updating node edge lookup */
  std::vector<internal_edge> edges;

  struct internal_edge {
    Node a_node;
    Node b_node;
  };

  struct internal_node {
    Point point;
    node_value_type value;
  };

  /** Adds Edge to lookup by both node indicies as the top level key*/
  void add_edge_to_lookup_by_node(const Node &a, const Node &b, const Edge &edge) {
    size_type a_index = a.index();
    size_type b_index = b.index();
    add_edge_to_node(a_index, b_index, edge);
    add_edge_to_node(b_index, a_index, edge);
  }

  /** Adds Edge to lookup by key of a's index and then using b as the index to the innermap */
  void add_edge_to_node(const size_type a_index, const size_type b_index, const Edge &edge) {
    bool has_node_index_a = node_edge_map.count(a_index);
    if (!has_node_index_a) {
      std::unordered_map<size_type, size_type> hash_map{};
      hash_map[b_index] = edge.index_;
      node_edge_map[a_index] = hash_map;
    } else {
      node_edge_map[a_index][b_index] = edge.index_;
    }
  }

  /** Complexity O(1) */
  size_type degree(const size_type a_index) const {
    return node_edge_map.count(a_index);
  }

  /** Complexity O(1) */
  map_iterator get_end_node_neighbor_iterator(const size_type a_index) const {
    return node_edge_map.at(a_index).end();
  }

  /** Complexity O(1) */
  map_iterator get_node_neighbor_iterator(const size_type a_index) const {
    return node_edge_map.at(a_index).begin();
  }

  /** Complexity O(1) */
  Edge getEdge(size_type edge_index, size_type node_index) const {
    internal_edge edge = edges.at(edge_index);
    if (edge.b_node.index() == node_index) {
      return Edge(this, edge_index, true);
    }
    return Edge(this, edge_index, false);
  }

};

#endif // CME212_GRAPH_HPP
