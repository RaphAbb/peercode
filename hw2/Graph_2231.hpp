#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <exception>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template<typename V, typename E>
class Graph {
private:
  struct internal_edge;
  struct internal_node;
  unsigned int node_uid = 0;
  unsigned int edge_uid = 0;

public:
  using graph_type = Graph;
  typedef V node_value_type;
  typedef E edge_value_type;

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
  using size_type = unsigned;
  using map_iterator = std::unordered_map<size_type, size_type>::const_iterator;

  /** Construct an empty graph. */
  Graph() = default;

  ~Graph() = default;

  class MissingNode : public std::exception {
    virtual const char *what() const throw() {
      return "Node is not currently part of Graph";
    }
  };

  class MissingEdge : public std::exception {
    virtual const char *what() const throw() {
      return "Edge is not currently part of Graph";
    }
  };

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
    Node() = default;

    /** @throws MissingNode if node is not part of graph*/
    Point &position() {
      return graph_->nodes[index()].point;
    }

    /** Return this node's position. */
    /** @throws MissingNode if node is not part of graph*/
    const Point &position() const {
      return graph_->nodes[index()].point;
    }

    /** Return this node's index, a number in the range [0, graph_size).
     * @throws MissingNode if node is not part of graph
     */
    size_type index() const {
      try {
        return graph_->nodes_uid_to_index.at(uid_);
      } catch (const std::out_of_range &oor) {
        throw MissingNode{};
      }
    }

    /** @throws MissingNode if node is not part of graph*/
    const node_value_type &value() const {
      return graph_->nodes[index()].value;
    }

    /** @throws MissingNode if node is not part of graph*/
    node_value_type &value() {
      return graph_->nodes[index()].value;
    }

    size_type degree() const {
      return graph_->degree(uid_);
    }

    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, uid_);
    }

    incident_iterator edge_end() const {
      return IncidentIterator(graph_, uid_, true);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node &n) const {
      return this->index() == n.index() && graph_->has_node(n);
    }

    /** Test whether this node is less than @a n in a global order.
     * Ordering is based on index and only is defined on Nodes from the
     * same graph.
     *
     * Checks first ordering the of graph and then the ordering of internal index when same graph
     */
    bool operator<(const Node &n) const {
      if (this->graph_->has_node(n)) {
        return this->index() < n.index();
      } else {
        return this->graph_ < n.graph_;
      }
    }

  private:
    friend class Graph;

    Graph *graph_;
    size_type uid_;

    Node(const Graph *graph, size_type uid)
        : graph_(const_cast<Graph *>(graph)), uid_(uid) {
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
    internal_node int_node{position, val, node_uid};
    nodes.push_back(int_node);

    nodes_uid_to_index[node_uid] = index;

    std::unordered_map<size_type, size_type> hash_map{};
    node_edge_map[node_uid] = hash_map;


    return {this, node_uid++};
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node &n) const {
    try {
      n.index();
    } catch (MissingNode &e) {
      return false;
    }

    return n.graph_ == this;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   * @throws MissingNode{}
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if (i >= num_nodes() || i < 0) {
      throw MissingNode{};
    }
    return {this, nodes[i].uid};
  }

  /** Removes Node by Node Reference
   * @pre 0 <= @a i < num_nodes()
   * @post  - Node is not in the graph
   *        - Calling graph related functions of Node will throw MissingNode{}
   *        - All edges incident to the node are removed and will throw MissingEdge{} if graph functions are used
   *        - Internal indexes of all other nodes may have changed
   *          (current implementation is just the last is moved to deleted index, but this might change in the future)
   *        - Adjacency List and Node index lookup updated
   *        - Edge and Incident iterators are invalid, as edges have been removed
   *        - Node iterator is still valid as long as it is not pointing at the last element or one past (
   *          (you will need to fetch a new end iterator and not advance the current iterator if you've just deleted
   *          a node)
   *
   * @returns 1 if the node was present and deleted
   *          0 if the node wasn't previously in the graph
   *
   * Complexity: O(num_neighbors).
  */
  size_type remove_node(const Node &node) {
    size_type i = 0;
    try {
      i = node.index();
    } catch (MissingNode &e) {
      return 0;
    }

    std::vector<size_type> edges;

    for (auto edge_i = node.edge_begin(); edge_i != node.edge_end(); ++edge_i) {
      edges.push_back((*edge_i).uid_);
    }

    for (size_type edge_index : edges) {
      remove_edge(Edge{node.graph_, edge_index});
    }

    // Move last element to deleted index
    nodes[i] = nodes[nodes.size() - 1];
    nodes_uid_to_index[nodes[i].uid] = i;
    nodes_uid_to_index.erase(node.uid_);
    node_edge_map.erase(node.uid_);
    nodes.pop_back();

    return 1;
  }

  /** @brief see remove_node(Node&)
  *   @pre iterator is not == the end_iterator */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return n_it;
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
    Edge() = default;

    /** Return a node of this Edge
     * Throws assertion error if edge has index -1 or no graph
     */
    Node node1() const {
      assert(uid_ >= 0 && graph_ != nullptr);
      if (relative_to_node_b_) {
        return graph_->edges.at(internal_index()).b_node;
      }
      return graph_->edges.at(internal_index()).a_node;
    }

    /** Return the other node of this Edge
     * @throws MissingEdge if edge is not part of graph
     */
    Node node2() const {
      if (relative_to_node_b_) {
        return graph_->edges.at(internal_index()).a_node;
      }
      return graph_->edges.at(internal_index()).b_node;
    }

    /**
     * @throws MissingEdge if edge is not part of graph
     */
    double length() const {
      return norm_2(node1().position() - node2().position());
    }

    /** @throws MissingEdge if edge is not part of graph*/
    edge_value_type &value() {
      return graph_->edges.at(internal_index()).value;
    }

    /** @throws MissingEdge if edge is not part of graph*/
    const edge_value_type &value() const {
      return graph_->edges.at(internal_index()).value;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge &e) const {
      return this->internal_index() == e.internal_index() && graph_->has_edge(e.node1(), e.node2());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * Checks first ordering the of graph and then the ordering of internal index when same graph
     */
    bool operator<(const Edge &e) const {
      if (this->graph_->has_edge(e.node1(), e.node2())) {
        return this->internal_index() < e.internal_index();
      } else {
        return this->graph_ < e.graph_;
      }
    }

  private:
    friend class Graph;

    Graph *graph_;
    size_type uid_;
    bool relative_to_node_b_ = false;

    Edge(const Graph *graph, size_type index)
        : graph_(const_cast<Graph *>(graph)), uid_(index) {
    }

    Edge(const Graph *graph, size_type index, bool relative_to_node_b)
        : graph_(const_cast<Graph *>(graph)), uid_(index), relative_to_node_b_(relative_to_node_b) {
    }

    /** @throws MissingEdge if edge is not part of graph*/
    size_type internal_index() const {
      try {
        return graph_->edges_uid_to_index.at(uid_);
      } catch (const std::out_of_range &oor) {
        throw MissingEdge{};
      }
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
   * @throws MissingEdge
   *
   * Complexity: O(1)
   */
  Edge edge(size_type i) const {
    if (i < 0 || i >= num_edges()) {
      throw MissingEdge{};
    }
    return Edge(this, edges[i].uid);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: O(1)
   */
  bool has_edge(const Node &a, const Node &b) const {
    return node_edge_map.count(a.uid_) > 0 && node_edge_map.at(a.uid_).count(b.uid_);
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
  Edge add_edge(const Node &a, const Node &b, const edge_value_type &val = edge_value_type()) {
    Edge edge;
    if (has_edge(a, b)) {
      edge = Edge(this, node_edge_map.at(a.uid_).at(b.uid_));
    } else {
      size_type index = num_edges();
      internal_edge int_edge{a, b, edge_uid, val};
      edges.push_back(int_edge);
      edges_uid_to_index[edge_uid] = index;
      edge = Edge(this, edge_uid++);
      add_edge_to_lookup_by_node(a, b, edge);
    }
    return edge;
  }

  /** Removes Edge by End Nodes
   * @pre   - Nodes must be valid
   *        - If Edge exists then Edge must also exist in node_edge_map
   * @post  - Edge is not in the graph
   *        - Edge will throw MissingEdge{} if graph functions are used
   *        - Internal indexes of all other edges may have changed
   *          (current implementation is just the last is moved to deleted index, but this might change in the future)
   *         - Adjacency List and Edge index lookup updated
   *        - Incident iterators are invalid, as edges have been removed
   *        - Edge iterator is still valid as long as it is not pointing at the last element or one past (
   *          (you will need to fetch a new end iterator and not advance the current iterator if you've just deleted
   *          a edge)
   *
   * @returns 1 if the edge was present and deleted
   *          0 if the edge wasn't previously in the graph
   *
   * Complexity: O(1).
 */
  size_type remove_edge(const Node &node_a, const Node &node_b) {
    size_type edge_uid;
    if (has_edge(node_a, node_b)) {
      edge_uid = node_edge_map.at(node_a.uid_).at(node_b.uid_);
    } else {
      return 0;
    }

    node_edge_map[node_a.uid_].erase(node_b.uid_);
    node_edge_map[node_b.uid_].erase(node_a.uid_);

    size_type index = edges_uid_to_index[edge_uid];
    edges[index] = edges[edges.size() - 1];
    edges_uid_to_index[edges[index].uid] = index;
    edges_uid_to_index.erase(edge_uid);
    edges.pop_back();
    return 1;
  }

  /** @brief See remove_edge(const Node &node_a, const Node &node_b)*/
  size_type remove_edge(const Edge &edge) {
    size_type index;
    try {
      index = edges_uid_to_index.at(edge.uid_);
    } catch (...) {
      return 0;
    }
    auto node_a = edges[index].a_node;
    auto node_b = edges[index].b_node;

    return remove_edge(node_a, node_b);
  }


  /** @brief See remove_edge(const Node &node_a, const Node &node_b)
   *  @pre iterator is not == the end_iterator */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return e_it;
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
    nodes_uid_to_index.clear();
    edges_uid_to_index.clear();
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
    NodeIterator() = default;

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

    Graph<V, E> *graph_;
    size_type iterator_index_;

    NodeIterator(const Graph<V, E> *graph, size_type node_index)
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
      return graph_->getEdge(edge_index, node_uid_);
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
      bool node_uid_equal = node_uid_ == other_iterator.node_uid_;
      bool map_it_equal = map_iterator_ == other_iterator.map_iterator_;

      return graphs_equal && node_uid_equal && map_it_equal;
    }

  private:
    friend class Graph;

    map_iterator map_iterator_;
    Graph<V, E> *graph_;
    size_type node_uid_;

    IncidentIterator(const Graph<V, E> *graph, const size_type node_uid, bool isEnd = false)
        : graph_(const_cast<Graph *>(graph)), node_uid_(node_uid) {
      if (isEnd) {
        map_iterator_ = graph->get_end_node_neighbor_iterator(node_uid);
      } else {
        map_iterator_ = graph->get_node_neighbor_iterator(node_uid);
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

    Graph<V, E> *graph_;
    size_type iterator_index_;

    EdgeIterator(const Graph<V, E> *graph, size_type index)
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
  std::unordered_map<size_type, size_type> nodes_uid_to_index;
  std::unordered_map<size_type, size_type> edges_uid_to_index;
  std::vector<internal_node> nodes;
  std::vector<internal_edge> edges;

  struct internal_edge {
    Node a_node;
    Node b_node;
    size_type uid;
    edge_value_type value;
  };

  struct internal_node {
    Point point;
    node_value_type value;
    size_type uid;
  };

  /** Adds Edge to lookup by both node indicies as the top level key*/
  void add_edge_to_lookup_by_node(const Node &a, const Node &b, const Edge &edge) {
    size_type a_uid = a.uid_;
    size_type b_uid = b.uid_;
    add_edge_to_node(a_uid, b_uid, edge);
    add_edge_to_node(b_uid, a_uid, edge);
  }

  /** Adds Edge to lookup by key of a's index and then using b as the index to the innermap */
  void add_edge_to_node(const size_type a_uid, const size_type b_uid, const Edge &edge) {
    bool has_node_uid_a = node_edge_map.count(a_uid);
    if (!has_node_uid_a) {
      std::unordered_map<size_type, size_type> hash_map{};
      hash_map[b_uid] = edge.uid_;
      node_edge_map[a_uid] = hash_map;
    } else {
      node_edge_map[a_uid][b_uid] = edge.uid_;
    }
  }

  /** Complexity O(1) */
  size_type degree(const size_type node_uid) const {
    if (node_edge_map.count(node_uid)) {
      return node_edge_map.at(node_uid).size();
    } else {
      return 0;
    }
  }

  /** Complexity O(1) */
  map_iterator get_end_node_neighbor_iterator(const size_type node_uid) const {
    return node_edge_map.at(node_uid).end();
  }

  /** Complexity O(1) */
  map_iterator get_node_neighbor_iterator(const size_type node_uid) const {
    return node_edge_map.at(node_uid).begin();
  }

  /** Complexity O(1)
   * @brief Returns new edge object by uid with the second argument as node1.
   * @throws MissingEdge
   **/
  Edge getEdge(size_type edge_uid, size_type node_uid) const {
    size_type index;
    try {
      index = edges_uid_to_index.at(edge_uid);
    } catch (...) {
      throw MissingEdge{};
    }

    internal_edge edge = edges.at(index);
    if (edge.b_node.uid_ == node_uid) {
      return Edge(this, edge_uid, true);
    }
    return Edge(this, edge_uid, false);
  }
};

#endif // CME212_GRAPH_HPP
