#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <set>
#include <unordered_map>

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

  struct NodeInternal;
  struct EdgeInternal;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V>;

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

  // Index of the true Node in the node_list_
  using NodeListIndex = size_type;

  // Index of true Edge in the edge_list_
  using EdgeListIndex = size_type;

  // ID of a node (not the index of underlying node_list_)
  using NodeID = size_type;

  // ID of an edge (not the index of underlying edge_list_)
  using EdgeID = size_type;

  using node_value_type = V;

  // Vector containing true nodes.
  using NodeList = std::vector<NodeInternal>;

  // Vector containing true edges.
  using EdgeList = std::vector<EdgeInternal>;

  // Stores the indices of the edges adjacent to the node.
  using AdjacentEdges = std::set<EdgeID>;

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
      node_id_ = std::numeric_limits<NodeID>::max();
      graph_ptr_ = nullptr;
    }

    Node(NodeID node_id, const Graph *graph_ptr)
      : node_id_(node_id), graph_ptr_(const_cast<Graph*>(graph_ptr)) {
    }

    /** Return this node's position. */
    const Point& position() const {
      // node is valid if index is within size of the node_list_
      assert((graph_ptr_->node_id_map_).count(node_id_));
      const NodeListIndex node_list_index = (graph_ptr_->node_id_map_).at(node_id_);
      return (graph_ptr_->node_list_)[node_list_index].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // NOTE: this function returns a NodeID, not a NodeListIndex
      assert((graph_ptr_->node_id_map_).count(node_id_));
      return node_id_;
    }

    // Redundant fxn to index() (preferring use of this one, to minimize confusion
    // between NodeID and NodeListIndex)
    NodeID get_node_id() const {
      return index();
    }

    node_value_type &value() {
      assert((graph_ptr_->node_id_map_).count(node_id_));
      const NodeListIndex node_list_index = (graph_ptr_->node_id_map_).at(node_id_);
      NodeInternal &node = (graph_ptr_->node_list_)[node_list_index];
      return node.node_value;
    }

    const node_value_type &value() const {
      assert((graph_ptr_->node_id_map_).count(node_id_));
      const NodeListIndex node_list_index = (graph_ptr_->node_id_map_).at(node_id_);
      const NodeInternal &node = (graph_ptr_->node_list_)[node_list_index];
      return node.node_value;
    }

    size_type degree() const {
      assert((graph_ptr_->node_id_map_).count(node_id_));
      const NodeListIndex node_list_index = (graph_ptr_->node_id_map_).at(node_id_);
      const NodeInternal &node = (graph_ptr_->node_list_)[node_list_index];
      return (node.adjacent_edges).size();
    }

    incident_iterator edge_begin() const {
      assert((graph_ptr_->node_id_map_).count(node_id_));
      const NodeListIndex node_list_index = (graph_ptr_->node_id_map_).at(node_id_);
      const NodeInternal &node = (graph_ptr_->node_list_)[node_list_index];
      return IncidentIterator((node.adjacent_edges).begin(), node_id_, graph_ptr_);
    }

    incident_iterator edge_end() const {
      assert((graph_ptr_->node_id_map_).count(node_id_));
      const NodeListIndex node_list_index = (graph_ptr_->node_id_map_).at(node_id_);
      const NodeInternal &node = (graph_ptr_->node_list_)[node_list_index];
      return IncidentIterator((node.adjacent_edges).end(), node_id_, graph_ptr_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      assert((graph_ptr_->node_id_map_).count(node_id_));
      return (graph_ptr_ == n.get_graph_pointer()) &&
              (node_id_ == n.get_node_id());
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
      assert((graph_ptr_->node_id_map_).count(node_id_));
      return (graph_ptr_ == n.get_graph_pointer()) &&
              (node_id_ < n.get_node_id());
    }

    // Access the pointer to the graph to which this node belongs.
    Graph *get_graph_pointer() const {
      return graph_ptr_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // ID of the node proxy. This is not necessarily the index in which
    // the corresponding NodeInternal exists
    NodeID node_id_;

    // Pointer to the graph to which this node belongs.
    Graph *graph_ptr_;

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_list_.size();
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

  Node add_node(const Point &position, const node_value_type &node_value = node_value_type()) {

    const NodeID new_node_id = node_list_.size();

    node_list_.emplace_back(new_node_id, position, node_value);

    node_id_map_[new_node_id] = new_node_id;

    return Node(new_node_id, this);
  }

  void update_node_value(const NodeID node_id, V node_value) {
    assert(node_id_map_.count(node_id));
    const NodeListIndex node_list_index = node_id_map_.at(node_id);
    NodeInternal &node = node_list_[node_list_index];
    node.node_value = node_value;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return n.get_graph_pointer() == this;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(node_id_map_.count(i));
    return Node(i, this);
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
    Edge() {
      edge_id_ = std::numeric_limits<EdgeID>::max();
      node1_id_ = std::numeric_limits<NodeID>::max();
      node2_id_ = std::numeric_limits<NodeID>::max();
      graph_ptr_ = nullptr;
    }

    Edge(EdgeID edge_id, NodeID node1_id, NodeID node2_id, const Graph *graph_ptr)
      : edge_id_(edge_id), node1_id_(node1_id), node2_id_(node2_id),
        graph_ptr_(const_cast<Graph*>(graph_ptr)) {
    }

    /** Return a node of this Edge */
    Node node1() const {
      assert((graph_ptr_->node_id_map_).count(node1_id_));
      return Node(node1_id_, graph_ptr_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      assert((graph_ptr_->node_id_map_).count(node2_id_));
      return Node(node2_id_, graph_ptr_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      assert((graph_ptr_->edge_id_map_).count(edge_id_));
      return (graph_ptr_ == e.get_graph_pointer()) &&
              (edge_id_ == e.get_edge_id());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      assert((graph_ptr_->edge_id_map_).count(edge_id_));
      return (graph_ptr_ == e.get_graph_pointer()) &&
              (edge_id_ < e.get_edge_id());
    }

    // Access the ID of the edge.
    // NOTE: this is EdgeID, not EdgeListID
    // Also NOTE: edge does not have to be valid for this fxn to be called
    EdgeID get_edge_id() const {
      return edge_id_;
    }

    // Access the pointer to the graph to which this edge belongs.
    Graph *get_graph_pointer() const {
      return graph_ptr_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // ID of the edge proxy. This is not necessarily the index in which
    // the corresponding EdgeInternal exists
    EdgeID edge_id_;

    // NOTE: these may not be in same order as stored in NodeInternal
    // (since node spawning incident iterator is set to be node1)
    NodeID node1_id_;
    NodeID node2_id_;

    // Pointer to the graph to which this node belongs.
    Graph *graph_ptr_;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_list_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(edge_id_map_.count(i));
    const EdgeListIndex edge_list_index = edge_id_map_.at(i);
    const EdgeInternal &edge = edge_list_[edge_list_index];
    return Edge(i, edge.node1_id, edge.node2_id, this);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    Edge edge = find_connecting_edge(a, b);
    const EdgeID edge_id = edge.get_edge_id();
    return (edge_id_map_.count(edge_id) && edge.get_graph_pointer() == this);
  }

  // Helper function to obtain the connecting edge.
  // Improves efficiency of add_edge() fxn. Otherwise would have to potentially
  // iterate through adjacent_edges twice.
  Edge find_connecting_edge(const Node &a, const Node &b) const {

    assert(node_id_map_.count(a.get_node_id()));
    const NodeListIndex node_list_index = node_id_map_.at(a.get_node_id());
    const NodeInternal &node = node_list_[node_list_index];

    for (const EdgeID edge_id : node.adjacent_edges) {
      assert(edge_id_map_.count(edge_id));
      const EdgeListIndex edge_list_index = edge_id_map_.at(edge_id);
      const EdgeInternal &edge = edge_list_[edge_list_index];
      if (edge.node1_id == b.get_node_id() || edge.node2_id == b.get_node_id()) {
        return Edge(edge_id, edge.node1_id, edge.node2_id, this);
      }
    }
    return Edge();
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
    Edge edge_proxy = find_connecting_edge(a, b);

    if (edge_id_map_.count(edge_proxy.get_edge_id()) &&
        edge_proxy.get_graph_pointer() == this) {
      return edge_proxy;
    } else {

      // create new EdgeID for this new edge
      const EdgeID new_edge_id = edge_list_.size();

      // add edge to graph
      NodeID id_a = a.get_node_id();
      NodeID id_b = b.get_node_id();
      edge_list_.emplace_back(new_edge_id, id_a, id_b);

      // add mapping from NodeID to NodeListIndex
      edge_id_map_[new_edge_id] = new_edge_id;

      // NOTE: do not specify these as const b/c need to add edges
      assert(node_id_map_.count(id_a));
      assert(node_id_map_.count(id_b));
      const NodeListIndex node_index1 = node_id_map_.at(id_a);
      const NodeListIndex node_index2 = node_id_map_.at(id_b);
      NodeInternal &node_a = node_list_[node_index1];
      NodeInternal &node_b = node_list_[node_index2];

      // add edge to nodes a and b
      (node_a.adjacent_edges).insert(new_edge_id);
      (node_b.adjacent_edges).insert(new_edge_id);

      // return edge proxy
      return Edge(new_edge_id, id_a, id_b, this);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    node_list_.clear();
    node_id_map_.clear();

    edge_list_.clear();
    edge_id_map_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {

   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    NodeIterator(typename NodeList::const_iterator node_list_iterator, const Graph *graph_ptr)
      : node_list_iterator_(node_list_iterator), graph_ptr_(const_cast<Graph*>(graph_ptr)) {
    }

    Node operator*() const {
      const NodeInternal &node = *node_list_iterator_;
      return Node(node.node_id, graph_ptr_);
    }

    NodeIterator& operator++() {
      node_list_iterator_++;
      return *this;
    }

    bool operator==(const NodeIterator &iterator) const {
      return (node_list_iterator_ == iterator.get_iterator()) &&
              (graph_ptr_ == iterator.get_graph_pointer());
    }

    typename NodeList::const_iterator get_iterator() const {
      return node_list_iterator_;
    }

    Graph *get_graph_pointer() const {
      return graph_ptr_;
    }

   private:
    friend class Graph;

    typename NodeList::const_iterator node_list_iterator_;

    Graph *graph_ptr_;

  };

  node_iterator node_begin() const {
    return NodeIterator(node_list_.begin(), this);
  }

  node_iterator node_end() const {
    return NodeIterator(node_list_.end(), this);
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    IncidentIterator(typename AdjacentEdges::const_iterator adj_edges_iterator,
                      const NodeID node_id, const Graph *graph_ptr)
      : adj_edges_iterator_(adj_edges_iterator), node_id_(node_id),
        graph_ptr_(const_cast<Graph*>(graph_ptr)) {
    }

    Edge operator*() const {
      const EdgeID edge_id = *adj_edges_iterator_;
      assert((graph_ptr_->edge_id_map_).count(edge_id));
      const EdgeListIndex edge_list_index = (graph_ptr_->edge_id_map_).at(edge_id);
      const EdgeInternal &edge = (graph_ptr_->edge_list_)[edge_list_index];

      // manipulate node1 and node2 so that node1 is node_id_
      if (node_id_ == edge.node1_id) {
        return Edge(edge_id, edge.node1_id, edge.node2_id, graph_ptr_);
      } else {
        return Edge(edge_id, edge.node2_id, edge.node1_id, graph_ptr_);
      }
    }

    IncidentIterator &operator++() {
      adj_edges_iterator_++;
      return *this;
    }

    bool operator==(const IncidentIterator &iit) const {
      return (adj_edges_iterator_ == iit.get_iterator()) &&
              (node_id_ == iit.get_node_id()) &&
              (graph_ptr_ == iit.get_graph_pointer());
    }

    typename AdjacentEdges::const_iterator get_iterator() const {
      return adj_edges_iterator_;
    }

    Graph *get_graph_pointer() const {
      return graph_ptr_;
    }

    NodeID get_node_id() const {
      return node_id_;
    }

   private:
    friend class Graph;

    typename AdjacentEdges::const_iterator adj_edges_iterator_;

    // ID of the node spawning the iterator (to ensure proper orientation)
    NodeID node_id_;

    Graph *graph_ptr_;
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    EdgeIterator(typename EdgeList::const_iterator edges_iterator, const Graph *graph_ptr)
      : edges_iterator_(edges_iterator), graph_ptr_(const_cast<Graph*>(graph_ptr)) {
    }

    Edge operator*() const {
      const EdgeInternal &edge = *edges_iterator_;
      return Edge(edge.edge_id, edge.node1_id, edge.node2_id, graph_ptr_);
    }

    EdgeIterator &operator++() {
      edges_iterator_++;
      return *this;
    }

    bool operator==(const EdgeIterator &eit) const {
      return (edges_iterator_ == eit.get_iterator()) &&
              (graph_ptr_ == eit.get_graph_pointer());
    }

    typename EdgeList::const_iterator get_iterator() const {
      return edges_iterator_;
    }

    Graph *get_graph_pointer() const {
      return graph_ptr_;
    }

   private:
    friend class Graph;

    typename EdgeList::const_iterator edges_iterator_;

    Graph *graph_ptr_;

  };

  edge_iterator edge_begin() const {
    return EdgeIterator(edge_list_.begin(), this);
  }

  edge_iterator edge_end() const {
    return EdgeIterator(edge_list_.end(), this);
  }

 private:

  // Container mapping Node IDs to index inside node_list_
  using NodeIDMap = std::unordered_map<NodeID, NodeListIndex>;

  // Container mapping Edge IDs to index inside edge_list_
  using EdgeIDMap = std::unordered_map<EdgeID, EdgeListIndex>;

  struct NodeInternal {
    NodeID node_id; // NOTE: not the same as position of node in node_list_
    Point position;
    AdjacentEdges adjacent_edges {};
    node_value_type node_value;

    NodeInternal(const NodeID node_id_, const Point &position_,
                  const node_value_type &node_value_)
      : node_id(node_id_), position(position_), node_value(node_value_) {
    }
  };

  struct EdgeInternal {
    EdgeID edge_id; // NOTE: not the same as position of edge in edge_list_
    NodeID node1_id;
    NodeID node2_id;
    EdgeInternal(const EdgeID edge_id_, const NodeID &node1_id_,
      const NodeID &node2_id_) : edge_id(edge_id_), node1_id(node1_id_), node2_id(node2_id_) {
    }
  };

  // Vector containing the true nodes.
  NodeList node_list_ {};

  // Container mapping Node IDs to index inside node_list_.
  NodeIDMap node_id_map_ {};

  // Vector containing the true edges.
  EdgeList edge_list_ {};

  // Container mapping Edge IDs to index inside edge_list_.
  EdgeIDMap edge_id_map_ {};

};

#endif // CME212_GRAPH_HPP
