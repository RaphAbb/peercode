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
template <typename V, typename E>
class Graph {
 private:

  struct NodeInternal;
  struct EdgeInternal;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V, E>;

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
  using edge_value_type = E;

  // Vector containing true nodes.
  using NodeList = std::vector<NodeInternal>;

  // Vector containing true edges.
  using EdgeList = std::vector<EdgeInternal>;

  // Stores the indices of the edges adjacent to the node.
  using AdjacentEdges = std::set<EdgeID>;

  // Container mapping Node IDs to index inside node_list_
  using NodeIDMap = std::unordered_map<NodeID, NodeListIndex>;

  // Container mapping Edge IDs to index inside edge_list_
  using EdgeIDMap = std::unordered_map<EdgeID, EdgeListIndex>;

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

    // Return modifiable, non-const reference to position
    Point &position() {
      assert(isValidNode());
      const NodeListIndex idx = (graph_ptr_->node_id_map_).at(node_id_);
      return (graph_ptr_->node_list_)[idx].position;
    }

    /** Return this node's position. */
    const Point &position() const {
      // node is valid if index is within size of the node_list_
      assert(isValidNode());
      const NodeListIndex idx = (graph_ptr_->node_id_map_).at(node_id_);
      return (graph_ptr_->node_list_)[idx].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // NOTE: this function returns a NodeID, not a NodeListIndex
      return node_id_;
    }

    // Redundant fxn to index() (preferring use of this one, to minimize confusion
    // between NodeID and NodeListIndex)
    NodeID get_node_id() const {
      return index();
    }

    node_value_type &value() {
      assert(isValidNode());
      const NodeListIndex idx = (graph_ptr_->node_id_map_).at(node_id_);
      NodeInternal &node = (graph_ptr_->node_list_)[idx];
      return node.node_value;
    }

    const node_value_type &value() const {
      assert(isValidNode());
      const NodeListIndex idx = (graph_ptr_->node_id_map_).at(node_id_);
      const NodeInternal &node = (graph_ptr_->node_list_)[idx];
      return node.node_value;
    }

    size_type degree() const {
      assert(isValidNode());
      const NodeListIndex idx = (graph_ptr_->node_id_map_).at(node_id_);
      const NodeInternal &node = (graph_ptr_->node_list_)[idx];
      return (node.adjacent_edges).size();
    }

    incident_iterator edge_begin() const {
      assert(isValidNode());
      const NodeListIndex idx = (graph_ptr_->node_id_map_).at(node_id_);
      const NodeInternal &node = (graph_ptr_->node_list_)[idx];
      return IncidentIterator((node.adjacent_edges).begin(), node_id_, graph_ptr_);
    }

    incident_iterator edge_end() const {
      assert(isValidNode());
      const NodeListIndex idx = (graph_ptr_->node_id_map_).at(node_id_);
      const NodeInternal &node = (graph_ptr_->node_list_)[idx];
      return IncidentIterator((node.adjacent_edges).end(), node_id_, graph_ptr_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      assert(isValidNode());
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
      assert(isValidNode());
      // what about ordering between two diff graphs but same node ID?
      return node_id_ < n.get_node_id();
    }

    bool isValidNode() const {
      return (graph_ptr_->node_id_map_).count(node_id_);
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
    return node_id_map_.size();
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

    const NodeID new_node_id = node_id_map_.size();
    const NodeListIndex new_node_idx = node_list_.size();

    node_list_.emplace_back(new_node_id, position, node_value);
    node_id_map_[new_node_id] = new_node_idx;

    return Node(new_node_id, this);
  }

  /** Remove node from the graph by deleting the node mapping from node_id_map_.
   *  Note that NodeInternal is not removed, only its NodeID --> NodeListIndex mapping.
   *
   *  @param[in]  n The Node to be removed from the graph.
   *  @result     Return 0 if no node found. Return 1 if node found and removed.
   *
   *  @pre
   *  @post       new node_id_map_.size() == old node_id_map_.size() - 1
   *              Max NodeID in node_id_map_ == new node_id_map_.size() - 1
   *              There exists a unique mapping for all IDs in the range
   *              [0, new node_id_map_.size()) for all valid nodes in the graph.
   */
   size_type remove_node(const Node &n) {

     // case where node not found
     if (!has_node(n)) {
       return 0;
     }

     const NodeID max_node_id = node_id_map_.size() - 1; // last ID in mapping

     // retrieve the node to be deleted
     const NodeID node_id = n.get_node_id();
     NodeInternal &node = node_list_[node_id_map_.at(node_id)];
     node.node_id = std::numeric_limits<NodeID>::max(); // ID set to out of range for safe measure

     // problem using incident iterator (iterating over set that I am deleting from)
     // instead, iterate through copy of adj edges and delete the edges
     AdjacentEdges adj_edges_copy = node_list_[node_id_map_.at(n.get_node_id())].adjacent_edges;
     for (EdgeID edge_id : adj_edges_copy) {
       Edge bad_edge = edge(edge_id);
       remove_edge(bad_edge);
     }

     // clear the actual adjacent edges set
     node_list_[node_id_map_.at(n.get_node_id())].adjacent_edges.clear();

     NodeListIndex &node_idx = node_id_map_.at(node_id);
     NodeListIndex &last_idx = node_id_map_.at(max_node_id);
     AdjacentEdges &adj_edges = node_list_[last_idx].adjacent_edges;

     // only update node info if not deleting last node
     if (max_node_id != node_id) {
       node_list_[last_idx].node_id = node_id; // update node ID of swapped (last) node
       update_swapped_node_edges(adj_edges, node_id, max_node_id);
     }

     // swap and pop
     std::swap(node_idx, last_idx);
     node_id_map_.erase(max_node_id);

     return 1;
   }

   // helper to update the node ID in the node swapped in for the deleted one
   void update_swapped_node_edges(AdjacentEdges &adj_edges, const NodeID node_id, const NodeID max_node_id) {
     for (auto it = adj_edges.begin(); it != adj_edges.end(); ++it) {
       EdgeInternal &edge = edge_list_[edge_id_map_.at(*it)];
       if (edge.node1_id == max_node_id) {
         edge.node1_id = node_id;
       } else if (edge.node2_id == max_node_id) {
         edge.node2_id = node_id;
       } else {
         assert(false); // catch error if adj edge does not have the right node
       }
     }
   }

   /** Remove node from the graph by deleting the node mapping from node_id_map_.
    *  Note that NodeInternal is not removed, only its NodeID --> NodeListIndex mapping.
    *
    *  @param[in]  n_it Iterator to the Node to be removed from the graph.
    *  @result     Return the same iterator, which will now point to what was
    *              previously the last element in old node_id_map_
    *
    *  @pre        has_node(*n_it) == true, i.e. the node must belong to this graph
    *  @post       new node_id_map_.size() == old node_id_map_.size() - 1
    *              Max NodeID in node_id_map_ == new node_id_map_.size() - 1
    *              There exists a unique mapping for all IDs in the range
    *              [0, new node_id_map_.size()) for all valid nodes in the graph.
    */
   node_iterator remove_node(node_iterator n_it) {
     const Node &n = *n_it;
     remove_node(n);
     return n_it;
   }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // if mapping doesn't exist, we've deleted the node or Node is invalid
    return node_id_map_.count(n.get_node_id());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if (node_id_map_.count(i)) {
      return Node(i, this);
    } else {
      return Node();
    }
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
      assert(isValidEdge());
      return Node(node1_id_, graph_ptr_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      assert(isValidEdge());
      return Node(node2_id_, graph_ptr_);
    }

    edge_value_type &value() {
      assert(isValidEdge());
      const EdgeListIndex idx = (graph_ptr_->edge_id_map_).at(edge_id_);
      return (graph_ptr_->edge_list_)[idx].edge_value;
    }

    const edge_value_type &value() const {
      assert(isValidEdge());
      const EdgeListIndex idx = (graph_ptr_->edge_id_map_).at(edge_id_);
      return (graph_ptr_->edge_list_)[idx].edge_value;
    }


    /** Obtain Eucidean length of edge */
    double length() const {
      assert(isValidEdge());
      const EdgeInternal &edge = (graph_ptr_->edge_list_)[(graph_ptr_->edge_id_map_).at(edge_id_)];
      const Point pos1 = (graph_ptr_->node_list_)[(graph_ptr_->node_id_map_).at(edge.node1_id)].position;
      const Point pos2 = (graph_ptr_->node_list_)[(graph_ptr_->node_id_map_).at(edge.node2_id)].position;
      return norm_2(pos1 - pos2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      assert(isValidEdge());
      return (graph_ptr_ == e.get_graph_pointer()) &&
              (edge_id_ == e.get_edge_id());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      assert(isValidEdge());
      return (edge_id_ < e.get_edge_id() || graph_ptr_ < e.get_graph_pointer());
    }

    // function for checking representation invariants
    bool isValidEdge() const {

      // check EdgeID --> EdgeListIndex mapping
      if (!(graph_ptr_->edge_id_map_).count(edge_id_)) {
        return false;
      }

      // check NodeID --> NodeListIndex mapping
      if (!(graph_ptr_->node_id_map_).count(node1_id_) || !(graph_ptr_->node_id_map_).count(node2_id_)) {
        return false;
      }
      return true;
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
    return edge_id_map_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    if (edge_id_map_.count(i)) {
      const EdgeListIndex idx = edge_id_map_.at(i);
      const EdgeInternal &edge = edge_list_[idx];
      return Edge(i, edge.node1_id, edge.node2_id, this);
    }
    else {
      return Edge();
    }
  }

  /** Remove edge from the graph, with an Edge as input.
   *  Note that EdgeInternal is not removed, only its EdgeID --> EdgeListIndex mapping.
   *
   *  @param[in]  edge The edge to be removed from the graph.
   *  @result     Return 0 if no edge found. Return 1 if edge found and removed.
   *
   *  @pre
   *  @post       new edge_id_map_.size() == old edge_id_map_.size() - 1
   *              Max EdgeID in edge_id_map_ == new edge_id_map_.size() - 1
   *              There exists a unique mapping for all IDs in the range
   *              [0, new edge_id_map_.size()) for all valid edges in the graph.
   */
  size_type remove_edge(const Edge &e) {

    const NodeID edge_id = e.get_edge_id();

    if (!edge_id_map_.count(edge_id)) {
      return 0;
    }

    const NodeID max_edge_id = edge_id_map_.size() - 1; // highest ID in mapping

    // retrieve edge to be deleted
    EdgeInternal &edge = edge_list_[edge_id_map_.at(edge_id)];
    edge.edge_id = std::numeric_limits<NodeID>::max(); // ID set to out of range for safe measure

    // delete this edge_id from the adjacent edges of both nodes
    NodeInternal &node_a = node_list_[node_id_map_.at(edge.node1_id)];
    NodeInternal &node_b = node_list_[node_id_map_.at(edge.node2_id)];
    (node_a.adjacent_edges).erase(edge_id);
    (node_b.adjacent_edges).erase(edge_id);

    EdgeListIndex &edge_idx = edge_id_map_.at(edge_id);
    EdgeListIndex &last_idx = edge_id_map_.at(max_edge_id);

    EdgeInternal &last_edge = edge_list_[last_idx];
    AdjacentEdges &adj_edges1 = node_list_[node_id_map_.at(last_edge.node1_id)].adjacent_edges;
    AdjacentEdges &adj_edges2 = node_list_[node_id_map_.at(last_edge.node2_id)].adjacent_edges;

    // also update the edge ID where it appears in adjacency list
    if (max_edge_id != edge_id) {
      last_edge.edge_id = edge_id; // update ID of edge swapped in for deleted one
      adj_edges1.erase(max_edge_id);
      adj_edges1.insert(edge_id);
      adj_edges2.erase(max_edge_id);
      adj_edges2.insert(edge_id);
    }

    // swap and pop
    std::swap(edge_idx, last_idx);
    edge_id_map_.erase(max_edge_id);

    return 1;
  }

  /** Remove edge from the graph, with its two connected nodes as input.
   *  Note that EdgeInternal is not removed, only its EdgeID --> EdgeListIndex mapping.
   *
   *  @param[in]  a The node at one end of the edge to be removed from the graph.
   *  @param[in]  b The node at the other end of the edge to be removed.
   *  @result     Return 0 if no edge found. Return 1 if edge found and removed.
   *
   *  @pre
   *  @post       new edge_id_map_.size() == old edge_id_map_.size() - 1
   *              Max EdgeID in edge_id_map_ == new edge_id_map_.size() - 1
   *              There exists a unique mapping for all IDs in the range
   *              [0, new edge_id_map_.size()) for all valid edges in the graph.
   */
  size_type remove_edge(const Node &a, const Node &b) {
     Edge edge = find_connecting_edge(a, b);
     return remove_edge(edge);
   }

  /** Remove edge from the graph, with an iterator to the edge as input.
   *  Note that EdgeInternal is not removed, only its EdgeID --> EdgeListIndex mapping.
   *
   *  @param[in]  e_it Iterator to the Edge to be removed from the graph.
   *  @result     Return the same iterator, which will now point to what was
   *              previously the last element in old node_id_map_.
   *
   *  @pre        has_edge(*e_it) == true, i.e. the edge must belong to this graph.
   *  @post       new edge_id_map_.size() == old edge_id_map_.size() - 1
   *              Max EdgeID in edge_id_map_ == new edge_id_map_.size() - 1
   *              There exists a unique mapping for all IDs in the range
   *              [0, new edge_id_map_.size()) for all valid edges in the graph.
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    const Edge &e = *e_it;
    remove_edge(e);
    return e_it;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    Edge edge = find_connecting_edge(a, b);

    if (!edge_id_map_.count(edge.get_edge_id())) {
      return false;
    }
    return true;
  }

  // Helper function to obtain the connecting edge.
  // Improves efficiency of add_edge() fxn. Otherwise would have to potentially
  // iterate through adjacent_edges twice.
  Edge find_connecting_edge(const Node &a, const Node &b) const {

    // no valid edge if either node is invalid
    if (!node_id_map_.count(a.get_node_id()) || !node_id_map_.count(b.get_node_id())) {
      return Edge();
    }

    // no valid edge if a and b are the same node
    if (a.get_node_id() == b.get_node_id()) {
      return Edge();
    }

    const NodeInternal &node = node_list_[node_id_map_.at(a.get_node_id())];

    for (const EdgeID edge_id : node.adjacent_edges) {
      if (!edge_id_map_.count(edge_id)) {
        continue; // pass if cannot find edge
      }
      const EdgeInternal &edge = edge_list_[edge_id_map_.at(edge_id)];
      // edge is only valid if node1_id and node2_id are different
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type &edge_value = edge_value_type()) {
    Edge edge_proxy = find_connecting_edge(a, b);

    if (edge_id_map_.count(edge_proxy.get_edge_id()) &&
        edge_proxy.get_graph_pointer() == this) {
      return edge_proxy;
    } else {

      // create new EdgeID for this new edge
      const EdgeID new_edge_id = edge_id_map_.size();
      const EdgeListIndex new_edge_idx = edge_list_.size();

      // add edge to graph
      NodeID id_a = a.get_node_id();
      NodeID id_b = b.get_node_id();
      edge_list_.emplace_back(new_edge_id, id_a, id_b, edge_value);

      // add mapping from NodeID to NodeListIndex
      edge_id_map_[new_edge_id] = new_edge_idx;

      assert(node_id_map_.count(id_a));
      assert(node_id_map_.count(id_b));
      NodeInternal &node_a = node_list_[node_id_map_.at(id_a)];
      NodeInternal &node_b = node_list_[node_id_map_.at(id_b)];

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

    NodeIterator(typename NodeIDMap::const_iterator node_map_iterator, const Graph *graph_ptr)
      : node_map_iterator_(node_map_iterator), graph_ptr_(const_cast<Graph*>(graph_ptr)) {
    }

    Node operator*() const {
      return Node(node_map_iterator_->first, graph_ptr_);
    }

    NodeIterator& operator++() {
      ++node_map_iterator_;
      return *this;
    }

    bool operator==(const NodeIterator &iterator) const {
      return (node_map_iterator_ == iterator.get_iterator()) &&
              (graph_ptr_ == iterator.get_graph_pointer());
    }

    typename NodeIDMap::const_iterator get_iterator() const {
      return node_map_iterator_;
    }

    Graph *get_graph_pointer() const {
      return graph_ptr_;
    }

   private:
    friend class Graph;

    typename NodeIDMap::const_iterator node_map_iterator_;

    Graph *graph_ptr_;

  };

  node_iterator node_begin() const {
    return NodeIterator(node_id_map_.begin(), this);
  }

  node_iterator node_end() const {
    return NodeIterator(node_id_map_.end(), this);
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
      ++adj_edges_iterator_;
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

    EdgeIterator(typename EdgeIDMap::const_iterator edge_map_iterator, const Graph *graph_ptr)
      : edge_map_iterator_(edge_map_iterator), graph_ptr_(const_cast<Graph*>(graph_ptr)) {
    }

    Edge operator*() const {
      const EdgeInternal &edge = (graph_ptr_->edge_list_)[edge_map_iterator_->second];
      return Edge(edge.edge_id, edge.node1_id, edge.node2_id, graph_ptr_);
    }

    EdgeIterator &operator++() {
      ++edge_map_iterator_;
      return *this;
    }

    bool operator==(const EdgeIterator &eit) const {
      return (edge_map_iterator_ == eit.get_iterator()) &&
              (graph_ptr_ == eit.get_graph_pointer());
    }

    typename EdgeIDMap::const_iterator get_iterator() const {
      return edge_map_iterator_;
    }

    Graph *get_graph_pointer() const {
      return graph_ptr_;
    }

   private:
    friend class Graph;

    typename EdgeIDMap::const_iterator edge_map_iterator_;

    Graph *graph_ptr_;

  };

  edge_iterator edge_begin() const {
    return EdgeIterator(edge_id_map_.begin(), this);
  }

  edge_iterator edge_end() const {
    return EdgeIterator(edge_id_map_.end(), this);
  }

 private:

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
    edge_value_type edge_value;
    EdgeInternal(const EdgeID edge_id_, const NodeID &node1_id_,
      const NodeID &node2_id_, const edge_value_type &edge_value_)
      : edge_id(edge_id_), node1_id(node1_id_),
        node2_id(node2_id_), edge_value(edge_value_) {
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
