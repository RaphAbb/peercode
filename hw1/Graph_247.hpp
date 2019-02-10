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

template <typename V>

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

  /** Type of value Node holds. */
  using node_value_type = V;

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
  class Node : private totally_ordered<Node> {
   public:
    /** @brief Construct an invalid node.*/
    Node() {
    }

    /** @brief the posiiton of this node
     *  @pre this Node is valid and is part of a graph
     *  @return a const Point object containing this node's position.
     */
    const Point& position() const {
      assert(graph_); // is *this a valid node?
      return graph_->nodes_[index_].position;
    }

    /** @brief the position of this node
     *  @pre this Node is valid and is part of a graph
     *  @return a Point object containing this node's position.
     */
    Point& position() {
      assert(graph_); // is *this a valid node?
      return graph_->nodes_[index_].position;
    }

    /** @brief the value of this node
     *  @pre this Node is valid and is part of a graph
     *  @return a mutable reference to this node's value.
     */
    node_value_type& value() {
      assert(graph_); // is *this a valid node?
      return graph_->nodes_[index_].value;
    }

    /** @brief the value of this node
     *  @pre this Node is valid and is part of a graph
     *  @return a const reference to this node's value.
     */
    const node_value_type& value() const {
      assert(graph_); // is *this a valid node?
      return graph_->nodes_[index_].value;
    }

    /** @brief the index of this Node
     *  @pre this Node is valid and is part of a graph
     *  @return this node's index, a number in the range [0, graph_size).
     */
    size_type index() const {
      assert(graph_); // is *this a valid node?
      return graph_->nodes_[index_].index;
    }

    /** @brief the degree of this Node
     *  @pre this Node is valid and is part of a graph
     *  @return the number of edges connected to this node
     *          this is equivalent to adjacency_vector.size()
     */
    size_type degree() const {
      assert(graph_); // is *this a valid node?
      return graph_->nodes_[index_].connections.size();
    }

    /** @brief an iterator to the beginning of this node's edges
     *  @pre this Node is vaild and is part of a graph
     *  @return an IncidentIterator to the beginning of this node's adjecency vector
     *          dereferences to Edge type
     */
    incident_iterator edge_begin() const {
      assert(graph_); // is *this a valid node?
      return incident_iterator(graph_, this, 0);
    }

    /** @brief an iterator to the end of this node's edges
     *  @pre this Node is valid and part of a graph
     *  @return an IncidentIterator to the one past the final edge of this node's
     *          adjacency vector. dereferences to Edge type
     */
    incident_iterator edge_end() const {
      assert(graph_);
      return incident_iterator(graph_, this, degree());
    }

    /** @brief Test whether this node and _n_ are equal.
     *  @pre this Node is valid, _n_ is a valid node in this graph
     *  @param[in] n Node to check equality with
     *  @return true if *this and _n_ are equal
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (!graph_ or !graph_->has_node(n)) return false; // are both valid?
      return this->index() == n.index();
    }

    /** @brief Test whether this node is less than _n_ in a global order.
     *  @pre this Node is valide, and _n_ is a valid node in this graph
     *  @param[in] n Node to compare this Node to
     *  @return true if *this < _n_ in a global order
     */
    bool operator<(const Node& n) const {
      assert(graph_); // is *this a valid node?
      assert(graph_->has_node(n)); // is n a valid node in this graph?
      return (this->index() < n.index());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    /** Pointer back to Graph */
    Graph* graph_ = nullptr;

    /** Index of node */
    size_type index_;

    /** @brief Construct a valid node
     *  @pre _index_ = _graph_.size()
     *  @param[in] graph a pointer to the parent graph
     *  @param[in] index _graph_.size() and index into graph's internal node data
     */
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
  Node add_node(const Point& pos, const node_value_type& val = node_value_type()) {
    // Initialize internal node struct
    size_type id = num_nodes();
    internal_node new_node;
    new_node.position = pos;
    new_node.value = val;
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {}

    /** @return this Edge's primary Node
     *  @pre this Edge is valid
     */
    const Node node1() const {
      return *node1_;
    }

    /** @return this Edge's secondary Node
     *  @pre this Edge is valid
     */
    const Node node2() const {
      return *node2_;
    }

    /** @brief Test whether this edge and _e_ are equal.
     *  @pre this Edge and _e_ are valid
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // Check if both nodes match
      return (node1() == e.node1() and node2() == e.node2())
          or (node2() == e.node1() and node1() == e.node2());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return index_ < e.index_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    /** Graph that this edge belongs to */
    Graph* graph_ = nullptr;

    /** Index */
    size_type index_;

    /** Node1 */
    const Node* node1_ = nullptr;

    /** Node2 */
    const Node* node2_ = nullptr;

    /** @brief Private constructor for a valid Edge
     *  @param[in] graph Pointer to parent graph
     *  @param[in] index Edge's index (has no meaning except within
     *                   the graph's internal data structure.
     *                   Equal to old num_edges() at time of construction
     *  @prarm[in] node1 Pointers to two nodes that this Edge
     *             node2 connects. Must be valid Nodes
     *
     *  Can only be called by the graph class through add_edge() or orient()
     */
    Edge(const Graph* graph, size_type index,
         const Node* node1, const Node* node2)
      : graph_(const_cast<Graph*>(graph)), index_(index),
        node1_(node1), node2_(node2) {
    }
  };

  /** @brief Return the total number of edges in the graph.
   *
   * Complexity: O(1)
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
    // Iterate through Node A's edges to see if it is connected with B
    for(auto e_it = a.edge_begin(); e_it != a.edge_end(); ++e_it) {
      if ((*e_it).node2() == b) return true;
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
    assert(a != b);                                  // distinct

    // Edge does not exist yet
    if(not has_edge(a,b)) {
      size_type id = this->num_edges();
      Edge new_edge(this,id,&a,&b);                      // Construct a valid Edge proxy
      edges_.emplace_back(a,b,new_edge);                 // Create the "real" internal struct
      nodes_[a.index()].connections.push_back(new_edge); // Add edge proxy to appropriate lists
      nodes_[b.index()].connections.push_back(new_edge);
      return new_edge;
    }

    // Edge exists - Iterate through Node A's edges to find it
    for(auto e_it = a.edge_begin(); e_it != a.edge_end(); ++e_it) {
      if ((*e_it).node2() == b) return *e_it;
    }

    return Edge();
  }

 private:

  /** @brief orients an edge _e_ such that the Node with _node1_index__ is _e_.node1()
   *  @pre _e_ is a valid Edge and has either node1 or node2 that has index _node_index__
   *  @param[in] e           The Edge to orient
   *  @param[in] node1_index The index representing the Node we wish to be e.node1()
   *  @return an Edge that is equal to _e_ with returned_edge.node1().index() = _node1_index_
   *
   *  @post old e.node1() == new e.node1()
   */
  Edge orient(const Edge& e, size_type node1_index) {
    if(e.node1().index() == node1_index) return e; // Edge already oriented correctly
    assert(e.node2().index() == node1_index);      // Make sure that edge is connected to node1_index
    Edge oriented_e(this,e.index_,e.node2_,e.node1_); // Create a new, flipped Edge (do not add to graph)
    return oriented_e;
  }

 public:

  /** @brief remove all nodes and edges from this graph.
   *  @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   *  @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** @brief Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    /** @brief dereferencing operator for a node_iterator
     *  @pre this iterator can be dereferenced, ie index_ < graph_->size_
     *  @return the Node that this node_iterator references
     */
    Node operator*() const {
      assert(index_ < graph_->size());
      return graph_->nodes_[index_].node;
    }

    /** @brief increment operator for a node_iterator
     *  @pre this iterator is not pointing to graph_->node_end(),
     *       ie index_ < graph_->size()
     *  @return a *this NodeIterator with its _index__ incremented
     */
    NodeIterator& operator++() {
      assert(index_ < graph_->size());
      index_++;
      return *this;
    }

    /** @brief the equality operator for a node_iterator
     *  @return true if these node iterators belong to the same graph and have
     *          the same index
     */
    bool operator==(const node_iterator& other_iterator) const {
      assert(graph_ == other_iterator.graph_);
      return index_ == other_iterator.index_;
    }

   private:
    friend class Graph;

    /** Pointer to parent Graph */
    Graph* graph_ = nullptr;

    /** Index of current internal_node */
    size_type index_;

    /** @brief Construct a valid NodeIterator
     *  @pre _index_ <= graph.size()
     *  @param[in] graph A pointer to the parent graph
     *  @param[in] index The index of the Node this Iterator references
     *
     *  Can only be called by the Graph class through node_begin() and node_end()
     */
    NodeIterator(const Graph* graph, size_type index)
      : graph_(const_cast<Graph*>(graph)), index_(index) { }
  };

  /** @brief Iterator to beginning of nodes_ vector
   *  @return  NodeIterator referencing the first Node in the graph
   *
   *  @post *(node_begin()) == graph.node(0)
   */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** @brief Iterator to end of nodes_ vector
   *  @return NodeIterator referencing one past the last Node in the graph
   *
   *  @post *(node_end()) is invalid, cannot be dereferenced
   */
  node_iterator node_end() const {
    return NodeIterator(this, this->size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   *  @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** @brief Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    /** @brief Dereferencing operator for IncidentIterator
     *  @pre   _index__ < (*node_).degree()
     *  @return the Edge that this IncidentIterator references
     */
    Edge operator*() const {
      auto incident_edges = this->incident_edges_();
      assert(index_ < incident_edges.size());
      Edge e = incident_edges[index_];
      return graph_->orient(e, node_index_);
    }

    /** @brief Increment operator for IncidentIterator
     *  @pre   _index__ < (*node_).degree()
     *  @return *this IncidentIterator with its _index__ incremented
     */
    IncidentIterator& operator++() {
      assert(index_ < this->incident_edges_().size());
      index_++;
      return *this;
    }

    /** @brief Equality operator for IncidentIterator
     *  @return true if *this and _other_it_ have the same Node and the same
     *          IncidentEdge index
     */
    bool operator==(const IncidentIterator& other_it) const {
      return (node_index_ == other_it.node_index_ &&
              index_ == other_it.index_);
    }

   private:
    friend class Graph;
    friend class Node;

    /** Pointer to parent Graph */
    Graph* graph_ = nullptr;

    /** Pointer to parent Node */
    size_type node_index_;

    /** Current index of incident edge */
    size_type index_;

    /** @brief return a reference to the vector of adjacent edges
     *  @pre the Node that _node_index_ references is valid,
     *       ie it is in a graph and its index < graph_->size()
     *  @return this Node's adjacency vector taken from the internal
     *          data structure
     *
     * This is used as a shortcut/helper to get a vector of the incident edges
     */
    const std::vector<edge_type>& incident_edges_() const {
      assert(graph_ && node_index_ < graph_->size());
      return graph_->nodes_[node_index_].connections;
    }

    /** @brief Construct a valid Incident Iterator
     *  @param[in] graph A pointer to the parent graph
     *  @param[in] node  A pointer to the parent node
     *  @param[in] index the index into the adjacency vector referenced by this IncidentIterator
     *  @pre the parent Node is a valid Node
     */
    IncidentIterator(const Graph* graph, const Node* node, size_type index)
      : graph_(const_cast<Graph*>(graph)), index_(index) {
        node_index_ = node->index();
      }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   *  @brief Iterator class for edges. A forward iterator. */
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

    /** @brief dereferencing operator for ae edge_iterator
     *  @pre this iterator can be dereferenced, ie index_ < graph_->num_edges()
     *  @return the Edge that this edge_iterator references
     */
    Edge operator*() const {
      assert(index_ < graph_->num_edges());
      return graph_->edges_[index_].edge;
    }

    /** @brief increment operator for an edge_iterator
     *  @pre this iterator is not pointing to graph_->edge_end(),
     *       ie index_ < graph_->num_edges()
     *  @return a *this EdgeIterator with its _index__ incremented
     */
    EdgeIterator& operator++() {
      assert(index_ < graph_->num_edges());
      index_++;
      return *this;
    }

    /** @brief the equality operator for n edge_iterator
     *  @return true if these edge iterators belong to the same graph and have
     *          the same index
     */
    bool operator==(const edge_iterator& other_iterator) const {
      assert(graph_ == other_iterator.graph_);
      return index_ == other_iterator.index_;
    }

   private:
    friend class Graph;

    /** Pointer to parent Graph */
    Graph* graph_ = nullptr;

    /** Index of current internal_node */
    size_type index_;

    /** @brief Construct a valid EdgeIterator
     *  @pre _index_ <= graph.num_edges()
     *  @param[in] graph A pointer to the parent graph
     *  @param[in] index The index of the Edge this Iterator references
      */
    EdgeIterator(const Graph* graph, size_type index)
      : graph_(const_cast<Graph*>(graph)), index_(index) { }

  };

  /** @brief Iterator to beginning of edges_ vector
   *  @return  EdgeIterator referencing the first Edge in the graph
   *
   *  @post *edge_begin() == edges_[0]
   */
  edge_iterator edge_begin() const {
    return edge_iterator(this, 0);
  }

  /** @brief Iterator to end of nodes_ vector
   *  @return  NodeIterator referencing one past the last Edge in the graph
   *
   *  @post *edge_end() is invalid, it cannot be derefereced
   */
  edge_iterator edge_end() const {
    return edge_iterator(this, num_edges());
  }

 private:

  /** Internal type for Node */
  struct internal_node {
    Point position;
    node_value_type value;
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
