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

template <typename V, typename E> // Graph is now a template for a certain value type V
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

  /** Type of this graph's values. */
  using node_value_type = V;

  /** Type of this graph's values. */
  using edge_value_type = E;

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

  /** Construct an empty graph, initializing node and edge vectors. */
  Graph()
    : nodes_(), i2u_(), edges_() {
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
  class Node : private totally_ordered<Node>  {
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

    /** @brief the posiiton of this node
     *  @pre this Node is valid and is part of a graph
     *  @return a const Point object containing this node's position.
     */
    const Point& position() const {
      assert(this->is_valid());
      return graph_->nodes_[uid_].position;
    }

    /** @brief the position of this node
     *  @pre this Node is valid and is part of a graph
     *  @return a Point object containing this node's position.
     */
    Point& position() {
      assert(this->is_valid());
      return graph_->nodes_[uid_].position;
    }

    /** @brief the value of this node
     *  @pre this Node is valid and is part of a graph
     *  @return a mutable reference to this node's value.
     */
    node_value_type& value() {
      assert(this->is_valid());
      return graph_->nodes_[uid_].value;
    }

    /** @brief the value of this node
     *  @pre this Node is valid and is part of a graph
     *  @return a const reference to this node's value.
     */
    const node_value_type& value() const {
      assert(this->is_valid());
      return graph_->nodes_[uid_].value;
    }

    /** @brief the index of this Node
     *  @pre this Node is valid and is part of a graph
     *  @return this node's index, a number in the range [0, graph_size).
     */
    size_type index() const {
      assert(this->is_valid());
      return graph_->nodes_[uid_].index;
    }

    /** @brief the degree of this Node
     *  @pre this Node is valid and is part of a graph
     *  @return the number of edges connected to this node
     *          this is equivalent to adjacency_vector.size()
     */
    size_type degree() const {
      assert(this->is_valid());
      return graph_->nodes_[uid_].connections.size();
    }

    /** @brief an iterator to the beginning of this node's edges
     *  @pre this Node is vaild and is part of a graph
     *  @return an IncidentIterator to the beginning of this node's adjecency vector
     *          dereferences to Edge type
     */
    incident_iterator edge_begin() const {
      assert(this->is_valid());
      return incident_iterator(graph_, this, 0);
    }

    /** @brief an iterator to the end of this node's edges
     *  @pre this Node is valid and part of a graph
     *  @return an IncidentIterator to the one past the final edge of this node's
     *          adjacency vector. dereferences to Edge type
     */
    incident_iterator edge_end() const {
      assert(this->is_valid());
      return incident_iterator(graph_, this, this->degree());
    }

    /** @brief Test whether this node and _n_ are equal.
     *  @pre this Node is valid, _n_ is a valid node in this graph
     *  @param[in] n Node to check equality with
     *  @return true if *this and _n_ are equal
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      assert(this->is_valid() and n.is_valid());
      return (graph_->has_node(n) and this->index() == n.index());
    }

    /** @brief Test whether this node is less than _n_ in a global order.
     *  @pre this Node is valide, and _n_ is a valid node in this graph
     *  @param[in] n Node to compare this Node to
     *  @return true if *this < _n_ in a global order
     */
    bool operator<(const Node& n) const {
      assert(this->is_valid() and n.is_valid());
      return (this->index() < n.index());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    /** Pointer back to Graph */
    graph_type* graph_ = nullptr;

    /** Unique ID of node */
    size_type uid_;

    /** @brief Construct a valid node
     *  @pre _uid_ = _graph_.size()
     *  @param[in] graph a pointer to the parent graph
     *  @param[in] uid _graph_.size() a unique ID into graph's internal node data
     */
    Node(const graph_type* graph, size_type uid)
      : graph_(const_cast<graph_type*>(graph)), uid_(uid) {
    }

    bool is_valid() const {
      return uid_>= 0 && uid_< graph_->nodes_.size()
          && graph_->nodes_[uid_].index < graph_->i2u_.size()
          && graph_->i2u_[graph_->nodes_[uid_].index] == uid_;
    }

  };

  /** @brief Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return i2u_.size();
  }

  /** @brief Synonym for size. */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] pos The new node's position
   * @param[in] val The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   * Parameter _val_ is default initialized if not provided
   */
  Node add_node(const Point& pos, const node_value_type& val = node_value_type () ) {

    size_type idx = i2u_.size();    // Next index in graph
    size_type uid = nodes_.size();  // Next uid in node container

    // Initialize internal node struct
    internal_node new_node;
    new_node.position = pos;
    new_node.value    = val;
    new_node.index    = idx;
    new_node.node     = Node(this, uid);

    // Add to vector of nodes and i2u containers
    i2u_.push_back(uid);
    nodes_.push_back(new_node);

    // return proxy
    return new_node.node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if _n_ is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (this == n.graph_ and n.is_valid());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i >= 0 and i < i2u_.size());
    return nodes_[i2u_[i]].node;
  }

  /** @brief Remove a Node _n_ from the graph, and all of its incident edges,
   *            returning a flag for operation success
   *  @param[in] n The Node to be removed
   *
   *  @post Node _n_ is not in the graph
   *  @post All edges incident to Node _n_ are not in the graph
   *  @post If Node _n_ was in the graph, then new num_nodes() = old num_nodes() - 1
   *  @post If Node _n_ was in the graph, then new num_edges() = old num_edges() - degree(n)
   *  @post If Node _n_ was in the graph, 
   *            then for all other nodes with old index _n_.index() < i < num_nodes(),
   *            new index() = old index() - 1
   *  @post If Node _n_ was in the graph, all NodeIterator's are invalidated
   *  @post If Node _n_ was in the graph, and degree(n) > 0, all EdgeIterator's are invalidated
   *  @post If Node _n_ was in the graph, and degree(n) > 0, 
   *            all IncidentIterator's for Node _n_ and connected nodes are invalidated
   *
   *  @return Binary flag 0 or 1, 
   *              0 indicating Node _n_ was not in the graph,
   *              1 indicating Node _n_ was successfully removed from the graph
   *
   * Complexity: O(k^2 + n) operations, where k is the degree of Node _n_, 
   *                and n is the number of nodes in the graph
   *
   */
  size_type remove_node(const Node& n) {

    if (has_node(n)) {
      // Go through adjacency list, delete edges of neighboring nodes
      auto n_edge_list = &(nodes_[i2u_[n.index()]].connections);
      for (auto niit = n_edge_list->begin(); niit != n_edge_list->end(); ++niit) {
        
        // Delete from attached node's adjacency lists
        auto n2 = niit->node2(); // Neighboring node is always the second one
        auto n2_edge_list = &(nodes_[i2u_[n2.index()]].connections);
        for (auto n2iit = n2_edge_list->begin(); n2iit != n2_edge_list->end(); ++n2iit) {
          if (n2iit->node2() == n) {
            n2_edge_list->erase(n2iit);
            break;        
          }
        }

        // Set the true edge container to be invalid
        edges_[niit->uid_].valid = false; 
      }
      n_edge_list->clear(); // Clear this node's adjacency list

      // Delete node and fix numbering
      auto it = i2u_.begin() + n.index();
      it = i2u_.erase(it);
      while (it != i2u_.end()) {
        nodes_[*it].index -= 1;
        ++it;
      }

      return 1;
    } else {  
      return 0;
    }
  }

  /** @brief Remove a Node referenced by *_nit_ and its incident edges, returning an iterator to the next node
   *  @param[in] nit An iterator referring to Node _n_ to be removed, *_nit_ == _n_
   *
   *  @post Node _n_ is not in the graph
   *  @post All edges incident to Node _n_ are not in the graph
   *  @post If Node _n_ was in the graph, then new num_nodes() = old num_nodes() - 1
   *  @post If Node _n_ was in the graph, then new num_edges() = old num_edges() - degree(n)
   *  @post If Node _n_ was in the graph, 
   *            then for all other nodes with old index _n_.index() < i < num_nodes(),
   *            new index() = old index() - 1
   *  @post If Node _n_ was in the graph, all NodeIterator's are invalidated
   *  @post If Node _n_ was in the graph, and degree(n) > 0, all EdgeIterator's are invalidated
   *  @post If Node _n_ was in the graph, and degree(n) > 0, 
   *            all IncidentIterator's for Node _n_ and connected nodes are invalidated
   *
   *  @return Binary flag 0 or 1, 
   *              0 indicating Node _n_ was not in the graph,
   *              1 indicating Node _n_ was successfully removed from the graph
   *
   * Complexity: O(k^2 + n) operations, where k is the degree of Node _n_, and n is the number of nodes in the graph
   *
   */
  node_iterator remove_node(node_iterator nit) {
    auto n = *nit;
    auto flag = this->remove_node(n);

    auto new_nit = this->node_begin();
    if (flag) { // If node removal occured, advance a new iterator to next position
      size_type count = 0;
      while (count < n.index()) {
        ++new_nit;
        ++count;
      }
    } else {   // If node removal did not occur, place new iterator at next position
      new_nit = ++nit;
    }

    return new_nit;
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
    /** @brief Construct an invalid Edge. */
    Edge() {}

    /** @return this Edge's primary Node
     *  @pre this Edge is valid
     */
    Node node1() const {
      assert(this->is_valid());
      if (swap_) {
        return graph_->edges_[uid_].node2;
      } else {
        return graph_->edges_[uid_].node1;
      }
    }

    /** @return this Edge's secondary Node
     *  @pre this Edge is valid
     */
    Node node2() const {
      assert(this->is_valid());
      if (swap_) {
        return graph_->edges_[uid_].node1;
      } else {
        return graph_->edges_[uid_].node2;
      }
    }

    /** Returns the distance between connected points */
    double length() const {
      return norm(this->node1().position() - this->node2().position());
    }

    /** @brief the value of this node
     *  @pre this Edge is valid and is part of a graph
     *  @return a mutable reference to this edge's value.
     */
    edge_value_type& value() {
      assert(this->is_valid());
      return graph_->edges_[uid_].value;
    }

    /** @brief the value of this node
     *  @pre this Edge is valid and is part of a graph
     *  @return a constant reference to this edge's value.
     */
    const edge_value_type& value() const {
      assert(this->is_valid());
      return graph_->edges_[uid_].value;
    }

    /** @brief Test whether this edge and _e_ are equal.
     *  @pre this Edge and _e_ are valid
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      node_type n1 = graph_->edges_[uid_].node1;
      node_type n2 = graph_->edges_[uid_].node2;
      // Check if both nodes match
      return (n1 == e.node1() and n2 == e.node2())
          or (n2 == e.node1() and n1 == e.node2());
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
    graph_type* graph_ = nullptr;

    /** Unique ID in graph container */
    size_type uid_;

    /** Flag for swap return of node 1 and 2 */
    bool swap_ = false;

    /** @brief Private constructor for a valid Edge
     *  @param[in] graph Pointer to parent graph
     *  @param[in] uid Edge's unique ID (has no meaning except within
     *                   the graph's internal data structure.
     *                   Equal to old num_edges() at time of construction
     *  @param[in] swap boolean flag for whether to swap return value of node1() and node2()
     *
     *  Can only be called by the graph class through add_edge()
     */
    Edge(graph_type* graph, size_type uid, bool swap = false)
      : graph_(const_cast<graph_type*>(graph)), uid_(uid), swap_(swap) {
    }

    /** Return valid flag for underlying node data */
    bool is_valid() const {
      return graph_->edges_[uid_].valid;
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    size_type count = 0;
    auto it = this->edge_begin();
    while (it != this->edge_end()) {
      ++it;
      ++count;
    }
    return count;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i >= 0 and i < this->num_edges());
    size_type count = 0;
    auto it = this->edge_begin();
    while (count < i) {
      ++it;
      ++count;
    }
    return *it;
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
    assert(this->has_node(a) and this->has_node(b));   // valid nodes of the graph
    assert(not (a == b));                              // distinct

    // Edge does not exist yet
    if(not has_edge(a,b)) {
      size_type uid = edges_.size(); // Next uid in edge container

      // Construct valid Edge proxies
      // Two edge proxies to ensure that each node connection has an Edge that returns itself as node1()
      Edge new_edge(this, uid);  
      Edge new_edge_swap(this, uid, true);

      edges_.emplace_back(a, b, new_edge);                 // Create the "real" internal struct
      nodes_[i2u_[a.index()]].connections.push_back(new_edge);    // Add edge proxy to appropriate lists
      nodes_[i2u_[b.index()]].connections.push_back(new_edge_swap);

      return new_edge;

    }

    // Edge exists - Iterate through Node A's edges to find it
    for(auto e_it = a.edge_begin(); e_it != a.edge_end(); ++e_it) {
      if ((*e_it).node2() == b) return *e_it;
    }

    return Edge();
  }

  /** @brief Remove an Edge _e_ from the graph
   *  @param[in] a The first Node to be removed
   *  @param[in] b The second Node to be removed
   *
   *  @post Edge _e_ is not in the graph
   *  @post If Edge _e_ was in the graph, then new num_edges() = old num_edges() - 1
   *  @post If Edge _e_ was in the graph, all EdgeIterator's are invalidated
   *  @post If Edge _e_ was in the graph, 
   *            NodeIterator's for Node _a_ and _b_ are invalidated
   *
   *  @return Binary flag 0 or 1, 
   *              0 indicating Edge _e_ was not in the graph,
   *              1 indicating Edge _e_ was successfully removed from the graph
   *
   * Complexity: O(k) operations, where k is the maximum degree of Node _a_ or _b_
   *
   */
  size_type remove_edge(const Node& a, const Node& b) {

    if (has_edge(a,b)) {
      // Edge proxy that is being deleted
      edge_type del_edge;

      // Delete from the adjacency list of node a
      auto a_edge_list = &(nodes_[i2u_[a.index()]].connections);
      auto ait = a_edge_list->begin();
      while (ait != a_edge_list->end()) {
        if (ait->node1() == b or ait->node2() == b) {
          // Delete the edge and end the loop
          del_edge = *ait;
          a_edge_list->erase(ait);
          break;
        } else {
          ++ait;
        }
      }

      // Delete from the adjacency list of node b
      auto b_edge_list = &(nodes_[i2u_[b.index()]].connections);
      auto bit = b_edge_list->begin();
      while (bit != b_edge_list->end()) {
        if (bit->node1() == a or bit->node2() == a) {
          // Delete the edge and end the loop
          b_edge_list->erase(bit);
          break;
        } else {
          ++bit;
        }
      }

      // Set the internal edge to an invalid status
      this->edges_[del_edge.uid_].valid = false;
      return 1;

    } else {
      return 0;
    }
  }

  /** @brief Remove an Edge _e_ from the graph
   *  @param[in] _e_ The Edge to be removed
   *
   *  @post Edge _e_ is not in the graph
   *  @post If Edge _e_ was in the graph, then new num_edges() = old num_edges() - 1
   *  @post If Edge _e_ was in the graph, all EdgeIterator's are invalidated
   *  @post If Edge _e_ was in the graph, 
   *            For the two connected nodes by _e_, Node _a_ and _b_,
   *            NodeIterator's for _a_ and _b_ are invalidated
   *
   *  @return Binary flag 0 or 1, 
   *              0 indicating Edge _e_ was not in the graph,
   *              1 indicating Edge _e_ was successfully removed from the graph
   *
   * Complexity: O(k) operations, where k is the maximum degree of Node _a_ or _b_
   *
   */
  size_type remove_edge(const Edge& e) {
    auto flag = remove_edge(e.node1(), e.node2());
    return flag;
  }

  /** @brief Remove an Edge _e_ from the graph
   *  @param[in] _eit_ An EdgeIterator referencing _e_, i.e. *_eit_ == _e_
   *
   *  @post Edge _e_ is not in the graph
   *  @post If Edge _e_ was in the graph, then new num_edges() = old num_edges() - 1
   *  @post If Edge _e_ was in the graph, all EdgeIterator's are invalidated
   *  @post If Edge _e_ was in the graph, 
   *            For the two connected nodes by _e_, Node _a_ and _b_,
   *            NodeIterator's for _a_ and _b_ are invalidated
   *
   *  @return A new EdgeIterator referring to the next valid Edge in the graph 
   *
   * Complexity: O(k) operations, where k is the maximum degree of Node _a_ or _b_
   *
   */
  edge_iterator remove_edge(edge_iterator eit) {
    auto e = *eit;
    remove_edge(e);
    return ++eit; // Just advance it since we don't delete from internal edge vector
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    i2u_.clear();
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
      assert(index_ < graph_->i2u_.size());
      return graph_->nodes_[graph_->i2u_[index_]].node;
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
    return NodeIterator(this, i2u_.size());
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
      return e; // Always oriented correctly for now
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
      return graph_->nodes_[graph_->i2u_[node_index_]].connections;
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
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator> {
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

    /** @brief dereferencing operator for the edge_iterator
     *  @pre this iterator can be dereferenced
     *  @return the Edge that this edge_iterator references
     */
    Edge operator*() const {
      assert(index_ < graph_->edges_.size());
      return graph_->edges_[index_].edge;
    }

    /** @brief increment operator for an edge_iterator
     *  @pre this iterator is not pointing to graph_->edge_end(),
     *  @return a *this EdgeIterator with its _index__ incremented
     */
    EdgeIterator& operator++() {
      assert(index_ < graph_->edges_.size());
      while (index_ < graph_->edges_.size()) {
        ++index_;
        if (graph_->edges_[index_].valid) { return *this; }
      }
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

    /** Index of current internal_edge */
    size_type index_;

    /** @brief Construct a valid EdgeIterator
     *  @param[in] graph A pointer to the parent graph
     *  @param[in] index The index of the Edge this Iterator references
      */
    EdgeIterator(const Graph* graph, size_type index)
      : graph_(const_cast<Graph*>(graph)), index_(index) { 
        bool valid = graph_->edges_[index_].valid;
        while (not valid and index_ < graph_->edges_.size()) {
          ++index_;
          valid = graph_->edges_[index_].valid;
        }
      }

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
    return edge_iterator(this, edges_.size());
  }

 private:

  /** Internal type for Node */
  struct internal_node {
    Point position;
    size_type index;
    node_value_type value;
    std::vector<edge_type> connections; // List of edges associated with this node
    Node node;
  };

  /** Internal type for Edge */
  struct internal_edge {
    Node node1;
    Node node2;
    Edge edge;
    edge_value_type value;
    bool valid;
    // Constructor for emplace_back
    internal_edge(Node a, Node b, Edge e, bool valid = true) 
    : node1(a), node2(b), edge(e), valid(valid) 
    {}
  };

  /** Vectors containing node structs (in order of assignment/indexing) */
  std::vector<internal_node> nodes_;  // Indexed by node uid
  std::vector<size_type> i2u_;        // Indexed by node index

  /** Vector containing edge structs (in order of assignment/indexing) */
  std::vector<internal_edge> edges_;  // Indexed by edge uid

};

#endif // CME212_GRAPH_HPP