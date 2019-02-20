#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <map>
#include <set>

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
  typedef V node_value_type;
  typedef E edge_value_type;

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

  /* Pair of nodes */
  typedef std::pair<size_type, size_type> NodePair;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() 
    : nodes_(), edges_(), idmap_(), edgemap_(),
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

    /**
     * @brief Get the value of the Node.
     * @return A reference to the Node's value.
     * @pre The Node belongs to a graph.
     */
    node_value_type& value() {
      return this->G_->nodes_[this->index()].value;
    }
    /**
     * @brief Get the value of the Node.
     * @return A constant reference to the Node's value.
     * @pre The Node is a valid node, i.e. it belongs to a graph.
     */
    const node_value_type& value() const {
      return this->G_->nodes_[this->index()].value;
    }
    /**
     * @brief Get the number of nodes adjacent to the current Node.
     * @return The number of adjacent Nodes.
     * @pre The Node is a valid node, i.e. it belongs to a graph.
     */
    size_type degree() const {
      return this->G_->nodes_[this->index()].adj_list.size();
    }
    /**
     * @brief Get the position of the Node.
     */
    Point& position() {
      return this->G_->nodes_[this->index()].point;
    }
    /**
     * @brief Create an iterator pointing to the first Edge incident 
     *        to the Node.
     * @return An IncidentIterator pointing to the first Edge
     *         incident to the Node.
     * @pre The Node is a valid node, i.e. it belongs to a graph.
     */
    incident_iterator edge_begin() const {
      IncidentIterator iit {this->G_, this, 0};
      return iit; 
    }
    /**
     * @brief Create an iterator pointing one element beyond the last 
     *        Edge incident to the Node.
     * @return An IncidentIterator pointing one element beyond the last
     *         Edge incident to the Node.
     * @pre The Node is a valid node, i.e. it belongs to a graph.
     */
    incident_iterator edge_end() const {
      IncidentIterator iit {this->G_, this, this->degree()};
      return iit; 
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
      // Node from another graph
      if (this->G_ < n.G_) {
        return true;
      }
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

  Node add_node(const Point& position, 
                const node_value_type& val = node_value_type()) {

    // Build internal node for storage
    internal_node in;
    in.point = position;
    in.value = val;
    in.adj_list = std::vector<size_type>();
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
  class Edge : private totally_ordered<Edge> {
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
        and std::min(nodepair_.first, nodepair_.second) == 
          std::min(e.nodepair_.first, e.nodepair_.second)
        and std::max(nodepair_.first, nodepair_.second) == 
          std::max(e.nodepair_.first, e.nodepair_.second)) { 
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
      // Edge from another graph
      if (this->G_ < e.G_) {
        return true;
      }
      // Min and max node IDs for each edge
      size_type min1, min2, max1, max2;
      min1 = std::min(this->nodepair_.first, this->nodepair_.second);
      min2 = std::min(e.nodepair_.first, e.nodepair_.second);
      max1 = std::max(this->nodepair_.first, this->nodepair_.second);
      max2 = std::max(e.nodepair_.first, e.nodepair_.second);

      if (min1 < min2) { 
        return true; 
      }
      else if (min1 == min2) {
        if (max1 < max2) { return true; }
      }
      return false;
    }
    
    /** Get the constant edge length */
    double length() const {
      return euclidean_dist(
        this->node1().position(), 
        this->node2().position());
    }

    /** Get the value of the Edge */
    edge_value_type& value() {
      return this->G_->edges_.at(this->index()).value;
    }

    const edge_value_type& value() const {
      return this->G_->edges_.at(this->index()).value;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Private data members
    Graph* G_;
    NodePair nodepair_;

    // Return index of this edge
    size_type index() {
      NodePair nodepair = std::make_pair(
        std::min(this->nodepair_.first, this->nodepair_.second),
        std::max(this->nodepair_.first, this->nodepair_.second));
      return this->G_->edgemap_.at(nodepair);
    }

    // Private constructor for Graph to construct a valid Edge object
    Edge(const Graph* G, NodePair nodepair)
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
    return Edge(this, this->edges_.at(i).nodepair);
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
    if (!has_node(a) or !has_node(b)) {
      return false;
    }

    // Iterate over nodes adjacent to a
    size_type idx_a = this->idmap_.at(a.uid_);
    for (size_type i = 0; i < a.degree(); ++i) {
      if (this->nodes_.at(idx_a).adj_list[i] == b.uid_) {
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
    NodePair nodepair = std::make_pair(
          std::min(a.uid_, b.uid_), 
          std::max(a.uid_, b.uid_));

    // Check if edge is in the graph (if so, return that edge)
    if (has_edge(a, b)) {

      // return Edge();
      return Edge(this, nodepair);
      
    }
    //// Add edge to graph
    internal_edge e;
    e.nodepair = nodepair;
    e.value = edge_value_type();

    // Add to list of edges
    this->edges_.insert(std::make_pair(this->num_edges_, e));
    // Add to the edge map
    this->edgemap_.insert(std::make_pair(nodepair, this->num_edges_));
    // Update number of edges
    this->num_edges_ += 1;
    // Update adjacency lists
    this->nodes_[a.uid_].adj_list.push_back(b.uid_);
    this->nodes_[b.uid_].adj_list.push_back(a.uid_);
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
    this->idmap_ = NIDMap();
    this->edgemap_ = EIDMap();
    this->next_uid_ = 0;
    this->size_ = 0;
    this->num_edges_ = 0;
  }

  /** 
   * @brief Remove a node from the graph
   * @param n Node to be removed
   *
   * @pre Node _n_ is a member of the graph
   * @post num_nodes() and size() have decreased by 1
   * @post has_node(n) returns false
   * @post All edges incident to node _n_ are removed from the graph
   *
   * Complexity: O(k * (num_edges() + num_edges()) where k is the 
   * maximum degree of a node in the graph
   */
  size_type remove_node (const Node& n) {
    if (!has_node(n)) {
      return 0;
    }

    // Remove all incident edges
    for (IncidentIterator iit = n.edge_begin(); iit != n.edge_end();) {
      remove_edge((*iit));
    }
    // Get node index
    size_type index = this->idmap_[n.uid_];

    // Remove from nodes_
    this->nodes_.erase(index);
    // Remove node n from idmap
    this->idmap_.erase(n.uid_);

    // Re-index the last node
    size_type last_idx = this->num_nodes() - 1;
    if (index < last_idx) {
      // Last node in the graph
      internal_node last_node = this->nodes_.at(last_idx);
      // Swap the last node with the current node
      this->nodes_.erase(last_idx);
      this->nodes_.insert(std::make_pair(index, last_node));
      this->idmap_.at(last_node.uid) = index;
    }

    // Decrement graph size
    this->size_ = this->size_ - 1;

    return 1;
  }
  /** 
   * @brief Remove a node from the graph
   * @param n_it node iterator indexed to the node to be removed
   *
   * @pre (*n_it) is a member of the graph
   * @post num_nodes() and size() have decreased by 1
   * @post All edges incident to the node are removed from the graph
   * @post _n_it_ points to the next node
   *
   * Complexity: O(k * (num_edges() + num_edges()) where k is the 
   * maximum degree of a node in the graph
   */
  node_iterator remove_node (node_iterator n_it ) {
    Node n = (*n_it);
    remove_node(n);
    return(n_it);
  }

  /** 
   * @brief Remove an edge from the graph.
   * @param a One node of the edge to be removed. 
   * @param b The other node of the edge to be removed. 
   *
   * @pre has_edge(a, b) returns true (i.e. the graph contains the edge)
   * @post num_edges() has decreased by 1
   * @post has_edge(a, b) returns false
   *
   * Complexity: Average case O(1), worst case O(num_edges).
   */
  size_type remove_edge (const Node& a, const Node& b) {

    if (!has_edge(a, b)) {
      return 0;
    } 
    // Create pair of node UIDs
    NodePair nodepair = std::make_pair(
          std::min(a.uid_, b.uid_), 
          std::max(a.uid_, b.uid_));

    // Get edge index
    size_type index = this->edgemap_.at(nodepair);

    // Remove edges from edges_
    this->edges_.erase(index);

    // Re-index the last edge
    size_type last_idx = this->num_edges_ - 1;
    if (index < last_idx) {
      // Node pair of the last edge in the graph
      internal_edge last_edge = this->edges_.at(last_idx);
      // Re-index
      this->edgemap_.at(last_edge.nodepair) = index;
      this->edges_.erase(last_idx);
      this->edges_.insert(std::make_pair(index, last_edge));
    }
    // Remove edge from edgemap_
    this->edgemap_.erase(nodepair);
    // Update number of edges
    this->num_edges_ -= 1;
    // Remove nodes from adjacency lists
    remove_adj(a, b);

    return 1;
  }

  // Remove nodes a and b from each others' adjacency lists
  void remove_adj (const Node& a, const Node& b) {
    std::vector<size_type>& adj1 = this->nodes_[this->idmap_[a.uid_]].adj_list;
    adj1.erase(std::remove(adj1.begin(), adj1.end(), b.uid_), adj1.end()); 
    std::vector<size_type>& adj2 = this->nodes_[this->idmap_[b.uid_]].adj_list;
    adj2.erase(std::remove(adj2.begin(), adj2.end(), a.uid_), adj2.end()); 
  }

  /** 
   * @brief Remove an edge from the graph.
   * @param e The edge to be removed.
   *
   * @pre has_edge(e.node1(), e.node2()) returns true (i.e. the graph 
   *  contains the edge.
   * @post num_edges() has decreased by 1
   * @post has_edge(e.node1(), e.node2()) returns false
   *
   * Complexity: Average case O(1), worst case O(num_edges)
   */
  size_type remove_edge (const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /** 
   * @brief Remove an edge from the graph.
   * @param e_it Edge iterator indexed to the edge to be removed.
   *
   * @pre _e_it_ is indexed to a valid edge
   * @post num_edges() has decreased by 1
   * @post _e_it_ is indexed to the next edge (unordered)
   *
   * Complexity: Average case O(1), worst case O(num_edges)
   */
  edge_iterator remove_edge (edge_iterator e_it ) {
    Edge e = (*e_it);
    remove_edge(e);
    return(e_it);
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
    NodeIterator() {}

    /**
     * @brief Dereference the NodeIterator.
     * @return The Node at the current position of the NodeIterator.
     *
     * @pre The NodeIterator is indexed to a valid Node; that is, 
     *      node_begin() <= index < node_end()
     */
    Node operator*() const {
      size_type uid = this->G_->nodes_.at(this->index_).uid;
      return Node(this->G_, uid);
    }
    /**
     * @brief Increment the NodeIterator.
     * @return A reference to a NodeIterator pointing to the next Node
     *         in the Graph.
     *
     * @pre The NodeIterator is indexed to a valid Node; that is, 
     *      node_begin() <= index < node_end().
     * @post The NodeIterator is indexed to a valid Node, or to the position
     *        beyond the last Node in the Graph:  
     *        node_begin() <= index <= node_end().
     */
    NodeIterator& operator++(){
      this->index_ += 1;
      return *this;
    }
    /**
     * @brief Check equality against another NodeIterator object.
     *
     * @param n A NodeIterator to compare against.
     * @return A boolean result of equality comparison.
     * @pre NodeIterators _this_ and _n_ must be valid NodeIterator objects.
     *
     * Equal NodeIterator objects have the same graph and the same index.
     */
    bool operator==(const NodeIterator& n) const {
      bool same_graph = ( this->G_ == n.G_ );
      bool same_index = ( this->index_ == n.index_ );
      return (same_index && same_graph);
    }

   private:
    friend class Graph;

    const Graph* G_;
    size_type index_;

    /** Construct a valid NodeIterator */
    NodeIterator(const Graph* G, size_type index) : 
      G_(const_cast<Graph*>(G)), index_(index) {}
  };

  /**
   * @brief Create an iterator pointing to the first Node in a Graph.
   * @return A NodeIterator object pointing to the Node with index 0.
   */
  node_iterator node_begin() const {
    return node_iterator(this, 0);
  }
  /**
   * @brief Create an iterator pointing past the last Node in a Graph.
   * @return A NodeIterator object pointing one position past the last
   *         Node in the graph.
   */
  node_iterator node_end() const {
    return node_iterator(this, this->size());
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
    IncidentIterator() : G_{nullptr}, node_{nullptr} {}

    /**
     * @brief Dereference the IncidentIterator.
     * @return Edge at the current index of the IncidentIterator.
     * @pre The IncidentIterator is indexed to a valid Edge; that is, 
     *      edge_begin() <= index < edge_end()
     */
    Edge operator*() const {
      // Get UID of current adjacent node
      size_type adj_uid = 
        this->G_->nodes_.at(this->node_->index()).adj_list[this->p_];
      // Build the edge
      NodePair nodepair = std::make_pair(this->node_->uid_, adj_uid);
      return Edge(this->G_, nodepair);
    }
    /**
     * @brief Increment the IncidentIterator.
     * @return A reference to a IncidentIterator pointing to the next Edge
     *         incident to the source Node.
     *
     * @pre The IncidentIterator is indexed to a valid Edge; that is, 
     *      edge_begin() <= index < edge_end().
     * @post The IncidentIterator is indexed to a valid Edge, or to the
     *       first position past the last edge incident Edge:  
     *       edge_begin() <= index <= edge_end().
     */
    IncidentIterator& operator++() {
      this->p_ += 1;
      return *this;
    }
    /**
     * @brief Check equality against another IncidentIterator object.
     *
     * @param iit A IncidentIterator to compare against.
     * @return A boolean result of equality comparison.
     * @pre NodeIterators _this_ and _iit_ must be valid IncidentIterator 
     *      objects.
     *
     * Equal IncidentIterator objects have the same graph, the same
     * Node, and the same position within that Node's adjacency list.
     */
    bool operator==(const IncidentIterator& iit) const {
      bool same_graph = ( this->G_ == iit.G_);
      bool same_node = ( this->node_ == iit.node_);
      bool same_index = ( this->p_ == iit.p_);
      return (same_graph && same_node && same_index);
    }

   private:
    friend class Graph;

    /** Private function to create a valid incident iterator */
    IncidentIterator(const Graph* G, const Node* node, size_type p)
        : G_{const_cast<Graph*>(G)}, node_{node}, p_{p} {}

    const Graph* G_; // Pointer to a graph
    const Node* node_; // Pointer to a node
    size_type p_; // Position in the adjacency list
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

    /**
     * @brief Dereference the EdgeIterator.
     * @return The Edge at the current position of the EdgeIterator.
     * @pre The EdgeIterator is indexed to a valid Edge; that is, 
     *      edge_begin() <= index < edge_end()
     */
    Edge operator*() const {
      Edge e {this->G_, this->G_->edges_.at(this->index_).nodepair};
      return e;
    }
    /**
     * @brief Increment the EdgeIterator.
     * @return A reference to a EdgeIterator pointing to the next Edge
     *         in the Graph.
     *
     * @pre The EdgeIterator is indexed to a valid Edge; that is, 
     *      edge_begin() <= index < edge_end().
     * @post The EdgeIterator is indexed to a valid Edge, or to the
     *       first position past the last edge Edge in the Graph:  
     *       edge_begin() <= index <= edge_end().
     */
    EdgeIterator& operator++() {
      this->index_ += 1;
      return *this;
    }
    /**
     * @brief Check equality against another EdgeIterator object.
     *
     * @param e A EdgeIterator to compare against.
     * @return A boolean result of equality comparison.
     * @pre NodeIterators _this_ and _e_ must be valid EdgeIterator 
     *      objects.
     *
     * Equal EdgeIterator objects have the same graph, and the same
     * index to an Edge.
     */
    bool operator==(const EdgeIterator& e) const {
      bool same_graph = ( this->G_ == e.G_ );
      bool same_index = ( this->index_ == e.index_ );
      return (same_index && same_graph);
    }

   private:
    friend class Graph;

    const Graph* G_;
    size_type index_;

    EdgeIterator(const Graph* G, size_type index) : 
      G_(const_cast<Graph*>(G)), index_(index) {}
  };

  /**
   * @brief Create an iterator pointing to the first Edge in the Graph.
   * @return An IncidentIterator pointing to the Edge with index 0.
   */
  edge_iterator edge_begin() const {
    return edge_iterator(this, 0);
  }
  /**
   * @brief Create an iterator pointing past the last Edge in a Graph.
   * @return An EdgeIterator object pointing one position past the last
   *         Edge in the graph.
   */ 
  edge_iterator edge_end() const {
    return edge_iterator(this, this->num_edges_);
  }

 private:

  /* Container for storing node data */
  typedef struct internal_node {
    Point point;   // The position of the node in (x, y, z)
    node_value_type value; // Value of the node
    size_type uid; // The unique identifcation for a node
    std::vector<size_type> adj_list; // List of all adjacent nodes
  } internal_node;

  /* Container for storing edge data */
  typedef struct internal_edge {
    NodePair nodepair; // Two nodes
    edge_value_type value; // Value of the edge 
  } internal_edge;


  /* Map from node index to node */
  typedef std::unordered_map<size_type, internal_node> Nodes;
  /* Map from node UID to node index */
  typedef std::unordered_map<size_type, size_type> NIDMap;
  /* Map from edge index to edge nodes */
  typedef std::unordered_map< size_type, internal_edge> Edges;
  /* Map from node pair to edge index */
  typedef std::map<NodePair, size_type> EIDMap;

  Nodes nodes_;
  Edges edges_;
  NIDMap idmap_;
  EIDMap edgemap_;

  // Next unique UID for a node
  size_type next_uid_;
  // Number of nodes in the graph
  size_type size_;
  // Number of edges in the graph
  size_type num_edges_;

};

#endif // CME212_GRAPH_HPP
