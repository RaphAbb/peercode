#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

#define NDEBUG // Ignore assert statements

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <map>
#include <stdexcept>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 *
 * The class is templated such that nodes can hold values of type V.
 */
template<typename V, typename E>
class Graph{
 private:
  // Declare types to be used in class Graph.
  struct InternalNode;
  struct InternalEdge;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V,E>;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
  /** Predeclaration of node value type. */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  /** Predeclaration of edge_value_type. */
  using edge_value_type = E;

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
  Graph(){}

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
    Node() {}

    /** Return this node's position. */
    const Point& position() const {
      assert(is_valid());
      return graph_->nodes_.at(this->uid_).p_;
    }

    Point& position(){
      assert(is_valid());
      return graph_->nodes_.at(this->uid_).p_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(is_valid());
      return graph_->nodes_.at(this->uid_).idx_;
    }

    /** Return this node's value and allows it to be set. */
    node_value_type& value(){
      assert(is_valid());
      return graph_->nodes_.at(this->uid_).val_;
    }

    /* Return this node's value. */
    const node_value_type& value() const{
      assert(is_valid());
      return graph_->nodes_.at(this->uid_).val_;
    }

    /* Return this node's degree, i.e. the number of adjacent nodes. */
    size_type degree() const{
      assert(is_valid());
      return graph_->adj_list_.at(this->uid_).size();
    }

    /* Return IncidentIterator begin instance. Iterates through incient nodes. */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, graph_->adj_list_.at(uid_).begin(), uid_);
    }

    /* Return IncidentIterator end instance. */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_, graph_->adj_list_.at(uid_).end(), uid_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (n.graph_ == this->graph_ and n.index() == this->index());
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
      assert(n.is_valid() && is_valid());
      if(graph_ != n.graph_) return graph_ < n.graph_;
      return index() < n.index();
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    bool is_valid() const {
      return graph_->nodes_.at(uid_).is_valid_;
    }

    Node(graph_type* graph, size_type uid) : graph_(graph), uid_(uid) {}

    graph_type* graph_;
    size_type uid_;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return i2u_n.size();
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
  Node add_node(const Point& p, const node_value_type& val = node_value_type()) {
    size_type idx = next_node_idx();
    size_type uid = next_node_uid();
    // Create and store internal node object
    nodes_.push_back(InternalNode(p, idx, uid, val));
    // Push back uid to i2u
    i2u_n.push_back(uid);
    // Add empty map of node uid -> edge uid to adj list for future use
    adj_list_.push_back(std::map<size_type, size_type>());
    return Node(this, uid); // Create & return node instance
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return this == n.graph_ && n.is_valid();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(0 <= i < num_nodes());
    return Node(const_cast<graph_type*>(this), i2u_n.at(i));
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

    /** Calculate and return the length of this Edge */
    double length() const {
      Point diff = node1().position() - node2().position();
      return norm_2(diff);
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, node2_);
    }

    /** Return the index of this Edge */
    size_type index() const {
      return graph_->edges_.at(uid_).idx_;
    }

    /** Return the value (modifiable) of this Edge */
    edge_value_type& value(){
      return graph_->edges_.at(uid_).val_;
    }

    /** Return the value of this Edge */
    const edge_value_type& value() const {
      return graph_->edges_.at(uid_).val_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes and
     * belong to the same graph.
     */
    bool operator==(const Edge& e) const {
      if(graph_ != e.graph_) return false;
      if(e.node1_ == node1_ and e.node2_ == node2_){
        return true;
      }
      else if(e.node1_ == node2_ and e.node2_ == node1_){
        return true;
      }
      else{
        return false;
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     *
     * If this edge and @a e come from different graphs the methods resorts to
     * operator< applied to graph pointers.
     */
    bool operator<(const Edge& e) const {
      if(graph_ != e.graph_) return graph_ < e.graph_;
      return index() < e.index();
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    bool is_valid() const {
      return graph_->edges_.at(uid_).is_valid_;
    }

    /** Construct valid edge object */
    Edge(size_type node1, size_type node2, size_type uid, graph_type* graph)
      : node1_(node1), node2_(node2), uid_(uid), graph_(graph) {}

    size_type node1_; // Node 1 uid
    size_type node2_; // Node 2 uid
    size_type uid_;   // This edge's uid
    graph_type* graph_;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1).
   */
  size_type num_edges() const {
    return i2u_e.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(1).
   */
  Edge edge(size_type i) const {
    assert(0 <= i < num_edges());
    size_type uid = i2u_e.at(i);
    return Edge(edges_.at(uid).node1_, edges_.at(uid).node2_, uid, const_cast<graph_type*>(this));
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert(has_node(a) and has_node(b));
    // Return whether edge exists.
    return (adj_list_.at(a.uid_).find(b.uid_) != adj_list_.at(a.uid_).end());
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val = edge_value_type()) {
    // Can't add edge to graph if nodes don't belong to the graph
    assert(!(a.graph_ != this || b.graph_ != this));
    assert(a != b); // Self loops are not allowed (actually, frowned upon)

    // If the edge exists loop over edges for node a, find index, return edge
    if (has_edge(a, b)){
      size_type uid;
      uid = adj_list_.at(a.uid_).at(b.uid_);
      return edge(edges_.at(uid).idx_);
    }
    // Else create new edge and store information in appropriate containers
    else{
      size_type idx = next_edge_idx();
      size_type uid = next_edge_uid();
      adj_list_.at(a.uid_)[b.uid_] = uid;
      adj_list_.at(b.uid_)[a.uid_] = uid;
      i2u_e.push_back(uid);
      edges_.push_back(InternalEdge(a.uid_, b.uid_, idx, uid, val));
      return Edge(a.uid_, b.uid_, uid, this);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 and num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // Clear graph containers
    nodes_.clear();
    adj_list_.clear();
    edges_.clear();
    i2u_e.clear();
    i2u_n.clear();
    // Check post condition
    assert(num_nodes() == 0 && num_edges() == 0);
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {}

    /** Return current node instance. */
    Node operator*() const{
      size_type uid = *it_;
      return Node(graph_, uid);
    }

    /** Get reference to next node instance. */
    NodeIterator& operator++(){
      ++it_;
      return *this;
    }

    /** Test whether this node iterator and @ other are equal. Wrapper for
     *  STL vector iterator operator==.
     */
    bool operator==(const NodeIterator& other) const{
      return it_ == other.it_;
    };

   private:
    friend class Graph; // Allow Graph to call private constructor.

    graph_type* graph_;
    typename std::vector<size_type>::const_iterator it_;

    NodeIterator(graph_type* graph, typename std::vector<size_type>::const_iterator it)
    : graph_(graph), it_(it) {}
  };

  /** Return begin instance of NodeIterator */
  node_iterator node_begin() const{
    return NodeIterator(const_cast<graph_type*>(this), i2u_n.begin());
  }
  /** Return end instance of NodeIterator */
  node_iterator node_end() const{
    return NodeIterator(const_cast<graph_type*>(this), i2u_n.end());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {}

    /* Return Edge instance of current incident edge to node. */
    Edge operator*() const{
      return Edge(node_, it_->first, it_->second, graph_);
    }

    /* Increment underlying iterator and return reference to current instance. */
    IncidentIterator& operator++(){
      it_++;
      return *this;
    }

    /* Wrapper to underlying operator== of STL map iterator. */
    bool operator==(const IncidentIterator& other) const{
      return it_ == other.it_;
    }

   private:
    friend class Graph; // Allow Graph to access private constructor

    graph_type* graph_;
    typename std::map<size_type, size_type>::const_iterator it_; // Underlying map iterator
    size_type node_;                                             // uid of base node

    // Private constructor called by Graph methods
    IncidentIterator(graph_type* graph, typename std::map<size_type, size_type>::const_iterator it, size_type node)
    : graph_(graph), it_(it), node_(node) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. Wrapper for STL
   *        vector container.
   */
  class EdgeIterator : private totally_ordered<EdgeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {}

    /* Dereference operator returns current Edge as Edge instance. */
    Edge operator*() const {
      size_type uid = *it_;
      InternalEdge ie = graph_->edges_.at(uid);
      return Edge(ie.node1_, ie.node2_, uid, graph_);
    }

    /* Increment underlying iterator and return this instance. */
    EdgeIterator& operator++(){
      it_++;
      return *this;
    }

    /* Test whether iterators are equal, wrapper for underlying iterator
     * operator==.
     */
    bool operator==(const EdgeIterator& other) const{
      return it_ == other.it_;
    }

   private:
    friend class Graph;

    graph_type* graph_;
    typename std::vector<size_type>::const_iterator it_; // Underlying iterator

    // Private constructor
    EdgeIterator(graph_type* graph, typename std::vector<size_type>::const_iterator it)
    : graph_(graph), it_(it) {}
  };

  /* Returns a begin instance of forward iterator EdgeIterator over edges. */
  edge_iterator edge_begin() const{
    return EdgeIterator(const_cast<graph_type*>(this), i2u_e.begin());
  }

  /* Returns a end instance of forward iterator EdgeIterator over edges. */
  edge_iterator edge_end() const{
    return EdgeIterator(const_cast<graph_type*>(this), i2u_e.end());
  }

  /** Remove an edge from the graph if it exists and is valid.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return a size_type denoting if the removal was successful in a boolean
   *         sense.
   * @post has_edge(@a a, @a b) == false
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges() - 1.
   *       Else,                        new num_edges() == old num_edges().
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: O(log(num_edges()))
   */
  size_type remove_edge(const Node& a, const Node& b){
    if(!has_edge(a, b)){
      return 0;
    }
    // Get unique id of deleted edge
    size_type del_edge_uid = adj_list_.at(a.uid_).at(b.uid_);

    // Update adjacency list
    adj_list_.at(a.uid_).erase(b.uid_);
    adj_list_.at(b.uid_).erase(a.uid_);

    // Get index of deleted node
    size_type del_edge_idx = edges_.at(del_edge_uid).idx_;

    edges_.at(del_edge_uid).is_valid_ = false; // Make deleted internal edge invalid

    // Swap and pop internal i2u entries
    size_type moved_uid = i2u_e.back();

    i2u_e.at(del_edge_idx) = i2u_e.back();    // Move back edge to replace deleted
    i2u_e.pop_back();                         // Remove the last edge
    edges_.at(moved_uid).idx_ = del_edge_idx; // Update index of moved edge

    assert(has_edge(a,b) == false);
    return 1;
  }

  /** Synonym for remove_edge(const Node& a, const Node& b) that takes an edge.
   * @pre @a e is a valid edge of this graph.
   */
  size_type remove_edge(const Edge& e){
    return remove_edge(e.node1(), e.node2());
  }

  /** Remove an edge from the graph pointed to by @a e_it if it exists and
   *  is valid.
   * @pre *(@a e_it) is a valid edge of this graph.
   * @return If @a e_it == edge_end(), edge_end()
             If !has_edge("*nodes"),   edge_begin()
             Else,                     Edge iterator pointing to the new edge
                                       with the same index as the removed edge.
   * @post has_edge("*nodes") == false
   * @post If old has_edge(),     new num_edges() == old num_edges() - 1.
   *       Else,                  new num_edges() == old num_edges().
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: Same complexity as remove_edge(Node& a, Node& b)
   */
  edge_iterator remove_edge(edge_iterator e_it){
    // If iterator is at end do nothing and return the iterator
    if(e_it == edge_end()){
      return e_it;
    }

    Edge edge = *e_it;
    size_type res = remove_edge(edge);
    // If remove edge removal fails return begin iterator
    if(!res){
      return ++e_it;
    }
    // Else return valid iterator to the new edge with old edge's index
    else{
      return e_it;
    }
  }

  /** Remove a node from the graph if it exists and is valid. Also removes
   *  all edges incident to the node.
   * @pre @a n is a valid node of this graph.
   * @return a size_type denoting if the removal was successful in a boolean
   *         sense.
   * @post has_node(@a n) == false.
   * @post If old has_edge(@a n), new num_nodes() == old num_nodes() - 1.
   *       Else,                  new num_nodes() == old num_nodes().
   *
   * Can invalidate node indexes -- in other words, old node(@a i) might not
   * equal new node(@a i). Must not invalidate outstanding Node objects.
   *
   * Complexity: O((@a n).degree())
   */
  size_type remove_node(const Node& n){
    if(!has_node(n)){
      return 0;
    }

    // Delete all edges incident to node
    for(incident_iterator it = n.edge_begin(); it != n.edge_end();){
      it = n.edge_begin();
      Edge edge = *it;
      remove_edge(edge);
    }
    // Swap and pop internal i2u entries.
    size_type del_node_idx = n.index();
    size_type del_node_uid = n.uid_;
    nodes_.at(del_node_uid).is_valid_ = false;
    size_type moved_uid = i2u_n.back();
    i2u_n.at(del_node_idx) = i2u_n.back();
    i2u_n.pop_back();
    nodes_.at(moved_uid).idx_ = del_node_idx;
    return 1;
  }

  /** Remove a node from the graph pointed to by @a n_it if it exists and
   *  is valid. Also removes all edges incident to the node.
   * @pre *(@a n_it) is a valid node of this graph.
   * @return If @a n_it == node_end(), node_end()
             If !has_node(*(@a n_it)), node_begin()
             Else,                     Node iterator pointing to the new node
                                       with the same index as the removed node.
   * @post has_node(*(@a n_it)) == false
   * @post If old has_edge(@a n), new num_nodes() == old num_nodes() - 1.
   *       Else,                  new num_nodes() == old num_nodes().
   *
   * Can invalidate node indexes -- in other words, old node(@a i) might not
   * equal new node(@a i). Must not invalidate outstanding Node objects.
   *
   * Complexity: Same complexity as remove_node(Node&)
   */
  node_iterator remove_node(node_iterator n_it){
    if(n_it == node_end()){
      return node_end();
    }
    Node n = *n_it;
    size_type res = remove_node(n);
    if(!res){
      return ++n_it;
    }
    else{
      return n_it;
    }
  }

 private:
   size_type next_node_idx(){
     return i2u_n.size();
   }
   size_type next_node_uid(){
     return nodes_.size();
   }
   size_type next_edge_uid(){
     return edges_.size();
   }
   size_type next_edge_idx(){
     return i2u_e.size();
   }

  // Structure to store an internal copy of node
  struct InternalNode{
    InternalNode(Point p, size_type idx, size_type uid, node_value_type val)
      : p_(p), idx_(idx), uid_(uid), val_(val) {}
    Point p_;
    size_type idx_;
    size_type uid_;
    node_value_type val_;
    bool is_valid_ = true;
  };

  // Store node indices and index of an edge, used to store internal Edge copies.
  struct InternalEdge{
    InternalEdge(size_type a, size_type b, size_type idx, size_type uid, edge_value_type val)
      : node1_(a), node2_(b), idx_(idx), uid_(uid), val_(val) {}
    size_type node1_;
    size_type node2_;
    size_type idx_;
    size_type uid_;
    edge_value_type val_;
    bool is_valid_ = true;
  };

  std::vector<InternalNode> nodes_;                      // node uid -> InternalNode
  std::vector<InternalEdge> edges_;                      // edge uid -> InternalEdge
  std::vector<std::map<size_type, size_type>> adj_list_; // For each node uid, map 2nd node uid to edge uid
  std::vector<size_type> i2u_n;                          // idx -> uid for nodes
  std::vector<size_type> i2u_e;                          // idx -> uid for edges
};

#endif // CME212_GRAPH_HPP
