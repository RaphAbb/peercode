#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <stdexcept>
#include <vector>
#include <map>

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
template<typename V>
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
  using graph_type = Graph<V>;

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
    Node() {

    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_.at(this->idx_).p_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return this->idx_;
    }

    /** Return this node's value and allows it to be set. */
    node_value_type& value(){
      return graph_->nodes_.at(this->idx_).val_;
    }

    /* Return this node's value. */
    const node_value_type& value() const{
      return graph_->nodes_.at(this->idx_).val_;
    }

    /* Return this node's degree, i.e. the number of adjacent nodes. */
    size_type degree() const{
      return graph_->adj_list_.at(idx_).size();
    }

    /* Return IncidentIterator begin instance. Iterates through incient nodes. */
    incident_iterator edge_begin() const{
      return IncidentIterator(const_cast<graph_type*>(graph_), graph_->adj_list_.at(idx_).begin(), idx_);
    }

    /* Return IncidentIterator end instance. */
    incident_iterator edge_end() const{
      return IncidentIterator(const_cast<graph_type*>(graph_), graph_->adj_list_.at(idx_).end(), idx_);
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
      return index() < n.index();
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    graph_type* graph_;
    size_type idx_;

    Node(graph_type* graph, size_type idx)
      : graph_(graph), idx_(idx) {

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
  Node add_node(const Point& p, const node_value_type& val = node_value_type()) {
    // Create and store internal node object
    nodes_.push_back(InternalNode(p, next_idx_, val));
    next_idx_++;                     // Increment index
    // Add empty map of node->edge to adj list for future use
    adj_list_.push_back(std::map<size_type, size_type>());
    return Node(this, next_idx_- 1); // Create & return node instance
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
    return Node(const_cast<graph_type*>(this), i);
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

    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, node2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
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
     */
    bool operator<(const Edge& e) const {
      return idx_ < e.idx_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    size_type node1_;
    size_type node2_;
    size_type idx_;
    graph_type* graph_;

    /** Construct valid edge object */
    Edge(size_type node1, size_type node2, size_type idx, graph_type* graph)
      : node1_(node1), node2_(node2), idx_(idx), graph_(graph) {

    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1).
   */
  size_type num_edges() const {
    return next_edge_idx_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(1).
   */
  Edge edge(size_type i) const {
    if(i > next_edge_idx_-1){
      throw std::out_of_range("Edge index out of range when calling edge.");
    }
    return Edge(edges_.at(i).node1_, edges_.at(i).node2_, i, const_cast<graph_type*>(this));
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Check that nodes belong to this graph
    if((has_node(a) and has_node(b)) == false){
      throw std::runtime_error("Node(s) do not belong to graph in has_edge()");
    }
    // Check if edge exists in both orders
    if(adj_list_.at(a.index()).count(b.index()) > 0){
      return true;
    }
    else if(adj_list_.at(b.index()).count(a.index()) > 0){
      return true;
    }
    else{
      return false;
    }
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
    // Check that a and b are distinct valid nodes of this graph
    if(a.graph_ != this || b.graph_ != this){
      throw std::runtime_error("One of the nodes do not belong to the graph.");
    }
    if(a == b){
      throw std::runtime_error("add_edge() received non-distinct nodes. Self-loops are frowned upon.");
    }

    // If the edge exists loop over edges for node a, find index, return edge
    if (has_edge(a, b)){
      size_type index;
      if(adj_list_.at(a.index()).count(b.index()) > 0){
        index = adj_list_.at(a.index()).at(b.index());
      }
      else{
        index = adj_list_.at(b.index()).at(a.index());
      }
      return edge(index);
    }
    // Else create new edge and store information in appropriate containers
    else{
      next_edge_idx_++;
      InternalEdge edge_ab(a.index(), b.index(), next_edge_idx_-1);
      adj_list_[a.index()][b.index()] = next_edge_idx_-1;
      adj_list_[b.index()][a.index()] = next_edge_idx_-1;
      edges_.push_back(InternalEdge(a.index(), b.index(), next_edge_idx_-1));
      return Edge(a.index(), b.index(), next_edge_idx_-1, this);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 and num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // Reset counters
    next_idx_      = 0;
    next_edge_idx_ = 0;
    // Clear graph containers
    nodes_.clear();
    adj_list_.clear();
    edges_.clear();
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
    NodeIterator() {
    }

    /** Return current node instance. */
    Node operator*() const{
      InternalNode int_node = *it_;
      return Node(graph_, int_node.idx_);
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
    typename std::vector<InternalNode>::const_iterator it_;
    NodeIterator(graph_type* graph, typename std::vector<InternalNode>::const_iterator it)
    : graph_(graph), it_(it){

    }
  };

  /** Return begin instance of NodeIterator */
  node_iterator node_begin() const{
    return NodeIterator(const_cast<graph_type*>(this), nodes_.begin());
  }
  /** Return end instance of NodeIterator */
  node_iterator node_end() const{
    return NodeIterator(const_cast<graph_type*>(this), nodes_.end());
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
    IncidentIterator() {
    }

    /* Return Edge instance of current incident edge to node. */
    Edge operator*() const{
      return Edge(node_, it_->first, it_->second, graph_);;
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
    size_type node_;
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
    EdgeIterator() {
    }

    /* Dereference operator returns current Edge as Edge instance. */
    Edge operator*() const {
      InternalEdge int_edge = *it_;
      return Edge(int_edge.node1_, int_edge.node2_, int_edge.idx_, graph_);
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
    typename std::vector<InternalEdge>::const_iterator it_; // Underlying iterator

    // Private constructor
    EdgeIterator(graph_type* graph, typename std::vector<InternalEdge>::const_iterator it)
    : graph_(graph), it_(it) {}
  };

  /* Returns a begin instance of forward iterator EdgeIterator over edges. */
  edge_iterator edge_begin() const{
    return EdgeIterator(const_cast<graph_type*>(this), edges_.begin());
  }

  /* Returns a end instance of forward iterator EdgeIterator over edges. */
  edge_iterator edge_end() const{
    return EdgeIterator(const_cast<graph_type*>(this), edges_.end());
  }

 private:
  // Structure to store an internal copy of node
  struct InternalNode{
    InternalNode(Point p, size_type idx, node_value_type val)
      : p_(p), idx_(idx), val_(val) {

    }
    Point p_;
    size_type idx_;
    node_value_type val_;
  };

  // Store node indices and index of an edge, used to store internal Edge copies.
  struct InternalEdge{
    InternalEdge(size_type a, size_type b, size_type idx)
      : node1_(a), node2_(b), idx_(idx){

    }
    size_type node1_;
    size_type node2_;
    size_type idx_;
  };

  std::vector<InternalNode> nodes_;                         // node idx -> Node
  std::vector<InternalEdge> edges_;                         // edge idx -> Edge
  std::vector<std::map<size_type, size_type>> adj_list_;    // For each node, map 2nd node idx to edge idx
  size_type next_idx_ = 0;
  size_type next_edge_idx_ = 0;
};

#endif // CME212_GRAPH_HPP
