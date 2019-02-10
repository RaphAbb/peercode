#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>

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
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /**  Synonym for V */
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

  /** Struct to hold node information */
  struct NodeInfo
  {
    Point pt;
    node_value_type val;
  };

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
    // HW0: YOUR CODE HERE
    : node_list_(), edge_list_(), adj_list_(), edge_map_() {
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
      // HW0: YOUR CODE HERE
      // Empty constructor creates an invalid Node
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return this->graph_->node_list_[this->uind_].pt;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return this->uind_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Return by reference the node's value.
     * @return by reference the node's value
     */
    node_value_type& value() {
      return this->graph_->node_list_[this->index()].val;
    }

    /** Return by const reference the node's value.
     * @return by const reference the node's value
     */
    const node_value_type& value() const {
      return this->graph_->node_list_[this->index()].val;
    }

    /** Return the degree of the node
     * @return the degree of the node
     */
    size_type degree() const {
      return this->graph_->adj_list_[this->uind_].size();
    }

    /** Return an incident iterator to the beginning of the node's edges
     * @return an incident iterator to the beginning of the node's edges
     */
    incident_iterator edge_begin() const{
      return IncidentIterator(this,0);
    }

    /** Return an incident iterator to the end of the node's edges
     * @return an incident iterator to the end of the node's edges
     */
    incident_iterator edge_end() const{
      return IncidentIterator(this,degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if ((index() == n.index()) and (this->graph_ == n.graph_)){
        return true;
      }
      else{
        return false;
      }
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
      // HW0: YOUR CODE HERE
      if (index() < n.index()){
        return true;
      } else {
        return false;
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_;
    size_type uind_;
    Node(const Graph* graph, size_type uind)
      : graph_(const_cast<Graph*>(graph)), uind_(uind){
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return this->node_list_.size();
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
  Node add_node(const Point& position, const node_value_type& = node_value_type()) {
    // HW0: YOUR CODE HERE
    // Add a vector to the adjacency matrix for each node
    adj_list_.push_back(std::vector<size_type>());
    NodeInfo new_node;
    new_node.pt = position;
    this->node_list_.push_back(new_node);
    return Node(this,this->num_nodes()-1);
  }



  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (this == n.graph_){
      return true;
    } else {
      return false;
    }
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(0 <= i and i < num_nodes());
    // HW0: YOUR CODE HERE
    return Node(this, i);
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
      // HW0: YOUR CODE HERE
      // Empty constructor constructs an invalid edge
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      size_type node_ind = this->graph_->edge_list_[this->uind_][0];
      return this->graph_->node(node_ind);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      size_type node_ind = this->graph_->edge_list_[this->uind_][1];
      return this->graph_->node(node_ind);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {

      size_type node_ind_0 = this->graph_->edge_list_[this->uind_][0];
      size_type node_ind_1 = this->graph_->edge_list_[this->uind_][1];

      size_type e_ind_0 = this->graph_->edge_list_[e.uind_][0];
      size_type e_ind_1 = this->graph_->edge_list_[e.uind_][1];

      if ((node_ind_0 == e_ind_0) and (node_ind_1 == e_ind_1) and (this->graph_ == e.graph_)){
        return true;
      } else if ((node_ind_1 == e_ind_0) and (node_ind_0 == e_ind_1) and (this->graph_ == e.graph_)){
        return true;
      } else {
        return false;
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (this->uind_ < e.uind_){
        return true;
      } else {
        return false;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type uind_;
    Edge(const Graph* graph, size_type uind)
      : graph_(const_cast<Graph*>(graph)), uind_(uind){
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return this->edge_list_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(0 <= i and i < num_edges());
    // HW0: YOUR CODE HERE
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    assert(has_node(a));
    assert(has_node(b));

    size_type a_ind = a.index();
    size_type b_ind = b.index();

    // go to adjacency list.  Iterate through adj_list[a_ind] and compare to b
    for (auto it = this->adj_list_[a_ind].begin(); it != this->adj_list_[a_ind].end(); ++it){
      if (b_ind == *it){
        return true;
      }
    }

    return false;
  }

  /** Find whether two nodes are connected by an edge and return the number.
   * @pre @a a and @a b are valid nodes of this graph and have an edge between them.
   * @return @a i, the index of edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type find_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    assert(has_node(a));
    assert(has_node(b));

    size_type a_ind = a.index();
    size_type b_ind = b.index();

    // go to adjacency list.  Iterate through adj_list[a_ind] and compare to b
    for (size_type i = 0; i < adj_list_[a_ind].size(); i++){
      if (b_ind == adj_list_[a_ind][i]){
        std::pair<size_type, size_type> edge_ends {a_ind, b_ind};
        return edge_map_.at(edge_ends);
      }
    }

    return -1;
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
    // HW0: YOUR CODE HERE
    // Make sure the adjacency list is the same size as the number of nodes
    assert(has_node(a));
    assert(has_node(b));

    size_type a_ind = a.index();
    size_type b_ind = b.index();
    assert(a_ind != b_ind);

    // Check if the edge is already in the edge list
    // This syntax is used instead of has_edge, because we want the index i
    if (has_edge(a,b)){
        return Edge(this,find_edge(a,b));
    }

    // If edge is not in already in the list
    // Add it to the edge list
    std::vector<size_type> new_edge;
    new_edge.push_back(a_ind);
    new_edge.push_back(b_ind);
    edge_list_.push_back(new_edge);
    // Add it in its two spots in the adjacency list and edge map
    adj_list_[a_ind].push_back(b_ind);
    adj_list_[b_ind].push_back(a_ind);
    std::pair<size_type, size_type> edge_ends1 {a_ind,b_ind};
    std::pair<size_type, size_type> edge_ends2 {b_ind,a_ind};
    edge_map_.insert(std::make_pair(edge_ends1,edge_list_.size()-1));
    edge_map_.insert(std::make_pair(edge_ends2,edge_list_.size()-1));
    // Add both orientations into the edge

    return Edge(this, edge_list_.size()-1);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    this->node_list_.clear();
    this->edge_list_.clear();
    this->adj_list_.clear();
    this->edge_map_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator> {
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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Returns the node referenced by the iterator
     * @return the node referenced by the iterator
     */
    Node operator*() const {
      return Node(this->graph_, this->node_ind_);
    }

    /** Increments the iterator and returns the new iterator
     * @return the new iterator
     */
    NodeIterator& operator++(){
      this->node_ind_++;
      return (*this);
    }

    /** Determines quality between NodeIterators
    * @return true or false if the iterators are equal
    */
    bool operator==(const NodeIterator& other) const {
      return (this->node_ind_ == other.node_ind_);
    }

   private:
    friend class Graph;

    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type node_ind_;
    NodeIterator(const Graph* graph, size_type ind)
      : graph_(const_cast<Graph*>(graph)), node_ind_(ind){
    };
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Return a node iterator to the beginning of the graph's nodes
   * @return a node iterator to the beginning of the graph's nodes
   */
  node_iterator node_begin() const{
    return NodeIterator(this,0);
  }

  /** Return a node iterator to the end of the graph's nodes
   * @return a node iterator to the end of the graph's nodes
   */
  node_iterator node_end() const{
    return NodeIterator(this,size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** Returns the edge referenced by the iterator
     * @return the edge referenced by the iterator
     */
    Edge operator*() const{
      size_type neighbor_ind = this->graph_->adj_list_[node_ind_][local_edge_ind_];
      size_type global_edge_ind = this->graph_->find_edge(Node(this->graph_,node_ind_),Node(this->graph_,neighbor_ind));
      return Edge(this->graph_, global_edge_ind);
    }

    /** Increments the iterator and returns the new iterator
     * @return the new iterator
     */
    IncidentIterator& operator++(){
      this->local_edge_ind_++;
      return (*this);
    }

    /** Determines quality between IncidentIterators
    * @return true or false if the iterators are equal
    */
    bool operator==(const IncidentIterator& other) const{
      return ((this->node_ind_ == other.node_ind_) && (this->local_edge_ind_ == other.local_edge_ind_));
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    Node* node_;
    size_type node_ind_;
    size_type local_edge_ind_; // Edge number within this row of the adjacency list
    IncidentIterator(const Node* node, size_type ind)
      : node_(const_cast<Node*>(node)), local_edge_ind_(ind){
        graph_ = node_->graph_;
        node_ind_ = node_->uind_;
    };
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator  : private equality_comparable<EdgeIterator> {
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /** Returns the edge referenced by the iterator
     * @return the edge referenced by the iterator
     */
    Edge operator*() const {
      return Edge(this->graph_, this->edge_ind_);
    }

    /** Increments the iterator and returns the new iterator
     * @return the new iterator
     */
    EdgeIterator& operator++(){
      this->edge_ind_++;
      return (*this);
    }

    /** Determines quality between EdgeIterators
    * @return true or false if the iterators are equal
    */
    bool operator==(const EdgeIterator& other) const {
      return (this->edge_ind_ == other.edge_ind_);
    }

   private:
    friend class Graph;

    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type edge_ind_;
    EdgeIterator(const Graph* graph, size_type ind)
      : graph_(const_cast<Graph*>(graph)), edge_ind_(ind){
    };
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Return a edge iterator to the beginning of the graph's edges
   * @return a edge iterator to the beginning of the graph's edges
   */
  edge_iterator edge_begin() const{
    return EdgeIterator(this,0);
  }

  /** Return a edge iterator to the beginning of the graph's edges
   * @return a edge iterator to the beginning of the graph's edges
   */
  edge_iterator edge_end() const{
    return EdgeIterator(this,this->edge_list_.size());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  std::vector<NodeInfo> node_list_; // List of node information
  std::vector<std::vector<size_type>> edge_list_; // List of edges and either node ids
  std::vector<std::vector<size_type>> adj_list_; // List of neighbors of each node
  std::map<std::pair<size_type, size_type>, size_type> edge_map_; // Map of edge nodes and numbers
};

#endif // CME212_GRAPH_HPP
