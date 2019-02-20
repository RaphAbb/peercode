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
template  <typename V, typename E>
class Graph : private totally_ordered<Graph<V, E>> {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  struct internal_node;
  struct internal_edge;

  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V, E>;

  /** Type of elements in graph **/
  using node_value_type = V;
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

  /** Construct an empty graph. */
  Graph() {
    next_node_uid = 0;
    next_edge_uid = 0;
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
    }

    /** Return this node's position. */
    const Point& position() const {
      return fetch().pos;
    }

    /** Return this node's position. */
    Point& position() {
      return fetch().pos;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Get the value of the node.
     * @return a (non-const) reference the Node's value.
     */
    node_value_type& value(){
      return fetch().value;
    }

    /** Get the Node's value.
     * @return a *const* reference to the Node's value.
     */
    const node_value_type & value () const{
      return const_cast<const node_value_type&>(fetch().value);
    }

    /** Get the number of edges adjacent to that Node.
     * @return the Node's degree.
     */
    size_type degree() const {
      return graph_->adjacency_[uid_].size();
    }

    /** Get a begin IncidentIterator for the edges pointing to that Node.
     * @return it, an IndicentIterator.
     * @post ((*it).node1() == *this) = true
     */
    incident_iterator edge_begin() const{
      return IncidentIterator(const_cast<Graph*>(graph_), uid_, 0);
    }

    /** Get an end IncidentIterator for the edges pointing to that node.
     * @return it, an IndicentIterator (will never be dereferenced as it
     *  points to an invalid final value).
     */
    incident_iterator edge_end() const{
      return IncidentIterator(const_cast<Graph*>(graph_), uid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
    	return (graph_ == n.graph_) && (index() == n.index());
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
      return graph_ < n.graph_ || n.index() < index();
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_;
    size_type uid_;

    internal_node& fetch() const {
      return graph_->nodes_[uid_];
    };

    Node(Graph* graph, size_type uid) : graph_(graph), uid_(uid) { };
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
   * @param[in] value (optional) The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() ==f old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()){
	  internal_node x = {
			const_cast<Point&>(position),
			next_node_uid,
      value
	  };
	  nodes_.push_back(x);
	  next_node_uid++;
    // Add new node to adjacency list.
    adjacency_.push_back(std::vector<size_type>());
	  return Node(this, next_node_uid-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return size() > n.index();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
	 assert(i < size());
	 return Node(const_cast<Graph*>(this), nodes_[i].uid);
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
  class Edge : private equality_comparable<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() { }

    /** Return a node of this Edge */
    Node node1() const {
      if(reverse){
          return graph_->node(fetch().node2_index);
      } else {
          return graph_->node(fetch().node1_index);
      }
    }

    /** Return the other node of this Edge */
    Node node2() const {
      if(reverse){
          return graph_->node(fetch().node1_index);
      } else {
          return graph_->node(fetch().node2_index);
      }
    }

    /** Return the length of this Edge in L2 norm */
    double length() const {
      return norm_2(node1().position() - node2().position());
    }

    /** Returns the UID of this Edge */
    size_type index() const {
      return uid_;
    }

    edge_value_type& value(){
      return fetch().value;
    }
    const edge_value_type& value() const{
      return fetch().value;
    };

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      bool condition1 = node1() == e.node1() && node2() == e.node2();
      bool condition2 = node2() == e.node1() && node1() == e.node2();
      return (condition1 || condition2);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return node1() < e.node1() || node2() < e.node2();
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type uid_;

    // Allow IncidentIterator to set the reverse flag.
    friend class IncidentIterator;
    bool reverse = false;

    internal_edge& fetch() const {
          return graph_->edges_.at(uid_);
        };

    Edge(const Graph* graph, size_type uid) : graph_(const_cast<Graph*>(graph)), uid_(uid) { }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
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
    assert(i < num_edges());
    return Edge(this, edges_[i].edge_uid);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   * HW2: Complexity reduced to O(max degree).
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert(has_node(a));
    assert(has_node(b));
    for(auto it = a.edge_begin(); it != a.edge_end(); ++it){
      if((*it).node2() == b) return true;
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
    if(has_edge(a,b)){
      return get_edge(a,b);
    }
    internal_edge edge;
    edge.edge_uid = next_edge_uid;
    edge.node1_index = a.index();
    edge.node2_index = b.index();
    edges_.push_back(edge);
    // Add edge to adjacency list.
    adjacency_[a.index()].push_back(edge.edge_uid);
    adjacency_[b.index()].push_back(edge.edge_uid);
    next_edge_uid++;
    return Edge(this, next_edge_uid-1);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
      next_edge_uid = 0;
      next_node_uid = 0;
      nodes_.clear();
      edges_.clear();
      adjacency_.clear();
  }


  /** Removes a node from the graph by reference.
   * @param[in] n - const reference to node to remove.
   * @return 0 if no node was removed (do not take post conditions into account), 1 if a node was removed.
   * @post new num_nodes() = old num_nodes() - 1
   * @post invalidates all existing nodes, NodeIterators and EdgeIterators (can be refined to nodes with index greater than n).
   * @post removes all edges incident to @a n
   * @post g.node(n.index()) == n for any node n constructed after function is used.
   * @post g.node(i).index() == i for all 0 <= i < new num_nodes().
   * @post new num_edges() == old num_edges() - n.degree()
   *
   * Complexity: O(max_degree*num_edges() + num_edges() + num_nodes())
   */
  size_type remove_node(const Node& n){
    if(!has_node(n)) return 0;
    // Remove all associated edges. Correct complexity assuming the graph is sparse.
    auto it = n.edge_begin();
    while(it != n.edge_end()){
      it = remove_edge(it);
    };
    size_type node_index = n.index();
    adjacency_.erase(adjacency_.begin() + node_index);
    nodes_.erase(nodes_.begin() + node_index);
    // Decrement all indices after removed node in nodes_.
    for(auto it = nodes_.begin() + node_index; it != nodes_.end(); ++it) (*it).uid--;
    // Decrement all indices after removed node in edges_.
    for(auto it = edges_.begin(); it != edges_.end(); ++it){
      if((*it).node1_index >= node_index) (*it).node1_index--;
      if((*it).node2_index >= node_index) (*it).node2_index--;
    }
    next_node_uid--;
    return 1;
  }

  /** Removes a node from the graph by iterator and returns a new iterator.
   * @param[in] n_it - iterator pointing to node to remove.
   * @pre n_it points to a valid node (ie it has not been invalidated by a removal operation).
   * @return NodeIterator to next node (no other node is skipped).
   * @post new num_nodes() = old num_nodes() - 1
   * @post invalidates all existing nodes, NodeIterators and EdgeIterators (can be refined to nodes with index greater than *n_it).
   * @post removes all edges incident to @a *n_it
   * @post g.node(n.index()) == n for any node n constructed after function is used.
   * @post g.node(i).index() == i for all 0 <= i < new num_nodes().
   * @post new num_edges() == old num_edges() - @a (*n_it).degree()
   *
   * Complexity: O(max_degree*num_edges() + num_edges() + num_nodes())
   */
  node_iterator remove_node(node_iterator n_it){
    size_type node_index = (*n_it).index();
    remove_node(*n_it);
    return NodeIterator(this, node_index);
  }

  /** Removes a node from the graph by iterator and returns a new iterator.
   * @param[in] a,b - Two nodes.
   * @pre @a a and @a b are valid nodes, and there is an edge between a and b (can be checked using has_edge).
   * @return i - index of the corresponding edge.
   * @post (edge(i).node1() == a && edge(i).node2() == b) || (edge(i).node2() == a && edge(i).node1() == b)
   *
   * Complexity: O(max_degree)
   */
  size_type get_edge_index(const Node& a, const Node& b){
    for(auto it = a.edge_begin(); it != a.edge_end(); ++it){
      if((*it).node2() == b) return (*it).index();
    }
    assert(false);
  }

  /** Removes an edge from the graph by end nodes.
   * @param[in] a,b - potential end nodes.
   * @return 0 if no edge was removed (ignore post conditions), 1 if an edge was removed.
   * @post new num_edges() = old num_edges() - 1
   * @post invalidates all existing edges, NodeIterators, EdgeIterators and IncidentIterators (can be refined to edges with index greater than e).
   * @post g.edge(e.index()) == e for any Edge e constructed after function is used.
   * @post there is no edge e such that (e.node1() == a && e.node2() == b) || (e.node1() == b && e.node2() == a)
   * @post new @a a.degree() = old @a a.degree() - 1
   * @post new @a b.degree() = old @a b.degree() - 1
   * @post g.edge(i).index() == i for all 0 <= i < new num_edges().
   *
   * Complexity: O(num_edges())
   */
  size_type remove_edge(const Node& a, const Node& b){
    if(!has_edge(a, b)) return 0;
    // If edge exists.
    size_type edge_index = get_edge_index(a, b);
    remove_edge(edge(edge_index));
    return 1;
  }

  /** Removes an edge from the graph by reference.
   * @param[in] e - reference to the Edge to remove.
   * @return 0 if no edge was removed (ignore post conditions), 1 if an edge was removed.
   * @post new num_edges() = old num_edges() - 1
   * @post invalidates all existing edges, NodeIterators, EdgeIterators and IncidentIterators (can be refined to edges with index greater than e).
   * @post there is no edge e such that (e.node1() == a && e.node2() == b) || (e.node1() == b && e.node2() == a)
   * @post g.edge(e.index()) == e for any Edge e constructed after function is used.
   * @post g.edge(i).index() == i for all 0 <= i < new num_edges().
   *
   * Complexity: O(num_edges())
   */
  size_type remove_edge(const Edge& e){
    size_type edge_index = e.index();
    // Remove edge from adjacency lists (both ends).
    adjacency_[e.node1().index()].erase(
      std::remove(adjacency_[e.node1().index()].begin(), adjacency_[e.node1().index()].end(), e.index()),
      adjacency_[e.node1().index()].end());
    adjacency_[e.node2().index()].erase(
      std::remove(adjacency_[e.node2().index()].begin(), adjacency_[e.node2().index()].end(), e.index()),
      adjacency_[e.node2().index()].end());
    // Remove edge from edge list.
    edges_.erase(edges_.begin() + edge_index);
    // Decrement edge uids in edges_
    for(auto it = edges_.begin() + edge_index; it != edges_.end(); ++it) (*it).edge_uid--;
    // Decrement edge uids in adjacency_
    for(auto outer_it = adjacency_.begin(); outer_it != adjacency_.end(); ++outer_it){
      for(auto inner_it = (*outer_it).begin(); inner_it != (*outer_it).end(); ++inner_it){
        if((*inner_it) >= edge_index) (*inner_it)--;
      }
    }
    next_edge_uid--;
    return 0;
  }

  /** Removes an edge from the graph by EdgeIterator.
   * @param[in] e_it - EdgeIterator pointing to edge to remove.
   * @return it, new EdgeIterator (only skips *e_it).
   * @post new num_edges() = old num_edges() - 1
   * @post invalidates all existing edges, NodeIterators, EdgeIterators and IncidentIterators (can be refined to edges with index greater than e).
   * @post there is no edge e such that (e.node1() == a && e.node2() == b) || (e.node1() == b && e.node2() == a)
   * @post g.edge(e.index()) == e for any Edge e constructed after function is used.
   * @post g.edge(i).index() == i for all 0 <= i < new num_edges().
   *
   * Complexity: O(num_edges())
   */
  edge_iterator remove_edge(edge_iterator e_it){
    size_type edge_index = (*e_it).index();
    remove_edge(*e_it);
    return EdgeIterator(this, edge_index);
  }

  /** Removes an edge from the graph by IncidentIterator (not asked in HW2, added for practical reasons).
   * @param[in] e_it - IncidentIterator pointing to edge to remove.
   * @return it, new IncidentIterator (only skips *e_it).
   * @post new num_edges() = old num_edges() - 1
   * @post invalidates all existing edges, NodeIterators, EdgeIterators and IncidentIterators (can be refined to edges with index greater than e).
   * @post there is no edge e such that (e.node1() == a && e.node2() == b) || (e.node1() == b && e.node2() == a)
   * @post g.edge(e.index()) == e for any Edge e constructed after function is used.
   * @post g.edge(i).index() == i for all 0 <= i < new num_edges().
   *
   * Complexity: O(num_edges())
   */
  incident_iterator remove_edge(incident_iterator e_it){
    size_type edge_index = e_it.index;
    remove_edge(*e_it);
    return IncidentIterator(this, e_it.node_uid, edge_index);
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator>{
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

    /** Returns the node the iterator is pointing to.
     */
    Node operator *() const{
      return graph->node(index);
    };

    /** Increments the iterator by one, to the next node in the list.
     */
    node_iterator& operator++(){
      ++index;
      return *this;
    };

    /** Compares the NodeIterator to another NodeIterator.
     * Both the index and the related graph are compared.
     */
    bool operator==(const node_iterator& other) const {
      return (index == other.index && graph == (other.graph));
    };

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    size_type index;
    const Graph* graph;

    NodeIterator(const Graph* g, size_type i){
        graph = g;
        index = i;
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Returns a node iterator that points to the first node in the graph.
   */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Returns a node iterator that points to the end of the graph's nodes.
   */
  node_iterator node_end() const {
    return NodeIterator(this, num_nodes());
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

    /** Returns the Edge the iterator is pointing to.
     * @return a valid Edge object edge
     * @post edge.node1() = graph->node(node_uid);
     * @post graph->has_edge(edge.node1(), edge.node2()) = true;
     */
    Edge operator*() const{
      Edge edge = graph->edge((graph->adjacency_)[node_uid][index]);
      if(edge.node1().uid_ != node_uid){
        edge.reverse = true;
      }
      return edge;
    }

    /** Increments the iterator to the next edge incident to node(node_uid).
     */
    IncidentIterator& operator++(){
      index++;
      return *this;
    }

    /** Compares the IncidentIterator to another IncidentIterator.
     * The index, the node and the graph are compared.
     */
    bool operator==(const IncidentIterator& other) const{
      return (index == other.index) && (graph == other.graph) && (node_uid == other.node_uid);
    }

   private:
    friend class Graph;
    const Graph* graph;
    size_type node_uid;
    size_type index;

    /** Private constructor, used by edge_begin() to construct an IncidentIterator to a given node.
     */
    IncidentIterator(const Graph* g, const size_type nuid, size_type i){
      graph = g;
      node_uid = nuid;
      index = i;
    }
    // HW1 #3: YOUR CODE HERE
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /** Returns the Edge the iterator is pointing to.
     */
    Edge operator*() const {
      return graph->edge(index);
    }

    /** Increments the iterator to the next edge in the graph.
     */
    EdgeIterator& operator++(){
      index++;
      return *this;
    }

    /** Compares the EdgeIterator to another EdgeIterator.
     * Both the index and the graph are compared.
     */
    bool operator==(const EdgeIterator& other) const {
      return (index == other.index) && (graph == other.graph);
    }

   private:
    friend class Graph;
    const Graph* graph;
    size_type index;

    /** Private constructor used to build and EdgeIterator pointing to a specific
     *  edge in a given graph.
     */
    EdgeIterator(const Graph* g, size_type i){
      graph = g;
      index = i;
    }
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Returns a begin edge iterator.
   * @return An EdgeIterator pointing to the first edge that was added to the Graph.
   *
   * Edges are enumerated in the order in which they were added to the graph.
   */
  edge_iterator edge_begin(){
    return EdgeIterator(const_cast<Graph*>(this), 0);
  }

  /** Returns the Edge the iterator is pointing to.
   * @return an EdgeIterator pointing to the end of the list of all edges added to the graph.
   */
  edge_iterator edge_end(){
    return EdgeIterator(const_cast<Graph*>(this), num_edges());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  std::vector<std::vector<size_type>> adjacency_;
  size_type next_node_uid;
  size_type next_edge_uid;

  struct internal_node {
    Point pos;
    size_type uid;
    node_value_type value;
  };

  struct internal_edge {
    size_type edge_uid;
    size_type node1_index;
    size_type node2_index;
    edge_value_type value;
  };

  // Helper function to get an EXISTING Edge from the graph.
  // HW1: modified to use iterators.
  Edge get_edge(const Node& a, const Node& b) const {
    // Only called when (a,b) is a valid edge!
    assert(has_node(a));
    assert(has_node(b));
    for(auto it = a.edge_begin(); it != a.edge_end(); ++it){
      if((*it).node2() == b){
        return *it;
      }
    }
    assert(false);
  }

};

#endif // CME212_GRAPH_HPP
