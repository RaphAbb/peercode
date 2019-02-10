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
 public:
  using node_value_type = V;    
  using size_type = unsigned;

 private:
  struct internal_node {
    Point position;
    size_type uid;
    node_value_type value;
  };

  struct internal_edge {
    size_type uid1;
    size_type uid2;
    size_type uid;
  };

  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;
  std::vector<std::vector<std::pair<size_type, bool>>> adj;
  size_type size_;
  size_type next_uid_;
  size_type num_edges_;
  size_type next_edge_uid_;
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
 public:
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

  /** Construct an empty graph. */
  Graph()
      : nodes_(), edges_(), size_(0), next_uid_(0), num_edges_(0), next_edge_uid_(0) {
  }

  ~Graph() = default;

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */

  class Node: private totally_ordered<Node> {
   public:
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return fetch().position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return uid_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    node_value_type& value (){
        return fetch().value;
    }
    const  node_value_type& value () const{
        return const_cast<const node_value_type&>(value());
    }
    // size_type degree() const;
    size_type degree() const {
      // HW0: YOUR CODE HERE
      return graph_->adj[uid_].size();
    }
    incident_iterator edge_begin() const {
        return IncidentIterator(graph_, index(), 0);
    }
    incident_iterator edge_end() const {
        return IncidentIterator(graph_, index(), degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
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
      const Point pos1 = position();
      const Point pos2 = n.position();
      if(pos1[0] < pos2[0]) return true;
      else if((pos1[0] == pos2[0]) && (pos2[1]< pos2[1])) return true;
      else if((pos1[1] == pos2[1]) && (pos2[2]< pos2[2])) return true;
      else return false;
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type uid_;
    /** Private Constructor */
    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }

    internal_node& fetch() const {
      for (size_type i = 0; i < graph_->size(); ++i)
        if (graph_->nodes_[i].uid == uid_)
          return graph_->nodes_[i];
      assert(false);
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return size_;
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    internal_node el;
    el.position = position;
    el.uid = next_uid_;
    el.value = value;
    adj.push_back(std::vector<std::pair<size_type,bool>>());
    nodes_.push_back(el);
    ++size_;
    ++next_uid_;

    return Node(this, next_uid_-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (this == n.graph_);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < size());
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
  class Edge: private totally_ordered<Edge> {
   public:

    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return (reverse ? graph_->node(fetch().uid2):graph_->node(fetch().uid1));
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return (reverse ? graph_->node(fetch().uid1):graph_->node(fetch().uid2));
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (node1() == e.node1() && node2() == e.node2()) || (node1() == e.node2() && node2() == e.node1());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      const Node n11 = node1();
      const Node n12 = node2();
      const Node n21 = e.node1();
      const Node n22 = e.node2();
      if ((n12 < n11) && (n22 < n21)){
          if (n21 < n11) return true;
          else if ((n11 == n21) && (n22 < n12)) return true;
          else return false;
      }
      else if ((n12 < n11) && (n21 < n22)){
          if (n22 < n11) return true;
          else if ((n11 == n22) && (n21 < n12)) return true;
          else return false;
      }
      else if ((n11 < n12) && (n22 < n21)){
          if (n21 < n12) return true;
          else if ((n12 == n21) && (n22 < n11)) return true;
          else return false;
      }
      else if ((n11 < n12) && (n21 < n22)){
          if (n22 < n12) return true;
          else if ((n12 == n22) && (n21 < n11)) return true;
          else return false;
      }
      return false;

    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type uid_;
    bool reverse;
    /** Private Constructor */
    Edge(const Graph* graph, size_type uid, bool invert)
        : graph_(const_cast<Graph*>(graph)), uid_(uid), reverse(invert) {
    }

    internal_edge& fetch() const {
      for (size_type i = 0; i < graph_->num_edges(); ++i)
        if (graph_->edges_[i].uid == uid_)
          return graph_->edges_[i];
      assert(false);
    }
    // Allow Graph to access Edge's private member data and functions.
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i, bool invert=false) const {
    // HW0: YOUR CODE HERE
    assert(i < num_edges());
    return Edge(this, i, invert);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    size_type uida(a.index());
    size_type uidb(b.index());
    for (size_type i = 0; i < num_edges(); ++i){
      if ((uida == edges_[i].uid1) && (uidb == edges_[i].uid2)) return true;
      if ((uida == edges_[i].uid2) && (uidb == edges_[i].uid1)) return true;
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
    // HW0: YOUR CODE HERE
    size_type uida(a.index());
    size_type uidb(b.index());
    for (size_type i = 0; i < num_edges(); ++i){
      if ((uida == edges_[i].uid1) && (uidb == edges_[i].uid2)) return edge(i);
      if ((uida == edges_[i].uid2) && (uidb == edges_[i].uid1)) return edge(i);
    }

    internal_edge el;
    el.uid1 = a.index();
    el.uid2 = b.index();
    el.uid = next_edge_uid_;
    adj[uida].push_back(std::make_pair(el.uid, false));
    adj[uidb].push_back(std::make_pair(el.uid, true));
    edges_.push_back(el);
    ++num_edges_;
    ++next_edge_uid_;

    return Edge(this, next_edge_uid_-1, false);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    edges_.clear();
    nodes_.clear();
    next_uid_ = 0;
    next_edge_uid_ = 0;
    size_ = 0;
    num_edges_ = 0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
  /** Return the ith node in the graph.
   *
   */
    Node operator*()  const{
        return graph->node(i);
    }
  /** Increments the iterator.
   * @post new (*it).index() = old (*it).index() + 1;
   *
   */
    node_iterator& operator++(){
        ++i;
        return *this;
    }
  /** Tests for equality.
   * @param[in] it: the other iterator
   *
   */
    bool  operator==(const  node_iterator& it)  const{
        return (it.i == i && graph == it.graph); 
    }

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }


   private:
    //Stores the nodes index
    Graph* graph;
    friend class Graph;
    size_type i;
    NodeIterator(const Graph* graph, size_type index)
        : graph(const_cast<Graph*>(graph)), i(index)  {
    }
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Returns first node in graph
   * @post (*node_begin()).index() = 0;
   */
  node_iterator node_begin () const{
      return NodeIterator(this, 0);
  }
  /** Returns interator of graph end
   * if (*it).index() == size() - 1 then ++it = node_end()
   */
  node_iterator node_end () const{
      return NodeIterator(this, size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator>{
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

   /** returns the current edge
   * (*it).node1().index() = it.node
   *
   */
    Edge operator*()  const{
        std::pair<size_type,bool> current(graph->adj[node][index]);
        return graph->edge(std::get<0>(current), std::get<1>(current));
    }
   /** Get the next adjacent edge
   *
   */
    incident_iterator& operator++(){
        ++index;
        return *this;
    }
   /** Check edge equality by checking all members
   *
   */
    bool operator==(const  incident_iterator& it)  const{
        return (it.index == index && graph == it.graph && node == it.node); 
    }

   private:
    size_type index;
    size_type node;
    Graph* graph;
    friend class Graph;
    IncidentIterator(const Graph* graph, size_type uid, size_type edge_index)
        : graph(const_cast<Graph*>(graph)), node(uid), index(edge_index)  {
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator> {
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
    Edge operator*()  const{
        return graph->edge(i);
    }
    edge_iterator& operator++(){
        ++i;
        return *this;
    }
    bool  operator==(const  edge_iterator& it)  const{
        return (it.i == i && graph == it.graph); 
    }

   private:
    size_type i;
    Graph* graph;
    friend class Graph;
    EdgeIterator(const Graph* graph, size_type index)
        : graph(const_cast<Graph*>(graph)), i(index)  {
    }

  };

  // HW1 #5: YOUR CODE HERE
  /** Returns first edge in graph
   */
  edge_iterator edge_begin () const{
      return EdgeIterator(this, 0);
  }
  /** Returns interator of graph edge end
   */
  edge_iterator edge_end () const{
      return EdgeIterator(this, num_edges());
  }
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
