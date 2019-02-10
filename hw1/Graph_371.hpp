#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <tuple>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

template <typename V>
/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Type of node values. */
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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
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
  class Node: private totally_ordered<Node> {
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
      // HW0: YOUR CODE HERE
      // Query the graph's list of nodes and retrieve the stored position.
      return this->graph_->nodes_[this->idx_].first;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return this->idx_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Return this node's internal value. */
    node_value_type& value() {
      // Query the graph's list of nodes and retrieve the stored value.
      return this->graph_->nodes_[this->idx_].second;
    }

    /** Return this node's internal value as an immutable constant. */
    const node_value_type& value() const{
      // Query the graph's list of nodes and retrieve the stored value.
      return this->graph_->nodes_[this->idx_].second;
    }

    /** Return the number of nodes that are adjacent to this node. */
    size_type degree() const {
      return (graph_->edges_[idx_]).size();
    }

    /**
     * @brief edge_begin
     * @return
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, idx_, 0);
    }

    //TODO
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, idx_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (this->graph_ == n.graph_ && this->idx_ == n.idx_);
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
      return (this->idx_ < n.idx_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Pointer back to Graph object to which the Node belongs.
    Graph* graph_;
    // Index of the Node in the Graph.
    size_type idx_ = 0;

    /** Private Constructor */
    Node(const Graph* graph, size_type idx)
        : graph_(const_cast<Graph*>(graph)), idx_(idx){
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return numNodes_;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] nodeValue The new node's internal value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& nodeValue = 0) {
    // HW0: YOUR CODE HERE
    nodes_.push_back(std::make_pair(position, nodeValue));
    numNodes_++;
    return Node(this, numNodes_ - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (n.graph_ == this && n.idx_ < numNodes_);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
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
  class Edge: totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(this->graph_, this->node1Idx_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(this->graph_, this->node2Idx_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((this->graph_ == e.graph_) &&
             ((this->node1Idx_ == e.node1Idx_ && this->node2Idx_ == e.node2Idx_) ||
             (this->node1Idx_ == e.node2Idx_ && this->node2Idx_ == e.node1Idx_)));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (this->graph_ == e.graph_ && this->edgeIdx_ < e.edgeIdx_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Pointer back to Graph object to which the Edge belongs.
    Graph* graph_;

    // Index of the Edge in the Graph.
    size_type edgeIdx_ = 0;

    // Indices of the two ends of the edge.
    size_type node1Idx_;
    size_type node2Idx_;

    /** Private Constructor */
    Edge(const Graph* graph, size_type idx, size_type node1Idx, size_type node2Idx)
        : graph_(const_cast<Graph*>(graph)), edgeIdx_(idx), node1Idx_(node1Idx), node2Idx_(node2Idx){
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return numEdges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(this, i, nodePairs_[i].first, nodePairs_[i].second);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    // A node is always adjacent to itself.
    if (a == b) return true;

    // Return false if @a a does not have adjacent nodes.
    if (edges_.find(a.idx_) == edges_.end()){
        return false;
    } else{
        for (auto const& p: edges_.at(a.idx_)){
          // Return true if @a a is adjacent to @a b.
          if (p.first == b.idx_) return true;
        }
    }
    // Return false if @a b was not in the list of @a a's adjacent nodes.
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
    // Return empty edge if @a a and @a b are same nodes.
    if (a == b) return Edge();

    // If there is an edge between @a a and @a b, find the edge index and return the corresponding edge.
    if (has_edge(a, b)){
        for (const auto& p: edges_[a.idx_]){
            if (p.first == b.idx_) return edge(p.second);
          }
    }

    // If there isn't an edge between @a a and @a b, create a new edge.
    edges_[a.idx_].push_back(std::make_pair(b.idx_, numEdges_));
    edges_[b.idx_].push_back(std::make_pair(a.idx_, numEdges_));
    if (a < b) nodePairs_.push_back(std::make_pair(a.idx_, b.idx_));
    else nodePairs_.push_back(std::make_pair(b.idx_, a.idx_));
    numEdges_++;
    return Edge(this, numEdges_ - 1, a.idx_, b.idx_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    edges_.clear();
    nodePairs_.clear();
    numNodes_ = 0;
    numEdges_ = 0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: equality_comparable<NodeIterator>{
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
    /** Dereferences iterator to return the current node.
     * @pre 0 <= idx_ <= size()
     * @return Node at index idx_
     */
    Node operator*() const {
      return graph_->node(idx_);
    }

    /** Increments the iterator to move to the next node.
     * @pre 0 <= idx_ <= size()
     * @return NodeIterator& with incremented idx_
     */
    node_iterator& operator++(){
      idx_++;
      return *this;
    }

    /** Check if two iterators are over the same graph and at the same node. */
    bool operator==(const node_iterator& nodeIter) const {
      return (graph_ == nodeIter.graph_ && idx_ == nodeIter.idx_);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_; // Graph to iterate.
    size_type idx_; // Index of node the iterator is currently at.

    /** Private Constructor */
    NodeIterator(const Graph* graph, size_type idx)
        : graph_(const_cast<Graph*>(graph)), idx_(idx){
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  /** Returns NodeIterator at the first node. */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Returns NodeIterator at one past the last node. */
  node_iterator node_end() const {
    return NodeIterator(this, numNodes_);
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: equality_comparable<IncidentIterator>{
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

    /** Returns the edge between the central node and the adjacent node.
     * @pre 0 <= adjIdx <= node1_.degree()
     * @return Edge with node1_ as the first node and the adjacent node as the second node
     */
    Edge operator*() const {
      // Find the index of the adjacent node and the index of the edge.
      size_type node2_ = graph_->edges_[node1_][adjIdx_].first;
      size_type edgeIdx = graph_->edges_[node1_][adjIdx_].second;
      return Edge(graph_, edgeIdx, node1_, node2_);
    }

    /** Increment the iterator to move to the next adjacent node.
     * @pre 0 <= adjIdx < node1_.degree()
     * @return IncidentIterator with incremented adjIdx
     */
    IncidentIterator& operator++() {
      adjIdx_++;
      return *this;
    }

    /** Check if two iterators are at the same central node and the same adjacent node of the same graph. */
    bool operator==(const IncidentIterator& incIter) const {
      return (graph_ == incIter.graph_) && (node1_ == incIter.node1_) && (adjIdx_ == incIter.adjIdx_);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_; // Graph to iterate.
    size_type node1_; // Index of central node.
    size_type adjIdx_; // Index of a node (node2) in the adjacency list of node1_.

    /** Private Constructor */
    IncidentIterator(const Graph* graph, size_type node1Idx, size_type adjIdx)
        : graph_(const_cast<Graph*>(graph)), node1_(node1Idx), adjIdx_(adjIdx){
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: equality_comparable<EdgeIterator>{
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

    /** Dereferences iterator to return edge corresponding to the current index.
     * @pre 0 <= edgeIdx_ <= numEdges_
     * @return Edge with corresponding index.
     */
    Edge operator*() const{
      return graph_->edge(edgeIdx_);
    }

    /** Increments the iterator to move to the next edge.
     * @pre 0 <= edgeIdx_ < numEdges_
     * @return EdgeIterator& with incremented edgeIdx_
     */
    EdgeIterator& operator++(){
      ++edgeIdx_;
      return *this;
    }

    /** Check if another EdgeIterator iterates over the same graph and is at the same edge. */
    bool operator==(const EdgeIterator& edgeIter) const{
      return graph_ == edgeIter.graph_ && edgeIdx_ == edgeIter.edgeIdx_;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_; // Graph to iterate over.
    size_type edgeIdx_; // Index of edge in nodePairs_ (vector of adjacent node pairs).

    /** Private Constructor */
    EdgeIterator(const Graph* graph, size_type edgeIdx)
        : graph_(const_cast<Graph*>(graph)), edgeIdx_(edgeIdx){
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  /** Returns EdgeIterator at the first edge. */
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0);
  }

  /** Returns EdgeIterator at one past the last edge. */
  edge_iterator edge_end() const{
    return EdgeIterator(this, numEdges_);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // Vector of Point and node value pairs.
  // e.g. {(node_a, node_a_value), (node_b, node_b_value), ... }
  std::vector<std::pair<Point, node_value_type>> nodes_;

  // Map with node index as keys, and vector of pairs with adjacent edge and edge index as values.
  // Note: Edge (a, b) gets stored twice in the map, one with key a_index and one with key b_index.
  // e.g. {node_a_index: (node_b_index, edge_index_1), (node_c_index, edge_index_2), ... }
  std::map<size_type, std::vector<std::pair<size_type, size_type>>> edges_;

  // Vector with adjecent nodes as pairs.
  // Note: Each edge(a, b) is stored only once. For each pair (a, b) in nodePairs_, a < b.
  std::vector<std::pair<size_type, size_type>> nodePairs_;

  size_type numNodes_ = 0; // Total number of nodes in the graph.
  size_type numEdges_ = 0; // Total number of edges in the graph.
};

#endif // CME212_GRAPH_HPP
