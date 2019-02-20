#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <map>
#include <vector>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

template <typename V, typename E = double>
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
  struct internal_node;
  struct internal_edge;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Type of node values. */
  using node_value_type = V;

  /** Type of edge values. */
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
      return graph_->nodes_[idx_].position;
    }

    /** Return this node's position in a modifiable state. */
    Point& position(){
      return graph_->nodes_[idx_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return idx_;
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
      return graph_->nodes_[idx_].value;
    }

    /** Return this node's internal value as an immutable constant. */
    const node_value_type& value() const{
      // Query the graph's list of nodes and retrieve the stored value.
      return graph_->nodes_[idx_].value;
    }

    /** Return the number of nodes that are adjacent to this node. */
    size_type degree() const {
      return (graph_->adjacencies_[idx_]).size();
    }

    /** Returns iterator for first edge.
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, idx_, 0);
    }

    /** Returns iterator for last edge.
     */
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
      return this->graph_ < n.graph_ ||
            (this->graph_ == n.graph_ && this->idx_ < n.idx_);
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
  Node add_node(const Point& position,
                const node_value_type& nodeValue = node_value_type()) {
    // HW0: YOUR CODE HERE
    // Make and store a new internal node object
    internal_node newNode = {numNodes_, position, nodeValue};
    nodes_.push_back(newNode);
    numNodes_++;
    return Node(this, numNodes_ - 1);
  }

  /** Remove a given node and all of its incident edges from the graph.
   * @param[in] n node to remove..
   * @post if has_node(@a n), new size() = old size() - 1
   * @post !has_node(n)
   * @post Edges incident to the node are removed, i.e. new num_edges() = old num_edges() - n.degree()
   * @return true(1) if the node existed in the graph and was successfully removed;
   *         false(0) if the node was not in the graph.
   *
   */
  size_type remove_node(const Node& n){
    if (!has_node(n)) return 0;

    // Remove all edges incident to the node
    auto ii = n.edge_begin();
    while (ii != n.edge_end()){
        ii = remove_edge(ii);
    }

    size_type oldIndex = size() - 1;
    size_type newIndex = n.idx_;

    // Update node idx in adjacency list
    for (auto& kv: adjacencies_){
      for (auto& p: kv.second){
          if (p.first == oldIndex) p.first = newIndex;
      }
    }
    adjacencies_[newIndex] = adjacencies_[oldIndex];
    adjacencies_.erase(oldIndex);

    // Update node idx in edge list
    for (auto& e: edges_){
        if (e.node1Idx == oldIndex) e.node1Idx = newIndex;
        else if (e.node2Idx == oldIndex) e.node2Idx = newIndex;
    }

    // Remove from node list by switching with last internal node
    nodes_[newIndex] = nodes_[oldIndex];
    nodes_[newIndex].idx = n.idx_;
    nodes_.pop_back();
    numNodes_--;

    return 1;
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
      return Node(graph_, node1Idx_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node2Idx_);
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
      return this->graph_ < e.graph_ ||
            (this->graph_ == e.graph_ && this->edgeIdx_ < e.edgeIdx_);
    }

    /** Return length of edge, i.e. Euclidean distance between the end-nodes.
     */
    double length() const{
      return norm(node1().position() - node2().position());
    }

    /** Return value of an edge, modifiable. */
    edge_value_type& value(){
      return (graph_->edges_[edgeIdx_]).value;
    }

    /** Return value of an edge as a const. */
    const edge_value_type& value() const{
      return (graph_->edges_[edgeIdx_]).value;
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
    return Edge(this, i, edges_[i].node1Idx, edges_[i].node2Idx);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    // Return false if @a a does not have adjacent nodes.
    if (adjacencies_.find(a.idx_) == adjacencies_.end()){
        return false;
    } else{
        for (auto const& p: adjacencies_.at(a.idx_)){
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& edgeValue = edge_value_type()) {
    // HW0: YOUR CODE HERE
    // Return empty edge if @a a and @a b are same nodes.
    if (a == b) return Edge();

    // If there is an edge between @a a and @a b, find the edge index and return the corresponding edge.
    if (has_edge(a, b)){
        for (const auto& p: adjacencies_[a.idx_]){
            if (p.first == b.idx_) return edge(p.second);
          }
    }

    // If there isn't an edge between @a a and @a b, create a new edge.
    adjacencies_[a.idx_].push_back(std::make_pair(b.idx_, numEdges_));
    adjacencies_[b.idx_].push_back(std::make_pair(a.idx_, numEdges_));
    internal_edge newEdge = {numEdges_, a.idx_, b.idx_, edgeValue};
    edges_.push_back(newEdge);
    numEdges_++;
    return Edge(this, numEdges_ - 1, a.idx_, b.idx_);
  }

  /** Remove an edge from the graph, given the two end nodes.
   * @param[in] a One end node of the edge to remove.
   * @param[in] b Other end node of the edge to remove.
   * @post if has_edge(a, b), new @a a.degree() = old @a a.degree() - 1;
   * @post if has_edge(a, b), new @a b.degree() = old @a b.degree() - 1;
   * @post if has_edge(a, b), new num_edges() = old num_edges() - 1;
   * @post !has_edge(a, b)
   * @post The edge is removed from the graph's list of internal edges and adjacencies.
   * @return true(1) if the edge was in the graph and was successfully removed;
   *         false(0) if the edge wasn't in the graph.
   */
  size_type remove_edge(const Node& a, const Node& b){
    if (!has_edge(a, b)) return 0;
    // Get edge index
    size_type edgeIdx;
    for (const auto& p: adjacencies_[a.idx_]){
        if (p.first == b.idx_){
          edgeIdx = p.second;
        }
      }

    // Switch with the last element
    edges_[edgeIdx] = edges_[num_edges() - 1];
    edges_[edgeIdx].idx = edgeIdx;
    edges_.pop_back();

    remove_adjacency(a.idx_, b.idx_);
    remove_adjacency(b.idx_, a.idx_);

    // Update edge that replaced the index of the removed edge.
    size_type node1Idx = edges_[edgeIdx].node1Idx;
    size_type node2Idx = edges_[edgeIdx].node2Idx;
    for (auto& p: adjacencies_[node1Idx]){
      if (p.first == node2Idx) p.second = edgeIdx;
    }
    for (auto& p: adjacencies_[node2Idx]){
      if (p.first == node1Idx) p.second = edgeIdx;
    }

    numEdges_--;

    return 1;
  }

  /** Remove a given edge from the graph.
   * @param[in] e Edge to remove.
   * @post if has_edge(e), new @a e.node1().degree() = old @a e.node1().degree() - 1;
   * @post if has_edge(e), new @a e.node2().degree() = old @a e.node2().degree() - 1;
   * @post if has_edge(e), new num_edges() = old num_edges() - 1;
   * @post The edge is removed from the graph's list of internal edges and adjacencies.
   * @return true(1) if the edge was in the graph and was successfully removed;
   *         false(0) if the edge wasn't in the graph.
   */
  size_type remove_edge(const Edge& e){
    if (!has_edge(e.node1(), e.node2())) return 0;
    remove_edge(e.node1(), e.node2());
    return 1;
  }

  /** Helper function to remove an adjacency from a center node to an adjacent node.
   *  Switches the pair of adjacent node and edge index with the last pair stored for the
   *  center node, then pops last element.
   * @param[in] centerNode Center node for the adjacency (i.e. key in the adjacency map)
   * @param[in] adjNode    Adjacent node for the adjacency
   */
  void remove_adjacency(const size_type centerNode, const size_type adjNode){
    size_type cNodeDegree = adjacencies_[centerNode].size();
    for (size_type i = 0; i < cNodeDegree; ++i){
        if (adjacencies_[centerNode][i].first == adjNode){
            adjacencies_[centerNode][i] = adjacencies_[centerNode][cNodeDegree - 1];
            adjacencies_[centerNode].pop_back();
            break;
        }
    }
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
    adjacencies_.clear();
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

  /** Removes node at a given NodeIterator.
   * @param[i] n_it NodeIterator pointing at the node to remove.
   * @pre n_it must be dereferenceable.
   * @return NodeIterator pointing at the previously next node.
   */
  node_iterator remove_node(node_iterator n_it){
    auto n = *n_it;
    remove_node(n);
    return n_it;
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
      size_type node2_ = graph_->adjacencies_[node1_][adjIdx_].first;
      size_type edgeIdx = graph_->adjacencies_[node1_][adjIdx_].second;
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

  /** Removes edge at a given EdgeIterator.
   * @param[in] e_it EdgeIterator pointing at the edge to remove.
   * @pre e_it must be dereferenceable.
   * @return EdgeIterator pointing at the previously next edge.
   */
  edge_iterator remove_edge(edge_iterator e_it){
    auto e = *e_it;
    remove_edge(e);
    return e_it;
  }

  /** Remove edge at a given IncidentIterator.
   * @param[in] i_it IncidentIterator pointing at the edge to remove.
   * @pre e_it must be dereferenceable.
   * @return IncidentIterator pointing at the previously next edge.
   */
  incident_iterator remove_edge(incident_iterator i_it){
    auto e = *i_it;
    remove_edge(e);
    return i_it;
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct internal_node {
    size_type idx;
    Point position;
    node_value_type value;
  };

  struct internal_edge {
    size_type idx;
    size_type node1Idx;
    size_type node2Idx;
    edge_value_type value;
  };

  // Vector of internal_node objects.
  std::vector<internal_node> nodes_;

  // Map with node index as keys, and vector of pairs with adjacent edge and edge index as values.
  // Note: Edge (a, b) gets stored twice in the map, one with key a_index and one with key b_index.
  // e.g. {node_a_index: (node_b_index, edge_index_1), (node_c_index, edge_index_2), ... }
  std::map<size_type, std::vector<std::pair<size_type, size_type>>> adjacencies_;

  // Vector with internal_edge objects.
  // Note: Each edge(a, b) is stored only once, unlike adjacencies_.
  std::vector<internal_edge> edges_;

  size_type numNodes_ = 0; // Total number of nodes in the graph.
  size_type numEdges_ = 0; // Total number of edges in the graph.
};

#endif // CME212_GRAPH_HPP
