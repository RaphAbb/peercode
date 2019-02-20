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
template <typename V, typename E>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  
  // Index start value
  const static unsigned FIRST_INDEX {0};
  // Invalid index value (assuming largest unsigned int is not needed)
  const static unsigned INVALID_INDEX {unsigned(-1)};

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
  
  /** Type of node value.
      Template type of Graph, return type of Graph::Node::value(), and type of
      second (optional) argument of Graph::add_node(Point&, node_value_type&) */
  using node_value_type = V;

  /** Type of edge value.
   *  Template type of graph, return type of Graph::Edge::value(), and type of
   *  second (optional) argument of Graph::add_edge(Node&, Node&,
   *  edge_value_type&)
   */
  using edge_value_type = E;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : next_nid_(FIRST_INDEX), first_valid_nid_(FIRST_INDEX),
            next_eid_(FIRST_INDEX), first_valid_eid_(FIRST_INDEX) {
    // HW0: YOUR CODE HERE
    // HW2: Use static local variable to give unique graph ID to each graph
    // If an edge or node with identical IDs are compared from different graphs,
    // this will be used to define the < operator (note that copy/assignment
    // would copy this as well; if this is a desired feature, something other
    // than a staic local variable would need to be used)
    static size_type gid = 0;
    gid_ = ++gid;
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
    Node() : graph_(nullptr), nid_(INVALID_INDEX) {
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      
      // Node and Point vectors use the same indexing; index() checks if the
      // node is valid
      return graph_->points_.at(this->index());
    }

    //HW2
    /** Return a reference to this node's position. */
    Point& position() {
      return graph_->points_.at(this->index());
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      
      this->assert_valid();
      
      return graph_->node_indices_.at(nid_);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
    
    /** Returns a reference to the node value
     *
     * @pre The node must be valid.
     */
    node_value_type& value() {
      return graph_->node_values_.at(this->index());
    } 

    /** Returns a const reference to the node value
     *
     * @pre The node must be valid.
     */
    const node_value_type& value() const {
      return graph_->node_values_.at(this->index());
    } 
    
    /** Return the degree of the node (i.e. number of incident edges) */
    size_type degree() const {
      return size_type(graph_->nodes_to_edges_.at(this->index()).size());
    }
    
    /** Return an iterator to the beginning of the list of incident edges to the
     * node.
     *
     * @pre The node must be valid.
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_,this->index(),FIRST_INDEX);
    }
    
    /** Return an iterator to the end of the list of incident edges to the node.
     *
     * @pre The node must be valid.
     */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_,this->index(),this->degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      
      // index() asserts the node is valid
      return (graph_->has_node(n) ? this->index() == n.index() : false);
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

      // HW2: changed ordering to be based on the nid_; added checks for
      // comparing nodes between graphs
      // == operator checks if node is valid
      if (*this == n) return false;

      // check if they belong to the same graph, or if the nid's are equal
      if (graph_->has_node(n) || nid_ != n.nid_) return nid_ < n.nid_;

      // if they do not belong to the same graph and nid's are equal, use graph ID
      return graph_->gid_ < n.graph_->gid_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    // Pointer to owning Graph
    Graph* graph_;
    // Unique node ID
    size_type nid_;
    
    /** Private constructor for valid nodes */
    Node(const Graph* graph, size_type nid)
        : graph_(const_cast<Graph*>(graph)), nid_(nid) {
    }
    
    /** Helper method for asserting that the node is valid */
    void assert_valid() const {
      // Check that the node belongs to a graph
      assert(graph_ != nullptr);
      // Check to see if the node was removed from the graph, e.g. via clear()
      assert(graph_->node_indices_.at(nid_) != INVALID_INDEX);
    }
    
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return size_type(nodes_.size());
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value The new node's value (default is the default constructor)
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, 
                const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    
    points_.push_back(position);
    node_indices_.push_back(size());  // Index is based on order in vector
    node_values_.push_back(value);
    nodes_.push_back(Node(this,next_nid_));
    nodes_to_edges_.push_back(std::vector<size_type>());//empty incident vector
    ++next_nid_;
    
    return nodes_.back();
  }

  /** Remove a node from the graph.
   *
   * @param[in] n The node to be removed.
   * @return 1 (true) if the node was a valid node of the graph and was removed,
   *         0 (false) if the node was not in the graph.
   *
   * @pre The edges in the incidence list of @a n are valid (i.e. if any edge
   *      incident to @a n is removed, it is also removed from the incidence list)
   * @post @a new @a num_nodes() == @a old @a num_nodes() - 1
   * @post @a new @a num_edges() == @a old @a num_edges() - @a n.degree()
   * @post @a new @a nj.degree() == @a old @a nj.degree() -1, where @a nj are
   *       the nodes connected to @a n by edges
   *
   * This operation invalidates the following:
   *    1) NodeIterators pointing to the @a old @a graph.node_end() and
   *       @a old @a nodes_.back()
   *    2) IncidenceIterators @a old @a edge_end() and @a old @a edge_end() - 1
   *       for all nodes @a nj connected to @a n
   *    3) All IncidenceIterators of @a n
   *    4) All copies of the node @a n
   *    5) All copies of the edges in the incidence list of @a n
   *    6) The last @a k = @a n.degree() + 1 EdgeIterators
   *
   * This operation changes to following:
   *    1) NodeIterators pointing to the removed node will now point to the
   *       @a old @a nodes_.back() (as long as it is not the last node)
   *    2) If the node to be removed is not the last edge, it is overwritten by
   *       the last node and the (duplicate) last node is removed. This changes
   *       the index of nodes with
   *       @a nid = @a old @a nodes_.back().nid_
   *    3) The relative node indexing order is not preserved
   *
   * The complexity of this operation is O(@a k1* @a k2), where
   *    @a k1 = @a n.degree()
   *    @a k2 = max(@a nj.degree()) is the maximum degree of the nodes incident
   *        to @a n
   * due to the iteration over the incidence lists in remove_edge(Edge&). Note
   * that we are always removing the first edge in the incidence list of @a n,
   * so only the iterations over the connecting nodes will be linear in the
   * degree. Since the graph is sparse, this should be much less than O(num_nodes())
   */
  size_type remove_node(const Node& n) {
    // Check if n is a valid node
    if (has_node(n)) {
      // Remove all incident edges
      auto it = n.edge_begin();
      for (; it != n.edge_end(); it = n.edge_begin()) {
        // If the edge has already been removed but is not a part of the
        // incidence list, it was improperly removed
        assert(remove_edge(*it));
      }
      // Get the index of the node and the last node
      size_type i = n.index();
      size_type last = num_nodes();
      --last;
      if (i != last) {
        // If i is not the last node, overwrite it with the last node
        size_type nid = nodes_[last].nid_;
        node_indices_[nid] = i;
        nodes_[i] = nodes_[last];
        points_[i] = points_[last];
        node_values_[i] = node_values_[last];
        nodes_to_edges_[i] = nodes_to_edges_[last];
      }
      // Remove the (duplicate) last node
      nodes_.pop_back();
      points_.pop_back();
      node_values_.pop_back();
      nodes_to_edges_.pop_back();
      // Invalidate the removed node
      node_indices_[n.nid_] = INVALID_INDEX;
      return size_type(true);
    }
    return size_type(false);
  }

  /** Remove a node from the graph.
   *
   * @param n_it A NodeIterator of the node @a n to be removed.
   * @return If @a n_it is in the graph:
   *             If @a n_it != @a old node_end() - 1 && @a old @a num_nodes() != 1
   *                 return a NodeIterator to node @a i = @a old @a (*n_it).index()
   *             Else if @a n_it == @a old @a node_end() - 1
   *                 return @a new @a node_end() - 1
   *             Else
   *                 return @a new @a node_end()
   *         Else
   *             return @a node_end()
   * @pre The edges in the incidence list of @a n are valid (i.e. if any edge
   *      incident to @a n is removed, it is also removed from the incidence list)
   * @post @a new @a num_nodes() == @a old @a num_nodes() - 1
   * @post @a new @a num_edges() == @a old @a num_edges() - @a n.degree()
   * @post @a new @a nj.degree() == @a old @a nj.degree() -1, where @a nj are
   *       the nodes connected to @a n by edges
   *
   * This operation invalidates the following:
   *    1) NodeIterators pointing to the @a old @a graph.node_end() and
   *       @a old @a nodes_.back()
   *    2) IncidenceIterators @a old @a edge_end() and @a old @a edge_end() - 1
   *       for all nodes @a nj connected to @a n
   *    3) All IncidenceIterators of @a n
   *    4) All copies of the node @a n
   *    5) All copies of the edges in the incidence list of @a n
   *    6) The last @a k = @a n.degree() + 1 EdgeIterators
   *
   * This operation changes to following:
   *    1) NodeIterators pointing to the removed node will now point to the
   *       @a old @a nodes_.back() (as long as it is not the last node)
   *    2) If the node to be removed is not the last edge, it is overwritten by
   *       the last node and the (duplicate) last node is removed. This changes
   *       the index of nodes with
   *       @a nid = @a old @a nodes_.back().nid_
   *    3) The relative node indexing order is not preserved
   *
   * The complexity of this operation is O(@a k1* @a k2), where
   *    @a k1 = @a n.degree()
   *    @a k2 = max(@a nj.degree()) is the maximum degree of the nodes incident
   *        to @a n
   * due to the iteration over the incidence lists in remove_edge(Edge&). Note
   * that we are always removing the first edge in the incidence list of @a n,
   * so only the iterations over the connecting nodes will be linear in the
   * degree. Since the graph is sparse, this should be much less than O(num_nodes())
   */
  node_iterator remove_node(node_iterator n_it) {
    // Obtain the node of the iterator and call remove_node(Node&)
    auto n = *n_it;
    // Remove the node if it is in the graph
    size_type is_removed = remove_node(n);

    if (is_removed) {
      // Check if the node was not the only node left in the graph
      if (num_nodes()) {
        size_type i = node_indices_[n.nid_];
        // If the node was the last node, return the new last node
        if (i == num_nodes()) return NodeIterator(this,--(nodes_.end()) );
        // Else, the current iterator points to the node now at index i
        return n_it;
      }
    }
    // Return node_end() if the node was not found or was the last node left
    return node_end();
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    
    return (this == n.graph_ && node_indices_.at(n.nid_) != INVALID_INDEX);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    
    assert(FIRST_INDEX <= i && i < size());
    return nodes_[i];
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
    Edge() : graph_(nullptr), node1_(INVALID_INDEX), node2_(INVALID_INDEX), 
             eid_(INVALID_INDEX) {
      // HW0: YOUR CODE HERE
    }

    // HW2: Added functions
    /** Return a reference to the edge value
     * @pre The edge must be valid.
     */
    edge_value_type& value() {
      return graph_->edge_values_.at(this->index());
    }

    /** Return a const reference to the edge value
     * @pre The edge must be valid.
     */
    const edge_value_type& value() const {
      return graph_->edge_values_.at(this->index());
    }

    /** Return the current length of the Edge */
    double length() const {
      return norm_2(node1().position()-node2().position());
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      
      // Assert edge is valid
      this->assert_valid();
      // Get node
      return graph_->node(graph_->node_indices_.at(node1_));
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      this->assert_valid();
      return graph_->node(graph_->node_indices_.at(node2_));
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // node1() and node2() functions check if edge is valid
      return  ( ( node1() == e.node1() && node2() == e.node2() ) || 
                ( node1() == e.node2() && node2() == e.node1() ) );
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // Assume this and e belong to the same graph

      // HW2: changed ordering to be based on eid_; added checks for comparing
      // edges between graphs
      // == operator checks if edge is valid
      if (*this == e) return false;

      // check if they belong to the same graph or if the eid's are equal
      if (graph_ == e.graph_ || eid_ != e.eid_)   return eid_ < e.eid_;

      // they do not belong to the same graph and eid's are equal, use graph ID
      return graph_->gid_ < e.graph_->gid_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    
    // Pointer to owning graph
    Graph* graph_;
    // Node IDs to nodes that define the edge
    size_type node1_;
    size_type node2_;
    // Unique edge ID
    size_type eid_;
    
    /** Private constructor for valid edges */
    Edge(const Graph* graph, size_type node1, size_type node2, size_type eid) :
           graph_(const_cast<Graph*>(graph)), node1_(node1), node2_(node2),
           eid_(eid) {
    }
    
    /** Helper method that checks if the edge is valid and returns the index */
    size_type index() const {
      this->assert_valid();
            
      return graph_->edge_indices_[eid_];
    }
    
    /** Helper method that checks if the edge is valid */
    void assert_valid() const {
      assert(graph_ != nullptr && node1_ != INVALID_INDEX && 
             node2_ != INVALID_INDEX);
      assert(graph_->edge_indices_.at(eid_) != INVALID_INDEX);
      // Check if either node is invalid
      assert(graph_->node_indices_.at(node1_) != INVALID_INDEX &&
             graph_->node_indices_.at(node2_) != INVALID_INDEX);
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return size_type(edges_.size());
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    
    assert(FIRST_INDEX <= i && i < num_edges());
    return edges_[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    
    // Assert precondition
    assert(has_node(a));
    assert(has_node(b));
    
    // eid_ of edge connecting a and b, if one exists
    size_type eid = get_edge_eid(a, b);
    // Check for case that a removed edge was not properly removed from the map
    if (eid != INVALID_INDEX) {
      assert(edge_indices_.at(eid) != INVALID_INDEX);
      return true;
    }
    return false;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @param[in] a One node defining the edge
   * @param[in] b The other node defining the edge
   * @param[in] value The new edge's value (default is the default constructor)
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
  Edge add_edge(const Node& a, const Node& b,
                const edge_value_type& value = edge_value_type()) {
    // HW0: YOUR CODE HERE
    
    // has_edge() will assert that a, b exist in the graph; assert that they are
    // distinct
    assert(!(a == b));
    bool edge_exists = has_edge(a, b);
    if ( !edge_exists ) {
      // Edge does not exist and will be added
      edge_indices_.push_back(num_edges());
      
      // Add edge to node incident vector
      nodes_to_edges_[a.index()].push_back(next_eid_);
      nodes_to_edges_[b.index()].push_back(next_eid_);
      
      // Add edge and value
      edges_.push_back( Edge(this, a.nid_, b.nid_, next_eid_) );
      edge_values_.push_back(value);
      ++next_eid_;
      
      return edges_.back();
      
    }
    size_type eid = get_edge_eid(a, b); // We know this edge exists in the graph
    return edges_[edge_indices_[eid]];
  }

  /** Remove an edge from the graph.
   *
   * @param[in] n1 A node defining the edge.
   * @param[in] n2 The other node defining the edge.
   * @return 1 (true) if the edge existed and was removed, 0 (false) if the edge
   *         was not found.
   *
   * @pre @a n1 and @a n2 are valid nodes of the graph.
   * @pre If the edge defined by @a n1 and @a n2 is in the graph, the edge is in
   *      their incidence lists.
   * @post @a new @a degree() == @a old @a degree() - 1 for both @a n1 and @a n2
   * @post @a new @a num_edges() == @a old @a num_edges() - 1
   *
   * This operation invalidates the following:
   *    1) EdgeIterators pointing to the @a old @a graph.edge_end() and
   *       @a old @a edges_.back()
   *    2) IncidenceIterators @a old @a edge_end() and @a old @a edge_end() - 1
   *       for both @a n1 and @a n2
   *    3) All copies of the removed edge.
   *
   * This operation changes the following:
   *    1) EdgeIterators pointing to the removed edge will now point to the
   *       @a old last edge (as long as it is not the last edge)
   *    2) If the edge to be removed is not the last edge (i.e. edges_.back()),
   *       it is overwritten by the last element and the duplicate last edge is
   *       then removed. This changes the index of edges with
   *       @a eid = @a old @a edges_.back().eid_
   *    3) IncidenceIterators pointing to the removed @a eid will point to the
   *       @a old last @a eid of the incidence list (unless that @a eid was
   *       removed).
   *    4) The relative edge indexing order is not preserved.
   *
   * The complexity of this operation is O( max(@a n1.degree(),@a n2.degree()) )
   * due to the iterations over the incidence lists.
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    // Check if edge exists in the graph
    if (has_edge(n1,n2)) {
      // Get edge ID
      size_type eid = get_edge_eid(n1,n2);

      // Remove edge from the incidence lists of n1 and n2; if the edge is not
      // in the lists, the edge was improperly removed earlier from the lists
      assert(remove_incidence_edge(n1,eid));
      assert(remove_incidence_edge(n2,eid));

      // Get current edge index and last index
      size_type i = edge_indices_[eid];
      size_type last = num_edges();
      --last;

      // If this is not the last edge, swap the edges
      if (i != last) {
        size_type eid2 = edges_[last].eid_;
        edge_indices_[eid2] = i;
        edges_[i] = edges_[last];
        edge_values_[i] = edge_values_[last];
      }
      // Remove the last edge
      edges_.pop_back();
      edge_values_.pop_back();
      edge_indices_[eid] = INVALID_INDEX;
      return size_type(true);
    }
    return size_type(false);
  }

  /** Remove an edge from the graph.
   *
   * @param[in] e The edge to be removed.
   * @return 1 (true) if the edge existed and was removed, 0 (false) if the edge
   *         was not found.
   * @pre The nodes @a n1 and @a n2 of @a e are valid nodes of the graph.
   * @pre If @a e is in the graph, the edge is in the incidence lists of @a n1
   *      and @a n2.
   * @post @a new @a degree() == @a old @a degree() - 1 for both @a n1 and @a n2
   * @post @a new @a num_edges() == @a old @a num_edges() - 1
   *
   * This operation invalidates the following:
   *    1) EdgeIterators pointing to the @a old @a graph.edge_end() and
   *       @a old @a edges_.back()
   *    2) IncidenceIterators @a old @a edge_end() and @a old @a edge_end() - 1
   *       for both @a n1 and @a n2
   *    3) All copies of the edge @a e.
   *
   * This operation changes the following:
   *    1) EdgeIterators pointing to the removed edge will now point to the
   *       @a old last edge (as long as it is not the last edge)
   *    2) If the edge to be removed is not the last edge (i.e. edges_.back()),
   *       it is overwritten by the last element and the duplicate last edge is
   *       then removed. This changes the index of edges with
   *       @a eid = @a old @a edges_.back().eid_
   *    3) IncidenceIterators pointing to the removed @a eid will point to the
   *       @a old last @a eid of the incidence list (unless that @a eid was
   *       removed).
   *    4) The relative edge indexing order is not preserved.
   *
   * The complexity of this operation is O( max(@a n1.degree(),@a n2.degree()) )
   * due to the iterations over the incidence lists.
   */
  size_type remove_edge(const Edge& e) {
    // Obtain the nodes of the edges and call remove_edge(Node&,Node&)
    Node n1 = e.node1();
    Node n2 = e.node2();
    return remove_edge(n1,n2);
  }

  /** Remove an edge from the graph.
   *
   * @param[in,out] e_it An EdgeIterator of the edge to be removed.
   * @return If @a e_it is in the graph:
   *             If @a e_it != @a old @a edge_end() - 1 && @a old @a num_edges() != 1
   *                 return an EdgeIterator to edge @a i = @a old @a (*e_it).index()
   *             Else if @a e_it == @a old @a edge_end() - 1
   *                 return @a new @a edge_end() - 1
   *             Else
   *                 return @a new @a edge_end()
   *         Else
   *             return @a edge_end()
   * @pre The nodes @a n1 and @a n2 of @a *e_it are valid nodes of the graph.
   * @pre If @a e_it is in the graph, the edge is in the incidence lists of @a n1
   *      and @a n2.
   * @post @a new @a degree() == @a old @a degree() - 1 for both @a n1 and @a n2
   * @post @a new @a num_edges() == @a old @a num_edges() - 1
   *
   * This operation invalidates the following:
   *    1) EdgeIterators pointing to the @a old @a graph.edge_end() and
   *       @a old @a edges_.back()
   *    2) IncidenceIterators @a old @a edge_end() and @a old @a edge_end() - 1
   *       for both @a n1 and @a n2
   *    3) All copies of the edge @a *e_it.
   *
   * This operation changes the following:
   *    1) EdgeIterators pointing to the removed edge will now point to the
   *       @a old last edge (as long as it is not the last edge)
   *    2) If the edge to be removed is not the last edge (i.e. edges_.back()),
   *       it is overwritten by the last element and the duplicate last edge is
   *       then removed. This changes the index of edges with
   *       @a eid = @a old @a edges_.back().eid_
   *    3) IncidenceIterators pointing to the removed @a eid will point to the
   *       @a old last @a eid of the incidence list (unless that @a eid was
   *       removed).
   *    4) The relative edge indexing order is not preserved.
   *
   * The complexity of this operation is O( max(@a n1.degree(),@a n2.degree()) )
   * due to the iterations over the incidence lists.
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    // Obtain the edge of the iterator and call remove_edge(Edge&)
    auto e = *e_it;
    // Remove the edge if it is in the graph; if the edge is not the last edge,
    // the iterator will still be valid and will point to the new edge with
    // index i
    size_type is_removed = remove_edge(e);
    if (is_removed) {
      // If the edge was the only edge of the graph, return edge_end()
      if ( num_edges() ) {
        size_type i = edge_indices_[e.eid_];
        // If the edge was the last edge, return the new last edge
        if (i == num_edges()) return EdgeIterator(this,--(edges_.end()) );
        // Else, the current iterator points to the edge now at index i
        return e_it;
      }
    }
    // Return edge_end() if the edge was not found or was the last edge left
    return edge_end();
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    points_.clear();
    node_values_.clear();
    edges_.clear();
    nodes_to_edges_.clear();
    edge_values_.clear();
    
    size_type i;
    // Starting from first_valid_nid_ avoids unnecessary iterating
    for (i = first_valid_nid_; i < next_nid_; ++i) {
      node_indices_[i] = INVALID_INDEX;
    }
    first_valid_nid_ = next_nid_;
    
    for (i = first_valid_eid_; i < next_eid_; ++i) {
      edge_indices_[i] = INVALID_INDEX;
    }
    first_valid_eid_ = next_eid_;
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
    NodeIterator() : graph_(nullptr), node_it_() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const
    
    /** Dereference the iterator to obtain the corresponding node.
     *
     * @pre The iterator must be constructed by the graph.
     */
    Node operator*() const {
      this->assert_valid();
      return *node_it_;
    }
    
    /** Forward advance the iterator.
     *
     * @pre The iterator must be constructed by the graph.
     */
    NodeIterator& operator++() {
      this->assert_valid();
      ++node_it_;
      return *this;
    }
    
    /** Test whether this iterator and @a nit are equal (i.e. do they point to
     * the same node)
     * @param[in] nit A Graph::NodeIterator to compare
     */
    bool operator==(const NodeIterator& nit) const {
      return node_it_ == nit.node_it_;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    
    /** Private constructor for valid node iterators
     *
     * @param[in] graph Pointer to the owning graph
     * @param[in] node_it Iterator to the underlying node vector
     *
     */
    NodeIterator(const Graph* graph,
                 typename std::vector<Node>::const_iterator node_it) : 
                      graph_(const_cast<Graph*>(graph)),
                      node_it_(node_it) {}
                      
    /** Helper method to assert the iterator is valid.
     * This will not catch iterators invalidated by adding or removing nodes.
     */
    void assert_valid() const {
      assert(graph_ != nullptr && 
             node_it_ != typename std::vector<Node>::const_iterator());
    }
    
    // Pointer to owning Graph
    Graph* graph_;
    // Iterator that points to the graph_->nodes_ vector
    typename std::vector<Node>::const_iterator node_it_;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  
  /** Return a node iterator to the beginning of the nodes */
  node_iterator node_begin() const {
    return NodeIterator(this,nodes_.begin());
  }
  
  /** Return a node iterator to one past the end of the nodes */
  node_iterator node_end() const {
    return NodeIterator(this,nodes_.end());
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
    IncidentIterator() : graph_(nullptr), node_index_(INVALID_INDEX), 
                         node_to_edges_it_() {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const
    
    /** Dereference the iterator to obtain the corresponding edge.
     *
     * @pre The iterator must be constructed by a valid node.
     * @post @a node1 of the returned edge should equal nodes_[node_index_]
     *
     */
    Edge operator*() const {     
      this->assert_valid(); 
      // Get the edge ID
      size_type eid = *node_to_edges_it_;
      
      // Obtain Edge and nodes
      Edge edge = graph_->edges_.at(graph_->edge_indices_.at(eid));
      Node node1 = edge.node1();
      Node node2 = edge.node2();
      
      // Return an edge such that edge.node1() is the owning node
      if ( node1 == graph_->nodes_.at(node_index_) ) {
        return Edge(graph_, node1.nid_, node2.nid_, eid);
      }
      return Edge(graph_, node2.nid_, node1.nid_, eid);
    }
    
    /** Forward advance the iterator.
     *
     * @pre The iterator must be constructed by a valid node.
     */
    IncidentIterator& operator++() {
    this->assert_valid();
      ++node_to_edges_it_;
      return *this;
    }
    
    /** Test whether this iterator and @a iit are equal (i.e. do they point to
     * the same edge of the same incident vector)
     * @param[in] iit A Graph::IncidentIterator to compare
     */
    bool operator==(const IncidentIterator& iit) const {
      return (node_index_ == iit.node_index_ &&
             node_to_edges_it_ == iit.node_to_edges_it_);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    
    /** Private constructor for valid incident iterators
     *
     * @param[in] graph Pointer to the owning graph
     * @param[in] node_index Index of the owning node
     * @param[in] it_index Index of iterator in the incident edge list
     *
     * @pre 0 <= @a node_index < num_nodes()
     * @pre 0 <= @a it_index <= @a size, where @a size is the number of incident
     * edges of node @a node_index
     *
     */
    IncidentIterator(const Graph* graph, size_type node_index,
                     size_type it_index) : graph_(const_cast<Graph*>(graph)),
                     node_index_(node_index) {
      
      // First precondition is will be asserted when the node index() function
      // is called
      assert(FIRST_INDEX <= it_index && 
             it_index <= graph_->nodes_to_edges_[node_index_].size());
             
      if (it_index == graph_->nodes_to_edges_[node_index_].size()) {
        node_to_edges_it_ = graph_->nodes_to_edges_[node_index_].end();
      }
      else {
        node_to_edges_it_ = graph_->nodes_to_edges_[node_index_].begin();
        node_to_edges_it_ += it_index;
      }
    }
    
    /** Helper method to assert that the iterator is valid.
     * This will not catch iterators invalidated by adding or removing edges.
     */
    void assert_valid() const {
      assert(graph_ != nullptr && node_index_ != INVALID_INDEX &&
             node_to_edges_it_ != std::vector<size_type>::const_iterator());
    }

    /** Helper method to return the edge ID of the edge to which the iterator
     * points.
     * @pre The iterator must be constructed by a valid node.
     */
    size_type get_eid() const {
      this->assert_valid();

      return *node_to_edges_it_;
    }
        
    // Pointer to owning graph
    Graph* graph_;
    // Index of owning node
    size_type node_index_;
    // Iterator that points to the graph_->nodes_to_edges_[node_index_] vector
    typename std::vector<size_type>::const_iterator node_to_edges_it_;
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : equality_comparable<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() : graph_(nullptr), edge_it_() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const
    
    /** Dereference the iterator to obtain the corresponding node.
     *
     * @pre The iterator must be constructed by the graph.
     */
    Edge operator*() const {
      this->assert_valid();
      return *edge_it_;
    }
    
    /** Forward advance the iterator
     *
     * @pre The iterator must be constructed by the graph.
     */
    EdgeIterator& operator++() {
      this->assert_valid();
      ++edge_it_;
      return *this;
    }
    
    /** Test whether this iterator and @a eit are equal (i.e. do they point to
     * the same edge).
     *
     * @param[in] eit A Graph::EdgeIterator to compare
     */
    bool operator==(const EdgeIterator& eit) const {
      return edge_it_ == eit.edge_it_;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    
    /** Private constructor for valid edge iterators.
     *
     * @param[in] graph Pointer to the owning graph
     * @param[in] edge_it Iterator to the underlying edge vector
     *
     */
    EdgeIterator(const Graph* graph,
                 typename std::vector<Edge>::const_iterator edge_it) :
                      graph_(const_cast<Graph*>(graph)),
                      edge_it_(edge_it) {}
    /** Helper method to assert that the iterator is valid.
     * This will not catch iterators invalidated by insertion or removal of
     * edges.
     */
    void assert_valid() const {
      assert(graph_ != nullptr && 
             edge_it_ != typename std::vector<Edge>::const_iterator());
    }
    
    // Pointer to owning graph
    Graph* graph_;
    // Iterator that points to the graph_->edges_ vector
    typename std::vector<Edge>::const_iterator edge_it_;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  
  /** Return an edge iterator to the beinning of the edges */
  edge_iterator edge_begin() const {
    return EdgeIterator(this,edges_.begin());
  }
  
  /** Return an edge iterator to one past the end of the edges */
  edge_iterator edge_end() const {
    return EdgeIterator(this,edges_.end());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  
  // Contains graph nodes
  std::vector<Node>  nodes_;
  // Contains node positions
  std::vector<Point> points_;
  // Vector that maps node IDs to node index
  std::vector<size_type> node_indices_;
  // Next node ID for assignment
  size_type  next_nid_;
  // First index for node_indices_ that corresonds to a valid node
  size_type first_valid_nid_;
  // Contains grah edges
  std::vector<Edge> edges_;
  // Vector that maps edge IDs to edge index
  std::vector<size_type> edge_indices_;
  // Next edge ID for assignment
  size_type next_eid_;
  // First index for edge_indices_ that corresponds to a valid node
  size_type first_valid_eid_;
  // Vector mapping a node to the edges connected to it; all edges will thus be
  // listed twice
  std::vector< std::vector<size_type> > nodes_to_edges_;
  // Contains node values
  std::vector<node_value_type> node_values_;
  // Contains edge values
  std::vector<edge_value_type> edge_values_;
  // Contains unique graph ID for Node and Edge comparisons between graphs;
  // since no custom copy/assignment functions are written, doing so will also
  // copy this ID; if we plan on copying graphs and want them to be treated as
  // separate for the case of matching nid_ and eid_ in two graphs, they will
  // need to be implemented
  size_type gid_;

  /** Helper method that returns the edge id if nodes @a a and @a b are
   *  connected by an edge, and @a INVALID_INDEX otherwise.
   */
  size_type get_edge_eid(const Node& a, const Node& b) const {
    
    for (auto it = a.edge_begin(); it != a.edge_end(); ++it) {
      if ((*it).node2() == b)
        return (*it).eid_;
    }
    return INVALID_INDEX;
  }

  /** Helper method to remove edges from the node incidence list.
   *
   * @param n The node owning the incidence list from which the edge is being
   *          removed
   * @param eid The edge ID of the edge being removed.
   * @return 1 (true) if the edge was removed, 0 (false) if the edge was not
   *         found.
   * @pre The edge @a eid is in the incidence list of @a n.
   * @post @a new @a n.degree() == @a old @a n.degree() - 1
   *
   * This operation invalidates IncidenceIterators @a old @a n.edge_end() and
   * iterators to the @a old last element (i.e. @a old @a n.edge_end()-1)
   *
   * This operation overwrites the element to be removed by the last element in
   * the list, thus changing the internal order of the edges in the incidence
   * list.
   */
  size_type remove_incidence_edge(const Node& n, const size_type eid) {
    // Search for edge in incidence list (don't use iterators so we can swap;
    // the const_iterator used in the IncidentIterators does not allow reassignemt)
    std::vector<size_type> list = nodes_to_edges_[n.index()];
    for (size_type i = 0; i < size_type(list.size()); ++i) {
      if (list[i] == eid) {
        // If in the list, check if it is the last element
        size_type last = list.back();
        if (last != eid) {
          // If not, swap the last element and the current element
          nodes_to_edges_[n.index()][i] = last;
        }
        // Remove the last element
        nodes_to_edges_[n.index()].pop_back();
        return size_type(true);
      }
    }
    return size_type(false);
  }

};

#endif // CME212_GRAPH_HPP
