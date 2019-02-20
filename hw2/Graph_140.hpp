#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <map>
#include <tuple>
#include <vector>

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

  /** Predeclaration of Internal Node type. */
  class Internal_Node;
  /** Predeclaration of Internal Edge type. */
  class Internal_Edge;
  
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V, E>;

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

  /** Type of node value */
  typedef V node_value_type;
  /** Type of edge value */
  typedef E edge_value_type;

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
      // HW0: YOUR CODE HERE
      return graph_->nodes_.at(this-> n_uid_).position_;
    }

    Point& position() {
    	return graph_->nodes_.at(this-> n_uid_).position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      //return this-> idx_;
    	return graph_->nodes_.at(this->n_uid_).n_idx_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;


    /** Retrieve a node's value, allows value to be set
    * 
    * Complexity: O(1) amortized operations. 
    */
    node_value_type& value() {
      return graph_->nodes_.at(this-> n_uid_).val_;
    }

    /** Retrieve a node's value
    *
    * Complexity: O(1)
    */
    const node_value_type& value() const {
      return graph_->nodes_.at(this-> n_uid_).val_;
    }

    /** Return node's degree i.e., the number of incident nodes
    * 
    * Complexity: O(1)
    */
    size_type degree() const {
      return graph_->connections_.at(this-> n_uid_).size();
    }

    /** Return incident iterator begin instance 
    *
    * Complexity: O(1)
    */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_->connections_.at(this-> n_uid_).begin(), this->n_uid_, this->graph_);

    }

    /** Return incident iterator end instance
    *
    * Complexity: O(1) 
    */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_->connections_.at(this-> n_uid_).end(), this->n_uid_, this->graph_);
    }


    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if ((n.graph_ == this->graph_) and (n.n_uid_ == this->n_uid_)) {
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
      // HW0: YOUR CODE HERE
      if (graph_ != n.graph_) {
        return graph_ < n.graph_;
      }
      else{
        return index() < n.index();
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Attribute to keep track of graph of node
    graph_type* graph_;

    // Attribute to keep track of index of node
    //size_type idx_;

    // Attribute to keep track of unique identifier of node
    size_type n_uid_;

    // Constructor for Node
    Node(graph_type* g, size_type n_uid)
      : graph_(g), n_uid_(n_uid) {
      }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    //return this->next_idx;
    return this->n_i2u_.size();
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
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    // HW0: YOUR CODE HERE

    // Add internal node to the Graph
    this->nodes_.push_back( Internal_Node(this, this-> n_i2u_.size() ,this->next_n_uid_, position, val) );

    // Add unique identifer to n_i2u_
    this->n_i2u_.push_back(next_n_uid_);

    // Increment next_n_uid_
    this->next_n_uid_++;

    // Add node and empty map to adjacency map
    this->connections_[next_n_uid_-1] = std::map<size_type, size_type>();

    // Create and return node instance
    return Node(this, next_n_uid_ - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return this == n.graph_;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE

    // Test for precondition
    // Note: no need to check 0 <= i, as 0 <= size_type always true
    if (i < this->num_nodes()) {
    	size_type uid = this->n_i2u_.at(i);
      return Node(const_cast<graph_type*>(this), uid);
    }
    // If precondition fails, return invalid node
    return Node();        // Invalid node
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, this->n1_uid_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, this->n2_uid_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if ((e.graph_ == this->graph_) and (e.e_uid_ == this->e_uid_)) {
        return true;
      }
      return false;
    }

    /** Return index of an edge
     *
     * Complexity: O(1)
     */
    size_type index() const {
    	return graph_->edges_.at(e_uid_).e_idx_;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
    	if (graph_ != e.graph_) {
    		return graph_ < e.graph_;
    	}
    	else{
    		return index() < e.index();
    	}
    }

    /** Returns length of edge, defined as euclidean instance between its two nodes
     *
     * Complexity: O(1)
     */
    double length() const {
    	return norm(node1().position() - node2().position());
    }

    /** Retrieve's an edges value, allows it to be set
     *
     * Complexity: O(1)
     */
    edge_value_type& value() {
    	return graph_->edges_.at(this-> e_uid_).val_;
    }

    /** Retrieve's an edges value
     * 
     * Complexit: O(1)
     */
    const edge_value_type& value() const {
    	return graph_->edges_.at(this-> e_uid_).val_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Attribute to keep track of graph of edge
    graph_type* graph_;

    // Attribute to keep track of index of edge
    size_type e_uid_;

    // Attribute to keep track of nodes of edge
    size_type n1_uid_;
    size_type n2_uid_;

    // Constructor for Edge
    Edge(graph_type* g, size_type uid, size_type n1, size_type n2)
      : graph_(g), e_uid_(uid), n1_uid_(n1), n2_uid_(n2) {
      }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    //return edges.size();
    return e_i2u_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(1)
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE

    // Check for precondition
    // Note: no need to check 0 <= i, as 0 <= size_type always true
    if (i < num_edges()) {
    	size_type uid = e_i2u_.at(i);
    	return Edge(const_cast<graph_type*>(this), uid, edges_.at(uid).n_uid1_, edges_.at(uid).n_uid2_);
    }
    // If precondition fails, return an invalid edge
    return Edge();        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully lessls
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE

    // Test for precondition
    if (has_node(a) and has_node(b)) {
      size_type uid1 = a.n_uid_;
      size_type uid2 = b.n_uid_;

      if (connections_.at(uid1).find(uid2) == connections_.at(uid1).end()) {
        return false;
      }
      else {
        return true;
      }
    }
    // If precondition fails, throw an exception
    else {
      throw std::invalid_argument( "Invalid nodes received" );
    }
    
    // Otherwise, if edge not in graph, return false
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val = edge_value_type()) {
    // HW0: YOUR CODE HERE

    //Test for precondition 
    if (a == b) {
      throw std::runtime_error("Self loop illegal");
    }
    if (!has_node(a) or !has_node(b)) {
      throw std::runtime_error("Invalid nodes for edge addition");
    }

    // If edge present in graph
    if (has_edge(a,b)) {
      size_type uid = connections_.at(a.n_uid_).at(b.n_uid_);
      return Edge(this, uid, a.n_uid_, b.n_uid_);
    }

    // Push back unqiue edge id to e_i2u_ vector
    this->e_i2u_.push_back(next_e_uid_);

    // Otherwise, add edge to graph, return edge object
    this->edges_.push_back( Internal_Edge(this, e_i2u_.size()-1, next_e_uid_, a.n_uid_, b.n_uid_, val) );

    edge_type new_edge =Edge(this, next_e_uid_, a.n_uid_, b.n_uid_);

    // Add edge to our adjacency map
    connections_[a.n_uid_][b.n_uid_] = next_e_uid_;
    connections_[b.n_uid_][a.n_uid_] = next_e_uid_;

    this-> next_e_uid_++;

    return new_edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE

    // Delete all nodes stored in graph
    nodes_.clear();

    // Reset next_idx (num_nodes) to 0
    next_n_uid_ = 0;

    // Delete all edges stored in adjacency list
    connections_.clear();

    // Delete all edges stored in vector
    edges_.clear();

    // Reset next_idy (num_edges) to 0
    next_e_uid_ = 0;

    // Delete all i2u_ vectors
    n_i2u_.clear();
    e_i2u_.clear();
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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) cons


    /** Return node corresponding to iterator position
    *
    * Complexity: O(1)
    */
    Node operator*() const {
    	size_type uid = *it_;
    	return Node(graph_, uid);
    }

    /** Increments value of node iterator
    * 
    * Complexity: O(1)
    */
    node_iterator& operator++(){
      it_++;
      return *this;
    }

    /** Test whether this node iterator and @a n are equal
    *
    * Equal node iterators have the same iterator and graph attribute
    * Complexity: O(1)
    */
    bool operator==(const NodeIterator& n) const {
      return (it_ == n.it_) & (graph_ == n.graph_);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE

    graph_type* graph_;
    typename std::vector<size_type>::const_iterator it_;

    NodeIterator(graph_type* graph, typename std::vector<size_type>::const_iterator it)
      : graph_(graph), it_(it) {
      }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Return node iterator begin instance
  *
  * Complexity: O(1)
  */
  node_iterator node_begin() const {
    node_iterator b = NodeIterator(const_cast<graph_type*>(this), n_i2u_.begin());
    return b;
  }

  /** Return node iterator end instance
  *
  * Complexity: O(1)
  */
  node_iterator node_end() const {
    node_iterator e = NodeIterator(const_cast<graph_type*>(this), n_i2u_.end());
    return e;
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
    IncidentIterator() {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** Retern Edge instance of current incident edge
    *
    * Complexity: O(1)
    */
    Edge operator*() const {
      size_type node2_uid = it_->first;
      size_type edge_uid = it_->second;
      return Edge(this->graph_, edge_uid, this->n_uid_, node2_uid);
    }

    /** Increments incident iterator instance
    *
    * Complexity: O(1)
    */
    IncidentIterator& operator++() {
      it_++;
      return *this;
    }

    /** Test whether this incident iterator and @a x are equal
    *
    * Equal node iterators have the same iterator, node id, and graph attribute
    * Complexity: O(1)
    */
    bool operator==(const IncidentIterator& x) const {
      if ( (x.it_ == this->it_) && (x.n_uid_ == this->n_uid_) && (x.graph_ == this->graph_) ){
        return true;
      }
      else{
        return false;
      }
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE

    typename std::map<size_type, size_type>::const_iterator it_;

    // Attribute to keep track of node index
    size_type n_uid_;

    // Attribute to keep track of graph of node
    graph_type* graph_;

    IncidentIterator(typename std::map<size_type, size_type>::const_iterator it, size_type n_uid, graph_type* g)
      : it_(it), n_uid_(n_uid), graph_(g) {
      }
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /** Return Edge instance of current edge iterator
    *
    * Complexity: O(1)
    */
    Edge operator*() const {
      size_type uid = *it_;
      Internal_Edge e = graph_->edges_.at(uid);
      return Edge(graph_, uid, e.n_uid1_, e.n_uid2_);
    }

    /** Increment edge iterator
    *
    * Complexity: O(1)
    */
    EdgeIterator& operator++() {
      it_++;
      return *this;
    }

    /** Test whether this edge iterator and @a e are equal
    *
    * Equal edge iterators have the same iterator attribute
    * Complexity: O(1)
    */
    bool operator==(const EdgeIterator& e) const {
      return (it_ == e.it_) & (graph_ == e.graph_);
    }


   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE

    graph_type* graph_;
    typename std::vector<size_type>::const_iterator it_;

    EdgeIterator(graph_type* graph, typename std::vector<size_type>::const_iterator it)
      : graph_(graph), it_(it) {

      }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Returns edge iterator begin instance
  *
  * Complexity: O(1)
  */
  edge_iterator edge_begin() const {
    edge_iterator b = EdgeIterator(const_cast<graph_type*>(this), e_i2u_.begin());
    return b;
  }

  /** Returns edge iterator end instance
  *
  * Complexity: O(1)
  */
  edge_iterator edge_end() const {
    edge_iterator e = EdgeIterator(const_cast<graph_type*>(this), e_i2u_.end());
    return e;
  }



  // HW2 



  /** Remove an edge to the graph, return the edge_begin() iterator
   * @param[in] e_it An iterator to an edge
   * @pre @a e_it is a valid edge iterator of the graph
   * @return edge_begin()
   * @post has_edge(@a *e_it.node1(), @a *e_it.node2()) == false
   * @post If old has_edge(@a *e_it.node1(), @a *e_it.node2()) false, new num_edges() == old num_edges().
   *       Else,                                                      new num_edges() == old num_edges() - 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: O(log(num_nodes()))
   */
  edge_iterator remove_edge(edge_iterator e_it) {
  	size_type rm = remove_edge(*e_it);
  	if (rm == 0) {
  		return edge_begin();
  	}
  	else {
  		if (e_it == edge_end()) {
  			return edge_begin();
  		}
  		else {
  			return edge_begin();
  		}
  	}
  }

  /** Remove an edge to the graph, return 0 if edge not in graph, 1 is edge successfully removed
   * @param[in] a One of the nodes of the edge
   * @param[in] b One of the nodes of the edge
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return size_type 0 if edge not in graph, size_type 1 if edge successfully removed
   * @post has_edge(@a a, @a b) == false
   * @post If old has_edge(@a a, @a b) false, new num_edges() == old num_edges().
   *       Else,                              new num_edges() == old num_edges() - 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: O(log(num_nodes()))
   */
  size_type remove_edge(const Node& a, const Node& b) {

  	// Check if edge is in graph
  	if (has_edge(a, b) == false) {
  		return 0;
  	}
  	else {
  		// Find unique identifier of edge
  		size_type uid = connections_.at(a.n_uid_).at(b.n_uid_);

  		// Remove edges from connections
  		connections_.at(a.n_uid_).erase(b.n_uid_);
  		connections_.at(b.n_uid_).erase(a.n_uid_);

        // Find internal edge
  		Internal_Edge int_e = edges_.at(uid);

  		// Keep track of this edges index
  		size_type tmp_index = int_e.e_idx_;

  		// Invalidate internal edge
  		edges_.at(uid).is_valid_ = false;

  		// Get back uid
  		size_type back_uid = e_i2u_.back();

  		// Swap with last edge in e_i2u
  		e_i2u_.at(tmp_index) = e_i2u_.back();

  		// Pop off last element
  		e_i2u_.pop_back();

  		// Change index of swapped edge
  		edges_.at(back_uid).e_idx_ = tmp_index;

  		return 1;
  	}
  }

  /** Remove an edge to the graph, return 0 if edge not in graph, 1 is edge successfully removed
   * @param[in] e An edge
   * @pre @a e.node1() and @a e.node2() are distinct valid nodes of this graph
   * @return size_type 0 if edge not in graph, size_type 1 if edge successfully removed
   * @post has_edge(@a e.node1(), @a e.node2()) == false
   * @post If old has_edge(@a e.node1(), @a e.node2()) false, new num_edges() == old num_edges().
   *       Else,                                              new num_edges() == old num_edges() - 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: O(log(num_nodes()))
   */
  size_type remove_edge(const Edge& e) {
  	return remove_edge(e.node1(), e.node2());
  }



  /** Remove a node to the graph, return 0 if node not in graph, 1 if node successfully removed
   * @param[in] n A node
   * @pre @a n is valid node of this graph
   * @return size_type 0 if node not in graph, size_type 1 if node successfully removed
   * @post has_node(@a n) == false
   * @post If old has_node(@a n) false, new num_nodes() == old num_nodes()
   *                                    new num_edges() == old num_edges().
   *       Else,                        new num_nodes() == old num_nodes() - 1                      
   *                                    new num_edges() == old num_edges() - old degree(@a n).
   *
   * Can invalidate node indexes -- in other words, old node(@a i) might not
   * equal new node(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: O(degree(@a n ))
   */
  size_type remove_node(const Node& n) {
  	if (has_node(n) == false) {
  		return 0;
  	}
  	else {

  		// Delete all edges from this node
  		for(auto it = n.edge_begin(); it != n.edge_end();) {
  			it = n.edge_begin();
  			Edge e = *it;
  			remove_edge(e);
  		}

  		size_type tmp_idx = n.index();

  		size_type tmp_uid = n.n_uid_;

  		nodes_.at(tmp_uid).is_valid_ = false;

  		size_type back_uid = n_i2u_.back();

  		n_i2u_.at(tmp_idx) = back_uid;

  		n_i2u_.pop_back();

  		nodes_.at(back_uid).n_idx_ = tmp_idx;

  		return 1;
  	}
  }


  /** Remove a node to the graph, return 0 if node not in graph, 1 if node successfully removed
   * @param[in] n_it An iterator to a node
   * @pre @a *n_it is valid node of this graph
   * @return size_type 0 if node not in graph, size_type 1 if node successfully removed
   * @post has_node(@a *n_it) == false
   * @post If old has_node(@a *n_it) false, new num_nodes() == old num_nodes()
   *                                        new num_edges() == old num_edges().
   *       Else,                            new num_nodes() == old num_nodes() - 1                      
   *                                        new num_edges() == old num_edges() - old degree(@a *n_it).
   *
   * Can invalidate node indexes -- in other words, old node(@a i) might not
   * equal new node(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: O(degree(@a *n_it ))
   */
  node_iterator remove_node(node_iterator n_it) {
  	return remove_node(*n_it);
  }


 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  class Internal_Node {

    friend class Graph;
    // Attribute to keep track of graph of node
    Graph* graph_;

    // Attribute to keep track of index of node
    size_type n_idx_;

    // Attribute to keep track of identifer of node
    size_type n_uid_;

    // Attribute to keep track of point of node
    Point position_;

    // Attribute to keep track of value of node
    node_value_type val_;

    // Attribute to keep track of the validity of internal node
    bool is_valid_ = true;

    // Default constructor
    Internal_Node() {
    }
    // Constructor for Internal_Node
    Internal_Node(const Graph* graph, size_type n_idx, size_type n_uid, Point position, node_value_type val)
      : graph_(const_cast<Graph*>(graph)), n_idx_(n_idx), n_uid_(n_uid), position_(position), val_(val) {
      }
  };

  class Internal_Edge {

    friend class Graph;
    // Attribute to keep track of graph of edge
    Graph* graph_;

    // Attribute to keep track of index of edge
    size_type e_idx_;

    // Attribute to keep track of identifier of edge
    size_type e_uid_;

    // Attribute to keep track of node identifiers of edge
    size_type n_uid1_;
    size_type n_uid2_;

    // Attribute to keep track of value of edge
    edge_value_type val_;

    // Attribute to keep track of the validity of internal edge
    bool is_valid_ = true;

    // Default constructor
    Internal_Edge() {
    }
    // Constructor for Internal_Node
    Internal_Edge(const Graph* graph, size_type e_idx, size_type e_uid, size_type n_uid1, size_type n_uid2, edge_value_type val)
      : graph_(const_cast<Graph*>(graph)), e_idx_(e_idx), e_uid_(e_uid), n_uid1_(n_uid1), n_uid2_(n_uid2), val_(val) {
      }
  };
  
  // Vector of internal nodes representing our Graphs nodes
  std::vector<Internal_Node> nodes_;

  // Unique Identifier for nodes
  size_type next_n_uid_ = 0;

  // Map of Maps representing an adjacency list
  std::map<size_type, std::map<size_type, size_type>> connections_;

  // Vector of tuples of node objects representing edges
  std::vector<Internal_Edge> edges_;

  // Unique identifier for edges
  size_type next_e_uid_ = 0;

  // Vector from node index to node unique identifer
  std::vector<size_type> n_i2u_;

  // Vector from edge index to edge unique identifier
  std::vector<size_type> e_i2u_;

};

#endif // CME212_GRAPH_HPP