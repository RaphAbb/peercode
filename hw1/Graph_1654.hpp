#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
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
template <typename V>
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
  using graph_type = Graph<V>;

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
  using node_value_type = V;

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
      return graph_->nodes.at(this-> idx_).position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return this-> idx_;
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
      return graph_->nodes.at(this-> idx_).val_;
    }

    /** Retrieve a node's value
    *
    * Complexity: O(1)
    */
    const node_value_type& value() const {
      return graph_->nodes.at(this-> idx_).val_;
    }

    /** Return node's degree i.e., the number of incident nodes
    * 
    * Complexity: O(1)
    */
    size_type degree() const {
      return graph_->connections.at(this-> idx_).size();
    }

    /** Return incident iterator begin instance 
    *
    * Complexity: O(1)
    */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_->connections.at(this-> idx_).begin(), this->idx_, this->graph_);

    }

    /** Return incident iterator end instance
    *
    * Complexity: O(1) 
    */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_->connections.at(this-> idx_).end(), this->idx_, this->graph_);
    }


    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if ((n.graph_ == this->graph_) and (n.idx_ == this->idx_)) {
      	return true;
      }

      (void) n;          // Quiet compiler warning
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
      if (this-> idx_ < n.idx_) {
      	return true;
      }
      (void) n;           // Quiet compiler warning
      return false;
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
    size_type idx_;

    // Constructor for Node
    Node(graph_type* g, size_type id)
      : graph_(g), idx_(id) {
      }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return this->next_idx;
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

    // Add internal node to the Graph
    this->nodes.push_back( Internal_Node(this, this->next_idx, position, node_value_type()) );

    // Increment next_idx
    this->next_idx++;

    // Add node and empty map to adjacency map
    this->connections[next_idx-1] = std::map<size_type, size_type>();

    // Create and return node instance
    return Node(this, next_idx - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return this == n.graph_;
    (void) n;            // Quiet compiler warning
    return false;
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
      return Node(const_cast<graph_type*>(this), i);
    }

    // If precondition fails, return invalid node
    (void) i;             // Quiet compiler warning
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
      return Node(graph_, this->node1_idx_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, this->node2_idx_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if ((e.graph_ == this->graph_) and (e.idy_ == this->idy_)) {
        return true;
      }

      (void) e;           // Quiet compiler warning
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (this-> idy_ < e.idy_) {
        return true;
      }
      (void) e;           // Quiet compiler warning
      return false;
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
    size_type idy_;

    // Attribute to keep track of nodes of edge
    size_type node1_idx_;
    size_type node2_idx_;

    // Constructor for Edge
    Edge(graph_type* g, size_type id, size_type node1_idx, size_type node2_idx)
      : graph_(g), idy_(id), node1_idx_(node1_idx), node2_idx_(node2_idx) {
      }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE

    // Check for precondition
    // Note: no need to check 0 <= i, as 0 <= size_type always true
    if (i < num_edges()) {

      return Edge(const_cast<graph_type*>(this), i, edges.at(i).n_id1_, edges.at(i).n_id2_);
    }
    
    // If precondition fails, return an invalid edge
    (void) i;             // Quiet compiler warning
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
      size_type idx1 = a.index();
      size_type idx2 = b.index();

      if (connections.at(idx1).find(idx2) == connections.at(idx1).end()) {
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
    (void) a; (void) b;   // Quiet compiler warning
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

    //Test for precondition 
    if (a == b) {
      throw std::runtime_error("Self loop illegal");
    }
    if (!has_node(a) or !has_node(b)) {
      throw std::runtime_error("Invalid nodes for edge addition");
    }

    // If edge present in graph
    if (has_edge(a,b)) {
      size_type id = connections.at(a.index()).at(b.index());
      return Edge(this, id, a.index(), b.index());
    }

    // Otherwise, add edge to graph, return edge object
    this->edges.push_back( Internal_Edge(this, next_idy, a.index(), b.index()) );

    edge_type new_edge =Edge(this, next_idy, a.index(), b.index());

    // Add edge to our adjacency map
    connections[a.index()][b.index()] = next_idy;
    connections[b.index()][a.index()] = next_idy;

    this->next_idy++;

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
    nodes.clear();

    // Reset next_idx (num_nodes) to 0
    next_idx = 0;

    // Delete all edges stored in adjacency list
    connections.clear();

    // Delete all edges stored in vector
    edges.clear();

    // Reset next_idy (num_edges) to 0
    next_idy = 0;
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
      Internal_Node internal_node = *it_;
      return Node(internal_node.graph_, internal_node.idx_);
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
    * Equal node iterators have the same iterator attribute
    * Complexity: O(1)
    */
    bool operator==(const NodeIterator& n) const {
      return it_ == n.it_;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE

    typename std::vector<Internal_Node>::const_iterator it_;
    NodeIterator(typename std::vector<Internal_Node>::const_iterator it)
      : it_(it) {
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
    node_iterator b = NodeIterator(nodes.begin());
    return b;
  }

  /** Return node iterator end instance
  *
  * Complexity: O(1)
  */
  node_iterator node_end() const {
    node_iterator e = NodeIterator(nodes.end());
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
      size_type node2_id = it_->first;
      size_type edge_id = it_->second;
      return Edge(this->graph_, edge_id, this->n_id_, node2_id);
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
      if ( (x.it_ == this->it_) && (x.n_id_ == this->n_id_) && (x.graph_ == this->graph_) ){
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
    size_type n_id_;

    // Attribute to keep track of graph of node
    graph_type* graph_;

    IncidentIterator(typename std::map<size_type, size_type>::const_iterator it, size_type n_id, graph_type* g)
      : it_(it), n_id_(n_id), graph_(g) {
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
      Internal_Edge e = *it_;
      return Edge(e.graph_, e.idx_, e.n_id1_, e.n_id2_);
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
      return it_ == e.it_;
    }


   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE

    typename std::vector<Internal_Edge>::const_iterator it_;
    EdgeIterator(typename std::vector<Internal_Edge>::const_iterator it)
      : it_(it) {

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
    edge_iterator b = EdgeIterator(edges.begin());
    return b;
  }

  /** Returns edge iterator end instance
  *
  * Complexity: O(1)
  */
  edge_iterator edge_end() const {
    edge_iterator e = EdgeIterator(edges.end());
    return e;
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
    size_type idx_;
    // Attribute to keep track of point of node
    Point position_;
    // Attribute to keep track of value of node
    node_value_type val_;

    // Default constructor
    Internal_Node() {
    }
    // Constructor for Internal_Node
    Internal_Node(const Graph* graph, size_type id, Point position, node_value_type val)
      : graph_(const_cast<Graph*>(graph)), idx_(id), position_(position), val_(val) {
      }
  };

  class Internal_Edge {

    friend class Graph;
    // Attribute to keep track of graph of edge
    Graph* graph_;
    // Attribute to keep track of index of edge
    size_type idx_;
    // Attribute to keep track of node indicies of edge
    size_type n_id1_;
    size_type n_id2_;

    // Default constructor
    Internal_Edge() {
    }
    // Constructor for Internal_Node
    Internal_Edge(const Graph* graph, size_type id, size_type n_id1, size_type n_id2)
      : graph_(const_cast<Graph*>(graph)), idx_(id), n_id1_(n_id1), n_id2_(n_id2) {
      }
  };
  
  // Vector of internal nodes representing our Graphs nodes
 	std::vector<Internal_Node> nodes;

  // Index for points
  size_type next_idx = 0;

  // Map of Maps representing an adjacency list
  std::map<size_type, std::map<size_type, size_type>> connections;

  // Vector of tuples of node objects representing edges
  std::vector<Internal_Edge> edges;

  // Index for edges
  size_type next_idy = 0;

};

#endif // CME212_GRAPH_HPP