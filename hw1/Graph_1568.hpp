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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : next_nid_(FIRST_INDEX), first_valid_nid_(FIRST_INDEX),
            next_eid_(FIRST_INDEX), first_valid_eid_(FIRST_INDEX) {
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
      
      // index() checks if node is valid
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
      
      // As per piazza post #35, assume there is no comparison between graphs
      // (i.e. assume this and n belong to the same graph)
      return (*this == n ? false : this->index() < n.index());
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
      return (*this == e ? false : this->index() < e.index());
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
      
      // Add edge
      edges_.push_back( Edge(this, a.nid_, b.nid_, next_eid_) );
      ++next_eid_;
      
      return edges_.back();
      
    }
    size_type eid = get_edge_eid(a, b); // We know this edge exists in the graph
    return edges_[edge_indices_[eid]];
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
  
  /** Reutrn an edge iterator to the beinning of the edges */
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
  
  /** Helper method that returns a pointer to the graph's copy of the node */
  Node* get_node_pointer(const Node& a) const {
    return const_cast<Node*>(&nodes_[a.index()]);
  }

};

#endif // CME212_GRAPH_HPP
