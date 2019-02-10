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

  /** Node value template type */
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
      // Return Point position using index and Graph pointer
      return ((this->graph_)->listOfNodes[this->index_].position_);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // Return this Node's index
      return this->index_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Node's value.
     * @return reference to Node's value
     *
     * Complexity: O(1) amortized operations.
     */
    node_value_type& value() {
      // Return Node value using index and Graph pointer
      return ((this->graph_)->listOfNodes[this->index_].node_value_);
    }

    /** Node's value.
     * @return constant reference to Node's value
     *
     * Complexity: O(1) amortized operations.
     */
    const node_value_type& value() const{
      // Return Node value using index and Graph pointer
      return ((this->graph_)->listOfNodes[this->index_].node_value_);
    }

    /** Number of Edges incident on a Node
     * @return number of Edges incident upon this Node
     *
     * Complexity: O(1) amortized operations.
     */
    size_type degree() const {
      return (((this->graph_)->listOfNodes[this->index_]).edge_indices_).size();
    }

    /** IncidentIterator pointing to the first Edge incident to Node
     * @return incident_iterator pointing to first Edge incident to this Node
     *
     * Complexity: O(1) amortized operations.
     */
    incident_iterator edge_begin() const {
      incident_iterator ni = IncidentIterator();
      ni.graph_ = this->graph_;
      ni.node_index_ = this->index_;
      ni.iter_index_ = 0;
      return ni;
    }

    /** Return IncidentIterator pointing to the "past-the-end" Edge incident to Node
     * @return incident_iterator pointing to the non-existent "past-the-end" Edge incident to this Node
     *
     * Complexity: O(1) amortized operations.
     */
    incident_iterator edge_end() const{
      incident_iterator ni = IncidentIterator();
      ni.graph_ = this->graph_;
      ni.node_index_ = this->index_;
      ni.iter_index_ = (((ni.graph_)->listOfNodes[ni.node_index_]).edge_indices_).size();
      return ni;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // Return whether Nodes have the same Graph and the same index.
      bool sameGraph = ( this->graph_ == n.graph_ ); // Check if Nodes have same Graph pointer
      bool sameIndex = ( this->index() == n.index() ); // Check if same index
      return (sameIndex && sameGraph);
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
      // Define order based on Node's index
      return ( this->index() < n.index() );
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Private Variables tracked by Node
    graph_type* graph_; // Pointer back to Graph containing Node
    size_type index_; // Index value associated with Node
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // Return number of Nodes
    return size_type(this->listOfNodes.size());
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
  Node add_node(const Point& position, const node_value_type& node_value = node_value_type()) {
    // Create a Node with pointer to Graph and index
    Node n = Node(); // Create (invalid) Node inside nodeInfo struct
    n.graph_ = this; // Assign this Graph to Node
    n.index_ = size_type(this->listOfNodes.size()); // Give Index to Node
    // Create a nodeInfo struct with Node and associated information
    nodeInfo n_info; // Create Empty nodeInfo struct
    n_info.position_ = position; // Give position to nodeInfo struct
    n_info.index_ = size_type(this->listOfNodes.size()); // Give Index to nodeInfo struct
    n_info.node_value_ = node_value; // Give Node Value to nodeInfo struct
    n_info.node_ = n; // Put Node inside nodeInfo struct
    // Add Node (with correct index) to Graph
    this->listOfNodes.push_back(n_info); // Assign index to Node
    return n; // Return Node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // Make sure that Node is inside Graph at index indicated by Node
    if (n.index_ < this->num_nodes()) { // Make sure index of Node is inside Graph
      return (this->listOfNodes[n.index_].node_ == n); // Match Nodes
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
    // Return Node with index i
    return this->listOfNodes[i].node_;
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
      // Return reference to dereferenced pointer
      return ((this->graph_)->node(this->node1_index_));
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // Return reference to dereferenced pointer
      return ((this->graph_)->node(this->node2_index_));
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // Equal Edges represent the same undirected edge between two Nodes
      bool nodesMatchExactly = ((e.node1() == this->node1()) && (e.node2() == this->node2()));
      bool nodesMatchSwitched = ((e.node1() == this->node2()) && (e.node2() == this->node1()));
      return (nodesMatchExactly || nodesMatchSwitched);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // Inequality based on Edge index
      return (this->index_ < e.index_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    graph_type* graph_; // Pointer back to Graph containing Edge
    size_type index_; // Index value associated with Edge
    size_type node1_index_; // Index of Node 1
    size_type node2_index_; // Index of Node 2
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // Return number of Edges in Graph
    return size_type(this->listOfEdges.size());
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // Return Edge from Graph
    return this->listOfEdges[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Find out how many edges are incident on each node
    size_type num_edges_a = ((this->listOfNodes[a.index_]).edge_indices_).size();
    size_type num_edges_b = ((this->listOfNodes[b.index_]).edge_indices_).size();
    // Pick Node with the fewest incident Edges
    if (num_edges_a < num_edges_b) {
      // Iterate through the Edges associated with this Node
      for (incident_iterator it = a.edge_begin(); it != a.edge_end(); ++it) {
        if ( (((*it).node1()==a) && ((*it).node2()==b)) || (((*it).node1()==b) && ((*it).node2()==a)) ) {
          return true; // Return true if this Edge connects Nodes
        }
      }
      return false; // Return false if no Edge connects those two Nodes
    } else {
      // Iterate through the Edges associated with this Node
      for (incident_iterator it = b.edge_begin(); it != b.edge_end(); ++it) {
        if ( (((*it).node1()==a) && ((*it).node2()==b)) || (((*it).node1()==b) && ((*it).node2()==a)) ) {
          return true; // Return true if this Edge connects Nodes
        }
      }
      return false; // Return false if no Edge connects those two Nodes
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
    // Create Edge with two associated Nodes, and add to Graph
    if (this->has_edge(a,b)) { // If Graph has Edge that connects these Nodes
      // Find out how many edges are incident on each node
      size_type num_edges_a = ((this->listOfNodes[a.index_]).edge_indices_).size();
      size_type num_edges_b = ((this->listOfNodes[b.index_]).edge_indices_).size();
      // Pick Node with the fewest incident Edges
      if (num_edges_a < num_edges_b) {
        // Iterate through the Edges associated with this Node
        for (incident_iterator it = a.edge_begin(); it != a.edge_end(); ++it) {
          if ( (((*it).node1()==a) && ((*it).node2()==b)) || (((*it).node1()==b) && ((*it).node2()==a)) ) {
            return (*it); // Return if this Edge connects Nodes
          }
        }
      } else {
        // Iterate through the Edges associated with this Node
        for (incident_iterator it = b.edge_begin(); it != b.edge_end(); ++it) {
          if ( (((*it).node1()==a) && ((*it).node2()==b)) || (((*it).node1()==b) && ((*it).node2()==a)) ) {
            return (*it); // Return if this Edge connects Nodes
          }
        }
      }
      assert(false); // Execution should not reach this point! If it does something is wrong
      return Edge(); // Gets rid of warning
    } else {
      // Create Edge
      edge_type e = Edge(); // Create (invalid) Edge
      e.graph_ = this; // Assign this Graph to Edge
      e.index_ = size_type(this->listOfEdges.size()); // Record index of Edge
      e.node1_index_ = a.index(); // First Node associated with Edge
      e.node2_index_ = b.index(); // Second Node associated with Edge
      (this->listOfEdges).push_back(e); // Add Edge to Graph
      // Add to the list of Edge indices for each Node
      (this->listOfNodes[a.index_].edge_indices_).push_back(e.index_);
      (this->listOfNodes[b.index_].edge_indices_).push_back(e.index_);
      return e; // Return edge
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // Remove all Nodes
    this->listOfNodes.clear();
    // Remove all Edges
    this->listOfEdges.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> { // REHMAN ADDED: " : private totally_ordered<NodeIterator> ""
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

    /** Define dereferencing operator for NodeIterator
     * @return Node that node_iterator "points" to
     *
     * Complexity: O(1) amortized operations.
     */
    Node operator*() const {
      // Return Node using index and Graph pointer
      return ((this->graph_)->node(this->index_));
    }

    /** Define increment operator for NodeIterator
     * @pre ((*this)->graph_)->node_begin() <= *(*this) < ((*this)->graph_)->node_end();
     * @return node_iterator that "points" to next Node on graph
     * @post ((*this)->graph_)->node_begin() < *(*this) <= ((*this)->graph_)->node_end();
     *
     * Complexity: O(1) amortized operations.
     */
    NodeIterator& operator++() {
      // Increment NodeIterator's index and return NodeIterator
      (this->index_)++;
      return (*this);
    }

    /** Define equality operator for two NodeIterators
     * @param ni Constant reference to node_iterator
     * @return boolean for whether "this" node_iterator and @a ni have the same graph and index
     *
     * Complexity: O(1) amortized operations.
     */
    bool operator==(const NodeIterator& ni) const {
      // Return whether NodeIterators have the same Graph and the same index.
      bool sameGraph = ( this->graph_ == ni.graph_ ); // Check if NodeIterators have same Graph pointer
      bool sameIndex = ( this->index_ == ni.index_ ); // Check if same index
      return (sameIndex && sameGraph);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE

    graph_type* graph_; // Pointer back to Graph containing Edge
    size_type index_; // Index value associated with Edge
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** NodeIterator pointing to the first Node in Graph
   * @return node_iterator pointing to the first Node in Graph
   *
   * Complexity: O(1) amortized operations.
   */
  node_iterator node_begin() const {
    node_iterator ni = NodeIterator();
    ni.graph_ = const_cast<graph_type*>(this); // "const_cast" removes constness from "this"
    ni.index_ = 0;
    return ni;
  }

  /** Return NodeIterator pointing to the "past-the-end" Node in Graph
   * @return node_iterator pointing to the nonexistent "past-the-end" Node in Graph
   *
   * Complexity: O(1) amortized operations.
   */
  node_iterator node_end() const {
    node_iterator ni = NodeIterator();
    ni.graph_ = const_cast<graph_type*>(this); // "const_cast" removes constness from "this"
    ni.index_ = this->size();
    return ni;
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> { // REHMAN ADDED: " : private totally_ordered<IncidentIterator> "
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

    /** Define dereferencing operator for IncidentIterator
     * @return Edge that this incident_iterator points to
     * @post node1() for returned Edge must be the Node that incident_iterator belongs to
     *
     * Complexity: O(1) amortized operations.
     */
    Edge operator*() const {
      // Return Node using Node index, Iterator Index and Graph pointer
      edge_type e = (this->graph_)->edge(((this->graph_)->listOfNodes[this->node_index_]).edge_indices_[this->iter_index_]);
      // Make sure that e.node1() is this Node
      if (e.node1() != (this->graph_)->node(this->node_index_)) {
        size_type n1_index = (e.node1()).index();
        size_type n2_index = (e.node2()).index();
        e.node1_index_ = n2_index;
        e.node2_index_ = n1_index;
        assert(e.node1() == (this->graph_)->node(this->node_index_));
      }
      return e;
    }

    /** Define increment operator for IncidentIterator
     * @return incident_iterator pointing to the next Edge incident to Node
     *
     * Complexity: O(1) amortized operations.
     */
    IncidentIterator& operator++() {
      // Increment IncidentIterator's index and return IncidentIterator
      (this->iter_index_)++;
      return (*this);
    }

    /** Define equality operator for two IncidentIterators
     * @param ii Constant reference to incident_iterator
     * @return Boolean that this incident_iterator and @a ii refer to the
     *    same Graph, Node Index, and Index into Incident Edges.
     *
     * Complexity: O(1) amortized operations.
     */
    bool operator==(const IncidentIterator& ii) const {
      // Return whether NodeIterators have the same Graph and the same index.
      bool sameGraph = ( this->graph_ == ii.graph_ ); // Check if NodeIterators have same Graph pointer
      bool sameNodeIndex = ( this->node_index_ == ii.node_index_ ); // Check if same Node index
      bool sameIteratorIndex = ( this->iter_index_ == ii.iter_index_ ); // Check if same Iterator index
      return (sameNodeIndex && sameGraph && sameIteratorIndex);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE

    graph_type* graph_; // Pointer back to Graph containing Node
    size_type node_index_; // Index value associated with Node
    size_type iter_index_; // Index value associated with Node Iterator
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> { // REHMAN ADDED: " : private totally_ordered<EdgeIterator> "
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

    /** Define dereferencing operator for EdgeIterator
     * @return Edge that this edge_iterator points to
     *
     * Complexity: O(1) amortized operations.
     */
    Edge operator*() const {
      // Return Edge using index and Graph pointer
      return ((this->graph_)->edge(this->index_));
    }

    /** Define increment operator for EdgeIterator
     * @return edge_iterator that points to next Edge in Graph.
     *
     * Complexity: O(1) amortized operations.
     */
    EdgeIterator& operator++() {
      // Increment EdgeIterator's index and return EdgeIterator
      (this->index_)++;
      return (*this);
    }

    /** Define equality operator for two EdgeIterators
     * @param ei Constant reference to edge_iterator
     * @return Boolean that "this" edge_iterator and ei have the same Graph and index
     *
     * Complexity: O(1) amortized operations.
     */
    bool operator==(const EdgeIterator& ei) const {
      // Return whether EdgeIterators have the same Graph and the same index.
      bool sameGraph = ( this->graph_ == ei.graph_ ); // Check if EdgeIterators have same Graph pointer
      bool sameIndex = ( this->index_ == ei.index_ ); // Check if same index
      return (sameIndex && sameGraph);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE

    graph_type* graph_; // Pointer back to Graph containing Edge
    size_type index_; // Index value associated with Edge
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** EdgeIterator pointing to the first Edge in Graph
   * @return edge_iterator pointing to the first Edge in Graph
   *
   * Complexity: O(1) amortized operations.
   */
  edge_iterator edge_begin() const {
    edge_iterator ei = EdgeIterator();
    ei.graph_ = const_cast<graph_type*>(this); // "const_cast" removes constness from "this"
    ei.index_ = 0;
    return ei;
  }

  /** Return EdgeIterator pointing to the "past-the-end" Edge in Graph
   * @return edge_iterator pointing to the nonexistent "past-the-end" Edge in Graph
   *
   * Complexity: O(1) amortized operations.
   */
  edge_iterator edge_end() const {
    edge_iterator ei = EdgeIterator();
    ei.graph_ = const_cast<graph_type*>(this); // "const_cast" removes constness from "this"
    ei.index_ = this->num_edges();
    return ei;
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // Structs for Holding Node and Edge information
  struct nodeInfo {
    Point position_; // position as a Point object
    size_type index_; // Index value associated with Node
    node_value_type node_value_; // Node Value Associated with Node
    node_type node_; // Node
    std::vector<size_type> edge_indices_; // Vector of Incident Edge Indices
  };

  // Node information within Graph
  std::vector<nodeInfo> listOfNodes; // Vector of nodeInfo structs

  // Edge information within graph
  std::vector<edge_type> listOfEdges; // Vector of Edges

};

#endif // CME212_GRAPH_HPP
