#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <tuple>
#include <vector>
#include <unordered_map>
#include <map>
#include <stdexcept>
#include <cassert>
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V> // Defines the node value type
class Graph {

// Publicly declare templated node value type 
public:
  using node_value_type = V;

 
// Declare private internal structures for their uses in nodevec_ and edgevec_ below
private:
  class InternalNode;
  class InternalEdge;

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  // Store the nodes using a vector of InternalNode objects. The i^th element in the vector 
  // contains an InternalNode object whose data attributes describe information about the Node
  // with nodeID == i.
  std::vector<InternalNode> nodevec_;

  // Store the edges with an adjacency list. adjlist_ maps nodeIDs to maps of
  // nodeIDs to edgeIDs. Thus if the edge with edgeID 4 connects the nodes with
  // nodeIDs 1 and 2, adjlist_ will contain the entries 1 -> (2 -> 4) and 2 -> (1 -> 4).
  // This allows O(1) lookup for both the presence AND the ID of an edge given two nodes, the
  // first of which is useful in has_edge() and the second of which is useful in add_edge().
  std::unordered_map<unsigned int, std::map<unsigned int,unsigned int>> adjlist_;

  // Store the edges using a vector of InternalEdge objects. The i^th element in the vector 
  // contains an InternalEdge object whose data attributes describe information about the Edge
  // with edgeID == i.
  std::vector<InternalEdge> edgevec_;

public:

 //
 // PUBLIC TYPE DEFINITIONS
 //

 /** Type of this graph. */
 using graph_type = Graph;

/** Predeclaration of Node and Edge types. */
 class Node;
 class Edge;

 /** Synonym for Node (following STL conventions). */
 using node_type = Node;

 /** Synonym for Edge (following STL conventions). */
 using edge_type = Edge;

 using internal_node_type = InternalNode;

 using internal_edge_type = InternalEdge;

 class NodeIterator;
 using node_iterator = NodeIterator;

 class EdgeIterator;
 using edge_iterator = EdgeIterator;

 class IncidentIterator;
 using incident_iterator = IncidentIterator;

 /** Type of indexes and sizes.
     Return type of Graph::Node::index(), Graph::num_nodes(),
     Graph::num_edges(), and argument type of Graph::node(size_type) */
 using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph():
  // HW0: YOUR CODE HERE
  nodevec_(), adjlist_(), edgevec_()
  {}

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
      return this->graph_->nodevec_[nodeid_].position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return this->nodeid_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

  /** Return this node's value.
   * @return an object of the templated node_value_type representing this node's value.
   * @pre *(this->graph_).node(this->nodeid_) == *this
   *
   * Complexity: O(1)
   */
    node_value_type& value(){
    return this->graph_->nodevec_[nodeid_].value_;
    }

  /** Return this node's value as an object which cannot be modified.
   * @return a const object of the templated node_value_type representing this node's value.
   * @pre *(this->graph_).node(this->nodeid_) == *this
   *
   * Complexity: O(1)
   */
    const node_value_type& value() const {
    return this->graph_->nodevec_[nodeid_].value_;
    }

  /** Return this node's degree.
   * @return a size_type object representing the number of edges incident to this node.
   * @pre *(this->graph_).node(this->nodeid_) == *this
   *
   * Complexity: O(1)
   */
    size_type degree() const {
      return this->graph_->adjlist_.at(this->nodeid_).size();
    }

  /** Creates an iterator to the first (ordered by increasing NodeID of the incident Node) Edge incident to this Node in this Graph.
   * @return  an incident_iterator which (for a non-isolated Node) dereferences to the Edge with minimal value of Edge.node2().index()
   * 
   * Intended to be functionally used exactly as the .begin() method of an STL vector.
   * 
   * Complexity: O(1)
   */
    incident_iterator edge_begin() const {
      return IncidentIterator((this->graph_->adjlist_.at(this->nodeid_)).begin(), *this);
    }

  /** Creates an iterator to one past the last (ordered by increasing NodeID of the incident Node) Edge incident to this Node in this Graph.
   * @return  an incident_iterator intended to signify that an iteration process over incident edges should terminate.
   * 
   * Intended to be functionally used exactly as the .end() method of an STL vector.
   * 
   * Complexity: O(1)
   */
    incident_iterator edge_end() const {
      return IncidentIterator(this->graph_->adjlist_.at(this->nodeid_).end(), *this);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (n.graph_ == this->graph_ && n.index() == this->index()){return true;}
      else{return false;}
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
      // Base the decision on node IDs.
      // As stated on Piazza, we can assume the nodes are in the same Graph
        if (this->nodeid_ >= n.index()){return false;}
        else {return true;}
      }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    Graph* graph_; // Pointer back to graph
    size_type nodeid_; // Unique node id
    // Private constructor
    Node(const Graph* graph, size_type id):
    graph_(const_cast<Graph*>(graph)), nodeid_(id) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodevec_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }


  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value    The new node's value (defaults to the default object of type node_value_type)
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    unsigned int numnodes_ = nodevec_.size();
    Node new_node(this, numnodes_); // Create new node
    nodevec_.push_back(InternalNode(position, value, numnodes_, this)); // Add the provided info of new node to the nodevec
    this->adjlist_[numnodes_]; // Add the node to the adjacency list
    return new_node;        // Invalid node
  }

/** Update the value of the Node with provided NodeID
   * @param[in] indx The NodeID of the Node to be updated
   * @tparam[in] value  the value (of templated type node_value_type) to be updated
   * @return none
   * @pre 0 <= indx < num_nodes()
   * 
   *
   * Complexity: O(1)
   */
  void update_value(size_type indx, const node_value_type& value){
    if (indx < nodevec_.size()){
    this->nodevec_[indx].newval(value);}
    return;
  }


  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (n.graph_ == this && nodevec_.size() > n.index()){return true;}
    else{return false;}
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i >= nodevec_.size()){
      // Out of bounds; throw error
    throw std::invalid_argument("Index out of range");
    }
    else{
    Node nodei = Node(this, i);
    return nodei;}       // Invalid node
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
      return Node(this->graph_, this->node1id_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(this->graph_, this->node2id_);
    }

    /** Return this edge's index, a number in the range [0, numedges_). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return this->edgeid_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */

    bool operator==(const Edge& e) const {
      if (e.graph_ == this->graph_){ // Verify the Graphs are the same
        // Verify that the pairs of nodes are identical
        if ((e.node1() == this->node1()) && (e.node2() == this->node2())){
          return true;
        }
        else if ((e.node1() == this->node2()) && (e.node2() == this->node1())){
          return true;
        }
        // Edge is between different pairs of nodes
        else{return false;}
      }
      // Edges are in different graphs
      else{return false;}
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */

    bool operator<(const Edge& e) const {
      // Again, can assume the edges are in the same graph, according to Piazza
      // Thus, we can just base this on edgeID.
      if (this->edgeid_ >= e.index()){return false;}
      else {return true;}
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_; // Pointer back to graph
    size_type edgeid_; // Unique edge id, the position in the edgevec_ vector
    size_type node1id_; // ID of one of the endpoints
    size_type node2id_; // ID of the other endpoint
    // Private constructor
    Edge(const Graph* graph, size_type id, size_type node1id, size_type node2id):
    graph_(const_cast<Graph*>(graph)), edgeid_(id), node1id_(node1id), node2id_(node2id) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return this->edgevec_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i >= edgevec_.size()){
      // Out of bounds; throw error
    throw std::invalid_argument("Index out of range");
    }
    else{
      // Return proxy. edgevec_ allows for O(1) lookup of the nodeIDs for this edgeID.
      size_type node1id = this->edgevec_[i].node1id_;
      size_type node2id = this->edgevec_[i].node2id_;
      Edge edgei = Edge(this, i, node1id, node2id);
      return edgei; }       // Invalid node
  }


  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const{
    // HW0: YOUR CODE HERE
    unsigned int ind1 = a.index();
    unsigned int ind2 = b.index();
    if (adjlist_.at(ind1).count(ind2) > 0){
        return true; // b is a neighbor of a
      }
      else{return false;} // b is not a neighbor of a
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
    if (has_edge(a,b)){
    // Return a proxy to the existing edge. First we need its index
    unsigned int indx = adjlist_[a.index()][b.index()];
    return Edge(this, indx, a.index(), b.index());
    }
    else{
    unsigned int numedges_ = edgevec_.size();
    // Now, the case where the edge was not already present in Graph
    Edge new_edge(this, numedges_, a.index(), b.index()); // Create new edge
    // Add edge information to the end of edgevec_
    edgevec_.push_back(InternalEdge(numedges_,a.index(),b.index(),this)); 
    // Update adjlist_
    adjlist_[a.index()][b.index()] = numedges_; // numedges_ is the new edgeID
    adjlist_[b.index()][a.index()] = numedges_; // numedges_ is the new edgeID
    return new_edge;
  }
}

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    this->nodevec_.clear();
    this->edgevec_.clear();
    this->adjlist_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : totally_ordered<NodeIterator> {
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

 /** Dereference the iterator to access the current node in an iteration.
   * @return the Node object at the current position of iteration
   * @pre this->curr_iter_ != nullptr
   *
   * Complexity: O(1)
   */
    Node operator*() const {
      InternalNode intnode = *curr_iter_;
      return Node(intnode.graph_, intnode.id_);
    }

  /** Increment the iterator to move to the next item in an iteration.
   * @return a reference to a NodeIterator pointing to the next position in an iteration
   * @pre this->curr_iter_ != nullptr
   *
   * Complexity: O(1)
   */
    NodeIterator& operator++() {
      curr_iter_++;
      return *this;
    }

  /** Tests equality with another NodeIterator object
   * @param[in] x      another NodeIterator
   * @return  True if x is identically equal to this NodeIterator
   * 
   * Two NodeIterators are equal if they point to the same place in memory.
   * 
   * Complexity: O(1)
   */
    bool operator==(const NodeIterator& x) const {
      return (this->curr_iter_ == x.curr_iter_);
    }


   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE

    // Define iterator attribute
    typename std::vector<InternalNode>::const_iterator curr_iter_;

    // Define iterator constructor
    NodeIterator(typename std::vector<InternalNode>::const_iterator it): curr_iter_(it){}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

/** Creates an iterator to the first (as ordered by NodeID) Node in this Graph.
   * @return  a node_iterator which (in a nonempty Graph) dereferences to the Node with NodeID == 0.
   * 
   * Intended to be functionally used exactly as the .begin() method of an STL vector.
   * 
   * Complexity: O(1)
   */
  node_iterator node_begin() const {
    return NodeIterator(this->nodevec_.begin());
  }

/** Creates an iterator to the last (as ordered by NodeID) Node in this Graph.
   * @return  a NodeIterator intended to signify that an iteration process over the Nodes in this Graph should terminate.
   * 
   * Intended to be functionally used exactly as the .end() method of an STL vector.
   * 
   * Complexity: O(1)
   */
  node_iterator node_end() const {
    return NodeIterator(this->nodevec_.end());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : totally_ordered<IncidentIterator> {
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

  /** Dereference the iterator to access the current Edge in an incidence iteration.
   * @return the Edge object at the current position of iteration
   * @pre this->curr_iter_ != nullptr
   *
   * Complexity: O(1)
   */
    Edge operator*() const {
      unsigned int adj_node_idx = curr_iter_->first;
      unsigned int edge_idx = curr_iter_->second;
      return Edge(this->graph_, edge_idx, this->node_idx_, adj_node_idx);
    }

  /** Increment the iterator to move to the next item in an iteration.
   * @return a reference to an IncidentIterator pointing to the next position in an iteration
   * @pre this->curr_iter_ != nullptr
   *
   * Complexity: O(1)
   */
    IncidentIterator& operator++() {
      curr_iter_++;
      return *this;
    }

  /** Tests equality with another IncidentIterator object
   * @param[in] x      another IncidentIterator
   * @return  True if x is identically equal to this IncidentIterator
   * 
   * Two IncidentIterators are equal if they point to the same place in memory.
   * 
   * Complexity: O(1)
   */
    bool operator==(const IncidentIterator& x) const {
      if((this->curr_iter_ == x.curr_iter_) && (this->graph_ == x.graph_) && (this->node_idx_ == x.node_idx_)){
        return true;}
        else{return false;}
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE

    // Define iterator attribute
    typename std::map<unsigned int, unsigned int>::const_iterator curr_iter_;

    // Store index of node over whose neighbors we are iterating
    const size_type node_idx_;

    // Store pointer back to graph for purpose of returning valid Edge objects
    const Graph* graph_;

    // Define iterator constructor
    IncidentIterator(typename std::map<unsigned int, unsigned int>::const_iterator it, const Node& node): 
    curr_iter_(it), node_idx_(node.nodeid_), graph_(const_cast<Graph*>(node.graph_)) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: totally_ordered<EdgeIterator> {
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

  /** Dereference the iterator to access the current edge in an iteration.
   * @return the Edge object at the current position of iteration
   * @pre this->curr_iter_ != nullptr
   *
   * Complexity: O(1)
   */
    Edge operator*() const {
      InternalEdge intedge = *curr_iter_;
      return Edge(intedge.graph_, intedge.edgeid_, intedge.node1id_, intedge.node2id_);
    }

  /** Increment the iterator to move to the next item in an iteration.
   * @return a reference to an EdgeIterator pointing to the next position in an iteration
   * @pre this->curr_iter_ != nullptr
   *
   * Complexity: O(1)
   */
    EdgeIterator& operator++() {
      curr_iter_++;
      return *this;
    }

  /** Tests equality with another EdgeIterator object
   * @param[in] x      another EdgeIterator
   * @return  True if x is identically equal to this EdgeIterator
   * 
   * Two EdgeIterators are equal if they point to the same place in memory.
   * 
   * Complexity: O(1)
   */
    bool operator==(const EdgeIterator& x) const {
      return (this->curr_iter_ == x.curr_iter_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    // Define iterator attribute
    typename std::vector<InternalEdge>::const_iterator curr_iter_;

    // Define iterator constructor
    EdgeIterator(typename std::vector<InternalEdge>::const_iterator it): curr_iter_(it){}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

/** Creates an iterator to the first (as ordered by EdgeID) Edge in this Graph.
   * @return  an edge_iterator which (in an edge-nonempty Graph) dereferences to the Edge with EdgeID == 0.
   * 
   * Intended to be functionally used exactly as the .begin() method of an STL vector.
   * 
   * Complexity: O(1)
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this->edgevec_.begin());
  }

/** Creates an iterator to the last (as ordered by EdgeID) Edge in this Graph.
   * @return  an EdgeIterator intended to signify that an iteration process over the Edges in this Graph should terminate.
   * 
   * Intended to be functionally used exactly as the .end() method of an STL vector.
   * 
   * Complexity: O(1)
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this->edgevec_.end());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // Define an internal class which stores various information about an Node. Assuming that
  // a user has "heavy"/memory-intensive Node information to be stored, this will facilitate
  // an overall lightweight Graph implementation using the proxy design pattern.
  class InternalNode {
    Point position_; // Node position
    node_value_type value_; // templated value
    size_type id_; // unique NodeID
    const Graph* graph_;
    friend class Graph;
    // Default private constructor
    InternalNode(){}
    // Initialize InternalNode with all information about Node.
    InternalNode(Point position, node_value_type value, size_type id, const Graph* graph): 
    position_(position), value_(value), id_(id), graph_(const_cast<Graph*>(graph))
    {}
    public:

  /** Update's the value of the corresponding Node.
   * @tparam[in] value, the value (of templated type node_value_type) to be updated
   * @return  none
   * 
   * Complexity: O(1)
   */
    void newval(const node_value_type& value){
      value_ = value;
      return;
    }
  };

  // Define an internal class which stores various information about an Edge. Assuming that
  // a user has "heavy"/memory-intensive Edge information to be stored, this will facilitate
  // an overall lightweight Graph implementation using the proxy design pattern.
  class InternalEdge {
    size_type edgeid_; // Unique edge id, the position in the edgevec_ vector
    size_type node1id_; // ID of one of the endpoints
    size_type node2id_; // ID of the other endpoint
    const Graph* graph_; // Pointer back to graph
    friend class Graph;
    // Default private constructor
    InternalEdge(){}
    // Initialize InternalEdge with all information about Edge.
    InternalEdge(size_type edgeid, size_type node1id, size_type node2id, const Graph* graph): 
    edgeid_(edgeid), node1id_(node1id), node2id_(node2id), graph_(const_cast<Graph*>(graph))
    {}
  };

};

#endif // CME212_GRAPH_HPP
