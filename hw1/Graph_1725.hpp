#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <set>
#include <map>

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

  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  struct pseudoedge;

  // pre-declare interior structures.
  struct interior_node;
  struct interior_edge;

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

  /*  Input variable type  */
  using node_value_type = V;





  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : numnodes_(0), numedges_(0){
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
    }

    /** Return this node's position. */
    const Point& position() const {
      return (this->pointer_to_graph->vecnodes_[unid_].inter_point);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(this->unid_ < this->pointer_to_graph->num_nodes());
      return this->unid_;
    }

    /**
     * @brief Returns the value for a node.
     *
     * @return A user-defined value for the given node
     *
     * O(1) runtime
     */
    node_value_type& value(){
      return(this->pointer_to_graph->vecnodes_[unid_].val);
    }
    const node_value_type& value() const{
      return(this->pointer_to_graph->vecnodes_[unid_].val);
    }

    /**
     * @brief Calculates the degree of a node - the number of adjacent edges to a node.
     *
     * @return A non-negative integer
     *
     * @post 0 <= result
     *
     * O(log(number of nodes))
     */
    size_type degree() const{
      size_type deg = 0;
      auto ptr = this->pointer_to_graph->mapedges_.find(unid_);
      if (ptr != this->pointer_to_graph->mapedges_.end()){
        deg = ptr->second.size();
      }
      return deg;
    }

    /**
     * @brief Returns a beginning edge incident iterator.
     *
     * @return an iterator to all the nodes adjacent to the node in question.
     *
     * @pre the called-on node is in the graph.
     *
     * Returns an iterator to the first of all the nodes that are adjacent to
     * the called node. O(log(number of nodes)) time.
     */
    incident_iterator edge_begin() const{
      auto ptr = this->pointer_to_graph->mapedges_.find(unid_);
      assert(ptr != this->pointer_to_graph->mapedges_.end());
      auto inptr = ptr->second.begin();
      return IncidentIterator(this->pointer_to_graph, inptr, this->unid_);
    }

    /**
     * @brief Returns an end edge incident iterator
     *
     * @return an iterator to one past the last node adjacent to the node in question.
     *
     * @pre the called-on node is in the graph.
     *
     * Returns an iterator to the end all the nodes that are adjacent to
     * the called node. O(log(number of nodes)) time.
     */
    incident_iterator edge_end() const{
      auto ptr = this->pointer_to_graph->mapedges_.find(unid_);
      assert(ptr != this->pointer_to_graph->mapedges_.end());
      auto inptr = ptr->second.end();
      return IncidentIterator(this->pointer_to_graph, inptr, this->unid_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if ((this->unid_ == n.unid_) & (this->pointer_to_graph == n.pointer_to_graph)){
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
      if ((this->unid_ < n.unid_) & (this->pointer_to_graph == n.pointer_to_graph)){
        return true;
      }
      return false;
    }

   private:
    friend class Graph;
    Graph* pointer_to_graph;
    size_type unid_;

    // private Node Constructor
    Node(const Graph* p_t_g, size_type unid_in)
      : pointer_to_graph(const_cast<Graph*>(p_t_g)), unid_(unid_in){
      }

  };


  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return numnodes_;
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
  Node add_node(const Point& position, const node_value_type& val = node_value_type()){
    this->vecnodes_.push_back({position, val, numnodes_});
    ++numnodes_;
    return Node(this, numnodes_-1);
  }





  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if ((n.pointer_to_graph == this) && (n.unid_ < this->numnodes_)){
      return true;
    }
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(0 <= i);
    assert(i < numnodes_);
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {}

    /** Return a node of this Edge */
    Node node1() const {
      return Node(p_t_g_edge_, node1_ID_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(p_t_g_edge_, node2_ID_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (((e.node1_ID_ == this->node1_ID_) && (e.node2_ID_ == this->node2_ID_)) ||
          ((e.node1_ID_ == this->node2_ID_) && (e.node2_ID_ == this->node1_ID_))){
            return true;
          }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (this->edge_ID_ < e.edge_ID_){
        return true;
      }
      return false;
    }

   private:
    friend class Graph;
    Graph* p_t_g_edge_;
    size_type edge_ID_;
    size_type node1_ID_;
    size_type node2_ID_;

    // private Edge constructor
    Edge(const Graph* p_t_g_in, size_type edge_ID_in, size_type nodeval1, size_type nodeval2)
      : p_t_g_edge_(const_cast<Graph*>(p_t_g_in)), edge_ID_(edge_ID_in), node1_ID_(nodeval1), node2_ID_(nodeval2){
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return numedges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(0 <= i);
    assert(i < num_edges());
    interior_edge ie = this->vecedges_[i];
    return Edge(this, i, ie.node1ID, ie.node2ID);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Check if @a a and @a b are valid nodes of this graph
    assert(this == a.pointer_to_graph);
    assert(this == b.pointer_to_graph);
    assert(a.unid_ < numnodes_);
    assert(b.unid_ < numnodes_);

    auto it = this->mapedges_.find(a.unid_);
      if (it != mapedges_.end()){
        auto itinner = it->second.find(b.unid_);
        if (itinner != it->second.end()){
          return true;
        }
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
    assert(this == a.pointer_to_graph);
    assert(this == b.pointer_to_graph);
    assert(a.index() < this->numnodes_);
    assert(b.index() < this->numnodes_);
    assert(a.index() != b.index());

    bool edgeexists = has_edge(a,b);
    if (!edgeexists){
      auto it = this->mapedges_.find(a.unid_);
      if (it != this->mapedges_.end()){
        // add the nodeID to the inner map
        it->second.insert({b.unid_, numedges_-1});
      }
      else{
        std::map<size_type, size_type> inner;
        inner.insert(std::make_pair(b.unid_,numedges_-1));
        this->mapedges_.insert({a.unid_, inner});
      }

      auto it2 = this->mapedges_.find(b.unid_);
      if (it2 != this->mapedges_.end()){
        // add the nodeID to the inner map
        it2->second.insert({a.unid_, numedges_-1});
      }
      else{
        std::map<size_type, size_type> inner2;
        inner2.insert(std::make_pair(a.unid_,numedges_-1));
        this->mapedges_.insert({b.unid_, inner2});
      }

      this->vecedges_.push_back(interior_edge(a.unid_, b.unid_, numedges_));
      numedges_++;
      return Edge(this, numedges_-1, a.unid_, b.unid_);
    }
    else{
      auto it = this->mapedges_.find(a.unid_);
      auto it2 = it->second.find(b.unid_);

      return Edge(this, it2->second, a.unid_, b.unid_);
    }
  }


  /**
   * @brief Returns the maximum distance from a root node in the graph.
   *
   * @return a non-negative integer
   *
   * @pre The graph has >1 node.
   *
   * Returns the maximum distance from a root note in the graph to all other nodes.
   * O(1) runtime.
   */
  size_type get_maxdist() const {
    return maxdist_;
  }

  /**
   * @brief Assigns the maximum distance from a root node in the graph.
   *
   * @param[in] in A non-negative integer calculated in the shortest path algorithm.
   *
   * Assigns the maximum distance private data member to the input variable.
   * O(1) runtime.
   */
  void set_maxdist(node_value_type in){
    this->maxdist_ = in;
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // remove all nodes from the graph
    vecnodes_.clear();
    assert(vecnodes_.size() == 0);
    numnodes_ = 0;

    // remove all edges from the graph
    mapedges_.clear();
    vecedges_.clear();
    numedges_ = 0;
    assert(mapedges_.size()==0);
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

    NodeIterator(const Graph* ptg_, typename std::vector<interior_node>::const_iterator veciterator_in) :
      NIptg_(const_cast<Graph*>(ptg_)), veciterator_(veciterator_in){

    }

    /**
     * @brief Returns a node accessed by the node iterator.
     *
     * @return A node in the graph.
     *
     * @pre The graph has at least 1 node.
     * @post 0 <= result.uniqueID <= the number of nodes.
     *
     * Derefernces the node iterator and returns the corresponding node.
     */
    Node operator*() const{
      interior_node temp = *veciterator_;
      return Node(NIptg_, temp.node_uniqueID);
    }

    /**
     * @brief Increments the node iterator.
     *
     * @return A node iterator
     *
     * @pre node iterator for a valid graph.
     * @post the node iterator no longer points to the node iterator that was called.
     *
     * Increments the node iterator.
     * O(1) runtime.
     */
    NodeIterator& operator++(){
      ++veciterator_;
      return *this;
    }

    /**
     * @brief Checks equivalence of two node iterators
     *
     * @param[in] ni Node Iterator
     * @return {true, false}
     *
     * @pre _ni_ and the node iterator that was called are from the same graph.
     * @post true if the iterators point to the same node.
     *
     * Checks equivalence of two node iterators.
     * O(1) runtime.
     */
    bool operator==(const NodeIterator& ni) const{
      return this->veciterator_ == ni.veciterator_;
    }

   private:
    friend class Graph;
    Graph* NIptg_;                                                    // Pointer to the graph
    typename std::vector<interior_node>::const_iterator veciterator_; // Iterator to internal node vector.
  };

  /**
   * @brief Returns a node iterator to the first node in the graph.
   *
   * @return a node iterator.
   *
   * @post if there are 0   nodes in the graph, result == node_end()
   *                    >=1 nodes in the graph, result != node_end()
   *
   * Returns a node iterator to the first node in the graph.
   * O(1) runtime.
   */
  node_iterator node_begin() const{
    return NodeIterator(this, vecnodes_.begin());
  }

  /**
   * @brief Returns a node iterator to (one past) the last node in the graph.
   *
   * @return a node iterator.
   *
   * @post if there are 0   nodes in the graph, result == node_begin()
   *                    >=1 nodes in the graph, result != node_begin()
   *
   * Returns a node iterator to (one past) the last node in the graph.
   * O(1) runtime.
   */
  node_iterator node_end() const{
    return NodeIterator(this, vecnodes_.end());
  }




  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>{
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

    // other Constructor:
    IncidentIterator(const Graph* IIptg_IN, std::map<size_type,size_type>::iterator mapiterIN, size_type IInode1IDIN) :
      IIptg_(const_cast<Graph*>(IIptg_IN)), mapiter_(mapiterIN), IInode1id(IInode1IDIN){

    }

    /**
     * @brief Dereferences an incident iterator and returns the corresponding Edge.
     *
     * @return An edge that the iterator was pointing to.
     *
     * @pre The graph has >= 1 edges in the graph.
     * @post The edge that was called on always is the first nodeID returned in _result_.
     *
     * Dereferences an incident iterator and returns the corresponding Edge.
     * O(1) runtime.
     */
    Edge operator*() const{
      auto i = this->mapiter_;
      return Edge(IIptg_, i->second, IInode1id, i->first);
    }

    /**
     * @brief Increments an Incident Iterator.
     *
     * @return An incident iterator one past the called-on Iterator.
     *
     * @post the incident iterator no longer points to the node iterator that was called.
     *
     * Increments an Incident Iterator.
     * O(1) runtime.
     */
    IncidentIterator& operator++(){
      this->mapiter_++;
      return *this;
    }


    /**
     * @brief Checks equality between an Incident Iterator and the Iterator called on.
     *
     * @param[in] iit Incident Iterator to be checked.
     * @return {true, false}
     *
     * @pre The Incident Iterators are from the same graph.
     * @post true if the Incident Iterators point to the same nodes.
     *
     * O(1) runtime.
     */
    bool operator==(const IncidentIterator& iit) const{
      if (this->IIptg_ == iit.IIptg_ && this->mapiter_ == iit.mapiter_){
        return true;
      }
      return false;
    }

   private:
    friend class Graph;
    Graph* IIptg_;
    std::map<size_type,size_type>::iterator mapiter_;
    size_type IInode1id;
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator(const Graph* EIptgIN, typename std::vector<interior_edge>::const_iterator veciterIN) :
      EIptg_(const_cast<Graph*>(EIptgIN)), veciter_(veciterIN){
    }

    /**
     * @brief Dereferences the Edge Iterator.
     *
     * @return the Edge to which the Iterator was pointing to.
     *
     * @pre the Edge Iterator is referencing a member of the appropriate graph.
     *
     * Dereferences the Edge Iterator. Returns an Edge with the corresponding
     * Edge unique ID, and both node IDs.
     * O(1) runtime.
     */
    Edge operator*() const{
      return Edge(EIptg_, (*veciter_).edge_UID, (*veciter_).node1ID, (*veciter_).node2ID);
    }

    /**
     * @brief Increments an Edge Iterator.
     *
     * @return An Edge Iterator one past the called-on Iterator.
     *
     * @post the Edge Iterator no longer points to the Edge from which it was called.
     *
     * O(1) runtime.
     */
    EdgeIterator& operator++(){
      ++veciter_;
      return *this;
    }

    /**
     * @brief Checks equality between Edge Iterators.
     *
     * @param[in] ei Edge Iterator to compare to.
     * @return {true, false}
     *
     * @pre The Edge Iterators are from the same graph.
     * @post <POSTCONDITION>
     *
     * O(1) runtime.
     */
    bool operator==(const EdgeIterator& ei) const{
      return (this->EIptg_ == ei.EIptg_) and (this->veciter_ == ei.veciter_);
    }

   private:
    friend class Graph;
    Graph* EIptg_;
    typename std::vector<interior_edge>::const_iterator veciter_;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /**
   * @brief Returns an Edge Iterator to the first Edge in the graph.
   *
   * @return Edge Iterator to the first Edge in the graph.
   *
   * Returns an Edge Iterator to the first Edge in the graph.
   * O(1) runtime.
   */
  edge_iterator edge_begin() const{
    return EdgeIterator(this, vecedges_.begin());
  }

  /**
   * @brief Returns an Edge Iterator to (one past) the last Edge in the graph.
   *
   * @return Edge Iterator to (one past) the last Edge in the graph.
   *
   * O(1) runtime.
   */
  edge_iterator edge_end() const{
    return EdgeIterator(this, vecedges_.end());
  }



 private:

   /** Interior Node
    * Holds the Point, a user-defined value, and a unique ID for each node.
    **/
  struct interior_node
  {
    Point inter_point;
    V val;
    size_type node_uniqueID;
  };

  /** Interior Edge
   * Holds the first and second nodeIDs, and a unique ID for each edge.
   **/
  struct interior_edge
  {
    interior_edge(size_type node1, size_type node2, size_type edge_UIDin) :
      node1ID(node1), node2ID(node2), edge_UID(edge_UIDin){};
    size_type node1ID;
    size_type node2ID;
    size_type edge_UID;
  };

  // Nodes:
  std::vector<interior_node> vecnodes_;
  size_type numnodes_;

  // Edges:
  std::vector<interior_edge> vecedges_;
  size_type numedges_;

  // Adjacency List
  //   outer map: key: node IDs, inner map key: node IDs, value: edgeID.
  std::map<size_type, std::map<size_type,size_type>> mapedges_;

  // maxdist_ is a value that holds the depth of a graph.
  node_value_type maxdist_;

};

#endif // CME212_GRAPH_HPP
