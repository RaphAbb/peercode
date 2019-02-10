#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <iostream>
#include <set>
#include <map>
#include <tuple>

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
  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

 private:

  //pre declare my internal types
  class NodeInside;
  class EdgeInside;

  //declare vector of Internal nodes
  std::vector<NodeInside> nodeVec;

  //declare  vector of Internal edges
  std::vector<EdgeInside> edgeVec;
  
  //declare map with
    //key tuple <size_type,size_type> node UIDs 
    //value size_type index
  std::map< std::tuple<size_type,size_type>, size_type > connection;
  
  //declare map of maps 
    // map1: key node_UID, value map2
    // map2: key node_UID, value sizetype index

  std::map< size_type, std::map<size_type, size_type> > nodenodeindex; 
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

  /** Synonym for tempalte class type value */
  using node_value_type = V;

  /* Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  //using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph<V>() {
    //defult behaviour is fine, my vectors get initalized to empty
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
      this->gptr_ = nullptr;
      this->id_ = -1;
    }

    /** Return this node's position. */
    const Point& position() const {
      if (gptr_ == nullptr){
      throw "invalid node called position()";
    } else {
      return (*gptr_).get_node_position(id_);
    }
  }

  /** Return this node's index, a number in the range [0, graph_size). */
  size_type index() const {
    if (gptr_ ==nullptr){
      return size_type(-1);
    } else{
      return (*gptr_).get_node_index(id_);
    }
  }

  // HW1: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_value_type& value();
  // const node_value_type& value() const;
  // size_type degree() const;
  // incident_iterator edge_begin() const;
  // incident_iterator edge_end() const;
  
  /** @brief return incident iterator at beginning
   *  @return IncidentIterator
   *  @post *edge_begin() == *(edge(nodenodeindex[id_].begin(),node(id_ ))
   */
  incident_iterator edge_begin() const{
    return (gptr_)->incident_start(id_ );
  }
  
  /** @brief return incident iterator at end
   *  @return IncidentIterator
   */
  incident_iterator edge_end() const {
    return (gptr_)->incident_end(id_ );
  }
  /**
   * @brief returns number of edges node has
   * @return size_type num of edges of node
   * calls to graph class
   * @post degree() == nodenodeindex[id_].size()
   */
  size_type degree() const{
    return gptr_->node_degree(id_);
  };
  
  /**
   * @brief returns a refernce to the nodes value of template type
   * @return returns a refernce of node_value_type
   */
  node_value_type& value(){
    return (gptr_)->node_value(id_);
    }

    /**
     * @brief returns a const refernce to the nodes value of template type
     * @return returns a cosnt refernce of node_value_type
     */
    const node_value_type& value() const{
      return gptr_->node_value_const(id_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE --DONE
      if ( (n.get_gptr() == gptr_) & (n.index() == this->index()) ) {
        return true;  
      } else {
        return false;
      }
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
      // HW0: YOUR CODE HERE ---DONE
      if ( n.index() > this->index() ) {
        return true;
      } else {
        return false;
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    
    //pointer to graph
    //edit to move const
    Graph* gptr_ = nullptr;

    //vector position/index as id
    size_type id_ = -1;

    /** constructor for Graph to call in add node function
     *  edit to add const_cast
     */
    Node(const Graph* ptr, const size_type id)
      : gptr_(const_cast<Graph*>(ptr)), id_(id) { 
    }
    
    /**method to return pointer */
    const Graph* get_gptr() const {
      return gptr_;
    }    

    /** method to return nodes id_ (for edge proxy) */
    size_type get_id() const {
      return id_;
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodeVec.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** @brief Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @tparam[in] val the new nodes value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    //get current number of nodes
    size_type graph_size = nodeVec.size();
    
    //construct a internal node inside the graph
    this->nodeVec.push_back( NodeInside( graph_size, position, val) );

    //construct a proxy and return it
    return Node(this, graph_size);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // check the pointer of node and then if the index less than size
    // will change if deletion occurs
    if ( (n.get_gptr() == this) & (n.index() < this->size()) ) {
      return true;
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
  Node node(size_type i) const { //remove const?
    //check predondition
    if ( i >= this->size()) {
      throw std::runtime_error("node index outside graph");
    } 

    // This will only work as long as the index = ID
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
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
      this->gptr_ = nullptr;
      this->edgeid_ = -1;
      this->genby_ = 'C';
    }

    /** @brief Return a node of this Edge 
     *  
     *  If genby_ == 'C' which node it returns is determined by the underlying storage
     *  if genby_ == 'A'|| genby_'B', node1() will return the node that spawned the edge
     */
    Node node1() const {
      // acess graph through pointer, feed graph edge ID and let it generate a node
      if (genby_ == 'C') {
        return (*gptr_).get_nodeA(edgeid_);
      } else if (genby_ == 'A') {
        return (*gptr_).get_nodeA(edgeid_);
      } else {
        return (*gptr_).get_nodeB(edgeid_);
      }
    }

    /** @brief Return a node of this Edge 
     *  
     *  If genby_ == 'C' which node it returns is determined by the underlying storage
     *  if genby_ == 'A'|| genby_'B', node1() will return the node that spawned the edge
     */
    Node node2() const {
      if (genby_ == 'C') {
        return (*gptr_).get_nodeB(edgeid_);
      } else if (genby_ == 'A') {
        return (*gptr_).get_nodeB(edgeid_);
      } else {
        return (*gptr_).get_nodeA(edgeid_);
      }
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // HW0: YOUR CODE HERE ---DONE
      //first check if graph is same by pointer, and then compare nodes
      if ( e.get_gptr() != gptr_) {
        return false;
      } else if ( ( (e.node1() == this->node1()) 
                  || (e.node1() == this->node2()) ) 
                  && ( (e.node2() == this->node1()) 
                  || (e.node2() == this->node2()) ) ) {
        return true;
      } else {
        return false;
      }
      //comparing edge IDs could be complicated because it is a 2x2 grid 

    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // HW0: YOUR CODE HERE ---DONE
      if ( e.get_index() > this->get_index()) {
        return true;
      } else {
        return false;  
      }
    
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    //pointer to graph
    const Graph* gptr_ = nullptr;

    //vector position/index as id
    size_type edgeid_ = -1;

    /** factor to determine which node to return
     *  If 2 arg constructor used, set to C and order or return doesnt mattter
     *  If 3 argument constructor use, will be set to A or B indicating which node should be returned by
     *    node1() and node2()
     */
    char genby_ = 'C';

    /**
     * @brief constructor where it doesnt matter order of node reutnr 
     */
    Edge( const Graph* ptr, size_type id) 
      : gptr_(ptr), edgeid_(id), genby_('C') {
    }
    
    /**
     * @brief constructor for Edge that will return nodes in a specific order
     */
    Edge( const Graph* ptr, size_type id, Node generator) 
      : gptr_(ptr), edgeid_(id), genby_('C') {
      if (this->node1() == generator) {
        genby_ = 'A';
      } else {
        genby_ = 'B';
      }
    }

    //method to return pointer 
    const Graph* get_gptr() const {
      return gptr_;
    }
    
    //get index of edge from EdgeInside
    size_type get_index() const {
      return (*gptr_).get_edge_index( edgeid_ );  
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edgeVec.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE --- DONE
    
    //check predondition
    if ( i >= this->num_edges() ){
      throw std::runtime_error("edge index outside graph");
    } 
    
    // This will only work as long as the index = ID
    return Edge(this, i);
  }
  
  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    //Check preconditons
    if ( not( (this->has_node(a)) && (this->has_node(b)) ) ) {
      throw std::runtime_error("nodes are not valid nodes for this graph");
    }
    
    // take input nodes and pull out IDs
    size_type AID, BID;
    AID = a.get_id();
    BID = b.get_id();
    
    //declare tuples
    std::tuple<size_type,size_type> tup1, tup2;
    //tup1
    std::get<0>(tup1) = AID; 
    std::get<1>(tup1) = BID; 
    //tup2
    std::get<0>(tup2) = BID; 
    std::get<1>(tup2) = AID; 
    
    //perform check
    if ( connection.find(tup1) != connection.end()){
      return true;
    }
    if ( connection.find(tup2) != connection.end()){
      return true;
    }
    else {return false;}    
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
    
    //Check preconditons
    if ( not( (this->has_node(a)) && (this->has_node(b)) && ((a==b)==false) ) ) {
      throw std::runtime_error("nodes are not distinct, valid nodes for this graph. Cannot create edge");
    }

    //NEED TO CHECK HAS EDGE
    if ( has_edge(a,b) == true) {
      return Edge( this,what_edge(a,b) );
    }
    // take input nodes and pull out IDs
    size_type AID, BID;
    AID = a.get_id();
    BID = b.get_id();
    
    //index and ID to use
    size_type edge_count  = edgeVec.size();

    //push back edge vec
    edgeVec.push_back( EdgeInside( edge_count, AID, BID) );

    //work with nodenodeindex for AID
    auto lvl1 = nodenodeindex.find(AID);
    //if AID not in map add it
    if ( lvl1 == nodenodeindex.end() ){     //add index AID with value of map
      nodenodeindex[AID][BID] = edge_count; //add BID->edge_count
    } else {                                // else check if B exists
      if (nodenodeindex[AID].find(BID) == nodenodeindex[AID].end()) {  
                                            //if BID does not exist then add
        nodenodeindex[AID][BID] = edge_count;
      }
    }
    
    //work with nodenodeindex for BID
    lvl1 = nodenodeindex.find(BID);
    //if BID not in map add it
    if ( lvl1 == nodenodeindex.end() ){    //add index BID with value of map
      nodenodeindex[BID][AID] = edge_count;//add AID->edge_count
    } else {                               // else check if A exists
      if (nodenodeindex[BID].find(AID) == nodenodeindex[BID].end()) {  
                                           //if AID does not exist then add
        nodenodeindex[BID][AID] = edge_count;
      }    
    }

    //is to maintain connection
    //create tuple key
    std::tuple<size_type,size_type> tup;
    std::get<0>(tup) = AID; 
    std::get<1>(tup) = BID; 
    
    //add new peice to map
    connection[tup] = edge_count;
    
    //construct a proxy and return it
    return Edge(this, edge_count);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    connection.clear();
    nodeVec.clear();
    edgeVec.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
    
    // type defintion for conveience
    using size_type         = unsigned;

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
      
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /**
     * @brief return Node object of the node at the current positon
     * by dereferncing the position we're getting a node object that referes
     * the underlying node at the position
     */
    Node operator*() const {
      return gptr_->node(pos_);
    }

    /**
     * @brief increment the position
     * @return a reference to the node iterator that called it 
     */
    NodeIterator& operator++(){
      pos_++;
      return *this;
    }

    /**
     * @brief evaluates true if the graph pointer and postion are the equal
     * @return boolean 
     */
    bool operator==(const NodeIterator& ni) const{
      if ( (ni.pos_ == this->pos_) && (ni.gptr_ == this->gptr_) ){
        return true;
      } else {
        return false;
      }
    }
    
    /**
     * @brief evaluates false if the graph pointer and postion are the equal
     * @return boolean 
     */
    bool operator!=(const NodeIterator& ni) const{
      if ( (ni.pos_ == this->pos_) && (ni.gptr_ == this->gptr_) ){
        return false;
      } else {
        return true;
      }
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    
    /** @brief constructor of valid NodeIterator 
     *  @param ptr pointer to the graph
     *  @param start_pos size_type value of the position for the iterator to start at
     *  
     *  receives a graph and starting value and constructs a iterator
     */
    NodeIterator(const Graph* ptr, size_type start_pos)
      : pos_(start_pos), gptr_(ptr)  {
    }
    
    /** position inside container */
    size_type pos_;
    
    /** pointer to graph */
    const Graph* gptr_;

    
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  
  /**
   *  @brief returns iterator with position 0
   *  @return node_iterator object
   */
  node_iterator node_begin() const{
    return NodeIterator(this,0);
  }

  /**
   *  @brief returns iterator with position num_nodess()
   *  @return node_iterator object
   *
   *  num_edges returns the size of the node vector, so it is one beyond the end 
   *  of the zero indexed vector
   */
  node_iterator node_end() const{
    return NodeIterator( this,this->num_nodes() );
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator {
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
    
    //working structure is to use the map operator to track and reuturn node Bs
    //when dereferenceing will return edge proxies cosntructed with the previous method 
    
    /**
     * @brief returns an Edge proxy for the current location
     * @return Edge
     * will return an Edge with node 1 as the node that spawned this Incident operator
     */
    Edge operator*() const {
      //get edge ID
      size_type edgeID = mapiter_->second; 
      //use special edge constructor
      return Edge( gptr_, edgeID, gptr_->node(nodeid_) );
    }

    /**
     * @brief increments iterator
     * @return reference to this indicent operator
     * increments underlying map iterator
     */
    IncidentIterator& operator++() {
      ++mapiter_;
      return *this;
    }

    /**
     * @breif returns true if all underlying elements are equal
     * @return bool
     * comarpe each data member
     */
    bool operator==(const IncidentIterator& it) const {
      if ( (it.gptr_ == gptr_) && (it.nodeid_ == nodeid_) && (it.mapiter_ == mapiter_) ){
        return true;
      } else {
        return false;
      }
    }
    
    /**
     * @breif returns true if all underlying elements are equal
     * @return bool
     * comarpe each data member
     */
    bool operator!=(const IncidentIterator& it) const {
      if ( (it.gptr_ == gptr_) 
            && (it.nodeid_ == nodeid_) 
            && (it.mapiter_ == mapiter_) ){
        return false;
      } else {
        return true;
      }
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE

    /** pointer to graph */
    const Graph* gptr_;
    
    /** uid of the node that spawned it */
    size_type nodeid_;

    /** map iterator to interact with */
    std::map<size_type,size_type>::const_iterator mapiter_;

    //made copy, then then made iterators. the equanlty check of the iterators rlys on memory location, 
    // so the problem is from the copys not having samem ememory location. 

    IncidentIterator(Graph* ptr, bool begin, size_type nid)
      : gptr_(ptr), nodeid_(nid) {
      
      if (begin == true){
        mapiter_ = ((gptr_->nodenodeindex.find(nodeid_))->second).begin();    
      } else {
        mapiter_ = ((gptr_->nodenodeindex.find(nodeid_))->second).end();    
      }
    }

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
    
    // type def for convenience
    using size_type = unsigned;

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const
    
    /**
     * @brief return Edge object of the edge at the current positon
     * by dereferncing the position we're getting a Edge object that referes the underlying edge at the position
     */
    Edge operator*() const {
      return gptr_->edge(pos_);
    }

    /**
     * @brief increment the position
     * @return a reference to the EdgeIterator that called it 
     */
    EdgeIterator& operator++(){
      pos_++;
      return *this;
    }

    /**
     * @brief evaluates true if the graph pointer and postion are the equal
     * @return boolean 
     */
    bool operator==(const EdgeIterator& ni) const{
      if ( (ni.pos_ == this->pos_) && (ni.gptr_ == this->gptr_) ){
        return true;
      } else {
        return false;
      }
    }
    
    /**
     * @brief evaluates false if the graph pointer and postion are the equal
     * @return boolean 
     */
    bool operator!=(const EdgeIterator& ni) const{
      if ( (ni.pos_ == this->pos_) && (ni.gptr_ == this->gptr_) ){
        return false;
      } else {
        return true;
      }
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    
    /** position inside container */
    size_type pos_;
    
    /** pointer to graph */
    const Graph* gptr_;

    /** @brief constructor of valid EdgeIterator 
     *  @param ptr pointer to the graph
     *  @param start_pos size_type value of the position for the iterator to start at
     *  
     *  receives a graph and starting value and constructs a iterator
     */
    EdgeIterator(const Graph* ptr, size_type start_pos)
      : pos_(start_pos),gptr_(ptr) {
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  
  /**
   *  @brief returns iterator with position 0
   *  @return EdgeIterator object
   *  @post *edge_begin() == gptr->edge(0)
   */
  edge_iterator edge_begin() const{
    return EdgeIterator(this,0);
  }

  /**
   *  @brief returns iterator with position num_edges()
   *  @return edge_iterator object
   *  @post *edge_end() == Vgptr->edge(0)
   *
   *  num_edges returns the size of the edge vector, so it is one beyond the end 
   *  of the zero indexed vector
   */
  edge_iterator edge_end() const{
    return EdgeIterator( this,this->num_edges() );
  }

 private:

  /**
   * This is the internal node class
   * contains
   *    a index 
   *    a point object
   *    a template type object
   */
  class NodeInside {
    
    using point_type = Point;
    using index_type = unsigned;
    using node_value_type = V;

    private:
      /** data members */
      index_type index_; /** index value */
      point_type position_; /** Point object */
      node_value_type value_;
      bool valid_; /** boolen value if valid or not, not yet used */ 
    
    public:
      /**
       * @brief constructs an invalid NodeInside
       * 
       * null initalizes index, position, value, and sets valid bool to false
       */
      NodeInside()
        : index_() , position_(), value_(), valid_(false) {
      }
      
      /**
       * @brief constructs an valid NodeInside with value defualt intalized
       * 
       * @param ind index of the node
       * @param pos reference to Point object with position of the node
       *
       * initalizes index, position with given params. sets valid bool to true
       * defualt initalized the value 
       */
      //constructor to use
      NodeInside( index_type ind, const point_type& pos )
        : index_(ind), position_(pos), value_(),  valid_(true) {
      }

      /**
       * @brief constructs an valid NodeInside with V value_ fed 
       * 
       * @param ind index of the node
       * @param pos reference to Point object with position of the node
       * @tparam val reference to node_value_type object 
       *
       * initalizes index, position, value  with given params. sets valid bool to true
       */
      NodeInside( index_type ind, const point_type& pos, const node_value_type& val )
        : index_(ind), position_(pos), value_(val),  valid_(true) {
      }
      
      /**
       * @brief returns the index
       * @return index_type index_ 
       */
      index_type get_index() const {
        return index_;
      }
      
      /**
       * @brief returns constant reference to the Point object
       * @return point_type& reference to the position point object 
       */
      const point_type& get_pointref() const {
        return position_;
      }

      /**
       * @brief returns const reference to value_ object
       * @return const node_value_type& reference
       */
      const node_value_type& get_valueref() const{
        return value_;
      }
      
      /**
       * @brief returns reference to value_ object
       * @return node_value_type& reference
       */
      node_value_type& get_valueref() {
        return value_;
      }

      /**
       * @brief changes the index value
       * @post this->get_index() = ind
       *
       * this isnt used yet
       */
      void set_index(const index_type ind){
        index_ = ind; //assign internal index to recived integer RIGHT WAY TO DO THIS?
      }

      /**
       * @brief changes the value of the node
       * 
       * Not implimented yet, will need to do so in HW1
       */
      void set_value(){
      } 
  };
  

  class EdgeInside {

    //type names
    using point_type = Point;
    using index_type = unsigned;
    private:
      //store index
      index_type index_;
      //internal set of points is a proxy Node
      index_type nodeAid_;
      index_type nodeBid_;
    public:
      //constructor to use
      EdgeInside(index_type ind, const index_type aa, const index_type  bb)
        : index_(ind), nodeAid_( aa ) , nodeBid_( bb ) {
      }

      //constructor for invalid?
      //destructor?
      
      //get functions
      index_type get_nodeAid() const {
        return nodeAid_;
      }

      index_type get_nodeBid() const {
        return nodeBid_;
      }
      
      index_type get_index() const {
        return index_;
      }

  };
  
  //
  //HW0 support functions
  //

  // get a position from a node
  const Point& get_node_position(size_type i) const {
    return nodeVec.at(i).get_pointref();
  }
  
  // get a position from a node
  size_type get_node_index(size_type i) const {
    return nodeVec.at(i).get_index();
  }
  
  //get nodeA from an edge
  Node get_nodeA( size_type i ) const {
    return Node(this, edgeVec.at(i).get_nodeAid()); 
  }
  
  //get nodeB from an edge
  Node get_nodeB( size_type i ) const {
    return Node(this, edgeVec.at(i).get_nodeBid()); 
  }
  
  //get egde index
  size_type get_edge_index(size_type i ) const{
    return edgeVec.at(i).get_index();
  }

/*
  // return the id/index of the exisiting node
  //basically has_edge but returns the index it finds
  size_type what_edge(const Node& a, const Node& b) const {
    if ( this->has_edge(a,b) == false ){
      throw "edge does not exist";
    }
    if ( connection.find(std::get<0>(get_tups(a,b))) != connection.end()){
      return connection.find(std::get<0>(get_tups(a,b)))->second;
    }
    if ( connection.find(std::get<1>(get_tups(a,b))) != connection.end()){
      return connection.find(std::get<1>(get_tups(a,b)))->second ;
    }
    else {
      throw "edge does not exist, dont call this function";
    }
  }
*/
 // return the id/index of the exisiting node
  //basically has_edge but returns the index it finds
  size_type what_edge(const Node& a, const Node& b) const {
    if ( this->has_edge(a,b) == false ){
      throw "edge does not exist";
    }

    // take input nodes and pull out IDs
    size_type AID, BID;
    AID = a.get_id();
    BID = b.get_id();

    //declare tuples
    std::tuple<size_type,size_type> tup1, tup2;
    //tup1
    std::get<0>(tup1) = AID;
    std::get<1>(tup1) = BID;
    //tup2
    std::get<0>(tup2) = BID;
    std::get<1>(tup2) = AID;

    if ( connection.find(tup1) != connection.end()){
      return connection.find(tup1)->second;
    }
    if ( connection.find(tup2) != connection.end()){
      return connection.find(tup1)->second;
    }
    else {
      throw "edge does not exist, dont call this function";
    }
}

  //
  //HW1 support functions
  //
  
  /** @brief get node value 
   *  @return reference to node value
  */
  node_value_type& node_value(size_type id){
    return nodeVec.at(id).get_valueref(); 
  }
  
  /** @brief get node value 
   *  @return const reference to node value
  */
  const node_value_type& node_value_const(size_type id) const{
    return nodeVec.at(id).get_valueref(); 
  }
  /** @brief returns number of edges of a node
   *  @return size_type 
  */
  size_type node_degree(size_type id) const {
    auto lvl1 = nodenodeindex.find(id);
    if (lvl1 == nodenodeindex.end()){
      throw "node doesnt exist";
    } else {
      return (lvl1->second).size();
    }
  }
  
  /**@brief generate a incedent operator for a specific node, at begining
   * @param nodeid unique ID of generating node 
   */
  incident_iterator incident_start(size_type nodeid){
    return IncidentIterator(this, true, nodeid);
  }
  
  /**@brief generate a incedent operator for a specific node, at end
   * @param nodeid unique ID of generating node 
   */
  incident_iterator incident_end(size_type nodeid){
    return IncidentIterator(this, false, nodeid);
  }
  

};

#endif // CME212_GRAPH_HPP
