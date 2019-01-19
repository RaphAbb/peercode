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
class Graph {
 public: 
  using size_type = unsigned;
 
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  
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

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  //using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
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
  class Node {
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
      // HW0: YOUR CODE HERE ---DONE
      this->gptr_ = nullptr;
      this->id_ = -1;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE ---DONE
      //need to access the node via the proxy index
      // need helper variable in Graph
      if (gptr_ == nullptr){
        std::cout << "invalid point returned" << std::endl;
        return Point();
      } else {
        return (*gptr_).get_node_position(id_);
      }
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE  ---DONE
      if (gptr_ ==nullptr){
        return size_type(-1);
      } else{
        return (*gptr_).get_node_index(id_);
      }
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    // LW- can access private function and attribute of the obejct that is calling this operator
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE --DONE
      if ( (n.get_gptr() == gptr_) & (n.index() == (*this).index()) ) {
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
    // LW- can access private function and attribute of the obejct that is calling this operator
    bool operator<(const Node& n) const {
      // HW0: YOUR CODE HERE ---DONE
      if ( n.index() > (*this).index() ) {
        return true;
      } else {
        return false;
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    //pointer to graph
    const Graph* gptr_ = nullptr;

    //vector position/index as id
    size_type id_ = -1;

    //constructor for Graph to call in add node function
    Node(Graph* ptr, const size_type id)
      : gptr_(ptr), id_(id) { 
    }
    
    // constuctor for const functions 
    Node(const Graph* ptr, const size_type id)
      : gptr_(ptr), id_(id) { 
    }

    //method to return pointer 
    const Graph* get_gptr() const {
      return gptr_;
    }    

    //method to return nodes id_ (for edge proxy)
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

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position) {
    //get current number of nodes
    size_type graph_size = nodeVec.size();
    
    //construct a internal node inside the graph
    this->nodeVec.push_back( NodeInside( graph_size, position) );

    //construct a proxy and return it
    return Node(this, graph_size);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE ---DONE
    // check the pointer of node and then if the index less than size
    // will change if deletion occurs
    if ( (n.get_gptr() == this) & (n.index() < (*this).size()) ) {
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
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE ---DONE
    // create a new node proxy
   
    //check predondition
    if ( i >= (*this).size()) {
      throw std::runtime_error("node index outside graph");
    } 

    // This will only work as long as the index = ID
    return Node(this, i);
    //return Node();        // Invalid node
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
    Edge() {
      this->gptr_ = nullptr;
      this->edgeid_ = -1;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE --DONE
      // acess graph through pointer, feed graph edge ID and let it generate a node
      return (*gptr_).get_nodeA(edgeid_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE ---DONE
      // acess graph through pointer, feed graph edge ID and let it generate a node
      return (*gptr_).get_nodeB(edgeid_);
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
      } else if ( (( e.node1() == (*this).node1()) || (e.node1() == (*this).node2() )) && (( e.node2() == (*this).node1()) || (e.node2() == (*this).node2() )) ) {
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
      if ( e.get_index() > (*this).get_index()) {
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

    //constructor for Graph to use
    Edge( const Graph* ptr, size_type id) 
      : gptr_(ptr), edgeid_(id) {
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
    // HW0: YOUR CODE HERE --- DONE 
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
    if ( i >= (*this).num_edges() ){
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
    // HW0: YOUR CODE HERE ---DONE
    
    //Check preconditons
    if ( not( ((*this).has_node(a)) && ((*this).has_node(b)) ) ) {
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
    std::get<0>(tup2) = AID; 
    std::get<1>(tup2) = BID; 
    
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
    // HW0: YOUR CODE HERE --- DONE
    
    //Check preconditons
    if ( not( ((*this).has_node(a)) && ((*this).has_node(b)) && ((a==b)==false) ) ) {
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
    // HW0: YOUR CODE HERE ---DONE
    connection.clear();
    nodeVec.clear();
    edgeVec.clear();
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  
  //create an internal node class that will be used inside the Graph datastructures
  class NodeInside {	
    //typenames
    using point_type = Point;
    using index_type = unsigned;
    private:
      //store index
      index_type index_;
      //store point
      point_type position_;
      //valid checker (becasue cant use -1 as index b/c unsigned)
      bool valid_;
    public:
      //constructor for invalid
      NodeInside()
        : index_() , position_(), valid_(false) {
      }
      //constructor to use
      NodeInside( index_type ind, const point_type& pos )
        : index_(ind), position_(pos), valid_(true) {
      }

      //destructor????
      
      //get functions
      index_type get_index() const {
        return index_;
      }

      const point_type& get_pointref() const {
        return position_;
      }

      //set functions
      void set_index(const index_type ind){
        index_ = ind; //assign internal index to recived integer RIGHT WAY TO DO THIS?
      }

      void set_position(const point_type pos){
        position_ = pos; //asign interal index to recived point object RIGHT WAY TO DO THIS?
      }
  };
  
  // I should store the node unique ID (same as node) becasue that is the unique position
  // of the nodes

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

  // return the id/index of the exisiting node
  //basically has_edge but returns the index it finds
  size_type what_edge(const Node& a, const Node& b) const {
    if ( (*this).has_edge(a,b) == false ){
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
    std::get<0>(tup2) = AID; 
    std::get<1>(tup2) = BID; 

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


};

#endif // CME212_GRAPH_HPP
