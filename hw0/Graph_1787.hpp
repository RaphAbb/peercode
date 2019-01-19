#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm> 
#include <map>
#include <unordered_map>
// For the moment I am using a convenient Boost hashing function, but no Boost containers.
// I was told that this will be OK for HW0 although I will likely change it going forward.
#include <boost/functional/hash.hpp>
#include <cassert>
#include <iterator>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
 private:

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
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    //Leaving 0 for invalid nodes and edges
    next_NodeID_=1;
    next_EdgeID_=1;
    deletedNodes_=0;
    deletedEdges_=0;
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
       // Invalid nodes have an ID of 0
    Node() : gr_(nullptr), uidN_(0) {}

    /** Return this node's position. */
    const Point& position() const {
      // Ensure that this is a valid node
      // I was unsure whether or not to 'trust' the user, so the assert
      // statements have been left in. They do not noticeably slow the program
      // even in the large case.
      assert(uidN_>0 && gr_->nodeMap.find(uidN_)!=gr_->nodeMap.end());
      return gr_->nodeMap.at(uidN_);
    }


    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(uidN_>0 );
      //Guaranteed to throw exception if out of bounds
      return gr_->nodeIndices.at(uidN_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
        assert(uidN_>0 && n.uidN_>0);
        // Two nodes with same uidN_ are guaranteed to have the same index
        if (n.uidN_==uidN_ && n.gr_==this->gr_)
              return true;
        else
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
         assert(n.uidN_>0 && uidN_>0);
    /** I chose to use the unique ID
     * of the node for ordering.
     */
         if(uidN_<n.uidN_)
             return true;
         else
             return false;
    }
   private:

    // No internal elements should change 
        const Graph* gr_;
        // Unique ID allows retrieval of Position and index
        const size_type uidN_;
        Node(const Graph* gr, size_type uid) : uidN_(uid) {
                gr_=const_cast<Graph*>(gr);
        }
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return (size_type) nodeMap.size();
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
      // Adds node information to each relevant unordered_map,
      // increments the ID counter, and returns the node.
    nodeMap[next_NodeID_]=position;
    // The deleted node counter makes sure the index is correct
    nodeIndicesRev[next_NodeID_-1-deletedNodes_]=next_NodeID_;
    nodeIndices[next_NodeID_]=next_NodeID_-1-deletedNodes_;
    next_NodeID_++;
    return Node(this, next_NodeID_-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // If a Node come from different Graph, or is not in the current Graph 
    // it returns false.
    if (this==n.gr_ && nodeMap.find(n.uidN_)!=nodeMap.end())
            return true;
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i<size());
    return Node(this,nodeIndicesRev.at(i));
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
    Edge() :  gr_(nullptr),n1_uid_(0),n2_uid_(0),uidE_(0){
      // Reserve positive ID  values for valid edges
    }

    /** Return a node of this Edge */
    Node node1() const {
      assert(n1_uid_>0 && gr_->nodeMap.find(n1_uid_)!=gr_->nodeMap.end());
      Node n1{gr_,n1_uid_};
      return n1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      assert(n2_uid_>0 && gr_->nodeMap.find(n2_uid_)!=gr_->nodeMap.end());
      Node n2{gr_,n2_uid_};
      return n2;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
        if((node1()==e.node1() && node2()==e.node2())|| 
        (node2()==e.node1() && node1()==e.node2())) 
            return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      assert (e.n1_uid_>0 && e.n2_uid_>0); 
      // Use the unique ID for ordering
      if (uidE_<e.uidE_)
        return true;
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Edge(const Graph* gr, const size_type uidE, const size_type n1_uid, const size_type n2_uid) 
    : gr_(gr),n1_uid_(n1_uid),n2_uid_(n2_uid),uidE_(uidE){
    }
    // No internal elements should change 
    const Graph* gr_;
    const size_type n1_uid_;
    const size_type n2_uid_;
    // Unique ID allows for retrieval of index and Nodes
    const size_type uidE_;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return (size_type) edgeMap.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i<num_edges() && edgeIndicesRev.find(i)!=edgeIndicesRev.end());
    
    size_type uid=edgeIndicesRev.at(i);
    return Edge(this,uid,edgeMap.at(uid).first,edgeMap.at(uid).second);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
      assert(has_node(a) && has_node(b));
      // If the unique hash is found already in the map, the Edge exists
      std::size_t hashE=hashNodes(a,b);
      if (uidEMap.find(hashE)!=uidEMap.end())
                return true;
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
    assert(!(a==b)&& has_node(a) && has_node(b));
    //Create a unique hash based on the passed Nodes
    std::size_t eHash=hashNodes(a,b);
    //Create an ordered pair of the Nodes
    std::pair<size_type,size_type> edgePair;
    if (a<b) {
     edgePair=std::make_pair(a.uidN_,b.uidN_);
    }
    else { 
        edgePair=std::make_pair(b.uidN_,a.uidN_);}
    //Attempt to insert the new Edge
    auto outVals = uidEMap.insert({eHash,next_EdgeID_});
    // Occurs if insert is successful
    if (outVals.second){
        // Save relevant values, used the number of deleted Edges
        // to insure that the index is correct
        edgeMap[next_EdgeID_]=edgePair;
        edgeIndicesRev[next_EdgeID_-1-deletedEdges_]=next_EdgeID_;
        next_EdgeID_++;
        return Edge(this,next_EdgeID_-1,edgePair.first,edgePair.second);
    }
    // Occurs if the Edge already exists
    else
       return Edge(this,uidEMap[eHash], edgePair.first,edgePair.second);
  }
  
  /* Hash the unique IDs of two Node objects for use as a key */
  // No Boost containers used
  // Data structures will be likely refined in later assignments
  // to avoid use of Boost here
  std::size_t hashNodes(const Node& a, const Node& b) const{
    std::size_t edgeHash=0;
    if (a.uidN_<b.uidN_){
        boost::hash_combine(edgeHash,a.uidN_);
        boost::hash_combine(edgeHash,b.uidN_);
    }
    else{
        boost::hash_combine(edgeHash,b.uidN_);
        boost::hash_combine(edgeHash,a.uidN_);
    }
    return edgeHash;
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
      // Update the number of deleted Nodes and Edges
       deletedNodes_+=(size_type) nodeMap.size();
       deletedEdges_+=(size_type) edgeMap.size();
       // All Nodes and Edges are invalidated since their 
       // unique IDs don't match the Graph's maps.
       nodeMap.clear(); 
       edgeMap.clear();
       uidEMap.clear();
       nodeIndicesRev.clear();
       nodeIndices.clear();
       edgeIndicesRev.clear();
  }


 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  //   ID counters
     size_type next_NodeID_;
     size_type next_EdgeID_;
    // Deletion counters
     size_type deletedNodes_;
     size_type deletedEdges_;
     // unordered_maps for finding needed information in O(1) time
     std::unordered_map<size_type,Point> nodeMap;
     std::unordered_map<size_type,std::pair<size_type,size_type>> edgeMap;
     std::unordered_map<std::size_t,size_type> uidEMap;
     std::unordered_map<size_type,size_type> nodeIndicesRev; 
     std::unordered_map<size_type,size_type> nodeIndices;
     std::unordered_map<size_type,size_type> edgeIndicesRev;
};

#endif // CME212_GRAPH_HPP
