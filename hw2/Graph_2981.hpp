#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <vector>
#include <unordered_map>
#include <algorithm> 
#include <cassert>
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
#include "filter_iterator.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V,typename E>
class Graph {
 private:
        // Declaring structs for internal storage
        struct fullNode;
        struct fullEdge;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  using node_value_type =V;
  using edge_value_type =E;
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
  using uid_type = size_type;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    next_NodeID_=0;
  }

  /** Default destructor */
  ~Graph() = default;
  
  // PREDICATES

  /** @class edgeExist
  *  @brief Predicate to check whether an Edge is still a valid Edge in the Graph
  *  @return Boolean checking existance of the Edge
  *
  *  Complexity: O(1)
  */
   struct edgeExist{ bool operator()(edge_type e){
           return(e.gr_->edgeVec[e.gr_->nodeVec[e.n1_uid_].otherNInd[e.n2_uid_]].exists);}
   };
  /** @class nodeExist
  *  @brief Predicate to check whether a Node is still a valid Node in the Graph
  *  @return Boolean checking existance of the Node
  *
  *  Complexity: O(1) 
  */
   struct nodeExist{
       bool operator()(node_type n){
           return(n.gr_->nodeVec[n.uidN_].exists);}
   };
    /** Helper function for constructing filter_iterators. This deduces the type of
     * the predicate function and the iterator so the user doesn't have to write it.
     * This also allows the use of lambda functions as predicates.
     *
     * Usage:
     * // Construct an iterator that filters odd values out and keeps even values.
     * std::vector<int> a = ...;
     * auto it = make_filtered(a.begin(), a.end(), [](int k) {return k % 2 == 0;});
     */
    template <typename Pred, typename Iter>
    filter_iterator<Pred,Iter> make_filtered(const Iter& it, const Iter& end,
                                             const Pred& p) const {
      return filter_iterator<Pred,Iter>(p, it, end);
    }
 
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
     * is occasionally useful to declare an invalid node, and assign a
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

    /** Return modifiable reference to this node's position. */
    Point& position(){ 
        assert(this->valid());
      return gr_->nodeVec[uidN_].pos;
    }
    /** Const function to return this node's position. */
    const Point& position() const {
        assert(valid());
     return gr_->nodeVec[uidN_].pos;
    }


    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
        assert(valid());
      return gr_->nodeVec[uidN_].ind;
    }

    /** Returns reference to stored Node value
     * @return the value, which can be changed
     */
    node_value_type& value(){
        assert(valid());
        return gr_->nodeVec[uidN_].nodeVal;}
    /** Returns const reference to stored node value
     * @return the value as const
     */
    const node_value_type& value() const{
        assert(valid());
        return gr_->nodeVec[uidN_].nodeVal;}
    /* Returns number of adjacent Nodes 
     * @return size_type number of adjacent Nodes 
     */
     size_type degree() const{
        assert(valid());
         return gr_->nodeVec[uidN_].assoc_edges.size();}
     /* Creates and returns filter_iterator of incident_iterators starting at this Node
      * @return filter_iterator of incident_iterators with its pointer to the first Edge from this Node
      */
     filter_iterator<edgeExist,incident_iterator> edge_begin() const{
        assert(valid());
        incident_iterator iit(gr_,gr_->nodeVec[uidN_].assoc_edges.begin());
        incident_iterator iitE(gr_,gr_->nodeVec[uidN_].assoc_edges.end());
        edgeExist ef{};
        while(iit!=iitE && !ef(*iit)){++iit;}
        return gr_->make_filtered(iit,iitE,ef);}
     /* Creates and returns filter_iterator of incident_iterators of this Node pointing to the end of its 
      * Edge list
      * @return filter_iterator of incident_iterators with its pointer to the end of the Edge list from this Node
      */
     filter_iterator<edgeExist,incident_iterator> edge_end() const{
        assert(valid());
        incident_iterator iit(gr_,gr_->nodeVec[uidN_].assoc_edges.end());
        incident_iterator iitE(gr_,gr_->nodeVec[uidN_].assoc_edges.end());
        edgeExist ef{};
        return gr_->make_filtered(iit,iitE,ef);}

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
        assert(n.valid() && this->valid());
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
    /** I chose to use the unique ID
     * of the node for ordering.
     */
        assert(n.valid() && valid());
         if(uidN_<n.uidN_)
             return true;
         else
             return false;
    }
   private:
        Graph* gr_;
        // Unique ID allows retrieval of Position and index
        size_type uidN_;
        Node(const Graph* gr, size_type uid) : uidN_(uid) {
                gr_=const_cast<Graph*>(gr);
        }
        /** Uses RI defined in class to check the validity of a Node object
        * @return Boolean specifying its validity
        *
        * Complexity: O(1)
        */
        bool valid() const {
            return (uidN_>=0 && uidN_<gr_->nodeVec.size() &&
                gr_->nodeVec[uidN_].ind<gr_->i2u.size() &&
                gr_->i2u[gr_->nodeVec[uidN_].ind]==uidN_);}

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
    return (size_type) i2u.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] nvt The new node's value, defaults to default value of type
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   * @post result_node.value() == @a nvt
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& nvt=node_value_type()) {
    fullNode newNode;
    Node ndAdd(this, next_NodeID_);
    newNode.pos=position;
    newNode.ind=(size_type) i2u.size();
    newNode.nd=ndAdd;
    newNode.nodeVal=nvt;
    nodeVec.push_back(newNode);
    i2u.push_back(next_NodeID_);
    next_NodeID_++;
    return newNode.nd;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // If a Node come from different Graph, or is not in the current Graph 
    // it returns false.
    assert(n.valid());
    if (this==n.gr_ && n.index()<size())
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
    return nodeVec[i2u[i]].nd;
  }
 /** Remove a Node from the Graph as well as all its incident Edges
  * Update index to unique ID vector and fullNode stored indices to ensure correct indexing
  * @return num_nodes() after deletion
  * @param[in] n Node to be removed
  * @pre All existing Nodes satisfy valid()
  * @post Any other copies of this Node object are invalidated
  * @post Any outstanding copies of removed incident Edges are invalidated
  * @post new num_nodes=old num_nodes()-1
  * @post new num_edges <= old num_edges()
  * @post Iterators on the i2u vector are invalidated
  * @post All indices are updated and remaining Nodes satisfy valid()
  *
  * Complexity: Up to O(N) (depending where the Node is in the fullNode vector)
 */
  size_type remove_node(const Node& n){
     assert(n.valid());
     //Remove incident Edges
     for (auto eit=n.edge_begin();eit!=n.edge_end();++eit)
         edgeVec[nodeVec[n.uidN_].otherNInd[(*eit).node2().uidN_]].exists=false;
         
     nodeVec[n.uidN_].exists=false; 
     //Update indices to unique ID vector
     auto it=i2u.erase(i2u.begin()+nodeVec[n.uidN_].ind); 
     //Update internal indices
     while (it!=i2u.end())
     {
         nodeVec[(*it)].ind--;
         ++it;
     }
     return num_nodes(); 
  }
 /** Remove a Node from the Graph as well as all its incident Edges using a node_iterator
  * Update index to unique ID vector and fullNode stored indices to ensure correct indexing
  * @return num_nodes() after deletion
  * @param[in] n_it node_iterator pointing to Node to be removed
  * @pre All existing Nodes satisfy valid()
  * @post Any other copies of this Node object are invalidated
  * @post Any outstanding copies of removed incident Edges are invalidated
  * @post new num_nodes=old num_nodes()-1
  * @post new num_edges <= old num_edges()
  * @post Iterators on the i2u vector are invalidated
  * @post All indices are updated and remaining Nodes satisfy valid()
  *
  * Complexity: Up to O(N) (depending where the Node is in the fullNode vector)
 */
  size_type remove_node(node_iterator n_it){
     return remove_node(*n_it); 
  }
 /** Remove an Edge from the Graph by its Nodes
  * @return num_edges() after deletion
  * @param[in] a, b Nodes that are on either end of the Edge
  * @pre @a a and @a b satisfy valid()
  * @pre has_edge(@a a,@a b) is true
  * @post Any outstanding copies of the removed Edge are invalidated
  * @post new num_edges() = old num_edges()-1
  *
  * Complexity: O(1)
 */
  size_type remove_edge(const Node& a, const Node& b){
     edgeVec[nodeVec[a.uidN_].otherNInd[(b.uidN_)]].exists=false;
     return num_edges(); 
  }
 /** Remove an Edge from the Graph
  * @return num_edges() after deletion
  * @param[in] e Edge to be removed
  * @pre @a e is a valid Edge 
  * @post Any outstanding copies of the removed Edge are invalidated
  * @post new num_edges() = old num_edges()-1
  *
  * Complexity: O(1)
 */
  size_type remove_edge(const Edge& e){
     return remove_edge(nodeVec[e.n1_uid_].nd,nodeVec[e.n2_uid_].nd);
  }
 /** Remove an Edge from the Graph using an edge_iterator
  * @return num_edges() after deletion
  * @param[in] e_it an edge_iterator pointing to the edge_iterator to be removed
  * @pre *(@a e_it) is a valid Edge 
  * @post Any outstanding copies of the removed Edge are invalidated
  * @post new num_edges() = old num_edges()-1
  *
  * Complexity: O(1)
 */
  size_type remove_edge(edge_iterator e_it){
     return remove_edge(*e_it);
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
    Edge() :  gr_(nullptr),n1_uid_(0),n2_uid_(0){
      // Reserve positive ID  values for valid edges
    }

    /** Returns reference to stored Edge value
     * @return the value, which can be changed
     */
    edge_value_type& value(){
       return gr_->edgeVec[gr_->nodeVec[n1_uid_].otherNInd[n2_uid_]].edgeVal;
      }
    /** Returns const reference to stored Edge value
     * @return the value as const
     */
    const edge_value_type& value() const {
       return gr_->edgeVec[gr_->nodeVec[n1_uid_].otherNInd[n2_uid_]].edgeVal;
      }

    /** Return length of Edge from Node positions */
    double length() const {
      return norm_2(gr_->nodeVec[n1_uid_].pos-gr_->nodeVec[n2_uid_].pos);
      }
    /** Return a node of this Edge */
    Node node1() const {
      assert(gr_->nodeVec[n1_uid_].nd.valid()); 
      return gr_->nodeVec[n1_uid_].nd;
    }

    /** Return a node of this Edge */
    Node node2() const {
      assert(gr_->nodeVec[n2_uid_].nd.valid()); 
      return gr_->nodeVec[n2_uid_].nd;
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
      // Use the node unique ID for ordering
      if (n1_uid_<e.n1_uid_ || (n1_uid_==e.n1_uid_ && n2_uid_<e.n2_uid_))
        return true;
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Edge(const Graph* gr, const size_type n1_uid, const size_type n2_uid) 
    : n1_uid_(n1_uid),n2_uid_(n2_uid){
                gr_=const_cast<Graph*>(gr);
    }
    Graph* gr_;
    size_type n1_uid_;
    size_type n2_uid_;
  };

  /** Return the total number of edges in the graph.
   * Could use an edge_iterator but complexity would be similar
   * @return Number of edges
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   * O(num_edges)
   */
  size_type num_edges() const {
    return (size_type) std::count_if(edgeVec.begin(),edgeVec.end(),[](fullEdge e){return e.exists;});
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i<num_edges());
    return edgeVec[i].edg;
  }



  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   * O(max degree() in Graph)
   */
  bool has_edge(const Node& a, const Node& b) const {
      assert(has_node(a) && has_node(b));
      Edge eTest(this,a.uidN_,b.uidN_);
      size_type smallInd=a.uidN_;
      // Search the incidence list of Node a to check for Edge
      for (auto eit=nodeVec[smallInd].nd.edge_begin();eit!=nodeVec[smallInd].nd.edge_end();++eit)
      {
          auto fnA=nodeVec[a.uidN_];
          auto eInd=fnA.otherNInd[b.uidN_];
          if (edgeVec[eInd].exists && eTest==(*eit))
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& evt=edge_value_type()) {
    assert(!(a==b)&& has_node(a) && has_node(b));
    Edge testEdge;
    bool aSmall;
    Edge edgeA(this,a.uidN_,b.uidN_);
    Edge edgeB(this,b.uidN_,a.uidN_);
    if (a<b) {   testEdge=Edge(this,a.uidN_,b.uidN_);
        aSmall=true;}
    else    {testEdge=Edge(this,b.uidN_,a.uidN_);
        aSmall=false;}
    // Occurs if the Edge already exists
    if (has_edge(a,b)) return testEdge;
    fullEdge newEdge;
    newEdge.edg=testEdge;
    newEdge.edgeVal=evt;
    if (aSmall) {newEdge.n1=a; newEdge.n2=b;}
    else {newEdge.n1=b; newEdge.n2=a;}
    // Add fullEdge to graph, and edges to each fullNode with that node as node1()
    edgeVec.push_back(newEdge);
    nodeVec[a.uidN_].assoc_edges.push_back(edgeA);
    nodeVec[a.uidN_].otherNInd[b.uidN_]=(size_type)edgeVec.size()-1;
    nodeVec[b.uidN_].assoc_edges.push_back(edgeB);
    nodeVec[b.uidN_].otherNInd[a.uidN_]=(size_type)edgeVec.size()-1;
    return testEdge;
  }
  

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
       i2u.clear();
       edgeVec.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator> {
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
/** Returns Node contained in the current fullNode struct pointed to by the node_iterator
 * @pre The node_iterator isn't pointing to the end of the Node list
 * @return Node contained in the current fullNode struct pointed to by the node_iterator
 */ 
      Node operator*() const {
          assert(nodeIt!=gr_->nodeVec.end());
            return (*nodeIt).nd;} 
      /** Increments internal iterator and returns this node_iterator
       * @return this node_iterator
       */
      NodeIterator& operator++() {
          nodeIt++;
          return *this;}
      /** Tests equality of this node_iterator with another based on Graph and internal iterator
       * @return whether they are equal
       */
      bool operator==(const node_iterator& ni) const {
          return (gr_==ni.gr_ && nodeIt==ni.nodeIt);} 
   private:
      // Use pointer to vector of fullNodes in the Graph
      typename std::vector<fullNode>::const_iterator nodeIt;
      Graph* gr_;
      NodeIterator(const Graph* gr, typename std::vector<fullNode>::const_iterator fNP){
                gr_=const_cast<Graph*>(gr);
                nodeIt=fNP;}
    friend class Graph;
  };
/** Returns filter_iterator of node_iterators pointing to first Node in the Graph
* Filter_iterator enables skipping of deleted Nodes
* @return filter_iterator of node_iterators pointing to first Node in the Graph
*/
     filter_iterator<nodeExist,node_iterator> node_begin() const{
          auto fnp=nodeVec.begin();
          node_iterator ni(this,fnp);
          node_iterator niE(this,nodeVec.end());
          nodeExist nf{};
          return make_filtered(ni,niE,nf);}
/** Returns filter_iterator of node_iterators pointing to the end of the fullNode vector in the Graph
* Filter_iterator enables skipping of deleted Nodes
* @return filter_iterator of node_iterators pointing to the end of the fullNode vector in the Graph
*/
     filter_iterator<nodeExist,node_iterator> node_end() const{
          node_iterator ni(this,nodeVec.end());
          node_iterator niE(this,nodeVec.end());
          nodeExist nf{};
          return make_filtered(ni,niE,nf);}

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


/** Returns Edge referenced by the internal iterator 
 * @pre The internal iterator isn't pointing to the end of the Edge list
 * @return Edge referenced by the internal iterator 
 */ 
     Edge operator*() const
     {
         return *edgeIt;}
     
     
      /** Increments internal iterator and returns this incident_iterator
       * @return this incident_iterator
       */
     IncidentIterator& operator++()
     { edgeIt++; return *this; }

      /** Tests equality of this incident_iterator with another based on Graph and internal iterator
       * @return whether they are equal
       */
     bool operator==(const IncidentIterator& IIT) const
     { return (gr_==IIT.gr_ && edgeIt==IIT.edgeIt) ;}


   private:
     // Use an internal iterator over edge_types contained in a fullNode struct
      typename std::vector<edge_type>::const_iterator edgeIt;
      Graph* gr_;
      IncidentIterator(const Graph* gr, typename std::vector<edge_type>::const_iterator iIt){
                gr_=const_cast<Graph*>(gr);
                edgeIt=iIt;}
    friend class Graph;
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : totally_ordered<EdgeIterator> {
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

    
    /** Return Edge contained in fullEdge pointed to by internal iterator
     * @return Edge contained in fullEdge pointed to by internal iterator
     */
     Edge operator*() const{
        return (*edgeIt).edg;}
     /** Increment internal iterator and return this EdgeIterator
      * @return this EdgeIterator
      */
     EdgeIterator& operator++(){ ++edgeIt; return *this;}
      /** Tests equality of this edge_iterator with another based on Graph and internal iterator
       * @return whether they are equal
       */
     bool operator==(const EdgeIterator& eit_) const{
        return (gr_==eit_.gr_ && edgeIt==eit_.edgeIt);}

   private:
    friend class Graph;
    // Use internal iterator over the Graph's 
    // Based on current Graph design, this is more efficient than methods using filter_iterator
      typename std::vector<fullEdge>::const_iterator edgeIt;
      Graph* gr_;
      EdgeIterator(const Graph* gr, typename std::vector<fullEdge>::const_iterator iIt){
                gr_=const_cast<Graph*>(gr);
                edgeIt=iIt;}
  };

        /** Returns filter_iterator of edge_iterators pointing to first Edge in the Graph
        * Filter_iterator enables skipping of deleted Edges
        * @return filter_iterator of edge_iterators pointing to first Edge in the Graph
        */
   filter_iterator<edgeExist,edge_iterator> edge_begin() const{
           edge_iterator eit{this,edgeVec.begin()};
           edge_iterator eit2{this,edgeVec.end()};
           edgeExist ef{};
           return filter_iterator<edgeExist,edge_iterator>(ef,eit,eit2);}
  /** Return filter_iterator of edge_iterators pointing to the end of the Graph's Edge list
   * @return filter_iterator of edge_iterators pointing to the end of the Graph's Edge list 
   */
   filter_iterator<edgeExist,edge_iterator> edge_end() const{
           edge_iterator eit{this,edgeVec.end()};
           edge_iterator eit2{this,edgeVec.end()};
           edgeExist ef{};
           return filter_iterator<edgeExist,edge_iterator>(ef,eit,eit2);}

 private:
//
    // Holds full information for a Node internally
     struct fullNode {
        Node nd;
        bool exists=true;
        size_type ind;
        Point pos;
        node_value_type nodeVal;
        // Map for constant time searching of incidence list of other Nodes connected to this Node
        std::unordered_map<size_type,size_type> otherNInd;
        std::vector<edge_type> assoc_edges;
        bool operator()(){return exists;}
     };
     // Holds full information for an Edge internally
     struct fullEdge {
        Edge edg;
        Node n1;
        Node n2;
        edge_value_type edgeVal;
        bool exists=true;
        bool operator()(){return exists;}
     };
     // Hold all Graph fullNodes in order
     std::vector<fullNode> nodeVec;
     // Hold all Graph fullEdges in order
     std::vector<fullEdge> edgeVec;
     // Get Node unique ID from index
     std::vector<uid_type> i2u;
     //   ID counters
     size_type next_NodeID_;
};

#endif // CME212_GRAPH_HPP
