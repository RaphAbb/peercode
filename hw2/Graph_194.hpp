#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
* @brief An undirected graph type
*/

#include <algorithm>
#include <vector>
#include <cassert>
#include <tuple>

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

  /**Predeclaration of internal struct of nodes. */
  struct InternalNode;

  /**Predeclaration of internal struct of edges. */
  struct InternalEdge;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {}

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
  using node_value_type = V;
  class Node : private totally_ordered<Node>{
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
    Node() {}

    /** Return this node's position. */
    const Point& position() const {return GetNode().Position_;}

    /** Return this node's position. */
    Point& position() {return GetNode().Position_;}

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {return GetNode().idx_;}

    /** Return this node's value
    * @pre Graph_->has_node(*this)
    * @return value of the current node
    */
    node_value_type& value() {return GetNode().Value_;}

    /** Return this node's value
    * @pre Graph_->has_node(*this)
    * @return value of the current node
    */
    const node_value_type& value() const {return GetNode().Value_;} //#TODO why this?

    /** Return the degree of the current node (number of adjacent vertexes)
    * @pre Graph_->has_node(*this)
    * @return degree of current node
    */
    size_type degree() const {return GetNode().AdjEdges_.size();};

    /** Set this node's value
    * @pre Graph_->has_node(*this)
    */
    void set_value(const node_value_type& value) {GetNode().Value_ = value;}

    /** Return the iterator of the first edge incident to a node
    * @pre Graph_->has_node(*this)
    * @post iterator to the first edge incident to a node
    */
    incident_iterator edge_begin() const {return IncidentIterator(Graph_, 0 ,uid_);}

    /** Return the iterator of the last edge incident to a node
    * @pre Graph_->has_node(*this)
    * @post iterator to the last edge incident to a node
    */
    incident_iterator edge_end() const {return IncidentIterator(Graph_, Graph_->Nodes_[uid_].AdjEdges_.size() ,uid_);}

    /** Test whether this node and @a n are equal.
    *
    * Equal nodes have the same graph and the same index.
    */
    bool operator==(const Node& n) const {
      if (this->Graph_ == n.Graph_ and this->uid_ == n.uid_) return true;
      else return false;
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
      if (Graph_ == n.Graph_) {return (uid_ < n.uid_);}
      else {return (Graph_ < n.Graph_);}
    }

  private:

    bool valid() const {return 1;}
    //   return uid_>= 0 && uid_< Graph_->Nodes_.size()
    //   && Graph_->Nodes_.at(uid_).idx_ < Graph_->i2u_Nodes_.size()
    //   && Graph_->i2u_Nodes_.at(Graph_->Nodes_.at(uid_).idx_) == uid_
    //   && Graph_->Nodes_[uid_].isActive_;
    // }
    void printDebug() const {}
    //   std::cout << "uid_" <<uid_<< std::endl;
    //   std::cout << "Checking Validity" << Graph_->Nodes_[uid_].isActive_<< std::endl;
    //   std::cout << "Graph_->Nodes_.size()" << Graph_->Nodes_.size()<<std::endl;
    //   std::cout << "Graph_->Nodes_.size()" << Graph_->Nodes_.size()<<std::endl;
    //   for (auto i = 0; i < Graph_->i2u_Nodes_.size(); ++i){std::cout<<"Printing " << Graph_->i2u_Nodes_[i]<<std::endl;}
    //   std::cout << "i2u_Nodes_.at(Graph_->Nodes_.at(uid_).idx_)"<<Graph_->i2u_Nodes_.at(Graph_->Nodes_.at(uid_).idx_) << std::endl;
    //   std::cout << "uid_"<< uid_ << std::endl;
    // }

    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    graph_type* Graph_;
    size_type uid_; // Unique ID
    Node(graph_type* Graph, size_type uid) : Graph_(Graph), uid_(uid) {}

    InternalNode& GetNode(){
      if(!valid()) {printDebug(); assert(1==0);}
      return Graph_->Nodes_[uid_];
    };

    const InternalNode& GetNode() const {
      if(!valid()) {printDebug(); assert(1==0);}
      return Graph_->Nodes_[uid_];
    };

  };

  /** Return the number of nodes in the graph.
  *
  * Complexity: O(1).
  */
  size_type size() const {return i2u_Nodes_.size();}

  /** Synonym for size(). */
  size_type num_nodes() const {return size();}

  /** Add a node to the graph, returning the added node.
  * @param[in] position The new node's position
  * @post new num_nodes() == old num_nodes() + 1
  * @post result_node.index() == old num_nodes()
  *
  * Complexity: O(1) amortized operations.
  */
  Node add_node(const Point& position,const node_value_type& Value = node_value_type()) {

    std::vector<size_type> EmptyVecNTE;
    auto idx = num_nodes();
    auto NewNode = InternalNode(position,EmptyVecNTE,Value,idx,true);
    Nodes_.push_back(NewNode);
    i2u_Nodes_.push_back(Nodes_.size()-1);
    return Node(this,idx); // zero based index
  }

  /** Determine if a Node belongs to this Graph
  * @return True if @a n is currently a Node of this Graph
  *
  * Complexity: O(1).
  */
  bool has_node(const Node& n) const {
    return (this == n.Graph_ and n.GetNode().isActive_); // check they point to the same graph
  }

  /** Return the node with index @a i.
  * @pre 0 <= @a i < num_nodes()
  * @post result_node.index() == i
  *
  * Complexity: O(1).
  */
  Node node(size_type i) const {
    assert(i <= num_nodes());
    return Node(const_cast<graph_type*>(this),i2u_Nodes_[i]);
  }

  /** Remove node @a n from the Graph.
  * @pre 0 <= @a n.index() < num_nodes()
  * @result 0 if @a n does not belong to the current graph, 1 otherwise
  * @post new num_nodes() == old num_nodes() - 1
  * @post all edges incident to @a n are eliminated
  * @post new num_edges() == old num_edges() - n.degree()
  *
  * Note: invalidates deleted node and its adjiancent edges
  * Complexity: O(@a n.degree() * d) where d is the maximum degree of the graph.
  */
  size_type remove_node(const Node& n) {
    // Check if the node is present
    if (!has_node(n)){return 0;}
    // needed for assertion
    auto old_numEdges = num_edges();
    // Copy Incident Edges for Removal
    std::vector<Edge> EdgesToRemove;
    for (auto iEdge = n.edge_begin(); iEdge != n.edge_end(); ++iEdge) {EdgesToRemove.push_back(*iEdge);}
    // now remove the edges and check operation successful
    for (size_type i = 0; i < EdgesToRemove.size(); ++i) {assert(remove_edge(EdgesToRemove[i]) == 1);}
    // check removed the correct number of edges
    assert(num_edges() == old_numEdges - EdgesToRemove.size());
    // copy index
    auto n_index = n.index();
    // swap this edge index to the end and delete it
    std::swap(i2u_Nodes_.at(n_index),i2u_Nodes_.at(i2u_Nodes_.size()-1));
    // Update last node index
    Nodes_.at(i2u_Nodes_.at(n_index)).idx_ = n_index;
    i2u_Nodes_.pop_back();
    // Invalidate Internal Representation of the Node
    Nodes_.at(n.uid_) = InternalNode(Point(),std::vector<size_type>(),node_value_type(),-1,false);
    return 1;
  }

  /** Remove edge with nodes @a node1 and @a node2 from the Graph.
  * @pre edge is an edge in the graph
  * @result 0 if edge does not belong to the current graph, 1 otherwise
  * @post new num_edges() == old num_edges() - 1 if preconditions are satisfies
  *
  * Note: invalidates removed edge
  * Complexity: O(@a e.node1().degree() * @a e.node2().degree())
  */
  size_type remove_edge(const Node& node1, const Node& node2){
    if (!(has_edge(node1,node2))) {return 0;}
    auto e = find_edge(node1,node2);
    return remove_edge(e);
  }

  /** Remove edge @a e from the Graph.
  * @pre @e is an edge in the graph
  * @result 0 if @a e does not belong to the current graph, 1 otherwise
  * @post new num_edges() == old num_edges() - 1 if preconditions are satisfies
  *
  * Note: invalidates removed edge
  * Complexity: O(@a e.node1().degree() * @a e.node2().degree())
  */
  size_type remove_edge(const Edge& e){
    // check the node belong to the graph
    if (!(has_edge(e.node1(),e.node2()))) {return 0;}

    size_type Edge_uid1 = -1; // only needed for assertion
    size_type Edge_uid2 = -1; // only needed for assertion

    // Scan adjacency list and remove edge uid
    for (size_type iEdge = 0; iEdge < e.node1().GetNode().AdjEdges_.size(); ++iEdge){
      // check uid
      if (e.node1().GetNode().AdjEdges_[iEdge] == e.uid_) {
        // mark for assertion
        Edge_uid1 = e.uid_;
        // swap uids in the adjacency list
        std::swap(e.node1().GetNode().AdjEdges_[iEdge],e.node1().GetNode().AdjEdges_[e.node1().GetNode().AdjEdges_.size()-1]);
        // eliminate current uid of the edge to be removed
        e.node1().GetNode().AdjEdges_.pop_back();
        break;
      }
    }
    // Scan adjacency list and remove edge uid
    for (size_type iEdge = 0; iEdge < e.node2().GetNode().AdjEdges_.size(); ++iEdge){
      // check uid
      if (e.node2().GetNode().AdjEdges_[iEdge] == e.uid_) {
        // mark for assertion
        Edge_uid2 = e.uid_;
        // swap uids in the adjacency list
        std::swap(e.node2().GetNode().AdjEdges_[iEdge],e.node2().GetNode().AdjEdges_[e.node2().GetNode().AdjEdges_.size()-1]);
        // eliminate current uid of the edge to be removed
        e.node2().GetNode().AdjEdges_.pop_back();
        break;
      }
    }
    // check you found the edge the edge
    assert((Edge_uid1 != (size_type)-1));
    assert((Edge_uid2 != (size_type)-1));
    // check you eliminated the same edge
    assert(Edge_uid1 == Edge_uid2);
    // check some internal consitency
    assert(Edges_[Edge_uid1].idx_ == e.index());
    // copy index
    auto e_index = e.index();
    // swap this edge index to the end and delete it
    std::swap(i2u_Edges_[e_index],i2u_Edges_[i2u_Edges_.size()-1]);
    i2u_Edges_.pop_back();
    Edges_[i2u_Edges_[e_index]].idx_ = e_index;
    // Invalidate Internal Representation of the Edge
    Edges_[Edge_uid1] = InternalEdge();
    // Success!
    return 1;
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

  using edge_value_type = E;
  class Edge : private totally_ordered<Edge>{
  public:
    /** Construct an invalid Edge. */
    Edge() {}

    /** Return a node of this Edge */
    Node node1() const {
      if (reversed_) return Node(Graph_,GetEdge().EndNode_);
      else return Node(Graph_,GetEdge().StartNode_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      if (reversed_) return Node(Graph_,GetEdge().StartNode_);
      else return Node(Graph_,GetEdge().EndNode_);
    }

    /** Test whether this edge and @a e are equal.
    *
    * Equal edges represent the same undirected edge between two nodes.
    */
    bool operator==(const Edge& e) const {
      if (this->Graph_ == e.Graph_ and this->node1() == e.node1() and this->node2() == e.node2())
      return true;
      else if (this->Graph_ == e.Graph_ and this->node1() == e.node2() and this->node2() == e.node1())
      return true;
      else return false;
    }

    /** Test whether this edge is less than @a e in a global order.
    *
    * This ordering function is useful for STL containers such as
    * std::map<>. It need not have any interpretive meaning.
    */
    bool operator<(const Edge& e) const {
      if (Graph_ == e.Graph_){return (GetEdge().idx_ < e.GetEdge().idx_);}
      else {return (Graph_ < e.Graph_);}
    }

    void flip_orientation(){
      if (reversed_ == true ) reversed_ = false;
      if (reversed_ == false) reversed_ = true ;
    }
    /** Return Edge Length*/
    double length() const {return GetEdge().Length_;}

    /** Return Edge Length*/
    double &length() {return GetEdge().Length_;}

    /** Return Edge Value*/
    edge_value_type& value() {return GetEdge().Value_;}

    /** Return Edge Value*/
    const edge_value_type& value() const {return GetEdge().Value_;}

    size_type index() const {return GetEdge().idx_;}

  private:

    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    graph_type* Graph_; // pointer to the graph
    size_type uid_;  // unique id_ of the graph
    bool reversed_; // true if it should return node2() when calling in node1() and viceversa;

    // Default Constructor
    Edge(graph_type* Graph, size_type uid, bool reversed = false)
    :  Graph_(Graph), uid_(uid), reversed_(reversed) {}

    InternalEdge& GetEdge(){
      assert(uid_ < Graph_->Edges_.size());
      return Graph_->Edges_[uid_];
    };

    const InternalEdge& GetEdge() const {
      assert(uid_ < Graph_->Edges_.size());
      return Graph_->Edges_[uid_];
    };
  };

  /** Return the total number of edges in the graph.
  *
  * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
  */
  size_type num_edges() const {
    return this->i2u_Edges_.size();
  }

  /** Return the edge with index @a i.
  * @pre 0 <= @a i < num_edges()
  *
  * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
  */
  Edge edge(size_type i) const {
    return Edge(const_cast<graph_type*>(this),i2u_Edges_[i]);        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
  * @pre @a a and @a b are valid nodes of this graph
  * @return True if for some @a i, edge(@a i) connects @a a and @a b.
  *
  * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
  */
  bool has_edge(const Node& a, const Node& b) const {
    for (auto iEdge = a.edge_begin(); iEdge != a.edge_end(); ++iEdge){
      if ((*iEdge).node2() == b) {return true;}
    }
    for (auto iEdge = b.edge_begin(); iEdge != b.edge_end(); ++iEdge){
      if ((*iEdge).node2() == a) {return true;}
    }
    return false;
  }

  Edge find_edge(const Node& a, const Node& b) {
    for (auto iEdge = a.edge_begin(); iEdge != a.edge_end(); ++iEdge){
      if ((*iEdge).node2() == b) {return (*iEdge);}
    }
    for (auto iEdge = b.edge_begin(); iEdge != b.edge_end(); ++iEdge){
      if ((*iEdge).node2() == a) {return (*iEdge);}
    }
    assert(1 != 0); // should not  reach this point
    return Edge(); // silence compiler warning
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {
    if (has_edge(a,b)){return find_edge(a,b);}
    //else
    assert(a.Graph_ == b.Graph_);
    auto NewEdge =  InternalEdge(a.uid_,b.uid_,value,0.0,num_edges(),value);
    Edges_.push_back(NewEdge);
    auto EdgeID = Edges_.size()-1;

    Nodes_[a.uid_].AdjEdges_.push_back(EdgeID);
    Nodes_[b.uid_].AdjEdges_.push_back(EdgeID);
    i2u_Edges_.push_back(EdgeID);
    return Edge(const_cast<graph_type*>(a.Graph_),EdgeID); // Correction for 0 index
  }
  /** Remove all nodes and edges from this graph.
  * @post num_nodes() == 0 && num_edges() == 0
  *
  * Invalidates all outstanding Node and Edge objects.
  */
  void clear() {
    Nodes_.clear();
    Edges_.clear();
    i2u_Nodes_.clear();
    i2u_Edges_.clear();
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

    /** Construct an invalid NodeIterator. */
    NodeIterator() {}

    /** Dereference operator
    * @pre valid iterator referring to a valid node
    * @return Node at the position of the iterator
    */
    Node operator*() const {
      assert(index_ < Graph_->i2u_Nodes_.size());
      return Node(const_cast<graph_type*>(Graph_),Graph_->i2u_Nodes_.at(index_)); }

      /** Increment operator
      * @pre valid iterator
      * @return an iterator referring to the next node
      */
      NodeIterator& operator++() {index_ += 1; return *this;}

      /** Equality comparison operator
      * param[in] NodeIterator _x_
      * @pre valid iterators
      * @return boolean true if iterators refer to the same node index of the same graph, false otherwise
      */
      bool operator==(const NodeIterator& x) const {return (index_ == x.index_ and Graph_ == x.Graph_);}

      /** Inequality comparison operator
      * param[in] NodeIterator _x_
      * @pre valid iterators
      * @return boolean false if iterators refer to the same node index of the same graph, true otherwise
      */
      bool operator!=(const NodeIterator& x) const {return !(index_ == x.index_ and Graph_ == x.Graph_);}
    private:
      friend class Graph;
      // Default Constructor
      NodeIterator(graph_type* Graph, size_type index)
      :  Graph_(Graph), index_(index) {}
      const graph_type* Graph_;// Parent Graph
      size_type index_; // node index corresponding to this iterator
    };

    /** Return node iterator pointing to the first node
    * @return iterator to first node
    */
    node_iterator node_begin() const {return NodeIterator(const_cast<graph_type*>(this),0);}

    /** Return node iterator pointing to the last node
    * @return iterator to last node
    */
    node_iterator node_end() const {return NodeIterator(const_cast<graph_type*>(this),i2u_Nodes_.size());}

    /** Remove node (and its incident edges) at the position pointed by the current node iterator
    *
    * @pre @a n_it is a valid node iterator whose dereference is an existing node
    * @post @a (*n_it) is deleted and new num_nodes() == old num_nodes()-1
    *       and new num_edges() == old num_edges()-n.degree() if preconditions are satisfied
    * @result an iterator pointing to the next element
    *
    * Note: invalidates deleted nodes and its adjiancent edges
    * Complexity O((*n_it).degree()*d) where d is the max degree in the graph
    */
    node_iterator remove_node(node_iterator n_it){
      assert(remove_node(*n_it) == 1);
      // After removal, nodes gets reindexed and n_it points
      // to the node whose index replaces the index of the deleted node.
      return n_it;
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

      /** Dereference operator, returns the edge corresponding to current iterator
      * @return Edge at iterator's position
      * @post (Edge.node1() == node) spawning the incident edges the iterator is looping over
      */
      Edge operator*() const {
        // Extract the undirected edge
        Edge EdgeToReturn = Edge(const_cast<graph_type*>(Graph_),Graph_->Nodes_[node_uid_].AdjEdges_[list_id_]);
        // make sure its the spawning node is numbered first
        // if (EdgeToReturn.node1().index() == node_id_) return EdgeToReturn;
        if (Graph_->i2u_Nodes_[EdgeToReturn.node1().index()] == node_uid_) return EdgeToReturn;
        else {
          EdgeToReturn.flip_orientation(); // reverse the orientation of the edge
          assert(Graph_->i2u_Nodes_[EdgeToReturn.node1().index()] == node_uid_); // make sure it actually worked
          return EdgeToReturn;
        }
      }

      /** Increment operator for the edge_iterator
      * @return Iterator with incremented position
      */
      IncidentIterator& operator++() {list_id_ += 1; return *this;}

      /** Equality comparison operator
      * @param[in] edge iterator
      * @return bool true if iterarors refer to same edge spawned by the same node of the same graph, false otherwise
      */
      bool operator==(const IncidentIterator& x) const {return (list_id_ == x.list_id_ and Graph_ == x.Graph_ and node_uid_ == x.node_uid_ );}

      /** Inequality comparison operator
      * @param[in] edge iterator
      * @return bool false if iterarors refer to same edge spawned by the same node of the same graph, true otherwise
      */
      bool operator!=(const IncidentIterator& x) const {return !(list_id_ == x.list_id_ and Graph_ == x.Graph_ and node_uid_ == x.node_uid_ );}

    private:
      friend class Graph;
      // Default constructor
      IncidentIterator(graph_type* Graph, size_type list_id_, size_type node_uid)
      :  Graph_(Graph), list_id_(list_id_), node_uid_(node_uid) {}

      const graph_type* Graph_; // parent Graph
      size_type list_id_; // position of the element in the _AdjEdges_ list
      size_type node_uid_; // Node spawning the adjiancent edges
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

      /** Construct an invalid EdgeIterator. */
      EdgeIterator() {}

      /** Dereference iterator for edge_iterator
      * @return edge pointed by the iterator
      */
      Edge operator*() const {return Edge(const_cast<graph_type*>(Graph_),Graph_->i2u_Edges_[index_]); }

      /** Increment iterator for edge_iterator
      * @return EdgeIterator incremented by one position
      */
      EdgeIterator& operator++() {index_++; return *this;}

      /** Equality comparator between EdgeIterator
      * @param[in] edge iterator
      * @return bool true if @a e and @a this refer to the edge of the same graph, false otherwise
      */
      bool operator==(const EdgeIterator& e) const {return (index_ == e.index_ and Graph_ == e.Graph_);}
      bool operator!=(const EdgeIterator& e) const {return !(index_ == e.index_ and Graph_ == e.Graph_);}

    private:
      friend class Graph;
      const graph_type* Graph_; // Parent Graph
      size_type index_; // Edge index
      // Default Constructor
      EdgeIterator(graph_type* g, size_type index)
      :  Graph_(g), index_(index) {}
    };

    /** Edge iterator at first edge of the graph
    * @return iterator at first edge
    */
    edge_iterator edge_begin() const {return EdgeIterator(const_cast<graph_type*>(this),0);}

    /** Edge iterator at last edge of the graph
    * @return iterator at last edge
    */
    edge_iterator edge_end()   const {return EdgeIterator(const_cast<graph_type*>(this),num_edges());}

    /** Remove edge at the position pointed by the current edge iterator
    *
    * @pre @a e_it is a valid edge iterator whose dereference is an existing edge
    * @post @a (*e_it) is deleted and new num_edges() == old num_edges()-1
    * @result an iterator pointing to the next element
    *
    * Note: invalidates deleted edge
    * Complexity O((*n_it).node1().degree()*(*n_it).node2().degree())
    */
    edge_iterator remove_edge(edge_iterator e_it){
      assert(remove_edge(*e_it) == 1);
      // After removal, edges gets reindexed and this point
      // to the edge whose index replaces the index of the deleted edge.
      return e_it;
    }

    struct InternalNode{
      Point Position_;
      std::vector<size_type> AdjEdges_;
      node_value_type Value_;
      size_type idx_;
      bool isActive_;
      InternalNode() {}
      InternalNode(Point Position,std::vector<size_type> AdjEdges,node_value_type Value, size_type idx, bool isActive) :
      Position_(Position), AdjEdges_(AdjEdges), Value_(Value), idx_(idx), isActive_(isActive) {}
    };

    struct InternalEdge{
      size_type StartNode_;
      size_type EndNode_;
      edge_value_type EdgeValues_;
      double Length_;
      size_type idx_;
      edge_value_type Value_;
      InternalEdge() {};
      InternalEdge(size_type StartNode, size_type EndNode, edge_value_type EdgeValues, double Length, size_type idx, edge_value_type Value) :
      StartNode_(StartNode), EndNode_(EndNode), EdgeValues_(EdgeValues), Length_(Length), idx_(idx), Value_(Value) {}
    };

  private:

    std::vector<InternalNode> Nodes_;
    std::vector<InternalEdge> Edges_;
    std::vector<size_type>    i2u_Nodes_;
    std::vector<size_type>    i2u_Edges_;
  };


  #endif // CME212_GRAPH_HPP
