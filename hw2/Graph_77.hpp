#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <iostream>
#include <map>
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief  A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Type of node value */
  using node_value_type = V;

  /** Type of edges value */
  using edge_value_type = E;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  //for an empty graph, we want the nodes_ and edges_ vectors to be empty
  Graph() 
      : nodes_(std::vector<internal_node>(0)),
	edges_(std::vector<internal_edge>(0)) {
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
    //invalid node constructor
    Node() {
    }

    /** Return reference to node's position, making it modifiable */
    Point& position() {
      return graph_->nodes_[node_id_].point;
    }

    /** Return this node's position. */
    const Point& position() const {
      //return the point field of the internal_node struct at the node index
      return graph_->nodes_[node_id_].point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return node_id_;
    }
    
    /** Return the value associated with this node */
    //return by reference: used for setting value
    node_value_type& value() {
      return graph_->nodes_[node_id_].val;
    }

    /** Return the value associated with this node */
    //return by constant reference: used for getting value
    const node_value_type& value() const {
      return graph_->nodes_[node_id_].val;
    }

    /** @brief Determines  the value of the degree of the node (ie. number of 
     * incident edges).
     *
     * @return the number of undirected edges incident to the node object
     */
    size_type degree() const {
      IncidentIterator ii(graph_,node_id_);
      return ii.edge_idx_list_.size();   
    }
    
    /** @brief Creates an iterator to the first adjacent edge to the node.
     *
     * @return incident_iterator iterator to first adjacent
     * edge of the node object.
     * @post If the node has no incident edges, returns nullptr
     */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_,node_id_);
    }
    
    /** @brif Creates an iterator to one past the last adjacent edge to
     * the node.
     *
     * @return incident_iterator iterator to one past the last adjacent
     * edge of the node object.
     * @post returned iterator is equavalent to nullptr (one past the last
     * is null)
     */
    incident_iterator edge_end() const{
      IncidentIterator ii(graph_,node_id_);
      ii.counter_ = ii.edge_idx_list_.size();
      return ii;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      //return true if nodes belong to same graph and have same index
      if (node_id_ == n.node_id_ and graph_ == n.graph_) return true;
      (void) n;          // Quiet compiler warning
      //return false otherwise
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
      /** Return true if the index of this node is less than that of n
       * and they belong to the same graph
       */
      if (node_id_ < n.node_id_ and graph_ == n.graph_) return true;
      else if (graph_ < n.graph_) return true;
      (void) n;           // Quiet compiler warning
      //return false otherwise
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    //pointer to the location of the Graph object of which Node is a node
    Graph* graph_;
    //current index of the node in the node_ vector 
    size_type node_id_;

    //valid Node constructor
    Node(const Graph* g, size_type idx) 
	: graph_(const_cast<Graph*>(g)), node_id_(idx) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    //return the size of the vector holding all nodes
    return nodes_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size(); //reuse the size() member function
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, 
                const node_value_type& value = node_value_type()) {
    //declare noew internal_node struct
    internal_node n;
    //initialize member fields
    n.node_id = num_nodes()-1;
    n.point = position;
    n.val = value;
    nodes_.push_back(n); //add to nodes_ vector
    return Node(this, num_nodes()-1); //return node object
    (void) position;      // Quiet compiler warning
    return Node();
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    //if n belongs to the same graph and has index less than # of graph nodes
    if (n.node_id_ < num_nodes() and n.graph_ == this) return true;
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
    //return a Node constructor for this graph object and specified point 
    return Node(this, i);
    (void) i;             // Quiet compiler warning
    return Node();        // Invalid node
  }

  //remove node from graph, as well as all edges connected to that node
  //returns 1 if successful, 0 if edge not present initially
  size_type remove_node(const Node& n) {
    if (!has_node(n)) return 0;
    else {
      nodes_.erase(nodes_.begin() + n.index());
      //update node_ids of all elements after removed element in nodes_
      for (auto it = nodes_.begin() + n.index(); it != nodes_.end(); ++it){
        (*it).node_id -= 1;
      }
      //remove all edges which are connected to that node
      for (auto it = edge_begin(); it != edge_end();){
	if ((*it).node1().index() == n.index() || (*it).node2().index() == n.index()) {
          it = remove_edge(it);
        }
        else ++it;
      }
      //updae the node1 and node2 member variables of internal edge
      for (auto it = edges_.begin(); it != edges_.end(); ++it) {
        if ((*it).node1.index() > n.index()) (*it).node1.node_id_ -= 1;
        if ((*it).node2.index() > n.index()) (*it).node2.node_id_ -= 1;
      }
    }
    return 1;
  }

  node_iterator remove_node(node_iterator ni){
    size_type old_node_id = (*ni).index();
    remove_node(*ni);
    //update iterator location
    /** since elements were shifted in place and ids updated, next element
     * now has same id as old_node_id: return iterator to that edge
     */
    return NodeIterator(this,old_node_id);
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
    //invalid Edge constructor
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      //return node1 object stored in internal_node struct
      return graph_->edges_[index()].node1;
      return Node();	  //Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      //return node2 object stored in internal_node struct
      return graph_->edges_[index()].node2;
      return Node();      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //if edges in the same graph and edge ids are the same 
      if (edge_id_ == e.edge_id_ and graph_ == e.graph_) return true; 
      (void) e;           // Quiet compiler warning
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //we will define Edge e1 < Edge e2 if e1.edge_id_ < e2.edge_id_
      if (edge_id_ < e.edge_id_ and graph_ == e.graph_) return true;
      (void) e;           // Quiet compiler warning
      return false;
    }

    /** Return the index that the edge_id_ maps to; one to one for now **/
    size_type index() const {
      return edge_id_;
      //return graph_->edge_indices_.at(edge_id_);
    }

    /** Returns the euclidean length of the edge */
    double length() const {
      return norm(node1().position() - node2().position());
    }

    /** Return the value associated with this edge */
    //return by reference: used for setting value
    edge_value_type& value() {
      return graph_->edges_[edge_id_].val;
    }

    /** Return the value associated with this edge */
    //return by constant reference: used for getting value
    const edge_value_type& value() const {
      return graph_->edges_[edge_id_].val;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    //pointer to Graph object to which the edge belongs
    Graph* graph_;
    //unique id of the edge 
    size_type edge_id_;

    //valid Edge contructor
    Edge(const Graph* g, size_type id)
	: graph_(const_cast<Graph*>(g)),
	  edge_id_(id) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    //return edge object constructed from current graph and the edge id
    return Edge(this,edges_[i].edge_id);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return edge(@ i).index() if for some @a i, edge(@a i) connects @a a and @a b
   *         -1 if no edge connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  int has_edge(const Node& a, const Node& b) const {
    int idx = -1;
    //if nodes belong to different graphs return false
    if (a.graph_ != b.graph_) return idx;
    auto p1 = std::make_pair(a,b);
    auto p2 = std::make_pair(b,a);
    //check whether either pair maps to an index in the nodes2edge_id_ map
    //if yes, set idx to equal the mapped index
    if (nodes2edge_id_.count(p1) == 1) idx = nodes2edge_id_.at(p1);
    else if (nodes2edge_id_.count(p2) == 1) idx = nodes2edge_id_.at(p2);
    (void) a; (void) b;   // Quiet compiler warning
    //otherwise edge does not exist
    return idx;
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
    //test if edge exists and if the edge does not aleady exist
    int test = has_edge(a,b);
    if (test == -1) {
      internal_edge e; //create new internal_edge struct
      //initialize struct member fields
      e.edge_id = edges_.size();
      e.node1 = a;
      e.node2 = b;
      //add internal_edge to the edges_ vector
      edges_.push_back(e);
      //make a pair consisting of nodes to be the key in node2edge_id_
      auto p = std::make_pair(a,b);
      //add to node2edge_id_ map
      nodes2edge_id_[p] = e.edge_id;
      //return edge object
      return Edge(this,e.edge_id);
    }
    //otherwise return edge that already exists
    //has_edge returns the edge index if it does already exist, so test = edge_id
    else {
       return Edge(this,test);
    }
    (void) a, (void) b;   // Quiet compiler warning
    return Edge();        // Invalid Edge
  }
  
  //removes edge connecting node 1 and node 2 from graph
  //returns 1 if successful, 0 if edge not present initially
  size_type remove_edge(const Node& n1, const Node& n2) {
    int idx = has_edge(n1,n2);
    if (idx == -1) return 0;
    else {
      unsigned int id = unsigned(idx);
      edges_.erase(edges_.begin()+id);
      //erase edge from node2edge map   
      for (auto it = nodes2edge_id_.begin(); it != nodes2edge_id_.end(); it++){
        if ((*it).second == id) {
          nodes2edge_id_.erase(it);
          break;
        }
      }
      //update edge_ids of all elements after removed element in edges_
      for (auto it = edges_.begin() + id; it != edges_.end(); ++it) {
        (*it).edge_id -= 1;
      }
      //make new map from node pair to edge id
      nodes2edge_id_.clear();
      for (auto it = edges_.begin(); it != edges_.end(); ++it) {
        auto p = std::make_pair((*it).node1,(*it).node2);
        nodes2edge_id_[p] = (*it).edge_id;
      }
    }
    return 1;
  }

  //removes edge (function argument) from graph
  //returns 1 if successful, 0 if edge not present initially
  size_type remove_edge(const Edge& e){
    //use previous remove_edge function (input two nodes)
    return remove_edge(e.node1(),e.node2());
  }
  
  //removes edge pointed to by the edge iterator
  //returns edge_iterator to edge with next edge_id
  edge_iterator remove_edge(edge_iterator ei){
    size_type old_edge_id = (*ei).index();
    remove_edge(*ei);
    //update iterator location
    /** since elements were shifted in place and ids updated, next element
     * now has same id as old_edge_id: return iterator to that edge
     */
    return EdgeIterator(this,old_edge_id);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    //make use of vector clear() function
    nodes_.clear();
    edges_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() : graph_(nullptr), idx_(0) {
    }

    /** @brief dereferencing operator defined to return a node
     * object of the specified node index. 
     *
     * @return Node object with the index given by NodeIterator.idx_
     * @post the index of the returned node < graph.size()
     */
    Node operator*() const {return Node(graph_,idx_);}
    
    /** @brief Operator to increment the NodeIterator
     *
     * @return The same NodeIterator object with the idx_ field incremented by 1
     * @post idx_ < graph.size()
     */
    NodeIterator& operator++(){
      idx_++;
      return *this;
    }

    NodeIterator& operator--(){
      idx_--;
      return *this;
    }
    
    /** @brief Binary operator to test equality of two node iterators
     *
     * @return True if node iterators are equal, False otherwise
     * @pre Node iterators are equal if they belong to the same graph
     *      and have their idx_ field is equal.
     */
    bool operator==(const NodeIterator& ni) const{
      return graph_ == ni.graph_ and idx_ == ni.idx_;
    }

   private:
    //graph can access fields
    friend class Graph;
    Graph* graph_; //pointer to graph to which node iterator belongs
    size_type idx_; //node index currently pointing to
    //valid constructor
    NodeIterator(const Graph* g, size_type idx)
                : graph_(const_cast<Graph*>(g)),
                  idx_(idx) {
    }
  };

  /** @brief Creates a node iterator to the first node in the graph object
   *
   * @return node_iterator object to the first node in the graph
   * @pre the first node in the graph has index 0
   * @pre all nodes in the graph have a unique index
   */
  node_iterator node_begin() const{
    return NodeIterator(this,0);
  }
  
  /** @brief Creates a node iterator to one past the last node in the graph
   *
   * @return node_iterator object to one past the last node
   * @post Returned node iterator should have index = graph.num_nodes()
   * @post There is no node with returned node iterator index
   */ 
  node_iterator node_end() const{
    return NodeIterator(this,num_nodes());
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
    IncidentIterator() : graph_(nullptr), node_id_(0), edge_idx_list_(0){
    }

    /** @brief dereferencing operator defined to return an Edge
     * object which is incident to the specified node
     *
     * @return Edge object with the index given in edge_idx_list_ member field
     * @post the index of the returned node < graph.num_edges()
     */
    Edge operator*() const{
      size_type edge_idx = edge_idx_list_[counter_];
      return Edge(graph_,graph_->edges_[edge_idx].edge_id);
    }

    /** @brief Operator to increment the IncidentIterator
     *
     * @return The same IncidentIterator object with the counter_ field 
     * incremented by 1
     * @post counter_ < node_id_.degree()
     */    
    IncidentIterator& operator++(){
      counter_++;
      return *this;      
    }

    /** @brief Binary operator to test equality of two incident iterators
     *
     * @return True if incident iterators are equal, False otherwise
     * @pre Incident iterators are equal if they belong to the same graph
     *      and have their node_id_ and counter_ fields equal.
     */
    bool operator==(const IncidentIterator& ii) const{
      return graph_ == ii.graph_ and node_id_ == ii.node_id_
             and counter_ == ii.counter_;
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type node_id_;
    size_type counter_;
    //vector to hold indices of incident edges to node object
    std::vector<int> edge_idx_list_;

    IncidentIterator(const Graph* g, size_type nid)
                    : graph_(const_cast<Graph*>(g)),
                      node_id_(nid), counter_(0) {
      edge_idx_list_ = {}; //initialize empty list for edge indices
      //loop through all graph edges
      for (size_type i = 0; i < graph_->edges_.size(); ++i){
        //if node1 field of edge matches node object, add index to list
        if (graph_->edges_[i].node1.node_id_ == nid)
          edge_idx_list_.push_back(i);
        //if node2 field of edge matches, switch node1 and node2 and add to lits
        else if (graph_->edges_[i].node2.node_id_ == nid)
          edge_idx_list_.push_back(i);
          std::swap(graph_->edges_[i].node1,graph_->edges_[i].node2);
      }
    }
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
    EdgeIterator() : graph_(nullptr), edge_id_(0){
    }

    /** @brief dereferencing operator defined to return an Edge
     * object to which the iterator points.
     *
     * @return Edge object with the index given in edge_id__ member field
     * @post the index of the returned node < graph.num_edges()
     */
    Edge operator*() const{
      return Edge(graph_,edge_id_);  
    }

    /** @brief Operator to increment the EdgeIterator
     *
     * @return The same EdgeIterator object with the edge_id_ field 
     * incremented by 1
     * @post edge_id_ < graph.num_edges()
     */
    EdgeIterator& operator++(){
      edge_id_++;
      return *this;
    }

    EdgeIterator& operator--(){
      edge_id_--;
      return *this;
    }

    /** @brief Binary operator to test equality of two edge iterators
     *
     * @return True if edge iterators are equal, False otherwise
     * @pre edge iterators are equal if they belong to the same graph
     *      and have their edge_id_ fields equal.
     */
    bool operator==(const EdgeIterator& ei) const{
      return graph_ == ei.graph_ and edge_id_ == ei.edge_id_;
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type edge_id_;
    //valid edge iterator constructor
    EdgeIterator(const Graph* g, size_type id)
                : graph_(const_cast<Graph*>(g)),
                  edge_id_(id) {
    }
  };

  /** @brief Creates a edge iterator to the first edge in the graph object
   *
   * @return edge_iterator object to the first edge in the graph
   * @pre the first edge in the graph has id 0
   * @pre all edge in the graph have a unique index
   */  
  edge_iterator edge_begin() const{
  return EdgeIterator(this,0);
}

  /** @brief Creates a edge iterator to one past the last edge in the graph
   *
   * @return edge_iterator object to one past the last edge
   * @post Returned edge iterator should have edge_id_ = graph.num_edge()
   * @post There is no edge with id of  returned edge iterator edge_id_
   */ 
  edge_iterator edge_end() const{
  return EdgeIterator(this,num_edges());
}

 private:
   /** struct for internal nodes */
   struct internal_node {
     size_type node_id; //id of node
     Point point; //associated point object specifying location
     node_value_type val; //value of node
   };

   /** stuct for internal endges */
   struct internal_edge {
     size_type edge_id; //id of idge
     node_type node1; //incident node 1
     node_type node2; //incident node 2
     edge_value_type val; //value of edge
   };
   /** Vector of internal nodes in graph */
   std::vector<internal_node> nodes_;

   /** Vector of internal_edges */
   std::vector<internal_edge> edges_;

   /** Map which maps node pair to edge id */
   std::map<std::pair<node_type,node_type>,size_type> nodes2edge_id_;
};

#endif // CME212_GRAPH_HPP
