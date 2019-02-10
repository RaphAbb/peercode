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

  /** Synonym for Node Value Type */
  using node_value_type = V;
  
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

  /** Type of pairs that designate adjacent edges
      first int is the destination node's index
      second int is index of corresponding edge */
  using pair_node_index_type = std::pair<size_type, size_type>;  
  
  /** Type of pairs that designate edges
      first int is the one node's index
      second int is the other node's index */
  using pair_nodes_type = std::pair<size_type, size_type>;  
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
      //All take default empty values
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
  class Node: private totally_ordered<Node> {
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
    Node(){
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      //return *pos;
      return graph_->points_[id_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return id_;
    }
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      (void) n;          // Quiet compiler warning
      return (graph_==n.graph_)&&(id_==n.id_);
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
      (void) n;           // Quiet compiler warning
      return id_<n.id_;
    }
   

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Return the value of this node.
     * @return the value of _this_   
     * @pre the node is valid
     *
     * Complexity: O(1).
     */
   node_value_type& value(){return graph_->values_[id_];}
   
    /** Return the value of this node.
     * @return the value of _this_   
     * @pre the node is valid
     *
     * Complexity: O(1).
     */
   const node_value_type& value() const {return const_cast<node_value_type&>(graph_->values_[id_]);}
    
   /** Return the beginning of an iterator to adjacent edges.
     * @return an incident_iterator of adjacent edges
     * @pre the node is valid
     *
     * Complexity: O(1).
     */
   incident_iterator edge_begin(){return incident_iterator(id_, 0, graph_);}
   
   /** Return the end iterator to adjacent edges.
     * @return an incident_iterator of adjacent edges
     * @pre the node is valid
     * @post _result_ is the end of _edge_begin()_  
     * Complexity: O(1).
     */
   incident_iterator edge_end(){return incident_iterator(id_, degree(), graph_);}
   
   /** Return the degree of this node.
     * @return number of adjacent edges
     * @pre the node is valid
     *  
     * Complexity: O(1).
     */
   size_type degree() const {return graph_->neigh_list_.at(id_).size();}
   
   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    size_type id_;
    graph_type* graph_;

    Node(const graph_type* g, size_type idx)
      : id_(idx), graph_(const_cast<graph_type*>(g)) {
    } 
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return this->n_nodes;
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
  Node add_node(const Point& position, const node_value_type& v = node_value_type()) {
    // HW0: YOUR CODE HERE
    Node n(this, n_nodes);
    points_.push_back(position);
    values_.push_back(v);
    std::vector<pair_node_index_type> emptyn;
    neigh_list_.push_back(emptyn);
    this->n_nodes++;
    //std::cout<<"n= "<<n_nodes<<std::endl;
    //std::cout<<"1adding "<<n.position()<<std::endl;
    //std::cout<<"2adding "<<node(n_nodes-1).position()<<std::endl;
    return n;        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    (void) n;            // Quiet compiler warning
    return (this==n.graph);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    (void) i;             // Quiet compiler warning
    //std::cout<<"Accessing node "<<i<<std::endl;
    //std::cout<<" at "<<points_[i]<<std::endl;
    return Node(this, i);        // Invalid node
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node1_);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node2_);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      bool orient1 = (node1_==e.node1_)&&(node2_==e.node2_);
      bool orient2 = (node2_==e.node1_)&&(node1_==e.node2_);
      return (orient1 || orient2);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      return id_<e.id_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    size_type node1_;
    size_type node2_;
    size_type id_;
    graph_type* graph_;

    /** Construct valid edge object */
    Edge(size_type node1, size_type node2, size_type idx, const graph_type* g)
      : node1_(node1), node2_(node2), id_(idx), graph_(const_cast<graph_type*>(g)) {

    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return n_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    (void) i;             // Quiet compiler warning
    return Edge(edges_.at(i).first, edges_.at(i).second, i, this);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    (void) a; (void) b;   // Quiet compiler warning
    for(auto c : neigh_list_.at(a.index())){
        if (c.first == b.index()){
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
    for(auto c : neigh_list_.at(a.index())){
        if (c.first == b.index()){
          return Edge(a.index(), b.index(), c.second, this);
        }
    }
    //Else create new edge
    Edge e = Edge(a.index(), b.index(), n_edges, this);;
     
    pair_nodes_type p(a.index(), b.index());
    edges_.push_back(p);
    pair_node_index_type p1(b.index(), n_edges);
    pair_node_index_type p2(a.index(), n_edges);
    neigh_list_.at(a.index()).push_back(p1);
    neigh_list_.at(b.index()).push_back(p2);
    this->n_edges++;
    (void) a, (void) b;   // Quiet compiler warning
    return e;        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    n_nodes = 0;
    values_.clear();
    neigh_list_.clear();
    points_.clear();
    edges_.clear();
    n_edges = 0;
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
        index_=0;
    }
    
    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const
   
   /** Return the current element of this iterator.
     * @return current Node
     * @pre the NodeIterator is valid
     * @pre @ index_ < @size()
     * @post _result_ is a valid node
     * Complexity: O(1).
     */
    value_type operator*() const{
        return Node(graph_, index_);
    }

   /** Increments iterator.
     * @return reference to new iterator
     * @pre the NodeIterator is valid
     * @pre @ index_ < @size()
     * Complexity: O(1).
     */
    NodeIterator& operator++(){
        index_++;
        return *this;
    }

   /** Tests equality of iterator.
     * @return boolean of equality
     * @pre this NodeIterator and @n_i are valid
     * @pre this NodeIterator and @n_i belong to the same graph 
     * Complexity: O(1).
     */
    bool operator==(const NodeIterator& n_i) const{
         return (index_==n_i.index_);
    }


   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    size_type index_;
    graph_type* graph_;
    NodeIterator(size_type id, const graph_type* g): index_(id), graph_(const_cast<graph_type*>(g)){}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  
   /** Return the beginning of an iterator to all nodes.
     * @return an iterator pointing to first node in graph 
     *
     * Complexity: O(1).
     */

   NodeIterator node_begin() const{
       return NodeIterator(0, this);
  }
   /** Return the end of an iterator to all nodes.
     * @return an iterator pointing to one past last node in graph 
     * @post when enough incremented @node_begin() == @node_end()
     * Complexity: O(1).
     */
  NodeIterator node_end() const{
       return NodeIterator(n_nodes, this); 
  }


  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator> {
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
   
    /** Return the current element of this iterator.
     * @return current Edge
     * @pre the node_ and IncidentIterator are valid
     * @pre @ index_ < @node_.degree()
     * @post _result_ is a valid edge
     * Complexity: O(1).
     */
    Edge operator*() const{
         std::pair<size_type, size_type> p = graph_->neigh_list_.at(node_).at(index_);
         return Edge(node_, p.first, p.second, graph_);
    }
   
    /** Increments iterator.
     * @return reference to new iterator
     * @pre the IncidentIterator is valid
     * @pre @ index_ < @node_.degree()
     * Complexity: O(1).
     */
    IncidentIterator& operator++(){
         index_++;
         return *this;
    }
   /** Tests equality of iterator.
     * @return boolean of equality
     * @pre this IncidentIterator and @iit are valid
     * @pre this IncidentIterator and @iit belong to the same graph 
     * Complexity: O(1).
     */
    bool operator==(const IncidentIterator& iit) const{
         return (node_==iit.node_)&& (index_ == iit.index_);
   }
   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    size_type node_;
    size_type index_;
    graph_type* graph_;
    IncidentIterator(size_type node, size_type index, const graph_type* g): node_(node), index_(index), graph_(const_cast<graph_type*>(g)){}
    
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator> {
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

    /** Return the current element of this iterator.
     * @return current Edge
     * @pre the EdgeIterator is valid
     * @pre @ index_ < @graph_->num_edges()
     * @post _result_ is a valid edge
     * Complexity: O(1).
     */
    Edge operator*() const{
        return Edge(graph_->edges_[index_].first, graph_->edges_[index_].second, index_, graph_);
    }
    
    /** Increments iterator.
     * @return reference to new iterator
     * @pre the EdgeIterator is valid
     * Complexity: O(1).
     */
    EdgeIterator& operator++(){
        index_++;
        return *this;
    }
     
   /** Tests equality of iterator.
     * @return boolean of equality
     * @pre this EdgeIterator and @ei are valid
     * @pre this EdgeIterator and @ei belong to the same graph 
     * Complexity: O(1).
     */
    bool operator==(const EdgeIterator& ei) const{
        return index_==ei.index_;
    }
    
   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    size_type index_;
    graph_type* graph_;
    EdgeIterator(size_type index, const graph_type* g): index_(index), graph_(const_cast<graph_type*>(g)){}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

   /** Return  an iterator to all edges.
     * @return an iterator pointing to one past last edge in Graph 
     * Complexity: O(1).
     */
  edge_iterator edge_begin() const{
     return EdgeIterator(0, this);
  }
  
   /** Return the end of an iterator to all edges.
     * @return an iterator pointing to one past last edge in graph 
     * @post when enough incremented @edge_begin() == @edge_end()
     * Complexity: O(1).
     */
  edge_iterator edge_end() const{
     return EdgeIterator(n_edges, this);
  }
 private:
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
     std::vector<Point> points_; 
     std::vector<node_value_type> values_; 
     std::vector<pair_nodes_type> edges_;
     std::vector<std::vector<pair_node_index_type>> neigh_list_; 
     
     size_type n_nodes;
     size_type n_edges;

};

#endif // CME212_GRAPH_HPP

