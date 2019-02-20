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
//template <typename V>
template <typename V, typename E>
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
  /** Synonym for Node Value Type */
  using edge_value_type = E;

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

    /** Return this node's position. */
    Point& position(){
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
      return (graph_<n.graph_)||(id_<n.id_);
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
   incident_iterator edge_begin() const {return incident_iterator(id_, 0, graph_);}
   
   /** Return the end iterator to adjacent edges.
     * @return an incident_iterator of adjacent edges
     * @pre the node is valid
     * @post _result_ is the end of _edge_begin()_  
     * Complexity: O(1).
     */
   incident_iterator edge_end() const {return incident_iterator(id_, degree(), graph_);}
   
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
    return n;        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
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

    size_type index() const{
      return id_;
    }

    double length() const{
      return graph_->lengths_[id_];
    }

    edge_value_type& value(){
      return graph_->edge_values_[id_];
    }
    const edge_value_type& value() const{
      return const_cast<edge_value_type&>(graph_->edge_values_[id_]);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (graph_==e.graph_)&&(e.id_ == id_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (graph_<e.graph_)||(id_<e.id_);
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
    
    //Adding length
    lengths_.push_back(norm(points_[a.index()]-points_[b.index()]));
    //Adding value
    edge_values_.push_back(100);
    //Adding edge
    pair_nodes_type p(a.index(), b.index());
    edges_.push_back(p);
    pair_node_index_type p1(b.index(), n_edges);
    pair_node_index_type p2(a.index(), n_edges);
    neigh_list_.at(a.index()).push_back(p1);
    neigh_list_.at(b.index()).push_back(p2);
    this->n_edges++;
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
         pair_node_index_type p = graph_->neigh_list_.at(node_).at(index_);
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
  // Removing

  /** Deletes an edge from the graph.
   
   * @param[in] e edge to delete
   * @return 0 if no exception
   * @pre the Edge is valid and in the graph
   * @post new _n_edges_ = old _n_edges_ - 1
   * @post if _e_ = (u, v) new v.degree() = old v.degree() - 1
   *                       new u.degree() = old u.degree() - 1
   * @post new edge(_e_.index()) = old edge(n_edges - 1)
   */
  size_type remove_edge(const Edge& e){
      //Delete info 
      edges_[e.index()] = edges_[n_edges-1];
      edges_.pop_back();
      lengths_[e.index()] = lengths_[n_edges-1];
      lengths_.pop_back();
      edge_values_[e.index()] = edge_values_[n_edges-1];
      edge_values_.pop_back();

      // As we changed the index of the previously last node 
      // we need to update its index
 
      //Update edge index in n1 incident edges
      std::vector<pair_node_index_type>& edges_1 = neigh_list_[edges_[e.index()].first];
      for (pair_node_index_type& pair:edges_1){
         if (pair.second == n_edges-1){
            pair.second = e.index();
            break;
         }
      }

      //Update edge index in n2 incident edges
      std::vector<pair_node_index_type>& edges_2 = neigh_list_[edges_[e.index()].second];
      for (pair_node_index_type& pair: edges_2){
         if (pair.second == n_edges-1){
            pair.second = e.index();
            break;
         }
      }

      //Delete edge in n1 incident edges
      std::vector<pair_node_index_type>& edges1 = neigh_list_[e.node1().index()];
      int count = 0;
      for (auto pair:edges1){
         if (pair.second == e.index()){
            break;
         }
         count+=1;
      }
      edges1[count] = edges1.back(); 
      edges1.pop_back();
      
      //Delete edge in n2 incident edges
      std::vector<pair_node_index_type>& edges2 = neigh_list_[e.node2().index()];
      count = 0;
      for (auto pair:edges2){
         if (pair.second == e.index()){
            break;
         }
         count+=1;
      }
      edges2[count] = edges2.back(); 
      edges2.pop_back();

      n_edges -= 1;
      return 0;
  }

  
  /** Delete an node and incident edges from the graph.
   * @param[in] n node to delete
   * @return 0 if no exception
   * @pre the Node is valid and in the graph
   * @post new _n_nodes_ = old _n_nodes_ - 1
   * @post new _n_edges_ = old _n_edges_ - _n_.degree()
   * @post new node(_n_.index()) = old node(n_nodes - 1)
   */
  size_type remove_node(const Node& n){
      // Delete incident edges
      auto eit = n.edge_begin();
      while (eit != n.edge_end()){remove_edge(*eit);}
      //Delete node info
      points_[n.index()] = points_[n_nodes-1];
      points_.pop_back();
      values_[n.index()] = values_[n_nodes-1];
      values_.pop_back();
      neigh_list_[n.index()] = neigh_list_[n_nodes-1];
      neigh_list_.pop_back();
      
      // update index of node id=(n_nodes - 1) in edges_ and neigh_list_
      // now n is the last node
      for (pair_node_index_type& pair : neigh_list_[n.index()]){
          for (pair_node_index_type& neigh : neigh_list_[pair.first]){
              if (neigh.first == n_nodes - 1)
                  neigh.first = n.index();
          }
          pair_nodes_type& edges = edges_[pair.second];
          if (edges.first == n_nodes - 1)
              edges.first = n.index();
          else 
              edges.second = n.index();
      }
      
      n_nodes -= 1;
      return 0;
  }
  
  /** Erase an edge from the edge_iterator.
   *
   * @param[in] e_it  	edge_iterator to erase from
   * @return new edge_iterator 
   *
   * @pre the EdgeIterator is valid and not empty
   * @post  if (++_e_it_ != end())
   *           *_result_ == (*(++_e_it_)
   *        else _result_ == end()
   */
  edge_iterator remove_edge(edge_iterator e_it){
     remove_edge(*e_it);
     return e_it;
  }

  /** Erase an node from the node_iterator.
   *
   * @param[in] n_it  	node_iterator to erase from
   * @return new node_iterator 
   *
   * @pre the NodeIterator is valid and not empty
   * @post  if (++_n_it_ != end())
   *           *_result_ == (*(++_n_it_)
   *        else _result_ == end()
   */
  node_iterator remove_node(node_iterator n_it){
     remove_node(*n_it);
     return n_it;
  }

  /** Delete an edge from the graph, if present.
   
   * @param[in] n1 first incident node
   * @param[in] n2 first incident node
   * @return 0 if no such edge or no exception
 
   * @pre the Node are valid and in the graph
   * @post if (n1, n2) in graph
   *          new _n_edges_ = old _n_edges_ - 1
   *          new n1.degree() = old n1.degree() - 1
   *          new n2.degree() = old n2.degree() - 1
   */
  size_type remove_edge(const Node& n1, const Node& n2){
    if (!has_edge(n1, n2)){
	return 0; 
     }
     std::vector<pair_node_index_type>& edges1 = neigh_list_[n1.index()];
     size_type e_id = 0;
     for (auto pair:edges1){
        if (pair.first == n2.index()){
           // Retrieve edge index
           e_id = pair.second;
           break;
        }
     }
     std::cout<<std::endl; //See Piazza
     remove_edge(Edge(n1.index(), n2.index(), e_id, this)); 
     return 0;
  }

 private:
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
     std::vector<Point> points_; 
     std::vector<node_value_type> values_; 
     std::vector<double> lengths_; 
     std::vector<edge_value_type> edge_values_; 
     std::vector<pair_nodes_type> edges_;
     std::vector<std::vector<pair_node_index_type>> neigh_list_; 
     
     size_type n_nodes;
     size_type n_edges;

};

#endif // CME212_GRAPH_HPP

