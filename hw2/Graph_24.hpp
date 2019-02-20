#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
#include <iostream>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = double,typename E = double>
class Graph {
 private:
  // HW0: YOUR CODE HERE
  struct internal_node;
  struct internal_edge;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  using node_value_type = V; // typedef V node_value_type;
  using edge_value_type = E; // typedef E edge_value_type;

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
  Graph() : 
    nodes_(), 
    idx_to_uid(),
    adjacency_list(),
    n_edges(0)
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
    Node() : graph_(nullptr), uid_(size_type(-1)) {}

    /** Return this node's position. */
    const Point& position() const {
      size_type i = index();
      return graph_->nodes_[i].position;
    }

    /** Return a reference to this node's position */
    Point& position() {
      size_type i = index();
      return graph_->nodes_[i].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // return graph_->idx_to_uid[uid_];
      if (uid_ > graph_->size()){
        return size_type(-1);
      }
      return graph_->nodes_[uid_].idx;
    }

    /** Return a reference to this node's value **/
    node_value_type& value() {
      // return graph_->nodes_[index()].val;
      return graph_->nodes_[uid_].val;
    }

    /** Return a const reference to this node's value **/
    const node_value_type& value() const {
      // return graph_->nodes_[index()].val;
      return graph_->nodes_[uid_].val;
    }

    /** Return the degree of this node **/
    size_type degree() const {
      return graph_->adjacency_list[index()].size();
    }

    /** Return the start incident_iterator for this node**/
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, uid_, graph_->adjacency_list[index()].begin());
      // return IncidentIterator(graph_, uid_, 0);
    }

    /** Return the end incident_iterator for this node**/
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, uid_, graph_->adjacency_list[index()].end());
      // return IncidentIterator(graph_, uid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      // std::cout << "Node::operator==(n)" << std::endl;
      // if (index() == n.index()) // check that index is the same
      if (uid_ == n.uid_) // check that index is the same
          if (graph_ == n.graph_) // check that graph is the same
              return true;
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
      // HW0: YOUR CODE HERE
      return (uid_ < n.uid_);
    }

    /** Test whether a node is valid */
    bool valid() const {
      // if (!((uid_) >= 0)) std::cout << "uid_ >= " << uid_ << " >= 0" << std::endl;
      // if (!(uid_ < graph_->nodes_.size())) std::cout << "uid_ < graph_->nodes_.size()" << std::endl;
      // if (!(graph_->idx_to_uid[index()] == uid_)) std::cout << "graph_->idx_to_uid[index()] == uid_" << std::endl;
      // if (!(graph_->nodes_[uid_].uid == uid_)) std::cout << "graph_->nodes_[uid_].uid == uid_" << std::endl;
      // if (!(graph_->nodes_[uid_].idx < graph_->idx_to_uid.size())) std::cout << "graph_->nodes_[uid_].idx < graph_->idx_to_uid.size()" << std::endl;
      return (uid_) >= 0 && (uid_ < graph_->nodes_.size())
             && (graph_->idx_to_uid[index()] == uid_) 
             && (graph_->nodes_[uid_].uid == uid_) 
             && (graph_->nodes_[uid_].idx < graph_->idx_to_uid.size())
             && true;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph_;          // pointer to graph
    size_type uid_;         // Unique identification number of element
    // Private constructor
    Node(const Graph* graph, size_type uid) : graph_(const_cast<graph_type*>(graph)), uid_(uid){}
    // Copy constructor
    // Node(const Node&);
    // Node& operator=(const Matrix&);

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // return idx_to_uid.size();
    return nodes_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    // return size();
    // return nodes_.size();
    return idx_to_uid.size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] val The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    // HW0: YOUR CODE HERE
    node_value_type myval = val;
    internal_node i_node(
      position,
      myval,
      size(), // uid
      num_nodes() // idx
      );
    // add uid of node to idx_to_uid container
    idx_to_uid.push_back(i_node.uid);
    // add i_node to this->nodes_
    nodes_.push_back(i_node); 
    // adjacency_list.push_back(std::vector<size_type>());
    adjacency_list.push_back(std::vector<internal_edge>());

    return Node(this,i_node.uid);
  }

  /** Remove Node from graph
   * @param[in] n     Node to be removed
   * 
   * @post if n is a valid node:
   *          graph no longer has node n or any edge to or from n
   *          new size() == old size() - 1
   *          new n_edges == old n_edges - n.degree()
   *          result = 1
   *       if n is not a valid node:
   *          no changes to graph
   *          result = 0
   * 
   * Complexity: O(n.degree()*D + index()) where D is the average 
   *          degree of all nodes connected to n
   */
  size_type remove_node(const Node& n){
    if (!has_node(n)) return 0;
    // remove all edges involving n
    while (n.degree() > 0){
      remove_edge(*(n.edge_begin()));
    }
    // remove from adjacency_list
    adjacency_list.erase(adjacency_list.begin()+n.index());
    idx_to_uid.erase(idx_to_uid.begin() + n.index());
    for (size_type id = n.uid_ + 1; id < size(); ++id){
      nodes_[id].idx -= 1;
    } 
    // remove n from this->idx_to_uid
    return 1; 
  }

  /** Remove Node from graph by iterator 
   * @param[in] n_it    iterator pointing to node to be deleted
   * 
   * @post result = n_it
   * @post all conditions of remove_node(*n_it) apply
   * 
   * Complexity: same as for remove_node(const Node& n)
  */
  node_iterator remove_node(node_iterator n_it){
    // if (remove_node(*n_it)) ++n_it;
    remove_node(*n_it);
    return n_it;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (this == n.graph_){
      return n.valid();
    } else return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i < num_nodes())
        // return Node(this,nodes_[idx_to_uid[i]].uid);
        return Node(this,idx_to_uid[i]);
    return Node();
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
    Edge() : node1_(Node()), node2_(Node()) {}

    /** Return a reference to the value of this Edge */
    edge_value_type& value() {
      for (auto it = node1().edge_begin(); it != node1().edge_end(); ++it){
        if (*this == *it){
          return *(it.it).val;
        }
      }
      return edge_value_type();
    }

    /** Return the value of this Edge */
    const edge_value_type& value() const {
      for (auto it = node1().edge_begin(); it != node1().edge_end(); ++it){
        if (*this == *it){
          return *(it.it).val;
        }
      }
      return edge_value_type();
    }

    /** Return the length of this Edge */
    double length() const {
      return norm(node2().position() - node1().position());
    }

    /** Return a node of this Edge */
    Node node1() const {
        return node1_;
    }

    /** Return the other node of this Edge */
    Node node2() const {
        return node2_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
        return (((node1() == e.node1()) && (node2() == e.node2())) 
                || ((node1() == e.node2()) && (node2() == e.node1())));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (node1() < e.node1()) return true;
      else if (node2() < e.node2()) return true;
      else if (node1().graph_ < e.node1().graph_) return true;
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Node node1_;
    Node node2_;
    // private constructor
    Edge(Node n1, Node n2) : node1_(n1), node2_(n2) {}
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
    // std::cout << "edge(i)" << std::endl;
    EdgeIterator it = edge_begin();
    for (size_type idx = 0; idx < i; ++idx){
      ++it;
    }
    return *it;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    if (!(has_node(a) && has_node(b))) return false;
    Edge e = Edge(a,b);
    for (auto it = a.edge_begin(); it != a.edge_end(); ++it){
      if (*it == e){
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
    // HW0: YOUR CODE HERE
    if (has_node(a) && has_node(b)){
      // return edge if it already exists
      if (has_edge(a,b)) return Edge(a,b); 
      // Otherwise add edge (both orientations) to adjacency_list
      adjacency_list[a.index()].push_back(internal_edge(a,b));
      adjacency_list[b.index()].push_back(internal_edge(b,a));
      // Increment edge counter
      n_edges += 1;
    }
    return Edge(a,b);
  }

  /** Remove edge from graph
   * @param[in] a   source node
   * @param[in] b   destination node
   * 
   * @pre src and dst are Nodes
   * @post if graph had an edge from a to b
   *            graph no longer has edge from a to b or from b to a
   *            new n_edges = old n_edges - 1
   *            result = 1
   *       if graph did not have an edge from a to b
   *            no changes to graph
   *            result = 0
   * 
   * Complexity: O(a.degree() + b.degree())
  */
  size_type remove_edge(const Node& a, const Node& b){
    if (!has_edge(a,b)){
      return 0;
    }
    // remove Edge(a,b) from adjacency_list[a.index()]
    adjacency_list[a.index()].erase(
      std::remove_if(
          adjacency_list[a.index()].begin(),
          adjacency_list[a.index()].end(),
          [this,a,b](internal_edge ie)->bool{return (this->edge(ie) == Edge(a,b));}),
      adjacency_list[a.index()].end());
    // remove Edge(a,b) from adjacency_list[b.index()]
    adjacency_list[b.index()].erase(
      std::remove_if(
          adjacency_list[b.index()].begin(),
          adjacency_list[b.index()].end(),
          [this,a,b](internal_edge ie)->bool{return (this->edge(ie) == Edge(b,a));}),
      adjacency_list[b.index()].end());
    // Decrement edge counter
    n_edges -= 1;
    return 1;
  }

  /** Remove edge e from graph 
   * @param[in] e   the edge to be removed
   * 
   * @post ( see conditions of remove_edge(e.node1(), e.node2()) )
   * Complexity: ( see conditions of remove_edge(e.node1(), e.node2()) )
  */
  size_type remove_edge(const Edge& e){
    return remove_edge(e.node1(),e.node2());
  }

  /** Erase edge correspondng to edge iterator e_it
   * @param[in,out] e_it    points to the edge to be removed
   * 
   * @pre @a e_it           is a valid edge_iterator
   * @post let L = adjacency_list[*(e_it.node_it).index()].size()
   *       if old @a e_it == (*e_it.node_it).edge_end() => new L == old L
   *       if old @a e_it != (*e_it.node_it).edge_end() => new L == old L - 1.
   *              new this->n_edges = old this->n_edges - 1
   * 
   * Complexity: O(1)
  */
  edge_iterator remove_edge(edge_iterator e_it){
    if (e_it.incident_iterator == *(e_it.node_it).edge_end())
      return e_it;
    adjacency_list[*(e_it.node_it)].erase(e_it.incident_it);
    n_edges -= 1;
    // e_it now points to the next element 
    return e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    idx_to_uid.clear();
    for (auto it = adjacency_list.begin(); it != adjacency_list.end(); ++it){
      (*it).clear(); // clear each individual vector within adjacency list
    }
    adjacency_list.clear();
    n_edges = 0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
      graph_ = nullptr;
      idx_ = 0;
    }

    /** Get Node corresponding to this iterator **/
    Node operator*() const {
      return graph_->node(idx_); // return a node matching the iterator idx
    }
    /** Increment iterator **/
    NodeIterator& operator++() {
      ++idx_;
      return *this;
    }
    /** Equality comparison between this iterator and another NodeIterator 
     * @param[in] iter    the other node_iterator 
     * **/
    bool operator==(const NodeIterator& iter) const {
      return (idx_ == iter.idx_);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    graph_type* graph_; // pointer to graph
    size_type idx_; // index into graph_->nodes_ (the graph's node vector)
    NodeIterator(const graph_type* graph, size_type idx) : graph_(const_cast<graph_type*>(graph)), idx_(idx) {}
  };

  /** Return start iterator pointing to the first node in the graph **/
  node_iterator node_begin() const {
    return NodeIterator(this,0); // 0 is the initial index
  }
  /** Return end iterator pointing just past the last node in the graph **/
  node_iterator node_end() const {
    // return NodeIterator(this,this->idx_to_uid.size()); // size of node vector
    return NodeIterator(this,this->num_nodes()); // size of node vector
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() : graph_(nullptr) {
    }

    /** Return Edge corresponding to this iterator **/
    edge_type operator*() const {
      return graph_->edge(*it);
    }

    /** Increment the iterator **/
    incident_iterator& operator++() {
      ++it;
      return *this;
    }

    /** Check equality with another incident_iterator @a iter 
     * @param[in] iter      the other incident_iterator
     * **/
    bool operator==(const incident_iterator& iter) const {
      return it == iter.it;
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type src_id;   // gives the uid_ of the node
    typename std::vector<internal_edge>::iterator it; // an iterator into graph->adjacency_list[src.index()]
    IncidentIterator(const graph_type* g, size_type s, typename std::vector<internal_edge>::iterator i) : 
      graph_(const_cast<graph_type*>(g)), src_id(s), it(i) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() : 
      graph_(nullptr), 
      node_it(),
      incident_it() 
      {}

    /** get edge element corresponding to this iterator 
     * **/
    edge_type operator*() const {
      return *incident_it;
    }
    /** Increment edge iterator 
     * 
     * @pre   @a this != graph_->edge_end() 
     * @pre   @a this->node_it != graph_->node_end() // implied by the previous condition
     * @pre   @a this->incident_it != (*(this->node_it)).edge_end()
     * **/
    edge_iterator& operator++(){
      // Iterate over the adjacency_list 
      ++incident_it;
      // Check that incident_it does not point past end of incident_edges
      if (incident_it == (*node_it).edge_end()){
        // Move to next node
        ++node_it;
        while (node_it != graph_->node_end()) {
          // Check if node has no edges
          if ((*node_it).degree() == 0){
            ++node_it;
          } else break;
        }
        // Check that node_it is not past end of nodes
        if (node_it != graph_->node_end()) 
          incident_it = (*node_it).edge_begin(); 
        else {
          incident_it = IncidentIterator();
          return *this;
        } 
      }
      // Check if we are on the upper triangle (j > i) of the adjacency matrix
      if ((*node_it) > (*incident_it).node2()) {
        // Skip this edge
        ++(*this);
      }
      return *this;
    }

    /** Check equality betwen edge_iterators **/
    bool operator==(const edge_iterator& iter) const {
      return (node_it == iter.node_it) && (incident_it == iter.incident_it);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type *graph_;
    NodeIterator node_it;
    IncidentIterator incident_it;
    EdgeIterator(const graph_type* g, NodeIterator n_it) : 
      graph_(const_cast<graph_type*>(g)), 
      node_it(n_it),
      incident_it((*n_it).edge_begin())
      {}
    EdgeIterator(const graph_type* g, NodeIterator n_it, IncidentIterator i_it) : 
      graph_(const_cast<graph_type*>(g)), 
      node_it(n_it),
      incident_it(i_it)
      {}

  };

  /** Return begin edge_iterator pointing to first edge 
   * this requires iterating through nodes until finding one with 
   * at least one edge. If no nodes have edges, edge_begin() == edge_end()
   * **/
  edge_iterator edge_begin() const{
    if (num_edges() == 0) return edge_end();
    for (auto it = node_begin(); it != node_end(); ++it){
      if ((*it).degree() != 0){
        return EdgeIterator(this,it,(*it).edge_begin());
      }
    }
    return edge_end();
  }
  /** Return end edge_iterator pointing past last edge */
  edge_iterator edge_end() const{
    return EdgeIterator(
            this,
            node_end(),
            IncidentIterator());
  }

  //
  // Helpers
  //

  /** verify that the number of edges matches the quantity of edges 
   * iterated through between edge_begin() and edge_end()
   * */
  bool verify_num_edges(){
    // std::cout << "VERIFYING NUM EDGES...\n";
    size_type n = 0;
    for (auto it = edge_begin(); it != edge_end(); ++it) 
        ++n;
    if (n == num_edges()){
      return true;
    } else {
      std::cout << "n = " << n << ", num_edges() = " << num_edges() << "\n";
    }
    return false;
  }

  /** verify that each edge between edge_begin() and edge_end() is part of the
   * graph and that each node exists too 
   * */
  bool verify_edges(){
    // std::cout << "VERIFYING EDGES...\n";
    for (auto it = edge_begin(); it != edge_end(); ++it){
      Edge e = *it;
      // std::cout << "(" << e.node1().index() << "," << e.node2().index() << "),";
      if (!has_node(e.node1())){
        std::cout << "!has_node(" << e.node1().index() << ")" << std::endl;
        return false;
      }
      if (!has_node(e.node2())){
        std::cout << "!has_node(" << e.node2().index() << ")" << std::endl;
        return false;
      }
      if (!has_edge(e.node1(),e.node2())){
        std::cout << "!has_edge(" << e.node1().index() << "," << e.node2().index() << ")" << std::endl;
        return false;
      }
    }
    return true;
  }

  // void check_move_assignable(){
  //   std::cout << "internal_node: " << std::is_move_assignable<internal_node>::value << "\n";
  //   std::cout << "internal_edge: " << std::is_move_assignable<internal_edge>::value << "\n";
  // }

 private:

  struct internal_node {
    Point position;
    node_value_type val;
    size_type uid;
    size_type idx;
    internal_node() : position(), val(), uid(), idx() {}
    internal_node(const Point& pos, size_type id) : position(pos), val(), uid(id), idx() {}
    internal_node(const Point& pos, node_value_type v, size_type id) : position(pos), val(v), uid(id), idx() {}
    internal_node(const Point& pos, node_value_type v, size_type id, size_type i) : position(pos), val(v), uid(id), idx(i) {}
  };
  struct internal_edge {
    size_type node1_uid;
    size_type node2_uid;
    edge_value_type val;
    internal_edge() : node1_uid(), node2_uid(), val() {}
    internal_edge(size_type u1, size_type u2) : node1_uid(u1), node2_uid(u2), val() {}
    internal_edge(Node n1, Node n2) : node1_uid(n1.uid_), node2_uid(n2.uid_), val() {}
    internal_edge(size_type u1, size_type u2, edge_value_type v) : 
      node1_uid(u1), node2_uid(u2), val(v) {}
  };
  // nodes of graph - nodes_[n.uid] = n
  std::vector<internal_node> nodes_;
  // vector of node uids, such that idx_to_uid[n.idx] = n.uid (way faster than a map!)
  std::vector<size_type> idx_to_uid;
  // Adjacency list - maps node.index() to a vector of node_ids
  std::vector<std::vector<internal_edge>> adjacency_list;
  // counter for number of edges
  int n_edges;
  /** Helper method to get an edge directly from node uids 
   * @param[in] src_id      uid of first node
   * @param[in] dst_id      uid of second node
   * **/
  Node node(internal_node i_node){
    return Node(this,i_node.uid);
  }
  Edge edge(size_type src_id, size_type dst_id) {
    return Edge(Node(this,src_id),Node(this,dst_id));
  }
  Edge edge(internal_edge i_edge) const {
    return Edge(Node(this,i_edge.node1_uid),Node(this,i_edge.node2_uid));
  }
};

#endif // CME212_GRAPH_HPP
