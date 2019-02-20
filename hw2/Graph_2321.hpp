#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

static const unsigned INVALID_VALUE= unsigned(-1);


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V= float, typename E= float>
class Graph {

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

  using node_value_type = V;
  using edge_value_type = E;
  /** Node_Information keeps all of the necessary information of a node. 
  */
 private:
  struct Node_Information {
    // The unique identiy of a node
    size_type num;
    // stores the position of anode
    Point point;
    // stores how many neighbors that have a smaller index than the current node. 
    size_type smaller_num;
    // a value of a node.
    node_value_type node_value;
    /** A constructor.*/
    Node_Information (size_type num_, const Point& point_ , node_value_type node_value_ , size_type smaller_num_ =0){
      num= num_ ;
      point= point_;
      smaller_num= smaller_num_;
      node_value= node_value_;
    }

  };
  /**
  * A structure that keeps necessary information of edges. 
  * neighbor_node is the index of the neighbor node (node2) of this edge.
  * edge_value keeps the value of some attributes of nodes. For application in mass_spring, this stores the distance of this edge.
  */
  class Edge_Information{
  public: 
    // A constructor 
    Edge_Information(const size_type& neighbor_node_= unsigned(-1), const edge_value_type& edge_value_ = unsigned(-1)){
      neighbor_node= neighbor_node_;
      edge_value= edge_value_;
    }
    /** Defien operators that is going to be of future use in comparisons
    * The comparison is based on the index of the neighbor node (node 2)
    */
    bool operator<(const Edge_Information& another_edge) const {
      return this-> neighbor_node< another_edge.neighbor_node;
    }

    bool operator== (const Edge_Information& another_edge) const {
      return this-> neighbor_node == another_edge.neighbor_node;
    }
  private:

    friend class Graph;
    /** The index of node2 of this edge.*/
    size_type neighbor_node;
    /** An attribute that is related to the edge.*/
    edge_value_type edge_value;    

  };

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //
 public:
  /** Construct an empty graph. */
  Graph() {
  	number_edges= 0;
  	nodes_index[INVALID_VALUE]= unsigned(-1);
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
    // Construct default and invalid node
    Node() {
    	Graph_N = NULL;
    	num= INVALID_VALUE; 
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return Graph_N->Nodes_Information[Graph_N->nodes_index.at(num)].point;
    }
    /** Return the position of a node */
    Point& position() {
      return Graph_N->Nodes_Information[Graph_N->nodes_index.at(num)].point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return Graph_N->nodes_index.at(num);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Return the attribute of a node.  */
    node_value_type& value(){
      return Graph_N->Nodes_Information[Graph_N->nodes_index.at(num)].node_value;
    }

    /** A rewrite of the function value() above that is const.*/
    const node_value_type& value() const{
      return Graph_N->Nodes_Information[Graph_N->nodes_index.at(num)].node_value;
    }
    /** Returns the degree of a node */
    size_type degree() const{
      return Graph_N->Edges_Information[Graph_N->nodes_index.at(num)].size();
    }
    /** return an incidend_iterator that points to the first edge of the node.*/
    incident_iterator edge_begin() const{
      return incident_iterator(Graph_N, index(), 0);
    }
    /** Return an incident iterator that points to the last edge of the node.*/
    incident_iterator edge_end() const{
      return incident_iterator(Graph_N, index(), this->degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
             // Quiet compiler warning
      return this->Graph_N== n.Graph_N && this->num == n.num;
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
      if (this->num<n.num){
      	return true;
      }         
      if (this->num>n.num){
      	return false;
      }
      return this->Graph_N< n.Graph_N;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    /** The graph that this node relates to */
    graph_type* Graph_N;
    /** The unique identity of this node, num*/
    size_type num;
    /** A constructor that takes the graph that this node is related to, and the identity of the node as argument.*/
    Node (const graph_type* Graph_N_, size_type num_){
    	Graph_N= const_cast<graph_type*>((Graph_N_)), num = num_;
    }
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    // std::cout<< Nodes_Information.size() <<std::endl;
    return Nodes_Information.size();
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
  Node add_node(const Point& position, const node_value_type& val_new= node_value_type()) {
    // HW0: YOUR CODE HERE
    size_type value_new= rand();
    while (nodes_index.count(value_new)){
    	value_new= rand();
    }
    this->nodes_index[value_new]= this->size();
    this->Nodes_Information.push_back(Node_Information(value_new, position, val_new));
    this->Edges_Information.push_back(std::vector<Edge_Information>(0));
    return Node(this, value_new); 
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
          // Quiet compiler warning
    return n.Graph_N== this && this->nodes_index.count(n.num) && n.num!= INVALID_VALUE;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node(this, Nodes_Information[i].num);        // Invalid node
  }

  /**
  * This function is for removing a node. If this node removed is not the last node of the graph 
  * in terms of index, then this function will move the last node to the location
  * of this removed code.
  * 
  * @param n is the node to be removed
  * @return The number of nodes in the graph after the function. 
  *
  * @pre n has to be already in the graph.
  * @post new num_edges() = old num_edges() - n.degree()
  * @post new num_nodes() = old num_nodes() - 
  * @post all information related to this removed node will get destroyed or move to the corresponding
  * end node.
  * Complexity: O(num_nodes())
  */

  size_type remove_node(const Node& n){
    size_type nidx= n.index();
    for (auto eit= n.edge_begin(); eit!= n.edge_end(); ++eit){
      size_type neighbor_idx= (*eit).node2().index();
      removeeleVec(Edges_Information[neighbor_idx], Edge_Information(nidx));
      if (nidx< neighbor_idx){
        Nodes_Information[neighbor_idx].smaller_num --; 
      }
    }
    number_edges -= n.degree();

    size_type end_idx= num_nodes() - 1; 

    if (nidx<end_idx){
      nodes_index[Nodes_Information[end_idx].num] = nidx;
      Nodes_Information[nidx] = Nodes_Information[end_idx];
      Edges_Information[nidx] = Edges_Information[end_idx];

      for (Edge_Information edgeinfo : Edges_Information[end_idx]){
        size_type neighbor_index = edgeinfo.neighbor_node;
        Edges_Information[neighbor_index].pop_back();
        edgeinfo.neighbor_node = nidx;
        addeleVec(Edges_Information[neighbor_index], edgeinfo);
        if (nidx< neighbor_index){
          Nodes_Information[neighbor_index].smaller_num ++;
          Nodes_Information[nidx].smaller_num -- ;
        }
      }
    }
    nodes_index.erase(n.num);
    // std::cout<<Nodes_Information.size()<<std::endl;
    Nodes_Information.pop_back();
    Edges_Information.pop_back();

    return num_nodes();
  }
  /** 
  * This function implements a node_iterator implementation of the same function remove_node(const Node& n)
  */
  node_iterator remove_node (node_iterator ni){
    remove_node(*ni);
    return ni;
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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    	Graph_e= NULL;
    	value_1= INVALID_VALUE;
    	value_2= INVALID_VALUE;
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(Graph_e, value_1);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(Graph_e, value_2);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return this->Graph_e == e.Graph_e && this->value_2== e.value_2 && this->value_1== e.value_1;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (this->value_1< e.value_1){
        return true;
      }
      if (this->value_1>e.value_1){
        return false;
      }
      if (this->value_2 < e.value_2){
        return true;
      }
      if (this->value_2 > e.value_2){
        return false;
      }
      return this-> Graph_e < e. Graph_e;
    }
    /**Returns the attribute of the edge. Makes use of the helper function true_value_fun below defined.*/
    edge_value_type& value(){
      return const_cast<edge_value_type&>(true_value_fun());
    }
    /**Returns the attribute of the edge. Makes use of the helper function true_value_fun below defined. */
    const edge_value_type& value() const {
      return true_value_fun();
    }
    /** Returns the length of an edge. */
    double length() const {
      return norm_2(node1().position()- node2().position());
    }


   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    /** A graph this edge relates to */
    graph_type* Graph_e;
    /** The nums of the two nodes that belongs to this edge*/
    size_type value_1, value_2;
    /**A legitmate edge that takes a graph, the nums of the two nodes of the edge. */
    Edge(const graph_type* Graph_e_, size_type value_1_, size_type value_2_){
    	Graph_e= const_cast<graph_type*>((Graph_e_)), value_1= value_1_, value_2= value_2_;
    }
    /** 
    * This is an internal helper function to let my code for value() be simpler.
    * This function is for getting the edge_value of an edge. 
    * @pre an edge to find its edge_value
    * @return return the edge's edge_value
    * @post find its corresponding Edge_Information, and get the edge_value.
    */
    const edge_value_type& true_value_fun() const{
      Node n1= node1();
      Node n2= node2();
      size_type node_idx1= n1.index();
      size_type node_idx2= n2.index();
      auto it= lower_bound(Graph_e->Edges_Information[node_idx1].begin(), Graph_e->Edges_Information[node_idx1].end() , Edge_Information(node_idx2));
      return it->edge_value;
    }




    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return number_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   * I know that because I have already implemented the smaller index of nodes, I
   *'m not going to traverse the same node twice but only consider nodes with a larger index.
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    size_type node_1= 0;             // Quiet compiler warning
    while (i>= Nodes_Information[node_1].smaller_num){    	
    	i-= Nodes_Information[node_1].smaller_num;
      node_1++;
    }
    return Edge(this, Nodes_Information[node_1].num, Nodes_Information[Edges_Information[node_1][i].neighbor_node].num);        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    size_type a_index= a.index() , b_index= b.index();
    auto it= lower_bound(Edges_Information[a_index].begin(), Edges_Information[a_index].end(), Edge_Information(b_index));
    if (it== Edges_Information[a_index].end()){
    	return false;
    }
    else{
    	return (*it).neighbor_node == b_index;
    }
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
	 assert (!(a== b));   // Make sure that a and b are not the same node.

	 size_type index_a= a.index(), index_b= b.index();
	 if (!this->has_edge(a, b)){
		  number_edges++;
      Nodes_Information[std::max(index_a, index_b)].smaller_num++;
		  auto it_a = lower_bound(Edges_Information[index_a].begin(), Edges_Information[index_a].end(), Edge_Information(index_b));
      Edges_Information[index_a].insert(it_a, Edge_Information(index_b));
      auto it_b= lower_bound(Edges_Information[index_b].begin(), Edges_Information[index_b].end(), Edge_Information(index_a));
      Edges_Information[index_b].insert(it_b, Edge_Information(index_a));
	 } 
    return Edge(this, Nodes_Information[index_a].num, Nodes_Information[index_b].num);        // Invalid Edge
  }
  /**
  * Remove an edge in the graph. If the edge is already there, then do nothing about it
  * @param n1, n2: the two nodes of the edge to be removed. 
  * @return the number of edges of the new graph that remove this edge.
  * @post new number_edges= number_edges -1 
  * @post new Nodes_Information and Edges_Information update all of the relevant information of nodes and edges. 
  */

  size_type remove_edge(const Node& n1, const Node& n2){
    if (!has_edge(n1, n2)){
      return number_edges;
    }
    size_type n1idx= n1.index();
    size_type n2idx= n2.index();
    removeeleVec(Edges_Information[n1idx], Edge_Information(n2idx)); 
    removeeleVec(Edges_Information[n2idx], Edge_Information(n1idx));
    number_edges = number_edges -1;
    Nodes_Information[std::max(n1idx, n2idx)].smaller_num --; 
    return number_edges;
  }
  /** Implements the same function as remove_edcge(const Node& n1, const Node& n2) in essence.
  * @param edge: The edge to be removed.
  * the other specifications are the same as remove_edge(const Node& n1, const Node& n2) 
  */
  size_type remove_edge(const Edge& edge){
    return remove_edge(edge.node_1(), edge.node_2());
  }
  /** Implements the same function as remove_edge(const Edge& edge) in essence.
  * @param eit: the edge_iterator pointing to the edge that the function wants to remove.
  * The other specifications are the same as remove_edge(const Edge& edge)
  */
  edge_iterator remove_edge(edge_iterator eit){
    Edge e= *eit;
    remove_edge(e); 
    return eit.update_edge();
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_index.clear();
    Edges_Information.clear();
    Nodes_Information.clear();
    number_edges= 0;
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
      Graph_NI= NULL;
      node_idx= unsigned(-1);
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const
    /** De-reference a node iterator return the node that it is pointing to */
    Node operator*() const {
      return Node (Graph_NI, Graph_NI->Nodes_Information[node_idx].num);
    }
    /** iterate to the next node.*/
    NodeIterator& operator++() {
      node_idx++;
      return *this;
    }
    /** Return equality between node_iterator only if they point to the same node of the same graph.*/
    bool operator == (const NodeIterator& another_nodeiterator) const{
      bool rule= this->Graph_NI== another_nodeiterator.Graph_NI && this->node_idx== another_nodeiterator.node_idx; 
      return rule;
    }

   private:
    friend class Graph;

    size_type node_idx;

    const graph_type * Graph_NI;

    NodeIterator(const graph_type* Graph_NI_, size_type node_idx_){
      Graph_NI= Graph_NI_; node_idx= node_idx_;
    }
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  /** This is a node iterator that points to the first node. */
  node_iterator node_begin() const{
    return node_iterator(this, 0);
  }
  /** This is a node iterator that points to the last node.*/
  node_iterator node_end() const{
    return node_iterator(this, this->size());
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
      Graph_II= NULL;
      node_idx= unsigned(-1);
      neighbor= unsigned(-1);
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const
    /** Derefrence an incident iterator*/
    Edge operator*() const{
      size_type neighbor_idx= Graph_II->Edges_Information[node_idx][neighbor].neighbor_node;
      return Edge(Graph_II, Graph_II->Nodes_Information[node_idx].num, Graph_II->Nodes_Information[neighbor_idx].num);
    }
    /** Point to the next neighbor of the node.*/
    IncidentIterator& operator++() {
      neighbor= neighbor+ 1;
      return *this;
    }
    /** Two incident iterators are the same only if they point to the same 
    neighbor node of same node of the same graph.*/
    bool operator==(const IncidentIterator& another_incidentiterator) const{
      bool rule= this->Graph_II== another_incidentiterator.Graph_II && this-> node_idx== another_incidentiterator.node_idx && this-> neighbor== another_incidentiterator.neighbor;
      return rule;
    }


   private:
    /** The graph of this incident_iterator, the node index, and the index of the neighbor that
    this iterator points to */
    friend class Graph;

    const graph_type * Graph_II;

    size_type node_idx;

    size_type neighbor;
    /**
    * A incident_iterator that takes a grpah, the index of the current node, and the index of 
    * the neighbor that this iterator is pointing to.
    */
    IncidentIterator (const graph_type* Graph_II_, size_type node_idx_, size_type neighbor_){
      Graph_II= Graph_II_;
      node_idx= node_idx_;
      neighbor= neighbor_;
    }


    // HW1 #3: YOUR CODE HERE
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
    EdgeIterator() {
      Graph_EI= NULL;
      node_idx= unsigned(-1);
      neighbor= unsigned(-1);
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /**
    * This is for for providing a reference to an edge that satisfied my constraints and definitions.
    */
    Edge operator*() const{
      return Edge(Graph_EI, Graph_EI->Nodes_Information[node_idx].num, Graph_EI->Nodes_Information[Graph_EI->Edges_Information[node_idx][neighbor].neighbor_node].num);
    }
    /**
    * Keep updating the edge_iterator until I have found a legitmate edge or I have found an edge that point to the end of the edge list.
    */
    EdgeIterator& operator++(){
      neighbor++;
      return this->update_edge();
      // while (node_idx< Graph_EI->size() && neighbor== Graph_EI->Nodes_Information[node_idx].smaller_num){
      //   neighbor= 0;
      //   node_idx+=1;
      // }
      // return *this;
    }
    /** 
    * Define if two edge_iterators are the same. They are the same only if they point to the same edge in the graph.
    */
    bool operator==(const EdgeIterator& another_EI) const{
      return this->Graph_EI== another_EI.Graph_EI && this-> node_idx== another_EI.node_idx && this->neighbor== another_EI.neighbor;
    }


   private:
    friend class Graph;
    /**
    * Parameters are: Graph, node1's index node_idx, and neighbor's index neighbor.
    */
    const graph_type * Graph_EI;
    size_type node_idx;
    size_type neighbor;

    /** 
    * Provide an edgeiterator constructor that takes a graph, the index of the first node, and the index
    * of the second node as input parameters.
    */
    EdgeIterator(const graph_type* Graph_EI_, size_type node_idx_, size_type neighbor_){
      Graph_EI= Graph_EI_;
      node_idx= node_idx_;
      neighbor= neighbor_;
    }
    /**
    * This function is for updating an edge because sometimes the edge returned may not be legitmate.
    * keep updating until one of the two conditions are satisfied: 
    * it point to a well define edge
    * or it has already reached the end of the edge list, and keep updating.
    */
    EdgeIterator& update_edge(){
      while (node_idx< Graph_EI->size() && neighbor== Graph_EI->Nodes_Information[node_idx].smaller_num){
        node_idx ++;
        neighbor= 0;
      }
      return *this;
    }

    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  // Return the first edge. Continue to 
  /** An edge_iterator that points to the first edge of the graph. If not legitmate, iterate until it is*/
  edge_iterator edge_begin() const{
    return EdgeIterator(this,0,0).update_edge();
  }
  /** Return an edge_iterator that points to the space after the last edge of the graph.
  This is pretty useful for loops and iterations. */
  edge_iterator edge_end() const{
    return EdgeIterator(this, this->size(), 0);
  }


 private:

  /**
  * Explanation for all of the variables: 
  * Each node has two attributes that I consider. First, there is index, which denotes a counting of the nodes. 
  * Second, there is num, which is unique so that I have an identity for each node. 
  * 
  * nodes_index is the correspondence relationship between num and index so that I can access the index of each node.
  *
  * Nodes_Information is a vector of Node_Information, where each Node_Information contains the necessary information of a node.
  *
  * Edge_Information in a vector of the Edge_Informations of a node. That is, the information of the edges incident to the node.  
  *
  * number_edges return the number of edges of the graph. 
  */
 	size_type number_edges;

 	std::unordered_map<size_type, size_type> nodes_index;

 	std::vector<Node_Information> Nodes_Information;

 	std::vector<std::vector<Edge_Information>> Edges_Information;
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals: 
  //   helper functions, data members, and so forth.
  /**
  * This function can be used for adding edges and etc.
  * @para vec: elements that already existed and are sorted so that I can use lower_bound
  * @para val: an element to be added.  
  * 
  * @pre: vec are sorted in increasing order
  * @post: a new vector that get inserted the element.
  */

  template<typename T>
  void removeeleVec(std::vector<T>& vec, T val){
    auto it = lower_bound(vec.begin(), vec.end(), val);
    vec.erase(it);
  }
  /**
  * This function can be used for adding edges and etc.
  * @para vec: elements that already existed and are sorted so that I can use lower_bound
  * @para val: an element that has already been in the vector and is going to get removed.
  * 
  * @pre: vec are sorted in increasing order and val is the same as one of the elements in vec.
  * @post: a new vector that has realdy removed the element val.
  */
  template<typename T>
  void addeleVec(std::vector<T>& vec, T val){
    auto it = lower_bound(vec.begin(), vec.end(), val);
    vec.insert(it, val);
  }

};

#endif // CME212_GRAPH_HPP
