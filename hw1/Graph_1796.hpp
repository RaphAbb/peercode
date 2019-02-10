#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_set>
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

  using node_value_type = V;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
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
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph->points[idx];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return idx;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
     
    node_value_type& value(){
     return graph->node_values[idx];
    }

    const node_value_type& value() const{
     return graph->node_values[idx];
    }

    size_type degree() const{
     return graph->edges[idx].size();
    }

    /*
    * @return edge iterator
    * @post iterator is positioned at the start of the node's
    * edge list
    */
    incident_iterator edge_begin() const {
     return IncidentIterator(graph, idx, 0);
    }

    /*
    * @return edge iterator
    * @post iterator is positioned at one past the end of the node's
    * edge list
    */
    incident_iterator edge_end() const {
     return IncidentIterator(graph, idx, graph->edges[idx].size());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return graph == n.graph && idx == n.idx;
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
      return graph == n.graph && idx < n.idx;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    graph_type* graph;
    size_type idx;
    Node(const graph_type* graph, size_type idx) : 
         graph(const_cast<graph_type*>(graph)), idx(idx){};
};

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return points.size();
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
  Node add_node(const Point& position, const node_value_type& data = node_value_type()) {
    // HW0: YOUR CODE HERE
    points.push_back(position);
    node_values.push_back(data);
    std::vector<size_type> list;
    edges.push_back(list);
    return Node(this, points.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return n.graph == this;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    if(i < size()){
     return Node(this, i);
    } else { 
     return Node(); //return invalid node if out of bounds
    }    
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
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return graph->node(node_1_idx);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph->node(node_2_idx);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return graph == e.graph && 
             ((node_1_idx == e.node_1_idx && node_2_idx == e.node_2_idx) ||
             (node_1_idx == e.node_2_idx && node_2_idx == e.node_1_idx));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const { 
      return (node_1_idx^node_2_idx) < (e.node_1_idx^e.node_2_idx);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects 
    graph_type* graph;
    size_type node_1_idx;
    size_type node_2_idx;
    Edge(const graph_type* graph, size_type node_1_idx, size_type node_2_idx) : 
     graph(const_cast<graph_type*>(graph)), node_1_idx(node_1_idx),
     node_2_idx(node_2_idx){};
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_cnt;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
  
    if(i < num_edges()){
      size_type n = 0;
      for(auto iter = edge_begin(); iter != edge_end(); ++iter){
       if(n == i){
        return *iter;
       }else {
        n++;
       }
      } 
    } else {
     return Edge(); //return invalid edge if out of bounds
    } 
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
  	if(a.graph != this || b.graph != this ||
       a.idx >= size() || b.idx >= size()){
      return false;
    } else {
      for(auto adjacent: edges[a.idx]){
       if(adjacent == b.idx){
         return true;
       }
      }
      return false;
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
    if(!has_edge(a,b)){
      edge_cnt++;
      edges[a.idx].push_back(b.idx);
      if(a.idx != b.idx){
       edges[b.idx].push_back(a.idx);
      }
      return Edge(this, a.idx, b.idx);
    } else {
	 return Edge();
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HEREi
    points.clear();
	edges.clear();
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
    NodeIterator() {
    }

    NodeIterator(const Graph* graph_, size_type pos_){
     graph = const_cast<Graph*>(graph_);
     pos = pos_;
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const
    
    /*
    * @return Node object
    * @post returned node is located at the iterator's
    * position in the node list
    */
    Node operator*() const{
     return Node(graph, pos);
    }

    /*
    * @return NodeIterator
    * @post returned NodeIterator is positioned at the 
    * next node in the node list
    */
    NodeIterator& operator++(){
     pos++;
     return *this;
    }

    /*
    * @return boolean
    * @post returns true if both iterators belong to the
    * same graph and have the same position
    */
    bool operator==(const NodeIterator& iter) const{
     return iter.graph == graph && iter.pos == pos;
    }


   private:
    friend class Graph;
    Graph* graph;
    size_type pos;
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /*
  * @return node iterator
  * @post iterator is positioned at the start of the graph's
  * node list
  */
  node_iterator node_begin() const{
   return NodeIterator(this, 0);
  }


  /*
  * @return node iterator
  * @post iterator is positioned at one past the end of the graph's
  * node list
  */
  node_iterator node_end() const{
   return NodeIterator(this, size());
  }
  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
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


    IncidentIterator(const Graph* graph_, size_type idx_, size_type pos_) {
     graph = const_cast<Graph*>(graph_);
     idx = idx_;
     pos = pos_;
    }
    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    
    /*
    * @return Edge object
    * @post returned edge is located at the iterator's
    * position in the node's edge list
    */
    Edge operator*() const{
     return Edge(graph, idx, graph->edges[idx][pos]);
    }

    /*
    * @return IncidentIterator
    * @post returned IncidentIterator is positioned at the 
    * next edge in the node's edge list
    */
    IncidentIterator& operator++(){
      pos++;
      return *this;
    }

    /*
    * @return boolean
    * @post returns true if both iterators belong to the same
    * graph and node, and their positions are the same
    */
    bool operator==(const IncidentIterator& iter) const{
     return iter.graph == graph && iter.idx == idx && iter.pos == pos;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph;
    size_type idx;
    size_type pos;
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
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

    EdgeIterator(const Graph* graph_, size_type node_, size_type adj_) {
     graph = const_cast<Graph*>(graph_);
     node = node_;
     adj = adj_;
     setpos();
    }
    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const
    
    /*
    * @return Edge object
    * @post returned edge is located at the iterator's
    * position in the graph's adjacency list
    */
    Edge operator*() const {
     return Edge(graph, node, graph->edges[node][adj]);
    }


    /*
    * @return EdgeIterator
    * @post returned EdgeIterator is positioned at the 
    * next edge in the graph's adjacency list
    */
    EdgeIterator& operator++(){
     adj++;
     setpos();
     return *this;
    }

    /*
    * @return boolean
    * @post returns true if both iterators belong to the same
    * graph and their positions in the adjacency list are the same
    */
    bool operator==(const EdgeIterator& iter) const{
     return iter.graph == graph && iter.node == node && iter.adj == adj;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    /*
    * Iterate through the adjacency list until the
    * next valid element is found.
    */
    void setpos(){
	while(node < graph->edges.size()){
	 while(adj < graph->edges[node].size()){
	  if(graph->edges[node][adj] < node){
		return;
	  }
	  adj++;
	 }
     node++;
     adj = 0;
     }   
	}

     Graph* graph;
     size_type node = 0;
     size_type adj = 0;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const


  /*
  * @return edge iterator
  * @post iterator is positioned at the start of the graph's
  * edge list
  */
  edge_iterator edge_begin() const {
   return EdgeIterator(this, 0 , 0);
  }
 
  /*
  * @return edge iterator
  * @post iterator is positioned at one past the end of the graph's
  * edges
  */
  edge_iterator edge_end() const {
   return EdgeIterator(this, edges.size()-1 , edges[edges.size()-1].size());

  }
 private:
  

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.i
  size_type edge_cnt = 0;  
  std::vector<node_value_type> node_values;
  std::vector<std::vector<size_type>> edges;
  std::vector<Point> points;

};

#endif // CME212_GRAPH_HPP
