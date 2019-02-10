#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <set>
#include <cassert>
#include <limits.h>

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
    
  /** Type of value in each node */
  using node_value_type = V;

  /** Type of this graph. */
  using graph_type = Graph<node_value_type>;

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
  static const size_type INVALID = std::numeric_limits<size_type>::max();

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
      idx = INVALID;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph->point_list[idx];
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
      
    /** return node value */
    node_value_type& value(){
        return graph->node_values[idx];
    }
      
    /** return node value const ref */
    const node_value_type& value() const{
        return graph->node_values[idx];
    }
    
    /** Number of neighbouring nodes
    */
    size_type degree() const{
        return graph->adjList[idx].size();
    }
      
    /** Iterator begin for edges incident on this node
    */
    incident_iterator edge_begin() const{
         return incident_iterator(graph, idx, graph->adjList[idx].begin());
    }
      
    /** Iterator begin for edges incident on this node
    */
    incident_iterator edge_end() const{
         return incident_iterator(graph, idx, graph->adjList[idx].end());
    }
      
     

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return  n.index() == index() && n.graph == graph;
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
      return index() < n.index() ;
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    size_type idx;  
    Graph * graph;
    Node(size_type i,Graph* g) :idx(i), graph(g)  {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return point_list.size();
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
  Node add_node(const Point& position, const node_value_type& nv = node_value_type () ) {
    // HW0: YOUR CODE HERE
    point_list.push_back(position);
    Node new_node(point_list.size()-1, this);
    node_list.push_back(new_node);
    adjList.push_back(Set());
    node_values.push_back(nv);
    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return this == n.graph;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    if( i > num_nodes() ){
        std::cout<<"Node index "<<i<<" out of range [0,"<<num_nodes()<<")"<<std::endl;
        return Node();
    }
    return node_list[i];     
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
      idx = INVALID;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return graph->node(n1);     
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph->node(n2);      
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return e.graph == graph && e.idx == idx;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return idx < e.idx;
    }
      
    /** Return this edge's index, a number in the range [0, num_edges()). */
    size_type index() const {
      return idx;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    size_type idx,n1,n2;
    const Graph* graph;
    
    Edge(size_type n1,size_type n2,size_type i, const Graph* g) : idx(i), n1(n1), n2(n2), graph(g)  {}
    
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_list.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    //(void) i;             // Quiet compiler warning
    if( i > num_edges() ){
        std::cout<<"Edge index "<<i<<" out of range [0,"<<num_edges()<<")"<<std::endl;
        return Edge();
    }
    return edge_list[i];       
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    //(void) a; (void) b;   // Quiet compiler warning
    Pair pr(a.index(),b.index());
    return edge_map.find(pr) != edge_map.end();
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
    //(void) a, (void) b;   // Quiet compiler warning
    size_type aidx = a.index();
    size_type bidx = b.index();
    if( aidx == INVALID || bidx == INVALID || aidx ==bidx){
        std::cout<<"Cant add edge between node indices : "<<aidx<<", "<<bidx<<std::endl;
        return Edge();
    }
    Pair pr1(aidx,bidx);
    Pair pr2(bidx,aidx);
    if(edge_map[pr1].index() != INVALID)
        return edge_map[pr1];
    Edge edge(aidx,bidx,edge_list.size(), this);
    edge_list.push_back(edge);
    edge_map[pr1]=edge;
    edge_map[pr2]=edge;
      
    adjList[aidx].insert(bidx);
    adjList[bidx].insert(aidx);
    return edge;       
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    point_list.clear();
    edge_list.clear();
    node_list.clear();
    edge_map.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : equality_comparable<NodeIterator>{
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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const
      
    /** dereference the iterator
    */
    Node operator*() const{
        return graph->node_list[idx];
    }
      
    /** increment the iterator
    */
    NodeIterator& operator++(){
        idx = (idx < graph->node_list.size())? idx+1 :graph->node_list.size() ;
        return *this;
    }
      
    /** check equality
    */
    bool operator==(const NodeIterator& other) const{
        return other.idx == idx;
    }
      

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    const graph_type* graph;
    size_type idx;
    
      
    NodeIterator(const graph_type* g, size_type init_idx) :graph(g) ,idx(init_idx) {}

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
    
  /** Returns begin of nodes list iterator
  */
  node_iterator node_begin() const {
      return node_iterator(this,0);
  }
   
  /** Returns end of nodes list iterator
  */
  node_iterator node_end() const {
      return node_iterator(this,node_list.size());
  }
  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : equality_comparable<IncidentIterator>{
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
      
    /** dereference the iterator
    */
    Edge operator*() const{
        Pair pr(node_idx,*adj_set_iter);
        return Edge(node_idx, *adj_set_iter, graph->edge_map.at(pr).index(),graph);
    }
      
    /** increment the iterator
    */
    IncidentIterator& operator++(){
        adj_set_iter++;
        return *this;
    }
      
    /** check equality
    */
    bool operator==(const IncidentIterator& other) const{
        return other.adj_set_iter == adj_set_iter;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    using SetIter = std::set<size_type>::iterator;
     
    const graph_type* graph;
    size_type node_idx;
    SetIter adj_set_iter;
    
      
    IncidentIterator(const graph_type* g, size_type node_idx, SetIter adj_set_iter) 
        :graph(g) ,node_idx(node_idx), adj_set_iter(adj_set_iter) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator  : equality_comparable<EdgeIterator>{
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
      
    Edge operator*() const{
        return graph->edge_list[idx];
    }
      
    /** increment the iterator
    */
    EdgeIterator& operator++(){
        idx = (idx < graph->edge_list.size())? idx+1 :graph->edge_list.size() ;
        return *this;
    }
      
    /** check equality
    */
    bool operator==(const EdgeIterator& other) const{
        return other.idx == idx;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
      
    const graph_type* graph;
    size_type idx;
    
    EdgeIterator(const graph_type* g, size_type init_idx) :graph(g) ,idx(init_idx) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
    
    /** Returns begin of nodes list iterator
  */
  edge_iterator edge_begin() const {
      return edge_iterator(this,0);
  }
   
  /** Returns end of nodes list iterator
  */
  edge_iterator edge_end() const {
      return edge_iterator(this,edge_list.size());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  using Vector = std::vector<size_type>;
  using Pair = std::pair<size_type,size_type>;
  using Set = std::set<size_type>;
  std::vector<Point> point_list;
  std::vector<Edge> edge_list;
    
  std::vector<Node> node_list;
  std::map<Pair,Edge> edge_map;
  std::vector<Set> adjList;
    
  std::vector<node_value_type> node_values;
};

#endif // CME212_GRAPH_HPP
