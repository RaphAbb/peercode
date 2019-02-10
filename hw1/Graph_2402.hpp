#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <map>
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
class Graph : private totally_ordered<Graph<V>> {
  struct internal_edge;
  struct internal_node;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V>;

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
  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //

  using node_value_type = V;

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node{
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
    }

    /** Return this node's position. */
    const Point& position() const {
      return (gph->point_list[this->idx]).point_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return size_type(this->idx);
    }

    /**
     * @breif Return this node's value, of type V, by reference
     * @return node value, type V&
     */
    node_value_type& value(){
      internal_node& n = gph->point_list[this->idx];
      return n.val_;
    }

  
    /**
     * @brief Return this node's value by constant reference
     * @return node value, type const V&
     */
    const node_value_type& value() const{
      return const_cast<node_value_type&>((gph->point_list[this->idx]).val_);
    }

    /** 
     * @brief Method to return the degree of this node, the number of
     *        edges incident to the node 
     * @return degree of this node, type size_type
     */
    size_type degree() const {
      size_type deg = 0;
      for (incident_iterator i = edge_begin(); i!=edge_end(); ++i)
        deg++;
      return deg;
    }

    /** 
     * @brief Method to return an iterator to the first edge incident to the 
     *        current node 
     * @return iterator pointing to first edge leaving this node
     */
    incident_iterator edge_begin() const{
      return incident_iterator(gph,idx,(gph->adjacency_list[idx]).begin());
    }

    /** 
     * @brief Method to return an iterator pointing to one past the end of the
     *        edges incident to the current node 
     * @return iterator pointing to last edge leaving this node
     */
    incident_iterator edge_end() const{
      return incident_iterator(gph, idx, (gph->adjacency_list[idx]).end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (n.index() == this->idx and n.gph == this->gph)
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
      internal_node a = gph->point_list[this->idx];
      if (n.index() < a.idx_)
        return true;
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // pointer back to node's graph
    graph_type* gph;
    // store node's index
    size_type idx;

    // private node constructor that sets node data
    Node(const graph_type* graph, int index) {
      this->gph = const_cast<graph_type*>(graph);
      this->idx = index;
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
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
  Node add_node(const Point& position, 
                const node_value_type& value=node_value_type()){
    Node new_node = Node(this, point_list.size());
    node_value_type v = value;
    internal_node a = internal_node(position, point_list.size(),v);
    point_list.push_back(a);

    // add empty map to graph adjacency list
    std::map<size_type, size_type> m;
    adjacency_list.push_back(m);

    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.index() < this->point_list.size())
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
    return Node(this, i);
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
    }

    /** Return node 1 of this Edge, or node 2 if switch_ == true */
    Node node1() const {
      if (switch_)
        return Node(this->gph, fetch().node2_idx);  
      return Node(this->gph, fetch().node1_idx);
    }

    /** Return node 2 of this Edge, or node 1 if switch_ == true */
    Node node2() const {
      if (switch_)
        return Node(this->gph, fetch().node1_idx);
      return Node(this->gph, fetch().node2_idx);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // get edge nodes to avoid looking up multiple times
      int e_node1 = e.node1().index();
      int e_node2 = e.node2().index();
      int t_node1 = fetch().node1_idx;
      int t_node2 = fetch().node2_idx;

      if ((e_node1==t_node1 and e_node2==t_node2) 
        or (e_node1==t_node2 and e_node2==t_node1))
        return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (this->idx < e.idx)
        return true;
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph; 
    // store pointer to graph the edge belongs to 
    const graph_type* gph;
    // store the index of the edge
    size_type idx;
    // switch node1() and node2() if necessary for IncidentIterator
    bool switch_;

    // private constructor that sets edge values
    Edge(const Graph* graph, int index, bool s = false){
      this->gph = graph;
      this->idx = index;
      this->switch_ = s;
    }

    // helper function to get internal edge object for this edge
    const internal_edge& fetch() const {
      return gph->edge_list[idx];
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_count;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if (adjacency_list[a.index()].count(b.index())!=0)
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
    // check if edge is already in graph
    if (this->has_edge(a,b)){
      size_type its_index = adjacency_list[a.index()][b.index()];
      return Edge(this,its_index);
    }

    // create new internal edge object
    internal_edge new_edge = internal_edge(a,b, edge_count);
    this->edge_list.push_back(new_edge);
    edge_count++;

    // add new edge to adjacency list
    adjacency_list[a.index()][b.index()] = edge_count-1;
    adjacency_list[b.index()][a.index()] = edge_count-1;

    return Edge(this, edge_count-1);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // clears edge list, point list, adjacency list, resets edge count
    this->edge_list.clear();
    this->point_list.clear();
    this->adjacency_list.clear();
    edge_count = 0;
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

    /**
     * @brief Method to access node object from iterator 
     * @return Node object for the node that the iterator currently 
     *         references
     */
    Node operator*() const {
      return Node(gph_, idx_);
    }

    /** 
     * @brief Method to increment iterator to the next node
     * @return a node iterator referencing the next element
     */
    NodeIterator& operator++(){
      idx_++;
      return *this;
    }

    /** 
     * @brief Method to compare if two iterators refer to the same node
     * @param[in] const NodeIterator& ni, the object to compare to the 
     *            current NodeIterator
     * @return true if iterators refer to the same graph node, false otherwise
     * 
     * Checks if the iterators refer to the same graph and the same node index
     */
    bool operator==(const NodeIterator& ni) const {
      if (idx_ == ni.idx_ and gph_ == ni.gph_)
        return true;
      return false;
    }

   private:
    friend class Graph;

    // pointer back to graph
    const graph_type* gph_;
    // current index that the iterator refers to 
    int idx_;

    /** Private constructor for a node iterator, instantiating values */
    NodeIterator(const graph_type* graph, int index1){
      gph_ = graph;
      idx_ = index1; 
    }
  };

  /** 
   * @brief Method to return an iterator pointing to the first node in the 
   *        graph
   * @return a node_iterator pointing to the first node in the graph 
   */
  node_iterator node_begin() const {
    return NodeIterator(this,0);
  }

  /**
   * @brief Method to return an interator pointing to one past the last node 
   *        in the graph
   * @return a node_iterator pointing to one past the last node in the graph
   */
  node_iterator node_end() const {
    return NodeIterator(this,this->num_nodes());
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
    }

    /** 
     * @brief Method to dereference the incident iterator
     * @return Edge e, an edge object matching the edge the iterator is
     *         pointing to
     * @post e.node1() returns the base node of the IncidentIterator
     */
    Edge operator*() const {
      Edge e = Edge(gph_, idx_);
      // if internal edge has nodes flipped, fix it
      if (node_idx_ == e.node2().index())
        return Edge(gph_, idx_,true);
      // otherwise just return the normal edge
      return e;
    }

    /** 
     * @brief Method to increment the incident iterator to the next edge 
     * @return an IncidentIterator that refers to the next incident edge
     * Increments the adjacency list map and updates the edge index accordingly
     */
    IncidentIterator& operator++(){
      map_iter++;
      idx_ = map_iter->second;
      return *this;
    }

    /** 
     * @brief Method to compare if two incident iterators refer to the same
     *        point. Checks if the base node is the same, and that the 
     *        iterator is pointing to the same edge leaving that node.
     * @param[in] const IncidentIterator iter, the object to compare to the
     *            current IncidentIterator
     * @return true if the iterators have the same base node and refer to the
     *         same edge, false otherwise.
     */
    bool operator==(const IncidentIterator& iter) const {
      if ((node_idx_ == iter.node_idx_) and (idx_ == iter.idx_))// and (map_iter == iter.map_iter))
        return true;
      return false;
    }

   private:
    friend class Graph;
    // pointer back to the original graph
    const graph_type* gph_;
    // index of the edge currently referenced by the iterator
    size_type idx_;
    // index of the node the desired edges are incident to 
    size_type node_idx_;
    // iterator for the adjacency list map
    std::map<size_type, size_type>::iterator map_iter; 

    /** Private constructor for the incident iterator */
    IncidentIterator(const graph_type* graph, int node_a, std::map<size_type, 
                      size_type>::iterator edge_iter)
    : gph_(graph), node_idx_(node_a), map_iter(edge_iter){
      idx_ = (*map_iter).second;
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
    EdgeIterator() {
    }

    /** 
     * @brief Method to dereference an edge_iterator 
     * @return an edge object matching the edge that the iterator is pointing 
     *         to. 
     */
    Edge operator*() const {
      return Edge(ei_gph_, ei_idx_);
    }

    /** 
     * @brief Method to increment an edge_iterator to the next edge 
     * @return an EdgeIterator that refers to the next edge
     */
    EdgeIterator& operator++(){
      ei_idx_++;
      return *this;
    }

    /**
     * @brief Method to check if two edge iterators are the same. Checks 
     *        if both iterators refer to both the same graph and the same
     *        edge. 
     * @param[in] const EdgeIterator& e, the EdgeIterator object to compare to
     *            the current EdgeIterator
     * @return true if both iterators refer to the same graph and edge,
     *         false otherwise
     */
    bool operator==(const EdgeIterator& e) const {
      if (ei_gph_==e.ei_gph_ and ei_idx_==e.ei_idx_)
        return true;
      return false;
    }

   private:
    friend class Graph;
    // pointer back to the graph
    const graph_type* ei_gph_;
    // index of the current edge
    int ei_idx_;


    /** Private constructor for an edge iterator */
    EdgeIterator(const graph_type* graph, int index){
      ei_gph_ = graph;
      ei_idx_ = index;
    }
  };

  /**
   * @brief Method to access an EdgeIterator that refers to the first edge
   *        in the graph.
   * @return an EdgeIterator object referring to the first edge in the graph
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /**
   * @brief Method to access an EdgeIterator that refers to the one past the
   *        last edge in in the graph.
   * @return an EdgeIterator object referring to one past the last edge in
   *         the graph
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, edge_count);
  }

 private:

  // adjacency map to store graph edge information
  std::vector<std::map<size_type, size_type>> adjacency_list;
  // keep a count of current number of edges
  size_type edge_count = 0;

  // internal edge object to store information about a point
  struct internal_edge {
    // store node indices and edge index
    size_type node1_idx;
    size_type node2_idx;
    size_type idx;
    // constructor to set internal edge values
    internal_edge(const Node& a, const Node& b, size_type index){
      node1_idx = a.index();
      node2_idx = b.index();
      idx = index;
    }
  };

  struct internal_node {
    // index, point, and value
    size_type idx_;
    const Point point_;
    V val_;
    internal_node(const Point& p, size_type i, V& v)
      : idx_(i), point_(p), val_(v) {}
  };

  // list of node data containers
  std::vector<internal_node> point_list;
  // list of edge data containers
  std::vector<internal_edge> edge_list;

};

Graph <double> graph1;
Graph <int> graph2;

#endif // CME212_GRAPH_HPP
