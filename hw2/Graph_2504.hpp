#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <set>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

/*
* Again, I have a fair amount of incomplete work on this assignment. 
* I have not yet worked on Problem 6 (node/edge removal 
* and the corresponding constraint), and have some bugs in the earlier
* questions too. Will give it a shot later and see if I can take a 
* late day and resubmit, and if not will push these to the next HW.
*/


/** @class Graph
 * @brief A template for 3D undirected graphs.
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
  using node_value_type = V;
  using edge_value_type = E;
  using graph_type = Graph<V, E>;

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

 private:

  // I really should have changed this...
  std::map<unsigned, const Point> node_map;
  std::map<unsigned, node_value_type> node_value_map;
  std::map<unsigned, std::pair<unsigned, unsigned>> edge_map;
  std::map<std::pair<unsigned, unsigned>, unsigned> inv_edge_map;
  std::map<unsigned, std::set<unsigned>> adj_list;    

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  public:
  Graph() {
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
    Graph* graph;
    size_type idx;

    Node() {
    }


    Node(const Graph* g, size_type i): idx(i)
    {
      graph = const_cast<Graph*> (g);
    }

    Node(const Node& n):idx(n.idx)
    {
      graph = n.graph;
    }

    /** Return this node's position. */
    Point& position() const {
      return graph->node_map.at(idx);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return idx;
    }

    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    node_value_type& value()
    {
      return graph->node_value_map.at(idx);
    }

    const node_value_type& value() const
    {
      return graph->node_value_map.at(idx);
    }

    size_type degree() const
    {
      return graph->adj_list.at(idx).size();
    }

    incident_iterator edge_begin const
    {
      return incident_iterator(this->graph, idx, graph->adj_list.at(idx).begin());
    }

    incident_iterator edge_end const
    {
      return incident_iterator(this->graph, idx, graph->adj_list.at(idx).end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (n.graph == graph) && (n.idx == idx);
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
      return idx < n.idx;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_map.size();
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

  Node add_node(const Point& position, const node_value_type val=node_value_type()) 
  {    
    const size_type old_size = node_map.size();
    node_map.insert({old_size, position});
    node_value_map.insert({old_size, val});
    return Node(this, old_size);
  }


  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.graph == this) && (n.idx < node_map.size());
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
  class Edge: private totally_ordered<Edge> {
   public:
    const Graph* graph;
    const size_type idx;
    edge_value_type val;
    Edge()
    {
    }
    
    /** Construct an invalid Edge. */    
    Edge(const Graph* g, size_type i, edge_value_type v): idx(i)
    {
      graph = g;
      val = v;
    }

    /* Constructor without specifying edge value */
    Edge(const Graph* g, size_type i): idx(i)
    {
      graph = g;
      val = edge_value_type(0);
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph, graph->edge_map.at(idx).first);  
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph, graph->edge_map.at(idx).second);  
    }

    edge_value_type& value() {
      return val;
    }

    const edge_value_type& value() const {
      return val;
    }

    /** Return the edge length */
    double length() const {
      return norm((graph->edge_map.at(idx).first).positon() - (graph->edge_map.at(idx).second).positon()); 
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((e.graph==graph) && (e.idx==idx));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return idx < e.idx;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    // const Graph* graph;
    // const size_type n_1;
    // const size_type n_2;
    // const size_type edge_idx;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return this->edge_map.size();
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
    const auto& alist = adj_list.find(a.idx);
    return (alist != adj_list.end()) && (alist->second.find(b.idx) != alist->second.end());
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
    if(has_edge(a, b))
    {
      return Edge(this,edge_map.at({a.idx, b.idx}));
    }

    adj_list[a.idx].insert(b.idx);
    adj_list[b.idx].insert(a.idx);
    auto const & orig_size = edge_map.size();
    edge_map[orig_size] = {a.idx, b.idx};
    inv_edge_map[{a.idx, b.idx}] = orig_size;
    inv_edge_map[{b.idx, a.idx}] = orig_size;
    return Edge(this, orig_size);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    node_map.clear();
    edge_map.clear();
    adj_list.clear();
  }


 // private:
  // Commented out the above line
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

    const Graph* graph;
    size_type idx;

    /** Construct an invalid NodeIterator. */
    NodeIterator(size_type i, const Graph* g) {
      idx = i;
      graph = g;
    }

    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const
    Node operator*() const
    {
      return Node(graph, idx);
    }

    NodeIterator& operator++()
    {
      ++idx;
      return *this;
    }

    bool operator==(const NodeIterator& iter) const
    {
      return iter.idx == idx;
    }

   private:
    friend class Graph;
  };

  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  node_iterator node_begin() const
  {
    return node_iterator(0, this);
  }

  node_iterator node_end() const
  {
    return node_iterator(node_map.size()-1, this);
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

    Graph* graph;
    size_type idx;
    std::set<unsigned>::iterator iter;

    /** Construct an invalid IncidentIterator. */
    IncidentIterator(const Graph* g, const size_type i, std::set<unsigned>::iterator it):idx(i) {
      graph = const_cast<Graph*>(g);
      iter = it;
    }


    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const
    Edge operator*() const
    {
      return Edge(graph, graph->inv_edge_map.at({idx, *iter}));
    }

    incident_iterator& operator++()
    {
      ++iter;
      return *this;
    }

    bool operator==(const IncidentIterator& newiter) const
    {
      return (newiter.idx == idx) && (newiter.iter == iter);
    }

   private:
    friend class Graph;
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

    const Graph* graph;
    size_type idx;

    /** Construct an invalid EdgeIterator. */
    EdgeIterator(size_type i, const Graph* g) {
      idx = i;
      graph = g;
    }

    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const
    Edge operator*() const
    {
      return Edge(graph, idx);
    }

    edge_iterator operator++()
    {
      ++idx;
      return *this;
    }

    bool operator==(const EdgeIterator& newiter) const
    {
      return newiter.idx == idx;
    }

   private:
    friend class Graph;
  };

  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  edge_iterator edge_begin() const
  {
    return edge_iterator(0, this);
  }

  edge_iterator edge_end() const
  {
    return edge_iterator(edge_map.size()-1, this);
  }

  //private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
