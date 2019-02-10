#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <iostream>
#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
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
  struct nodeinfo;
  
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** type of Node value. */
  using node_value_type = V ;

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
    nnodes = 0;
    nedges = 0;
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
      graph_n = nullptr;
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_n->nodes[node_index].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return node_index;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_n->has_node(n) && node_index == n.index());
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
      return (node_index < n.index());
    }

    /** Return this node's value. */
    node_value_type& value() {
      return graph_n->nodes[node_index].node_value;
    }

    /** Return this node's value. */
    const node_value_type& value() const {
      return graph_n->nodes[node_index].node_value;
    }

    /** Return this node's degree. */
    size_type degree () const {
      return graph_n->adj_idx[node_index].size();
    }

    /** Return the incident_iterator, as the first edge. */
    incident_iterator edge_begin () const {
      return incident_iterator(graph_n, node_index, 0);
    }

    /** Return the incident_iterator, as the (last+1) edge. */
    incident_iterator edge_end () const {
      return incident_iterator(graph_n, node_index, 
                               graph_n->adj_idx[node_index].size());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    graph_type* graph_n;
    size_type node_index;
    Node(const graph_type* graph, size_type i)
      : graph_n(const_cast<graph_type*>(graph)), node_index(i) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nnodes;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return nnodes;
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& node_value = node_value_type()) {
    node_value_type n_value = node_value;
    nodeinfo new_node{position, nnodes, n_value};
    nodes.push_back(new_node);
    nnodes += 1;
    adj_idx.push_back(std::vector<size_type> ());
    return Node(this, nnodes - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (this == n.graph_n);
  }

  /** Return the node ith index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < nnodes);
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      edge_node1 = nullptr;
      edge_node2 = nullptr;
      graph_e = nullptr;
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_e, edge_node1->idx);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_e, edge_node2->idx);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const edge_type& e) const {
      if ((node1() == e.node1() && node2() == e.node2()) ||
          (node1() == e.node2() && node2() == e.node1()))
        return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const edge_type& e) const {  
      return ((edge_node1->idx   + edge_node2->idx) <
              (e.node1().index() + e.node1().index()));
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    const nodeinfo*  edge_node1;
    const nodeinfo*  edge_node2;
    graph_type* graph_e;
    Edge(const graph_type* graph, const size_type n1, const size_type n2)
      : graph_e(const_cast<graph_type*>(graph)) {
      assert(n1 < graph_e->nnodes && n2 < graph_e->nnodes);
      edge_node1 = &graph_e->nodes[n1];
      edge_node2 = &graph_e->nodes[n2];
    }
  };
  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return nedges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < nedges);
    return edges[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert (has_node(a) && has_node(b));
    size_type a_ind = a.index();
    size_type b_ind = b.index();
    std::vector<size_type> a_edge = adj_idx[a_ind];
    return(std::find(a_edge.begin(), a_edge.end(), b_ind) != a_edge.end());
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
  Edge add_edge(const node_type& a, const node_type& b) {
    assert (has_node(a) && has_node(b));
    assert (! (a == b));
    size_type a_ind = a.index();
    size_type b_ind = b.index();
    if (has_edge(a, b))
      return Edge(this, a_ind, b_ind);

    edge_type new_edge(this, a_ind, b_ind);
    edges.push_back(new_edge);
    
    adj_idx[a_ind].push_back(b_ind);
    adj_idx[b_ind].push_back(a_ind);
  
    nedges += 1;
    return new_edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    adj_idx.clear();
    edges.clear();
    nodes.clear();
    nnodes = 0;
    nedges = 0;
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
      graph_ni = nullptr;
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() const {
      return Node(graph_ni, curr_idx);
    }

    NodeIterator& operator++() {
      curr_idx++;
      return *this;
    }
 
    bool operator==(const NodeIterator& ni) const {
      return (graph_ni == ni.graph_ni && curr_idx == ni.curr_idx);
    }

   private:
    friend class Graph;
    graph_type* graph_ni;
    size_type curr_idx;
    NodeIterator(const graph_type* graph, size_type i)
      : graph_ni(const_cast<graph_type*>(graph)), curr_idx(i) {
    }
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const {
    return node_iterator(this, 0);
  }

  node_iterator node_end() const {
    return node_iterator(this, nnodes);
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
      graph_ii = nullptr;
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {
      return Edge(graph_ii, node_idx, graph_ii->adj_idx[node_idx][curr_edge]);
    }

    IncidentIterator& operator++() {
      curr_edge++;
      return *this;
    }

    bool operator==(const IncidentIterator& ii) const {
      return (graph_ii  == ii.graph_ii &&
              node_idx  == ii.node_idx &&
              curr_edge == ii.curr_edge);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    graph_type* graph_ii;
    size_type node_idx;
    size_type curr_edge;
    IncidentIterator(const graph_type* graph, size_type idx, size_type i)
      : graph_ii(const_cast<graph_type*>(graph)), node_idx(idx), curr_edge(i){
    }
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {
      return graph_ei->edges[edge_idx];
    }

    EdgeIterator& operator++() {
      edge_idx++;
      return *this;
    }

    bool operator==(const EdgeIterator& ei) const {
      return (graph_ei == ei.graph_ei &&
              edge_idx == ei.edge_idx);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_ei;
    size_type   edge_idx;
    EdgeIterator(const graph_type* graph, size_type idx)
      : graph_ei(const_cast<graph_type*>(graph)), edge_idx(idx){
    }

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  edge_iterator edge_end() const {
    return EdgeIterator(this, nedges);
  }

 private:

  struct nodeinfo {
    Point                                     position;
    size_type                                 idx;
    node_value_type                           node_value;
  };


  std::vector<edge_type>                      edges;
  std::vector<nodeinfo>                       nodes;
  std::vector<std::vector<size_type> >        adj_idx; // check edge existence
  size_type                                   nnodes;
  size_type                                   nedges;

};

#endif // CME212_GRAPH_HPP
