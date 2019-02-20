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
template <typename V, typename E>
class Graph {
 private:
  struct NodeInfo;
  struct EdgeInfo;
  
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** type of Node value. */
  using node_value_type = V ;

  using edge_value_type = E ;

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
      return graph_n->nodes[n_uid].position;
    }

    /** Return this node's position, modifiable. */
    Point& position() {
      return graph_n->nodes[n_uid].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_n->nodes[n_uid].idx;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_n->has_node(n) && index() == n.index());
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
      if (graph_n != n.graph_n)
        return true;
      return (n_uid < n.n_uid);
    }

    /** Return this node's value. */
    node_value_type& value() {
      return graph_n->nodes[n_uid].node_value;
    }

    /** Return this node's value. */
    const node_value_type& value() const {
      return graph_n->nodes[n_uid].node_value;
    }

    /** Return this node's degree. */
    size_type degree () const {
      return graph_n->adj_uid[n_uid].size();
    }

    /** Return the incident_iterator, as the first edge. */
    incident_iterator edge_begin () const {
      return incident_iterator(graph_n, index(), 0);
    }

    /** Return the incident_iterator, as the (last+1) edge. */
    incident_iterator edge_end () const {
      return incident_iterator(graph_n, index(), 
                               graph_n->adj_uid[n_uid].size());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    graph_type* graph_n;
    size_type n_uid;
    Node(const graph_type* graph, size_type i)
      : graph_n(const_cast<graph_type*>(graph)), n_uid(i) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return ni2nu.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return ni2nu.size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, 
                const node_value_type& node_value = node_value_type()) {
    NodeInfo new_node{position, (unsigned int)ni2nu.size(), node_value};
    nodes.push_back(new_node);
    ni2nu.push_back(nodes.size() - 1);
    adj_uid.push_back(std::vector<id_pair> ());
    return Node(this, nodes.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (this == n.graph_n &&
            n.n_uid == ni2nu[n.index()]);
  }

  /** Return the node ith index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < ni2nu.size());
    return Node(this, ni2nu[i]);
  }

  /** Remove the given node and all the edges that are incident to it.
   * @param[in] n The target node for removal.
   * @result 0 or 1 (0 means @a n is not in the graph, 
   *                 1 means @a n is in the graph).
   * @post if @a n is not in the graph, 0 is returned;
   *       if @a n is in the graph, ni2nu[@a n.index()]=ni2nu[ni2nu.size()-1],
   *              the last element of ni2nu is pop_backed,
   *              new ni2nu size = old ni2nu size - 1,
   *              all the edges that are incident to @a n are removed,
   *              1 is returned.
   * Complexity: O(num_local_edges * num_edges()).
   * Can invalidate parts of containers related to nodes.
   * for ni2nu, [n.index(), end()) would be invalidated.
   * for nodes, the element with position ni2nu[n.index()] is modified
   * 
   */
  size_type remove_node(const Node& n) {
    if (!has_node(n))
      return 0;
    size_type n_uid = ni2nu[n.index()];
    std::vector<id_pair> n_edges = adj_uid[n_uid];

    // remove edges that are incident to the node
    for (unsigned int i = 0; i < n_edges.size(); ++i) {
      auto no_use = remove_edge(n, Node(this, n_edges[i].n_uid));
      (void) no_use;
    }
   
    ni2nu[n.index()] = ni2nu[ni2nu.size() - 1];
    ni2nu.pop_back();
    nodes[ni2nu[n.index()]].idx = n.index();
    return 1;
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
      graph_e = nullptr;
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_e, n1_uid);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_e, n2_uid);
    }

    double length() const {
      return norm(node1().position() - node2().position());
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
      if (graph_e != e.graph_e)
        return true;
      return ((node1().index() + node2().index()) <
              (e.node1().index() + e.node1().index()));
    }

    edge_value_type& value() {
      return graph_e->edges[e_uid].edge_value;
    }

    const edge_value_type& value() const {
      return graph_e->edges[e_uid].edge_value;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    size_type  n1_uid;
    size_type  n2_uid;
    size_type  e_uid;
    graph_type* graph_e;
    Edge(const graph_type* graph, size_type n1, size_type n2, size_type i)
      : graph_e(const_cast<graph_type*>(graph)) {
      assert(n1 < graph_e->nodes.size() && n2 < graph_e->nodes.size());
      n1_uid = n1;
      n2_uid = n2;
      e_uid  = i;
    }
  };
  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return ei2eu.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < ei2eu.size());
    return Edge(this, edges[ei2eu[i]].n1_uid, edges[ei2eu[i]].n2_uid, ei2eu[i]);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert (has_node(a) && has_node(b));
    size_type a_uid = ni2nu[a.index()];
    size_type b_uid = ni2nu[b.index()];
    std::vector<id_pair> a_edge = adj_uid[a_uid];

    for (unsigned int i = 0; i < a_edge.size(); ++i) {
      if (a_edge[i].n_uid == b_uid)
        return true;
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
  Edge add_edge(const node_type& a, const node_type& b) {
    assert (has_node(a) && has_node(b));
    assert (! (a == b));
    size_type a_uid = ni2nu[a.index()];
    size_type b_uid = ni2nu[b.index()];
    std::vector<id_pair> a_edge = adj_uid[a_uid];
    for (unsigned int i = 0; i < a_edge.size(); ++i) {
      if (a_edge[i].n_uid == b_uid)
        return Edge(this, a_uid, b_uid, a_edge[i].e_uid);
    }

    edge_type new_edge(this, a_uid, b_uid, edges.size());
    edges.push_back(EdgeInfo{a_uid, b_uid, edge_value_type()});
    ei2eu.push_back(edges.size() - 1);
    
    adj_uid[a_uid].push_back(id_pair{b_uid, (unsigned int)edges.size() - 1});
    adj_uid[b_uid].push_back(id_pair{a_uid, (unsigned int)edges.size() - 1});
  
    return new_edge;
  }

 /** Remove an edge with 2 given nodes in the graph
  * @param[in] n1 One node of the target edge
  * @param[in] n2 The other node of the target edge
  * @result 0 or 1 (0 means the edge is not in the graph
  *                 1 means the edge is successfully removed)
  * @post if @a n1 or @a n2 or the edge constructed by @a n1 and @a n2 is not
  *          in the graph, 0 is returned;
  *       if the edge constructed by @a n1 and @a n2 is in the graph, 
  *          1 is returned,
  *          the connection of nodes of @a e is deleted in adj_uid,
  *          ei2eu[k] = e.e_uid, 
  *             for any i, 0 <= i < k, ei2eu[i] is unchanged,
  *             for any j, k <= j < ei2eu.size()-1, ei2eu[i] = ei2eu[i+1]
  *          ei2eu[ei2eu.size()-1] is deleted.
  * @post for 0 <= i < j < ei2eu.size(), ei2eu[i] < ei2eu[j].
  * Complexity: O(num_edges())
  * Can invalidate edge indices of the index of @a e and all after that.
  */
  size_type remove_edge(const Node& n1, const Node& n2) {
    if (!has_edge(n1, n2))
      return 0;
    return remove_edge(add_edge(n1, n2));
  }

 /** Remove a given edge in the graph
  * @param[in] e The edge that should be removed
  * @result 0 or 1 (0 means the edge is not in the graph
  *                 1 means the edge is successfully removed)
  * @post if @a e is not in the graph, 0 is returned;
  *       if @a e is in the graph, 1 is returned,
  *          the connection of nodes of @a e is deleted in adj_uid,
  *          ei2eu[k] = e.e_uid, 
  *             for any i, 0 <= i < k, ei2eu[i] is unchanged,
  *             for any j, k <= j < ei2eu.size()-1, ei2eu[i] = ei2eu[i+1]
  *          ei2eu[ei2eu.size()-1] is deleted.
  * @post for 0 <= i < j < ei2eu.size(), ei2eu[i] < ei2eu[j].
  * Complexity: O(num_edges())
  * Can invalidate edge indices of the index of @a e and all after that.
  */
  size_type remove_edge(const Edge& e) {
    size_type e_uid = e.e_uid;
    size_type n1_uid = ni2nu[e.node1().index()];
    size_type n2_uid = ni2nu[e.node2().index()];
    auto nedges = adj_uid[n1_uid].size();
    for (unsigned int i = 0; i < nedges; ++i) {
      if (adj_uid[n1_uid][i].n_uid==n2_uid && adj_uid[n1_uid][i].e_uid==e_uid){
        adj_uid[n1_uid][i].n_uid = adj_uid[n1_uid][nedges-1].n_uid;
        adj_uid[n1_uid][i].e_uid = adj_uid[n1_uid][nedges-1].e_uid;
        adj_uid[n1_uid].pop_back();
      }
    }
    if (nedges == adj_uid[n1_uid].size()) // e is not in the graph
      return 0; 
    auto nedges2 = adj_uid[n2_uid].size();
    for (unsigned int i = 0; i < nedges2; ++i) {
      if (adj_uid[n2_uid][i].n_uid==n1_uid && adj_uid[n2_uid][i].e_uid==e_uid){
        adj_uid[n2_uid][i].n_uid = adj_uid[n2_uid][nedges2-1].n_uid;
        adj_uid[n2_uid][i].e_uid = adj_uid[n2_uid][nedges2-1].e_uid;
        adj_uid[n2_uid].pop_back();
      }
    }

    for (unsigned int i = 0; i < ei2eu.size() - 1; ++i) {
      if (ei2eu[i] >= e.e_uid)
        ei2eu[i] = ei2eu[i+1];
    }
    ei2eu.pop_back();
    return 1;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    adj_uid.clear();
    edges.clear();
    nodes.clear();
    ei2eu.clear();
    ni2nu.clear();
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

    /** Construct the node that the current iterator points to. */
    Node operator*() const {
      return Node(graph_ni, graph_ni->ni2nu[curr_idx]);
    }

    /** Construct a node iterator points to the next one. */
    NodeIterator& operator++() {
      curr_idx++;
      return *this;
    }
 
    /** Equivalence check. */
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
  };


  /** construct a node iterator points to the node with index of 0. */
  node_iterator node_begin() const {
    return node_iterator(this, 0);
  }

  /** construct a node iterator points to the (last + 1) node */
  node_iterator node_end() const {
    return node_iterator(this, ni2nu.size());
  }

  /** Remove the node cooreponding to the given node iterator
   *        and all the edges that are incident to it.
   * @param[in] n_it The node iterator of target node for removal.
   * @result 0 or 1 (0 means @a (*n_it) is not in the graph, 
   *                 1 means @a (*n_it) is in the graph).
   * @post if @a n is not in the graph, 0 is returned;
   *       if @a n is in the graph, ni2nu[@a n.index()]=ni2nu[ni2nu.size()-1],
   *              the last element of ni2nu is pop_backed,
   *              new ni2nu size = old ni2nu size - 1,
   *              all the edges that are incident to @a n are removed,
   *              1 is returned.
   * Complexity: O(num_local_edges * num_edges()).
   * Can invalidate iterators [n.index(), end()).
   * With the implementation of remove_node(const Node&), only iterators at 
   * n.index() and end() - 1 would be invalidated.
   * 
   */
  node_iterator remove_node(node_iterator n_it) {
    auto node = *n_it;
    auto idx = node.index();
    remove_node(node);
    return node_iterator(this, idx);
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

    /** Construct the edge that the current iterator points to. */
    Edge operator*() const {
      size_type n1 = graph_ii->ni2nu[node_idx];
      return Edge(graph_ii, n1, graph_ii->adj_uid[n1][curr_edge].n_uid, 
                  graph_ii->adj_uid[n1][curr_edge].e_uid);
    }

    /** Construct an incident iterator points to the next one. */
    IncidentIterator& operator++() {
      curr_edge++;
      return *this;
    }

    /** Equaivalence check. */
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

    /** Construct the edge that current iterator points to. */
    Edge operator*() const {
      size_type n1 = graph_ei->edges[graph_ei->ei2eu[edge_idx]].n1_uid;
      size_type n2 = graph_ei->edges[graph_ei->ei2eu[edge_idx]].n2_uid;
      return Edge(graph_ei, n1, n2, graph_ei->ei2eu[edge_idx]);
    }

    /** Construct a iterator that points to the next one. */
    EdgeIterator& operator++() {
      edge_idx++;
      return *this;
    }

    /** Equivalence check. */
    bool operator==(const EdgeIterator& ei) const {
      return (graph_ei == ei.graph_ei &&
              edge_idx == ei.edge_idx);
    }

   private:
    friend class Graph;
    graph_type* graph_ei;
    size_type   edge_idx;
    EdgeIterator(const graph_type* graph, size_type idx)
      : graph_ei(const_cast<graph_type*>(graph)), edge_idx(idx){
    }

  };

  /** Construct an edge iterator points to the edge with index 0. */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /** Construct an edge iterator points to the (last + 1) edge. */
  edge_iterator edge_end() const {
    return EdgeIterator(this, ei2eu.size());
  }

 /** Remove an edge coorespoping to the given edge iterator in the graph.
  * @param[in] e_it The edge iterator of the target edge for removal.
  * @result edge iterator of this graph with the same edge index.
  * @pre @a e_it is the iterator in this graph.
  * @pre @a e_it.edge_idx < ei2eu.size().
  * @post new ei2eu.size() <= old ei2eu.size()
  * Complexity: O(num_edges()), depends on remove_edge(*e_it).
  * Can invalidate edge indices of the index of @a e_it and all after that.
  */
  edge_iterator remove_edge(edge_iterator e_it) {
    assert(this == e_it.graph_ei);
    assert(e_it.edge_idx < ei2eu.size());
    size_type idx = e_it.edge_idx;
    remove_edge(*e_it);
    return edge_iterator(this, idx);
  }

 private:

  struct NodeInfo {
    Point                                     position;
    size_type                                 idx;
    node_value_type                           node_value;
  };

  struct EdgeInfo {
    size_type                                 n1_uid;
    size_type                                 n2_uid;
    edge_value_type                           edge_value;
  };

  struct id_pair {
    size_type                                 n_uid;
    size_type                                 e_uid;
  };

  std::vector<EdgeInfo>                       edges;   // order in e_uid
  std::vector<NodeInfo>                       nodes;   // order in n_uid
  std::vector<size_type>                      ei2eu;   // e_idx to e_uid
  std::vector<size_type>                      ni2nu;   // n_idx to n_uid
  std::vector<std::vector<id_pair> >        adj_uid; // check edge existence
};

#endif // CME212_GRAPH_HPP
