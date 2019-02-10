#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

using namespace std;

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph {
private:
  // Predeclare the internal struct
  struct internal_node;
  
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
  
  /** Type of node value for this graph. */
  using node_value_type = V;
  
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
      graph_ = nullptr;
      node_ind_ = size_type(-1);
    }
    
    /** Return this node's position. */
    const Point& position() const {
      return fetch().position;
    }
    
    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return node_ind_;
    }
    
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Set this node's value. */
    V& value() {
      return fetch().value;
    }
    
    /** Return this node's value. */
    const V& value() const {
      return fetch().value;
    }
    
    /** Return this node's degree (number of edges incident to it). */
    size_type degree() const {
      return graph_->node_neighbors_[node_ind_].size();
    }
    
    /** Return an incident iterator that points to the first edge in our
     * node neighbors vector for this node.
     *
     * Complexity: O(1).
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, node_ind_,
                              graph_->node_neighbors_[node_ind_].begin());
    }
    
    /** Return an incident iterator that points to one past the last edge
     * in our node neighbors vector.
     *
     * Complexity: O(1).
     */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, node_ind_,
                              graph_->node_neighbors_[node_ind_].end());
    }
    
    
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return ((index() == n.index()) && (graph_ == n.graph_));
    }
    
    
    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     *
     * Here, we define the ordering by node indices.
     */
    bool operator<(const Node& n) const {
      return (index() < n.index());
    }
    
  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    
    Graph* graph_;
    size_type node_ind_;
    
    /** Private Constructor for Node */
    Node(const Graph* graph, size_type node_ind)
    : graph_(const_cast<Graph*>(graph)), node_ind_(node_ind) {
    }
    
    /** Helper method to return the appropriate node. */
    internal_node& fetch() const {
      return graph_->nodes_[index()];
    }
  };
  
  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes_.size();
  }
  
  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }
  
  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const V& value = V()) {
    size_type new_index = num_nodes();
    internal_node new_node = {};
    new_node.index = new_index;
    new_node.position = position;
    new_node.value = value;
    nodes_.push_back(new_node);
    node_neighbors_.push_back(vector<size_type>());
    return Node(this, new_index);
  }
  
  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.index() < size() && n.graph_ == this) {
      return true;
    } else {
      return false;
    }
  }
  
  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < size());
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
      graph_ = nullptr;
      node1_ind_ = size_type(-1);
      node2_ind_ = size_type(-1);
      edge_ind_ = size_type(-1);
    }
    
    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, node1_ind_);
    }
    
    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, node2_ind_);
    }
    
    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((edge_ind_ == e.edge_ind_) && (graph_ == e.graph_));
    }
    
    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     *
     * Here, ordering is defined by the edge indicies.
     */
    bool operator<(const Edge& e) const {
      return (edge_ind_ < e.edge_ind_);
    }
    
  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    
    Graph* graph_;
    size_type node1_ind_;
    size_type node2_ind_;
    size_type edge_ind_;
    
    /** Private Constructor for Edge */
    Edge(const Graph* graph, size_type node1_ind, size_type node2_ind,
         size_type edge_ind) : graph_(const_cast<Graph*>(graph)),
    node1_ind_(node1_ind), node2_ind_(node2_ind), edge_ind_(edge_ind) {
    }
  };
  
  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_map_.size();
  }
  
  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());
    
    // get node_pair value for current edge index i
    set<size_type> node_pair = edges_map_flipped_.at(i);
    
    // use iterator to access elements in set
    set<size_type>::iterator it = node_pair.begin();
    advance(it, 0);
    size_type nodea_index = *it;
    advance(it, 1);
    size_type nodeb_index = *it;
    
    return Edge(this, nodea_index, nodeb_index, i);
  }
  
  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert(a.graph_ == this && b.graph_ == this);
    
    size_type nodea_index = a.index();
    size_type nodeb_index = b.index();
    
    set<size_type> node_pair;  // create set with node a and b indicies
    node_pair.insert(nodea_index);
    node_pair.insert(nodeb_index);
    
    // check if this node_pair is a key in edges_map_
    if (edges_map_.find(node_pair) == edges_map_.end()) {
      return false;   // key not found
    } else {
      return true;    // key found
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
    assert((a == b) == false);
    assert(a.graph_ == this && b.graph_ == this);
    
    size_type nodea_index = a.index();
    size_type nodeb_index = b.index();
    
    set<size_type> node_pair;  // create set with node a and b indicies
    node_pair.insert(nodea_index);
    node_pair.insert(nodeb_index);
    
    if (has_edge(a, b) == true) {
      // return the current edge if it already exists
      return Edge(this, nodea_index, nodeb_index, edges_map_[node_pair]);
    } else {
      size_type edge_index = edges_map_.size();
      edges_map_[node_pair] = edge_index;  // add key, val pair for new edge
      edges_map_flipped_[edge_index] = node_pair; // add flipped info to copy map
      // add node b to neighbor list of a and vise vera
      node_neighbors_[nodea_index].push_back(nodeb_index);
      node_neighbors_[nodeb_index].push_back(nodea_index);
      return Edge(this, nodea_index, nodeb_index, edge_index);
    }
  }
  
  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_ = {};
    node_neighbors_ = {};
    edges_map_.clear();
    edges_map_flipped_.clear();
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
      graph_ = nullptr;
      node_it_ = nodes_.begin();
    }
   
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Return the node that this node iterator points to */
    Node operator*() const {
      return Node(graph_, (*node_it_).index);
    }
    
    /** Return a node iterator that points to the node one past the node that
     * this node iterator points to */
    NodeIterator& operator++() {
      ++node_it_;
      return *this;
    }
    
    /** Test whether this node iterator and @a nit are equal.
     *
     * Equal node iterators represent the same pointer to a node.
     */
    bool operator==(const NodeIterator& nit) const {
      return ((node_it_ == nit.node_it_) && (graph_ == nit.graph_));
    }
    
  private:
    friend class Graph;
    
    Graph* graph_;
    typename vector<internal_node>::const_iterator node_it_;
    
    /** Private Constructor for NodeIterator */
    NodeIterator(const Graph* graph,
                 typename vector<internal_node>::const_iterator node_it) :
    graph_(const_cast<Graph*>(graph)), node_it_(node_it) {
    }
  };
  
  // Supply definitions AND SPECIFICATIONS for:
  
  /** Return a node iterator that points to the first node in our node vector.
   *
   * Complexity: O(1).
   */
  node_iterator node_begin() const {
    return NodeIterator(this, nodes_.begin());
  }
  
  /** Return a node iterator that points to one past the last node in our node vector.
   *
   * Complexity: O(1).
   */
  node_iterator node_end() const {
    return NodeIterator(this, nodes_.end());
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
      graph_ = nullptr;
      node1_ind_ = size_type(-1);
      inc_it_ = nodes_.begin();
    }
    
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Return the edge that this incident iterator points to */
    Edge operator*() const {
      set<size_type> node_pair;  // create set with node a and b indicies
      node_pair.insert(node1_ind_);
      node_pair.insert(*inc_it_);
      
      size_type edge_index = graph_->edges_map_[node_pair];
    
      return Edge(graph_, node1_ind_, *inc_it_, edge_index);
    }
    
    /** Return an incident iterator that points to the edge one past the edge that
     * this incident iterator points to */
    IncidentIterator& operator++() {
      ++inc_it_;
      return *this;
    }
    
    /** Test whether this incident iterator and @a iit are equal.
     *
     * Equal incident iterators represent the same pointer to an edge.
     */
    bool operator==(const IncidentIterator& iit) const {
      return ((inc_it_ == iit.inc_it_) && (graph_ == iit.graph_));
    }
    
  private:
    friend class Graph;

    Graph* graph_;
    size_type node1_ind_;
    typename vector<size_type>::const_iterator inc_it_;
    
    /** Private Constructor for IncidentIterator */
    IncidentIterator(const Graph* graph, size_type node1_ind,
                     typename vector<size_type>::const_iterator inc_it):
    graph_(const_cast<Graph*>(graph)), node1_ind_(node1_ind), inc_it_(inc_it) {
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
      graph_ = nullptr;
      edge_it_ = edges_map_.begin();
    }
    
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Return the edge that this edge iterator points to */
    Edge operator*() const {
      set<size_type> node_pair = edge_it_->first;
      size_type edge_index = edge_it_->second;
      
      // use iterator to access elements in set
      set<size_type>::iterator it = node_pair.begin();
      advance(it, 0);
      size_type nodea_index = *it;
      advance(it, 1);
      size_type nodeb_index = *it;
      
      return Edge(graph_, nodea_index, nodeb_index, edge_index);
    }
    
    /** Return an edge iterator that points to the edge one past the edge that
     * this edge iterator points to */
    EdgeIterator& operator++() {
      ++edge_it_;
      return *this;
    }
    
    /** Test whether this edge iterator and @a eit are equal.
     *
     * Equal edge iterators represent the same pointer to an edge.
     */
    bool operator==(const EdgeIterator& eit) const {
      return ((edge_it_ == eit.edge_it_) && (graph_ == eit.graph_));
    }
    
  private:
    friend class Graph;
    
    Graph* graph_;
    map<set<size_type>, size_type>::const_iterator edge_it_;
    
    /** Private Constructor for EdgeIterator */

    EdgeIterator(const Graph* graph, map<set<size_type>,
                 size_type>::const_iterator edge_it) :
    graph_(const_cast<Graph*>(graph)), edge_it_(edge_it) {
    }
  };
  
  // Supply definitions AND SPECIFICATIONS for:
  
  /** Return an edge iterator that points to the first edge in our edge map.
   *
   * Complexity: O(1).
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, edges_map_.begin());
  }
  
  /** Return an edge iterator that points to one past the last edge in our edge map.
   *
   * Complexity: O(1).
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, edges_map_.end());
  }
  
private:
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  
  // internal type for graph nodes
  struct internal_node {
    size_type index;
    Point position;
    V value;
  };
  
  vector<internal_node> nodes_;
  // the 0th vector in node_neighbors_ holds all the adjacent nodes to node 0
  vector<vector<size_type>> node_neighbors_;
  //edges_map_ contains: key = pair of node a and b indicies, val = unique edge index
  map<set<size_type>, size_type> edges_map_;
  map<size_type, set<size_type>> edges_map_flipped_;  // keys and vals flipped here
};

#endif // CME212_GRAPH_HPP
