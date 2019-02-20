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
 *
 * Citation: Modified data structures in Graph class insprired by peer code #2889
 */
template <typename V, typename E>
class Graph {
private:
  // Predeclare the internal struct
  struct internal_node;
  struct internal_edge;
  
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
  
  /** Type of edge value for this graph. */
  using edge_value_type = E;
  
  
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //
  
  /** Construct an empty graph. */
  Graph() {
  }
  
  /** Default destructor */
  ~Graph() = default;
  
  // Supply definitions AND SPECIFICATIONS for:   // TO DO

  //
  /** Remove node @a n from this graph.
  * @param[in] n     A node to remove.
  * @return 1 if graph had node and was successfuly removed, else
  *         0 if node did not exist in this graph
  *
  * @pre @a n is a valid node.
  * @post has_node(@a n) == false
  * @post If old has_node(@a n), new num_nodes() == old num_nodes() - 1.
  *                              new num_edges() == old num_edges() - n.degree().
  *                              For all nodes i that were adjacent to @a n,
  *                              new i.degree() == old i.degree() - 1 and
  *                              has_edge(@a n, i) == false.
  *       Else,                  new num_nodes() == old num_nodes().
  *                              new num_nodes() == old num_edges().
  *                              new degree() == old degree() for all nodes.
  *
  * Complexity: No more than O(num_nodes()) + complexity of remove_edge function
  *
  */
  size_type remove_node(const Node& n) {

    size_type new_index = n.index();   // new index for replacement node
    size_type n_uid = index2uid_[new_index];
    
    if (has_node(n)) {

      // Iterate through incident edges of n and remove these edges
      vector<Edge> edges_to_remove;
      for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
        Edge e = *ei;
        edges_to_remove.push_back(e);
      }

      for (size_type i = 0; i < edges_to_remove.size(); i++) {
        Edge e = edges_to_remove[i];
        Node n2 = e.node2();   //neighboring node
        remove_edge(n, n2);   // remove edge between n and all its neighbors
      }

      // Swap last node with n
      size_type last_node_uid = index2uid_[num_nodes()-1];
      swap(index2uid_[num_nodes()-1], index2uid_[new_index]);
      uid2index_[last_node_uid] = new_index;   // update index of replacement node
      
      // Remove node n from internal structures
      index2uid_.pop_back();
      nodes_.erase(n_uid);
      uid2index_.erase(n_uid);

      // Iterate through all neighbors of n and remove it as neighhbor in node_neightbors
      vector<size_type> neighbors = node_neighbors_[n_uid];
      for (size_type i = 0; i < neighbors.size(); i++) {
        size_type next_uid = neighbors[i];
        vector<size_type> next_uid_neighbors = node_neighbors_[next_uid];
        vector<size_type>::iterator it = find(next_uid_neighbors.begin(),
                                              next_uid_neighbors.end(), n_uid);
        next_uid_neighbors.erase(it);
        node_neighbors_[next_uid] = next_uid_neighbors;
      }
      node_neighbors_.erase(n_uid);
      
      return 1;  // graph had node and was successfuly removed
    } else {
      return 0;  // node did not exist in this graph
    }
  }


  /** Remove edge between nodes @a n1 and @a n2 from this graph.
   * @param[in] n1     A node of edge to remove.
   * @param[in] n2     The other node of edge to remove.
   * @return 1 if graph had node and was successfuly removed, else
   *         0 if node did not exist in this graph
   *
   * @pre @a n1 and @a n1 are valid nodes.
   *
   * @post has_edge(@a n1, @a n2) == false
   * @post If old has_edge(@a n1, @a n2), new num_edges() == old num_edges() - 1.
   *                                      new n1.degree() == old n1.degree() - 1.
   *                                      new n2.degree() == old n2.degree() - 1.
   *       Else,                          new num_edges() == old num_edges().
   *                                      new n1.degree() == old n1.degree().
   *                                      new n2.degree() == old n2.degree().
   *
   * Complexity: No more than O(num_nodes() + num_edges())
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    if (has_edge(n1, n2)) {
      
      // Find uid of n1 and n2
      size_type n1_uid = index2uid_[n1.index()];
      size_type n2_uid = index2uid_[n2.index()];

      // Find edge index of n1 and n2
      set<size_type> node_pair = {};
      node_pair.insert(n1_uid);
      node_pair.insert(n2_uid);
      size_type edge_index = uid2edge_index_[node_pair];
      
      // Keep track of last edge (to swap with edge from n1 to n2)
      internal_edge last_edge = edges_[num_edges()-1];
      set<size_type> last_node_pair = {};
      last_node_pair.insert(last_edge.node1_uid);
      last_node_pair.insert(last_edge.node2_uid);
      
      // Swap last edge with edge to remove
      swap(edges_[edge_index], edges_[num_edges()-1]);
      uid2edge_index_[last_node_pair] = edge_index;  // update index
      
      // Remove n1 as neighbor of n2
      vector<size_type> n1_neighbors = node_neighbors_[n1_uid];
      vector<size_type>::iterator it1 = find(n1_neighbors.begin(), n1_neighbors.end(), n2_uid);
      n1_neighbors.erase(it1);
      node_neighbors_[n1_uid] = n1_neighbors;

      // Remove n2 as neighbor of n1
      vector<size_type> n2_neighbors = node_neighbors_[n2_uid];
      vector<size_type>::iterator it2 = find(n2_neighbors.begin(), n2_neighbors.end(), n1_uid);
      n2_neighbors.erase(it2);
      node_neighbors_[n2_uid] = n2_neighbors;
      
      // Remove edge from internal structures
      edges_.pop_back();
      uid2edge_index_.erase(node_pair);

      return 1;  // graph had edge and was successfuly removed
    } else {
      return 0;  // node did not exist in this graph
    }
  }


  /** Remove node pointed to by @a n_it from this graph.
   * @param[in] n_it     A node iterator that points to the node to remove.
   * @return A new node iterator that points to the beginning of our nodes_ vector;
   *
   * @pre @a n_it is a valid node iterator.
   *
   * @post All old iterators are invalidated.
   * @post For node n that the @a n_it points to: has_node(@a n) == false
   * @post If old has_node(@a n), new num_nodes() == old num_nodes() - 1.
   *                              new num_edges() == old num_edges() - n.degree().
   *                              For all nodes i that were adjacent to @a n,
   *                              new i.degree() == old i.degree() - 1 and
   *                              has_edge(@a n, i) == false.
   *       Else,                  new num_nodes() == old num_nodes().
   *                              new num_nodes() == old num_edges().
   *                              new degree() == old degree() for all nodes.
   *
   * Complexity: No more than O(num_nodes()) + complexity of remove_edge function
   *
   */
  node_iterator remove_node(node_iterator n_it) {
    Node n = *n_it;
    if (n_it != node_end() && node_begin() != node_end()) {
      remove_node(n);
    }
    node_iterator new_it = node_begin();
    return new_it;
  }


  /** Remove edge @a e from this graph.
   * @param[in] e      An edge to remove.
   * @return 1 if graph had node and was successfuly removed, else
   *         0 if node did not exist in this graph
   *
   * @pre @a e is a valid edge.
   *
   * @post For nodes n1 and n2, previously incident to edge @e, has_edge(@a n1, @a n2) == false
   * @post If old has_edge(@a n1, @a n2), new num_edges() == old num_edges() - 1.
   *                                      new n1.degree() == old n1.degree() - 1.
   *                                      new n2.degree() == old n2.degree() - 1.
   *       Else,                          new num_edges() == old num_edges().
   *                                      new n1.degree() == old n1.degree().
   *                                      new n2.degree() == old n2.degree().
   *
   * Complexity: No more than O(num_nodes() + num_edges())
   */
  size_type remove_edge(const Edge& e) {
    Node n1 = e.node1();
    Node n2 = e.node2();
    if (has_edge(n1, n2)) {
      remove_edge(n1, n2);
      return 1;  // graph had edge and was successfuly removed
    } else {
      return 0;  // node did not exist in this graph
    }
  }


  /** Remove edge pointed to by @a e_it from this graph.
   * @param[in] e_it   Edge iterator that points to the edge to remove.
   * @return A new edge iterator that points to the beginning of our edges_ vector.
   *
   * @pre @a e_it is a valid edge iterator.
   *
   * @post All old iterators are invalidated.
   * @post For nodes n1 and n2 that are incident to the edge @a e_it points to:
            has_edge(@a n1, @a n2) == false
   * @post If old has_edge(@a n1, @a n2), new num_edges() == old num_edges() - 1.
   *                                      new n1.degree() == old n1.degree() - 1.
   *                                      new n2.degree() == old n2.degree() - 1.
   *       Else,                          new num_edges() == old num_edges().
   *                                      new n1.degree() == old n1.degree().
   *                                      new n2.degree() == old n2.degree().
   *
   * Complexity: No more than O(num_nodes() + num_edges())
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    Edge e = *e_it;
    Node n1 = e.node1();
    Node n2 = e.node2();
    if (e_it != edge_end() && edge_begin() != edge_end()) {
      remove_edge(e);
    }
    edge_iterator new_it = edge_begin();
    return new_it;
  }
  
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
      uid_ = size_type(-1);
    }
    
    /** Set this node's position. */
    Point& position() {
      return fetch().position;
    }
    
    /** Return this node's position. */
    const Point& position() const {
      return fetch().position;
    }
    
    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->uid2index_[uid_];
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
      return graph_->node_neighbors_[uid_].size();
    }
    
    /** Return an incident iterator that points to the first edge in our
     * neighbors map for this node.
     *
     * Complexity: O(1).
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_,
                              uid_,
                              graph_->node_neighbors_[uid_].begin());
    }
    
    /** Return an incident iterator that points to one past the last edge
     * in our neighbors map for this node.
     *
     * Complexity: O(1).
     */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_,
                              uid_,
                              graph_->node_neighbors_[uid_].end());
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
     * Here, if nodes are in different graphs, ordering is defined
     * by the graph addresses. If nodes are in the same graph,
     * ordering is defined by the node indicies.
     */
    bool operator<(const Node& n) const {
      if (graph_ != n.graph_) {
        return (&graph_ < &n.graph_);
      } else {
        return (index() < n.index());
      }
    }
    
  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    
    Graph* graph_;
    size_type uid_;
    
    /** Private Constructor for Node */
    Node(const Graph* graph, size_type uid)
    : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
    
    /** Helper method to return the appropriate node. */
    internal_node& fetch() const {
      return graph_->nodes_[uid_];
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
    new_node.position = position;
    new_node.value = value;
    new_node.uid = next_uid_;
    nodes_[next_uid_] = new_node;
    
    uid2index_[next_uid_] = new_index;
    index2uid_.push_back(next_uid_);

    node_neighbors_[next_uid_] = vector<size_type>();
    next_uid_++;

    return Node(this, next_uid_-1);
  }
  
  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.index() < size() && n.graph_ == this);
  }
  
  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < size());
    size_type uid = index2uid_[i];
    return Node(this, uid);
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
      node1_uid_ = size_type(-1);
      node2_uid_ = size_type(-1);
      edge_ind_ = size_type(-1);
    }
    
    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, node1_uid_);
    }
    
    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, node2_uid_);
    }
    
    /** Set this edge's value. */
    E& value() {
      return graph_->edges_[edge_ind_].value;
    }
    
    /** Return this edge's value. */
    const E& value() const {
      return graph_->edges_[edge_ind_].value;
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
     * Here, if edges are in different graphs, ordering is defined
     * by the graph addresses. If edges are in the same graph,
     * ordering is defined by the edge indicies.
     */
    bool operator<(const Edge& e) const {
      if (graph_ != e.graph_) {
        return (&graph_ < &e.graph_);
      } else {
        return ((edge_ind_ < e.edge_ind_));
      }
    }
    
  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    
    Graph* graph_;
    size_type node1_uid_;
    size_type node2_uid_;
    size_type edge_ind_;
    
    /** Private Constructor for Edge */
    Edge(const Graph* graph, size_type node1_uid, size_type node2_uid,
         size_type edge_ind) : graph_(const_cast<Graph*>(graph)),
    node1_uid_(node1_uid), node2_uid_(node2_uid), edge_ind_(edge_ind) {
    }
  };
  
  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_.size();
  }
  
  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());
    return Edge(this, edges_[i].node1_uid, edges_[i].node2_uid, i);
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
    size_type nodea_uid = index2uid_[nodea_index];
    size_type nodeb_uid = index2uid_[nodeb_index];
    // check if this node b is a key in node a's neighbor vector
    vector<size_type> nodea_neighbors = node_neighbors_.at(nodea_uid);
    if (find(nodea_neighbors.begin(), nodea_neighbors.end(),
             nodeb_uid) != nodea_neighbors.end()) {
      return true;  //key found
    } else {
      return false;  //key not found
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
    
    size_type nodea_uid = index2uid_[nodea_index];
    size_type nodeb_uid = index2uid_[nodeb_index];
    
    if (has_edge(a, b) == true) {
      set<size_type> node_pair = {};
      node_pair.insert(nodea_uid);
      node_pair.insert(nodeb_uid);
      size_type edge_index = uid2edge_index_[node_pair];
      return Edge(this, nodea_uid, nodeb_uid, edge_index);
    } else {
      size_type new_index = num_edges();
      internal_edge new_edge = {};
      new_edge.node1_uid = nodea_uid;
      new_edge.node2_uid = nodeb_uid;
      new_edge.edge_index = new_index;
      edges_.push_back(new_edge);
      
      set<size_type> node_pair = {};
      node_pair.insert(nodea_uid);
      node_pair.insert(nodeb_uid);
      uid2edge_index_[node_pair] = new_index;
      
      // add node b to neighbor list of a and vise vera
      node_neighbors_[nodea_uid].push_back(nodeb_uid);
      node_neighbors_[nodeb_uid].push_back(nodea_uid);

      return Edge(this, nodea_uid, nodeb_uid, new_index);
    }
  }
  
  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    next_uid_ = 0;
    nodes_.clear();
    uid2index_.clear();
    index2uid_ = {};
    edges_ = {};
    uid2edge_index_.clear();
    node_neighbors_.clear();
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
      node_it_ = index2uid_.begin();
    }
    
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Return the node that this node iterator points to */
    Node operator*() const {
      return Node(graph_, *node_it_);
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
    typename vector<size_type>::const_iterator node_it_;
    
    /** Private Constructor for NodeIterator */
    NodeIterator(const Graph* graph,
                 typename vector<size_type>::const_iterator node_it) :
    graph_(const_cast<Graph*>(graph)), node_it_(node_it) {
    }
  };
  
  // Supply definitions AND SPECIFICATIONS for:
  
  /** Return a node iterator that points to the first node in our node vector.
   *
   * Complexity: O(1).
   */
  node_iterator node_begin() const {
    return NodeIterator(this, index2uid_.begin());
  }
  
  /** Return a node iterator that points to one past the last node in our node vector.
   *
   * Complexity: O(1).
   */
  node_iterator node_end() const {
    return NodeIterator(this, index2uid_.end());
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
      node1_uid_ = size_type(-1);
    }
    
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Return the edge that this incident iterator points to */
    Edge operator*() const {
      size_type node2_uid = *inc_it_;
      set<size_type> node_pair = {};
      node_pair.insert(node1_uid_);
      node_pair.insert(node2_uid);
      size_type edge_index = graph_->uid2edge_index_[node_pair];
      return Edge(graph_, node1_uid_, node2_uid, edge_index);
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
    size_type node1_uid_;
    typename vector<size_type>::const_iterator inc_it_;
    
    /** Private Constructor for IncidentIterator */
    IncidentIterator(const Graph* graph, size_type node1_uid,
                    typename vector<size_type>::const_iterator inc_it):
    graph_(const_cast<Graph*>(graph)), node1_uid_(node1_uid), inc_it_(inc_it) {
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
      edge_it_ = edges_.begin();
    }
    
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Return the edge that this edge iterator points to */
    Edge operator*() const {
      return Edge(graph_, (*edge_it_).node1_uid, (*edge_it_).node2_uid,
                  (*edge_it_).edge_index);
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
    typename vector<internal_edge>::const_iterator edge_it_;
    
    /** Private Constructor for EdgeIterator */
    EdgeIterator(const Graph* graph,
                 typename vector<internal_edge>::const_iterator edge_it) :
    graph_(const_cast<Graph*>(graph)), edge_it_(edge_it) {
    }
  };
  
  // Supply definitions AND SPECIFICATIONS for:
  
  /** Return an edge iterator that points to the first edge in our edge vector.
   *
   * Complexity: O(1).
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, edges_.begin());
  }
  
  /** Return an edge iterator that points to one past the last edge in our edge vector.
   *
   * Complexity: O(1).
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, edges_.end());
  }
  
private:
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  
  // internal type for graph nodes
  struct internal_node {
    size_type uid;
    Point position;
    V value;
  };
  
  
  // internal type for graph edges
  struct internal_edge {
    size_type node1_uid;
    size_type node2_uid;
    size_type edge_index;
    E value;
  };
  
  size_type next_uid_ = 0;  // counter for uid
  map<size_type, internal_node> nodes_;
  
  
  map<size_type, size_type> uid2index_;   // maps uid to index: key = uid, val = index
  vector<size_type> index2uid_;   // index to uid: position at index i is uid for node i
  
  vector<internal_edge> edges_;
  // maps uid to edge index
  map<set<size_type>, size_type> uid2edge_index_; // key = pair of uids of edge, val = edge index
  
  // key = node uid, val = vector of adj nodes to that node
  map<size_type, vector<size_type>> node_neighbors_;
};

#endif // CME212_GRAPH_HPP
