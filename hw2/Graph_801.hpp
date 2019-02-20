#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

// added by me
#include <map>


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = int, typename E = int>
class Graph {
 public:
  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Type of the nodes' value. */
  using node_value_type = V;

  /** Type of the edges' value. */
  using edge_value_type = E;

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
  /** Type of unique identifiers. Same as the indexes type. */
  using uid_type  = unsigned;

 /** Forward declaration of two structures containing node and edge data. */
 private:
  // private struct containing node info
  struct node_info {
    size_type idx;
    V value;
    Point position;
  };

  // private struct containing edge info
  struct edge_info {
    size_type idx;
    E value;
    uid_type node1_uid;
    uid_type node2_uid;
  };

 public:
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {}

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
  class Node : private totally_ordered<Node>{
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
    Node() {} // this invalid constructor does not do anything

    /** Non const version of position(). */
    Point& position() {
      assert(is_valid());
      return graph_->nodes_[uid_].position;
    }

    /** Return this node's position. */
    const Point& position() const {
      assert(is_valid());
      return graph_->nodes_[uid_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(is_valid());
      return graph_->nodes_[uid_].idx;
    }

    /** Return the uinique identifier of the node. */
    uid_type uid() const {
      assert(is_valid());
      return uid_;
    }

    /**
     * @brief Return the value associated with the node
     * @return The value of the node
     */
    node_value_type& value() {
      assert(is_valid());
      return graph_->nodes_[uid_].value;
    }

    /**
     * @brief Return the value associated with the node
     * @return The value of the node
     * const version
     */
    const node_value_type& value() const {
      assert(is_valid());
      return graph_->nodes_[uid_].value;
    }

    /**
     * @brief Return the degree of the node
     * @return The degree of the node
     *
     * Look at the adjacency list and return the number of edges
     * for which one of the end is the node.
     */
    size_type degree() const {
      std::map<uid_type, std::vector<uid_type>>& adj_list = graph_->adj_list_;
      if (adj_list.find(uid_) != adj_list.end()) return adj_list.at(uid_).size();
      else return 0;
    }

    /**
     * @brief Valid IncidentIterator constructor
     * @return A valid IncidentIterator with iedge_ member variable initialized
     *         to zero.
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, uid_, 0);
    }

    /**
     * @brief Valid IncidentIterator constructor
     * @return A valid IncidentIterator with iedge_ member variable initialized
     *         to the degree of the node.
     */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, uid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return ((n.uid_ == uid_) && (n.graph_ == graph_));
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
      return (uid_ < n.uid_);
    }

   private:
    friend class Graph;
    Graph* graph_;
    uid_type uid_;

    // private constructor accessible in Graph class
    Node(const Graph* graph, const uid_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {}

    // check validity of node
    bool is_valid() const {
      uid_type max_uid = graph_->nodes_.size();
      size_type max_idx = graph_->ni2u_.size();
      return uid_ >=0 && uid_ < max_uid &&
             graph_->nodes_[uid_].idx >= 0 &&
             graph_->nodes_[uid_].idx < max_idx &&
             graph_->ni2u_[graph_->nodes_[uid_].idx] == uid_;
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return ni2u_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return ni2u_.size();
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return n.is_valid();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < size());
    return Node(this, ni2u_[i]);
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = {}) {
    size_type idx = size();
    uid_type  uid = nodes_.size();
    Node new_node = Node(this, uid);
    struct node_info ni = {idx, value, position};
    nodes_.push_back(ni);
    ni2u_.push_back(uid);
    return new_node;
  }

  /**
   * @brief Remove a node
   *
   * @param n A node
   * @return The number of nodes after remove
   *
   * @pre  @ n must be a valid node
   * @post Remove all the edges incident to @ n. Refer to the
   *       remove_edge(const Node&, const Node&) for details on this.
   *       Once those edges are renoved, @ n itself is removed from the
   *       graph. All Edge objects, Node object, Edge iteraators, Node
   *       iterators, Incident iterators refering to the deleted edges
   *       or node are invalidated.
   *
   * The complextity is O(d*n_e + n_n) where d is the degree of the node,
   * n_e is the number of edges in the graph and n_n is the number of nodes.
   * This is for the worst case scenario.
   */
  size_type remove_node(const Node& n) {
    assert(has_node(n));
    // remove all incident edges
    for (auto it=n.edge_begin(); it!=n.edge_end();) {
      remove_edge(*it);
    }
    // invalidate the node
    size_type idx = nodes_[n.uid_].idx;
    auto it = ni2u_.erase(ni2u_.begin() + idx);
    // update all subsequent nodes
    for (; it!=ni2u_.end(); ++it) nodes_[*it].idx--;
    return size();
  }

  /**
   * @brief Remove a node
   *
   * Refer the remove_node(const Node&) for details.
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return n_it;
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
    Edge() {}

    /** Return a node of this Edge */
    Node node1() const {
      assert(is_valid());
      return Node(graph_, node1_uid_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      assert(is_valid());
      return Node(graph_, node2_uid_);
    }

    /**
     * @brief Return the value associated with the edge
     * @return The value of the edge
     */
    edge_value_type& value() {
      assert(is_valid());
      return graph_->edges_[uid_].value;
    }

    /**
     * @brief Return the value associated with the edge
     * @return The value of the edge
     * const version
     */
    const edge_value_type& value() const {
      assert(is_valid());
      return graph_->edges_[uid_].value;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (e.uid_ == uid_) && (e.graph_ == graph_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_ != e.graph_) return uid_ <= e.uid_;
      return (uid_ < e.uid_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph_;
    uid_type uid_;
    uid_type node1_uid_;
    uid_type node2_uid_;

    Edge(const Graph* graph, const  uid_type uid,
         const uid_type node1_uid, const uid_type node2_uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid),
          node1_uid_(node1_uid), node2_uid_(node2_uid) {}

    // check validity of edge
    bool is_valid() const {
      uid_type max_uid = graph_->edges_.size();
      size_type max_idx = graph_->ei2u_.size();
      return uid_ >=0 && uid_ < max_uid &&
             graph_->edges_[uid_].idx >= 0 &&
             graph_->edges_[uid_].idx < max_idx &&
             graph_->ei2u_[graph_->edges_[uid_].idx] == uid_;
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return ei2u_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());
    uid_type uid = ei2u_[i];
    return Edge(this, uid, edges_[uid].node1_uid, edges_[uid].node2_uid);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert(has_node(a) && has_node(b));
    if (a == b) return false;
    uid_type nuid_min = std::min(a.uid_, b.uid_);
    uid_type nuid_max = std::max(a.uid_, b.uid_);
    // check if an edge might exist in the first place
    if (adj_list_.find(nuid_min) != adj_list_.end()) {
      for (uid_type i : adj_list_.at(nuid_min))
        if (edges_[i].node2_uid == nuid_max) return true;
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
  Edge add_edge(const Node& a, const Node& b,
      const edge_value_type value = {}) {
    assert(has_node(a) && has_node(b) && !(a == b));
    uid_type nuid_min = std::min(a.uid_, b.uid_);
    uid_type nuid_max = std::max(a.uid_, b.uid_);
    uid_type uid      = find_edge(a, b);

    if (uid == edges_.size()){
      // the edge does not exist, create it
      size_type idx = num_edges();
      struct edge_info ei = {idx, value, nuid_min, nuid_max};
      edges_.push_back(ei);
      ei2u_.push_back(uid);
      adj_list_[nuid_min].push_back(uid);
      adj_list_[nuid_max].push_back(uid);
      return Edge(this, uid, nuid_min, nuid_max);
    }
    // the edge exists, return it
    else return Edge(this, uid, nuid_min, nuid_max);
  }

  /**
   * @brief Remove an edge from the graph
   *
   * @param a A node
   * @param b A node
   * @return The number of edges in the graph
   *
   * @post If @ a or @ b is not valid, or if they are equal,
   *       or if there isn't an edge between them, the graph is
   *       left unchanged. Otherwise, it is removed and all Edge
   *       objects representing this edge or Edge iterators refering
   *       to this edge are invalidated.
   *
   * The complexity of the function is O(num_edges).
   */
  size_type remove_edge(const Node& a, const Node& b) {
    // either the arguments are not valid, or they are identical; return immediately
    if (!has_node(a) || !has_node(b) || a==b) return num_edges();

    uid_type nuid_min = std::min(a.uid_, b.uid_);
    uid_type nuid_max = std::max(a.uid_, b.uid_);
    if (adj_list_.find(nuid_min) != adj_list_.end()){
      std::vector<uid_type>& neighbors1 = adj_list_.at(nuid_min);

      for (auto it=neighbors1.begin(); it!=neighbors1.end(); ++it)
        if (edges_[*it].node2_uid == nuid_max) {
          // the edge was found, fetch its index
          size_type idx = edges_[*it].idx;
          // invalidate the edge
          auto iter = ei2u_.erase(ei2u_.begin() + idx);
          // update the index of the subsequent edges
          for (; iter!=ei2u_.end(); ++iter) edges_[*iter].idx--;
          // update the adjacency list; 2 entries to delete
          neighbors1.erase(it);
          // delete symmetric entry
          std::vector<uid_type>& neighbors2 = adj_list_.at(nuid_max);
          for (auto it2=neighbors2.begin(); it2!=neighbors2.end(); ++it2)
            if (edges_[*it2].node1_uid == nuid_min) {
              neighbors2.erase(it2);
              break;
            }
          break;
        }
    }
    return num_edges();
  }

  /**
   * @brief Remove an edge of the graph
   *
   * Refer to the function remove_edge(const Node&, const Node&)
   * for details.
   */
  size_type remove_edge(const Edge& e) {
    Node n1 = Node(this, e.node1_uid_);
    Node n2 = Node(this, e.node2_uid_);
    return remove_edge(n1, n2);
  }

  /**
   * @brief Remove an edge of the graph
   *
   * Refer to the function remove_edge(const Node&, const Node&)
   * for details.
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    ni2u_.clear();
    edges_.clear();
    ei2u_.clear();
    adj_list_.clear();
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
    NodeIterator() {}

    /**
     * @brief Dereferencing operator for the node iterator.
     * @return The node at which the iterator points to
     */
    Node operator*() const {
      return Node(graph_, graph_->ni2u_[idx_]);
    }

    /**
     * @brief Increment operator for the node iterator.
     * @return Itself after incrementing
     *
     * This function increments the node_id_ member variable.
     */
    node_iterator& operator++() {
      idx_++;
      return *this;
    }

    /**
     * @brief Equality comparator
     * @return True if both iterators point to the same node in the same graph
     *         False otherwise
     */
    bool operator==(const node_iterator& ni) const {
      return (graph_ == ni.graph_) && (idx_ == ni.idx_);
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type idx_;

    NodeIterator(const Graph* graph, const size_type idx)
        : graph_(const_cast<Graph*>(graph)), idx_(idx) {}
  };

  /**
   * @brief Valid NodeIterator constructor
   * @return A NodeIterator object with node_id_ initialized to 0
   */
  node_iterator node_begin() const {
    return node_iterator(this, 0);
  }

  /**
   * @brief Valid NodeIterator constructor
   * @return A NodeIterator object with node_id_ initialized to the size of the graph
   */
  node_iterator node_end() const {
    return node_iterator(this, size());
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
    IncidentIterator() {}

    /**
     * @brief Dereferencing operator for IncidentIterator
     * @return The edge it points to
     *
     * The function makes sure that node1() of the returned edge always
     * is the node corresponding to the node_id_ of the object.
     */
    Edge operator*() const{
      uid_type edge_uid  = graph_->adj_list_.at(node_uid_)[iedge_];
      uid_type node1_uid = graph_->edges_[edge_uid].node1_uid;
      uid_type node2_uid = graph_->edges_[edge_uid].node2_uid;
      if (node_uid_ == node1_uid) return Edge(graph_, edge_uid, node1_uid, node2_uid);
      else return Edge(graph_, edge_uid, node2_uid, node1_uid);
    }

    /**
     * @brief Increment operator
     * @return An incremented IncidentIterator
     */
    incident_iterator& operator++() {
      iedge_++;
      return *this;
    }

    /**
     * @brief Equality comparator
     * @return True if the object points to the same graph, the same node and the
     *         same edge incident to the node
     *         False otherwise
     */
    bool operator==(const incident_iterator& iit) const {
      return (graph_ == iit.graph_) && (node_uid_ == iit.node_uid_)
                                    && (iedge_ == iit.iedge_);
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type node_uid_;
    size_type iedge_;

    IncidentIterator(const Graph* graph, const size_type node_uid,
        const size_type iedge)
        : graph_(const_cast<Graph*>(graph)), node_uid_(node_uid), iedge_(iedge) {}
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
    EdgeIterator() {}

    /**
     * @brief Dereferencing operator for the edge iterator.
     * @return The edge at which the iterator points to
     */
    Edge operator*() const {
      uid_type uid       = graph_->ei2u_[idx_];
      uid_type node1_uid = graph_->edges_[uid].node1_uid;
      uid_type node2_uid = graph_->edges_[uid].node2_uid;
      return Edge(graph_, uid, node1_uid, node2_uid);
    }

    /**
     * @brief Increment operator for the node iterator.
     * @details Itself after incrementing
     *
     * This function increments the edge_id_ member variable.
     */
    EdgeIterator& operator++() {
      idx_++;
      return *this;
    }

    /**
     * @brief Equality comparator
     * @return True if both iterators point to the sane edge in the same graph
     *         False otherwise
     */
    bool operator==(const EdgeIterator& ei) const {
      return (graph_ == ei.graph_) && (idx_ == ei.idx_);
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type idx_;

    EdgeIterator(const Graph* graph, const size_type idx)
        : graph_(const_cast<Graph*>(graph)), idx_(idx) {}
  };

  /**
   * @brief Valid EdgeIterator constructor
   * @return An EdgeIterator object with edge_id_ initialized to 0
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /**
   * @brief Valid EdgeIterator constructor
   * @return An EdgeIterator object with edge_id_ initialized to the number of edges
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

  /** A function used for debugging. Check the validity of the entire graph. */
  /*
  bool check_val(){
    for (auto it = node_begin(); it!=node_end(); ++it)
      if (!((*it).is_valid())) return false;
    for (auto it = edge_begin(); it!=edge_end(); ++it)
      if (!((*it).is_valid())) return false;
    return true;
  }
  */

 private:

  /**
   * @ nodes_ vector of all nodes created so far
   * @ ni2u_  vector mapping indexes of current nodes to node uid
   * @ edges_ vector of all edges created so far
   * @ ei2u_  vector mapping indexes of current edges to edge uid
   * @ adj_list_ map from node uid to vector of edge uid of
   *             edges adjacent to the node
   */
  std::vector<struct node_info> nodes_;
  std::vector<uid_type>  ni2u_;
  std::vector<struct edge_info> edges_;
  std::vector<uid_type>  ei2u_;
  std::map<uid_type, std::vector<uid_type>> adj_list_;

  /** Find an edge and return its id
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return either the id of the edge if it exists,
   *         otherwise the number of edges in the graph.
   */
  uid_type find_edge(const Node& a, const Node& b){
    assert(has_node(a) && has_node(b) && !(a == b));
    uid_type nuid_min = std::min(a.uid_, b.uid_);
    uid_type nuid_max = std::max(a.uid_, b.uid_);
    if (adj_list_.find(nuid_min) != adj_list_.end())
      for (uid_type i : adj_list_.at(nuid_min))
        if (edges_[i].node2_uid == nuid_max) return i;
    return edges_.size();
  }
};

#endif // CME212_GRAPH_HPP
