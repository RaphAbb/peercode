#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <utility>

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

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Type of Node value **/
  using node_value_type = V;

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

  typedef std::pair<size_type, size_type> pair;

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
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->points[index_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return index_;
    }

    /**
     * @brief returns the value of a node in the graph
     *
     * @return node_value_type of the node
     *
     * @post value = values[index_]
     *
     * Runtime: O(1), also allows us to assign a new value to a node through Node object
     */
    node_value_type& value() {
      return graph_->values[index_];
    }

    /**
     * @brief returns reference to the value of a node in the graph
     *
     * @return node_value_type& of the node's value
     *
     * @post value = values[index_]
     *
     * Runtime: O(1)
     */
    const node_value_type& value() const {
      return graph_->values[index_];
    }

    /**
     * @brief calculates the degree of a node
     *
     * @return degree of node with type size_type
     *
     * @post deg = edges[deg].size()
     *
     * Runtime: O(deg n), uses edge iterator to calculate degree
     */
    size_type degree() const {
      size_type deg = 0;
      for (auto it = this->edge_begin(); it != this->edge_end(); ++it) {
        deg++;
      }
      return deg;
    }

    /**
     * @brief returns a node iterator to the beginning of the nodes (list)
     *
     * @return incident iterator object
     *
     * @post only returns valid IncidentIterator if index_ is a key of the edges map (has edges in the graph)
     *
     * Runtime O(1) since searching a hash map is average O(1) time
     */
    incident_iterator edge_begin() const {
      if (graph_->edges.find(index_) != graph_->edges.end()) {
        return IncidentIterator(index_, graph_->edges.at(index_).begin(), graph_);  
      }
      return IncidentIterator();
    }

    /**
     * @brief returns a node iterator to (one past) the end of the nodes (list)
     *
     * return incident iterator object
     *
     * @post returns valid IncidentIterator if index_ is a key of the edges map (has edges in the graph)
     *
     * Runtime O(1) since searching a hash map is average O(1) tune
     */
    incident_iterator edge_end() const {
      if (graph_->edges.find(index_) != graph_->edges.end()) {
        return IncidentIterator(index_, graph_->edges.at(index_).end(), graph_);
      }
      return IncidentIterator();
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if ((graph_ == n.graph_) and (index_ == n.index_)) {
        return true;
      }
      (void) n;          // Quiet compiler warning
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
      if (index_ < n.index_) {
        return true;
      }
      (void) n;           // Quiet compiler warning
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    size_type index_;
    graph_type* graph_;

    //Private constructor for use in the Graph class
    Node(size_type idx, const graph_type* g)
      : index_(idx), graph_(const_cast<graph_type*>(g)) { }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return points.size();
  }

  /** Synonym for size()make . */
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
  Node add_node(const Point& position, const node_value_type& val= node_value_type(-1)) {
    points.push_back(position);
    values.push_back(val);
    (void) position;      // Quiet compiler warning
    return Node(this->num_nodes() - 1, this);        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if ((n.graph_ == this) and (n.index_ < this->num_nodes())) {
      return true;
    }
    (void) n;            // Quiet compiler warning
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if (i < this->num_nodes()) {
      return Node(i, this);
    }
    (void) i;             // Quiet compiler warning
    return Node();        // Invalid node
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->node(node1_);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(node2_);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if ((node1_ == e.node1_ and node2_ == e.node2_) or (node1_ == e.node2_ and node2_ == e.node1_)) {
        return true; 
      }
      (void) e;           // Quiet compiler warning
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (index_ < e.index_) {
        return true;
      }
      (void) e;           // Quiet compiler warning
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    size_type node1_;
    size_type node2_;
    size_type index_;
    graph_type* graph_;

    //Private constructor for Edge
    Edge(size_type n1, size_type n2, size_type idx, const graph_type* g)
      : node1_(n1), node2_(n2), index_(idx), graph_(const_cast<graph_type*>(g))  { }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_idxs.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    if (i < this->num_edges()) {
      return Edge(edge_idxs[i].first, edge_idxs[i].second, i, this);
    }
    (void) i;             // Quiet compiler warning
    return Edge();        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // only need to check this case because in add_node we actually add a, b and b, a cases
    // in analogy with adjacency list
    if (edges.find(a.index_) != edges.end()) {
      if (edges.at(a.index_).find(b.index_) != edges.at(a.index_).end()) {
        return true;
      }
    }
    (void) a; (void) b;   // Quiet compiler warning
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
    // check if edge already exists 
    if (edges.find(a.index_) != edges.end()) {
      if (edges.at(a.index_).find(b.index_) != edges.at(a.index_).end()) {
        return Edge(a.index_, b.index_, edges.at(a.index_).at(b.index_), this);
      }
    }
    size_type curr_size = this->num_edges();
    //create map if it doesn't yet exist in edges (i.e., node a has no neighbors yet)
    if (edges.find(a.index_) == edges.end()) {
      std::unordered_map<size_type, size_type> a_nbhd;
      a_nbhd.insert({b.index_, curr_size});
      edges.insert({a.index_, a_nbhd});
    } else {
      edges.at(a.index_).insert({b.index_, curr_size});
    }
    //add a (b, a) edge (similar to adjacency list)
    if (edges.find(b.index_) == edges.end()) {
      std::unordered_map<size_type, size_type> b_nbhd;
      b_nbhd.insert({a.index_, curr_size});
      edges.insert({b.index_, b_nbhd});
    } else {
      edges.at(b.index_).insert({a.index_, curr_size});
    }
    std::pair<size_type, size_type> e (a.index_, b.index_);
    edge_idxs.push_back(e);
    (void) a, (void) b;   // Quiet compiler warning
    return Edge(a.index_, b.index_, curr_size, this);        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    points.clear();
    values.clear();
    edges.clear();
    edge_idxs.clear();
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
     * @brief dereferences node iterator to get the underlying node
     *
     * @post node referred to by NodeIterator
     *
     * Runtime: O(1)
     */
    Node operator*() const {
      return graph_->node(niter_);
    }

    /**
     * @brief increments node iterator
     *
     * @post reference to NodeIterator that's "next"
     *
     * Runtime O(1)
     */
    node_iterator& operator++() {
      niter_++;
      return *this;
    }

    /** 
     * @brief checks equality of two node iterators
     *
     * @pre any NodeIterator
     * @post True if two iterators are equal (i.e. have same graph and are at the same "position"), false otherwise
     *
     * Runtime O(1)
     */
    bool operator==(const node_iterator& node_iter) const {
      if (graph_ == node_iter.graph_ and niter_ == node_iter.niter_) {
        return true;
      }
      return false;
    }

   private:
    friend class Graph;

    size_type niter_;
    graph_type* graph_;

    NodeIterator(size_type idx, const graph_type* g)
      : niter_(idx), graph_(const_cast<graph_type*>(g)) { }
  };

  /**
   * @brief returns a NodeIterator the the "first" node (i.e., referring the beginning of the Node list)
   *
   * @post NodeIterator which dereferences to the first node in the graph (in the nodes list)
   *
   * Runtime O(1)
   */
  node_iterator node_begin() const {
    return NodeIterator(0, this);
  }

  /**
   * @brief returns a NodeIterator the the "last" node (i.e., referring to one past the end of the Node list)
   *
   * @post NodeIterator which refers to the non-existent node "one past the end" of our Nodes
   *
   * Runtime O(1)
   */
  node_iterator node_end() const {
    return NodeIterator(this->num_nodes(), this);
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

    /**
     * @brief dereferences IncidentIterator to get an edge in the graph incident to the original node
     *
     * @post returns a valid edge in our graph with the correct orientation, incident to _original__
     *
     * Runtime O(1)
     */
    Edge operator*() const {
      return Edge(original_, iiter_->first, iiter_->second, graph_);
    }

    /**
     * @brief increments the IncidentIterator to the "next" edge
     *
     * @post returns a valid IncidentIterator that points to the "next" edge incident to the original node
     *
     * Runtime O(1)
     */
    incident_iterator& operator++() {
      iiter_++;
      return *this;
    }

    /**
     * @brief checks if two IncidentIterators are equal
     *
     * @pre any IncidentIterator
     * @post true if this and _iit_ are equal (i.e., refer to the same graph and the same index_ "position" and have same incident node)
     *
     * Runtime O(1)
     */
    bool operator==(const incident_iterator& iit) const {
      if (graph_ == iit.graph_ and original_ == iit.original_ and iiter_ == iit.iiter_) {
        return true;
      }
      return false;
    }
    

   private:
    friend class Graph;
    
    size_type original_;
    std::unordered_map<size_type, size_type>::iterator iiter_;
    graph_type* graph_;

    IncidentIterator(size_type o, std::unordered_map<size_type, size_type>::iterator iit, const graph_type* g)
      : original_(o), iiter_(iit), graph_(const_cast<graph_type*>(g)) { }
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

    /**
     * @brief dereferences the EdgeIterator to get a valid edge in our graph
     *
     * @post valid Edge in our graph
     *
     * Runtime O(1)
     */
    Edge operator*() const {
      return graph_->edge(eiter_);
    }

    /**
     * @brief Increments the EdgeIterator to get the "next" valid edge in our graph
     *
     * @post EdgeIterator referring to the "next" edge in our graph
     *
     * Runtime O(1)
     */
    edge_iterator& operator++() {
      eiter_++;
      return *this;
    }

    /**
     * @brief Checks if two EdgeIterator are equal
     *
     * @pre any EdgeIterator
     * @post true iff the this and _eit_ are the same (i.e., refer to the same graph and same "position")
     *
     * Runtime O(1)
     */
    bool operator==(const edge_iterator& eit) const {
      if (graph_ == eit.graph_ and eiter_ == eit.eiter_) {
        return true;
      }
      return false;
    }

   private:
    friend class Graph;

    size_type eiter_;
    graph_type* graph_;

    EdgeIterator(size_type eit, const graph_type* g)
      : eiter_(eit), graph_(const_cast<graph_type*>(g)) { }
  };

  /**
   * @brief Returns an EdgeIterator the the "first" edge in our graph (that is, the first edge in our list of edges)
   *
   * @post EdgeIterator referring the first edge in edge_indxs
   *
   * Runtime O(1)
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(0, this);
  }

  /**
   * @brief Returns an EdgeIterator the the "last" edge in our graph
   *
   * @post EdgeIterator referring to "one past the end" of the list of valid edges
   *
   * Runtime O(1)
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this->num_edges(), this);
  }

 private:

  std::vector<Point> points;
  std::unordered_map<size_type, std::unordered_map<size_type, size_type>> edges;
  std::vector<std::pair<size_type, size_type>> edge_idxs;
  std::vector<node_value_type> values;

};

#endif // CME212_GRAPH_HPP
