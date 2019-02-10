#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <set>
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
class Graph {
 private:

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

  /** Type of value contained in nodes */
  using node_value_type = V;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() 
    : nodes_(), edges_(), num_nodes_(0), num_edges_(0) {
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
      // Do Nothing
    }

    /** Return this node's position. */
    const Point& position() const {
      return fetch().position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(graph_ != nullptr);
      return index_;
    }

    /** 
     * @brief Return this node's value.
     * 
     * @return A node_value_type value. 
     *
     * @pre This node is valid.
     * @post Result is the value associated with this node. Can be used
     *   for assignment.
     *
     * Complexity: O(1).
     */
    node_value_type& value() { return fetch().value_; }

       /** 
     * @brief Return this (const) node's value.
     * 
     * @return A node_value_type value. 
     *
     * @pre This node is valid.
     * @post Result is the value associated with this node.
     *
     * Complexity: O(1).
     */
    const node_value_type& value() const { return fetch().value_; }

    /**
     * @brief Get the degree of this node.
     *
     * @return A size_type value.
     *
     * @pre This node is valid.
     * @post Result is the degree of this node.
     *
     * Complexity: O(1).
     */
    size_type degree() const {
      return graph_->adj_list_[index_].size();
    }

    /**
     * @brief Get an iterator to the first edge incident to this node.
     *
     * @return An iterator.
     *
     * @pre This node is valid.
     * @post Result points to the first edge incident to this node.
     *
     * Complexity: O(1).
     */
    incident_iterator edge_begin() const {
      assert(graph_ != nullptr);
      return IncidentIterator(graph_, *this
        , graph_->adj_list_[index_].begin());
    }

    /**
     * @brief Get an iterator to one past the last edge incident to this node.
     *
     * @return An iterator.
     *
     * @pre This node is valid. 
     * @post Result is the end iterator for edges incident to this node.
     * 
     * Complexity: O(1).
     */
    incident_iterator edge_end() const {
      assert(graph_ != nullptr);
      return IncidentIterator(graph_, *this
        , graph_->adj_list_[index_].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return graph_ == n.graph_ && index_ == n.index_;
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
      assert(graph_ == n.graph_);
      return index_ < n.index_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    
    Graph* graph_;
    size_type index_;

    /** Private constructor */
    Node(const Graph* graph, size_type index)
      : graph_(const_cast<Graph*>(graph)), index_(index) {
    }
    
    /** Find the corresponding internal_node in the graph
     *  by looping over all nodes and matching by index
     */
    internal_node& fetch() const {
      assert(graph_ != nullptr);
      return graph_->nodes_[index_];
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return num_nodes_;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value    The new nodes value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   * @post result node has specified @a position and @a value.
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, 
                const node_value_type& value = node_value_type()) {
    size_type index = num_nodes_;

    internal_node new_node {position, value};

    nodes_.push_back(new_node);
    adj_list_.push_back(std::map<size_type, size_type>());
    num_nodes_++;
    return node(index);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return n.graph_ == this && n.index_ < size();
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

    /** Return a node of this Edge */
    Node node1() const {
      return swap_ ? fetch().node2_ : fetch().node1_;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return swap_ ? fetch().node1_ : fetch().node2_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (node1() == e.node1())
        return node2() == e.node2();
      if (node1() == e.node2())
        return node2() == e.node1();
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      assert(graph_ == e.graph_);
      return edgeid_ < e.edgeid_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    friend class IncidentIterator;

    Graph* graph_;
    size_type edgeid_;
    bool swap_;
    
    /** Private constructor */
    Edge(const Graph* graph, size_type edgeid) 
      : graph_(const_cast<Graph*>(graph)), edgeid_(edgeid), swap_(false) {
    }

    /** Find the corresponding edge in the graph by edgeid */
    internal_edge& fetch() const {
      assert(graph_ != nullptr);
      return graph_->edges_[edgeid_];
    }

    /**
     * @brief Set whether the output of node1() and node2() are swapped.
     *
     * @param[in] b The new desired value of swap_
     * @return  None
     *
     * @pre This edge is valid.
     * @post If @a b is false, then node1() and node2() will return the
     *   endpoints of the edge in the same order as when the edge was added
     *   to the graph. Otherwise, they will be swapped.
     */
    void set_swap(bool b) { swap_ = b; }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return num_edges_;
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
    return get_edgeid(a, b) < num_edges_;
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
    size_type edgeid = get_edgeid(a,b);
    if (edgeid == num_edges_) {
      internal_edge new_edge {a, b};

      edges_.push_back(new_edge);
      adj_list_[a.index_][b.index_] = edgeid;
      adj_list_[b.index_][a.index_] = edgeid;
      num_edges_++;
    }
    return edge(edgeid);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    adj_list_.clear();
    num_nodes_ = 0;
    num_edges_ = 0;
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
     * @brief Dereference operator.
     *
     * @return A node.
     *
     * @pre This iterator is valid, and is not the one-past-the-end iterator.
     * @post Result is the node pointed to by this iterator.
     *
     * Complexity: O(1).
     */
    Node operator*() const { 
      assert(idx_ < graph_->size());
      return graph_->node(idx_); 
    }

    /**
     * @brief Increment operator.
     *
     * @return This NodeIterator.
     *
     * @pre This iterator is valid, and is not the one-past-the-end iterator.
     * @post If this iterator initially pointed to the last node, then it
     *   becomes the one-past-the-end iterator. Otherwise, it points to the
     *   next node, ordered by index.
     *
     * Complexity: O(1).
     */
    NodeIterator& operator++() { 
      assert(idx_ < graph_->size());
      idx_++;
      return *this;
    }

    /**
     * @brief Equality comparison operator.
     *
     * @param[in] iter Another iterator to compare with.
     * @return True or false.
     *
     * @pre Both iterators are valid.
     * @post Result is true iff this and @a iter point to the same node,
     *   or both are the end iterator.
     * 
     * Complexity: O(1).
     */
    bool operator==(const NodeIterator& iter) const {
      return graph_ == iter.graph_ && idx_ == iter.idx_;
    }

   private:
    friend class Graph;
    const Graph* graph_;
    size_type idx_;

    /**
     * @brief Private constructor for NodeIterator.
     *
     * @param[in] g   The graph containing the nodes.
     * @param[in] idx The index of the node to point to.
     * @return A NodeIterator.
     *
     * @pre 0 <= @a idx < @a g.num_nodes()
     * @post This iterator points to the specified node.
     */
    NodeIterator(const Graph* g, size_type idx) 
      : graph_(const_cast<Graph*>(g)), idx_(idx) {}
  };

  /**
   * @brief Get a NodeIterator to the first node.
   *
   * @return A NodeIterator.
   *
   * @post Result points to the first node of the graph.
   *
   * Complexity: O(1).
   */
  node_iterator node_begin() const { return NodeIterator(this, 0); }

  /**
   * @brief Get the one-past-the-end NodeIterator.
   *
   * @return A NodeIterator.
   *
   * @post Result is the end iterator for nodes.
   * 
   * Complexity: O(1).
   */
  node_iterator node_end() const { return NodeIterator(this, num_nodes_); }

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
     * @brief Dereference operator.
     *
     * @return An edge.
     *
     * @pre This iterator is valid, and is not the one-past-the-end iterator.
     * @post Result is the edge pointed to by this.
     * 
     * Complexity: O(1).
     */
    Edge operator*() const {
      assert(iter_ != node1_.edge_end().iter_);
      Edge e = graph_->edge(iter_->second);
      if (e.node1() != node1_) e.set_swap(true);
      return e;
    }

    /**
     * @brief Increment operator.
     *
     * @return This iterator.
     *
     * @pre This iterator is valid, and is not the one-past-the-end iterator.
     * @post If this iterator pointed to the last incident edge, it now points
     *   to one-past-the-end. Otherwise, it now points to the next incident
     *   edge, ordered by the other endpoint's index.
     *
     * Complexity: O(1) amortized.
     */
    IncidentIterator& operator++() {
      assert(iter_ != node1_.edge_end().iter_);
      iter_++;
      return *this;
    }

    /**
     * @brief Equality comparison operator.
     *
     * @param[in] it Another IncidentIterator to compare to.
     * @return True or false.
     *
     * @pre Both iterators are valid.
     *
     * @post Result is true iff they are iterators for edges
     *   incident to the same node, and if those edges are the same.
     */
    bool operator==(const IncidentIterator& it) const {
      return graph_ == it.graph_ && iter_ == it.iter_;
    }

   private:
    friend class Graph;
    /** Type of adjacency list iterator */
    using internal_iter_type = std::map<size_type,size_type>::iterator;

    const Graph* graph_;
    const Node node1_;
    internal_iter_type iter_;

    /**
     * @brief Private constructor.
     *
     * @param[in] g  The graph we're looking at.
     * @param[in] n  The node whose incident edges we want to iterate through.
     * @param[in] it Refers to where the edge is stored internally.
     * @return An IncidentIterator.
     *
     * @pre @a n is a node in @a g, and @a it is a valid iterator.
     * @post This is an IncidentIterator that points to @a *it, and 
     *   iterates through the edges incident to @a n.
     */
    IncidentIterator(const Graph* g, const Node& n, internal_iter_type it)
      : graph_(const_cast<Graph*>(g)), node1_(n), iter_(it) {}
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
     * @brief Dereference operator.
     *
     * @return An edge.
     *
     * @pre This iterator is valid, and is not the one-past-the-end iterator.
     * @post Result is the edge that this iterator points to.
     *
     * Complexity: O(1).
     */
    Edge operator*() const { 
      assert(edgeid_ < graph_->num_edges());
      return graph_->edge(edgeid_); 
    }

    /**
     * @brief Increment operator.
     *
     * @return This iterator.
     *
     * @pre This iterator is valid, and is not the one-past-the-end iterator.
     * @post If this iterator used to point to the last edge, it now points
     *   to one-past-the-end. Otherwise, it now points to the next edge,
     *   ordered by edgeid.
     *
     * Complexity: O(1).
     */
    EdgeIterator& operator++() {
      assert(edgeid_ < graph_->num_edges());
      edgeid_++;
      return *this;
    }

    /**
     * @brief Equality comparison operator.
     *
     * @param[in] iter Another EdgeIterator.
     * @return True or false.
     *
     * @pre Both iterators are valid.
     * @post Result is true iff *this and @a iter point to the same edge,
     *   or both are the end iterator.
     */
    bool operator==(const EdgeIterator& iter) const {
      return graph_ == iter.graph_ && edgeid_ == iter.edgeid_;
    }

   private:
    friend class Graph;

    const Graph* graph_;
    size_type edgeid_;

    /**
     * @brief Private constructor.
     *
     * @param[in] g      Pointer to a graph.
     * @param[in] edgeid The id of the edge we want to point at.
     * @return An edge iterator.
     *
     * @pre 0 <= @a edgeid < @a g.num_edges()
     * @post Result points to the edge in @a *g with given @a edgeid.
     */
    EdgeIterator(const Graph* g, size_type edgeid)
      : graph_(const_cast<Graph*>(g)), edgeid_(edgeid) {}
  };

  /**
   * @brief Get an iterator to the first edge of the graph.
   *
   * @return An iterator for edges in the graph.
   *
   * @post Result points to the first edge of the graph.
   *
   * Complexity: O(1).
   */
  edge_iterator edge_begin() const { return EdgeIterator(this, 0); }

  /**
   * @brief Get the one-past-the-end EdgeIterator.
   *
   * @return An iterator for edges in the graph.
   *
   * @post Result is the one-past-the-end iterator for edges.
   *
   * Complexity: O(1).
   */
  edge_iterator edge_end()   const { return EdgeIterator(this, num_edges_); }


 private:

  struct internal_node {
    Point position_;
    node_value_type value_;
  };

  struct internal_edge {
    Node node1_;
    Node node2_;
  };

  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;
  std::vector<std::map<size_type, size_type>> adj_list_;
  size_type num_nodes_;
  size_type num_edges_;

  /** Looks for an edge in the graph and returns its index
   * @param[in] a one node of the edge to find
   * @param[in] b the other node of the edge find
   * @return An index into the array of edges.
   * 
   * @pre @a a and @b b are valid nodes in the graph
   *
   * @post If (a,b) is an edge in the graph, result is it's id. Otherwise,
   *   return is the number of edges in the graph.
   */
  size_type get_edgeid(const Node& a, const Node& b) const{
    auto search = adj_list_[a.index_].find(b.index_);
    if (search == adj_list_[a.index_].end())
      return num_edges_;
    else
      return search->second;
  }
};

#endif // CME212_GRAPH_HPP
