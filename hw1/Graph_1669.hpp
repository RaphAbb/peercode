#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

// Added by me
#include <map>

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
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Type of the nodes' value. */
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

    /** Return this node's position. */
    const Point& position() const {
      // proxy design: let the graph object handles everything
      assert(id_ < graph_->size());
      return graph_->point_list_[id_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return this->id_;
    }

    /**
     * @brief Return the value associated with the node
     * @return The value of the node
     */
    node_value_type& value(){
        // proxy design
        assert(id_ < graph_->size());
        return graph_->node_value_list_[id_];
    }

    /**
     * @brief Return the value associated with the node
     * @return The value of the node
     * const version
     */
    const node_value_type& value() const{
        // proxy design; first check if the node exists
        assert(id_ < graph_->size());
        return graph_->node_value_list_[id_];
    }

    /**
     * @brief Return the degree of the node
     * @return The degree of the node
     *
     * Look at the adjacency list and return the number of edges
     * for which one of the end is the node.
     */
    size_type degree() const {
      std::map<size_type, std::vector<size_type>>& adj_list = graph_->adj_list_;
      if (adj_list.find(id_) != adj_list.end()) return adj_list.at(id_).size();
      else return 0;
    }

    /**
     * @brief Valid IncidentIterator constructor
     * @return A valid IncidentIterator with iedge_ member variable initialized
     *         to zero.
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, id_, 0);
    }

    /**
     * @brief Valid IncidentIterator constructor
     * @return A valid IncidentIterator with iedge_ member variable initialized
     *         to the degree of the node.
     */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, id_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return ((n.id_ == id_) && (n.graph_ == graph_));
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
      return (id_ < n.id_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type id_;
    // private constructor accessible in Graph class
    Node(const Graph* graph, const size_type id)
        : graph_(const_cast<Graph*>(graph)), id_(id) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return point_list_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return point_list_.size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position,
        const node_value_type& value = {}) {
    size_type s = size();
    Node new_node = Node(this, s);
    point_list_.push_back(position);
    node_value_list_.push_back(value);
    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return n.id_ < size();
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
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {}

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, node1_id_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, node2_id_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (e.id_ == id_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (id_ < e.id_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph_;
    size_type node1_id_; // id of node 1
    size_type node2_id_; // id of node 2
    size_type id_;       // id of the edge
    Edge(const Graph* graph, const size_type node1_id, const size_type node2_id,
        const size_type id)
        : graph_(const_cast<Graph*>(graph)), node1_id_(node1_id), node2_id_(node2_id),
          id_(id) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_list_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());
    return edge_list_[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert(has_node(a) && has_node(b));
    size_type id_min = std::min(a.id_, b.id_);
    size_type id_max = std::max(a.id_, b.id_);
    // check if an edge might exist in the first place
    if (adj_list_.find(id_min) != adj_list_.end()){
      for (size_type i : adj_list_.at(id_min))
        if (edge_list_[i].node2().id_ == id_max) return true;
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
  Edge add_edge(const Node& a, const Node& b) {
    assert(has_node(a) && has_node(b) && !(a == b));
    size_type id_min = std::min(a.id_, b.id_);
    size_type id_max = std::max(a.id_, b.id_);
    size_type i {find_edge(a, b)};

    if (i == num_edges()){
      // the edge does not exist, create it
      Edge new_edge = Edge(this, id_min, id_max, i);
      edge_list_.push_back(new_edge);
      adj_list_[id_min].push_back(i);
      adj_list_[id_max].push_back(i);
      return new_edge;
    }
    // the edge exists, return it
    else return edge_list_[i];
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    point_list_.clear();
    node_value_list_.clear();
    edge_list_.clear();
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
      return Node(graph_, node_id_);
    }

    /**
     * @brief Increment operator for the node iterator.
     * @return Itself after incrementing
     *
     * This function increments the node_id_ member variable.
     */
    node_iterator& operator++() {
      node_id_++;
      return *this;
    }

    /**
     * @brief Equality comparator
     * @return True if both iterators point to the same node in the same graph
     *         False otherwise
     */
    bool operator==(const node_iterator& ni) const {
      return (graph_ == ni.graph_) && (node_id_ == ni.node_id_);
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type node_id_;
    NodeIterator(const Graph* graph, const size_type node_id)
        : graph_(const_cast<Graph*>(graph)), node_id_(node_id) {}
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
      Edge& e = graph_->edge_list_[graph_->adj_list_.at(node_id_)[iedge_]];
      if (e.node1_id_ == node_id_) return e;
      else return Edge(graph_, e.node2_id_, e.node1_id_, e.id_);
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
      return (graph_ ==iit.graph_) && (node_id_ == iit.node_id_)
                                   && (iedge_ == iit.iedge_);
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type node_id_;
    size_type iedge_;
    IncidentIterator(const Graph* graph, const size_type node_id,
        const size_type iedge)
        : graph_(const_cast<Graph*>(graph)), node_id_(node_id), iedge_(iedge) {}
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
      return graph_->edge_list_[edge_id_];
    }

    /**
     * @brief Increment operator for the node iterator.
     * @details Itself after incrementing
     *
     * This function increments the edge_id_ member variable.
     */
    EdgeIterator& operator++() {
      edge_id_++;
      return *this;
    }

    /**
     * @brief Equality comparator
     * @return True if both iterators point to the sane edge in the same graph
     *         False otherwise
     */
    bool operator==(const EdgeIterator& ei) const {
      return (graph_ == ei.graph_) && (edge_id_ == ei.edge_id_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type edge_id_;
    EdgeIterator(const Graph* graph, const size_type edge_id)
        : graph_(const_cast<Graph*>(graph)), edge_id_(edge_id) {}
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

 private:

  /** For now, the graph is represented with an node and edge vectors. */
  std::vector<Point> point_list_;
  std::vector<V> node_value_list_;
  std::vector<edge_type> edge_list_;
  /** Kind of adjacency list
   * It relate a node with a vector of edge ids' of the strictly greater
   * nodes connected to it. Don't know yet if it is a good data structure,
   * will see in future assignments if I have to modify it.
   */
  std::map<size_type, std::vector<size_type>> adj_list_;


  /** Find an edge and return its id
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return either the id of the edge if it exists,
   *         otherwise the number of edges in the graph.
   */
  size_type find_edge(const Node& a, const Node& b){
    assert(has_node(a) && has_node(b) && !(a == b));
    size_type id_min = std::min(a.id_, b.id_);
    size_type id_max = std::max(a.id_, b.id_);
    if (adj_list_.find(id_min) != adj_list_.end()){
      for (size_type i : adj_list_.at(id_min))
        if (edge_list_[i].node2().id_ == id_max) return i;
    }
    return num_edges();
  }

};

#endif // CME212_GRAPH_HPP
