#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
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
 private:

  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  class Internal_Node;
  class Internal_Edge;

  // For storage of Nodes
  std::vector<Internal_Node> nodes_;

  unsigned next_index_;
  unsigned next_edge_index_;

  // For storage of Edges
  std::vector<Internal_Edge> edges_;
  std::map<unsigned, std::map<unsigned, unsigned>> adjacency_list;

 public:
  using node_value_type = V;
  //
  // PUBLIC TYPE DEFINITIONS
  //
  using size_type = unsigned;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph(): next_index_(0), next_edge_index_(0){
    // HW0: YOUR CODE HERE
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

    /** Return this node's position.*/
    const Point& position() const {
      return (graph_->nodes_[index_]).point_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {

      // With the current functionality the index will never be greater than graph_size
      // since we are not enabling the user to delete any nodes.
      unsigned index = index_;

      return index;
    }

    // HW1: YOUR CODE HERE
    // @Todo:Supply definitions AND SPECIFICATIONS for:

    /** Return the value of the node instance.
     * @return value of the node
     * */
    node_value_type& value() {
      return (graph_->nodes_[index_]).value_;
    }

    /** Return the value of the node instance.
     *
     * @return value of the node
     * */
    const node_value_type& value() const {
      return (graph_->nodes_[index_]).value_;
    }

    /** Return the degree of the node instance.
     *
     * @return degree of the node
     */
    size_type degree() const {
      size_type node_id = index_;
      return (graph_->adjacency_list.at(node_id)).size();
    }

    /** Return an iterator that points to the first incident edge of node.
     *
     * @return IncidentIterator
     */
    incident_iterator edge_begin() const {
      return IncidentIterator((graph_->adjacency_list.at(index_)).begin(), graph_, index_);
    }

    /**
     * Return an iterator that points to the last incident edge of node.
     *
     * @return IncidentIterator
     */
    incident_iterator edge_end() const {
      return IncidentIterator((graph_->adjacency_list.at(index_)).end(), graph_, index_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {

      bool equals = (graph_ == n.graph_) and (index_ == n.index_);
      return equals;
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

      return index_ < n.index();
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    Graph* graph_;
    size_type index_;

    /** Private Constructor **/
    Node(const Graph* graph, size_type index)
              : graph_(const_cast<Graph*>(graph)), index_(index){
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
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // Add point to graph data structure
    Internal_Node internal_node = Internal_Node(this, next_index_, position, value);
    nodes_.push_back(const_cast<Internal_Node&>(internal_node));

    // Make proxy node instance
    Node n = Node(this, next_index_);

    // Add node to adjacency list
    adjacency_list[next_index_] = std::map<unsigned, unsigned>();

    // Update graph info
    next_index_ = next_index_ + 1;

    return n;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {

    return ((nodes_[n.index()]).point_ == n.position()) and (this == n.graph_);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < num_nodes());
    Node n = Node(this, i);
    return n;
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
      unsigned node1_id = node1_id_;
      Node node = graph_->node(node1_id);

      return node;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      unsigned node2_id = node2_id_;
      Node node = graph_->node(node2_id);

      return node;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {

      bool equal_same_order = (node1_id_ == e.node1_id_) and (node2_id_ == e.node2_id_);
      bool equal_opposite_order = (node1_id_ == e.node2_id_) and (node2_id_ == e.node1_id_);

      bool equal = equal_same_order or equal_opposite_order;

      return equal;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {

      bool less = index_ < e.index_;

      return less;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    Graph* graph_;
    unsigned index_;

    unsigned node1_id_;
    unsigned node2_id_;

    /** Constructor for Edge */
    Edge(const Graph* graph, unsigned index, unsigned node1_id, unsigned node2_id):
      graph_(const_cast<Graph*>(graph)), index_(index), node1_id_(node1_id),
      node2_id_(node2_id) {
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
    Internal_Edge internal_edge = edges_[i];
    unsigned node_a_index = internal_edge.node1_id_;
    unsigned node_b_index = internal_edge.node2_id_;

    Edge e = Edge(this, i, node_a_index, node_b_index);

    return e;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {

    // Check if a and b are in this graph
    assert(a.graph_ == b.graph_);

    unsigned node_a_index = a.index();
    unsigned node_b_index = b.index();

    // Check if the edge is represented in the adjacency list
    if (adjacency_list.count(node_a_index)) {
      if ((adjacency_list.at(node_a_index)).count(node_b_index)) {
        return true;
      }
      else {
        return false;
      }
    }
    else {
      return false;
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

    assert(not(a == b));

    unsigned node_a_id = a.index();
    unsigned node_b_id = b.index();

    // If edge exists, find the edge index and return the edge
    if (has_edge(a, b)) {
      unsigned index = adjacency_list.at(a.index()).at(b.index());
      return Edge(this, index, a.index_, b.index_);
    }

    // Add entry for key_value of a in the adjacency list
    adjacency_list[node_a_id][node_b_id] = next_edge_index_;
    // Add entry for key_value of b in the adjacency list
    adjacency_list[node_b_id][node_a_id] = next_edge_index_;

    // Add edge to internal edge_vector
    Internal_Edge internal_edge = Internal_Edge(this, next_edge_index_, node_a_id, node_b_id);
    edges_.push_back(internal_edge);

    // Make new edge to return
    Edge e = Edge(this, next_edge_index_ , node_a_id, node_b_id);

    // Increment edge counter
    next_edge_index_++;

    return e;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // Set counters to 0
    next_index_ = 0;
    next_edge_index_ = 0;

    // Free storage of nodes and edges
    nodes_ = std::vector<Internal_Node>();
    edges_ = std::vector<Internal_Edge>();
    adjacency_list = std::map<unsigned, std::map<unsigned, unsigned>>();
  }

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

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    /**
     * Dereferences the node pointed to by the NodeIterator.
     * @return Node
     */
    Node operator*() const {
      Internal_Node internal_node = *iter_;
      return Node(internal_node.graph_, internal_node.index_);
    }

    /**
     * Increments the NodeIterator.
     * @return A reference to the NodeIterator after incrementing.
     */
    NodeIterator& operator++() {
      iter_++;
      return *this;
    }

    /**
     * Checks if two NodeIterators are equal.
     * @param x NodeIterator to compare to.
     * @return True if equal, False otherwise.
     */
    bool operator==(const NodeIterator& x) const {
      return (iter_ == x.iter_);
    }

   private:
    friend class Graph;
    typename std::vector<Internal_Node>::const_iterator iter_;

    // Private constructor
    NodeIterator(typename std::vector<Internal_Node>::const_iterator iter) :
    iter_(iter) {
    }

  };

  /**
   * Return an iterator to the beginning of the nodes in the graph.
   */
  node_iterator node_begin() const {
    NodeIterator beginIterator = NodeIterator(nodes_.begin());
    return beginIterator;
  }

  /**
   * Return an iterator to the end of hte nodes in the graph.
   */
  node_iterator node_end() const {
    NodeIterator endIterator = NodeIterator(nodes_.end());
    return endIterator;
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /**
     * Dereferences the IncidentIterator
     * @return Edge pointed to by the IncidentIterator
     */
    Edge operator*() const {
      // Fixme: iter is a map iterator it does not contain an internal edge.
      unsigned node_b_id = iter_->first;
      unsigned edge_id = iter_->second;
      return Edge(graph_, edge_id, node_a_id_, node_b_id);
    }

    /**
     * Increments the IncidentIterator
     * @return a reference to the IncidentIterator.
     */
    IncidentIterator& operator++() {
      iter_++;
      return *this;
    }

    /**
     * Checks if two IncidentIterators are equal.
     * @param other_iter : The IncidentIterator to compare to.
     * @return True if the IncidentIterators are equal, False otherwise.
     */
    bool operator==(const IncidentIterator& other_iter) const {
      return (iter_ == other_iter.iter_);
    }

   private:
    friend class Graph;
    typename std::map<unsigned, unsigned>::const_iterator iter_;
    Graph* graph_;
    unsigned node_a_id_;

    //Private constructor
    IncidentIterator(typename std::map<unsigned, unsigned>::const_iterator iter, Graph* graph,
            unsigned node_a_id) :
    iter_(iter), graph_(graph), node_a_id_(node_a_id){
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>  {
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
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /**
     * Dereferences the EdgeIterator.
     * @return Edge pointed to by the EdgeIterator.
     */
    Edge operator*() const {
      Internal_Edge internal_edge = *iter_;
      return Edge(internal_edge.graph_, internal_edge.index_, internal_edge.node1_id_,
              internal_edge.node2_id_);
    }

    /**
     * Increments the EdgeIterator.
     * @return Reference to the incremented EdgeIterator.
     */
    EdgeIterator& operator++() {
      iter_++;
      return *this;
    }

    /**
     * Checks if two EdgeIterators are equal.
     * @param other_iter : The other iterator to compare to.
     * @return True if the EdgeIterators are equal, False otherwise.
     */
    bool operator==(const EdgeIterator& other_iter) const {
      return iter_ == other_iter.iter_;
    }

   private:
    friend class Graph;
    const Graph* graph_;
    typename std::vector<Internal_Edge>::const_iterator iter_;

    // Private constructor
    EdgeIterator(const Graph* graph, typename std::vector<Internal_Edge>::const_iterator iter) :
      graph_(graph), iter_(iter) {
      }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /**
   * Makes EdgeIterator that points to the first edge.
   * @return EdgeIterator
   */
  EdgeIterator edge_begin() const {
    return EdgeIterator(this, edges_.begin());
  }

  /**
   * Makes an EdgeIterator that points to the last edge.
   * @return EdgeIterator
   */
  EdgeIterator edge_end() const {
    return EdgeIterator(this, edges_.end());
  }

 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  class Internal_Node {

      friend class Graph;

      const Graph *graph_;
      unsigned index_;

      Point point_;
      node_value_type value_;

      // Invalid constructor
      Internal_Node() {
      }

      Internal_Node(const Graph* graph, unsigned index, Point point, node_value_type value) :
          graph_(const_cast<Graph*>(graph)), index_(index), point_(point), value_(value) {
      };
  };

  class Internal_Edge {

      friend class Graph;

      const Graph *graph_;
      unsigned index_;

      unsigned node1_id_;
      unsigned node2_id_;

      // Invalid constructor
      Internal_Edge() {
      }

      Internal_Edge(const Graph* graph, unsigned index, unsigned node1_id, unsigned node2_id):
              graph_(const_cast<Graph*>(graph)), index_(index), node1_id_(node1_id),
              node2_id_(node2_id) {
      }
  };

};

#endif // CME212_GRAPH_HPP
