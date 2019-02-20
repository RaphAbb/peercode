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
template <typename V, typename E>
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
  std::vector<unsigned> node_itou; //maps indices to uids

  // For storage of Edges
  std::vector<Internal_Edge> edges_;
  std::map<unsigned, std::map<unsigned, unsigned>> adjacency_list;
  std::vector<unsigned> edge_itou; //maps indices to uids

 public:
  using node_value_type = V;
  using edge_value_type = E;

  //
  // PUBLIC TYPE DEFINITIONS
  //
  using size_type = unsigned;

  /** Type of this graph. */
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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph(){
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
    Point& position(){
      return (graph_->nodes_[uid_]).point_;
    }

    /** Return this node's position.*/
    const Point& position() const {
        return (graph_->nodes_[uid_]).point_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      Internal_Node internal_node = graph_->nodes_[uid_];
      return internal_node.index_;
    }

    /** Return the value of the node instance.
     *
     * @return value of the node
     * */
    node_value_type& value() {
      return (graph_->nodes_[uid_]).value_;
    }

    /** Return the value of the node instance.
     *
     * @return value of the node
     * */
    const node_value_type& value() const {
      return (graph_->nodes_[uid_]).value_;
    }

    /** Return the degree of the node instance.
     *
     * @return degree of the node
     */
    size_type degree() const {
      size_type node_id = uid_;
      return (graph_->adjacency_list.at(node_id)).size();
    }

    /** Return an iterator that points to the first incident edge of node.
     *
     * @return IncidentIterator
     */
    incident_iterator edge_begin() const {
      return IncidentIterator((graph_->adjacency_list.at(uid_)).begin(), graph_, uid_);
    }

    /**
     * Return an iterator that points to the last incident edge of node.
     *
     * @return IncidentIterator
     */
    incident_iterator edge_end() const {
      return IncidentIterator((graph_->adjacency_list.at(uid_)).end(), graph_, uid_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {

      bool equals = (graph_ == n.graph_) and (uid_ == n.uid_);
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
      return uid_ < n.uid_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    Graph* graph_;
    size_type uid_;

    /** Private Constructor **/
    Node(const Graph* graph, size_type uid)
              : graph_(const_cast<Graph*>(graph)), uid_(uid){
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_itou.size();
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
    unsigned next_index = node_itou.size();
    unsigned next_uid = nodes_.size();

    Internal_Node internal_node = Internal_Node(this, next_uid, position, value,next_index, true);

    // Add node to internal data structures
    nodes_.push_back(const_cast<Internal_Node&>(internal_node));
    node_itou.push_back(next_uid);

    // Make proxy node instance
    Node n = Node(this, next_uid);

    // Add node to adjacency list
    adjacency_list[next_uid] = std::map<unsigned, unsigned>();

    return n;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    Internal_Node internal_node = nodes_[n.uid_];
    // Check that the node is valid and that it belongs to the graph.
    return (internal_node.valid_ and (this == n.graph_));
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < num_nodes());
    unsigned uid = node_itou[i];
    Node n = Node(this, uid);
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
      Node node = Node(graph_, node1_id);

      return node;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      unsigned node2_id = node2_id_;
      Node node = Node(graph_, node2_id);

      return node;
    }

    /** Return the length of an edge */
    double length() const {
        return norm(node1().position() - node2().position());
    }

    /** Return a reference to a value of an edge.*/
    edge_value_type& value() {
        return (graph_->edges_[uid_].value_);
    }

    /** Return the value of an edge*/
    const edge_value_type value() const {
        return (graph_->edges_[uid_].value_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {

      if (graph_ != e.graph_) {
          return false;
      }

      bool equal_same_order = (node1_id_ == e.node1_id_) and (node2_id_ == e.node2_id_);
      bool equal_opposite_order = (node1_id_ == e.node2_id_) and (node2_id_ == e.node1_id_);

      return equal_same_order or equal_opposite_order;;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {

        // If the edges are not in the same graph compare the graphs.
      if (graph_ != e.graph_) {
          return graph_ < e.graph_;
      }
      return uid_ < e.uid_;;
    }

    /** Return the index of an edge
     *
     * Complexity: O(1)
     * */
    unsigned index() {
        Internal_Edge internal_edge = graph_->edges_[uid_];
        return internal_edge.index_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    Graph* graph_;
    unsigned uid_; //Unique ID over all created Edges.

    unsigned node1_id_;
    unsigned node2_id_;

    /** Private Constructor for Edge */
    Edge(const Graph* graph, unsigned index, unsigned node1_id, unsigned node2_id):
      graph_(const_cast<Graph*>(graph)), uid_(index), node1_id_(node1_id),
      node2_id_(node2_id) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    return edge_itou.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(1)
   */
  Edge edge(size_type i) const {

    assert(i < num_edges());

    unsigned uid = edge_itou[i];    // get index from index to uid vector
    Internal_Edge internal_edge = edges_[uid];
    unsigned node_a_index = internal_edge.node1_id_;
    unsigned node_b_index = internal_edge.node2_id_;

    Edge e = Edge(this, uid, node_a_index, node_b_index);

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

    unsigned node_a_uid = a.uid_;
    unsigned node_b_uid = b.uid_;

    // Check if the edge is represented in the adjacency list
    if (adjacency_list.count(node_a_uid)) {
        return (adjacency_list.at(node_a_uid)).count(node_b_uid);
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
   * Complexity: No more than O(num_nodes() + num_edges())
   */
  Edge add_edge(const Node& a, const Node& b) {

    assert(not(a == b));

    unsigned node_a_id = a.uid_;
    unsigned node_b_id = b.uid_;

    unsigned next_edge_uid = edges_.size();
    unsigned next_index = edge_itou.size();

    // If edge exists, find the edge index and return the edge
    if (has_edge(a, b)) {
      unsigned edge_uid = adjacency_list.at(a.uid_).at(b.uid_);
      return Edge(this, edge_uid, a.uid_, b.uid_);
    }

    // Add entry for key_value of a in the adjacency list
    adjacency_list[node_a_id][node_b_id] = next_edge_uid;
    // Add entry for key_value of b in the adjacency list
    adjacency_list[node_b_id][node_a_id] = next_edge_uid;

    // Add edge to internal edge_vector
    Internal_Edge internal_edge = Internal_Edge(this, next_edge_uid, node_a_id, node_b_id, next_index, true);
    edges_.push_back(internal_edge);
    edge_itou.push_back(next_edge_uid);

    // Make new edge to return
    Edge e = Edge(this, next_edge_uid , node_a_id, node_b_id);

    return e;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // Free storage of nodes and edges
    nodes_ = std::vector<Internal_Node>();
    edges_ = std::vector<Internal_Edge>();
    adjacency_list = std::map<unsigned, std::map<unsigned, unsigned>>();

    node_itou = std::vector<unsigned>();
    edge_itou = std::vector<unsigned>();
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
      unsigned uid = *iter_;
      return Node(graph_, uid);
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
    typename std::vector<unsigned>::const_iterator iter_;
    const Graph* graph_;

    // Private constructor
    NodeIterator(const Graph* graph, typename std::vector<size_type>::const_iterator iter) :
    iter_(iter),
    graph_(graph){
    }

  };

  /**
   * Return an iterator to the beginning of the nodes in the graph.
   */
  node_iterator node_begin() const {
    NodeIterator beginIterator = NodeIterator(this, node_itou.begin());
    return beginIterator;
  }

  /**
   * Return an iterator to the end of hte nodes in the graph.
   */
  node_iterator node_end() const {
    NodeIterator endIterator = NodeIterator(this, node_itou.end());
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

    /**
     * Dereferences the IncidentIterator
     * @return Edge pointed to by the IncidentIterator
     */
    Edge operator*() const {
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
      unsigned uid = *iter_;
      Internal_Edge internal_edge = graph_->edges_[uid];
      return Edge(internal_edge.graph_, internal_edge.uid_, internal_edge.node1_id_,
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
    typename std::vector<unsigned>::const_iterator iter_;

    // Private constructor
    EdgeIterator(const Graph* graph, typename std::vector<unsigned>::const_iterator iter) :
      graph_(graph), iter_(iter) {
      }
  };

  /**
   * Makes EdgeIterator that points to the first edge.
   * @return EdgeIterator
   */
  EdgeIterator edge_begin() const {
    return EdgeIterator(this, edge_itou.begin());
  }

  /**
   * Makes an EdgeIterator that points to the last edge.
   * @return EdgeIterator
   */
  EdgeIterator edge_end() const {
    return EdgeIterator(this, edge_itou.end());
  }

  /**
   * Removes node from the graph.
   * @param node
   * @post if the node was in the graph num_nodes() = num_nodes() - 1
   * @return 1 if the node was removed successfully and 0 otherwise.
   *
   * Complexity: O(log(num_nodes) + num_edges())
   *            The log factor is from look up in the adjacency list for the incident iterator.
   */
  unsigned remove_node(const Node& node) {
      // Check if node is in graph
      if (not has_node(node)){
          return 0;
      }

      // Remove incident edges
      while (node.edge_begin() != node.edge_end()) {
          Edge edge = *node.edge_begin();
          remove_edge(edge);
      }

      // Swap and pop
      unsigned node_uid = node.uid_;
      unsigned node_index = node.index();

      // Invalidate the node
      nodes_[node_uid].valid_ = false;

      // Swap the node with the last node
      node_itou[node_index] = node_itou[node_itou.size()-1];

      // Reset the internal node index of the swapped edge to the index of the deleted node
      nodes_[node_itou[node_index]].index_ = node_index;

      // Pop the end of the itou
      node_itou.pop_back();

      return 1;
  }

  /**
   * Removes a node given an iterator.
   * @param n_it
   * @post if the node was in the graph num_nodes() = num_nodes() - 1
   * @return a valid edge iterator.
   *
   * Complexity: O(log(num_nodes) + num_edges())
   *            The log factor is from look up in the adjacency list for the incident iterator.
   */
  node_iterator remove_node(node_iterator n_it){

      Node node = *n_it;
      remove_node(node);

      return node_begin();
  }

  /**
   * Remove an edge between two nodes.
   * @param node1
   * @param node2
   * @post if the edge was in the graph num_edges() = num_edges() - 1
   * @return 1 if the node was removed successfully and 0 otherwise.
   *
   * Complexity: O(log(num_nodes))
   */
  unsigned remove_edge(const Node& node1, const Node& node2) {

      if (!has_edge(node1, node2)){
          return 0;
      }

      // Look up edge uid to get the index
      unsigned edge_uid = (adjacency_list.at(node1.uid_)).at(node2.uid_);
      Edge edge = Edge(this, edge_uid, node1.uid_, node2.uid_);
      unsigned edge_index = edge.index();

      // Invalidate the edge
      edges_[edge_uid].valid_ = false;

      // Swap the edge with the last edge
      edge_itou[edge_index] = edge_itou[edge_itou.size()-1];

      // Reset the internal edge index of the swapped edge to the index
      // of the deleted edge
      edges_[edge_itou[edge_index]].index_ = edge_index;

      // pop the end of the itou
      edge_itou.pop_back();

      // Update adjacency list
      adjacency_list.at(node1.uid_).erase(node2.uid_);
      adjacency_list.at(node2.uid_).erase(node1.uid_);

      return 1;

  }

   /**
   * Remove an edge.
   * @param edge
   * @post if the edge was in the graph num_edges() = num_edges() - 1
   * @return 1 if the node was removed successfully and 0 otherwise.
   *
   * Complexity: O(log(num_nodes))
   */
  unsigned remove_edge(const Edge& edge) {
      Node node1 = edge.node1();
      Node node2 = edge.node2();
      return remove_edge(node1, node2);
  }

   /**
   * Remove an edge given an edge iterator.
   * @param node1
   * @param node2
   * @post if the edge was in the graph num_edges() = num_edges() - 1
   * @return a valid edge iterator.
   *
   * Complexity: O(log(num_nodes))
   */
  edge_iterator remove_edge(edge_iterator e_it) {
      Edge edge = *e_it;
      remove_edge(edge);
      return edge_begin();
  }

 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  /**
   * Class to hold information associated with nodes.
   */
  class Internal_Node {

      friend class Graph;

      const Graph *graph_;
      unsigned uid_;
      unsigned index_;

      bool valid_;

      Point point_;
      node_value_type value_;

      // Invalid constructor
      Internal_Node() {
      }

      Internal_Node(const Graph* graph, unsigned uid, Point point, node_value_type value, unsigned index,
              bool valid) :
          graph_(const_cast<Graph*>(graph)),
          uid_(uid),
          index_(index),
          valid_(valid),
          point_(point),
          value_(value) {
      };
  };

  /**
   * Class to hold information associated with edges.
   */
  class Internal_Edge {

      friend class Graph;

      const Graph *graph_;
      unsigned uid_;
      unsigned index_;

      bool valid_;

      unsigned node1_id_;
      unsigned node2_id_;

      edge_value_type value_;

      // Invalid constructor
      Internal_Edge() {
      }

      Internal_Edge(const Graph* graph, unsigned uid, unsigned node1_id, unsigned node2_id, unsigned index,
              bool valid):
              graph_(const_cast<Graph*>(graph)),
              uid_(uid),
              index_(index),
              valid_(valid),
              node1_id_(node1_id),
              node2_id_(node2_id){
      }
  };
};

#endif // CME212_GRAPH_HPP
