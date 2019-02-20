#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <map>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E=int>
class Graph {
 private:
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

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
  /** Added in HW2-3. */
  typedef V node_value_type;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  /** Added in HW2-3 */
  typedef E edge_value_type;

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
  Graph() : nodes_data_(), edges_data_(), size_(0), edges_(), adjacency_(), ind_to_edges_(), n_edges_(0) {

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
      return graph_->nodes_data_[node_id_].position;
    }

    /* Return's modifiable node position */
    Point& position() {
      return graph_->nodes_data_[node_id_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return this->node_id_;
    }

    /**
     * @brief Node value type specified by the user, of type V.
     *
     * @return Value of the current node.
     *
     * @pre Current node_id_ is in the range [0, graph_size).
     *
     */
    node_value_type& value() {
      return graph_->nodes_data_[node_id_].value;
    }

    /**
     * @brief Node value type specified by the user, of type V, that cannot be modified.
     *
     * @return Value of the current node.
     *
     * @pre Current node_id_ is in the range [0, graph_size).
     * @post Output cannot be modified.
     *
     */
    const node_value_type& value() const {
      return graph_->nodes_data_[node_id_].value;
    }

    /**
     * @brief Returns the number of incident edges to a node.
     *
     * @return Number of nodes incident to the current node.
     *
     * @pre Current node_id_ is valid.
     * @post Output is in the range [0, graph_size).
     *
     */
    size_type degree() const {
      return graph_->adjacency_[node_id_].size();
    }

    /**
     * @brief Pointer to the first edge incident to a given node.
     *
     * @return Pointer to the iterator at the first indicent node.
     *
     * O(1) complexity
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, node_id_, 0);
    }

    /**
     * @brief Pointer to the element after the last edge incident to a given node.
     *
     * @return Pointer to the iterator at the element after the last incident node.
     *
     * O(1) complexity
     */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, node_id_, degree());
    }

    /**
     * @brief Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // Check that indices and graphs are equal
      if ((this->graph_ == n.graph_) and (this->node_id_ == n.node_id_)) {
        return true;
      }

      // Return false otherwise
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
      // Check if ordering is verified
      if (node_id_ < n.node_id_) {
        return true;
      }

      // Return false otherwise
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Pointer back to the graph
    Graph* graph_;

    // This node's identifcation number
    size_type node_id_;

    // Private constructor
    Node(const Graph* graph, size_type node_id)
      : graph_(const_cast<Graph*>(graph)), node_id_(node_id) {
      }

    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return size_;
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
    // New node id
    size_type new_node_id = size_;

    // New internal node
    struct internal_node new_internal_node = {new_node_id, position, value};

    // Push new node
    nodes_data_.push_back(new_internal_node);

    // Initialize adjacency list
    adjacency_.push_back(std::vector<size_type> ());

    // Initialize dictionary of incident_nodes : edge_id associated to the new node
    edges_.push_back(std::map<size_type, size_type> ());

    // Increment size
    size_ += 1;

    // Return added node
    return Node(this, new_node_id);
  }

  /**
   * @brief Removes node _n_ (if present).
   *
   * @param n Valid Node object.
   * @return 1 if the node removal was successful (node was indeed present in the graph),
   *         0 if no node removal occurred (no such node in the graph).
   *
   * @pre node _n_ may or may not exist in the graph.
   * @post if a node _n_ of degree _d_ existed, then size_[after] = size_[before] - 1
   *                                             and n_edges_[after] = n_edges_[before] - d
   *                              otherwise, we have size_[after] = size_[before]
   *                                             and n_edges_[after] = n_edges_[before]
   * @post has_node(n) == false.
   * @post if a node _n_ existed, then it is erased from nodes_data_, edges_,
   *       adjacency_, ind_to_edges_.
   * @post if a node _n_ existed, then all nodes of index greater than it are decreased by 1.
   * @post Node iterator associated to node _n_ (if applicable) is invalidated.
   * @post Edge and incident iterators associated to node _n_ (if applicable) are invalidated.
   *
   * The complexity of the remove node algorithm is O(num_nodes + max_degree * complexity(remove_edge)).
   */
  size_type remove_node(const Node& n) {
    // Check if node exists
    if (has_node(n) == false) {
      return 0;
    }

    // Parameters
    size_type node_id = n.index();

    // Delete associated data
    nodes_data_.erase(nodes_data_.begin() + node_id);

    // Remove incident edges to current node
    while (adjacency_[node_id].size() > 0) {
      // Retrieve corresponding edge ID
      size_type incident_node_id = adjacency_[node_id][0];

      // Remove corresponding edge
      remove_edge(node(node_id), node(incident_node_id));
    }

    // Delete entry associated to current node in edges_
    edges_.erase(edges_.begin() + node_id);

    // Decrement id of nodes with higher id in edges_
    for (size_type i = 0; i < edges_.size(); i++) {
      // Create temporary map
      std::map<size_type, size_type> temp;

      // Loop over previous map
      for (auto const& x : edges_[i]) {
        // Local parameters
        size_type j = x.first;
        size_type local_edge_id = x.second;

        // Decrease other node id if applicable
        if (j > node_id) {
          j -= 1;
        }

        // Store information in updated map
        temp[j] = local_edge_id;
      }

      // Replace old map by updated map
      edges_[i] = temp;
    }

    // Delete entry associated to current node in adjacency_
    adjacency_.erase(adjacency_.begin() + node_id);

    // Decrement id of nodes with higher id in adjacency_
    for (size_type i = 0; i < adjacency_.size(); i++) {
      for (size_type j = 0; j < adjacency_[i].size(); j++) {
        // Decrement id if applicable
        if (adjacency_[i][j] > node_id) {
          adjacency_[i][j] -= 1;
        }
      }
    }

    // Decrement id of nodes with higher id in ind_to_edges_
    for (size_type local_edge_id = 0; local_edge_id < ind_to_edges_.size(); local_edge_id++) {
      // Decrement id of first element if applicable
      if (ind_to_edges_[local_edge_id].first > node_id) {
        ind_to_edges_[local_edge_id].first -= 1;
      }

      // Decrement id of second element if applicable
      if (ind_to_edges_[local_edge_id].second > node_id) {
        ind_to_edges_[local_edge_id].second -= 1;
      }
    }

    // Decrease number of nodes
    size_ -= 1;

    // Return true
    return 1;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // Check if node belongs to same graph
    if (n.graph_ == this) {
      return true;
    }

    // Return false otherwise
    else {
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
    // Verify input
    assert((0 <= i) and (i < this->size_));

    // Return corresponding node
    return Node(this, i);
  }

////////////////////////////////////////////////////////////////////////////////

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
      return Node(graph_, node1_id_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, node2_id_);
    }

    /** Return the length of the edge (L2 norm of between the two nodes) */
    double length() const {
      return norm(node1().position() - node2().position());
    }

    /** Access edge's modifiable value */
    edge_value_type& value() {
      return graph_->edges_data_[edge_id_].value;
    }

    /** Access edge's value */
    const edge_value_type& value() const {
      return graph_->edges_data_[edge_id_].value;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
     // Case where both edges are part of the same graph
     if (this->graph_ == e.graph_) {
       // Check if edges are the same
       if ((node1_id_ == e.node1_id_ and node2_id_ == e.node2_id_) or
           (node1_id_ == e.node2_id_ and node2_id_ == e.node1_id_)) {
             return true;
           }
     }

     // All other situations
     return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // Case both belong to the same graph
      if (graph_ == e.graph_) {
        return (this->edge_id_ < e.edge_id_);
      }

      // Case do not belong to same graph
      else {
        return graph_ < e.graph_;
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Pointer back to the graph
    Graph* graph_;

    // This edge's identifcation number
    size_type edge_id_;

    // The two ends of the edge
    size_type node1_id_;
    size_type node2_id_;

    // Private constructor
    Edge(const Graph* graph, size_type edge_id, size_type node1_id, size_type node2_id)
      : graph_(const_cast<Graph*>(graph)), edge_id_(edge_id), node1_id_(node1_id), node2_id_(node2_id) {
      }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return n_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // Check if requirement satisfied
    assert((0 <= i) and (i < this->num_edges()));

    // Return corresponding edge
    return Edge(this, i, ind_to_edges_[i].first, ind_to_edges_[i].second);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Check if assertion correct
    assert(a.graph_ == b.graph_);

    // Select the two end nodes
    size_type node1_id = std::min(a.node_id_, b.node_id_);
    size_type node2_id = std::max(a.node_id_, b.node_id_);

    // Edge not found
    if (edges_[node1_id].find(node2_id) == edges_[node1_id].end()) {
      return false;
    }

    // Key found otherwise
    return true;
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {
    // Select the two end nodes
    size_type node1_id = std::min(a.node_id_, b.node_id_);
    size_type node2_id = std::max(a.node_id_, b.node_id_);

    // Do nothing if edge already exists
    if (this->has_edge(a, b)) {

    }

    // Otherwise, create it
    else {
      // Otherwise, create new edge
      edges_[node1_id][node2_id] = n_edges_;
      edges_[node2_id][node1_id] = n_edges_;

      // Store adjacency information in both directions
      adjacency_[node1_id].push_back(node2_id);
      adjacency_[node2_id].push_back(node1_id);

      // Map index to edge
      ind_to_edges_.push_back({node1_id, node2_id});

      // New internal edge
      struct internal_edge new_internal_edge = {n_edges_, value};

      // Push new node
      edges_data_.push_back(new_internal_edge);

      // Increment number of edges
      n_edges_ += 1;
    }

    // Return edge
    return Edge(this, edges_[node1_id].at(node2_id), node1_id, node2_id);
  }

  /**
   * @brief Removes edge of extremities _a_ and _b_ (if present).
   *
   * @param a Valid Node object.
   * @param b Valid Node object.
   * @return 1 if the edge removal was successful (edge was indeed present in the graph),
   *         0 if no edge removal occurred (no such edge in the graph).
   *
   * @pre edge of extremities _a_ and _b_ may or may not exist in the graph.
   * @post if an edge of extremities _a_ and _b_ existed, then n_edges_[after] = n_edges_[before] - 1
   *                                        otherwise, we have n_edges_[after] = n_edges_[before]
   * @post has_edge(a, b) == false.
   * @post if an edge of extremities _a_ and _b_ existed, then it is erased from edges_data_, edges_,
   *       adjacency_, ind_to_edges_.
   * @post Edge iterator associated to the edge linking _a_ and _b_ (if applicable) is invalidated.
   * @post Incident iterators associated to nodes _a_ and _b_ (if applicable) are invalidated.
   *
   * The complexity of the remove edge algorithm is O(num_edges).
   */
  size_type remove_edge(const Node& a, const Node& b) {
    // Return false if edge does not exist
    if (has_edge(a, b) == false) {
      return 0;
    }

    // Parameters
    size_type node1_id = a.index();
    size_type node2_id = b.index();
    size_type edge_id = edges_[node1_id][node2_id];

    // Erase edge data from edges_data_
    edges_data_.erase(edges_data_.begin() + edge_id);

    // Erase edge index from edges_, both ways
    edges_[node1_id].erase(node2_id);
    edges_[node2_id].erase(node1_id);

    // Erase adjacency information from adjacency_, both ways
    adjacency_[node1_id].erase(std::remove(adjacency_[node1_id].begin(), adjacency_[node1_id].end(), node2_id), adjacency_[node1_id].end());
    adjacency_[node2_id].erase(std::remove(adjacency_[node2_id].begin(), adjacency_[node2_id].end(), node1_id), adjacency_[node2_id].end());

    // Remove edge entry from ind_to_edges_
    ind_to_edges_.erase(ind_to_edges_.begin() + edge_id);

    // Decrement id of subsequent edges by 1 from edges_
    for (size_type i = 0; i < edges_.size(); i++) {
      // Loop over previous map
      for (auto const& x : edges_[i]) {
        // Local parameters
        size_type j = x.first;
        size_type local_edge_id = x.second;

        // Decrease id if applicable
        if (local_edge_id > edge_id) {
          local_edge_id -= 1;
        }

        // Store updated id
        edges_[i][j] = local_edge_id;
      }
    }

    // Decrease number of edges
    n_edges_ -= 1;

    // Return true
    return 1;
  }


  /**
   * @brief Removes edge _e_ from the graph (if present).
   *
   * @param _e_, edge object.
   * @return 1 if the edge removal was successful (edge was indeed present in the graph),
   *         0 if no edge removal occurred (no such edge in the graph).
   *
   * @pre _e_ is a valid edge object.
   * @post please see above.
   *
   * For extended documentation about operations, complexity and @post conditions
   * entailed by this function, please refer to the 'size_type remove_edge' function
   * above.
   */
  size_type remove_edge(const Edge& e) {
    // Retrieve nodes associated to edge
    Node a = e.node1();
    Node b = e.node2();

    // Remove corresponding edge
    return remove_edge(a, b);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // Clear all nodes
    nodes_data_.clear();

    // Clear all edges
    edges_data_.clear();

    // Update graph size to 0
    size_ = 0;

    // Clear all edges
    edges_.clear();

    // Clear adjacency list
    adjacency_.clear();

    // Clear edge indices to edges
    ind_to_edges_.clear();

    // Update number of edges to 0
    n_edges_ = 0;
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
     * @brief Retrieve pointer to the Node of the current node index.
     *
     * O(1) complexity.
     */
    Node operator*() const {
      return graph_->node(index_);
    }

    /**
     * @brief Increment index of node and return corresponding node iterator.
     *
     * O(1) complexity.
     */
    NodeIterator& operator++() {
      // Increment index
      index_ += 1;

      // Return corresponding node iterator
      return *this;
    }

    /**
     * @brief Checks equality between two node iterators.
     *
     * O(1) complexity.
     */
    bool operator==(const NodeIterator& other_node_iter) const {
      // Checks if belong to same graph
      if (graph_ != other_node_iter.graph_) {
        return false;
      }

      // Checks if point to same index
      if (index_ != other_node_iter.index_) {
        return false;
      }

      // Two quantities are equal otherwise
      return true;
    }


   private:
    friend class Graph;

    // Pointer back to the graph
    Graph* graph_;

    // The iterator's index value
    size_type index_;

    // Private constructor
    NodeIterator(const Graph* graph, size_type index)
      : graph_(const_cast<Graph*>(graph)), index_(index) {
      }
  };


  /**
   * @brief Pointer to the first node of the graph.
   *
   * O(1) complexity.
   */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /**
   * @brief Pointer to the element after the last one of the graph.
   *
   * O(1) complexity.
   */
  node_iterator node_end() const {
    return NodeIterator(this, num_nodes());
  }

  /**
   * @brief Pointer to the node to be removed.
   *
   * @param _n_it_, node iterator pointing to the Node object to be removed.
   * @return node iterator on removed node.
   *
   * @pre valid node iterator pointing to a node that may or may not belong to
   * the graph.
   *
   * For extended documentation about operations, complexity and @post conditions
   * entailed by this function, please refer to the 'size_type remove_node' function.
   */
  node_iterator remove_node(node_iterator n_it) {
    // Extract node to remove
    Node n = *n_it;

    // Remove node
    remove_node(n);

    // Return iterator
    return n_it;
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
     * @brief Pointer to the edge corresponding to the current incident node.
     *
     * O(1) complexity.
     */
    Edge operator*() const {
      // Find incident neighbor
      size_type neighbor_index = graph_->adjacency_[root_index_][counter_];

      // Find edge ID
      size_type edge_id = graph_->edges_[root_index_][neighbor_index];

      // Return corresponding edge
      return Edge(graph_, edge_id, root_index_, neighbor_index);
    }

    /**
     * @brief Increment counter of incident nodes to the root.
     *
     * @pre counter_ in [0, number of incident nodes).
     * @post counter_ in [1, number of incident nodes].
     *
     * O(1) complexity.
     */
    IncidentIterator& operator++() {
      // Increment counter
      counter_ += 1;

      // Return corresponding iterator
      return *this;
    }


    /**
     * @brief Checks whether two incident iterators are equal.
     *
     * O(1) complexity.
     */
    bool operator==(const IncidentIterator& other_incident_iterator) const {
      // Check if belong to the same graph
      if (graph_ != other_incident_iterator.graph_) {
        return false;
      }

      // Check if iteration on neighbors of a same node
      if (root_index_ != other_incident_iterator.root_index_) {
        return false;
      }

      // Check if currently pointing at the same incident node
      if (counter_ != other_incident_iterator.counter_) {
        return false;
      }

      // Equality if passes all checks
      return true;
    }


   private:
    friend class Graph;

    // Pointer back to the graph
    Graph* graph_;

    // The iterator's root node index value
    size_type root_index_;

    // Edge counter
    size_type counter_;

    // Private constructor
    IncidentIterator(const Graph* graph, size_type root_index, size_type counter)
      : graph_(const_cast<Graph*>(graph)), root_index_(root_index), counter_(counter) {
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

    /**
     * @brief Returns pointer to the current edge.
     *
     * @pre counter_ between [0, number of edges in the graph).
     *
     * O(1) complexity.
     */
    Edge operator*() const {
      // Retrieve information relative to edge
      std::pair<size_type, size_type> nodes = graph_->ind_to_edges_[counter_];

      // Get node 1
      size_type node_1 = std::get<0>(nodes);

      // Get node 2
      size_type node_2 = std::get<1>(nodes);;

      // Return corresponding edge
      return Edge(graph_, counter_, node_1, node_2);
    }

    EdgeIterator& operator++() {
      // Increment counter
      counter_ += 1;

      // Return self
      return *this;
    }

    /**
     * @brief Checks whether two edge iterators are the same.
     *
     * O(1) complexity.
     */
    bool operator==(const EdgeIterator& other_edge_iterator) const {
      // Check if same graph
      if (graph_ != other_edge_iterator.graph_) {
        return false;
      }

      // Check if same index
      if (counter_ != other_edge_iterator.counter_) {
        return false;
      }

      // Passed all sanity checks
      return true;
    }


   private:
    friend class Graph;

    // Pointer back to the graph
    Graph* graph_;

    // Current edge ID
    size_type counter_;

    // Initialization
    EdgeIterator(const Graph* graph, size_type counter)
      : graph_(const_cast<Graph*>(graph)), counter_(counter) {
      }
  };


  /**
   * @brief Return iterator pointing on beginning edge.
   *
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /**
   * @brief Return iterator pointing on edge after last.
   *
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

  /**
   * @brief Pointer to the edge to be removed.
   *
   * @param _e_it_, Edge iterator pointing to the Edge object to be removed.
   * @return edge iterator on removed node.
   *
   * @pre valid edge iterator pointing to a node that may or may not belong to
   * the graph.
   *
   * For extended documentation about operations, complexity and @post conditions
   * entailed by this function, please refer to the 'size_type remove_edge' function.
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    // Extract edge to remove
    Edge e = *e_it;

    // Remove edge
    remove_edge(e);

    // Return iterator
    return e_it;
  }


 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // Internal type for node elements
  struct internal_node {
    // Node id
    size_type node_id;

    // Position of node
    Point position;

    // Node value
    node_value_type value;
  };

  // Internal type for edge elements
  struct internal_edge {
    // Edge id
    size_type edge_id;

    // Edge value
    edge_value_type value;
  };

  // List of nodes with associated data
  std::vector<internal_node> nodes_data_;

  // List of edges with associated data
  std::vector<internal_edge> edges_data_;

  // Size of graph
  size_type size_;

  // Edges as a list of dictionaries. Each index node1 of the list is a
  // dictionary of keys node2 with associated values edge_index.
  // E.g. edges_[node1_id][node2_id] = edge_id
  std::vector<std::map<size_type, size_type>> edges_;

  // Edges as a list of lists. Each key node lists the nodes adjacent to
  // it.
  std::vector<std::vector<size_type>> adjacency_;

  // Edge indices to edges
  std::vector<std::pair<size_type, size_type>> ind_to_edges_;

  // Number of edges
  size_type n_edges_;
};

#endif // CME212_GRAPH_HPP
