#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <string>
#include <set>
#include <unordered_map>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  // Predeclare the internal node and edge storage containers
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

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    add_nid_ = 1;
    add_eid_ = 1;
    num_nodes_ = 0;
    num_edges_ = 0;
  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES SECTION
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node {
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
      // HW0: YOUR CODE HERE
      graph_ = nullptr;
      nid_ = 0;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return fetch_node().position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return fetch_index();
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (graph_ == n.graph_ && this->index() == n.index());
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
      // HW0: YOUR CODE HERE
      return this->index() < n.index();
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    graph_type* graph_;     // Pointer reference to the parent Graph it exists within
    size_type nid_;         // Unique ID for Node

    // // Implement a private constructor for the node
    Node(const graph_type* graph, size_type nid)
      : graph_(const_cast<graph_type*>(graph)), nid_(nid) 
      {
      }

    // Retrieve the index of node with appropriate unique ID
    size_type fetch_index() const {
      return graph_->nid_to_idx_[nid_];
    }

    // Retrieve the node in container with appropriate unique ID
    internal_node& fetch_node() const {
      size_type idx = fetch_index();
      return graph_->nodes_[idx];
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return num_nodes_;
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
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    
    nid_to_idx_[add_nid_] = num_nodes_; // ID of created Node is mapped to end index in Graph
    internal_node n = {add_nid_, position};
    nodes_.push_back(n);

    ++num_nodes_;
    ++add_nid_;

    return Node(this, add_nid_ - 1);

  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE\
    // Node ID is always positively incremented, so check if it is less than current counter
    return (this == n.graph_ && n.nid_ < add_nid_);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(0 < i < num_nodes());
    size_type nid = nodes_[i].nid;
    return Node(this, nid);
  }

  //
  // EDGES SECTION
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
      graph_ = nullptr;
      eid_ = 0;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      size_type idx = fetch_index();
      std::set<size_type> pair = graph_->edges_[idx].pair;
      size_type nid1 = *(pair.begin());
      return Node(graph_, nid1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      size_type idx = fetch_index();
      std::set<size_type> pair = graph_->edges_[idx].pair;
      size_type nid2 = *(pair.rbegin());
      return Node(graph_, nid2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      node_type n1 = this->node1();
      node_type n2 = this->node2();
      if (n1 == e.node1() && n2 == e.node2()) {
        return true;
      } else {
        return false;
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return this->fetch_index() < e.fetch_index();
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    graph_type* graph_;     // Pointer reference to the parent Graph it exists within
    size_type eid_;         // Unique ID of Edge

    // Implement a private constructor for the edge
    Edge(const graph_type* graph, size_type eid)
      : graph_(const_cast<graph_type*>(graph)), eid_(eid)
      {
      }

    /** Returns the index of this Edge in teh Graph
     *  @return @a i corresponding to index of Edge in Graph container
     */
    size_type fetch_index() const {
      // HW0: YOUR CODE HERE
      return graph_->eid_to_idx_[eid_];
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(0 < i < num_edges());
    size_type eid = edges_[i].eid;
    return Edge(this, eid);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    std::string pair_key = pair_string(a, b);
    auto search = pair_to_eid_.find(pair_key);
    return search != pair_to_eid_.end() ? true : false;
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
    // HW0: YOUR CODE HERE

    std::string pair_key = pair_string(a, b);
    std::set<size_type> pair = {a.nid_, b.nid_};

    if (has_edge(a, b)) {

      size_type eid = pair_to_eid_[pair_key];
      return Edge(this, eid);

    } else {

      internal_edge e = {add_eid_, pair};
      edges_.push_back(e);

      eid_to_idx_[add_eid_] = num_edges_;
      pair_to_eid_[pair_key] = add_eid_;

      ++num_edges_;
      ++add_eid_;

      return Edge(this, add_eid_ - 1); // Returns Node proxy object with ID reference

    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    edges_.clear();

    nid_to_idx_.clear();
    eid_to_idx_.clear();
    pair_to_eid_.clear();

    // Reset counters numbering items
    add_nid_ = 1;
    add_eid_ = 1;
    num_nodes_ = 0;
    num_edges_ = 0;
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // Holder for actual data of the node
  struct internal_node {
    size_type nid;
    Point position;
  };

  // Holder for actual data of the node
  struct internal_edge {
    size_type eid;
    std::set<size_type> pair; // Ordered pair of Node IDs
  };

  // Counters for the Node/Edge IDs and element sizes
  size_type add_nid_;
  size_type add_eid_;
  size_type num_nodes_;
  size_type num_edges_;

  // Vectors holder actual data objects
  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;

  // Hash table for unique ID of node -> index in Graph storage container;
  std::unordered_map<size_type, size_type> nid_to_idx_;
  // Hash table for unique ID of edge -> index in Graph storage container;
  std::unordered_map<size_type, size_type> eid_to_idx_;
  // Hash table of pair of node IDs -> edge ID
  std::unordered_map<std::string, size_type> pair_to_eid_;

   /** Create a unique ID string for a pair of nodes that form an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return Pair string containing @a a.nid_ "-" + @a b.nid_
   */
  std::string pair_string(const Node& a, const Node& b) const {
    size_type nid1 = a.nid_;
    size_type nid2 = b.nid_;
    if (nid1 < nid2) {
      return std::to_string(nid1) + "-" + std::to_string(nid2);
    } else {
      return std::to_string(nid2) + "-" + std::to_string(nid1);
    }
  }

};

#endif // CME212_GRAPH_HPP
