#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
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

  // predeclaring the internal node struct
  struct internal_node;

 private:

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

  Graph()
      : intnodes_(), size_(0) {
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
    // empty constructor for invalid node
    Node() {
    
    }

    /** Return this node's position. */
    const Point& position() const {

      // grab the id'th index of the vector of internal nodes
      // and return the point value of that internal node
      return graph_->intnodes_[uid_node_].point_int_node;

    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {

      // return the ID of that node
      return uid_node_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */

    bool operator==(const Node& n) const {

      // check that the index of this node and the index of the input
      // node are the same, and that they refer to the same graph
      if (uid_node_ == n.uid_node_ and graph_ == n.graph_) {
        return true;
      }
      else {
        return false;
      }
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

      // check that the two indices are the same and that the graph
      // is the same
      if (index() < n.index() and graph_ == n.graph_) {
        return true;
      }
      else {
        return false;
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    /** Private Constructor */
    Node(const Graph* g, size_type i)
        : graph_(const_cast<Graph*>(g)), uid_node_(i) {
    }

    // Pointer back to the graph container
    const Graph* graph_;

    // Node should take in a uid always
    size_type uid_node_;

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {

    // return the size of the internal node vector
    return intnodes_.size();
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

    // push back a new node to the vector, with the 
    // position and size attributes
    // NOTE: the reason that this works is because in the internal
    // node struct, the order of the attributes is the same as the
    // order that the arguments are passed in the push back
    intnodes_.push_back({position, size_});

    size_++;

    return Node(this, size_-1);

  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {

    // if the index of the passed node is between 0 and 
    // the size of the internal node vector, it is a valid query
    if (n.uid_node_ < size_ and n.uid_node_ >= 0) {
      return true;
    }
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

    // make sure that the index is lower than the size of the
    // vector, and if so, return that node.
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

  class Edge {
   public:

    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      
      // at function finds the value associated with key
      // return the FIRST node in the vector associated with key
      return Node(graph_, graph_->umap_byedgeid.at(uid_edge_)[0]);

    }

    /** Return the other node of this Edge */
    Node node2() const {

      // at function finds the value associated with key
      // return the SECOND node in the vector associated with key
      return Node(graph_, graph_->umap_byedgeid.at(uid_edge_)[1]);

    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {

      // check that the ids are the same and that the graphs are the same
      if (uid_edge_ == e.uid_edge_ and graph_ == e.graph_) {
        return true;
      }
      else {
        return false;
      }

    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {

      // check that the id is less than input and that the graphs are the same
      if (uid_edge_ < e.uid_edge_ and graph_ == e.graph_) {
        return true;
      }
      else {
        return false;
      }

    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    /** Private Constructor */
    // i is the edge ID
    Edge(const Graph* g, size_type i)
        : graph_(const_cast<Graph*>(g)), uid_edge_(i) {
    }

    // Pointer back to the graph container
    const Graph* graph_;

    // Node should take in a uid always
    size_type uid_edge_;

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {

    // return the size of the umap
    return umap_byedgeid.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {

    // instantiate new edge object with this graph and with the index i
    return Edge(this, i);

  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {

    // Check if the key of node a's ID has an associated ID for node b
    // in the nested umap
    if (umap_bynodeid.count(a.uid_node_) > 0
        and umap_bynodeid.at(a.uid_node_).count(b.uid_node_) > 0) {
      return true;
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

    size_type a_id = a.uid_node_;
    size_type b_id = b.uid_node_;

    if (has_edge(a, b) == true) {
        // if the edge already exists, return that edge
        size_type edge_id = umap_bynodeid.at(a_id).at(b_id);
        return Edge(this, edge_id);
    }
    else {
        // add to umap of edges by id
        std::vector<size_type> node_vector{a_id, b_id};
        umap_byedgeid[num_edges()] = node_vector;

        // add to umap of edges by node
        std::unordered_map<size_type, size_type> temp_map1;
        temp_map1[a_id] = num_edges();

        std::unordered_map<size_type, size_type> temp_map2;
        temp_map2[b_id] = num_edges();

        umap_bynodeid[b_id] = temp_map1;
        umap_bynodeid[a_id] = temp_map2;

        return Edge(this, num_edges());
    }

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {

    // clearing all containers
    intnodes_.clear();
    umap_byedgeid.clear();
    umap_bynodeid.clear(); 

  }


 private:

  // Internal type for set elements
  struct InternalNode {
    Point point_int_node;   // The text held by an element
    size_type uid_int_node;      // The unique identifcation for an element
  };

  std::vector<InternalNode> intnodes_;

  // FOR EDGES
  // first unordered map: key is the edge id, values are a set
  // with the two nodes
  std::unordered_map<size_type, std::vector<size_type>> umap_byedgeid;
  // second unordered map: key is a NODE id, value is ANOTHER map where the key
  // is the SECOND node id, and the value is the edge id
  std::unordered_map<size_type, 
  std::unordered_map<size_type, size_type>> 
  umap_bynodeid;

  // size of the graph
  size_type size_;

};

#endif // CME212_GRAPH_HPP