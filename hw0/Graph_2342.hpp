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
class Graph {
 private:

  /** Predeclaration of the internal_node class. */
  struct internal_node;

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
    : nodes_vector() {
  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Lightweight proxy class to the graph's internal_node class.
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
    // Default public constructor. Creates an invalid Node
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      // dereference the graph pointer to access the vector of internal nodes. 
      // Then get the element with this Node's index and return the point/position
      return graph_->nodes_vector[uid_].point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // return the uid_ because currently the uid == the index
        // (note: assumption will likely change in later assignments)
      return uid_; 
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index (unique id in this case).
     */
    bool operator==(const Node& n) const {
      // check that the 2 nodes have pointers to the same graph and the same uid
      if (graph_ == n.graph_ && uid_ == n.uid_) {
        return true;
      }
      // if the two nodes are from different graphs or have different uid's
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
      // check that the 2 nodes have pointers to the same graph and that this
        // node's uid is smaller than the uid of the specified node
      if (graph_ == n.graph_ && uid_ < n.uid_) {
        return true;
      }
      else {
        return false;
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Pointer back to the Graph container
    Graph* graph_;
    // The node's unique identification number
    size_type uid_;
    // Private Constructor
    Node(const Graph* graph, size_type uid) 
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // get the current length of the vector containing internal nodes
    return nodes_vector.size();
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
    // get the index for a new node
    size_type idx = num_nodes();
    // create a new instance of internal node using the position specified
      // and add it to our internal node vector
    nodes_vector.push_back(internal_node(position, idx));
    // return a proxy Node that points to the newly created internal node
    return Node(this, idx);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // check if the specified Node object has a pointer to this graph and its
      // unique id exists in this graph
    if (n.graph_ == this && n.uid_ < num_nodes()) {
      return true;
    }
    // if the Node object has a pointer to another graph or its unique id doesn't
      // currently exist in the graph
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
    // check that the specified index exists in this graph
    assert(i < num_nodes());
    // return the node with the specified index
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
    /** Default public constructor. Constructs an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      // Dereference the graph pointer to access the unordered map containing pairs
        // of an edge uid & a vector of two node uid's. Get the first node uid. 
      size_type node1_uid = graph_->edge_uid_umap.at(uid_)[0];
      // return the node corresponding to node1_uid
      return Node(graph_, node1_uid);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // Dereference the graph pointer to access the unordered map containing pairs
        // of an edge uid & a vector of two node uid's. Get the second node uid. 
      size_type node2_uid = graph_->edge_uid_umap.at(uid_)[1];
      // return the node corresponding to node2_uid
      return Node(graph_, node2_uid);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // check that the 2 edges have pointers to the same graph and the same uid
        // (this check suffices b/c the Graph class's add_edges() functions only
        // allows for one edge uid between every pair of nodes)
      if (graph_ == e.graph_ && uid_ == e.uid_) {
        return true;
      }
      // if the two edges are from different graphs or have different uid's
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
      // check that the 2 edges have pointers to the same graph and that this
        // edge's uid is smaller than the uid of the specified edge
      if (graph_ == e.graph_ && uid_ < e.uid_) {
        return true;
      }
      else {
        return false;
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Pointer back to the Graph container
    Graph* graph_;
    // The edge's unique identification number
    size_type uid_;
    // Private Constructor
    Edge(const Graph* graph, size_type uid) 
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // get the current number of elements in the unordered map where the key is
      // edge uid
    return edge_uid_umap.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // check that the specified edge index exists in this graph
    assert(i < num_edges());
    // return the edge with the specified index
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // check that the two nodes both have pointers to this graph
    assert(a.graph_ == this && b.graph_ == this);
    // check that the two nodes have uids that exist in this graph
    assert(a.uid_ < num_nodes() && b.uid_ < num_nodes());
    // since the order of the nodes doesn't matter, create two pairs with the
      // opposite ordering of node a's uid and node b's uid
    std::pair<size_type, size_type> pair1 {a.uid_, b.uid_};
    std::pair<size_type, size_type> pair2 {b.uid_, a.uid_}; 
    // check if one of the two orderings exists in the map where keys are unique
      // node pairs
    if (node_uid_map.count(pair1) > 0 || node_uid_map.count(pair2) > 0) {
      return true;
    }
    // if neither of the orderings exists, return false
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
    // use the has_edge() method to see if 2 nodes are already connected by an edge
    if (has_edge(a, b)) {
      // if so, define a variable where we will store the existing edge's uid
      size_type edge_uid;
      // tuple containing one of the two possible orderings of the nodes' uid's
      std::pair<size_type, size_type> pair1 {a.uid_, b.uid_};
      // if this first ordering exists as a key in our map, return the edge uid
      if (node_uid_map.count(pair1) > 0) {
        edge_uid = node_uid_map.at(pair1);
      }
      // if the first ordering isn't the key, create a tuple containing the other
        // possible ordering and use it to fetch the edge uid
      else {
        std::pair<size_type, size_type> pair2 {b.uid_, a.uid_};
        edge_uid = node_uid_map.at(pair2);
      }
      // return the Edge corresponding to this edge uid
      return Edge(this, edge_uid);
    }
    // if the two specified nodes aren't already connected by an edge
    else {
      // get the index for a new edge
      size_type idx = num_edges();
      // update the unordered map that has edge uid's as keys with a new edge
      edge_uid_umap.insert(make_pair(idx, std::vector<size_type>{a.uid_, b.uid_}));
      // update the map that has the 2 node's uids as the key with a new edge
      std::pair<size_type, size_type> pair {a.uid_, b.uid_}; 
      node_uid_map.insert(make_pair(pair, idx));
      return Edge(this, idx);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // clear the vector storing the graph's internal nodes
    nodes_vector.clear();
    // clear the 2 maps storing the relationship between edge uid's and node
      // uid's for this graph
    edge_uid_umap.clear();
    node_uid_map.clear();
  }

 private:

  // Internal type for node elements
  struct internal_node {
    // Constructor
    internal_node(Point p, size_type index)
      : point(p), 
        uid(index)
    {}

    Point point; // the node's position
    size_type uid; // the node's unique id
  };

  // a vector containing our internal nodes
  std::vector<internal_node> nodes_vector;

  // an unordered map where the key is an edge's unique id and the value is
    // a vector containing the node 1's unique id and node's 2 unique id
  std::unordered_map<size_type, std::vector<size_type>> edge_uid_umap;

  // an map where the key is a tuple containing node 1's unique id & node 2's
    // unique id, and the corresponding value is the edge uid
  std::map<std::pair<size_type, size_type>, size_type> node_uid_map;
};

#endif // CME212_GRAPH_HPP
