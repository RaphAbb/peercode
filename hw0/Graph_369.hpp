#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
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
class Graph {

  struct internal_edge;

 private:

  // list of node data containers
  std::vector<Point> point_list;
  // list of edge data containers
  std::vector<internal_edge> edge_list;

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
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      return gph->point_list[this->idx];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return size_type(this->idx);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (n.index() == this->idx and n.gph == this->gph)
        return true;
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
      if (n.index() < this->idx)
        return true;
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // pointer back to node's graph
    const Graph* gph;
    // store node's index
    size_type idx;

    // private node constructor that sets node data
    Node(const Graph* graph, int index) {
      this->gph = graph;
      this->idx = index;
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return point_list.size();
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
    Node new_node = Node(this, point_list.size());
    point_list.push_back(position);

    // add empty map to graph adjacency list
    std::map<size_type, size_type> m;
    adjacency_list.push_back(m);

    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.index() < this->point_list.size())
      return true;
    return false;
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(this->gph, fetch().node1_idx);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(this->gph, fetch().node2_idx);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // get edge nodes to avoid looking up multiple times
      int e_node1 = e.node1().index();
      int e_node2 = e.node2().index();
      int t_node1 = fetch().node1_idx;
      int t_node2 = fetch().node2_idx;

      if ((e_node1==t_node1 and e_node2==t_node2) 
        or (e_node1==t_node2 and e_node2==t_node1))
        return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (this->idx < e.idx)
        return true;
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph; 
    // store pointer to graph the edge belongs to 
    const Graph* gph;
    // store the index of the edge
    size_type idx;

    // private constructor that sets edge values
    Edge(const Graph* graph, int index){
      this->gph = graph;
      this->idx = index;
    }

    // helper function to get internal edge object for this edge
    const internal_edge& fetch() const {
      return gph->edge_list[idx];
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_count;
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
    if (adjacency_list[a.index()].count(b.index())!=0)
      return true;
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
    // check if edge is already in graph
    if (this->has_edge(a,b)){
      size_type its_index = adjacency_list[a.index()][b.index()];
      return Edge(this,its_index);
    }

    // create new internal edge object
    internal_edge new_edge = internal_edge(a,b, edge_count);
    this->edge_list.push_back(new_edge);
    edge_count++;

    // add new edge to adjacency list
    adjacency_list[a.index()][b.index()] = edge_count-1;
    adjacency_list[b.index()][a.index()] = edge_count-1;

    return Edge(this, edge_count-1);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // clears edge list, point list, adjacency list, resets edge count
    this->edge_list.clear();
    this->point_list.clear();
    this->adjacency_list.clear();
    edge_count = 0;
  }

 private:

  // adjacency map to store graph edge information
  std::vector<std::map<size_type, size_type>> adjacency_list;
  // keep a count of current number of edges
  size_type edge_count = 0;

  // internal edge object to store information about a point
  struct internal_edge {
    // store node indices and edge index
    size_type node1_idx;
    size_type node2_idx;
    size_type idx;
    // constructor to set internal edge values
    internal_edge(const Node& a, const Node& b, size_type index){
      node1_idx = a.index();
      node2_idx = b.index();
      idx = index;
    }
  };

};

#endif // CME212_GRAPH_HPP
