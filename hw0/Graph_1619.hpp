#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <set>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
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

 private:
  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  /* Explanation of below data structures
     - node_to_position_map: (key,value) = (node index, Point object)
     - edge_to_node_map: (key,value) = (edge index, vector of node indexes)
     - node_to_edge_map: (key,value) = (set of node index, edge index)
  */
  std::map<size_type,Point> node_to_position_map;
  std::map<size_type,std::vector<size_type>> edge_to_node_map;
  std::map<std::set<size_type>,size_type> node_to_edge_map;
  size_type number_nodes; // number of nodes
  size_type number_edges; // number of edges

 public:

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    number_nodes = 0;
    number_edges = 0;
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
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return owner_graph->node_to_position_map.at(node_id);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
     return node_id;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
     if ((owner_graph == n.owner_graph) && (node_id == n.node_id)) {return true;}
     else {return false;}
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
      // I'm assuming we'll only compare nodes that are in same graph
     if (node_id < n.node_id) {return true;}
     else {return false;}
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    const Graph *owner_graph; // pointer to a constant object of type Graph
    size_type node_id;  // index of node

    // Constructor that will be used by Graph::add_node()
    Node(const Graph *graph, size_type id)
        : owner_graph(graph), node_id(id) { }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return number_nodes;
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

    /*
    This code could be made more robust by first checking that there doesn't
    already exist a node with position @param[in]. This could not be implemented
    in O(1) time if we searched through class member node_to_position_map. To do so in
    O(1) time would require creating a new data structure that maps
    Point objects to indexes of nodes to which they correspond. It seems
    like we can assume all Point objects are unique... so that's what
    I'm doing here.
    */
    // Note: index of new node is current value of number_nodes
    node_to_position_map.insert({number_nodes, position});
    number_nodes++; // just made new node, so increment number of nodes
    return node(number_nodes-1);  // -1 to get back to node just created
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    /* For graph to have node n, two things must be true:
        (1) n is associated with current graph instance
        (2) graph instance must contain a node with same index as n */
    if ((n.owner_graph == this) && (n.node_id < number_nodes)) {return true;}
    else {return false;}
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    // Ensure function can only return nodes that already exist in graph
    assert( (0 <= i) && (i < number_nodes) );
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
     return owner_graph->node(owner_graph->edge_to_node_map.at(edge_id)[0]);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
     return owner_graph->node(owner_graph->edge_to_node_map.at(edge_id)[1]);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
     if ( (owner_graph == e.owner_graph) && (edge_id == e.edge_id) ) {return true;}
     else {return false;}
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
     if (edge_id < e.edge_id) {return true;}
     else {return false;}
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    const Graph *owner_graph; // pointer to a constant object of type Graph
    size_type edge_id;  // index of edge

    // Constructor that will be used by Graph::add_edge()
    Edge(const Graph *graph, size_type id)
        : owner_graph(graph), edge_id(id) { }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return number_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    // Ensure function can only return edges that already exist in graph
    assert( (0 <= i) && (i < number_edges) );
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
     std::set<size_type> node_set = {a.index(),b.index()};
     // Check if node_set is in node_to_edge_map
     if (node_to_edge_map.find(node_set) == node_to_edge_map.end()) {return false;}
     else {return true;}
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
     size_type edge_1 = a.index();
     size_type edge_2 = b.index();
     std::set<size_type> node_set = {edge_1, edge_2};
     std::vector<size_type> node_vec = {edge_1, edge_2};
    if (has_edge(a,b) == false)
    {
     // Add new edge
     // Note: index of new edge is current value of number_edges
     node_to_edge_map.insert({node_set, number_edges});
     edge_to_node_map.insert({number_edges, node_vec});
     number_edges++;  // just made new edge, so increment number of edges
     return edge(number_edges-1);  // -1 to get back to edge just created
    }
    else
    {
      // Note: index of existing edge is value of node_to_edge_map
      return edge(node_to_edge_map.at(node_set));
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    // Empty all maps
    node_to_position_map.clear();
    node_to_edge_map.clear();
    edge_to_node_map.clear();
  }

};

#endif // CME212_GRAPH_HPP
