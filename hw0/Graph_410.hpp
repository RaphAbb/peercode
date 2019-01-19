#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
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

  std::vector<Point> nodeLocations;
  std::vector<Point*> edgeNode1;
  std::vector<Point*> edgeNode2;

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
      // Return a default initialized node.
    }

    /** Return this node's position. */
    const Point& position() const {
      return this->graphContainer->nodeLocations[this->nodeIndex];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return this->nodeIndex;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (this->nodeIndex == n.nodeIndex);
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
      return (this->nodeIndex < n.nodeIndex);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Allow Edge to also access Node's private information.
    friend class Edge;

    // On 64-bit architecture, size_type and pointers generally use 8 bytes
    // each, so a Node object will take about 16 bytes total.
    size_type nodeIndex;
    const Graph* graphContainer;

    // A legitimate constructor for Nodes that is private. Only the graph
    // class can call this object.
    Node(const Graph* graphPointer, size_type index) :
      nodeIndex(index), graphContainer(graphPointer) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return this->nodeLocations.size();
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
    nodeLocations.push_back(position);
    return Node(this, this->num_nodes()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (this == n.graphContainer);
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
      // Extract the address location of node 1's position.
      Point* positionPointer = graphContainer->edgeNode1[this->edgeIndex];

      // Then, we subtract it off the location of the 0th Point we have
      // stored to extract the index.
      size_type index = positionPointer - &(graphContainer->nodeLocations[0]);

      return Node(this->graphContainer, index);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // Extract the address location of node 1's position.
      Point* positionPointer = graphContainer->edgeNode2[this->edgeIndex];

      // Then, we subtract it off the location of the 0th Point we have
      // stored to extract the index.
      size_type index = positionPointer - &(graphContainer->nodeLocations[0]);

      return Node(this->graphContainer, index);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {

      // If both the graph containers and the edges match, we return true.
      // Otherwise false.
      if (this->graphContainer == e.graphContainer) {
        return (this->edgeIndex == e.edgeIndex);
      }
      return false;

    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // Assuming that this and the provided edge are in the same graph.
      return (this->edgeIndex < e.edgeIndex);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    
    Edge(const Graph* graphPointer, size_type indexOfEdge) :
     graphContainer(graphPointer), edgeIndex(indexOfEdge) {
    }

    const Graph* graphContainer;
    size_type edgeIndex;
    
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edgeNode1.size();
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
    size_type index = findEdge(a,b);
    return (index < this->num_edges());
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

    size_type index = findEdge(a, b);

    if (index == this->num_edges()) {  
      edgeNode1.push_back(&(this->nodeLocations[a.index()]));
      edgeNode2.push_back(&(this->nodeLocations[b.index()]));
    }

    return Edge(this, index);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodeLocations.clear();
    edgeNode1.clear();
    edgeNode2.clear();
  }

 private:
  /** Search for the index of an edge in the graph given the end nodes.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return the index of the corresponding edge if it is in the graph
   *         or the number of edges in the graph otherwise.
   */ 
  size_type findEdge(const Node& a, const Node& b) const {

    const Point* pointA = &(this->nodeLocations[a.index()]);
    const Point* pointB = &(this->nodeLocations[b.index()]);

    for (size_type i = 0; i<this->num_edges(); i++) {
      if (pointA == this->edgeNode1[i] and pointB == this->edgeNode2[i]) {
        return i;
      }
      if (pointA == this->edgeNode2[i] and pointB == this->edgeNode1[i]) {
        return i;
      }

    }

    return this->num_edges();

  }

};

#endif // CME212_GRAPH_HPP
