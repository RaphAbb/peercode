#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
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
  //
  // Internal type declarations using public aliased size types
  //
  /** Orders a pair of node indices in lexicographic ordering. */
  static void orderEdge(size_type& nodeIndex1, size_type& nodeIndex2){
    if (nodeIndex1 < nodeIndex2){
      std::swap(nodeIndex1, nodeIndex2);
    }
  }
  /** Internal node representation (16-byte) */
  struct NodeInternal{
    const Point point;
    const size_type id;
    NodeInternal(Point p, size_type i) : point{p}, id{i} {}
  };
  /** Internal edge representation (16-byte) */
  struct EdgeInternal {
    size_type node_id_first;
    size_type node_id_second;
    /** Constructs internal node using first and second indices. */
    EdgeInternal(size_type i_first,
                 size_type i_second) : node_id_first{i_first},
                                       node_id_second{i_second} {
      orderEdge(node_id_first, node_id_second);
    }
  };

  /** Private collection of nodes (containing ground truth). */
  std::vector<NodeInternal> nodesInternal;
  /** Private collection of edges (containing ground truth). */
  std::vector<EdgeInternal> edgesInternal;
  /** Redundant fast edge-check map mapping ordered pairs to edge index. */
  std::map<std::pair<size_type, size_type>, size_type> adjacencyMap;
  /** ID shifting scheme to keep proxies to re-built graphs invalid. **/
  size_type node_id_minimum{0};
  size_type edge_id_minimum{0};

 public:

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() { }

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
    Node() { }

    /** Return this node's position. */
    const Point& position() const {
      return parent_graph->pointAt(proxy_id);
    }

    /** Return this node's index, a number in the range [0, graph_size).
     * Index checking is delegated to the parent graph.
     */
    size_type index() const {
      return parent_graph->indexAt(proxy_id);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return parent_graph->indexAt(n.proxy_id) == 
             parent_graph->indexAt(proxy_id);
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
      return parent_graph->indexAt(n.proxy_id) < 
             parent_graph->indexAt(proxy_id);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // ID known to the proxy (guaranteed == ID in Graph)
    size_type proxy_id;
    // Pointer to graph (allows reassignment for invalid nodes but not mutate)
    const Graph* parent_graph;

    /** External Node constructor used by Graph class. */
    Node(size_type new_id, const Graph* new_graph): proxy_id{new_id},
                                              parent_graph{new_graph}{ }
    /** Reassigns properties of proxy node */
    void set_new(size_type new_id, Graph* new_graph){
      proxy_id = new_id;
      parent_graph = new_graph;
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodesInternal.size();
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
    nodesInternal.emplace_back(position, size());
    return Node(size()-1, this);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // Shifts indices by the minimum index to prevent expired proxy access
    return has_node_index(n.proxy_id);
  }

  /** Return the node with index @a i. If none, returns invalid node.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if (has_node_index(i)){
      return Node(i, this);
    }
    return Node();        // Invalid node
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
    Edge() { }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(parent_graph->nodeIndexOfEdge(edge_proxy_id, 0),
                  parent_graph);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(parent_graph->nodeIndexOfEdge(edge_proxy_id, 1),
                  parent_graph);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // Unpack edge node indices
      size_type e1_n1, e1_n2, e2_n1, e2_n2;
      e1_n1 = parent_graph->nodeIndexOfEdge(edge_proxy_id, 0);
      e1_n2 = parent_graph->nodeIndexOfEdge(edge_proxy_id, 1);
      e2_n1 = parent_graph->nodeIndexOfEdge(e.edge_proxy_id, 0);
      e2_n2 = parent_graph->nodeIndexOfEdge(e.edge_proxy_id, 1);
      // Return equality, exhausting combinatorics
      return (e1_n1 == e2_n1 && e1_n2 == e2_n2) ||
        (e1_n1 == e2_n2 && e1_n2 == e2_n1);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // Unpack edge node indices
      size_type e1_n1, e1_n2, e2_n1, e2_n2;
      e1_n1 = parent_graph->nodeIndexOfEdge(edge_proxy_id, 0);
      e1_n2 = parent_graph->nodeIndexOfEdge(edge_proxy_id, 1);
      e2_n1 = parent_graph->nodeIndexOfEdge(e.edge_proxy_id, 0);
      e2_n2 = parent_graph->nodeIndexOfEdge(e.edge_proxy_id, 1);
      // Return lexicographic order
      return e1_n1 != e2_n1 ?
        e1_n1 < e2_n1 : (e1_n2 != e2_n2 ? e1_n2 < e2_n2 : false);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    
    // ID known to the proxy (guaranteed == ID in Graph)
    size_type edge_proxy_id;
    // Pointer to graph (allows reassignment for invalid nodes but not mutate)
    const Graph* parent_graph;

    /** Reassigns properties of proxy edge */
    void set_new(size_type new_id, const Graph* new_graph){
      edge_proxy_id = new_id;
      parent_graph = new_graph;
    }

    /** Parametrized constructor exposed to Graph class */
    Edge(size_type edge_id,
         const Graph* new_graph): edge_proxy_id{edge_id},
                                  parent_graph{new_graph}{
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    return edgesInternal.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges() + num_edges_deleted
   *
   * Complexity: O(1) due to vector backing
   */
  Edge edge(size_type i) const {
    if (has_edge_index(i)){
      return Edge(i, this);
    }
    return Edge();
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid dinstinct nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: O(log(num_edges())) due to map
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Order the edge (a,b) so that a >= b
    size_type nodeIndex1 = indexAt(a.proxy_id);
    size_type nodeIndex2 = indexAt(b.proxy_id);
    orderEdge(nodeIndex1, nodeIndex2);
    // Check membership of target edge index in map
    auto target = std::pair<size_type, size_type>(nodeIndex1, nodeIndex2);
    auto it = adjacencyMap.find(target);
    if (it != adjacencyMap.end()){
      return true;
    }
    return false;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b or
   *   vice versa
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: O(log(num_edges())) due to map
   */
  Edge add_edge(const Node& a, const Node& b) {
    // Order the edge (a,b) so that a >= b
    size_type nodeIndex1 = indexAt(a.proxy_id);
    size_type nodeIndex2 = indexAt(b.proxy_id);
    orderEdge(nodeIndex1, nodeIndex2);
    // Find target edge in map
    auto target = std::pair<size_type, size_type>(nodeIndex1, nodeIndex2);
    auto it = adjacencyMap.find(target);
    if (it != adjacencyMap.end()){ // Edge found
      return Edge(it->second, this);
    }
    // Create new edge: update internal representation
    edgesInternal.emplace_back(a.proxy_id, b.proxy_id);
    adjacencyMap.insert(std::pair<std::pair<size_type, size_type>, size_type>(
      std::make_pair(a.proxy_id, b.proxy_id), edgesInternal.size()));
    return Edge(edgesInternal.size()-1, this);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // Shift indices to keep old references invalid
    node_id_minimum += nodesInternal.size();
    edge_id_minimum += edgesInternal.size();
    nodesInternal.clear();
    edgesInternal.clear();
    adjacencyMap.clear();
  }

 private:
  //
  // Access helpers
  //

  /** Checks if the graph has node @a i after shifting.
   * Shifting ensures old proxies stay invalid after clearing the
   * graph and then adding new nodes and edges. */
  bool has_node_index(size_type i) const {
    return i - node_id_minimum < size();
  }

  /** Checks if the graph has edge @a i after shifting.
   * Shifting ensures old proxies stay invalid after clearing the
   * graph and then adding new nodes and edges. */
  bool has_edge_index(size_type i) const {
    return i - edge_id_minimum < edgesInternal.size();
  }

  /** Returns point (spatial) at @a proxy_id.
   * @pre @a proxy_id may or may not be valid
   * @post point for node id is returned, or program is aborted if invalid
  */
  const Point& pointAt(size_type proxy_id) const{
    assert(has_node_index(proxy_id));
    return nodesInternal.at(proxy_id).point;
  }

  /** Returns point (spatial) at @a proxy_id.
   * @pre @a proxy_id may or may not be valid
   * @post index for node id is returned, or program is aborted if invalid
  */
  size_type indexAt(size_type proxy_id) const{
    assert(has_node_index(proxy_id));
    return nodesInternal.at(proxy_id).id;
  }

  /** Returns the index of one of the nodes of specified edge.
   * @param[in] edge_proxy_index index of the edge reckoned by the proxy
   * @param[in] nodeNumber number (0 or 1) of node attached to edge
  */
  size_type nodeIndexOfEdge(size_type edge_proxy_index,
                            size_type nodeNumber) const{
    assert(has_edge_index(edge_proxy_index));
    assert(nodeNumber == 0 || nodeNumber == 1);
    if (nodeNumber == 0)
      return edgesInternal.at(edge_proxy_index).node_id_first;
    else
      return edgesInternal.at(edge_proxy_index).node_id_second;
  };
};

#endif // CME212_GRAPH_HPP
