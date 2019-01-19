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
using namespace std;

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
 private:
    // HW0 Addition
    struct NodeElement;
    struct EdgeElement;
    vector<NodeElement> node_vector;
    vector<EdgeElement> edge_vector;

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
  Graph() : node_vector(), edge_vector() {}

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
    Node() {}

    /** Return this node's position. */
    const Point& position() const {
      // HW0 Addition
      return graph->node_vector.at(id).pos;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0 Addition
      return id;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0 Addition
      if (n.id == id && n.graph == graph) {
        return true;
      }
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
      // HW0 Addition
      if (n.id < id) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0 Addition
    Graph* graph;
    size_type id;
    Node(const Graph* graph_, size_type id_) : graph(const_cast<Graph*>(graph_)), id(id_) {}
  };
  //*****END NODE DEF*****

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0 Addition
    return node_vector.size();
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
    // HW0 Addition
    node_vector.push_back(NodeElement(position));
    return Node(this, size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0 Addition
    if (node(n.id) == n) {
      return true;
    }
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0 Addition
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
    Edge() {}

    /** Return a node of this Edge */
    Node node1() const {
      // HW0 Addition
      return Node(graph, graph->edge_vector.at(edge_id).node1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0 Addition
      return Node(graph, graph->edge_vector.at(edge_id).node2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // HW0 Addition
      if ((e.node1() == node1() && e.node2() == node2()) || 
	  (e.node1() == node2() && e.node2() == node1())) {
        return true;
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // HW0 Addition
      if (edge_id < e.edge_id) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0 Addition
    Graph* graph;
    size_type edge_id;
    Edge(const Graph* graph_, size_type edge_id_) : graph(const_cast<Graph*>(graph_)), edge_id(edge_id_) {}
  };
  //*****END EDGE DEF*****

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0 Addition
    return edge_vector.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0 Addition
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0 Addition
    // Iterate through adjacency list for node a to check for node b 
    vector<size_type> edge_nodes = node_vector.at(a.id).node_map;
    if (!edge_nodes.empty() && a.graph == b.graph) {
      vector<size_type>::iterator it;
      it = find(edge_nodes.begin(), edge_nodes.end(), b.id);
      if (it != edge_nodes.end()) {
        return true;
      }
    }
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
    // HW0 Addition
    // Iterate through adjacency list for node a to check for node b 
    vector<size_type>& edge_nodes = node_vector.at(a.id).node_map;
    if (!edge_nodes.empty() && a.graph == b.graph) {
      vector<size_type>::iterator it;
      it = find(edge_nodes.begin(), edge_nodes.end(), b.id);
      if (it != edge_nodes.end()) {
        size_type edge_loc = distance(edge_nodes.begin(), it);
        return edge(node_vector.at(a.id).edge_map.at(edge_loc));
      }
    }
    // Add new edge if not existing
    size_type edge_id = num_edges();
    edge_vector.push_back(EdgeElement(a.id, b.id));
    // Add edge and edge_ids to adjacency lists
    node_vector.at(a.id).node_map.push_back(b.id);
    node_vector.at(b.id).node_map.push_back(a.id);
    node_vector.at(a.id).edge_map.push_back(edge_id);
    node_vector.at(b.id).edge_map.push_back(edge_id);
    return edge(edge_id);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0 Addition
    node_vector.clear(); 
    edge_vector.clear();
  }

 private:
  // HW0 Addition
  struct NodeElement {
    NodeElement(const Point p) : pos(p), node_map(), edge_map() {}
    const Point pos;       			//node position
    vector<size_type> node_map;			//adjacency list: nodeid -> nodeid
    vector<size_type> edge_map;  		//adjacency list: nodeid -> edgeid
  };
  struct EdgeElement {
    EdgeElement(size_type n1, size_type n2) : node1(n1), node2(n2) {}
    size_type node1;				//node1 id
    size_type node2;				//node2 id
  };

};

#endif // CME212_GRAPH_HPP
