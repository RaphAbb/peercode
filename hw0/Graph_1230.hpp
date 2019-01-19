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
#include <set>

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
  private:
    std::vector<Point> internal_nodes; // vector storing all of the nodes' positions
    std::vector<std::set<unsigned int>> internal_edges; // vector storing all of edges information
    unsigned int graph_size;
    
    
    
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
  using size_type = unsigned int;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
      // HW0: YOUR CODE HERE
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
        // assert(this->node_idx);
        std::vector<Point> vec = *this->node_ptr;
        // begin debug
        std::cout << std::endl;
        return vec[this->node_idx];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
        // HW0: YOUR CODE HERE
        //assert(this->node_idx);
        return this->node_idx;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
        // HW0: YOUR CODE HERE
        if (this->node_idx == n.node_idx && this->parent == n.parent) {
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
        // HW0: YOUR CODE HERE
        // assuming the nodes being compared belong to the same graph
        if (this->node_idx < n.node_idx) {
            return true;
        }
      return false;
    }

   private:
      // Allow Graph to access Node's private member data and functions.
      friend class Graph;
      // HW0: YOUR CODE HERE
      size_type node_idx;
      const Graph* parent;
      const std::vector<Point>* node_ptr;
      
      // set up node
      void Node_setup(size_type node_idx, const Graph* parent, const std::vector<Point>* node_ptr){
          this->node_idx = node_idx;
          this->parent = parent;
          this->node_ptr = node_ptr;
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
      if (graph_size) {
          return graph_size;
      }
      return 0;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
      return internal_nodes.size();
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
      // add node if node doesn't exist
      // push back the node
      this->internal_nodes.push_back(position);
      
      // set up the Node
      Node to_add;
      to_add.Node_setup(internal_nodes.size(), this, &(this->internal_nodes));
      //node_vec.push_back(to_add);
      return to_add;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
      // HW0: YOUR CODE HERE
      // Recover the position from the node
      Point pos = n.position();
      
      // retrieve the idx num of the node
      size_type idx = n.index();
      
      // see if the node belongs to the graph
      if (internal_nodes[idx] == pos) {
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
      // HW0: YOUR CODE HERE
      size_type idx_range = num_nodes();
      if (i < idx_range) {
          Node tmp_node;
          const Graph* this_graph = this;
          tmp_node.Node_setup(i, this_graph, &(this->internal_nodes));
          const Node to_return = tmp_node;
          return to_return;
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
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
        // HW0: YOUR CODE HERE
        return this->node_1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
        // HW0: YOUR CODE HERE
        return this->node_2;

    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
        // 1-1 2-2 correspondence
        if (this->node_1 == e.node_1 && this->node_2 == e.node_2) {
            return true;
        }
        // 1-2 2-1 correspondence
        if (this->node_2 == e.node_1 && this->node_1 == e.node_2) {
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
        if (this->edge_idx < e.edge_idx) {
            return true;
        }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
      friend class Graph;
      Node node_1;
      Node node_2;
      size_type edge_idx;
      
      Edge edge_setup(Node node1, Node node2, size_type edge_idx) {
          this->node_1 = node1;
          this->node_2 = node2;
          this->edge_idx = edge_idx;
          return *this;
      }

   
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return this->internal_edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
      // HW0: YOUR CODE HERE
      // set up and return
      Edge to_return;
      Node node1 = this->node(*(this->internal_edges[i].begin()));
      Node node2 = this->node(*(this->internal_edges[i].end()));
      to_return.edge_setup(node1, node2, i);
    return to_return;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
      // HW0: YOUR CODE HERE
      std::set<unsigned int> to_check;
      to_check = {a.index(), b.index()};
      for (unsigned int i = 0; i < this->internal_edges.size(); i++) {
          if (internal_edges[i] == to_check) {
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
    // HW0: YOUR CODE HERE
      Edge to_add;
      this->internal_edges.push_back({a.node_idx, b.node_idx});
      return to_add.edge_setup(a, b, this->num_edges());
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
      // HW0: YOUR CODE HERE
      this->internal_edges = {};
      this->internal_nodes = {};
      this->graph_size = 0;
  }

 private:


};

#endif // CME212_GRAPH_HPP
