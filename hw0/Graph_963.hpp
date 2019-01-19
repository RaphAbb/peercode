#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <tuple>
#include <utility>
#include <functional>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

/*Hashing Function for unordered maps of pairs*/
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1,T2> &pair) const {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

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

  // Structs for nodes and edges of the graph
  struct graph_node;
  struct graph_edge;

  // Vectors containing nodes and edges based on their indices
  std::vector<graph_node> nodes;
  std::vector<graph_edge> edges;

  // Maps of uid to the indices for O(1) search 
  std::unordered_map<std::size_t,std::size_t> node_uid_index;
  std::unordered_map<std::size_t,std::size_t> edge_uid_index;

  // Map of nodes to edges
  std::unordered_map<std::pair<std::size_t,std::size_t>,
                std::pair<std::size_t,std::size_t>,pair_hash> node_edge_map;

  // Uids of the next node and edge
  std::size_t next_node_uid_;
  std::size_t next_edge_uid_;

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
    nodes = std::vector<graph_node>();
    edges = std::vector<graph_edge>(); 
    next_node_uid_=0; next_edge_uid_=0;
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
      return fetch().point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return fetch().index;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if ((this->graph_ == n.graph_) && (this->index()==n.index()))
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
      // HW0: YOUR CODE HERE
      if (this->index() < n.index()) return true;
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_;
    size_type uid_;

    /* Fetch the node of the graph with its uid == this->uid_ */
    graph_node& fetch() const {
        auto node_index = this->graph_->node_uid_index[this->uid_];
        return this->graph_->nodes[(int)node_index];
    }

    Node(const Graph* graph, size_type uid) : graph_(const_cast<Graph*>
                                         (graph)),uid_(uid) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes.size();
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
    
    // Initialize a node of the graph with its appropriate member attributes
    graph_node n_node = { num_nodes(), position, next_node_uid_ };

    // Update the map of uid vs index
    node_uid_index[next_node_uid_] = num_nodes();

    // Add the node to the vector of graph nodes
    nodes.push_back(n_node);

    ++next_node_uid_;
    
    return Node(this,next_node_uid_-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (nodes[n.index()].uid == n.uid_) return true;
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
    assert((i>=0)&&(i<num_nodes()));
    return Node(this,nodes[i].uid);
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
      return fetch().node_1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return fetch().node_2;      
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if ((this->node1() == e.node1()) && (this->node2() == e.node2()))
          return true;
      else if ((this->node1() == e.node2()) && (this->node2() == e.node1()))
          return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (this->fetch().index < e.fetch().index) return true;
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type uid_;
    
    // Fetch the edge of the graph with uid == this->uid_
    graph_edge& fetch() const {
        auto edge_index = this->graph_->edge_uid_index[this->uid_];
        return this->graph_->edges[edge_index];
    }

    Edge(const Graph* graph, size_type uid) : graph_(const_cast<Graph*>
                                         (graph)),uid_(uid) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert((i>=0)&&(i<num_edges()));
    return Edge(this,edges[i].uid);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    assert(has_node(a) && has_node(b));
    
    std::pair<std::size_t,std::size_t> node_pair,node_pair_swap;
    node_pair = std::make_pair(a.uid_,b.uid_);
    node_pair_swap = std::make_pair(b.uid_,a.uid_);

    if ((node_edge_map.find(node_pair) == node_edge_map.end()) &&
        (node_edge_map.find(node_pair_swap) == node_edge_map.end()))
        return false;

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
  Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    assert(not(a==b));
    assert(has_node(a) && has_node(b));

    std::pair<std::size_t,std::size_t> node_pair,edge_details;

    // If the graph already has the edge, return it
    if (has_edge(a,b)) {
        node_pair = std::make_pair(a.uid_,b.uid_);
        std::pair<std::size_t,std::size_t> node_pair_swap(b.uid_,a.uid_);

        if (node_edge_map.find(node_pair) != node_edge_map.end()) {
            auto edge_uid = node_edge_map[node_pair].first;
            return Edge(this,edge_uid);
        }
        auto edge_uid = node_edge_map[node_pair_swap].first;
        return Edge(this,edge_uid);
    }    
    // If the edge doesn't exist, create a new one           
    graph_edge n_edge = { num_edges(), a, b, next_edge_uid_ };
    
    // Update edge uid vs index map
    edge_uid_index[next_edge_uid_] = num_edges();

    // Update node vs edge map
    node_pair = std::make_pair(a.uid_,b.uid_);
    edge_details = std::make_pair(next_edge_uid_,num_edges());
    node_edge_map[node_pair] = edge_details;
    
    // Add to the edges of the graph vector
    edges.push_back(n_edge);

    ++next_edge_uid_;

    return Edge(this,next_edge_uid_-1);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();
    edges.clear();
    node_uid_index.clear();
    edge_uid_index.clear();
    return;
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct graph_node {
      std::size_t index;
      Point point;
      std::size_t uid;
  };

  struct graph_edge {
      std::size_t index;
      Node node_1;
      Node node_2;
      std::size_t uid;
  };

};

#endif // CME212_GRAPH_HPP
