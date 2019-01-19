#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <set> 
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
  struct node_element;  
  struct edge_element;

 public:
  class Node;
  class Edge;
 
  using graph_type = Graph;
  using node_type = Node;
  using edge_type = Edge;
  using size_type = unsigned;

  /** Constructor: Empty Graph. */
  Graph(): curr_id_(0) {
    // Set current id counter to zero
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
    Node(): graph_(NULL) {
      // Invalid nodes have a null pointer for graph_
    }

    /** Return this node's position. */
    const Point& position() const {
      // Call a private method, might be helpful later for encapsulation?
      return getPoint();
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return getIndex();
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      assert(graph_ != NULL);
      if (graph_ == n.graph_ && node_id_ == n.node_id_ ) {
        return true; 
      } else {
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
      if (node_id_ < n.node_id_) { 
        return true;
      } else { 
        return false;
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    const graph_type* graph_; // Pointer to graph this node is part of
    size_type node_id_;       // ID number of node (NOT necessarily the index)

    /** Private Constructor for Node - should only be called by Graph */
    Node(const graph_type* graph, size_type nid) 
      : graph_(graph), node_id_(nid) { 
    } 

    /** Helper Method - Reads Graph object to get the most current Point() */ 
    const Point& getPoint() const { 
      assert(graph_ != NULL);
      // could use binary search instead of brute force
      // node_id_ is guaranteed to be increasing
      for (size_type i = 0; i < graph_->size(); i++) {
        if (graph_->nodeList_[i].id == node_id_) {
          return graph_->nodeList_[i].p;
        }
      }
      assert(false); // If the loop exits without return, then node is invalid
    } 

    /** Helper Method - Reads Graph object to get the most current index */
    size_type getIndex() const { 
      assert(graph_ != NULL);
      // could use binary search instead of brute force, see above
      for (size_type i = 0; i < graph_->size(); i++) {
        if (graph_->nodeList_[i].id == node_id_) {
          return i;
        }
      }
      assert(false);
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodeList_.size();
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
    node_element newNode; 
    newNode.p = position; // This calls the copy constructor for Point()
    newNode.id = curr_id_++; 
    nodeList_.push_back( newNode );
 
    return Node(this, newNode.id);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.graph_ == this) {
      return true;
    } else { 
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
    assert( i >= 0 && i < num_nodes());
    return Node(this, nodeList_[i].id);
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
    Edge(): graph_(NULL) {
      // Set graph_ to null for invalid edges
    }

    /** Return a node of this Edge */
    Node node1() const {
      assert(graph_ != NULL); 
      if (flipped_) { 
        return Node(graph_, id2_);
      } else { 
        return Node(graph_, id1_); 
      } 
    }

    /** Return the other node of this Edge */
    Node node2() const {
      assert(graph_ != NULL);
      if (flipped_) { 
        return Node(graph_, id1_);
      } else { 
        return Node(graph_, id2_);  
      }
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      assert(graph_ != NULL); 
      // Private constructor enforces id1 < id2, so don't need to track order
      if (graph_ == e.graph_ && id1_ == e.id1_ && id2_ == e.id2_) {
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
      assert(graph_ != NULL);
      return std::tie(id1_, id2_) < std::tie(e.id1_, e.id2_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    const graph_type* graph_; // Pointer to graph this node is part of
    size_type id1_;           // ID number of node 1 
    size_type id2_;           // ID number of node 2
    bool flipped_;            // True if a.id < b.id, False otherwise
    
    /** Private Edge Constructor, for use by Graph */
    Edge(const graph_type* graph, const Node& a, const Node& b) 
      : graph_(graph), id1_(b.node_id_), id2_(a.node_id_), flipped_(true) {
      if (a.node_id_ < b.node_id_) {
        id1_ = a.node_id_;
        id2_ = b.node_id_;
        flipped_ = false; 
      } 
    }
  
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    assert(edgeList_.size() == edgeSet_.size()); 
    return edgeList_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return edgeList_[i]; 
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert( has_node(a) && has_node(b) ); 
    return edgeSet_.count( Edge(this,a,b) );
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
    // Return current edge if it already exists 
    if (has_edge(a,b)) {
      return Edge(this,a,b); 
    } else { 
      Edge new_Edge = Edge(this,a,b); 
      edgeList_.push_back(new_Edge); 
      edgeSet_.insert(new_Edge);

      return new_Edge;  
    }  
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodeList_.clear();
    edgeList_.clear(); 
    edgeSet_.clear(); 
  }

 private:
  // As with proxy_example.cpp, internal type for a nodeList element
  struct node_element{ 
    size_type id;
    Point p;
  };

  // Private Members
  size_type curr_id_;                  // The current node ID (not index)
  std::vector<node_element> nodeList_; // Contains (Point, node_id) pairs
  std::vector<Edge> edgeList_;         // std::vector of edges
  std::set<Edge> edgeSet_;             // std::set of edges

};

#endif // CME212_GRAPH_HPP
