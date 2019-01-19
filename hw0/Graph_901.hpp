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

// HW0: Add unordered_map
#include <unordered_map>

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
  
  // Index start value
  const static unsigned FIRST_INDEX {0};
  // Invalid index value (assuming largest unsigned int is not needed)
  const static unsigned INVALID_INDEX {unsigned(-1)};

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
  Graph() : next_nid_(FIRST_INDEX), first_valid_nid_(FIRST_INDEX),
            next_eid_(FIRST_INDEX), first_valid_eid_(FIRST_INDEX) {
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
    Node() : graph_(nullptr), nid_(INVALID_INDEX) {
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      
      // Node and Point vectors use the same indexing; index() checks if the
      // node is valid
      return graph_->points_.at(this->index());
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      
      this->assert_valid();
      
      return graph_->node_indices_.at(nid_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      
      // index() checks if node is valid
      return (graph_->has_node(n) ? this->index() == n.index() : false);
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
      
      // As per piazza post #35, assume there is no comparison between graphs
      // (i.e. assume this and n belong to the same graph)
      return (*this == n ? false : this->index() < n.index());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    // Pointer to owning Graph
    Graph* graph_;
    // Unique node ID
    size_type nid_;
    
    /** Private constructor for valid nodes */
    Node(const Graph* graph, size_type nid)
        : graph_(const_cast<Graph*>(graph)), nid_(nid) {
    }
    
    /** Helper method for asserting that the node is valid */
    void assert_valid() const {
      // Check that the node belongs to a graph
      assert(graph_ != nullptr);
      // Check to see if the node was removed from the graph, e.g. via clear()
      assert(graph_->node_indices_.at(nid_) != INVALID_INDEX);
    }
    
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return size_type(nodes_.size());
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
    
    points_.push_back(position);
    node_indices_.push_back(size());  // Index is based on order in vector
    nodes_.push_back(Node(this,next_nid_));
    ++next_nid_;
    
    return nodes_.back();
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    
    return (this == n.graph_ && node_indices_.at(n.nid_) != INVALID_INDEX);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    
    assert(FIRST_INDEX <= i && i < size());
    return nodes_[i];
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
    Edge() : graph_(nullptr), node1_(nullptr), node2_(nullptr), 
             eid_(INVALID_INDEX) {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      
      // Access proper edge (index checks for validity) and dereference the
      // pointer to node1
      return *(graph_->edge(this->index()).node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return *(graph_->edge(this->index()).node2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // node1() and node2() functions check if edge is valid
      return  ( ( node1() == e.node1() && node2() == e.node2() ) || 
                ( node1() == e.node2() && node2() == e.node1() ) );
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // Assume this and e belong to the same graph
      return (*this == e ? false : this->index() < e.index());
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    
    // Pointer to owning graph
    Graph* graph_;
    // Pointers to nodes that define the edge
    Node* node1_;
    Node* node2_;
    // Unique edge ID
    size_type eid_;
    
    /** Private constructor for valid edges */
    Edge(const Graph* graph, const Node* node1, const Node* node2, 
        size_type eid) : graph_(const_cast<Graph*>(graph)), 
        node1_(const_cast<Node*>(node1)), node2_(const_cast<Node*>(node2)),
        eid_(eid) {
    }
    
    /** Helper method that checks if the edge is valid and returns the index */
    size_type index() const {
      this->assert_valid();
            
      return graph_->edge_indices_[eid_];
    }
    
    /** Helper method that checks if the edge is valid */
    void assert_valid() const {
      assert(graph_ != nullptr && node1_ != nullptr && node2_ != nullptr);
      assert(graph_->edge_indices_.at(eid_) != INVALID_INDEX);
      // Check if either node is invalid
      assert(node1_->index() != INVALID_INDEX && node2_->index() != INVALID_INDEX);
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return size_type(edges_.size());
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    
    assert(FIRST_INDEX <= i && i < num_edges());
    return edges_[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    
    // Assert precondition
    assert(has_node(a));
    assert(has_node(b));
    
    // eid_ of edge connecting a and b, if one exists
    size_type eid = get_edge_eid(a, b);
    // Check for case that a removed edge was not properly removed from the map
    if (eid != INVALID_INDEX) {
      assert(edge_indices_.at(eid) != INVALID_INDEX);
      return true;
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
    
    // has_edge() will assert that a, b exist in the graph; assert that they are
    // distinct
    assert(!(a == b));
    bool edge_exists = has_edge(a, b);
    if ( !edge_exists ) {
      // Edge does not exist and will be added
      edge_indices_.push_back(num_edges());
      
      // Add map elements (one for each node to speed up the search
      nodes_to_edge_.emplace(a.nid_, std::vector<size_type>({b.nid_, next_eid_}));
      nodes_to_edge_.emplace(b.nid_, std::vector<size_type>({a.nid_, next_eid_}));
      
      // Add edge
      Node* node1 = get_node_pointer(a);
      Node* node2 = get_node_pointer(b);
      edges_.push_back( Edge(this, node1, node2, next_eid_) );
      ++next_eid_;
      
      return edges_.back();
      
    }
    size_type eid = get_edge_eid(a, b); // We know this edge exists in the graph
    return edges_[edge_indices_[eid]];
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    points_.clear();
    edges_.clear();
    nodes_to_edge_.clear();
    
    size_type i;
    // Starting from first_valid_nid_ avoids unnecessary iterating
    for (i = first_valid_nid_; i < next_nid_; ++i) {
      node_indices_[i] = INVALID_INDEX;
    }
    first_valid_nid_ = next_nid_;
    
    for (i = first_valid_eid_; i < next_eid_; ++i) {
      edge_indices_[i] = INVALID_INDEX;
    }
    first_valid_eid_ = next_eid_;
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  
  // Contains graph nodes
  std::vector<Node>  nodes_;
  // Contains node positions
  std::vector<Point> points_;
  // Vector that maps node IDs to node index
  std::vector<size_type> node_indices_;
  // Next node ID for assignment
  size_type  next_nid_;
  // First index for node_indices_ that corresonds to a valid node
  size_type first_valid_nid_;
  // Contains grah edges
  std::vector<Edge> edges_;
  // Vector that maps edge IDs to edge index
  std::vector<size_type> edge_indices_;
  // Next edge ID for assignment
  size_type next_eid_;
  // First index for edge_indices_ that corresponds to a valid node
  size_type first_valid_eid_;
  // Maps a node to another node and the edge joining them; key is the nid_ of
  // the first node, value is a vector with nid_ of the second node and the eid_
  // Will contain two entries per edge to be order invariant when searching
  std::unordered_multimap< size_type, std::vector<size_type> > nodes_to_edge_;
  
  /** Helper method that returns an iterator to the nodes_to_edge_ map if an
   *  edge exists between them; otherwise the end() iterator is returned
   *
   *  Checking if nodes @a a and @a b are valid is done externally
   */
  size_type get_edge_eid(const Node& a, const Node& b) const {
    auto node1_edges = nodes_to_edge_.equal_range(a.nid_);
    
    for (auto it = node1_edges.first; it != node1_edges.second; ++it) {
      if (it->second[FIRST_INDEX] == b.nid_)
        return it->second.back();
    }
    return INVALID_INDEX;
  }
  
  /** Helper method that returns a pointer to the graph's copy of the node */
  Node* get_node_pointer(const Node& a) const {
    return const_cast<Node*>(&nodes_[a.index()]);
  }

};

#endif // CME212_GRAPH_HPP
