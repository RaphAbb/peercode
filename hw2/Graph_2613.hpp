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
//HW2
template <typename V, typename E = double>
class Graph {
  //HW2
  typedef V node_value_type;
  typedef E edge_value_type;
 public:
  
  // PUBLIC TYPE DEFINITIONS

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

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  using node_iterator = NodeIterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  using incident_iterator = IncidentIterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;
 
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.

  struct internal_node {
    Point point;
    size_type uid;   
    node_value_type value;
    // adjacency list contains uid of edges
    std::vector<size_type> adj_list;
  };

  std::vector<internal_node> nodes_;
  size_type node_size_;
  size_type next_node_uid_;

  struct internal_edge {
    size_type node1_uid;
    size_type node2_uid;
    size_type uid;
    edge_value_type value;
  };

  std::vector<internal_edge> edges_;
  size_type edge_size_;
  size_type next_edge_uid_;

 public:

  // CONSTRUCTORS AND DESTRUCTOR

  /** Construct an empty graph. */
  
  Graph()
    // HW0: YOUR CODE HERE
      : nodes_(), node_size_(0), next_node_uid_(0), edges_(), edge_size_(0), next_edge_uid_(0) {
  }

  /** Default destructor */
  ~Graph() = default;

  // NODES

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node : private totally_ordered<Node> {
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
      return graph_->nodes_[uid_].point;
    }

    Point& position() {
      // HW2
      return graph_->nodes_[uid_].point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return uid_; 
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    // Return a reference to this node's value
    node_value_type& value() {return graph_->nodes_[uid_].value;}
    
    // Return a constant reference to this node's constant value
    const node_value_type& value() const {return graph_->nodes_[uid_].value;}
    
    // Return degree of node
    size_type degree() const {return graph_->nodes_[uid_].adj_list.size();}
   
    // Iterator begin
    incident_iterator edge_begin() const {return IncidentIterator(graph_, 0, uid_);} 
   
    // Iterator end 
    incident_iterator edge_end() const {return IncidentIterator(graph_, this->degree(), uid_);} 

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (graph_ == n.graph_) and (index() == n.index());
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
      const Point p1 = position();
      const Point p2 = n.position();
      if (p1[0] < p2[0]) {return true;}
      else if (p1[0] > p2[0]) {return false;}
      else 
      {
        if (p1[1] < p2[1]) {return true;}
        else if (p1[1] > p2[1]) {return false;}
        else
          {
            if (p1[2] < p2[2]) {return true;}
          }
      }
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

    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return node_size_;
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
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    // HW1
    internal_node new_node;
    new_node.point = position;
    new_node.uid = next_node_uid_;
    new_node.value = val;
    nodes_.push_back(new_node);
    ++node_size_;
    ++next_node_uid_;
    return Node(this, next_node_uid_-1); 
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (this == n.graph_);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < size());
    return Node(this, i);
  }

  // EDGES

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    // Return a reference to this node's value
    edge_value_type& value() {return graph_->edges_[uid_].value;}
    
    // Return a constant reference to this node's constant value
    const edge_value_type& value() const {return graph_->edges_[uid_].value;}
    
    // HW2
    double length() const {
      Point p1 = node1().position();
      Point p2 = node2().position();
      double dist = (p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y) + (p1.z-p2.z)*(p1.z-p2.z);
      return std::sqrt(dist);
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, graph_->edges_[uid_].node1_uid);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, graph_->edges_[uid_].node2_uid);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      bool case1 = (node1() == e.node1()) and (node2() == e.node2());
      bool case2 = (node1() == e.node2()) and (node2() == e.node1());
      return (case1 or case2);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (node1() < node2()) {
        if (e.node1() < e.node2()) {
          if (node1() < e.node1()) {return true;}
          else if (e.node1() < node1()) {return false;}
          else {
            if (node2() < e.node2()) {return true;}
          }
          return false;
        }
        else {
          if (node1() < e.node2()) {return true;}
          else if (e.node2() < node1()) {return false;}
          else
          {
            if (node2() < e.node1()) {return true;}
          }
          return false;
        }
      }
      else {      
        if (e.node1() < e.node2()) {
          if (node2() < e.node1()) {return true;}
          else if (e.node1() < node2()) {return false;}
          else {
            if (node1() < e.node2()) {return true;}
          }
          return false;
        }
        else {
          if (node2() < e.node2()) {return true;}
          else if (e.node2() < node2()) {return false;}
          else
          {
            if (node1() < e.node1()) {return true;}
          }
          return false;
        }
      }

    }
    
    /** Return this edge's index, a number in the range [0, edge_size). */
    size_type index() const {
      return uid_; 
    }

    void switch_nodes() {
      size_type temp = graph_->edges_[uid_].node1_uid;
      graph_->edges_[uid_].node1_uid = graph_->edges_[uid_].node2_uid;
      graph_->edges_[uid_].node2_uid = temp;
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

    Edge(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_size_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < num_edges());
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
    internal_node internal_node_a = nodes_[a.index()];
    std::vector<size_type> adj_list_a = internal_node_a.adj_list;
    for (size_type i = 0; i < adj_list_a.size(); i++) {
      if (edges_[adj_list_a[i]].node1_uid == b.index()) {return true;}
      if (edges_[adj_list_a[i]].node2_uid == b.index()) {return true;}
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val = edge_value_type()) {
    // HW0: YOUR CODE HERE
    internal_node internal_node_a = nodes_[a.index()];
    std::vector<size_type> adj_list_a = internal_node_a.adj_list;
    for (size_type i = 0; i < adj_list_a.size(); i++) {
      if (edges_[adj_list_a[i]].node1_uid == b.index()) {
        return Edge(this, adj_list_a[i]);
      }
      if (edges_[adj_list_a[i]].node2_uid == b.index()) {
        return Edge(this, adj_list_a[i]);
      }
    }
    
    internal_edge new_edge;
    new_edge.node1_uid = a.index();
    new_edge.node2_uid = b.index();
    new_edge.uid = next_edge_uid_;
    new_edge.value = val;
    edges_.push_back(new_edge);
    ++edge_size_;
    nodes_[a.index()].adj_list.push_back(next_edge_uid_);
    nodes_[b.index()].adj_list.push_back(next_edge_uid_);
    ++next_edge_uid_;
    return Edge(this, next_edge_uid_-1); 
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    node_size_ = 0;
    next_node_uid_ = 0;
    edges_.clear();
    edge_size_ = 0;
    next_edge_uid_ = 0;
  }

  // Node Iterator

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    // Operator *
    //Node operator*() const {return *pos_;}
    Node operator*() const {
      Node n = graph_->node(pos_);
      return n;
    }
    
    // Operator ++
    NodeIterator& operator++() {
      pos_ ++;
      return *this;
    }

    // Operator ==
    bool operator==(const NodeIterator& it) const {return pos_ == it.pos_;}

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    NodeIterator (const Graph* graph, size_type uid) 
        : graph_(const_cast<Graph*>(graph)), pos_(uid) {
    }
    Graph* graph_;
    size_type pos_;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  
  // Iterator begin
  node_iterator node_begin() const {return NodeIterator(this, 0);} 
  
  // Iterator end
  node_iterator node_end() const {return NodeIterator(this, num_nodes());} 

  // Incident Iterator

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    // Operator *
    Edge operator*() const {
      std::vector<size_type> adj_list = graph_->nodes_[node_uid_].adj_list;
      Edge e = graph_->edge(adj_list[pos_]);
      if (e.node1().index() == node_uid_) {return e;}
      else {
        e.switch_nodes();
        return e;
      } 
    }
    
    // Operator ++ 
    IncidentIterator& operator++() {
      pos_++;
      return *this;
    }
    
    // Operator == 
    bool operator==(const IncidentIterator& it) const {
      return (pos_ == it.pos_ and node_uid_ == it.node_uid_);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    IncidentIterator (const Graph* graph, size_type uid, size_type node_uid) 
        : graph_(const_cast<Graph*>(graph)), pos_(uid), node_uid_(node_uid) {
    }
    Graph* graph_;
    size_type pos_;
    size_type node_uid_;
  };

  // Edge Iterator

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for: 
    Edge operator*() const {
      Edge e = graph_->edge(pos_);
      return e;
    }
   
    EdgeIterator& operator++() {
      pos_++;
      return *this;
    }

    bool operator==(const EdgeIterator& it) const {return pos_ == it.pos_;}

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type pos_;
    EdgeIterator (const Graph* graph, size_type uid) 
        : graph_(const_cast<Graph*>(graph)), pos_(uid) {
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const {return EdgeIterator(this, 0);}
  
  edge_iterator edge_end() const {return EdgeIterator(this, num_edges());}

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
