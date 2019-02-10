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

using namespace std;

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */

template <typename V>
class Graph {
 private:
  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Type of argument */
  using node_value_type = V;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    size_ = 0;
    edgeNum_ = 0;
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

    Node() {}

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->nodePositionVector_[index_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return index_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Return this node's value, @a V stored in a vector, with an opportunity to change it later*/
    node_value_type& value() {
      return graph_->nodeValueVector_[index_];
    }

    /** Return this node's value, @a V stored in a vector, as a const*/
    const node_value_type& value() const {
      return graph_->nodeValueVector_[index_];
    }

    /** Return this node's degree, the number of edges incident to it. Loops count as two edges. This @a with signature size_type is the size of the myEdges_ vector*/
    size_type degree() const {
      return myEdges_.size();     
    }

    /** Return an iterator pointer to the first edge of this Node. */
    incident_iterator edge_begin() const {
      return IncidentIterator(this, 0); 
    }

    /** Return an iterator pointer to the last edge of this Node. */
    incident_iterator edge_end() const {
      return IncidentIterator(this, myEdges_.size());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if(n.graph_ == graph_ && n.index() == index_) {
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
      // HW0: YOUR CODE HERE
      if (index_ < n.index()) {
        return true;
      } else {
        return false;
      } 
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    size_type index_;
    node_value_type value_;
    vector<Edge> myEdges_;
    Graph* graph_;

    Node(size_type ind, const Graph* graph)
      : index_(ind),
        graph_(const_cast<Graph*>(graph))
    {
      // HW0: YOUR CODE HERE
    }

    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return size_;
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
  Node add_node(const Point& position, const node_value_type& valueName = node_value_type()) {
    // HW0: YOUR CODE HERE
    nodePositionVector_.push_back(position);
    nodeValueVector_.push_back(valueName);

    size_ = nodeVector_.size();
    nodeVector_.push_back(Node(this, size_));
   
    ++size_;

    return nodeVector_[size_ - 1];
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (n.index() < size_ && n.graph_ == this) {
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
    // HW0: YOUR CODE HERE
    assert(i < size());
    return nodeVector_[i];
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return node1_;      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return node2_;      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (e.node1() == node1_ && e.node2() == node2_) {
        return true;
      } else if (e.node1() == node2_ && e.node2() == node1_) {
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
      Node tmp;

      if (node1_ < node2_) {
        tmp = node1_;
      } else {
        tmp = node2_;
      }

      if (e.node1() < tmp) {
        return false;
      } else if (e.node2() < tmp) {
        return false;
      } else {
        return true;
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    Node node1_;
    Node node2_;
    
    Edge(Node nodeA, Node nodeB)
      : node1_(nodeA),
        node2_(nodeB)
    {
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
    return edgeNum_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < edgeNum_);
    return edgeVector_[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    for (size_type j = 0; j < edgeNum_; ++j) {
      if (edgeVector_[j].node1() == a && edgeVector_[j].node2() == b) {
        return true;
      } else if (edgeVector_[j].node2() == a && edgeVector_[j].node1() == b) {
        return true;
      } else {
        return false;
      }
    }
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
    if (has_edge(a, b)) {
      for (size_type j = 0; j < edgeNum_; ++j) {
        if (edgeVector_[j].node1() == a && edgeVector_[j].node2() == b) {
        return edgeVector_[j];
      } else if (edgeVector_[j].node2() == a && edgeVector_[j].node1() == b) {
        (edgeVector_[j]).node1_ = a;
        (edgeVector_[j]).node2_ = b;
        return edgeVector_[j];
      }
      }
    } else {
      edgeNum_ = edgeVector_.size();
      edgeVector_.push_back(Edge(a, b));
      ++edgeNum_;

      a.myEdges_.push_back(Edge(a, b));
      b.myEdges_.push_back(Edge(a, b));

      return edgeVector_[edgeNum_ - 1];
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodeVector_.clear();
    size_ = 0;

    edgeVector_.clear();
    edgeNum_ = 0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator {
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
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Return the Node this iterator is pointing to.  */
    Node operator*() const {
      return graph_.nodeVector_[uid_];
    }

    /** Move NodeIterator to point to next Node in graph.
        @post: size_type uid = uid_old + 1

        return new NodeIterator
    */
    NodeIterator& operator++() {
      uid_ = uid_ + 1; 
      return this;
    }

    /** @brief check if this NodeIterator points to the same node as @a NodeIterator a
        @post If this == @a a, this.uid_ == @a a.uid_ && this.graph_ == @a a.graph_
              else             this.uid_ != @a a.uid_ || this.graph_ != @a a.graph_ 
        return bool
    */
    bool operator==(const NodeIterator& a) const {
      if (uid_ == a.uid_ && graph_ == a.graph_) {
        return true;
      } else {
        return false;
      }  
    }

    /** @brief check if this NodeIterator points to the same node as @a NodeIterator a
        @post If this != @a a, this.uid_ != @a a.uid_ || this.graph_ != @a a.graph_
              else             this.uid_ == @a a.uid_ && this.graph_ == @a a.graph_ 
        return bool
    */
    bool operator!=(const NodeIterator& a) const {
      if (uid_ == a.uid_ && graph_ == a.graph_) {
        return false;
      } else {
        return true;
      }
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type uid_;

    NodeIterator(const Graph* graph, size_type iter)
      : graph_(const_cast<Graph*>(graph)),
        uid_(iter)
    {
      
    }

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Return NodeIterator to this Graph's first node. */
  NodeIterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Return NodeIterator to this Graph's last node. */
  NodeIterator node_end() const {
    return NodeIterator(this, nodeVector_.size());
  } 

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator {
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
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** Return the Edge that this IncidentIterator is pointing to.  */
    Edge operator*() const {
      return node_.myEdges_[uid_]; 
    }

    /** Move IncidentIterator to point to next Edge around Node node_.
        @post: size_type uid = uid_old + 1

        return new IncidentIterator
    */
    IncidentIterator& operator++() {
      uid_ = uid_ + 1;
      return this;
    }

    /** @brief check if this IncidentIterator points to the same edge as @a IncidentIterator iter
        @param[in] IncidentIterator iter to be compared against
        @post If this == @a iter, this.uid_ == @a iter.uid_ && this.node_ == @a a.node_
              else             this.uid_ != @a a.uid_ || this.node_ != @a a.node_ 
        return bool
    */
    bool operator==(const IncidentIterator& iter) const {
      if (uid_ == iter.uid_ && node_ == iter.node_) {
        return true;
      } else {
        return false;
      } 
    }

    /** @brief check if this IncidentIterator points to the same edge as @a IncidentIterator iter
        @param[in] IncidentIterator iter to be compared against
        @post If this != @a iter, this.uid_ != @a iter.uid_ || this.node_ != @a a.node_
              else             this.uid_ == @a a.uid_ && this.node_ == @a a.node_ 
        return bool
    */
    bool operator!=(const IncidentIterator& iter) const {
      if (uid_ == iter.uid_ && node_ == iter.node_) {
        return false;
      } else {
        return true;
      }
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Node* node_;
    size_type uid_;

    IncidentIterator(const Node* node, size_type itern)
      : node_(const_cast<Node*>(node)),
        uid_(itern)
    {
      
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator {
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
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /** Return the Edge that this IncidentIterator is pointing to.  */    
    Edge operator*() const {
      return graph_.edgeVector[uid_];
    }

    /** Move EdgeIterator to point to next Edge in Graph graph_.
        @post: size_type new uid = uid_old + 1

        return new EdgeIterator
    */
    EdgeIterator& operator++() {
      uid_ = uid_ + 1;
      return this;
    }

    /** @brief check if this EdgeIterator points to the same edge as @a EdgeIterator eit
        @param[in] EdgeIterator eit to be compared against
        @post If this == @a eit, this.uid_ == @a iter.uid_ && this.graph_ == @a a.graph_
              else             this.uid_ != @a a.uid_ || this.graph_ != @a a.graph_ 
        return bool
    */
    bool operator==(const EdgeIterator& eit) const {
      if (uid_ == eit.uid_ && graph_ == eit.graph_) {
        return true;
      } else {
        return false;
      }
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    size_type uid_;
    Graph graph_;

    EdgeIterator(const Graph* graph, size_type iter)
      : graph_(const_cast<Graph*>(graph)),
        uid_(iter)
    {

    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  EdgeIterator edge_begin() const {
    return EdgeIterator(this, (size_type) 0);
  }

  EdgeIterator edge_end() const {
    return EdgeIterator(this, edgeNum_);
  }

 private:
  size_type size_;
  size_type edgeNum_;

  vector<Node> nodeVector_;
  vector<Point> nodePositionVector_;
  vector<node_value_type> nodeValueVector_;
  vector<Edge> edgeVector_;


  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
