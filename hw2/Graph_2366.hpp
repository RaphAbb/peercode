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

template <typename V, typename E>
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
  using edge_value_type = E;

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

    Point& position() {
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
      return graph_->myEdges_[index_].size();     
    }

    /** Return an iterator pointer to the first edge of this Node. */
    incident_iterator edge_begin() const {
      return IncidentIterator(this, 0);
    }

    /** Return an iterator pointer to the last edge of this Node. */
    incident_iterator edge_end() const {
      return IncidentIterator(this, this->degree());
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
    nodeVector_.push_back(Node(size_, this));
    myEdges_.push_back({});  
 
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
      return graph_->nodeVector_[nodeID1_];   
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph_->nodeVector_[nodeID1_];
    }

    /** Return the value of this Edge */
    edge_value_type& value() {
      return graph_->edgeValueVector_[eidx_];
    }

    /** Return the value of this Edge */
    const edge_value_type& value() const {
      return graph_->edgeValueVector_[eidx_];
    }

    /** Return the Eucledian distance between positions of this Edge's Nodes */
    double length() const {
      Node node1_ = graph_->nodeVector_[nodeID1_];
      Node node2_ = graph_->nodeVector_[nodeID2_];

      Point n1 = node1_.position();
      Point n2 = node2_.position();

      double x1 = n1[0];
      double y1 = n1[1];
      double z1 = n1[2];

      double x2 = n2[0];
      double y2 = n2[1];
      double z2 = n2[2];

      double edge_length = sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2) + pow((z1 - z2), 2));

      return edge_length;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (e.node1() == graph_->nodeVector_[nodeID1_] && e.node2() == graph_->nodeVector_[nodeID2_]) {
        return true;
      } else if (e.node1() == graph_->nodeVector_[nodeID2_] && e.node2() == graph_->nodeVector_[nodeID1_]) {
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

      if (graph_->nodeVector_[nodeID1_] < graph_->nodeVector_[nodeID2_]) {
        tmp = graph_->nodeVector_[nodeID1_];
      } else {
        tmp = graph_->nodeVector_[nodeID2_];
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
    Graph* graph_;
    size_type eidx_;
    size_type nodeID1_;
    size_type nodeID2_;   
 
    Edge(Graph* graph, size_type eidx, size_type nodeIDA, size_type nodeIDB)
      : graph_(graph),
        eidx_(eidx),
        nodeID1_(nodeIDA),
        nodeID2_(nodeIDB)
    {} 

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
    bool os;
    for (size_type j = 0; j < edgeNum_; ++j) { // Go through all edges
      os = false; // Assume false

      if (edgeVector_[j].node1() == a && edgeVector_[j].node2() == b) { // If jth edge has nodes a and b
        os = true; // Set true
        break;
      } else if (edgeVector_[j].node2() == a && edgeVector_[j].node1() == b) { // If jth edge has nodes b and a
        os = true; // Set true
        break;
      }
    }
    return os;
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& valueName = edge_value_type()) {
    // HW0: YOUR CODE HERE
    if (has_edge(a, b)) {
      for (size_type j = 0; j < edgeNum_; ++j) {
        if (edgeVector_[j].node1() == a && edgeVector_[j].node2() == b) {
        return edgeVector_[j];
        break;
        } else if (edgeVector_[j].node2() == a && edgeVector_[j].node1() == b) {
        (edgeVector_[j]).node1() = a;
        (edgeVector_[j]).node2() = b;
        return edgeVector_[j];
        break;
        }
      }
    }

    edgeNum_ = edgeVector_.size();
    Edge new_edge = Edge(this, edgeNum_, a.index(), b.index());
    edgeVector_.push_back(new_edge);
    edgeValueVector_.push_back(valueName);
    ++edgeNum_;

    myEdges_[a.index()].push_back(new_edge);
    myEdges_[b.index()].push_back(new_edge);
      
    return edgeVector_[new_edge.eidx_];
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodeVector_.clear();
    nodePositionVector_.clear();
    nodeValueVector_.clear();
    size_ = 0;

    edgeVector_.clear();
    edgeValueVector_.clear();
    edgeNum_ = 0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
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
      return *it_;
    }

    /** Move NodeIterator to point to next Node in graph.
        @post: size_type uid = uid_old + 1

        return new NodeIterator
    */
    NodeIterator& operator++() {
      it_++; 
      return *this;
    }

    /** @brief check if this NodeIterator points to the same node as @a NodeIterator a
        @post If this == @a a, this.it_ == @a a.it_ && this.graph_ == @a a.graph_
              else             this.it_ != @a a.it_ || this.graph_ != @a a.graph_ 
        return bool
    */
    bool operator==(const NodeIterator& a) const {
      if (it_ == a.it_ && graph_ == a.graph_) {
        return true;
      } else {
        return false;
      }  
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    typename vector<Node>::const_iterator it_;

    NodeIterator(const Graph* graph, typename vector<Node>::const_iterator iter)
      : graph_(const_cast<Graph*>(graph)),
        it_(iter)
    {}
    NodeIterator(const Graph* graph, typename vector<Node>::iterator iter)
      : graph_(const_cast<Graph*>(graph)),
        it_(iter)
    {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Return NodeIterator to this Graph's first node. */
  NodeIterator node_begin() const {
    return NodeIterator(this, nodeVector_.begin());
  }

  /** Return NodeIterator to this Graph's last node. */
  NodeIterator node_end() const {
    return NodeIterator(this, nodeVector_.end());
  } 

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
      node_ = nullptr;
      uid_ = 0;
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** Return the Edge that this IncidentIterator is pointing to.  */
    Edge operator*() const {
      return ((*node_).graph_)->myEdges_[(*node_).index()][uid_]; 
    }

    /** Move IncidentIterator to point to next Edge around Node node_.
        @post: size_type new uid = old uid + 1

        return new IncidentIterator
    */
    IncidentIterator& operator++() {
      uid_++;
      return *this;
    }

    /** @brief check if this IncidentIterator points to the same edge as @a IncidentIterator iter
        @param[in] IncidentIterator iter to be compared against
        @post If this == @a iter, this.uid_ == @a iter.uid_ && this.node_ == @a a.node_
              else                this.uid_ != @a a.uid_ || this.node_ != @a a.node_ 
        return bool
    */
    bool operator==(const IncidentIterator& iter) const {
      if (uid_ == iter.uid_ && node_ == iter.node_) {
        return true;
      } else {
        return false;
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
    {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
      it_ = nullptr;
      graph_ = nullptr;
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /** Return the Edge that this IncidentIterator is pointing to.  */    
    Edge operator*() const {
      return *it_;
    }

    /** Move EdgeIterator to point to next Edge in Graph graph_.
        @post: size_type new uid = old uid + 1

        return new EdgeIterator
    */
    EdgeIterator& operator++() {
      it_++;
      return *this;
    }

    /** @brief check if this EdgeIterator points to the same edge as @a EdgeIterator eit
        @param[in] EdgeIterator eit to be compared against
        @post If this == @a eit, this.it_ == @a iter.it_ && this.graph_ == @a a.graph_
              else               this.it_ != @a a.it_ || this.graph_ != @a a.graph_ 
        return bool
    */
    bool operator==(const EdgeIterator& eit) const {
      if (it_ == eit.it_ && graph_ == eit.graph_) {
        return true;
      } else {
        return false;
      }
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const Graph* graph_;
    typename vector<Edge>::const_iterator it_;

    EdgeIterator(const Graph* graph, typename vector<Edge>::const_iterator iter)
      : graph_(const_cast<Graph*>(graph)),
        it_(iter)
    {}

    EdgeIterator(const Graph* graph, typename vector<Edge>::iterator iter)
      : graph_(const_cast<Graph*>(graph)),
        it_(iter)
    {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Return EdgeIterator to first Edge in Graph */
  EdgeIterator edge_begin() const {
    return EdgeIterator(this, edgeVector_.begin());
  }

  /** Return EdgeIterator to last Edge in Graph */
  EdgeIterator edge_end() const {
    return EdgeIterator(this, edgeVector_.end());
  }

  /** @brief Method to remove Node from Graph
   *
   * param[in] a Node to be removed
   * return    1 or 0
   *
   * @pre this.has_node(a)
   *
   * @post !(this.has_node(a))
   * @post new (this.size()) = old (this.size()) - 1
   * @post new this.nodeVector_[a.index()] = old this.nodeVector_[old this.size() - 1] != a
   * @post new this.nodePositionVector_[a.index()] = old this.nodePositionVector_[old this.size() - 1] != a.position()
   * @post new this.nodeValueVector_[a.index()] = old this.nodeValueVector_[old this.size() - 1] != a.value()
   * @post new this.myEdges_[a.index()] = old this.myEdges_[old this.size() - 1] != old this.myEdges_[a.index()]
   * @post for all b such that old this.has_edge(a, b), new (this.myEdges_[b.index()]).size()  = old (this.myEdges_[b.index()]).size() - 1
   * @post for all b such that this.has_node(b), !(this.has_edge(a, b))
   * @post new edgeNum_ <= old edgeNum_
   *
   * Complexity: O(a.degree() + num_edges())
   */
  size_type remove_node(const Node& a) {
    if (has_node(a)) { // Checks if Node is in Graph
      size_type this_idx = a.index();
      nodePositionVector_[this_idx] = nodePositionVector_[size_ - 1]; // Swap Position with last value in vector
      nodePositionVector_.pop_back(); // Pop last value
      nodeValueVector_[this_idx] = nodeValueVector_[size_ - 1]; // Swap Value with last value in vector
      nodeValueVector_.pop_back(); // Pop last value

      for (auto it = a.edge_begin(); it != a.edge_end(); ++it) { // For all incident edges
        remove_edge(*it); // remove edge
       
        if (myEdges_[this_idx].size() == 0) {break;} // if no incident edges get out of loop
      }

      for (auto it = (nodeVector_[size() - 1]).edge_begin(); it != (nodeVector_[size() - 1]).edge_end(); ++it) {
        edgeVector_[(*it).eidx_].nodeID1_ = this_idx;
      }

      myEdges_[this_idx] = myEdges_[size_-1]; // Swap empty vector of incident Edges with last vector
      myEdges_.pop_back(); // Pop last vector
      nodeVector_[this_idx] = nodeVector_[size_ - 1]; // Swap Node with last Node in vector
      (nodeVector_[this_idx]).index_ = this_idx; // Reassign index
      nodeVector_.pop_back(); // Pop last value in vector
      size_--; // Reduce size

      return (size_type) 1;
    } else {
      return (size_type) 0;
    }
  }

  /** @brief Method to remove Node from Graph
   *
   * param[in] n_it      NodeIterator to Node that should be removed
   * return    new n_it  NodeIterator to new Node
   *
   * @pre this.has_node(*(old n_it))
   *
   * @post new (this.size()) = old (this.size()) - 1
   * @post new this.nodeVector_[(*(new n_it)).index()] = old this.nodeVector_[old this.size() - 1] != *(old n_it)
   * @post new this.nodePositionVector_[*(new n_it).index()] = old this.nodePositionVector_[old this.size() - 1] != *(old n_it).position()
   * @post new this.nodeValueVector_[*(new n_it).index()] = old this.nodeValueVector_[old this.size() - 1] != *(old n_it).value()
   * @post new this.myEdges_[*(new n_it).index()] = old this.myEdges_[old this.size() - 1] != old this.myEdges_[*(old n_it).index()]
   * @post for all Nodes b such that old this.has_edge(*(old n_it), b), new (this.myEdges_[b.index()]).size()  = old (this.myEdges_[b.index()]).size() - 1
   * @post for all Nodes b such that this.has_node(b), !(this.has_edge(*(old n_it), b))
   * @post new edgeNum_ <= old edgeNum_
   *
   * Complexity: O(*(old n_it).degree() + num_edges())
   */
  node_iterator remove_node(node_iterator n_it) {
    Node this_n = *n_it;
    remove_node(this_n); // Call remove_node(Node& a) function
    return n_it;
  }

  /** @brief Method to remove Edge from Graph
   *
   * param[in] a, b  Nodes forming the Edge which should be removed
   * return    1 or 0
   * 
   * @pre has_edge(a, b)
   *
   * @post new num_edges() = old num_edges() - 1
   * @post new size() = old size()
   * @post has_node(a) && has_node(b)
   * @post !has_edge(a, b)
   * @post For i such that old edge(i) = e && e.node1() = a && e.node2() = b,   new edge(i) = old edge(old num_edges()) != old edge(i)
   *
   * Note: It is an invalid operation to request the value() of a removed edge. This memory is released when removing the edge
   *
   * Complexity(a.degree() + num_edges())
   */
  size_type remove_edge(const Node& a, const Node& b) {
    if (has_edge(a, b)) { // Check if it has Edge
      for (auto it = a.edge_begin(); it != a.edge_end(); ++it) { // Go through all incident Edges of a 
        if ((*it).node2() == b) { // If Edge's other vertex is b
          size_type tORf = remove_edge(*it); // Call remove_edge(Edge& e)
          (void) tORf;
          break;
        }
      }
      return (size_type) 1;
    } else {
      return (size_type) 0;
    }
  }

  /** @brief Method to remove Edge from Graph
   *
   * param[in] e      Edge which should be removed
   * return    1 or 0
   * 
   * @pre For a = e.node1() & b = e.node2(),  has_edge(a, b)
   *
   * @post new num_edges() = old num_edges() - 1
   * @post new size() = old size()
   * @post For a = e.node1() and b = e.node2(), has_node(a) && has_node(b)
   * @post For a = e.node1() and b = e.node2(), !has_edge(a, b)
   * @post For i such that old edge(i) = e,   new edge(i) = old edge(old num_edges()) != old edge(i)
   *
   * Note: It is an invalid operation to request the value() of a removed edge. This memory is released when removing the
 edge
   *
   * Complexity(num_edges())
   */
  size_type remove_edge(const Edge& e) {
    size_type this_idx = e.eidx_; // Extract information from Edge
    size_type a_idx = e.nodeID1_;
    size_type b_idx = e.nodeID2_;

    Node a = nodeVector_[a_idx]; // Find Edge's Nodes
    Node b = nodeVector_[b_idx];

    if(has_edge(a, b)) { // Check if Graph has Edge
      if (edgeNum_ == 1) { // If it's only Edge
        edgeValueVector_.pop_back(); // Simply remove it from everywhere
        myEdges_[a_idx].pop_back();
        myEdges_[b_idx].pop_back();
        edgeVector_.pop_back();
      } else {
        edgeValueVector_[this_idx] = edgeValueVector_[edgeNum_-1]; // Swap with last Edge's value
        edgeValueVector_.pop_back(); // Pop last element

        if (a.degree() == 1) { // If only Edge of Node a
          (myEdges_[a_idx]).pop_back(); // Simply remove it
        } else {
          for (auto it = a.edge_begin(); it != a.edge_end(); ++it) { // Go through all Edges of a
            if ((*it) == e) { // Find Edge e
              myEdges_[a_idx][it.uid_] = myEdges_[a_idx][a.degree()-1]; // Swap with last Edge in vector
              (myEdges_[a_idx]).pop_back(); // Pop last element
              break;
            }
          }
        }

        if (b.degree() == 1) { // If only Edge of Node b
          (myEdges_[b_idx]).pop_back(); // Simply remove it
        } else {
          for (auto it = b.edge_begin(); it != b.edge_end(); ++it) { // Go through all Edges of b
            if ((*it) == e) { // Find Edge e
              myEdges_[b_idx][it.uid_] = myEdges_[b_idx][b.degree()-1]; // Swap with last Edge in vector
              (myEdges_[b_idx]).pop_back(); // Pop last element
              break;
            }
          }
        }
        edgeVector_[this_idx] = edgeVector_[edgeNum_ - 1]; // Swap with last element in general Edge vector
        (edgeVector_[this_idx]).eidx_ = this_idx; // Reassign Edge index
        edgeVector_.pop_back(); // Pop last element
      }

      edgeNum_--; // Reduce Edge counter
      return (size_type) 1;
    } else {
      return (size_type) 0;
    }
  }

  /** @brief Method to remove Edge from Graph
   *
   * param[in] e_it   EdgeIterator to Edge which should be removed
   * return    new e_it
   * 
   * @pre For a = old (*e_it).node1() & b = old (*e_it).node2(),  has_edge(a, b)
   *
   * @post new num_edges() = old num_edges() - 1
   * @post new size() = old size()
   * @post For a = old (*e_it).node1() and b = old (*e_it).node2(), has_node(a) && has_node(b)
   * @post For a = old (*e_it).node1() and b = old (*e_it).node2(), !has_edge(a, b)
   * @post For i such that old edge(i) = old (*e_it),   new edge(i) = old edge(old num_edges()) != old edge(i)
   *
   * Note: It is an invalid operation to request the value() of a removed edge. This memory is released when removing the
 edge
   *
   * Complexity(num_edges())
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    Edge this_e = *e_it;
    remove(this_e);
    return e_it;
  }

 private:
  size_type size_;
  size_type edgeNum_;

  vector<Node> nodeVector_;
  vector<Point> nodePositionVector_;
  vector<vector<Edge>> myEdges_;
  vector<node_value_type> nodeValueVector_;
  vector<Edge> edgeVector_;
  vector<edge_value_type> edgeValueVector_;

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
