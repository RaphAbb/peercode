#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <map>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph : private totally_ordered<Graph<V>> {
 struct internal_node;
 struct internal_edge;
 
 private: 

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V>;
  
  /** Synonym for V. */
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
  class Node: private totally_ordered<Node> {
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
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_[index_].position; //access position via proxy vector of internal_node's
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return index_;
    }

    /** Return this node's value.
    * @return a reference to the Node's value.
    */
    V& value() {
      //access value via proxy vector of internal_node's
      return graph_->nodes_[index_].value; 
    }

    /** Return this node's value.
    * @return a copy of the Node's value
    */
    const V& value() const {
      //access value via proxy vector of internal_node's
      return graph_->nodes_[index_].value; 
    }

    /** Return this node's degree.
    * @return the number of edges incident to this node.
    */
    size_type degree() const {
      return (graph_->adjacency_[index_]).size(); 
    }    

    /** Return iterator to the first edge incident to this node.
    * @return  An iterator to the first edge incident to this node.
    */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, (graph_->adjacency_[index_]).begin());
    }    

    /** Return iterator to the last edge incident to this node.
    * @return An iterator to the last edge incident to this node
    */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, (graph_->adjacency_[index_]).end());
    }

    /** Test whether this node and @a n are equal.
     * @param n   The node the compare to.
     * @return    True if the nodes are equal, meaning
     *            they have the same graph and the same index,
     *            and false otherwise.
     */
    bool operator==(const Node& n) const {
      return ((graph_ == n.graph_)&&(index_ == n.index_));
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
      return (index_ < n.index_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    friend class NodeIterator;
    friend class IncidentIterator;
    
    // Pointer back to the Graph
    Graph* graph_;
    
    //index of node in graph.nodes
    size_type index_;

    // value of node in graph.nodes
    node_value_type value_;

    /** Private Constructor */
    Node(const Graph* graph, size_type index, node_value_type value = V()):graph_(const_cast<Graph*>(graph)), index_(index), value_(value) {};

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes_.size(); //accessing member of vector is O(1)
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
  Node add_node(const Point& position, const V& value = V()) {
    size_type index = num_nodes();
    //add internal_node proxy; inserting is O(1) amortized
    nodes_.push_back({position, index, value}); 

    //map in adjacency vector will store <neighbor, edge> pairs
    adjacency_.push_back(std::map<size_type, size_type>()); 
    return Node(this, index);
  }
  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.index() >= num_nodes())  
      return false; //invalid index for graph
    else {
      internal_node n2 = nodes_[n.index()]; //internal_node proxy
      return ((n2.position==n.position())&&(n2.index==n.index())); 
    }
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i >= 0);
    assert(i < num_nodes());
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //equal edges should be in the same graph and connect the same nodes (regardless)
      //of the order of the connection They will also have the same index since each
      //edge is added only once.
      return ((graph_ == e.graph_)&&(((node1_ == e.node1_)&&(node2_ == e.node2_)) ||
          ((node1_ == e.node2_)&&(node2_ == e.node1_))));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (index_ < e.index_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph; 
    // Pointer back to the Graph
    Graph* graph_;

    // index of edge in graph.edges_
    size_type index_;

    // index of node1 in graph.nodes_
    size_type node1_;

    // index of node2 in graph.nodes_
    size_type node2_;

    /** Private Constructor */
    Edge(const Graph* graph, size_type index, size_type node1, size_type node2):graph_(const_cast<Graph*>(graph)), index_(index), node1_(node1), node2_(node2) {};
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_.size(); //O(1) to get size of vector

  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i >= 0);
    assert(i < num_edges());
    //determine edge's nodes via proxy vector of internal_edge's; O(1) for vector operator []
    return Edge(this, i, (edges_[i]).node1, (edges_[i]).node2); 
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    //maps neighbors of a to index of connecting edge in proxy edges_
    std::map<size_type, size_type> a_edges = adjacency_[a.index()]; 

    //edge exists iff b in neighbors of a
    return (a_edges.find(b.index()) != a_edges.end()); 
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
    size_type i; //index of connecting edge in proxy edges_
    if (has_edge(a, b))
      i = (adjacency_[a.index()])[b.index()];  
    else {
      i = num_edges(); 
      //add that b is a neighbor of a, connected by edge i; O(num_edges()) to insert into map of size O(num_edges())
      (adjacency_[a.index()]).insert(std::pair<size_type, size_type>(b.index(), i)); 

      //add that a is a neighbor of b, connected by edge i; O(num_edges()) as well
      (adjacency_[b.index()]).insert(std::pair<size_type, size_type>(a.index(), i)); 

      //add internal_edge i to edges_ as proxy for edge i; O(1) for vector push_back
      edges_.push_back({a.index(), b.index(), i}); 
    }
    return Edge(this, i, a.index(), b.index());
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    adjacency_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator> {
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

    /** Get the Node to which this iterator points.
    *@return    The Node to which this iterator points.
    */
    Node operator*() const {
      return Node(graph_, (*internal_iterator_).index, (*internal_iterator_).value);
    }

    /** Set this iterator to point to the next Node in the graph.
    *@return    An iterator to the next Node in the graph.
    */
    NodeIterator& operator++() {
      internal_iterator_++;
      return *this;
    }    

    /** Check if this iterator equals a different iterator
    *@param other    The iterator to which to compare.
    *@return         True if this iterator and @a other are equal,
    *                meaning that they point to the same object in
    *                the same graph, and false otherwise..
    */
    bool operator==(const NodeIterator& other) const {
      return ((graph_==other.graph_)&&(internal_iterator_==other.internal_iterator_));
    }

    private:
     friend class Graph;
    /** Private Constructor */
    NodeIterator(const Graph* graph, typename std::vector<internal_node>::const_iterator internal_iterator):graph_(const_cast<Graph*>(graph)), internal_iterator_(internal_iterator) {};
    Graph* graph_;
    typename std::vector<internal_node>::const_iterator internal_iterator_;
  };

  /** Return an iterator to the first node in the Graph.
  * @return  An iterator to the first node in the Graph.
  */
  node_iterator node_begin() const {
    return NodeIterator(this, nodes_.begin());
  }

  /** Return an iterator to the last node in the Graph.
  * @return  An iterator to the last node in the Graph.
  */
  node_iterator node_end() const {
    return NodeIterator(this, nodes_.end());
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
    }

     /** Get the Edge to which this iterator points.
    *@return    The Edge to which this iterator points,
    *           where node1() of the edge is the node 
    *           that spawned the incident iterator and
    *           node2() returns the adjacent node.
    */
    Edge operator*() const {
      size_type node2 = (*internal_iterator_).first;
      size_type edge_index = (*internal_iterator_).second;

      //get the representation of this edge from the proxy vector
      internal_edge edge = graph_->edges_[edge_index];

      //determine which node of the edge should be returned as
      //node1() and which node should be returned as node2()
      size_type nodea = edge.node1;
      size_type nodeb = edge.node2;
      size_type node1;
      if (nodea != node2) {
        //nodea is the node that spanwed the iterator
        node1 = nodea;
      }
      else {
        //nodeb is the node that spanwed the iterator
        node1 = nodeb;
      }
      return Edge(graph_, edge_index, node1, node2);
    }

    /** Set this iterator to point to the next Edge in the graph
    *   that is incident to this node.
    *@return    An iterator to the next Edge in the graph
    *           that is incident to the node.
    */
    IncidentIterator& operator++() {
      internal_iterator_++;
      return *this;
    }


    /** Check if this iterator equals a different iterator
    *@param other    The iterator to which to compare.
    *@return         True if this iterator and @a other are equal,
    *                meaning that they point to the same object in
    *                the same graph, and false otherwise..
    */
    bool operator==(const IncidentIterator& other) const {
      return ((graph_==other.graph_)&&(internal_iterator_ == other.internal_iterator_));
    }

   private:
    friend class Graph;
    
    Graph* graph_;
    typename std::map<size_type, size_type>::const_iterator internal_iterator_;

    /** Private Constructor */
    IncidentIterator(const Graph* graph, typename std::map<size_type, size_type>::const_iterator internal_iterator):graph_(const_cast<Graph*>(graph)), internal_iterator_(internal_iterator) {};
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator> {
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

    /** Get the Edge to which this iterator points.
    *@return    The Edge to which this iterator points.
    */
    Edge operator*() const {
      return Edge(graph_, (*internal_iterator_).index, (*internal_iterator_).node1, (*internal_iterator_).node2);
    }

    /** Set this iterator to point to the next Edge in the graph.
    *@return    An iterator to the next Edge in the graph, where
    *           incrementing the iterator from beginning to end
    *           will iterate over each undirected edge exactly once.
    */
    EdgeIterator& operator++() {
      internal_iterator_++;
      return *this;
    }    

    /** Check if this iterator equals a different iterator
    *@param other    The iterator to which to compare.
    *@return         True if this iterator and @a other are equal,
    *                meaning that they point to the same object in
    *                the same graph.
    */
    bool operator==(const EdgeIterator& other) const {
      return ((graph_==other.graph_)&&(internal_iterator_==other.internal_iterator_));
    }

   private:
    friend class Graph;
    
    /** Private Constructor */
    EdgeIterator(const Graph* graph, typename std::vector<internal_edge>::const_iterator internal_iterator):graph_(const_cast<Graph*>(graph)), internal_iterator_(internal_iterator) {};
    Graph* graph_; 
    typename std::vector<internal_edge>::const_iterator internal_iterator_;
  };

  /** Return iterator to the first edge in the graph.
  * @return  An iterator to the first edge in the graph.
  */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, edges_.begin());
  }

  /** Return iterator to the last edge in the graph..
  * @return  An iterator to the last edge in the graph.
  */
  edge_iterator edge_end() const {
    return EdgeIterator(this, edges_.end());
  }

 private:
  struct internal_node {
    const Point position;
    const size_type index; //index in proxy vector nodes_
    node_value_type value; //value in proxy vector nodes_
  };
  struct internal_edge {
    const size_type node1; //index in proxy vector nodes_
    const size_type node2; //index in proxy vector nodes_
    const size_type index; //index in proxy vector edges_
  };

  std::vector<internal_edge> edges_; //proxy for Edge's
  //adjacency vector, map in index a with entry <b, i> denotes edge i between nodes a and b
  std::vector<std::map<size_type, size_type>> adjacency_; 
  std::vector<internal_node> nodes_; //proxy for Node's
};

#endif // CME212_GRAPH_HPP
