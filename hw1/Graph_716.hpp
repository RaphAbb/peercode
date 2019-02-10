#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
//#include <iostream>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


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

  // Node Point data (index is node id)
  std::vector<Point> Points;
  // Edge Node id's data (index is edge id)
  std::vector<unsigned> EdgeA_ids;
  std::vector<unsigned> EdgeB_ids;
  // adjacency list: Node id -> connected Node id's
  std::vector<std::vector<unsigned>> Nadj_list;
  // adjacency list: Node id -> corresponding Edge id's to the Node id's above
  std::vector<std::vector<unsigned>> Eadj_list;
  // storage for node value type
  std::vector<V> Values;
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

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Synonym for Node Value type. */
  using node_value_type = V;

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
    Node() {
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->Points[id_];
    }

    /** Return this node's index, a number in the range [0, num_Nodes). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return id_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if(graph_ == n.graph_)
        return id_ == n.index();
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
      return id_ < n.index();
    }

    /** returns the reference to the node value */
    node_value_type& value() {
      return graph_->Values[id_];
    }

    /** returns the constant reference to the node value */
    const node_value_type& value() const {
      return graph_->Values[id_];
    }

    /** returns the degree of the node, the number of edges incident to it */
    size_type degree() const {
      return graph_->Eadj_list[id_].size();
    }

    /** returns an iterator pointing to the first edge that is incident to this node.
    *  @return A IncidentIterator object
    */
    IncidentIterator edge_begin() const {
      return IncidentIterator(graph_, id_, 0);
    }

    /** returns an iterator pointing to the last edge that is incident to this node.
    *  @return A IncidentIterator object
    */
    IncidentIterator edge_end() const {
      return IncidentIterator(graph_, id_, degree());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // internal Node variables
    Graph* graph_;
    size_type id_;
    /** Private Constructor */
    Node(const Graph* g_in, size_type id_in) : graph_(const_cast<Graph*>(g_in)), id_(id_in) {}
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
    return num_Nodes;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] val The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    // HW0: YOUR CODE HERE
    // updates node information
    Points.push_back(position);
    Values.push_back(val);

    // increases size of adjacency lists
    std::vector<unsigned> newvec;
    Nadj_list.push_back(newvec);
    Eadj_list.push_back(newvec);

    num_Nodes++;
    return Node(this,num_Nodes-1);
  }


  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    // also makes sure if Node id is not too large
    if(n.graph_ == this)
      return n.id_ < size();
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
    assert(i < size());
    return Node(this,i);
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

    /** Return this edge's index, a number in the range [0, num_Edges). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return id_;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      // retrieves Node id
      size_type A_i = graph_->EdgeA_ids[id_];
      return graph_->node(A_i);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      // retrieves Node id
      size_type B_i = graph_->EdgeB_ids[id_];
      return graph_->node(B_i);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // compares node id's in both directions
      if(graph_->EdgeA_ids[id_] == graph_->EdgeA_ids[e.index()])
        if(graph_->EdgeB_ids[id_] == graph_->EdgeB_ids[e.index()])
          return true;
      if(graph_->EdgeA_ids[id_] == graph_->EdgeB_ids[e.index()])
        return (graph_->EdgeB_ids[id_] == graph_->EdgeA_ids[e.index()]);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return id_ < e.index();
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // internal Edge variables
    Graph* graph_;
    size_type id_;
    /** Private Constructor */
    Edge(const Graph* g_in, size_type id_in) : graph_(const_cast<Graph*>(g_in)), id_(id_in) {}
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
    return num_Edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i<num_edges());
    return Edge(this,i);        // Valid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    // Uses the node adjacency list to determine if edge already exists
    for(size_type i = 0; i < Nadj_list[a.index()].size(); i++)
    {
      if(Nadj_list[a.index()][i] == b.index())
      {
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
    // Uses the node adjacency list to determine if edge exists
    for(size_type i = 0; i < Nadj_list[a.index()].size(); i++)
    {
      if(Nadj_list[a.index()][i] == b.index())
      {
        // determines the edge id corresponding to that position on the node adjacency list
        size_type E_ind = Eadj_list[a.index()][i];
        // returns existing edge
        if(EdgeA_ids[E_ind] == a.index())
          if(EdgeB_ids[E_ind] == b.index())
          {
            return Edge(this,E_ind);
          }
        // switches a,b and returns existing edge with nodes switched
        if(EdgeB_ids[E_ind] == a.index())
          if(EdgeA_ids[E_ind] == b.index())
          {
            EdgeA_ids[E_ind] = a.index();
            EdgeB_ids[E_ind] = b.index();
            return Edge(this,E_ind);
          } 
      }
    }

    // else, creates and returns new edge
    // updates adjacency lists
    Nadj_list[a.index()].push_back(b.index());
    Nadj_list[b.index()].push_back(a.index());
    Eadj_list[a.index()].push_back(num_Edges);
    Eadj_list[b.index()].push_back(num_Edges);
    // adds edge id data
    EdgeA_ids.push_back(a.index());
    EdgeB_ids.push_back(b.index());
    num_Edges++;
    return Edge(this,num_Edges-1);  
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    Points.clear();
    EdgeA_ids.clear();
    EdgeB_ids.clear();
    Nadj_list.clear();
    Eadj_list.clear();
    num_Nodes = 0;
    num_Edges = 0;
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

    /** Constructor for a NodeIterator
      @param g_in Pointer to the current graph
      @param Niter_in A number representing the current position of the iterator
    */
    NodeIterator(const Graph* g_in, size_type Niter_in) : graph_(const_cast<Graph*>(g_in)), Niter_(Niter_in) {}

    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Returns the Node that this iterator is pointing to.
      @return A Node object
    */
    Node operator*() const {
      return Node(graph_,Niter_);
    }

    /** Increments the iterator. 
      @return A NodeIterator object
    */
    NodeIterator& operator++() {
      Niter_++;
      return *this;
    }

    /** Checks for equality of this NodeIterator and @a NI_in. */
    bool operator==(const NodeIterator& NI_in) const {
      return Niter_==NI_in.Niter_;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type Niter_;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Returns an iterator pointing to the first node in the graph.
  *  @return A NodeIterator object
  */
  NodeIterator node_begin() const {
    return NodeIterator(this,0);
  }

  /** Returns an iterator pointing to the last node in the graph.
  *  @return A NodeIterator object
  */
  NodeIterator node_end() const {
    return NodeIterator(this,size());
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

    // HW1 #3: YOUR CODE HERE

    /** Constructor for a IncidentIterator
      @param g_in Pointer to the current graph
      @param Nid_in The node ID of the current node
      @param Niter_in A number representing the current position of the iterator
    */
    IncidentIterator(const Graph* g_in, size_type Nid_in, size_type Niter_in) :
           graph_(const_cast<Graph*>(g_in)), Nid_(Nid_in), Niter_(Niter_in) {}

    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** Returns the Edge that this iterator is pointing to.
      @return An Edge object
    */
    Edge operator*() const {
      return Edge(graph_,graph_->Eadj_list[Nid_][Niter_]);
    }

    /** Increments the iterator. 
      @return An IncidentIterator object
    */
    IncidentIterator& operator++() {
      Niter_++;
      return *this;
    }

    /** Checks for equality of this IncidentIterator and @a II_in. */
    bool operator==(const IncidentIterator& II_in) const {
      return Niter_==II_in.Niter_;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type Nid_;
    size_type Niter_;
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
    }

    // HW1 #5: YOUR CODE HERE

    /** Constructor for an EdgeIterator
      @param g_in Pointer to the current graph
      @param Niter_in A number representing the current position of the iterator
    */
    EdgeIterator(const Graph* g_in, size_type Niter_in) : graph_(const_cast<Graph*>(g_in)), Niter_(Niter_in) {}

    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

     /** Returns the Edge that this iterator is pointing to.
      @return An Edge object
    */
    Edge operator*() const {
      return Edge(graph_,Niter_);
    }

    /** Increments the iterator. 
      @return An EdgeIterator object
    */
    EdgeIterator& operator++() {
      Niter_++;
      return *this;
    }

    /** Checks for equality of this EdgeIterator and @a EI_in. */
    bool operator==(const EdgeIterator& EI_in) const {
      return Niter_==EI_in.Niter_;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type Niter_;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Returns an iterator pointing to the first edge in the graph.
  *  @return An EdgeIterator object
  */
  EdgeIterator edge_begin() const {
    return EdgeIterator(this,0);
  }

  /** Returns an iterator pointing to the last edge in the graph.
  *  @return An EdgeIterator object
  */
  EdgeIterator edge_end() const {
    return EdgeIterator(this,num_edges());
  }

 private:

  // HW0: YOUR CODE HERE
  size_type num_Nodes = 0;
  size_type num_Edges = 0;
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
