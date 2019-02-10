#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <tuple>

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
  using node_value_type = V;
  class Node : private totally_ordered<Node>{
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
    const Point& position() const {return Graph_->Positions_[id_];}

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {return this->id_;}


    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Return this node's value
      * @pre Graph_->has_node(*this)
      * @return value of the current node
      */
    node_value_type& value() {return Graph_->NodeValues_[id_];}

    /** Return this node's value
      * @pre Graph_->has_node(*this)
      * @return value of the current node
      */
    const node_value_type& value() const {return Graph_->NodeValues_[id_];} //#TODO why this?

    /** Return the degree of the current node (number of adjacent vertexes)
      * @pre Graph_->has_node(*this)
      * @return degree of current node
      */
    size_type degree() const {return Adjacency_[id_].size();};

    /** Set this node's value
      * @pre Graph_->has_node(*this)
      */
    void set_value(const node_value_type& value) {Graph_->NodeValues_[id_] = value;}

    /** Return the iterator of the first edge incident to a node
      * @pre Graph_->has_node(*this)
      * @post iterator to the first edge incident to a node
      */
    incident_iterator edge_begin() const {return IncidentIterator(Graph_, 0 ,id_);}

    /** Return the iterator of the last edge incident to a node
      * @pre Graph_->has_node(*this)
      * @post iterator to the last edge incident to a node
      */
    incident_iterator edge_end() const {return IncidentIterator(Graph_, Graph_->NodesToEdges_[id_].size() ,id_);}

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (this->Graph_ == n.Graph_ and this->id_ == n.id_) return true;
      else return false;
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
      return (id_ < n.id_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    graph_type* Graph_;
    size_type id_;
    Node(graph_type* Graph, size_type id) : Graph_(Graph), id_(id) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {return Positions_.size();}

  /** Synonym for size(). */
  size_type num_nodes() const {return size();}

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position,const node_value_type& Value = node_value_type()) {
    Positions_.push_back(position);
    std::vector<size_type> EmptyVecNTE;
    std::vector<size_type> EmptyVecADJ;

    NodesToEdges_.push_back(EmptyVecNTE);
    Adjacency_.push_back(EmptyVecADJ);
    NodeValues_.push_back(Value);

    // nNodes_++;
    return Node(this,num_nodes()-1); // zero based index
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return this == n.Graph_; // check they point to the same graph
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i <= num_nodes());
    return Node(const_cast<graph_type*>(this),i); // tell the compiler it's a constant
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
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {}

    /** Return a node of this Edge */
    Node node1() const {
      if (reversed_) return Node(Graph_,Graph_->EdgeEndNode_[id_]);
      else return Node(Graph_,Graph_->EdgeStartNode_[id_]);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      if (reversed_) return Node(Graph_,Graph_->EdgeStartNode_[id_]);
      else return Node(Graph_,Graph_->EdgeEndNode_[id_]);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (this->Graph_ == e.Graph_ and this->node1() == e.node1() and this->node2() == e.node2())
      return true;
      else if (this->Graph_ == e.Graph_ and this->node1() == e.node2() and this->node2() == e.node1())
      return true;
      else return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (this->id_ < e.id_);
    }

    void flip_orientation(){
      if (reversed_ == true ) reversed_ = false;
      if (reversed_ == false) reversed_ = true ;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    graph_type* Graph_; // pointer to the graph
    const size_type id_;  // unique id_ of the graph
    bool reversed_; // true if it should return node2() when calling in node1() and viceversa;

    // Default Constructor
    Edge(graph_type* Graph, size_type id, bool reversed = false)
      :  Graph_(Graph), id_(id), reversed_(reversed) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return this->EdgeStartNode_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(const_cast<graph_type*>(this),i);        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (size_type i = 0; i < Adjacency_[a.id_].size(); i += 1)
    {
      if (Adjacency_[a.id_][i] == b.id_) return true;
    }
    for (size_type i = 0; i < Adjacency_[b.id_].size(); i += 1)
    {
      if (Adjacency_[b.id_][i] == a.id_) return true;
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
    if (has_edge(a,b)){
      for (size_type i = 0; i < NodesToEdges_[a.id_].size(); i += 1)
      {
        for (size_type j = 0; j < NodesToEdges_[b.id_].size(); j += 1)
        {
          if (NodesToEdges_[a.id_][i] == NodesToEdges_[b.id_][j]) return Edge(const_cast<graph_type*>(this), NodesToEdges_[a.id_][i]);
        }
      }
      assert(1 != 0); // should not  reach this point
    }
    //else
    assert(a.Graph_ == b.Graph_);
    EdgeStartNode_.push_back(a.id_);
    EdgeEndNode_.push_back(b.id_);
    NodesToEdges_[a.id_].push_back(EdgeStartNode_.size() - 1);
    NodesToEdges_[b.id_].push_back(EdgeStartNode_.size() - 1);
    Adjacency_[a.id_].push_back(b.id_);
    Adjacency_[b.id_].push_back(a.id_);
    return Edge(const_cast<graph_type*>(a.Graph_),EdgeStartNode_.size() - 1); // Correction for 0 index
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    Positions_.clear();
    EdgeStartNode_.clear();
    EdgeEndNode_.clear();
    NodesToEdges_.clear();
    Adjacency_.clear();
    NodeValues_.clear();
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

    /** Dereference operator
      * @pre valid iterator referring to a valid node
      * @return Node at the position of the iterator
      */
    Node operator*() const {return Node(const_cast<graph_type*>(Graph_),id_); }

    /** Increment operator
      * @pre valid iterator
      * @return an iterator referring to the next node
      */
    NodeIterator& operator++() {id_ += 1; return *this;}

    /** Equality comparison operator
      * param[in] NodeIterator _x_
      * @pre valid iterators
      * @return boolean true if iterators refer to the same node index of the same graph, false otherwise
      */
    bool operator==(const NodeIterator& x) const {return (id_ == x.id_ and Graph_ == x.Graph_);}

    /** Inequality comparison operator
      * param[in] NodeIterator _x_
      * @pre valid iterators
      * @return boolean false if iterators refer to the same node index of the same graph, true otherwise
      */
    bool operator!=(const NodeIterator& x) const {return !(id_ == x.id_ and Graph_ == x.Graph_);}
   private:
    friend class Graph;
    // Default Constructor
    NodeIterator(graph_type* Graph, size_type id)
    :  Graph_(Graph), id_(id) {}
    const graph_type* Graph_; // Parent Graph
    size_type id_; // node id corresponding to this iterator
  };

  /** Return node iterator pointing to the first node
  * @return iterator to first node
  */
  node_iterator node_begin() const {return NodeIterator(const_cast<graph_type*>(this),0);}

  /** Return node iterator pointing to the last node
  * @return iterator to last node
  */
  node_iterator node_end() const {return NodeIterator(const_cast<graph_type*>(this),Positions_.size());}

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

    /** Dereference operator, returns the edge corresponding to current iterator
      * @return Edge at iterator's position
      * @post (Edge.node1() == node) spawning the incident edges the iterator is looping over
      */
    Edge operator*() const {
      // Extract the undirected edge
      Edge EdgeToReturn = Edge(const_cast<graph_type*>(Graph_),Graph_->NodesToEdges_[node_id_][list_id_]);
      // make sure its the spawning node is numbered first
      if (EdgeToReturn.node1().index() == node_id_) return EdgeToReturn;
      else {
        EdgeToReturn.flip_orientation(); // reverse the orientation of the edge
        assert(EdgeToReturn.node1().index() == node_id_); // make sure it actually worked
        return EdgeToReturn;
      }
    }

    /** Increment operator for the edge_iterator
      * @return Iterator with incremented position
      */
    IncidentIterator& operator++() {list_id_ += 1; return *this;}

    /** Equality comparison operator
      * @param[in] edge iterator
      * @return bool true if iterarors refer to same edge spawned by the same node of the same graph, false otherwise
      */
    bool operator==(const IncidentIterator& x) const {return (list_id_ == x.list_id_ and Graph_ == x.Graph_ and node_id_ == x.node_id_ );}

    /** Inequality comparison operator
      * @param[in] edge iterator
      * @return bool false if iterarors refer to same edge spawned by the same node of the same graph, true otherwise
      */
    bool operator!=(const IncidentIterator& x) const {return !(list_id_ == x.list_id_ and Graph_ == x.Graph_ and node_id_ == x.node_id_ );}

   private:
    friend class Graph;
    // Default constructor
    IncidentIterator(graph_type* Graph, size_type list_id_, size_type node_id)
    :  Graph_(Graph), list_id_(list_id_), node_id_(node_id) {}

    const graph_type* Graph_; // parent Graph
    size_type list_id_; // position of the element in the _NodesToEdges_ list
    const size_type node_id_; // Node spawning the adjiancent edges
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

    /** Dereference iterator for edge_iterator
      * @return edge pointed by the iterator
      */
    Edge operator*() const {return Edge(const_cast<graph_type*>(Graph_),id_); }

    /** Increment iterator for edge_iterator
      * @return EdgeIterator incremented by one position
      */
    EdgeIterator& operator++() {id_++; return *this;}

    /** Equality comparator between EdgeIterator
      * @param[in] edge iterator
      * @return bool true if @a e and @a this refer to the edge of the same graph, false otherwise
      */
    bool operator==(const EdgeIterator& e) const {return (id_ == e.id_ and Graph_ == e.Graph_);}
    bool operator!=(const EdgeIterator& e) const {return !(id_ == e.id_ and Graph_ == e.Graph_);}

   private:
    friend class Graph;
    const graph_type* Graph_; // Parent Graph
    size_type id_; // Edge ID
    // Default Constructor
    EdgeIterator(graph_type* g, size_type id)
    :  Graph_(g), id_(id) {}
  };

  /** Edge iterator at first edge of the graph
    * @return iterator at first edge
    */
  edge_iterator edge_begin() const {return EdgeIterator(const_cast<graph_type*>(this),0);}

  /** Edge iterator at last edge of the graph
    * @return iterator at last edge
    */
  edge_iterator edge_end()   const {return EdgeIterator(const_cast<graph_type*>(this),EdgeStartNode_.size());}

 private:

  std::vector<Point> Positions_; // Nodes Positions
  std::vector<size_type> EdgeStartNode_; // Nodes IDs Defining the Edges
  std::vector<size_type> EdgeEndNode_; // Nodes IDs Defining the Edges
  std::vector<std::vector<size_type>> NodesToEdges_;
  std::vector<std::vector<size_type>> Adjacency_;
  std::vector<node_value_type> NodeValues_;
};

#endif // CME212_GRAPH_HPP
