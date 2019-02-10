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
  /** Type of node value. */
  using node_value_type = V;

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
    // HW0: YOUR CODE HERE (keep it empty)
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
      // HW0: YOUR CODE HERE (DONE)
      this->graph = nullptr;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE (DONE)
      return this->graph->coords[this->nid];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE (DONE)
      return this->nid;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Get value of this node. Return by reference. */
    node_value_type& value() {
        return this->v;
    }

    /** Get read-only value of this node. */
    const node_value_type& value() const {
        return this->v;
    }

    /** Get number of incident edges. */
    size_type degree() const {return this->graph->n2n[this->nid].size();}

    /** Get begin point of incident_iterator. */
    incident_iterator edge_begin() const {return incident_iterator(*this, 0);}

    /** Get end point of edge_itorator. */
    incident_iterator edge_end() const {
        return incident_iterator(*this, this->degree());
    }



    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE (DONE)
      return this->graph == n.graph && this->nid == n.nid;
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
      // HW0: YOUR CODE HERE (DONE)
      return (this->graph == n.graph ? this->nid < n.nid : this->graph < n.graph);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    size_type nid;   // node index
    Graph* graph;    // graph pointer
    node_value_type v; // value
  };


 
 public:

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE (DONE)
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
  Node add_node(const Point& position, const node_value_type& = node_value_type()) {
    // HW0: YOUR CODE HERE (DONE)
    Node n;
    coords.push_back(position);
    n.graph = this;
    n.nid = nodes.size();
    nodes.push_back(n);
    n2n.push_back(std::vector<size_type>());
    return n;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE (DONE)
    return (this == n.graph);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node& node(size_type i) {
    // HW0: YOUR CODE HERE (DONE)
    return nodes[i];
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
      // HW0: YOUR CODE HERE (DONE)
      n1 = nullptr;
      n2 = nullptr;
    }

    /** Construct a valid Edge. */
    Edge(const Node& node_1, const Node& node_2) : n1(&node_1), n2(&node_2) {} 

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE (DONE)
      if (n1 != nullptr) {return *n1;}
      return Node();
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE (DONE)
      if (n2 != nullptr) {return *n2;}
      return Node();
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const { // (Not sure I should do this in HW0)
      return (e.n1==this->n1 && e.n2==this->n2) || (e.n1==this->n2 && e.n2==this->n1);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return this < &e;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    const Node* n1;
    const Node* n2;
    Graph* graph;
  };

 public:

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE (DONE)
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE (DONE)
    return edges[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE (DONE)
    size_type ia = a.nid;
    size_type ib = b.nid;
    const std::vector<size_type>& dict = this->n2n[ia];
    for (size_type i=0; i<dict.size(); i++) {
        if (dict[i] == ib) {return true;}
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
    // HW0: YOUR CODE HERE (DONE)
    Edge e;
    if (!this->has_edge(a, b)) {
        e.graph = this;
        e.n1 = &a;
        e.n2 = &b;
        edges.push_back(e);
        n2n[a.nid].push_back(b.nid);
        n2n[b.nid].push_back(a.nid);
    }
    return e;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE (DONE)
    n2n.clear();
    edges.clear();
    nodes.clear();
    coords.clear();
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

    typedef std::vector<value_type> NodeContainer;      // Node container used in Graph class

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    /** Construct an valide NodeIterator. */
    NodeIterator(const Graph* g, size_type id=0) : node_container(g->nodes), nid(id) {}


    // HW1 #2: YOUR CODE HERE (DONE)
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** De-reference Node object. (required for an iterator)*/
    value_type operator*() const {return node_container[nid];}

    /** Forward move to the next node address. (required for an iterator)*/
    NodeIterator& operator++() {
        nid++;
        return *this;
    }
    
    /** Check equality. (required for an iterator)*/
    bool operator==(const NodeIterator& other) const {
        return this->nid == other.nid;
    }

    /** Overload copy assignment operator. */
    NodeIterator& operator=(const NodeIterator&) {}

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    const NodeContainer& node_container;
    size_type nid;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Get starting point of node_iterator. */
  node_iterator node_begin() const {return node_iterator(this, 0);}

  /** Get end point of node_iterator. */
  node_iterator node_end() const {return node_iterator(this, nodes.size());}

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

    typedef std::vector<size_type> Dict;
    typedef std::vector<Node> NodeContainer;
    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {}


    /** Construct a valid IncidentIterator. */
    IncidentIterator(const Node& n, size_type r=0) :
        n1(n), rank(r), dict(n.graph->n2n[n.nid]), all_nodes(n.graph->nodes) {}

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const


    /** De-reference Edge. */
    Edge operator*() const {
        return Edge(n1, all_nodes[dict[rank]]);
    }

    /** Increament of iterator. */
    IncidentIterator& operator++() {
        rank++;
        return *this;
    }


    /** Check equality.*/
    bool operator==(const IncidentIterator& other) const {
        return this->n1 == other.n1 && this->rank == other.rank;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    const Node& n1;
    size_type rank;   // An index of an incident node relative to this node.
    const Dict& dict; // reference vector that stores the connected nodes.
    const NodeContainer& all_nodes; // reference vector of all nodes stored in graph.
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
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /** Construct an valid node iterator. */
    EdgeIterator(const Graph& g, const size_type i) : graph(g), eid(i) {}

    /** De-reference an edge iterator. */
    Edge operator*() const {
      //Node& n1 = graph.node(graph.edges[eid].node1.index());
      //Node& n2 = graph.node(graph.edges[eid].node2.index());
      return graph.edges[eid];
    }

    /** Forward moving. */
    EdgeIterator& operator++() {
        eid++;
        return *this;
    }

    /** Compare equality. */
    bool operator==(const EdgeIterator& other) const {
        return this->eid == other.eid && &(this->graph)==&(other.graph);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const Graph& graph;
    size_type eid;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Return the start point of the edge iterator. */
  edge_iterator edge_begin() const {
    return edge_iterator(*this, 0);
  }

  /** Return the end point of the edge iterator. */
  edge_iterator edge_end() const {
    return edge_iterator(*this, this->edges.size());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::vector<Node> nodes; // use std::vector to contain each Node
  std::vector<Edge> edges; // use std::vector to contain each edge
  std::vector<std::vector<size_type>> n2n; // map of connected nodes assuming edges cannot cross graphs
  std::vector<Point> coords; // positions
};

#endif // CME212_GRAPH_HPP
