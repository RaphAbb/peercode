#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <set> // used in HW0
#include <vector>
#include <map> // HW2

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
 public:
  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V,E>;

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
  using node_value_type = V; // type of objects held in each node
  using edge_value_type = E; // type of objects held in each edge

 private:
  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  /* Explanation of below data structures
     Nodes:
     - node_to_position_map: (key,value) = (node index, Point object)
     - node_to_value_map: (key,value) = (node index, node_value_type object)
     - node_to_incident_map: (key,value) = (node index, set of edge indexes)
     Edges:
     - edge_to_nodeVec_map: (key,value) = (edge index, vector of node indexes)
     - nodeSet_to_edge_map: (key,value) = (set of node indexes, edge index)
     - edge_to_value_map: (key,value) = (edge index, edge_value_type object)
  */

  // Node data containers
  std::map<size_type, Point> node_to_position_map;
  std::map<size_type, node_value_type> node_to_value_map;
  std::map<size_type, std::set<size_type>> node_to_incident_map;

  // Edge data containers
  std::map<size_type, std::vector<size_type>> edge_to_nodeVec_map;  // used by edge1() and edge2()
  std::map<std::set<size_type>, size_type> nodeSet_to_edge_map;  // used by has_edge()
  std::map<size_type, edge_value_type> edge_to_value_map;

  size_type number_nodes; // number of nodes
  size_type number_edges; // number of edges

 public:

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    number_nodes = 0;
    number_edges = 0;
  }

  /** Default destructor */
  ~Graph() = default;

  // ====================================================================
  // NODES
  // ====================================================================

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node : private totally_ordered<Node> {
   public:
    Graph *owner_graph; // pointer to a constant object of type Graph
    size_type node_id;  // index of node

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

    /** Return this node's position, a Point object */
    Point& position() {
      // HW0: YOUR CODE HERE
      return (owner_graph->node_to_position_map).at(node_id);
    }

    /** Return this node's position, a Point object */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return (owner_graph->node_to_position_map).at(node_id);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
     return node_id;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Return the value stored in this node */
    node_value_type& value() {
      return (owner_graph->node_to_value_map).at(node_id);
    }

    /** Return the value stored in this node */
    const node_value_type& value() const {
      return (owner_graph->node_to_value_map).at(node_id);
    }

    /* Return the number of edges connected to this node */
    size_type degree() const {
      return (owner_graph->node_to_incident_map).at(node_id).size();
      // be carefule about type conversion here
      // set::size returns an unsigned int
    }

    /* Returns the beginning iterator which will iterate over all edges 
     * connected to this node */
    incident_iterator edge_begin() const {
      auto set_iter = owner_graph->node_to_incident_map.at(node_id).begin();
      return IncidentIterator(owner_graph, node_id, set_iter);
    }

    /* Returns the ending iterator which will iterate over all edges 
     * connected to this node */
    incident_iterator edge_end() const {
      auto set_iter = owner_graph->node_to_incident_map.at(node_id).end();
      return IncidentIterator(owner_graph, node_id, set_iter);
    }


    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
     if ((owner_graph == n.owner_graph) && (node_id == n.node_id)) {return true;}
     else {return false;}
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
      // I'm assuming we'll only compare nodes that are in same graph
     if (node_id < n.node_id) {return true;}
     else {return false;}
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    Node(const Graph *graph, size_type id) : node_id(id) { 
      owner_graph = const_cast<Graph*>(graph);  // allow nodes to modify graph's data
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return number_nodes;
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
  Node add_node(const Point& position, 
                const node_value_type& value = node_value_type()) {
    // Note: index of new node is current value of number_nodes
    node_to_position_map.insert({number_nodes, position});
    node_to_value_map.insert({number_nodes, value});
    node_to_incident_map.insert({number_nodes, std::set<size_type>()}); // add empty set
    number_nodes++; // just made new node, so increment number of nodes
    return node(number_nodes-1);  // -1 to get back to node just created
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    /* For graph to have node n, two things must be true:
        (1) n is associated with current graph instance
        (2) graph instance must contain a node with same index as n */
    if ((n.owner_graph == this) && (n.node_id < number_nodes)) {return true;}
    else {return false;}
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    // Ensure function can only return nodes that already exist in graph
    assert( (0 <= i) && (i < number_nodes) );
    return Node(this, i);
  }

  // ====================================================================
  // EDGES
  // ====================================================================

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Edge> {
   public:
    Graph *owner_graph; // pointer to a constant object of type Graph
    size_type edge_id;  // index of edge

    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
     return owner_graph->node((owner_graph->edge_to_nodeVec_map).at(edge_id)[0]);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
     return owner_graph->node((owner_graph->edge_to_nodeVec_map).at(edge_id)[1]);
    }

    /* Return the length of the edge, which equals the Euclidian distance
       between the two nodes to which the edge is incident */
    double length() const {
      return norm(this->node1().position() - this->node2().position());
    }

    /** Return the value stored in this edge */
    edge_value_type& value() {
      return (owner_graph->edge_to_value_map).at(edge_id);
    }

    /** Return the value stored in this edge */
    const edge_value_type& value() const {
      return (owner_graph->edge_to_value_map).at(edge_id);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
     if ( (owner_graph == e.owner_graph) && (edge_id == e.edge_id) ) {return true;}
     else {return false;}
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
     if (owner_graph != e.owner_graph) {
      if (owner_graph < e.owner_graph) {return false;}
      else {return true;}
     }
     if (edge_id < e.edge_id) {return true;}
     else {return false;}
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Constructor that will be used by Graph::add_edge()
    // Edge(const Graph *graph, size_type id)
    //     : owner_graph(graph), edge_id(id) { }
    Edge(const Graph *graph, size_type id) : edge_id(id) {
      owner_graph = const_cast<Graph*>(graph);  // allow edges to modify graph's data
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return number_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    // Ensure function can only return edges that already exist in graph
    assert( (0 <= i) && (i < number_edges) );
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
    std::set<size_type> node_set = {a.index(), b.index()};
    // Check if node_set is in nodeSet_to_edge_map
    if (nodeSet_to_edge_map.find(node_set) == nodeSet_to_edge_map.end()) {return false;}
    else {return true;}
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
    /*
    Give this function the ability to flip the order of elements
    in edge_to_node_vec_map. If edge (a,b) does not exist yet, make
    edge_to_node_vec_map = [a,b]. If you then call add_edge(b,a), has_edge()
    will return true, but now flip edge_to_node_vec_map from [a,b] to [b,a]

    TLDR: you can use this function to return an existing edge but with the
    orientation encoded in edge_to_node_vec_map flipped
    */
    // HW0: YOUR CODE HERE
    assert(this->has_node(a) && this->has_node(b)); // ensure nodes are valid
    assert(a != b); // ensure nodes are distinct
    size_type a_id = a.index();
    size_type b_id = b.index();
    std::set<size_type> node_set = {a_id, b_id};
    std::vector<size_type> node_vec = {a_id, b_id};
    if (has_edge(a,b) == false)
    {
     size_type new_edge_id = number_edges;
     nodeSet_to_edge_map.insert({node_set, new_edge_id});
     edge_to_nodeVec_map.insert({new_edge_id, node_vec});
     // default initialize value of edges based on edge_value_type
     edge_to_value_map.insert({new_edge_id, edge_value_type()});
     // Add incident edges
     node_to_incident_map[a_id].insert({new_edge_id});
     node_to_incident_map[b_id].insert({new_edge_id});
     number_edges++;  // just made new edge, so increment number of edges
     return edge(new_edge_id);
    }
    else
    {
      size_type edge_id = nodeSet_to_edge_map.at(node_set);
      // If node pair is NOT stored in correct order, flip the order
      if ( edge_to_nodeVec_map.at(edge_id).at(0) != a_id ) {
        edge_to_nodeVec_map.at(edge_id).at(0) = a_id;
        edge_to_nodeVec_map.at(edge_id).at(1) = b_id;
      }
      return edge(edge_id);
    }
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    // Empty all maps
    node_to_value_map.clear();
    node_to_position_map.clear();
    edge_to_nodeVec_map.clear();
    nodeSet_to_edge_map.clear();
    node_to_incident_map.clear();
    number_nodes = 0;
    number_edges = 0;
  }


  // ====================================================================
  // NODE ITERATOR
  // ====================================================================

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
   public:
    Graph *owner_graph;
    size_type index;

    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }
    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Return the node object corresponding to current NodeIterator object */
    Node operator*() const { 
      return owner_graph->node(index);
    }

    /** Increment index of NodeIterator object
     * @pre this.index != node_end()
     * @post new index = old index + 1
     */
    node_iterator& operator++() {
      index++;
      return *this;
    }

    /** Compare two NodeIterator objects for equality
     * @pre Both NodeIterator objects must belong to the same graph instance
     */
    bool operator==(const node_iterator& other_iterator) const {
      assert(owner_graph == other_iterator.owner_graph);
      if (index == other_iterator.index){return true;}
      else {return false;}
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    // Add private constructor that a Graph instance will use.
    // The node_begin() and node_end() functions will call this constructor
    // and return a NodeIterator object
    NodeIterator(const Graph *graph, size_type i) {
      owner_graph = const_cast<Graph*>(graph);
      index = i;
    }
  };

  // HW1 #3: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Returns the node iterator object corresponding to the
   * first node in the graph
   * @post 0 <= @a i < number_nodes, where @a i is 
   *     the index of returned NodeIterator object
   */
  node_iterator node_begin() const {
    return node_iterator(this, 0); // first node has index 0
  }

  /** Returns the node iterator object corresponding to on past the
   * last node in the graph
   * @post @a i == number_nodes, where @a i is 
   *     the index of returned NodeIterator object
   */
  node_iterator node_end() const {
    return node_iterator(this, number_nodes); 
    // last node has index number_nodes-1
    // ==> one past last node has index number_nodes
  }


  // ====================================================================
  // INCIDENT ITERATOR
  // ====================================================================

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
   public:
    Graph *owner_graph;
    size_type node_id;
    std::set<size_type>::iterator set_iterator;  // iterator for set of edge indexes

    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    // HW1 #4: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** Return edge object corresponding to current IncidentIterator object */
    Edge operator*() const {
      size_type edge_id = *set_iterator;
      size_type one = owner_graph->edge_to_nodeVec_map.at(edge_id).at(0);
      size_type two = owner_graph->edge_to_nodeVec_map.at(edge_id).at(1);
      // If edge.node1() = node_id, as we want it to
      if (one == node_id){ return owner_graph->edge(edge_id);}
      // If edge.node1() != node_id, use add_edge to flip nodes
      else { return owner_graph->add_edge( owner_graph->node(two), 
                                           owner_graph->node(one));}
    }

    /** Increment index of IncidentIterator object
     * @pre this.set_iterator != this.set_iterator.end()
     */
    IncidentIterator& operator++() {
      set_iterator++;
      return *this;
    }

    /** Compare two IncidentIterator objects for equality
     * @pre Both IncidentIterator objects must belong to same graph instance
     */
    bool operator==(const IncidentIterator& other_iterator) const {
      assert(owner_graph == other_iterator.owner_graph);
      if (set_iterator == other_iterator.set_iterator) {return true;}
      else {return false;}
    }

   private:
    friend class Graph;
    // HW1 #4: YOUR CODE HERE
    IncidentIterator(Graph *graph, 
                     size_type node, 
                     std::set<size_type>::iterator iter)
      : owner_graph(graph), node_id(node), set_iterator(iter) {}

  };

  // ====================================================================
  // EDGE ITERATOR
  // ====================================================================

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
   public:
    Graph *owner_graph;
    size_type index;

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

    /** Return edge object corresponding to current EdgeIterator object */
    Edge operator*() const {
      return owner_graph->edge(index);
    }

    /** Increment index of EdgeIterator object
     * @pre this.index != edge_end()
     * @post new index = old index + 1
     */
    EdgeIterator& operator++() {
      index++;
      return *this;
    }

    /** Compare two EdgeIterator objects for equality
     * @pre Both EdgeIterator objects must belong to the same graph instance
     */
    bool operator==(const EdgeIterator& other_iterator) const {
      assert(owner_graph == other_iterator.owner_graph);
      if (index == other_iterator.index) {return true;}
      else {return false;}
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    // EdgeIterator(const Graph *graph, size_type i)
    //   : owner_graph(graph), index(i) { }
    EdgeIterator(const Graph *graph, size_type i) {
      owner_graph = const_cast<Graph*>(graph);
      index = i;
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Returns the EdgeIterator object corresponding to the
   * first edge in the graph
   * @post 0 <= @a i < number_edges, where @a i is 
   *     the index of returned EdgeIterator object
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this,0);
  }

  /** Returns the EdgeIterator object corresponding to one
   * past the last edge in the graph
   * @post @a i == number_edges, where @a i is 
   *     the index of returned EdgeIterator object
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this,number_edges);
    // last edge has index number_edges-1
    // ==> one past last node has index number_edges
  }


  size_type remove_edge(const Edge& e) {
    // Edge to remove from graph
    // -----------------------------------------------
    size_type id_prime = e.edge_id;  // id of edge to remove
    Node node_1_prime = e.node1();
    Node node_2_prime = e.node2();
    std::set<size_type> incident_set_prime = {node_1_prime.index(), node_2_prime.index()};
    if (has_edge(node_1_prime, node_2_prime) == false) {return 0;} // Nothing to do if graph doesn't have edge

    // Edge with last index, i.e. edge_id equal to number_edges - 1
    // -----------------------------------------------
    size_type id_last = number_edges-1;
    Edge e_last = this->edge(id_last);
    Node node_1_last = e_last.node1();
    Node node_2_last = e_last.node2();
    std::set<size_type> incident_set_last = {node_1_last.index(), node_2_last.index()};
    std::vector<size_type> incident_vec_last = {node_1_last.index(), node_2_last.index()};


    // Graph modifications
    // -----------------------------------------------
    // Remove edge
    nodeSet_to_edge_map.erase(incident_set_prime);
    node_to_incident_map.at(node_1_prime.index()).erase(id_prime);
    node_to_incident_map.at(node_2_prime.index()).erase(id_prime);

    if (id_prime != id_last) {
      // Only execute this block if edge to remove isn't the last edge
      // Reassign data corresponding to last edge
      nodeSet_to_edge_map.at(incident_set_last) = id_prime;
      node_to_incident_map.at(node_1_last.index()).erase(id_last);
      node_to_incident_map.at(node_1_last.index()).insert(id_prime);
      node_to_incident_map.at(node_2_last.index()).erase(id_last);
      node_to_incident_map.at(node_2_last.index()).insert(id_prime);

      // Reassign values of edge to "remove"
      edge_to_nodeVec_map.at(id_prime) = incident_vec_last;
      edge_to_value_map.at(id_prime) = e_last.value();
    }
    // Below code will execute as intended if (id_prime == id_last) or (id_prime != id_last)
    edge_to_nodeVec_map.erase(id_last);
    edge_to_value_map.erase(id_last);

    number_edges--; // decrement number of edges by 1

    return number_edges;
  }


  size_type remove_node(const Node& n) {
    // Node to remove from graph
    // ---------------------------------------------------------------
    size_type id_prime = n.index();
    if (has_node(n) == false) {return 0;} // if graph doesn't have node, there's nothing to do

    // Node with last index, i.e. with node_id equal to number_nodes-1
    // ---------------------------------------------------------------
    size_type id_last = number_nodes-1;
    Node n_last = this->node(id_last);

    // Remove all edges incident to node n
    auto it = n.edge_begin();
    while (it != n.edge_end()) {
      this->remove_edge(*it);
      it = n.edge_begin();
    }


    if (id_prime != id_last) {

      node_to_incident_map.at(id_prime).clear(); // should already be empty but just to make sure

      // Loop through edges incident to last node
      for (auto it = n_last.edge_begin(); it != n_last.edge_end(); ++it) {
        Edge e = *it;
        Node node_b = e.node2();

        std::vector<size_type> old_incident_vec = {id_last, node_b.index()};
        std::vector<size_type> new_incident_vec = {id_prime, node_b.index()};

        std::set<size_type> old_incident_set = {id_last, node_b.index()};
        std::set<size_type> new_incident_set = {id_prime, node_b.index()};

        // Graph modifications
        edge_to_nodeVec_map.at(e.edge_id) = new_incident_vec;
        nodeSet_to_edge_map.erase(old_incident_set);
        nodeSet_to_edge_map.insert({new_incident_set,e.edge_id});
        node_to_incident_map.at(id_prime).insert({e.edge_id});
      }

      node_to_position_map.at(id_prime) = n_last.position();
      node_to_value_map.at(id_prime) = n_last.value();

    }

    // Below code will execute as intended if (id_prime == id_last) or (id_prime != id_last)
    node_to_position_map.erase(id_last);
    node_to_value_map.erase(id_last);
    node_to_incident_map.erase(id_last);

    number_nodes--; // decrement number of nodes by 1
    return number_nodes;
  }

  node_iterator remove_node(node_iterator n_it) {
    Node n = *n_it;
    this->remove_node(n);
    return n_it;
  }

  size_type remove_edge(const Node& node_a, const Node& node_b) {
    if (has_edge(node_a, node_b) == false) {return 0;} // Nothing to do if graph doesn't have edge
    std::set<size_type> node_set = {node_a.index(), node_b.index()};
    size_type e_id = nodeSet_to_edge_map.at(node_set);
    Edge e = this->edge(e_id);
    return this->remove_edge(e);
  }

  edge_iterator remove_edge(edge_iterator e_it) {
    Edge e = *e_it;
    this->remove_edge(e);
    return e_it;
  }

};

#endif // CME212_GRAPH_HPP





