#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>

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
  //structures to store node & edge info in Graph 
  struct GraphNode;
  struct GraphEdge;

 public:
  //
  // PUBLIC TYPE DEFINITIONS
  //
  /** define node_value_type */ 
  using node_value_type = V ;
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
      _graph = nullptr;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      //return position of this Node's id
      return _graph->_nodes[_id].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return _id;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /**function return current node's value by reference 
     * we can set node's value through this function
     * @return node's value by reference 
     **/
    node_value_type& value() {
      return _graph->_nodes.at(_id).val;
    }
    /**function return current node's value by reference 
     * we can read node's value through this function
     * @return node's value by reference in constant 
     **/
    const node_value_type& value() const {
      return _graph->_nodes[_id].val;
    }
    /**function to get number of edges node connected to
     * @return number of of edges connected the current node
     */
    size_type degree() const {
      if (_graph->_adj_mat.find(_id) == _graph->_adj_mat.end()) return 0;
      return (size_type)_graph->_adj_mat.at(_id).size();
    }
    /** function return begin of the adjacency edge list
     *  @pre current node have adjacent edges
     *  @return incident_iterator pointing to start of adjacency edge hashmap
     */
    incident_iterator edge_begin() const {
      return incident_iterator(_graph->_adj_mat.at(_id).begin(), _graph);
    }
    /** function return end of the adjacency edge list
     *  @pre current node have adjacent edges
     *  @return incident_iterator pointing to end of adjacency edge hashmap
     */
    incident_iterator edge_end() const {
      return incident_iterator(_graph->_adj_mat.at(_id).end(), _graph);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return _id == n._id && _graph == n._graph;
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
      return _id < n._id;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    //private constuctor, support functions that return nodes in Graph class
    Node(size_type id, const graph_type * graph) : _id(id), _graph(const_cast<graph_type*>(graph)){}
    size_type _id; //id of this node
    graph_type *_graph; //pointer to Graph instance containing this node
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return (size_type) _nodes.size();
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
  Node add_node(const Point& position, const node_value_type& node_val = node_value_type()) {
    // HW0: YOUR CODE HERE
    //create a new GraphNode, add it to _nodes,
    //return new Node with correct id
    GraphNode new_gn;
    new_gn.position = position;
    new_gn.id = size();
    new_gn.val = node_val;
    _nodes.push_back(new_gn);
    return Node(new_gn.id, this);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    //check if node is in graoh by comparing n with node in graph with same id
    if (n._id < _nodes.size()) {
      if (n == node(n._id)) return true;
    }
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
    return Node(i, this); 
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
      _graph = nullptr;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return _graph->_edges[_id].a;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return _graph->_edges[_id].b;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //check if two nodes in e are corresponding to two nodes in current edge
      return (node1() == e.node1() && node2() == e.node2()) || \
              (node2() == e.node1() && node1() == e.node2());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //WARNING: We do not handle the case that two edges are from different graph
      // since the behavior is undefined
      return _id < e._id;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    size_type _id; //id of current edge
    graph_type *_graph; //pointer to graph containing this edge
    //private constructor for graoh member functions returning edge 
    Edge(size_type id, const Graph *graph) : _id(id), _graph(const_cast<graph_type*>(graph)){}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return (size_type) _edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(i, this);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    return find_edge(a, b) >= 0;
    //loop through edges and compare nodes in edge with a and b
    // for (auto& edge: _edges) {
    //   if ((edge.a == a && edge.b == b) || (edge.b == a && edge.a == b))
    //     return true;
    // }
    // return false;
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
    int edge_index = find_edge(a, b);
    //if return value >=0, return the found edge
    if (edge_index >= 0) {
      return Edge((size_type)edge_index, this);
    } else { //otherwise create a new GraphEdge and return an edge
      GraphEdge edge;
      edge.a = a;
      edge.b = b;
      size_type siz = (size_type)_edges.size();
      _edges.push_back(edge);
      _adj_mat[a._id].emplace(b._id, siz);
      _adj_mat[b._id].emplace(a._id, siz);
      return Edge(siz, this);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    _nodes.clear();
    _edges.clear();
    _adj_mat.clear();
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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** deference operator to return the node the iterator point to
     *  @pre _iter_index is a valid index
     *  @return Node object _iter_index refer
     */
    Node operator*() const {
      return Node(_iter_index, _graph);
    }
    /** increment operator of NodeIterator
     *  @post _iter_index is a valid index after increment
     *  @return updated NodeIterator by reference
     */
    NodeIterator& operator++() {
      ++_iter_index;
      return *this;
    }
    /** equality operator determine if two NodeIterators are the same
     *  @pre @param iter is a valid iterator
     *  @return the two iterators are same if they have same graph and index 
     */
    bool operator==(const NodeIterator& iter) const {
      return (_graph == iter._graph && _iter_index == iter._iter_index);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    //_iter_index keep track current position of iterator
    size_type _iter_index;
    graph_type *_graph;
    //NodeIterator constructor
    NodeIterator(size_type index, const graph_type* graph)
        : _iter_index(index), _graph(const_cast<graph_type*>(graph)){
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** function to return node_iterator pointing to the start of node list
   *  graph contains
   */
  node_iterator node_begin() const {
    return NodeIterator(0, this);
  }
  /** function to return node_iterator pointing to the end of node list
   */
  node_iterator node_end() const {
    return NodeIterator(_nodes.size(), this);
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>{
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
    /** deference operator of IncidentIterator
     *  @return adjacent edge the iterator currently points to
     *  @post node @param a of the edge is the node own's the current iterator
     *  @post node @param b of the edge is the adjacent node of the current node and edge  
     */
    Edge operator*() const {
      size_type edge_id = _map_iter->second;
      Edge edge = Edge(edge_id, _graph);
      //map_iter first should be dst
      //need to swap nodes if @param a is the adjacent node
      if (edge.node1()._id == _map_iter->first) {
        Node temp = _graph->_edges[edge_id].a;
        _graph->_edges[edge_id].a = _graph->_edges[edge_id].b;
        _graph->_edges[edge_id].b = temp;
      }
      return edge;
    }
    /** increment operator of IncidentIterator
     *  @post _iter_index is a valid index after increment
     *  @return updated IncidentIterator by reference
     */
    IncidentIterator& operator++() {
      ++_map_iter;
      return *this;
    }
    /** equality operator determine if two IncidentIterator are the same
     *  @pre @param inc_iter is a valid IncidentIterator
     *  @return the two iterators are same if they have same graph and index 
     */
    bool operator==(const IncidentIterator& inc_iter) const {
      if (_graph == inc_iter._graph && _map_iter == inc_iter._map_iter) return true;
      return false;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    unordered_map<size_type, size_type>::iterator _map_iter;
    Graph * _graph; 
    /**@pre iter points to the start of hashmap containing adjacent edges to current node
     */
    IncidentIterator(unordered_map<size_type, size_type>::iterator iter, const graph_type * graph)
      :_map_iter(iter), _graph(const_cast<graph_type*>(graph)) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>{
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
    /** deference operator of EdgeIterator
     *  @return edge the iterator currently points to
     */
    Edge operator*() const {
      return Edge(_iter_index, _graph);
    }
    /** increment operator of EdgeIterator
     *  @post _iter_index is a valid index after increment
     *  @return updated EdgeIterator by reference
     */
    EdgeIterator& operator++() {
      ++_iter_index;
      return *this;
    }
    /** equality operator determine if two EdgeIterator are the same
     *  @pre @param edge_iter is a valid EdgeIterator
     *  @return the two iterators are same if they have same graph and index 
     */
    bool operator==(const EdgeIterator& edge_iter) const {
      return (_graph == edge_iter._graph && _iter_index == edge_iter._iter_index);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    size_type _iter_index;
    graph_type* _graph;
    /**@pre 0 <= index <= size of graph's edges
     */
    EdgeIterator(size_type index, const graph_type* graph)
        : _iter_index(index), _graph(const_cast<graph_type*>(graph)){
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** function return the begin of the graph's edge list
   *  @pre graph have edges
   *  @return edge_iterator pointing to start of edge list
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(0, this);
  }
  /** function return the end of the graph's edge list
   *  @pre graph have edges
   *  @return edge_iterator pointing to the end of edge list
   */
  edge_iterator edge_end() const {
    return EdgeIterator(_edges.size(), this);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct GraphNode {
    Point position;
    size_type id;
    node_value_type val;
  };
  struct GraphEdge {
    Node a;
    Node b;
  };

  //using vectors to record nodes and edges in graph
  vector<GraphNode> _nodes;
  vector<GraphEdge> _edges;
  //node a -> (node b, edge id)
  /** keys of _adj_mat are nodes having adjacent edges
   *  value of _adj_mat are hash maps containing nodes adjacents to corresponding key 
   *  and the id of the edge connecting the two nodes
   */
  unordered_map<size_type, unordered_map<size_type, size_type>> _adj_mat;

  //internal helper function
  // @post edge id if find corresponding edge, otherwise -1
  int find_edge(const Node& a, const Node& b) {
    if (_adj_mat.find(a._id) == _adj_mat.end()) return -1;
    unordered_map<size_type, size_type>::const_iterator res = _adj_mat[a._id].find(b._id);
    if (res != _adj_mat[a._id].end()) return (int)res->second;
    return (int) -1;
  }


};

#endif // CME212_GRAPH_HPP
