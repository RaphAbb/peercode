#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */
#include <algorithm>
#include <vector>
#include <cassert>
#include <utility>
#include <map>

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

 private:
  //
  // Internal type declarations using public aliased size types
  //
  /** Internal node representation
   * @details Contains node's data and proxy id
  */
  struct NodeInternal{
    Point point;
    size_type id;
    node_value_type value;
    size_type degree;
    NodeInternal(Point p, size_type i, node_value_type v)
      : point{p}, id{i}, value{v}, degree{0} {}
  };
  /** Internal edge representation
   * @details Contains internal node indices and this edge's proxy id
  */
  struct EdgeInternal {
    size_type node_id_first;
    size_type node_id_second;
    size_type id;
    edge_value_type value;
    /** Constructs internal node using first and second indices. */
    EdgeInternal(size_type i_first,
                 size_type i_second,
                 size_type edge_id) : node_id_first{i_first},
                                      node_id_second{i_second},
                                      id{edge_id} {}
  };

  /** Private collection of nodes (containing ground truth)
   * Size is at least the number of valid nodes in graph
   */
  std::vector<NodeInternal> nodes_internal;
  /** Edge container (containing ground truth) for O(1) random access. */
  std::vector<EdgeInternal> edges_internal;
  /** Edge container for O(log(N)) edge checking (via RBT backing).
   * @details For some edge from tail to head, use the head node index as key
   * to return a <key, map> pair. Then use the tail node index as key on
   * to return a <key, edge_index> pair
  */
  using inner_map_type = std::map<size_type, size_type>;
  std::map<size_type, inner_map_type> adjacency_map;
  /** Index redirection layer
   * @details Vectors whose indices map to the indices of nodes_internal and
   *   edges_internal, i.e., lookups using proxy indices return the internal
   *   index.
  */
  std::vector<size_type> node_index_mapping;
  std::vector<size_type> edge_index_mapping;

 public:

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() { }

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
    Node() { }

    /** Return this node's position. */
    const Point& position() const {
      return parent_graph->point_at(proxy_id);
    }

    /** Return this node's position. */
    Point& position() {
      return parent_graph->point_at(proxy_id);
    }

    /** Return this node's proxy index
     * @pre This node is valid
     * @return A proxy node index in the range [0, graph_size).
     */
    size_type index() const {
      return proxy_id;
    }

    /** Return mutable reference to node's value.
     * @return mutable reference to node's internal value
     * @pre 0 <= @a this->proxy_id < num_nodes
     */
    node_value_type& value() {
      return parent_graph->value_at(proxy_id);
    }

    /** Return immutable reference to node's value.
     * @return immutable reference to node's internal value
     * @pre node is assigned a valid value
     * @pre 0 <= @a this->proxy_id < num_nodes
     * @post internal value of node is not changed
     */
    const node_value_type& value() const {
      // return parent_graph->valueAt(proxy_id);
      return parent_graph->value_at(proxy_id);
    }

    /** Queries degree of this node using parent graph's degree.
     * @return degree as recorded by graph's internal representation
     * @pre 0 <= @a this->proxy_id < num_nodes
     */
    size_type degree() const {
      return parent_graph->degree_at(proxy_id);
    }

    /** Returns iterator pointing to start of collection of incident edges.
     * @return Beginning iterator of the container at its current state
     * @pre 0 <= @a this->proxy_id < num_nodes
     * @post result is valid until deletion of nodes or edges
     * @details Method is agnostic to internal data representation
     *   in the parent graph.
     */
    incident_iterator edge_begin() const {
      // Delegate data representation matters to parent graph
      return parent_graph->incident_edge_begin(proxy_id);
    }

    /** Returns iterator pointing to end of collection of incident edges.
     * @return End iterator of the container at its current state
     * @pre 0 <= @a this->proxy_id < num_nodes
     * @post result is valid until insertion or deletion of nodes or edges
     * @details Method is agnostic to internal data representation
     *   in the parent graph.
     */
    incident_iterator edge_end() const {
      // Delegate data representation matters to parent graph
      return parent_graph->incident_edge_end(proxy_id);
    }

    /** Test whether this node and @a n are equal.
     * @param[in] node to be compared to
     * 
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return parent_graph == n.parent_graph &&
             proxy_id == n.proxy_id;
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
      return parent_graph != n.parent_graph ?
               parent_graph < n.parent_graph :
               proxy_id < n.proxy_id;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Id known to the proxy (guaranteed == Id in Graph)
    size_type proxy_id;
    // Pointer to graph (allows reassignment for invalid nodes but not mutate)
    const Graph* parent_graph;

    /** Parametrized constructor for node proxy exposed to friend Graph class.
     * @brief Constructs fully specified node using node id and graph pointer
     */
    Node(size_type new_id, const Graph* graph_ptr) : proxy_id{new_id},
                                                     parent_graph{graph_ptr} {}

    /** Reassigns properties of proxy node
     * @brief Overwrite proxy id and pointer to graph pointer. 
    */
    void set_new(size_type new_id, Graph* graph_ptr){
      proxy_id = new_id;
      parent_graph = graph_ptr;
    }
  };

  /** Return the number of nodes in the graph.
   * @return Number of valid accessible nodes
   * 
   * Complexity: O(1).
   * @details Returns number of valid nodes in the container, and not the
   * actual quantity of data written, but it is guaranteed that the number
   * of valid nodes <= number of data written.
   */
  size_type size() const {
    return node_index_mapping.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position,
                const node_value_type& value=node_value_type{}) {
    // Create new node in random access structure
    auto stack_top = node_index_mapping.size();
    nodes_internal.emplace_back(position, stack_top, value);
    // Create empty edge structure in map corresponding to internal index
    adjacency_map.insert(std::pair<size_type, inner_map_type> (
      nodes_internal.size()-1, inner_map_type{}));
    // Add new mapping from new proxy id to new internal id
    node_index_mapping.push_back(nodes_internal.size()-1);
    // Return node with id equal to end of node vector
    return Node(stack_top, this);
  }

  /** Determine if a Node belongs to this Graph
   * @return True iff @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return has_node_index(n.proxy_id);
  }

  /** Return the node with index @a i. If none, returns invalid node.
   * @pre 0 <= @a proxy_id < num_nodes()
   * @return Node with result.index() == i referring to this graph
   *
   * Complexity: O(1).
   */
  Node node(size_type proxy_id) const {
    if (has_node_index(proxy_id)){
      return Node(proxy_id, this);
    }
    return Node(); // Invalid node
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
    Edge() { }

    /** Return first node of this Edge.
     * @pre 0 <= @a this->edge_proxy_id < num_edges
     * @return The source node of edge if directed information was provided
     * at construction. Otherwise, returns first node in internal
     * (arbitrarily ordered) representation. Self-loops ok.
     */
    Node node1() const {
      // Query graph for first internal node index
      auto node_id = parent_graph->node_int_of_edge(edge_proxy_id, 0);
      // Return node if it is the source, else return the other
      if (node_id == source_node_internal_id)
        return Node(parent_graph->node_proxy_at_internal(node_id),
                    parent_graph);
      else {
        return Node(parent_graph->node_proxy_at_internal(
          parent_graph->node_int_of_edge(edge_proxy_id, 1)),
          parent_graph);
      }
    }

    /** Return second node of this Edge.
     * @pre 0 <= @a this->edge_proxy_id < num_edges
     * @return The non-source node of edge if directed information was provided
     * at construction. Otherwise, returns first node in internal
     * (arbitrarily ordered) representation. Self-loops ok.
     */
    Node node2() const {
      // Query graph for first node of internal representation
      auto node_id = parent_graph->node_int_of_edge(edge_proxy_id, 0);
      // Return node if it is not the source, else return the other
      if (node_id != source_node_internal_id)
        return Node(parent_graph->node_proxy_at_internal(node_id),
                    parent_graph);
      else {
        return Node(parent_graph->node_proxy_at_internal(
          parent_graph->node_int_of_edge(edge_proxy_id, 1)),
          parent_graph);
      }
    }

    /** Return mutable reference to edge's value.
     * @return mutable reference to edge's internal value
     * @pre 0 <= @a this->edge_proxy_id < num_edges
     */
    edge_value_type& value() {
      return parent_graph->value_edge(edge_proxy_id);
    }

    /** Return immutable reference to edge's value.
     * @return immutable reference to edge's internal value
     * @pre edge is assigned a valid value
     * @pre 0 <= @a this->edge_proxy_id < num_edges
     * @post internal value of edge is not changed
     */
    const edge_value_type& value() const {
      return parent_graph->value_edge(edge_proxy_id);
    }

    /** Test whether this edge and @a e are equal.
     * @param[in] e Other edge to compare to
     * 
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // Unpack edge node indices
      size_type e1_n1, e1_n2, e2_n1, e2_n2;
      e1_n1 = parent_graph->node_int_of_edge(edge_proxy_id, 0);
      e1_n2 = parent_graph->node_int_of_edge(edge_proxy_id, 1);
      e2_n1 = parent_graph->node_int_of_edge(e.edge_proxy_id, 0);
      e2_n2 = parent_graph->node_int_of_edge(e.edge_proxy_id, 1);
      // Return equality, exhausting combinatorics of nodes
      return parent_graph == e.parent_graph && (
        (e1_n1 == e2_n1 && e1_n2 == e2_n2) ||
        (e1_n1 == e2_n2 && e1_n2 == e2_n1));
    }

    /** Test whether this edge is less than @a e in a global order.
     * @param[in] e Other edge to compare to
     * 
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // Unpack edge node indices
      size_type e1_n1, e1_n2, e2_n1, e2_n2;
      e1_n1 = parent_graph->node_int_of_edge(edge_proxy_id, 0);
      e1_n2 = parent_graph->node_int_of_edge(edge_proxy_id, 1);
      e2_n1 = parent_graph->node_int_of_edge(e.edge_proxy_id, 0);
      e2_n2 = parent_graph->node_int_of_edge(e.edge_proxy_id, 1);
      // Return lexicographic order
      return parent_graph != e.parent_graph ? parent_graph < e.parent_graph :
        (e1_n1 != e2_n1 ? e1_n1 < e2_n1 :
          (e1_n2 != e2_n2 ? e1_n2 < e2_n2 : false));
    }

    /** Return length of this edge in physical space
     * @return Euclidean distance between end nodes
     * @pre @a this is a valid edge
     */
    double length() const {
      return norm(node1().position() - node2().position());
    }

    /** Return proxy index of this edge **/
    size_type index() const {
      return edge_proxy_id;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    
    // Edge Id known to the proxy (guaranteed == ID in Graph if valid in Graph)
    size_type edge_proxy_id;
    // Pointer to graph (allows reassignment for invalid nodes but not mutate)
    const Graph* parent_graph;
    // Node Id of source node (in context of incident iterator)
    size_type source_node_internal_id;

    /** Param constructor for directed edge proxy exposed to friend Graph.
     * @brief Constructs fully specified edge using edge id, graph pointer and
     *   an internal node id for the source of the directed edge
     */
    Edge(size_type edge_id,
         const Graph* graph_ptr,
         size_type new_source_node_id):
            edge_proxy_id{edge_id},
            parent_graph{graph_ptr},
            source_node_internal_id{new_source_node_id} {}

    /** Param constructor for directed edge proxy exposed to friend Graph.
     * @brief Constructs fully specified edge using edge id and graph pointer.
     *   Takes first node of edge in internal representation as source node.
     */
    Edge(size_type edge_id,
         const Graph* graph_ptr): 
            Edge(edge_id,
                 graph_ptr,
                 graph_ptr->node_int_of_edge(edge_id, 0)) {}

    /** Reassigns properties of proxy edge
     * @brief Overwrite proxy id, pointer to parent graph, and
     *   source node internal id. 
    */
    void set_new(size_type new_edge_id,
                 const Graph* graph_ptr,
                 size_type new_source_node_id){
      edge_proxy_id = new_edge_id;
      parent_graph = graph_ptr;
      source_node_internal_id = new_source_node_id;
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    return edge_index_mapping.size();
  }

  /** Return the edge with index @a proxy_edge_index.
   * @pre 0 <= @a proxy_edge_index < num_edges()
   *
   * Complexity: O(1) due to vector backing
   */
  Edge edge(size_type proxy_edge_index) const {
    if (has_edge_index(proxy_edge_index)){
      return Edge(proxy_edge_index, this);
    }
    return Edge();
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid dinstinct nodes of this graph
   * @return True if for some i, edge i connects @a a and @a b.
   *
   * Complexity: O(log(num_nodes())) due to map
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Check membership of target edge index in adjacency map
    auto it_outer = adjacency_map.find(index_at(a.proxy_id));
    if (it_outer != adjacency_map.end()){
      if (it_outer->second.find(index_at(b.proxy_id))
          != it_outer->second.end()){
        return true;
      }
    }
    return false;
  }

  /** Add an edge to the graph, or return the edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b or
   *   vice versa
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: O(log(num_nodes())) due to RBT-backed map
   */
  Edge add_edge(const Node& a, const Node& b) {
    // Find undirected edge {a,b}; if it exists, return it
    auto it_outer = adjacency_map.find(index_at(a.proxy_id));
    if (it_outer != adjacency_map.end()){
      auto it_inner = it_outer->second.find(index_at(b.proxy_id));
      if (it_inner != it_outer->second.end()){
        // Return edge constructed with edge proxy id
        return Edge(edges_internal[it_inner->second].id, this);
      }
    }

    // Save index to top of the array (location for new element)
    auto array_top = edges_internal.size();
    // Create new edge: update internal representation in array
    edges_internal.emplace_back(index_at(a.proxy_id),
                                index_at(b.proxy_id),
                                edge_index_mapping.size());
    // Insert (a,b) directed edge into adjacency map with internal edge id
    it_outer->second.insert(std::pair<size_type, size_type>(
      index_at(b.proxy_id), array_top));

    // Reassign outer iterator to reverse edge
    it_outer = adjacency_map.find(index_at(b.proxy_id));
    assert(it_outer != adjacency_map.end());
    // Insert (b,a) directed edge into adjacency map with internal edge id
    it_outer->second.insert(std::pair<size_type, size_type>(
      index_at(a.proxy_id), array_top));

    // Reflect edge addition in index mapping layer
    edge_index_mapping.push_back(array_top);

    // Count degree of each node a and b
    ++nodes_internal[index_at(a.proxy_id)].degree;
    ++nodes_internal[index_at(b.proxy_id)].degree;

    return Edge(edge_index_mapping.size()-1, this);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_internal.clear();
    edges_internal.clear();
    adjacency_map.clear();
    node_index_mapping.clear();
    edge_index_mapping.clear();
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
    NodeIterator() { }

    /** Returns proxy of iterated node.
     * @return Proxy of iterated node returned; if iterator refers to end,
     *         returns default constructed object instead
     * @pre 0 <= @a this->index <= num_nodes */
    Node operator*() const {
      // Return invalid node if end of collection reached
      if (index == parent_graph->node_end().index)
        return Node();
      return Node(index, parent_graph);
    }

    /** Increment operator.
     * @return Reference to @a this
     * @pre 0 <= @a this->index <= num_nodes
     * @post 0 <= @a this->index <= num_nodes
     * 
     * @details Increments the underlying iterator if the end of the
     * parent graph has not been reached; otherwise does nothing.
     */
    NodeIterator& operator++() {
      if (index != parent_graph->node_end().index)
        ++index;
      return *this;
    }

    /** Equality comparison operator. 
     * @param[in] other Iterator object @a this is compared to
     * @return Equality of underlying iterator and parent
     * @pre @a this and @a other are constructed with valid underlying
     *   iterator and parent graph
     */
    bool operator==(const NodeIterator& other) const {
      return parent_graph == other.parent_graph &&
             index == other.index;
    }

   private:
    friend class Graph;

    // Array index (proxy) as underlying iterator
    size_type index;
    // Parent graph pointer
    const Graph* parent_graph;

    /** Graph-accessible parametrized constructor.
     * @brief Constructs edge iterator using a proxy index and the parent
     * graph's pointer.
     */
    NodeIterator(size_type start_index, const Graph* graph_ptr)
      : index{start_index}, parent_graph{graph_ptr} {}
  };

  /** Returns iterator referring to beginning of node collection.
   * @return Beginning iterator valid for graph lifetime.
   */
  node_iterator node_begin() const {
    return node_iterator(0, this);
  }

  /** Returns iterator referring to end of node collection.
   * @return End iterator valid until node set size changes.
   */
  node_iterator node_end() const {
    return node_iterator(node_index_mapping.size(), this);
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
    IncidentIterator() { }

    /** Returns proxy of iterated object.
     * @return Proxy of iterated incident edge; if iterator refers to end,
     *         returns (invalid) default constructed edge instead
     * @pre iterator refers to a valid incident edge or end of incident set
     */
    Edge operator*() const {
      // Get end iterator
      auto end_iter = parent_graph->incident_edge_end(
        parent_graph->node_proxy_at_internal(
        source_node_internal_id));
      // End-check the underlying iterator
      if (wrapped_iterator == end_iter.wrapped_iterator){
          return Edge();
      }
      return Edge(parent_graph->edge_proxy_at_internal(wrapped_iterator->second),
                  parent_graph, this->source_node_internal_id);
    }

    /** Increment operator.
     * @return Reference to @a this
     * @pre @a this refers to a valid incident edge or end of incident set
     * @post @a this refers to a valid incident edge or end of incident set
     * 
     * @details Increments the underlying iterator if the end of the
     * incident set of nodes has not been reached; otherwise does nothing.
     */
    IncidentIterator& operator++() {
      // Get end iterator
      auto end_iter = parent_graph->incident_edge_end(
        parent_graph->node_proxy_at_internal(
        source_node_internal_id));
      // End-check the underlying iterator
      if (wrapped_iterator != end_iter.wrapped_iterator){
          // Delegate increment to wrapped iterator type
          ++wrapped_iterator;
      }
      return *this;
    }

    /** Equality comparison operator. 
     * @param[in] other Iterator object @a this is compared to
     * @return Equality of underlying iterator and parent
     * @pre @a this and @a other are constructed with valid underlying
     *   iterator and parent graph
     * @details Delegates equality of iterator to the underlying
     *   STL iterator class.
     */
    bool operator==(const IncidentIterator& other) const {
      return parent_graph == other.parent_graph &&
             wrapped_iterator == other.wrapped_iterator;
    }

   private:
    friend class Graph;

    // Alias for the wrapped inner map of adjacency map of maps
    using wrapped_container_type = std::map<size_type,size_type>;
    // Underlying iterator to inner map
    wrapped_container_type::const_iterator wrapped_iterator;
    // Parent graph pointer
    const Graph* parent_graph;
    // Unique id of node that edges are incident to
    size_type source_node_internal_id;

    /** Graph-accessible parametrized constructor.
     * @brief Constructs incident iterator using an iterator to
     *   the inner map structure of the adjacency map, a pointer to
     *   the graph and the unique id of the source node.
    */
    IncidentIterator(const wrapped_container_type::const_iterator &it,
                     const Graph* graph_ptr,
                     size_type new_source_node_id)
                     : wrapped_iterator{it},
                       parent_graph{graph_ptr},
                       source_node_internal_id {new_source_node_id} {}
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
    EdgeIterator() { }

    /** Returns proxy of iterated edge.
     * @return Proxy of iterated edge returned; if iterator refers to end,
     *         returns default constructed object instead
     * @pre 0 <= @a this->index <= num_edges */
    Edge operator*() const {
      // Return invalid node if end of collection reached
      if (index == parent_graph->edge_end().index)
        return Edge();
      return Edge(index, parent_graph);
    }

    /** Increment operator.
     * @return Reference to @a this
     * @pre 0 <= @a this->index <= num_edges
     * @post 0 <= @a this->index <= num_edges
     * @details Increments the underlying iterator if the end of the
     * parent graph has not been reached; otherwise does nothing.
     */
    EdgeIterator& operator++() {
      if (index != parent_graph->edge_end().index)
        ++index;
      return *this;
    }

    /** Equality comparison operator. 
     * @param[in] other Iterator object @a this is compared to
     * @return Equality of underlying iterator and parent
     * @pre @a this and @a other are constructed with valid underlying
     *   iterator and parent graph
     */
    bool operator==(const EdgeIterator& other) const {
      return parent_graph == other.parent_graph &&
             index == other.index;
    }

   private:
    friend class Graph;

    // Array-index backing as underlying iterator
    size_type index;
    // Parent graph pointer
    const Graph* parent_graph;

    /** Graph-accessible parametrized constructor.
     * @brief Constructs edge iterator using a proxy edge index and the parent
     * graph's pointer.
     */
    EdgeIterator(size_type start_index, const Graph* graph_ptr)
      : index{start_index}, parent_graph{graph_ptr} {}
  };

  /** Returns iterator referring to beginning of edge collection.
   * @return Beginning iterator valid for graph lifetime.
   */
  edge_iterator edge_begin() const {
    return edge_iterator(0, this);
  }

  /** Returns iterator referring to end of edge collection.
   * @return End iterator valid until edge set size changes.
   */
  edge_iterator edge_end() const {
    return edge_iterator(edge_index_mapping.size(), this);
  }

  //
  // Deletion methods
  //

  /** Removes node specified by proxy from graph
   * @param[in] n Proxy node object to delete
   * @return Index to next node, or if none, to last possible node
   * @pre @a n is a valid proxy node
   * @post If old @a n was a valid node with index s.t.
   *   0 <= index < old num_nodes, then:
   *     1. the corresponding internal node is inaccessible
   *     2. proxy nodes are invalidated
   *     3. all incident undirected edges are deleted
   *     4. node iterators are valid but out of order
   *     5. all edge and incident iterators are invalidated.
   *   Otherwise, nothing happens.
   * @details Hides internal data until graph is cleared, and redirects
   *   requests using proxy indices via an index mapping layer.
   * 
   * Leaves fragmented internal node data for aggressive runtime optimization.
   * 
   * Runtime: O(deg*log(num_edges)) where deg* is the max degree of any vertex;
   *   essentially O(log(num_nodes^2)) = O(log(num_nodes)) in sparse graphs.
  */
  size_type remove_node(const Node& n){
    // Safeguard for double deletion
    if(!has_node_index(n.proxy_id)){
      return n.index();
    }
    // Remove all neighbouring edges
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
      remove_edge(*it);
    }
    // Re-assign proxy index corresponding to stack-top
    nodes_internal[node_index_mapping[node_index_mapping.size()-1]]
      .id = n.index();
    // Reflect deletion in index mapping layer using swap and pop
    std::swap(node_index_mapping[n.proxy_id],
              node_index_mapping[node_index_mapping.size()-1]);
    node_index_mapping.pop_back();
    // Return index, clamped to num_nodes
    return std::min(n.index(), num_nodes());
  }

  /** Removes node specified by node iterator from graph
   * @param[in] n_it Node iterator to delete
   * @return iterator to next node, or if none, to end
   * @pre @a n is a valid node iterator
   * @post If old @a *n_it was a valid node with index s.t.
   *   0 <= index < old num_nodes, then:
   *     1. the corresponding internal node is inaccessible
   *     2. proxy nodes are invalidated
   *     3. all incident undirected edges are deleted
   *     4. node iterators are valid but out of order
   *     5. all edge and incident iterators are invalidated.
   *   Otherwise, nothing happens.
   * @details Hides internal data until graph is cleared, and redirects
   *   requests using proxy indices via an index mapping layer.
   * 
   * Leaves fragmented internal node data for aggressive runtime optimization.
   * 
   * Runtime: O(deg*log(num_edges)) where deg* is the max degree of any vertex;
   *   essentially O(log(num_nodes^2)) = O(log(num_nodes)) in sparse graphs.
  */
  node_iterator remove_node(node_iterator n_it){
    remove_node(*n_it);
    return n_it;
  }

   /** Removes edge specified by endpoints from graph
   * @param[in] n1, n2 nodes defining the edge
   * @return 1 if deletion was successful; 0 otherwise
   * @pre @a n1, n2 are valid proxy node objects
   * @post If old edge was in the graph, then:
   *     1. the corresponding internal edge is inaccessible
   *     2. proxy edges are invalidated
   *     3. edge iterators are valid but out of order
   *     4. all incident iterators are invalidated.
   *   Otherwise, nothing happens.
   * @details Hides internal data until graph is cleared, and redirects
   *   requests using proxy indices via an index mapping layer. Adjacency
   *   map representation is deleted for incident iterators to work.
   * 
   * Runtime: O(log(num_edges))
  */
  size_type remove_edge(const Node& n1, const Node& n2){
    // Abort if non-existent edge
    if (!has_edge(n1, n2)){
      return 0;
    }
    // Clear adjacency map and return internal edge index
    auto internal_edge_id = clear_adj_map_entry(index_at(n1.index()),
                                       index_at(n2.index()));
    // Clear edge at internal edge index
    clear_edge_array_entry(edge_proxy_at_internal(internal_edge_id));
    return 1;
  }

  /** Removes edge from graph corresponding to given proxy edge
   * @param[in] e Proxy edge to delete
   * @return 1 if deletion was successful; 0 otherwise
   * @pre @a e is a valid proxy edge object
   * @post If old edge was in the graph, then:
   *     1. the corresponding internal edge is inaccessible
   *     2. proxy edges are invalidated
   *     3. edge iterators are valid but out of order
   *     4. all incident iterators are invalidated.
   *   Otherwise, nothing happens.
   * @details Hides internal data until graph is cleared, and redirects
   *   requests using proxy indices via an index mapping layer. Adjacency
   *   map representation is deleted for incident iterators to work.
   * 
   * Runtime: O(log(num_edges))
   */
  size_type remove_edge(const Edge& e){
    // Abort if non-existent edge
    if (!has_edge_index(e.index())){
      return 0;
    }
    // Clear corresponding edge in adjacency map
    clear_adj_map_entry(node_int_of_edge(e.index(), 0),
                        node_int_of_edge(e.index(), 1));
    // Clear edge at proxy edge index and get internal node ids of endpoints
    clear_edge_array_entry(e.index());
    return 1;
  }

  /** Removes edge from graph corresponding to edge iterator position
   * @param[in] e_it Iterator to edge to delete
   * @return Edge iterator to next element, or if none, to end
   * @pre @a e is a valid proxy edge object
   * @post If old edge was in the graph, then:
   *     1. the corresponding internal edge is inaccessible
   *     2. proxy edges are invalidated
   *     3. edge iterators are valid but out of order
   *     4. all incident iterators are invalidated.
   *   Otherwise, nothing happens.
   * @details Hides internal data until graph is cleared, and redirects
   *   requests using proxy indices via an index mapping layer. Adjacency
   *   map representation is deleted for incident iterators to work.
   * 
   * Runtime: O(log(num_edges))
  */
  edge_iterator remove_edge(edge_iterator e_it){
    // Abort if end iterator
    if (e_it == edge_end()) { return 0; }
    // Delegate to Edge-parametrized method
    remove_edge(*e_it);
    // Array index should stay at same position
    return e_it;
  }
  
  //
  // Access helpers
  //
  // Private access functions are provided to Node and Edge
  // so the latter can be remain agnostic to the data representation
  // internal to the graph.
  //

  /** Returns internal index corresponding to proxy id.
   * @param[in] Id of proxy node
   * @return index corresponding to input proxy id
   * @pre @a 0 <= @a proxy_id < num_nodes
   * @details Redirects proxy id to internal data location
  */
  inline size_type index_at(size_type proxy_id) const{
    assert(has_node_index(proxy_id));
    return node_index_mapping.at(proxy_id);
  }

  /** Returns internal index corresponding to proxy id.
   * @param[in] Id of proxy node
   * @return index corresponding to input proxy id
   * @pre @a 0 <= @a proxy_id < num_edges
   * @details Redirects proxy id to internal data location
  */
  inline size_type edge_index_at(size_type proxy_id) const{
    assert(has_edge_index(proxy_id));
    return edge_index_mapping.at(proxy_id);
  }

  /** Returns proxy node index corresponding to internal node id.
   * @param[in] internal_id Internal node id
   * @return Corresponding proxy node index 
   * @pre @a 0 <= @a internal_id < @a this->nodes_internal.size()
   * @details Inverse mapping from internal node id to proxy.
   */
  inline size_type node_proxy_at_internal(size_type internal_id) const{
    return nodes_internal[internal_id].id;
  }

  /** Returns proxy edge index corresponding to internal edge id.
   * @param[in] internal_id Internal edge id
   * @return Corresponding proxy edge index 
   * @pre @a 0 <= @a internal_id < @a this->edges_internal.size()
   * @details Inverse mapping from internal edge id to proxy.
   */
  inline size_type edge_proxy_at_internal(size_type internal_id) const{
    return edges_internal[internal_id].id;
  }

 private:

  /** Returns true iff the graph has node @a i. */
  inline bool has_node_index(size_type i) const {
    return i < node_index_mapping.size();
  }

  /** Returns true iff the graph has edge @a i. */
  inline bool has_edge_index(size_type i) const {
    return i < edge_index_mapping.size();
  }

  /** Returns point (spatial representation) at @a proxy_id.
   * @param[in] Id of proxy node
   * @pre @a 0 <= @a proxy_id < num_nodes
   * @return point for node id is returned
   * @details Provides an extra layer of abstraction to the
   *   internal container for nodes.
  */
  Point& point_at(size_type proxy_id) const{
    return const_cast<Point&>(
      nodes_internal.at(index_at(proxy_id)).point);
  }

  /** Returns mutable reference to value of node corresponding to proxy id
   * @param[in] Id of proxy node
   * @pre @a 0 <= @a proxy_id < num_nodes
   * @post index corresponding to  input id is returned
   * @details Provides an extra layer of abstraction to get
   *   references of node values via const graph pointers.
  */
  node_value_type& value_at(size_type proxy_id) const{
    return const_cast<node_value_type&>(
      nodes_internal.at(index_at(proxy_id)).value);
  }

  /** Returns mutable reference to value of edge corresponding to proxy id
   * @param[in] Id of proxy edge
   * @pre @a 0 <= @a proxy_edge_id < this->edges_internal.size()
   * @post index corresponding to  input id is returned
   * @details Provides an extra layer of abstraction to get
   *   references of node values via const graph pointers.
  */
  edge_value_type& value_edge(size_type edge_proxy_id) const{
    assert(has_edge_index(edge_proxy_id));
    return const_cast<edge_value_type&>(
      edges_internal.at(edge_index_at(edge_proxy_id)).value);
  }

  /** Returns the internal index of one of the nodes of specified edge.
   * @param[in] edge_proxy_index index of the edge reckoned by the proxy
   * @param[in] nodeNumber number (0 or 1) of node attached to edge
   * @return Internal node id, such that 0 <= result < num_nodes
   * @pre 0 <= @a edge_proxy_index < num_edges
   * @pre @a nodeNumber in {0,1}
   * @details Returns internal node id corresponding to edge id.
   *   Nodes are not ordered in internal representation.
  */
  size_type node_int_of_edge(size_type edge_proxy_index,
                            size_type nodeNumber) const{
    assert(has_edge_index(edge_proxy_index));
    assert(nodeNumber == 0 || nodeNumber == 1);
    if (nodeNumber == 0)
      return edges_internal.at(edge_index_at(edge_proxy_index)).node_id_first;
    else
      return edges_internal.at(edge_index_at(edge_proxy_index)).node_id_second;
  };

  /** Returns degree of node corresponding to proxy id
   * @param[in] Id of proxy node
   * @pre @a 0 <= @a proxy_id < num_nodes
   * @return degree corresponding to node such that result >= 0
   * @details Provides an extra layer of abstraction to
   *   interface with the internal representation of node degrees.
  */
  size_type degree_at(size_type proxy_id) const {
    return nodes_internal[index_at(proxy_id)].degree;
  }
  
  /** Returns iterator pointing to the beginning of this node's incident set.
   * @param[in] proxy_id Unique id of the proxy node
   * @return Iterator referring to beginning of incident set
   * @pre 0 <= @a proxy_id < num_nodes
   * @post result is valid until deletion of nodes or edges
   * @details Provides a layer of abstraction for interfacing with the adjacency
   *   map, i.e., so that Node does not need to know implementation of
   *   adjacency map
   */
  incident_iterator incident_edge_begin(const size_type proxy_id) const {
    return incident_iterator(adjacency_map.find(index_at(proxy_id))->
      second.begin(), this, index_at(proxy_id));
  }

  /** Returns iterator pointing to the end of this node's incident set.
   * @param[in] proxy_id Unique id of the proxy node
   * @return Iterator referring to end of incident set
   * @pre 0 <= @a proxy_id < num_nodes
   * @post result is valid until insertion or deletion of nodes or edges
   * @details Provides a layer of abstraction for interfacing with the adjacency
   *   map, i.e., so that Node does not need to know implementation of
   *   adjacency map
   */
  incident_iterator incident_edge_end(const size_type proxy_id) const {
    return incident_iterator(adjacency_map.find(index_at(proxy_id))->
      second.end(), this, index_at(proxy_id));
  }

  /** Deletes undirected adjacency map entry.
   * @return Internal edge_id corresponding to the edge
   * @param[in] id_node1 Internal id of node 1
   * @param[in] id_node2 Internal id of node 2
   * 
   * @details Must be called in conjunction with clear_edge_array_entry(...),
   *   before or after this method executes to cleanup the array.
   * 
   * Complexity: O(log (num_edges))
   */
  size_type clear_adj_map_entry(const size_type id_node1,
                                const size_type id_node2) {
    size_type edge_id = -1;
    // Delete edge node 1 -> node 2
    auto it_outer = adjacency_map.find(id_node1);
    if (it_outer != adjacency_map.end()){
      auto it_inner = it_outer->second.find(id_node2);
      if (it_inner != it_outer->second.end()){
        // Grab edge id before deleting
        edge_id = it_inner->second;
        it_outer->second.erase(it_inner);
      }
      else {
        std::cerr << "Illegal delete." << std::endl;
      }
    }
    else {
      std::cerr << "Illegal delete." << std::endl;
    }
    // Delete edge node 2 -> node 1
    it_outer = adjacency_map.find(id_node2);
    if (it_outer != adjacency_map.end()){ 
      auto it_inner = it_outer->second.find(id_node1);
      if (it_inner != it_outer->second.end()){
        it_outer->second.erase(it_inner);
      }
      else {
        std::cerr << "Illegal delete." << std::endl;
      }
    }
    else {
      std::cerr << "Illegal delete." << std::endl;
    }
    // Update degree representation
    --nodes_internal[id_node1].degree;
    --nodes_internal[id_node2].degree;
    return edge_id;
  }

  /** Marks edge id as invalid, "deleting" it from the outside world.
   * @param[in] proxy_id Internal edge id
   * @details Swap and pop in the mapping layer, but leave the internal
   *   edge corresponding to the delete orphaned to reduce computation
   * 
   * Complexity: O(1)
   */
  void clear_edge_array_entry(size_type proxy_id) {
    // Re-assign proxy index corresponding to stack-top
    edges_internal[edge_index_mapping[edge_index_mapping.size()-1]]
      .id = proxy_id;
    // Reflect deletion in index mapping layer using swap and pop
    std::swap(edge_index_mapping[proxy_id],
              edge_index_mapping[edge_index_mapping.size()-1]);
    edge_index_mapping.pop_back();
  }
};

#endif // CME212_GRAPH_HPP
