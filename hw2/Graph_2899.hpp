#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <utility>
#include <cstdint>

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
	// Node and edge data is stored in vectors of structs
  struct NodeElement;                 // struct containing node data
  struct EdgeElement;                 // struct containing edge data
  vector<NodeElement> node_vector;    // vector containing valid node information
  vector<EdgeElement> edge_vector;    // vector containing valid edge information
  vector<int> node_idx;               // vector that maps node id to node index
  vector<int> edge_idx;               // vector that maps edge id to edge index
  const int invalid = -1;  // flag for invalidating nodes and edges 

 public:
  //
  // PUBLIC TYPE DEFINITIONS
  //
  
  /** Type of this graph. */
  using graph_type = Graph;
  /** Type of value. **/
  using node_value_type = V;  
  /** Type of edge value. **/
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
  Graph() : node_vector(), edge_vector() {}

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
   
    /** Return a reference to this node's position, can be used to reset node position. */
    Point& position() {
      return graph->node_vector[graph->node_idx[uid]].pos;
    }
    
    /** Return this node's position, not to be changed. */
    const Point& position() const {
      return graph->node_vector[graph->node_idx[uid]].pos;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph->node_idx[uid];
    }

    /** Return a reference to this node's value, can be used to reset node value. */
    node_value_type& value() {
      return graph->node_vector[graph->node_idx[uid]].value;   
    } 
     
    /** Return this node's value, not to be changed. */
    const node_value_type& value() const {
      return const_cast<node_value_type>(graph->node_vector[graph->node_idx[uid]].value);  
    }
    
    /** Return the number of edges that currently connect to this node. **/	
    size_type degree() const {
      return graph->node_vector[graph->node_idx[uid]].node_map.size();
    }
    
    /** Return iterator start point to traverse connected edges. **/	
    incident_iterator edge_begin() const {
      return IncidentIterator(graph, uid, 0);  
    }
    
    /** Return iterator end point to traverse connected edges. **/	
    incident_iterator edge_end() const {
      return IncidentIterator(graph, uid, degree()); 
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (n.uid == uid && n.graph == graph);
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
      return ((((uintptr_t) graph) + uid) < (((uintptr_t) n.graph) + n.uid));
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph;
    size_type uid;
    Node(const Graph* graph_, size_type id_) : graph(const_cast<Graph*>(graph_)), uid(id_) {}
  };
  //*****END NODE DEF*****

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_vector.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value    The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    size_type uid = node_idx.size();
    node_idx.push_back(size());
    node_vector.push_back(NodeElement(uid, position, value));
    return Node(this, uid);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {  
    if (n.uid >= node_idx.size() || node_idx[n.uid] == invalid) {
      return false;
    }
    return (node(node_idx[n.uid]) == n);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, node_vector[i].uid);
  }
  
  /** Removes node from the graph. Removes all edges connected to node. Returns true if successfully removed.
   * @param[in] n     Node to be deleted
   * @return true/false indicating whether node was successfully removed.
   *
   * @pre 0 <= @a i < old num_nodes()
   * @post 0 <= @a i < old num_nodes()-1
   * @post new num_nodes() == old num_nodes() - 1
   * @post new num_edges() == old num_edges() - n.degree()
   * @post Node @a n is invalidated.
   * @post old index @a i = old num_nodes()-1 is invalidated.
   * @post All edges connecting to node @a n are invalidated.
   * @post All node_iterators are invalidated.
   * @post All edge_iterators are invalidated.
   *
   * Complexity: O(n1.degree())
   */
  size_type remove_node(const Node& n) {
    // Check if node exists
    if (!has_node(n)) {
      return false;
    }
    // Remove all edges
    Node n_copy = n;
    incident_iterator ii = n_copy.edge_begin();
    while (ii != n_copy.edge_end()) {
      remove_edge(*ii);
    }
    // Reassign mapping and remove node
    size_type cur_id = n.uid;   
    size_type end_id = node_vector.back().uid;
    if (cur_id != end_id) {
      swap(node_vector[node_idx[cur_id]], node_vector[node_idx[end_id]]);
      node_idx[end_id] = node_idx[cur_id];
    } 
    node_idx[cur_id] = invalid;
    node_vector.pop_back();
    return true;
  }
  
  
  /** Removes node from the graph. Removes all edges connected to node. Returns iterator pointing to an unvisited node.
   * @param[in] n_it    node_iterator pointing to the node @a n to be deleted
   * @return node_iterator pointing to a node that has not been iterated over.
   *
   * @pre 0 <= @a i < old num_nodes()
   * @post 0 <= @a i < old num_nodes()-1
   * @post new num_nodes() == old num_nodes() - 1
   * @post new num_edges() == old num_edges() - n.degree()
   * @post Node @a n is invalidated.
   * @post old index @a i = old num_nodes()-1 is invalidated.
   * @post All edges connecting to node @a n are invalidated.
   * @post @a old node_end() is invalidated. 
   *
   * Complexity: O(n.degree())
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return n_it;
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
      if (do_swap) {
        return Node(graph, graph->edge_vector[graph->edge_idx[edge_uid]].node2);
      }
      return Node(graph, graph->edge_vector[graph->edge_idx[edge_uid]].node1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      if (do_swap) {
        return Node(graph, graph->edge_vector[graph->edge_idx[edge_uid]].node1);
      }
      return Node(graph, graph->edge_vector[graph->edge_idx[edge_uid]].node2);
    }
    
    /** Return value of the Edge for modification */
    edge_value_type& value() {
      return graph->edge_vector[graph->edge_idx[edge_uid]].value;
    }
    
    /** Return value of the Edge, not to be modified */
    const edge_value_type& value() const {
      return graph->edge_vector[graph->edge_idx[edge_uid]].value;
    }

    /** Return length of the Edge */
    double length() const {
      return norm(node1().position()-node2().position());
    }
    
    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((e.node1() == node1() && e.node2() == node2()) || 
              (e.node1() == node2() && e.node2() == node1())) ;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return ((((uintptr_t) graph) + edge_uid) < (((uintptr_t) e.graph) + e.edge_uid));
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph;
    size_type edge_uid;
    bool do_swap;
    Edge(const Graph* graph_, size_type edge_id_) : graph(const_cast<Graph*>(graph_)), edge_uid(edge_id_), do_swap(false) {}
    Edge(const Graph* graph_, size_type edge_id_, bool swap_) : graph(const_cast<Graph*>(graph_)), edge_uid(edge_id_), do_swap(swap_) {}
  };
  //*****END EDGE DEF*****

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_vector.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, edge_vector[i].edge_uid);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Check if nodes are in the graph
    if (!has_node(a) || !has_node(b)) {
      return false;
    }
    // Iterate through adjacency list for node a to check for node b 
    vector<size_type> n_map = node_vector[node_idx[a.uid]].node_map;
    if (!n_map.empty() && a.graph == b.graph) {
      vector<size_type>::iterator it = find(n_map.begin(), n_map.end(), b.uid);
      if (it != n_map.end()) {
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
    // Iterate through adjacency list for node a to check for node b 
    vector<size_type>& n_map = node_vector[node_idx[a.uid]].node_map;
    if (!n_map.empty() && a.graph == b.graph) {
      vector<size_type>::iterator it = find(n_map.begin(), n_map.end(), b.uid);
      if (it != n_map.end()) {
        size_type edge_loc = distance(n_map.begin(), it);
        return Edge(this, node_vector[node_idx[a.uid]].edge_map[edge_loc]);
      }
    }
    // Add new edge if not existing
    size_type edge_uid = edge_idx.size();
    edge_idx.push_back(num_edges());
    edge_vector.push_back(EdgeElement(edge_uid, a.uid, b.uid));
    // Add edge and edge_ids to adjacency lists
    node_vector[node_idx[a.uid]].node_map.push_back(b.uid);
    node_vector[node_idx[b.uid]].node_map.push_back(a.uid);
    node_vector[node_idx[a.uid]].edge_map.push_back(edge_uid);
    node_vector[node_idx[b.uid]].edge_map.push_back(edge_uid);
    return edge(edge_idx[edge_uid]);
  }

  /** Remove edge from the graph. Returns true if successfully removed.
   * @param[in] n1, n2    nodes that the edge connects.
   * @return true/false indicating whether edge was successfully removed.
   *
   * @pre 0 <= @a i < old num_edges()
   * @post 0 <= @a i < old edges()-1
   * @post new num_edges() == old edges() - 1
   * @post Edge @a e is invalidated.
   * @post old index @a i = old num_edges()-1 is invalidated.
   * @post All edge_iterators are invalidated.
   *
   * Complexity: O(n1.degree())
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    // Check if edge exists
    if (!has_edge(n1, n2)) {
      return false;
    }
    // Find edge uid, then remove using uid.
    vector<size_type> n_map = node_vector[node_idx[n1.uid]].node_map;
    vector<size_type>::iterator it = find(n_map.begin(), n_map.end(), n2.uid);
    size_type loc = distance(n_map.begin(), it);
    size_type e_uid = node_vector[node_idx[n1.uid]].edge_map[loc];
    return remove_edge(edge(edge_idx[e_uid]));
  }
  
  /** Remove edge from the graph. Returns true if successfully removed.
   * @param[in] e    edge to be removed.
   * @return true/false indicating whether edge was successfully removed.
   *
   * @pre 0 <= @a i < old num_edges()
   * @post 0 <= @a i < old edges()-1
   * @post new num_edges() == old edges() - 1
   * @post Edge @a e is invalidated.
   * @post old index @a i = old num_edges()-1 is invalidated.
   * @post All edge_iterators are invalidated.
   *
   * Complexity: O(n1.degree())
   */
  size_type remove_edge(const Edge& e) {
    // Check if edge exists
    if (!has_edge(e.node1(), e.node2())) {
      return false;
    }
    // Remove current edge from node1 adjacency lists
    size_type n_uid = e.node1().uid;  
    vector<size_type> e_map = node_vector[node_idx[n_uid]].edge_map; 
    vector<size_type>::iterator it = find(e_map.begin(), e_map.end(), e.edge_uid);
    size_type loc = distance(e_map.begin(), it);
    if (loc != e_map.size() - 1) {
      node_vector[node_idx[n_uid]].edge_map[loc] = node_vector[node_idx[n_uid]].edge_map.back();
      node_vector[node_idx[n_uid]].node_map[loc] = node_vector[node_idx[n_uid]].node_map.back();
    }
    node_vector[node_idx[n_uid]].edge_map.pop_back();
    node_vector[node_idx[n_uid]].node_map.pop_back();
    // Remove current edge from node2 adjacency lists
    n_uid = e.node2().uid;  
    e_map = node_vector[node_idx[n_uid]].edge_map; 
    it = find(e_map.begin(), e_map.end(), e.edge_uid);
    loc = distance(e_map.begin(), it);
    if (loc != e_map.size() - 1) {
      node_vector[node_idx[n_uid]].edge_map[loc] = node_vector[node_idx[n_uid]].edge_map.back();
      node_vector[node_idx[n_uid]].node_map[loc] = node_vector[node_idx[n_uid]].node_map.back();
    }
    node_vector[node_idx[n_uid]].edge_map.pop_back();
    node_vector[node_idx[n_uid]].node_map.pop_back();
    // Reassign mapping and remove edge
    size_type cur_id = e.edge_uid;
    size_type end_id = edge_vector.back().edge_uid;    
    if (cur_id != end_id) {
      swap(edge_vector[edge_idx[cur_id]], edge_vector[edge_idx[end_id]]);
      edge_idx[end_id] = edge_idx[cur_id];
    }
    edge_idx[cur_id] = invalid;
    edge_vector.pop_back();
    return true;
  }
  
  /** Remove edge from the graph.  Returns iterator pointing to an unvisited edge.
   * @param[in] e_it    edge_iterator pointint to edge to be removed.
   * @return edge_iterator pointing to an edge that has not been iterated over.
   *
   * @pre 0 <= @a i < old num_edges()
   * @post 0 <= @a i < old edges()-1
   * @post new num_edges() == old edges() - 1
   * @post Edge @a e is invalidated.
   * @post old index @a i = old num_edges()-1 is invalidated.
   * @post old @a edge_end() is invalidated. 
   *
   * Complexity: O(n1.degree())
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return e_it;
  }  
  
  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    node_vector.clear(); 
    edge_vector.clear();
    fill(node_idx.begin(), node_idx.end(), invalid);
    fill(edge_idx.begin(), edge_idx.end(), invalid);
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator>{
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
    
    /** Returns the dereferenced node currently pointed to by the iterator. */
    Node operator*() const {
      return Node(graph, graph->node_vector[current_id].uid);
    }
  
    /** Increments the NodeIterator to point to the next node. */
    NodeIterator& operator++() {
      current_id++;
      return *this;
    }
    
	  /** Checks if two NodeIterators are currently pointing to the same node during iterator traversal. */
    bool operator==(const NodeIterator& ni) const {
      return ((current_id == ni.current_id) && (graph == ni.graph));
    }
  
   private:
    friend class Graph;
    Graph* graph;
    size_type current_id;
    NodeIterator(const Graph* graph_, size_type id_) : graph(const_cast<Graph*>(graph_)), current_id(id_) {}    
  };

  /** Return iterator start point to traverse all nodes in a graph. */	
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }
  
  /** Return iterator end point to traverse all nodes in a graph. */	
  node_iterator node_end() const {
    return NodeIterator(this, size());
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
    
    /** Returns the dereferenced edge currently pointed to by the iterator. Swaps order of
     *  the edge so that node1 points to the root node and node2 points to the connected node. */
    Edge operator*() const {
      size_type edge_uid = graph->node_vector[graph->node_idx[node_uid]].edge_map[map_id];
      if (graph->edge_vector[graph->edge_idx[edge_uid]].node1 != node_uid) {
        return Edge(graph, edge_uid, true); 
      }
      return Edge(graph, edge_uid);
    }
	
    /** Increments the IncidentIterator to point to the next edge. */
    IncidentIterator& operator++() {
      map_id++;
      return *this;
    }
    
    /** Checks if two IncidentIterators are currently pointing to the same edge during iterator traversal. */
    bool operator==(const IncidentIterator& ii) const {
      return ((node_uid == ii.node_uid) && (map_id == ii.map_id) && (graph == ii.graph));
    }
      
   private:
    friend class Graph;
    Graph* graph;
    size_type node_uid;
    size_type map_id;
    IncidentIterator(const Graph* graph_, size_type node_id_, size_type map_id_) : 
      graph(const_cast<Graph*>(graph_)), node_uid(node_id_), map_id(map_id_) {}
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

    /** Returns the dereferenced edge currently pointed to by the iterator. */
    Edge operator*() const {
      return Edge(graph, graph->edge_vector[current_id].edge_uid);
    }
    
	  /** Increments the EdgeIterator to point to the next edge. */
    EdgeIterator& operator++() {
      current_id++;
      return *this;
    }
    
	  /** Checks if two EdgeIterators are currently pointing to the same edge during iterator traversal. */
    bool operator==(const EdgeIterator& ei) const {
      return ((current_id == ei.current_id) && (graph == ei.graph));
    }

   private:
    friend class Graph;
    Graph* graph;
    size_type current_id;
    EdgeIterator(const Graph* graph_, size_type id_) : graph(const_cast<Graph*>(graph_)), current_id(id_) {} 
  };

  /** Return iterator start point to traverse all edges in a graph. */	
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
  
  /** Return iterator end point to traverse all edges in a graph. */	
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

 private:
  /** NodeElement structure to store all node data. Contains unique id, position, value, and adjacency lists. */
  struct NodeElement {
    NodeElement(const size_type id_, const Point p) : uid(id_), pos(p), node_map(), edge_map() {}
    NodeElement(const size_type id_, const Point p, node_value_type v) : uid(id_), pos(p), value(v), node_map(), edge_map() {}
    size_type uid;                      //node unique id
    Point pos;                          //node position
    node_value_type value;              //node value
    vector<size_type> node_map;         //adjacency list: node_idx -> node_uid
    vector<size_type> edge_map;         //adjacency list: node_idx -> edge_uid
  };
  
  /** EdgeElement structure to store all edge data. Contains unique id, position, value, and two connected nodes. */
  struct EdgeElement {
    EdgeElement(size_type id_, size_type n1, size_type n2) : edge_uid(id_), node1(n1), node2(n2) {}
    size_type edge_uid;     //edge unique id
    edge_value_type value;  //edge value
    size_type node1;        //node1 id
    size_type node2;        //node2 id
  };

};

#endif // CME212_GRAPH_HPP

