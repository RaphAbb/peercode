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

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template<typename V, typename E>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  // Struct of graph nodes and edges
  struct graph_node;
  struct graph_edge;

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

  /** Type of values stored in the node and edge */
  using node_value_type = V;
  using edge_value_type = E;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    nodes = std::vector<graph_node>();
    edges = std::vector<graph_edge>(); 
    next_node_uid_=0; next_edge_uid_=0;
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
    Point& position() {
      return fetch().point;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return fetch().point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return graph_->node_uid_index[uid_];
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Get the value associated with the node */
    node_value_type& value() {
      return fetch().value;
    }

    const node_value_type& value() const {
      return fetch().value;
    }

    /** Get the degree (number of edges incident) of the node */
    size_type degree() const {
      return graph_->node_edge_list[this->index()].size();
    }

    /** Get an incident_iterator to start traversing the incident edges of
     *  this node */
   incident_iterator edge_begin() const {
      // Get the node_id of first element stored in node_edge_vecmap[node]
      size_type node_index = this->index();
      auto it = graph_->node_edge_vecmap[node_index].begin();

      return IncidentIterator(graph_,uid_,it->first);
    }

    /** Get an incident_iterator to the end of traversing the incident edges
     *  of this node */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_,uid_,uid_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if ((this->graph_ == n.graph_) && (this->index()==n.index()))
          return true;
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
      if (this->index() < n.index()) return true;
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_;
    size_type uid_;

    /* Fetch the node of the graph with its uid == this->uid_ */
    graph_node& fetch() const {
        auto node_index = this->graph_->node_uid_index[this->uid_];
        return this->graph_->nodes[node_index];
    }

    Node(const Graph* graph, size_type uid) : graph_(const_cast<Graph*>
                                         (graph)),uid_(uid) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
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
  Node add_node(const Point& position, 
                const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    
    // Initialize a node of the graph with its appropriate member attributes
    graph_node n_node = { position, next_node_uid_,value };

    // Update the map of uid vs index
    node_uid_index[next_node_uid_] = num_nodes();

    // Add the node to the vector of graph nodes
    nodes.push_back(n_node);

    ++next_node_uid_;

    // Push an empty element for node_edge_vecmap[a.index()]
    node_edge_vecmap.push_back({});
    
    return Node(this,next_node_uid_-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (nodes[n.index()].uid == n.uid_) return true;
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
    assert((i>=0)&&(i<num_nodes()));
    return Node(this,nodes[i].uid);
  }

  /** Remove the node @a n.
   *  @param[in] n Node to be removed
   *  @return 0 if @a n is not present in the graph
   *          1 if @a n is successfully removed
   *
   *  @post new num_nodes() == old num_nodes() - 1
   *  @post has_node(@a n) == false
   *  @post new node(i) == old node(i) if i < @a n.index()
   *        new node(i) == old node(i+1) if i >= @a n.index()
   *  @post new num_edges() == old num_edges - @a n.degree()
   *
   *  The node @a n is invalidated. All the edges connected to it are also
   *  invalidated. So, the current node_iterator would be invalidated. But
   *  by construction, ++iterators would function as expected.
   *
   *  Complexity : O(nn+k*ne) where nn - num nodes, ne - num edges, k - degree
   */
  size_type remove_node(const Node& n) {
    if (!has_node(n)) 
      return 0;

    // Id and index of the node to be removed - O(1)
    size_type n_id = n.uid_;
    size_type n_index = node_uid_index[n_id];
    
    // Get the map containing the nodes and edges it is connected to - O(1)
    auto conn_nodes = node_edge_vecmap[n_index];

    // Loop over all the incident edges / nodes - O(degree)
    for (auto& v : conn_nodes) {
      size_type e_id = v.second;   // Edge id

      // Remove the edge - O(num_edge)
      remove_edge(Edge(this,e_id,n_id));
    }
    // Erase the vector component at index n_index - O(num_nodes)
    nodes.erase(nodes.begin()+n_index);

    // Update the node_uid_index map - O(num_nodes)
    update_node_uid_index(n_id);

    // Erase the vector element containing the incidents - O(num_nodes)
    node_edge_vecmap.erase(node_edge_vecmap.begin()+n_index);
    return 1;
  }

  /** Remove the node @a n pointed by the iterator @a n_it
   *  @param[in] n_it Node iterator to the node to be removed
   *  @return Node iterator to the next node
   *
   *  @param n Node pointed by the node iterator
   *
   *  @post new num_nodes() == old num_nodes() - 1
   *  @post has_node(@a n) == false
   *  @post new node(i) == old node(i) if i < @a n.index()
   *        new node(i) == old node(i+1) if i >= @a n.index()
   *  @post new num_edges() == old num_edges - @a n.degree()
   *
   *  The node @a n is invalidated. All the edges connected to it are also
   *  invalidated. So, the current node_iterator would be invalidated. But
   *  by construction, ++iterators would function as expected.
   *
   *  Complexity : O(nn+k*ne) where nn - num nodes, ne - num edges, k - degree
   */
  node_iterator remove_node(node_iterator n_it) {
    auto node = *n_it;  // Node pointed by the iterator
    ++n_it;             // Make the iterator point to the next node
    remove_node(node);  // Remove the node
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      // Check which of node_1 or 2 corresponds to node_1_id_ stored
      if (fetch().node_1.uid_==node_1_id_)
        return fetch().node_1;
      return fetch().node_2;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      if (fetch().node_1.uid_==node_1_id_)
        return fetch().node_2;
      return fetch().node_1;      
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if ((this->node1() == e.node1()) && (this->node2() == e.node2()))
          return true;
      else if ((this->node1() == e.node2()) && (this->node2() == e.node1()))
          return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_ == e.graph_) {
        if (graph_->edge_uid_index[uid_] < e.graph_->edge_uid_index[e.uid_]) 
          return true;
        return false;
      }
      else {
       if (graph_ < e.graph_)
         return true;        
       return false;
      }
    }

    /** Get the length of the edge */
    double length() const {
      return norm(node1().position()-node2().position());
    }

    /** Get the value stored in the edge */
    edge_value_type& value() {
      return fetch().value;
    }

    const edge_value_type& value() const {
      return fetch().value;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type uid_;
    size_type node_1_id_;
    
    // Fetch the edge of the graph with uid == this->uid_
    graph_edge& fetch() const {
        auto edge_index = this->graph_->edge_uid_index[this->uid_];
        return this->graph_->edges[edge_index];
    }

    Edge(const Graph* graph, size_type uid, size_type node_1_id) : 
        graph_(const_cast<Graph*>(graph)), uid_(uid), node_1_id_(node_1_id) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert((i>=0)&&(i<num_edges()));
    return Edge(this,edges[i].uid,edges[i].node_1.uid_);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    assert(has_node(a) && has_node(b));

    if (node_edge_vecmap[a.index()].find(b.uid_) != 
                           node_edge_vecmap[a.index()].end())
      return true;

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
  Edge add_edge(const Node& a, const Node& b,
                const edge_value_type& value = edge_value_type()) {
    // HW0: YOUR CODE HERE
    assert(not(a==b));
    assert(has_node(a) && has_node(b));

    // If the graph already has the edge, return it
    if (has_edge(a,b)) {
      size_type edge_id = node_edge_vecmap[a.index()][b.uid_];
      return Edge(this,edge_id,a.uid_);
    }    
    // If the edge doesn't exist, create a new one           
    graph_edge n_edge = { a, b, next_edge_uid_, value };
    
    // Update edge uid vs index map
    edge_uid_index[next_edge_uid_] = num_edges();

    // Update node_edge_vecmap
    node_edge_vecmap[a.index()][b.uid_] = next_edge_uid_;
    node_edge_vecmap[b.index()][a.uid_] = next_edge_uid_;

    // Add to the edges of the graph vector
    edges.push_back(n_edge);

    ++next_edge_uid_;

    return Edge(this,next_edge_uid_-1,a.uid_);
  }

  /** Remove the edge @a e connected to nodes @a n1 and @a n2
   *  @param[in] n1,n2 Nodes to which the edge to be removed is connected
   *  @return 0 if the edge doesn't exist in the graph
   *          1 if the edge is removed successfully
   *
   *  @param e Edge connecting @a n1 and @a n2
   *
   *  @post new num_edges() == old num_edges() - 1
   *  @post has_edge(@a e) == false
   *  @post new edge(i) == old edge(i) if i < @a e.index()
   *        new edge(i) == old edge(i+1) if i >= @a e.index()
   *  @post new num_nodes() == old num_nodes()
   *
   *  The edge @a e is invalidated. But no nodes are invalidated. The
   *  edge_iterator pointing to the current edge would be invalidated. But
   *  by construction, ++iterators all work perfectly.
   *
   *  Complexity : O(ne) where ne - num edges
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    // Get the indices of the nodes n1, n2 - O(1)
    size_type n1_index = node_uid_index[n1.uid_];
    size_type n2_index = node_uid_index[n2.uid_];

    // Check if the edge exists - O(1)
    auto search = node_edge_vecmap[n2_index].find(n1.uid_);
    if (search == node_edge_vecmap[n2_index].end())
      return 0;

    // Get the id and index of the edge to be removed - O(1)
    size_type edge_id = node_edge_vecmap[n2_index][n1.uid_];
    size_type edge_index = edge_uid_index[edge_id];

    // Remove the edge from the edges vector - O(num_edges)
    edges.erase(edges.begin()+edge_index);

    // Update the edge_uid_index map - O(num_edges)
    update_edge_uid_index(edge_id);

    // Remove the edge entry from the adjacency vecmap - O(1)
    node_edge_vecmap[n2_index].erase(n1.uid_);
    node_edge_vecmap[n1_index].erase(n2.uid_);
    return 1;
  }
   
  /** Remove the edge @a e
   *  @param[in] e Edge to be removed
   *  @return 0 if the edge doesn't exist in the graph
   *          1 if the edge is removed successfully
   *
   *  @post new num_edges() == old num_edges() - 1
   *  @post has_edge(@a e) == false
   *  @post new edge(i) == old edge(i) if i < @a e.index()
   *        new edge(i) == old edge(i+1) if i >= @a e.index()
   *  @post new num_nodes() == old num_nodes()
   *
   *  The edge @a e is invalidated. But no nodes are invalidated. The
   *  edge_iterator pointing to the current edge would be invalidated. But
   *  by construction, ++iterators all work perfectly.
   *
   *  Complexity : O(ne) where ne - num edges
   */
  size_type remove_edge(const Edge& e) {
    // Get the nodes n1 and n2 connected by the edge
    auto n1 = e.node1();
    auto n2 = e.node2();
    return remove_edge(n1,n2);
  }

  /** Remove the edge @a e pointed by the edge iterator @a eit
   *  @param[in] eit Edge iterator pointing to the edge to be removed
   *  @return Edge iterator pointing to the next edge
   *
   *  @post new num_edges() == old num_edges() - 1
   *  @post has_edge(@a e) == false
   *  @post new edge(i) == old edge(i) if i < @a e.index()
   *        new edge(i) == old edge(i+1) if i >= @a e.index()
   *  @post new num_nodes() == old num_nodes()
   *
   *  The edge @a e is invalidated. But no nodes are invalidated. The
   *  edge_iterator pointing to the current edge would be invalidated. But
   *  by construction, ++iterators all work perfectly.
   *
   *  Complexity : O(ne) where ne - num edges
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    auto edge = *e_it;   // Edge pointed to by the edge iterator
    ++e_it;              // Make the iterator to point to the next edge
    remove_edge(edge);   // Remove the edge
    return e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();
    edges.clear();
    node_uid_index.clear();
    edge_uid_index.clear();
    node_edge_vecmap.clear();
    return;
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
    
   /** Dereferencing the iterator and returning the node pointed by
    *  the node-iterator */
    Node operator* () const {
      return Node(graph_,node_id_);
    }

   /** Moving to the next node according to its index */
    NodeIterator& operator++() {
      size_type node_index = graph_->node_uid_index[node_id_];
      node_index++;
      if (node_index < graph_->num_nodes()) {
        node_id_ = graph_->nodes[node_index].uid;
        return *this;
      }
      node_id_++;
      return *this;
    }

   /** Compare two node iterators */
    bool operator==(const NodeIterator& ni) const {
      if ((this->graph_ == ni.graph_) && (this->node_id_ == ni.node_id_))
        return true;
      return false;
    }

    /** Constructor for the class */
    NodeIterator(const Graph* graph, size_type id) : graph_(const_cast<Graph*>
                                                 (graph)),node_id_(id) {}

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type node_id_;

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Get the beginning of the node vector to be traversed.
   *
   *  Check the case of no nodes in the class.
   */
  node_iterator node_begin() const {
    if (num_nodes() > 0)
      return node_iterator(this,nodes[0].uid);
    return node_iterator(this,0);
  }

  /** Get the ending of the node vector to be traversed.
   *
   *  End is defined using end_node_uid + 1. Check for the case of no nodes.
   */
  node_iterator node_end() const {
    if (num_nodes() > 0)
      return node_iterator(this,nodes[num_nodes()-1].uid+1);
    return node_iterator(this,0);
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
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** Dereferencing of the iterator to give the current edge */
    Edge operator*() const {
      // Set the node_id_1 for edge class to be the node id which spawned it
      size_type node_index_1 = graph_->node_uid_index[node_id_1_];
      size_type edge_id = graph_->node_edge_vecmap[node_index_1][node_id_2_];

      return Edge(graph_,edge_id,node_id_1_);
    }

    /** Increment the iterator to move to the next edge.
     *  
     *  Check the case where the it++ goes to map.end().
     *  End iter is defined using node_id_1 = node_id_2.
     */
    IncidentIterator& operator++() {
      size_type node_index_1 = graph_->node_uid_index[node_id_1_];
      auto it = graph_->node_edge_vecmap[node_index_1].find(node_id_2_);
      it++;
      if (it == graph_->node_edge_vecmap[node_index_1].end()) {
        node_id_2_ = node_id_1_;
        return *this;
      }
      node_id_2_ = it->first;
      return *this;
    }

    /** Compare two incident iterators for equality.
     *  
     *  Check if the graphs they point to, node_id_1 and node_id_2 are same.
     */
    bool operator==(const IncidentIterator& iit) const {
      if ((graph_ == iit.graph_) && (node_id_1_ == iit.node_id_1_) &&
                                       (node_id_2_ == iit.node_id_2_))
        return true;

      return false;
    }
      
    /** Constructor for the class */
    IncidentIterator(const Graph* graph, size_type node_id_1, 
        size_type node_id_2) : graph_(const_cast<Graph*>(graph)),
              node_id_1_(node_id_1), node_id_2_(node_id_2) {}

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type node_id_1_;
    size_type node_id_2_;
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

    /** Dereferencing the edge iterator to get the edge it is pointing to */
    Edge operator*() const {
      // Choose node_1 id to be node_1's id (arbitrarily)
      size_type edge_index = graph_->edge_uid_index[edge_id_];
      size_type node_id_1 = graph_->edges[edge_index].node_1.uid_;

      return Edge(graph_,edge_id_,node_id_1);
    }

    /** Increment the edge iterator to go to the edge with next index.
     *
     *  Check the case where edge_index++ flows out of edge vector bounds.
     *  End edge iter is defined by using last_edge_id + 1.
     */
    EdgeIterator& operator++() {
      size_type edge_index = graph_->edge_uid_index[edge_id_];
      edge_index++;
      if (edge_index < graph_->num_edges()) {
        edge_id_ = graph_->edges[edge_index].uid;
        return *this;
      }
      edge_id_++;
      return *this;
    }

    /** Compare two edge iterators for equality.
     *
     *  Check if the graphs they point to and the edge_ids are same.
     */
    bool operator==(const EdgeIterator& eit) const {
      if ((graph_ == eit.graph_) && (edge_id_ == eit.edge_id_))
        return true;
      return false;
    }

    /** Constructor for the edge iterators */
    EdgeIterator(const Graph* graph, size_type edge_id) : 
                graph_(const_cast<Graph*>(graph)), edge_id_(edge_id) {}

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type edge_id_;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Get the edge iterator to the beginning of the edge list to be traversed.
   *
   *  Check for the case of no edges in the graph.
   */
  edge_iterator edge_begin() const {
    if (num_edges() > 0) 
      return EdgeIterator(this,edges[0].uid);
    return EdgeIterator(this,0);
  }

  /** Get the edge iterator to the end of the edge list to be traversed.
   *
   *  Check for the case of no edges in the graph.
   *  End edge iter is defined using last_edge_id + 1.
   */
  edge_iterator edge_end() const {
    if (num_edges() > 0)
      return EdgeIterator(this,edges[num_edges()-1].uid+1);
    return EdgeIterator(this,0);
  }

 private:
  /** Function to update the node_uid_index map - O(num_nodes)*/
  void update_node_uid_index(const size_type& id) {
    // Loop over all the node ids in the map
    for (auto it = node_uid_index.begin(); it != node_uid_index.end();) {
      // If the node_id matches the id removed, erase the map element
      if (it->first == id)
        it = node_uid_index.erase(it);
      // Else if the id is greater than id removed, reduce its index by 1
      else if (it->first > id) {
        it->second -= 1;
        ++it;
      }
      // Else if the id is lesser than id removed, just move ahead
      else
        ++it;
    }
    return;
  }

  /** Function to update the edge_uid_index map - O(num_edges */
  void update_edge_uid_index(const size_type& id) {
    // Loop over all the edge ids in the map
    for (auto it = edge_uid_index.begin(); it != edge_uid_index.end();) {
      // If the id is the id removed, erase the map element
      if (it->first == id)
        it = edge_uid_index.erase(it);
      // Else if the id is greater than id removed, reduce the index by 1
      else if (it->first > id) {
        it->second -= 1;
        ++it;
      }
      // Else if the id is lesser, just move ahead
      else
        ++it;
    }
    return;
  }

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct graph_node {
      Point point;
      size_type uid;
      node_value_type value;
  };

  struct graph_edge {
      Node node_1;
      Node node_2;
      size_type uid;
      edge_value_type value;
  };

  // Vectors containing nodes and edges based on their indices
  std::vector<graph_node> nodes;
  std::vector<graph_edge> edges;

  // Maps of uid to the indices for O(1) search 
  std::unordered_map<size_type,size_type> node_uid_index;
  std::unordered_map<size_type,size_type> edge_uid_index;

  // Vector of maps of node_2 vs edge ids for each node_1
  // node_edge_vecmap[node_1.index][node_2.id] = edge_id
  std::vector<std::unordered_map<size_type,size_type>> node_edge_vecmap;

  // Uids of the next node and edge
  size_type next_node_uid_;
  size_type next_edge_uid_;
};

#endif // CME212_GRAPH_HPP


