#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <boost/functional/hash.hpp> 
#include <iostream>
#include <cstdlib>


#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename NV, typename EV>
class Graph {

 private:

  struct node_element;
  struct edge_element;

 public:

  using graph_type = Graph<NV, EV>;

  // Type of indexes and sizes.
  //    Return type of Graph::Node::index(), Graph::num_nodes(),
  //    Graph::num_edges(), and argument type of Graph::node(size_type) 
  using size_type     = unsigned;
  using node_map_type = std::unordered_map<size_type, node_element>;
  using edge_map_type = std::unordered_map<size_type, edge_element>;

  // Predeclaration of Node class and associated types/classes 
  class Node;
    class NodeIterator;
    using node_type       = Node;
    using node_value_type = NV;
    using node_iterator   = NodeIterator;

  // Predeclaration of Edge class and associated types/classes 
  class Edge;
    class EdgeIterator;
    using edge_type       = Edge;
    using edge_value_type = EV;
    using edge_iterator   = EdgeIterator;

  // Predeclaration of IncidentIterator 
  class IncidentIterator;
    using incident_iterator = IncidentIterator;


  // Default constructor
  Graph() 
    : node_map_(), edge_map_(), nodes_to_edge_(), n_nodes_(0), n_edges_(0),
      n_nodes_generated_(0), n_edges_generated_(0) {
      g_idx_ = std::rand();
  }

  // Default destructor 
  ~Graph() = default;

  //  ---------------------------------- Nodes ----------------------------- //

  // @class Graph::Node. Proxy object to get info about Graph's nodes.
  class Node : private totally_ordered<Node> {
   public:

    // Invalid constructor. Valid ctor is private.
    Node() {
    }

    // Grab node with Node ID = id from Graph.
    node_element& get() const {
      return graph_->node_map_[this->uid_];
    }

    // Return this node's position.
    const Point& position() const {
      return get().position;
    }

    // Return reference to edit this node's position.
    Point& position() {
      return get().position;
    }

    // Return Node's index (idx), a number in the range [0, i2u.size()).
    // Used to take in data from external-facing methods, like graph.node()
    size_type index() const { 
      return get().idx;
    }

    // Return Node's index (idx), a number in the range [0, i2u.size()).
    // Used to take in data from external-facing methods, like graph.node()
    size_type index() { 
      return get().idx;
    }

    // Return Node's unique ID, a number in the range [0, node_map_.size()).
    // Used to reference Node's data container and as an identifier for edges.
    size_type uid() const {
      return get().uid;
    }

    // Return this node's container, a pointer to a Graph.
    const Graph* graph() const {
      return this->graph_;
    }

    // Return reference to edit this node's value.
    node_value_type& value() {
      return get().value;
    }

    // Return this node's value.
    const node_value_type& value() const {
      return get().value;
    }

    // Return this node's degree: the number of Nodes it's connected to.
    size_type degree() const {
      return get().connections.size();
    }

    // Return a mutable reference to this node's connections.
    std::vector<edge_type>& connections() const {
      return get().connections;
    }

    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, this, 0); 
    }

    incident_iterator edge_end() const {
      return IncidentIterator(graph_, this, degree());
    }

    // Test whether this node and @a n are equal.
    // Equal nodes have both the same graph and the same idx.
    bool operator==(const Node& n) const {
      return (this->graph() == n.graph() && this->index() == n.index()); // idx √
    }

    // Compares indices of two nodes, not necessarily of the same graph.
    // It has no geometric meaning.
    bool operator<(const Node& n) const {
      assert(graph_);    
      assert(n.graph());
      return (this->index() < n.index()); // idx √ 
    }

   private:

    friend class Graph; // Proxy design pattern necessitates Graph access.

    Graph* graph_;      // Pointer to parent Graph class (container)
    size_type uid_;     // Unique id. Key for graph_'s unordered_map.

    // The only valid constructor. Accessible only through graph_.
    Node(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
    
  };

  // -------------------------- Graph Node Methods ------------------------- //

  // Return the number of nodes in the graph.
  size_type size() const {
    return n_nodes_;
  }
  
  // Synonym of size().
  size_type num_nodes() const {
    return size();
  }

  // Add a new Node to Graph. 
  Node add_node(const Point& position, 
                  const node_value_type& value = node_value_type()) {
   
    // Create and populate the associated node_element data structure.
    node_element node_el;
    node_el.node     = Node(this, n_nodes_generated_);  // use the unique ID
    node_el.uid      = n_nodes_generated_;              // unique ID
    node_el.idx      = n_nodes_;                        // position w/in graph
    node_el.position = position;
    node_el.value    = value;
    
    // Leave node_el.connections unassigned. Add_edge updates this attribute.

    // Store node in Graph's map of nodes.
    node_map_[n_nodes_generated_] = node_el;
    i2u.push_back(n_nodes_generated_);

    // Update node counters
    ++n_nodes_;
    ++n_nodes_generated_;

    return node_el.node;
  }

  // Determine if a Node is a member of this Graph, comparing pointers.
  bool has_node(const Node& n) const {
    
    return this == n.graph();
  }

  // Return the Node with index @a i, for 0 =< i < num_nodes().
  Node node(size_type i) const {
    assert(i < size());
    return Node(this, i2u[i]);  // Valid Node.
  }

  //  ---------------------------------- Edges ----------------------------- //

  /** @class Graph::Edge
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Edge> { 

   public:

    // Invalid constructor. Valid ctor is private.
    Edge() {
    }

    // Return this Edge's index.
    size_type index() const {
      return eid_;
    }

    // Return the first Node of this Edge
    Node node1() const {
      size_type uid = graph_->edge_map_[eid_].node1;
      return graph_->node_map_[uid].node;
    }

    // Return the second Node of this Edge 
    Node node2() const {
      size_type uid = graph_->edge_map_[eid_].node2;
      return graph_->node_map_[uid].node;
    }

    // Return an immuatable copy of Edge's underlying data type
    const edge_value_type& value() const {
      return graph_->edge_map_[eid_].value;
    }

    // Return a muatable copy of Edge's underlying data type
    edge_value_type& value() {
      return graph_->edge_map_[eid_].value;
    }

    // Return this edge's length
    double length() const {
      Point dist = node1().position() - node2().position();
      return norm(dist);
    }

    // Test whether this edge and @a e are equal (i.e. equivalent).
    // Equal edges represent the same undirected edge between two nodes.
    bool operator==(const Edge& e) const {
      if  (     this->graph_ == e.graph_   and
              ((this->node1() == e.node1() and this->node2() == e.node2())
            or (this->node1() == e.node2() and this->node2() == e.node1()))) {
          return true;
        } else {
          return false;
        }
    }

    // Test whether this edge is less than @a e in a global order.
    // Compares Edge IDs.
    // Result has no interpretive meaning.
    bool operator<(const Edge& e) const {
      if (this->graph_->g_idx_ == e.graph_->g_idx_) {
        return this->index() < e.index();
      } else {
        return this->graph_->g_idx_ < e.graph_->g_idx_;
      }
    }

   private:

    friend class Graph; // Proxy design pattern necessitates Graph access.

    Graph* graph_;      // Pointer to parent Graph class (container)
    size_type eid_;     // Unique id. Key for graph_'s unordered_map.

    // The only valid constructor. Accessible only through graph_.
    Edge(const Graph* graph, size_type eid)
      : graph_(const_cast<Graph*>(graph)), eid_(eid) {
    }
    
    // Grab Edge with Edge ID = eid_ from Graph.
    edge_element& get() const {
      return graph_->edge_map_[this->eid_];
    } 

    // Swap unique node IDs of edge. Allows nid1 > nid2 for IncidentIterator.
    // Grabbing an edge from edge_map_ should still be done with sorted keys.
    void swap_node_indices() {
      
      size_type nid1 = graph_->edge_map_[eid_].node1;
      size_type nid2 = graph_->edge_map_[eid_].node2;
      graph_->edge_map_[eid_].node1 = nid2;
      graph_->edge_map_[eid_].node2 = nid1;

      return;
    }

  };

  // -------------------------- Graph Edge Methods ------------------------- //

  // Return the total number of edges in the graph.
  size_type num_edges() const {
    return n_edges_;
  }

  // Return the Edge with index i.
  Edge edge(size_type i) const {
    if (i >= num_edges()) {
      return Edge();         // Invalid Edge
    } else {
      return Edge(this, i);  // Valid Edge
    }
  }

  // Test whether two nodes are connected by an edge.
  // Complexity: O(num_edges()).
  bool has_edge(const Node& a, const Node& b) const {

    assert(has_node(a) and has_node(b));

    if (a != b and a.degree() > 0 and b.degree() > 0) {
      for (auto it = a.edge_begin(); it != a.edge_end(); ++it) {
        auto connection = *it;
        if (connection.node2() == b) {
          // Found it: an edge in @a connecting it to @b.
          return true;
        }
      }
    }

    // No edge exists between Nodes
//     std::cout << "edge not found between " << a.uid() << " and " << b.uid() << std::endl;
    return false;
  }

  // Add an edge to the graph or return the current edge if it already exists.
  // Invalidates neither Edge indices nor outstanding Edge objects.
  // Complexity: No more than O(num_nodes() + num_edges()), hopefully less.
  Edge add_edge(const Node& a, const Node& b) {

    if (has_edge(a, b)) {
      // return this Edge if it already exists.
      // note: has_edge contains assertion that graph has nodes a and b.
      return Edge(this, nodes_to_edge_[node_hash(a, b)]);

    } else {

      // Create the new Edge, with unique id = n_edges_generated_
      Edge new_edge(this, n_edges_generated_);

      // Create and populate the associated edge_element data structure.
      edge_element    edge_el;
      edge_el.edge  = new_edge;

      size_type uid1 = std::min(a.uid(), b.uid()); 
      size_type uid2 = std::max(a.uid(), b.uid());
      edge_el.node1  = uid1;
      edge_el.node2  = uid2;
      edge_el.value  = edge_value_type();

      // Store Edge in Graph's map of Edges.
      edge_map_[n_edges_generated_]   = edge_el;
      nodes_to_edge_[node_hash(a, b)] = n_edges_generated_;

      // Add this edge to "connections" containers of the nodes it connects.
      node_map_[uid1].connections.push_back(new_edge);
      node_map_[uid2].connections.push_back(new_edge);

      // Increment n_edges_ and n_edges_generated_ (latter is unique ID)
      ++n_edges_;
      ++n_edges_generated_;

      return new_edge;  // Return a new, valid Edge.
    }

  }

  // Remove all nodes and edges from this graph.
  // Invalidates all outstanding Node and Edge objects.
  void clear() {
    // Reset Counters
    n_nodes_ = 0;
    n_edges_ = 0;
    n_nodes_generated_ = 0;
    n_edges_generated_ = 0;

    // Reset Maps
    node_map_.clear();
    edge_map_.clear();
    i2u.clear();
    nodes_to_edge_.clear();
  }

  // Search for an Edge connecting the two given Nodes, a and b.
  // Both has_edge and add_edge need the value of this search operation.
  std::unordered_map<std::size_t, size_type>::const_iterator
        node_connection_search(const Node& a, const Node& b) const {

    std::size_t hash = node_hash(a, b);
    std::unordered_map<std::size_t, size_type>::const_iterator search;
    search = nodes_to_edge_.find(hash);
    return search;
  }

  // ----------------- Iterator Members, Associated Aliases ----------------- //

  // ------------------------------- NodeIterator --------------------------- //
  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. Wrapper around std::unordered_map iterator.*/
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

    // Return current node object.
    Node operator*() const {
      assert(idx < graph->size());
      return graph->node_map_[graph->i2u[idx]].node;
    }

    // Advance one element in the i2u vector.
    NodeIterator& operator++() {
      assert(idx < graph->num_nodes());
      idx++;
      return *this;
    }

    // Assess equality of current node with another iterator's node.
    bool operator==(const node_iterator& it) const {
      assert(graph == it.graph);
      return idx == it.idx;
    }

   private:

    friend class Graph;
//     typename std::vector<size_type>::const_iterator it_;
    
    // Initialize variables
    Graph*    graph = nullptr; // Pointer to Graph
    size_type idx;             // Index of current node

    // Construct the iterator, instantiating a begin() pointer.
    NodeIterator(const Graph* graph_, size_type idx_)
      : graph(const_cast<Graph*>(graph_)), idx(idx_) { }

  };

  // Type alias for NodeIterator, returning the begin() pointer.
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  };
  
  // Type alias for NodeIterator, returning the end() pointer.
  node_iterator node_end() const {
    return NodeIterator(this, this->num_nodes());
  };

  // ----------------------------- IncidentIterator ------------------------- //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : totally_ordered<IncidentIterator> {
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

    // Return the current Edge the iterator is on. 
    // If necessary, underlying Edge indices will be adjusted so
    //    base node will be current_edge.node1().
    Edge operator*() const {

      auto edge_list = this->incident_edges_();
      assert(it_ < edge_list.size());
      Edge e = edge_list[it_]; 

      // Swap node indices if necessary to return current node first.
      if (id_ != e.node1().index()) { 
        e.swap_node_indices();
        assert(id_ == e.node1().index()); 
      }

      return e;
    }

    // Move on to the next element.
    IncidentIterator& operator++() {
      assert(it_ < this->incident_edges_().size());
      ++it_;
      return *this;
    }

    // Compare two iterators, checking node indices and position in n.connections().
    bool operator==(const IncidentIterator& x) const {
      return (id_ == x.id_ and it_ == x.it_);
    }

   private:

    friend class Graph;

    // Initialize variables
    Graph* graph_ = nullptr;  // Pointer to parent Graph
    size_type id_;            // index of parent Node
    size_type it_;            // position of iterator

    // Return an immutable reference to node's internal data
    const std::vector<edge_type>& incident_edges_() const {
      assert(graph_ && id_ < graph_->size());
      return graph_->node_map_[graph_->i2u[id_]].connections;
//       return graph_->node_map_[uid_].connections;
    }

    // Return a mutable reference to node's internal data
    std::vector<edge_type>& incident_edges_() {
      assert(graph_);
      assert(id_ < graph_->size());
      return graph_->node_map_[graph_->i2u[id_]].connections;
    }

    // Construct a valid IncidentIterator.
    IncidentIterator(const Graph* graph, const Node* node, size_type index)
      : graph_(const_cast<Graph*>(graph)), it_(index) {
        id_ = node->index();
    }

  };

  // ------------------------------- EdgeIterator --------------------------- //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : totally_ordered<EdgeIterator> {
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

    // Return current edge object.
    Edge operator*() const {
      return Edge(graph_, it_->second.edge.index());
    };

    // Advance one element in the map.
    EdgeIterator& operator++() {
      ++it_;
      return *this;
    };

    // Assess equality of current iterator with another.
    bool operator==(const EdgeIterator& x) const {
      return this->it_ == x.it_;
    };

   private:
    friend class Graph;
    const Graph* graph_;
    typename edge_map_type::const_iterator it_;

    // Construct the iterator, instantiating a begin() pointer.
    EdgeIterator(const Graph* g, typename edge_map_type::const_iterator it)
      : graph_(g), it_(it) {
    }
  };

  // Type alias for EdgeIterator, returning the begin() pointer.
  edge_iterator edge_begin() const {
    return EdgeIterator(this, edge_map_.begin());
  };
  
  // Type alias for EdgeIterator, returning the end() pointer.
  edge_iterator edge_end() const {
    return EdgeIterator(this, edge_map_.end());
  };

 // -------------------------- Edge and Node Removal ------------------------ //
 
 public:
 

  /* @brief Removes a Node from Graph.
   * @param[in] n node to be deleted        
   * @result    number of nodes remaining in Graph.
   *
   * @pre Node is in graph: 
   *         0 =< n.index() < g.num_nodes()
   *         for i s.t. i2u[i] == n.uid(), node_map[i2u[i]].index() == i
   *
   * @post Node @n is not in graph:
   *         there exists no i s.t. i2u[i] == n.uid()
   *
   * @post Graph size is reduced by one: 
   *        new g.num_nodes() == old g.num_nodes() - 1;
   *        new i2u.size()    == old i2u.size() - 1
   * 
   * @post Graph maintains a list of valid nodes:
   *        for all 0 =< i < new g.num_nodes(), g.node_map_[new i2u[i]].index()] == i
   *        for all 0 =< i < n.index(),
   *          old i2u[i]   == new i2u[i] 
   *        for all i =< n.index() < new g.num_nodes(),
   *          old i2u[i-1] == new i2u[i]
   *
   * @post for all i s.t i \in n.connections(0 =< j < n.degree()):
   *        new i.degree() == old i.degree() - 1
   *
   * @post new g.num_edges() == old g.num_edges() - n.degree();
   *        
   * @post Containers indexed by unique indices remain unaffected:
   *        new g.node_map_().size()   == old g.node_map_().size()
   *        new g.edge_map_().size()   == old g.edge_map_().size()
   *        new g.nodes_to_edge.size() == old nodes_to_edge.size()
   *        
   * Complexity is O(degree^2) for degree << n_nodes \sim n_edges
   * Invalidates no other nodes besides Node @n.
   * Invalidates no other edges besides Edge @e \in n.connections().
   * Invalidates i2u iterators for all i in i2u s.t. i =< n.index() < new g.num_nodes(),
   *
   */
  size_type remove_node(const Node& n) {
   
    assert(has_node(n));

    // Remove each Edge connected to this node
    for (auto it = n.edge_begin(); it != n.edge_end(); ) {

      auto e = *it;

      size_type s10 = n.degree();
      size_type s20 = e.node2().degree();

      // Remove @e from the following:
      assert(has_edge(n, e.node2()));
      delete_connection(e.node2(), e);  // @e.node2()'s connections;
      delete_connection(n, e);          // @e.node2()'s connections;

      size_type s11 = n.degree();
      size_type s21 = e.node2().degree();
      assert(s11 == s10 - 1);
      assert(s21 == s20 - 1);

      // Decrement number of edges
      --n_edges_;
    }

    // Remove node from i2u
    auto it = i2u.begin() + n.index(); 
    i2u.erase(std::remove(it, i2u.end(), *it), i2u.end()); // erase it
    
    // Adjust remaining indices in i2u
    for (; it != i2u.end(); ++it) {
      --(node_map_[*it].idx);
    }

    // decrement number of nodes
    --n_nodes_;
    
    return num_nodes();
  }

  /* @brief Removes a Node from Graph.
   * @param[in] node_iterator n_it, pointing to Node @n to be deleted        
   * @result    number of nodes remaining in Graph.
   *
   * @pre Node is in graph: 
   *         0 =< n.index() < g.num_nodes()
   *         for i s.t. i2u[i] == n.uid(), node_map[i2u[i]].index() == i
   *
   * @post Node @n is not in graph:
   *         there exists no i s.t. i2u[i] == n.uid()
   *
   * @post Graph size is reduced by one: 
   *        new g.num_nodes() == old g.num_nodes() - 1;
   *        new i2u.size()    == old i2u.size() - 1
   * 
   * @post Graph maintains a list of valid nodes:
   *        for all 0 =< i < new g.num_nodes(), g.node_map_[new i2u[i]].index()] == i
   *        for all 0 =< i < n.index(),
   *          old i2u[i]   == new i2u[i] 
   *        for all i =< n.index() < new g.num_nodes(),
   *          old i2u[i-1] == new i2u[i]
   *
   * @post for all i s.t i \in n.connections(0 =< j < n.degree()):
   *        new i.degree() == old i.degree() - 1
   *
   * @post new g.num_edges() == old g.num_edges() - n.degree();
   *        
   * @post Containers indexed by unique indices remain unaffected:
   *        new g.node_map_().size()   == old g.node_map_().size()
   *        new g.edge_map_().size()   == old g.edge_map_().size()
   *        new g.nodes_to_edge.size() == old nodes_to_edge.size()
   *        
   * Complexity is O(degree^2) for degree << n_nodes \sim n_edges
   * Invalidates no other nodes besides Node @n.
   * Invalidates no other edges besides Edge @e \in n.connections().
   * Invalidates i2u iterators for all i in i2u s.t. i =< n.index() < new g.num_nodes(),
   *
   */
  node_iterator remove_node(node_iterator n_it) {

    if (n_it != node_end()) {
      
      // Because node_map_ is std::unordered_map, erasing Node n = *n_it
      //    invalidates no other iterators. 
      //    Advancing before deletion avoids invalidation.

      Node n = *n_it;                      // Get the node pointed to by n_it
      ++n_it;                              // Advance iterator, avoiding invalidation.
      size_type n_nodes = remove_node(n);  // Delete from node_map_.

      (void) n_nodes;
      return n_it;
       
    } else {
      throw std::invalid_argument(
              "Cannot remove: node_map_.end() is not an Node");
    }
    
  }

  /* @brief Removes an Edge from Graph.
   * @param[in] n1, one end of Edge @e to be deleted
   * @param[in] n2, other end of Edge @e to be deleted
   * @result        number of edges remaining in Graph
   *
   * @pre @n1 and @n2 are in graph: 
   *         0 =< n{1,2}.index() < g.num_nodes()
   *         for i s.t. i2u[i] == n{1,2}.uid(), node_map[i2u[i]].index() == i
   *
   * @post Edge @e is not in graph:
   *         n{1,2}.uid() \not\in n{2,1}.connections().node2().uid()
   *
   * @post Graph size remains constant: new g.num_nodes() == old g.num_nodes()
   * @post if g.has_edge(@e), new n{1,2}.degree() == old n{1,2}.degree() - 1,
   *                     else new n{1,2}.degree() == old n{1,2}.degree().
   *
   * @post new g.num_edges()   == old g.num_edges() - 1.
   * @post if g.has_edge(@e), new g.num_edges() == old g.num_edges() - 1
   *                     else new g.num_edges() == old g.num_edges().
   *        
   * @post Containers indexed by unique indices remain unaffected:
   *        new g.edge_map_().size()   == old g.edge_map_().size()
   *        new g.nodes_to_edge.size() == old nodes_to_edge.size()
   *        
   * Complexity is O(degree^2) for degree << n_nodes \sim n_edges
   * Invalidates no nodes.
   * Invalidates no other edges besides Edge @e.
   * Invalidates iterators on n{1,2}.connections() that point to elements after @e.
   *
   */
  size_type remove_edge(const Node& n1, const Node& n2) {

    // If edge exists, create a proxy and delete all asosciated data.
    if (has_edge(n1, n2)) {
      Edge e = Edge(this, nodes_to_edge_[node_hash(n1, n2)]);
      delete_edge_data(e, n1, n2);
    } else {
      std::cout << "no edge between " << n1.uid() << " and  " << n2.uid() << std::endl;
    }

    return num_edges();
  }

  /* @brief Removes an Edge from Graph.
   * @param[in] e, the Edge to be removed. connects Nodes @n1 and @n2.
   * @result       number of edges remaining in Graph
   *
   * @pre @n1 and @n2 are in graph: 
   *         0 =< n{1,2}.index() < g.num_nodes()
   *         for i s.t. i2u[i] == n{1,2}.uid(), node_map[i2u[i]].index() == i
   *
   * @post Edge @e is not in graph:
   *         n{1,2}.uid() \not\in n{2,1}.connections().node2().uid()
   *
   * @post Graph size remains constant: new g.num_nodes() == old g.num_nodes()
   * @post if g.has_edge(@e), new n{1,2}.degree() == old n{1,2}.degree() - 1,
   *                     else new n{1,2}.degree() == old n{1,2}.degree().
   *
   * @post if g.has_edge(@e), new g.num_edges() == old g.num_edges() - 1
   *                     else new g.num_edges() == old g.num_edges().
   *        
   * @post Containers indexed by unique indices remain unaffected:
   *        new g.edge_map_().size()   == old g.edge_map_().size()
   *        new g.nodes_to_edge.size() == old nodes_to_edge.size()
   *        
   * Complexity is O(degree^2) for degree << n_nodes \sim n_edges
   * Invalidates no nodes.
   * Invalidates no other edges besides Edge @e.
   * Invalidates iterators on n{1,2}.connections() that point to elements after @e.
   *
   */
  size_type remove_edge(const Edge& e) {

    Node n1 = e.node1();
    Node n2 = e.node2();
    assert(has_edge(e.node1(), e.node2()));
    delete_edge_data(e, n1, n2);
    return num_edges();
  }

  /* @brief Removes an Edge from Graph.
   * @param[in] Edge iterator e_it pointing to @e, the Edge to be removed, 
   *               which connects Nodes @n1 and @n2.
   * @result    number of edges remaining in Graph
   *
   * @pre @n1 and @n2 are in graph: 
   *         0 =< n{1,2}.index() < g.num_nodes()
   *         for i s.t. i2u[i] == n{1,2}.uid(), node_map[i2u[i]].index() == i
   *
   * @post Edge @e is not in graph:
   *         n{1,2}.uid() \not\in n{2,1}.connections().node2().uid()
   *
   * @post Graph size remains constant: new g.num_nodes() == old g.num_nodes()
   * @post if g.has_edge(@e), new n{1,2}.degree() == old n{1,2}.degree() - 1,
   *                     else new n{1,2}.degree() == old n{1,2}.degree().
   *
   * @post if g.has_edge(@e), new g.num_edges() == old g.num_edges() - 1
   *                     else new g.num_edges() == old g.num_edges().
   *        
   * @post Containers indexed by unique indices remain unaffected:
   *        new g.edge_map_().size()   == old g.edge_map_().size()
   *        new g.nodes_to_edge.size() == old nodes_to_edge.size()
   *        
   * Complexity is O(degree^2) for degree << n_nodes \sim n_edges
   * Invalidates no nodes.
   * Invalidates no other edges besides Edge @e.
   * Invalidates iterators on n{1,2}.connections() that point to elements after @e.
   *
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    if (e_it != edge_end()) {
      
      // Because edge_map_ is std::unordered_map, erasing Edge e = *e_it
      //    invalidates no other iterators. 
      //    Advancing before deletion avoids invalidation.

      Edge e = *e_it;                      // Get the edge pointed to by e_it
      ++e_it;                              // Advance iterator, avoiding invalidation.
      size_type n_edges = remove_edge(e);  // Delete from edge_map_.

      (void) n_edges;
      return e_it;
       
    } else {
      throw std::invalid_argument(
              "Cannot remove: edge_map_.end() is not an Edge");
    }
  }

  /* @brief Removes an Edge from Graph.
   * @param[in] e  the Edge to be removed, 
   * @param[in] n1 first connection of @e
   * @param[in] n2 second connection of @e
   * @result    (void)
   *
   * @pre @n1 and @n2 are in graph: 
   *         0 =< n{1,2}.index() < g.num_nodes()
   *         for i s.t. i2u[i] == n{1,2}.uid(), node_map[i2u[i]].index() == i
   *
   * @pre @e connects @n1 and @n2:
   *       (e.node1 == n1.uid() && e.node2 == n2.uid() or
   *        e.node2 == n1.uid() && e.node1 == n1.uid())
   *       &&
   *        n1.uid() == e.node2().connections(i).node2().uid() for some 0 <= i < n2.degree() or
   *        n2.uid() == e.node1().connections(i).node1().uid() for some 0 <= i < n1.degree().
   *
   * @post Edge @e is not in graph:
   *         n{1,2}.uid() \not\in n{2,1}.connections().node2().uid()
   *
   * @post Graph size remains constant: new g.num_nodes() == old g.num_nodes()
   * @post new g.num_edges() == old g.num_edges() - 1
   * @post new n{1,2}.degree() == old n{1,2}.degree() - 1
   *
   * @post Containers indexed by unique indices remain unaffected:
   *        new g.edge_map_().size()   == old g.edge_map_().size()
   *        new g.nodes_to_edge.size() == old nodes_to_edge.size()
   *        
   * Complexity is O(degree^2) for degree << n_nodes \sim n_edges
   * Invalidates no nodes.
   * Invalidates no other edges besides Edge @e.
   * Invalidates iterators on n{1,2}.connections() that point to elements after @e.
   *
   */
  void delete_edge_data(const Edge& e, const Node& n1, const Node& n2) {
    
    size_type s10 = n1.degree();
    size_type s20 = n2.degree();
    delete_connection(n1, e); // Remove @e from @n1's connections list
    delete_connection(n2, e); // Remove @e from @n2's connections list
    size_type s11 = n1.degree();
    size_type s21 = n2.degree();
    assert(s11 == s10 - 1);
    assert(s21 == s20 - 1);

    // Key for edge_map_ is unique edge ID.
    edge_map_.erase(e.index());                 

    // Key for nodes_to_edge is hash of node indices.
    nodes_to_edge_.erase(node_hash(n1, n2));

    // Decrement number of edges.
    --n_edges_;                    
  }

  // deletes @e from @n.connections().
  // @pre (has_node(n) and has_edge(e)).
  void delete_connection(const Node& n, const Edge& e) {


    size_type s0 = n.degree();
    for (auto it = n.connections().begin(); it != n.connections().end(); ) {
      auto e_it = *it;
      if (e == e_it) {
        n.connections().erase(std::remove(it, n.connections().end(), *it),
                                n.connections().end());
        break;
      } else {
        ++it;
      }
    }
    size_type s1 = n.degree();
    assert(s1 == s0 - 1);
  }

 // ------------------------- Miscellaneous Functions ----------------------- //

  // Hashes two Node IDs into a unique key for nodes_to_edge_ map.
  std::size_t node_hash(const Node& n1, const Node& n2) const {

      // Make sure the IDs are sorted before we query for their hash.
      size_type id1 = std::min(n1.uid(), n2.uid()); // uid √
      size_type id2 = std::max(n1.uid(), n2.uid()); // uid √

      std::size_t seed = 0;
      boost::hash_combine(seed, id1);
      boost::hash_combine(seed, id2);

      return seed;
  }

 // -------------------------- Private Graph Methods ------------------------ //

 private:

  // Data structures for nodes and edges
  struct node_element {
    Node node;        // A node object this container is connected to.
    size_type uid;    // ID unique to this node, immutable
    size_type idx;    // Position of node in graph's node container, mutable

    std::vector<edge_type> connections;  // Edges connecting Node to others
    Point                  position;     // position in space
    node_value_type        value;        // Can be arbitrary data structure

  };

  struct edge_element {
    Edge            edge;  // Pointer back to the Edge.
    size_type       node1; // uid of Node 1
    size_type       node2; // uid of Node 2
    edge_value_type value; // struct containing physical parameters

  };

  // Containers for nodes and edges. 
  node_map_type node_map_;    // map unique index to a node_element
  edge_map_type edge_map_;    // map edge_index to edge_element
  std::vector<size_type> i2u; // map graph index to unique index

  // Map Node IDs to Edge ID. 
  // Allows for quick check of whether an Edge exists (i.e. for has_edge).
  std::unordered_map<std::size_t, size_type> nodes_to_edge_;

  // count of each element.
  size_type n_nodes_  = 0;           // accessible via num_nodes()
  size_type n_edges_  = 0;           // accessible via num_edges()
  size_type n_nodes_generated_ = 0;  // used to assign unique IDs to new nodes
  size_type n_edges_generated_ = 0;  // used to assign unique IDs to new edges

  // Forbid copy or assignment of Graphs
  Graph(const Graph&) = delete;
  Graph& operator=(const Graph&) = delete;

  // Create a graph index for comparison of edges between graphs
  int g_idx_;

  // Create global ordering of graph
  bool operator<(const Graph& g)  {
    return this->g_idx_ < g.g_idx_;
  }


};


#endif // CME212_GRAPH_HPP
