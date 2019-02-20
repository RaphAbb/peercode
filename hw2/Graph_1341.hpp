#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <set>
#include <string>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = double, typename E = double>
class Graph {

 private:

  //
  // PRIVATE PREDECLARATIONS
  //

  struct internal_node;
  struct internal_edge;

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
  /** Type of node value */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  /** Type of edge value */
  using edge_value_type = E;

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
  Graph()
    : nodes_(), size_(0), edges_(), edge_size_(0) {
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
    Node()
      : graph_(nullptr), id_() {
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_[id_].position;
    }

    /** Returns modifiable node's position. */
    Point& position() {
      return graph_->nodes_[id_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return id_;
    } 
    /** Returns the value of the node */
    node_value_type& value() {
      return graph_->nodes_[id_].value;
    }

    /** Overloading const version */
    const node_value_type& value() const {
      return graph_->nodes_[id_].value;
    }

    // Number of nodes this node is connected to O(1)
    size_type degree() const {
      return graph_->nodes_[id_].connectivity.size();
    };

    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, id_, 0);
    };

    incident_iterator edge_end() const {
      return IncidentIterator(graph_, id_, degree());
    };


    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (id_ == n.id_ && graph_ == n.graph_)
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
      if (id_ < n.id_)
        return true;
      return false;
    }


   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Pointer back to the Graph container
    Graph* graph_;
    // This node's index number
    size_type id_;

    // Private Constructor
    Node(const Graph* graph, size_type id)
      : graph_(const_cast<Graph*>(graph)), id_(id) {
    }

  };

  //
  // END NODES
  //

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return size_;
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type(),
      const std::vector<size_type> connectivity = std::vector<size_type>()) {
    internal_node temp_node { position, value, connectivity};
    nodes_.push_back(temp_node);
    ++size_;
    return Node(this, size_-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.id_ < size_ && n.graph_ == this)
      return true;
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < size_);
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
    Edge()
      : graph_(nullptr), id_(0) {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->edges_[id_].node1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->edges_[id_].node2;
    }

    /** Returns the length between 2 nodes */
    double length() const {
      Node n1 = graph_->edges_[id_].node1;
      Node n2 =  graph_->edges_[id_].node2;
      return norm(n2.position() - n1.position());
    }

    /** Returns the index of Edge */
    size_type index() const {
      return id_;
    }

    /** Return value of edge */
    edge_value_type& value() {
      return graph_->edges_[id_].value;
    }

    /** Const version of returning value of edge*/
    const edge_value_type& value() const {
      return graph_->edges_[id_].value;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (id_ == e.id_ && graph_ == e.graph_)
        return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (id_ < e.id_)
        return true;
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    
    Edge(const Graph* graph, size_type id)
      : graph_(const_cast<Graph*>(graph)), id_(id) {
    }

    // Pointer back to graph container
    Graph* graph_;
    // ID of the edge
    size_type id_;

  };

  // 
  // END EDGES
  //

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_size_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert(a.index() < num_nodes() && b.index() < num_nodes());
    // Return false if nodes are invalid
    if (a.index() == b.index()) {
      return false;
    }

    std::vector<size_type> connected_edge_ids = nodes_[a.index()].connectivity;

    for (auto edge_id : connected_edge_ids) {
      Edge e = Edge(this, edge_id);
      if (e.node1() == b or e.node2() == b)
        return true;
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
    // Check if nodes are valid
    assert(a.index() < num_nodes() && b.index() < num_nodes());
    assert(a.index() != b.index());

    // Check if edge already exists
    if ( has_edge(a, b) ) {
      std::vector<size_type> connected_edge_ids = nodes_[a.index()].connectivity;
      for (auto edge_id : connected_edge_ids) {
        Edge e = Edge(this, edge_id);
        if (e.node1() == b or e.node2() == b)
          return e;
      }
    }

    // Create new edge
    edges_.emplace_back(internal_edge(a, b));

    // Update connectivity
    nodes_[a.index()].connectivity.emplace_back(edge_size_);
    nodes_[b.index()].connectivity.emplace_back(edge_size_);

    // Update num edges
    ++edge_size_;
    return Edge(this, edge_size_-1);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    edge_size_ = 0;
    size_ = 0;
  }


  // 
  // Removal Functions
  //
  
  /** Remove edge e from this graph
   * @param[in] e     The edge that we want to remove from the graph
   * returns true if edge exists and removed, false if edge doesn't exist
   *
   * note: Replaces the removed edge with the edge with the largest id
   * note: define @a e_last as the edge with the largest index
   *
   * @pre @a e.index() < num_edges()
   * @post new num_edges() = old num_edges() - 1
   * @post new edges_.size() = old edges_.size() - 1
   * @post new edges_[@a e.index()] = old edges_[@a e_last.index()]
   * @post for i, 0 <= i < num_nodes()
   *        if nodes_[i].connectivity contains edge @a e, delete it
   *        if nodes_[i].connectivity contains edge @a e_last, rename the index
   *          of edge to be that of edge @a e.
   *          
   */
  bool remove_edge(const Edge& e) {
    // Checks whether edge is valid
    if (e.index() >= num_edges()){
      return false;
    }

    size_type removed_edge_id = e.index();
    size_type replaced_edge_id = num_edges()-1;

    // Update nodes_ by
    // 1. Removing all mentions of the removed edge
    // 2. Updating the new edge
    for (auto ni = node_begin(); ni != node_end(); ++ni) {
      Node n = *ni;
      std::vector<size_type>& edge_ids = nodes_[n.index()].connectivity;
      for (auto ei = edge_ids.begin(); ei != edge_ids.end(); ) {
        if (*ei == removed_edge_id) { // Removing mention of the removed edge
          ei = edge_ids.erase(ei);
          continue;
        } else if (*ei == replaced_edge_id) { // Renaming the replaced edge
          *ei = removed_edge_id;
        }
        ++ei;
      }
    } 

    // Remove from edges_
    edges_[removed_edge_id] = edges_[replaced_edge_id]; // Replace with last one first
    edges_.pop_back(); 

    // Reduce edge_size_
    edge_size_--;
    return true;
  }

  /** Remove edge from this graph
   * @param[in] e_it    Iterator of edge to remove
   * returns Iterator of the edge that replaced the removed edge
   *
   * note: see remove_edge(Edge e) for complete specification
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return e_it;
  }

  /** Remove edge from this graph
   * @param[in] e_it    Iterator of edge to remove
   * returns Iterator of the edge that replaced the removed edge
   *
   * note: see remove_edge(Edge e) for complete specification
   */
  incident_iterator remove_edge(incident_iterator e_it) {
    remove_edge(*e_it);
    return e_it;
  }


  /** Remove edge between node @a a and node @a b from this graph
   * @param[in] a       Node A
   * @param[in] b       Node B
   * returns true if edge exists and removed, false if edge doesn't exist
   *
   * note: see remove_edge(Edge e) for complete specification
   */
  bool remove_edge(const Node& a, const Node& b) {
    // Find edge associated with the nodes
    Edge edge;
    std::vector<size_type> connected_edge_ids = nodes_[a.index()].connectivity;
    for (auto edge_id : connected_edge_ids) {
      Edge e = Edge(this, edge_id);
      if (e.node1() == b or e.node2() == b)
        edge = e;
    }
    if ( edge == Edge() )
      return false;

    return remove_edge(edge);
  }


  /** Removes node n from this graph
   * @param[in] n       Node to be removed
   * returns true if node exists and is removed, false if node doesn't exist
   *
   * note: Replaces the removed node with the node with the largest id
   * note: define @a n_last as the node with the largest index
   * note: deletes all edges incident to the node
   *
   * @pre @a n.index() < num_nodes()
   * @post new num_nodes() = old num_nodes() - 1
   * @post new nodes_.size() = old nodes_.size() - 1
   * @post new nodes_[@a n.index()] = old nodes_[@a n_last.index()]
   * @post for i, 0 <= i < old nodes_[@a n.index()].connectivity.size()
   *          delete edge with index i
   * @post for i, 0<= i < old nodes_[@a n_last.index()]
   *          update edges_[i] to replace the node @a n_last with @a n
   *
   */
  size_type remove_node(const Node& n) {
    // Checks whether node is valid
    size_type removed_node_id = n.index();
    if ( !(removed_node_id < num_nodes()) ) {
      return false;
    }
    size_type replaced_node_id = num_nodes()-1;


    // Removing all edges which are adjacent to this node
    int countess = 0;
    for (auto ei = n.edge_begin(); ei != n.edge_end();) {
      ei = remove_edge(ei);
      countess++;
    }

    // Updating the edge indices of the replaced node
    auto edge_ids = nodes_[replaced_node_id].connectivity;
    for (auto edge_id: edge_ids) {
      Edge e = Edge(this, edge_id);
      if (e.node1().index() == replaced_node_id) {
        edges_[edge_id].node1 = n;
      } else if (e.node2().index() == replaced_node_id) {
        edges_[edge_id].node2 = n;
      }
    }

    // Update nodes_
    nodes_[removed_node_id] = nodes_[replaced_node_id]; // Replace with last
    nodes_.pop_back();

    // Reduce the number of nodes
    size_--;
    return removed_node_id;
  }

  /** Removes node from this graph
   * @param[in] n_int       Iterator of node to remove
   * returns iterator of the replacement node
   *
   * see remove_node(Node n) for the complete specification
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return n_it;
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

    Node operator*() const {
      return Node(graph_, p_);
    }

    NodeIterator& operator++() {
      p_++;
      return *this;
    }

    bool operator==(const node_iterator& x) const {
      return p_ == x.p_;
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type p_;
    NodeIterator(const Graph* graph, size_type p) 
      : graph_(const_cast<Graph*>(graph)), p_(p) {
    }
  };

  //
  // END NODE ITERATOR
  //

  node_iterator node_begin() const { return NodeIterator(this, 0); }
  node_iterator node_end() const { return NodeIterator(this, num_nodes()); }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator> {
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

    Edge operator*() const {
      // Gets the adjacent node id from the iterator
      size_type edge_id = graph_->nodes_[h_node_id_].connectivity[edge_iter_];
      Edge e = Edge(graph_, edge_id);

      // If Node1 isn't the "home node", swap the ordering such that it is
      if (e.node1().index() != h_node_id_) {
        Node tmp_node = graph_->edges_[edge_id].node1;
        graph_->edges_[edge_id].node1 = graph_->edges_[edge_id].node2;
        graph_->edges_[edge_id].node2 = tmp_node;
      }

      return e;
    }

    IncidentIterator& operator++() {
      edge_iter_++;
      return *this;
    }

    bool operator==(const IncidentIterator& x) const {
      return edge_iter_ == x.edge_iter_;
    }


   private:
    friend class Graph;

    Graph* graph_;
    size_type h_node_id_; // ID of the central node
    size_type edge_iter_; // Indexing into nodes_[h_node_id_].connectivity

    // Private constructor to make IncidetIterator tied to a certain graph and home node
    IncidentIterator(const Graph* graph, size_type h_node_id, size_type edge_iter) 
      : graph_(const_cast<Graph*>(graph)), h_node_id_(h_node_id), 
        edge_iter_(edge_iter) {
    }
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

    Edge operator*() const {
      return Edge(graph_, edge_id_);
    };

    EdgeIterator& operator++() {
      edge_id_++;
      return *this;
    };

    bool operator==(const EdgeIterator& x) const {
      return edge_id_ == x.edge_id_;
    };

   private:
    friend class Graph;
    Graph* graph_;
    size_type edge_id_;

    EdgeIterator(const Graph* graph, size_type edge_id)
      : graph_(const_cast<Graph*>(graph)), edge_id_(edge_id) {
      };
  };

  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  };

  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  };

 private:

  /* 
   * Node helper functions, data members, etc.
   */

  struct internal_node {
    Point position;
    node_value_type value;
    std::vector<size_type> connectivity; // Edge ids connected to this node
  };

  std::vector<internal_node> nodes_;

  size_type size_;

  /* 
   * Edge helper functions, data members, etc.
   */


  struct internal_edge {
    Node node1;
    Node node2;
    edge_value_type value;
    internal_edge(Node a, Node b): 
      node1(a), node2(b), value() {
    }
  };
  std::vector<internal_edge> edges_;

  size_type edge_size_;



};


#endif // CME212_GRAPH_HPP

