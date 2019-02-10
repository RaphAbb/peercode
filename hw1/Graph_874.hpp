#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
#include <iostream>

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

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_value_type = V;
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
  Graph() : 
    nodes_(), 
    edges_(), 
    next_uid_(0), 
    node_indices(), 
    edge_map(),
    adjacency_list() 
    { 
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
    Node() : graph_(nullptr), uid_(size_type(-1)) {
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      size_type i = index();
      return graph_->nodes_[i].position;
      // if (graph_ != nullptr)
      //     if (i < graph_->size())
      //         return graph_->nodes_[i].position;
      // return Point(); // okay, how to change?
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      // TODO: speed up
      return graph_->node_indices.at(uid_); // just let it break if node is invalid
      // if (graph_ != nullptr)
      //   return graph_->node_indices.at(uid_);
      // return size_type(-1);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // [x] node_value_type& value();
    // [x] const node_value_type& value() const;
    // [x] size_type degree() const;
    // [ ] incident_iterator edge_begin() const;
    // [ ] incident_iterator edge_end() const;

    /** Return a reference to this node's value **/
    node_value_type& value() {
      return graph_->nodes_[index()].val;
    }

    /** Return a const reference to this node's value **/
    const node_value_type& value() const {
      return graph_->nodes_[index()].val;
    }

    /** Return the degree of this node **/
    size_type degree() const {
      return graph_->adjacency_list[index()].size();
    }

    /** Return the start incident_iterator for this node**/
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, uid_, graph_->adjacency_list[index()].begin());
    }

    /** Return the end incident_iterator for this node**/
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, uid_, graph_->adjacency_list[index()].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (index() == n.index()) // check that index is the same
          if (graph_ == n.graph_) // check that graph is the same
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
      return (uid_ < n.uid_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    Graph* graph_;    // points to graph
    size_type uid_;         // Unique identification number of element
    // Private constructor
    Node(const Graph* graph, size_type uid) : graph_(const_cast<graph_type*>(graph)), uid_(uid){}

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes_.size();
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
    node_value_type myval = val;
    internal_node i_node(position,myval,next_uid_);

    nodes_.push_back(i_node); // add to nodes_
    edge_map[i_node.uid] = std::map<size_type,size_type>();
    // adjacency_map[i_node.uid] = std::vector<size_type>();
    adjacency_list.push_back(std::vector<size_type>());
    node_indices[i_node.uid] = num_nodes() - 1; // add to node_indices
    Node node(this,i_node.uid);
    
    ++next_uid_;
    return node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (this == n.graph_);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i < nodes_.size())
        return Node(this,nodes_[i].uid);
    return Node();
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
    Edge() : node1_(Node()), node2_(Node()) {}

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
        return node1_;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
        return node2_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
        return (((node1_== e.node1_) && (node2_ == e.node2_)) || ((node1_== e.node2_) && (node2_ == e.node1_)));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (node1_ < e.node1_) return true;
      else if (node2_ < e.node2_) return true;
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Node node1_;
    Node node2_;
    Edge(Node n1, Node n2) : node1_(n1), node2_(n2) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(Node(this,edges_[i].node1_uid),Node(this,edges_[i].node2_uid));
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    auto& vec = adjacency_list[a.index()]; // vector of edges emanating from a
    for (auto it = vec.begin(); it != vec.end(); it++){
      if (*it == b.uid_) return true;
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
    if (has_node(a) && has_node(b)){
      if (has_edge(a,b)) return Edge(a,b); // return edge if it already exists
      // Otherwise add edge
      edges_.push_back(internal_edge(a.uid_,b.uid_));
      // add edge (both orientations) to edge_map (they share the same index)
      edge_map[a.uid_][b.uid_] = edges_.size() - 1;
      edge_map[b.uid_][a.uid_] = edges_.size() - 1;
      // add edge (both orientations) to adjacency_list
      adjacency_list[a.index()].push_back(b.uid_);
      adjacency_list[b.index()].push_back(a.uid_);
    }
    return Edge(a,b);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    edges_.clear();
    node_indices.clear();
    edge_map.clear();
    adjacency_list.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
      graph_ == nullptr;
      idx_ = 0;
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Get Node corresponding to this iterator **/
    Node operator*() const {
      return graph_->node(idx_); // return a node matching the iterator idx
    }
    /** Increment iterator **/
    NodeIterator& operator++() {
      ++idx_;
      return *this;
    }
    /** Equality comparison between this iterator and another NodeIterator 
     * @param[in] iter    the other node_iterator 
     * **/
    bool operator==(const NodeIterator& iter) const {
      return (idx_ == iter.idx_);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    graph_type* graph_; // pointer to graph
    size_type idx_; // index into graph_->nodes_ (the graph's node vector)
    NodeIterator(const graph_type* graph, size_type idx) : graph_(const_cast<graph_type*>(graph)), idx_(idx) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Return start iterator pointing to the first node in the graph **/
  node_iterator node_begin() const {
    return NodeIterator(this,0); // 0 is the initial index
  }
  /** Return end iterator pointing just past the last node in the graph **/
  node_iterator node_end() const {
    return NodeIterator(this,this->size()); // size of node vector
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() : graph(nullptr) {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** Return Edge corresponding to this iterator **/
    edge_type operator*() const {
      size_type dst_id = *it;
      return graph->edge(src_id,dst_id); // uses private graph helper method
    }

    /** Increment the iterator **/
    incident_iterator& operator++() {
      ++it;
      return *this;
    }

    /** Check equality with another incident_iterator @a iter 
     * @param[in] iter      the other incident_iterator
     * **/
    bool operator==(const incident_iterator& iter) const {
      return it == iter.it;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    graph_type* graph;
    size_type src_id;   // gives the uid_ of the node
    std::vector<size_type>::iterator it; // an iterator into graph->adjacency_list[src.index()]
    IncidentIterator(const graph_type* g, size_type s, std::vector<size_type>::iterator i) : 
      graph(const_cast<graph_type*>(g)), src_id(s), it(i) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() : graph_(nullptr) {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /** get edge element corresponding to this iterator 
     * **/
    edge_type operator*() const {
      // return graph_->edges_[idx_];
      return graph_->edge(idx_);
    }
    /** Increment edge iterator 
     * Since only unique edges are stored in graph_->edges, we can iterate straight through
     * **/
    edge_iterator& operator++(){
      ++idx_;
      return *this;
    }
    /** Check equality betwen edge_iterators **/
    bool operator==(const edge_iterator& iter) const {
      return idx_ == iter.idx_;
    }


   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type *graph_;
    size_type idx_;
    EdgeIterator(const graph_type* g, size_type i) : graph_(const_cast<graph_type*>(g)), idx_(i) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Return begin edge_iterator pointing to first edge **/
  edge_iterator edge_begin() const{
    return EdgeIterator(this,0);
  }
  /** Return end edge_iterator pointing past last edge **/
  edge_iterator edge_end() const{
    return EdgeIterator(this,edges_.size());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct internal_node {
    const Point position;
    node_value_type val;
    size_type uid;
    internal_node() : position(), val(), uid() {}
    internal_node(const Point& pos, size_type id) : position(pos), val(), uid(id) {}
    internal_node(const Point& pos, node_value_type v, size_type id) : position(pos), val(v), uid(id) {}
  };
  struct internal_edge {
    size_type node1_uid;
    size_type node2_uid;
    internal_edge() : node1_uid(), node2_uid() {}
    internal_edge(size_type u1, size_type u2) : node1_uid(u1), node2_uid(u2) {}
  };
  // nodes of graph
  std::vector<internal_node> nodes_;
  // edges of graph
  std::vector<internal_edge> edges_;
  // next unique id to be assigned to a node
  size_type next_uid_;
  // map of node indices (a.uid_ => a.index)
  std::map<size_type,size_type> node_indices; 
  // map of edges indices
  std::map<size_type,std::map<size_type,size_type>> edge_map; 
  // Adjacency list
  std::vector<std::vector<size_type>> adjacency_list;
  /** Helper method to get an edge directly from node uids 
   * @param[in] src_id      uid of first node
   * @param[in] dst_id      uid of second node
   * **/
  Edge edge(size_type src_id, size_type dst_id) {
    return Edge(Node(this,src_id),Node(this,dst_id));
  }
};

#endif // CME212_GRAPH_HPP
