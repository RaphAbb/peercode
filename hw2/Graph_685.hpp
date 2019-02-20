#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file  Graph.hpp
 *  @brief Undirected graph type. 
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_set>
#include <cstdint>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 *  @brief A template for 3D undirected graphs.
 *
 *  Users can add and retrieve nodes and edges. 
 *  Edges are unique (there is at most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
 private:
  // Pre-declarations for internal types.
  struct node_element; 
  struct edge_element;
  struct edge_hash; 

 public:
  // Pre-declarations for public types.
  class Node;               // Node type 
  class NodeIterator;       // Iterator for graph nodes
  class Edge;               // Edge type
  class EdgeIterator;       // Iterator for graph edges
  class IncidentIterator;   // Iterator over incident edges to a node

  // Type synonyms.
  using graph_type = Graph<V,E>;
  using size_type = unsigned;

  using node_type = Node;
  using node_value_type = V;
  using node_iterator = NodeIterator; 

  using edge_type = Edge;
  using edge_value_type = E; 
  using edge_iterator = EdgeIterator;
  using incident_iterator = IncidentIterator; 

  /** Construct an empty graph (default constructor)  
   *  @return  A new Graph object with zero nodes and zero edges.
   */ 
  Graph() {} 

  /** Graph uses default destructor. 
   *  @post  All of this graph's nodes and edges are invalidated.
   */
  ~Graph() = default;

  /** @class Graph::Node
   *  @brief Class representing the graph's nodes.
   *
   *  Node objects are used to access information about the Graph's nodes.
   *  Each Node object is lightweight, storing only:
   *     @a graph_ : A pointer to the Graph which contains it
   *     @a uid_   ; A unique node ID number
   */
  class Node : private totally_ordered<Node>{
   public:
    /** Construct an invalid node. */
    Node(): graph_(NULL) { } 

    /** Return this node's index, a number in the range [0, graph_->size()). 
     *  @return  The index of this node in the graph. 
     *  @pre     Node is valid, i.e. this.valid() == true.
     *  
     *  Complexity: O(1)
     */
    size_type index() const {
      assert(valid()); 
      return graph_->nodeList_[uid_].idx;
    }

    /** Read/Write access to the position of this node.
     *  @pre     Node is valid (this.valid() == true) 
     *  @return  Reference to a Point object with this node's position. 
     *
     *  Complexity: O(1) 
     */
    Point& position() {
      assert(valid());
      return (const_cast<Graph*>(graph_))->nodeList_[uid_].p;
    }

    /** Read-only access to the position of this node. 
     *  @pre     Node is valid (this.valid() == true) 
     *  @return  Const reference to a Point object with this node's position. 
     *  
     *  Complexity: O(1)
     */
    const Point& position() const {
      assert(valid()); 
      return graph_->nodeList_[uid_].p;
    }

    /** Read/Write access to value of this node. 
     *  @pre     Node is valid (this.valid() == true) 
     *  @return Reference to the value of this node. 
     *
     *  Complexity: O(1)  
     */
    node_value_type& value() { 
      assert(valid()); 
      return (const_cast<Graph*>(graph_))->nodeList_[uid_].val;
    }
    
    /** Read-only access to value of this node. 
     *  @pre     Node is valid (this.valid() == true) 
     *  @return Const reference to the value of this node. 
     *
     *  Complexity: O(1)
     */
    const node_value_type& value() const { 
      assert(valid()); 
      return graph_->nodeList_[uid_].val;
    }

    /** Return the degree (num. of incident edges) to this node.
     *  @pre     Node is valid (this.valid() == true) 
     *  @return  The number of distinct edges in graph_ which contain this node.
     *
     *  Complexity: O(1)
     */
    size_type degree() const {
      assert(valid()); 
      return graph_->nodeList_[uid_].adj_list.size();
    } 
 
    /** Return an iterator to the incident edges of this node.
     *  @return An iterator pointing to the first element in 
     *          this node's adjacency list.  
     */
    incident_iterator edge_begin() const {
      return IncidentIterator( graph_->nodeList_[uid_].adj_list.begin() );
    }
 
    /** Return the end iterator to the incident edges of this node. 
     *  @return An iterator pointing to one past the last element
     *          in this node's adjacency list. 
     */
    incident_iterator edge_end() const {
      return IncidentIterator( graph_->nodeList_[uid_].adj_list.end() );
    }

    /** Test whether this node and @a n are equal.
     *  @param[in] n   Node to compare this to. 
     *  @return        True if @a n and this have the same
     *                 private members graph_ and uid_. 
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_ && uid_ == n.uid_ );
    }

    /** Test whether this node is less than @a n in a global order.
     *  @param[in] n   Node to compare this to.
     *
     *  @post For any two nodes x and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const { 
      if (graph_ != n.graph_) { 
        auto g1 = reinterpret_cast<std::uintptr_t>(graph_);
        auto g2 = reinterpret_cast<std::uintptr_t>(n.graph_); 
        return g1 < g2;
      } else { 
        return uid_ < n.uid_; 
      } 
    }

   public:
    friend class Graph;        // Allow Graph to access Node's private members/functions.
    const graph_type* graph_;  // Pointer to graph this node is part of.
    size_type uid_;            // Unique ID number of this node.

    /** Private Node constructor, should only be called by Graph 
     *  @param graph  Pointer to the graph which contains the new node
     *  @param nid    Internal node id assigned by Graph.  
     *  @return       A new Node object, associated with @a graph and @a nid. 
     */
    Node(const graph_type* graph, size_type id) 
      : graph_(graph), uid_(id) { 
    }

    /** Checks if this node object satisfies its representation invariants
     *  @return   True if the node is valid, False if the RI has been violated
     *
     *  Checks:
     *     1. uid_ is in bounds of of nodeList_[uid_] 
     *     2. The idx_ associated with this uid_ is less than graph_->num_nodes()
     *     3. idx2uid_ correctly maps this node's uid_ to its idx_ in nodeList_. 
     */
    bool valid() const { 
      return uid_ >= 0 && uid_ < graph_->nodeList_.size()         
          && graph_->nodeList_[uid_].idx < graph_->idx2uid_.size()   
          && graph_->idx2uid_[graph_->nodeList_[uid_].idx] == uid_;              
    } 

  };

  /** Return the number of nodes in the graph.
   *
   *  Complexity: O(1).
   */
  size_type size() const {
    return idx2uid_.size();
  }

  /** Synonym for size(). 
   *
   *  Complexity: O(1).
   */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   *  @param[in] position  The new node's position. 
   *  @param[in] value     The value associated with the new node. 
   *  @return              The new node.
   * 
   *  @post new num_nodes() == old num_nodes() + 1 
   *  @post For 0 <= i < size(), nodeList_[idx2uid_[i]].idx == i
   *  
   *  Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // Initialize the new node
    node_element newNode; 
    newNode.idx = size();
    newNode.p = position;
    newNode.val = value;
    size_type new_uid = nodeList_.size();  
    newNode.n = Node(this, new_uid);  

    // Update Graph's data structures
    nodeList_.push_back( newNode ); // Copies newNode into nodeList_
    idx2uid_.push_back( new_uid );  // Update index-uid mapping
 
    return newNode.n; // Returns a copy of newly constructed node
  }

  /** Determine if a Node belongs to this Graph
   *  @return True if @a n is currently a Node of this Graph
   *
   *  Complexity: O(1).
   */
  bool has_node(const Node& n) const { 
    return (n.graph_ == this && n.valid()); 
  }

  /** Return the node with index @a i.
   *  @param i  The index at which to retrieve the node. 
   *  @return   The node at index @a i. 
   * 
   *  @pre 0 <= @a i < num_nodes()
   *  @post @ return.index() == i
   *
   *  Complexity: O(1).
   */
  Node node(size_type i) const {
    assert( i >= 0 && i < num_nodes());
    return nodeList_[idx2uid_[i]].n;
  }

  /** Removes this node and all of its edges from the graph.
   *  @param[in] n  The node to remove
   *  @return       1 if node is successfully removed, 0 otherwise
   * 
   *  @pre   @a n is a valid node. 
   *  @post  All nodes m, where m == @a n, are invalidated.
   *  @post  All edges e, where (e.node1() == @a n || e.node2() == @a n),
   *         are invalidated. 
   *
   *  Complexity: O(n), n = graph.size().
   */
  size_type remove_node(const Node& n) { 
    // Check if node is in graph:  O(1)
    if ( !has_node(n) ) { std::cout << "a" << std::endl; return 0; } 

    // Remove edges incident to node: O(deg(n)) 
    while( n.degree() > 0 ) {
      remove_edge( *(n.edge_begin()) );
    }    

    // Invalidate the removed node: O(1) 
    size_type ind = n.index();
    nodeList_[n.uid_].n.graph_ = NULL;

    // Update node indexing: O(n)
    idx2uid_.erase( idx2uid_.begin() + ind );
    size_type count = 0;
    for (auto it = idx2uid_.begin(); it != idx2uid_.end(); ++it) { 
      nodeList_[*it].idx = count++;
    }
    
    return 1;
  } 

  /** Iterator form of remove node.
   *  @param[in] n_it  Iterator to the node to be removed.
   *  @return          The 'next' iterator ( ++@a n_it) 
   *
   *  @pre   @n_it points to a valid node.
   *  @post  All nodes m, where m == @a n, are invalidated.
   *  @post  All edges e, where (e.node1() == @a n || e.node2() == @a n),
   *         are invalidated. 
   * 
   *  Complexity: O(n), n = graph.size() .
   */
  node_iterator remove_node(node_iterator n_it) { 
    remove_node( *(n_it) ); 
    return ++n_it;
  } 

  /** @class Graph::Edge
   *  @brief Class representing the graph's edges.
   *
   *  Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   *  are considered equal if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Edge> {
   public:
    /** Create an invalid edge */ 
    Edge(): graph_(NULL) {}

    /** Return this first node of this edge 
     */ 
    Node node1() const {
      assert(valid());
      if (flipped_) { 
        return graph_->edgeList_[uid_].n2; 
      } else {
        return graph_->edgeList_[uid_].n1; 
      } 
    }

    /** Return the other node of this Edge 
     */
    Node node2() const {
      assert(valid()); 
      if (flipped_) { 
        return graph_->edgeList_[uid_].n1; 
      } else { 
        return graph_->edgeList_[uid_].n2; 
      }
    }

    /** Read/Write access to value of this edge. 
     *  @return Reference to the value of this edge. 
     *
     *  Complexity: O(1)
     */
    edge_value_type& value() {
      assert(valid()); 
      return (const_cast<Graph*>(graph_))->edgeList_[uid_].val;
    }
    
    /** Read-only access to value of this edge. 
     *  @return Const reference to the value of this edge. 
     *  @post   The value of this edge is unchanged.
     *
     *  Complexity: O(1)
     */
    const edge_value_type& value() const { 
      assert(valid()); 
      return graph_->edgeList_[uid_].val;
    }

    /** Test whether this edge and @a e are equal.
     *
     *  Note that two edges are equal even if their nodes are 'flipped'.
     *  The order of nodes doesn't matter (undirected edges). 
     */
    bool operator==(const Edge& e) const {
      assert(valid()); 
      return (graph_ == e.graph_ && uid_ == e.uid_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     *  This ordering function is useful for STL containers such as
     *  std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_ != e.graph_) { 
        auto g1 = reinterpret_cast<std::uintptr_t>(graph_);
        auto g2 = reinterpret_cast<std::uintptr_t>(e.graph_); 
        return g1 < g2;
      } else { 
        return uid_ < e.uid_; 
      } 
    }

   private:
    friend class Graph;       // Allow Graph to access Edge's private members/functions
    const graph_type* graph_; // Pointer to the graph this edge is a part of. 
    size_type uid_;           // Unique ID number of this edge object
    bool flipped_;            // True if node1() > node2() (tracks order) 

    /** Private Edge Constructor for use by Graph. 
     *  @param graph   Pointer to the graph to which this edge belongs.
     *  @param uid     Unique ID number for new edge. 
     *  @param n1      The first node in the edge.
     *  @param n2      The second node in the edge.
     *  @return        A new Edge object in @a graph. 
     * 
     *  @post this.node1() = @a n1, this.node2() = @a n2.
     */ 
    Edge(const graph_type* graph, size_type uid, const Node& n1, const Node& n2) 
      : graph_(graph), uid_(uid), flipped_(n1 > n2) {} 

    /** Checks if this edge object satisfies its representation invariants
     *  @return   True if the edge is valid, False if the RI has been violated
     *
     *  Checks:
     *     1. uid_ is in bounds of of edgeList_[uid_]
     *
     *  Note the validity checks for edge are less stringent than for nodes,
     *  because we do not impose the index invariant on edge. 
     */ 
    bool valid() const { 
      return graph_!=NULL && uid_ >= 0 && uid_ < graph_->edgeList_.size();         
    } 
  
  };

  /** Return the total number of edges in the graph.
   *
   *  Complexity: O(1)
   */
  size_type num_edges() const {
    return edgeMap_.size();
  }

  /** Return the edge with index @a i.
   *  @pre 0 <= @a i < num_edges()
   *
   *  Complexity: O(E), E is the number of edges in the graph.
   *  I have sacrificed run time on this function
   *  to speed up remove operations, since the EdgeIterator
   *  can be used for efficient traversal (e.g. printing all edges).
   */
  Edge edge(size_type i) const {
    assert( i >= 0 && i < num_edges() ); 
    auto it = edge_begin();
    for (size_type ind = 0; ind < i; ++ind) { ++it; } 
    return *it;
  }

  /** Test whether two nodes are connected by an edge.
   *  @pre @a a and @a b are distinct valid nodes of this graph
   *  @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   *  Complexity: O( a.degree() )
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert( has_node(a) && has_node(b) ); 
    // Search node a's adjacency list for node b
    for (auto it = a.edge_begin(); it != a.edge_end(); ++it) { 
      if ((*it).node2() == b) { return true; } 
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
   * Complexity: Same as has_edge, O( a.degree() )
   */
  Edge add_edge(const Node& a, const Node& b) {
    assert( has_node(a) && has_node(b) ); 
    // Return current edge if it already exists 
    for (auto it = a.edge_begin(); it != a.edge_end(); ++it) { 
      if ((*it).node2() == b) { 
        return *it; 
      }
    } 
 
    // No existing edge was found. Make a new one. 
    size_type new_uid = edgeList_.size(); 

    // Update the Adjacency Lists
    nodeList_[a.uid_].adj_list.insert( Edge(this,new_uid,a,b) );
    nodeList_[b.uid_].adj_list.insert( Edge(this,new_uid,b,a) );

    // Make a new edge_element, and add to edgeList_ 
    Edge new_edge = Edge(this,new_uid,a,b); 
    edgeList_.push_back( edge_element(edge_value_type(), a, b, new_edge) );  
    edgeMap_.insert( new_uid ); 
    
    return new_edge;
  }

  /** Removes edge (a,b) from the graph, if it exists.
   *  @param[in] a   One node in edge to remove.
   *  @param[in] b   The other node in the edge to remove.
   *  @return        1 if edge is successfully removed, 0 otherwise. 
   *  
   *  @post  All edges e with nodes (@a a, @a b) are invalidated.
   *
   *  Complexity: O(k), k is the maximum degree of a node in graph.
   */
  size_type remove_edge(const Node& a, const Node& b) { 
    // Return 0 if either node is not in graph.
    if ( !has_node(a) || !has_node(b) ) { return 0; }
    
    // Check if edge is in graph: O(k), k = a.degree() 
    Edge e;  
    for (auto it = a.edge_begin(); it != a.edge_end(); ++it) { 
      if ((*it).node2() == b) { return remove_edge(*it); }
    } 
   
    return 0; // Edge wasn't found
  } 

  /** Removes edge e from the graph.
   *  @param[in] e   The edge to remove from the graph.
   *  @return        1 if edge is successfully removed, 0 otherwise. 
   * 
   *  @pre   @a e is a valid edge. 
   *  @post  All edges f, where f == @a e, are invalidated.
   *
   *  Complexity: Average O(1) 
   */
  size_type remove_edge(const Edge& e) {
    // Check if edge is in graph (std::unordered_map edgeMap_, i.e. hash table)
    auto in_graph = edgeMap_.find(e.uid_); // Average O(1) 
    if ( e.graph_ == NULL || in_graph == edgeMap_.end() ) { return 0; }

    // Update adjacency lists (std::undordered_map i.e. hash table) 
    nodeList_[e.node1().uid_].adj_list.erase(e); // Average O(1)
    nodeList_[e.node2().uid_].adj_list.erase(e); // Average O(1)

    // Remove from edgeMap_ (actually a hashtable) 
    edgeMap_.erase( e.uid_ ); // Average O(1)
 
    return 1; 
} 

  /** Removes edge e from the graph.
   *  @param[in] e_it  Iterator to edge to remove from the graph.
   *  @return          1 if edge is successfully removed, 0 otherwise. 
   * 
   *  @pre   @a e_it points to a valid edge. 
   *  @post  All edges f, where f == *(@a e_it), are invalidated.
   *
   *  Complexity: Average O(1) 
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge( *e_it ); 
    return ++e_it;
  } 

  /** Remove all nodes and edges from this graph.
   *  @post num_nodes() == 0 && num_edges() == 0
   *  @post All existing Node and Edge objects are invalidated
   *
   *  Complexity: O(V+E), V is the number of nodes, E the number of edges in graph. 
   */
  void clear() {
    nodeList_.clear();
    idx2uid_.clear(); 
    edgeList_.clear();
    edgeMap_.clear(); 
  }

  /** @class Graph::NodeIterator
   *  @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. 
     */
    NodeIterator() {} 

    /** Dereference this NodeIterator. 
     *  @return  The Node object that this is pointing to. 
     *  @pre     This node iterator is in range [begin, end) 
     */
    Node operator*() const {
      return graph_->nodeList_[graph_->idx2uid_[idx_]].n;
    }
 
    /** Increment this NodeIterator. 
     *  @pre     This node iterator is in range [begin, end) 
     *  @return  The incremented iterator           
     */ 
    node_iterator& operator++() { 
      idx_ += 1;
      return *this; 
    } 

    /** Check if this iterator is equal to @a iter2
     *  @pre     Both iterators are in range [begin, end) 
     */ 
    bool operator==(const node_iterator& iter2) const {
      return (graph_ == iter2.graph_ && idx_ == iter2.idx_);
    } 

   private:
    friend class Graph;       // Allows Graph to access private members. 
    const graph_type* graph_; // Pointer to the graph this node is part of
    size_type idx_;           // Index of this iterator in the Graph

    /** Private constructor for node iterator.
     *  @pre   0 <= iter_index < graph.size() 
     */
    NodeIterator(const graph_type* graph, size_type iter_index) : 
      graph_(graph), idx_(iter_index) { 
    } 
  
   };

  /** Return iterator pointing to front of the node list. 
   */
  node_iterator node_begin() const {
    return node_iterator(this, 0); 
  } 

  /** Return iterator pointing to one past the end of the node list.  
   */
  node_iterator node_end() const {
    return node_iterator(this, size());
  } 

  /** @class Graph::IncidentIterator
   *  @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
    using iter_type         = typename std::unordered_set<edge_type,edge_hash>::const_iterator;

    /** Construct an invalid IncidentIterator. 
     */
    IncidentIterator() {}

    /** Dereference this IncidentIterator. 
     *  @return  The Edge object that this is pointing to. 
     *  @pre     This IncidentIterator is in range [begin, end) 
     */
    Edge operator*() const {
      return (*it_); 
    }
 
    /** Increment this IncidentIterator. 
     *  @pre     This IncidentIterator is in range [begin, end) 
     *  @return  The incremented iterator           
     */ 
    IncidentIterator& operator++() { 
      ++it_;
      return *this; 
    } 

    /** Check if this iterator is equal to @a iter2
     *  @pre     Both iterators are in range [begin, end) 
     */ 
    bool operator==(const IncidentIterator& iter2) const {
      return (it_ == iter2.it_);
    } 

   private:
    friend class Graph;  // Graph can access private members/functions
    iter_type it_;       // Iterator to adjacency list representation

    /** Private constructor for IncidentIterator
     *  @pre it_ must be in range.
     */
    IncidentIterator(iter_type it) : it_(it) {}  
  
  };

  /** @class Graph::EdgeIterator
   *  @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
    using iter_type = std::unordered_set<size_type>::const_iterator;

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {}

    /** Dereference this EdgeIterator. 
     *  @return  The Edge object that this is pointing to. 
     *  @pre     This EdgeIterator is in range [begin, end) 
     */
    Edge operator*() const {
      return graph_->edgeList_[(*it_)].edge;
    }
 
    /** Increment this EdgeIterator. 
     *  @pre     This EdgeIterator is in range [begin, end) 
     *  @return  The incremented iterator           
     */ 
    EdgeIterator& operator++() { 
      ++it_;
      return *this; 
    } 

    /** Check if this iterator is equal to @a iter2
     *  @pre     Both iterators are in range [begin, end) 
     */ 
    bool operator==(const EdgeIterator& iter2) const {
      return (graph_ == iter2.graph_ && it_ == iter2.it_);
    } 

   private:
    friend class Graph;        // Gives Graph access to private functions/members
    const graph_type* graph_;  // Keeps a pointer to the graph
    iter_type it_;             // Iterator to edgeMap_

    /** Private constructor for Edge_Iterator
     *  @pre it_ must be in range. 
     */
    EdgeIterator(const graph_type* graph, iter_type it) : graph_(graph), it_(it) {}

  };
  
  /** Iterator pointing to the first element in the edge list.
   */
  edge_iterator edge_begin() const {
    return edge_iterator(this, edgeMap_.begin()); 
  } 

  /** Iterator pointing to one past the end of the edge list.
   */
  edge_iterator edge_end() const {
    return edge_iterator(this, edgeMap_.end());
  } 


 private:
  /** Hash function for edges, using their unique ID number */
  struct edge_hash{ 
    std::size_t operator()(const edge_type& e) const { 
      return std::hash<size_type>()(e.uid_);
    }
  }; 

  // Internal type for node objects
  struct node_element{ 
    size_type idx;                   // Node index (changes with deletes) 
    Point p;                         // Node position 
    node_value_type val;             // Node Value (of type V) 
    node_type n;                     // Node object
    std::unordered_set<edge_type,edge_hash> adj_list; // Adjacency list 
  };

  // Internal type for edge objects
  struct edge_element{
    edge_value_type val;           // Edge value (of type E) 
    node_type n1;                  // Node 1 in Edge
    node_type n2;                  // Node 2 in Edge
    edge_type edge;                // Edge object

    /** Helpful constructor for edge_element, helps make sure (a,b) == (b,a). 
     *  
     *  @post @a n1 < @a n2
     */
    edge_element(edge_value_type value, node_type node1, node_type node2, edge_type e) 
               : val(value), n1(node1), n2(node2), edge(e) {
      if (node2 < node1) { 
        n1 = node2; 
        n2 = node1;
      }; 
    } 
  }; 

  // Private Members
  std::vector<node_element>    nodeList_;  // Contains node elements
  std::vector<size_type>        idx2uid_;  // Maps an index (idx) to a node id (uid)
  std::vector<edge_element>    edgeList_;  // Contains edge elements
  std::unordered_set<size_type> edgeMap_;  // Keeps indices of existing edges.

};

#endif // CME212_GRAPH_HPP
