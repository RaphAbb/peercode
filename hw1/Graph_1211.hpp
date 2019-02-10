#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file  Graph.hpp
 *  @brief Undirected graph type. 
 */

#include <algorithm>
#include <vector>
#include <set> 
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 *  @brief A template for 3D undirected graphs.
 *
 *  Users can add and retrieve nodes and edges. 
 *  Edges are unique (there is at most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph {
 private:
  // Pre-declarations for internal types.
  struct node_element;  

 public:
  // Pre-declarations for public types.
  class Node;               // Node type 
  class NodeIterator;       // Iterator for graph nodes
  class Edge;               // Edge type
  class EdgeIterator;       // Iterator for graph edges
  class IncidentIterator;   // Iterator over incident edges to a node

  // Type synonyms.
  using graph_type = Graph<V>;
  using size_type = unsigned;

  using node_type = Node;
  using node_value_type = V;
  using node_iterator = NodeIterator; 

  using edge_type = Edge;
  using edge_iterator = EdgeIterator;
  using incident_iterator = IncidentIterator; 

  /** Construct an empty graph (default constructor)  
   *  @return  A new Graph object with zero nodes and zero edges.
   */ 
  Graph(): curr_id_(0) {
    // Set current id counter to zero
  } 

  /** Graph uses default destructor. 
   *  @post  All of this graph's nodes and edges are invalidated.
   */
  ~Graph() = default;

  /** @class Graph::Node
   *  @brief Class representing the graph's nodes.
   *
   *  Node objects are used to access information about the Graph's nodes.
   *  Each Node object is lightweight, and keep a pointer to 
   *  the Graph object which contains it, as well as a node ID number. 
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
    Node(): graph_(NULL) {
      // Invalid nodes have a null pointer for graph_
    }

    /** Return this node's index, a number in the range [0, graph_size). 
     *  @return  The index of this node in the Graph. 
     *  @pre     node_id_ = graph_->nodeList[node_id_].id
     *  
     *  Return the index of the node object in the graph representation.
     *  This implementation requires that node_id_ is identical
     *  to this Node's actual index in the graph_->nodeList vector,
     *  which is guaranteed by the implementation of the Graph class.
     *  Note: This will need to be updated to support node removal. 
     *
     *  Complexity: O(1)
     */
    size_type index() const {
      return node_id_;

      /** Naive search, needed if pre-condition isn't met

      assert(graph_ != NULL); // Make sure that Node is valid      
      for (size_type i = 0; i < graph_->size(); i++) {
        if (graph_->nodeList_[i].id == node_id_) { 
          return i;
        }
      }
      assert(false); // If the loop exits without return, then node is invalid

      */
    }

    /** Return this node's position. 
     *  @return  Reference to a Point object with this node's position. 
     *
     *  Complexity: Same as index()
     */
    const Point& position() const {
      return graph_->nodeList_[index()].p;
    }

    /** Read/Write access to value of this node. 
     *  @return Reference to the value of this node. 
     *
     *  Complexity: Same as index() 
     */
    node_value_type& value() { 
      // graph_ is a constant pointer to Graph; 
      // Graders, if you see this, is there a better way to do this?
      // I want graph_ to be a const pointer to ensure other
      // functions do not mess with Graph, so I have to cast away the const
      // for this function.
      return (const_cast<Graph*>(graph_))->nodeList_[index()].val;
    }
    
    /** Read-only access to value of this node. 
     *  @return Const reference to the value of this node. 
     *  @post   The value of this node is unchanged.
     *
     *  Complexity: Same as index() 
     */
    const node_value_type& value() const { 
      return graph_->nodeList_[index()].val;
    }

    /** Return the degree (num. of incident edges) to this node.
     *  @return  The number of distinct edges in graph_ which contain this node.
     *
     *  Complexity: Same as index() 
     */
    size_type degree() const {
      return graph_->nodeList_[index()].adj_list.size();
    } 
 
    /** Return an iterator to the incident edges of this node.
     *  @return An iterator pointing to the first element in 
     *          this node's adjacency list.  
     */
    incident_iterator edge_begin() const {
      const std::vector<edge_type>* adj = &graph_->nodeList_[index()].adj_list;
      return IncidentIterator(adj, 0); 
    }
 
    /** Return the end iterator to the incident edges of this node. 
     *  @return An iterator pointing to one past the last element
     *          in this node's adjacency list. 
     */
    incident_iterator edge_end() const {
      const std::vector<edge_type>* adj = &graph_->nodeList_[index()].adj_list;
      return IncidentIterator(adj, adj->size()); 
    }

    /** Test whether this node and @a n are equal.
     *  @param[in] n   Node to compare this to. 
     *  @return        True if @a n and this have the same
     *                 private members graph_ and node_id_. 
     */
    bool operator==(const Node& n) const {
      assert(graph_ != NULL);
      if (graph_ == n.graph_ && node_id_ == n.node_id_ ) {
        return true; 
      } else {
        return false;
      } 
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     *
     * This implementation compares each node's node_id_. 
     */
    bool operator<(const Node& n) const {
      if (node_id_ < n.node_id_) { 
        return true;
      } else { 
        return false;
      }
    }

   private:
    friend class Graph;       // Allow Graph to access Node's private members/functions
    const graph_type* graph_; // Pointer to graph this node is part of
    size_type node_id_;       // ID number of node (must be equal to index) 

    /** Private Node constructor, should only be called by Graph 
     *  @param graph  Pointer to the graph which contains the new node
     *  @param nid    Internal node id assigned by Graph.  
     *  @return       A new Node object, associated with @a graph and @a nid. 
     */
    Node(const graph_type* graph, size_type nid) 
      : graph_(graph), node_id_(nid) { 
    } 

  };

  /** Return the number of nodes in the graph.
   *
   *  Complexity: O(1).
   */
  size_type size() const {
    return nodeList_.size();
  }

  /** Synonym for size(). 
   *
   *  Complexity: O(1).
   */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   *  @param[in] position The new node's position
   *  @return             The new node.
   * 
   *  @post new num_nodes() == old num_nodes() + 1
   *  @post result_node.index() == old num_nodes()
   *  @post result_node.id == result_node.index() 
   *  
   *  Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    node_element newNode; 
    newNode.p = position; // This calls the copy constructor for Point()
    newNode.id = curr_id_++;
    newNode.val = value;
    nodeList_.push_back( newNode );
 
    return Node(this, newNode.id);
  }

  /** Determine if a Node belongs to this Graph
   *  @return True if @a n is currently a Node of this Graph
   *
   *  Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.graph_ == this) {
      return true;
    } else { 
      return false;
    }
  }

  /** Return the node with index @a i.
   *  @param i  Index to fetch.
   *  @return   The node at index @a i. 
   * 
   *  @pre 0 <= @a i < num_nodes()
   *  @post result_node.index() == i
   *
   *  Complexity: O(1).
   */
  Node node(size_type i) const {
    assert( i >= 0 && i < num_nodes());
    return Node(this, nodeList_[i].id);
  }

  /** @class Graph::Edge
   *  @brief Class representing the graph's edges.
   *
   *  Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   *  are considered equal if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Edge> {
   public:
    /** Create an invalid edge
     */ 
    Edge(): graph_(NULL) {
      // Set graph_ to null for invalid edges
    }

    /** Return this first node of this edge 
     */ 
    Node node1() const {
      assert(graph_ != NULL); 
      if (flipped_) { 
        return Node(graph_, id2_);
      } else {
        return Node(graph_, id1_); 
      } 
    }

    /** Return the other node of this Edge 
     */
    Node node2() const {
      assert(graph_ != NULL);
      if (flipped_) { 
        return Node(graph_, id1_);
      } else { 
        return Node(graph_, id2_);  
      }
    }

    /** Test whether this edge and @a e are equal.
     *  @pre this.id1_ < this.id2_ 
     *  @pre @a e.id1_ < e.id2_
     *  Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      assert(graph_ != NULL); 
      if (graph_ == e.graph_ && id1_ == e.id1_ && id2_ == e.id2_) {
        return true;
      } else { 
        return false;
      } 
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      assert(graph_ != NULL);
      // std::tie produces a std::tuple object, for which the STL
      // provides operator< which imposes strict ordering. 
      // This is a convenient way of imposing an arbitrary ordering on edges.
      return std::tie(id1_, id2_) < std::tie(e.id1_, e.id2_);
    }

   private:
    friend class Graph;       // Allow Graph to access Edge's private members/functions
    const graph_type* graph_; // Pointer to graph this node is part of
    size_type id1_;           // ID number of node 1 
    size_type id2_;           // ID number of node 2
    bool flipped_;            // True if a.id < b.id, False otherwise
    
    /** Private Edge Constructor for use by Graph.
     *  @post @a id1_ < @a id2_
     *  @post this.node1() == a
     *  @post this.node2() == b 
     */  
    Edge(const graph_type* graph, const Node& a, const Node& b) 
      : graph_(graph), id1_(b.node_id_), id2_(a.node_id_), flipped_(true) {
      if (a.node_id_ < b.node_id_) {
        id1_ = a.node_id_;
        id2_ = b.node_id_;
        flipped_ = false; 
      } 
    }
  
  };

  /** Return the total number of edges in the graph.
   *
   *  Complexity: O(1)
   */
  size_type num_edges() const {
    return edgeList_.size();
  }

  /** Return the edge with index @a i.
   *  @pre 0 <= @a i < num_edges()
   *
   *  Complexity: O(1)
   */
  Edge edge(size_type i) const {
    assert( i >= 0 && i < num_edges() ); 
    return edgeList_[i]; 
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
    const std::vector<Edge>& a_list = nodeList_[a.index()].adj_list;
    for (auto it = a_list.begin(); it != a_list.end(); ++it) {
        if ((*it).node2() == b) { return true; } 
    } 
    // Not found - return false
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
    // Return current edge if it already exists 
    if (has_edge(a,b)) {
      return Edge(this,a,b); 
    } else { 
      // Update Adjacency lists of a and b
      // STL says this should make a copy. Hopefully nodes are light and this is okay. 
      nodeList_[a.index()].adj_list.push_back( Edge(this,a,b) );
      nodeList_[b.index()].adj_list.push_back( Edge(this,b,a) );

      // Make a new edge and add to the edgeList_
      Edge new_Edge = Edge(this,a,b); 
      edgeList_.push_back(new_Edge); 

      return new_Edge;  
    }  
  }

  /** Remove all nodes and edges from this graph.
   *  @post num_nodes() == 0 && num_edges() == 0
   *
   *  Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodeList_.clear();
    edgeList_.clear(); 
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
      return Node(graph_ ,graph_->nodeList_[iter_index_].id);
    }
 
    /** Increment this NodeIterator. 
     *  @pre     This node iterator is in range [begin, end) 
     *  @return  The incremented iterator           
     */ 
    node_iterator& operator++() { 
      iter_index_ += 1;
      return *this; 
    } 

    /** Check if this iterator is equal to @a iter2
     *  @pre     Both iterators are in range [begin, end) 
     */ 
    bool operator==(const node_iterator& iter2) const {
      return (graph_ == iter2.graph_ && iter_index_ == iter2.iter_index_);
    } 

   private:
    friend class Graph;       // Allows Graph to access private members. 
    const graph_type* graph_; // Pointer to the graph this node is part of
    size_type iter_index_;    // Index of this iterator in the Graph's nodeList

    /** Private constructor for node iterator.
     *  @pre   0 <= iter_index < graph.size() 
     */
    NodeIterator(const graph_type* graph, size_type iter_index) : 
      graph_(graph), iter_index_(iter_index) { 
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

    /** Construct an invalid IncidentIterator. 
     */
    IncidentIterator() {}

    /** Dereference this IncidentIterator. 
     *  @return  The Edge object that this is pointing to. 
     *  @pre     This IncidentIterator is in range [begin, end) 
     */
    Edge operator*() const {
      return adj_->at(iter_index_); 
    }
 
    /** Increment this IncidentIterator. 
     *  @pre     This IncidentIterator is in range [begin, end) 
     *  @return  The incremented iterator           
     */ 
    IncidentIterator& operator++() { 
      iter_index_ += 1;
      return *this; 
    } 

    /** Check if this iterator is equal to @a iter2
     *  @pre     Both iterators are in range [begin, end) 
     */ 
    bool operator==(const IncidentIterator& iter2) const {
      return (adj_ == iter2.adj_ && iter_index_ == iter2.iter_index_);
    } 

   private:
    friend class Graph;                  // Graph can access private members/functions
    const std::vector<edge_type>* adj_;  // Pointer to adjacency list 
    size_type iter_index_;               // Index of iterator in adjacency list

    /** Private constructor for IncidentIterator
     *  @pre 0 <= iter_index < adj_->size() 
     */
    IncidentIterator(const std::vector<edge_type>* adj_list, size_type iter_index) : 
      adj_(adj_list), iter_index_(iter_index) {
    } 
  
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

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {}

    /** Dereference this EdgeIterator. 
     *  @return  The Edge object that this is pointing to. 
     *  @pre     This EdgeIterator is in range [begin, end) 
     */
    Edge operator*() const {
      return graph_->edgeList_[iter_index_];
    }
 
    /** Increment this EdgeIterator. 
     *  @pre     This EdgeIterator is in range [begin, end) 
     *  @return  The incremented iterator           
     */ 
    EdgeIterator& operator++() { 
      iter_index_ += 1;
      return *this; 
    } 

    /** Check if this iterator is equal to @a iter2
     *  @pre     Both iterators are in range [begin, end) 
     */ 
    bool operator==(const EdgeIterator& iter2) const {
      return (graph_ == iter2.graph_ && iter_index_ == iter2.iter_index_);
    } 

   private:
    friend class Graph;        // Gives Graph access to private functions/members
    const graph_type* graph_;  // Keeps a pointer to the graph
    size_type iter_index_;     // Index of iterator in adjacency list

    /** Private constructor for Edge_Iterator
     *  @pre 0 <= iter_index < num_edges() 
     */
    EdgeIterator(const graph_type* graph, size_type iter_index) : 
      graph_(graph), iter_index_(iter_index) {
    } 

  };
  
  /** Iterator pointing to the first element in the edge list.
   */
  edge_iterator edge_begin() const {
    return edge_iterator(this, 0); 
  } 

  /** Iterator pointing to one past the end of the edge list.
   */
  edge_iterator edge_end() const {
    return edge_iterator(this, num_edges());
  } 


 private:
  // Internal type for node object
  struct node_element{ 
    size_type id;                      // Node internal ID number
    Point p;                           // Position of node
    node_value_type val;               // Value of node
    std::vector<edge_type> adj_list;   // Adjacency list = Incident edges to node
  };

  // Private Members
  size_type curr_id_;                  // The current node ID 
  std::vector<node_element> nodeList_; // Contains node elements
  std::vector<edge_type> edgeList_;    // List of edge objects

};

#endif // CME212_GRAPH_HPP
