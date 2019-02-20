#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

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
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  using node_value_type = V;
  using edge_value_type = E;
  
  /** Type of this graph. */
  using graph_type = Graph<V, E>;
  //  using GraphType = graph_type;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
  //  using NodeType = node_type;

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

  
  /** custom info for internal node (not used this time) */
  struct internal_edge {
    size_type nid_a_;
    size_type nid_b_;
    
    internal_edge(size_type nid_a, size_type nid_b):
      nid_a_(nid_a), nid_b_(nid_b) {}
  };

  struct edge_info {
    size_type uid_;
    edge_value_type edge_value_;

    edge_info(size_type uid, const edge_value_type &edge_value = edge_value_type()):
      uid_(uid) {}
  };
    
  /** custom info for internal node */
  struct internal_node {
    Point p_;
    size_type index_;
    node_value_type node_value_;
    std::vector<edge_info> adj_;
    internal_node(const Point& pos, size_type index, const node_value_type &node_value = node_value_type()):
      p_(pos), index_(index), node_value_(node_value), adj_() {}
  };
  
  
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    num_of_edges = 0;
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
      assert(valid());
      return graph_->node_list[uid_].p_;
    }
    
    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      assert(valid());
      return graph_->node_list[uid_].p_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      assert(valid());
      return graph_->node_list[uid_].index_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
    
    /** Return the value of the node
     * @return node's value
     */
    node_value_type& value() {
      assert(valid());
      return graph_->node_list[uid_].node_value_;
    }

    /** 
     * @return a const reference of this node's value
     */
    const node_value_type& value() const {
      assert(valid());
      return graph_->node_list[uid_].node_value_;
    }

    /**
     * @return the degree of this node
     */
    size_type degree() const {
      return graph_->node_list[uid_].adj_.size();
    }

    /**
     * @return the beginning of incident iterator
     */
    incident_iterator edge_begin() const {
      return incident_iterator(graph_, uid_, 0);
    }

    /**
     * @return the end of incident iterator
     */
    incident_iterator edge_end() const {
      return incident_iterator(graph_, uid_, graph_->node_list[uid_].adj_.size());
    }
    
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (graph_ == n.graph_) && (uid_ == n.uid_);
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
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    Graph* graph_;
    size_type uid_;

    /** Check if the node is valid */
    bool valid() const {
      return 0 <= uid_ && uid_ < graph_->node_list.size() && 
      graph_->node_list[uid_].index_ < graph_->i2u_.size() &&
      graph_->i2u_[graph_->node_list[uid_].index_] == uid_;
    }
    
    
    Node(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)), uid_(uid) {}
    
  };

  
  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return i2u_.size();
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
  Node add_node(const Point& position, const node_value_type& node_value = node_value_type()) {
    // HW0: YOUR CODE HERE
    node_list.push_back(internal_node(position, num_nodes(), node_value));
    i2u_.push_back(node_list.size() - 1);
    
    return Node(this, i2u_[i2u_.size() - 1]);   
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    
    return (n.graph_ == this) && n.valid();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    
    return Node(this, i2u_[i]);
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
      return Node(graph_, nid_a_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, nid_b_);
    }

    //bool valid_edge const {
    //  return  Node(this, nid_a_).valid() && Node(this, nid_b_).valid();
    //}
		   
    /** Return the internal value of the edge */
    edge_value_type& value() {
      size_type left = std::min(nid_a_, nid_b_);
      size_type right = std::max(nid_a_, nid_b_);
      
      for (size_type i = 0; i < graph_->node_list[left].adj_.size(); ++i) {
	if (graph_->node_list[left].adj_[i].uid_ == right) {
	  return graph_->node_list[left].adj_[i].edge_value_;
	}
      }
    }

    /** Calculate the length of the edge */
    double length() const {
      return norm(graph_->node_list[nid_a_].p_ - graph_->node_list[nid_b_].p_);
    }	

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (this->graph_ != e.graph_) return false;
      return (nid_a_ == e.nid_a_ && nid_b_ == e.nid_b_) ||
	(nid_a_ == e.nid_b_ && nid_b_ == e.nid_a_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      size_type sum_this = nid_a_ + nid_b_;
      size_type sum_e = e.nid_a_ + e.nid_b_;
      return sum_this < sum_e;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge object

    Graph * graph_;
    size_type nid_a_;
    size_type nid_b_;

    Edge(const Graph* graph, size_type nid_a, size_type nid_b)
      : graph_(const_cast<Graph*> (graph)), nid_a_(nid_a), nid_b_(nid_b) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return num_of_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    for (auto e = this->edge_begin(); e != this->edge_end(); ++e) {
      if (i == 0) {
	return *e;
      }
      else {
	i --;
      }
    }
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE

    for (auto e : node_list[a.uid_].adj_) {
      if (e.uid_ == b.uid_) {
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
    // HW0: YOUR CODE HERE

    if (has_edge(a, b)) {
      return Edge(this, a.uid_, b.uid_);
    }
    else {
      node_list[a.uid_].adj_.push_back(edge_info(b.uid_));
      node_list[b.uid_].adj_.push_back(edge_info(a.uid_));

      size_type left = std::min(a.uid_, b.uid_);
      size_type right = std::max(a.uid_, b.uid_);
      edge_list.push_back(internal_edge(left, right));
			  
      num_of_edges ++;
    }
    
    return Edge(this, a.uid_, b.uid_);  
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    node_list.clear();
    edge_list.clear();
    i2u_.clear();
    num_of_edges = 0;
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
    
    NodeIterator(const Graph* graph, size_type index) {
      graph_ = const_cast<Graph*> (graph);
      index_ = index;
    }
    
    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const
    
    /** Dereference operator of NodeIterator
     * @return the Node object the iterator is referenced to
     */
    Node operator*() const {
      return Node(graph_, graph_->i2u_[index_]);
    }

    /** Incrementor of NodeIterator
     * @pre 0 <= index_ < number of nodes
     * @post index_ increment by 1
     * @return the current iterator
     */
    NodeIterator& operator++() {
      index_ ++;
      return *this;
    }

    /** Compare if the current NodeIterator equals to another iterator
     * @param other another NodeIterator
     * @return true if they are equal, false otherwise
     */
    bool operator==(const NodeIterator& other) const {
      return graph_ == other.graph_ && index_ == other.index_;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type index_;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /**
   * @return the beginning of the iterator
   */
  node_iterator node_begin() {
    return NodeIterator(this, size_type(0));
  }

  /**
   * @return the beginning of the iterator
   */
  node_iterator node_begin() const {
    return NodeIterator(this, size_type(0));
  }

  /**
   * @return the end of the iterator
   */
  node_iterator node_end() {
    return NodeIterator(this, size_type(i2u_.size()));
  }

  /**
   * @return the end of the iterator
   */
  node_iterator node_end() const {
    return NodeIterator(this, size_type(i2u_.size()));
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

    IncidentIterator(const Graph* graph, size_type node_index, size_type degree_count):
      graph_(const_cast<Graph*> (graph)), node_index_(node_index), degree_count_(degree_count) {}
  
    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** Return the dereferenced edge of the incident iterator
     * @return the reference edge 
     */
    Edge operator*() const {
      return Edge(graph_, node_index_, graph_->node_list[node_index_].adj_[degree_count_].uid_);
    }

    /**
     * Incrementor of incident iterator
     * @pre 0 <= degree_count < number of incident of current node
     * @post degree_count increment by 1
     * @return the current incident iterator
     */
    IncidentIterator& operator++() {
      degree_count_++;
      return *this;
    }

    /**
     * Compare two IncidentIterators
     * @param other another IncidentIterator
     * @return true if the two iterators are equal, false otherwise
     */
    bool operator==(const IncidentIterator& other) const {
      return graph_ == other.graph_ && node_index_ == other.node_index_ && degree_count_ == other.degree_count_;
    }
  
   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;

    /** uid of the node that incident iterator is iterating on */
    size_type node_index_;
  
    /** degree of the node */
    size_type degree_count_;
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

    /** Construct an EdgeIterator */
    EdgeIterator(const Graph* graph, size_type index) :
      graph_(const_cast<Graph*> (graph)), index_(index) {}

    
    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

  
    /** Dereference operator of the EdgeIterator
     * @return the referenced edge
     */  
    Edge operator*() const {
      size_type a = graph_->edge_list[index_].nid_a_;
      size_type b = graph_->edge_list[index_].nid_b_;
      return Edge(graph_, a, b);
    }

   /**
    * Incrementor of EdgeIterator
    * @return the current iterator
    */
    EdgeIterator& operator++() {
      index_ ++;
      return *this;
    }

   /** Compare two EdgeIterator
    * @param other another EdgeIterator to compare with
    * @return true if two iterators are equal
    */  
    bool operator==(const EdgeIterator& other) const {
      return graph_ == other.graph_ && index_ == other.index_;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type index_;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Return the beginning of the edge_iterator
   *
   * @return the beginning of the iterator
   */
   edge_iterator edge_begin() const {
     return EdgeIterator(this, size_type(0));
   }  



  /** Return the end of the edge_iterator
   * @pre num_of_edge_ is valid
   * @return the end of the iterator
   */
   edge_iterator edge_end() const {
     return EdgeIterator(this, num_of_edges);
   }


  /** Remove an edge between two given nodes
   *  @return Return 1 if the edge is a valid one and has been removed. Return 0 if the edge was invalid
   *  Complexity: O(v+e), v is the maxium degree of nodes in the graph, e is the number of edges
   */
   size_type remove_edge(const Node& n1, const Node& n2) {

  	if (!has_edge(n1, n2)) {
  		return 0;
  	}
 
      size_type left = std::min(n1.uid_, n2.uid_);
      size_type right = std::max(n1.uid_, n2.uid_);
      
      for (size_type i = 0; i < edge_list.size(); ++ i) {
	if (edge_list[i].nid_a_ == left && edge_list[i].nid_b_ == right) {
	  edge_list.erase(edge_list.begin() + i);
	}
      }
      num_of_edges --;
      
      for (size_type i = 0; i < node_list[n1.uid_].adj_.size(); ++i) {
	if (node_list[n1.uid_].adj_[i].uid_ == n2.uid_) {
	  node_list[n1.uid_].adj_.erase(node_list[n1.uid_].adj_.begin() + i);
	  break;
        }
      }
      for (size_type i = 0; i < node_list[n2.uid_].adj_.size(); ++i) {
	if (node_list[n2.uid_].adj_[i].uid_ == n1.uid_) {
	  node_list[n2.uid_].adj_.erase(node_list[n2.uid_].adj_.begin() + i);
	  break;
        }
      }
		  
      return 1;
  }

  /** Remove a given edge.
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  } 

  /** Remove the edge that the edge iterator points to. Same as the previous function.
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(Node(this, edge_list[e_it.index_].nid_a_), Node(this, edge_list[e_it.index_].nid_b_));
    ++ e_it;
    return e_it;
   }

  /** Remove a given node.
   * @param n given node
   */  
  size_type remove_node(const Node& n) {
    if (!has_node(n)) {
      return 0;
    }
    size_type index = n.index();
    // Adjust index
    for (size_type i = n.index() + 1; i < num_nodes(); ++i) {
    	node_list[i2u_[i]].index_ --;
    }
    // Remove edges that link to this node
    std::vector<edge_type> edges_to_remove;
    std::copy(n.edge_begin(), n.edge_end(), std::back_inserter(edges_to_remove));
    std::for_each(edges_to_remove.begin(), edges_to_remove.end(), [this] (edge_type &e) { remove_edge(e); });
    // Remove the nodes 
    i2u_.erase(i2u_.begin() + index);

    for (size_type l = 0; l < edge_list.size(); ++ l) {
      if (edge_list[l].nid_a_ == n.uid_ || edge_list[l].nid_b_ == n.uid_) {
	edge_list.erase(edge_list.begin() + l);
	num_of_edges --;
      }
    }
    return 1;
  } 


  /** Remove the given node and the edges that link to it
   *  @pre @a n_it is a valid interator and @a n_it != node_end().
   *  @return A node_iterator new n_it such that new n_it == ++n_it.
   *  Complexity: O(num_nodes()) given a sparse graph.
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return n_it;
  } 

  
 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

    
  std::vector<internal_node> node_list;
  std::vector<internal_edge> edge_list;
  size_type num_of_edges;
  std::vector<size_type> i2u_;

};

#endif // CME212_GRAPH_HPP
