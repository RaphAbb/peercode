#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

# define print_debug 0

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
	
	struct internal_node;
  struct internal_edge;


 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node ;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
  
  /** Templated value types for node and edge */
  using node_value_type = V;
  using edge_value_type = E;

  /** Predeclaration of Edge type. */
  class Edge ;
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
  Graph() : nodes_(), size_(0), edges_(), num_edges_(0){
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
  class Node : private totally_ordered<Node>{
   public:
    /** Construct an invalid node.
     */
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      return this->graph_->nodes_[uid_]->P;
    }
    
    /** Return modifiable position. */
    Point& position() {
      return this->graph_->nodes_[uid_]->P;
    }
    
    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return this->graph_->nodes_[uid_]->uid;
    }
		
		//** Return value at node (modify-able) */
		node_value_type& value() {
			return this->graph_->nodes_[uid_]->value;
		}
		
	  /** Return value at node */
		const node_value_type& value() const {
			return this->graph_->nodes_[uid_]->value;
		}
		
		//** Change value at node... */
		void value(node_value_type v) {
			this->graph_->nodes_[uid_]->value = v;
		}
		
		/** Return how many edges eminate from this node */
		size_type degree() const {
			return this->graph_->adj_[index()].size();
		}
		
    /** Method to begin an iterator over connected edges */
		incident_iterator edge_begin() const {
      size_type begin = 0;
      return IncidentIterator(this->graph_, index(), begin);
		}
		
		incident_iterator edge_end() const {
      return IncidentIterator(this->graph_, index(), this->degree());
		}

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_) and (index() == n.index());
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
			if (index() < n.index()) return true;
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
		
		// Allow Graph to access pointer to Graph container
		Graph* graph_;
		
		// This node's uid
		size_type uid_;
		
		/** Private Constructor */
		Node(const Graph* set, size_type uid)
			 : graph_(const_cast<Graph*>(set)), uid_(uid){
		}
  };

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
	Node add_node(const Point& position) {
		// Overload the add_node function
		node_value_type val;
		return add_node(position, val);
	}
	
  Node add_node(const Point& position, const node_value_type& val) {

    internal_node* new_node = new internal_node;
		
		// add new node
		new_node->P = position;
		new_node->uid = size_;
		new_node->value = val;
    nodes_.push_back(new_node);
		++size_;
    
    // initialize node adjacency
    std::vector<size_type> temp_vec {};
    adj_.push_back(temp_vec);
		adj_edge_.push_back(temp_vec);
    
    return Node(this, size_-1);        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    (void) n;            // Quiet compiler warning
    return (n.index() < size() && this == n.graph_);
  }

  /** Return the node with index @a i. (get_node)
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
		if (i >= num_nodes()) {
			std::cout << "Inputted index: " << i << " and num_nodes: " << num_nodes() << std::endl;
			assert(i < num_nodes());
		}
    return Node(this, i);        // returns proxy node
  }
	
	
 /** Remove @a node from graph and erase connected edges.
  * @pre @a node has been added to the @a g graph.
	* @return 1 if successful.
	*
  * @post g.num_nodes() = old g.num_nodes() - 1
	* @post g.node(i).index() == i for all i with 0 ≤ i < g.num_nodes()
	* @post Edges connected to @a node are invalid.
	* @post g.num_edges() = old g.num_edges - @a node.degree()
	*
	* @note Invalidates all node iterators and edge iterators.
  *
  * Complexity: much less than O(num_nodes) if graph is sparse.
  */
	size_type remove_node(const Node& n) {
		/** idea is to remove node index i from all relevent structures, 
		* replace with last node, then pop off (now duplicated) last node */
		if (!has_node(n)) return 0;
		
		size_type ind = n.index();
		
		// remove edges associated with node
		std::vector<size_type> neighbors = adj_[ind];
		for (size_type k = 0; k < neighbors.size(); ++k){
			size_type neighbor = neighbors[k];
			remove_edge(n, node(neighbor));
			adj_[ind][k] = adj_[ind].back();
			adj_[ind].pop_back();
			
		}
		
		// update neighbors of last node, which now has index i
		internal_node* end_node = nodes_.back();
		size_type i_last = end_node->uid;
		for (size_type k = 0; k < adj_[i_last].size(); ++k){
			size_type neighbor = adj_[i_last][k];
			for (size_type j = 0; j < adj_[neighbor].size(); ++j) {
				if (adj_[neighbor][j] == i_last) {
					adj_[neighbor][j] = ind;
				}
			}
		}
		
		// remove node's adjacency list
		adj_[ind] = adj_[i_last];
		adj_.pop_back();
		adj_edge_[ind] = adj_edge_[i_last];
		adj_edge_.pop_back();
		
		// remove node
		nodes_[ind] = nodes_.back();
		nodes_.pop_back();
		--size_;
		
		// change uid of last node
		nodes_[ind]->uid = ind;
		
		return 1;
	}
	
  /** Remove node pointed at by @a n_it node iterator.
	 * @param[in] n_it  Node iterator
   * @pre @a n_it points at a valid node in the graph
   * @return iterator that points to next node
   *
   * @post g.num_nodes() = old g.num_nodes() - 1
   * @post g.node(i).index() == i for all i with 0 ≤ i < g.num_nodes()
   * @post Edges connected to @a node are invalid.
   * @post g.num_edges() = old g.num_edges - @a node.degree()
	 *
	 * @note Invalidates all node iterators and edge iterators.
   *
   * Complexity: much less than O(num_nodes) if graph is sparse.
   */
	node_iterator remove_node(node_iterator n_it) {
		assert(n_it != node_end());
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
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return n1;      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return n2;      // Invalid Node
    }
    
    /** Return the index this Edge */
		size_type index() const {
			return this->g_->edges_[uid_]->uid;
		}
		
    /** Return the length of this Edge */
    double length() const {
      return norm(n1.position() - n2.position());
    }
    
    /** Return value of edge (defined via template) - modifyable */
    edge_value_type& value() {
      return this->g_->edges_[uid_]->value;
    }
    
    /** Return value of edge (defined via template) */
    const edge_value_type& value() const {
      return this->g_->edges_[uid_]->value;
    }
    
    // edit which node is node1 and node2
    // this is unused currently but leaving here just in case
    void spawn_node(const Node& a) const {
      if (n1 != a) {
        Node b = g_->edges_[uid_]->n1;
        if (print_debug) {
          std::cout << "Spawn node uid: " << a.index() << "Adj node: " << n1.index() << std::endl;
        }
        n1 = a;
        n2 = b;
        if (print_debug) {
          std::cout << "Edge node 1: " << n1.index() << "Edge node 2: " << n2.index() << std::endl;
        }
      }
    }
    
    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((node1() == e.node1()) and (node2() == e.node2())) or ((node2() == e.node1()) and (node1() == e.node2()));
    }

    /** Test whether this edge is less than @a e in a global order.*/
    bool operator<(const Edge& e) const {
			if (index() < e.index()) return true;
			return false;
    }

   private:
    friend class Graph;
		
		Graph* g_;
		
		// This edge's uid
		size_type uid_;
		
		node_type n1;
		node_type n2;
		
    /** Private Constructor */
    Edge(const Graph* g, size_type uid, node_type n1_, node_type n2_)
        : g_(const_cast<Graph*>(g)), uid_(uid), n1(n1_), n2(n2_) {
    }
		
		
		
  };
  
  
  
  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type uid, size_type n1_, size_type n2_) const {
		assert(uid < num_edges());
		
    return Edge(this, uid, node(n1_), node(n2_));        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    
    for( size_type i = 0; i < adj_[a.index()].size(); ++i ) {
      if (adj_[a.index()][i] == b.index()) return true;
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
		
		// check if edge already exists
		if (has_edge(a, b)) {
			return Edge();
		}
		// check that nodes are in graph
		if (!has_node(a) or !has_node(b)) {
			return Edge();
		}

    internal_edge* new_edge = new internal_edge;
    
		
		// Set the rest of the internal_edge structure and add edge
		new_edge->n1 = a;
		new_edge->n2 = b;
		new_edge->uid = num_edges_; // may need to differentiate later from size attribute
		
    edges_.push_back(new_edge);
		++num_edges_;
		
    // fill out adjacency list
    // Node a is connected to Node b, and vice versa:
    if (print_debug) {
      std::cout << " Adjacency: Node a uid: " << a.index() << " and Node b uid: " << b.index() << std::endl; }
    
    adj_[a.index()].push_back(b.index());
    adj_[b.index()].push_back(a.index());
		
		adj_edge_[a.index()].push_back(num_edges_-1);
		adj_edge_[b.index()].push_back(num_edges_-1);
    
    return Edge(this, num_edges_-1, a, b);        // Invalid Edge
  }


  /** Remove an edge from a graph. Does not remove nodes.
   * @pre @a n1 and @a n2 are distinct valid nodes of this graph.
	 * @pre g.has_edge(@a n1, @a n2) returns true.
   * @return 1 if successful
	 *
   * @post has_edge(@a n1, @a n2) == false
   * @post g.num_edges() = old g.num_edges() -1
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
	size_type remove_edge(const Node& n1, const Node& n2) {
		//assert(has_edge(n1, n2));
		
		// find index of edge that connects node1 and node2
		size_type n1_id = n1.index();
		size_type n2_id = n2.index();
		
		for (size_type i = 0; i < adj_edge_[n1_id].size(); ++i) {
			size_type edge1 = adj_edge_[n1_id][i];
			
			for (size_type j = 0; j < adj_edge_[n2_id].size(); ++j) {
				size_type edge2 = adj_edge_[n2_id][j];
				
				if (edge1 == edge2) {
					//std::cout << "removing edge connected to nodes " << n1_id << " and " << n2_id << std::endl;
					//std::cout << edge1 << " edge id, but num_edges = " << num_edges_ << std::endl;
					return remove_edge(edge(edge1, n1_id, n2_id));
				}
			}
		}
		
		return 0; // if the assert passes, this shouldn't be needed
	}

  /** Remove an edge from a graph. Does not remove nodes.
   * @pre @a e is a valid edge of this graph.
   * @pre g.has_edge(@a e.n1, @a e.n2) returns true.
   * @return 1 if successful
   *
   * @post has_edge(@a e.n1, @a e.n2) == false
   * @post g.num_edges() = old g.num_edges() -1
	 *
	 * @note All iterators in the range [0, num_edges] are invalidated.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
	size_type remove_edge(const Edge& e) {
		assert(has_edge(e.node1(), e.node2()));
		
		size_type n1_id = e.node1().index();
		size_type n2_id = e.node2().index();
		size_type e_id = e.index();
		
		// remove edge uid from nodes' edge adjacency list
		for (size_type i = 0; i < adj_edge_[n1_id].size(); ++i) {
			if (e_id == adj_edge_[n1_id][i]) {
				adj_edge_[n1_id][i] = adj_edge_[n1_id].back();
				adj_edge_[n1_id].pop_back();
			}
		}
		for (size_type i = 0; i < adj_edge_[n2_id].size(); ++i) {
			if (e_id == adj_edge_[n2_id][i]) {
				adj_edge_[n2_id][i] = adj_edge_[n2_id].back();
				adj_edge_[n2_id].pop_back();
			}
		}
		
		// remove adjacency of nodes 
		for (size_type k = 0; k < adj_[n1_id].size(); ++k){
			if (adj_[n1_id][k] == n2_id) {
				adj_[n1_id][k] = adj_[n1_id].back();
				adj_[n1_id].pop_back();
			}
		}
		for (size_type k = 0; k < adj_[n2_id].size(); ++k){
			if (adj_[n2_id][k] == n1_id) {
				adj_[n2_id][k] = adj_[n2_id].back();
				adj_[n2_id].pop_back();
			}
		}
		
		// update edge id of last edge
		internal_edge* end_edge = edges_.back();
		size_type end_n1 = end_edge->n1.index();
		size_type end_n2 = end_edge->n2.index();
		for (size_type k = 0; k < adj_edge_[end_n1].size(); ++k){
			if (end_edge->uid == adj_edge_[end_n1][k]) adj_edge_[end_n1][k] = e_id;
		}
		for (size_type k = 0; k < adj_edge_[end_n2].size(); ++k){
			if (end_edge->uid == adj_edge_[end_n2][k]) adj_edge_[end_n2][k] = e_id;
		}
		
		edges_[e_id] = edges_.back();
		edges_.pop_back();
		edges_[e_id]->uid = e_id;
		
		--num_edges_;
		return 1;
	}
	
	
 /** Remove an edge pointed at by @a e_it edge iterator. Does not remove nodes.
  * @pre @a *e_it is a valid edge of this graph.
  * @pre g.has_edge(@a e.n1, @a e.n2) returns true.
	* @pre @a e_it is not the end iterator.
  * @return edge iterator to next edge.
  *
  * @post has_edge(@a e.n1, @a e.n2) == false
  * @post g.num_edges() = old g.num_edges() -1
  *
	* @note All iterators in the range [0, num_edges] are invalidated.
	*
  * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
  */
	edge_iterator remove_edge(edge_iterator e_it) {
		assert(e_it != edge_end());
		remove_edge(*e_it);
		return e_it;
	}
	
	
	
  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
		size_ = 0;
    edges_.clear();
		num_edges_ = 0;
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
		
		// p = uid for a Node. Used "p" since I was following the lecture notes
		size_type p; 
		
		// Dereference operator
		Node operator*() const {
			return graph_->node(p);
		}
		
		// Increment operator
		node_iterator& operator++() {
			++p;
			return *this;
		}
		
		// Equivalency operator
    bool operator==(const node_iterator& node_iter) const {
      return (p == node_iter.p);
    }

   private:
    friend class Graph;
		
		// allow pointer to graph (to use node() function)
		Graph* graph_;
		
		// private constructor
		NodeIterator(const Graph* g_, size_type p_)
			 : p(p_), graph_(const_cast<Graph*>(g_)){
		}
		
  };

  // Begin and end definitions for iterating over node list
  node_iterator node_begin() const {
		size_type p = 0;
  	return NodeIterator(this, p);
  }
	
	node_iterator node_end() const {
		size_type p = num_nodes();
		return NodeIterator(this, p);
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
		
		// spawning node uid:
    size_type node_uid;
		
		// counter for edges connected to spawn node:
		size_type edge_number;
		
		// Dereference operator
		Edge operator*() const {
			size_type edge_uid = graph_->adj_edge_[node_uid][edge_number];
			size_type node2_uid = graph_->adj_[node_uid][edge_number];
			
			return graph_->edge(edge_uid, node_uid, node2_uid);
		}
		
		// Increment operator
		incident_iterator& operator++() {
			++edge_number;
			return *this;
		}
		
		// Equivalency operator
    bool operator==(const incident_iterator& iit) const {
      return (node_uid == iit.node_uid) and (edge_number == iit.edge_number);
    }
		
   private:
    friend class Graph;
    
    // allow pointer to graph
    Graph* graph_;
    
    // private constructor
    IncidentIterator(const Graph* g_, size_type node_uid_, size_type edge_number_) :
    node_uid(node_uid_), edge_number(edge_number_), graph_(const_cast<Graph*>(g_)) {
    }
    
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
		
		size_type edge_uid;
		
    
		Edge operator*() const {
			size_type n1_uid = graph_->edges_[edge_uid]->n1.index();
			size_type n2_uid = graph_->edges_[edge_uid]->n2.index();
			
			return graph_->edge(edge_uid, n1_uid, n2_uid);
		}
		
		EdgeIterator& operator++(){
			++edge_uid;
			return *this;
		}
		
		bool operator==(const EdgeIterator& eit) const {
			return (edge_uid == eit.edge_uid);
		} 

   private:
    friend class Graph;
		
		Graph* graph_;
		
    // private constructor
    EdgeIterator(const Graph* g_, size_type edge_uid_) :
    edge_uid(edge_uid_), graph_(const_cast<Graph*>(g_)) {
    }
  };

	
	edge_iterator edge_begin() const {
		return EdgeIterator(this, 0);
	}
	
	edge_iterator edge_end() const {
		return EdgeIterator(this, num_edges());
	}

 private:
	 
	 // internal type for set nodes
	 struct internal_node {
	 	 Point P;
	 	 size_type uid; 	// index of node
		 node_value_type value; // can be temp, distance, weight, etc
	 };
	 
   std::vector<internal_node*> nodes_;
	 size_type size_;
  
   // internal type for set edges
	 struct internal_edge {
		 Node n1;
		 Node n2;
		 size_type uid;
     edge_value_type value;
	 };
	 
   std::vector<internal_edge*> edges_;
	 size_type num_edges_;
  
   // vector of (vector of node uid's), for adjacency
   std::vector<std::vector<size_type>> adj_;
   // adj_[1] should return an array of node uids adjacent to node 1
   // adjacent nodes are connected by an edge
	 
	 // vector of connected edges per node
	 std::vector<std::vector<size_type>> adj_edge_;
	 // adj_edge_[1] returns list of edge uids connected to node 1
  
};

#endif // CME212_GRAPH_HPP
