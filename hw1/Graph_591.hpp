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
template <typename V>
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
  
  using node_value_type = V;

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

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return this->uid_;
    }
		
		//** Return value at node */
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
		
		incident_iterator edge_begin() const {
      size_type begin = 0;
      return IncidentIterator(this->graph_, index(), begin);
		}
		
		incident_iterator edge_end() const {
      return IncidentIterator(this->graph_, index(), this->degree());
		}
		
    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value(); - DONE
    // const node_value_type& value() const; - DONE
    // size_type degree() const; - DONE
    // incident_iterator edge_begin() const; - DONE
    // incident_iterator edge_end() const; - DONE

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
		node_value_type val = 1;
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
    return (this == n.graph_);
  }

  /** Return the node with index @a i. (get_node)
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
		assert(i < num_nodes());
    return Node(this, i);        // returns proxy node
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

		size_type index() const {
			return uid_;
		}
		
    // edit which node is node1 and node2
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
    
    for( int i = 0; i < adj_[a.index()].size(); i++ ) {
      if (adj_[a.index()][i] == b.index()) return true;
    }
    //for( int i = 0; i < adj_[b.index()].size(); i++ ) {
    //  if (adj_[b.index()][i] == a.index()) return true;
    //}
    
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
		// Would like to use a template struct with uid and nxt, in the future
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
			 : graph_(const_cast<Graph*>(g_)), p(p_){
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
			
			// std::cout << "Edge uid, node1, node2 uid:  " << edge_uid << ", " << node_uid << ", " << node2_uid << std::endl;
			
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
    graph_(const_cast<Graph*>(g_)),node_uid(node_uid_), edge_number(edge_number_){
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
		
    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const
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
    // HW1 #5: YOUR CODE HERE
		
		Graph* graph_;
		
    // private constructor
    EdgeIterator(const Graph* g_, size_type edge_uid_) :
    graph_(const_cast<Graph*>(g_)), edge_uid(edge_uid_) {
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
	
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
		 size_type next_uid; // next node
	 };
	 
   std::vector<internal_node*> nodes_;
	 size_type size_;
  
   // internal type for set edges
	 struct internal_edge {
		 Node n1;
		 Node n2;
		 size_type uid;
		 
	 };
	 
   std::vector<internal_edge*> edges_;
	 size_type num_edges_;
  
   // vector of (vector of node uid's), for adjacency
   std::vector<std::vector<size_type>> adj_;
   // adj_[1] should return an array of nodes adjacent to node 1
   // adjacent nodes are connected by an edge
	 
	 // vector of connected edges per node
	 std::vector<std::vector<size_type>> adj_edge_;
	 // adj_edge_[1] returns list of edge uids connected to node 1
	 
	 
	 // 
	 // Iterator Helper Functions (had trouble implementing)
	 //
	 
	 // template list of elements
	 //template <typename T>
	 //struct list_elem {
	 //	 T uid;
	 //	 //list_elem<T>* nxt;
	 //};
	 
	 //// template list for 
	 //template <typename T>
	 //class list {
	 //	 list_iterator<T> 
	 //};
};

#endif // CME212_GRAPH_HPP
