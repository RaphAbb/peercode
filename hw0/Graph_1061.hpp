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


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
	
	struct internal_node;
	struct internal_edge;
	
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

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : nodes_(), size_(0), edges_(), num_edges_(0){
    // HW0: YOUR CODE HERE - DONE
  }

  /** Default destructor */
  ~Graph() {
  	delete[] nodes_;
  }

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node {
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
      // HW0: YOUR CODE HERE - DONE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE - DONE
      return fetch().P;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE - DONE
      return fetch().uid;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE - DONE
			if (fetch().uid == n.index()) return true;
      (void) n;          // Quiet compiler warning
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
      // HW0: YOUR CODE HERE - DONE
			if (fetch().uid < n.index()) return true;
			
      (void) n;           // Quiet compiler warning
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
		
		// Allow Graph to access pointer to Graph container
		Graph* set_;
		
		// This node's uid
		size_type uid_;
		
		/** Private Constructor */
		Node(const Graph* set, size_type uid)
			 : set_(const_cast<Graph*>(set)), uid_(uid){
		}
		
		// helper function to return the appropriate node
		// loops over nodes until it finds nodes with correct uid
    internal_node& fetch() const {
      for (size_type i = 0; i < set_->size(); ++i) 
        if (set_->nodes_[i].uid == uid_)
          return set_->nodes_[i];
      assert(false);
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE - DONE
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
    // HW0: YOUR CODE HERE - DONE
		
		internal_node* new_node = new internal_node[size_+1];
		
		// Copy current nodes to new map
		for (size_type i = 0; i < size_; ++i) {
			new_node[i] = nodes_[i];
		}
		
		// add new node
		new_node[size_].P = position;
		new_node[size_].uid = size_;
		
		// Delete old nodes and reassign its value
		delete[] nodes_;
		nodes_ = new_node;
		++size_;
		
    (void) position;      // Quiet compiler warning
    return Node(this, size_-1);        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    (void) n;            // Quiet compiler warning
    return false;
  }

  /** Return the node with index @a i. (get_node)
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE - DONE
		assert(i < num_nodes());
    (void) i;             // Quiet compiler warning
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return fetch().n1;      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return fetch().n2;      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE - DONE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
		
		Graph* g_;
		size_type uid_;
		
    /** Private Constructor */
    Edge(const Graph* g, size_type uid)
        : g_(const_cast<Graph*>(g)), uid_(uid) {
    }
		
    internal_edge& fetch() const {
      for (size_type i = 0; i < g_->num_edges(); ++i)
        if (g_->edges_[i].uid == uid_)
          return g_->edges_[i];
      assert(false);
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE - DONE
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE - DONE
		assert(i < num_edges());
    (void) i;             // Quiet compiler warning
    return Edge(this, i);        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE - DONE
		for (size_type i = 0; i < num_edges_; ++i) {
			if (edges_[i].n1 == a) 
				if (edges_[i].n2 == b) return true;
			
			if (edges_[i].n2 == a) 
				if (edges_[i].n1 == b) return true;
			
		}
    (void) a; (void) b;   // Quiet compiler warning
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
    // HW0: YOUR CODE HERE - DONE
		
		// check if edge already exists
		if (has_edge(a, b)) {
			std::cout << "Edge already exists" << std::endl;
			return Edge();
		}
		
		internal_edge* new_edges = new internal_edge[num_edges_ + 1];
		
		// Copy current edges to a new array
		for (size_type i = 0; i < num_edges_; ++i)
			new_edges[i] = edges_[i];
		
		// Set the rest of the internal_edge structure
		new_edges[num_edges_].n1 = a;
		new_edges[num_edges_].n2 = b;
		new_edges[num_edges_].uid = num_edges_;
		
		// Delete old edges and reassign value
		delete[] edges_;
		edges_ = new_edges;
		++num_edges_;
		
    (void) a, (void) b;   // Quiet compiler warning
    return Edge(this, num_edges_-1);        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE - DONE?
		delete[] nodes_;
		nodes_ = new internal_node[0];
		size_ = 0;
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
	 
	 // internal type for set nodes
	 struct internal_node {
	 	 Point P;
	 	 size_type uid; 	// index of node
	 };
	 
	 internal_node* nodes_;
	 size_type size_;
	 
	 struct internal_edge {
		 Node n1;
		 Node n2;
		 size_type uid;
	 };
	 
	 internal_edge* edges_;
	 size_type num_edges_;
};

#endif // CME212_GRAPH_HPP
