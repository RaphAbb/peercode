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

using namespace std;

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
 private:
	vector<Point> point_list_;//this is a vector of point (coordinates) from the input .nodes file
	vector<vector<unsigned>> edge_list_;//vector of vectors, internal vectors are of size 2. 1st dimension is the index of the edge. Second dimension gives the indices for the nodes.
	vector<vector<unsigned>> neighbor_list_;//vector of vectors, this is an adjacency list. If nodes a and b share an edge, then in "row" a, "b" will appear
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
  Graph() {
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
    Node() {//nothing is needed here
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE

	return graph_->point_list_[this->index_];//return's the position using the index_ private variable

//      return Point();
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
	return this->index_;//this is defined later on in Node's private scope


//      return size_type(-1);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE

	if ((index_ == n.index()) and (graph_ == n.graph_)) {//checking the graph_ object is one way of ensuring invalid nodes can't be input
		  return true;
	  }else{
		  return false;
	  }

//      (void) n;          // Quiet compiler warning
//      return false;
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
	if (graph_ == n.graph_)
	{
		if (index_ < n.index()) {//check the node indices
			  return true;
	 	 }else{
			  return false;
	 	 }
	}
	else
	{
		printf("Invalid node, not in correct Graph object");
		return false;
	}

      //(void) n;           // Quiet compiler warning
      //return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
	size_type index_;
	Graph* graph_;
	
	//Default Node constructor
	Node(const Graph* gra, size_type index)
		: graph_(const_cast<Graph*>(gra)), index_(index){
	}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
   	 return this->point_list_.size();
	//return 0;
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
    // HW0: YOUR CODE HERE
	
//	Node n = Node(this,position)

//	if (has_node(n) == true)
//	{
//		return Node();
//	}else{ 

	vector<unsigned> nothing; 
	this->point_list_.push_back(position);
	this->neighbor_list_.push_back(nothing);//push back an empty row onto the neighbor_list_ in order to ensure it has the same number of size in the first dimension as the number of nodes - this avoids a segmentation fault
	return Node(this, num_nodes()-1);  
//	}

//(void) position;      // Quiet compiler warning
//    return Node();        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE

	if (this == n.graph_){
		return true;
	}else{
		return false;
	}

   // (void) n;            // Quiet compiler warning
   // return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE

	if (i >= num_nodes())
	{
		printf("Node index out of bounds.");
		return Node();
	}

	return Node(this,i);//just like get_element in proxy_example

//    (void) i;             // Quiet compiler warning
//    return Node();        // Invalid node
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

	size_type node_index = graph_->edge_list_[this->index_][0];//node_index
	return graph_->node(node_index);

 //	 return Node();      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE

	size_type node_index = graph_->edge_list_[this->index_][1];
	  return graph_->node(node_index);

 //     return Node();      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
	//compare the indices of 2 nodes; I assume the graph_ object will the be the same
	size_type index_node_0 = this->graph_->edge_list_[index_][0];
	size_type index_node_1 = this->graph_->edge_list_[index_][1];

	size_type index_e_0 = this->graph_->edge_list_[e.index_][0];
	size_type index_e_1 = this->graph_->edge_list_[e.index_][1];

	
	if (e.graph_ == graph_)
	{
      		if ((index_node_0 == index_e_0) and (index_node_1 == index_e_1)) 
		{	
		 	 return true;
		 }
		else if ((index_node_0 == index_e_1) and (index_node_1 == index_e_0))
		{
			return true;
		}
		else
		{
		  return false;
		 }
	}
	else
	{	printf("Invalid edge, not in this Graph.");
		return false;
	}

//      (void) e;           // Quiet compiler warning
//      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
	if (e.graph_ == graph_)
	{
		if (index_ < e.index_) {
			  return true;
	 	 }else{
			  return false;
	 	 }
	}else{
		printf("Invalid edge, not in this Graph");
		return false;
	}
//      (void) e;           // Quiet compiler warning
//      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
	
	size_type index_;
	Graph* graph_;
	
	//make a constructor for Edge
	
	Edge(const Graph* gra, size_type index)	//same implementation as Node constructor
		: graph_(const_cast<Graph*>(gra)), index_(index){
	}

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE

   return edge_list_.size();

    //return 0;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE

   if (i>= num_edges())//check if index out of bounds
	{
		printf("Edge index out of bounds");
		return Edge();
	}

   return Edge(this,i);//beacuse this is a constructor of a proxy class, this doesn't create a new instance of the object

 //   (void) i;             // Quiet compiler warning
 //   return Edge();        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE

	
	if ((a.graph_ != this) or (b.graph_ != this))//check validity
	{
		printf("Invalid node input");
		return false;
	}	
	
	//Get indices
	size_type a_index = a.index();
	size_type b_index = b.index();
	
	//Check for out of bounds error
	if ((a_index >= num_nodes()) or (b_index >= num_nodes()))
	{
		printf("Invalid node input");
		return false;
	}

	for (auto i = this->neighbor_list_[a_index].begin(); i != this->neighbor_list_[a_index].end(); ++i)//or may use for loop to check if this edge exists
	{
		if(b_index == *i)
		{
			return true;		
		}
	}
    return false;

//    (void) a; (void) b;   // Quiet compiler warning
//  return false;
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
	
	if ((a.graph_ != this) or (b.graph_ != this))//check validity of a and b
		{
		printf("Invalid node input");
		return Edge();
		}
	
	//Get indices of a and b
	size_type a_index = a.index();
	size_type b_index = b.index();

	for (unsigned i = 0; i != this->edge_list_.size(); ++i)//chose to use a for loop to check if a and b already comprise an edge; return if true
	{
		size_type index_0 = edge_list_[i][0];
		size_type index_1 = edge_list_[i][1];

		if ((a_index == index_0) and (b_index == index_1))
		{
			return Edge(this,i);
		}
		else if ((a_index == index_1) and (b_index == index_0))
		{
			return Edge(this,i);
		}
	}
			
	vector<unsigned> inserted_new_edge;//create a new edge vector to end of edge_list_
	inserted_new_edge.push_back(a_index);
	inserted_new_edge.push_back(b_index);

	edge_list_.push_back(inserted_new_edge);

	neighbor_list_[a_index].push_back(b_index);
	neighbor_list_[b_index].push_back(a_index);

	return Edge(this,edge_list_.size()-1);
	

//   (void) a, (void) b;   // Quiet compiler warning
//   return Edge();        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
   // HW0: YOUR CODE HERE
	point_list_.clear();
	edge_list_.clear();
	neighbor_list_.clear();

  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
 

 //nothing additional needed here
};

#endif // CME212_GRAPH_HPP
