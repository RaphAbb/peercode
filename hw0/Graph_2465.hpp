#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

static const unsigned INVALID_VALUE= unsigned(-1);


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {

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

 private:
  struct Node_Information {
    size_type value;

    Point point;

    size_type smaller_num;

    Node_Information (size_type value_, const Point& point_, size_type smaller_num_ =0){
      value= value_ ;
      point= point_;
      smaller_num= smaller_num_;
    }

  };

  class Edge_Information{
  public:
    Edge_Information(const size_type& neighbor_node_= unsigned(-1)){
      neighbor_node= neighbor_node_;
    }

    bool operator<(const Edge_Information& another_edge) const {
      return this-> neighbor_node< another_edge.neighbor_node;
    }

    bool operator== (const Edge_Information& another_edge) const {
      return this-> neighbor_node == another_edge.neighbor_node;
    }
  private:
    friend class Graph;

    size_type neighbor_node;

  };

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //
 public:
  /** Construct an empty graph. */
  Graph() {
  	number_edges= 0;
  	nodes_index[INVALID_VALUE]= unsigned(-1);
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
    Node() {
    	Graph_N = NULL;
    	value= INVALID_VALUE;
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return Graph_N->Nodes_Information[Graph_N->nodes_index.at(value)].point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return Graph_N->nodes_index.at(value);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
             // Quiet compiler warning
      return this->Graph_N== n.Graph_N && this->value == n.value;
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
      if (this->value<n.value){
      	return true;
      }         
      if (this->value>n.value){
      	return false;
      }

      return this->Graph_N< n.Graph_N;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    graph_type* Graph_N;

    size_type value;

    Node (const graph_type* Graph_N_, size_type value_){
    	Graph_N= const_cast<graph_type*>((Graph_N_)), value = value_;
    }
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return Nodes_Information.size();
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
    size_type value_new= rand();
    while (nodes_index.count(value_new)){
    	value_new= rand();
    }
    this->nodes_index[value_new]= this->size();
    this->Nodes_Information.push_back(Node_Information(value_new, position));
    this->Edges_Information.push_back(std::vector<Edge_Information>(0));
    return Node(this, value_new); 
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
          // Quiet compiler warning
    return n.Graph_N== this && this->nodes_index.count(n.value) && n.value!= INVALID_VALUE;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node(this, Nodes_Information[i].value);        // Invalid node
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
    	Graph_e= NULL;
    	value_1= INVALID_VALUE;
    	value_2= INVALID_VALUE;
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(Graph_e, value_1);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(Graph_e, value_2);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return this->Graph_e == e.Graph_e && this->value_2== e.value_2 && this->value_1== e.value_1;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (this->value_1< e.value_1){
        return true;
      }
      if (this->value_1>e.value_1){
        return false;
      }
      if (this->value_2 < e.value_2){
        return true;
      }
      if (this->value_2 > e.value_2){
        return false;
      }
      return this-> Graph_e < e. Graph_e;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    graph_type* Graph_e;

    size_type value_1, value_2;

    Edge(const graph_type* Graph_e_, size_type value_1_, size_type value_2_){
    	Graph_e= const_cast<graph_type*>((Graph_e_)), value_1= value_1_, value_2= value_2_;
    }




    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return number_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    size_type node_1= 0;             // Quiet compiler warning
    while (i>= Nodes_Information[node_1].smaller_num){
    	node_1++;
    	i-= Nodes_Information[node_1].smaller_num;
    }
    return Edge(this, Nodes_Information[node_1].value, Nodes_Information[Edges_Information[node_1][i].neighbor_node].value);        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    size_type a_index= a.index() , b_index= b.index();
    auto it= lower_bound(Edges_Information[a_index].begin(), Edges_Information[a_index].end(), Edge_Information(b_index));
    if (it== Edges_Information[a_index].end()){
    	return false;
    }
    else{
    	return true;
    }
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
	assert (!(a== b));   // Make sure that a and b are not the same node.

	size_type index_a= a.index(), index_b= b.index();
	if (!this->has_edge(a, b)){
		number_edges++;
    Nodes_Information[std::max(index_a, index_b)].smaller_num++;
		auto it_a = lower_bound(Edges_Information[index_a].begin(), Edges_Information[index_a].end(), Edge_Information(index_a));
    Edges_Information[index_a].insert(it_a, Edge_Information(index_a));
    auto it_b= lower_bound(Edges_Information[index_b].begin(), Edges_Information[index_b].end(), Edge_Information(index_b));
    Edges_Information[index_b].insert(it_b, Edge_Information(index_b));
	}
    return Edge(this, Nodes_Information[index_a].value, Nodes_Information[index_b].value);        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_index.clear();
    Edges_Information.clear();
    Nodes_Information.clear();
    number_edges= 0;
  }

 private:
 	size_type number_edges;

 	std::unordered_map<size_type, size_type> nodes_index;

 	std::vector<Node_Information> Nodes_Information;

 	std::vector<std::vector<Edge_Information>> Edges_Information;
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
