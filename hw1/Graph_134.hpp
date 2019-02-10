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

template <typename V>
class Graph {
 private:
//Public Type Definitions
 public:
	vector<Point> point_list_;//vec of points (coordinates) from input files
	vector<vector<unsigned>> edge_list_;//1st dimension is the index of the edge. Second dimension gives the indices for the nodes.
	vector<vector<unsigned>> neighbor_list_;//adjacency list
	vector<V> value_list_;//vector of node_value_type
	vector<vector<unsigned>> edge_node_ref_;//vec of vecs, 1st dimension is index of node a, 2nd dimension is index of edge
//Types, Predeclarations
  /** Type of this graph. */
  using graph_type = Graph<V>;
  using node_value_type = V;

  /** Predeclaration of Node type. */
  class Node; /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Predeclaration of Edge type. */
  class Edge; /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator; /** Synonym for NodeIterator */
  using node_iterator = NodeIterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator; /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator; /** Synonym for IncidentIterator */
  using incident_iterator = IncidentIterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
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
    }

/** @brief Return this node's position.
*		@return (Point) position <- xyz coordinates at current index_
*
*		@pre valid Graph object and Points
*		@post result has x,y,z, spatial coordinates
*
*		Complexity: O(1)
*/
    const Point& position() const {
	return graph_->point_list_[this->index_];//get position
    }

/** @brief Return this node's index, a number in the range [0, graph_size). 
*		@return size_type index_ of current node
*
*		@pre valid Graph object, Points, and Nodes
*		@post 0<= result < num_nodes()
*
*		Complexity: O(1)
*/
    size_type index() const {
	return this->index_;
    }

/**@brief Return number of edges incident to a node 
*		@return size_type that specifies number of edges incident to node
*
*		@pre valid Node
*		@post result.size() = number of edges incident to edge
*
*		Complexity: O(1) (calls size())
*/
    size_type degree() const{
	return neighbor_list_[this->index()].size();
	}

/**@brief Return first initial starting point for incident_iterator object
*
*		@return starting point for incident_iterator class iterations
*
*		@pre valid incident_iterator object
*		@post first iterator point of incident_iterator object
*/
    incident_iterator edge_begin() const{
	return incident_iterator(graph_,this,0);
	}

/**@brief Return end point for incident_iterator object
*
*		@return ending point for incident_iterator class iterations
*
*		@pre valid incident_iterator object
*		@post last iterator point of incident_iterator object
*/
     incident_iterator edge_end() const{
	int n = edge_node_ref_[this->index()].size()-1;
	return incident_iterator(graph_,this,n);
	}

/**@brief Test whether this node and @a n are equal.
*
*Note: Equal nodes have the same graph and the same index.
*
*	@return True if @a n is currently a node of this graph
*/
    bool operator==(const Node& n) const {
	if ((index_ == n.index()) and (graph_ == n.graph_)) {
		  return true;
	  }else{
		  return false;
	  }
    }

/**@brief Test whether this node is less than @a n in a global order.
*
* @return True if index of current node is less than @a n
*
* This ordering function is useful for STL containers such as
* std::map<>. It need not have any geometric meaning.
*
* The node ordering relation must obey trichotomy: For any two nodes 
* and y, exactly one of x == y, x < y, and y < x is true.
*/
    bool operator<(const Node& n) const {
		if (graph_ == n.graph_){
				if (index_ < n.index()) {
						return true;
				}else{
						return false;
						}
		}else{
				printf("Invalid node, not in correct Graph object");
				return false;
				}
		}

/*@brief returns value stored by current node, at current index
* 
* return node_value_type object <- value at current index_
*/
    node_value_type& value(){	
		return graph_->value_list_[this->index_];
		}
		
/*@brief same as above, but returns the value as a constant
*
*return node_value_type object (constant) <- value at current index_
*/
    const node_value_type& value() const{
		return graph_->value_list_[this->index_];
		}

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
	  size_type index_;
	  Graph* graph_;
	
/*@brief Default node constructor
*
* @param[in] gra, valid graph object
* @param[in] index, 0 <= @a index < num_nodes()
*
* @return valid Node object
*/
	Node(const Graph* gra, size_type index)
		: index_(index), graph_(const_cast<Graph*>(gra)) {
	}
  };
/**@brief Return the number of nodes in the graph.
*
* @pre valid point_list_
* @pre size of point_list_
*
* @result size_type_ size of point_list_
*
* Complexity: O(1).
*/
  size_type size() const {
   	 return this->point_list_.size();
  }

  /** Synonym for size() -> see size() for documentation */
  size_type num_nodes() const {
    return size();
  }

/**@brief Add a node to the graph, returning the added node.
* @param[in] position The new node's position
* @post new num_nodes() == old num_nodes() + 1
* @post result_node.index() == old num_nodes()
*
* Complexity: O(1) amortized operations.
*/
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
	vector<unsigned> nothing; 
	this->point_list_.push_back(position);
	this->neighbor_list_.push_back(nothing);// avoids a segmentation fault
	this->edge_node_ref_.push_back(nothing);// avoid a segmentation fault
	this->value_list_.push_back(val);
	return Node(this, num_nodes()-1);  
  }

/**@brief  Determine if a Node belongs to this Graph
* @param[in] n, The input node
* @return True if @a n is currently a Node of this Graph
*
* Complexity: O(1).
*/
  bool has_node(const Node& n) const {
	if (this == n.graph_){
		return true;
	}else{
		return false;
	}
  }

/**@brief Return the node with index @a i.
* @pre 0 <= @a i < num_nodes()
* @post result_node.index() == i
*
* Complexity: O(1).
*/
  Node node(size_type i) const {
	if (i >= num_nodes())
	{
		printf("Node index out of bounds.");
		return Node();
	}

	return Node(this,i);//just like get_element in proxy_example
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
  class Edge: private totally_ordered<Edge> {
   public:
    /**@brief Construct an invalid Edge. */
    Edge() {
    }

    /**@brief Return a node of this Edge*/
    Node node1() const {
	size_type node_index = graph_->edge_list_[this->index_][0];//node_index
	return graph_->node(node_index);
    }

    /**@brief Return the other node of this Edge */
    Node node2() const {
	size_type node_index = graph_->edge_list_[this->index_][1];
	  return graph_->node(node_index);
    }

/**@brief Test whether this edge and @a e are equal.
* @param[in] e, Edge object to be tested
*
* @return True if edge is equal to current edge
*
* Equal edges represent the same undirected edge between two nodes.
* Tests edge by comparing index of node1 and node2 w/ that of current edge
* 
* Complexity: O(1)
*/
    bool operator==(const Edge& e) const {
	size_type index_node_0 = this->graph_->edge_list_[index_][0];
	size_type index_node_1 = this->graph_->edge_list_[index_][1];

	size_type index_e_0 = this->graph_->edge_list_[e.index_][0];
	size_type index_e_1 = this->graph_->edge_list_[e.index_][1];
	
	if (e.graph_ == graph_)	{
      		if ((index_node_0 == index_e_0) 
			and (index_node_1 == index_e_1)) {	
		 	 return true;
		 }else if ((index_node_0 == index_e_1) 
			and (index_node_1 == index_e_0)){
			return true;
		}
		else{
		  return false;
		 }
	}
	else{	printf("Invalid edge, not in this Graph.");
		return false;
	}
    }

    /**@brief Test whether this edge is less than @a e in a global order.
     *
		 * @param[in] e, Edge object to be tested
		 *
		 * 
		 * @return True if current edge index is less than that of tested edge
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
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
	size_type index_;
	Graph* graph_;
	
		/*@brief Default edge constructor
		*
		* @param[in] gra, valid graph object 
		* @param[in] index, 0 <= @a index <= num_edges()
		* @post Edge constructor 
		*
		* @return valid Edge object
		*/
			
	Edge(const Graph* gra, size_type index)	//same impl. as Node const
		: index_(index), graph_(const_cast<Graph*>(gra)) {
	}

  };

  /**@brief Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
   return edge_list_.size();
  }

  /**@brief Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
   if (i>= num_edges())//check if index out of bounds
	{
		printf("Edge index out of bounds");
		return Edge();
	}
   return Edge(this,i);
  }

  /** @brief Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {

	
	if ((a.graph_ != this) or (b.graph_ != this))
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
  }

  /**@briefAdd an edge to the graph, or return current edge if already exists.
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

	//return current edge if edge already exists
	if(has_edge(a,b) or has_edge(b,a)){
		for (size_type i = 0; i < edge_node_ref_[a_index].size(); i++){
				if (edge(edge_node_ref_[a_index][i]).node2() == b)
						return edge(edge_node_ref_[a_index][i]);
		}
		for (size_type i = 0; i < edge_node_ref_[b_index].size(); i++){
				if(edge(edge_node_ref_[b_index][i]).node2()==a)
			  //return Edge(this,i);
						return edge(edge_node_ref_[b_index][i]);
		}
	}
	
	//create a new edge vector, push to end of edge_list_		
	vector<unsigned> inserted_new_edge;
	inserted_new_edge.push_back(a_index);
	inserted_new_edge.push_back(b_index);
	edge_list_.push_back(inserted_new_edge);
	
	//add to edge/vec reference
	edge_node_ref_[a_index].push_back(edge_list_.size()-1);
	
	//add entry to adjacency list
	neighbor_list_[a_index].push_back(b_index);
	neighbor_list_[b_index].push_back(a_index);//for symmetry

	return Edge(this,edge_list_.size()-1);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
	point_list_.clear();
	edge_list_.clear();
	neighbor_list_.clear();
  }

  //
  // Node Iterator
  //
  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator {
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
		
		/** @brief dereferencing operator, takes in NodeIterator, returns Node
		*
		*@return Node at current index
		*/
     Node operator*() const{
			return Node(graph_,u_ind);
		}
		/*@brief increment operator, goes to iterator pointing to next node
		 *
		 *@return NodeIterator& pointing to next node
		 *
		 */
     NodeIterator& operator++(){
				++u_ind;
				return *this;
				}
		 /*Check equality between input node and current NodeIterator's node */
		 bool operator==(const NodeIterator& n) const{
						return (u_ind==n.u_ind) and (graph_==n.graph_);
						}
		 /*Check inequality between input node and current NodeIterator's node*/
		 bool operator!=(const NodeIterator& n) const{
						return u_ind!=n.u_ind;
						}

   private:
    friend class Graph;
    Graph* graph_;
    size_type u_ind;

		/*@brief NodeIterator constructor
		* 
		* @param[in] gra, valid graph object
		* @param[in] u, index where 0<= @a u < num_nodes()
		*
		* @return valid NodeIterator object
		*/
    NodeIterator(const Graph* gra, size_type u) :
		 graph_(const_cast<Graph*>(gra)), u_ind(u) {
		}
  };

   /*@brief return iterator associated with first node of graph
		*
		* @return valid NodeIterator object, associated with last node
		*used as initial condition in loop
		*
		* Complexity: O(1)
		 */
   node_iterator node_begin() const{
		return NodeIterator(this,0);
			}
   /*@brief return iterator associated with last node of graph
		*
		* @return valid NodeIterator object, associated with last node
		* used as termination condition in loop
		*
		* Complexity: O(1)
		 */
   node_iterator node_end() const{
		return NodeIterator(this,num_nodes());
		}	

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator {
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
     Edge operator*() const{
				return gra_->edge(gra_->edge_node_ref_[(*node_).index()][u_ind]);
				}

		 /** Increment IncidentIterator to next valid instance */
     IncidentIterator& operator++(){
			++u_ind;
		  return *this;
		  }

		/*Check equality between input Incidentiterator iit, and current instance */
     bool operator==(const IncidentIterator& iit) const{
				return (iit.u_ind==u_ind) and (iit.node_==node_) and (iit.gra_==gra_);
				}

   private:
    friend class Graph;
    friend class Node;
    Graph* gra_;
    Node* node_;
    size_type u_ind;//index of Edge

		/*@brief IncidentIterator constructor, iterates over edges incident to Node
		* 
		* @param[in] gra, valid graph object
		* @param[in] nod, valid Node object
		* @param[in] index of Edge
		* 
		* @return valid IncidentIterator object
		*/
    IncidentIterator(const Graph* gra, const Node* nod, size_type ind) :
	 gra_(const_cast<Graph*>(gra)), node_(nod), u_ind(ind){}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator {
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
		/*dereference EdgeIterator to return the Edge at this instance */
     Edge operator*() const{
		return Edge(graph_,u_ind);
		}
		/*increment EdgeIterator to point to the next Edge*/
     EdgeIterator& operator++(){
		++u_ind;
		return *this;
		}

		/*check equality between input EdgeIterator object, e, and current */
     bool operator==(const EdgeIterator& e) const{
		return u_ind==e.u_ind and graph_==e.graph_;
		}

		/*check inequality between input EdgeIterator object, e, and current*/
		bool operator!=(const EdgeIterator& e) const{
		 return u_ind!=e.u_ind;
		}

   private:
    friend class Graph;
    Graph* graph_;
    size_type u_ind;

/*@brief EdgeIterator consturctor
*
* @param[in] gra, valid graph object
* @param u, index of Edge
* @post valid EdgeIterator constructor
* 
* @return valid Edge object
*/
    EdgeIterator(const Graph* gra, size_type u) :
	 graph_(const_cast<Graph*>(gra)), u_ind(u) {
		}

   
  };
/*@brief return iterator associated with first edge of graph
*
* @return valid EdgeIterator object, associated with first edge
* used as initial condition in loop
*
* Complexity: O(1)
*/
   edge_iterator edge_begin() const{
				return EdgeIterator(this,0);
		}
/*@brief return iterator associated with last edge of graph
*
* @return valid EdgeIterator object, associated with last edge
* used as termination condition
*
* Complexity: O(1)
*/
   edge_iterator edge_end() const{
				return EdgeIterator(this,num_edges());
		}
 private:
 //nothing additional needed here (from HW0, Graph class internals)
};

#endif // CME212_GRAPH_HPP
