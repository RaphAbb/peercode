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

using namespace std;
/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */

template <typename V, typename E>
class Graph {
 private:
  public:
//Public Type Definitions

int node_counter_uid=0;//uid counter, used for add_nodes and remove_nodes
int node_counter_id=0;//id counter, used for add_nodes and remove_nodes

	vector<Point> point_list_;//vec of points (coordinates) from input files
	vector<vector<unsigned>> edge_list_;//1st dimension is the index of the edge. Second dimension gives the indices for the nodes.
	vector<vector<unsigned>> neighbor_list_;//adjacency list
	vector<V> node_value_list_;//vector of node_value_type
  vector<E> edge_value_list_;//vector of edge_value_type
	vector<vector<unsigned>> edge_node_ref_;//vec of vecs, 1st dimension is index of node a, 2nd dimension is index of edges connected to it

//HW2P5
 unordered_map<unsigned,unsigned> u2i_;//unordered map, key = uid, val = id
 unordered_map<unsigned,unsigned> i2u_;//unordered map, key = id, val = uid

 unordered_map<unsigned,unsigned> u2i_edge_;//unordered map, key = uid, val = id
 unordered_map<unsigned,unsigned> i2u_edge_;//unordered map, key = id, val = uid


//Types, Predeclarations
  /** Type of this graph. */
  using graph_type = Graph<V,E>;
  using node_value_type = V;
  using edge_value_type = E;

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

/* @brief Return this node's position
* 
* return position <- value at current index_ (uid_), using l_value call
*/
    Point& position(){
				//assert(this->valid());
				return graph_->point_list_[this->uid_];
		}
/** @brief Return this node's position.
*		@return (Point) position <- xyz coordinates at current index_(uid)
*
*		@pre valid Graph object and Points
*		@post result has x,y,z, spatial coordinates
*
*		Complexity: O(1)
*/
    const Point& position() const {
				//assert(this->valid());
				return graph_->point_list_[this->uid_];//get position
    }

/** @brief Return this node's index, a number in the range [0, graph_size). 
*		@return size_type index_ of current node (id)
*
*		@pre valid Graph object, Points, and Nodes
*		@post 0<= result < num_nodes()
*
*		Complexity: O(1)
*/
    size_type index() const {
				return graph_->u2i_.find(this->uid_)->second;//from uid to id
    }

/**@brief Return number of edges incident to a node 
*		@return size_type that specifies number of edges incident to node
*
*		@pre valid Node
*		@post result.size() = number of edges incident to edge, deleted or not
*
* Note: this current degree() function returns the total edges, past and current, that are attached to a node; Further post-processing is required downstream
*
*		Complexity: O(constant)
*/
    size_type degree_total() const{
				return graph_->neighbor_list_[this->uid_].size();
	}

/*@ brief: true degree() function; returns total nodes currently incident to node
*
*see degree_total() for specs; same except
*
* Complexity: O(num_nodes)
*/
		size_type degree() const{
				int counter = 0;
				for (auto i = 0; i < graph_->neighbor_list_[this->uid_].size(); ++i){
						if (Node(graph_,graph_->neighbor_list_[this->uid_][i]).valid())
								counter++;
				}
				return counter;
		}



/**@brief Return first initial starting point for incident_iterator object
*
*		@return starting point for incident_iterator class iterations
*
*		@pre valid incident_iterator object
*		@post first iterator point of incident_iterator object
*
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
	return incident_iterator(graph_,this,degree_total());
	}

/**@brief Test whether this node and @a n are equal.
*
*Note: Equal nodes have the same graph and the same index (uid)
*
*	@return True if @a n is currently a node of this graph
*/
    bool operator==(const Node& n) const {
		if((uid_==n.uid_) and (graph_ == n.graph_))
		  return true;
	  else
		  return false;
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
				//if (graph_->u2i_.find(this->uid_)->second < n.index()) {
				if (graph_ < n.graph_)
						return true;
				if((uid_<n.uid_))
						return true;
				else
						return false;
						
				return false;
		}

/*@brief returns value stored by current node, at current uid
* 
* return node_value_type object <- value at current index_ (uid)
*/
    node_value_type& value(){
		return graph_->node_value_list_[this->uid_];
		}
		
/*@brief same as above, but returns the value as a constant
*
*return node_value_type object (constant) <- value at current uid
*/
    const node_value_type& value() const{
		return graph_->node_value_list_[this->uid_];
		}

/*@brief returns the indices of edges that are conncted to this node
*/
		vector<size_type> get_adj_edges() const{
		return graph_->edge_node_ref_[this->uid_];
		}


   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
	  size_type uid_;//set at creation, this is a uid
	  Graph* graph_;

	
/*@brief Default node constructor
*
* @param[in] gra, valid graph object
* @param[in] index, 0 <= @a index < num_nodes_ever_created
*
* @return valid Node object
*/
	Node(const Graph* gra, size_type index)
		: uid_(index), graph_(const_cast<Graph*>(gra)) {
	}

/*

*/
  public:
/*@brief: checks validity of node
* 
* Checks whether id of node exists, given uid
* @return bool value for whether node is valid
*/
  bool valid() const{
		if (graph_->u2i_.find(uid_) == graph_->u2i_.end())
				return false;
		return true;
		}
/*get uid, used for debugging*/
  size_type get_uid() const{
		return uid_;
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
		 return u2i_.size();
  }

  /** Synonym for size() -> see size() for documentation */
  size_type num_nodes() const {
    return size();
  }

 /* print statement used for debugging, same specs as size()*/
  void print_num_nodes() const {
		std::cout<<"print i2u_"<<std::endl;
		for(auto i = 0;i<i2u_.size();i++){
				auto pair = i2u_.find(i);
				std::cout<< pair->first<<":"<< pair->second <<std::endl;
		}
		std::cout<<"\n num_nodes() "<<num_nodes()<<std::endl;
  }

 /*return node_counter_uid, used for debugging*/
 int get_node_counter(){
				return node_counter_uid;
		}

/**@brief Add a node to the graph, returning the added node.
* @param[in] position The new node's position
* @post new num_nodes() == old num_nodes() + 1
* @post result_node.index() == old num_nodes()
*
* Note: also add uid,id pair to unordered maps for conversion to and from id/uid
*
* Complexity: O(1) amortized operations.
*/
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
  //increment
  node_counter_uid++;
  node_counter_id++;

	vector<unsigned> nothing; 
	this->point_list_.push_back(position);
	this->neighbor_list_.push_back(nothing);// avoids a segmentation fault
	this->edge_node_ref_.push_back(nothing);// avoid a segmentation fault
	this->node_value_list_.push_back(val);//push back node value

  //HW2P6; add (uid,id) or (id, uid) pair to unordered_map
  u2i_.insert(make_pair(node_counter_uid-1,node_counter_id-1));
  i2u_.insert(make_pair(node_counter_id-1,node_counter_uid-1));

  return Node(this,node_counter_uid-1);//create new, unique uid
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
  Node node(size_type i) const {//i is an inputted id
		auto uid = i2u_.find(i)->second;//from id to uid
		return Node(this,uid);
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

/* @brief Return length of Edge between 2 nodes
*  @return double, length of Edge 
*
*  @pre: valid Edge object
*
*  Complexity: O(1)
*/
		double length() const {
		Point p1 = node1().position();
		Point p2 = node2().position();

		return sqrt(pow(p1[0]-p2[0],2) + 
								pow(p1[1]-p2[1],2) + 
								pow(p1[2]-p2[2],2));
		}

    /**@brief Return a node of this Edge, using the uid
		*
		* @pre: valid Edge object
		*
		* @return Node object for one of nodes of this edge
		*/
		Node node1() const {
		auto node_index = graph_->edge_list_[this->uid_][0];//get node_uid
		return Node(graph_,node_index);
    }

    /**@brief Return the other node of this Edge, see node1() for specs */
    Node node2() const {
	  auto node_index = graph_->edge_list_[this->uid_][1];//get node_uid
		return Node(graph_,node_index);
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
	size_type index_node_0 = this->graph_->edge_list_[this->uid_][0];
	size_type index_node_1 = this->graph_->edge_list_[this->uid_][1];

	size_type index_e_0 = this->graph_->edge_list_[e.uid_][0];
	size_type index_e_1 = this->graph_->edge_list_[e.uid_][1];
	
		if (e.graph_ == graph_)	{
      		if ((index_node_0 == index_e_0) and (index_node_1 == index_e_1)) {	
						return true;
				}else if ((index_node_0 == index_e_1) and (index_node_1 == index_e_0)){
						return true;
		}else
		  return false;
		 
	}
	else{	printf("Invalid edge, not in this Graph.");
		return false;
	}
}

/**@brief Test whether this edge is less than @a e in a global order.
*
* @param[in] e, Edge object to be tested
* Note: in the case of comparison between edges in 2 graphs, one graph's edges will be assumed to be all strictly greater in magnitude than the other
* 
* @return True if current edge index is less than that of tested edge
*
* This ordering function is useful for STL containers such as
* std::map<>. It need not have any interpretive meaning.
*/
    bool operator<(const Edge& e) const {
		if (graph_ < e.graph_)
				return true;
		if (uid_ < e.uid_) {
			  return true;
	 	 }else{
			  return false;
	 	 }
    }

/*@brief Return index of edge, after conversion from uid to id
*
* @pre: valid Edge object
* @post 0<= result < num_edges()
* @result size_type id of Edge
*/
  size_type edge_index() const{
				return graph_->u2i_edge_.find(this->uid_)->second;//get id from uid
		}
/*@brief Return uid of this edge
* @pre: valid Edge object
* @post: 0 <= result < num_edges()
* @return: size_type uid of Edge
*
* Note!: edge uids are not necessarily unique. please see add_edge and remove_edge for more details
*/
  size_type edge_uid() const{
				return uid_;
		}

/*@brief returns value stored by current edge, at current index
* 
* return edge_value_type object <- value at current index_
*/
    edge_value_type& value(){	
		return graph_->edge_value_list_[this->uid_];
		}
		
/*@brief same as above, but returns the value as a constant
*
*return edge_value_type object (constant) <- value at current index_
*/
    const edge_value_type& value() const{
		return graph_->edge_value_list_[this->uid_];
		}


   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
	size_type uid_;
	Graph* graph_;
			
		/*@brief Default edge constructor
		*
		* @param[in] gra, valid graph object 
		* @param[in] index, 0 <= @a index <= num_edges_ever_created
		* @post Edge constructor 
		*
		* @return valid Edge object
		*/
			
	Edge(const Graph* gra, size_type index)	//same impl. as Node ctr
		: uid_(index), graph_(const_cast<Graph*>(gra)) {
	}

  /*

*/
  public:
  /*@brief: return validity of edge, first by checking validity of nodes, then of edge
  * see Node class valid() function for specs
  */
  bool valid() const{
		if(!node1().valid() or !node2().valid())
				return false;
  	if(graph_->u2i_edge_.find(uid_)==graph_->u2i_edge_.end())
				return false;
		return true;
		}

  };

  /**@brief Return the total number of edges in the graph.
   *
   * @pre: valid Graph object
   * @post: size_Type number of edges in current Graph, accounting for deletions
   * Complexity: No more than O(num_edges());
   */
  size_type num_edges() const {
	 return u2i_edge_.size();
  }

  /**@brief Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   * @return Edge edge object
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {//inputs a id
		auto uid = i2u_edge_.find(i)->second;//convert id to uid
    return Edge(this,uid);
  }
  /** @brief Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   * (also takes into account validity, in return statement)
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
	for (auto uid : edge_node_ref_[a.uid_]){
				auto e = Edge(this,uid);
				if (e.node2()==b)
						return e.valid();
		}
  for (auto uid : edge_node_ref_[b.uid_]){
				auto e = Edge(this,uid);
				if (e.node2()==a)
						return e.valid();
		}	
		return false;
}

  /**@brief Add an edge to the graph, or return current edge if already exists.
   * @pre @a a, @a b are distinct nodes; validity is checked in the first lines
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Also pushes back relevant uid/id onto the mappings
   *
   * Note: by design of the coder, edge uid = edge id always, upon construction
   * This is different from node uid and node ids
   *
   * Complexity: No more than O(num_edges())
   */
	
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val = edge_value_type()) {
		assert(a.valid());
		assert(b.valid());

		if (a == b)
				return Edge();//if nodes are equal, edge cannot exist

		if (has_edge(a,b)){//check if edge exists
				for (auto uid : edge_node_ref_[a.uid_]){
						auto e = Edge(this,uid);
				if (e.node2()==b){
						return e;
						}
				}
				for (auto uid : edge_node_ref_[b.uid_]){
						auto e = Edge(this,uid);
				if (e.node2()==a){
						return e;
						}
				}
		}

  //If edge not found, implement procedure to add edges
  size_type edge_counter = num_edges();
	
	//create a new edge vector, push to end of edge_list_		
	vector<unsigned> inserted_new_edge;
	inserted_new_edge.push_back(a.uid_);
  inserted_new_edge.push_back(b.uid_);
	std::vector<unsigned> ends {a.uid_,b.uid_};
  if(edge_counter < edge_value_list_.size()){
				edge_list_[edge_counter] = ends;
				edge_value_list_[edge_counter] = val;
		}else{
				edge_list_.push_back(inserted_new_edge);
				this->edge_value_list_.push_back(val);//avoid segmentation fault
		}

		//add to edge/vec reference
	edge_node_ref_[a.uid_].push_back(num_edges());
  //edge_node_ref_[b_index].push_back(edge_counter_uid-1);

	//add entry to adjacency list
	neighbor_list_[a.uid_].push_back(b.uid_);
	neighbor_list_[b.uid_].push_back(a.uid_);//for symmetry
	
  //HW2P6: add to maps
	u2i_edge_.insert(make_pair(edge_counter,edge_counter));
  i2u_edge_.insert(make_pair(edge_counter,edge_counter));
  
   return Edge(this,edge_counter);//return edge object, using uid
		//return edge(u2i_edge_.find(edge_counter_uid-1)->second);
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
  node_value_list_.clear();
  edge_value_list_.clear();
  edge_node_ref_.clear();

  u2i_.clear();
  i2u_.clear();

  u2i_edge_.clear();
  i2u_edge_.clear();

 node_counter_id = 0;
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
		*Note: due to assert statement, validity of node not necessary as @pre
		*@return Node at current index
		*/
     Node operator*() const{
			  assert(graph_->node(u_ind).valid());	//make sure node is valid
				return graph_->node(u_ind);
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
					return (u_ind!=n.u_ind) or (graph_!=n.graph_);
					}

   private:
    friend class Graph;
    Graph* graph_;
    size_type u_ind;//this is an id

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
		* @return valid NodeIterator object, associated with first node
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
	
		/** @brief dereferencing operator, takes in IncidentIterator, returns Edge
		*
		*@return next adjacent Edge object
		*/

     Edge operator*() const{
				//get current node uid
				auto n_uid = (*node_).uid_;
				//get adjacent node uid
				auto adj_n_uid = gra_->neighbor_list_[n_uid][u_ind];
				
				//find uid of edge between these 2 nodes
				size_type edge_uid;
				for (auto uid : gra_->edge_node_ref_[n_uid]){
						auto e = Edge(gra_,uid);
						if (e.node2()==Node(gra_,adj_n_uid)){
								edge_uid = e.uid_;
						}
				}
				for (auto uid : gra_->edge_node_ref_[adj_n_uid]){
						auto e = Edge(gra_,uid);
						if (e.node2()==Node(gra_,n_uid)){
								edge_uid = e.uid_;
						}
				}
						//return associated Edge object
						return Edge(gra_,edge_uid);
				}


		 /*Move IncidentIterator to next valid index, validity checked downstream */
     IncidentIterator& operator++(){
			++u_ind;
			return *this;
		  }

		/*Check equality between input Incidentiterator iit, and current instance */
     bool operator==(const IncidentIterator& iit) const{
				return (iit.u_ind==u_ind) and (iit.node_==node_) and (iit.gra_==gra_);
				}

	/*Check inequality between input Incidentiterator iit, and current instance */
     bool operator!=(const IncidentIterator& iit) const{
				return (iit.u_ind!=u_ind) or (iit.node_!=node_) or (iit.gra_!=gra_);
	}
  
	   private:
	    friend class Graph;
    friend class Node;
    Graph* gra_;
    Node* node_;
    size_type u_ind;//index of Edge, this is a uid

		/*@brief IncidentIterator constructor, iterates over edges incident to Node
		* 
		* @param[in] gra, valid graph object
		* @param[in] nod, valid Node object
		* @param[in] index of Edge
		* 
		* @return valid IncidentIterator object
		*/
    IncidentIterator(const Graph* gra, const Node* nod, size_type ind) :
	 gra_(const_cast<Graph*>(gra)), node_(const_cast<Node*>(nod)), u_ind(ind){}
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
    size_type u_ind;//index of edge (uid)

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
//HW2P6, additional functions

/* @brief: remove Node object, rearrange id/uid maps in order to account for removal
* @param[in]: e, Node& object to be removed
* @return: 1 for successful removal of node, 0 for unsuccesful removal of node
*
* Note: removal of nodes requires removal of adjacent edges
* 
* n need not be a valid node; this program will check for that; thus, no @pres
*
* Note: because edge and node maps are constructed differently, their removals are implemented differently
*
* @post: new num_nodes() = old num_nodes() - 1
* @post: for all i with 0<= i<num_nodes(), 0<=edge(i).id < num_nodes()
* @post: n.valid() == false
* @post: n.degree() == 0
*
* Notes: No objects are invalidated via deletion; validity statements and maps are used in lieu of deleting objects
*
* Complexity: O(num_edges), due to edge removal algorithm
*/
size_type remove_node(const Node& n){
	  if(n.valid() == false)
				return 0;	
		
		//remove adjacent edges
		for (auto itr = n.edge_begin(); itr != n.edge_end(); ++itr)	{
				auto ed = *itr;
				if(ed.valid()){
						size_type id_remove_e = ed.edge_index();//this is an id value
						 auto uid_remove_e = i2u_edge_.find(id_remove_e)->second;
						u2i_edge_.erase(uid_remove_e);
						i2u_edge_.erase(id_remove_e);

						for( auto z = id_remove_e+1;z<num_edges()+1;z++){
								auto pair = i2u_edge_.find(z);
						id_remove_e = pair->first;
						uid_remove_e = pair->second;

						i2u_edge_.erase(id_remove_e);
						i2u_edge_.insert(make_pair(id_remove_e-1,uid_remove_e));
						u2i_edge_.erase(uid_remove_e);
						u2i_edge_.insert(make_pair(uid_remove_e,id_remove_e-1));
				
						}
				}
			}
		size_type id_remove = n.index();//this is an id value
		auto uid_remove = i2u_.find(n.index())->second;

		auto uid_keep = i2u_.find((num_nodes()-1))->second;

	  //remove from inverse iterator (key = id, value = uid)
		i2u_.find(id_remove)->second = i2u_.find(num_nodes()-1)->second;
		i2u_.erase(num_nodes()-1);//erase last pair
		//remove from forward iterator (key = uid, value = id)
		u2i_.find(uid_keep)->second = id_remove;//set val= id in last pair
		u2i_.erase(uid_remove);//erase pair that should be erased		
		node_counter_id--;
	
		return 1;
		
}

/*@brief: Node Iterator version input of remove_node; see remove_node(const Node&) for specifications
*/
node_iterator remove_node(node_iterator n_it){
		//return NodeIterator(this,i2u_.find(remove_node(*n_it))->second);
		return remove_node((*n_it));
}


/*@brief: 2-node-input (Nodes input are the end nodes of the edge) version of remove_edge; see remove_edge(const Edge&) for specifications
*/
size_type remove_edge(const Node& a, const Node& b){
		if(has_edge(a,b)==false){
				printf("No such edge exists");
				return 0;
		}
		size_type e_uid = add_edge(a,b).uid_;
		return remove_edge(Edge(this,e_uid));
		
}	

/* @brief: remove Edge object, rearrange id/uid maps in order to account for removal
* @param[in]: e, Edge& object to be removed
* @return: 1 for successful removal of edge, 0 for unsuccesful removal of edge
* 
* e need not be a valid edge; this program will check for that; thus, no @pres
*
* Note: because edge and node maps are constructed differently, their removals are implemented differently
*
* @post: new num_edges() = old num_edges() - 1
* @post: for all i with 0<= i<num_edges(), 0<=edge(i).id < num_edges()
* @post: e.valid() == false
*
* Notes: No objects are invalidated via deletion; validity statements and maps are used in lieu of deleting objects
*
* Complexity: O(num_edges), due to removal algorithm
*/

size_type remove_edge(const Edge& e){
		//assert(e.valid());
		if (!e.valid())//check validity
				return 0;

		size_type id_remove = e.edge_index();//this is an id value
	  auto uid_remove = i2u_edge_.find(id_remove)->second;

		//remove corresponding elements from each map
		u2i_edge_.erase(uid_remove);
		i2u_edge_.erase(id_remove);
		if(id_remove == num_edges())
				return 1;

		for( auto z = id_remove+1;z<num_edges()+1;z++){//removal algorithm
				auto pair = i2u_edge_.find(z);
				id_remove = pair->first;
				uid_remove = pair->second;

				//deincrement id in each map, for all instances whos id > id_remove
				
				//remove from forward edge map (key = id, val = uid)
				i2u_edge_.erase(id_remove);
				i2u_edge_.insert(make_pair(id_remove-1,uid_remove));

				//remove from reverse iterator (key = uid, val = id)
				u2i_edge_.erase(uid_remove);
				u2i_edge_.insert(make_pair(uid_remove,id_remove-1));
				}
		
		return 1;
		
}
/* @brief remove_edge function where input is edge_iterator, see remove_edge(Edge&) for specs
*/
edge_iterator remove_edge(edge_iterator e_it){
		return EdgeIterator(this,i2u_edge_.find(remove_edge(*e_it))->second);
}

 private:
 //nothing additional needed here (from HW0, Graph class internals)
};

#endif // CME212_GRAPH_HPP
