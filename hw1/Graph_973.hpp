#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <map>
#include <cassert>
#include <queue>
#include <set>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph {
 private:

	 struct internal_edges;
	 struct internal_nodes;
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
  Graph() : node_index(), edge_index(), graph_nodes(), graph_edges(), next_nid(0), next_eid(0)  {

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
    
    
    
    
    class Node: private totally_ordered<Node> {
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
        //Create an invalid node with nullptr as graph pointer
		graph_ = nullptr;
      // HW0: YOUR CODE HERE
    }

        
        //node_value_type& value();
    /** Return this node's position. */
    const Point& position() const {
        //Test if the node is valid and return the corresponding point
		assert(graph_->node_index.count(node_id)!=0);
		return graph_->graph_nodes[graph_->node_index.at(node_id)].point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
        //Test if the node is valid
    	assert(graph_!=nullptr);
    	assert(graph_->node_index.count(node_id)!=0);
      return graph_->node_index.at(node_id);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
      /** Find the node's value.
       * @return A reference to the value of this node.
       */
        node_value_type& value() {
            //test if the node is valid
            assert(graph_->node_index.count(node_id)!=0);
            return graph_->graph_nodes[graph_->node_index.at(node_id)].data;
        }
        /** Find the node's value.
         * @return A constant reference to the value of this node.
         */
        const node_value_type& value() const {
            //test if the node is valid
            assert(graph_->node_index.count(node_id)!=0);
            return graph_->graph_nodes[graph_->node_index.at(node_id)].data;
        }
        /** Find the node's degree.
         * @return The number of edges incident to this node.
         */
        size_type degree() const {
            return graph_->graph_nodes[this->index()].incident_edges.size();
        }
        /** Find the beginning from which the iterator goes through all the incident edges
         * @return The iterator pointing to the first incident edge to this node.
         *           this->edge_end() if this->degree()==0.
         */
        incident_iterator edge_begin() {
            return IncidentIterator(graph_,this, 0);
        }
        /** Find the end which the iterator reaches after it goes through all the incident edges
         * @return The iterator incremented by one after it reaches the last incident edge
         */
        incident_iterator edge_end() {
            return IncidentIterator(graph_,this, this->degree());
        }
        
        

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
       //Compare the unique id of two nodes
        if (this->graph_ != n.graph_) {
            return false;
        }
        
		return (this->node_id == n.node_id);
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
    //Compare the unique id of two nodes
      return ((this->node_id < n.node_id));
    }


        


        


   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Allow the node the access the graph
	graph_type* graph_;
    // the unique internal id for this node
	size_type node_id;
	// Node constructor that initializes the graph and its internal node id
	Node(const Graph* graph, size_type id) : graph_(const_cast<Graph*>(graph)), node_id(id) {}
    };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {

    return this->graph_nodes.size();
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
 /* Node add_node(const Point& position) {
      //update the map between node id and node index and the corresponding vector of points and node ids
      //Also update id for next node
	  node_index[this->next_nid] = graph_nodes.size();
	  next_nid++;
	  internal_nodes t;
	  t.nid=next_nid-1;
	  t.point=position;
	  graph_nodes.push_back(t);

    return Node(this, next_nid-1);        // Invalid node
  }*/

    
    Node add_node(const Point& position, const node_value_type& nval = node_value_type())  {
        node_index[this->next_nid] = graph_nodes.size();
        next_nid++;
        internal_nodes t;
        t.nid=next_nid-1;
        t.point=position;
        t.data=nval;
        graph_nodes.push_back(t);
        
        return Node(this, next_nid-1);
        
    }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    //Test if the id of the node is contained in the map node_index
    return (node_index.count(n.node_id)==1);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // Find the node id of the corresponding index and return the node
  	assert(i<graph_nodes.size());
    return Node(this, graph_nodes[i].nid);       
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
    /** Construct an invalid Edge. */
    Edge() {
        
        //create an invalid edge with graph pointer as the nullptr
      graph_=nullptr;
    }

    /** Return a node of this Edge */
    Node node1() const {
      //Use the index of the edge to access the first end point in graph
        if (orientation) {
            return graph_->graph_edges[index()].firstendpoint;
        }
        else {
            return graph_->graph_edges[index()].secondendpoint;
        }
    }

    /** Return the other node of this Edge */
    Node node2() const {
      //Use the index of the edge to access the second end point in graph
        if (orientation) {
            return graph_->graph_edges[index()].secondendpoint;
        }
        else {
            return graph_->graph_edges[index()].firstendpoint;
        }
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
        // Test if two edges have equal ids
      return (this->edge_id==e.edge_id);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
        // Compare the edge id of two edges
      return (this->edge_id<e.edge_id);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // The graph pointer that allows the edge to access the graph
    graph_type* graph_;
    // The unique id for this edge
	size_type  edge_id;
    bool orientation;
    // Returns the index of the edge in the graph
	size_type index() const   {
		assert(graph_!=nullptr);
		assert(graph_->edge_index.count(edge_id)==1);
		return graph_->edge_index.at(edge_id);
	}

      // Initialize the valid edge with the corresponding graph pointer and edge id
	Edge(const Graph* graph, size_type eid, bool os=true): graph_(const_cast<Graph*>(graph)), edge_id(eid), orientation(os) {}


  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    
    return graph_edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // Test if the index is valid
    assert(i<graph_edges.size());
    return Edge(this, graph_edges[i].eid);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Make sure that a and b are valid
    assert((a.graph_!=nullptr) && (b.graph_!=nullptr));
    assert(has_node(a) && has_node(b));
      
    // Find the edge connecting the two nodes if possible
    for (unsigned int i=0; i<graph_edges.size(); i++) {
    	node_type fpt=graph_edges[i].firstendpoint;
    	node_type spt=graph_edges[i].secondendpoint;
    	if ((a==fpt && b==spt) || (b==fpt && a==spt)) {
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
    // Update the map between edge id and edge index and the corresponding vector of node pairs and edge ids
    // Also update the next id for edge
      for (unsigned int i=0; i<graph_edges.size(); i++) {
          node_type fpt=graph_edges[i].firstendpoint;
          node_type spt=graph_edges[i].secondendpoint;
          if ((a==fpt && b==spt) || (b==fpt && a==spt)) {
              return Edge(this, graph_edges[i].eid);
          }
          
      }
    internal_edges t;
    t.eid=next_eid;
    t.firstendpoint=a;
    t.secondendpoint=b;
      
      for (unsigned int i=0; i<graph_nodes.size(); i++) {
          node_type fpt=node(i);
          if (a==fpt) {
              graph_nodes[i].incident_edges.push_back(std::pair<size_type,bool>(t.eid, true));
          }
          if (b==fpt) {
              graph_nodes[i].incident_edges.push_back(std::pair<size_type,bool>(t.eid, false));
          }
      }
      
    edge_index[next_eid]=graph_edges.size();
    graph_edges.push_back(t);
    next_eid++;
    return Edge(this, next_eid-1);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // Clear all the maps and vectors
  	node_index.clear();
  	edge_index.clear();
  	graph_nodes.clear();
  	graph_edges.clear();
    // HW0: YOUR CODE HERE
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
    class NodeIterator: private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
        it=nullptr;
        graph_=nullptr;
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
        
        
    /** Find the node the iterator is pointing to.
     * @return The node pointed to by the iterator.
     */
      Node operator*() const {
          
          //Check if the iterator is valid.
          assert(graph_ != nullptr);
          return Node(graph_, it->nid);
      }
      
    /** Return the iterator incremented by one. */
      NodeIterator& operator++() {
          it++;
          return *this;
      }
    // bool operator==(const NodeIterator&) const
    /** Evaluate whether or not two iterators are equal.
     * @return True if the private iterators it of the two are equal.
     *           False if the two are unequal.
     */
      bool operator==(const NodeIterator& other) const {
          return (it == other.it);
      }

   private:
    friend class Graph;
      // Pointer to its graph
      Graph* graph_;
      // Iterator that goes through its graph's vectors containing data about nodes.
      typename std::vector<internal_nodes>::const_iterator it;
      /** Constructor */
      NodeIterator(const Graph* graph, typename std::vector<internal_nodes>::iterator t): graph_(const_cast<Graph*>(graph)), it(t) {}
      /** Constructor */
      NodeIterator(const Graph* graph, typename std::vector<internal_nodes>::const_iterator t): graph_(const_cast<Graph*>(graph)), it(t) {}
      
      
      
    // HW1 #2: YOUR CODE HERE
  };
    
    

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
    
    /** Return the beginning iterator that goes through all the nodes*/
    node_iterator node_begin() const {
        return NodeIterator(this, graph_nodes.begin());
    }
  // node_iterator node_end() const
    
    /** Return the iterator incremented by one after it reaches the last node */
    node_iterator node_end() const {
        typename std::vector<internal_nodes>::const_iterator ss=graph_nodes.end();
        return NodeIterator(this, ss);
    }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
    class IncidentIterator: private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
        graph_=nullptr;
        t_node=nullptr;
        ind = 0;
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
        

    /** Return the edge which the iterator is pointing to */
      Edge operator*() const {
          //check if the iterator is valid
          assert(graph_ != nullptr);
          std::pair<size_type, bool> curedge=graph_->graph_nodes[t_node->index()].incident_edges[ind];
          return graph_->edge(curedge.first, curedge.second);
      }
    // IncidentIterator& operator++()
    /** Return the incremented iterator*/
      IncidentIterator& operator++() {
          ind++;
          return *this;
      }
    // bool operator==(const IncidentIterator&) const
    /** Check if two iterators are equal
     * @return True if the iterators have the same graph, nodes and corresponding incident edges.
     *           False otherwise.
     */
      bool operator==(const IncidentIterator& other) const {
          
          return ((graph_ == other.graph_) && (t_node == other.t_node) && (ind == other.ind));
      }

   private:
    friend class Graph;
      //The node to which the edges are incident
      Node* t_node;
      //The graph that contains the node
      Graph* graph_;
      //Index of the current edge
      size_type ind;
      /** Constructor */
      IncidentIterator(const Graph* graph, const Node* node, size_type index): graph_(const_cast<Graph*>(graph)), t_node(const_cast<Node*>(node)), ind(index) {}
      
    // HW1 #3: YOUR CODE HERE
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
    class EdgeIterator: private totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
        it=nullptr;
        graph_=nullptr;
    }

    /** Find the edge the iterator is pointing to.
     * @return The edge pointed to by the iterator.
     */

      Edge operator*() const {
          
          return Edge(graph_, it->eid);
      }
      
    /** Return the iterator incremented by one */
      EdgeIterator& operator++() {
          it++;
          return *this;
      }
    /** Evaluate whether or not two iterators are equal.
     * @return True if the private iterators it of the two are equal.
     *           False if the two are unequal.
     */
      
      bool operator==(const EdgeIterator& other) const {
          return (it == other.it);
      }
    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
      // Pointer to its graph
      Graph* graph_;
      // The iterator that goes through graph's vector containing data of the edges
      typename std::vector<internal_edges>::const_iterator it;
      /** Constructor */
      EdgeIterator(const Graph* graph, typename std::vector<internal_edges>::iterator t): graph_(const_cast<Graph*>(graph)), it(t) {}
      /** Constructor */
      EdgeIterator(const Graph* graph, typename std::vector<internal_edges>::const_iterator t): graph_(const_cast<Graph*>(graph)), it(t) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
    
  // edge_iterator edge_end() const
    /** Return the beiginning iterator that goes through all the edges of the graph */
    edge_iterator edge_begin() const {
        return EdgeIterator(this, graph_edges.begin());
    }
    // node_iterator node_end() const
    /** Return the iterator incremented by one after it reaches the last edge */
    edge_iterator edge_end() const {
        typename std::vector<internal_edges>::const_iterator ss=graph_edges.end();
        return EdgeIterator(this, ss);
    }

 private:
    // internal_nodes hold the information for nodes, including the corresponding node id and point
 	 struct internal_nodes {
 	 	size_type nid;
 	 	Point point;
        node_value_type data;
         std::vector<std::pair<size_type, bool>> incident_edges;
 	 };
    // internal_edges hold the information for edges, including the corresponding edge id and nodes
	 struct internal_edges {
	 	 size_type eid;
		 node_type firstendpoint;
		 node_type secondendpoint;
	 };
    
    /** Return the edge with the desired orientation and designated edge internal id */
    Edge edge(size_type id, bool ot) {
        return Edge(this, id, ot);
    }
    
    
	 // node_index records the correspondence between node internal id and node's index
	 // the same applies to edges
	 // the node's and edge's index are then used to access points and nodes in the corresponding vectors
	 std::map<size_type, size_type> node_index;
	 std::map<size_type, size_type> edge_index;
	 std::vector<internal_nodes> graph_nodes;
	 std::vector<internal_edges> graph_edges;
     // next_nid and next_eid determine the next node id and the next edge id respecrtively
	 size_type next_nid;
	 size_type next_eid;



};

#endif // CME212_GRAPH_HPP
