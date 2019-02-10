#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP
/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <limits.h>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


template <typename V>
/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
 private:

    /** Pre-declaration for internal type for node */
    struct Internal_Node;
    
    /** Pre-declaration for internal type for edge */
    struct Internal_Edge;
    
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

  /** synonym of value type V **/
  using node_value_type = V;
    
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph, initialize node and edge vectors and the sizes of these vectors*/
    Graph()
    : vec_internal_edges(), vec_internal_nodes(),  size_nodes(0),  size_edges(0)
    {      
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
     *
     * Invalid Node with nullpointer graph and maxiumm int for ID
     * We use max int for a invalid ID so that it does not correspond to any internal node
     */
    Node():graph_container(nullptr),Node_ID(INT_MAX) {
    }
      
    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return get_IN().coordinate;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return get_IN().IN_id;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Returns the value associated with this node*/
    node_value_type& value(){
        return get_IN().value;
    }
    /** Returns value associated with the node that cannot be changed */
    const node_value_type& value() const{
        return get_IN().value;
    }

    /** Returns the number of nodes connected to this node
     * 
     * Complexity O(1)
     * 
     */
    size_type degree() const {
      return get_IN().associated_edges.size();
    }

    /** Instantiate a incident_iterator that points to the beginning of the vector of edges associated with this node
     * @return IncidentIterator that points to the first element of associated_edges
     */
    incident_iterator edge_begin() const {
        return IncidentIterator(this,0);
    }

    /** Instantiate a incident_iterator that points to the end of the vector of edges associated with this node
     * @return IncidentIterator that points to the past the last element of associated_edges
     */
    incident_iterator edge_end() const {
        return IncidentIterator(this,degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
        return (this->index() == n.index() && graph_container==n.graph_container);
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
      // Ordering nodes by their index
      assert(graph_container->has_node(n));
      return (this->index() < n.index());

    }
      

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
      
      //Private attributes of Node
      graph_type* graph_container;   // Pointer back to graph container
      size_type Node_ID;             // The Node's unique identification number

      
      /** Private Constructor**/
      Node(const graph_type* this_graph, size_type uid)
      :graph_container(const_cast<graph_type*>(this_graph)),Node_ID(uid){
      }
      
      /*Helper method to return internal node that corresponds to Node
       *The corresponding internal node should be at index [Node ID]
       */
      Internal_Node& get_IN() const {
          if (Node_ID<graph_container->size())
              return graph_container->vec_internal_nodes[Node_ID];
          assert(false);
      }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
      // can get size of vector of nodes:
      // return vec_internal_nodes.size();
    return size_nodes;
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
  Node add_node(const Point& position, const node_value_type& node_value= node_value_type()) {
    // HW0: YOUR CODE HERE
    // Create a new internal node and set the position and ID of that internal node
    Internal_Node temp_node; 
    temp_node.coordinate = position;
    temp_node.IN_id = size_nodes;  // *or* temp_node.IN_id = num_nodes();
    temp_node.node = Node(this,size_nodes);
    temp_node.value = node_value;
      
    //Add that internal node to the end of the internal nodes vector
    vec_internal_nodes.push_back(temp_node);

    //Increment the size of the graph
    ++size_nodes;
  
    return temp_node.node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
   //See if this the node's container is this graph. If it does this graph has this node
      return this == n.graph_container;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    //Make sure the graph has this node at index i
    assert (i>=0 and i<size());
    return vec_internal_nodes[i].node;
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
    Edge() :graph_container(nullptr),Edge_ID(INT_MAX){
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return get_IE().N1;      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return get_IE().N2;      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
        // Find the nodes that correspond to this edge
        node_type node_1 = graph_container->vec_internal_edges[Edge_ID].N1;
        node_type node_2 = graph_container->vec_internal_edges[Edge_ID].N2;
        
        //compare them to the nodes of Edge e
        return (node_1 == e.node1() and node_2 = e.node2())
           or (node_2 == e.node1() and node_1 ==e.node2());;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // No interpretive meanting to the ordering of edge. Order edges by indices
        return (Edge_ID < e.Edge_ID);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

     //Private attributes of Edge
      graph_type* graph_container;   // Pointer back to graph container
      size_type Edge_ID;        // The Edge's unique identification number
      
      /** Private Constructor**/
      Edge(const graph_type* this_graph, size_type uid)
      :graph_container(const_cast<graph_type*>(this_graph)),Edge_ID(uid){
      }
      
      /**Helper method to return internal node that corresponds to Node
       *@return internal node at index [Node ID]
       *
       *Complexity O(1)
       */
      Internal_Edge& get_IE() const {
          if (Edge_ID<graph_container->num_edges())
              return graph_container->vec_internal_edges[Edge_ID];
          assert(false);
      }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return size_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    //Make sure the graph has this edge at index i
    assert (i>= 0 and i<num_edges());
    return vec_internal_edges[i].edge;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
      
    // Make sure Node a and Node b are valid nodes of the graph
      assert (this->has_node(a) and this->has_node(b));
      
    //Find the vector of edges associated with node A
      std::vector<edge_type> Edges = vec_internal_nodes[a.index()].associated_edges;
      // For all the edges in this list, see if node B exists
      for (const auto& edge : Edges){
          if (edge.node1() == b or edge.node2() ==b) return true;
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
        // Make sure Node a and Node b are valid nodes of the
        // No self reference
        assert (this->has_node(a) and this->has_node(b));
        assert(not(a==b));
      
      //If the graph does not have this edge, must create new Internal Edge and add this edge to the two containers that keep track of edges
      // Create a new internal edge and set the position and ID of that internal edge
      if (not has_edge(a,b)) {
        Internal_Edge temp_internal_edge;
        temp_internal_edge.N1 = a; 
        temp_internal_edge.N2 = b; 
        temp_internal_edge.IE_id = size_edges;
        temp_internal_edge.edge = Edge(this, size_edges);

        //Add that internal edge to the end of the internal edge vector that orders by index
        vec_internal_edges.push_back(temp_internal_edge);

        //Add this edge to the vectors of edges of each nodes to keep track of connections
          vec_internal_nodes[a.index()].associated_edges.push_back(temp_internal_edge.edge);
          vec_internal_nodes[b.index()].associated_edges.push_back(temp_internal_edge.edge);
          
        //Increment the size of the graph
        ++size_edges;
      
        return temp_internal_edge.edge;
      }
      
      /*Test if the graph already has this edge*/
      if (has_edge(a,b)){
          //Find the vector of edges associated with node a
          std::vector<edge_type> Edges = vec_internal_nodes[a.index()].associated_edges;
          // For all the edges in this list, see if node b exists
          for (const auto& edge : Edges){
              if (edge.node1() == b or edge.node2() ==b) return edge;
          }
      }
          
      //Return invalid edge if nodes do not exist
      return Edge();
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {

    size_edges = 0;
    size_nodes = 0;
    vec_internal_edges.clear();
    vec_internal_nodes.clear();
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
        }
        
        // HW1 #2: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:

        /** Dereferences the Node associated with this node iterator
         * @return A node associated with this node iterator
         */
        Node operator*() const {
            return graph_container->vec_internal_nodes[iterator_index].node;
        }
        
        /** Advance iterator
         * @return A Node iterator with that points to the next Node
         * 
         * @post iterator index increase by 1
         */
        NodeIterator& operator++() {
            iterator_index++;
            return *this;
        }
        
        /** Test whether iterators are equivalent
         *@return True if @a NI.graph_container is the same as this NodeIterator's
         *        and @a NI.iterator_index is the same as this Nodeiterators'
         */
        bool operator==(const NodeIterator& NI) const {
            return (graph_container == NI.graph_container and iterator_index == NI.iterator_index);
        }
        
    private:
        friend class Graph;
        // HW1 #2: YOUR CODE HERE
        graph_type* graph_container;   // Pointer back to graph container
        size_type iterator_index;      //Index into the vector of Internal Nodes
        
        //Private Node iterator constructor
        NodeIterator(const graph_type* this_graph, size_type index) 
        :graph_container(const_cast<graph_type*>(this_graph)),iterator_index(index){
        }
    };
    
    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Instantiate a node_iterator that points to the beginning of the vector of internal nodes
     * @return NodeIterator that points to the first element of vec_internal_nodes
     */
    node_iterator node_begin() const {
        return NodeIterator(this, 0);
    }
    
    /** Instantiate a node_iterator that points to the end of the vector of internal nodes
     * @return NodeIterator that points to the element past the end of vec_internal_nodes
     */
    node_iterator node_end() const {
        return NodeIterator(this, vec_internal_nodes.size());
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /**Dereference the edge associated with this IncidentIterator
     * @return A Edge object this_edge associated with this IncidentIterator 
     * 
     * @post         this_edge.node1() == spawning node
     * @post         this_edge.node2() == adjacent node
     * */
    Edge operator*() const {
        edge_type& this_edge = node_container->get_IN().associated_edges[iterator_index];
        node_type spawning_node = *node_container;
        
        //See if node1 of the edge is the spawning node. If it is return the edge. 
        if (this_edge.node1() == spawning_node) return this_edge;
        //If node1 is not the spawning node, change the orders of the nodes
        this_edge.get_IE().N2 = this_edge.get_IE().N1;
        this_edge.get_IE().N1 = spawning_node;
        return this_edge;
    }

    /** Advance Incident Iterator
     *@return A Incident iterator with that points to the next Edge associated with this node
     *
     *@post iterator index increase by 1
     */
    IncidentIterator& operator++() {
        iterator_index++;
        return *this;
        }

    /** Test if iterators are equivalent
    *@return True if @a II.node_container is the same as this IncidentIterator's
    *        and @a II.iterator_index is the same as this Incidentiterators'
    */
    bool operator==(const IncidentIterator& II) const {
      return (node_container == II.node_container and iterator_index == II.iterator_index);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE

    node_type* node_container;   // Pointer back to graph container
    size_type iterator_index;      //Index into the vector of associated nodes
    
    //Private Node iterator constructor
    IncidentIterator(const Node* this_node, size_type index) 
    :node_container(const_cast<node_type*>(this_node)),iterator_index(index){
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Dereferences the Edge associated with this edge iterator
     * @return A edge associated with this edge iterator
     */
    Edge operator*() const {
        return graph_container->vec_internal_edges[iterator_index].edge;
    }
    
    /** Advance iterator
     * @return A Edge iterator with that points to the next Edge
     * 
     * @post iterator_index increase by 1
     */
    EdgeIterator& operator++() {
        iterator_index++;
        return *this;
    }
    
    /** Test whether iterators are equivalent
     *@return True if @a EI.graph_container is the same as this EdgeIterator's
    *        and @a EI.iterator_index is the same as this Edgeiterators'
    */
    bool operator==(const EdgeIterator& EI) const {
        return (graph_container == EI.graph_container and iterator_index == EI.iterator_index);
    }
    
   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_container;   // Pointer back to graph container
    size_type iterator_index;      //Index into the vector of Internal Edges
        
    //Private Node iterator constructor
    EdgeIterator(const graph_type* this_graph, size_type index) 
    :graph_container(const_cast<graph_type*>(this_graph)),iterator_index(index){
    }

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

    /** Instantiate a edge_iterator that points to the beginning of the vector of internal edges
     * @return EdgeIterator that points to the first element of vec_internal_edges
     */
    edge_iterator edge_begin() const {
        return EdgeIterator(this, 0);
    }

    /** Instantiate a edge_iterator that points to the beginning of the vector of internal edges
     * @return EdgeIterator that points to the element past the end of vec_internal_edges
     */
    edge_iterator edge_end() const {
        return EdgeIterator(this, vec_internal_edges.size());
    }
 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
    
    //Internal Nodes of the graph. Contains the actual point data
    struct Internal_Node {
        Point coordinate;       // A Node's cartesian coordinate
        size_type IN_id;        // The unique identification for an Internal Node
        Node node;              // Proxy node that corresponds to this Internal Node
        std::vector<edge_type> associated_edges; // Vector of edges associated with this node
        node_value_type value; // The value associated with this Node
    };
    
    //Internal Edges of the graph. Contains the nodes of each edge and an edge ID
    struct Internal_Edge {
      Node N1;          //First Node
      Node N2;          //Second Node
      size_type IE_id;  //The unique identification for an Internal Edge
      Edge edge;        // Proxy Edge that corresponds to this Internal Edge
      
    };
  
    // Attributes of the graph class
    
    std::vector<Internal_Edge> vec_internal_edges;
    std::vector<Internal_Node> vec_internal_nodes;
    size_type size_nodes;
    size_type size_edges;
};


#endif // CME212_GRAPH_HPP
