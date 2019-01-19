#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <limits.h>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

    struct Internal_Node;
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

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
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
    Node():graph_container(nullptr),Node_ID(INT_MAX) {
      //Invalid Node with nullpointer graph and maxiumm int for ID
      // We use max int for a invalid ID so that it does not correspond to any internal node
      //HW0: YOUR CODE HERE
    }
      
    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return get_IN().coordinate;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return Node_ID;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (Node_ID == n.index() && graph_container==n.graph_container)
          return true;
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
      // HW0: YOUR CODE HERE
      // Ordering nodes by their index
      if (Node_ID < n.index())
          return true;
      return false;
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
      size_type Node_ID;        // The Node's unique identification number
      
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
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    // Create a new internal node and set the position and ID of that internal node
    Internal_Node temp_node; 
    temp_node.coordinate = position;
    temp_node.IN_id = size_nodes; 

    //Add that internal node to the end of the internal nodes vector
    vec_internal_nodes.push_back(temp_node);

    //Increment the size of the graph
    ++size_nodes;
  
    return Node(this,size_nodes-1);       
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    //Since no individual node is deleted, the ID of this node must be smaller than the number of nodes of the graph. 
    if (n.index()<size())
      return true;
    return false;
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
    assert (i<size());
    return Node(this,i);   
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
      //Edges that span between the same two nodes should have the same ID
      //There are no duplicate Edges
      if (Edge_ID == e.Edge_ID)
        return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // No interpretive meanting to the ordering of edge. Order edges by indices
      if (Edge_ID < e.Edge_ID)
        return true;
      return false;
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
       *The corresponding internal node should be at index [Node ID]
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
    assert (i<num_edges());
    return Edge(this,i);
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
    if (has_node(a) and has_node(b)){
        // Decide order of the nodes. Node 1 has the smaller node by index.
        Node SmallNode;
        Node LargeNode;
        if (a<b) {
          SmallNode = a;
          LargeNode = b;
        } else {
          SmallNode = b;
          LargeNode = a;
        }

        // See if the first node is in the map at all. If it is, then go the vector of edges that correspond to that node. Iterate through those edges to see if those edge contains the index of the other node.
        if (Node_to_Edge_map.find(SmallNode.index())!=Node_to_Edge_map.end()){
          std::vector<size_type> Edges(Node_to_Edge_map.at(SmallNode.index()));
          for (size_type i=0; i<Edges.size(); i++){
            if (vec_internal_edges[Edges[i]].N2 == LargeNode){
              return true;
            }
          }
        }
        return false;
    }else{
        return false;
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
      // Make sure Node a and Node b are valid nodes of the graph
      if (has_node(a) and has_node(b)){
      
        // Decide order of the nodes. Node 1 has the smaller node by index.
        Node SmallNode;
        Node LargeNode;
        if (a<b) {
          SmallNode = a;
          LargeNode = b;
        } else {
          SmallNode = b;
          LargeNode = a;
        }
      
        /*See if the graph already has this edge, if it does, simply return the edge by looking it up in the map.
        Find the value in the map associated with the first Node index.
        Loop through the map value and see if the Edge associated with that index has the second matching Node*/
        if (has_edge(SmallNode,LargeNode)){
          std::vector<size_type> Edges(Node_to_Edge_map[SmallNode.index()]);
          for (size_type i=0; i<Edges.size(); i++){
            if (vec_internal_edges[Edges[i]].N2 == LargeNode){
              return Edge(this,Edges[i]);
            }
          }
        }
        
        //If the graph does not have this edge, must create new Internal Edge and add this edge to the two containers that keep track of edges
        // Create a new internal edge and set the position and ID of that internal edge
        Internal_Edge temp_internal_edge;
        temp_internal_edge.N1 = SmallNode; // First node should be the smaller node
        temp_internal_edge.N2 = LargeNode; // Second node should be the larger node
        temp_internal_edge.IE_id = size_edges;

        //Add that internal edge to the end of the internal edge vector that orders by index
        vec_internal_edges.push_back(temp_internal_edge);

        //Add this edge to the map that keeps track of nodes to edges
        Node_to_Edge_map[SmallNode.index()].push_back(temp_internal_edge.IE_id);

        //Increment the size of the graph
        ++size_edges;
      
        return Edge(this,size_edges-1);
      }
      else{
          //Return invalid edge if nodes do not exist
          return Edge();
      }
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
    Node_to_Edge_map.clear();
    // HW0: YOUR CODE HERE
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
    
    //Internal Nodes of the graph. Contains the actual point data
    struct Internal_Node {
        Point coordinate;     // A Node's cartesian coordinate
        size_type IN_id;  // The unique identification for an Internal Node
    };
    
    //Internal Edges of the graph. Contains the nodes of each edge and an edge ID
    struct Internal_Edge {
      Node N1; //First Node, should be the smaller node by ID
      Node N2; //Second Node, should be the larger node by ID
      size_type IE_id; //The unique identification for an Internal Edge
    };
  
    // Attributes of the graph class
    
    std::vector<Internal_Edge> vec_internal_edges;
    std::vector<Internal_Node> vec_internal_nodes;
    std::unordered_map<size_type, std::vector<size_type>> Node_to_Edge_map;
    size_type size_nodes;
    size_type size_edges;
};

#endif // CME212_GRAPH_HPP
