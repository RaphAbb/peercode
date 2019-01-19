// CME 212 Winter 2019 HW0: Graph.hpp
// Brian Ryu
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

  // Note: Nothing Added Here

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

    // Note: Nothing Added Here
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
      // HW0: YOUR CODE HERE

      // Note: Nothing Added Here
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      // Nodes and Points are index-matched.
      // Return the Point with the same index in the Same graph
      return (*homeGraphPtr).pointVector[nodeIndex];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      // Node knows its index
      return nodeIndex;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      // If (Belong to same graph && same index)
      if ( homeGraphPtr == n.homeGraphPtr   && index() == n.index()) {
        // Then return true
        return true;
      }
      // Otherwise,return false
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
      // My rule:
      // First compare indices, then x,y,and z positions in that order

      // First check index
      if (index() < n.index()) {
        return true;
      } else if (index() > n.index()) {
        return false;
      }

      // If index is the same for some reason:
      // Compare x:
      if (position().x < n.position().x) {
        return true;
      } else if (position().x > n.position().x) {
        return false;
      }
      // If x is the same, compare y.
      if (position().y < n.position().y) {
        return true;
      } else if (position().y > n.position().y) {
        return false;
      }
      // If y is the same, compare z.
      if (position().z < n.position().z) {
        return true;
      } else if (position().z > n.position().z) {
        return false;
      }

      // False if x, y, z, and index are the same
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Max of 16 bytes
    int nodeIndex; // Index of the node; 8 bytes
    Graph* homeGraphPtr; // Pointer to the graph it belongs to; 8 bytes
  }; 

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    // We have a vector of nodes. Use size()
    return nodeVector.size(); // O(1)
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

    // Push position to position vector
    pointVector.push_back(position); // O(1)

    // Add "empty" node, set index, and set position
    nodeVector.push_back(Node()); // O(1)
    nodeVector[nodeVector.size()-1].nodeIndex = nodeVector.size()-1; // O(1)
    nodeVector[nodeVector.size()-1].homeGraphPtr = this; // O(1)

    return nodeVector[nodeVector.size()-1];   
    // Total O(1) complexity
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (nodeVector[n.index()] == n) {
      // If node with n.index is the same as n
      // vector index access is O(1), == is O(1)
      return true;
    }
    // Otherwise false
    return false;
    // Total O(1) complexity
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    // nodeVector index access. Total O(1) complexity
    return nodeVector[i];
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

      // Note: Nothing Added Here
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE

      // Return beginning node
      return begNode;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE

      // Return ending node
      return endNode;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // If first nodes are the same and second nodes are the same
      // i.e. "Same" direction, though undirected
      if (node1() == e.node1() && node2() == e.node2()) {
        // Return true
        return true;
      }
      // If first and end, and end and first are the same
      // i.e. "Opposite" direction
      if (node1() == e.node2() && node2() == e.node1()) {
        // Return true
        return true;
      }
      // Otherwise, false
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // First compare node 1
      if (node1() < e.node1()) {
        return true;
      } else if (e.node1() < node1()) {
        return false;
      }

      // If the same, compare node 2
      if (node2() < e.node2()) {
        return true;
      } else if (e.node2() < node2()) {
        return false;
      }
      // If all nodes the same, return false
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Max of 32 bytes
    Node begNode; // Beginning node, or node1(). 16 bytes
    Node endNode; // Ending node, or node2(). 16 bytes

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    
    // We iterate through the adjacency list, see how many edges a node has
    // Compute the cumulative sum
    size_t counter = 0;
    for (size_t i = 0 ; i < adjList.size() ; i++) {
        counter+= adjList[i].size();
    }
    return counter;
    // For loop iterating over all nodes: O(num_nodes)
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    
    // What I do here is iterate through all nodes, compute the cumulative sum
    // by calling .size() of each sub-adjacency list and incrementing
    // If the cumulative sum exceeds the index, I go back one step, and find
    // Where the element is.

    // Use two counters for this
    size_t counter1 = 0;
    size_t counter2 = 0;

    // Iterate through nodes O(num_nodes)
    for (size_t j = 0 ; j < adjList.size(); j++) {

      counter1 += adjList[j].size(); // Increment cumulative sum O(1)
      // std::cout<< "counter1: " << counter1 << std::endl;
      if (i < counter1) {
        // If we went further than the input index
        // See where our Edge is in the sub-adjacency list
        size_t difference = i - counter2; 
        // Return that
        return adjList[j][difference].second;
      }
      counter2 += adjList[j].size();  // Increment cumulative sum O(1)
    }

    // If a valid index is called, we should not come here, but
    // Return an invalid edge, if an invalid index is inputted.
    // This will result in a runtime error under normal circumstances
    return Edge();

    // Total complexity: Worst case O(num_nodes()). 
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    
    // Create an edge, such that the smaller node is the begNode.
    Edge tempEdge = Edge();
    tempEdge.begNode = Node();
    tempEdge.endNode = Node();
    size_t lowIdx, highIdx;

    if (a.nodeIndex < b.nodeIndex) {
      lowIdx = a.nodeIndex;
      tempEdge.begNode.nodeIndex = lowIdx;
      tempEdge.begNode.homeGraphPtr = a.homeGraphPtr;
      
      highIdx = b.nodeIndex;
      tempEdge.endNode.nodeIndex = highIdx;
      tempEdge.endNode.homeGraphPtr = b.homeGraphPtr;

    } else {
      lowIdx = b.nodeIndex;
      tempEdge.begNode.nodeIndex = lowIdx;
      tempEdge.begNode.homeGraphPtr = b.homeGraphPtr;
      
      highIdx = a.nodeIndex;
      tempEdge.endNode.nodeIndex = highIdx;
      tempEdge.endNode.homeGraphPtr = a.homeGraphPtr;
    }
    
    if ((adjList.size()-1) < lowIdx) {
      // If the index of last element is less than lowIdx
      // Of adjacency list. Then edge should not be there
      return false;
    }
    // Iterate through the elements of a sub-adjacency list.
    // Worst case: O(num_edges).
    // If edges are evenly connected to nodes (uniformly):
    // We expect about O(num_edges/num_nodes)
    for (size_t i = 0 ; i < adjList[lowIdx].size(); i++) {
      //Iterate through sub-adjacency list
      if (adjList[lowIdx][i].second==tempEdge) { //O(1)
        // Found, return true
        return true;
      }
    }
    // Not found
    return false;
    // Worst case complexity: O(num_edges)
    // However, this happens ONLY when all edges connected one node.
    // If edges are balanced (uniform distribution of edges 
    // connected to nodes), expected runtime should be 
    // O((num_edges)/(num_nodes)).

  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges () == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    
    // first resize adjacency list if it is too small
    if (adjList.size() != nodeVector.size()) {
      adjList.resize(nodeVector.size());
    }
    
    // Create temporary edge so we can put in
    Edge tempEdge = Edge();
    tempEdge.begNode = Node();
    tempEdge.endNode = Node();
    size_t lowIdx, highIdx;

    if (a.nodeIndex < b.nodeIndex) {
      lowIdx = a.nodeIndex;
      tempEdge.begNode.nodeIndex = lowIdx;
      tempEdge.begNode.homeGraphPtr = a.homeGraphPtr;
      
      highIdx = b.nodeIndex;
      tempEdge.endNode.nodeIndex = highIdx;
      tempEdge.endNode.homeGraphPtr = b.homeGraphPtr;

    } else {
      lowIdx = b.nodeIndex;
      tempEdge.begNode.nodeIndex = lowIdx;
      tempEdge.begNode.homeGraphPtr = b.homeGraphPtr;
      
      highIdx = a.nodeIndex;
      tempEdge.endNode.nodeIndex = highIdx;
      tempEdge.endNode.homeGraphPtr = a.homeGraphPtr;
    }

    // Check if edge is already there
    // Very worst case: O(N_e) = O(num_edges) complexity
    // This worst case is extreme, and should be better typically
    if (has_edge(a, b) == true) {
      return tempEdge; // O(1)
    }

    // Create a pair to add to adjacency list
    std::pair<int, Edge> toAdd;
    toAdd.first = highIdx;
    toAdd.second = tempEdge;
    adjList[lowIdx].push_back(toAdd);
    return tempEdge;

    // Absolute worst case complexity: O(N_e) = O(num_edges)
    // Should typically be better
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    // Simply clear vectors
    nodeVector.clear();
    adjList.clear();
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  // Two vectors, containing Nodes, and P oints.
  std::vector<Node> nodeVector;
  std::vector<Point> pointVector;

  // One adjacency list
  std::vector<std::vector<std::pair<int, Edge>>> adjList;

};

#endif // CME212_GRAPH_HPP
