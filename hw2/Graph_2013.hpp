// CME 212 Winter 2019 HW2: Graph.hpp
// Brian Ryu
#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <vector>
#include <unordered_map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
 private:
  // HW0: YOUR CODE HERE

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
  // Now class template with type V
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  // Now class template with node type E
  // Using using instead of typedef to be consistent with other definitions
  using edge_value_type = E;

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
  class Node : private totally_ordered<Node> {
    // For homework 1: Node is totally ordered now
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

    // Added for HW2
    /** Return the node's position
     *  @return The position of current node
     *
     * Complexity: O(1) 
     */ 
    Point& position() {
      return (*homeGraphPtr).pointVector[nodeId()];
    }

    /** Return the node's position
     *  @return The position of current node
     *
     * Complexity: O(1) 
     */ 
    const Point& position() const {
      // HW0: YOUR CODE HERE
      // Nodes and Points are index-matched.
      // Return the Point with the same index in the Same graph
      return (*homeGraphPtr).pointVector[nodeId()];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return nodeIndex;
    }

    // HW1: YOUR CODE HERE
    /** Return the node value contained in a node as a reference.
     * @return The reference to the current node value
     *
     * @tparam node_value_type Type of node value variable used in this graph
     *
     * Complexity: O(1)
     */
    node_value_type& value() {
      // Access current graph's vector of nodes at the current index.
      return (*homeGraphPtr).nodeValueVector[nodeId()];
    }

    /** Return the node value contained in a node as a constant reference.
     * @return The reference to the current node value
     *
     * @tparam node_value_type Type of node value variable used in this graph
     *
     * Complexity: O(1)
     */
    const node_value_type& value() const {
      // Access current graph's vector of nodes at the current index.
      return (*homeGraphPtr).nodeValueVector[nodeId()];
    }

    /** Return the number of incident edges to the current node.
     *  @return the number of edges connected to current node
     *
     *  Complexity: O(1)
     */
    size_type degree() const {
      return (*homeGraphPtr).unsortedAdjList[nodeIndex].size();
    }

    /** Create a IncidentIterator that points to the first incident edge
     *   This is a constant function that does not modify the graph instance 
     *  @return IncidentIterator points to the first incident edge
     *
     *  Complexity: O(1)
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(nodeIndex, 0, homeGraphPtr);
    }

    /** Create a IncidentIterator that points to once space after 
     *   last incident edge. This is a constant function that 
     *   does not modify the graph instance 
     *  @return IncidentIterator points to one space after the last
     *          incident edge.
     *
     *  Complexity: O(1)
     */
    incident_iterator edge_end() const {
      int endIdx = (*homeGraphPtr).unsortedAdjList[nodeIndex].size();
      return IncidentIterator(nodeIndex, endIdx, homeGraphPtr);
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
      // Compare indices, assuming that nodes belong to the same graph
      // Return true if < is true for index
      if (index() < n.index()) {
        return true;
      } 
      if (homeGraphPtr < n.homeGraphPtr){
        return true;
      }
      return false; // False otherwise
    }

    /** Return this node's unique nodeID, a number in the range [0, N_add)
     *  where N_add is the total number of nodes that has been added to 
     *  the current graph */
    int nodeId() const {
      return (*homeGraphPtr).nodeIds[index()];
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
   * @param[in] inputNodeValue The contained value of new node
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   * @post result_node.value() == @a inputNodeValue
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& inputNodeValue = node_value_type()) {
    // HW0: YOUR CODE HERE
    // Push position and value to position and value vectors
    nodeIds.push_back(nodeIDCount); // O(1)
    pointVector.resize(nodeIDCount+1);
    pointVector[nodeIDCount] = position; // O(1)

    nodeValueVector.resize(nodeIDCount+1);
    nodeValueVector[nodeIDCount] = inputNodeValue; // O(1)


    // Add "empty" node, set index, and set position
    nodeVector.push_back(Node()); // O(1)
    nodeVector[nodeVector.size()-1].nodeIndex = nodeVector.size()-1; // O(1)
    nodeVector[nodeVector.size()-1].homeGraphPtr = this; // O(1)

    // Added for HW1 Problem2:
    nodeIDCount++;

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
  class Edge : private totally_ordered<Edge> {
    // For homework 1: Edge is totally ordered now
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // Note: Nothing Added Here
    }

    /** Return a node of this Edge */
    Node node1() const {
      return begNode;
    }

    /** Return the other node of this Edge */
    Node node2() const {
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

    // Added for HW2 Problem 2
    /** Return the Euclidean length of the edge
     *  @return the Euclidean length of current edge
     *
     *  Complexity: O(1)
     */
    double length() const {
      return norm(node1().position()-node2().position());
    }

    /* Return the node value contained in an edge as a reference.
     * @return The reference to the current edge value
     *
     * @tparam edge_value_type Type of edge value variable used in this graph
     *
     * Complexity: O(1)
     */
    edge_value_type& value() {
      // Get node indices
      int n1Idx = node1().nodeIndex;
      int n2Idx = node2().nodeIndex;
      // and then grab ID's that are tied to edge values: O(1)
      int lowId, highId;
      if ((*(node1().homeGraphPtr)).nodeIds[n1Idx] 
          < (*(node1().homeGraphPtr)).nodeIds[n2Idx]){
        lowId = (*(node1().homeGraphPtr)).nodeIds[n1Idx];
        highId = (*(node1().homeGraphPtr)).nodeIds[n2Idx];
      } else {
        highId = (*(node1().homeGraphPtr)).nodeIds[n1Idx];
        lowId = (*(node1().homeGraphPtr)).nodeIds[n2Idx];
      }
      // Retrun Value from hash map: O(1) expected
      return (*(node1().homeGraphPtr)).edgeValues[lowId][highId];
    }

    /* Return the node value contained in an edge as a reference.
     * @return The reference to the current edge value
     *
     * @tparam edge_value_type Type of edge value variable used in this graph
     *
     * Complexity: O(1)
     */
    const edge_value_type& value() const {
      // Get node indices
      int n1Idx = node1().nodeIndex;
      int n2Idx = node2().nodeIndex;

      // and then grab ID's that are tied to edge values: O(1)
      int lowId, highId;
      if ((*(node1().homeGraphPtr)).nodeIds[n1Idx] 
          < (*(node1().homeGraphPtr)).nodeIds[n2Idx]){
        lowId = (*(node1().homeGraphPtr)).nodeIds[n1Idx];
        highId = (*(node1().homeGraphPtr)).nodeIds[n2Idx];
      } else {
        highId = (*(node1().homeGraphPtr)).nodeIds[n1Idx];
        lowId = (*(node1().homeGraphPtr)).nodeIds[n2Idx];
      }
      // Retrun Value from hash map: O(1) expected
      return (*(node1().homeGraphPtr)).edgeValues[lowId][highId];
    }
   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

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
        // std::cout<<counter<<std::endl;
    }
    return counter;
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

    if (a == b) {
      // No loops allowed. Fast check!
      return false;
    }
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

    
    if ((unsortedAdjList.size()-1) < lowIdx) {
      // If the index of last element is less than lowIdx
      // Of adjacency list. Then edge should not be there
      return false;
    }
    // Iterate through the elements of a sub-adjacency list.
    // Worst case: O(num_edges).
    // If edges are evenly connected to nodes (uniformly):
    // We expect about O(num_edges/num_nodes)
    for (size_t i = 0 ; i < unsortedAdjList[lowIdx].size(); i++) {
      //Iterate through sub-adjacency list
      if (unsortedAdjList[lowIdx][i].second==tempEdge) { //O(1)
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
    resizeAdjList();
    
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

    std::pair<int, Edge> revAdd;
    revAdd.first = lowIdx;
    revAdd.second = tempEdge;
    unsortedAdjList[highIdx].push_back(revAdd);
    unsortedAdjList[lowIdx].push_back(toAdd);
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
    // Simply clear all data
    nodeVector.clear();
    nodeValueVector.clear();
    pointVector.clear();
    adjList.clear();
    unsortedAdjList.clear();
    edgeValues.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
    // For homework 1: NodeIterator is totally ordered
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
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Returns the Node the current iterator is pointing at.
     *  @return The Node of the current iterator is pointing at
     *      as a constant copy
     *
     *  Complexity: O(1)
     */
    Node operator*() const {
      // Access current graph's vector of nodes at current index.
      return (*homeGraphPtr).nodeVector[nodeIdx];
    }

    /** Increments the node index of current iterator by 1
     *  @return The iterator itself, advnaced by one node index
     * 
     *  @post old nodeIdx + 1 == new nodeIdx
     *
     *  Complexity: O(1)
     */
    NodeIterator& operator++() {
      nodeIdx++;
      return *this;
    }

    /** Compares whether two node Iterators point at the same node
     *    Implicitly assume that iterators are from same graph instance.
     *  @return True if at same position, false otherwise
     *
     *  Complexity: O(1)
     */
    bool operator==(const NodeIterator& inputItr) const {
      // Compare node index only.
      return nodeIdx == inputItr.nodeIdx;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    int nodeIdx; // Keep tracks of index
    const Graph* homeGraphPtr; // The graph iterator belongs to

    /** Constructor for NodeIterator based on which graph and index
     *  @return the NodeIterator of requested graph and index
     */
    NodeIterator(int x, const Graph* inHomeGraphPtr){
      nodeIdx = x; 
      homeGraphPtr = inHomeGraphPtr;
    }
  };

  // HW1 #2: YOUR CODE HERE
  /** Create a NodeIterator that points to the first node
   *   This is a constant function that does not modify the graph instance 
   *  @return NodeIterator points to the first node
   *
   *  Complexity: O(1)
   */
  node_iterator node_begin() const {
    return NodeIterator(0, this);
  }

  /** Create a NodeIterator that points to one space after the last node
  *   This is a constant function that does not modify the graph instance
   *  @return NodeIterator points to one space after the last node
   *
   *  Complexity: O(1)
   */
  node_iterator node_end() const {
    return NodeIterator(nodeVector.size(), this);
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>  {
    // For homework 1: IncidentIterator is totally ordered
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
    /** Returns the Edge the current iterator is pointing at, 
     *  while setting the current node as node1() and connected
     *  node as node2()
     *
     *  @return The Edge of the current iterator is pointing at
     *      as a constant copy, with node1 and node2 set such that
     *      the current node is node1
     *
     *  Complexity: O(1)
     */
    Edge operator*() const {
      // first grab edge current iterator is pointing at
      Edge currEdge= (*homeGraphPtr).unsortedAdjList[nodeIdx][edgeIdx].second;

      Edge toReturn;
      // Check and flip nodes 1 and 2, if necessary
      if (currEdge.node1().nodeIndex == nodeIdx) {
        toReturn.begNode = currEdge.node1();
        toReturn.endNode = currEdge.node2();
      } else {
        toReturn.begNode = currEdge.node2();
        toReturn.endNode = currEdge.node1();
      }

      return toReturn;
    }

    /** Increments the edge index of current iterator by 1
     *  @return The iterator itself, advnaced by one edge index
     * 
     *  @post old edgeIdx + 1 == new edgeIdx
     *
     *  Complexity: O(1)
     */
    IncidentIterator& operator++() {
      edgeIdx++;
      return *this;
    }

    /** Compares whether two incident edge Iterators point at the same edge
     *  Checks for current node, edge index, and whether the iterator
     *  belongs to the same graph.
     *  @return True if at same position and graph, false otherwise
     *
     *  Complexity: O(1)
     */
    bool operator==(const IncidentIterator& inputItr) const {
      // Sequentially check three conditions, returning false if
      // any fails.
      if (edgeIdx != inputItr.edgeIdx) return false;
      if (nodeIdx != inputItr.nodeIdx) return false;
      if (homeGraphPtr != inputItr.homeGraphPtr) return false;

      return true;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
      int nodeIdx; // Index of current node
      int edgeIdx; // Index of currently pointing edge
      const Graph* homeGraphPtr; // The graph iterator belongs to

      /** Constructor for IncidentIterator based on which graph and indices
       *  @return the IncidentIterator of requested graph and index
       */
      IncidentIterator(int nIdx, int eIdx, const Graph* inHomeGraphPtr) {
        nodeIdx = nIdx;
        edgeIdx = eIdx;
        homeGraphPtr = inHomeGraphPtr;
      }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator> {
    // For homework 1: EgeIterator is totally ordered
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
    /** Returns the Edge the current iterator is pointing at.
     *  @return The Edge of the current iterator is pointing at
     *      as a constant copy
     *
     *  Complexity: O(1)
     */
    Edge operator*() const {
      return (*homeGraphPtr).edge(edgeIdx);
    }

    /** Increments the edge index of current iterator by 1
     *  @return The iterator itself, advnaced by one edge index
     * 
     *  @post old edgeIdx + 1 == new edgeIdx
     *
     *  Complexity: O(1)
     */
    EdgeIterator& operator++() {
      edgeIdx++;
      return *this;
    }

    /** Compares whether two edge Iterators point at the same node
     *    Implicitly assume that iterators are from same graph instance.
     *  @return True if at same position, false otherwise
     *
     *  Complexity: O(1)
     */
    bool operator==(const EdgeIterator& inputItr) const {
      if((homeGraphPtr == inputItr.homeGraphPtr) 
          && (edgeIdx == inputItr.edgeIdx)) {
        return true;
      }
      return false;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const Graph* homeGraphPtr; // Which graph the iterator belongs to
    int edgeIdx; // Index of edge

    /** Constructor for EdgeIterator based on which graph and index
     *  @return the EdgeIterator of requested graph and index
     */
    EdgeIterator(int x, const Graph* inHomeGraphPtr) {
      edgeIdx = x;
      homeGraphPtr = inHomeGraphPtr;
    }
  };

  // HW1 #5: YOUR CODE HERE
  /** Create a EdgeIterator that points to the first Edge
   *   This is a constant function that does not modify the graph instance 
   *  @return EdgeIterator points to the first edge
   *
   *  Complexity: O(1)
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(0, this);
  }

  /** Create a EdgeIterator that points to one space after the last edge
  *   This is a constant function that does not modify the graph instance
   *  @return EdgeIterator points to one space after the last edge
   *
   *  Complexity: O(1)
   */
  edge_iterator edge_end() const {
    return EdgeIterator(num_edges(), this);
  }  

 
 public:
  /** Remove a node from the current graph. All edges incident
   *   to that node are removed/invalidated
   *  @param[in] n Node to remove from graph
   *  @return 1 if successfully removed, 0 if node did not exist
   *  @post new size() = old size() - 1 
   *  @post nodeId of removed node is invalidated
   *  @post new node(i) = old node(size(i)-1) where i is 
   *   the index of the removed node. That is, the node that had the last
   *   index prior to deletion of node i, now occupies index i
   *  @post All incident iterators pertaining to the removed node
   *   are invalidated
   *
   *  Complexity: O(num_nodes) Assuming a sparse graph
   *              If not sparse, O(num_nodes*degree of node)
   */
  size_type remove_node(const Node& n) {
    if (has_node(n)) {
      // Resize to prevent error
      resizeAdjList();

      // Remove all incident edges
      auto begItr = n.edge_begin();
      auto endItr = n.edge_end();
      std::vector<int> incidentNodes;
      for (; begItr != endItr ; ++begItr) {
          incidentNodes.push_back((*begItr).node2().index());
      }
      // Now know all incident edge indices Remove!
      while(n.edge_begin() != n.edge_end()) {
        remove_edge(*(n.edge_begin()));
        
      }

      // Done removing incident edges. Now initiate swapping
      int lastIndex = size() - 1; // Remember index of last element
      int currIndex = n.index();

      // Swap with very last element.
      // To Swap: nodeVector, adjList, unsortedAdjList
      std::swap(nodeVector[lastIndex], nodeVector[currIndex]);
      std::swap(adjList[lastIndex], adjList[currIndex]);
      std::swap(unsortedAdjList[lastIndex], unsortedAdjList[currIndex]);
      std::swap(nodeIds[lastIndex], nodeIds[currIndex]);

      // Delete last of the nodes
      nodeVector.pop_back();
      adjList.pop_back();
      unsortedAdjList.pop_back();
      nodeIds.pop_back();
      
      // Updating index and edge information
      //  Index of previously last node
      nodeVector[currIndex].nodeIndex = currIndex;
      //  Iterate through incident edges of previously last node, update index.
      std::pair<int, Edge> currentEdge;

      for (size_t i = 0 ; i < adjList[currIndex].size(); i++) {
        //No need to change first element, because it is the one this edge is connected to
        currentEdge = adjList[currIndex][i];
        int connectedIdx = currentEdge.first;
        // Change in this sublist
        Edge toChange;
        
        toChange.begNode = nodeVector[currIndex];
        toChange.endNode = nodeVector[connectedIdx];
        adjList[currIndex][i].second = toChange; 
      }

      for (size_t i = 0 ; i < unsortedAdjList[currIndex].size(); i++) {
        //No need to change first element, because it is the one this edge is connected to
        currentEdge = unsortedAdjList[currIndex][i];
        int connectedIdx = currentEdge.first;
        // Change in this sublist
        Edge toChange;
        toChange.begNode = nodeVector[currIndex];
        toChange.endNode = nodeVector[connectedIdx];
        unsortedAdjList[currIndex][i].second = toChange; 

        // We know which one this is connected to. Change
        for (size_t j = 0 ; j < unsortedAdjList[connectedIdx].size(); j++) {
          if (unsortedAdjList[connectedIdx][j].first == lastIndex) {
            unsortedAdjList[connectedIdx][j].first = currIndex;
            unsortedAdjList[connectedIdx][j].second = toChange; 
          }
        }
        for (size_t j = 0 ; j < adjList[connectedIdx].size(); j++) {
          if (adjList[connectedIdx][j].first == lastIndex) {

            adjList[connectedIdx][j].first = currIndex;
            adjList[connectedIdx][j].second = toChange; 
          }
        }

      }
      return 1;
    }
    return 0;
  }

  /** Remove a node from the current graph. All edges incident
   *   to that node are removed/invalidated
   *  @param[in] n_it Node iterator poining to the node to be deleted
   *  @return node iterator pointing to the node at the 
   *   same index after deletion
   *  @post new size() = old size() - 1 
   *  @post nodeId of removed node is invalidated
   *  @post new node(i) = old node(size(i)-1) where i is 
   *   the index of the removed node. That is, the node that had the last
   *   index prior to deletion of node i, now occupies index i
   *  @post All incident iterators pertaining to the removed node
   *   are invalidated
   *
   *  Complexity: O(num_nodes) Assuming a sparse graph
   *              If not sparse, O(num_nodes*degree of node)
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return n_it;
  }

  /** Remove an edge from the current graph
   *  @param[in] n1 A node connected to the edge to be deleted
   *  @param[in] n2 The other node connected to the edge to be deleted
   *  @return 1 if successfully removed, 0 if edge did not exist
   *  @post new num_edges() = old num_edges() - 1 
   *  @post edge value of the removed edge is invalidated
   *  @post All edge iterators pertaining to the removed edge
   *   are invalidated
   *
   *  Complexity: O(degree of nodes connected to removed edge)
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    if (has_edge(n1, n2)) {
      // Grab indices of the two nodes
      int lowIdx, highIdx;
      if (n1.index() < n2.index()) {
        lowIdx = n1.index();
        highIdx = n2.index();
      } else {
        lowIdx = n2.index();
        highIdx = n1.index();
      }

      // Remove Edge from adjacency list
      for (size_t i = 0 ; i < adjList[lowIdx].size(); i++) {
        if (adjList[lowIdx][i].first == highIdx) {
          adjList[lowIdx].erase(adjList[lowIdx].begin() + i);
          break;
        }
      }
      for (size_t i = 0 ; i < adjList[highIdx].size(); i++) {
        if (adjList[highIdx][i].first == lowIdx) {
          adjList[highIdx].erase(adjList[highIdx].begin() + i);
          break;
        }
      }
      // Remove Edge from unsorted adjacency list
      for (size_t i = 0 ; i < unsortedAdjList[lowIdx].size(); i++) {
        if (unsortedAdjList[lowIdx][i].first == highIdx) {
          unsortedAdjList[lowIdx].erase(unsortedAdjList[lowIdx].begin() + i);
          break;
        }
      }
      for (size_t i = 0 ; i < unsortedAdjList[highIdx].size(); i++) {
        if (unsortedAdjList[highIdx][i].first == lowIdx) {
          unsortedAdjList[highIdx].erase(unsortedAdjList[highIdx].begin() + i);
          break;
        }
      }
      return 1;
    }
    return 0;
  }

  /** Remove an edge from the current graph
   *  @param[in] e edge to be removed
   *  @return 1 if successfully removed, 0 if edge did not exist
   *  @post new num_edges() = old num_edges() - 1 
   *  @post edge value of the removed edge is invalidated
   *  @post All edge iterators pertaining to the removed edge
   *   are invalidated
   *
   *  Complexity: O(degree of nodes connected to removed edge)
   */
  size_type remove_edge(const Edge& e) {
    const Node n1 = e.node1();
    const Node n2 = e.node2();
    return remove_edge(n1, n2);
  }

  /** Remove an edge from the current graph
   *  @param[in] e_it Edge iterator poining to the node to be deleted
   *  @return 1 if successfully removed, 0 if edge did not exist
   *  @post new num_edges() = old num_edges() - 1 
   *  @post edge value of the removed edge is invalidated
   *  @post All edge iterators pertaining to the removed edge
   *   are invalidated
   *
   *  Complexity: O(degree of nodes connected to removed edge)
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return e_it;
  }
  
private:
  std::vector<Node> nodeVector;
  std::vector<node_value_type> nodeValueVector;
  std::vector<Point> pointVector;
  std::vector<int> nodeIds;
  int nodeIDCount= 0;

  /** Private method that resize adjacency list to prevent segfault
   *  error. Only resizes if the list is smaller than the size of
   *  vector containing nodes
   */
  void resizeAdjList() {
    if (adjList.size() != nodeVector.size()) {
        adjList.resize(nodeVector.size());
        unsortedAdjList.resize(nodeVector.size());
    }
    return;
  }

  // One adjacency list for fast index access
  std::vector<std::vector<std::pair<int, Edge>>> adjList;
  // Reverse-adjacency list for traversing incident edges
  std::vector<std::vector<std::pair<int, Edge>>> unsortedAdjList;

  // Nested hashmap for edge values
  std::unordered_map<int, std::unordered_map<int, edge_value_type>> edgeValues;
};

#endif // CME212_GRAPH_HPP