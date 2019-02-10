#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <deque>
#include <vector>


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

  std::vector<Point> nodeLocations;
  std::vector<V> nodeValues;

  unsigned numEdges;
  std::vector<std::deque<unsigned>> incidentArray;
  std::vector<unsigned> smallerNodesConnected;
 

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
  /** Declaration of what type of data the node represents. */
  using node_value_type = V;

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
  Graph() : numEdges(0) {
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
      // Return a default initialized node.
    }

    /** Return this node's position. */
    const Point& position() const {
      return this->graphContainer->nodeLocations[this->nodeIndex];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return this->nodeIndex;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return a non-const reference to the value of this node.
     *  
     *  @pre The node that this is called upon must be valid. (The node is
     *       contained in some graph.)
     */
    node_value_type& value() {
      this->assert_valid();
      // We temporarily suppress the const modifier we placed on the Graph
      // pointer so that we can appropriately modify the value.
      return const_cast<Graph*>(this->graphContainer)->nodeValues[this->nodeIndex];
    }

    /** Return a const reference to the value associated with this node.
     *
     *  @pre The node that this is called upon must be valid. (The node is
     *       contained in some graph.)
     */
    const node_value_type& value() const {
      this->assert_valid();
      return this->graphContainer->nodeValues[this->nodeIndex];
    }

    /** Return the degree of this node. (The number of nodes adjacent to it.)
     *
     *  @pre The node that this is called upon must be valid. (The node is
     *       contained in some graph.)
     */
    size_type degree() const {
      this->assert_valid();
      return this->graphContainer->incidentArray[this->nodeIndex].size();
    }

    /** Return the starting incident_iterator for this node.
     *  
     *  @pre The node that this is called upon must be valid. (The node is
     *       contained in some graph.)
     */
    incident_iterator edge_begin() const {
      this->assert_valid();
      return IncidentIterator(this->graphContainer, this->nodeIndex, 0);
    }
  
    /** Return the ending incident_iterator for this node.
     *  
     *  @pre The node that this is called upon must be valid. (The node is
     *       contained in some graph.)
     */
    incident_iterator edge_end() const {
      this->assert_valid();
      return IncidentIterator(this->graphContainer, this->nodeIndex, this->degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes are contained in the same graph and have the same index.
     */
    bool operator==(const Node& n) const {
      return (this->nodeIndex == n.nodeIndex);
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
      return (this->nodeIndex < n.nodeIndex);
    }

   private:
    // Allow Graph, Edge, and NodeIterator to access Node's private member data
    //  and functions.
    friend class Graph;
    friend class Edge;
    friend class NodeIterator;

    size_type nodeIndex;
    const Graph* graphContainer;

    // A legitimate constructor for Nodes that is private. Only the Graph,
    // Edge, and NodeIterator classes can call this object.
    Node(const Graph* graphPointer, size_type index) :
      nodeIndex(index), graphContainer(graphPointer) {
    }

    /** A function that asserts that the IncidentIterator is valid for
     *  operations concerning incrementing and dereferencing.
     *
     *  Checks that:
     *    - The node is associated with a legitimate graph
     *    - The node index is somewhere between [0, number of nodes in the
     *        graph).
     */
    void assert_valid() const {
      assert (this->graphContainer != nullptr);
      assert (this->nodeIndex >= 0);
      assert (this->nodeIndex < this->graphContainer->num_nodes());
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return this->nodeLocations.size();
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
  Node add_node(const Point& position, 
                const node_value_type& value = node_value_type()) {
    nodeLocations.push_back(position);
    nodeValues.push_back(value);

    incidentArray.push_back(std::deque<size_type>());
    smallerNodesConnected.push_back(0);

    return Node(this, this->num_nodes()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (this == n.graphContainer);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {

    // Check to see that the node index is appropriate.
    assert (i < this->num_nodes());

    return Node(this, i);
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
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(this->graphContainer, this->edgeNode1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      assert (this->graphContainer != nullptr);
      return Node(this->graphContainer, this->edgeNode2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {

      if (this->graphContainer != e.graphContainer) return false;

      size_type thisNode1 = this->edgeNode1;
      size_type thisNode2 = this->edgeNode2;
      size_type thatNode1 = e.edgeNode1;
      size_type thatNode2 = e.edgeNode2;      

      if (thisNode1 == thatNode1 and thisNode2 == thatNode2) return true;
      if (thisNode1 == thatNode2 and thisNode2 == thatNode1) return true;
      return false;

    }

    /** Test whether this edge and @a e are inequal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator!=(const Edge& e) const {
      return !(*this==e);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      size_type thisSmaller = std::min(this->edgeNode1, this->edgeNode2);
      size_type thisBigger  = std::max(this->edgeNode1, this->edgeNode2);
      size_type thatSmaller = std::min(e.edgeNode1, e.edgeNode2);
      size_type thatBigger  = std::max(e.edgeNode1, e.edgeNode2);

      if (thisSmaller < thatSmaller) return true;
      if (thisSmaller > thatSmaller) return false;
      if (thisBigger  < thatBigger ) return true;
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
   
    // Allow EdgeIterator to access Edge's private data and functions.
    friend class EdgeIterator;
 
    Edge(const Graph* graphPointer, size_type node1, size_type node2) :
     graphContainer(graphPointer), edgeNode1(node1), edgeNode2(node2) {
    }

    const Graph* graphContainer;
    size_type edgeNode1;
    size_type edgeNode2;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return numEdges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   * Currently complexity O(num_nodes() + d) where d is the largest degree of
   * nodes in the graph.
   */
  Edge edge(size_type i) const {

    bool validIndex = (i >= 0 and i < this->num_edges());
    assert (validIndex);
    
    // We iterate through the nodes to find the ith edge from a node with a
    // larger index to a node with a smaller index.
    for (size_type nodeIndex = 0; nodeIndex<this->num_nodes(); nodeIndex++) {

      // If the number of edges to smaller index nodes is greater than i, then
      // return the respective edge.
      size_type numNodesSmaller = this->smallerNodesConnected[nodeIndex];
      if (i < numNodesSmaller) {
        size_type otherNode = numNodesSmaller;
        return Edge(this, nodeIndex, otherNode);
      }

      // Otherwise, decrease the provided index by the number of edges.
      i -= numNodesSmaller;
    }

    assert (false); // With proper behavior, we should never get here.
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert (this->has_node(a));
    assert (this->has_node(b));

    size_type aIndex = a.index();
    size_type bIndex = b.index();

    std::vector<size_type>& incidentNodes = this->incidentArray[aIndex];

    for (size_type i = 0; i < incidentNodes.size(); i++) {
      if (incidentNodes[i] == bIndex) return true;
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

    assert (this->has_node(a));
    assert (this->has_node(b));
    assert (a.index() != b.index());

    size_type aIndex = a.index();
    size_type bIndex = b.index();

    if (aIndex < bIndex) {
      this->smallerNodesConnected[bIndex]++;
      this->incidentArray[aIndex].push_back(bIndex);
      this->incidentArray[bIndex].push_front(aIndex);
    } else {
      this->smallerNodesConnected[aIndex]++;
      this->incidentArray[aIndex].push_front(bIndex);
      this->incidentArray[bIndex].push_back(aIndex);
    }

    this->numEdges++;

    return Edge(this, aIndex, bIndex);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    numEdges = 0;
    nodeLocations.clear();
    incidentArray.clear();
    smallerNodesConnected.clear();
    nodeValues.clear();
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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Dereference operator for the node iterator.
     *
     * @pre the node must refer to a valid node within the graph.
     */
    Node operator*() const{
      this->assert_valid();
      return Node(this->graphContainer, this->index);
    }

    // Increment operator for the node iterator.
    NodeIterator& operator++(){
      this->assert_valid();
      this->index++;
      return *this;
    }

    // Equality operator for node iterators.
    bool operator==(const NodeIterator& otherIter) const {
      bool sameIndex = (this->index == otherIter.index);
      bool sameGraph = (this->graphContainer == otherIter.graphContainer);
      return (sameIndex and sameGraph);
    }

    // Inequality operator for node iterators.
    bool operator!=(const NodeIterator& otherIter) const {
      return !(*this == otherIter);
    }

   private:
    friend class Graph;

    // HW1 #2: YOUR CODE HERE
    size_type index;
    const Graph* graphContainer;

    // Main constructor for the NodeIterator class.
    NodeIterator(size_type index, const Graph* graphContainer) : index(index),
        graphContainer(graphContainer) {
    }

    /** A helper function to assert that this NodeIterator is valid.
     *
     * A NodeIterator is considered valid if:
     *   - It contains a pointer to a valid graph.
     *   - The node it refers to is valid.
     */
    void assert_valid() const {
      assert (this->graphContainer != nullptr);
      assert (this->index >= 0);
      assert (this->index < this->graphContainer->num_nodes());
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  // Initializer for the starting NodeIterator.
  node_iterator node_begin() const {
    return NodeIterator(0, this);
  }

  // Initializer for the ending NodeIterator.
  node_iterator node_end() const {
    return NodeIterator(this->num_nodes(), this);
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Dereference operator for this IncidentIterator.
     *
     * @pre the IncidentIterator must refer to a valid node.
     */
    Edge operator*() const {
      this->assert_valid();
      const size_type& cruxNode = this->cruxNode;
      const Graph& graph = *(this->graphContainer);
      size_type otherNode = graph.incidentArray[cruxNode][this->adjNodeIndex];
      return Edge(this->graphContainer, cruxNode, otherNode);
    }

    /** Increment operator for this IncidentIterator.
     *
     * @pre the IncidentIterator must refer to a valid node.
     */
    IncidentIterator& operator++() {
      this-> assert_valid();
      this->adjNodeIndex++;
      return *this;
    }

    // Equality operator for the IncidentIterator.
    bool operator==(const IncidentIterator& other) const {
      bool sameGraph = (this->graphContainer == other.graphContainer);
      bool sameCrux = (this->cruxNode == other.cruxNode);
      bool sameAdj = (this->adjNodeIndex == other.adjNodeIndex);
      return (sameGraph and sameCrux and sameAdj);
    }

    // Inequality operator for the IncidentIterator.
    bool operator!=(const IncidentIterator& other) const {
      return !(*this == other);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    const Graph* graphContainer;
    const size_type cruxNode;
    size_type adjNodeIndex;

    // Main constructor for the IncidentIterator class.
    IncidentIterator(const Graph* graphContainer, const size_type cruxNode,
        size_type adjNodeIndex) : graphContainer(graphContainer), 
        cruxNode(cruxNode), adjNodeIndex(adjNodeIndex) {
    }

    /** A function that asserts that the iterator is valid.
     *
     * Makes sure that the iterator this is called upon satisfies the following
     *   - The Incident Iterator is associated with a graph.
     *   - The adjacent index is somewhere in [0, degree of crux vertex).
     */
    void assert_valid() const {
      assert (this->graphContainer != nullptr);
      assert (this->cruxNode >= 0);
      assert (this->cruxNode < this->graphContainer->num_nodes());
      assert (this->adjNodeIndex >= 0);

      const Graph& graphRef = *(this->graphContainer);
      size_type numAdjacent = graphRef.incidentArray[this->cruxNode].size();
      assert (this->adjNodeIndex < numAdjacent);
    }
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Dereference operator for valid EdgeIterators.
     *
     * @pre this EdgeIterator must refer to a valid edge.
     */
    Edge operator*() const {
      this->assert_valid();
      const Graph& graph = *(this->graphContainer);

      size_type node1 = this->firstNode;
      size_type node2 = graph.incidentArray[node1][this->adjacencyIndex];
      return Edge(this->graphContainer, node1, node2);
    }

    /** Increment operator for valid edgeIterators.
     *
     * @pre this EdgeIterator must refer to a valid edge.
     */
    EdgeIterator& operator++() {
      this->assert_valid();
      const Graph& graph = *(this->graphContainer);
      
      // Increase the adjacency index of the first node.
      this->adjacencyIndex++;

      // If we adjacency index surpasses the number of smaller nodes, then
      // continue onto the next node with smaller nodes connected to it.
      size_type numNodesSmaller = graph.smallerNodesConnected[this->firstNode];
      if (this->adjacencyIndex == numNodesSmaller) {
        this->adjacencyIndex = 0;
        this->firstNode++;

        numNodesSmaller = graph.smallerNodesConnected[this->firstNode];
        while (numNodesSmaller==0 and this->firstNode!=graph.num_nodes()) {
          this->firstNode++;
          numNodesSmaller = graph.smallerNodesConnected[this->firstNode];
        }

      }

      return *this;
    }

    // Equality operator for EdgeIterators.
    bool operator==(const EdgeIterator& other) const {
      bool sameFirstIndex = (this->firstNode == other.firstNode);
      bool sameOtherIndex = (this->adjacencyIndex == other.adjacencyIndex);
      bool sameGraph = (this->graphContainer == other.graphContainer);
      return (sameFirstIndex and sameOtherIndex and sameGraph);
    }

    // Inequality operator for EdgeIterators.
    bool operator!=(const EdgeIterator& other) const {
      return !(*this == other);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE

    const Graph* graphContainer;
    size_type firstNode;
    size_type adjacencyIndex;

    // Default constructor for the EdgeIterator class.
    EdgeIterator(const Graph* graphContainer, size_type firstNode, 
      size_type adjacencyIndex) : graphContainer(graphContainer), 
      firstNode(firstNode), adjacencyIndex(adjacencyIndex) {
    }

    // A function that asserts if the edgeIterator refers to a valid edge.
    void assert_valid() const {

      assert (this->graphContainer != nullptr);
      assert (this->firstNode >= 0);
      assert (this->adjacencyIndex >= 0);

      const Graph& graph = *(this->graphContainer);
      
      assert (this->firstNode <= graph.num_nodes());

      if (this->firstNode != graph.num_nodes()) {
        size_type maxAdj = graph.smallerNodesConnected[this->firstNode];
        assert (this->adjacencyIndex < maxAdj);
      } else {
        assert (this->adjacencyIndex == 0);
      }
      
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  // The constructor for the starting EdgeIterator.
  edge_iterator edge_begin() const {
    
    // Find the first node that is connected to another node with a smaller
    // index.
    size_type firstNode = 0;
    while (this->smallerNodesConnected[firstNode] == 0) {
      firstNode++;
    }

    return EdgeIterator(this, firstNode, 0);
  }
  
  // The constructor for the ending EdgeIterator.
  edge_iterator edge_end() const {
    return EdgeIterator(this, this->num_nodes(), 0);
  }

 private:
  /** A function that prints the data from the incidentArray and the
   *  smallerNodesConnected class variables. Mainly used for personal
   *  debugging on smaller graphs.
   */
  void printIncidentData() const {

    std::cout << "Incident edge information:" << std::endl;
    for (size_type i = 0; i < this->num_nodes(); i++) {
      std::cout << i << ":";
      for (size_type j = 0; j < this->incidentArray[i].size(); j++) {
        std::cout << " " << this->incidentArray[i][j];
      }
      std::cout << std::endl;
    }

    std::cout << "Smaller Incident Nodes Count:" << std::endl;
    for (size_type i = 0; i < this->num_nodes(); i++) {
      std::cout << this->smallerNodesConnected[i] << " ";
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }

};

#endif // CME212_GRAPH_HPP
