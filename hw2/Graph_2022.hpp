#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <set>
#include <map>
#include <iostream>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph
{
private:
  /** Predeclaration of InternalEdge, InternalNode and EdgeMapper classes. */
  class InternalNode;
  class InternalEdge;
  class EdgeMapper;

public:
  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V, E>;
  using node_value_type = V;
  using edge_value_type = E;

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
  Graph() : edgeCounter(0), nodeCounter(0), node_info(std::vector<InternalNode>()), edge_info(std::vector<InternalEdge>())
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
  class Node : private totally_ordered<Node>
  {

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

    /** Node class that constructs an invalid node. */
    Node() : graph_(NULL), u_id_node_(-1)
    {
    }

    /** Return this node's value, which is user specified, e.g. shortest distance 
     * to a node. Calling this function allows the modification of the node 
     * value because of the reference type.
     * 
     * Complexity: O(1)
     * 
     */

    node_value_type &value()
    {
      return graph_->node_info[index()].node_value_;
    }
    /** Return this node's value, which is user specified, e.g. shortest distance 
     * to a node. Calling this function does not allow the modification of the node 
     * value because of the use of keyword const. 
     * 
     * Complexity: O(1)
     * 
     */

    const node_value_type &value() const
    {
      return graph_->node_info[index()].node_value_;
    }

    /** Return the number of edges incident to this node 
     * 
     * Complexity: O(1)
     * 
     */

    size_type degree() const
    {
      size_type degree = graph_->adjacency_map[u_id_node_].size();
      return degree;
    }

    /** Return a pointer to the beginning of the adjacency list of this node
     * which is a map of the adjacent node's index to the index of the edge
     * between those two nodes
     * 
     * @return an IncidentIterator object with a pointer intialised to the start of 
     * adjacency_map[u_id_node_]
     * 
     * Complexity: O(1)
     * 
     */
    IncidentIterator edge_begin() const
    {
      return IncidentIterator(graph_, graph_->adjacency_map.at(u_id_node_).begin(), u_id_node_);
    }

    /** Return a pointer to one past the end of the adjacency list of this node
     * which is a map of the adjacent node's index to the index of the edge
     * between those two nodes
     * 
     * @return an IncidentIterator object with a pointer intialised to one past the end of
     * adjacency_map[u_id_node_]
     * 
     * Complexity: O(1)
     * 
     */
    IncidentIterator edge_end() const
    {
      return IncidentIterator(graph_, graph_->adjacency_map.at(u_id_node_).end(), u_id_node_);
    }

    /** Return the Point object (coordinates) of this node 
     * 
     * Complexity: O(1)
     * 
     */

    const Point &position() const
    {
      return graph_->node_info[index()].point_;
    }

    /** Return the Point object (coordinates) of this node as a reference
     * 
     * Complexity: O(1)
     * 
     */
    Point &position()
    {
      return graph_->node_info[index()].point_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const
    {
      return graph_->Node_ID_to_index_[u_id_node_];
    }
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node &n) const
    {
      return (graph_ == n.graph_ && index() == n.index());
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node &n) const
    {
      return ((graph_ == n.graph_ && u_id_node_ < n.u_id_node_) || (graph_ != n.graph_ && u_id_node_ < n.u_id_node_));
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph *graph_;
    size_type u_id_node_;
    // Node constructor taking a pointer to the graph and an index
    Node(const Graph *graph, size_type u_id_node) : graph_(const_cast<Graph *>(graph)), u_id_node_(u_id_node)
    {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const
  {
    return node_info.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const
  {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position and value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   * @post result_node.value() == node_value_type()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point &position, const node_value_type &value = node_value_type())
  {
    Node_ID_to_index_.push_back(size());
    InternalNode internalNode(nodeCounter, position, value);
    node_info.push_back(internalNode);
    index_to_Node_ID.push_back(nodeCounter);
    std::map<size_type, size_type> innerMap;
    adjacency_map.insert(std::pair<size_type, std::map<size_type, size_type>>(nodeCounter, innerMap));
    nodeCounter += 1;
    return Node(this, nodeCounter - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node &n) const
  {
    return (this == n.graph_ && Node_ID_to_index_[n.u_id_node_] != -1);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const
  {
    return Node(this, index_to_Node_ID[i]);
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
  class Edge : private totally_ordered<Edge>
  {
  public:
    /** Construct an invalid Edge. */
    Edge() : graph_(NULL), u_id_edge_(-1)
    {
    }

    /** Return this edge's value, which is user specified, e.g. spring constant.
     * Calling this function allows the modification of the edge 
     * value because of the reference type.
     * 
     * Complexity: O(1)
     * 
     */

    edge_value_type &value()
    {
      size_type currentIndex = graph_->Edge_ID_to_index_[u_id_edge_].edge_index_;
      return graph_->edge_info[currentIndex].spring_constant_;
    }

    /** Return this edge's value, which is user specified, e.g. spring constant
     * Calling this function does not allow the modification of the edge 
     * value because of the use of keyword const. 
     * 
     * Complexity: O(1)
     * 
     */

    const edge_value_type &value() const
    {
      size_type currentIndex = graph_->Edge_ID_to_index_[u_id_edge_].edge_index_;
      return graph_->edge_info[currentIndex].spring_constant_;
    }
    /** Return a node of this Edge */
    Node node1() const
    {
      return graph_->node(node1_ID_);
    }

    /** Return the other node of this Edge */
    Node node2() const
    {
      return graph_->node(node2_ID_);
    }

    /** Returns the length between two nodes */
    double length() const
    {
      Point point1 = node1().position();
      Point point2 = node2().position();
      double dist = Helper::EuclideanDistance(point1, point2);
      return dist;
    }

    /** Returns the initial between two nodes (the initial length of an edge before potential
     * movement of nodes) */
    double get_initial_length() const
    {
      size_type currentIndex = graph_->Edge_ID_to_index_[u_id_edge_].edge_index_;
      return graph_->edge_info[currentIndex].initial_length_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */

    bool operator==(const Edge &e) const
    {
      return ((e.graph_ == graph_ && ((e.node1_ID_ == node1_ID_ && e.node2_ID_ == node2_ID_) || (e.node1_ID_ == node2_ID_ && e.node2_ID_ == node1_ID_))));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge &e) const
    {
      return ((e.graph_ == graph_ && u_id_edge_ < e.u_id_edge_) || (e.graph_ != graph_ && u_id_edge_ <= e.u_id_edge_));
    }

  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph *graph_;
    size_type u_id_edge_;
    size_type node1_ID_;
    size_type node2_ID_;
    // Edge constructor taking a pointer to the graph,an index, and the two node indices
    // of the nodes connected by this edge
    Edge(const Graph *graph, size_type u_id_edge, size_type node1_ID, size_type node2_ID) : graph_(const_cast<Graph *>(graph)), u_id_edge_(u_id_edge), node1_ID_(node1_ID), node2_ID_(node2_ID)
    {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const
  {
    return edge_info.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   * std::vector.at is O(1)
   */
  Edge edge(size_type i) const
  {
    if (i < num_edges())
    {
      InternalEdge internalEdge = edge_info.at(i);
      return Edge(this, internalEdge.u_id_edge_, internalEdge.node1_ID_, internalEdge.node2_ID_);
    }
    throw "Invalid index";
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node &a, const Node &b) const
  {
    assert(a.graph_ == this && b.graph_ == this);
    if ((has_node(a) && has_node(b)))
    {
      if (edge_info.size() == 0)
      {
        return false;
      }
      for (unsigned i = 0; i < Edge_ID_to_index_.size(); i++)
      {
        size_type node1_ID = Edge_ID_to_index_[i].node1_ID_;
        size_type node2_ID = Edge_ID_to_index_[i].node2_ID_;
        if ((a.u_id_node_ == node1_ID && b.u_id_node_ == node2_ID) || (a.u_id_node_ == node2_ID && b.u_id_node_ == node1_ID))
        {
          if (Edge_ID_to_index_[i].edge_index_ == -1)
          {
            return false;
          }
          else
          {
            return true;
          }
        }
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
  Edge add_edge(const Node &a, const Node &b)
  {
    assert(a.graph_ == this && b.graph_ == this);
    assert(has_node(a) && has_node(b));
    size_type smaller_ID = std::min(a.u_id_node_, b.u_id_node_);
    size_type larger_ID = std::max(a.u_id_node_, b.u_id_node_);
    if (has_edge(a, b))
    {
      std::map<size_type, size_type> firstNodePair = adjacency_map[smaller_ID];
      return Edge(this, firstNodePair.at(larger_ID), a.u_id_node_, b.u_id_node_);
    }
    adjacency_map[smaller_ID][larger_ID] = edgeCounter;
    adjacency_map[larger_ID][smaller_ID] = edgeCounter;
    edge_value_type length = Helper::EuclideanDistance(a.position(), b.position());
    InternalEdge internalEdge(edgeCounter, a.u_id_node_, b.u_id_node_, length);
    edge_info.push_back(internalEdge);
    checkEdgeMapper(a, b);
    edgeCounter += 1;
    return Edge(this, (edgeCounter - 1), a.u_id_node_, b.u_id_node_);
  }

  /*Checking if an edge has been added and erased before in order
  to keep track of the correct indices and refrain from adding 
  an edge twice to Edge_ID_to_index_ */
  void checkEdgeMapper(const Node &a, const Node &b)
  {
    for (unsigned i = 0; i < Edge_ID_to_index_.size(); i++)
    {
      size_type node1_ID = Edge_ID_to_index_[i].node1_ID_;
      size_type node2_ID = Edge_ID_to_index_[i].node2_ID_;
      if ((a.u_id_node_ == node1_ID && b.u_id_node_ == node2_ID) || (a.u_id_node_ == node2_ID && b.u_id_node_ == node1_ID))
      {
        if (Edge_ID_to_index_[i].edge_index_ == -1)
        {
          Edge_ID_to_index_[i].edge_index_ = num_edges();
          break;
        }
      }
    }
    EdgeMapper edgeMapper = EdgeMapper(num_edges() - 1, a.u_id_node_, b.u_id_node_);
    Edge_ID_to_index_.push_back(edgeMapper);
    
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear()
  {
    node_info.clear();
    Node_ID_to_index_.clear();
    index_to_Node_ID.clear();
    Edge_ID_to_index_.clear();
    adjacency_map.clear();
    edge_info.clear();
    edgeCounter = 0;
    nodeCounter = 0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */

  class NodeIterator : private totally_ordered<NodeIterator>
  {

  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Node;
    // Element type

    using pointer = Node *;
    // Pointers to elements

    using reference = Node &;
    // Reference to elements

    using difference_type = std::ptrdiff_t;
    // Signed difference

    using iterator_category = std::input_iterator_tag;
    // Weak Category, Proxy

    /** Construct an invalid NodeIterator */
    NodeIterator() : graph_(NULL), ptr_(NULL)
    {
    }
    /** Construct a valid NodeIterator, containing a pointer to the graph
     * and a pointer to the vector of Internal Node objects.
     */
    NodeIterator(const Graph *graph, typename std::vector<InternalNode>::const_iterator ptr) : graph_(graph), ptr_(ptr)
    {
    }

    /** Dereference a pointer pointing to an Internal Node object.
   * @pre ptr_ is pointing to an Internal Node object
   * @return Node object corresponding to the data in the Internal Node
   *
   * Complexity: O(1)
   */
    Node operator*() const
    {
      return Node(graph_, (*ptr_).u_id_node_);
    }

    /** Increment a pointer pointing to an Internal Node object.
   * @pre ptr_ is pointing to an Internal Node object
   * @post ptr_ is pointing to the next Internal Node object in the vector
   * of internal nodes
   * @return a NodeIterator object with the updated ptr_
   * 
   * Complexity: O(1)
   */

    NodeIterator operator++()
    {
      ptr_++;
      return *this;
    }

    /** Test whether two NodeIterators are equal.
   * @return True if the two pointers point to the same Internal Node object
   * and if the two graph pointers are the same
   *
   * Complexity: O(1)
   */

    bool operator==(const NodeIterator &n) const
    {
      return (n.graph_ == graph_ && n.ptr_ == ptr_);
    }

  private:
    friend class Graph;
    const Graph *graph_;
    typename std::vector<InternalNode>::const_iterator ptr_;
  };

  /** Return a pointer to the first element of the vector of Internal Nodes
  * @return a NodeIterator object with a pointer intialised to the start of
  * node_info.
  * 
  * Complexity: O(1)
  * 
  */
  NodeIterator node_begin() const
  {
    return NodeIterator(this, node_info.begin());
  }

  /** Return a pointer to one past the end of the vector of Internal Nodes
   * @return a NodeIterator object with a pointer intialised to one past the end of
   * node_info. 
   *     
   * Complexity: O(1)
   * 
   */

  NodeIterator node_end() const
  {
    return NodeIterator(this, node_info.end());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>
  {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() : graph_(NULL), ptr_(std::map<size_type, size_type>::const_iterator())
    {
    }

    /** Construct a valid IncidentIterator, containing a pointer to the graph, a
     * pointer to a map of size_type to size_type elements, and the node ID 
     * corresponding to the node whose adjacency list we are exploring.
     */
    IncidentIterator(const Graph *graph, typename std::map<size_type, size_type>::const_iterator ptr, size_type u_id_node) : graph_(graph), ptr_(ptr), u_id_node_(u_id_node)
    {
    }

    /** Dereference a pointer pointing to a map of size_type to size_type.
   * @pre ptr_ is pointing to an std::pair object (adjacent nodeID, edge index)
   * @return Edge object corresponding to the data in the pair in the map
   *
   * Complexity: O(1)
   */

    Edge operator*() const
    {
      return Edge(graph_, ptr_->second, u_id_node_, ptr_->first);
    }

    /** Increment a pointer pointing to a map object.
   * @pre ptr_ is pointing to an std::pair object (adjacent nodeID, edge index)
   * @post ptr_ is pointing to the next std::pair object (adjacent nodeID, edge index)
   * @return an IncidentIterator object with the updated ptr_
   * 
   * Complexity: O(1)
   */

    IncidentIterator &operator++()
    {
      ptr_++;
      return *this;
    }

    /** Test whether two IncidentIterators are equal.
   * @return True if the two pointers point to the same std::pair object
   * and if the two graph pointers are the same, and if the two node IDs are 
   * the same.
   *
   * Complexity: O(1)
   */

    bool operator==(const IncidentIterator &i) const
    {
      return (i.graph_ == graph_ && i.ptr_ == ptr_ && i.u_id_node_ == u_id_node_);
    }

  private:
    friend class Graph;
    const Graph *graph_;
    typename std::map<size_type, size_type>::const_iterator ptr_;
    size_type u_id_node_;
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>
  {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() : graph_(NULL), ptr_(NULL)
    {
    }
    /** Construct a valid EdgeIterator, containing a pointer to the graph
     * and a pointer to the vector of Internal Edge objects.
     */
    EdgeIterator(const Graph *graph, typename std::vector<InternalEdge>::const_iterator ptr) : graph_(graph), ptr_(ptr)
    {
    }

    /** Dereference a pointer pointing to an Internal Edge object.
   * @pre ptr_ is pointing to an Internal Edge object
   * @return Edge object corresponding to the data in the Internal Edge
   *
   * Complexity: O(1)
   */
    Edge operator*() const
    {
      return Edge(graph_, (*ptr_).u_id_edge_, (*ptr_).node1_ID_, (*ptr_).node2_ID_);
    }

    /** Increment a pointer pointing to an Internal Edge object.
   * @pre ptr_ is pointing to an Internal Edge object
   * @post ptr_ is pointing to the next Internal Edge object in the vector
   * of internal edges
   * @return an EdgeIterator object with the updated ptr_
   * 
   * Complexity: O(1)
   */

    EdgeIterator &operator++()
    {
      ptr_++;
      return *this;
    }

    /** Test whether two EdgeIterators are equal.
   * @return True if the two pointers point to the same Internal Edge object
   * and if the two graph pointers are the same
   *
   * Complexity: O(1)
   */

    bool operator==(const EdgeIterator &e) const
    {
      return (e.graph_ == graph_ && e.ptr_ == ptr_);
    }

  private:
    friend class Graph;
    const Graph *graph_;
    typename std::vector<InternalEdge>::const_iterator ptr_;
  };

  /** Return a pointer to the first element of the vector of Internal Edges
  * @return an EdgeIterator object with a pointer intialised to the start of
  * edge_info.
  * 
  * Complexity: O(1)
  * 
  */

  EdgeIterator edge_begin() const
  {
    return EdgeIterator(this, edge_info.begin());
  }

  /** Return a pointer to one past the end of the vector of Internal Edges
   * @return an EdgeIterator object with a pointer intialised to one past the end of
   * edge_info.
   * 
   * Complexity: O(1)
   * 
   */

  EdgeIterator edge_end() const
  {
    return EdgeIterator(this, edge_info.end());
  }

  /** Remove a node from the graph.
   * @pre @a nodes is a valid node of this graph
   * @return an int (1) indicating the successful removal of a node
   * @post has_node(@a node) == false
   * @post If old has_node(@a node) == false, new num_nodes() == old num_nodes().
   *       Else,                              new num_nodes() == old num_nodes() - 1.
   * @post adjacency_map.size() decreases by the degree of the node * 2
   * @post new num_edges() is now decreased by the degree of the node
   * Complexity: At most O(num nodes) since the graph is sparse and the degree of the node is
   * far less than the number of nodes.
   */

  size_type remove_node(const Node &n)
  {
    int currentNodeIndex = n.index();
    size_type currentNodeID = n.u_id_node_;
    InternalNode lastNode = node_info[node_info.size() - 1];
    size_type lastNodeID = lastNode.u_id_node_;
    if (currentNodeIndex != -1)
    {
      Node_ID_to_index_[currentNodeID] = -1;
      if (currentNodeIndex != (int)node_info.size() - 1)
      {
        node_info[currentNodeIndex] = lastNode;
        Node_ID_to_index_[lastNodeID] = currentNodeIndex;
        index_to_Node_ID[currentNodeIndex] = lastNodeID;
      }
      node_info.pop_back();
      index_to_Node_ID.pop_back();
      for (IncidentIterator it = n.edge_begin(); it != n.edge_end(); ++it)
      {
        size_type adjNodeID = (*it).node2_ID_;
        int adjNodeIndex = Node_ID_to_index_[adjNodeID];
        if (adjNodeIndex > -1 && adjNodeIndex < (int)size())
        {
          Node adjNode = node(adjNodeIndex);
          remove_edge(adjNode, n);
        }
      }
    }
    return 1;
  }

  /** Remove a node from the graph.
   * @pre @a nodes is a valid node of this graph
   * @return a node iterator pointing to the node currently at index i in @a node_info
   * @post has_node(@a node) == false
   * @post If old has_node(@a node) == false, new num_nodes() == old num_node().
   *       Else,                              new num_nodes() == old num_nodes() - 1.
   * @post If old has_node(@a node) == false, new index_to_Node_ID.size() = old index_to_Node_ID.size()
   *       Else,                              new index_to_Node_ID.size() = old index_to_Node_ID.size() - 1
   * @post adjacency_map.size() decreases by the degree of the node * 2
   * * @post new num_edges() is now decreased by the degree of the node
   * Complexity: At most O(num nodes) since the graph is sparse and the degree of the node is
   * far less than the number of nodes.
   */
  node_iterator remove_node(node_iterator n_it)
  {
    Node node = *n_it;
    remove_node(node);
    return n_it;
  }

  /** Remove an edge from the graph 
   * @pre @a a and @a b are distinct valid nodes of this graph 
   * @return an int (1) indicating the successful removal of an edge
   * @post has_edge(@a a, @a b) == false
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() - 1.
   *
   * @post new adjacency_map.size() = old adjacency_map.size() - 2
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type remove_edge(const Node &node1, const Node &node2)
  {
    if (edge_info.size() == 0)
    {
      return 0;
    }
    size_type currentEdgeID = adjacency_map[node1.u_id_node_][node2.u_id_node_];
    int currentEdgeIndex = Edge_ID_to_index_[currentEdgeID].edge_index_;
    InternalEdge lastEdge = edge_info[edge_info.size() - 1];
    size_type lastEdgeID = lastEdge.u_id_edge_;
    if (currentEdgeIndex != -1)
    {
      Edge_ID_to_index_[currentEdgeID].edge_index_ = -1;
      //make sure we are not swapping the last edge with itself
      if (currentEdgeIndex != (int)edge_info.size() - 1)
      {
        edge_info[currentEdgeIndex] = lastEdge;
        Edge_ID_to_index_[lastEdgeID].edge_index_ = currentEdgeIndex;
      }
      std::map<size_type, size_type> innerMap1 = adjacency_map.at(node1.u_id_node_);
      innerMap1.erase(node2.u_id_node_);
      std::map<size_type, size_type> innerMap2 = adjacency_map.at(node2.u_id_node_);
      innerMap2.erase(node1.u_id_node_);
      edge_info.pop_back();
      return 1;
    }
    return 0;
  }

  /** Remove an edge from the graph 
   * @pre @a edge is a valid edge of this graph 
   * @return an int (1) indicating the successful removal of an edge
   * @post has_edge(@a a, @a b) == false
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() - 1.
   *
   * @post new adjacency_map.size() = old adjacency_map.size() - 2
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */

  size_type remove_edge(const Edge &edge)
  {
    Node node1 = edge.node1();
    Node node2 = edge.node2();
    return remove_edge(node1, node2);
  }

  /** Remove an edge from the graph 
   * @pre @a e_it is a valid edge iterator of this graph 
   * @return an edge iterator pointing to the edge currently at index i in @a edge_info
   * @post has_edge(@a a, @a b) == false
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() - 1.
   *
   * @post new adjacency_map.size() = old adjacency_map.size() - 2
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */

  edge_iterator remove_edge(edge_iterator e_it)
  {
    Edge edge = *e_it;
    Node node1 = edge.node1();
    Node node2 = edge.node2();
    remove_edge(node1, node2);
    return e_it;
  }

private:
  /** Internal Node class containing private attributes
   * of each node. These are the Point object, the node ID, and the node value.
   */
  class InternalNode
  {
  private:
    InternalNode(size_type u_id_node, Point point, node_value_type node_value) : u_id_node_(u_id_node), point_(point), node_value_(node_value)
    {
    }
    friend class Graph;
    size_type u_id_node_;
    Point point_;
    node_value_type node_value_;
  };

  /** Internal Edge class containing private attributes
  * of each edge. These are the edge ID, and the two IDs of the nodes
  * connected by this edge, the spring constant for that edge, and the intial length.
  */
  class InternalEdge
  {
  private:
    InternalEdge(size_type u_id_edge, size_type node1_ID, size_type node2_ID, edge_value_type initial_length) : u_id_edge_(u_id_edge), node1_ID_(node1_ID), node2_ID_(node2_ID), initial_length_(initial_length)
    {
    }
    friend class Graph;
    size_type u_id_edge_;
    size_type node1_ID_;
    size_type node2_ID_;
    edge_value_type spring_constant_;
    edge_value_type initial_length_;
  };

  /** EdgeMapper class containing private attributes of edges to keep track of 
   * deleted edges. These are the edge index (the position of the internal Edge object
   * in edge_info), and the node IDs of the two nodes that are or were connected by that edge.
   */
  class EdgeMapper
  {
  private:
    EdgeMapper(int edge_index, size_type node1_ID, size_type node2_ID) : edge_index_(edge_index), node1_ID_(node1_ID), node2_ID_(node2_ID)
    {
    }
    friend class Graph;
    int edge_index_;
    size_type node1_ID_;
    size_type node2_ID_;
  };

  /** Private attributes of the Graph class. These are an edgeCounter and nodeCounter
   * that increment with each edge/node added and act as unique identifiers,
   * a vector of Interal Node objects, a vector of Internal Edge objects, a vector
   * mapping Node IDs to indices, a vector mapping indices to Node IDs, and a 
   * vector of EdgeMapper object mapping Edge ID to corresponding edge information
   * a map that contains the adjacency list for each node, whereby each node
   * maps to a map of the adjacent nodes to the corresponding edge ID. 
  */

  size_type edgeCounter;
  size_type nodeCounter;

  std::vector<InternalNode> node_info;

  std::vector<InternalEdge> edge_info;
  std::map<size_type, std::map<size_type, size_type>> adjacency_map;

  std::vector<int> index_to_Node_ID;
  std::vector<int> Node_ID_to_index_;
  std::vector<EdgeMapper> Edge_ID_to_index_;
};
#endif // CME212_GRAPH_HPP
