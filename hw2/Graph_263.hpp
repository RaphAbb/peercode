#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
//#include <iostream>

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

  /** Note that the data in Points, NValues, EValues, EdgeA_ids, and EdgeB_ids are never deleted.
      The data in Nadj_list and Eadj_list are selectively deleted when nodes and edges are deleted.
      This is to speed up the delete methods. */
  // Node Point data (index is unique node id)
  std::vector<Point> Points;
  // Edge Node id's data (index is unique edge id)
  std::vector<unsigned> EdgeA_ids;
  std::vector<unsigned> EdgeB_ids;
  // adjacency list: unique Node id -> connected unique Node id's
  std::vector<std::vector<unsigned>> Nadj_list;
  // adjacency list: Node id -> corresponding unique Edge id's to the unique Node id's above
  std::vector<std::vector<unsigned>> Eadj_list;
  // storage for node value type (index is unique node id)
  std::vector<V> NValues;
  // storage for edge value type (index is unique edge id)
  std::vector<E> EValues;
  // maps node unique ids to current node index
  // NOTE values at indices corresponding to deleted node unique ids are garbage
  std::vector<unsigned> Nu2i;
  // maps current node index to node unique ids
  std::vector<unsigned> Ni2u;
  // maps edge unique ids to current edge index 
  // NOTE values at indices corresponding to deleted edge unique ids are garbage
  std::vector<unsigned> Eu2i;
  // maps current edge index to edge unique ids
  std::vector<unsigned> Ei2u;
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

  /** Synonym for Node Value type. */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Synonym for Edge Value type. */
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
      // HW0: YOUR CODE HERE
    }

#if 0
    /** Returns constant reference to this node's position. (old function) */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->Points[id_];
    }
#endif

    /** Returns reference to this node's position based on index. */
    Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->Points[graph_->Ni2u[id_]];
    }

    /** Return this node's index, a number in the range [0, Ni2u.size()). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return id_;
    }

    /** Return this node's unique id, a number in the range [0, num_uNodes). */
    size_type uid() const {
      // HW0: YOUR CODE HERE
      return graph_->Ni2u[id_];
    }


    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if(graph_ == n.graph_)
        return id_ == n.index();
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
      if(graph_ != n.graph_)
        return graph_ < n.graph_;
      return id_ < n.index();
    }

    /** returns the reference to the node value */
    node_value_type& value() {
      return graph_->NValues[graph_->Ni2u[id_]];
    }

    /** returns the constant reference to the node value */
    const node_value_type& value() const {
      return graph_->NValues[graph_->Ni2u[id_]];
    }

    /** returns the degree of the node, the number of edges incident to it */
    size_type degree() const {
      return graph_->Eadj_list[graph_->Ni2u[id_]].size();
    }

    /** returns an iterator pointing to the first edge that is incident to this node.
    *  @return A IncidentIterator object
    */
    IncidentIterator edge_begin() const {
      return IncidentIterator(graph_, id_, 0);
    }

    /** returns an iterator pointing to the last edge that is incident to this node.
    *  @return A IncidentIterator object
    */
    IncidentIterator edge_end() const {
      return IncidentIterator(graph_, id_, degree());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // internal Node variables
    Graph* graph_;
    size_type id_;
    /** Private Constructor */
    Node(const Graph* g_in, size_type id_in) : graph_(const_cast<Graph*>(g_in)), id_(id_in) {}
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
    return Ni2u.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Returns reference to a node's position based on its unique id. 
      Can be executed on deleted nodes, but will have no real effect. */
  Point& uid_position(size_type uid) {
    return Points[uid];
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] val The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    // HW0: YOUR CODE HERE
    // updates node information
    Points.push_back(position);
    NValues.push_back(val);

    // increases size of adjacency lists
    std::vector<unsigned> newvec;
    Nadj_list.push_back(newvec);
    Eadj_list.push_back(newvec);

    // updates graph information
    Ni2u.push_back(num_uNodes);
    Nu2i.push_back(Ni2u.size()-1);
    num_uNodes++;
    return Node(this,Ni2u.size()-1);
  }

  /** Remove a node and its associated edges from the graph.
   * @param[in] n The node to be deleted.
   * @return 1 if node exists in graph and was successfully removed
   *         0 if node did not exist in the graph and nothing was done
   * @post the Node @a n and its adjacent Edges are invalidated
   * @post new num_nodes() == old num_nodes() - 1
   * @post new num_edges() == old num_edges() - @a n.degree()
   * @post all existing iterators are now considered invalid
   *
   * Complexity: maximum O(Nnodes+Nedges*degree), usually less
   */
  size_type remove_node(const Node& n) {
    // if node exists
    if(n.index()<size())
    {
      // creates list of unique edge ids to delete, then deletes those edges, max order O(Nedges*degree)
      std::vector<size_type> edgestodelete = Eadj_list[Ni2u[n.index()]];
      for(auto eit = edgestodelete.begin() ; eit != edgestodelete.end() ; ++eit) {
        size_type indx = Eu2i[*eit];
        remove_edge(edge(indx));
      }
      // deletes the node, max order O(Nnodes)
      std::vector<size_type>::iterator it = Ni2u.erase(Ni2u.begin() + n.index());
      // updates the indices in Nu2i, max order O(Nnodes)
      for(; it!=Ni2u.end() ; ++it)
      {
        Nu2i[*it] = Nu2i[*it]-1;
      }
      return 1;
    }
    return 0;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    // also makes sure if Node id is not too large
    if(n.graph_ == this)
      return n.id_ < size();
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
    assert(i < size());
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return this edge's index, a number in the range [0, num_Edges). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return id_;
    }

    /** returns the reference to the edge value */
    edge_value_type& value() {
      return graph_->EValues[graph_->Ei2u[id_]];
    }

    /** returns the constant reference to the edge value */
    const edge_value_type& value() const {
      return graph_->EValues[graph_->Ei2u[id_]];
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      // retrieves unique Node id
      size_type uA_i = graph_->EdgeA_ids[graph_->Ei2u[id_]];
      return graph_->node(graph_->Nu2i[uA_i]);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      // retrieves unique Node id
      size_type uB_i = graph_->EdgeB_ids[graph_->Ei2u[id_]];
      return graph_->node(graph_->Nu2i[uB_i]);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // compares node id's in both directions
      if(graph_ != e.graph_)
        return false;
      if(graph_->EdgeA_ids[graph_->Ei2u[id_]] == graph_->EdgeA_ids[graph_->Ei2u[e.index()]])
        if(graph_->EdgeB_ids[graph_->Ei2u[id_]] == graph_->EdgeB_ids[graph_->Ei2u[e.index()]])
          return true;
      if(graph_->EdgeA_ids[graph_->Ei2u[id_]] == graph_->EdgeB_ids[graph_->Ei2u[e.index()]])
        return (graph_->EdgeB_ids[graph_->Ei2u[id_]] == graph_->EdgeA_ids[graph_->Ei2u[e.index()]]);
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if(graph_ != e.graph_)
        return graph_ < e.graph_;
      return id_ < e.index();
    }

    /** Returns the length of this Edge. */
    double length() const {
      return norm(node1().position()-node2().position());
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // internal Edge variables
    Graph* graph_;
    size_type id_;
    /** Private Constructor */
    Edge(const Graph* g_in, size_type id_in) : graph_(const_cast<Graph*>(g_in)), id_(id_in) {}
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
    return Ei2u.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i<num_edges());
    return Edge(this,i);        // Valid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    if(num_edges() == 0)
      return false;
    // Uses the node adjacency list to determine if edge already exists
    for(size_type i = 0; i < Nadj_list[a.uid()].size(); i++)
    {
      if(Nadj_list[a.uid()][i] == b.uid())
      {
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val = edge_value_type()) {
    // HW0: YOUR CODE HERE
    // Uses the node adjacency list to determine if edge exists
    for(size_type i = 0; i < Nadj_list[a.uid()].size(); i++)
    {
      if(Nadj_list[a.uid()][i] == b.uid())
      {
        // determines the edge id corresponding to that position on the node adjacency list
        size_type uE_ind = Eadj_list[a.uid()][i];
        // returns existing edge
        if(EdgeA_ids[uE_ind] == a.uid())
          if(EdgeB_ids[uE_ind] == b.uid())
          {
            return Edge(this,Eu2i[uE_ind]);
          }
        // switches a,b and returns existing edge with nodes switched
        if(EdgeB_ids[uE_ind] == a.uid())
          if(EdgeA_ids[uE_ind] == b.uid())
          {
            EdgeA_ids[uE_ind] = a.uid();
            EdgeB_ids[uE_ind] = b.uid();
            return Edge(this,Eu2i[uE_ind]);
          } 
      }
    }

    // else, creates and returns new edge
    // updates adjacency lists
    Nadj_list[a.uid()].push_back(b.uid());
    Nadj_list[b.uid()].push_back(a.uid());
    Eadj_list[a.uid()].push_back(num_uEdges);
    Eadj_list[b.uid()].push_back(num_uEdges);
    // adds edge data
    EdgeA_ids.push_back(a.uid());
    EdgeB_ids.push_back(b.uid());
    EValues.push_back(val);
    // adds graph data
    Ei2u.push_back(num_uEdges);
    Eu2i.push_back(Ei2u.size()-1);
    num_uEdges++;
    return Edge(this,Ei2u.size()-1);  
  }

  /** Removes an edge from the graph.
   * @param[in] a,b The Nodes at both ends of the edge to be deleted.
   * @return 1 if edge exists in graph and was successfully removed
   *         0 if edge did not exist in the graph and nothing was done
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @post the Edge @a e is invalidated
   * @post new num_edges() == old num_edges() - 1
   * @post all existing edge iterators are now considered invalid
   *
   * Complexity: maximum O(Nedges), usually less
   */
  size_type remove_edge(const Node& a, const Node& b) {
    // first check validity of nodes in current graph
    if(a.index() >= size() || b.index() >= size())
      return 0;
    // then check if edge exists and return its unique id, searching a small vector, so order <<O(Nedges)
    int Euid = -1;
    for(size_type i = 0; i < Nadj_list[a.uid()].size(); i++)
    {
      if(Nadj_list[a.uid()][i] == b.uid())
      {
        Euid = Eadj_list[a.uid()][i];
      }
    }
    // if edge exists
    if(Euid>=0)
    {
      // finds current index from unique index
      unsigned Eindx = Eu2i[Euid];
      // deletes data from adjacency lists - these are finds and deletes on small vectors, so order <<O(Nedges)
      Eadj_list[a.uid()].erase(std::remove(Eadj_list[a.uid()].begin(), Eadj_list[a.uid()].end(), Euid), Eadj_list[a.uid()].end());
      Eadj_list[b.uid()].erase(std::remove(Eadj_list[b.uid()].begin(), Eadj_list[b.uid()].end(), Euid), Eadj_list[b.uid()].end());
      Nadj_list[a.uid()].erase(std::remove(Nadj_list[a.uid()].begin(), Nadj_list[a.uid()].end(), b.uid()), Nadj_list[a.uid()].end());
      Nadj_list[b.uid()].erase(std::remove(Nadj_list[b.uid()].begin(), Nadj_list[b.uid()].end(), a.uid()), Nadj_list[b.uid()].end());
      // deletes the Edge in Ei2u, max order O(Nedges)
      std::vector<size_type>::iterator it = Ei2u.erase(Ei2u.begin() + Eindx);
      // updates the indices in Eu2i, max order O(Nedges)
      for(; it!=Ei2u.end() ; ++it)
      {
        Eu2i[*it] = Eu2i[*it]-1;
      }
      return 1;
    }
    return 0;
  }

  /** Removes an edge from the graph.
   * @param[in] e The edge to be deleted.
   * @return 1 if edge exists in graph and was successfully removed
   *         0 if edge did not exist in the graph and nothing was done
   * @pre The edge @a e must be a valid edge in this graph. 
   * @post the Edge @a e is invalidated
   * @post new num_edges() == old num_edges() - 1
   * @post all existing edge iterators are now considered invalid
   *
   * Complexity: maximum O(Nedges), usually less
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(),e.node2());
  }



  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    Points.clear();
    EdgeA_ids.clear();
    EdgeB_ids.clear();
    Nadj_list.clear();
    Eadj_list.clear();
    NValues.clear();
    EValues.clear();
    Ni2u.clear();
    Nu2i.clear();
    Ei2u.clear();
    Eu2i.clear();
    num_uNodes = 0;
    num_uEdges = 0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
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

    /** Constructor for a NodeIterator
      @param g_in Pointer to the current graph
      @param Niter_in A number representing the current position of the iterator
    */
    NodeIterator(const Graph* g_in, size_type Niter_in) : graph_(const_cast<Graph*>(g_in)), Niter_(Niter_in) {}

    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Returns the Node that this iterator is pointing to.
      @return A Node object
    */
    Node operator*() const {
      return Node(graph_,Niter_);
    }

    /** Increments the iterator. 
      @return A NodeIterator object
    */
    NodeIterator& operator++() {
      Niter_++;
      return *this;
    }

    /** Checks for equality of this NodeIterator and @a NI_in. */
    bool operator==(const NodeIterator& NI_in) const {
      if(graph_ != NI_in.graph_)
        return false;
      return Niter_==NI_in.Niter_;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type Niter_;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Returns an iterator pointing to the first node in the graph.
  *  @return A NodeIterator object
  */
  NodeIterator node_begin() const {
    return NodeIterator(this,0);
  }

  /** Returns an iterator pointing to the last node in the graph.
  *  @return A NodeIterator object
  */
  NodeIterator node_end() const {
    return NodeIterator(this,size());
  }

  /** Remove a node and its associated edges from the graph using an iterator.
   * @param[in] n_it A NodeIterator pointing to the node to be deleted.
   * @return A NodeIterator pointing to the next node in the graph.
   * @post the Node n=*(@a n_it) and its adjacent Edges are invalidated
   * @post new num_nodes() == old num_nodes() - 1
   * @post new num_edges() == old num_edges() - n.degree()
   * @post all existing iterators are now considered invalid
   *
   * Complexity: maximum O(Nnodes+Nedges*degree), usually less
   */
  NodeIterator remove_node(NodeIterator n_it) {
    size_type indx = (*n_it).index();
    remove_node(node(indx));
    return NodeIterator(this,indx);
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

    /** Constructor for a IncidentIterator
      @param g_in Pointer to the current graph
      @param Nid_in The node ID of the current node
      @param Niter_in A number representing the current position of the iterator
    */
    IncidentIterator(const Graph* g_in, size_type Nid_in, size_type Niter_in) :
           graph_(const_cast<Graph*>(g_in)), Nid_(Nid_in), Niter_(Niter_in) {}

    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** Returns the Edge that this iterator is pointing to.
      @return An Edge object
    */
    Edge operator*() const {
      size_type uedge_id = graph_->Eadj_list[graph_->Ni2u[Nid_]][Niter_];
      // checks if node1() of Edge is spawning node
      // if not, switches orientation of edge
      if(graph_->EdgeA_ids[uedge_id] != graph_->Ni2u[Nid_])
      {
        graph_->EdgeB_ids[uedge_id] = graph_->EdgeA_ids[uedge_id];
        graph_->EdgeA_ids[uedge_id] = graph_->Ni2u[Nid_];
      }
      return Edge(graph_,graph_->Eu2i[uedge_id]);
    }

    /** Increments the iterator. 
      @return An IncidentIterator object
    */
    IncidentIterator& operator++() {
      Niter_++;
      return *this;
    }

    /** Checks for equality of this IncidentIterator and @a II_in. */
    bool operator==(const IncidentIterator& II_in) const {
      if(graph_ != II_in.graph_)
        return false;
      if(Nid_ != II_in.Nid_)
        return false;
      return Niter_==II_in.Niter_;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type Nid_;
    size_type Niter_;
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
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

    /** Constructor for an EdgeIterator
      @param g_in Pointer to the current graph
      @param Niter_in A number representing the current position of the iterator
    */
    EdgeIterator(const Graph* g_in, size_type Niter_in) : graph_(const_cast<Graph*>(g_in)), Niter_(Niter_in) {}

    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

     /** Returns the Edge that this iterator is pointing to.
      @return An Edge object
    */
    Edge operator*() const {
      return Edge(graph_,Niter_);
    }

    /** Increments the iterator. 
      @return An EdgeIterator object
    */
    EdgeIterator& operator++() {
      Niter_++;
      return *this;
    }

    /** Checks for equality of this EdgeIterator and @a EI_in. */
    bool operator==(const EdgeIterator& EI_in) const {
      if(graph_ != EI_in.graph_)
        return false;
      return Niter_==EI_in.Niter_;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type Niter_;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Returns an iterator pointing to the first edge in the graph.
  *  @return An EdgeIterator object
  */
  EdgeIterator edge_begin() const {
    return EdgeIterator(this,0);
  }

  /** Returns an iterator pointing to the last edge in the graph.
  *  @return An EdgeIterator object
  */
  EdgeIterator edge_end() const {
    return EdgeIterator(this,num_edges());
  }

  /** Removes an edge from the graph using an iterator.
   * @param[in] e_it An EdgeIterator pointing to the edge to be deleted.
   * @return An EdgeIterator pointing to the next edge in the graph.
   * @post the Edge *(@a e_it) is invalidated
   * @post new num_edges() == old num_edges() - 1
   * @post all other existing edge iterators are now considered invalid
   *
   * Complexity: maximum O(Nedges), usually less
   */
  EdgeIterator remove_edge(EdgeIterator e_it) {
    size_type indx = (*e_it).index();
    remove_edge(edge(indx));
    return EdgeIterator(this,indx);
  }

 private:

  // HW0: YOUR CODE HERE
  size_type num_uNodes = 0;
  size_type num_uEdges = 0;
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
