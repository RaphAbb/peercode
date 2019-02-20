#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type containing Node and Edge objects.
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <utility>

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

  /** Predeclaration of the internal_node class and internal edge class. */
  struct internal_node;
  struct internal_edge;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Allows the Node class to support user-specified values */
  using node_value_type = V;

  /** Allows the Edge class to support user-specified values */
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
      Return type of Graph::Node::index(), Graph::Node::degree(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) and
      Graph::edge(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Default destructor 
   * @post Deletion of the object's private variables. 
   */
  ~Graph() = default;

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Lightweight proxy class to the graph's internal_node class.
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
     * @return An invalid node whose member variables are default initialized.
     */
    // Default public constructor. Creates an invalid Node.
    Node() {
    }

    /** Return this node's position. 
     * @return A point object containing the node's position in 3D space. Position
     *            is a const object in this method, so we can inspect the position
     *            that is returned but cannot modify it.
     * 
     * Complexity: O(1).
     */
    const Point& position() const {
      // dereference the graph pointer to access the vector of internal nodes. 
      // Then get the element with this Node's uid and return the point/position
      return graph_->nodes_vector[uid_].point;
    }

    /** Return this node's position.
     * @return A point object containing the node's position in 3D space. Position
     *            is a mutable object in this method, allowing us to modify it.
     * 
     * Complexity: O(1).
     */
    Point& position() {
      // dereference the graph pointer to access the vector of internal nodes. 
      // Then get the element with this Node's uid and return the point/position
      return graph_->nodes_vector[uid_].point;
    }

    /** Return this node's index, a number in the range [0, num_nodes()). 
     * @return An integer specifying the node's index, which can range from 0 
     *    to the total number of nodes - 1. 
     * 
     * Complexity: average case O(1).
     */
    size_type index() const {
      // return the node's index in the vector of valid nodes
        // (note: uid != index after nodes have been removed)
      return graph_->uid2id[uid_]; 
    }

    /** Return this node's associated value. 
     * @return The node's value whose type is user-specified. Value is a mutable 
     *    object in this method, allowing us to modify it. Value can be the 
     *    default value, node_value_type(), if no value was passed to `add_node`
     *    during the node's creation. 
     * 
     * Complexity: O(1).
     */
    node_value_type& value() {
      // dereference the graph pointer to access the vector of internal nodes.
      // Then get the element with this Node's index and return its value
      return graph_->nodes_vector[uid_].value;
    }

    /** Return this node's associated value. 
     * @return The node's value whose type is user-specified. Value is a const 
     *    object in this method, so we can inspect the object that is returned 
     *    but cannot modify it. Value can be the default value, node_value_type(), 
     *    if no value was passed to `add_node` during the node's creation.
     * 
     * Complexity: O(1).
     */
    const node_value_type& value() const {
        // dereference the graph pointer to access the vector of internal nodes.
        // Then get the element with this Node's index and return its value
        return graph_->nodes_vector[uid_].value;
    }

    /** Return the number of edges incident to this node. 
     * @return An integer specifying the number of edges incident to this node. 
     * 
     * Complexity: Average case O(1).
     */
    size_type degree() const {
      // dereference the graph pointer to access the unordered map where keys
        // are node 1's uid. Get the number of values (in this case incident edges)
        // associated with this node's uid.
      return graph_->node_uid_umap.at(uid_).size();
    }

    /** Return an iterator to the first incident edge.
     * @return An incident iterator pointing to the first edge that is incident
     *    to this node (as recorded in the unordered map where the keys are node
     *    1's uid.)
     *
     * Complexity: Average case O(1).
     */
    incident_iterator edge_begin() const {
      // create an incident iterator that points to the beginning of the nested 
        // unordered map, which contains this node's incident edges 
      return IncidentIterator(graph_, uid_, graph_->node_uid_umap.at(uid_).begin());
    }

    /** Return an iterator to the last incident edge.
     * @return An incident iterator pointing to the last edge that is incident
     *    to this node) as recorded in the unordered map where the keys are node
     *    1's uid.)
     *
     * Complexity: Average case O(1).
     */
    incident_iterator edge_end() const {
      // create an incident iterator that points to the end of the nested unordered
        // map, which contains this node's incident edges
      return IncidentIterator(graph_, uid_, graph_->node_uid_umap.at(uid_).end());
    }

    /** Test whether this node and @a n are equal.
     * @param n Node object being compared to the current node.
     * @return True if @a n and the current node have the same graph and the same 
     *    unique id, and false if either the graph or the index are different.
     * 
     * Complexity: O(1).
     */
    bool operator==(const Node& n) const {
      // check that the 2 nodes have pointers to the same graph and the same uid
      if (graph_ == n.graph_ && uid_ == n.uid_) {
        return true;
      }
      // if the two nodes are from different graphs or have different uid's
      else {
        return false;
      }
    }

    /** Test whether this node is less than @a n in a global order.
     * This ordering function is useful for STL containers such as std::map<>.
     * It need not have any geometric meaning.
     * @param n Node object being compared to the current node.
     * @return True if @a n and the current node have the same graph and the current
     *    node's unique id < @a n's uid. False if @a n's graph is different from the
     *    current node's graph, the current node's uid == @a n's uid, or the
     *    current node's uid > @a n's uid.
     * 
     * Complexity: O(1).
     */
    bool operator<(const Node& n) const {
      // check that the 2 nodes have pointers to the same graph and that this
        // node's uid is smaller than the uid of the specified node
      if (graph_ < n.graph_) {
        return true;
      }
      else if (graph_ == n.graph_ && uid_ < n.uid_) {
        return true;
      }
      else {
        return false;
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Pointer back to the Graph container
    Graph* graph_;
    // The node's unique identification number
    size_type uid_;
    // Private Constructor
    Node(const Graph* graph, size_type uid) 
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {}
  };

  /** Return the number of nodes in the graph.
   * @return An integer specifying this graph's total number of nodes. 
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // get the current length of the vector containing internal nodes
    return id2uid.size();
  }

  /** Synonym for size(). 
   * @return An integer specifying this graph's total number of nodes.
   * 
   * Complexity: O(1).
   */
  size_type num_nodes() const {
    return size();
  }

  /** Determine if a Node belongs to this Graph
   * @param n A Node object
   * @return True if @a n is currently a Node of this Graph. False if it is not.
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // check if the specified Node object has a pointer to this graph and if its
      // unique id exists in this graph
    if ((n.graph_ == this) && (uid2id.find(n.uid_) != uid2id.end())) {
      return true;
    }
    // if the Node object has a pointer to another graph or its unique id doesn't
      // currently exist in the graph
    else {
      return false;
    }
  }

  /** Return the node with index @a i.
   * @param i An integer specifying the index of the node to be returned.
   * @return The Node with index @a i.
   * 
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // check that the specified index exists in this graph
    assert(i < num_nodes());
    // get the uid associated with the node at the specified index
    size_type uid = id2uid[i]; 
    // return the node with the relevant uid
    return Node(this, uid);
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value The new node's value
   * @return A Node object for the new node
   *
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   * @post (result_node.value() == @a value) or (result_node.value() == the default
   *          value node_value_type() if no @a value was passed in.)
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, 
      const node_value_type& value = node_value_type()) {
    // get a unique id for a new node
    size_type uid = nodes_vector.size();
    // create a new instance of internal node using the position specified
      // and add it to our internal node vector
    nodes_vector.push_back(internal_node(position, uid, value));
    
    // count how many valid nodes are currently in the graph
    size_type idx = num_nodes();
    // add the node's uid to the vector of valid node uid's 
    id2uid.push_back(uid);    
    // add a mapping between the node's uid and its index in the vector of valid nodes
    uid2id[uid] = idx; 
    // add the node uid to our unordered map of edges (keys in this map are one
      // of the node's uid's)
    node_uid_umap[uid];
    // return a proxy Node that points to the newly created internal node
    return Node(this, uid);
  }

  /** Remove node from a graph, returning 1 if successful and 0 if not.
   * If the specified node is found in the graph, this method invalidates that one
   * node along with all of its incident edges.
   * @param[in] n The Node object that we want to remove from the graph
   * @return The integer 1 if we were able to remove the specified node and its
   *          adjacent edges, and a 0 if we were not able to find the node.
   *
   * @pre @a n is a valid node from an existing graph.
   * @pre num_nodes() >= 1
   * @post If old has_node(@a n) == true, new num_nodes() = old num_nodes() - 1.
   *        Else, new num_nodes() = old num_nodes().
   * @post new num_edges <= old num num_edges()
   * @post has_node (@a n) == false
   *
   * Complexity: Average case O(1). Worst case: O(n).
   */
  size_type remove_node(const Node& n) {
    if (has_node(n)) {
      // get the index of the node that is going to removed
      size_type idx = uid2id.at(n.uid_);
      // overwrite the uid at that index with the last uid currently in the vector
      size_type last_uid = *(id2uid.end() - 1);
      id2uid[idx] = last_uid;
      // remove the last element of the vector since it's now been duplicated
      id2uid.pop_back();

      // remove the overwritten node from the umap linking node id with node uid
      uid2id.erase(n.uid_);
      // in the same unordered map, if the overwritten node was not the last node (i.e.,
        // the one we just removed), update the formerly last node with its new index
      if (n.uid_ != last_uid) {
        uid2id.at(last_uid) = idx;
      }

      // remove edges incident to the node while they still exist in the graph
      while (n.degree() != 0) {
        remove_edge(*n.edge_begin());
      }
      
      // return 1 to indicate that we were able to remove a node
      return 1;
    } // end of if clause
    else {
      // return 0 to indicate that we were unable to remove a node
      return 0;
    } // end of else clause
  }

  /** Remove the node that the specified node iterator points to from a graph.
   * In addition, this method shifts later nodes down and returns an iterator to
   * the formerly next element and invalidates all iterators in the range
   * [@a n_it, end()].
   * @param[in] n_it Node iterator pointing to an element in the node vector
   * @return         The same iterator (i.e., @a n_it), but now pointing to old
   *                  *(@a n_it + 1) (i.e., the formerly next element).
   * 
   * @pre @a n_it points to a node in a node vector
   * @pre num_nodes() >= 1
   * @post If old has_node(@a n) == true, new num_nodes() = old num_nodes() - 1.
   *        Else, new num_nodes() = old num_nodes().
   * @post has_node (*(@a n_it)) == false
   * @post for all iterators in the range [@a n_it, old_end() - 1]:
              new *@a n_it = old *(@a n_it + 1)
   *
   * Complexity: Average case O(1). Worst case: O(n).
   */
  node_iterator remove_node(node_iterator n_it) {
    // use the `remove_node` function that takes a node by dereferencing the
      // specified node iterator
    remove_node(*n_it);
    // return the iterator, which nows points to a different node (i.e., formerly
      // the last node).
    return n_it;
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
    /** Default public constructor that creates an invalid edge.
     * @return An invalid edge whose member variables are default initialized.
     */
    Edge() {
    }

    /** Return a node of this Edge 
    * @return One of the nodes associated with this edge.
    *
    * Complexity: O(1)
    */
    Node node1() const {
      // return the node corresponding to node1's uid
      return Node(graph_, n1_uid_);
    }

    /** Return the other node of this Edge 
     * @return The other node associated with this edge. 
     *
     * Complexity: O(1)
     */
    Node node2() const { 
      // return the node corresponding to node2's uid
      return Node(graph_, n2_uid_);
    }

    /** Test whether this edge and @a e are equal.
     * @param e An Edge to be compared with the current edge
     * @return True if @a e and the current edge are part of the same graph and
     *    connect the same nodes (order does not matter as edges are undirected).
     *    False if either their graphs are different or they connect different nodes.
     *
     * Complexity: O(1)
     */
    bool operator==(const Edge& e) const {
      // check that the 2 edges have pointers to the same graph and the same uid
        // (this check suffices b/c the Graph class's add_edges() functions only
        // allows for one edge uid between every pair of nodes)
      if (graph_ == e.graph_ && edge_uid_ == e.edge_uid_) {
        return true;
      }
      // if the two edges are from different graphs or have different uid's
      else {
        return false;
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     *
     * @param e An Edge to be compared with the current edge
     * @return True if @a e and the current edge are part of the same graph and
     *    the current edge's unique id < @a e's uid. False if @a e's graph is
     *    different from the current edge's graph, the current edge's uid ==
     *    @a e's uid, or the current edge's uid > @a e's uid.
     *
     * Complexity: O(1)
     */
    bool operator<(const Edge& e) const {
      // check that the 2 edges have pointers to the same graph and that this
        // edge's uid is smaller than the uid of the specified edge
      if (graph_ < e.graph_) {
        return true;
      } 
      else if (graph_ == e.graph_ && edge_uid_ < e.edge_uid_) {
        return true;
      }
      else {
        return false;
      }
    }

    /** Return this edge's associated value. 
     * @return The edge's value whose type is user-specified. Value is a mutable 
     *    object in this method, allowing us to modify it. Value can be the 
     *    default value, edge_value_type(), if no value was passed to `add_edge`
     *    during the edge's creation. 
     *
     * Complexity: O(1)
     */
    edge_value_type& value() {
        // dereference the graph pointer to access the vector of internal nodes.
        // Then get the element with this Node's index and return its value
        return graph_->edges_vector[edge_uid_].value;
    }

    /** Return this edge's associated value. 
     * @return The edge's value whose type is user-specified. Value is a const 
     *    object in this method, so we can inspect the object that is returned 
     *    but cannot modify it. Value can be the default value, edge_value_type(), 
     *    if no value was passed to `edge_node` during the node's creation. 
     *
     * Complexity: O(1)
     */
    const edge_value_type& value() const {
        // dereference the graph pointer to access the vector of internal nodes.
        // Then get the element with this Node's index and return its value
        return graph_->edges_vector[edge_uid_].value;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Pointer back to the Graph container
    Graph* graph_;
    // The edge's unique identification number
    size_type edge_uid_;
    // Node 1's unique identification number
    size_type n1_uid_;
    // Node 2's unique identification number
    size_type n2_uid_;
    // Private Constructor
    Edge(const Graph* graph, size_type edge_uid, size_type n1_uid, size_type n2_uid) 
        : graph_(const_cast<Graph*>(graph)), edge_uid_(edge_uid), n1_uid_(n1_uid),
          n2_uid_(n2_uid) {}
  };

  /** Return the total number of edges in the graph.
   * @return An integer specifying this graph's total number of edges.
   * 
   * Complexity: O(1).
   */
  size_type num_edges() const {
    // get the current number of elements in the unordered map where the key is
      // edge uid
    return edges_vector.size();
  }

  /** Return the edge with index @a i.
   * @param i An integer specifying the index of the edge to be returned.
   * @return The Edge with index @a i.
   * 
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(1).
   */
  Edge edge(size_type i) const {
    // check that the specified edge index exists in this graph
    assert(i < num_edges());
    // get the uid's for the corresponding node 1 and node 2
    size_type n1_uid = edges_vector[i].n1_uid;
    size_type n2_uid = edges_vector[i].n2_uid;
    // return the edge with the specified index
    return Edge(this, i, n1_uid, n2_uid);
  }

 /** Test whether two nodes are connected by an edge.
   * @param a A Node object
   * @param b Another Node object
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * @pre @a a and @a b are distinct nodes
   *
   * Complexity: Average case O(1). Worst case O(2 * num_nodes())
   */
  bool has_edge(const Node& a, const Node& b) const {
    // check that the two nodes both have pointers to this graph
    assert(a.graph_ == this && b.graph_ == this);
    // check if a's uid exists in the unordered map where keys are one of the edge's
      // node uid's. If so, see if the nested umap contains b's uid as the key.
    if (node_uid_umap.count(a.uid_)
      && node_uid_umap.at(a.uid_).count(b.uid_)) {
    return true;
  }
    // there are no edges connecting node a and node b
    else {
      return false;
    }
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @param a A Node object
   * @param b Another Node object
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * 
   * @pre @a a and @a b are distinct nodes
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: Average case O(1). Worst case O(2 * num_nodes()).
   */
  Edge add_edge(const Node& a, const Node& b, 
                const edge_value_type& value = edge_value_type()) {
    // use the has_edge() method to see if 2 nodes are already connected by an edge
    if (has_edge(a, b)) {
      // get the uid of the edge that already exists between the 2 nodes
      size_type edge_uid = node_uid_umap.at(a.uid_).at(b.uid_);
      // return the Edge corresponding to this edge uid
      return Edge(this, edge_uid, a.uid_, b.uid_);
    }
    // if the two specified nodes aren't already connected by an edge
    else {
      // get the uid for a new edge
      size_type edge_uid = num_edges();
      // update the vector of internal edges with a new internal edge
      edges_vector.push_back(internal_edge(edge_uid, a.uid_, b.uid_, value));
      // update the unordered map that has node 1's uid as the key with a new edge
        // first use a's uid as the key
      node_uid_umap[a.uid_][b.uid_] = edge_uid;
        // then use b's uid as the key since node order doesn't matter for edges
      node_uid_umap[b.uid_][a.uid_] = edge_uid;
      // return the Edge corresponding to this new edge uid
      return Edge(this, edge_uid, a.uid_, b.uid_);
    }
  }

  /** Remove the edge connecting two specified nodes from a graph. Return 1 if
   * successful and 0 if not. If the specified edge is found in the graph, this
   * method invalidates that specific edge. 
   * @param[in] a A Node object 
   * @param[in] b Another Node object 
   * @return      The integer 1 if an edge was found between @a a and @a b and
   *                invalidated and the integer 0 if an edge was not found.
   *
   * @pre @a a and @a b are distinct nodes
   * @pre num_edges() >= 1
   * @post If old has_edge(@a a, @a b) == true, new num_edges() = old num_edges() - 1.
   *        Else, new num_edges() = old num_edges().
   * @post has_edge(@a a, @a b) == false
   *
   * Complexity: Average case O(1). Worst case O(2 * num_nodes()).
   */
  size_type remove_edge(const Node& a, const Node& b) {
    if (has_edge(a, b)) {
      // get the uid of the edge that already exists between the 2 nodes
      size_type edge_uid = node_uid_umap.at(a.uid_).at(b.uid_);
      
      // create temp variables to store info. from the internal edge that's
        // currently last in the edges vector
      internal_edge last_edge = *(edges_vector.end() - 1);
      size_type last_n1_uid = last_edge.n1_uid;
      size_type last_n2_uid = last_edge.n2_uid;
      edge_value_type last_value = last_edge.value;

      // overwrite the edge we want to remove with data from the internal 
        // edge that's currently last in the vector
      edges_vector[edge_uid].n1_uid = last_n1_uid;
      edges_vector[edge_uid].n2_uid = last_n2_uid;
      edges_vector[edge_uid].value = last_value;

      // remove the last internal edge since it's data is now duplicated
      edges_vector.pop_back();

      // remove the overwritten edge from the unordered map of edges
        // need to erase 2x (once with node a as key and once with node b as key)
      node_uid_umap[a.uid_].erase(b.uid_);
      node_uid_umap[b.uid_].erase(a.uid_);

      // in the edge unordered map, updated the formerly last edge with its new
        // edge uid if we didn't just delete the last edge in the step above
      if ((edge_uid != last_edge.edge_uid) && edges_vector.size() != 0) {
        node_uid_umap[last_n1_uid][last_n2_uid] = edge_uid;
        node_uid_umap[last_n2_uid][last_n1_uid] = edge_uid;
      }

      // return 1 to indicate that we were able to remove an edge
      return 1;
    } // end of if clause
    else {
      // return 0 to indicate that we were unable to remove an edge
      return 0;
    } // end of return clause
  }

  /** Remove the specified edge from a graph, returning 1 if successful and 0 if not.
   * If the specified edge is found in the graph, this method invalidates that edge.
   * @param[in] e An Edge object
   * @return      The integer 1 if @a e was found and invalidated and the integer
   *                 0 if the edge was not found.
   *
   * @pre @a e must be an edge in a graph
   * @pre num_edges() >= 1
   * @post If old has_edge(@a e.node1(), @a e.node2()) == true, 
   *        new num_edges() = old num_edges() - 1. Else, new num_edges() = old num_edges().
   * @post has_edge(@a e.node1(), @a e.node2()) == false
   *
   * Complexity: Average case O(1). Worst case O(2 * num_nodes())
   */
  size_type remove_edge(const Edge& e) {
    // get the two nodes affiliated with this edge
    Node a = e.node1();
    Node b = e.node2();

    // use the `remove_edge` function that takes two nodes
    return remove_edge(a, b);
  }

  /** Remove the edge that the specified edge iterator points to from a graph.
   * In addition, this method shifts later edges down and returns an iterator to the
   * formerly next edge. It also invalidates all edges in the range [@a e_it, end()].
   * @param[in] e_it Edge iterator pointing to an element in the edge vector
   * @return         The same iterator (i.e., @a e_it), but now pointing to old
   *                  *(@a e_it + 1) (i.e., the formerly next edge).
   *
   * @pre @a e_it must be a valid iterator in this graph
   * @pre num_edges() >= 1
   * @post If old has_edge((@a *e_it).node1(), (@a *e_it).node2()) == true, 
   *        new num_edges() = old num_edges() - 1. Else, new num_edges() = old num_edges().
   * @post has_edge((@a *e_it).node1(), (@a *e_it).node2()) == false
   *
   * Complexity: Average case O(1). Worst case O(2 * num_nodes())
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    // use the `remove_edge` function that takes an edge by dereferencing the
      // specified edge iterator
    remove_edge(*e_it);
    // return the iterator, which is now pointing to a different edge (i.e., the
      // formerly last edge of the edge vector)
    return e_it;
  }

  /** Remove all nodes and edges from this graph.
   * Invalidates all outstanding Node and Edge objects as well as all iterators,
   *  pointers, and references.
   *
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Complexity: O(num_edges())
   */
  void clear() {
    // clear the vectors storing the graph's internal nodes and internal vectors
    nodes_vector.clear();
    edges_vector.clear();
    // clear the vector and unordered map linking node index to node uid
    id2uid.clear();
    uid2id.clear();
    // clear the unordered map storing the relationship between edge uid's and node
      // uid's for this graph
    node_uid_umap.clear();
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
    using internal_node_type= internal_node;
    /** Default public constructor. Creates an invalid NodeIterator. 
     * @return An invalid node iterator whose member variables are default initialized.
     */
    NodeIterator() {
    }

    /** Dereference a node iterator and return the associated Node. 
     * @return The Node object that the current node iterator is pointing to.
     *
     * Complexity: O(1).
     */
    Node operator*() const {
      // dereference the iterator to get the current node it's pointing to
      // create a new node using the graph pointer and the node's uid
      return Node(graph_, *iter_);
    }

    /** Increment the node iterator. 
     * @return old node iterator + 1
     * 
     * Complexity: O(1).
     */
    node_iterator& operator++() {
      ++iter_;
      return *this;
    }

    /** Test whether this node iterator and the node iterator being read in are equal.
     * @param ni A node iterator object to be compared to the current node iterator
     * @return True if *(@a ni) == *(current node iterator).
     *
     * Complexity: O(1).
     */
    bool operator==(const node_iterator& ni) const {
      // check that the 2 node iterators have pointers to the same graph and that
        // they point to the same element in the internal node vector
      if (graph_ == ni.graph_ && iter_ == ni.iter_) {
        return true;
      }
      // if the two iterators are from different graphs or point to different
        // elements in the internal node vector
      else {
        return false;
      }
    }

   private:
    friend class Graph;

    // Pointer back to the Graph container
    Graph* graph_;
    // an iterator to a vector of valid node uid's
    std::vector<size_type>::const_iterator iter_;
    // Private Constructor
    NodeIterator(const Graph* graph, 
      std::vector<size_type>::const_iterator iter) 
      : graph_(const_cast<Graph*>(graph)), iter_(iter) {}
  };

  /** Return a node iterator pointing to the first element of the node vector
   * @return node iterator pointing to the first element of the node vector
   *
   * Complexity: O(1).
   */
  node_iterator node_begin() const {
    // create a node iterator pointing to the 1st element of the internal node vector
    return NodeIterator(this, id2uid.begin()); 
  }

  /** Return a node iterator pointing to the last element of the node vector.
   * @return node iterator pointing to the last element of the node vector
   *
   * Complexity: O(1).
   */
  node_iterator node_end() const {
    // create a node iterator pointing to the last element of the internal node vector
    return NodeIterator(this, id2uid.end()); 
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

    /** Construct an invalid IncidentIterator. 
     * @return An invalid incident iterator whose member variables are default 
                  initialized.
     */
    IncidentIterator() {
    }

    /** Dereference an incident iterator and return the associated edge. 
     * @return The Edge object that the current incident iterator is pointing to.
     *
     * Complexity: O(1).
     */
    Edge operator*() const {
      // create a new edge by getting the graph pointer, edge id, node 1 id, and
        // and node 2 id from the class's private members
      return Edge(graph_, iter_->second, node_id_, iter_->first);
    }

    /** Increment the incident iterator. 
     * @return old incident iterator + 1
     * 
     * Complexity: O(1).
     */
    incident_iterator& operator++() {
      ++iter_;
      return *this;
    }

    /** Test whether this incident iterator and the incident iterator being read
     *  in are equal.
     * @param iit A incident iterator to be compared to the current incident iterator
     * @return True if *(@a iit) == *(current node iterator), meaning both iterators
     *             point to the same Edge object.
     *
     * Complexity: O(1).
     */
    bool operator==(const incident_iterator& iit) const {
      // check that the 2 incident iterators point to the same graph, point to the
        // same element in the unordered map, and have the same spawning node
      if (graph_ == iit.graph_ && iter_ == iit.iter_ && node_id_ == iit.node_id_) {
        return true;
      }
      // if the two iterators are from different graphs, point to different incident
        // edges, or are associated with different spawning nodes
      else {
        return false;
      }
    }

   private:
    friend class Graph;

    // Pointer back to the Graph container
    Graph* graph_;
    // The id of the node for which we are looking for incident edges
    size_type node_id_;
    // an iterator for an unordered map whose keys are the id's of incident nodes
      // and values are the id's of the connecting edges
    std::unordered_map<size_type, size_type>::const_iterator iter_;
    // Private Constructor
    IncidentIterator(const Graph* graph, size_type node_id,
      std::unordered_map<size_type, size_type>::const_iterator iter) 
      : graph_(const_cast<Graph*>(graph)), node_id_(node_id), 
        iter_(iter) {}
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
    using internal_edge_type= internal_edge;

    /** Construct an invalid EdgeIterator. 
     * @return An invalid edge iterator whose member variables are default initialized.
     */
    EdgeIterator() {
    }

    /** Dereference an edge iterator and return the associated Edge. 
     * @return The Edge object that the current edge iterator is pointing to.
     *
     * Complexity: O(1).
     */
    Edge operator*() const {
      // create a new edge by dereferencing iterators for the edge and node uid's
      return Edge(graph_, iter_->edge_uid, iter_->n1_uid, iter_->n2_uid);
    }

    /** Increment the edge iterator. 
     * @return old edge iterator + 1
     * 
     * Complexity: O(1).
     */
    edge_iterator& operator++() {
      ++iter_;
      return *this;
    }

    /** Test whether this edge iterator and the edge iterator being read in are equal.
     * @param eit A edge iterator to be compared to the current edge iterator 
     * @return True if *(@a eit) == *(current edge iterator), meaning both iterators
     *             point to the same Edge object.
     *
     * Complexity: O(1).
     */
    bool operator==(const edge_iterator& eit) const {
      // check that the 2 edge iterators have pointers to the same graph and that
        // they point to the same element in the edge id unordered map
      if (graph_ == eit.graph_ && iter_ == eit.iter_) {
        return true;
      }
      // if the two iterators are from different graphs or point to different
        // elements in unordered map
      else {
        return false;
      }
    }

   private:
    friend class Graph;

    // Pointer back to the Graph container
    Graph* graph_;
    // an iterator to the vector containing internal edges
    typename std::vector<internal_edge_type>::const_iterator iter_;
    // Private Constructor
    EdgeIterator(const Graph* graph, 
      typename std::vector<internal_edge_type>::const_iterator iter)
      : graph_(const_cast<Graph*>(graph)), iter_(iter) {}
  };

  /** Return an edge iterator pointing to the first edge of the edge vector.
   * @return edge iterator pointing to the first element of the vector
   *          containing internal edges.
   *
   * Complexity: O(1).
   */
  edge_iterator edge_begin() const {
    // create an edge iterator pointing to the 1st element of the edge vector
    return EdgeIterator(this, edges_vector.begin());
  }

  /** Return a edge iterator pointing to the last element of the edge vector
   * @return edge iterator pointing to the last element of the vector
   *          containing internal edges.
   *
   * Complexity: O(1).
   */
  edge_iterator edge_end() const {
    // create an edge iterator pointing to the last element of the edge vector
    return EdgeIterator(this, edges_vector.end()); 
  }

 private:

  /** @struct Graph::internal_node
   * @brief Internal node class for storing data about nodes (i.e., their position,
   *          unique id, and value). 
   */
  struct internal_node {
    /** Construct an internal node object. 
     * @param p Point object containing the point's position in 3D space
     * @param unique_id The unique id of the node
     * @param val The value of the node.
     * @return An internal node object whose member variables are initialized
     *           with the specified values.
     */
    internal_node(Point p, size_type unique_id, node_value_type val)
      : point(p), 
        uid(unique_id),
        value(val)
    {}

    Point point; // the node's position
    size_type uid; // the node's unique id
    node_value_type value; // the node's value
  };

  /** @struct Graph::internal_edge
   * @brief Internal edge class for storing data about edges (i.e., the uid's of
   *          their two adjacent nodes, the edge's uid, and its value).
   */
  struct internal_edge {
    /** Construct an internal edge object. 
     * @param[in] e_id The unique id of the edge
     * @param[in] n1_id The unique id of one of its adjacent nodes
     * @param[in] n2_id The unique id of its other adjacent
     * @param[in] val The value of the edge
     * @return An internal edge object whose member variables are initialized with
     *            the specified values.
     */
    internal_edge(size_type e_id, size_type n1_id, size_type n2_id, edge_value_type val)
      : edge_uid(e_id),
        n1_uid(n1_id), 
        n2_uid(n2_id),
        value(val)
    {}

    size_type edge_uid; // the uid of the edge 
    size_type n1_uid; // the uid of one of the edge's nodes
    size_type n2_uid; // the uid of the other one of the edge's nodes
    edge_value_type value; // the edge's value
  };

  // a vector containing our internal nodes
  std::vector<internal_node> nodes_vector;

  // a vector tracking the uid's of nodes that haven't been removed
  std::vector<size_type> id2uid;

  // unordered map where keys are a node's uid and values are the node's index
    // in the id2uid vector
  std::unordered_map<size_type, size_type> uid2id;

  // a vector containing our internal edges
  std::vector<internal_edge> edges_vector;

  // an unordered map where the key is node 1's unique id and the value is a
    // nested unordered map where the key is node 2's uid and value is the edge id
  std::unordered_map<size_type, std::unordered_map<size_type, size_type>> 
    node_uid_umap;

};

#endif // CME212_GRAPH_HPP
