#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <functional>
#include <iterator>
#include <memory>
#include <list>
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template  <typename V, typename E>
class Graph {

 public:
	 using  node_value_type = V;
	 using  edge_value_type = E;
	 class Node;
 private:
	 
  // HW0: YOUR CODE HERE
	 struct internal_node_element {
		 Point position;   // The text held by an element
		 // unsigned uid;      // The unique identifcation for an element
		 node_value_type val;
		 unsigned index;
		 std::shared_ptr<std::unordered_set<unsigned>> edges;
	 };

	 struct internal_edge_element {
		 std::shared_ptr<Node> node_1;
		 std::shared_ptr<Node> node_2;
		 // unsigned uid;
		 std::list<unsigned>::iterator index;
		 // unsigned index;
		 edge_value_type val;
	 };


	 // std::vector<internal_node_element> nodes;
	 std::vector<unsigned> nodes;
	 // std::vector<internal_edge_element> edges;
	 std::list<unsigned> edges;

	 unsigned next_node_uid;
	 unsigned next_edge_uid;
	 // std::unordered_map<unsigned, unsigned> nodes_uid_map;
	 std::unordered_map<unsigned, internal_node_element> nodes_uid_map;
	 //std::unordered_map<unsigned, unsigned> edges_uid_map;	 
	 std::unordered_map<unsigned, internal_edge_element> edges_uid_map;
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
  Graph() : next_node_uid(0), next_edge_uid(0) {
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

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
		return fetch_node().position;
    }

	/** Return this node's position. */
	Point& position() {
		return fetch_node().position;
	}

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
		// return ref_graph->nodes_uid_map.at(node_uid);
		return ref_graph->nodes_uid_map.at(node_uid).index;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
	/** return the value for a given node.
	* @return a reference to the value of a given node
	*/
	node_value_type& value() {
		return fetch_node().val;
	}
    
	/** Return the value for a given node.
	* @return a constant reference to the value of a given node
	*/
	const node_value_type& value() const {
		return fetch_node().val;
	}
    
	/** Return the the number of edges incident to a given node.
	* @return a value equal to the number of edges incident to a given node.
	*/
	size_type degree() const {
		return fetch_node().edges->size();
	}
    
	/** Return the the number of edges incident to a given node.
	* @return an iterator pointing to the beginning of a collection of edges
	* for a given node.
	*/
	incident_iterator edge_begin() const {
		return IncidentIterator(this, (*(fetch_node().edges)).begin());
	}
    
	/** Return the the number of edges incident to a given node.
	* @return an iterator pointing to the end of a collection of edges
	* for a given node.
	*/
	incident_iterator edge_end() const {
		return IncidentIterator(this, fetch_node().edges->end());
	}

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
		if (ref_graph != n.ref_graph) return false;

		return index() == n.index();
      //(void) n;          // Quiet compiler warning
      //return false;
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
		if (ref_graph == n.ref_graph) return node_uid < n.node_uid;
		std::less<graph_type*> compare;
		return compare(ref_graph, n.ref_graph);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
	size_type node_uid;
	graph_type* ref_graph;
	Node(const graph_type* reference_graph, size_type uid)
		: node_uid(uid), ref_graph(const_cast<graph_type*>(reference_graph)) {
	}
	internal_node_element& fetch_node() const {
		// return ref_graph->nodes[ref_graph->nodes_uid_map.at(node_uid)];
		return ref_graph->nodes_uid_map.at(node_uid);
	}
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
    return nodes.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] v The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const  node_value_type& v = node_value_type()) {
    // HW0: YOUR CODE HERE
	/*size_type n_uid = next_node_uid;
	next_node_uid++;
	Node n(this, n_uid);
	internal_node_element n_el;
	n_el.position = position;
	n_el.uid = n_uid;
	n_el.edges = std::make_shared<std::unordered_set<unsigned>>();
	n_el.val = v;
	nodes.push_back(n_el);
	nodes_uid_map[n_uid] = n.ref_graph->size()-1;
	return n;  */ 

	size_type n_uid = next_node_uid;
	next_node_uid++;
	Node n(this, n_uid);
	nodes.push_back(n_uid);
	internal_node_element n_el;
	n_el.position = position;
	n_el.index = nodes.size() - 1;
	n_el.edges = std::make_shared<std::unordered_set<unsigned>>();
	n_el.val = v;	
	nodes_uid_map[n_uid] = n_el;
	return n;  
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
	if (this != n.ref_graph) return false;
	
	auto it = nodes_uid_map.find(n.node_uid);
	if (it != nodes_uid_map.end()) return true;
    //(void) n;            // Quiet compiler warning
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
	  // return Node(this, nodes[i].uid);
	  return Node(this, nodes[i]);
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

	/** return the value for a given node.
	* @return a reference to the value of a given edge
	*/
	edge_value_type& value() {
		return fetch_edge().val;
	}
		
	/** return the value for a given node.
	* @return a reference to the value of a given edge
	*/
	const edge_value_type& value() const {
		return fetch_edge().val;
	}

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE		
		Node x = *(fetch_edge().node_1);
		if (ref_node == nullptr) {
			return x;
		}
		assert(*ref_node == *(fetch_edge().node_1) || *ref_node == *(fetch_edge().node_2));
		if (x == *ref_node) {
			return x;
		}
		return *(fetch_edge().node_2);
    }

    /** Return the other node of this Edge */
    Node node2() const {
		Node x = *(fetch_edge().node_2);
		if (ref_node == nullptr) {
			return x;
		}
		assert(*ref_node == *(fetch_edge().node_1) || *ref_node == *(fetch_edge().node_2));
		if (x == *ref_node) {
			return *(fetch_edge().node_1);
		}
		return x;
		
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
		return ((node1() == e.node1()) && (node2() == e.node2())) ||
			((node1() == e.node2()) && (node2() == e.node1()));
      //(void) e;           // Quiet compiler warning
      //return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
		// This function compares the minimum of the two nodes for each edge.
		// If the minimum node for both edges are equal, it then compares
		// the second node for each edge.
		if (ref_graph == e.ref_graph) return edge_uid < e.edge_uid;
		std::less<graph_type*> compare;
		return compare(ref_graph, e.ref_graph);
		return ref_graph < e.ref_graph;
		/*if (operator==(e)) return false;
		if (node1() <= node2()) {
			if (e.node1() <= e.node2()) {
				return (node1() < e.node1()) || 
					((node1() == e.node1()) && (node2() < e.node2()));
			}
			return (node1() < e.node2()) ||
				((node1() == e.node2()) && (node2() < e.node1()));
		}

		if (e.node1() <= e.node2()) {
			return (node2() < e.node1()) ||
				((node2() == e.node1()) && (node1() < e.node2()));
		}
		return (node2() < e.node2()) ||
			((node2() == e.node2()) && (node1() < e.node1()));*/
    }

	double length() const {
		return norm(node2().position() - node1().position());
	}

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
	// Use this space to declare private data members and methods for Edge
	// that will not be visible to users, but may be useful within Graph.
	// i.e. Graph needs a way to construct valid Edge objects

	size_type edge_uid;
	graph_type* ref_graph;
	node_type* ref_node;
	Edge(const graph_type* reference_graph, size_type uid, const node_type* reference_node = nullptr)
		: edge_uid(uid), ref_graph(const_cast<graph_type*>(reference_graph)), ref_node(const_cast<node_type*>(reference_node)) {
	}

	internal_edge_element& fetch_edge() const {
		// return ref_graph->edges[ref_graph->edges_uid_map.at(edge_uid)];
		return ref_graph->edges_uid_map.at(edge_uid);
	}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    // return edges.size();
	  return edges_uid_map.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
	  // return Edge(this, edges[i].uid);
	  // return Edge(this, edges[i]);
	  auto it = edges.begin();
	  std::advance(it, i);
	  return Edge(this, *it);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
	if (!has_node(a) || !has_node(b)) return false;
	for (unsigned u : *((a.fetch_node().edges))) {
		auto it = b.fetch_node().edges->find(u);
		if (it != b.fetch_node().edges->end()) return true;
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

	  if (has_edge(a, b))
		  for (unsigned u : *((a.fetch_node().edges))) {
			  auto it = b.fetch_node().edges->find(u);
			  if (it != b.fetch_node().edges->end()) return Edge(this, *it);
		  }
	  

	  //size_type e_uid = next_edge_uid;
	  //next_edge_uid++;
	  //Edge e(this, e_uid);
	  //internal_edge_element e_el;
	  //e_el.node_1 = std::make_shared<Node>(a);
	  //e_el.node_2 = std::make_shared<Node>(b);
	  //e_el.uid = e_uid;
	  //// edges.push_back(e_el);
	  //a.fetch_node().edges->insert(e_uid);
	  //b.fetch_node().edges->insert(e_uid);
	  //// edges_uid_map[e_uid] = e.ref_graph->num_edges() - 1;
	  //edges_uid_map[e_uid] = e_el;

	  //return e;   

	  size_type e_uid = next_edge_uid;
	  next_edge_uid++;
	  Edge e(this, e_uid);
	  edges.push_back(e_uid);
	  internal_edge_element e_el;
	  e_el.node_1 = std::make_shared<Node>(a);
	  e_el.node_2 = std::make_shared<Node>(b);
	  // e_el.index = edges.size() - 1;
	  e_el.index = --(edges.end());
	  a.fetch_node().edges->insert(e_uid);
	  b.fetch_node().edges->insert(e_uid);
	  // edges_uid_map[e_uid] = e.ref_graph->num_edges() - 1;
	  edges_uid_map[e_uid] = e_el;

	  return e;
  }

  /** Remove an edge from graph if it exists.
  * @pre @a e is a distince edge of the graph.
  * @return
  * @post has_edge(@a a, @a b) == false, where a and b were the two nodes
  * of the edge.
  * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges() - 1.
  *       Else,                        new num_edges() == old num_edges().
  *
  * Can invalidate edge indexes -- in other words, old edge(@a i) might not
  * equal new edge(@a i). Must not invalidate outstanding Edge objects.
  *
  * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
  */
  bool remove_edge(const  Edge& e) {
	  if (edges_uid_map.find(e.edge_uid) == edges_uid_map.end()) return false;

	  edges.erase(edges_uid_map[e.edge_uid].index);
	  e.node1().fetch_node().edges->erase(e.edge_uid);
	  e.node2().fetch_node().edges->erase(e.edge_uid);
	  edges_uid_map.erase(e.edge_uid);
	  return true;
  }

  /** Remove an edge from graph if it exists.
  * @pre @a a and @a b are distinct valid nodes of this graph
  * @return a boolean value indicating whether the edge was removed.
  * @post has_edge(@a a, @a b) == false
  * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges() - 1.
  *       Else,                        new num_edges() == old num_edges().
  *
  * Can invalidate edge indexes -- in other words, old edge(@a i) might not
  * equal new edge(@a i). Must not invalidate outstanding Edge objects.
  *
  * Complexity: O(1)
  */
  bool remove_edge(const  Node& a, const  Node& b) {
	  Edge e;
	  if (has_edge(a, b))
		  for (unsigned u : *((a.fetch_node().edges))) {
			  auto it = b.fetch_node().edges->find(u);
			  if (it != b.fetch_node().edges->end()) e = Edge(this, *it);
		  }
	  else return false;

	  return remove_edge(e);
  }

  /** Remove a node from the graph.
  * @param[in] n The node to remove
  * @return a boolean value indicating whether the node was removed.
  * @post new num_nodes() == old num_nodes() - 1
  * @post All nodes with indexes past the index of @a n are shifted down by 1.
  * @post If old has_node(@a n), new num_nodes() == old num_nodes() - 1.
  *       Else,                        new num_nodes() == old num_nodes().
  * Complexity: O(m), where m is the number of edges incident to the node.
  */
  bool remove_node(const  Node & n) {
	  if (!has_node(n)) return false;
	  // Erase all edges connected to node
	  auto it = n.fetch_node().edges->begin();
	  while (true) {
		  if (it == n.fetch_node().edges->end()) break;
		  assert(edges_uid_map.find(*it) != edges_uid_map.end());
		  Edge e(this, *it, &n);
		  e.node2().fetch_node().edges->erase(*it);
		  edges.erase(edges_uid_map[*it].index);
		  edges_uid_map.erase(*it);
		  it = n.fetch_node().edges->erase(it);	  
	  }

	  // Now erase node
	  size_type index = n.index();
	  nodes[index] = nodes[nodes.size()-1];
	  nodes_uid_map[nodes[index]].index = index;
	  nodes_uid_map.erase(n.node_uid);
	  nodes.erase(--(nodes.end()));
	  return true;

	  /*size_type index = n.index();
	  auto vec_it = nodes.begin() + index;
	  nodes.erase(vec_it);
	  for (; index < nodes.size(); index++) {
		  nodes_uid_map[nodes[index]].index = index;
	  }
	  nodes_uid_map.erase(n.node_uid);
	  return true;*/
  }

  /** Remove an edge from graph if it exists.
  * @param[in] e_it The edge_iterator object pointing to the edge to remove
  * @pre @a e_it is a valid edge iterator to the graph
  * @return an edge_iterator pointing to the next edge in the graph.
  * @post has_edge(e) == false, where @a e was the edge to be removed.
  * @post If old has_edge(@a e), new num_edges() == old num_edges() - 1.
  *       Else,                        new num_edges() == old num_edges().
  *
  * Can invalidate edge indexes -- in other words, old edge(@a i) might not
  * equal new edge(@a i). Must not invalidate outstanding Edge objects.
  *
  * Complexity: O(1)
  */
  edge_iterator remove_edge(edge_iterator e_it) {
	  Edge e = *e_it;
	  auto it = ++(e.fetch_edge().index);
	  assert(remove_edge(e));
	  return EdgeIterator(this, it);
  }

  /** Remove a node from the graph.
  * @param[in] n_it The edge_iterator object pointing to the node to remove
  * @pre @a n_it is a valid node iterator to the graph
  * @return a node_iterator pointing to the next node in the graph.
  * @post If old has_node(@a n), new num_nodes() == old num_nodes() - 1.
  *       Else,                        new num_nodes() == old num_nodes().
  * Complexity: O(m), where m is the number of edges incident to the node.
  */
  node_iterator remove_node(node_iterator n_it) {
	  Node n = *n_it;
	  size_type index = n.index();
	  assert(remove_node(n));
	  return NodeIterator(this, index);
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
	  // nodes = std::vector<internal_node_element>();
	  // edges = std::vector<internal_edge_element>();
	  nodes = std::vector<unsigned>();
	  edges = std::list<unsigned>();
	  next_node_uid = 0;
	  next_edge_uid = 0;
	  nodes_uid_map = std::unordered_map<unsigned, internal_node_element>();
	  // edges_uid_map = std::unordered_map<unsigned, unsigned>();
	  edges_uid_map = std::unordered_map<unsigned, internal_edge_element>();
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
	/** Reference the Node stored in stored in the iterator.
	* @pre the iterator object is a valid iterator and is in a valid state
	* (meaning it's not referencing 1 past the last element in the collection)
	* @pre nodes haven't been removed such that the index of the iterator object isn't
	* larger than one minus the number of nodes.
	* 
	* @return a Node object n
	*
	*/
	Node operator*() const {
		// return Node(ref_graph, ref_graph->nodes[index].uid);
		return Node(ref_graph, ref_graph->nodes[index]);
	}

	/** Increment the iterator to point to the next node object.
	* @pre the iterator object is a valid iterator and is in a valid state
	* (meaning it's not referencing 1 past the last element in the collection)
	* @pre the no nodes have been added or removed since the iterator was constructed.
	* 
	* @post the index 
	* @return the iterator iter with its index incremented by 1.
	*/
	NodeIterator& operator++() {
		index++;
		return *this;
	}
    
	/** Compare this iterator with another NodeIterator object for equality,
	* meaning that the reference graph and the current index of the iterator are
	* the same.
	* @param[in] ni the other node iterator object to compare with for equality.
	* @pre the iterator objects are both in a valid state
	* (meaning neither has an index larger than the size of the collection of nodes, and neither
	* have indexes less than 0.
	* @pre the iterator parameter is a valid graph.
	* 
	* @return bool indicating whether the iterator objects are equal.
	*/
	bool operator==(const NodeIterator& ni) const {
		return (ref_graph == ni.ref_graph) && (index == ni.index);
	}

	bool operator!=(const NodeIterator& ni) const {
		return !(*this == ni);
	}

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
	unsigned index;
	graph_type* ref_graph;

	NodeIterator(const graph_type* reference_graph,unsigned i) : index(i), ref_graph(const_cast<graph_type*>(reference_graph)) { }
  };

  // HW1 #2: YOUR CODE HERE
  /** Return a node iterator pointing to the beginning of a collection of nodes
  * the same.
  * @pre the underlying graph has at least 1 node in the graph.
  * @pre the iterator parameter is a valid graph.
  *
  * @return a node_iterator object
  */
  node_iterator node_begin() const {
	  return NodeIterator(this,0);
  }

  /** Return a node iterator pointing to the end of a collection of nodes
  * the same.
  * @pre the iterator parameter is a valid graph.
  * @return a node_iterator object
  */
  node_iterator node_end() const {
	  return NodeIterator(this, size());
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
	/** Reference the Edge stored in stored in the iterator.
	* @pre the iterator object is a valid iterator and is in a valid state
	* (the edge list for a node hasn't changed since the iterator was initialized)
	*
	* @return a Edge object
	*/
    Edge operator*() const {
		return Edge(ref_node->ref_graph, *it, ref_node);
	}
    
	/** Increment the iterator to the next Edge in the collection of edges
	* incident to a given node.
	* @pre the iterator object is a valid iterator and is in a valid state
	* (the edge list for a node hasn't changed since the iterator was initialized, and
	* the node is still valid)
	*
	* @return the iterator pointing to the next edge in the collection.
	*/
	IncidentIterator& operator++() {
		it++;
		return *this;
	}
    
	/** Compare this iterator with another Incidence Iterator object for equality,
	* meaning that the reference graph and the current index of the iterator are
	* the same.
	* @param[in] ii the other incident iterator object to compare with for equality.
	*
	* @pre the iterator parameter is a valid graph.
	*
	* @return bool indicating whether the iterator objects are equal.
	*/
	bool operator==(const IncidentIterator& ii) const {
		return (*ref_node == *(ii.ref_node)) && (it == ii.it);
	}

	bool operator!=(const IncidentIterator& ii) const {
		return !(*this == ii);
	}

   private:
    friend class Graph;
	std::unordered_set<unsigned>::iterator it;
	node_type* ref_node;
	IncidentIterator(const node_type* reference_node, const std::unordered_set<unsigned>::iterator& iter) : 
		it(iter), ref_node(const_cast<node_type*>(reference_node)) { }
    // HW1 #3: YOUR CODE HERE
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
	/** Reference the Edge stored in stored in the iterator.
	* @pre the iterator object is a valid iterator and is in a valid state
	* (the edge list for hasn't changed since the iterator was initialized)
	*
	* @return a Edge object
	*/
	Edge operator*() const {
		// return Edge(ref_graph, ref_graph->edges[index].uid);
		return Edge(ref_graph, *index);
	}


	/** Increment the iterator to the next Edge in the collection of edges
	* for the entire graph.
	* @pre the iterator object is a valid iterator and is in a valid state
	* (the edge list for a node hasn't changed since the iterator was initialized)
	*
	* @return the iterator pointing to the next edge in the collection.
	*/
	EdgeIterator& operator++() {
		++index;
		return *this;
	}


	/** Compare this iterator with another Edge Iterator object for equality,
	* meaning that the reference graph and the current index of the iterator are
	* the same.
	* @param[in] ei the other incident iterator object to compare with for equality.
	*
	* @pre the iterator parameter is a valid graph.
	*
	* @return bool indicating whether the iterator objects are equal.
	*/
	bool operator==(const EdgeIterator& ei) const {
		return (ref_graph == ei.ref_graph) && (index == ei.index);
	}

	bool operator!=(const EdgeIterator& ei) const {
		return !(*this == ei);
	}

   private:
    friend class Graph;
	// unsigned index;
	std::list<unsigned>::iterator index;
	graph_type* ref_graph;
	EdgeIterator(const graph_type* reference_graph, const std::list<unsigned>::iterator& i) :
		index(i), ref_graph(const_cast<graph_type*>(reference_graph)) { }
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Return an edge iterator pointing to the beginning of a collection of edges
  * for the graph.
  * @pre the iterator parameter is a valid graph.
  *
  * @return an edge_iterator object
  */
  edge_iterator edge_begin() {
	  std::list<unsigned>::iterator it = edges.begin();
	  return EdgeIterator(this, it);
  }

  /** Return an edge iterator pointing to the end of a collection of edges
  * for the graph.
  * @pre the iterator parameter is a valid graph.
  * @return an edge_iterator object
  */
  edge_iterator edge_end() {
	  std::list<unsigned>::iterator it = edges.end();
	  return EdgeIterator(this, it);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
