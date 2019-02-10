//I have worked with Benoit Zhou and Charbel Saad for this assignement

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
template<typename V>
class Graph{
 private:

 public:

  // PUBLIC TYPE DEFINITIONS

  using node_value_type = V;
  /** Type of this graph. */
  using graph_type = Graph<V>;

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

  // CONSTRUCTORS AND DESTRUCTOR
  /** Construct an empty graph. */
  
  Graph() : size_(0), next_uid_(0), nodes_(), nb_edges_(0), next_edge_uid_(0), edges_(){}

  /** Default destructor */
  ~Graph() = default;

  struct internal_nodes;

  struct internal_position{	
    Point position;
    size_type uid;
    
    //vector of incident edges
    std::vector<internal_nodes*> in_edges;
    size_type degree;

    //used for BFS
    bool mark = false;
    node_value_type value;
  };

  size_type size_;
  size_type next_uid_;
  //Vector of Points
  std::vector<internal_position> nodes_;


  size_type nb_edges_;
  size_type next_edge_uid_;
  //Vector of Nodes
  std::vector<internal_nodes> edges_;


  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node : private totally_ordered<Node> {

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    
    //Pointer back to the Graph containing the node
    Graph* graph_;
    //Node unique identifier
    size_type uid_;

   public:
   Node() {
    }

    const Point& position() const {
      return fetch().position;
    }


    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }


    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value() { return fetch().value; }
    const node_value_type& value() const { return fetch().value; }

    size_type degree() const{ return fetch().degree; }
    incident_iterator edge_begin() const { return IncidentIterator(graph_, 0, this);}
    incident_iterator edge_end() const { return IncidentIterator(graph_, degree(), this); }


    bool operator==(const Node& n) const {
      if((uid_ == n.uid_) and (graph_ == n.graph_))
        return true;
      return false;
    }

    bool operator<(const Node& n) const {
      if(fetch().position.x < n.fetch().position.x){return true;}
      else if (fetch().position.x > n.fetch().position.x){return false;}
      else if (fetch().position.y < n.fetch().position.y){return true;}
      else if (fetch().position.y > n.fetch().position.y){return false;}
      else if (fetch().position.y < n.fetch().position.y){return true;}
      else if (fetch().position.z > n.fetch().position.z){return false;}

      return false;
    }

    internal_position& fetch() const {
          return graph_->nodes_[uid_];

      assert(false);
    }

   private:
    /** Private Constructor */
    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }


  };
  //End Node Class

  
/** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return size_;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  Node add_node(const Point& position, const node_value_type& = node_value_type()) {
    // Create a new vector of internal elemts
    internal_position new_node;
    new_node.position = position;
    new_node.uid = next_uid_;
    nodes_.push_back(new_node);

    ++size_;
    ++next_uid_;

    return Node(this, next_uid_-1);
  }

  /** Determine if a Node belongs to this Graph
   */
  bool has_node(const Node& n) const {
    if(this == n.graph_)
      return true;
    return false;
  }

  /** Return the node with index @a i.
   */
  Node node(size_type i) const {
    return Node(this, i);
  }


 public:
  struct internal_nodes{	
    Node* node1_;
    Node* node2_;
    size_type uid_;
  };


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

    //This is a dummy variable that allows us to have node1 the spawning node
    //for incident_iterator
    int trick = 0;

    /** Return a node of this Edge */
    Node node1() const {
      if(trick == 0) { return *(fetch_nodes().node1_); }
      else{ return *(fetch_nodes().node2_); }
    }

    /** Return the other node of this Edge */
    Node node2() const {
      if(trick == 0) { return *(fetch_nodes().node2_); }
      else{ return *(fetch_nodes().node1_); }
    }

   private:
    Graph* graph_;
    size_type edge_uid_;

    /** Private Constructor */
    Edge(const Graph* graph, size_type edge_uid)
        : graph_(const_cast<Graph*>(graph)), edge_uid_(edge_uid) {
    }

    internal_nodes& fetch_nodes() const {
      return graph_->edges_[edge_uid_];

      assert(false);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
   
   public:
    bool operator==(const Edge& e) const {
      if( ((this->node1() == e.node1()) and (this->node2() == e.node2())) 
        or ((this->node1() == e.node2()) and (this->node2() == e.node1())) ) 
        return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if(this->node1() < e.node1())
        return true;
      else if(e.node1() < this->node1())
        return false;
      //in this case the two minimal nodes are equal
      else{
        if(this->node2() < e.node2())
          return true;
      }

      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
  };


  // Return the total number of edges in the graph.
  size_type num_edges() const {
    return nb_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (size_type i = 0; i < this->num_edges(); ++i)
      if( ((this->edges_[i].node1_ == &a) and (this->edges_[i].node2_ == &b)) 
         or ((this->edges_[i].node1_ == &b) and (this->edges_[i].node2_ == &a)) ) 
        return true;
   
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

  // Edge add_edge(const Node& a, const Node& b) {
  //Implementation detail: we build this method so as to always have node1 < node2
  // for the order relation defined for nodes
  Edge add_edge(Node& a, Node& b) {
    //be sure there is no self-loop
    assert((a == b) == false);



    internal_nodes new_edge;
    
    if(a < b){
      new_edge.node1_ = &a;
      new_edge.node2_ = &b;
    }
    else{
      new_edge.node1_ = &b;
      new_edge.node2_ = &a;
    }

    new_edge.uid_ = next_edge_uid_;
    edges_.push_back(new_edge);
   
    //adding incident edges to nodes 
    nodes_[a.uid_].in_edges.push_back(&new_edge);
    ++nodes_[a.uid_].degree;
    nodes_[b.uid_].in_edges.push_back(&new_edge);
    ++nodes_[b.uid_].degree;
   
     
    ++nb_edges_;
    ++next_edge_uid_; 
    return Edge(this, next_edge_uid_ -1);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    size_ = 0;
    next_uid_ = 0;
    nb_edges_ = 0;
    next_edge_uid_ = 0;

    nodes_.clear();
    edges_.clear();
  }


  class NodeIterator {
   private:
    Graph* graph;
    size_type node_uid;
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const
    NodeIterator(const Graph* graph, size_type uid) : 
         graph(const_cast<Graph*>(graph)), node_uid(uid) {}

    Node operator*() const{ return Node(graph, node_uid); }
    NodeIterator& operator++(){
      ++node_uid;
      return *this;
    }
    bool operator==(const NodeIterator& node_iter) const { return node_uid == node_iter.node_uid; }
    bool operator!=(const NodeIterator& node_iter) const { return node_uid != node_iter.node_uid; }    

   private:
    friend class Graph;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  node_iterator node_begin() const { return NodeIterator(this, 0); }
  node_iterator node_end() const { return NodeIterator(this, size_); }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator {
   private:
    Node* n;
    Graph* graph;
    size_type index;
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    IncidentIterator(const Graph* graph, size_type uid, const Node* n) :
         graph(const_cast<Graph*>(graph)), index(uid), n(const_cast<Node*>(n)) {}

    Edge operator*() const{
      Edge e = Edge(graph, index);
      if(e.node1().uid_ == (*n).uid_) {}
      else{ e.trick = 1; } 
      return e;
    }
    incident_iterator& operator++(){
      ++index;
      return *this;
    }
    bool operator==(const incident_iterator& iit) const { return index == iit.index; }
    bool operator!=(const incident_iterator& iit) const { return index != iit.index; }
    
   private:
    friend class Graph;
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator {
   private:
    Graph* graph;
    size_type edge_uid;

   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator(const Graph* graph, size_type uid) :
         graph(const_cast<Graph*>(graph)), edge_uid(uid) {}

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const { return Edge(graph, edge_uid); }
    EdgeIterator& operator++() {
      ++edge_uid;
      return *this;
    }
    bool operator==(const EdgeIterator& edge_iter) const { return edge_uid == edge_iter.edge_uid; }
    bool operator!=(const EdgeIterator& edge_iter) const { return edge_uid != edge_iter.edge_uid; }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const { return EdgeIterator(this, 0); }
  edge_iterator edge_end() const { return EdgeIterator(this, nb_edges_); }
  

};

#endif // CME212_GRAPH_HPP
