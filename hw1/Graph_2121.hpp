#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <set>

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
  public:
  /** Type of this graph. */
  using graph_type = Graph;

  /** Type of node value. */
  using node_value_type = V;

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


  private:
  size_type graph_size;
  size_type edge_size;
  std::map<size_type,Point> index_to_point;
  std::map<size_type,node_value_type> index_to_value;
  std::map<size_type,std::pair<size_type,size_type> > index_to_edge;
  std::map<size_type,std::set<size_type> > edge_map;
  std::map<std::pair<size_type,size_type>,size_type> edge_to_index;

  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //
  public:
  /** Construct an empty graph. */
  Graph()
  {
    graph_size=0;
    edge_size=0;
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
  class Node: private totally_ordered<Node> {
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

     Graph* owner;
     size_type index_no;

    Node(){}

    Node(const Graph* o,const size_type i): index_no(i)
    {
      owner=const_cast<Graph*>(o);
    }

    Node(const Node& n):index_no(n.index_no)
    {
      owner=n.owner;
    }

    /*Node operator=(const Node n)
    {
      index_no=n.index_no;
      owner=n.owner;
      return *this;
    }*/

    /** Return this node's position. */
    const Point& position() const
    {
      return owner->index_to_point.at(index_no);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const
    {
      return index_no;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    node_value_type& value()
    {
      return owner->index_to_value.at(index_no);
    }

    const node_value_type& value() const
    {
      return owner->index_to_value.at(index_no);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const
    {
      if((owner==n.owner)&&(index_no==n.index_no))
      {
        return true;
      }
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
    bool operator<(const Node& n) const
    {
      if(owner<n.owner)
      {
        return true;
      }
      if((owner==n.owner)&&(index_no<n.index_no))
      {
        return true;
      }
      return false;
    }

    size_type degree() const
    {
      return owner->edge_map.at(index_no).size();
    }

    incident_iterator edge_begin() const
    {
      return incident_iterator(this->owner,index_no,owner->edge_map.at(index_no).begin());
    }

    incident_iterator edge_end() const
    {
      return incident_iterator(this->owner,index_no,owner->edge_map.at(index_no).end());
    }


   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const
  {
    return graph_size;
  }

  /** Synonym for size(). */
  size_type num_nodes() const
  {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position,const node_value_type& val=node_value_type())
  {
    index_to_point.insert({graph_size,position});
    index_to_value.insert({graph_size,val});
    graph_size++;
    return Node(this,graph_size-1);        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const
  {
    if((n.owner==this)&&(n.index_no<graph_size))
    {
      return true;
    }
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const
  {
    assert((i<graph_size)&&(i>=0));  //Precondition
    return Node(this,i);        // Invalid node
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
  class Edge: private totally_ordered<Edge> {
   public:
     const Graph * owner;
     const size_type index_no;

    Edge(){}

    /** Construct an invalid Edge. */
    Edge(const Graph* o,size_type index): index_no(index)
    {
      owner=o;
    }

    /** Return a node of this Edge */
    Node node1() const
    {
      return Node(owner,owner->index_to_edge.at(index_no).first);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const
    {
      return Node(owner,owner->index_to_edge.at(index_no).second);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const
    {
      if((owner==e.owner)&&(index_no==e.index_no))
      {
        return true;
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const
    {
      if(owner<e.owner)
      {
        return true;
      }
      if((owner==e.owner)&&(index_no<e.index_no))
      {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const
  {
    return edge_size;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const
  {
    assert((i<edge_size)&&(i>=0)); //Precondition
    return Edge(this,i);        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const
  {
    assert(has_node(a)&&has_node(b)); //Precondition
    //std::cout<<"Has "<<a.index_no<<" "<<b.index_no<<" ";
    if((a.owner!=this)||(b.owner!=this))
    {
      return false;
    }
    const auto& it = edge_map.find(a.index_no);
    bool detect = it!=edge_map.end() && it->second.find(b.index_no)!=it->second.end();
    return detect;
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
    assert(has_node(a) && has_node(b) && (!(a==b)));
    if (edge_map[a.index_no].find(b.index_no)!=edge_map[a.index_no].end())
    {
      return Edge(this,edge_to_index.at({a.index_no,b.index_no}));
    }
    edge_map[a.index_no].insert(b.index_no);
    edge_map[b.index_no].insert(a.index_no);
    index_to_edge[edge_size] = {a.index_no, b.index_no};
    edge_to_index[{a.index_no,b.index_no}]=edge_size;
    edge_to_index[{b.index_no,a.index_no}]=edge_size;
    edge_size++;
    //std::cout<<"Adding "<<a.index_no<<" "<<b.index_no<<std::endl;
    return Edge(this,edge_size-1);        // Invalid Edge
  }

  //Helper function
  void displaygraph() const
  {
    for(auto it = edge_map.cbegin(); it != edge_map.cend(); ++it)
    {
        std::cout << it->first  <<  "\n";
        for(auto iter = it->second.begin(); iter != it->second.end(); iter++)
        {
          std::cout << *iter;
        }
        std::cout <<  "\n";
    }
    std::cout<<"Listing Edges "<<std::endl;
    for(auto it = index_to_edge.cbegin(); it != index_to_edge.cend(); ++it)
    {
        std::cout << it->first  <<  "  " << it->second.first << "  " << it->second.second <<"\n";
    }
    std::cout<<"Number of edges "<< edge_size <<std::endl;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    edge_map.clear();
    index_to_edge.clear();
    index_to_point.clear();
    graph_size=0;
    edge_size=0;
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
    size_type index;
    const Graph * owner;

    /** Construct an invalid NodeIterator. */
    NodeIterator(size_type ind,const Graph *own)
    {
      owner=own;
      index=ind;
    }

    Node operator*() const
    {
      return Node(owner,index);
    }

    node_iterator& operator++(int)
    {
      index++;
      return *this;
    }

    node_iterator& operator++()
    {
      ++index;
      return *this;
    }

    bool operator==(const node_iterator& other_iterator) const
    {
        assert(owner==other_iterator.owner);
        if(index==other_iterator.index)
        {
          return true;
        }
        return false;
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
  };


node_iterator node_begin() const
{
  return node_iterator(0,this);
}

node_iterator node_end() const
{
  return node_iterator(graph_size-1,this);
}

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator>  {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
    size_type index;
    Graph * owner;
    std::set<size_type>::iterator it;

    /** Construct an invalid IncidentIterator. */
    IncidentIterator(const Graph* o,const size_type in,std::set<size_type>::iterator i):index(in)
    {
      owner=const_cast<Graph*>(o);
      it=i;

    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    Edge operator*() const
    {
      return Edge(owner,owner->edge_to_index.at({index,*it}));
    }

    incident_iterator& operator++()
    {
      ++it;
      return *this;
    }

    incident_iterator& operator++(int)
    {
      it++;
      return *this;
    }

    bool operator==(const incident_iterator& iit) const
    {
      if((it==iit.it)&&(index==iit.index)&&(owner==iit.owner))
      {
        return true;
      }
      return false;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
    size_type index;
    const Graph * owner;

    /** Construct an invalid EdgeIterator. */
    EdgeIterator()
    {
    }

    EdgeIterator(size_type ind,const Graph *own)
    {
      owner=own;
      index=ind;
    }

    Edge operator*() const
    {
      return Edge(owner,index);
    }

    edge_iterator& operator++(int)
    {
      index++;
      return *this;
    }

    edge_iterator& operator++()
    {
      ++index;
      return *this;
    }

    bool operator==(const edge_iterator& other_iterator) const
    {
        assert(owner==other_iterator.owner);
        if(index==other_iterator.index)
        {
          return true;
        }
        return false;
    }


    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  edge_iterator edge_begin() const
  {
    return edge_iterator(0,this);
  }

  edge_iterator edge_end() const
  {
    return edge_iterator(edge_size-1,this);
  }

};

#endif // CME212_GRAPH_HPP
