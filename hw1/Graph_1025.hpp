#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <stdexcept>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <cassert>
#include <iterator>
#include <string>


#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

using namespace std;

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template<typename V>
class Graph {
 private:
  struct internal_node;
  // HW0: YOUR CODE HERE
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
  Graph(){
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
        graph_ = nullptr;
        node_ind_ = size_type(-1);
    }


    /** Return this node's position. */
    const Point& position() const {
        return fetch().position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
        return node_ind_;
    }
        
        
    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
    
    node_value_type& value ()
    {
        return fetch().value;
    }
    
    const node_value_type& value () const
    {
        return fetch().value;
    }
    
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
     bool operator==(const Node& n) const {
            if ((index() == n.index()) && (graph_ == n.graph_)) {
                return true;
            } else {
                return false;
            }
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
            if (index() < n.index()) {
                return true;
            } else {
                return false;
            }
        }

     size_type degree () const
     {
     	auto it = graph_->adj_list.find(node_ind_);
        size_type degree = (it->second).size();
     	return degree;
     }
     incident_iterator edge_begin () const
     {
     	auto it = graph_->adj_list.find(node_ind_);  
     	return incident_iterator(graph_, (it->second).begin(), node_ind_);
     }
     incident_iterator edge_end () const
     {
     	auto it = graph_->adj_list.find(node_ind_); 
     	return incident_iterator(graph_, (it->second).end());
     }
    

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

        Graph* graph_;
        size_type node_ind_;

        /** Private Constructor for Node */
        Node(const Graph* graph, size_type node_ind)
        : graph_(const_cast<Graph*>(graph)), node_ind_(node_ind) {
        }
        
        /** Helper method to return the appropriate node. */
        internal_node& fetch() const {
            return graph_->nodes_[index()];
        }    
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
    size_type size() const {
        return nodes_.size();
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
    Node add_node(const Point& position, const node_value_type &value_ = node_value_type ()) {
        size_type new_index = num_nodes();
        internal_node new_node = {};
        new_node.index = new_index;
        new_node.position = position;
        new_node.value = value_;
        nodes_.push_back(new_node);
        return Node(this, new_index);
    }
  
  
  

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
    bool has_node(const Node& n) const {
        if (n.index() < size() && n.graph_ == this) {
            return true;
        } else {
            return false;
        }
    }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
   Node node(size_type i) const {
        //assert(i < size());
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
       Edge() {
            graph_ = nullptr;
            node1_ind_ = size_type(-1);
            node2_ind_ = size_type(-1);
            edge_ind_ = size_type(-1);
        }

    /** Return a node of this Edge */
        Node node1() const {
            return Node(graph_, node1_ind_);
        }

    /** Return the other node of this Edge */
        Node node2() const {
            return Node(graph_, node2_ind_);
        }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
        bool operator==(const Edge& e) const {
            if ((edge_ind_ == e.edge_ind_) && (graph_ == e.graph_)) {
                return true;
            } else {
                return false;
            }
        }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
        bool operator<(const Edge& e) const {
            if (edge_ind_ < e.edge_ind_) {
                return true;
            } else {
                return false;
            }
        }
        
        
   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

        Graph* graph_;
        size_type node1_ind_;
        size_type node2_ind_;
        size_type edge_ind_;
        

        /** Private Constructor for Edge */
        Edge(const Graph* graph, size_type node1_ind, size_type node2_ind,
             size_type edge_ind) : graph_(const_cast<Graph*>(graph)),
             node1_ind_(node1_ind), node2_ind_(node2_ind), edge_ind_(edge_ind) {
        }
    };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
    size_type num_edges() const {
        return edges_map_.size();
    }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
        Edge edge(size_type i) const {
        assert(i < num_edges());
        
        // get node_pair value for current edge index i
        set<size_type> node_pair = edges_map_flipped_.at(i);
        
        // use iterator to access elements in set
        set<size_type>::iterator it = node_pair.begin();
        advance(it, 0);
        size_type nodea_index = *it;
        advance(it, 1);
        size_type nodeb_index = *it;
        
        return Edge(this, nodea_index, nodeb_index, i);
    }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
        assert(a.graph_ == this && b.graph_ == this);
        
        size_type nodea_index = a.index();
        size_type nodeb_index = b.index();
        
        set<size_type> node_pair;  // create set with node a and b indicies
        node_pair.insert(nodea_index);
        node_pair.insert(nodeb_index);
        
        // check if this node_pair is a key in edges_map_
        if (edges_map_.find(node_pair) == edges_map_.end()) {
            return false;   // key not found
        } else {
            return true;    // key found
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
        assert((a == b) == false);
        assert(a.graph_ == this && b.graph_ == this);
        
        size_type nodea_index = a.index();
        size_type nodeb_index = b.index();
        
        set<size_type> node_pair;  // create set with node a and b indicies
        node_pair.insert(nodea_index);
        node_pair.insert(nodeb_index);
        

        if (has_edge(a, b) == true) {
            // return the current edge if it already exists
            return Edge(this, nodea_index, nodeb_index, edges_map_[node_pair] );
        } else {
            size_type edge_index = edges_map_.size();
            edges_map_[node_pair] = edge_index;  // add key, val pair for new edge
            edges_map_flipped_[edge_index] = node_pair; // add flipped info to copy map

            if (adj_list.count(nodea_index) == 0)
            {
               vector<size_type> v;
               v.push_back(nodeb_index);
               adj_list[nodea_index] = v;
            }

            else
            	adj_list[nodea_index].push_back(nodeb_index);


            return Edge(this, nodea_index, nodeb_index, edge_index);
        }
    }
    
    
        /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        nodes_ = {};
        edges_map_.clear();
        edges_map_flipped_.clear();
    }



  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() 
    {
    }
    
    // HW1 #2: YsOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
     Node operator*() const
     {
         //cout<< itr->index <<endl;
         //return Node(graph_, itr->index);
     	//cout<<itr->position[0]<<endl;
         return graph_->node(itr->index);
     }
     NodeIterator& operator++()
    {
    	/*
    	for(auto iter = graph_->adj_list.cbegin(); iter != graph_->adj_list.cend(); ++iter)
       { 
            cout<< "key: " << iter-> first  << endl;
            for (size_type i = 0; i < (iter->second).size(); i++)
            cout<< " value: " << iter->second[i] << " "<< endl;    
       }
       */

         itr++;
         return *this;
         
    }
     bool operator==(const NodeIterator& i) const
     {
        // cout<< to_string(itr == i.itr) <<endl;
         return itr == i.itr;
     }

   private:
  
    friend class Graph;
    Graph* graph_;
    typename vector<internal_node>::const_iterator itr;
    
    /** Private Constructor for NodeIterator */
    NodeIterator(const Graph* graph,typename vector<internal_node>::const_iterator itr_)
    : graph_(const_cast<Graph*>(graph)), itr(itr_) {}
    

    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
   node_iterator node_begin() const
   {
       return node_iterator(this, nodes_.begin());
   }
   node_iterator node_end() const
   {
       return node_iterator(this, nodes_.end());
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
    IncidentIterator() {}

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
     
     Edge operator*() const
     {
     	set<size_type> pair_ab = {target_node_idx, *i_itr};
        map<set<size_type>, size_type>::iterator iter = graph_->edges_map_.find(pair_ab);
     	return graph_->edge(iter->second);//index of edge
     }
     IncidentIterator& operator++()
     {
     	 i_itr++;
         return *this;
     }
     bool operator==(const IncidentIterator& iit) const
     {
     	return i_itr == iit.i_itr;
     }
     

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE]
    Graph* graph_;
    vector<size_type>::iterator i_itr;
    size_type target_node_idx;
    

    /** Private Constructor for IncidnentIterator */
    IncidentIterator(const Graph* graph, vector<size_type>::iterator i_itr_)
    : graph_(const_cast<Graph*>(graph)), i_itr(i_itr_) {}

  

     /** Another private Constructor for IncidentIterator */
    IncidentIterator(const Graph* graph, vector<size_type>::iterator i_itr_, size_type idx_)
    : graph_(const_cast<Graph*>(graph)), i_itr(i_itr_), target_node_idx(idx_) {}

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>{
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
     Edge operator*() const
     {
          return graph_->edge(e_itr->second);
     }

     EdgeIterator& operator++()
     {
         e_itr++;
         return *this;
     }
     bool operator==(const EdgeIterator& iit) const
     {
       return e_itr == iit.e_itr;
     }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    map<set<size_type>, size_type>::const_iterator e_itr;

 /** Private Constructor for IncidnentIterator */
    EdgeIterator(const Graph* graph, map<set<size_type>, size_type>::const_iterator e_itr_)
    : graph_(const_cast<Graph*>(graph)), e_itr(e_itr_) {}

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
   edge_iterator edge_begin() const 
   {
       return edge_iterator(this, edges_map_.begin());
   }
   edge_iterator edge_end() const
   {
      return edge_iterator(this, edges_map_.end());
   }

 private:
  // Structure to store a pair of node indices
  
    // internal type for graph nodes
    struct internal_node {
        size_type index;
        Point position;
        node_value_type value;
    };
    
    vector<internal_node> nodes_;
    //edges_map_ contains: key = pair of node a and b indicies, val = unique edge index
    map<set<size_type>, size_type> edges_map_;
    map<size_type, set<size_type>> edges_map_flipped_;  // keys and vals flipped here
    map<size_type, vector<size_type>> adj_list; //adjacency list 
};

#endif // CME212_GRAPH_HPP