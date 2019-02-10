#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

//--documentation_0
//--Documentation is sufficient
//--END

//--style_0
//--Looks good
//--END

//--functionality_0
//--Code passed all unit tests, nice work!
//--END



/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph {
private:
 
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
    Graph() : nodes_(), edges_(), neighbors(){
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
    class Node : private totally_ordered<Node>{
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
        }

        
        /** Return this node's position. */
        const Point& position() const {
            return graph_->nodes_[uid_].first;
        }
        
        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
            return uid_;
        }
        
        // HW1: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        // node_value_type& value();
        // const node_value_type& value() const;
        // size_type degree() const;
        // incident_iterator edge_begin() const;
        // incident_iterator edge_end() const;
        
        /** Return this node's value.
         * @post result == value of the node.
         */
        node_value_type& value(){
            return graph_->nodes_[uid_].second;
        }
        
        /** Return this node's value, without any change
         * @post result == value of the node.
         */
        const node_value_type& value() const{
            return graph_->nodes_[uid_].second;
        }
        
        /** Return this node's degree
         * @post result == # of this node's neighbors.
         */
        size_type degree() const{
            return graph_->neighbors[uid_].size();
        }
        
        /** Return the first incident iterator*/
        incident_iterator edge_begin() const{
            return IncidentIterator(graph_, uid_, 0);
        }
        
        /** Return the last incident iterator*/
        incident_iterator edge_end() const{
            return IncidentIterator(graph_, uid_, this->degree());
        }
        
        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        bool operator==(const Node& n) const {
            return ((graph_ == n.graph_)&&(uid_ ==  n.uid_) );
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
            if ((graph_ == n.graph_)&&(uid_ <  n.uid_)){
                return true;
            }
            if (graph_ < n.graph_){
                return true;
            }
            return false;
        }
        
    private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;
        Graph* graph_;
        size_type uid_;
        Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
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
        return nodes_.size();
    }
    
    /** Synonym for size(). */
    size_type num_nodes() const {
        return size();
    }
    
    /** Add a node to the graph, returning the added node.
     * @param[in] position The new node's position
     * @param[in] value    The new node's value
     * @post new num_nodes() == old num_nodes() + 1
     * @post result_node.index() == old num_nodes()
     *
     * Complexity: O(1) amortized operations.
     */
    Node add_node( const Point& position, const node_value_type& value = node_value_type() ) {
        nodes_.push_back(std::make_pair(position, value));
        neighbors.push_back(std::vector<std::pair<size_type, size_type>> ());
        return Node(this, num_nodes()-1);
    }
    

    
    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
        return ((n.graph_ == this) && (n.uid_<num_nodes()));
    }
    
    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
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
            return Node(graph_, uid_1);
        }
        
        /** Return the other node of this Edge */
        Node node2() const {
            return Node(graph_, uid_2);
        }
        
        
        
        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge& e) const {
            if (graph_ == e.graph_){
                if ((node1() == e.node1() && node2() == e.node2()) ||
                    (node1() == e.node2() && node2() == e.node1())){
                    return true;
                }
            }
            return false;
        }
        
        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge& e) const {
            if (graph_ == e.graph_){
                if (uid_ < e.uid_)
                    return true;
            }
            else{
                if (graph_ < e.graph_)
                    return true;
            }
            return false;
        }
        
    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;
        Graph* graph_;
        size_type uid_;
        size_type uid_1;
        size_type uid_2;
        Edge(const Graph* graph, size_type uid, size_type uid1, size_type uid2)
        : graph_(const_cast<Graph*>(graph)), uid_(uid), uid_1(uid1), uid_2(uid2){
        }
        // Use this space to declare private data members and methods for Edge
        // that will not be visible to users, but may be useful within Graph.
        // i.e. Graph needs a way to construct valid Edge objects
    };
    
    /** Return the total number of edges in the graph.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    size_type num_edges() const {
        return edges_.size();
    }
    
    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
        return Edge(this, i, edges_[i].first, edges_[i].second);
    }
    
    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node& a, const Node& b) const {
        for (unsigned int i=0; i<neighbors[a.uid_].size(); i++)
            if (b.uid_ == neighbors[a.uid_][i].first)
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
    Edge add_edge(const Node& a, const Node& b) {
        for (unsigned int i=0; i < neighbors[a.uid_].size(); i++){
            if (b.uid_ == neighbors[a.uid_][i].first){
                return Edge(this, neighbors[a.uid_][i].second, a.uid_, b.uid_);}
        }
        if (a<b){
            edges_.push_back(std::make_pair(a.uid_,b.uid_));
        }else{
            edges_.push_back(std::make_pair(b.uid_,a.uid_));
        }
        neighbors[a.uid_].push_back(std::make_pair(b.uid_,edges_.size()-1));
        neighbors[b.uid_].push_back(std::make_pair(a.uid_,edges_.size()-1));
        return Edge(this,edges_.size()-1,a.uid_,b.uid_);
    }
    
    /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        nodes_.clear();
        edges_.clear();
        neighbors.clear();
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
        NodeIterator() {
        }
        
        // HW1 #2: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        // Node operator*() const
        // NodeIterator& operator++()
        // bool operator==(const NodeIterator&) const
        
        /** Return the value of the iterator
         * @post result = node n with n.uid_ == index_ and n.graph_ == graph_.
         * Complexity: O(1).
         */
        Node operator*() const{
            return graph_->node(index_);
        }
        
        /** Return the next iterator
         * @post result == index_ + 1.
         */
        NodeIterator& operator++(){
            ++index_;
            return *this;
        }
        
        /** Test equivalence '=='
        * Two node_iterators are equal if both graphs and indices match
        */
        bool operator==(const NodeIterator& ni) const{
            return ((graph_ == ni.graph_ )&&(index_ == ni.index_));
        }
    private:
        friend class Graph;
        // HW1 #2: YOUR CODE HERE
        Graph* graph_;
        size_type index_;
        NodeIterator(const Graph* graph, size_type index)
        : graph_(const_cast<Graph*>(graph)), index_(index){
        }
    };
    
    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_iterator node_begin() const
    // node_iterator node_end() const
    
    /** Return the first node iterator*/
    node_iterator node_begin() const{
        return NodeIterator(this, 0);
    }
    
    /** Return the last node iterator*/
    node_iterator node_end() const{
        return NodeIterator(this, this->num_nodes());
    }
    
    //
    // Incident Iterator
    //
    
    /** @class Graph::IncidentIterator
     * @brief Iterator class for edges incident to a node. A forward iterator. */
    class IncidentIterator : private totally_ordered<IncidentIterator>{
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
        // Edge operator*() const
        // IncidentIterator& operator++()
        // bool operator==(const IncidentIterator&) const
        
        
        /** Return the edge that the iterator pointing to.
         * @post result = node n with n.uid_ == index_ and n.graph_ == graph_.
         * Complexity: O(1).
         */
        Edge operator*() const{
            return Edge(graph_, graph_->neighbors[central_id_][index_].second, central_id_, graph_->neighbors[central_id_][index_].first);
        }
        
        /** Return the next iterator
         * @post result == index_ + 1.
         */
        IncidentIterator& operator++(){
            ++index_;
            return *this;
        }
        
        /** Test equivalence '=='
         * Two incident_iterators are equal if graphs, entral nodes, and indices of iterators match
         */
        bool operator==(const IncidentIterator& ii) const{
            return ((graph_ == ii.graph_)&&(central_id_ == ii.central_id_)&&(index_ == ii.index_));
        }
    private:
        friend class Graph;
        // HW1 #3: YOUR CODE HERE
        Graph* graph_;
        size_type central_id_; // The id of the central node
        size_type index_; // Iterating
        IncidentIterator(const Graph* graph, size_type central_id, size_type index)
        : graph_(const_cast<Graph*>(graph)), central_id_(central_id), index_(index){
        }
        
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
        // Supply definitions AND SPECIFICATIONS for:
        // Edge operator*() const
        // EdgeIterator& operator++()
        // bool operator==(const EdgeIterator&) const
        Edge operator*() const{
            return graph_->edge(index_);
        }
        
        /** Return the next iterator
         * @post result == index_ + 1.
         */
        EdgeIterator& operator++(){
            ++index_;
            return *this;
        }
        /** Test equivalence '=='
         * Two edge_iterators are equal if both graphs and indices match
         */
        bool operator==(const EdgeIterator& ei) const{
            return ((graph_ == ei.graph_)&&(index_ == ei.index_));
        }
    private:
        friend class Graph;
        // HW1 #5: YOUR CODE HERE
        Graph* graph_;
        size_type index_;
        EdgeIterator(const Graph* graph, size_type index) : graph_(const_cast<Graph*>(graph)), index_(index){
        }
    };
    
    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // edge_iterator edge_begin() const
    // edge_iterator edge_end() const
    
    /** Return the first edge iterator*/
    edge_iterator edge_begin() const{
        return EdgeIterator(this, 0);
    }
    /** Return the last edge iterator*/
    edge_iterator edge_end() const{
        return EdgeIterator(this, this->num_edges());
    }
private:
    
    //--design_0
    //--This storage makes sense for now, but one thing to think about over the next few weeks
    //--is whether you really need a pointer back to the edge_uid in your adjacency list
    //--(The nested structure of vector<vector<pair>> is a bit complicated)
    //--START
    
    
    // A vector of nodes with their points and values stored
    std::vector<std::pair<Point, node_value_type>> nodes_;
    // A vector of edges, with each edge's both end nodes' identities stored
    std::vector<std::pair<size_type, size_type>> edges_;
    // This stores each node's neighbor nodes, which means that they have an
    // edge exist in between. The second vector, neighbor nodes vector, has pair
    // struture, <id of neighbor node, id of the connected edge>
    std::vector<std::vector<std::pair<size_type, size_type>>> neighbors;
    //--END
};

#endif // CME212_GRAPH_HPP
