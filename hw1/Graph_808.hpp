#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template<typename V>
class Graph: private totally_ordered<Graph<V>> {
private:

public:
    //
    // PUBLIC TYPE DEFINITIONS
    //
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
    using size_type = unsigned int;
    
    //
    // CONSTRUCTORS AND DESTRUCTOR
    //
    
    /** Construct an empty graph. */
    Graph(){
        // construct an empty graph
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
        
        /** Constructing an empty node. */
        Node() {
            
        }
        
        
        /** Return this node's position. */
        const Point& position() const {
            return graph_->points_.at(idx_);
        }
        
        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
            return this->idx_;
        }
        
        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        // HW1: YOUR CODE HERE

        // node_value_type& value();
        /** return node the reference to the node value */
        node_value_type& value() {
            return const_cast<graph_type*>(graph_)->values_.at(idx_);
        };
        
        // const node_value_type& value() const;
        /* return the node value, as a default constant */
        const node_value_type& value() const{
            return node_value_type();
        }
        
        // size_type degree() const;
        /** Return num of edges to a node */
        size_type degree() const{
            return graph_->adj_list_[idx_].size();
        }
        
        // incident_iterator edge_begin() const;
        /** Return the first incident iterator when looping over all edges around the node */
        incident_iterator edge_begin() const{
            return IncidentIterator(const_cast<graph_type*>(graph_), idx_, 0);
        }
        
        // incident_iterator edge_end() const;
        /** Return the final incident iterator when looping over all edges around the node */
        incident_iterator edge_end() const{
            size_type adj_idx_end = graph_->adj_list_[idx_].size();
            return IncidentIterator(const_cast<graph_type*>(graph_), idx_, adj_idx_end);
        }
        
        
        bool operator==(const Node& n) const {
            if(n.graph_ == this->graph_ and n.index() == this->index()){
                return true;
            }
            else{
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
            return index() < n.index();
        }
        

    private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;
        friend class IncidentIterator;
        
        const graph_type* graph_;
        size_type idx_;
        
        //constructor of node
        Node(graph_type* graph, size_type idx)
        : graph_(graph), idx_(idx){}
    };
    
    /** Return the number of nodes in the graph.
     *
     * Complexity: O(1).
     */
    size_type size() const {
        return points_.size();
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
    Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
        points_.push_back(position);     // Add node's point to graph
        values_.push_back(value);
        next_idx_++;                     // Increment index
        adj_list_.push_back(std::vector<InternalEdge>());
        
        // debug
        //std::cout << "debug addnote" << std::endl;
        //std::cout << value << std::endl;
        //node_value_ = node_value_type();
        
        return Node(this, next_idx_- 1); // Create & return node instance
    }
    
    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
        return this == n.graph_;
    }
    
    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
        return Node(const_cast<graph_type*>(this), i);
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
        /** Construct an invalid Edge. */
        Edge() {
        }
        
        /** Return a node of this Edge */
        Node node1() const {
            return Node(graph_, node1_);
        }
        
        /** Return the other node of this Edge */
        Node node2() const {
            return Node(graph_, node2_);
        }
        
        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge& e) const {
            if(e.node1_ == node1_ and e.node2_ == node2_){
                return true;
            }
            else if(e.node1_ == node2_ and e.node2_ == node1_){
                return true;
            }
            else{
                return false;
            }
        }
        
        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge& e) const {
            return idx_ < e.idx_;
        }
        
    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;
        
        size_type node1_;
        size_type node2_;
        size_type idx_;
        graph_type* graph_;
        
        /** Construct valid edge object */
        Edge(size_type node1, size_type node2, size_type idx, graph_type* graph)
        : node1_(node1), node2_(node2), idx_(idx), graph_(graph) {}
    };
    
    /** Return the total number of edges in the graph.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    size_type num_edges() const {
        return next_edge_idx_;
    }
    
    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
        if(i < 0 || i > next_edge_idx_-1){
            throw std::out_of_range("Edge index out of range when calling edge.");
        }
        return Edge(edges_.at(i).node1_, edges_.at(i).node2_, i, const_cast<graph_type*>(this));
    }
    
    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node& a, const Node& b) const {
        // Check that nodes belong to this graph
        if((has_node(a) and has_node(b)) == false){
            throw std::runtime_error("Node(s) do not belong to graph in has_edge()");
        }
        // Loop through edges that indclude a and look for (a,b) edge
        for(auto edge : adj_list_.at(a.index())){
            if (edge.node1_ == a.index() and edge.node2_ == b.index()){
                return true;
            }
            else if (edge.node1_ == b.index() and edge.node2_ == a.index()){
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
    Edge add_edge(const Node& a, const Node& b) {
        // Check that a and b are distinct valid nodes of this graph
        if(a.graph_ != this || b.graph_ != this){
            throw std::runtime_error("One of the nodes does not belong to the graph.");
        }
        if(a == b){
            throw std::runtime_error("add_edge() received non-distinct nodes. Self-loops are frowned upon.");
        }
        
        // If the edge exists loop over edges for node a, find index, return edge
        if (has_edge(a, b)){
            size_type index;
            for(auto edge_ : adj_list_.at(a.index())){
                if(edge_.node1_ == a.index() and edge_.node2_ == b.index()){
                    index = edge_.idx_;
                    break;
                }
                else if (edge_.node1_ == b.index() and edge_.node2_ == a.index()){
                    index = edge_.idx_;
                    break;
                }
            }
            return edge(index);
        }
        // Else create new edge and store information in appropriate containers
        else{
            next_edge_idx_++;
            InternalEdge edge_ab(a.index(), b.index(), next_edge_idx_-1);
            adj_list_[a.index()].push_back(edge_ab);
            adj_list_[b.index()].push_back(edge_ab);
            edges_.push_back(NodePair(a.index(), b.index()));
            return Edge(a.index(), b.index(), next_edge_idx_-1, this);
        }
    }
    
    /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 and num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        // Reset counters
        next_idx_      = 0;
        next_edge_idx_ = 0;
        // Clear graph containers
        points_.clear();
        values_.clear();
        adj_list_.clear();
        edges_.clear();
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
        
        /** Construct an invalid NodeIterator. */
        NodeIterator() {
        }
        
        /** NodeIterator Operator *() redefine, return the pointer to the node */
        Node operator*() const {
            //std::cout << "*" << std::endl;
            return graph_->node(idx_);
        }

        // NodeIterator& operator++()
        /** NodeIterator Operator ++ redefine */
        node_iterator& operator++() {
            idx_ = idx_ + 1;
            return *this;
        }
        
        // bool operator==(const NodeIterator&) const
        /** NodeIterator Operator == redefine */
        bool operator==(const NodeIterator& node_iterator) const {
            //std::cout << "==" << std::endl;
            if (node_iterator.graph_ == this->graph_ and node_iterator.idx_ == this->idx_){
                return true;
            } else {
                return false;
            }
        }
        
    private:
        friend class Graph;
        const graph_type* graph_;
        size_type idx_;
        
        // constructor for NodeIterator
        NodeIterator(graph_type* graph, size_type idx)
        : graph_(graph), idx_(idx){}
    };
    
    // HW1 #2: YOUR CODE HERE
    // node_iterator node_begin() const
    node_iterator node_begin() const{
        return NodeIterator(const_cast<graph_type*>(this), 0);
    }
    
    // node_iterator node_end() const
    node_iterator node_end() const{
        return NodeIterator(const_cast<graph_type*>(this), size());
    }
    //
    // Incident Iterator
    //
    
    /** @class Graph::IncidentIterator
     * @brief Iterator class for edges incident to a node. A forward iterator. */
    class IncidentIterator: private totally_ordered<IncidentIterator> {
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
        /** Return the edges the incident operator is currently referring to */
        Edge operator*() const{
            // find the adjacent node number
            size_type node_idx2_ {};
            if (node_idx_ != graph_->adj_list_[node_idx_][adj_idx_].node2_) {
                node_idx2_ = graph_->adj_list_[node_idx_][adj_idx_].node2_;
            } else {
                node_idx2_ = graph_->adj_list_[node_idx_][adj_idx_].node1_;
            }
            
            // find the edge idx
            size_type edge_idx_ = graph_->adj_list_[node_idx_][adj_idx_].idx_;
            // return edge
            //return Edge(node_idx_, node_idx2_, edge_idx_, const_cast<graph_type*>(graph_));
            return graph_->edge(edge_idx_);
        }
        
        // IncidentIterator& operator++()
        /** Define ++ operator */
        IncidentIterator& operator++() {
            adj_idx_ = adj_idx_ + 1;
            return *(this);
        }
        
        // bool operator==(const IncidentIterator&) const
        /** Define operator == */
        bool operator==(const incident_iterator& iit) const{
            if (graph_ == iit.graph_ and node_idx_ == iit.node_idx_ and adj_idx_ == iit.adj_idx_) {
                return true;
            }
            return false;
        }
        
    private:
        friend class Graph;
        friend class Edge;
        const graph_type* graph_;
        size_type node_idx_; // idx of the incident node
        size_type adj_idx_; // idx of the node the edge is in the adj_list,
      
        // private constructor
        IncidentIterator(graph_type* graph, size_type node_idx, size_type adj_idx): graph_(graph), node_idx_(node_idx), adj_idx_(adj_idx) {};
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
        
        /** Construct an invalid EdgeIterator. */
        EdgeIterator() {
        }
        
        // HW1 #5: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        // Edge operator*() const
         /** Define operator * for dereferencing */
        Edge operator*() const {
            return graph_->edge(edge_idx_);
        }
        
        // EdgeIterator& operator++()
         /** Define operator ++ */
        EdgeIterator& operator++() {
            edge_idx_ = edge_idx_ + 1;
            return *this;
        }
        
        // bool operator==(const EdgeIterator&) const
         /** Define operator == */
        bool operator==(const EdgeIterator& ei) const {
            if (graph_ == ei.graph_ and edge_idx_ == ei.edge_idx_) {
                return true;
            }
            return false;
        }
        
    private:
        friend class Graph;
        friend class Edge;
        const graph_type* graph_;
        size_type edge_idx_; // idx of the incident node
        
        // private constructor
        EdgeIterator(graph_type* graph, size_type edge_idx): graph_(graph), edge_idx_(edge_idx){};
    };
    
    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    /** call the first iterator for iterating through the edges */
    // edge_iterator edge_begin() const
    edge_iterator edge_begin() const{
        return EdgeIterator(const_cast<graph_type*>(this), 0);
    }
    
    /** call the last iterator for iterating through the edges */
    // edge_iterator edge_end() const
    edge_iterator edge_end() const{
        return EdgeIterator(const_cast<graph_type*>(this), next_edge_idx_);
    }
        
private:
    // Structure to store a pair of node indices
    // friend class EdgeIterator;
    struct NodePair{
        NodePair(size_type a, size_type b)
        : node1_(a), node2_(b) {
            
        }
        size_type node1_;
        size_type node2_;
    };
    
    // Structure to store the node indices and index of an edge
    struct InternalEdge{
        InternalEdge(size_type a, size_type b, size_type idx)
        : node1_(a), node2_(b), idx_(idx){
            
        }
        size_type node1_;
        size_type node2_;
        size_type idx_; // edge idx
    };
    
    std::vector<Point> points_;                       // Node idx -> Point
    std::vector<node_value_type> values_;              // Node idx -> values
    std::vector<NodePair> edges_;                     // Map node idxs -> edge idx
    std::vector<std::vector<InternalEdge>> adj_list_; // Map node idx -> Edges
    size_type next_idx_ = 0;
    size_type next_edge_idx_ = 0;
};

#endif // CME212_GRAPH_HPP

