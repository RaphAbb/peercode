#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */
#include <algorithm>
#include <vector>
#include <unordered_map> 
#include <cassert>
#include <utility>      // std::pair, std::make_pair

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

    // predeclaring the internal node struct
    struct InternalNode;

 public:

    //
    // PUBLIC TYPE DEFINITIONS
    //

    /** Type of this graph. */
    using graph_type = Graph;

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

    //
    // CONSTRUCTORS AND DESTRUCTOR
    //

    /** Construct an empty graph. 
    * @return empty graph class
    */
    Graph()
        : intnodes_(), size_(0) {
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
     * Node objects only store unique id and pointer to the graph. Internal nodes
     * store actual node data (point and value).
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
        // empty constructor for invalid node
        Node() {
        
        }

        /** Return this node's position. 
        * @return node position
        */
        const Point& position() const {
            // grab the id'th index of the vector of internal nodes
            // and return the point value of that internal node
            return graph_->intnodes_[uid_node_].point_int_node;
        }

        /** Return this node's index, a number in the range [0, graph_size). 
        * @return unique node id / index
        */
        size_type index() const {
            // return the ID of that node
            return uid_node_;
        }

        /** Return this node's value. User can also change the value.
        * @return node value
        */
        node_value_type& value() {
            // go into the vector of internal nodes at the id of this node,
            // then grab the value of that internal node
            return graph_->intnodes_[uid_node_].val_int_node;
        }

        /** Return this node's value.
        * @return node value
        */
            const node_value_type& value() const {
            // go into the vector of internal nodes at the id of this node,
            // then grab the value of that internal node
            return graph_->intnodes_[uid_node_].val_int_node;
        }

        /** Return the number of edges connected to this node object.
        * @return number of edges connected to node.
        */        
        size_type degree() const {
            // return the size of the nested unordered map which is the number
            // of edges that are connected to it
            return graph_->umap_bynodeid[uid_node_].size();
        }

        /** Return an incident iterator to this node's first edge.
        * @return incident iterator of node.
        */ 
        incident_iterator edge_begin() const { 
            /* create a new incident iterator passing in the graph, the id of
            * this node, of this node, and an iterator to the end of
            * the NESTED unordered map map which is located in the parent umap
            * where the key is the id of this node
            */
            return IncidentIterator(graph_, uid_node_, 
                graph_->umap_bynodeid.at(uid_node_).begin());
        }

        /** Return an incident iterator to this node's last edge.
        * @return incident iterator of node.
        */         
        incident_iterator edge_end() const {
            /* create a new incident iterator passing in the graph, the id of
            * this node, of this node, and an iterator to the start of
            * the NESTED unordered map map which is located in the parent umap
            * where the key is the id of this node
            */
            return IncidentIterator(graph_, uid_node_, 
                graph_->umap_bynodeid.at(uid_node_).end());
        }

        /** Test whether this node and @a n are equal. Equal nodes have the same
         *  graph and the same index.
         *
         * @param[in] n    Node object
         * @return         True if equal, false if not equal
         */
        bool operator==(const Node& n) const {
            // check that the index of this node and the index of the input
            // node are the same, and that they refer to the same graph
            return (uid_node_ == n.uid_node_ and graph_ == n.graph_);
        }

        /** Test whether this node is less than @a n in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any geometric meaning.
         *
         * The node ordering relation must obey trichotomy: For any two nodes x
         * and y, exactly one of x == y, x < y, and y < x is true.
         *
         * @param[in] n    Node object
         * @return         True if this node less than @a n, false if not equal
         */
        bool operator<(const Node& n) const {
            // check that the two indices are the same and that the graph
            // is the same
            return (index() < n.index() and graph_ == n.graph_);
        }

     private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;

        /** Private Constructor of node that takes in graph and id
         *
         * @param[in] g    Graph object
         * @param[in] i    Unique identifier of node
         * @post           Assign graph_ and uid_node_ members to inputs
         */
        Node(const Graph* g, size_type i)
                : graph_(const_cast<Graph*>(g)), uid_node_(i) {
        }

        /** Pointer to the graph that this node belongs to */
        Graph* graph_;

        /** Unique identifier for this node */
        size_type uid_node_;
    };

    /** Return the number of nodes in the graph.
     *
     * @return    Total number of nodes in the graph.
     * Complexity: O(1).
     */
    size_type size() const {
        // return the size of the internal node vector
        return intnodes_.size();
    }

    /** Synonym for size(). */
    size_type num_nodes() const {
        return size();
    }

    /** Add a node to the graph, returning the added node.
     * @param[in] position   The new node's position
     * @param[in] value      The new node's value
     *
     * @post new num_nodes() == old num_nodes() + 1
     * @post result_node.value() == value
     * @post result_node.index() == old num_nodes()
     *
     * @return new node object.
     *
     * Complexity: O(1) amortized operations.
     */
    Node add_node(const Point& position, const node_value_type& val = 
        node_value_type()) {

        // push back a new node to the vector, with the 
        intnodes_.push_back({position, val, size_});

        // increment the size counter
        size_++;

        // return a new node object with the ID of the size - 1
        return Node(this, size_-1);
    }

    /** Determine if a Node belongs to this Graph
     * @param[in] n    Node object
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
        // if the index of the passed node is between 0 and 
        // the size of the internal node vector, it is a valid query
        return (n.uid_node_ < size_);
    }

    /** Return the node with index @a i.
     * @param[in] i Index/unique id of node
     *
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * @return New node object with index @a i
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
        // make sure that the index is lower than the size of the
        // vector, and if so, return that node.
        assert(i < size());
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

        /** Return a node of this Edge 
        * @return node object with ID of node 1 member
        */
        Node node1() const {
            // return new edge object, passing in node1 id
            return Node(graph_, node1_id);
        }

        /** Return a node of this Edge 
        * @return node object with ID of node 2 member
        */        
        Node node2() const {
            // return new edge object, passing in node2 id
            return Node(graph_, node2_id);
        }

        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         *
         * @param[in] e    Edge object
         * @return         True if @a e id and graph is equal to this edge's
         *                  id and graph
         */
        bool operator==(const Edge& e) const {
            // check that the ids are the same and that the graphs are the same
            return (uid_edge_ == e.uid_edge_ and graph_ == e.graph_);
        }

        /** Test whether this node is less than @a n in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any geometric meaning.
         *
         * The edge ordering relation must obey trichotomy: For any two nodes x
         * and y, exactly one of x == y, x < y, and y < x is true.
         *
         * @param[in] e    Edge object
         * @return         True if this edge less than @a e, false if not equal
         */
        bool operator<(const Edge& e) const {
            // check that the id is less than input and that the graphs are equal
            return (uid_edge_ < e.uid_edge_ and graph_ == e.graph_);
        }

     private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;

        /** Private Constructor of edge that takes in graph, id, and 2 node ids
         *
         * @param[in] g    Graph object
         * @param[in] i    Unique identifier of edge
         * @param[in] n1   Unique identifier of first node
         * @param[in] n2   Unique identifier of second node 
         * @post           Assign graph_ and uid_node_ members to inputs
         */
        Edge(const Graph* g, size_type i, size_type n1_id, size_type n2_id)
                : graph_(const_cast<Graph*>(g)), uid_edge_(i), 
                node1_id(n1_id), node2_id(n2_id) {
        }

        /** Pointer to the graph that this edge belongs to */
        const Graph* graph_;

        /** Unique identifier for this edge */
        size_type uid_edge_;

        /** Unique identifier for the first node of this edge */
        size_type node1_id;

        /** Unique identifier for the second node of this edge */
        size_type node2_id;

    };

    /** Return the total number of edges in the graph.
     *
     * @return    Total number of edges in the graph.
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    size_type num_edges() const {
        // return the size of the umap
        return umap_byedgeid.size();
    }

    /** Return the edge with index @a i.
     * @param[in] i Index/unique id of edge
     *
     * @pre 0 <= @a i < num_edges()
     * @post result_edge.index() = i
     *
     * @return New edge object with index @a i
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
        // make sure that the value is within the number of edges
        assert(i < num_edges());

        // instantiate new edge object with this graph and with the index i
        return Edge(this, i, umap_byedgeid.at(i).first, 
            umap_byedgeid.at(i).second);
    }

    /** Test whether two nodes are connected by an edge.
     *
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node& a, const Node& b) const {
        // Check if the key of node a's ID has an associated ID for node b
        // in the nested umap
        return (umap_bynodeid.count(a.uid_node_) > 0 
            and umap_bynodeid.at(a.uid_node_).count(b.uid_node_) > 0);
    }

    /** Add an edge to the graph, or return the current edge if it already exists.
     *
     * @param[in] n1    First node object
     * @param[in] n2    Second node object
     *
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

        // save the ID numbers for each node
        size_type a_id = a.uid_node_;
        size_type b_id = b.uid_node_;

        // check if the edge already exists; if so, return that edge
        if (has_edge(a, b) == true) {
            // create new edge object, passing in the graph, the edge id,
            // the node1 id and the node2 id
            return Edge(this, umap_bynodeid.at(a_id).at(b_id), a_id, b_id);
        }
        // if the edge doesn't exist, add it to the underlying containers and
        // return that new edge
        else {
            // calculate the number of existing edges
            size_type edge_id = num_edges();

            // add to umap of edges by id
            umap_byedgeid[edge_id] = std::make_pair(a_id, b_id);

            // add to umap of edges by node id
            umap_bynodeid[b_id][a_id] = edge_id;
            umap_bynodeid[a_id][b_id] = edge_id;

            return Edge(this, edge_id, a_id, b_id);
        }
    }

    /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     * @post clear containers (vector and 2 unordered maps)
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        // clearing all containers
        intnodes_.clear();
        umap_byedgeid.clear();
        umap_bynodeid.clear(); 
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

        /** Dereference the NodeIterator. 
         * @return node object
         */
        Node operator*() const {
            // dereference the internal node iterator, call the id of the 
            // internal node, then pass the graph and that id to make a new
            // node and return
            return Node(graph_, (*node_it).uid_int_node);
        }
        
        /** Increment the NodeIterator. 
         * @return node iterator object
         */
        NodeIterator& operator++() {
            // incrememnt the underlying iterator, return this
            ++node_it;
            return *this;
        }

        /** Test whether this node iterator and @a node iterator are equal.
         *
         * @param[in] ni2    Node iterator object
         * @return    True if @a ni2 equals this in underlying iterator and graph
         *
         * Equal node iterators have the same graph and the same underlying
         * iterator.
         */
        bool operator==(const NodeIterator& node_it_2) const {
            return (node_it == node_it_2.node_it and graph_ == node_it_2.graph_);
        }

     private:
        friend class Graph;

        /** Private Constructor of node iterator that takes in graph and 
         * an iterator of an underlying container of graph
         *
         * @param[in] g    Graph object
         * @param[in] it   Underlying iterator of internal node container
         * @post           Assign graph_ and node_it members to inputs
         */
        NodeIterator(const Graph* g, typename std::vector<InternalNode>::const_iterator i)
                : graph_(const_cast<Graph*>(g)), node_it(i) {
        }

        /** Pointer to the graph that this node iterator belongs to */
        Graph* graph_;

        /** Iterator of data structure in the underlying graph */
        typename std::vector<InternalNode>::const_iterator node_it;
    };

    /** Return a node iterator to this graph's first associated node 
     * @return node iterator object
     */ 
    node_iterator node_begin() const {
        // return a new node iterator belonging to this graph, and setting the 
        // iterator to the beginning of the vector of internal nodes
        return NodeIterator(this, intnodes_.begin());
    }
    
    /** Return a node iterator to this graph's last associated node 
     * @return node iterator object
     */
    node_iterator node_end() const {
        // return a new node iterator belonging to this graph, and setting the
        // iterator to the end of the vector of internal nodes
        return NodeIterator(this, intnodes_.end());
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

        /** Dereference the IncidentIterator
        * @return an edge object
        */
        Edge operator*() const {
            // create new edge passing in the graph, the edge id, the spawning node id,
            // and the other node id
            return Edge(graph_, umap_it->second, node_id, umap_it->first);
        }

        /** Increment the IncidentIterator.
        * @post Increment the underlying iterator of unordered map 
        * @return an incident iterator object
        */
        IncidentIterator& operator++() {
            // incrememnt the underlying iterator, then dereference and return this
            ++umap_it;
            return *this;
        }

        /** Test whether this incident iterator and @a incident iterator are equal.
         *
         * Equal incident iterators have the same graph and the same underlying
         * iterator.
         *
         * @param[in] ii  Incident iterator object
         * @return        True if this incident iter less than @a ii, false if not equal
         */
        bool operator==(const IncidentIterator& inc_it_2) const {
            // check that the umap iterators of the incident iterators are the same
            return (umap_it == inc_it_2.umap_it and graph_ == inc_it_2.graph_);
        }

     private:
        friend class Graph;
 
        /** Private Constructor of incident iterator that takes in graph, a node,  
         * and an iterator of an underlying container of graph
         *
         * @param[in] g    Graph object
         * @param[in] i    Unique identifier of node
         * @param[in] it   Underlying iterator specific to node with id @a i
         * @post           Assign graph_ and uid_node_ members to inputs
         */
        IncidentIterator(const Graph* g, size_type n_id, typename 
            std::unordered_map<size_type, size_type>::const_iterator i)
                : graph_(const_cast<Graph*>(g)), node_id(n_id), umap_it(i) {
        }

        /** Pointer to the graph that this incident iterator belongs to */
        const Graph* graph_;

        /** Pointer to the node that this incident iterator belongs to */
        size_type node_id;

        /** Iterator of data structure in the underlying graph */
        typename std::unordered_map<size_type, size_type>::const_iterator umap_it;
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

        /** Dereference the EdgeIterator
        * @return an edge object
        */        
        Edge operator*() const {
            // return a new Edge object and pass in
            // 1) the graph, 2) the first element of the underlying iterator, which is the edge
            // id, 3) the first element OF the second element, which is the node 1 id, and 
            // 4) the second element OF the second elemnt, which is the node 2 id
            return Edge(graph_, edge_it->first, edge_it->second.first, edge_it->second.second);
        }

        /** Increment the EdgeIterator. 
        * @return an edge iterator object
        */
        EdgeIterator& operator++() {
            // increment the underlying iterator and return this
            ++edge_it;
            return *this;
        }

        /** Test whether this edge iterator and @a edge iterator are equal.
         *
         * Equal edge iterators have the same graph and the same underlying
         * iterator.
         *
         * @param[in] ei  Edge iterator object
         * @return        True if this edge iter less than @a ei, false if not equal
         */
        bool operator==(const EdgeIterator& edge_it_2) const {
            // check that the umap iterators of the incident iterators are the same, as well as
            // the graphs of both objects
            return (edge_it == edge_it_2.edge_it and graph_ == edge_it_2.graph_);
        }

     private:
        friend class Graph;
        /** Private Constructor of edge iterator that takes in graph and   
         * and an iterator of an underlying container of graph in which they keys
         * are the edge ids and the values are pairs of node ids
         *
         * @param[in] g    Graph object
         * @param[in] i    Unique identifier of edge
         * @param[in] it   Underlying iterator of unordered map where keys are edge id
         * @post           Assign graph_ and uid_node_ members to inputs
         */
        EdgeIterator(const Graph* g, typename 
            std::unordered_map<size_type, std::pair<size_type, size_type>>::const_iterator i)
                : graph_(const_cast<Graph*>(g)), edge_it(i) {
        }

        /** Pointer to the graph that this edge iterator belongs to */
        const Graph* graph_;

        /** Iterator of data structure in the underlying graph */
        typename std::unordered_map<size_type, std::pair<size_type, size_type>>::const_iterator edge_it;

    };

    /** Return an edge iterator to this graph's first associated edge 
     * @return iterator of all edge.
     */
    edge_iterator edge_begin() const { 
        return EdgeIterator(this, umap_byedgeid.begin());
    }

    /** Return an edge iterator to this graph's last associated edge 
     * @return iterator of all edge.
     */
    edge_iterator edge_end() const {
        return EdgeIterator(this, umap_byedgeid.end());
    }

 private:
    /** @struct Graph::InternalNode
     * @brief struct representing the graph's internal nodes.
     *
     * Internal node objects are used to store information about the Graph's nodes.
     * Information includes a 3-d point, a value, and a unique id.
     */
    struct InternalNode {
        Point point_int_node;   // The text held by an internal node
        node_value_type val_int_node;   // The value held by an internal node
        size_type uid_int_node;      // The unique identification by an internal node
    };

    /** Vector of internal node objects */
    std::vector<InternalNode> intnodes_;

    /** Unordered map 1 to store edges
    * Key is the edge id, values are a pair with two node ids
    */
    std::unordered_map<size_type, std::pair<size_type, size_type>> umap_byedgeid;

    /** Unordered map 2 to store edges
    * Key is the node id, values are a nested unordered map, with each node that
    * the key node is connected to, and the edge id of the edge that connects them
    */
    std::unordered_map<size_type, std::unordered_map<size_type, 
        size_type>> umap_bynodeid;

    /** Number of nodes in the graph */
    size_type size_;

};

#endif // CME212_GRAPH_HPP