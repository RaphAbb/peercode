#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP


/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>

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

private:

    // Pre declared the internal struct for the nodes
    struct internal_node;

    // Pre declared the internal struct for the edges
    struct internal_edge;

public:

    //
    // PUBLIC TYPE DEFINITIONS
    //

    /** Type of this graph. */
    using graph_type = Graph;

    /** Pre declaration of Node type. */
    class Node;
    /** Synonym for Node (following STL conventions). */
    using node_type = Node;

    /** Pre declaration of Edge type. */
    class Edge;
    /** Synonym for Edge (following STL conventions). */
    using edge_type = Edge;

    /** Type of indexes and sizes.
     * Return type of Graph::Node::index(), Graph::num_nodes(),
     * Graph::num_edges(), and argument type of Graph::node(size_type) */
    using size_type = unsigned;

    /** Pre declaration of the node iterator type */
    class node_iterator;

    /** Pre declaration of the edge iterator type */
    class edge_iterator;

    /** Pre declaration of the incident iterator type */
    class incident_iterator;


    //
    // CONSTRUCTORS AND DESTRUCTOR
    //

    /** Construct an empty graph. */
    Graph()
        : index_to_nodes_(),
          index_to_edges_(),
          edges_to_index_(),
          node_list_(),
          edge_list_(),
          size_(0),
          next_node_index_(0),
          num_edge_(0),
          next_edge_index_(0){
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
    using node_value_type = V;

    class Node : private totally_ordered<Node> {
    public:
        /** @brief Construct an invalid node.
         *
         * Valid nodes are obtained from the Graph class, but it
         * is occasionally useful to declare an @i invalid node, and
         * assign a valid node to it later. For example:
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
        Node()
            : graph_(), index_(0) {
        }

        /** @brief Construct a valid node.
         *
         * @param[in] graph A pointer to the graph representation.
         * @param[in] index the index of the node to create.
         *
         * @pre index >= 0
         * @post Create a valid of the graph with the given index.
         */
        Node(Graph<V> *graph, size_type index) {
            graph_ = graph;
            index_ = index;
        }

        /** @brief Return this node's position.
         *
         * The latest information about the node are stored
         * in a proxy design pattern. The node_fetch() method aims to look for
         * the corresponding node information in O(1).
         */
        const Point &position() const {
            return node_fetch().point;
        }

        /** @biref Return this node's index.
         *
         * @post The index of the node which is a a number in the
         * range [0, graph_size).
         */
        size_type index() const {
            return node_fetch().index;
        }

        /** @brief Test whether this node and @a n are equal.
         *
         * @pram[in] n A node to compare this note with.
         *
         * @return True if and only if the nodes have the same graph
         * and the same index.
         *
         * Equal nodes have the same graph and the same index.
         */
        bool operator==(const Node &n) const {
            bool equal = (graph_ == n.graph_ && index_ == n.index());
            return equal;
        }

        /** @brief Test whether this node is less than @a n in a global order.
         *
         * @parm[in] n A node to compare this note with.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any geometric meaning.
         *
         * The node ordering relation must obey trichotomy: For any two
         * nodes x and y, exactly one of x == y, x < y, and y < x is
         * true.
         */
        bool operator<(const Node &n) const {
            return (node_fetch().index < n.index_);
        }

        /** @brief Return the value of this node.
         *
         * @return The value fo this node.
         */
        node_value_type &value() {
            return node_fetch().value;
        }

        /** @brief Return the value of this node.
         *
         * @return A constant value of the value this node.
         */
        const node_value_type &value() const {
            const node_value_type& val = node_fetch().value;
            return val;
        }

        /** @brief Return the degree of this node.
         *
         * @return The number of incident edges i.e. the number of edges
         * with this node as first node.
         *
         * @post degree >= 0
         */
        size_type degree() const {
            size_type degree = node_fetch().degree;
            return degree;
        }

        /** @brief Return the first pointer to the list of edges incident to
         * this node.
         *
         * @pre The degree of this node is at least 1
         *
         * @return The first pointer to the list of edges incident to
         * this node.
         */
        incident_iterator edge_begin () const {
            return incident_iterator(&node_fetch().incident_edges[0]);
        }

        /** @brief Return the first pointer after the last element of the list
         * of edges incident to this node.
         *
         * @pre The degree of this node is at least 1
         *
         * @return The first pointer after the last element of the list
         * of edges incident to this node.
         */
        incident_iterator edge_end () const {
            return incident_iterator(&(node_fetch().incident_edges).back()+1);
        }

    private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;

        // Pointer back to the Graph container.
        Graph<V> *graph_;

        // The node unique identification number.
        size_type index_;

        /** @brief Private constructor to a valid node.
         *
         * @param[in] graph A pointer to the graph representation.
         * @param[in] index the index of the node to create.
         *
         * @pre index >= 0
         * @post Create a valid of the graph with the given index.
         */
        Node(const Graph<V> *graph, size_type index)
                : graph_(const_cast<Graph<V> *>(graph)), index_(index) {
        }

        /** @brief Helper method to return the appropriate node.
         * This finds the node with the correct index using an
         * unordered map.
         */
        internal_node &node_fetch() const {
            return graph_->index_to_nodes_[index_];
        }
    };


    /** @class Graph::node_iterator
     *
    * @brief Class representing the iterator objects for the graph nodes.
    *
    * node_iterator objects are used to loop on the nodes of the graph.
    */
    class node_iterator : private totally_ordered<node_iterator>{

    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type        = Node;                     // Element type
        using pointer           = Node*;                    // Pointers to elements
        using reference         = Node&;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

        /** @brief Default constructor of node_iterator. */
        node_iterator()
            : p{nullptr} {
        }

        /** @brief Overriding the operator*()
         *
         * @return The Node value of the pointer
         */
        Node operator *() const{
            return *p ;
        }

        /** @brief Overriding the operator==()
         *
         * @param[in] node_itr to compare this node iterator with
         *
         * @return True if the two node_iterator are equal
         */
        bool operator ==( const node_iterator &node_itr) const {
            return p == node_itr.p;
        }

        /** @brief Overriding the operator++()
         *
         * @return The next node_iterator
         */
        node_iterator & operator ++() {
            p++;
            return *this;
        }

    private :
        // The pointer to the node of the iterator.
        const Node* p;

        // Pointer back to the Graph container
        friend class Graph<V>;

        /** Private Constructor for node_iterator
         *
         * @param p the pointer to the Node
         */
        explicit node_iterator (const Node* p)
            : p{p} {
        }
    };

    /** @brief Method to access the pointer of the first node of this graph
     *
     * @return The pointer of the first node of this graph
     */
    node_iterator node_begin() const {
        return node_iterator(&node_list_[0]);
    };

    /** @brief Method to access the pointer after the last node of this graph
     *
     * @return The pointer after the last node of this graph
     */
    node_iterator node_end() const {
        return node_iterator(&node_list_.back()+1);
    };


    /** @class Graph::incident_iterator
     *
     * @brief Class representing the iterator objects for the graph nodes
     * that are incident to a given node.
     *
     * incident_iterator objects are used to loop on the nodes of the graph
     * that are incident to a given node.
     */
    class incident_iterator : private totally_ordered<incident_iterator>{

    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type        = Node;                     // Element type
        using pointer           = Node*;                    // Pointers to elements
        using reference         = Node&;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

        /** @brief Default constructor of incident_iterator. */
        incident_iterator()
            : p{nullptr} {
        }

        /** @brief Overriding the operator*()
         *
         * @return The incident Edge value of this incident_iterator
         */
        Edge operator*() const{
            return *p;
        };

        /** @brief Overriding the operator++()
         *
         * @return The next incident_iterator
         */
        incident_iterator &operator++(){
            p++;
            return *this;
        };

        /** @brief Overriding the operator==()
         *
         * @param[in] iit to compare this incident iterator with
         *
         * @return True if the two incident_iterator are equal
         */
        bool operator==(const incident_iterator &iit) const{
            return p == iit.p;
        };

    private:
        // The pointer to the incident edge of the iterator
        const Edge* p;

        // Allow Graph to access Node's private member data and functions.
        friend class Graph<V>;

        /** Private Constructor for incident_iterator
         *
         * @param p the pointer to the incident Edge
         */
        explicit incident_iterator (const Edge* p)
            : p{p} {
        }
    };


    /** @class Graph::edge_iterator
     *
     * @brief Class representing the edge objects for the graph.
     *
     * edge_iterator objects are used to loop on the edges of the graph.
     */
    class edge_iterator : private totally_ordered<edge_iterator>{

    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type        = Edge;                     // Element type
        using pointer           = Edge*;                    // Pointers to elements
        using reference         = Edge&;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

        /** @brief Default constructor of edge_iterator. */
        edge_iterator ()
            : p{nullptr} {
        }

        /** @brief Overriding the operator++()
         *
         * @return The next edge_iterator
         */
        Edge operator *() const{
            return *p ;
        }

        /** @brief Overriding the operator==()
         *
         * @param[in] eit to compare this edge iterator with
         *
         * @return True if the two edge_iterator are equal
         */
        bool operator ==( const edge_iterator &eit) const {
            return p == eit.p;
        }

        /** @brief Overriding the operator++()
         *
         * @return The next edge_iterator
         */
        edge_iterator & operator ++() {
            p++;
            return *this;
        }

    private :
        // The pointer to the edge of the iterator
        const Edge* p;

        // Allow Graph to access Node's private member data and functions.
        friend class Graph<V>;

        /** Private Constructor for edge_iterator
         *
         * @param p the pointer to the Edge
         */
        explicit edge_iterator (const Edge* p)
            : p{p} {
        }
    };


    /** @brief Method to access the pointer of the first edge of this graph
     *
     * @return The pointer of the first edge of this graph
     */
    edge_iterator edge_begin() const {
        return edge_iterator (&edge_list_[0]);
    };


    /** @brief Method to access the pointer after the last edge of this graph
     *
     * @return The pointer after the last edge of this graph
     */
    edge_iterator edge_end() const {
        return edge_iterator (&edge_list_.back()+1);
    };


    /** Return the number of nodes in the graph.
     *
     * Complexity: O(1).
     */
    size_type size() const {
        return size_;
    };


    /** Synonym for size(). */
    size_type num_nodes() const {
        return size();
    };


    /** Add a node to the graph, returning the added node.
     *
     * @param[in] position The new node's position
     *
     * @post new num_nodes() == old num_nodes() + 1
     * @post result_node.index() == old num_nodes()
     *
     * Complexity: O(1) amortized operations.
     */
    Node add_node(const Point &position) {
        // Take the position of the node
        Point p = position;

        // Initialize node value o default
        node_value_type val = node_value_type();

        // Create the new Node to add
        Node new_node = Node(this, next_node_index_);

        // Add node.
        index_to_nodes_[next_node_index_] = {p, next_node_index_, val, 0, {}};
        node_list_.push_back(new_node);

        // Increment number of nodes and latest index given
        ++next_node_index_;
        ++size_;

        return {this, next_node_index_ - 1};
    };


    /** Add a node to the graph, returning the added node.
     *
     * @param[in] position The new node's position
     * @param[in] value The new node's value
     *
     * @post new num_nodes() == old num_nodes() + 1
     * @post result_node.index() == old num_nodes()
     *
     * Complexity: O(1) amortized operations.
     */
    Node add_node(const Point &position, const node_value_type &val) {
        // Take the position of the node
        Point p = position;

        // Create the new Node to add
        Node new_node = Node(this, next_node_index_);

        // Add node.
        index_to_nodes_[next_node_index_] = {p, next_node_index_, val, 0, {}};
        node_list_.push_back(new_node);

        // Increment number of nodes and latest index given
        ++next_node_index_;
        ++size_;

        return {this, next_node_index_ - 1};
    };


    /** Determine if a Node belongs to this Graph
     *
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node &n) const {
        return ((n.graph_ == this)
                & (n.index() < this->size())
                & (0 <= n.index()));
    };


    /** Return the node with index @a i.
     *
     * @pre 0 <= @a i < num_nodes()
     *
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
        return Node(this, i);
    };


    /** Increase the degree of a given node.
     *
     * @pre @a a is a node of this graph
     *
     * @param[in] a The node to update
     *
     * @post The degree of the nodde is increased by one
     *
     * Complexity: O(1).
     */
    void incr_degree(Node a){
        // Extract current information of the node
        size_type a_index = a.index();
        internal_node node = index_to_nodes_[a_index];
        node_value_type val = node.value;
        Point point = node.point;

        // Increment the degree
        size_type degree = node.degree + 1;

        // Store the node with its new degree
        std::vector<Edge> edges = node.incident_edges;
        index_to_nodes_[a_index] = {point, a_index, val, degree, edges};
    };


    /** Update the list of incident edges of a node.
     *
     * @pre @a a is a node of this graph
     * @pre @a edge is an edge incident to @a a
     *
     * @param[in] a The node to update
     * @param[in] edge The incident edge to the node
     *
     * @post The edge is added to the list of incident edges of the node
     *
     * Complexity: O(1).
     */
    void add_incident_edge(Node a, Edge edge) {
        size_type a_index = a.index();
        internal_node node = index_to_nodes_[a_index];
        node.incident_edges.push_back(edge);
        index_to_nodes_[a_index] = node;
    };


    /** Update the value of a node.
     *
     * @pre @a a is a node of this graph
     *
     * @param[in] a The node to update
     * @param[in] value The new value of the node
     *
     * @post The value of the node is set to value
     *
     * Complexity: O(1).
     */
    void set_node_value(Node a, node_value_type val){
        size_type a_index = a.index();
        internal_node node = index_to_nodes_[a_index];
        node.value = val;
        index_to_nodes_[a_index] = node;
    };


    //
    // EDGES
    //

    /** @class Graph::Edge
     *
     * @brief Class representing the graph's edges.
     *
     * Edges are order-insensitive pairs of nodes. Two Edges with the same
     * nodes are considered equal if they connect the same nodes, in either
     * order.
     */
    class Edge : private totally_ordered<Edge>{

    public:
        /** Construct an invalid Edge. */
        Edge()
            : graph_(), index_(0) {
        }

        /** Return a node of this Edge */
        Node node1() const {
            return edge_fetch().node1;
        }

        /** Return the other node of this Edge */
        Node node2() const {
            return edge_fetch().node2;
        }

        /** Return the index of this Edge */
        size_type index() const {
            return edge_fetch().index;
        }

        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two
         * nodes.
         */
        bool operator==(const Edge& e) const {
            // Test if the nodes are equal in one of the two directions
            return ((edge_fetch().node1 == e.node1()
                     && edge_fetch().node2 == e.node2())
                    || (edge_fetch().node2 == e.node1()
                        && edge_fetch().node1 == e.node2()));
        }

        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge& e) const {
            return edge_fetch().index < e.index_;
        }

    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph<V>;

        // Pointer back to the Graph container
        Graph<V>* graph_;

        // The node unique identification number
        size_type index_;

        /** Private Constructor */
        Edge(const Graph* graph, size_type index)
            : graph_(const_cast<Graph*>(graph)), index_(index) {
        }

        /** Helper method to return the appropriate edge.
         * This finds the edge with the correct index using an
         * unordered map.
         */
        internal_edge& edge_fetch() const {
            return graph_->index_to_edges_[index_];
        }
    };


    /** Return the total number of edges in the graph.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    size_type num_edges() const {
        return num_edge_;
    }


    /** Return the edge with index @a i.
     *
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
        return Edge(this, i);
    }


    /** Test whether two nodes are connected by an edge.
     *
     * @pre @a a and @a b are valid nodes of this graph
     *
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node& a, const Node& b) const {
        // Result of if the two edges are equal
        bool edge_in_graph = false;

        // Extract the index of the two nodes a and b
        size_type a_index = a.index();
        size_type b_index = b.index();

        // Declare the unordered map used to map the node indexes connected
        // to a given node and the corresponding edges
        std::unordered_map<size_type, internal_edge> connected_nodes;

        // Look for the two nodes of the edges in the map of internal edges
        if(edges_to_index_.find(a_index) != edges_to_index_.end()){
            connected_nodes = edges_to_index_.at(a_index);
            if(connected_nodes.find(b_index) != connected_nodes.end()){
                edge_in_graph = true;
            }
        }
        return edge_in_graph;
    };


    /** Add an edge to the graph, or return the current edge if it already exists.
     *
     * @pre @a a and @a b are distinct valid nodes of this graph
     *
     * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
     *
     * @post has_edge(@a a, @a b) == true
     *
     * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
     *       Else,                        new num_edges() == old num_edges() + 1.
     *
     * Can invalidate edge indexes -- in other words, old edge(@a i) might not
     * equal new edge(@a i). Must not invalidate outstanding Edge objects.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge add_edge(const Node& a, const Node& b) {
        // Declare the index of the new edge
        size_type new_edge_index = next_edge_index_;;

        // Extract the index of the two nodes a and b
        size_type a_index = a.index();
        size_type b_index = b.index();

        // Declare the unordered map used to map the node indexes connected
        // to a given node and the corresponding edges
        std::unordered_map<size_type, internal_edge> connected_nodes;

        // Test if it is a new edge
        bool is_in_edges = has_edge(a, b);

        // If the edge is not new then it look for its index
        if(is_in_edges){

            if(edges_to_index_.find(a_index) != edges_to_index_.end()){
                connected_nodes = edges_to_index_.at(a_index);
                new_edge_index = connected_nodes.at(b_index).index;

            } else if(edges_to_index_.find(b_index) != edges_to_index_.end()){
                connected_nodes = edges_to_index_.at(b_index);
                new_edge_index = connected_nodes.at(a_index).index;
            }

        //If the edge is new then it is added to the map of edges
        } else {

            // Add the new incident edge to a
            add_incident_edge(a, {this, new_edge_index});
            incr_degree(a);

            // Add the new edge
            edge_list_.push_back({this, new_edge_index});
            index_to_edges_[new_edge_index] = {a, b, new_edge_index};

            if (edges_to_index_.find(a_index) != edges_to_index_.end()) {
                std::pair<size_type, internal_edge> new_connection
                        (b_index, {a, b, new_edge_index});
                edges_to_index_[a_index].insert(new_connection);
            } else {
                std::pair<size_type, internal_edge> new_connection
                        (b_index, {a, b, new_edge_index});
                connected_nodes.insert(new_connection);
                std::pair<size_type, std::unordered_map<size_type, internal_edge>>
                        new_edge(a_index, connected_nodes);
                edges_to_index_.insert(new_edge);
            }

            ++next_edge_index_;
            ++num_edge_;

        }
        return Edge(this, new_edge_index);
    };


    /** Remove all nodes and edges from this graph.
    * @post num_nodes() == 0 && num_edges() == 0
    *
    * Invalidates all outstanding Node and Edge objects.
    */
    void clear() {
        // Clear Node data
        size_ = 0;
        next_node_index_ = 0;
        index_to_nodes_ = {};
        node_list_ = {};

        // Clear Edge data
        num_edge_ = 0;
        next_edge_index_ = 0;
        index_to_edges_ = {};
        edges_to_index_ = {};
        edge_list_ = {};
    };

private:
    // Internal type for node
    struct internal_node {
        Point point;
        size_type index;
        node_value_type value;

        // Degree of the node
        size_type degree;

        std::vector<Edge> incident_edges;
    };

    // Internal type for edge
    struct internal_edge {
        Node node1;
        Node node2;
        size_type index;

        // Test whether this internal_edge and edge are equal
        bool operator==(const internal_edge &edge) const{
            return (node1 == edge.node1 && node2 == edge.node2);
        }

    };

    // Map of internal nodes with the node index as key
    std::unordered_map<size_type, internal_node> index_to_nodes_;

    // Map of internal edges with the edge index as key
    std::unordered_map<size_type, internal_edge> index_to_edges_;

    // Map of internal edges with the edge as key
    std::unordered_map<size_type, std::unordered_map<size_type, internal_edge>> edges_to_index_;

    //Node list
    std::vector<Node> node_list_;

    //Node list
    std::vector<Edge> edge_list_;

    // Number of nodes in the graph
    size_type size_;

    // Latest index given to a node in the graph
    size_type next_node_index_;

    // Number of distinct edges in the graph
    size_type num_edge_;

    // Latest index given to an edge in the graph
    size_type next_edge_index_;

};

#endif // CME212_GRAPH_HPP