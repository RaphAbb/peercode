#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP


/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 *
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
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
    Graph():
            nodes_info_(),
            edges_info_(),
            nodes_to_edges_(),
            graph_nodes_(),
            graph_edges_(),
            num_nodes_(0),
            next_node_index_(0),
            next_node_id_(0),
            num_edges_(0),
            next_edge_index_(0),
            next_edge_id_(0){
    }


    /** Default destructor */
    ~Graph() = default;


    //
    // NODES
    //

    /** @class Graph::Node
     *
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
            : graph_(), node_id_(){
        }

        /** @brief Return this node's position.
         *
         * The latest information about the node are stored
         * in a proxy design pattern. The node_fetch() method aims to look for
         * the corresponding node information in O(1).
         */
        Point& position(){
            return node_fetch().point_;
        }

        /** @brief Return a const to this node's position.
         *
         * The latest information about the node are stored
         * in a proxy design pattern. The node_fetch() method aims to look for
         * the corresponding node information in O(1).
         */
        const Point& position() const{
            return node_fetch().point_;
        }

        /** @biref Return this node's index.
         *
         * @post The index of the node which is a a number in the
         * range [0, graph_size).
         */
        size_type index() const {
            return graph_->nodes_info_[node_id_].index_;
        }

        /** @brief Return the value of this node.
         *
         * @return The value fo this node.
         */
        node_value_type &value() {
            return node_fetch().value_;
        }

        /** @brief Return the value of this node.
         *
         * @return A constant value of the value this node.
         */
        const node_value_type &value() const {
            const node_value_type& value = node_fetch().value_;
            return value;
        }

        /** @brief Return the degree of this node.
         *
         * @return The number of incident edges i.e. the number of edges
         * with this node as first node.
         *
         * @post degree >= 0.
         */
        size_type degree() const {
            size_type degree = node_fetch().degree_;
            return degree;
        }

        /** @brief Test whether this node and @a n are equal.
         *
         * @pram[in] n the node to compare this note with.
         *
         * @return True if and only if the nodes have the same graph
         * and the same index.
         *
         * Equal nodes have the same graph and the same index.
         */
        bool operator==(const Node &n) const {
            return (graph_ == n.graph_ && this->index() == n.index());
        }

        /** @brief Test whether this node is less than @a n in a global order.
         *
         * @parm[in] n the node to compare this note with.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any geometric meaning.
         *
         * The node ordering relation must obey trichotomy: For any two
         * nodes x and y, exactly one of x == y, x < y, and y < x is
         * true.
         */
        bool operator<(const Node &n) const {
            if (graph_ == n.graph_){
                return (this->index() < n.index());
            } else {
                auto graph_id = reinterpret_cast<intptr_t>(graph_);
                auto n_graph_id = reinterpret_cast<intptr_t>(n.graph_);
                return (graph_id < n_graph_id);
            }
        }

        /** @brief Return the first pointer to the list of incident edges of
         * this node.
         *
         * @pre The degree of this node is at least 1.
         *
         * @return The first pointer to the list of edges incident to
         * this node.
         */
        incident_iterator edge_begin() const {
            const Edge* p = &graph_->nodes_info_[node_id_].incident_edges_.front();
            return incident_iterator(this, p);
        }

        /** @brief Return the first pointer after the last element of the list
         * of incident edges to this node.
         *
         * @pre The degree of this node is at least 1
         *
         * @return The first pointer after the last element of the list
         * of edges incident to this node.
         */
        incident_iterator edge_end() const {
            const Edge* p = &(graph_->nodes_info_[node_id_].incident_edges_.back());
            ++p;
            return incident_iterator(this, p);
        }

    private:
        // Allow Graph to access class's private member data and functions
        friend class Graph;

        // Pointer back to the Graph container
        Graph<V, E>* graph_;

        // The node unique identification number
        unsigned int node_id_;

        /** @brief Construct a valid node.
         *
         * @param[in] graph the pointer to the graph representation.
         * @param[in] node_id the node id of the node to create.
         *
         * @pre node_id >= 0
         * @post Create a valid of the graph with the given index.
         */
        Node(const Graph<V, E>* graph, unsigned int node_id):
                graph_(const_cast<Graph<V, E> *>(graph)),
                node_id_(node_id) {
        }

        /** @brief Helper method to return the corresponding internal node.
         *
         * This finds the internal node with the correct id using an
         * unordered map.
         *
         * @return the corresponding internal node of this node.
         *
         * Complexity: O(1)
         */
        internal_node &node_fetch() const {
            return graph_->nodes_info_[node_id_];
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
        using value_type        = Node;
        using pointer           = Node*;
        using reference         = Node&;
        using difference_type   = std::ptrdiff_t;
        using iterator_category = std::input_iterator_tag;

        /** @brief Default constructor of node_iterator. */
        node_iterator()
            : graph_(nullptr), pointer_(nullptr) {
        }

        /** @brief Overriding the operator*()
         *
         * @return The Node value of the pointer
         */
        Node operator *() const{
            return {graph_, *pointer_};
        }

        /** @brief Overriding the operator==()
         *
         * @param[in] node_itr to compare this node iterator with
         *
         * @return True if the two node_iterator are equal
         */
        bool operator==(const node_iterator &node_itr) const {
            return pointer_ == node_itr.pointer_;
        }

        /** @brief Overriding the operator++()
         *
         * @return The next node_iterator
         */
        node_iterator & operator ++() {
            ++pointer_;
            return *this;
        }

    private :
        // Pointer back to the Graph container.
        Graph<V, E>* graph_;

        // The pointer to the node id of the iterator.
        const unsigned int* pointer_;

        // Allow Graph to access class's private member data and functions.
        friend class Graph<V, E>;

        /** Private Constructor for node_iterator
         *
         * @param[in] graph Pointer back to the Graph container
         * @param[in] p the pointer to the Node ID
         *
         */
        node_iterator(const Graph<V, E>* graph, const unsigned int* p)
            : graph_(const_cast<Graph<V, E> *>(graph)), pointer_(p) {
        }
    };


    /** @brief Method to access the pointer of the first node of this graph.
     *
     * @return The pointer of the first node of this graph.
     */
    node_iterator node_begin() const {
        const unsigned int* p = &graph_nodes_.front();
        return node_iterator(this, p);
    };


    /** @brief Method to access the pointer after the last node of this graph.
     *
     * @return The pointer after the last node of this graph.
     */
    node_iterator node_end() const {
        const unsigned int* p = &graph_nodes_.back();
        ++p;
        return node_iterator(this, p);
    };


    /** @class Graph::incident_iterator
     *
     * @brief Class representing the iterator objects for the graph nodes
     * that are incident to a given node.
     *
     * incident_iterator objects are used to loop on the edges of the graph
     * that are incident to a given node.
     */
    class incident_iterator : private totally_ordered<incident_iterator>{
    public:
        using value_type        = Edge;
        using pointer           = Edge*;
        using reference         = Edge&;
        using difference_type   = std::ptrdiff_t;
        using iterator_category = std::input_iterator_tag;

        /** @brief Default constructor of incident_iterator.
         *
         * Construct an invalid incident_iterator object.
         */
        incident_iterator()
            : node_(nullptr), pointer_(nullptr) {
        }

        /** @brief Overriding the operator*()
         *
         * @return The incident Edge value of this incident_iterator.
         */
        Edge operator*() const{
            return *pointer_;
        };

        /** @brief Overriding the operator++()
         *
         * @return The next incident_iterator.
         */
        incident_iterator &operator++(){
            ++pointer_;
            return *this;
        };

        /** @brief Overriding the operator==()
         *
         * @param[in] iit to compare this incident iterator with.
         *
         * @return True if the two incident_iterator are equal.
         */
        bool operator==(const incident_iterator &iit) const{
            return pointer_ == iit.pointer_;
        };

    private:
        // Allow Graph to access class's private member data and functions.
        friend class Graph<V, E>;

        // Allow Node to access class's private member data and functions.
        friend class Node;

        // Pointer back to the Graph container.
        Graph<V, E>::Node* node_;

        // The pointer to the incident edge of the iterator
        const Edge* pointer_;

        /** Private Constructor for incident_iterator
         *
         * @param[in] node the pointer to the node representation.
         * @param[in] p the pointer to the incident Edge
         *
         */
        incident_iterator(const Graph<V, E>::Node* node, const Edge* p)
            : node_(const_cast<Graph<V, E>::Node*>(node)), pointer_(p) {
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
        using value_type        = Edge;
        using pointer           = Edge*;
        using reference         = Edge&;
        using difference_type   = std::ptrdiff_t;
        using iterator_category = std::input_iterator_tag;

        /** @brief Default constructor of edge_iterator. */
        edge_iterator()
            : graph_(nullptr), pointer_(nullptr) {
        }

        /** @brief Overriding the operator++()
         *
         * @return The next edge_iterator.
         */
        Edge operator*() const{
            unsigned int edge_id = *pointer_;
            return {graph_, edge_id, graph_->edges_info_[edge_id].node1_, graph_->edges_info_[edge_id].node2_};
        }

        /** @brief Overriding the operator==()
         *
         * @param[in] eit to compare this edge iterator with.
         *
         * @return True if the two edge_iterator are equal
         */
        bool operator==(const edge_iterator &eit) const {
            return pointer_ == eit.pointer_;
        }

        /** @brief Overriding the operator++()
         *
         * @return The next edge_iterator.
         */
        edge_iterator & operator ++() {
            ++pointer_;
            return *this;
        }

    private :
        // Pointer back to the Graph container.
        Graph<V, E>* graph_;

        // The pointer to the edge of the iterator
        const unsigned int* pointer_;

        // Allow Graph to access class's private member data and functions.
        friend class Graph<V, E>;

        /** Private Constructor for edge_iterator
         *
         * @param[in] graph the pointer to the graph representation.
         * @param[in] p the pointer to the Edge ID
         */
        edge_iterator(const Graph<V, E>* graph, const unsigned int* p)
            : graph_(const_cast<Graph<V, E> *>(graph)), pointer_(p) {
        }
    };


    /** @brief Method to access the pointer of the first edge of this graph
     *
     * @return The pointer of the first edge of this graph
     */
    edge_iterator edge_begin() const {
        const unsigned int* p = &graph_edges_.front();
        return edge_iterator(this, p);
    };


    /** @brief Method to access the pointer after the last edge of this graph
     *
     * @return The pointer after the last edge of this graph.
     */
    edge_iterator edge_end() const {
        const unsigned int* p = &(graph_edges_.back());
        ++p;
        return edge_iterator(this, p);
    };


    /** @brief Return the number of nodes in the graph.
     *
     * Complexity: O(1)
     */
    size_type size() const {
        return num_nodes_;
    };


    /** Synonym for size().
     *
     * Complexity: O(1)
     */
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

        // Initialize node value to default
        node_value_type value = node_value_type();

        // Create the Node to add with its node id
        Node node = Node(this, next_node_id_);

        // Add new node to the list of node
        graph_nodes_.push_back(next_node_id_);

        // Add corresponding internal node
        internal_node intern_node = internal_node();
        intern_node.point_ = position;
        intern_node.index_ = next_node_index_;
        intern_node.value_ = value;
        nodes_info_.insert(std::pair<unsigned int, internal_node>(next_node_id_, intern_node));

        // Increment latest node index and node id
        ++next_node_index_;
        ++next_node_id_;

        // Update total number of nodes
        ++num_nodes_;

        return node;
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
    Node add_node(const Point &position, const node_value_type &value) {

        // Create the Node to add with its node id
        Node node = Node(this, next_node_id_);

        // Add new node to the list of node
        graph_nodes_.push_back(next_node_id_);

        // Add corresponding internal node
        internal_node intern_node = internal_node();
        intern_node.point_ = position;
        intern_node.index_ = next_node_index_;
        intern_node.value_ = value;
        nodes_info_.insert(std::pair<unsigned int, internal_node>(next_node_id_, intern_node));

        // Increment latest node index and node id
        ++next_node_index_;
        ++next_node_id_;

        // Update total number of nodes
        ++num_nodes_;

        return node;
    };


    /** Determine if a Node belongs to this Graph
     *
     * @param[in] n the node to cjeck if it is in this graph.
     *
     * @return True if @a n is currently a Node of this Graph.
     *
     * Complexity: O(1)
     */
    bool has_node(const Node &n) const {
        size_type node_index = n.index();
        return ((n.graph_ == this)
                & (node_index < this->size())
                & (0 <= node_index));
    };


    /** Return the node with index @a i.
     *
     * @pre 0 <= @a i < num_nodes()
     *
     * @post result_node.index() == i
     *
     * Complexity: O(1)
     */
    Node node(size_type i) const {
        assert((0 <= i) && (i < this->size()));
        unsigned int node_id = graph_nodes_[i];
        return {this, node_id};
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
     * Complexity: O(1)
     */
    void update_node_value(Node &n, node_value_type value){
        assert(has_node(n));
        unsigned int node_id = n.node_id_;
        nodes_info_[node_id].value_ = value;
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
    using edge_value_type = E;

    class Edge : private totally_ordered<Edge>{
    public:
        /** Construct an invalid Edge. */
        Edge()
            : graph_(nullptr), edge_id_(), node1_id_(), node2_id_() {
        }

        /** Return a node of this Edge */
        Node node1() const {
            return {graph_, node1_id_};
        }

        /** Return the other node of this Edge */
        Node node2() const {
            return {graph_, node2_id_};
        }

        /** Return the index of this Edge */
        size_type index() const {
            return edge_fetch().index_;
        }

        /** Return the distance between the two nodes of this Edge */
        double length() const{
            Point p = this->node2().position() - this->node1().position();
            return norm(p);
        }

        /** @brief Return the value of this node.
         *
         * @return The value fo this node.
         */
        edge_value_type &value() {
            return edge_fetch().value_;
        }

        /** @brief Return the value of this node.
         *
         * @return A constant value of the value this node.
         */
        const edge_value_type &value() const {
            const edge_value_type& value = edge_fetch().value_;
            return value;
        }

        /** Test whether this edge and @a e are equal.
         *
         * @param[in] e the edge to compare this edge with.
         *
         * Equal edges represent the same undirected edge between two
         * nodes.
         *
         * Complexity: O(1)
         */
        bool operator==(const Edge& e) const {
            // Test if the nodes are equal in one of the two directions
            return ((this->node1() == e.node1()
                     && this->node2() == e.node2())
                    || (this->node2() == e.node1()
                        && this->node1() == e.node2()));
        }

        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge& e) const {
            if(*this == e){
                return false;
            } else{
                return (this->node1() < e.node1()
                    || ((this->node1() == e.node1()) && this->node2() < e.node2()));
            }

        }

    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph<V, E>;

        // Pointer back to the Graph container
        Graph<V, E>* graph_;

        // The edge unique identification number
        unsigned int edge_id_;

        // The first node ID delimiting the edge
        unsigned int node1_id_;

        // The second node ID delimiting the edge
        unsigned int node2_id_;

        /** Private constructor for Edge object
         *
         * @param graph the pointer to the graph representation.
         * @param edge_id the edge unique identification number of this edge.
         * @param node1 the first node ID of the edge.
         * @param node2 the second node ID of the edge.
         */
        Edge(const Graph<V, E>* graph, unsigned int edge_id, const unsigned int node1_id, const unsigned int node2_id):
                graph_(const_cast<Graph<V, E> *>(graph)),
                edge_id_(edge_id),
                node1_id_(node1_id),
                node2_id_(node2_id){
        }

        /** Private constructor for Edge object
         *
         * @param graph the pointer to the graph representation.
         * @param edge_id the edge unique identification number of this edge.
         * @param node1 the first node of the edge.
         * @param node2 the second node of the edge.
         */
        Edge(const Graph<V, E>* graph, unsigned int edge_id, const Node& node1, const Node& node2):
                graph_(const_cast<Graph<V, E> *>(graph)),
                edge_id_(edge_id),
                node1_id_(node1.node_id_),
                node2_id_(node2.node_id_){
        }

        /** @brief Helper method to return the corresponding internal edge.
          *
          * This finds the internal edge with the correct id using an
          * unordered map.
          *
          * @return the corresponding internal edge of this edge.
          *
          * Complexity: O(1)
          */
        internal_edge& edge_fetch() const {
            return graph_->edges_info_[edge_id_];
        }
    };


    /** Return the total number of edges in the graph.
     *
     * Complexity: O(1)
     */
    size_type num_edges() const {
        return num_edges_;
    };


    /** Return the edge with index @a i.
     *
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: O(1)
     */
    Edge edge(size_type i) const {
        unsigned int edge_id = graph_edges_[i];
        internal_edge intern_edge = const_cast<Graph<V, E>*>(this)->edges_info_[edge_id];
        return {this, edge_id, intern_edge.node1_, intern_edge.node2_};
    };


    /** Test whether two nodes are connected by an edge.
     *
     * @pre @a a and @a b are valid nodes of this graph
     *
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: O(1) amortized operations.
     */
    bool has_edge(const Node& a, const Node& b) const {
        unsigned int node_a_id = a.node_id_;
        unsigned int node_b_id = b.node_id_;
        std::pair<unsigned int, unsigned int> node_pair_ab(node_a_id, node_b_id);
        std::pair<unsigned int, unsigned int> node_pair_ba(node_b_id, node_a_id);

        if(nodes_to_edges_.find(node_pair_ab) != nodes_to_edges_.end()){
            return true;
        } else return nodes_to_edges_.find(node_pair_ba) != nodes_to_edges_.end();
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
     * Complexity: O(1) amortized operations.
     */
    Edge add_edge(const Node& a, const Node& b) {
        unsigned int node_a_id = a.node_id_;
        unsigned int node_b_id = b.node_id_;
        std::pair<unsigned int, unsigned int> node_pair_ab(node_a_id, node_b_id);
        std::pair<unsigned int, unsigned int> node_pair_ba(node_b_id, node_a_id);

        // If the edge is already in the graph then it returns it
        if(nodes_to_edges_.find(node_pair_ab) != nodes_to_edges_.end()){
            unsigned int edge_id = nodes_to_edges_[node_pair_ab];
            return {this, edge_id, a, b};

        } else if (nodes_to_edges_.find(node_pair_ba) != nodes_to_edges_.end()){
            unsigned int edge_id = nodes_to_edges_[node_pair_ba];
            return {this, edge_id, a, b};
        }
        // If the edge is new it is added to the graph
        else {
            // Create new edge
            Edge edge_ab = {this, next_edge_id_, a, b};
            Edge edge_ba = {this, next_edge_id_, b, a};

            // Add corresponding internal edge
            edge_value_type value = edge_value_type();
            edges_info_.insert(std::pair<unsigned int, internal_edge>(next_edge_id_, {a, b, next_edge_index_, value}));

            // Add new edge
            graph_edges_.push_back(next_edge_id_);

            // Update incident edges of node a
            size_type degree_a = nodes_info_[node_a_id].degree_;
            nodes_info_[node_a_id].incident_edges_.push_back(edge_ab);
            nodes_info_[node_a_id].incident_indexes_.insert(std::pair<unsigned int, size_type>(next_edge_id_, degree_a));
            nodes_info_[node_a_id].degree_ = degree_a + 1;

            // Update incident edges of node b
            size_type degree_b = nodes_info_[node_b_id].degree_;
            nodes_info_[node_b_id].incident_edges_.push_back(edge_ba);
            nodes_info_[node_b_id].incident_indexes_.insert(std::pair<unsigned int, size_type>(next_edge_id_, degree_b));
            nodes_info_[node_b_id].degree_ = degree_b + 1;

            // Add new edge to mapping node pair to edge
            nodes_to_edges_.insert(std::pair<std::pair<unsigned int, unsigned int>, unsigned int>(node_pair_ab, next_edge_id_));
            nodes_to_edges_.insert(std::pair<std::pair<unsigned int, unsigned int>, unsigned int>(node_pair_ba, next_edge_id_));

            // Increment latest edge index and edge id
            ++next_edge_index_;
            ++next_edge_id_;

            // Update total number of edges
            ++num_edges_;

            return edge_ab;
        }

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
     * Complexity: O(1) amortized operations.
     */
    Edge add_edge(const Node& a, const Node& b, edge_value_type value) {
        unsigned int node_a_id = a.node_id_;
        unsigned int node_b_id = b.node_id_;
        std::pair<unsigned int, unsigned int> node_pair_ab(node_a_id, node_b_id);
        std::pair<unsigned int, unsigned int> node_pair_ba(node_b_id, node_a_id);

        // If the edge is already in the graph then it returns it
        if(nodes_to_edges_.find(node_pair_ab) != nodes_to_edges_.end()){
            unsigned int edge_id = nodes_to_edges_[node_pair_ab];
            return {this, edge_id, a, b};

        } else if (nodes_to_edges_.find(node_pair_ba) != nodes_to_edges_.end()){
            unsigned int edge_id = nodes_to_edges_[node_pair_ba];
            return {this, edge_id, a, b};
        }
            // If the edge is new it is added to the graph
        else {
            // Create new edge
            Edge edge_ab = {this, next_edge_id_, a, b};
            Edge edge_ba = {this, next_edge_id_, b, a};

            // Add corresponding internal edge
            edges_info_.insert(std::pair<unsigned int, internal_edge>(next_edge_id_, {a, b, next_edge_index_, value}));

            // Add new edge
            graph_edges_.push_back(next_edge_id_);

            // Update incident edges of node a
            size_type degree_a = nodes_info_[node_a_id].degree_;
            nodes_info_[node_a_id].incident_edges_.push_back(edge_ab);
            nodes_info_[node_a_id].incident_indexes_.insert(std::pair<unsigned int, size_type>(next_edge_id_, degree_a));
            nodes_info_[node_a_id].degree_ = degree_a + 1;

            // Update incident edges of node b
            size_type degree_b = nodes_info_[node_b_id].degree_;
            nodes_info_[node_b_id].incident_edges_.push_back(edge_ba);
            nodes_info_[node_b_id].incident_indexes_.insert(std::pair<unsigned int, size_type>(next_edge_id_, degree_b));
            nodes_info_[node_b_id].degree_ = degree_b + 1;

            // Add new edge to mapping node pair to edge
            nodes_to_edges_.insert(std::pair<std::pair<unsigned int, unsigned int>, unsigned int>(node_pair_ab, next_edge_id_));
            nodes_to_edges_.insert(std::pair<std::pair<unsigned int, unsigned int>, unsigned int>(node_pair_ba, next_edge_id_));

            // Increment latest edge index and edge id
            ++next_edge_index_;
            ++next_edge_id_;

            // Update total number of edges
            ++num_edges_;

            return edge_ab;
        }

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
     * Complexity: O(1)
     */
    void update_edge_value(Edge e, edge_value_type value){
        assert(has_edge(e.node1(), e.node2()));
        unsigned int edge_id = e.edge_id_;
        edges_info_[edge_id].value_ = value;
    };


    /** Remove a node.
     *
     * @pre @a n is a node of this graph
     *
     * @param[in] n The node to update
     *
     * @post The node @n and all the connected edges are removed
     *
     * Complexity: O(degree of @n).
     */
    size_type remove_node(const Node &n){
        if (has_node(n)){

            // Extract the infos of the node to be removed
            size_type node_index = n.index();
            unsigned int node_id = n.node_id_;

            // Remove the incident edges of the node to be removed
            if(n.degree() > 0){
                for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {

                    // Extract the infos of the edge to be removed
                    Edge e = *ei;
                    unsigned int edge_id = e.edge_id_;
                    size_type edge_index = e.index();

                    unsigned int node1_id = e.node1().node_id_;
                    unsigned int node2_id = e.node2().node_id_;

                    // Remove the edge from the incident edges of node 2
                    std::pair<unsigned int, unsigned int> node_pair_12(node1_id, node2_id);
                    std::pair<unsigned int, unsigned int> node_pair_21(node2_id, node1_id);
                    nodes_to_edges_.erase(node_pair_12);
                    nodes_to_edges_.erase(node_pair_21);

                    // Remove edge from incident edge
                    Edge back_incident_edge = nodes_info_[node2_id].incident_edges_.back();
                    unsigned int back_incident_edge_id = back_incident_edge.edge_id_;
                    size_type incident_index = nodes_info_[node2_id].incident_indexes_[edge_id];

                    if (e == back_incident_edge){
                        // Remove the incident edge from the vector of edge id
                        nodes_info_[node2_id].incident_edges_.pop_back();
                        nodes_info_[node2_id].incident_indexes_.erase(edge_id);

                    } else {
                        // Remove the incident edge from the vector of edge id by switching
                        // it with the last edges of the vector and then removing
                        // the last edge of the vector
                        nodes_info_[node2_id].incident_edges_[incident_index] = back_incident_edge;
                        nodes_info_[node2_id].incident_edges_.pop_back();

                        // Update incident_indexes_
                        nodes_info_[node2_id].incident_indexes_.erase(edge_id);
                        nodes_info_[node2_id].incident_indexes_[back_incident_edge_id] = incident_index;
                    }

                    // Update degree
                    nodes_info_[node2_id].degree_ -= 1;

                    // Remove the edge from the list of edges
                    unsigned int back_edge_id = graph_edges_.back();

                    if(edge_id == back_edge_id){
                        // Remove the node from the list of nodes
                        graph_edges_.pop_back();

                    } else {
                        // Remove the node from the list of nodes by switching it with
                        // the last node of the list and then removing the last node of
                        // the list
                        graph_edges_[edge_index] = back_edge_id;
                        graph_edges_.pop_back();

                        // Update the index of the remaining nodes
                        edges_info_[back_edge_id].index_ = edge_index;
                    }

                    // Remove the internal node corresponding to the node to be
                    // removed
                    edges_info_.erase(edge_id);

                    // Update the number of edges of the graph
                    --num_edges_;
                }
            }

            // Extract the node id which is at the end
            unsigned int back_node_id = graph_nodes_.back();

            if (node_id == back_node_id){
                // Remove the node from the list of nodes
                graph_nodes_.pop_back();

            } else {
                // Remove the node from the list of nodes by switching it with
                // the last node of the list and then removing the last node of
                // the list
                graph_nodes_[node_index] = back_node_id;
                graph_nodes_.pop_back();

                // Update the index of the remaining nodes
                nodes_info_[back_node_id].index_ = node_index;

            }

            // Remove the internal node corresponding to the node to be
            // removed
            nodes_info_.erase(node_id);

            // Update the size of the graph
            --num_nodes_;

            // Update the latest node index used in the graph
            --next_node_index_;

            return 1;

        } else {
            return 0;
        }

    };


    /** Remove a node based on its node iterator.
     *
     * @pre @a n_it is a node_iterator of this graph
     *
     * @param[in] n_it is a node_iterator to remove
     *
     * @post The node associated to @n_it and all
     * the connected edges are removed
     *
     * @return the next valid node iterator
     *
     * Complexity: O(degree of n associated to n_it).
     */
    node_iterator remove_node(node_iterator n_it){
        Node n = *n_it;
        size_type node_index = n.index();
        remove_node(n);
        if (node_index == num_nodes_){
            return node_iterator(&graph_nodes_.back());

        } else{
            return node_iterator(&graph_nodes_[node_index]);
        }

    };


    /** Remove an edge based on the two nodes that define it.
     *
     * @pre @a a is a node of this graph
     * @pre @a b is a node of this graph
     *
     * @param[in] a the first node of the edge
     * @param[in] b the second node of the edge
     *
     * @post The edge @a - @b is removed
     *
     * @return 1 or 0 if the edge is indeed removed
     *
     * Complexity: O(1) amortized operations.
     */
    size_type remove_edge(const Node &a, const Node &b){

        unsigned int node_a_id = a.node_id_;
        unsigned int node_b_id = b.node_id_;
        std::pair<unsigned int, unsigned int> node_pair_ab(node_a_id, node_b_id);
        std::pair<unsigned int, unsigned int> node_pair_ba(node_b_id, node_a_id);

        if(nodes_to_edges_.find(node_pair_ab) != nodes_to_edges_.end()){
            Edge e = {this, nodes_to_edges_[node_pair_ab], a, b};
            size_type rst = this->remove_edge(e);
            return rst;

        } else if (nodes_to_edges_.find(node_pair_ba) != nodes_to_edges_.end()){
            Edge e = {this, nodes_to_edges_[node_pair_ba], b, a};
            size_type rst = this->remove_edge(e);
            return rst;

        } else{
            return 0;
        }

    };


    /** Remove an edge.
     *
     * @pre @a e is an edge of this graph
     *
     * @param[in] e the edge to remove
     *
     * @post The edge @e is removed
     *
     * @return 1 or 0 if the edge is indeed removed
     *
     * Complexity: O(1) amortized operations.
     */
    size_type remove_edge(const Edge &e){

        if (has_edge(e.node1(), e.node2())){

            // Extract the infos of the edge to be removed
            unsigned int edge_id = e.edge_id_;
            size_type  edge_index = e.index();
            unsigned int node1_id = e.node1().node_id_;
            unsigned int node2_id = e.node2().node_id_;

            // Remove the edge from nodes_to_edges_
            std::pair<unsigned int, unsigned int> node_pair_12(node1_id, node2_id);
            std::pair<unsigned int, unsigned int> node_pair_21(node2_id, node1_id);
            nodes_to_edges_.erase(node_pair_12);
            nodes_to_edges_.erase(node_pair_21);

            // Remove edge from incident edges of node 1
            Edge back_incident_edge1 = nodes_info_[node1_id].incident_edges_.back();
            unsigned int back_incident_edge1_id = back_incident_edge1.edge_id_;
            size_type incident_index1 = nodes_info_[node1_id].incident_indexes_[edge_id];

            if (e == back_incident_edge1){
                // Remove the incident edge from the vector of edge id
                nodes_info_[node1_id].incident_edges_.pop_back();
                nodes_info_[node1_id].incident_indexes_.erase(edge_id);

            } else {
                // Remove the incident edge from the vector of edge id by switching
                // it with the last edges of the vector and then removing
                // the last edge of the vector
                nodes_info_[node1_id].incident_edges_[incident_index1] = back_incident_edge1;
                nodes_info_[node1_id].incident_edges_.pop_back();

                // Update incident_indexes_
                nodes_info_[node1_id].incident_indexes_.erase(edge_id);
                nodes_info_[node1_id].incident_indexes_[back_incident_edge1_id] = incident_index1;
            }

            // Remove edge from incident edges of node 2
            Edge back_incident_edge2 = nodes_info_[node2_id].incident_edges_.back();
            unsigned int back_incident_edge2_id = back_incident_edge2.edge_id_;
            size_type incident_index2 = nodes_info_[node2_id].incident_indexes_[edge_id];

            if (e == back_incident_edge2){
                // Remove the incident edge from the vector of edge id
                nodes_info_[node2_id].incident_edges_.pop_back();
                nodes_info_[node2_id].incident_indexes_.erase(edge_id);

            } else {
                // Remove the incident edge from the vector of edge id by switching
                // it with the last edges of the vector and then removing
                // the last edge of the vector
                nodes_info_[node2_id].incident_edges_[incident_index2] = back_incident_edge2;
                nodes_info_[node2_id].incident_edges_.pop_back();

                // Update incident_indexes_
                nodes_info_[node2_id].incident_indexes_.erase(edge_id);
                nodes_info_[node2_id].incident_indexes_[back_incident_edge2_id] = incident_index2;
            }

            // Update degree
            nodes_info_[node1_id].degree_ -= 1;
            nodes_info_[node2_id].degree_ -= 1;

            // Remove the node from the list of nodes by switching it with
            // the last node of the list and then removing the last node of
            // the list
            unsigned int back_edge_id = graph_edges_.back();

            if (edge_id == back_edge_id){
                graph_edges_.pop_back();
            } else {
                graph_edges_[edge_index] = back_edge_id;
                graph_edges_.pop_back();

                // Update the index of the remaining nodes
                edges_info_[back_edge_id].index_ = edge_index;
            }

            // Remove the internal node corresponding to the node to be
            // removed
            edges_info_.erase(edge_id);

            // Update the number of edges of the graph
            --num_edges_;

            // Update the latest edge index used in the graph
            --next_edge_index_;

            return 1;

        }else {
            return 0;
        }

    };


    /** Remove an edge based on its node iterator.
     *
     * @pre @a e_it is a edge_iterator of this graph
     *
     * @param[in] e_it is a edge_iterator to remove
     *
     * @post The edge associated to @e_it is removed
     *
     * @return the next valid edge iterator
     *
     * Complexity: O(1) amortized operations.
     */
    edge_iterator remove_edge(edge_iterator e_it){
        Edge e = *e_it;
        size_type edge_index = e.index();
        this->remove_edge(e);

        if (edge_index == num_edges_){
            return edge_iterator(&graph_edges_.back());
        } else{
            return edge_iterator(&graph_edges_.at(edge_index));
        }
    };


    /** Remove all nodes and edges from this graph.
     *
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        nodes_info_ = {};
        edges_info_ = {};
        nodes_to_edges_ = {};
        graph_nodes_ = {};
        graph_edges_ = {};
        num_nodes_ = 0;
        next_node_index_ = 0;
        next_node_id_ = 0;
        num_edges_ = 0;
        next_edge_index_ = 0;
        next_edge_id_ = 0;
    };

private:
    // Map of internal nodes with the node index as key
    std::unordered_map<unsigned int, internal_node> nodes_info_;

    // Map of internal edges with the edge index as key
    std::unordered_map<unsigned int, internal_edge> edges_info_;

    // Map of internal edges with the edge as key
    std::map<std::pair<unsigned int, unsigned int>, unsigned int> nodes_to_edges_;

    // Vector of nodes in the graph
    std::vector<unsigned int> graph_nodes_;

    // Vector of edges in the graph
    std::vector<unsigned int> graph_edges_;

    // Number of nodes in the graph
    size_type num_nodes_;

    // Index to give to the next node to be added in to the graph
    size_type next_node_index_;

    // Id to give to the next node to be added in to the graph
    unsigned int next_node_id_;

    // Number of distinct edges in the graph
    size_type num_edges_;

    // Index to give to the next edge to be added in to the graph
    size_type next_edge_index_;

    // Id to give to the next edge to be added in to the graph
    unsigned int next_edge_id_;

    // Internal node containing the information about its corresponding node
    struct internal_node {
        // Position of the corresponding node
        Point point_;

        // Index of the corresponding node
        size_type index_;

        // Value of the corresponding node
        node_value_type value_;

        // Degree of the corresponding node
        size_type degree_;

        // Incident edges of the corresponding node
        std::vector<Edge> incident_edges_;

        // Mapping of the edge id with the position of an edge in the vector
        // of incident edges
        std::unordered_map<unsigned int, size_type> incident_indexes_;

        /** @brief Default constructor of internal node. */
        internal_node(){
            point_ = Point();
            index_ = size_type();
            value_ = node_value_type();
            degree_ = 0;
            incident_edges_ = std::vector<Edge>();
            incident_indexes_ = std::unordered_map<unsigned int, size_type>();
        }

        /** @brief Constructor of internal node. */
        internal_node(Point point,
                      size_type index,
                      node_value_type value,
                      size_type degree,
                      std::vector<Edge> incident_edges,
                      std::unordered_map<unsigned int, size_type> incident_indexes){
            point_ = point;
            index_ = index;
            value_ = value;
            degree_ = degree;
            incident_edges_ = incident_edges;
            incident_indexes_ = incident_indexes;
        }

    };

    // Internal edge containing the information about its corresponding edge
    struct internal_edge {
        // First node of the corresponding edge
        Node node1_;

        // Second node of the corresponding edge
        Node node2_;

        // Index of the corresponding edge
        size_type index_;

        // Value of the corresponding edge
        edge_value_type value_;

        /** @brief Default constructor of internal edge. */
        internal_edge():
                node1_(),
                node2_(),
                index_(0),
                value_(){
        }

        /** @brief Constructor of internal node. */
        internal_edge(Node node1,
                      Node node2,
                      size_type index,
                      edge_value_type value)
                : node1_(node1),
                  node2_(node2),
                  index_(index),
                  value_(value){
        }

    };

};

#endif // CME212_GRAPH_HPP