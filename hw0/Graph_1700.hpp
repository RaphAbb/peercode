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
class Graph {
    private:

        // Predeclare the internal struct for the nodes
        struct internal_node;

        // Predeclare the internal struct for the edges
        struct internal_edge;

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

        /** Predeclaration of Edge type. */
        class Edge;
        /** Synonym for Edge (following STL conventions). */
        using edge_type = Edge;

        /** Type of indexes and sizes.
          Return type of Graph::Node::index(), Graph::num_nodes(),
          Graph::num_edges(), and argument type of Graph::node(size_type) */
        using size_type = unsigned;

        //
        // CONSTRUCTORS AND DESTRUCTOR
        //

        /** Construct an empty graph. */
        Graph() 
            : index_to_nodes_(), 
            index_to_edges_(),
            edges_to_index_(),
            size_(0), 
            next_node_index_(0), 
            num_edge_(0), 
            next_edge_index_(0) {
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

                /** Return this node's position. */
                const Point& position() const {
                    return node_fetch().point;
                }

                /** Return this node's index, a number in the 
                * range [0, graph_size). 
                */
                size_type index() const {
                    return node_fetch().index;
                }

                /** Test whether this node and @a n are equal.
                 *
                 * Equal nodes have the same graph and the same index.
                 */
                bool operator==(const Node& n) const {
                    bool equal = (graph_ == n.graph_ && index_ == n.index());
                    return equal;
                }

                /** Test whether this node is less than @a n in a global order.
                 *
                 * This ordering function is useful for STL containers such as
                 * std::map<>. It need not have any geometric meaning.
                 *
                 * The node ordering relation must obey trichotomy: For any two 
                 * nodes x and y, exactly one of x == y, x < y, and y < x is 
                 * true.
                 */
                bool operator<(const Node& n) const {
                    if (node_fetch().index < n.index_){
                        return true;
                    } else {
                        return false;
                    }
                }

            private:
                // Allow Graph to access Node's private member data and 
                // functions.
                friend class Graph;

                // Pointer back to the Graph container
                Graph* graph_;

                // The node unique identification number
                size_type index_;

                /** Private Constructor */
                Node(const Graph* graph, size_type index)
                    : graph_(const_cast<Graph*>(graph)), index_(index) {
                }

                /** Helper method to return the appropriate node.
                 * This finds the node with the correct index using an 
                 * unordered map.
                 */
                internal_node& node_fetch() const {
                    return graph_->index_to_nodes_[index_];
                    assert(false);
                }
        };


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

        /** Add a node to the graph, returning the added node.
        * @param[in] position The new node's position
        * @post new num_nodes() == old num_nodes() + 1
        * @post result_node.index() == old num_nodes()
        *
        * Complexity: O(1) amortized operations.
        */
        Node add_node(const Point& position) {
            Point p = position;

            // Add node
            index_to_nodes_[next_node_index_] = {p, next_node_index_};

            // Increment number of nodes and latest index given
            ++next_node_index_;
            ++size_;
            return Node(this, next_node_index_-1);
        }

        /** Determine if a Node belongs to this Graph
        * @return True if @a n is currently a Node of this Graph
        *
        * Complexity: O(1).
        */
        bool has_node(const Node& n) const {
            return (n.graph_ == this);
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
        * Edges are order-insensitive pairs of nodes. Two Edges with the same 
        * nodes are considered equal if they connect the same nodes, in either 
        * order.
        */
        class Edge {
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
                    if (edge_fetch().index < e.index_){
                        return true;
                    } else {
                        return false;
                    }
                }

            private:
                // Allow Graph to access Edge's private member data and 
                // functions.
                friend class Graph;

                // Pointer back to the Graph container
                Graph* graph_;

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
                    assert(false);
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
            // Result of if the two edges are equal
            bool edge_in_graph = false;

            // Extract the index of the two nodes a and b
            size_type a_index = a.index();
            size_type b_index = b.index();

            // Declare the unordered map used to map the node indeces connected
            // to a given node and the corresponnding edges
            std::unordered_map<size_type, internal_edge> connected_nodes;

            // Look for the two nodes of the edges in the map of internal edges
            if(edges_to_index_.find(a_index) != edges_to_index_.end()){
                connected_nodes = edges_to_index_.at(a_index);
                if(connected_nodes.find(b_index) != connected_nodes.end()){
                    edge_in_graph = true;
                }
            } else if(edges_to_index_.find(b_index) != edges_to_index_.end()){
                connected_nodes = edges_to_index_.at(b_index);
                if(connected_nodes.find(a_index) != connected_nodes.end()){
                    edge_in_graph = true;
                }
            }

            return edge_in_graph;
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
            // Declare the index of the new edge
            size_type new_edge_index;

            // Extract the index of the two nodes a and b
            size_type a_index = a.index();
            size_type b_index = b.index();

            // Declare the unordered map used to map the node indeces connected
            // to a given node and the corresponnding edges
            std::unordered_map<size_type, internal_edge> connected_nodes;

            // Test if it is a new edge. If the edge is not new then it looks 
            // for its index. If the edge is new then it is added to the map
            // of edges
            bool is_in_edges = has_edge(a, b);
            if(is_in_edges){
                if(edges_to_index_.find(a_index) != edges_to_index_.end()){
                    connected_nodes = edges_to_index_.at(a_index);
                    new_edge_index = connected_nodes.at(b_index).index;
                    
                } else if(edges_to_index_.find(b_index) != edges_to_index_.end()){
                    connected_nodes = edges_to_index_.at(b_index);
                    new_edge_index = connected_nodes.at(a_index).index;
                }
            } else {
                new_edge_index = next_edge_index_;
                index_to_edges_[new_edge_index] = {a, b, new_edge_index};

                if(edges_to_index_.find(a_index) != edges_to_index_.end()){
                    std::pair<size_type, internal_edge> new_connection 
                        (b_index, {a, b, new_edge_index});
                    edges_to_index_[a_index].insert(new_connection);
                } else {
                    std::pair<size_type, internal_edge> new_connection 
                        (b_index, {a, b, new_edge_index});
                    connected_nodes.insert(new_connection);
                    std::pair<size_type, std::unordered_map<size_type, internal_edge>> 
                        new_edge (a_index, connected_nodes);
                    edges_to_index_.insert(new_edge);
                }
                ++next_edge_index_;
                ++num_edge_;
            }

            return Edge(this, new_edge_index);
        }

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

            // Clear Edge data
            num_edge_ = 0;
            next_edge_index_ = 0;
            index_to_edges_ = {};
            edges_to_index_ = {};
            
        }

 private:

    // Internal type for node
    struct internal_node {
        Point point;  
        size_type index;    
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

    // Number of nodes in the graph
    size_type size_;

    // Lastest index given to a node in the graph
    size_type next_node_index_;

    // Number of distinct edges in the graph
    size_type num_edge_;

    // Lastest index given to an edge in the graph
    size_type next_edge_index_;


};

#endif // CME212_GRAPH_HPP
