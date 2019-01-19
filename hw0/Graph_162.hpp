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
class Graph {
private:
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
    Graph() {
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

        }

        /** Return this node's position. */
        const Point& position() const {
            return this->graph->pointVector[indexNode];
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
            return indexNode;
        }

        /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
        bool operator==(const Node& n) const {
            if (this->graph==n.graph && this->indexNode == n.indexNode){
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
        bool operator<(const Node& n) const {
            if (this->index()<n.index()){
                return true;
            }
            return false;
        }

    private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;
        // Use this space to declare private data members and methods for Node
        // that will not be visible to users, but may be useful within Graph.
        // i.e. Graph needs a way to construct valid Node objects
        Graph* graph;
        int indexNode;
        // private constructor used to add node in the graph
        Node(const Graph* graphInput,int index) : graph(const_cast<Graph*>(graphInput)), indexNode(index) {
        }
    };

    /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
    size_type size() const {
        return pointVector.size();
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
        // create new node and add to graph
        Node newNode(this,pointVector.size());
        pointVector.push_back(position);
        return newNode;
    }

    /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
    bool has_node(const Node& n) const {
        if (this == n.graph){
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
    Node node(size_type i) const {
        Node newNode(this,i);
        return newNode;
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
        }

        /** Return a node of this Edge */
        Node node1() const {
            // look at the edge variable and get the node index to fetch the corresponding node
            std::vector<size_type> ev = *edge;
            int node1Index=ev[0];
            Node newNode(graph,node1Index);
            return newNode;
        }

        /** Return the other node of this Edge */
        Node node2() const {
            // look at the edge variable and get the node index to fetch the corresponding node
            std::vector<size_type> ev = *edge;
            int node2Index=ev[1];
            Node newNode(graph,node2Index);
            return newNode;
        }

        /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
        bool operator==(const Edge& e) const {
            // check if the edge has the same pair of node, the order of the pair does not have to be the same
            int firstIndex = (*this->edge)[0];
            int secondIndex = (*this->edge)[1];
            int firstIndexE = (*e.edge)[0];
            int secondIndexE = (*e.edge)[1];
            bool sameOrderEqual = ((this->graph->pointVector)[firstIndex] == (e.graph->pointVector)[firstIndexE]) && ((this->graph->pointVector)[secondIndex] == (e.graph->pointVector)[secondIndexE]);
            bool flipOrderEqual = ((this->graph->pointVector)[firstIndex] == (e.graph->pointVector)[secondIndexE]) && ((this->graph->pointVector)[secondIndex] == (e.graph->pointVector)[firstIndexE]);
            if (sameOrderEqual || flipOrderEqual){
                return true;
            }
            return false;
        }

        /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
        bool operator<(const Edge& e) const {
            if (this->indexEdge < e.indexEdge){
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
        Graph* graph;
        std::vector<size_type>* edge;
        int indexEdge;
        // private constructor used to add edge to graph
        Edge(const Graph* graphInput,const std::vector<size_type>* edgeInput, int index) : graph(const_cast<Graph*>(graphInput)),edge(const_cast<std::vector<size_type>*>(edgeInput)),indexEdge(index) {
        }
    };

    /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
    size_type num_edges() const {
        return edgeVector.size();
    }

    /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
    Edge edge(size_type i) const {
        Edge newEdge(this,&(edgeVector[i]),i);
        return newEdge;
    }

    /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
    bool has_edge(const Node& a, const Node& b) const {
        // check if both point to the same graph
        assert(a.graph==b.graph);
        // check if the pair is in the map
        std::vector<size_type> edge1={a.index(),b.index()};
        std::vector<size_type> edge2={b.index(),a.index()};
        if (edgeMap.find(edge1) != edgeMap.end() || edgeMap.find(edge2) != edgeMap.end()){
            return true;
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
        std::vector<size_type> edge={a.index(),b.index()};
        if (has_edge(a,b)){
            int index=edgeMap[edge];
            return Edge(this,&edge,index);
        } else {
            int newSize=edgeVector.size();
            edgeMap[edge]=newSize;
            edgeVector.push_back(edge);
            return Edge(this,&edge,newSize);
        }
    }

    /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
    void clear() {
        pointVector.clear();
        edgeVector.clear();
        edgeMap.clear();
    }

private:
    // Use this space for your Graph class's internals:
    //   helper functions, data members, and so forth.
    std::vector<Point> pointVector{};
    std::vector<std::vector<size_type>> edgeVector{};
    // construct a map from a pair of nodes to edge indes in order to check if the edge already existed in the Graph class
    std::map<std::vector<size_type>, size_type> edgeMap{};
};

#endif // CME212_GRAPH_HPP
