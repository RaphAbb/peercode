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
#include "CME212/Color.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph {

private:

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

        Node() {

        }


        /** Return this node's position. */
        const Point& position() const {
            return graph->internalNodeVector[indexNode].point;
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
            return indexNode;
        }

        /** return the value associated to the called node */
        node_value_type& value(){
            return graph->internalNodeVector[indexNode].value;
        }

        /** return the value associated to the called node */ //???????????????????????????
        const node_value_type& value() const{
            return graph->internalNodeVector[indexNode].value;
        }

        /** return the number of edge that is incidenting into the called node */
        size_type degree() const{
            return graph->nodeAdjacencyMap[indexNode].size();
        }

        /** return the first edge that is incidenting into the called node */
        incident_iterator edge_begin() const{
            IncidentIterator newiter = IncidentIterator(graph,indexNode,0);
            return newiter;
        }

        /** return the last edge that is incidenting into the called node */
        incident_iterator edge_end() const{
            IncidentIterator newiter = IncidentIterator(graph,indexNode,this->degree());
            return newiter;
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
        return internalNodeVector.size();
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
    
    Node add_node ( const Point & position, const node_value_type & valueInput = node_value_type ()){
        // create new node and add to graph
        Node newNode(this,internalNodeVector.size());
        internalNode newInternalNode({position,valueInput});
        internalNodeVector.push_back(newInternalNode);
        return newNode;
    }

    /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
    bool has_node(const Node& n) const {
        assert(n.index() < this->size());
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
    class Edge : private totally_ordered<Edge>{
    public:
        /** Construct an invalid Edge. */
        Edge() {
        }

        /** Return a node of this Edge */
        Node node1() const {
            // look at the edge variable and get the node index to fetch the corresponding node
            size_type node1Index=graph->edgeVector[indexEdge][0];
            Node newNode(graph,node1Index);
            return newNode;
        }

        /** Return the other node of this Edge */
        Node node2() const {
            // look at the edge variable and get the node index to fetch the corresponding node
            size_type node2Index=graph->edgeVector[indexEdge][1];
            Node newNode(graph,node2Index);
            return newNode;
        }

        /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
        bool operator==(const Edge& e) const {
            // check if the edge has the same pair of node, the order of the pair does not have to be the same
            size_type firstIndex=this->graph->edgeVector[indexEdge][0];
            size_type secondIndex=this->graph->edgeVector[indexEdge][1];
            size_type firstIndexE=e.graph->edgeVector[indexEdge][0];
            size_type secondIndexE=e.graph->edgeVector[indexEdge][1];
            bool sameOrderEqual = ((this->graph->internalNodeVector)[firstIndex].point == (e.graph->internalNodeVector)[firstIndexE].point) && ((this->graph->internalNodeVector)[secondIndex].point == (e.graph->internalNodeVector)[secondIndexE].point);
            bool flipOrderEqual = ((this->graph->internalNodeVector)[firstIndex].point == (e.graph->internalNodeVector)[secondIndexE].point) && ((this->graph->internalNodeVector)[secondIndex].point == (e.graph->internalNodeVector)[firstIndexE].point);
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
        Graph* graph;
        int indexEdge;
        Edge(const Graph* graphInput, int index) : graph(const_cast<Graph*>(graphInput)),indexEdge(index) {
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
        Edge newEdge(this,i);
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
        // sort before insert
        std::vector<size_type> edge, edge2;
        edge={a.index(),b.index()};
        edge2={b.index(),a.index()};
        if (has_edge(a,b)){
            int index=edgeMap[edge];
            return Edge(this,index);
        } else {
            // add node in the adjacency list
            nodeAdjacencyMap[a.index()].push_back(b.index());
            nodeAdjacencyMap[b.index()].push_back(a.index());
            int newSize=edgeVector.size();
            edgeMap[edge]=newSize;
            edgeMap[edge2]=newSize+1;//*
            edgeVector.push_back(edge);
            edgeVector.push_back(edge2);//*

            return Edge(this,newSize);
            //return Edge(this,&edge,newSize);
        }
    }

    /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */

    void clear() {
        internalNodeVector.clear();
        edgeVector.clear();
        edgeMap.clear();
        nodeAdjacencyMap.clear();
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

        /**
         * @brief dereference nodeiterator and return the corresponding node
         * @pre indexNodeIterator is less than the total graph size otherwise invalid node
         */
        // check for derefernecing it hs to be in the range, or just do assert
        Node operator *() const{
            //int index=graph->internalNodeVector[indexNodeIterator].value;
            //Node newNode = Node(graph,index);
            assert(indexNodeIterator < graph->size());
            Node newNode = Node(graph,indexNodeIterator);
            return newNode;
        }

        /**
         * @brief increment node iterator to go to the next node
         * @post if it is not at the end of the list, new indexNodeIterator = old indexNodeIterator + 1
         */
        NodeIterator& operator++() {
            if (indexNodeIterator<graph->size()){
                indexNodeIterator=indexNodeIterator+1;
            }
            return *this;
        }

        /** overloading == operator to compare two nodeIterators */
        bool operator == (const NodeIterator& nodeIterator_) const {
            if ((this->graph == nodeIterator_.graph) && (this->indexNodeIterator == nodeIterator_.indexNodeIterator)){
                return true;
            }
            return false;
        }

    private:
        friend class Graph;
        friend class Node;
        Graph* graph;
        size_type indexNodeIterator;
        // private constructor used to to used in node_begin() and node_end()
        NodeIterator(const Graph* graphInput, int index) : graph(const_cast<Graph*>(graphInput)),indexNodeIterator(index) {
        }
    };

    /** return nodeIterator associated to the first node*/
    NodeIterator node_begin() const {
        NodeIterator n = NodeIterator(this,0);
        return n;
    };

    /** return nodeIterator associated to the last node*/
    NodeIterator node_end() const {
        NodeIterator n = NodeIterator(this,this->size());
        return n;
    };

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_iterator node_begin() const
    // node_iterator node_end() const

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

        /**
         * @brief dereference current egde associated to incident iterator
         * @pre indexNodeIterator is less than the total graph size otherwise invalid node
         * @pre connectedNodeIndexIncidentIterator is less than the total size of the total incident edge of the parent node
         */
        /** dereference current egde associated to incident iterator*/
        Edge operator*() const{
            assert(parentNodeIndexIncidentIterator < graph->size() && connectedNodeIndexIncidentIterator < graph->nodeAdjacencyMap[parentNodeIndexIncidentIterator].size());
            // construct the nodePair from both parent node index and child node index
            size_type node1=parentNodeIndexIncidentIterator;
            size_type node2=graph->nodeAdjacencyMap[parentNodeIndexIncidentIterator][connectedNodeIndexIncidentIterator];
            std::vector<size_type> nodePair;
            nodePair={node1,node2};
            size_type indexEdge=(graph->edgeMap)[nodePair];
            Edge newEdge = Edge(graph, indexEdge);
            return newEdge;
        }

        /**
         * @brief increment incident iterator iterator to go to the next edge
         * @post if it is not at the end of the list, new connectedNodeIndexIncidentIterator = old connectedNodeIndexIncidentIterator + 1
         */
        IncidentIterator& operator++(){
            if (connectedNodeIndexIncidentIterator<graph->nodeAdjacencyMap[parentNodeIndexIncidentIterator].size()){
                this->connectedNodeIndexIncidentIterator++;
            }
            return *this;
        }

        /** Test whether this IncidentIterator and incidentIteratorInput are equal.
        *
        * Equal incidentIterator have the same graph and the same index for both child node and parent node.
        */
        bool operator==(const IncidentIterator& incidentIteratorInput) const{
            if ((this->graph==incidentIteratorInput.graph) &&
                    (this->parentNodeIndexIncidentIterator==incidentIteratorInput.parentNodeIndexIncidentIterator) &&
                    (this->connectedNodeIndexIncidentIterator==incidentIteratorInput.connectedNodeIndexIncidentIterator)){
                return true;
            }
            return false;
        }

    private:
        friend class Graph;
        Graph* graph;
        size_type parentNodeIndexIncidentIterator;
        size_type connectedNodeIndexIncidentIterator;
        // private constructor
        IncidentIterator(const Graph* graphInput,
                         size_type parentNodeIndexIncidentIteratorInput,
                         size_type connectedNodeIndexIncidentIteratorInput) :
            graph(const_cast<Graph*>(graphInput)),
            parentNodeIndexIncidentIterator(parentNodeIndexIncidentIteratorInput)
          ,connectedNodeIndexIncidentIterator(connectedNodeIndexIncidentIteratorInput) {
        }
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

        /**
         * @brief dereference current egde associated to edge iterator
         * @pre indexEdgeIterator is less than the total number of edges in the graph
         */
        Edge operator *() const{
            assert(indexEdgeIterator<graph->edgeVector.size());
            Edge newEdge = Edge(graph,indexEdgeIterator);
            return newEdge;
        }

        /**
         * @brief increment edge iterator iterator to go to the next edge
         * @post if it is not at the end of the list, new indexEdgeIterator = old indexEdgeIterator + 1
         */
        EdgeIterator& operator++() {
            indexEdgeIterator=indexEdgeIterator+1;
            return *this;
        }

        /** Test whether this EdgeIterator and edgeIterator_ are equal.
        *
        * Equal EdgeIterator have the same graph and the same index for the edge.
        */
        bool operator == (const EdgeIterator& edgeIterator_) const {
            if ((this->graph == edgeIterator_.graph) && (this->indexEdgeIterator == edgeIterator_.indexEdgeIterator)){
                return true;
            }
            return false;
        }

    private:
        friend class Graph;
        friend class Edge;
        Graph* graph;
        size_type indexEdgeIterator;
        // private constructor
        EdgeIterator(const Graph* graphInput, int index) : graph(const_cast<Graph*>(graphInput)),indexEdgeIterator(index) {
        }
    };

    /** return EdgeIterator associated to the first edge*/
    EdgeIterator edge_begin() const {
        EdgeIterator n = EdgeIterator(this,0);
        return n;
    }

    /** return EdgeIterator associated to the last edge*/
    EdgeIterator edge_end() const {
        EdgeIterator n = EdgeIterator(this,this->num_edges());
        return n;
    }

private:
    // internal node store node and value associated to each node
    struct internalNode {
        Point point;
        node_value_type value;
    };

    std::vector<internalNode> internalNodeVector{};
    std::vector<std::vector<size_type>> edgeVector{};
    // construct a map from a pair of nodes to edge indes in order to check if the edge already existed in the Graph class
    std::map<std::vector<size_type>, size_type> edgeMap{};
    // create adjacency list to know which edge comes out of a specific node
    std::map<size_type,std::vector<size_type>> nodeAdjacencyMap{};
};

#endif // CME212_GRAPH_HPP
