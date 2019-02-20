#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
    
private:
    class internalNode;
    class internalEdge;
    
    
public:
    
    //
    // PUBLIC TYPE DEFINITIONS
    //
    
    typedef V node_value_type;
    
    typedef E edge_value_type;
    
    
    /** Type of this graph. */
    using graph_type = Graph<V, E>;
    
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
        
        //invalid node
        Node()
        {
            
            
        }
        
        
        /** Return this node's position. */
        
        Point& position()
        {
            //We get an internalNode object from the vector vectorOfNodes
            //and we get the Point object.
            return (graph_ -> vectorOfNodes[identificationNumberNode].pointOfTheNode);
        }
        
        /** Return this node's position. */
        
        const Point& position() const {
            
            //We get an internalNode object from the vector vectorOfNodes
            //and we get the Point object.
            return (graph_ -> vectorOfNodes[identificationNumberNode].pointOfTheNode);
        }
        
        
        /** Return this node's index, a number in the range [0, graph_size). */
        
        size_type index() const
        {
            internalNode intNode = graph_ -> vectorOfNodes[identificationNumberNode];
            return intNode.nodeIndex;
        }
        
        /** This function returns the value of a node.
         * @return value of a node.
         *
         * Complexity: O(1).
         */
        
        
        node_value_type& value()
        {
            return (graph_->vectorOfNodes[identificationNumberNode].nodeValue);
        }
        
        /** This function returns the value of a node.
         * @return value of a node.
         *
         * Complexity: O(1).
         */
        
        const node_value_type& value() const
        {
            return (graph_->vectorOfNodes[identificationNumberNode].nodeValue);
        }
        
        /** This function returns the degree of a node.
         * @return degree of node.
         *
         * Complexity: O(1).
         */
        
        size_type degree() const
        {
            return graph_-> adjacencyList.at(identificationNumberNode).size();
        }
        
        /** This function returns an iterator that points to the first incident edge of a node.
         * @return incident_iterator first incident edge of a node.
         *
         * Complexity: O(1).
         */
        
        incident_iterator edge_begin() const
        {
            return IncidentIterator(graph_, (graph_ -> adjacencyList.at(identificationNumberNode)).begin(), identificationNumberNode);
        }
        
        /** This function returns an iterator that points to the last incident edge of a node.
         * @return incident_iterator last incident edge of a node.
         *
         * Complexity: O(1).
         */
        
        incident_iterator edge_end() const
        {
            return IncidentIterator(graph_, (graph_ -> adjacencyList.at(identificationNumberNode)).end(), identificationNumberNode);
        }
        
        
        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        
        bool operator==(const Node& n) const {
            
            return ((graph_ == n.graph_) && (identificationNumberNode == n.identificationNumberNode));
            
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
            
            //We compare the index of these 2 nodes to check which one is smaller.
            return (identificationNumberNode < n.identificationNumberNode);
        }
        
    private:
        
        // Allow Graph to access Node's private member data and functions.
        
        //We define a pointer to the graph and an ID for the node that represents
        // the index for now.
        //We also define a private constructor that takes as parameter the pointer
        // to the graph and the ID of the nodes to be able to define a Node object.
        
        friend class Graph;
        
        Graph* graph_;
        size_type identificationNumberNode;
        
        /** Private Constructor */
        Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), identificationNumberNode(uid) {
        }
        
    };
    
    /** Return the number of nodes in the graph.
     *
     * Complexity: O(1).
     */
    
    size_type size() const {
        
        return indexToID_Node.size();
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
    
    Node add_node(const Point& position, const node_value_type& value = node_value_type())
    {
        //Creating new index and ID.
        size_type nextIndexNode = indexToID_Node.size();
        size_type nextUidNode = vectorOfNodes.size();
        
        // Creating a new Node pointing to the graph with the new ID.
        Node newNode = Node(this, nextUidNode);
        
        // Creating a new Internal Node with the same attributes and we push it to the vector
        // of internal nodes.
        internalNode newInternalNode = internalNode(this, position, value, nextUidNode, nextIndexNode, true);
        
        vectorOfNodes.push_back(const_cast<internalNode&>(newInternalNode));
        indexToID_Node.push_back(nextUidNode);
        
        // We add an empty map to ensure that the indexing is correct.
        adjacencyList[nextUidNode] = std::map<size_type, size_type>();
        
        return newNode;
    }
    
    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    
    bool has_node(const Node& n) const {
        
        return (this == n.graph_);
    }
    
    
    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    
    Node node(size_type i) const
    {
        if (i < num_nodes())
        {
            size_type idNode = indexToID_Node[i];
            return Node(this, idNode);
        }
        else
        {
            throw std::invalid_argument("This is an invalid node");
        }
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
        Edge()
        {
            
        }
        
        /** Return the length of an Edge using the Euclidean norm */
        
        double length() const
        {
            return norm(node1().position() - node2().position());
        }
        
        /** Return a node of this Edge */
        
        Node node1() const {
            
            return Node(graph_, idFirstNode);
        }
        
        /** Return the other node of this Edge */
        
        Node node2() const {

            return Node(graph_, idSecondNode);
        }
        
        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        
        bool operator==(const Edge& e) const
        {
            //I'm checking whether the first node of the first edge is
            // the same as the first node of the second edge and whether
            // the second node of the first edge is the same as the second node
            // of the second edge (and in the opposite order).
            
            bool firstCond = (node1() == e.node1()) && (node2() == e.node2());
            bool secondCond = (node1() == e.node2()) && (node2() == e.node1());
            return (firstCond or secondCond);
            
        }
        
        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        
        
        
        bool operator<(const Edge& e) const
        {
            // We compare the graphs of two edges
            if (graph_ != e.graph_)
            {
                return (graph_ < e.graph_);
            }
            // We just compare the indices of these 2 edges and check which one
            // is larger.
            return (identificationNumberEdge < e.identificationNumberEdge);
            
        }
        
        /** Return the value assigned by the edge */
        
        edge_value_type& value()
        {
            return (graph_->vectorOfEdges[identificationNumberEdge].edgeValue);
        }
        
        /** Return the value assigned by the edge */
        
        const edge_value_type& value() const
        {
            return (graph_->vectorOfEdges[identificationNumberEdge].edgeValue);
        }
        
        /** Return the index of the edge, a number in the range [0, num_edges())*/
        
        size_type index()
        {
            internalEdge intEdge = graph_->vectorOfEdges[identificationNumberEdge];
            return intEdge.edgeIndex;
        }
    private:
        // Allow Graph to access Edge's private member data and functions.
        
        // We define the private data members that we'll be using which are a graph
        // pointer, the ID of the edge, the ID of the first Node and the ID of
        // the second Edge.
        // We also define a constructor that takes the private data members that I
        // defined as parameters.
        
        friend class Graph;
        
        Graph* graph_;
        size_type identificationNumberEdge;
        size_type idFirstNode;
        size_type idSecondNode;
        
        /** Private Constructor */
        Edge(const Graph* graph, size_type uid, size_type node1, size_type node2)
        : graph_(const_cast<Graph*>(graph)), identificationNumberEdge(uid), idFirstNode(node1), idSecondNode(node2) {
        }
    };
    
    /** Return the total number of edges in the graph.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    
    size_type num_edges() const {
        
        return indexToID_Edge.size();
    }
    
    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */

    
    Edge edge(size_type i) const
    {
        if (i < num_edges())
        {
            
            size_type idOfEdge = indexToID_Edge[i];
            
            // Get internal Edge from the vector of Edges at position idOfEdge.
            internalEdge nodesOfThisEdge = vectorOfEdges[idOfEdge];
            
            // Create an edge using the index idOfEdge and the indices of nodes which
            // we got from the internal edge.
            
            Edge e (this, idOfEdge, nodesOfThisEdge.iD_Of_firstNode, nodesOfThisEdge.iD_Of_secondNode);
            return e;
        }
        else
        {
            throw std::invalid_argument("This is an invalid edge.");
        }
    }
    
    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */

    
    bool has_edge(const Node& a, const Node& b) const
    {
        // We first check that nodes a and b are valid nodes of this graph.
        
        if (has_node(a) && has_node(b))
        {
            size_type iDa = a.identificationNumberNode;
            size_type iDb = b.identificationNumberNode;
            
            // We check if th edge formed by nodes a and b exists, i.e it is in the
            // adjacency list.
            
            if (adjacencyList.count(iDa))
            {
                return (adjacencyList.at(iDa).count(iDb));
            }
            else
            {
                return false;
            }
        }
        else
        {
            throw std::invalid_argument("These nodes are not valid");
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

    
    Edge add_edge(const Node& a, const Node& b)
    {
        
        size_type nextIndexEdge = indexToID_Edge.size();
        size_type nextUidEdge = vectorOfEdges.size();
        
        // We first check that a and b are distinct valid nodes of this graph.
        if (has_node(a) && has_node(b) && (not (a == b)))
        {
            // if we have a new edge
            if (has_edge(a, b) == false)
            {
                //when updating the adjacency list, we need to fill both
                //edges (a,b) and (b.a).
                adjacencyList[a.identificationNumberNode][b.identificationNumberNode] = nextUidEdge;
                adjacencyList[b.identificationNumberNode][a.identificationNumberNode] = nextUidEdge;
                
                // we push this edge to the vector of edges.
                internalEdge interEdge = internalEdge(this, nextUidEdge, a.identificationNumberNode, b.identificationNumberNode, nextIndexEdge, true);
                vectorOfEdges.push_back(interEdge);
                indexToID_Edge.push_back(nextUidEdge);
                
                Edge e = Edge(this, nextUidEdge, a.index(), b.index());
                return e;
            }
            // edge already existed
            else
            {
                size_type edgeID = adjacencyList.at(a.identificationNumberNode).at(b.identificationNumberNode);
                return Edge(this, edgeID, a.identificationNumberNode, b.identificationNumberNode);
            }
        }
        else
        {
            throw std::invalid_argument("These nodes are either the same or not valid");
        }
    }
    
    /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    
    void clear()
    {
        vectorOfNodes = std::vector<internalNode>();
        adjacencyList = std::map<size_type, std::map<size_type, size_type>>();
        vectorOfEdges = std::vector<internalEdge>();
        indexToID_Node = std::vector<size_type>();
        indexToID_Edge = std::vector<size_type>();
        
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
        using pointer           = Node*;            // Pointers to elements
        using reference         = Node&;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
        
        /** Construct an invalid NodeIterator. */
        NodeIterator()
        {
            
        }
        
        
        /** This function returns a node dereferenced by the node iterator.
         * @return Node for the specific node iterator.
         *
         * Complexity: O(1).
         */
        
        Node operator*() const
        {
            size_type idnode = *iterator;
            return Node(graph_, idnode);
        }
        
        /** This function returns a reference to the node iterator that was incremented by 1.
         * @return NodeIterator incremented by 1.
         *
         * Complexity: O(1).
         */
        
        NodeIterator& operator++()
        {
            iterator++;
            return *this;
        }
        
        /** This function checks whether two iterators are the same.
         * @param[in] newNodeIterator This is another node iterator that will be
         * compared to the one defined in the private attributes of the class NodeIterator.
         * @return true if the iterators are the same and false if they aren't.
         *
         * Complexity: O(1).
         */
        
        bool operator==(const NodeIterator& newNodeIterator) const
        {
            return (iterator == newNodeIterator.iterator);
        }

    private:
        friend class Graph;
        typename std::vector<size_type>::const_iterator iterator;
        const Graph* graph_;
        
        //Private constructor
        NodeIterator(typename std::vector<size_type>::const_iterator iter, const Graph* graph)
        :  iterator(iter), graph_(graph) {
        }
    };

    /** This function returns an iterator of the first node.
     * @return node_iterator first node iterator.
     *
     * Complexity: O(1).
     */

    node_iterator node_begin() const
    {
        NodeIterator iteratorBegin = NodeIterator(indexToID_Node.begin(), this);
        return iteratorBegin;
    }

    /** This function returns an iterator of the last node.
     * @return node_iterator last node iterator.
     *
     * Complexity: O(1).
     */

    node_iterator node_end() const
    {
        return NodeIterator(indexToID_Node.end(), this);
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
        IncidentIterator()
        {
            
        }
        
        /** This function returns an Edge dereferenced by the incident iterator, it is
         * the edge incident to a node.
         * @return Edge for the specific node iterator.
         *
         * Complexity: O(1).
         */
        
        Edge operator*() const
        {
            return Edge(graph_, iterator->second, nodeA, iterator->first);
        }
        
        /** This function returns a reference to the incident iterator that was incremented by 1.
         * @return incident_iterator incremented by 1.
         *
         * Complexity: O(1).
         */
        
        incident_iterator& operator++()
        {
            iterator++;
            return *this;
        }
        
        /** This function checks whether two iterators are the same.
         * @param[in] newIncidentiter This is another incident iterator that will be
         * compared to the one defined in the private attributes of the class IncidentIterator.
         * @return true if the iterators are the same and false if they aren't.
         *
         * Complexity: O(1).
         */
        
        bool operator==(const incident_iterator& newIncidentiter) const
        {
            return (iterator == newIncidentiter.iterator);
        }
        
    private:
        friend class Graph ;
        // HW1 #3: YOUR CODE HERE
        Graph* graph_;
        typename std::map<size_type, size_type>::const_iterator iterator;
        size_type nodeA;
        
        //Private Constructor
        IncidentIterator(const Graph* graph, typename std::map<size_type, size_type>::const_iterator iter, size_type nodeOfA)
        : graph_(const_cast<Graph*>(graph)), iterator(iter), nodeA(nodeOfA){
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
        
        /** This function returns an Edge by getting the internal Edge
         * at the specific iD.
         * @return Edge for the specific edge iterator.
         *
         * Complexity: O(1).
         */
        
        Edge operator*() const
        {
            size_type idedge = *iteratorEdge;
            internalEdge iterEdge = graph_ -> vectorOfEdges[idedge];
            return Edge(graph_, iterEdge.iD_Of_Edge, iterEdge.iD_Of_firstNode, iterEdge.iD_Of_secondNode);
        }
        
        /** This function returns a reference to the Edge iterator that was incremented by 1.
         * @return Edgeiterator incremented by 1.
         *
         * Complexity: O(1).
         */
        
        EdgeIterator& operator++()
        {
            iteratorEdge++;
            return *this;
        }
        
        /** This function checks whether two iterators are the same.
         * @param[in] newEdgeIterator This is another incident iterator that will be
         * compared to the one defined in the private attributes of the class EdgeIterator.
         * @return true if the iterators are the same and false if they aren't.
         *
         * Complexity: O(1).
         */
        
        bool operator==(const EdgeIterator& newEdgeIterator) const
        {
            return (iteratorEdge == newEdgeIterator.iteratorEdge);
        }
        
    private:
        friend class Graph;
        const Graph* graph_;
        typename std::vector<size_type>::const_iterator iteratorEdge;
        
        EdgeIterator(const Graph* graph, typename std::vector<size_type>::const_iterator iter)
        : graph_(const_cast<Graph*>(graph)), iteratorEdge(iter) {
        }
    };
    
    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    /** This function returns an iterator of the first edge.
     * @return node_iterator first edge iterator.
     *
     * Complexity: O(1).
     */
    
    edge_iterator edge_begin() const
    {
        return EdgeIterator(this, indexToID_Edge.begin());
    }
    
    /** This function returns an iterator of the last edge.
     * @return edge_iterator last edge iterator.
     *
     * Complexity: O(1).
     */
    
    edge_iterator edge_end() const
    {
        return EdgeIterator(this, indexToID_Edge.end());
    }
    
    
    /** Remove a node from the graph.
     * @param[in] A Node object passed by reference called node which is a node of this graph.
     * @return This function returns 0 if the node doesn't exist in the graph
     * and 1 if it does and got removed from the graph.
     * @post num_nodes() has decreased by 1.
     *
     * Complexity: O(num_nodes() + num_edges())
     */
    
    
    size_type remove_node(const Node& node)
    {
        if (not has_node(node))
        {
            return 0;
        }
        
        // When removing a node, some edges won't exist anymore
        // so we have to remove the edges that are incident to
        // this node.
        while (node.edge_begin() != node.edge_end())
        {
            Edge edge = *node.edge_begin();
            remove_edge(edge);
            
        }
            
        //This node is invalid.
        vectorOfNodes[node.identificationNumberNode].validInternalNode = false;
            
        //Swap the node we want to remove with the last node in the vector.
        indexToID_Node[node.index()] = indexToID_Node[indexToID_Node.size() - 1];
            
        // We set the index of the new node at the position of the deleted node to
        // the index of the deleted node.
        vectorOfNodes[indexToID_Node[node.index()]].nodeIndex = node.index();
            
        // We pop the last element of the vector which is the one we want to remove.
        indexToID_Node.pop_back();
            
        return 1;
            
    }
    
    
    /** Remove a node from the graph.
     * @param[in] A node_iterator n_it which is a node iterator of the node
     * we want to remove from the graph.
     * @return This function returns a node_iterator pointing to the first Node.
     * @post num_nodes() has decreased by 1.
     *
     * Complexity: O(num_nodes() + num_edges())
     */
    
    node_iterator remove_node(node_iterator n_it)
    {
        Node node = *n_it;
        remove_node(node);
        
        return node_begin();

    }
    
    
    /** Remove an edge from the graph.
     * @param[in] A Node object passed by reference called node1 which is the first node
     * forming the edge.
     * @param[in] A Node object passed by reference called node2 which is the first node
     * forming the edge.
     * @return This function returns 0 if the edge doesn't exist in the graph
     * and 1 if it does and got removed from the graph.
     * @post num_edges() has decreased by 1.
     *
     * Complexity: O(num_nodes() + num_edges())
     */
    
    
    size_type remove_edge(const Node& node1, const Node& node2)
    {
        // The edge doesn't exist.
        if (!has_edge(node1, node2))
        {
            return 0;
        }
        size_type edgeID = (adjacencyList.at(node1.identificationNumberNode)).at(node2.identificationNumberNode);
        Edge e = Edge(this, edgeID, node1.identificationNumberNode, node2.identificationNumberNode);
        
        //The edge that we want to remove is now invalid.
        vectorOfEdges[edgeID].validInternalEdge = false;
        
        //Swapping the edge that we want to remove with the last edge in the vector.
        indexToID_Edge[e.index()] = indexToID_Edge[indexToID_Edge.size() - 1];
        
        // We set the index of the new edge at the position of the deleted edge to
        // the index of the deleted edge.
        vectorOfEdges[indexToID_Edge[e.index()]].edgeIndex = e.index();
        
        // We pop the last element of the vector
        indexToID_Edge.pop_back();
        
        //We update our adjacency list by erasing the adjacent nodes for every node.
        
        adjacencyList.at(node1.identificationNumberNode).erase(node2.identificationNumberNode);
        adjacencyList.at(node2.identificationNumberNode).erase(node1.identificationNumberNode);
        
        return 1;
        
    }
    
    
    /** Remove an edge from the graph.
     * @param[in] A Edge object passed by reference called e which is the
     * first node forming the edge.
     * @return This function returns 0 if the edge doesn't exist in the graph
     * and 1 if it does and got removed from the graph.
     * @post num_edges() has decreased by 1.
     *
     * Complexity: O(num_nodes() + num_edges())
     */
    
    size_type remove_edge(const Edge& e)
    {
        Node firstNode = e.node1();
        Node secondNode = e.node2();
        
        return remove_edge(firstNode, secondNode);
    }
    
    
    /** Remove an edge from the graph.
     * @param[in] An edge_iterator n_it which is a edge iterator of the edge
     * we want to remove from the graph.
     * @return This function returns an edge_iterator pointing to the first Edge.
     * @post num_edges() has decreased by 1.
     *
     * Complexity: O(num_nodes() + num_edges())
     */
    
    
    edge_iterator remove_edge(edge_iterator e_it)
    {
        Edge e = *e_it;
        remove_edge(e);
        return edge_begin();
    }
    
    
private:
    
    class internalNode
    {
        friend class Graph;
        
        // Parameters of class internalNode : Point of the Node, value of the node
        // and the index of the node.
        Graph* graph_;
        Point pointOfTheNode;
        node_value_type nodeValue;
        size_type nodeID;
        size_type nodeIndex;
        bool validInternalNode;
        
        /** Construct an invalid internalNode. */
        internalNode()
        {
            
        }
        
        //Private Constructor
        internalNode(const Graph* graph, Point pointNode, node_value_type valueOfNode, size_type idOfNode, size_type indexNode, bool nodeValid)
        : graph_(const_cast<Graph*>(graph)), pointOfTheNode(pointNode), nodeValue(valueOfNode), nodeID(idOfNode), nodeIndex(indexNode), validInternalNode(nodeValid) {
        }
    };
    
    // This class has the object internalEdge which is defined by the
    // index of the edge, the index of the first node and the index of
    // the second node.
    // We also create an invalid internalEdge constructor and a valid
    // one which takes the index of the edge and the index of the nodes
    // as parameters.
    
    class internalEdge
    {
        friend class Graph;
        
        // Parameters of class internalEdge : index of the Edge, index of the first node
        // and index of the second node.
        Graph* graph_;
        size_type iD_Of_Edge;
        size_type iD_Of_firstNode;
        size_type iD_Of_secondNode;
        edge_value_type edgeValue;
        size_type edgeIndex;
        bool validInternalEdge;
        
        /** Construct an invalid internalEdge. */
        internalEdge()
        {
            
        }
        //Private Constructor
        internalEdge(const Graph* graph, size_type edgeId, size_type idNode1, size_type idNode2, size_type indexEdge, bool edgeValid)
        : graph_(const_cast<Graph*>(graph)), iD_Of_Edge(edgeId), iD_Of_firstNode(idNode1), iD_Of_secondNode(idNode2), edgeIndex(indexEdge), validInternalEdge(edgeValid) {
        }
        
        
    };
    
    // vector of internal Nodes
    std::vector<internalNode> vectorOfNodes;
    
    // vector of internal Edges
    std::vector<internalEdge> vectorOfEdges;
    
    // vector of incides to uid's for Nodes
    std::vector<size_type> indexToID_Node;
    
    // vector of incides to uid's for Edges
    std::vector<size_type> indexToID_Edge;
    
    
    // map that has the indices of two nodes and the index of the edges
    // joins them
    std::map<size_type, std::map<size_type, size_type>> adjacencyList;
};

#endif // CME212_GRAPH_HPP
