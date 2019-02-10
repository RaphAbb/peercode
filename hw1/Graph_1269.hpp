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
template <typename V>
class Graph {
    
private:
    class internalNode;
    class internalEdge;
    
    
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
    using size_type = unsigned;
    
    //
    // CONSTRUCTORS AND DESTRUCTOR
    //
    
    /** Construct an empty graph. */
    Graph() {
        next_uid_Node = 0;
        next_uid_Edge = 0;
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
        
        const Point& position() const {
            
            //We get an internalNode object from the vector nodes_In_The_Graph
            //and we get the Point object.
            return graph_ -> nodes_In_The_Graph[identificationNumberNode].pointOfTheNode;
        }
        
        /** Return this node's index, a number in the range [0, graph_size). */
        
        
        size_type index() const
        {
            return identificationNumberNode;
        }
        
        /** This function returns the value of a node.
         * @return value of a node.
         *
         * Complexity: O(1).
         */
        
        
        node_value_type& value()
        {
            return (graph_->nodes_In_The_Graph[identificationNumberNode].nodeValue);
        }
        
        /** This function returns the value of a node.
         * @return value of a node.
         *
         * Complexity: O(1).
         */
        
        const node_value_type& value() const
        {
            return (graph_->nodes_In_The_Graph[identificationNumberNode].nodeValue);
        }
        
        /** This function returns the degree of a node.
         * @return degree of node.
         *
         * Complexity: O(1).
         */
        
        size_type degree() const
        {
            return graph_->adjacencyList.at(identificationNumberNode).size();
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
            
            if ((graph_ == n.graph_) && (identificationNumberNode == n.index()))
            {
                return true;
            }
            else
            {
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
            
            //We compare the index of these 2 nodes to check which one is smaller.
            if (identificationNumberNode < n.index())
            {
                return true;
            }
            else
            {
                return false;
            }
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
        
        return nodes_In_The_Graph.size();
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
        // Creating a new Node pointing to the graph with the new ID.
        Node newNode = Node(this, next_uid_Node);
        
        // Creating a new Internal Node with the same attributes and we push it to the vector
        // of internal nodes.
        internalNode newInternalNode = internalNode(position, value, next_uid_Node);
        nodes_In_The_Graph.push_back(newInternalNode);
        
        // We add an empty map to ensure that the indexing is correct.
        adjacencyList[next_uid_Node] = std::map<size_type, size_type>();
        
        //We then increment the index to get valid results when adding the next node.
        next_uid_Node = next_uid_Node + 1;
        
        return newNode;
    }
    
    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    
    bool has_node(const Node& n) const {
        internalNode internalNodeOfn = nodes_In_The_Graph[n.index()];
        if (n.position() == internalNodeOfn.pointOfTheNode)
        {
            return true;
        }
        else
        {
            return false;
        }
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
            return Node(this, i);
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
        
        /** Return a node of this Edge */
        
        Node node1() const {
            
            return Node(graph_, firstNode);
        }
        
        /** Return the other node of this Edge */
        
        Node node2() const {
            
            return Node(graph_, secondNode);
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
            // of the second edge.
            
            if (((firstNode == (e.node1()).index()) && (secondNode == (e.node2()).index())) or ((firstNode == (e.node2()).index()) && (secondNode == (e.node1()).index())))
            {
                return true;
            }
            else
            {
                return false;
            }
            
        }
        
        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        
        
        
        bool operator<(const Edge& e) const
        {
            // We just compare the indices of these 2 edges and check which one
            // is larger.
            
            if (identificationNumberEdge < e.identificationNumberEdge)
            {
                return true;
            }
            else
            {
                return false;
            }
            
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
        size_type firstNode;
        size_type secondNode;
        
        /** Private Constructor */
        Edge(const Graph* graph, size_type uid, size_type node1, size_type node2)
        : graph_(const_cast<Graph*>(graph)), identificationNumberEdge(uid), firstNode(node1), secondNode(node2) {
        }
    };
    
    /** Return the total number of edges in the graph.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    
    size_type num_edges() const {
        
        return vectorOfEdges.size();
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
            // Get internal Edge from the vector of Edges at position i.
            internalEdge nodesOfThisEdge = vectorOfEdges[i];
            
            // Create an edge using the index i and the indices of nodes which
            // we got from the internal edge.
            
            Edge e (this, i, nodesOfThisEdge.index_Of_firstNode, nodesOfThisEdge.index_Of_secondNode);
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
            size_type indexOfa = a.index();
            size_type indexOfb = b.index();
            
            // We check if th edge formed by nodes a and b exists, i.e it is in the
            // adjacency list.
            
            if (adjacencyList.count(indexOfa) > 0)
            {
                if (adjacencyList.at(indexOfa).count(indexOfb) > 0)
                {
                    return true;
                }
                else
                {
                    return false;
                }
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
        // We first check that a and b are distinct valid nodes of this graph.
        
        if (has_node(a) && has_node(b) && (not (a == b)))
        {
            // if we have a new edge
            if (has_edge(a, b) == false)
            {
                //when updating the adjacency list, we need to fill both
                //edges (a,b) and (b.a).
                adjacencyList[a.index()][b.index()] = next_uid_Edge;
                adjacencyList[b.index()][a.index()] = next_uid_Edge;
                
                // we push this edge to the vector of edges.
                vectorOfEdges.push_back(internalEdge(next_uid_Edge, a.index(), b.index()));
                Edge e = Edge(this, next_uid_Edge, a.index(), b.index());
                next_uid_Edge = next_uid_Edge + 1;
                return e;
            }
            // edge already existed
            else
            {
                size_type edgeID = adjacencyList.at(a.index()).at(b.index());
                return Edge(this, edgeID, a.index(), b.index());
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
        next_uid_Node = 0;
        next_uid_Edge = 0;
        nodes_In_The_Graph = std::vector<internalNode>();
        adjacencyList = std::map<size_type, std::map<size_type, size_type>>();
        vectorOfEdges = std::vector<internalEdge>() ;
        
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
            return Node(graph_, (*iteratorNode).nodeID);
        }
        
        /** This function returns a reference to the node iterator that was incremented by 1.
         * @return NodeIterator incremented by 1.
         *
         * Complexity: O(1).
         */
        
        NodeIterator& operator++()
        {
            iteratorNode++;
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
            if ((graph_ == newNodeIterator.graph_) && (iteratorNode == newNodeIterator.iteratorNode))
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        
    private:
        friend class Graph;
        
        Graph* graph_;
        typename std::vector<internalNode>::const_iterator iteratorNode;
        
        //Private constructor
        NodeIterator(const Graph* graph, typename std::vector<internalNode>::const_iterator pointerNodeInternal)
        : graph_(const_cast<Graph*>(graph)), iteratorNode(pointerNodeInternal) {
        }
    };

    /** This function returns an iterator of the first node.
     * @return node_iterator first node iterator.
     *
     * Complexity: O(1).
     */

    node_iterator node_begin() const
    {
        return NodeIterator(this, nodes_In_The_Graph.begin());
    }

    /** This function returns an iterator of the last node.
     * @return node_iterator last node iterator.
     *
     * Complexity: O(1).
     */

    node_iterator node_end() const
    {
        return NodeIterator(this, nodes_In_The_Graph.end());
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
            if ((graph_ == newIncidentiter.graph_) && (iterator == newIncidentiter.iterator))
            {
                return true;
            }
            else
            {
                return false;
            }
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
        
        /** This function returns an Edge dereferenced by the edge iterator, it is
         * @return Edge for the specific edge iterator.
         *
         * Complexity: O(1).
         */
        
        Edge operator*() const
        {
            return Edge(graph_, (*iteratorEdge).index_Of_Edge, (*iteratorEdge).index_Of_firstNode, (*iteratorEdge).index_Of_secondNode );
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
            if ((graph_ == newEdgeIterator.graph_) && (iteratorEdge == newEdgeIterator.iteratorEdge))
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        
    private:
        friend class Graph;
        Graph* graph_;
        typename std::vector<internalEdge>::const_iterator iteratorEdge;
        
        EdgeIterator(const Graph* graph, typename std::vector<internalEdge>::const_iterator pointerEdgeInternal)
        : graph_(const_cast<Graph*>(graph)), iteratorEdge(pointerEdgeInternal) {
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
        return EdgeIterator(this, vectorOfEdges.begin());
    }
    
    /** This function returns an iterator of the last edge.
     * @return edge_iterator last edge iterator.
     *
     * Complexity: O(1).
     */
    
    edge_iterator edge_end() const
    {
        return EdgeIterator(this, vectorOfEdges.end());
    }
    
    
private:
    
    class internalNode
    {
        friend class Graph;
        
        // Parameters of class internalNode : Point of the Node, value of the node
        // and the index of the node.
        Point pointOfTheNode;
        node_value_type nodeValue;
        size_type nodeID;
        
        /** Construct an invalid internalNode. */
        internalNode()
        {
            
        }
        
        //Private Constructor
        internalNode(Point pointNode, node_value_type valueOfNode, size_type idOfNode)
        : pointOfTheNode(pointNode), nodeValue(valueOfNode), nodeID(idOfNode) {
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
        size_type index_Of_Edge;
        size_type index_Of_firstNode;
        size_type index_Of_secondNode;
        
        /** Construct an invalid internalEdge. */
        internalEdge()
        {
            
        }
        //Private Constructor
        internalEdge(size_type edgeId, size_type idNode1, size_type idNode2)
        : index_Of_Edge(edgeId), index_Of_firstNode(idNode1), index_Of_secondNode(idNode2) {
        }
        
        
    };
    
    // vector of internal Nodes
    std::vector<internalNode> nodes_In_The_Graph;
    
    // new indices of Nodes and Edges
    size_type next_uid_Node;
    size_type next_uid_Edge;
    
    // vector of internal Edges
    std::vector<internalEdge> vectorOfEdges;
    
    // map that has the indices of two nodes and the index of the edges
    // joins them
    std::map<size_type, std::map<size_type, size_type>> adjacencyList;
};

#endif // CME212_GRAPH_HPP
