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
        number_Of_Edges = 0;
        number_Of_Nodes = 0;
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
        
        //invalid node
        Node()
        {
            
    
        }
        
        /** Return this node's position. */
        
        //This function returns an object of type Point which refers
        // to the position of the node. We return the point at the
        // index's position in the vector of points.
        
        const Point& position() const {

            return graph_ -> nodes_In_The_Graph[identificationNumberNode];
        }
        
        /** Return this node's index, a number in the range [0, graph_size). */
        
        //If the node's index is >= 0 or < number of nodes in the graph then
        // this node is valid and we'll return it's index. However, if the node's
        // index is not within that interval then we'll output an error.
        
        size_type index() const
        {

            if ((identificationNumberNode >= 0) && (identificationNumberNode < graph_ -> number_Of_Nodes))
            {
                return identificationNumberNode;
            }
            else
            {
                throw std::invalid_argument("This is an invalid node");
            }
        }
        
        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        
        // We compare whether these nodes are in the same graph and have the same
        // index.
        
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
        
        //We compare the index of these 2 nodes to check which one is smaller.
        
        bool operator<(const Node& n) const {

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
        
        const Graph* graph_;
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
    
    // We return the number of nodes.
    size_type size() const {
       
        return number_Of_Nodes;
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
    
    //We create a newNode pointing to the graph with the new Id and push it
    // into the vector of nodes. We then increment the number of Nodes and
    // the index to get valid results when adding the next node.
    // We add an empty vector to ensure that the indexing is correct.
    
    Node add_node(const Point& position)
    {
        Node newNode = Node(this, next_uid_Node);
        nodes_In_The_Graph.push_back(position);
        number_Of_Nodes = number_Of_Nodes + 1;
        next_uid_Node = next_uid_Node + 1;
        adjacencyList.push_back(std::vector<internalEdge>());
        
        return newNode;
    }
    
    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    
    // We check if the Point characterizing the Node n exists in the vector of
    // of nodes at the index of n.
    
    bool has_node(const Node& n) const {
        
        if (n.position() == nodes_In_The_Graph[n.index()])
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
    
    //We return the Node at index i.
    
    Node node(size_type i) const
    {
        if ((0 <= i) && (i < number_Of_Nodes))
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
    class Edge {
    public:
        /** Construct an invalid Edge. */
        Edge()
        {
    
        }
        
        /** Return a node of this Edge */
        
        //We return the first Node of this Edge.
        
        Node node1() const {

            return Node(graph_, firstNode);
        }
        
        /** Return the other node of this Edge */
        
        //We return the second Node of this Edge.
        
        Node node2() const {
        
            return Node(graph_, secondNode);
        }
        
        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        
        // First, I'm checking whether these edges are in the same graph.
        // Second, I'm checking whether the first node of the first edge is
        // the same as the first node of the second edge and whether
        // the second node of the first edge is the same as the second node
        // of the second edge.
        // We also need to check the reverse order of the nodes since the edge
        // (1,2) and the edge (2,1) are the same.
        
        bool operator==(const Edge& e) const
        {
            if (graph_ == e.graph_)
            {
                if (((firstNode == (e.node1()).index()) && (secondNode == (e.node2()).index())) or ((firstNode == (e.node2()).index()) && (secondNode == (e.node1()).index())))
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
        
        //COMPARE EDGES HOW
        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        
        // We just compare the indices of these 2 edges and check which one
        // is larger.
        
        bool operator<(const Edge& e) const
        {
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
        
        const Graph* graph_;
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
    
    // We return the number of Edges in the Graph.
    
    size_type num_edges() const {
        
        return number_Of_Edges;
    }
    
    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    
    // We define an object pairOfNodes (discussed below) which is a pair of indices
    // of Nodes. The vectorOfEdges has pairOfNodes as elements. So, we create an
    // Edge of index i, index of Node1 and index of node2 which we get using the
    // struct pairOfNodes.
    
    Edge edge(size_type i) const
    {
        if ((0 <= i) && (i < number_Of_Edges))
        {
            pairOfNodes nodesOfThisEdge = vectorOfEdges[i];
            
            Edge e (this, i, nodesOfThisEdge.node1, nodesOfThisEdge.node2);
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
    
    //This function checks whether two nodes are connected by an edge.
    // We created a vector of vectors (adjacency list) of internal edges.(discussed below)
    // We loop over the vector of internal edges and compare the nodes of the internal
    // edges in the vector and the nodes a and b and check whether they have the same
    // indices. We also take care of the fact that the edge (a,b) and (b,a) are the
    // same.
    
    bool has_edge(const Node& a, const Node& b) const
    {
        if (has_node(a) && has_node(b))
        {
            std::vector<internalEdge> edgesConnectedToa = adjacencyList[a.index()];
            
            for (auto e : edgesConnectedToa)
            {
                if ((e.index_Of_firstNode == a.index()) && (e.index_Of_secondNode == b.index()))
                {
                    return true;
                }
                
                if ((e.index_Of_firstNode == b.index()) && (e.index_Of_secondNode == a.index()))
                {
                    return true;
                }
            }
            
            return false;
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
     ASK HERE
     * Can invalidate edge indexes -- in other words, old edge(@a i) might not
     * equal new edge(@a i). Must not invalidate outstanding Edge objects.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    
    // We first check that a and b are distinct valid nodes of this graph.
    // Then, we need to check if the edge (a,b) exists in the graph.
    // If it doesn't, we create a new Edge with nodes a and b. While doing
    // this, we should update vectorOfEdges and the adjacencyList (both
    // discussed below). Indeed, when updating the adjacency list, we need
    // to fill both edges (a,b) and (b.a). We then increment the number of
    // edges and th index of the edge.
    // If the edge already existed,we loop over the vector of internal edges
    // to dinf the right internalEdge to get its index and finally return the
    // existing edge.
    
    
    Edge add_edge(const Node& a, const Node& b)
    {
        if (has_node(a) && has_node(b) && (not (a == b)))
        {
            if (has_edge(a, b) == false)
            {
                internalEdge edgeOfab = internalEdge(next_uid_Edge, a.index(), b.index());
                adjacencyList[a.index()].push_back(edgeOfab);
                adjacencyList[b.index()].push_back(edgeOfab);
                vectorOfEdges.push_back(pairOfNodes(a.index(),b.index()));
                Edge e = Edge(this, next_uid_Edge, a.index(), b.index());
                
                number_Of_Edges = number_Of_Edges + 1;
                next_uid_Edge = next_uid_Edge + 1;
                return e;
            }
            else
            {
                std::vector<internalEdge> edgesConnectedToa = adjacencyList[a.index()];
                internalEdge edgeOfab;
                
                for (auto e : edgesConnectedToa)
                {
                    if ((e.index_Of_firstNode == a.index()) && (e.index_Of_secondNode == b.index()))
                    {
                        edgeOfab = e;
                    }
                    
                    if ((e.index_Of_firstNode == b.index()) && (e.index_Of_secondNode == a.index()))
                    {
                        edgeOfab = e;
                    }
                }
                Edge existingEdge = Edge(this, edgeOfab.index_Of_Edge, a.index(), b.index());
                
                return existingEdge;
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
    
    // We set the number of nodes and edges to 0.
    // We set the vectors used to empty vectors.
    
    void clear()
    {
        number_Of_Edges = 0;
        number_Of_Nodes = 0;
        nodes_In_The_Graph = std::vector<Point>();
        adjacencyList = std::vector<std::vector<internalEdge>>();
        vectorOfEdges = std::vector<pairOfNodes>() ;
        
    }
    
private:
    
    // This class has the object internalEdge which is defined by the
    // index of the edge, the index of the first node and the index of
    // the second node.
    // We also create an invalid internalEdge constructor and a valid
    // one which takes the index of the edge and the index of the nodes
    // as parameters.
    
    class internalEdge
    {
        friend class Graph;
        
        size_type index_Of_Edge;
        
        size_type index_Of_firstNode;
        size_type index_Of_secondNode;
        
        internalEdge()
        {
            
        }
        
        internalEdge(size_type edgeId, size_type idNode1, size_type idNode2)
        : index_Of_Edge(edgeId), index_Of_firstNode(idNode1), index_Of_secondNode(idNode2) {
        }
        
        
    };
    
    // This struct pairOfNodes has the object pairOfNodes which has the index
    // of the first node and the index of the second node.
    
    struct pairOfNodes
    {
        pairOfNodes(size_type firstNode, size_type secondNode)
        {
            node1 = firstNode;
            node2 = secondNode;
        }
        size_type node1;
        size_type node2;
        
    };
    
    // These are the private data attributes and data structures I used.
    // We have the number of nodes in the graph, number of edges in the graph, a
    // vector of Point object that represents the nodes in the graph, the index of
    // the Node, the index of the Edge, a vector of pair nodes which represents a
    // vector of edges and finally a vector of vectors of internal Edge. The
    // vectorOfEdges has a vector of internalEdge which represents the edges containing
    // the node i.
    
    unsigned int number_Of_Nodes;
    unsigned int number_Of_Edges;
    std::vector<Point> nodes_In_The_Graph;
    size_type next_uid_Node;
    size_type next_uid_Edge;
    std::vector<pairOfNodes> vectorOfEdges;
    std::vector<std::vector<internalEdge>> adjacencyList;
    
};

#endif // CME212_GRAPH_HPP
