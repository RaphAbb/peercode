#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

#include <algorithm>
#include <vector>
#include <cassert>
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


# define print_debug 0 


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

    // Declare classes
    class Node;
    class Edge;
    class NodeIterator;
    class EdgeIterator;
    class IncidentIterator;

    // Synonyms
    using node_value_type = V;
    using graph_type = Graph<V>;
    using node_type = Node;
    using edge_type = Edge;
    using node_iterator = NodeIterator;
    using edge_iterator = EdgeIterator;
    using incident_iterator = IncidentIterator;
    using size_type = unsigned;

    // Globally accessible variables
    std::vector<Node*> nodes;   // vec of edge pointers for convenience
    std::vector<Edge*> edges;   // vec of edge pointers for convenience
    std::vector<Point> coords;  // vec of positions of each node
    std::vector<V> values;      // vec of values at each node
    std::vector<std::vector<size_type>>ebuddies;//ebuddies[node]: edges attached w/ node1=node
    std::vector<std::vector<size_type>>nbuddies;//nbuddies[node]: nodes attached to node 
    std::vector<std::vector<size_type>> ecaps;  //ecaps[edge]:    nodes capping the edge

    // Constructor init with zero nodes and edges 
    Graph() {}


    // Default destructor
    ~Graph() = default;
  

    /** @class Graph::Node
     * @brief Class representing the graph's nodes.
     * Node objects are used to access information about the Graph's nodes.
     */
    class Node : private totally_ordered<Node> {
        private:
            // Allow Graph to access Node's private member data and functions.
            friend class Graph;
            size_type idx_;
            Graph* graph_;
        public:
            Node(const Graph* graph_, size_type idx_)
                    : idx_(idx_), graph_(const_cast<Graph*>(graph_)) {}
            Node() {}

            // Return this node's position.
            const Point& position() const {
                return this->graph_->coords[this->idx_];
            }

            // Return this node's index, a number in the range [0, graph_size)
            size_type index() const {
                return this->idx_;
            }

            // Test whether this node and n are equal.
            bool operator==(const Node& n) const {
                if ((idx_==n.index()) and (graph_==n.graph_)){return true;}
                else{return false;}
            }

            // Test whether this node is less than @a n in a global order.
            // The node ordering relation must obey trichotomy: For any two nodes x
            // and y, exactly one of x == y, x < y, and y < x is true.
            bool operator<(const Node& n) const {
                return idx_< n.index(); 
            }
            
            node_value_type& value(){
                return this->graph_->values[this->idx_];
            };

            const node_value_type& value() const {
                return this->graph_->values[this->idx_];
            };

            size_type degree() const{
                return (graph_->nbuddies[this->index()]).size();
            };
            
            incident_iterator edge_begin() const{
                return IncidentIterator(graph_,this,0);
            }
            
            incident_iterator edge_end() const{
                int N = ebuddies[this->index()].size()-1;
                return IncidentIterator(graph_,this,N);
            }
    };


    /** Add a node to the graph, returning the added node.
     * @param[in] position The new node's position
     * @post new num_nodes() == old num_nodes() + 1
     * @post result_node.index() == old num_nodes()
     * Complexity: O(1) amortized operations.
     */
    Node add_node(const Point& position, const node_value_type& val = node_value_type()){
        // Make a new node, but don't push into the global vector yet.
        Node* n = new Node(this,nodes.size());
            
        // put into global vectors
        this->nodes.push_back(n);
        this->coords.push_back(position);
        this->values.push_back(val);
        //(*n).value() = val;
        
        // make space to list node and edge buddies to this new node
        std::vector<unsigned>temp_vec {};
        nbuddies.push_back(temp_vec);
        ebuddies.push_back(temp_vec);
        
        if (print_debug){
            std::cout << "adding node " << (*n).index();
            std::cout << " at " << (*n).position();
            std::cout << " address: " << n;
            std::cout<<std::endl;
        }
        return *n; 
    }


    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
        if (this==n.graph_){return true;}
        else{return false;}
    }


    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     * Complexity: O(1).
     */
    Node node(size_type i) const {return Node(this,i); }


    /** @class Graph::Edge
     * @brief Class representing the graph's edges.
     * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
     * are considered equal if they connect the same nodes, in either order.
     */
    class Edge : private totally_ordered<Edge> {
        private:
            // Allow Graph to access Edge's private member data and functions.
            friend class Graph;
            Graph* graph_;
            size_type idx_;
        public:
            // Constructor
            Edge() {}
            Edge(const Graph* graph_, size_type idx_) : 
                 graph_(const_cast<Graph*>(graph_)), idx_(idx_){}

            /** Return a node of this Edge */
            Node node1() const {
                int i = graph_->ecaps[idx_][0];   
                return graph_->node(i);
            }

            /** Return the other node of this Edge */
            Node node2() const {
                int i = graph_->ecaps[idx_][1];   
                return graph_->node(i);
            }

            // Test whether this edge and @a e are equal.
            //Equal edges represent the same undirected edge between two nodes.
            bool operator==(const Edge& e) const {
                // test n1=n1 and n2=n2 (same direction)
                if (node1()==e.node1()){
                    if (node2()==e.node2()){return true;}
                }
                // test n2=n1 and n1=n2 (opp direction)
                else if (node1()==e.node2()){
                    if (node2()==e.node1()){return true;}
                }
                return false;
            }

            // Test whether this edge is less than @a e in a global order.
            bool operator<(const Edge& e) const {
                return (idx_<e.idx_);
            }
    };
   

    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     * Complexity: O(1) 
     */
    Edge edge(size_type i) const {
        return Edge(this,i); // Invalid Edge
    }


    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     * Complexity: approx O(1)     
     */
    bool has_edge(const Node& a, const Node& b) const {
        // check if b is a neighbor to a
        for (auto idx : nbuddies[a.index()]){
            if (idx==b.index()){return true;}
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
        // Check that a and b are valid nodes of this graph
        if ( !( has_node(a) and has_node(b) ) ){
            if (print_debug){
                std::cout << "Node(s) at "<< a.position() ;
                std::cout <<" or "<< b.position() << " invalid " ;
                std::cout << std::endl;
            }
            return Edge();
        }
       
        // Check that a and b are distinct
        if ( a==b ){
            if (print_debug){
                std::cout << "Node(s) at "<< a.position() ;
                std::cout <<" or "<< b.position() << " duplicate " ;
                std::cout << std::endl;
            }
            return Edge();
        }

        // Return a new edge object with corresponding nodes
        int idx = edges.size();
        Edge* e = new Edge(this,idx);

        // Check if edge already exists
        if (has_edge(a,b) or has_edge(b,a)){
            if (print_debug){
                std::cout << "Edge between "<< a.position() ;
                std::cout <<" and "<< b.position() << " already exists" ;
                std::cout << std::endl;
            }
            // Check if b is one of the nodes attached to a
            for ( int idx : ebuddies[a.index()] ){
                if ( (edge(idx)).node2() == b ){
                    return edge(idx);
                }
            }
            // Check if a is one of the nodes attached to b
            for ( int idx : ebuddies[b.index()] ){
                if ( (edge(idx)).node2() == a ){
                    return edge(idx);
                }
            }
        }

        // Put into global vec
        edges.push_back(e);

        // Keep track of the caps on this edge
        std::vector<size_type> caps {a.index(),b.index()};
        ecaps.push_back(caps);
        
        // Record this edge as a buddy to node1
        ebuddies[((*e).node1()).index()].push_back(idx);

        // Record b as a node buddy to a and vice versa
        nbuddies[a.index()].push_back(b.index());
        nbuddies[b.index()].push_back(a.index());

        if (print_debug){
            std::cout << "adding edge between nodes at ";
            std::cout << a.position() << " & " << b.position();
            std::cout << std::endl;
        }
        return *e;
    }


    /** @class Graph::NodeIterator
     * @brief Iterator class for nodes. A forward iterator. Just an index we 
     * move forward when requested to access the graph's node or coordinate vec.
     */
    class NodeIterator {
        public:
            // These type definitions let us use STL's iterator_traits.
            using value_type      = Node;                     // Element type
            using pointer         = Node*;                    // Pointers to elements
            using reference       = Node&;                    // Reference to elements
            using difference_type = std::ptrdiff_t;           // Signed difference
            using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
            
            /** @brief Construct an invalid NodeIterator.
             */
            NodeIterator(const Graph* graph, size_type itr ):
                graph_(const_cast<Graph*>(graph)), itr_(itr) {} 
            NodeIterator(){}
           
            /** @brief dereferencer
             * @pre     num_nodes > 0
             */
            Node operator*() const{
                return graph_->node(itr_);
            };

            /** @brief Iterates the iterator by one step 
             * @pre !(itr_==end) 
             */
            NodeIterator& operator++(){
                ++itr_;
                return *this;
            };
            
            /** @brief Compares two predicates based off current and end iterators
            */
            bool operator!=(const NodeIterator& n) const{
                return !(itr_==n.itr_) or !(graph_==n.graph_);

            }

            /** @brief Compares two predicates based off current and end iterators
            */
            bool operator==(const NodeIterator& n) const{
                return (itr_==n.itr_) and (graph_==n.graph_);
            }
        private:
            friend class Graph;
            Graph* graph_;
            size_type itr_;
    };


    /** Returns an iterator beginning at the first node of the graph. 
     * @pre num_nodes>0 
     * Complexity: O(1).
     * */
    node_iterator node_begin() const{
        return NodeIterator(this,0);
    };


    /** Returns an iterator with position at one past the end node of the graph. 
     * @pre num_nodes>0 
     * Complexity: O(1).
     * */
    node_iterator node_end() const{
        return NodeIterator(this,num_nodes());
    };
   

    /** @class Graph::IncidentIterator
     * @brief Iterator class for edges incident to a node. A forward iterator. */
    class IncidentIterator {
        public: 
            // These type definitions let us use STL's iterator_traits.
            using value_type        = Edge;                     // Element type
            using pointer           = Edge*;                    // Pointers to elements
            using reference         = Edge&;                    // Reference to elements
            using difference_type   = std::ptrdiff_t;           // Signed difference
            using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
            
            /** @brief Construct a valid or an invalid IncidentIterator. 
             */
            IncidentIterator(const Graph* graph, const Node* node, size_type itr )
                : graph_(const_cast<Graph*>(graph)), 
                  node_(const_cast<Node*>(node)), itr_(itr) {} 
            IncidentIterator(){}

            /** @brief dereferencer first fetches the index of the edge
             * connected to this node listed in ebuddies, then returns that edge.
             * @pre     num_nodes > 0
             * @pre     num_edges > 0
             */
            Edge operator*() const{
                size_type idx = graph_->ebuddies[(*node_).index()][itr_];
                return graph_->edge(idx);    
            }

            /** @brief Iterates on which position in 
             * the row ebuddies[node index] we point to
             */
            IncidentIterator& operator++(){
                ++itr_;
                return *this;
            }

            /** @brief Compares two predicates based off current and end iterators
            */
            bool operator==(const IncidentIterator& i) const{
                return (itr_==i.itr_) and (graph_==i.graph_) and (node_==i.node_);
            }
        private:
            friend class Graph;
            friend class Node;
            Graph* graph_;
            Node* node_;
            size_type itr_;
    };
   

    /** @class Graph::EdgeIterator
     * @brief Iterator class for edges. A forward iterator. 
     */
    class EdgeIterator {
        public:
            // These type definitions let us use STL's iterator_traits.
            using value_type      = Edge;                     // Element type
            using pointer         = Edge*;                    // Pointers to elements
            using reference       = Edge&;                    // Reference to elements
            using difference_type = std::ptrdiff_t;           // Signed difference
            using iterator_category = std::input_iterator_tag;// Weak Category, Proxy
            
            /** @brief Construct a valid or an invalid IncidentIterator. 
             */
            EdgeIterator(){}
            EdgeIterator(const Graph* graph, size_type itr )
                : graph_(const_cast<Graph*>(graph)), itr_(itr)  {} 

            /** @brief dereferencer
             * @pre     num_edges > 0
             */
            Edge operator*() const{
                return graph_->edge(itr_);
            };

            /** @brief Iterates the iterator by one step 
             * @pre !(itr_==end) 
             */
            EdgeIterator& operator++(){
                ++itr_;
                return *this;
            };
           
            /** @brief Compares two predicates based off private vars 
            */
            bool operator!=(const EdgeIterator& e) const{
                return (graph_!=e.graph_) or (itr_!=e.itr_);
            }

            /** @brief Compares two predicates based off private vars 
            */
            bool operator==(const EdgeIterator& e) const{
                return (graph_==e.graph_) and (itr_==e.itr_);
            }
        private:
            friend class Graph;
            Graph* graph_;
            size_type itr_;
    };

    
    /** Returns an iterator beginning at the first edge of the graph. 
     * @pre num_nodes>0 
     * Complexity: O(1).
     */
    edge_iterator edge_begin() const{
        return EdgeIterator(this,0);
    };
    
    
    /** Returns an iterator at one past the last  edge of the graph. 
     * @pre num_nodes>0 
     * Complexity: O(1).
     */
    edge_iterator edge_end() const{
        return EdgeIterator(this,num_edges());
    };
   

    /** @brief Remove all items from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        for (Edge* e : edges){ delete e;}
        for (Node* n : nodes){ delete n;}
        edges.clear();
        nodes.clear();
        values.clear();
        coords.clear();
        ebuddies.clear();
        nbuddies.clear();
        ecaps.clear();
    }


    /** Return the number of nodes in the graph.
     * Complexity: O(1).
     */
    size_type size() const {return nodes.size();}

    
    /** Return the number of nodes in the graph.
     * Complexity: O(1).
     */
    size_type num_nodes() const {return size();}


    /** Return the number of eddes in the graph.
     * Complexity: O(1).
     */
    size_type num_edges() const {return edges.size();}
};//Graph

#endif // CME212_GRAPH_HPP
