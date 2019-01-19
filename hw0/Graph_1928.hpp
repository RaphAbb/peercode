#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

#include <algorithm>
#include <vector>
#include <cassert>
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

# define print_debug 0

class Graph {

private:
public:

    class Node;
    class Edge;
    using graph_type = Graph;
    using node_type = Node;
    using edge_type = Edge;
    using size_type = unsigned;
    std::vector<Node*> nodes;
    std::vector<Edge*> edges;
    std::vector<Point> coords;
    std::vector<std::vector<unsigned>> ebuddies;
    std::vector<std::vector<unsigned>> nbuddies;

    // constructor init with zero nodes and edges 
    Graph() {}


    // Default destructor
    ~Graph() = default;
   

    class Node {
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
    };


    class Edge {
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
                int i = graph_->ebuddies[idx_][0];   
                return graph_->node(i);
            }

            /** Return the other node of this Edge */
            Node node2() const {
                int i = graph_->ebuddies[idx_][1];   
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
    


    //Add a node to the graph, returning the added node. 
    Node add_node(const Point& position) {
        // Make a new node, but don't push into the global vector yet.
        Node* n = new Node(this,nodes.size());

        // Check if node already exists:
        if (false){//(has_node(*n)){
            return Node();
        }else{
            // put into global vectors
            nodes.push_back(n);
            coords.push_back(position);
            
            // make space to list node buddies
            std::vector<unsigned>temp_vec {};
            nbuddies.push_back(temp_vec);
           
            if (print_debug){
                std::cout << "adding node " << (*n).index();
                std::cout << " at " << (*n).position();
                std::cout << " address: " << n;
                std::cout<<std::endl;
            }
            return *n; 
        }
    }


    // Determine if a Node belongs to this Graph
    bool has_node(const Node& n) const {
        if (this==n.graph_){return true;}
        else{return false;}
    }


    // Return the node with index i<nnodes
    Node node(size_type i) const {return Node(this,i); }


    // Return the edge with index i<num_edges
    Edge edge(size_type i) const {
        return Edge(this,i); // Invalid Edge
    }


    // Test whether two nodes are connected by an edge.
    bool has_edge(const Node& a, const Node& b) const {
        // check if b is a neighbor to a
        for (unsigned idx : nbuddies[a.index()]){
            if (idx==b.index()){return true;}
        }
        return false;
    }


    // Add an edge to the graph, or return the current edge if it already exists.
    // Can invalidate edge indexes -- in other words, old edge(i) might not
    // equal new edge(i). Must not invalidate outstanding Edge objects.
    // Complexity: No more than O(num_nodes() + num_edges()), hopefully less
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
            return Edge();
        }
        else{
            // Put into global vec
            edges.push_back(e);
            std::vector<unsigned>temp_vec {a.index(),b.index()};
            ebuddies.push_back(temp_vec);

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
    }


    // Remove all nodes and edges from this graph.
    //Invalidates all outstanding Node and Edge objects.
    void clear() {
        for (Edge* e : edges){ delete e;}
        for (Node* n : nodes){ delete n;}
        edges.clear();
        nodes.clear();
        ebuddies.clear();
        nbuddies.clear();
    }


    // Return the number of nodes in the graph. Complexity: O(1).
    size_type size() const {return nodes.size();}

    
    // Size of graph = number of nodes
    size_type num_nodes() const {return size();}


    // Return the total number of edges in the graph.
    size_type num_edges() const {return edges.size();}


};//Graph

#endif // CME212_GRAPH_HPP
