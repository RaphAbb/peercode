#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

#include <stdio.h>     
#include <assert.h>      
#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

#define debug_g 0 // general graph debug 
#define debug_d 0 // debug node/edge deletion

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
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
    using edge_value_type = E;
    using graph_type = Graph<V,E>;
    using node_type = Node;
    using edge_type = Edge;
    using node_iterator = NodeIterator;
    using edge_iterator = EdgeIterator;
    using incident_iterator = IncidentIterator;
    using size_type = unsigned;

    // Globally accessible variables
    std::vector<Point> coords;  // vec of positions of each node
    std::vector<V> nvalues;     // vec of values at each node
    std::vector<E> evalues;     // vec of values at each edge
    std::vector<unsigned> empty_vec {}; // Just use this to allocate space later
    std::vector<std::vector<size_type>> ebuddies;
    std::vector<std::vector<size_type>> nbuddies;
    std::vector<std::vector<size_type>> ecaps;  
    //ebuddies[node]: edges attached w/ node1=node
    //nbuddies[node]: nodes attached to node     
    //ecaps[edge]:    nodes capping the edge     
    std::unordered_map<unsigned,unsigned> nID;     // [originalID,currentID] for nodes 
    std::unordered_map<unsigned,unsigned> eID;     // [originalID,currentID] for nodes
    std::unordered_map<unsigned,unsigned> nID_inv; // [currentID,originalID] for edges
    std::unordered_map<unsigned,unsigned> eID_inv; // [currentID,originalID] for edges


    // Constructor
    Graph() {
        assert (num_nodes()==0);
        assert (num_edges()==0);
    }

    // Default destructor
    ~Graph() = default;
//-----------------------------------------------------------------------
    /** @class Graph::Node
     * @brief Class representing the graph's nodes.
     * Node objects are used to access information about the Graph's nodes.
     */
    class Node : private totally_ordered<Node> {
        public:
            Node(const Graph* graph, size_type uid)
                : uid_(uid), graph_(const_cast<Graph*>(graph)) {}
            Node() {}

            /** @brief Return this node's position. 
             * Searches in the coords/position library for the
             * position in the uid row. Can be lvalue as well.
             * Complexity: O(1)
             */
            const Point& position() const{return graph_->coords.at(uid_);}
            Point& position(){return graph_->coords.at(uid_);}

            /** @brief Return this node's @param id, a number in the range [0, graph_size)
             * finds @param uid_ in the @param nIDmap, and returns the id pair
             * @pre num_nodes()>0
             * Complexity: O(1)
             */
            size_type index() const {
                unsigned id = graph_->nID.find(uid_)->second;
                assert (id <= graph_->num_nodes());
                return id;
            }

            /** @brief Does the node really exist?
             *  Get the id:uid pair in both maps, basically 
             *  checking if the uid is still in the node map.
             * Complexity: O(1)
             */
            bool valid() const {
                if ( (graph_->nID.find(uid_)==graph_->nID.end()) )
                    return false;
                return true; 
            }

            /** @brief Test whether this node and n are equal.
             * Tests the @param uid_ and @param graph_
             * Complexity: O(1)
             */
            bool operator==(const Node& n) const {
                if ((uid_==n.uid_) and (graph_==n.graph_))
                    return true;
                return false;
            }

            /** @brief Test whether this node is less than @a n in a global order.
             * The node ordering relation must obey trichotomy: For any two nodes x
             * and y, exactly one of x == y, x < y, and y < x is true.
             * Complexity: O(1)
             */
            bool operator<(const Node& n)const{return uid_ < n.uid_;}
            
            /** @brief Return this node's value 
             * Searches in @param values vector for the
             * position in the uid row. Can be lvalue as well.
             * Complexity: O(1)
             */
            node_value_type& value(){return this->graph_->nvalues.at(uid_);}
            const node_value_type& value()const{return this->graph_->nvalues.at(uid_);}
            
            /** @brief @param nbuddies contains @param uid's,
             * so degree returns the size of the @param nbuddies row vector
             * Note: this is not technically correct once nodes are removed,
             * but we like to check the validity of each neighboring node outside 
             * this function.
             * @pre num_edges > 0
             * @pre num_nodes > 0
             * Complexity: O(1)
             */
            size_type degree()const{return (graph_->nbuddies.at(uid_)).size();}

            /** @brief Return this node's edge iterators 
             * to select all the edges, valid or not, 
             * which are attached to this node. We check validity elsewhere.
             * Complexity: O(1)
             */
            incident_iterator edge_begin() const{
                return IncidentIterator(graph_,this,0);
            }
            incident_iterator edge_end() const{
                return IncidentIterator(graph_,this,degree());
            }
        private:
            // Allow Graph to access Node's private member data and functions.
            friend class Graph;
            size_type uid_;
            Graph* graph_;
    };


    /** Add a node to the graph, returning the added node.
     * @param[in] position The new node's position
     * @post new num_nodes() == old num_nodes() + 1
     * @post result_node.index() == old num_nodes()
     * Complexity: O(1) 
     */
    Node add_node(const Point& position, const node_value_type& val = node_value_type()){
        // Make a new node, but don't push into the global vector yet.
        unsigned uid = num_nodes();
        Node n(this,uid);

        // put indices in
        nID.insert(std::make_pair(uid,uid));
        nID_inv.insert(std::make_pair(uid,uid));

        // Make space to list node and edge buddies to this new node
        // If there is already space allocated, then empty the rows
        assert (nbuddies.size()==ebuddies.size());
        assert (coords.size()==nvalues.size());
        if (uid < nbuddies.size()){ 
            nbuddies.at(uid) = empty_vec;
            ebuddies.at(uid) = empty_vec;
            coords.at(uid) = position;
            nvalues.at(uid) = val;
        }else{ // allocate some space, we haven't added this many nodes/edges yet
            nbuddies.push_back(empty_vec);
            ebuddies.push_back(empty_vec);
            coords.push_back(position);
            nvalues.push_back(val);
        }
        assert (n.position()==position);
        //assert (n.value()==val);

        if (debug_g){
            std::cout << "adding node(" << uid << ")";
            std::cout << " at " << Point(coords.at(uid));
            //std::cout << "\twith value " << nvalues[uid];
            std::cout<<std::endl;
        }
        return n;
    }


    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     * Complexity: O(1)
     */
    bool has_node(const Node& n)const{return this==n.graph_;}


    /** Return the node with current index @a i.
     * Essentially reconstructs the node with corresponding uid.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     * Complexity: O(1).
     */
    Node node(size_type id) const {
        assert (id <= num_nodes());
        auto uid = nID_inv.find(id)->second;
        return Node(this,uid); 
    }
//-----------------------------------------------------------------------
    /** @class Graph::Edge
     * @brief Class representing the graph's edges.
     * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
     * are considered equal if they connect the same nodes, in either order.
     */
    class Edge : private totally_ordered<Edge> {
        public:
            // Constructor for valid or invalid edges
            Edge(const Graph* graph, size_type uid)  
                 : graph_(const_cast<Graph*>(graph)), uid_(uid){}
            Edge() {}

            /** @brief Return the nodes of this Edge 
             * by accessing the graph's @param ecaps data structure
             * Complexity: O(1)
             */
            Node node1() const{
                unsigned uid = graph_->ecaps.at(uid_)[0];
                return Node(graph_,uid);
            }
            Node node2() const{
                unsigned uid = graph_->ecaps.at(uid_)[1];
                return Node(graph_,uid);
            }

            /** @brief Does this edge really exist in the graph? 
             * we check if both nodes attached to this edge are valid
             * and if the edge's @a uid_ is still in the edge map @param eID
             * Complexity: O(1)
             */
            bool valid() const{
                // both end nodes must be valid
                if (!node1().valid() or !node2().valid()){return false;}

                // Check if the uid is still in the edge map
                if (graph_->eID.find(uid_)==graph_->eID.end()){return false;}

                return true;
            }

            /** @brief Returns the Euclidean distance betw node1 and node2
             * Complexity: O(1)
             */
            double length() const{
                Point coord1 = node1().position();
                Point coord2 = node2().position();
                return norm(coord1-coord2); 
            }
            
            /** @brief Return or set the value for this edge 
             * by accessing the graph's @param evalues data structure.
             * Complexity: O(1)
             */
            edge_value_type& value(){return graph_->evalues.at(uid_);}
            const edge_value_type& value() const{return graph_->evalues.at(uid_);}

            /** @brief Test whether this edge and @a e are equal.
            * Equal edges represent the same undirected edge between two nodesi.
            * Complexity: O(1)
            */
            bool operator==(const Edge& e) const{
                // Graphs must be the same
                if (graph_!=e.graph_){return false;}
                
                // test n1=n1 and n2=n2 (same direction)
                if (node1()==e.node1()){
                    if (node2()==e.node2()){
                        return true;
                    }
                }
                // test n2=n1 and n1=n2 (opp direction)
                else if (node1()==e.node2()){
                    if (node2()==e.node1()){
                        return true;
                    }
                }
                return false;
            }
            
            /** @brief Test whether this edge is less than @a e in a global order.
            * Complexity: O(1)
            */
            bool operator<(const Edge& e) const {
                return graph_<e.graph_;
            }
        private:
            // Allow Graph to access Edge's private member data and functions.
            friend class Graph;
            Graph* graph_;
            size_type uid_;
    };
   

    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     * Complexity: O(1) 
     */
    Edge edge(size_type id) const {
        unsigned uid = eID_inv.find(id)->second;
        return Edge(this,uid); 
    }
   

    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     * Complexity: approx O(1)     
     */
    bool has_edge(const Node& a, const Node& b) const {
        //check if b is a neighbor to a (compares the uids)
        for (auto uid : ebuddies.at(a.uid_)){
            auto e = Edge(this,uid);
            if (e.node2()==b)
                return e.valid();
        }
        for (auto uid : ebuddies.at(b.uid_)){
            auto e = Edge(this,uid);
            if (e.node2()==a)
                return e.valid();
        }
        //unsigned uid = get_edge_uid(a,b);
        return false;//Edge(this,uid).valid();
    }
    
    
    /** @brief Fetches the ORIGINAL edge index @a uid between nodes a and b by checking
    *  each of the rows in ebuddies[a] and ebuddies[b]. Can be given either
    *  the nodes or the nodes' @a uid's
    */
    unsigned get_edge_uid(const Node& a, const Node& b) const{
        // Check if b is one of the nodes attached to a
        for ( unsigned uid : ebuddies.at(a.uid_) ){
            if ( (Edge(this,uid)).node2() == b )
                return uid;
        }
        // Check if a is one of the nodes attached to b
        for ( unsigned uid : ebuddies.at(b.uid_) ){
            if ( (Edge(this,uid)).node2() == a )
                return uid;
        }
        std::cout << "Error getting edge index\n";
        return -1;
    }
    unsigned get_edge_uid(const unsigned a_uid, const unsigned b_uid) const{
        // Check if b is one of the nodes attached to a
        for ( unsigned uid : ebuddies.at(a_uid) ){
            if ( (Edge(this,uid)).node2() == Node(this,b_uid) )
                return uid;
        }
        // Check if a is one of the nodes attached to b
        for ( unsigned uid : ebuddies.at(b_uid) ){
            if ( (Edge(this,uid)).node2() == Node(this,a_uid) )
                return uid;
        }
        std::cout << "Error getting edge index\n";
        return -1;
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
    Edge add_edge(const Node& a, const Node& b,
        const edge_value_type& val = edge_value_type()) {
        assert (a.valid());
        assert (b.valid());

        // Check that a and b are valid nodes of this graph
        if ( !( has_node(a) and has_node(b) ) ){
            if (debug_g){
                std::cout << "Node(s) at "<< a.position() ;
                std::cout <<" or "<< b.position() << " invalid " ;
                std::cout << std::endl;
            }
            return Edge();
        }
       
        // Check that a and b are distinct
        if ( a==b ){
            if (debug_g){
                std::cout << "Node(s) at "<< a.position() ;
                std::cout <<" or "<< b.position() << " duplicate " ;
                std::cout << std::endl;
            }
            return Edge();
        }

        // Check if edge already exists
        if (has_edge(a,b)){
            if (debug_g){
                std::cout << "Edge between "<< a.position() ;
                std::cout <<" and "<< b.position() << " already exists" ;
                std::cout << std::endl;
            }
            unsigned uid = get_edge_uid(a,b);
            return Edge(this,uid);
        }

        // What is the uid of this edge?
        unsigned uid = num_edges();

        // Make a vector of node caps for this edge (node1=a,node2=b)
        std::vector<unsigned> caps {a.uid_, b.uid_}; 

        // Record the original index as the current one
        eID.insert(std::make_pair(uid,uid));
        eID_inv.insert(std::make_pair(uid,uid));

        // Record node caps for this edge
        // Record values for this edge
        // If there is already space allocated, then overwrite 
        assert (evalues.size()==ecaps.size());
        if (uid < evalues.size()){ 
            ecaps.at(uid)   = caps;
            evalues.at(uid) = val;
        }else{ // if there isn't space allocated, then pushback
            ecaps.push_back(caps);
            evalues.push_back(val);
        }
         
        ebuddies.at(a.uid_).push_back(uid);    // Record this edge as a buddy to node1
        nbuddies.at(a.uid_).push_back(b.uid_); // Record b as a node buddy to a
        nbuddies.at(b.uid_).push_back(a.uid_); // Record a as a node buddy to b
        assert (a.degree()>0);
        assert (b.degree()>0);
        if (debug_d){
            std::cout << "Added edge uid = "<< uid << std::endl;;
            std::cout << "\tNode a degree: "<<a.degree() << std::endl;
            std::cout << "\tNode b degree: "<<b.degree() << std::endl;
        }
        auto e = Edge(this,uid);
        assert (e.node1()==a);
        assert (e.node2()==b);
        //assert (e.value()==val);

        if (debug_g){
            std::cout << "adding edge(" << uid << ")";
            std::cout << " between node("<< a.uid_<< ")";
            std::cout << " at " << a.position();
            std::cout << " & node("<< b.uid_<< ")";
            std::cout << " at " << b.position();
            std::cout << std::endl;
        }
        return e; 
    }
//-----------------------------------------------------------------------
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
                assert (graph_->node(itr_).valid());
                return graph_->node(itr_);
            };

            /** @brief Iterates the iterator by one step 
             * @pre     !(itr_==end) 
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
//-----------------------------------------------------------------------
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
           
            // for each node attached: 
            // check ebuddies row[node1] and row[node2] checks (add_edge) 
            // gets edge index
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
                // get the paired node
                unsigned uid_  = (*node_).uid_; 
                unsigned n_uid = graph_->nbuddies.at(uid_).at(itr_);
                unsigned e_uid = graph_->get_edge_uid( uid_, n_uid );
                return Edge(graph_,e_uid);    
            }

            /** @brief Iterates on which position in 
             * the row nbuddies[node index] we point to
             * i.e. which node buddy are we pointing to?
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
            bool operator!=(const IncidentIterator& i) const{
                return (itr_!=i.itr_) or (graph_!=i.graph_) or (node_!=i.node_);
            }
        private:
            friend class Graph;
            friend class Node;
            Graph* graph_;
            Node* node_;
            size_type itr_;
    };
//-----------------------------------------------------------------------
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
//-----------------------------------------------------------------------
    /** @brief Removes a node and the edges attached to it. 
     * Returns false if the node has already been removed.
     * Returns true is the node has sucessfully been removed.
     * Uses the incident iterator to obtain the @a uid's for this node, then calls
     * remove_edge().
     * Note: All the removal methods call remove_node(n).
     *
     * @pre num_nodes > 0
     * @post updated Nnodes = Nnodes-1
     * @post updated Nedges < Nedges
     * Complexity: worst case, O(Nnodes+Nedges)
     */ 
    bool remove_node(const Node& n){
        // Make sure this node is valid
        if (!n.valid()){return false;} 

        // Let's save this id:uid pair for later use
        unsigned id = n.index();
        unsigned uid = n.uid_;
        assert (n.degree()==nbuddies.at(uid).size());
        if (debug_d){
            std::cout << "\nErasing node with id:uid = ";
            std::cout << id << ":" << uid << "\n";
            std::cout << "\tDegree = "<< n.degree() << " max edges to remove\n";
        }

        // Remove the attached edges (only the valid ones!)
        bool success=true;
        if (debug_d){std::cout << "\tRemoving edges...\n";}
        for (auto itr=n.edge_begin(); itr!=n.edge_end(); ++itr){
            auto e = *itr;
            if ( e.valid() ){
                if (debug_d){
                    std::cout << "\tRemoving attached edge #" << itr.itr_;
                    std::cout << ", uid = " << e.uid_<<"\n";
                }
                success = success and remove_edge(e);
            }else{ 
                if (debug_d){
                    std::cout << "\tAttached edge #" << itr.itr_;
                    std::cout << ", uid = " << e.uid_ <<" already removed\n";
                }
            }
        }
        if (!success){return false;}

        // remove the key-value pairs from our nodeID maps
        nID.erase(uid);
        nID_inv.erase(id); 
        
        // check if we removed the maps' last node
        if (id==num_nodes()){return true;}

        // Otherwise, update the remaining higher id's
        if (debug_d){
            std::cout << "\tUpdating node indices ";
            std::cout << id+1 <<" to "<<num_nodes() << "...\n";
        }
        for (auto j=id+1; j<num_nodes()+1; j++){
            auto pair = nID_inv.find(j);
            id  = pair->first;  // this will be decremented
            uid = pair->second; // this remains the same
            
            // decrement the id in the nodeID maps:
            // remove [id:uid] and insert [id-1:uid] since keys are immutable
            nID_inv.erase(id);
            nID_inv.insert(std::make_pair(id-1,uid));
            nID.erase(uid);
            nID.insert(std::make_pair(uid,id-1));
        }

        if (debug_d)
            std::cout << "\tDone updating node ids, Nnodes = "<< num_nodes() <<"\n";
        assert (!n.valid());
        return true; 
    }
    node_iterator remove_node(node_iterator n_it){
        // if correctly removed, return iterator to the updated node (same id)
        unsigned id = (*n_it).index();
        if (remove_node(*n_it)){
            return NodeIterator(this,id);
        }
        // else, return the end iterator
        return node_end();
    } 
    
    
    /** Removes an edge by removing key-values pairs from edgeID
     * maps. Again, we update the indices here that are higher than the @id removed.
     * Returns false if the edge is already removed.
     * Returns true if the edge is successfully removed.
     * Note: All the removal methods call remove_edge(e).
     *
     * @pre Nnodes >= 2
     * @post updated Nedges = Nedges-1
     * Complexity: O(Nedges), worst case;
     */ 
    bool remove_edge(const Node& n1, const Node& n2){
        if (!has_edge(n1,n2)){return false;}
        unsigned uid = get_edge_uid(n1,n2);
        return remove_edge(Edge(this,uid));
    } 
    bool remove_edge(const Edge& e){
        if (!e.valid()){return false;}
        
        // remove these key-value pairs from our edgeID maps
        unsigned uid = e.uid_;
        unsigned id = eID.find(uid)->second;
        eID.erase(uid);
        eID_inv.erase(id);
        if (debug_d){
            std::cout << "\t\tErasing edge with id:uid = ";
            std::cout << id << ":" << uid;
            std::cout << ", new Nedges = "<< num_edges() <<"\n";
        }

        // check if we removed the maps' last node
        if (id>=num_edges()){
            if (debug_d) {std::cout << "\tNo indicies to update\n";}    
            return true;
        }
        if (debug_d){
            std::cout << "\t\tUpdating edge indices ";
            std::cout << id+1 <<" to "<<num_edges() << "\n";
            //std::cout << "\tuid\tCurrent id\tUpdated id\n";
        }
        
        // update the remaining indices:
        for (auto j=id+1; j!=num_edges()+1; j++){
            auto pair = eID_inv.find(j);
            id  = pair->first;  // this will be decremented
            uid = pair->second; // this remains the same
            if (0)
                std::cout << "\t\t" << uid << "\t" << id << "\t\t";
            
            // decrement the id in the nodeID maps:
            // remove [id:uid] and insert [id-1:uid] since keys are immutable
            eID_inv.erase(id);
            eID_inv.insert(std::make_pair(id-1,uid));
            eID.erase(uid);
            eID.insert(std::make_pair(uid,id-1));
            if (0)
                std::cout << eID.find(uid)->second << "\n";
        }
        if (debug_d){std::cout << "\t\tDone updating edge indices\n";}
        return true; 
    }
    edge_iterator remove_edge(edge_iterator e_it){
        // if correctly removed, return iterator to the updated edge (same id)
        unsigned uid = (*e_it).uid_;
        unsigned id  = eID.find(uid)->second;
        if (remove_edge(*e_it)){
            return EdgeIterator(this,id);
        }
        // else, return the end iterator
        return edge_end();
    }
//-----------------------------------------------------------------------
    /** @brief Remove all items from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        nvalues.clear();
        evalues.clear();
        coords.clear();
        ecaps.clear();
        ebuddies.clear();
        nbuddies.clear();
        nID.clear();
        eID.clear();
        nID_inv.clear();
        eID_inv.clear();
    }


    /** Return the number of nodes in the graph.
     * Complexity: O(1).
     */
    size_type size() const {return nID.size();}

    
    /** Return the number of nodes in the graph.
     * Complexity: O(1).
     */
    size_type num_nodes() const {return size();}


    /** Return the number of eddes in the graph.
     * Complexity: O(1).
     */
    size_type num_edges() const {return eID.size();}
};//Graph

#endif // CME212_GRAPH_HPP
