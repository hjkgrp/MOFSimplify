
/* This file contains all of the data structures used when converting a
 * VORONOI_NETWORK into a DIJKSTRA_NETWORK.
 */

//#include "network.h"
#include <cstdlib>
#include <cmath>

#include "graphstorage.h"

using namespace std;

/* Begin of DELTA_POS functionality*/
/** Simple class used to store and manipulate set of three integers. */
DELTA_POS::DELTA_POS(int myX, int myY, int myZ) : x(myX), y(myY), z(myZ){}
const bool DELTA_POS::equals(DELTA_POS other) const
{
    return (other.x == x) && (other.y == y) && (other.z == z);
}
const DELTA_POS DELTA_POS::operator*(const int& factor) const
{
	return DELTA_POS(x*factor, y*factor, z*factor);
}
const DELTA_POS DELTA_POS::operator*(const DELTA_POS& other) const
{
	return DELTA_POS(x*other.x, y*other.y, z*other.z);
}
const DELTA_POS DELTA_POS::operator+(const DELTA_POS& other) const
{
	return DELTA_POS(x+other.x, y+other.y, z+other.z);
}
const DELTA_POS DELTA_POS::operator-(const DELTA_POS& other) const
{
	return DELTA_POS(x-other.x, y-other.y, z-other.z);
}
const DELTA_POS DELTA_POS::absoluteValue() const
{
    return DELTA_POS(abs(x), abs(y), abs(z));
}
const bool DELTA_POS::isZero() const
{
    return (x == 0) && (y == 0) && (z == 0);
}
const double DELTA_POS::magnitude() const
{
    return sqrt(x*x + y*y + z*z);
}
void DELTA_POS::print(ostream &out) const
{
    out << x << " " << y << " " << z << "\n";
}


bool deltaPosLessThan (DELTA_POS p1, DELTA_POS p2) {
    if(p1.x != p2.x)
        return p1.x < p2.x;
    else if(p1.y != p2.y)
        return p1.y < p2.y;
    else if(p1.z != p2.z)
        return p1.z < p2.z;
    else
        return false;
}
/* End of DELTA_POS functionality*/




/* Begin of CONN functionality */
/** Class used to represent a connection between DIJKSTRA_NODES. */

/* Create a CONN using the provided parameters */
CONN::CONN(int myFrom, int myTo, double len, double maxR, DELTA_POS deltaP){
    from       = myFrom;
    to         = myTo;
    length     = len;
    max_radius = maxR;
    deltaPos   = deltaP;
}

/* Output information about the connection to the provided output stream*/
void CONN::print(ostream &out) const{
    out << from << "->" << to << "   Length:" << length
    << "   Max radius:" << max_radius << "   Change in Unit Cell: ("
    << deltaPos.x << "," << deltaPos.y << "," << deltaPos.z << ")" << "\n";
}

/* End of CONN functionality*/




/* Begin of DIJKSTRA_NODE functionality */
/* Represents a node in a graph data structure created from a Voronoi network. Stores
 * information about the connections that stem from it as well as its position.
 */

/** Create a node using the provided parameters. */
DIJKSTRA_NODE::DIJKSTRA_NODE(int myID, double myX, double myY, double myZ, double maxR, bool active_flag){
    id = myID; x = myX; y = myY; z = myZ;
    connections = vector<CONN> ();
    max_radius = maxR;
    active = active_flag;
}

/** Output information about the node to the provided output stream. */
void DIJKSTRA_NODE::print(ostream &out) const {
    out << " Node info:" << "\n"
    << "    #: " << id << "    X: " << x << "    Y: " << y << "    Z:" << z << "\n"
    << "   Connections:" << "\n";
    for(unsigned int i = 0; i<connections.size(); i++){
        out << "     ";
        connections.at(i).print();
    }
}

/* End of DIJKSTRA_NODE functionality*/





/* Begin of DIJKSTRA_NETWORK functionality*/
/** Data structure that represents a set of nodes and edges created from a
 *  corresponding Voronoi network. Represents a single unit cell of the material.
 */

DIJKSTRA_NETWORK::DIJKSTRA_NETWORK(){
    nodes = vector<DIJKSTRA_NODE>();
}

/* Output information about the network to the provided output stream. */
void DIJKSTRA_NETWORK::print(ostream &out){
    for(unsigned int i = 0; i < nodes.size(); i++)
        nodes.at(i).print(out);
}

/* Build a DIJKSTRA_NETWORK from the provided VORONOI_NETWORK
 *  and store it using the provided pointer.
 */
void DIJKSTRA_NETWORK::buildDijkstraNetwork(const VORONOI_NETWORK *vornet, DIJKSTRA_NETWORK *dnet){
    vector<VOR_NODE> ::const_iterator niter = vornet->nodes.begin();
    int i = 0;
    dnet->nodes.clear();
    
    // Add copies of all nodes to the network
    while(niter != vornet->nodes.end()){
        DIJKSTRA_NODE node = DIJKSTRA_NODE(i, niter->x, niter->y, niter->z, niter->rad_stat_sphere, niter->active);
        i++;
        niter++;
        dnet->nodes.push_back(node);
    }
    
    // For each edge, store it in the DIJKSTRA node's list of connections that the connection
    // stems from
    vector<VOR_EDGE> ::const_iterator eiter = vornet->edges.begin();
    while(eiter != vornet->edges.end()){
        DELTA_POS pos = DELTA_POS(eiter->delta_uc_x, eiter->delta_uc_y, eiter->delta_uc_z);
        CONN conn = CONN(eiter->from, eiter->to, eiter->length, eiter->rad_moving_sphere,pos);
        dnet->nodes.at(conn.from).connections.push_back(conn);
        eiter++;
    }
    
    // Copy the unit cell vectors into the DIJKSTRA_NETWORK
    dnet->v_a = vornet->v_a;
    dnet->v_b = vornet->v_b;
    dnet->v_c = vornet->v_c;
}

/* Creates a new DIJKSTRA_NETWORK containing only the edges in which both terminuses are contained
 within the provided vector of ids.
 */
void DIJKSTRA_NETWORK::filterDnetEdges(vector<int> nodeIDs,
                                       VORONOI_NETWORK *origVornet,
                                       DIJKSTRA_NETWORK *newDnet){
    //VORONOI_NETWORK newVorNet;
    //VORONOI_NETWORK::filterVornetEdges(nodeIDs, origVornet, &newVorNet);
    VORONOI_NETWORK newVorNet = origVornet->filterEdges(nodeIDs);
    DIJKSTRA_NETWORK::buildDijkstraNetwork(&newVorNet, newDnet);
}

/* End of DIJKSTRA_NETWORK functionality*/





SEGCONN::SEGCONN(){}


