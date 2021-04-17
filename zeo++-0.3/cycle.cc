/* Cycles for each voronoi node are computed
 * Useful in material science where the voronoi nodes
 * are taken as tetrahedra interstitial sites
 * An octahedral interstitial site is computed as centroid of
 * voronoi nodes present in a cycle of length 4.
 * Author: Bharat Medasani
 * Date: Jan 02, 2013
 */
#include <vector>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cassert>

#include "mindist.h"
#include "geometry.h"
#include "networkstorage.h"
#include "graphstorage.h"
#include "sphere_approx.h"
#include "network.h"
#include "cycle.h"
#include "voronoicell.h"
#include "network.h"


using namespace std;

typedef vector<VOR_NODE>::iterator VnIt;
typedef vector<VOR_NODE>::const_iterator VnCit;
typedef vector<DIJKSTRA_NODE>::iterator DnIt;
typedef vector<DIJKSTRA_NODE>::const_iterator DnCit;
typedef vector<CONN>::iterator CnIt;


CYCLE::CYCLE() {length = 0;}

/* Function to compute girth of voronoi network with integer weights
 * Args:
 *  Input: 
 *      inp_vor: Pointer to input voronoi network
 *      range_flag: Denotes wether the weights have value other than 1
 *      range: Max. value of weights. weights = {1,2,...,range}
 *  Output:
 *      Returns length of smallest cycle
 */
int girth(const VORONOI_NETWORK* inp_vor, bool weight_flag, int wegith_range)
{
    return 0;

}

void edge_finder(const VORONOI_NETWORK* vornet, int node_index, vector<int>* edge_indices)
{
    typedef vector<VOR_EDGE>::const_iterator VeCit;
    int i = 0;
    VeCit it = vornet->edges.begin();
    for (; it != vornet->edges.end(); ++i, ++it){
        if (it->from == node_index || it->to == node_index) 
            edge_indices->push_back(i);
    }
    return;
}

/* At present I am assuming the least possible GIRTH to be 4 in a voronoi node
 * of regular crystals.
 */
static const int GIRTH = 4;

/* Function to compute a cycle of given length for each voronoi node
 * Args:
 *  Input:
 *      dnet: Pointer to input Dijkstra network
 *      cycle_length: Length of cycle to be computed
 *                    If 0, minimum length cycle is computed
 *      weight_flag: Denotes whether the weights have value other than 1
 *                   Ignored at present. Weight 1 is always used.
 *      weight_range: Max. value of weights. weights = {1,2,...,range}
 *                    Ignored at present. Range 1 is used.
 *  Output:
 *      Returns false if a cycle of given length cannot be found
 */
    /* Algorithm (Shortest cycle):
     * 1) For each Dijkstra node, identifiy all the edges from the node.
     * 2) For each edge from the node remove the edge and compute the
     * shortest path between the pair of nodes in the removed edge
     * using the Dijkstra's algorithm.
     */
/*
bool compute_cycle(DIJKSTRA_NETWORK* dnet, vector<CYCLE> *cycles,
                   int cycle_length, bool weight_flag, int weight_range)
{
    // First clear the active flag for each Dijkstra node
    // Set the max_radius of each connection to 1. Used in enabling or
    // disabling a connection without deleting the connection
    // Adds an additional check
    DnIt it = dnet->nodes.begin();
    for (; it != dnet->nodes.end(); ++it) {
        it->active = false;
        CnIt cit = it->connections.begin();
        for (; cit != it->connections.end(); ++cit) cit->max_radius = 1;
    }
    
    if (!cycle_length) {        // Shortest length cycle
        int girth = 0;
        for (int i = 0; i < dnet->nodes.size(); ++i){
            int short_cycle_length = 0;  //Shortest cycle length for this node
            
            DIJKSTRA_NODE origNode = dnet->nodes.at(i);
            CnIt it = origNode.connections.begin();
            int conn_ind = 0;
            for (; it != origNode.connections.end(); ++it, ++conn_ind){
                int path_length = 0;
                int orgNodeId = it->from;
                int endNodeId = it->to;
                
                map<int, int> visitedNodeDisplacement;
                vector<pair<int,int> > stack;
                stack.push_back(pair<int, int> (i, 0));
                visitedNodeDisplacement.insert(pair<int,int>(i,0));
                
                while (stack.size() != 0) {
                    // Remove the top-most node in stack
                    pair<int,int> nodeInfo = stack.back();
                    DIJKSTRA_NODE currNode = dnet->nodes.at(nodeInfo.first);
                    stack.pop_back();
                    
                    CnIt cit = currNode.connections.begin();
                    int conn_ind = 0;
                    for (; cit != currNode.connections.end(); ++cit, ++conn_ind){
                        //Disable the connection and the reverse connection
                        cit->max_radius = 0;
                        // Find the index of connection from to-node to from-node
                        CONN reverse_conn = *it;
                        reverse_conn.from = it->to;
                        reverse_conn.to = it->from;
                        
                        CnIt rev_cit;
                        rev_cit = find(dnet->nodes[it->to].connections.begin(),
                                       dnet->nodes[it->to].connections.end(),
                                       reverse_conn);
                        assert (rev_cit != dnet->nodes[it->to].connections.end());
                        rev_cit->max_radius = 0;
                        
                        //Find the shortest path between the nodes of connection
                        
                        
                        
                    }
                }
            }
        }
    }
    else{
        cout << "Functionality with given cycle length not implemented" << endl;
        return false;
    }
    return false;
}
*/

/* Function to compute a cycle of given length for each voronoi node
 * Args:
 *  Input: 
 *      vornet: Pointer to input voronoi network
 *      cycle_length: Length of cycle to be computed
 *                    If 0, minimum length cycle is computed
 *      weight_flag: Denotes wether the weights have value other than 1
 *      weight_range: Max. value of weights. weights = {1,2,...,range}
 *  Output:
 *      Returns false if the cycle of given length cannot be found
 */
/*
bool compute_cycle(const VORONOI_NETWORK* vornet, vector<CYCLE> *cycles,
                   int cycle_length, bool weight_flag, int weight_range)
{
    //Build graph data structure
    DIJKSTRA_NETWORK dnet;
    DIJKSTRA_NETWORK::buildDijkstraNetwork(vornet, &dnet);
    return compute_cycle(&dnet, cycles, cycle_length, weight_flag, weight_range);
}
*/

/* Function to compute a cycle of given length for each voronoi node
 * Args:
 *  Input:
 *      dnet: Pointer to input Dijkstra network
 *      cycle_length: Length of cycle to be computed
 *                    If 0, minimum length cycle is computed
 *      weight_flag: Denotes whether the weights have value other than 1
 *                   Ignored at present. Weight 1 is always used.
 *      weight_range: Max. value of weights. weights = {1,2,...,range}
 *                    Ignored at present. Range 1 is used.
 *  Output:
 *      Returns false if a cycle of given length cannot be found
 */
bool compute_4cycle(DIJKSTRA_NETWORK* dnet, vector<CYCLE> *cycles,
                    bool weight_flag, int weight_range)
{
    /* Algorithm (Shortest cycle):
     * 1) For each Dijkstra node, identifiy all the edges from the node.
     * 2) For each edge from the node remove the edge and compute the
     * shortest path between the pair of nodes in the removed edge
     * using the Dijkstra's algorithm.
     */
    // First clear the active flag for each Dijkstra node
    // Set the max_radius of each connection to 1. Used in enabling or
    // disabling a connection without deleting the connection
    // Adds an additional check
    DnIt it = dnet->nodes.begin();
    int count = 0;
    for (; it != dnet->nodes.end(); ++it, ++count) {
        it->active = false;
        CnIt cit = it->connections.begin();
        for (; cit != it->connections.end(); ++cit) cit->max_radius = 1;
    }
    cout << "Length of vornet: " << count << endl;

    it = dnet->nodes.begin();
    for (; it != dnet->nodes.end(); ++it) {
        //cout << "x: " << it->x << " y: " << it->y << " z: " <<  it->z << endl;
        ;
    }
    
    for (int i = 0; i < dnet->nodes.size(); ++i){
        DIJKSTRA_NODE org_node = dnet->nodes.at(i);
        if (i == 0) cout << "org_node: "<< org_node.x << " " << org_node.y << " " << org_node.z << endl;
        CnIt it = org_node.connections.begin();
        int conn_ind = 0;
        for (; it != org_node.connections.end(); ++it, ++conn_ind){
            int path_length = 0;
            int org_node_id = it->from;
            int end_node_id = it->to;
            if (i == 0) cout << "end_node_id: " << end_node_id << endl;
            DIJKSTRA_NODE end_node = dnet->nodes.at(end_node_id);
            
            if (end_node.active) continue;

            CnIt oit = it + 1;       //Check the other connections of begin node
            int oconn_ind=conn_ind + 1;
            for (; oit != org_node.connections.end(); ++oit, ++oconn_ind) {
                int inter_node_id = oit->to;
                if (inter_node_id==org_node_id || inter_node_id==end_node_id) continue;
		if (i==0) cout << "    inter_node_id: " << inter_node_id << endl;
                DIJKSTRA_NODE inter_node = dnet->nodes.at(inter_node_id);
                if (inter_node.active) continue;
                for (CnIt iit = inter_node.connections.begin(); iit != inter_node.connections.end(); ++iit) 
                    if (i==0) cout << "        to: "<< iit->to << endl;
                CnIt eit = end_node.connections.begin();
                int econn_ind=0;
                for (; eit != end_node.connections.end(); ++eit, ++econn_ind) {
                    int final_node_id = eit->to;
                    if (final_node_id==end_node_id || final_node_id==org_node_id) continue;
                    if (i==0) cout << "    final_node_id: " << final_node_id << endl;
                    DIJKSTRA_NODE final_node = dnet->nodes.at(final_node_id);
                    if (final_node.active) continue;

                    CnIt iit = inter_node.connections.begin();
                    for (; iit != inter_node.connections.end(); ++iit) {
                        if (iit->to == final_node_id){  //Cycle found
                            CYCLE new_cycle;
                            new_cycle.length=4;
                            new_cycle.nodes.push_back(org_node);
                            new_cycle.nodes.push_back(end_node);
                            new_cycle.nodes.push_back(inter_node);
                            new_cycle.nodes.push_back(dnet->nodes.at(final_node_id));
                            cycles->push_back(new_cycle);
                        }
                    }
                }
            }
        }
        org_node.active = true;
    }
    cout << "No. of cycles: " << cycles->size() << endl;
    if (!cycles->size()) return false;
    return true;
}

/* Function to compute a cycle of given length for each voronoi node
 * Args:
 *  Input: 
 *      vornet: Pointer to input voronoi network
 *      cycle_length: Length of cycle to be computed
 *                    If 0, minimum length cycle is computed
 *      weight_flag: Denotes wether the weights have value other than 1
 *      weight_range: Max. value of weights. weights = {1,2,...,range}
 *  Output:
 *      Returns false if the cycle of given length cannot be found
 */
bool compute_4cycle(VORONOI_NETWORK* vornet, vector<CYCLE> *cycles,
                    bool weight_flag, int weight_range)
{
    //Build graph data structure
    DIJKSTRA_NETWORK dnet;
    DIJKSTRA_NETWORK::buildDijkstraNetwork(vornet, &dnet);
    return compute_4cycle(&dnet, cycles, weight_flag, weight_range);
}

void centroid(const CYCLE* const cycle, XYZ* center, vector<int>* node_ids)
{
    DnCit it = cycle->nodes.begin();
    //cout << "Length of nodes in the cycle: " << cycle->nodes.size() << endl;

    *center = XYZ(0,0,0);
    node_ids->clear();
    //cout << "Length of node_ids vector: " << node_ids->size() << endl;
    for (; it != cycle->nodes.end(); ++it) {
        XYZ new_pt (it->x, it->y, it->z);
        *center = *center + new_pt;
	node_ids->push_back(it->id);
    }
    //cout << "Length of node_ids vector: " << node_ids->size() << endl;

    *center = center->scale(1.0/cycle->length);

}


/* Function to compute the voronoi face centers
 * Arguments:
 *  atmnet (Input) : Unit cell containing atoms/ions in ATOM_NETWORK class
 *  face_centers: (Output), pointer to vector of XYZ objects holding the 
 *                coordinates of centroids of voronoi faces
 */
void face_center(ATOM_NETWORK* atmnet, vector<XYZ>* face_centers)
{
    VORONOI_NETWORK vornet;
    vector<VOR_CELL> vcell;
    vector<BASIC_VCELL> bvcell;
    performVoronoiDecomp(true, atmnet, &vornet, &vcell, true, &bvcell);

    vector<VOR_CELL>::iterator vc_it;
    vector<VOR_FACE>::iterator vf_it;

    int vor_face_count = 0;
    int vor_cell_count = 0;

    for (vc_it = vcell.begin(); vc_it != vcell.end(); ++vc_it){
        for (vf_it = vc_it->faces.begin(); vf_it != vc_it->faces.end(); ++vf_it){
            ++vor_face_count;
            vector<int>::iterator iit;
            if (vf_it->nodeIDs.size() > 4) {
                continue;
            }
            vector<Point>::iterator ip_it;
            cout << "Orderd vertices in the face: " << endl;
            for (ip_it = vf_it->orderedVertices.begin(); ip_it != vf_it->orderedVertices.end(); ++ip_it){;
                //cout << "(";
                //ip_it->print(cout);
                //cout << ")  " ;
            }
            //cout <<  endl;
        }
    }
    cout << "VOR_FACE_COUNT " << vor_face_count << endl;
    vcell.clear();
    bvcell.clear();
}
