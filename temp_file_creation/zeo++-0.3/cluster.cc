/* Clusters of voronoi nodes are analyzed and then compacted to a single 
 * voronoi node for each cluster. 
 * Different methods of pruning. 
 * Useful in material science where the high symmetry voronoi nodes
 * are taken as interstitial sites
 * Author: Bharat Medasani
 * Date: Dec 09, 2013
 */
#include <vector>
#include <set>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <cmath>

#include "mindist.h"
#include "geometry.h"
#include "networkstorage.h"
#include "graphstorage.h"
#include "sphere_approx.h"
#include "network.h"
#include "cluster.h"


using namespace std;

typedef vector<ATOM>::iterator AtmIt;
typedef vector<VOR_NODE>::iterator VnIt;
typedef vector<VOR_EDGE>::iterator VeIt;
typedef vector<VOR_NODE>::const_iterator VnCit;
typedef vector<DIJKSTRA_NODE>::iterator DnIt;
//typedef vector<DIJKSTRA_NODE>::const_iterator DnCit;

/* conn_comp 
 * Compares two objects of CONN class defined in graphstorage.h
 * Ensure that from has same value in both conn objects for proper operation
 */
bool edge_comp(const CONN& lhs, const CONN& rhs) 
{
    if (lhs.from != rhs.from) throw 11;  // from nodes should be equal
    return lhs.length < rhs.length;
}


/* Sort the dijkstra network such that for each node all the connections
 * are arranged in order based on distance 
 */
void dijkstra_sort(DIJKSTRA_NETWORK* dnet)
{
    vector<DIJKSTRA_NODE>::iterator it = dnet->nodes.begin();
    for ( ; it != dnet->nodes.end(); ++it) 
        sort(it->connections.begin(), it->connections.end(), edge_comp);
    return;
}


pair<DnIt, DnIt> get_first_closer_nonassigned_node(const XYZ pt, const ATOM_NETWORK& atmnet,
        vector<DIJKSTRA_NODE> nodes, const float cutoff_dist)
{
    DnIt it = nodes.begin(); 
    for (; it != nodes.end(); ++it){
        double period_dist = atmnet.calcDistanceXYZ(it->x, it->y, it->z, pt.x, pt.y, pt.z);
        if (!it->active && period_dist < cutoff_dist) break;
    }
    return make_pair(it, nodes.end());
}


/* Function to partition high accuracy voronoi graph into clusters that are close to 
 * voronoi nodes in low accuracy cluster.
 * Args:
 *  Input: 
 *      atmnet: ATOM_NETWORK
 *      cutoff_dist: Maximum distance from regular voronoi network node that belong to a cluster
 *  Output:
 *      Returns clusters as vector of VORONOI_NODES
 */
vector< vector<XYZ> > cluster_partition(ATOM_NETWORK* atmnet, const float cutoff_dist)
{
    //Algorithm:
    
    // Partition a high accuracy voronoi net into clusters.
    ATOM_NETWORK ha_atmnet = *atmnet;
    string acc_setting = "S30";
    setupHighAccuracyAtomNetwork(&ha_atmnet, acc_setting);

    vector<VOR_CELL> vcell;
    vector<BASIC_VCELL> bvcell;
    VORONOI_NETWORK vornet, ha_vornet;
    performVoronoiDecomp(true, atmnet, &vornet, &vcell, false, &bvcell);
    vcell.clear();
    bvcell.clear();
    performVoronoiDecomp(true, &ha_atmnet, &ha_vornet, &vcell, false, &bvcell);
    DIJKSTRA_NETWORK dnet;
    DIJKSTRA_NETWORK::buildDijkstraNetwork(&ha_vornet, &dnet);

    for (DnIt it = dnet.nodes.begin(); it != dnet.nodes.end(); ++it)
        it->active = false;
    
    for (VnIt it = vornet.nodes.begin(); it != vornet.nodes.end(); ++it)
        cout << it->x << " " << it->y << " " << it->z << endl;

    vector <vector<XYZ> > clusters;
    
    for (VnIt nodeit = vornet.nodes.begin(); nodeit != vornet.nodes.end(); ++nodeit) {
        // Get the first non assigned high accuracy voronoi node within the 
        // cutoff distance to current low accuracy voronoi node
        XYZ low_acc_node_loc(nodeit->x, nodeit->y, nodeit->z);
        pair<DnIt, DnIt> p = get_first_closer_nonassigned_node(   // Highaccuracy dijkstra
                low_acc_node_loc, *atmnet, dnet.nodes, cutoff_dist);

        // Initialize cluster and tovisit_nodes
        vector<XYZ> cluster;
        vector<int> tovisit_ids;
        vector<int> visited_ids;

        // Add the current high accuracy node to tovisit ids
        tovisit_ids.push_back(p.first->id);

        // Populate the cluster
        while (!tovisit_ids.empty()){
            //Work on the current node
            int curr_id = tovisit_ids.back();
            DIJKSTRA_NODE curr_nd = dnet.nodes.at(curr_id);
            cluster.push_back(XYZ(curr_nd.x, curr_nd.y, curr_nd.z));
            visited_ids.push_back(curr_nd.id);
            tovisit_ids.pop_back();

            //If any connected nodes are unvisited and are within cutoff distance,
            //push them to tovisit_ids
            vector<CONN>::iterator connit = curr_nd.connections.begin();
            for (; connit != curr_nd.connections.end(); ++connit){
                int to_nodeid = connit->to;
                DIJKSTRA_NODE to_node = dnet.nodes.at(to_nodeid);
                double dist_to_low_acc_node = atmnet->calcDistanceXYZ(
                        nodeit->x, nodeit->y, nodeit->z,
                        to_node.x, to_node.y, to_node.z);
                if (find(visited_ids.begin(), visited_ids.end(), to_nodeid) == visited_ids.end() 
                        && dist_to_low_acc_node < cutoff_dist)
                    tovisit_ids.push_back(to_nodeid); 
            }
        }
        
        clusters.push_back(cluster);
        tovisit_ids.clear();
        visited_ids.clear();
        cluster.clear();
        
    }
    // assert (clusters.size() == vornet.nodes.size());

    return clusters;
}


/* cluster_aggregate: Reduces a set of points to a single point.
 * The reduction is done thr' computing the centroid of the points */
vector<XYZ> cluster_aggregate(const vector< vector<XYZ> >& clusters, const ATOM_NETWORK* atmnet)
{
    MIN_PER_DISTANCE distcalc = atmnet->getDistCalc();
    vector<XYZ> reduced_clusters;
    vector <vector<XYZ> >::const_iterator clust_it = clusters.begin();


    for (; clust_it != clusters.end(); ++clust_it) {
        vector<XYZ>::const_iterator it = clust_it->begin();
        XYZ centroid(it->x, it->y, it->z);  // Assign the centroid as first point
        Point c_pt = atmnet->xyz_to_abc(centroid); // Fractional coords
        int count = 1;  //measures how many points are added in c_pt
        ++it;
        for (; it != clust_it->end(); ++it) {
            // Find the closest periodic of *it to centroid and add it to it.
            Point curr_pt = atmnet->xyz_to_abc(it->x, it->y, it->z);
            double tmp_x, tmp_y, tmp_z;    // C.P.I. of curr_pt w.r.t. centroid
            distcalc.closest_periodic_image(c_pt.vals[0], c_pt.vals[1], c_pt.vals[2], 
                                 curr_pt.vals[0], curr_pt.vals[1], curr_pt.vals[2], 
                                 tmp_x, tmp_y, tmp_z);
            c_pt = c_pt + Point(tmp_x, tmp_y, tmp_z);
            ++count;
        }
        c_pt = atmnet->abc_to_xyz(c_pt);     // Convert to absolute coordinates
        c_pt = c_pt.scale(1.0/double(count));
        reduced_clusters.push_back(XYZ(c_pt.vals[0], c_pt.vals[1], c_pt.vals[2]));
    }
    return reduced_clusters;
}

/* cluster_aggregate: Reduces a set of points to a single point.
 * The reduction is done thr' computing the centroid of the points */
void cluster_aggregate(const vector< vector<XYZ> >& clusters, const ATOM_NETWORK* atmnet, vector<XYZ>* reduced_clusters)
{
    MIN_PER_DISTANCE distcalc = atmnet->getDistCalc();
    //vector<XYZ> reduced_clusters;
    vector <vector<XYZ> >::const_iterator clust_it = clusters.begin();


    for (; clust_it != clusters.end(); ++clust_it) {
        vector<XYZ>::const_iterator it = clust_it->begin();
        XYZ centroid(it->x, it->y, it->z);  // Assign the centroid as first point
        Point c_pt = atmnet->xyz_to_abc(centroid); // Fractional coords
        int count = 1;  //measures how many points are added in c_pt
        ++it;
        for (; it != clust_it->end(); ++it) {
            // Find the closest periodic of *it to centroid and add it to it.
            Point curr_pt = atmnet->xyz_to_abc(it->x, it->y, it->z);
            double tmp_x, tmp_y, tmp_z;    // C.P.I. of curr_pt w.r.t. centroid
            distcalc.closest_periodic_image(c_pt.vals[0], c_pt.vals[1], c_pt.vals[2], 
                                 curr_pt.vals[0], curr_pt.vals[1], curr_pt.vals[2], 
                                 tmp_x, tmp_y, tmp_z);
            c_pt = c_pt + Point(tmp_x, tmp_y, tmp_z);
            ++count;
        }
        c_pt = atmnet->abc_to_xyz(c_pt);     // Convert to absolute coordinates
        c_pt = c_pt.scale(1.0/double(count));
        reduced_clusters->push_back(XYZ(c_pt.vals[0], c_pt.vals[1], c_pt.vals[2]));
    }
    return;
}

void simplify_ha_vornet(ATOM_NETWORK* atmnet)
{
    vector <vector<XYZ> > clusters =  cluster_partition(atmnet, 0.2);
    vector<XYZ> reduced_clusters = cluster_aggregate(clusters, atmnet);
    
    vector<XYZ>::iterator it = reduced_clusters.begin();
    for (; it != reduced_clusters.end(); ++it) {
        it->print(cout);
        ;
    }
    return;
}

void high_accuracy_vornodes_reduction(ATOM_NETWORK* atmnet, Vector_XYZ* vornodes)
{
    vector <vector<XYZ> > clusters =  cluster_partition(atmnet, 0.2);
    cluster_aggregate(clusters, atmnet, &(vornodes->nodes));
}

void high_accuracy_vornodes_reduction(ATOM_NETWORK* atmnet, vector<XYZ>* vornodes)
{
    vector <vector<XYZ> > clusters =  cluster_partition(atmnet, 0.2);
    cluster_aggregate(clusters, atmnet, vornodes);
}

/* Function to compact voronoi network
 * Args:
 *  Input: 
 *      inp_vor: Pointer to input voronoi network
 *      cluster_rad: Maximum distance between nodes that belong to a cluster
 *                   Default is 0.5 Angstrom    // May be too high. Confirm
 *  Output:
 *      Returns compacted voronoi network
 */
/*
VORONOI_NETWORK cluster_reduce(const VORONOI_NETWORK* vornet, const float cutoff_dist)
{
    // Algorithm:
    // Identify the clusters.
    // For each cluster find the voronoi node with maximum circle of radius
    // Eliminate all other nodes in the cluster
    // Generate new edges linking the remaining voronoi node in the cluster to other clusters
    
    // Partition the dijkstra net into clusters.
    DIJKSTRA_NETWORK dnet;
    DIJKSTRA_NETWORK::buildDijkstraNetwork(vornet, &dnet);
    dijkstra_sort(&dnet);

    vector<DIJKSTRA_NODE>::iterator nodeiter = dnet.nodes.begin();
    vector< vector<int> > clusters;
    while (nodeiter != dnet.nodes.end()) {
        if (!nodeiter->active){ // Node still not assigned to any cluster
            // Check if the node belongs to any existing cluster by first checking if a 
            // connecting node belongs to cluster and then test the distance to connecting node
            for (vector<CONN>::const_iterator conniter = nodeiter->connections.begin();
                    conniter != nodeiter->connections.end(); ++conniter){
                int conn_node = conniter->to;
                for (vector< vector<int> >::iterator clusteriter = clusters.begin();
                        clusteriter != clusters.end(); ++clusteriter) {
                    //If a connecting node belongs to cluster check otherwise skip
                    vector<int>::iterator it = find(clusteriter->begin(), clusteriter->end(), conn_node);
                    if (it != clusteriter->end() && conniter->length <= cutoff_dist){
                        clusteriter->push_back(nodeiter->id);
                        nodeiter->active = true;
                        break;
                    }
                }
                if (nodeiter->active) break;
            }
            if (!nodeiter->active) { // Create a new cluster and assign node to new cluster
                vector<int> cluster;
                cluster.push_back(nodeiter->id);
                clusters.push_back(cluster);
            }
        }
        ++nodeiter;
    }

    VORONOI_NETWORK net;
    return net;

}
*/

/* Function to prune high accuracy voronoi network.
 * Removes the voronoi nodes within the higher radius atoms
 * Args:
 *  Input: 
 *      low_atm_net: Pointer to original atom network
 *      ha_atm_net: Pointer to high accuracy atom network
 *  Input/Output:
 *      ha_vor: Pointer to high accuracy voronoi network
 */
void prune_high_accuracy_voronoi_network( VORONOI_NETWORK* ha_vor, 
                            ATOM_NETWORK* low_atm_net, ATOM_NETWORK* high_atm_net,
                            double delta, bool print)
{
    double minr, maxr;

    /* Analyze the original atom network*/
    minr = 100000;      // A very big number as initial min. value
    maxr = -100000;     // " " " " "  negative number as init max. value
    for (AtmIt it = low_atm_net->atoms.begin(); it != low_atm_net->atoms.end(); ++it){
        if (it->radius < minr) minr = it->radius;
        if (it->radius > maxr) maxr = it->radius;
    }
    /* Print the high accuracy atom network for debug */
    //vector<int>::iterator id_it;
    //for (id_it = high_atm_net->IDmapping.begin(); id_it != high_atm_net->IDmapping.end(); ++id_it){
    //    cout << *id_it << endl;
    //}

    /* Print the sizes of high accuracy and regular atom networks for debug */
    cout << "Size of regular atom network " << low_atm_net->atoms.size() << endl;
    cout << "Size of high accuracy atom network " << high_atm_net->atoms.size() << endl;

    if (print) {
        cout << "Radii analysis:" << endl;
        cout << "the smallest atom r = " << minr << ", while the largest atom r = " << maxr << endl;

        cout << "Length of vornet nodes before pruning: " << ha_vor->nodes.size() << endl;
        cout << "Length of vornet edges before pruning: " << ha_vor->edges.size() << endl;
    }
    MIN_PER_DISTANCE distcalc = low_atm_net->getDistCalc();
    for (AtmIt it = low_atm_net->atoms.begin(); it != low_atm_net->atoms.end(); ++it){
        double atm_rad = it->radius;
        if (atm_rad == minr) continue;
        if (print) {
            cout << "Atom radius" <<  atm_rad << endl;
            cout << "Atoms locations: " << endl;
            cout << it->x <<  ", " <<  it->y << ", " <<  it->z << ", " <<  endl;
        }

        double a1 = it->a_coord; double b1 = it->b_coord; double c1 = it->c_coord;
        for (VnIt vnit = ha_vor->nodes.begin(); vnit != ha_vor->nodes.end();) {
            double x = vnit->x; double y = vnit->y; double z = vnit->z;
            double dist = low_atm_net->calcDistanceXYZABC(x, y, z, a1, b1, c1);
            if (dist <= atm_rad-delta){       // Delete the voronoi node and associated edges 
                int node_id = vnit - ha_vor->nodes.begin();
                // Delete the edges with the node_id and decrement the node index when its
                // greater than the id of deleted node
                for (VeIt veit = ha_vor->edges.begin(); veit != ha_vor->edges.end();) {
                    int from = veit->from;
                    int to = veit->to;
                    if (from == node_id || to == node_id) {
                        veit = ha_vor->edges.erase(veit);
                        //cout << "current veit index" << veit - ha_vor->edges.begin();
                    }
                    else{
                        if (from > node_id) --(veit->from);
                        if (to > node_id) --(veit->to);
                        ++veit;
                        //cout << "current veit index" << veit - ha_vor->edges.begin();
                    }
                }
                //cout << "Length of vornet edges after deleting node: " << ha_vor->edges.size() << endl;
                vnit = ha_vor->nodes.erase(vnit);
            }
            else
                ++vnit;
        }
        if (print) {
            cout << "Length of vornet nodes after pruning: " << ha_vor->nodes.size() << endl;
            cout << "Length of vornet edges after pruning: " << ha_vor->edges.size() << endl;
        }
    }
    return;
}

/* Function to identify the nearest high accuracy voronoi node that has a high radius.
 * Removes the voronoi nodes within the higher radius atoms
 * Args:
 *  Input: 
 *      ha_vornet: Pointer to high accuracy voronoi network
 *      vornet: Pointer to original voronoi network
 *      atm_net: Pointer to atom network
 *      cutoff: Radius of cutoff sphere from each voronoi ode
 *  Output:
 *      red_vornet: Pointer to reduced voronoi network
 */
void nearest_largest_diameter_ha_vornet(VORONOI_NETWORK* ha_vornet, VORONOI_NETWORK* vornet, 
        ATOM_NETWORK* atmnet, VORONOI_NETWORK* red_vornet, float cutoff)
{
    cout << "vornet size " << ha_vornet->nodes.size() << endl;
    for (VnIt vit = vornet->nodes.begin(); vit != vornet->nodes.end(); ++vit){
        VOR_NODE* biggest_near_vornode = NULL;
        double x2 = vit->x;
        double y2 = vit->y;
        double z2 = vit->z;
        for (VnIt vh_it = ha_vornet->nodes.begin(); vh_it != ha_vornet->nodes.end(); ++vh_it){
            double x1 = vh_it->x;
            double y1 = vh_it->y;
            double z1 = vh_it->z;
            double dist = atmnet->calcDistanceXYZ(x1,y1,z1,x2,y2,z2);
            //cout << dist << endl;
            if (dist <= cutoff){
                if (!biggest_near_vornode)
                    biggest_near_vornode = &*vh_it;
                else{
                    if (vh_it->rad_stat_sphere > biggest_near_vornode->rad_stat_sphere) {
                        biggest_near_vornode = &*vh_it;
                    }
                }
                //cout << "Size of vornode " << sizeof(biggest_near_vornode) << endl;
            }
        }
        //cout << "Reduced vornet nodes size " << red_vornet->nodes.size() << endl;
        //cout << "Size of vornode " << sizeof(biggest_near_vornode) << endl;
        if (biggest_near_vornode){
            red_vornet->nodes.push_back(*biggest_near_vornode);
        }
        else{
            cout << "Not able to find closest ha node" << endl;
            //throw 10;
        }
        //cout << "Size of vornode " << sizeof(biggest_near_vornode) << endl;
    }

    /* Ignoring edges at the moment as they are not needed */
    // Now that the voronoi nodes in reduced voronoi network are in the order 
    // of low accuracy vornet, generate the edges connecting the nodes in the order corresponding
    // to that in the low accuracy voronoi network. 
    // Adjust parameters of edges such as length, radius etc
}

/* Function to simplify the high accuracy voronoi network such that 
 * the small spheres coordinating the vornodes originate from different
 * atoms in the original atom network.  
 * Args:
 *  Input: 
 *      ha_vornet: Pointer to high accuracy voronoi network
 *      ha_atment: Pointer to high accuracy atom network
 *  Output:
 *      red_vornet: Pointer to reduced voronoi network
 */
void simplify_high_accuracy_vornet(VORONOI_NETWORK* ha_vornet, 
                                   ATOM_NETWORK* ha_atmnet, 
                                   VORONOI_NETWORK* red_vornet)
{
    // Identify nodes that belong to unitcell crossing edges
    vector<int> nodes_with_boundary_crossing_edges;
    for (VeIt vit = ha_vornet->edges.begin(); vit != ha_vornet->edges.end(); ++vit){
        if (vit->delta_uc_x || vit->delta_uc_y || vit->delta_uc_y){
            nodes_with_boundary_crossing_edges.push_back(vit->from);
            nodes_with_boundary_crossing_edges.push_back(vit->to);
        }
    }
    // Sort and remove duplicates
    vector<int>::iterator vbeg = nodes_with_boundary_crossing_edges.begin(); 
    vector<int>::iterator vend = nodes_with_boundary_crossing_edges.end(); 
    int size = nodes_with_boundary_crossing_edges.size(); 
    sort(vbeg, vbeg+size);
    nodes_with_boundary_crossing_edges.erase( unique( vbeg, vend ), vend );

    vector<int>* id_map = &(ha_atmnet->IDmapping);
    int i = 0;
    for (VnIt vit = ha_vornet->nodes.begin(); vit != ha_vornet->nodes.end(); ++vit, ++i){

        if (find(vbeg, vend, i) != vend){
            red_vornet->nodes.push_back(*vit);
        }
        else{

            set<int> original_atmnet_ids;
            for (vector<int>::iterator it = vit->atomIDs.begin(); it != vit->atomIDs.end(); ++it){
                original_atmnet_ids.insert(id_map->at(*it));
            }
            
            if (original_atmnet_ids.size() >= 4) {
                red_vornet->nodes.push_back(*vit);
            }
        }
    }
}


/* Function to prune the high accuracy voronoi network based on geometry
 * such that within a 0.1Ang^2 grid only one voronoi node is retained. 
 * Implemented that such that no two nodes are with 0.1 Ang distance 
 * Args:
 *  Input: 
 *      ha_vornet: Pointer to high accuracy voronoi network
 *      ha_atment: Pointer to high accuracy atom network
 *  Output:
 *      red_vornet: Pointer to reduced voronoi network
 */
void geometry_pruning(VORONOI_NETWORK* ha_vornet, ATOM_NETWORK* atmnet, 
                      float cutoff, VORONOI_NETWORK* red_vornet)
{
    for (VnIt vh_it = ha_vornet->nodes.begin(); vh_it != ha_vornet->nodes.end(); ++vh_it){
        if (red_vornet->nodes.size() == 0) {
            red_vornet->nodes.push_back(*vh_it);
        }
        else {
            vector<double> dists;
            double x1 = vh_it->x;
            double y1 = vh_it->y;
            double z1 = vh_it->z;
            for (VnIt vit = red_vornet->nodes.begin(); vit != red_vornet->nodes.end(); ++vit){
                double x2 = vit->x;
                double y2 = vit->y;
                double z2 = vit->z;
                double dist = atmnet->calcDistanceXYZ(x1,y1,z1,x2,y2,z2); 
                dists.push_back(dist);
            }
            sort(dists.begin(), dists.end());
            if (dists[0] > cutoff) {
                red_vornet->nodes.push_back(*vh_it);
            }
        }
    }

    cout << "size of reduced vornet " << red_vornet->nodes.size() << endl;

}

/* Function to prune the high accuracy voronoi network such that nodes within 
 * the original atoms are pruned. Nodes within 0.1Ang^2 from surface are 
 * retained. 
 * Args:
 *  Input: 
 *      ha_vornet: Pointer to high accuracy voronoi network
 *      atment: Pointer to original atom network
 *  Output:
 *      red_vornet: Pointer to reduced voronoi network
 */
void ha_prune_within_atom(VORONOI_NETWORK* ha_vornet, ATOM_NETWORK* atmnet, 
                      float cutoff, VORONOI_NETWORK* red_vornet)
{
    for (VnIt vh_it = ha_vornet->nodes.begin(); vh_it != ha_vornet->nodes.end(); ++vh_it){
        double x2 = vh_it->x;
        double y2 = vh_it->y;
        double z2 = vh_it->z;
        bool near_flag = false;
        for (AtmIt ait = atmnet->atoms.begin(); ait != atmnet->atoms.end(); ++ait) {
            double x1 = ait->x;
            double y1 = ait->y;
            double z1 = ait->z;
            double dist = atmnet->calcDistanceXYZ(x1,y1,z1,x2,y2,z2); 
            if (dist < ait->radius-cutoff) {
                near_flag = true;
                break;
            }
        }
        if (not near_flag) 
            red_vornet->nodes.push_back(*vh_it);
    }

    cout << "size of reduced vornet " << red_vornet->nodes.size() << endl;

}
