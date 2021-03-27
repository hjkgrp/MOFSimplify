#ifndef CHANNEL_H
#define CHANNEL_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <utility>

#include "graphstorage.h"
#include "networkstorage.h"
#include "heap.h"


/* Classes for topological analysis of pores */

/* The main class is PORE class, which will be used to create subsequent classes:
   CHANNEL class - inherit PORE plus adds functions to analyze channels
   POCKET class - inherit PORE plus adds functions to analyze inaccessible pockets

*/

class PORE {
public:
    std::map<int,int> idMappings;        // (node id, new index) pairs
    std::map<int,int> reverseIDMappings; // (new index, node id) pairs
    std::vector<DIJKSTRA_NODE> nodes;    // List of nodes in PORE indexed by new index
    std::vector<CONN> connections;       // List of connections in PORE

    std::vector<DELTA_POS> unitCells;    // List of displacements of unitcells in single PORE unit
    std::vector< std::vector<int> > ucNodes;  // List of nodes contained in each unit cell
    XYZ v_a, v_b, v_c;              // Unit cell vectors

    int dimensionality;             // Dimensionality of pore
    int basis [3][3];               // Basis vectors of pore


    /* Reconstructs the PORE by propagating paths until all nodes have been accessed. 
     * Stores each node in the unit cell in which it is encountered.
     * Attempts to minimize the number of unit cells required for reconstruction by 
     * favoring nodes in already-accessed unit cells over nodes in new unit cells. */
    void reconstruct();

    /* Constructs a pore from the provided  nodes in the DIJKSTRA_NETWORK. 
     * Reconstructs the pore by trying to minimize the number of unit cells
     * required to show a single pore unit.
     */
    PORE(std::vector<int> nodeIDs, DIJKSTRA_NETWORK *dnet, int dim, int basisVecs[3][3]);
  
    /* Create a PORE that does not contain any nodes or connections
     * and spans 0 unit cells.*/
    PORE();

    void storeNetwork(DIJKSTRA_NETWORK *);

    /* Converts the current channel into a DIJKSTRA_NETWORK */
    void buildDijkstraNetwork(DIJKSTRA_NETWORK *); //(unused)

    /* provides a vector with IDs of original voronoi nodes that correspond to the pore*/
    std::vector <int> nodeIds();

    /** Nodes within the provided DIJKSTRA_NETWORK can be classified as either accessible or inaccessible, 
     *  where sets of accessible  nodes constitute a CHANNEL and inaccessible POCKET. This function identifies
     *  the POREs that exist for the provided particle diameter and stores them using the provided pointer 
     *  to a vector of POREs. CHANNEL and PORE can be distunguished by dimentionality. 
     *  In addition, the pointer to the vector of bools is used to a store bool
     *  for each VORONOI_NODE, where infoStorage[i] is true iff node #i is accessible. */
    static void findChannelsAndPockets(DIJKSTRA_NETWORK *, std::vector<bool> *, std::vector<PORE> *);

    /** Nodes within the provided VORONOI_NETWORK can be classified as either accessible or inaccessible.
     *  This function does the same as the other findChannelsAndPockets but operates on VORONOI_NETWORK,
     *  the one provided as argument */
    static void findChannelsAndPockets(VORONOI_NETWORK *, double, std::vector<bool> *, std::vector<PORE> *);

    /** Prints all pore information (nodes, their positions and radii) to a file 
     *  Atom_network is needed to convert xyz coordinates to abc **/

    void printPoreSummary(std::ostream &out, ATOM_NETWORK *atmNet);

    /** Calculates center of mass and internal void radii (distance between the center of mass and its nearest atom)
     *  This function attempts pore reconstruction in cases where pores cross the cell boundaries. 
     */
    pair <XYZ, double> getCenterOfMass();

    /** Attempt to reconstruct a pore (dim==0) to get structue without PBC 
     */ 
    vector< pair <int,XYZ> >  getReconstructedPore();

    /**  Attempts to reconstructure a pore (dim==0) to get structure without PBC,
         for pores that cross cell boundry, multile copies of molecules are saved for each
         corresponding periodic image
     */

    vector<  vector< pair <int,XYZ> > >  getReconstructredPoresWithCrossBoundryCopies();

    /** Return the largest free sphere diameter **/
    double getIncludedSphereDiameter();

    /** Return the largest free and included along free sphere path for a path between two nodes **/

    pair <double,double> getFreeIncludedSphereDiameterforNodePair(int node1, int node2); 

    // function that fills a vector with data on pocket
    // Di, coordinates of Di, and radii that encapsulates the pocket
    void getSimplifiedPocketInfo(ATOM_NETWORK *, std::vector <double> *);


    /** PORE ENVELOPE CALCULATION **
     *  Functions that segment pore into sub fragments and calculate Df between sub fragments
     **/

     void getRestrictingDiameters(int nSegments, vector<int> vorNetID, vector< vector<double> > *PLDtable, vector< vector< pair<int,int> > > *PLDEdgeTable, 
                                 vector<double> *segmentDi, vector<int> *segmentDiNodeID, vector <double> *segmentDiFinal, vector<int> *segmentDiFinalNodeID);

};

class CHANNEL : public PORE {
public:

    CHANNEL(std::vector<int> nodeIDs, DIJKSTRA_NETWORK *dnet, int dim, int basisVecs[3][3])
           : PORE(nodeIDs, dnet, dim, basisVecs){};
  
    CHANNEL() : PORE(){};

    CHANNEL(PORE *);

    /* Prints information about the CHANNEL to the provided output stream, including
     * the number of nodes, unitcells, and the nodes located in each unit cell. Additional
     * node information is outputted if requested.*/
    void print(std::ostream &out, bool dispNodeInfo);

    /* Prints information about the CHANNEL to the standard output stream, including
     * the number of nodes, unitcells, and the nodes located in each unit cell. Additional
     * node information is outputted if requested.*/
    void print(bool dispNodeInfo);

    /** Write the commands necessary to draw the CHANNEL in ZeoVis 
     *  to the provided output stream. */
    void writeToVMD(int n, std::fstream &output);

    /** Write the commands necessary to draw the CHANNEL in ZeoVis 
     *  to the provided output stream. Includes a type because features and segments
     *  are drawn using the same command.*/
    void writeToVMD(std::string type, int n, std::fstream &output);

    /** Nodes within the provided DIJKSTRA_NETWORK can be classified as either accessible or inaccessible, 
     *  where sets of accessible  nodes constitute a CHANNEL. This function identifies
     *  the CHANNELs that exist for the provided particle diameter and stores them using the provided pointer 
     *  to a vector of channels. In addition, the pointer to the vector of bools is used to a store bool
     *  for each VORONOI_NODE, where infoStorage[i] is true iff node #i is accessible. */
    static void findChannels(DIJKSTRA_NETWORK *, std::vector<bool> *, std::vector<CHANNEL> *);

    /** Nodes within the provided VORONOI_NETWORK can be classified as either accessible or inaccessible, 
     *  where sets of accessible  nodes constitute a CHANNEL. This function identifies
     *  the CHANNELs that exist for the provided particle diameter and stores them using the provided pointer 
     *  to a vector of channels. In addition, the pointer to the vector of bools is used to a store bool
     *  for each VORONOI_NODE, where infoStorage[i] is true iff node #i is accessible. */
    static void findChannels(VORONOI_NETWORK *, double, std::vector<bool> *, std::vector<CHANNEL> *);
 
  
    /* Stores the ids of all atoms that bound this channel using the provided vector reference. An atom is considered
    *  to bound a channel if a node in the channel is a member of the atom's Voronoi cell. */
    void findBoundingAtoms(ATOM_NETWORK *, std::vector<BASIC_VCELL> &, std::vector<int> &);

    /* Returns true iff the CHANNEL unit can be depicted within one unit cell*/
    bool isUnicellular();

    /* Returns the largest free sphere diameter for the current channel */
    std::pair<double, std::pair <double,double> >  findFreeIncludedSphereDiameter();

    /* Return the largest free sphere starting from a partiular node 
       It uses a pair<pair> object that storage the current 'record' Di/df/dif
       This is to speed-up by early discart of paths with more restriction than the current path */
    std::pair<double, std::pair<double,double> >  findFreeIncludedSphereDiameterforNode(int, std::pair<double, std::pair<double,double> >);

};


/** Class POCKET is handling inaccessible pockets and have different functions than CHANNEL 
 *
 * */

class POCKET : public PORE {
public:

    POCKET(std::vector<int> nodeIDs, DIJKSTRA_NETWORK *dnet, int dim, int basisVecs[3][3])
           : PORE(nodeIDs, dnet, dim, basisVecs){};

    POCKET() : PORE(){};

// moved to PORE class    double getIncludedSphereDiameter();


};

 
/** Special class used to compare pair<int,DELTA_POS> instances when rebuilding
 *  CHANNEL instances. Used in conjunction with a HEAP, this class ensures that nodes 
 *  located in the current unit cell or in past unit cells are given preference over 
 *  other nodes when reconstructing the channel. */
class ReconstructorComparator{
    DELTA_POS currentPos;
    std::set<DELTA_POS, bool(*)(DELTA_POS,DELTA_POS)> positions;

public:
    /* Reset the visited positions and set the current position to the origin. */
    ReconstructorComparator();
  
    /* Set the current position and store the old position. */
    void setPosition(DELTA_POS p);

    bool compare(std::pair<int,DELTA_POS> p1, std::pair<int,DELTA_POS> p2); 
};





bool compareNodes(std::pair<int,DELTA_POS> p1, std::pair<int,DELTA_POS> p2);


#endif
