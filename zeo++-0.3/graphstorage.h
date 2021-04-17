#ifndef GRAPHSTORAGE_H
#define GRAPHSTORAGE_H

/* This file contains all of the data structures used when converting a
 * VORONOI_NETWORK into a similar DIJKSTRA_NETWORK.
 */

#include <iostream>
#include <vector>

#include "geometry.h"
#include "voronoicell.h"
#include "networkstorage.h"

/** Simple class used to store and manipulate set of three integers. */
class DELTA_POS{
public:
    int x, y, z;
    
    DELTA_POS(int myX = 0, int myY = 0, int myZ = 0);
    const bool equals(DELTA_POS other) const;
    const DELTA_POS operator*(const int& factor) const;
    const DELTA_POS operator*(const DELTA_POS& other) const;
    const DELTA_POS operator+(const DELTA_POS& other) const;
    const DELTA_POS operator-(const DELTA_POS& other) const;
    const DELTA_POS absoluteValue() const;
    const bool isZero() const;
    const double magnitude() const;
    void print(std::ostream &out = std::cout) const;
    //void print() const;
};

bool deltaPosLessThan (DELTA_POS p1, DELTA_POS p2);

/** Class used to represent a connection between DIJKSTRA_NODES. */
class CONN{
public:
    int from, to;       // IDs of nodes involved in connectioon
    double length;      // Length of edge
    double max_radius;  // Radius of largest spherical probe that can travel along edge
    DELTA_POS deltaPos; // Change in unit cell
    
    /* Create a CONN using the provided parameters */
    CONN(int myFrom, int myTo, double len, double maxR, DELTA_POS deltaP);
    
    /* Output information about the connection to the provided output stream*/
    void print(std::ostream& out = std::cout) const;
    //void print();
};

/* Represents a node in a graph data structure created from a Voronoi network. Stores
 * information about the connections that stem from it as well as its position.
 */
class DIJKSTRA_NODE {
public:
    int id;                    // Node id
    double x,y,z;              // (x,y,z) coordinate
    std::vector <CONN> connections; // List of connections that lead from the node to other nodes
    double max_radius;         // Maximum radius of particle that can sit at the node
    bool active;               // flag used to disactivate a node
    /* Creates a node with an invalid id, no connections and a radius of 0.*/
    /** Create a node using the provided parameters. */
    //DIJKSTRA_NODE();
    
    /* Creates a node. Default is with an invalid id, no connections and a radius of 0.*/
    DIJKSTRA_NODE(int myID = -1, double myX = 0, double myY = 0, double myZ = 0,
                  double maxR = 0, bool active_flag = true);
    
    /** Output information about the node to the provided output stream. Default is std::cout*/
    void print(std::ostream &out = std::cout) const;
};

/** Data structure that represents a set of nodes and edges created from a
 *  corresponding Voronoi network. Represents a single unit cell of the material.
 */
class DIJKSTRA_NETWORK {
public:
    std::vector<DIJKSTRA_NODE> nodes; // List of nodes
    XYZ v_a, v_b, v_c;           // Unit cell vectors of original unit cell
    
    DIJKSTRA_NETWORK();
    
    /* Output information about the network to the provided output stream. Default is std::cout*/
    void print(std::ostream &out = std::cout);
    
    static void buildDijkstraNetwork(const VORONOI_NETWORK *vornet, DIJKSTRA_NETWORK *dnet);
    static void filterDnetEdges(std::vector<int> nodeIDs, VORONOI_NETWORK *vornet, DIJKSTRA_NETWORK *dnet);
};


/* added by Maciek to handle segments */
class SEGCONN {
public:
    int from, to;         // nodes involved in connection
    int from_seg, to_seg; // segments that are connected
    double max_radius;
    double length;        // connection lenght
    int merged;           // connection to be merged
    DELTA_POS deltaPos;   // Change in unit cell
    
    SEGCONN();
};

#endif
