#ifndef NETWORKANALYSIS_H
#define NETWORKANALYSIS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <utility>

#include "graphstorage.h"

/** Data structure used to represent the path taken when traversing
 *  a DIJKSTRA_NETWORK. 
 */
class PATH {
 public:
     std::vector<DIJKSTRA_NODE> visitedNodes; // List of visited nodes
  DIJKSTRA_NODE currentNode;          // Next node to be visited
  double max_radius;                  // Maximum radius of sphere that can travel along PATH
  double max_inc_radius;              // Maximum included sphere along PATH
  double length;                      // Current length of PATH
  bool toSink;                        // True iff next node to be visited is a sink node
  std::vector<int> visitedIDs;             // IDs of all visited nodes
  std::set<int> visitedSourceIDs;          // IDs of visited source nodes
 
  /* Create an empty PATH with invalid length and radii. */
  PATH();

  /* Create a PATH with the provided current node, radius and length.*/
  PATH(DIJKSTRA_NODE cNode, double maxRadius, double maxIncRadius, double len);

  /* Prints the information about the edges in the DIJKSTRA_NETWORK that were
  *  used when following the PATH. */
  void printPathEdges (std::vector<int> nodeList, DIJKSTRA_NETWORK *dnet);

  /* Print information about the PATH to the provided output stream.
   * If none given, default is the standard output stream. */
  void print(std::ostream &output = std::cout);
};


/* Returns true iff the first PATH has a lower maximu radius than the second PATH.*/
bool hasLowerMaxR(PATH n1, PATH n2);


/* Data structure used when studying paths across a DIJKSTRA_NETWORK in 
 * one of the unit cell directions. Source nodes are identified as the 
 * terminus of a connection with a -1 in the direction of interest. Connections
 * are also classified in a similar manner as source connections (-1 change in uc in d.o.i.), 
 * sink connections (+1 in d.o.i.) or a regular (0 in d.o.i.).
 */
class TRAVERSAL_NETWORK{
 public:
     std::vector<int> sourceNodeIDs;
     std::vector< std::vector<CONN> > regConnections;
     std::vector< std::vector<CONN> > connectToSource;
     std::vector< std::vector<CONN> > connectToSink;
  DELTA_POS direction;
  DIJKSTRA_NETWORK *dnet;

  /** Construct a traversal network in the direction (dx,dy,dz) using the provided source ids, sink
   *  ids, and classified connections. */
  TRAVERSAL_NETWORK(int dx, int dy, int dz, std::vector<int> srcIDs, std::vector< std::vector<CONN> > regConns, 
		    std::vector< std::vector<CONN> > srcConns, std::vector< std::vector<CONN> > sinkConns, DIJKSTRA_NETWORK * dnetwork);

  /** Construct a TRAVERSAL_ NETWORK in the direction (dx,dy,dz) relative to the unit cell vectors.
   *  Identify all source nodes and classify the connections based on their component change in unit cell
   *  in the direction of interest. */
  TRAVERSAL_NETWORK(int dx, int dy, int dz, DIJKSTRA_NETWORK * dnetwork);

  /* Print information about the network to the provided output
  *  stream. Information includes the source node ids, and the connections
  *  of type source, regular and sink.*/
  void print(std::ostream &out);

  /* Returns a pair containing a bool and the maximum sized PATH a particle
   * can take to traverse this portion of the network. The bool is true iff the PATH is 
   * periodically viable (node #N to node #N).*/
  std::pair<bool,PATH> findMaxFreeSphere(std::map<int,int> *idAliases, std::set<int> *sourceNodes);
  
  /* Returns a pair containing the maximum radius of a free sphere that can
   *  traverse the portion of the network and the corresponding optimal PATH. No 
   *  check is made as to whether it is periodically viable (node #N to node #N).*/
  PATH findMaxFreeSphere();
};


/* Returns the maximum diameter of a sphere that can be inside of the
   VORONOI_NETWORK .*/ 
double findMaxIncludedSphere(VORONOI_NETWORK *vornet);


bool betterPath(std::pair<int, std::pair<DELTA_POS, double> > p1, std::pair<int, std::pair<DELTA_POS, double> > p2);

double calculateNodalFreeSphere(DIJKSTRA_NETWORK *dnet, DELTA_POS filter);

bool compareNodeData(std::pair<int, double> d1, std::pair<int, double> d2);



/* Calculates the fraction of the node's atomic radii that overlap
   and returns the result in the form of a 2D vector.
 */
std::vector< std::vector<double> > calculateNodeOverlap(ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet);
  

  /* Divide the set of node ids into subsets representing individual cages.
   * Returns the results in the form of a vector of sets,
   * where each set represents the node ids that form a single cage. MIN_FRACTION
   * represents the minimum amount of sphere-sphere overlap required for two
   * nodes to belong to the same cluster
   */
std::vector< std::set<int> > clusterElements(std::set<int> cageNodes, std::vector< std::vector<double> > overlapFractions, double MIN_FRACTION);

class CAGE {
    std::map<int,int> nodeIndices;
  std::vector<int> nodeIDs;
  std::vector<DELTA_POS> nodeOffsets;
  double cx, cy, cz, radius;

  /* Compute the center and radius of the cage formed by the provided set of node ids in
   * the Voronoi network. Store the resulting center and radius as well as the node ids
   * and displacements.
   */
  void reconstructCage(std::set<int> ids, ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet);

public:
  CAGE(std::set<int> ids, ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet);

  DELTA_POS getNodeOffset(int nodeID);

  double getRadius();
  Point  getCenter();
  std::vector<int> getNodeIDs();
  void writeToVMD(int cageIndex, std::fstream &output);
}; 

 // Defined in area_and_volume.cc
 double calcDensity(ATOM_NETWORK *);

 
 /* Identify the cages present in Voronoi network for the given probe size. If the
  * visualization option is specified, commands necessary to visualize the cages in 
  * ZeoVis is written to the provided file stream. Otherwise, the cages are classifed 
  * according to radius and the relevant cage sizes and volume is output using the file stream.
  */
void identifyCages(ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet, DIJKSTRA_NETWORK *dnet, double probeRadius, bool visualize, std::fstream &output, std::vector<CAGE> &cages); 



void simplifyCageNetwork(ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet, DIJKSTRA_NETWORK *dnet, double probe_radius);
 

#endif
