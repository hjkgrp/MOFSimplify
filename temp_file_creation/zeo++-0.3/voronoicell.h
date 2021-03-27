/* This file contains two defintions of Voronoi cells: one that stores
 * information about the faces, edges and nodes in each cell and one that only
 * stores information about the node coordinates. The first definition, VOR_CELL, is used
 * when visualization is important. The other, BASIC_VCELL, is used in surface area/volume
 * calculations to decrease overhead time. 
 *
 * The file also contains all of the functions required to visualize an ATOM_NETWORK
 * and VORONOI_NETWORK using the ZeoVis tool.
 */


/** Class used to store the node ids and coordinates of the vertices 
 *  that comprise a face of a Voronoi cell. The class orders the vertices
 *  in clockwise or counter-clockwise fashion and keeps track of the edges
 *  that outline it.
 */

#ifndef VORONOICELL_H
#define VORONOICELL_H

#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <utility>

#include "zeo_consts.h"
#include "geometry.h"
#include "networkstorage.h"

const unsigned int randSeed = 994879221;
//const double pi = 3.1415926535;
//const double AVOGRADOS_NUMBER = 6.0221415e23;

class VOR_FACE{
public:
    std::vector<Point> orderedVertices;
    std::vector<int> nodeIDs;

  /* Store the provided vertices and their ids.   */
  VOR_FACE(std::vector<Point> vertices, ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet);

  /* Store the provided vertices and their ids.  */
  VOR_FACE(std::vector<Point> vertices, std::vector<int> vertexIDs);
  
  /* Returns a vector of pairs of integers and Points, where each entry represents
  *  the pair (node id, node coordinates).*/
  std::vector< std::pair<int,Point> > getNodes();

  /* Returns a vector of pairs of Points representing the start and end coordinates
   * of each edge that outlines the face.   */
  std::vector< std::pair<Point, Point> > getEdgeCoords();

  /* Using the provided ouput stream, write the geometric objects that would fill the Voronoi 
   * face when drawn.   */
  void writeVMDFilled(std::fstream &output);
};

/* Class used to reconstruct Voronoi cells from their constituent faces.  Used to write configuration files for 
 * ZeoVis visualization tool. 
 */
class VOR_CELL {
public:
    std::vector<VOR_FACE> faces;                             // List of faces in cell
  int numVertices;                                    // Number of vertices
  std::map<Point, int, bool(*)(Point,Point)> vertexIDs;    // Vertex number associated with coordinate
  std::map<int,int> idMappings;                            // (Vertex ID, Node ID) pairs
  std::map<int,std::vector<int> > reverseIDMappings;            // (Node ID, List of Vertex IDs) pairs
  std::map<int, Point> vertexCoords;                       // Coordinates associated with each vertex ID
  std::vector< std::set<int> > edgeConnections;                 // List of vertex id's connected to each vertex

  /*Constructs a VOR_CELL that does not initially have any faces or vertices.*/
  VOR_CELL();

  /* Add the provided coordinate and its corresponding node id if
  *  no such vertex has been previously added.*/
  void addNode(int nodeID, Point coord);
  
  /* Add the edge that spans the two points if it has not yet been added.*/
  void addEdge(Point from, Point to);

  /* Add the face to the VOR_CELL. Adds all edges and vertices that have not yet been added 
  *  to the cell.*/
  void addFace(VOR_FACE face);

  /* Returns a vector containing all of the coordinates in the VOR_CELL that
  *  whose node id is as provided.*/
  std::vector<Point> getNodeCoords(int nodeID);

  /* Using the underlying list of faces, write the set of commands necessary
  *  to fill the VOR_CELL's exterior to the provided output stream. Labels 
  *  the corresponding commands as faces(n).*/
  void writeVMDFilled(std::fstream &output, int n);
 
  /* Using the set of nodes and edges, write the set of commands necessary to 
  *  draw the outline of the VOR_CELL to the provided output stream. Labels 
  *  the corresponding commands as vorcells(n)*/
  void writeVMDOutlined(std::fstream &output, int n);

  /*  Using the provided list of VOR_FACEs, reconstruct each VOR_CELL
   *  and store it using the provided vector pointer. */
  static void getVoronoiCells(std::vector<VOR_CELL> *cells, std::vector< std::vector<VOR_FACE> > faceInfo,  ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet);
};


/* Simple class used to store the positions and ids 
** of the nodes that comprise a Voronoi cell.
*/
class BASIC_VCELL{
    std::vector<Point> nodeCoords;
    std::vector<int> nodeIDs;
public:
  BASIC_VCELL();
  
  BASIC_VCELL (std::vector<Point> myNodeCoords, std::vector<int> myNodeIDs);
  int getNumNodes();

  Point getNodeCoord(int index);
  
  int getNodeID(int index);

  /* Removes the nodes from the cell which would overlap with a sphere center around
   the provided atom with a radius of r_probe + r_atom*/
  void removeOverlappedNodes(int atomID, ATOM_NETWORK *atmnet, double r_probe);

  void writeToVMD(std::fstream &output, int n);
};

/* Writes the commands to the provided output stream necessary to approriately
 * display the unit cell in the ZeoVis visualization tool. Labels the unitcell
 * as unitcells(0).
 */
void writeVMDUC(std::fstream &output, ATOM_NETWORK *atmnet);

/* Writes the commands to the provided output stream necessary to display 
 * the entire Voronoi network. Labels the Voronoi network as vornets(0). 
 */
void writeVornet(std::fstream &output, ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet);


/* Writes the commands to the provided output stream necessary to display all atoms and
 * Voronoi nodes in ZeoVis. Labels each atom as atoms(i) and each node as nodes(j).
 */
void writeVMDAtomsAndNodes(std::fstream &output, ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet);

/* Writes the commands to the provided output stream for ZeoVis
 * that establishes the number of each component as well as 
 * unitcell vector information.
 */
void writeVMDEnvVars(std::fstream &output,  ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet);

/* Writes all of the information necessary to properly visualize the Voronoi and atom network in ZeoVis to the 
 * file referred to by filename. Information includes atoms, nodes, the unitcell, the voronoi network, 
 *  voronoi cells (outlined and filled) and environment variables.
 */
void writeZeoVisFile(char *filename, std::vector<VOR_CELL> *cells, ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet);

void writeSpecialZeoVisFile(char *filename, std::vector<VOR_CELL> *cells, ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet, std::vector<BASIC_VCELL> &vcells);

#endif
