
//#include "network.h"
#include <cstdlib>

#include "voronoicell.h"
#include "graphstorage.h"

using namespace std;

/* This file contains two definitions of Voronoi cells: one that stores
 * information about the faces, edges and nodes in each cell and one that only
 * stores information about the node coordinates. The first definition, VOR_CELL, is used
 * when visualization is important. The other, BASIC_VCELL, is used in surface area/volume
 * calculations to decrease overhead time. 
 *
 * The file also contains all of the functions required to visualize an ATOM_NETWORK
 * and VORONOI_NETWORK using the ZeoVis tool.
 */



/* Store the provided vertices and their ids.   */
VOR_FACE::VOR_FACE(vector<Point> vertices, ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet){
  orderedVertices = vertices;
  for(unsigned int i = 0; i < orderedVertices.size(); i++){
    nodeIDs.push_back(getNodeID(orderedVertices.at(i), atmnet, vornet));
  }
}

/* Store the provided vertices and their ids.  */
VOR_FACE::VOR_FACE(vector<Point> vertices, vector<int> vertexIDs) {
  orderedVertices = vertices;
  nodeIDs = vertexIDs;
}

/* Returns a vector of pairs of integers and Points, where each entry represents
 *  the pair (node id, node coordinates).*/
vector< pair<int,Point> > VOR_FACE::getNodes(){
  vector< pair<int,Point> > results = vector< pair<int,Point> > ();
  for(unsigned int i = 0; i < orderedVertices.size(); i++){
    results.push_back(pair<int,Point> (nodeIDs[i], orderedVertices[i]));
  }
  return results;
}

/* Returns a vector of pairs of Points representing the start and end coordinates
 * of each edge that outlines the face.   */
vector< pair<Point, Point> > VOR_FACE::getEdgeCoords() {
  vector< pair<Point, Point> > results = vector< pair<Point,Point> > ();
  for(unsigned int i = 0; i < orderedVertices.size() - 1; i++){
    Point curPt = orderedVertices[i];
    Point nextPt = orderedVertices[i+1];
    results.push_back(pair<Point,Point> (curPt,nextPt));
  }
  results.push_back(pair<Point,Point> (orderedVertices[orderedVertices.size() - 1],orderedVertices[0]));
  return results;
}

/* Using the provided ouput stream, write the geometric objects that would fill the Voronoi 
 * face when drawn.   */
void VOR_FACE::writeVMDFilled(fstream &output){
  Point p1 = orderedVertices[0];
  unsigned int n2 = 1;
  unsigned int n3 = 2;
  while(n3 < orderedVertices.size()){
    Point p2 = orderedVertices.at(n2);
    Point p3 = orderedVertices.at(n3);
    output << "{triangle {" 
	   << p1[0] << " " << p1[1] << " " << p1[2] << "} {" 
	   << p2[0] << " " << p2[1] << " " << p2[2] << "} {" 
	   << p3[0] << " " << p3[1] << " " << p3[2] << "} }" << "\n";
    n2++;
    n3++;
  }
}



/*Constructs a VOR_CELL that does not initially have any faces or vertices.*/
VOR_CELL::VOR_CELL(){
  faces = vector<VOR_FACE> ();
  numVertices = 0;
  vertexIDs = map<Point, int, bool(*)(Point,Point)> (pointIsLess);
  idMappings = map<int,int> ();
  reverseIDMappings = map<int, vector<int> > ();
  vertexCoords = map<int,Point> ();
  edgeConnections = vector< set<int> > ();
}

/* Add the provided coordinate and its corresponding node id if
 *  no such vertex has been previously added.*/
void VOR_CELL::addNode(int nodeID, Point coord){
  if(vertexIDs.find(coord) == vertexIDs.end()){
    idMappings.insert(pair<int,int> (numVertices, nodeID));
    
    map<int, vector<int> >::iterator revMap = reverseIDMappings.find(nodeID);      
    if(revMap == reverseIDMappings.end()){
      vector<int> newList = vector<int> ();
      newList.push_back(numVertices);
      reverseIDMappings.insert(pair<int, vector<int> > (nodeID, newList));
    }
    else {
      revMap->second.push_back(numVertices);
    }
    vertexIDs.insert(pair<Point,int> (coord, numVertices));
    vertexCoords.insert(pair<int,Point> (numVertices, coord));
    edgeConnections.push_back(set<int>());
    numVertices++;
  }
}
  
/* Add the edge that spans the two points if it has not yet been added.*/
void VOR_CELL::addEdge(Point from, Point to){
  map<Point,int>::iterator iter1 = vertexIDs.find(from);
  map<Point,int>::iterator iter2 = vertexIDs.find(to);
  if((iter1 == vertexIDs.end()) || (iter2 == vertexIDs.end())){
    cerr << "Unable to add edge because nodes have not been added." << "\n"
	 << "Point 1: (" << from[0] <<  ", " << from[1] << ", " << from[2] << ")" << "\n"
	 << "Point 2: (" << to[0]   <<  ", " << to[1]   << ", " << to[2]   << ")" << "\n"
	 << "Exiting..." << "\n";
    exit(1);
  }
  
  if(edgeConnections[iter2->second].find(iter1->second) == edgeConnections[iter2->second].end())
    edgeConnections[iter1->second].insert(iter2->second);
}

/* Add the face to the VOR_CELL. Adds all edges and vertices that have not yet been added 
 *  to the cell.*/
void VOR_CELL::addFace(VOR_FACE face) {
  faces.push_back(face);
  vector< pair<int,Point> > nodes = face.getNodes();
  for(unsigned int i = 0; i < nodes.size(); i++){
    pair<int, Point> node = nodes[i];
    addNode(node.first, node.second);
  }
  vector< pair<Point, Point> > edges= face.getEdgeCoords();
  for(unsigned int i = 0; i < edges.size(); i++){
    pair<Point,Point> coord = edges[i];
    addEdge(coord.first, coord.second);
  }
}

/* Returns a vector containing all of the coordinates in the VOR_CELL that
 *  whose node id is as provided.*/
vector<Point> VOR_CELL::getNodeCoords(int nodeID){
  map<int, vector<int> >::iterator iter = reverseIDMappings.find(nodeID);
  if(iter == reverseIDMappings.end()){
    cerr << "Error: Node #" << nodeID << " isn't in this Voronoi cell." << "\n";
    cerr << "Cell contains these nodes: "; 
    map<int, vector<int> >::iterator nIter = reverseIDMappings.begin();
    while(nIter != reverseIDMappings.end()){
      cerr << nIter->first << " ";
      nIter++;
    }
    cerr << "\n";
    cerr << "Exiting..." << "\n";
    exit(1);
  }
  vector<int> vertexIDs = iter->second;
  vector<Point> coords = vector<Point> ();
  for(unsigned int i = 0; i < vertexIDs.size(); i++){
    coords.push_back(vertexCoords.find(vertexIDs[i])->second);
  }
  return coords;
}

/* Using the underlying list of faces, write the set of commands necessary
 *  to fill the VOR_CELL's exterior to the provided output stream. Labels 
 *  the corresponding commands as faces(n).*/
void VOR_CELL::writeVMDFilled(fstream &output, int n){
  output << "set faces(" << n << ") {" << "\n"
	 << "{color $faceColors(" << n << ") }" << "\n";
  for(unsigned int i = 0; i < faces.size(); i++){
    faces[i].writeVMDFilled(output);
  }
  output << "}" << "\n";
}

/* Using the set of nodes and edges, write the set of commands necessary to 
 *  draw the outline of the VOR_CELL to the provided output stream. Labels 
 *  the corresponding commands as vorcells(n)*/
void VOR_CELL::writeVMDOutlined(fstream &output, int n){
  output << "set vorcells(" << n <<") {" << "\n";
  
  // Iterate over all nodes in VOR_CELL
  for(int i = 0; i < numVertices; i++){
    Point curPoint = vertexCoords.find(i)->second;
    int nodeID = idMappings.find(i)->second;
    output << "{color $nodeColors(" << nodeID << ") }" << "\n"; 
    output << "{sphere {" << curPoint[0] << " " << curPoint[1] << " " << curPoint[2] << "} radius $nodeRadii(" << nodeID << ") resolution $sphere_resolution}" << "\n";
  }
  
  // Iterate over all edges in VOR_CELL
  output << "{color $vorcellColors(" << n << ") }" << "\n";
  for(int i = 0; i < numVertices; i++){
    Point p1 = vertexCoords[i];
    set<int>::iterator vertexIter = edgeConnections[i].begin();
    while(vertexIter != edgeConnections[i].end()){
      Point p2 = vertexCoords[*vertexIter];
	output << "{line {" 
	       << p1[0] << " " << p1[1] << " " << p1[2] << "} {" 
	       << p2[0] << " " << p2[1] << " " << p2[2] << "} width 1}" << "\n";
	vertexIter++;
    }
  }
  output << "}" << "\n";
}

/*  Using the provided list of VOR_FACEs, reconstruct each VOR_CELL
 *  and store it using the provided vector pointer. */
void VOR_CELL::getVoronoiCells(vector<VOR_CELL> *cells, vector< vector<VOR_FACE> > faceInfo,
  ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet){
  cells->clear();
  for(unsigned int i = 0; i < faceInfo.size(); i++){
    VOR_CELL newCell;
    for(unsigned int j = 0; j < faceInfo[i].size(); j++)
      newCell.addFace(faceInfo[i][j]);
    cells->push_back(newCell);
  }
}






BASIC_VCELL::BASIC_VCELL(){
  nodeCoords = vector<Point> ();
  nodeIDs = vector<int> ();
}

BASIC_VCELL::BASIC_VCELL (vector<Point> myNodeCoords, vector<int> myNodeIDs){
  nodeCoords = myNodeCoords;
  nodeIDs = myNodeIDs;
}

int BASIC_VCELL::getNumNodes(){
  return nodeCoords.size();
}

Point BASIC_VCELL::getNodeCoord(int index){
  return nodeCoords[index];
}
  
int BASIC_VCELL::getNodeID(int index){
  return nodeIDs[index];
}

/* Removes the nodes from the cell which would overlap with a sphere center around
   the provided atom with a radius of r_probe + r_atom*/
void BASIC_VCELL::removeOverlappedNodes(int atomID, ATOM_NETWORK *atmnet, double r_probe){
  vector<int> newIDs;
  vector<Point> newCoords;
  ATOM center = atmnet->atoms[atomID];
  for(unsigned int i = 0; i < nodeCoords.size(); i++){
    Point node  = nodeCoords[i];
    if(calcEuclideanDistance(node[0], node[1], node[2], center.x, center.y, center.z) >= 
		(center.radius + r_probe)){
      newIDs.push_back(nodeIDs[i]);
      newCoords.push_back(nodeCoords[i]);
    }
  }
  nodeIDs = newIDs; 
  nodeCoords = newCoords;
}

void BASIC_VCELL::writeToVMD(fstream &output, int n){
  output << "set nodecells(" << n <<") {" << "\n";
  
  // Iterate over all nodes in the cell
  for(unsigned int i = 0; i < nodeCoords.size(); i++){
    Point curPoint = nodeCoords[i];
    int nodeID = nodeIDs[i];
    output << "{color $nodeColors(" << nodeID << ") }" << "\n"; 
    output << "{sphere {" << curPoint[0] << " " << curPoint[1] << " " << curPoint[2]
           << "} radius $nodeRadii(" << nodeID << ") resolution $sphere_resolution}" << "\n";
  }
  output << "}" << "\n";
}




/* Writes the commands to the provided output stream necessary to approriately
 * display the unit cell in the ZeoVis visualization tool. Labels the unitcell
 * as unitcells(0).
 */
void writeVMDUC(fstream &output, ATOM_NETWORK *atmnet){
  XYZ v_a = atmnet->v_a;
  XYZ v_b = atmnet->v_b;
  XYZ v_c = atmnet->v_c;
  
  output << "set unitcells(0) {" << "\n"
	 << "{color $unitcellColors(0)}" << "\n";
  DELTA_POS directions [3] = {DELTA_POS(1,0,0), DELTA_POS(0,1,0), DELTA_POS(0,0,1)};
  DELTA_POS limits [3] = {DELTA_POS(0,1,1), DELTA_POS(1,0,1), DELTA_POS(1,1,0)};
  for(unsigned int i = 0; i < 3; i++){
    DELTA_POS direction = directions[i];
    DELTA_POS limit = limits[i];
    for(int a = 0; a < 2; a++){
      for(int b = 0; b < 2; b++){
	for(int c = 0; c < 2; c++){
	  if((limit.x < a) || (limit.y < b) || (limit.z < c))
	    continue;
	  
	  // Calculate starting coordinate
	  double x1 = v_a.x *a + v_b.x*b + v_c.x*c;
	  double y1 = v_a.y *a + v_b.y*b + v_c.y*c;
	  double z1 = v_a.z *a + v_b.z*b + v_c.z*c;
	  
	  // Calculate ending coordinate
	  double x2 = x1 + v_a.x*direction.x + v_b.x*direction.y + v_c.x*direction.z;
	  double y2 = y1 + v_a.y*direction.x + v_b.y*direction.y + v_c.y*direction.z;
	  double z2 = z1 + v_a.z*direction.x + v_b.z*direction.y + v_c.z*direction.z;

	  output << "{line "
		 << "{" << x1 << " " << y1 << " " << z1 << "} "
		 << "{" << x2 << " " << y2 << " " << z2 << "} }" << "\n";
	}
      }
    }
  }
  output << "}" << "\n";
}

/* Writes the commands to the provided output stream necessary to display 
 * the entire Voronoi network. Labels the Voronoi network as vornets(0). 
 */
void writeVornet(fstream &output, ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet){
  output << "set vornets(0) {" << "\n";

  // Iterate over all Voronoi nodes
  for(unsigned int i = 0; i < vornet->nodes.size(); i++){
    VOR_NODE curNode = vornet->nodes.at(i);
    output << "{color $nodeColors("<< i << ") }" << "\n"
	       << "{sphere {" << curNode.x << " " << curNode.y << " " << curNode.z 
	       << "} radius $nodeRadii(" << i << ") resolution $sphere_resolution}" << "\n";
  }

  // Iterate over all Voronoi edges
  output << "{color $vornetColors(0)}" << "\n";
  for(unsigned int i = 0; i < vornet->edges.size(); i++){
    VOR_EDGE curEdge = vornet->edges.at(i);
    VOR_NODE startNode = vornet->nodes.at(curEdge.from);
    Point start = Point(startNode.x, startNode.y, startNode.z);
    VOR_NODE endNode = vornet->nodes.at(curEdge.to);
    Point end = Point(endNode.x, endNode.y, endNode.z);
    atmnet->translatePoint(&end, curEdge.delta_uc_x, curEdge.delta_uc_y, curEdge.delta_uc_z);
    output << "{line {" 
	   << start[0] << " " << start[1] << " " << start[2] << "} " 
	   << "{" << end[0] << " " << end[1] << " " << end[2] << "}" 
	   << "}" << "\n";
  }
  output << "}" << "\n";
}

/* Writes the commands to the provided output stream necessary to display all atoms and
 * Voronoi nodes in ZeoVis. Labels each atom as atoms(i) and each node as nodes(j).
 */
void writeVMDAtomsAndNodes(fstream &output, ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet){
  // Iterate over all atoms. Set the radius of the atom to its actual value
  for(unsigned int i = 0; i < atmnet->atoms.size(); i++){
    ATOM curAtom = atmnet->atoms.at(i);
    output << "set atoms(" << i << ") {" << "\n"
	       << "{color $atomColors(" << i << ") }" << "\n"
	       << "{sphere {" << curAtom.x << " " << curAtom.y << " " << curAtom.z 
	       << "} radius $atomRadii(" << i << ") resolution $sphere_resolution}" 
	       << "\n" << "}" << "\n";
    output << "set atomRadii(" << i << ") " << curAtom.radius << "\n";
  }

  // Iterate over all nodes. Set the radius of the node to that of the maximum static sphere
  for(unsigned int i = 0; i < vornet->nodes.size(); i++){
    VOR_NODE curNode = vornet->nodes.at(i);
    output << "set nodes(" << i << ") {" << "\n"
	   << "{color $nodeColors("<< i << ") }" << "\n"
	   << "{sphere {" << curNode.x << " " << curNode.y << " " << curNode.z 
	   << "} radius $nodeRadii(" << i <<  ") resolution $sphere_resolution}" << "\n"
	   << "}" << "\n";
    output << "set nodeRadii(" << i << ") " << curNode.rad_stat_sphere << "\n";
  }    
}

/* Writes the commands to the provided output stream for ZeoVis
 * that establishes the number of each component as well as 
 * unitcell vector information.
 */
void writeVMDEnvVars(fstream &output,  ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet){
  // Write information about the number of each component
  output << "set num_vorcells " << atmnet->numAtoms << "\n";
  output << "set num_faces " << atmnet->numAtoms << "\n";
  output << "set num_vornets 1" << "\n";
  output << "set num_nodes " << vornet->nodes.size() << "\n";
  output << "set num_atoms " << atmnet->numAtoms << "\n";
  output << "set num_unitcells 1" << "\n";
  output << "set num_channels 0" << "\n";
  
  // Write unitcell vector information
  output << "set uc_a_vector {" << atmnet->v_a.x << " " << atmnet->v_a.y 
         << " " << atmnet->v_a.z << "}" << "\n";
  output << "set uc_b_vector {" << atmnet->v_b.x << " " << atmnet->v_b.y
         << " " << atmnet->v_b.z << "}" << "\n";
  output << "set uc_c_vector {" << atmnet->v_c.x << " " << atmnet->v_c.y
         << " " << atmnet->v_c.z << "}" << "\n";
  output << "set sphere_resolution 100" << "\n";
}



/* Writes all of the information necessary to properly visualize the Voronoi
 * and atom network in ZeoVis to the file referred to by filename. Information 
 * includes atoms, nodes, the unitcell, the voronoi network, 
 * voronoi cells (outlined and filled) and environment variables.
 */
void writeZeoVisFile(char *filename, vector<VOR_CELL> *cells,
  ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet)
{
  fstream output;
  output.open(filename, fstream::out);
  if(!output.is_open()){
    cout << "Error: Failed to open output file for ZeoVis settings" << filename;
    cout << "Exiting ..." << "\n";
    exit(1);
  }
  else{
    cout << "Writing ZeoVis information to " << filename << "\n";
    
    writeVMDEnvVars(output, atmnet, vornet);
    writeVMDAtomsAndNodes(output, atmnet, vornet);
    writeVornet(output, atmnet, vornet);
    writeVMDUC(output, atmnet);

    for(unsigned int i = 0; i < cells->size(); i++){
      cells->at(i).writeVMDOutlined(output, i);
      cells->at(i).writeVMDFilled(output, i);
    }
    output << "set num_faces " << cells->size() << "\n"
	   << "set num_channels " << 0 << "\n"
	   << "set num_features " << 0 << "\n"
	   << "set num_segments " << 0 << "\n"
	   << "set num_cages "    << 0 << "\n";
  }
  output.close();
}

/* Identical to writeZeoVisFile except that information about the basic voronoi 
 * cells is also outputted. */
void writeSpecialZeoVisFile(char *filename, vector<VOR_CELL> *cells, 
  ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet, vector<BASIC_VCELL> &vcells)
{
  fstream output;
  output.open(filename, fstream::out);
  if(!output.is_open()){
    cout << "Error: Failed to open output file for ZeoVis settings" << filename;
    cout << "Exiting ..." << "\n";
    exit(1);
  }
  else{
    cout << "Writing ZeoVis information to " << filename << "\n";
    
    writeVMDEnvVars(output, atmnet, vornet);
    writeVMDAtomsAndNodes(output, atmnet, vornet);
    writeVornet(output, atmnet, vornet);
    writeVMDUC(output, atmnet);

    for(unsigned int i = 0; i < cells->size(); i++){
      cells->at(i).writeVMDOutlined(output, i);
      cells->at(i).writeVMDFilled(output, i);
    }
    output << "set num_faces " << cells->size() << "\n";
    output << "set num_channels " << 0 << "\n";
    
    for(unsigned int i = 0; i < vcells.size(); i++){
      vcells[i].writeToVMD(output, i);
    }

  }
  output.close();
}

