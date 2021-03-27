//#include "network.h"
#include <cstdlib>
#include <cmath>
#include <float.h>
//#include <algorithm>

#include "networkanalysis.h"
#include "heap.h"

using namespace std;

/* Create an empty PATH with invalid length and radii. */
PATH::PATH(){
  max_radius = max_inc_radius = length = -1;
  visitedNodes = vector<DIJKSTRA_NODE> ();
  toSink = false;
  visitedSourceIDs = set<int> ();
}

/* Create a PATH with the provided current node, radius and length.*/
PATH::PATH(DIJKSTRA_NODE cNode, double maxRadius, double maxIncRadius, double len){
  currentNode = cNode;
  max_radius = maxRadius;
  max_inc_radius = maxIncRadius;
  length = len;
  visitedNodes = vector<DIJKSTRA_NODE> ();
  toSink = false;
  visitedSourceIDs = set<int> ();
}

/* Prints the information about the edges in the DIJKSTRA_NETWORK that were
   used when following the PATH. */
void PATH::printPathEdges (vector<int> nodeList, DIJKSTRA_NETWORK *dnet){
  if(nodeList.size() < 2)
    return;
  
  unsigned int index = 0;
  int from, to;
  
  while(index < nodeList.size()-1){
    from = nodeList.at(index);
    to = nodeList.at(index+1);
    
    DIJKSTRA_NODE node = dnet->nodes.at(from);
    for(unsigned int i = 0; i < node.connections.size(); i++){
      CONN conn = node.connections.at(i);
      if(conn.to == to){
	conn.print();
      }
    }
    index++;
  }
}

/* Print information about the PATH to the provided output stream.*/
void PATH::print(ostream &output){
  output << "  "<< "Node lineage: ";
  for(unsigned int i=0; i<visitedNodes.size()-1; i++){
    output << "   " << visitedNodes.at(i).id << " -> ";
  }
  output << "   " << visitedNodes.back().id << "\n"
	   << "  Current node: "   << currentNode.id << "\n"
	 << "  Maximum radius: " << max_radius     << "\n";
}


/* Returns true iff the first PATH has a lower maximu radius than the second PATH.*/
bool hasLowerMaxR(PATH n1, PATH n2){
  return n1.max_radius < n2.max_radius;
}



/** Construct a traversal network in the direction (dx,dy,dz) using the provided source ids, sink
 *  ids, and classified connections. */
TRAVERSAL_NETWORK::TRAVERSAL_NETWORK(int dx, int dy, int dz, vector<int> srcIDs, vector< vector<CONN> > regConns, 
				     vector< vector<CONN> > srcConns, vector< vector<CONN> > sinkConns, DIJKSTRA_NETWORK * dnetwork){
  sourceNodeIDs = srcIDs;
  regConnections = regConns;
  connectToSource = srcConns;
  connectToSink = sinkConns;
  direction = DELTA_POS(dx,dy,dz);
  dnet = dnetwork;
}

/** Construct a TRAVERSAL_ NETWORK in the direction (dx,dy,dz) relative to the unit cell vectors.
 *  Identify all source nodes and classify the connections based on their component change in unit cell
 *  in the direction of interest. */
TRAVERSAL_NETWORK::TRAVERSAL_NETWORK(int dx, int dy, int dz, DIJKSTRA_NETWORK * dnetwork){
  direction = DELTA_POS(dx,dy,dz);
  sourceNodeIDs = vector<int> ();
  regConnections = vector< vector<CONN> > ();
  connectToSource = vector< vector<CONN> > ();
  connectToSink = vector< vector<CONN> > ();
  dnet = dnetwork;
    
  map<int,bool> sourceIDs;
  
  // Iterate over all nodes
  for(unsigned int i = 0; i < dnet->nodes.size(); i++){
    DIJKSTRA_NODE currentNode = dnet->nodes.at(i);
    vector<CONN> regConns = vector <CONN> ();
    vector<CONN> srcConns = vector <CONN> ();
    vector<CONN> drainConns = vector<CONN> ();
    
    // Iterate over all connections
    for(unsigned int j = 0; j < currentNode.connections.size(); j++){
      CONN currentConn = currentNode.connections.at(j);
      DELTA_POS component = currentConn.deltaPos*direction;
      
      int directionComp;
      int otherComp;
      
      if(!component.isZero()){
	if(direction.x != 0){
	  directionComp = direction.x;
	  otherComp = component.x;
	}
	else if (direction.y != 0){
	  directionComp = direction.y;
	  otherComp = component.y;
	}
	else if (direction.z != 0){
	  directionComp = direction.z;
	  otherComp = component.z;
	}
	else{
	  cerr << "Invalid argument reached when building TRAVERSAL_NETWORK. Please contact the source code provider with your input" 
	       << "\n" << "Exiting..." << "\n";
	  exit(1);
	}
	
	if(directionComp == otherComp){
	  //Connection is a drain connection and receiver is a drain node
	  drainConns.push_back(currentConn);
	}
	else if (directionComp == -1*otherComp){
	  //Connection is source connection and receiver is a source node
	  srcConns.push_back(currentConn);
	  sourceIDs.insert(pair<int,bool> (currentConn.to,true));
	}
	else{
	  // There's a problem with the logic
	  cerr << "Invalid argument reached when building TRAVERSAL_NETWORK. Please contact the source code provider with your input" 
	       << "\n" << "Exiting..." << "\n";
	  exit(1);
	}
      }
      else{
	// Connection is a drain connection
	regConns.push_back(currentConn);
      }
    }      
    
    // Store connections according to type
      regConnections.push_back(regConns);
      connectToSource.push_back(srcConns);
      connectToSink.push_back(drainConns);
  }
  
  //Copy from map to id list to prevent duplicate id's
  map<int,bool>::iterator iter = sourceIDs.begin();
  while(iter != sourceIDs.end()){
    sourceNodeIDs.push_back(iter->first);
    iter++;
  }
}

/* Print information about the network to the provided output
   stream. Information includes the source node ids, and the connections
   of type source, regular and sink.*/
void TRAVERSAL_NETWORK::print(ostream &out){
  out << "Source nodes ids:  ";
  for(unsigned int i = 0; i < sourceNodeIDs.size(); i++)
    out << sourceNodeIDs.at(i) << "  ";
  out << "\n";
  out << "Regular connections:" << "\n";
  for(unsigned int i = 0; i <regConnections.size();i++){
    vector<CONN> nodeConns = regConnections.at(i);
    if(nodeConns.size() == 0)
      continue;
    out << "From #"<<i << "   To: ";
    for(unsigned int j = 0; j< nodeConns.size(); j++){
      out << nodeConns.at(j).to << "  ";
    }
    out << "\n";
  }
  out << "Connections to source node:" << "\n";
  for(unsigned int i = 0; i <connectToSource.size();i++){
    vector<CONN> nodeConns = connectToSource.at(i);
    if(nodeConns.size() == 0)
      continue;
    out << "From #"<< i << "   To:";
    for(unsigned int j = 0; j< nodeConns.size(); j++){
      cout << nodeConns.at(j).to << "  ";
    }
    cout << "\n";
  }
    cout << "Connections to sink node:" << "\n";
    for(unsigned int i = 0; i <connectToSink.size();i++){
      vector<CONN> nodeConns = connectToSink.at(i);
      if(nodeConns.size() == 0)
	continue;
      cout << "From #"<<i << "   To:";
      for(unsigned int j = 0; j< nodeConns.size(); j++){
	cout << nodeConns.at(j).to << "  ";
      }
      cout << "\n";
    }
    cout << "\n" << "\n";
}

/* Returns a pair containing a bool and the maximum sized PATH a particle
 * can take to traverse this portion of the network. The bool is true iff the PATH is 
 * periodically viable (node #N to node #N).*/
pair<bool,PATH> TRAVERSAL_NETWORK::findMaxFreeSphere(map<int,int> *idAliases, set<int> *sourceNodes){
  PATH bestPath; bestPath.max_radius = -1;
  vector<bool> haveVisited = vector<bool> (dnet->nodes.size(), false);
  HEAP<PATH> heap (hasLowerMaxR);
  bool NtoN = false;
  
  // Initialize the stack with all nodes connected to the source nodes
  for(unsigned int i = 0; i < sourceNodeIDs.size(); i++){
    DIJKSTRA_NODE sourceNode = dnet->nodes.at(sourceNodeIDs.at(i));
    vector<CONN> viableConns = connectToSink.at(sourceNode.id);
    for(unsigned int j = 0; j <  viableConns.size(); j++){
      CONN nextConn = viableConns.at(j);
      DIJKSTRA_NODE nextNode = dnet->nodes.at(nextConn.to);
      PATH newPath = PATH(nextNode, nextConn.max_radius, max(sourceNode.max_radius, nextNode.max_radius), nextConn.length);
      newPath.visitedIDs.push_back(sourceNode.id);
      newPath.visitedSourceIDs.insert(idAliases->find(sourceNode.id)->second);
      heap.insert(newPath);
    }
  }
  
  while(heap.size() != 0){
    PATH best = heap.pop();
    
    //Done searching if max radius path leads to sink node
    if(best.toSink){
      best.visitedIDs.push_back(best.currentNode.id);
      bestPath = best;
      break;
    }
    
    // Mark nodes that have been visited or abandon current path if node
    // already visited
    if(haveVisited.at(best.currentNode.id))
      continue;
    else
      haveVisited.at(best.currentNode.id) = true;
    
    if(sourceNodes->find(best.currentNode.id) != sourceNodes->end()){
      int origID = idAliases->find(best.currentNode.id)->second;
      if(best.visitedSourceIDs.find(origID) != best.visitedSourceIDs.end()){
	best.visitedIDs.push_back(best.currentNode.id);
	bestPath = best;
	NtoN = true;
	break;
      }
      else{
	best.visitedSourceIDs.insert(origID);
      }	
    }
    best.visitedIDs.push_back(best.currentNode.id);
    
    // Add all nodes that are connected within the unit cell to the stack
    vector<CONN> regConns = regConnections.at(best.currentNode.id);
    for(unsigned int i = 0; i < regConns.size(); i++){
      CONN nextConn = regConns.at(i);
      if(!haveVisited.at(nextConn.to)){
	DIJKSTRA_NODE nextNode = dnet->nodes[nextConn.to];
	PATH nextPath = PATH(nextNode,min(best.max_radius,nextConn.max_radius), max(best.max_inc_radius, nextNode.max_radius), best.length + nextConn.length);
	nextPath.visitedSourceIDs = best.visitedSourceIDs;
	nextPath.visitedIDs = best.visitedIDs;
	heap.insert(nextPath);
      }
    }

    // Add all nodes that are connected in the sink to the stack
    vector<CONN> sinkConns = connectToSink.at(best.currentNode.id);
    for(unsigned int i = 0; i < sinkConns.size(); i++){
      CONN nextConn = sinkConns.at(i);
      DIJKSTRA_NODE nextNode = dnet->nodes[nextConn.to];
      PATH nextPath = PATH(nextNode, min(best.max_radius, nextConn.max_radius), max(best.max_inc_radius, nextNode.max_radius), best.length + nextConn.length);
      nextPath.toSink = true;
      nextPath.visitedIDs = best.visitedIDs;
      nextPath.visitedSourceIDs = best.visitedSourceIDs;
      heap.insert(nextPath);
    }
    
  }
  return pair<bool,PATH> (NtoN,bestPath);
}
 
/* Returns a pair containing the maximum radius of a free sphere that can
 *  traverse the portion of the network and the corresponding optimal PATH. No 
 *  check is made as to whether it is periodically viable (node #N to node #N).*/
PATH TRAVERSAL_NETWORK::findMaxFreeSphere(){
  double maxR = -1;
  PATH bestPath;
  vector<bool> haveVisited = vector<bool> (dnet->nodes.size(), false);
  HEAP<PATH> heap (hasLowerMaxR);
  
  // Initialize the stack with all nodes connected to the source nodes
  for(unsigned int i = 0; i < sourceNodeIDs.size(); i++){
    DIJKSTRA_NODE sourceNode = dnet->nodes.at(sourceNodeIDs.at(i));
    vector<CONN> viableConns = connectToSink.at(sourceNode.id);
    for(unsigned int j = 0; j <  viableConns.size(); j++){
      CONN nextConn = viableConns.at(j);
      DIJKSTRA_NODE nextNode = dnet->nodes.at(nextConn.to);
      PATH newPath = PATH(nextNode, nextConn.max_radius, max(sourceNode.max_radius, nextNode.max_radius), nextConn.length);
      newPath.visitedIDs.push_back(sourceNode.id);
      heap.insert(newPath);
    }
  }

  while(heap.size() != 0){
    PATH best = heap.pop();
    best.visitedIDs.push_back(best.currentNode.id);
    
    //Done searching if max radius path leads to sink node
    if(best.toSink){
      maxR = best.max_radius;
      bestPath = best;
      break;
    }
    
    // Mark nodes that have been visited or abandon current path if node
    // already visited
    if(haveVisited.at(best.currentNode.id))
      continue;
    else
      haveVisited.at(best.currentNode.id) = true;
    
    // Add all nodes that are connected within the unit cell to the stack
    vector<CONN> regConns = regConnections.at(best.currentNode.id);
    for(unsigned int i = 0; i < regConns.size(); i++){
      CONN nextConn = regConns.at(i);
      if(!haveVisited.at(nextConn.to)){
	DIJKSTRA_NODE nextNode = dnet->nodes.at(nextConn.to);
	PATH nextPath = PATH(nextNode,min(best.max_radius,nextConn.max_radius), max(best.max_inc_radius, nextNode.max_radius), best.length + nextConn.length);
	nextPath.visitedIDs = best.visitedIDs;
	heap.insert(nextPath);
      }
    }

    // Add all nodes that are connected in the sink to the stack
    vector<CONN> sinkConns = connectToSink.at(best.currentNode.id);
    for(unsigned int i = 0; i < sinkConns.size(); i++){
      CONN nextConn = sinkConns.at(i);
      DIJKSTRA_NODE nextNode = dnet->nodes.at(nextConn.to);
      PATH nextPath = PATH(nextNode, min(best.max_radius, nextConn.max_radius), max(best.max_inc_radius, nextNode.max_radius), best.length + nextConn.length);
      nextPath.toSink = true;
      nextPath.visitedIDs = best.visitedIDs;
      heap.insert(nextPath);
    }
  }
  return bestPath;
}


/* Returns the maximum diameter of a sphere that can be inside of the
   VORONOI_NETWORK .*/ 
double findMaxIncludedSphere(VORONOI_NETWORK *vornet){
  double max = 0;
  vector<VOR_NODE> ::iterator iter = vornet->nodes.begin();
  while(iter != vornet->nodes.end()){
      if(iter->rad_stat_sphere > max)
	max = iter->rad_stat_sphere;
      iter++;
  }
  return max;
}

bool betterPath(pair<int, pair<DELTA_POS, double> > p1, pair<int, pair<DELTA_POS, double> > p2){
  return p1.second.second <= p2.second.second;
}

double calculateNodalFreeSphere(DIJKSTRA_NETWORK *dnet, DELTA_POS filter){
  vector<bool> nodeProcessed = vector<bool>(dnet->nodes.size(), false);
  vector<double> bestRadii   = vector<double>(dnet->nodes.size(), -1);
  unsigned int nodeIndex = 0;

  double maxRadius = 0.0;

  // Iterate over all nodes
  while(nodeIndex < nodeProcessed.size()){
   
    // Place starting node on stack with (0,0,0) displacement
    DELTA_POS displacement = DELTA_POS(0,0,0);
    map<int, pair<DELTA_POS, double> > visitedNodes;
    
    HEAP<pair<int, pair<DELTA_POS, double> > > heap (betterPath);
    heap.insert(pair<int, pair<DELTA_POS, double> > (nodeIndex, pair<DELTA_POS, double> (displacement, dnet->nodes[nodeIndex].max_radius)));
  
    while(heap.size() != 0){
      // Remove best remaining path
      pair<int, pair<DELTA_POS, double> > pathInfo = heap.pop();
      map<int, pair<DELTA_POS, double> >::iterator nodeInfo = visitedNodes.find(pathInfo.first);
      
      if(pathInfo.second.second < maxRadius)
	break;

      // Ending on a preprocessed node
      if(pathInfo.first == -1){
	bestRadii[nodeIndex] = pathInfo.first;
	break;
      }

      if(nodeInfo != visitedNodes.end()){
	if(nodeInfo->second.first.equals(pathInfo.second.first)){
	  // Circling back to previous node so terminate path
	  continue;
	}
	else {
	  // Found the optimal path through the start node
	  bestRadii[nodeIndex] = min(pathInfo.second.second, nodeInfo->second.second);
	  maxRadius = max(maxRadius, bestRadii[nodeIndex]);
	  break;
	}
      }
      else{
	visitedNodes.insert(pathInfo);
	DIJKSTRA_NODE currentNode = dnet->nodes[pathInfo.first];
	vector<CONN>::iterator connIter = currentNode.connections.begin();
	// Follow all edges leading to adjoining nodes
      
	while(connIter != currentNode.connections.end()){
	  int to = connIter->to;
	  if(nodeProcessed[to]){
	    // Know the best radii for this path because the next node has already been processed
	    // Change its node id to -1 to denote this
	    heap.insert(pair<int, pair<DELTA_POS, double> > (
	      -1, pair<DELTA_POS, double> (DELTA_POS(0,0,0), min(
	      connIter->max_radius, min(pathInfo.second.second, bestRadii[to])))));
	  }
	  else{
	    DELTA_POS newDisplacement = pathInfo.second.first + (
	      connIter->deltaPos*filter);
	    heap.insert(pair<int, pair<DELTA_POS, double> > (
	      to, pair<DELTA_POS, double> (newDisplacement, min(
	      dnet->nodes[to].max_radius, min(pathInfo.second.second,
	      connIter->max_radius)))));
	  }
	  connIter++;
	}
      }
    }
    nodeProcessed[nodeIndex] = true;
    nodeIndex++;
  }
  return maxRadius;
}




bool compareNodeData(pair<int, double> d1, pair<int, double> d2){
  return d1.second > d2.second;
}



/* Calculates the fraction of the node's atomic radii that overlap
   and returns the result in the form of a 2D vector.
 */
vector< vector<double> > calculateNodeOverlap(ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet){
  unsigned int numNodes = vornet->nodes.size();
  vector< vector<double> > overlapFractions = vector< vector<double> > (numNodes, vector<double> (numNodes, 0));
  for(unsigned int i = 0; i < numNodes; i++){
    VOR_NODE n1 = vornet->nodes[i];
    for(unsigned int j = i + 1; j < numNodes; j++){
	VOR_NODE n2 = vornet->nodes[j];
	double sumRadii =  n1.rad_stat_sphere + n2.rad_stat_sphere;
	overlapFractions[i][j] = overlapFractions[j][i] = max(0.0, (sumRadii - atmnet->calcDistanceXYZ(n1.x, n1.y, n1.z, n2.x, n2.y, n2.z))/sumRadii);
      }
   }	    
    return overlapFractions;
}
  

/* Divide the set of node ids into subsets representing individual cages.
 * Returns the results in the form of a vector of sets,
 * where each set represents the node ids that form a single cage. MIN_FRACTION
 * represents the minimum amount of sphere-sphere overlap required for two
 * nodes to belong to the same cluster
 */
vector< set<int> > clusterElements(set<int> cageNodes, vector< vector<double> > overlapFractions, double MIN_FRACTION){
  vector< set<int> > clusters;
  set<int>::iterator nodeIter = cageNodes.begin();
  
  while(nodeIter != cageNodes.end()) {
    int elementID = *nodeIter;
    bool merged = false;
    vector<int> clusterIDs;
    for(unsigned int i = 0; i < clusters.size(); i++){
      set<int>::iterator clusterIter = clusters[i].begin();
      while(clusterIter != clusters[i].end()){
	if(overlapFractions[elementID][*clusterIter] > MIN_FRACTION){
	  clusters[i].insert(elementID);
	  merged = true;
	  clusterIDs.push_back(i);
	  break;
	}
	clusterIter++;
      }
    }
    if(!merged){
      // Node not added to any clusters so form a new cluster
      set<int> newCluster; newCluster.insert(elementID);
      clusters.push_back(newCluster);
    }
    else if (clusterIDs.size() > 1){
      // Node was added to more than one cluster. Therefore,
      // join the resulting clusters into a single cluster

      for(unsigned int i = 1; i < clusterIDs.size(); i++){
	clusters[clusterIDs[0]].insert(clusters[clusterIDs[i]].begin(), clusters[clusterIDs[i]].end());
      }
      int removeCount = 0;
      for(unsigned int i = 1; i < clusterIDs.size(); i++){
	clusters.erase(clusters.begin() + clusterIDs[i] - removeCount);
	removeCount++;
      }
    }
    nodeIter++;
  }
  return clusters;
}


void CAGE::reconstructCage(set<int> ids, ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet){
  int idCount = 0;
  set<int>::iterator iter = ids.begin();
  Point pivot = atmnet->xyz_to_abc(vornet->nodes[*iter].x, vornet->nodes[*iter].y, vornet->nodes[*iter].z);
  vector<Point> nodeXYZCoords;   nodeXYZCoords.push_back(Point(vornet->nodes[*iter].x, vornet->nodes[*iter].y, vornet->nodes[*iter].z));
  nodeIDs.push_back(*iter); 
  nodeOffsets.push_back(DELTA_POS(0,0,0));
  nodeIndices.insert(pair<int,int>(*iter, idCount));
  iter++;
  idCount++;
  
  // Calculate and store the closest positions of other points to pivot point
  while(iter != ids.end()) {
    double minDa = DBL_MAX, minDb = DBL_MAX, minDc = DBL_MAX;
    Point other = atmnet->xyz_to_abc(vornet->nodes[*iter].x,  vornet->nodes[*iter].y, vornet->nodes[*iter].z);
    atmnet->getDistCalc().minimum_periodic_distance(pivot[0], pivot[1], pivot[2], other[0], other[1], other[2], minDa, minDb, minDc);
    Point shifted = Point(minDa, minDb, minDc).add(pivot);
    nodeXYZCoords.push_back(atmnet->abc_to_xyz(shifted));
    
    Point offset = shifted.subtract(other);
    nodeIDs.push_back(*iter);  
    nodeOffsets.push_back(DELTA_POS((int)(floor(offset[0] + 0.5)), (int)(floor(offset[1] + 0.5)), (int)(floor(offset[2] + 0.5))));
    nodeIndices.insert(pair<int,int>(*iter, idCount));
    iter++;
    idCount++;
  }
  
  // Compute the center point of the cage
  cx = cy = cz = 0;
  for(unsigned int i = 0; i < nodeXYZCoords.size(); i++){
    cx += nodeXYZCoords[i][0]; cy += nodeXYZCoords[i][1]; cz += nodeXYZCoords[i][2];
  }
  cx /= nodeXYZCoords.size(); cy /= nodeXYZCoords.size(); cz /= nodeXYZCoords.size();
  
  // Compute the radius of the cage
  iter = ids.begin();
  radius = DBL_MAX;
  for(unsigned int i = 0; i < nodeXYZCoords.size(); i++){
    radius = min(radius, calcEuclideanDistance(cx, cy, cz, nodeXYZCoords[i][0], nodeXYZCoords[i][1], nodeXYZCoords[i][2]) + vornet->nodes[*iter].rad_stat_sphere);
    iter++;
  }
}

CAGE::CAGE(set<int> ids, ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet){
  reconstructCage(ids, atmnet, vornet);
}

DELTA_POS CAGE::getNodeOffset(int nodeID){
  map<int,int>::iterator iter = nodeIndices.find(nodeID);
  if(iter == nodeIndices.end()){
    cerr << "Error: Node #" << nodeID << " not found in cage." << "\n"
	 << "Exiting..." << "\n";
    exit(1);
  }
  else
    return nodeOffsets[iter->second];
}

double CAGE::getRadius() { return radius; }
Point  CAGE::getCenter() { return Point(cx, cy, cz); }
vector<int> CAGE::getNodeIDs() { return nodeIDs; }
void CAGE::writeToVMD(int cageIndex, fstream &output){
  output << "set cages(" << cageIndex << ") {" << "\n"
	 << "{color $cageColors(" << cageIndex << ")}" << "\n"
	 << "{sphere {" << cx << " " << cy << " " << cz << "} radius " << radius << " resolution 100 }" << "\n"
	 <<"}" << "\n";
}


 // Defined in area_and_volume.cc
 double calcDensity(ATOM_NETWORK *);

 
/* Identify the cages present in Voronoi network for the given probe size. If the
 * visualization option is specified, commands necessary to visualize the cages in 
 * ZeoVis is written to the provided file stream. Otherwise, the cages are classifed 
 * according to radius and the relevant cage sizes and volume is output using the file stream.
 */
void identifyCages(ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet, DIJKSTRA_NETWORK *dnet, double probeRadius, bool visualize, fstream &output, vector<CAGE> &cages) {
  vector<bool> nodeProcessed = vector<bool>(dnet->nodes.size(), false);
  vector<double> bestRadii   = vector<double>(dnet->nodes.size(), -1);
  unsigned int nodeIndex = 0;

  double maxRadius = 0.0;
  DELTA_POS filter(1,1,1);

  // Iterate over all nodes
  while(nodeIndex < nodeProcessed.size()){
   
    // Place starting node on stack with (0,0,0) displacement
    DELTA_POS displacement(0,0,0);
    map<int, pair<DELTA_POS, double> > visitedNodes;
    
    HEAP<pair<int, pair<DELTA_POS, double> > > heap (betterPath);
    heap.insert(pair<int, pair<DELTA_POS, double> > (nodeIndex, pair<DELTA_POS, double> (displacement, dnet->nodes[nodeIndex].max_radius)));
  
    while(heap.size() != 0){
      // Remove best remaining path
      pair<int, pair<DELTA_POS, double> > pathInfo = heap.pop();
      map<int, pair<DELTA_POS, double> >::iterator nodeInfo = visitedNodes.find(pathInfo.first);
      
      // Ending on a preprocessed node
      if(pathInfo.first == -1){
	bestRadii[nodeIndex] = pathInfo.second.second;
	break;
      }

      if(nodeInfo != visitedNodes.end()){
	if(nodeInfo->second.first.equals(pathInfo.second.first)){
	  // Circling back to previous node so terminate path
	  continue;
	}
	else {
	  // Found the optimal path through the start node
	  bestRadii[nodeIndex] = min(pathInfo.second.second, nodeInfo->second.second);
	  maxRadius = max(maxRadius, bestRadii[nodeIndex]);
	  break;
	}
      }
      else{
	visitedNodes.insert(pathInfo);
	DIJKSTRA_NODE currentNode = dnet->nodes[pathInfo.first];
	vector<CONN>::iterator connIter = currentNode.connections.begin();
	// Follow all edges leading to adjoining nodes
      
	while(connIter != currentNode.connections.end()){
	  int to = connIter->to;
	  if(nodeProcessed[to]){
	    // Know the best radii for this path because the next node has already been processed
	    // Change the path's node id to -1 to denote this
	    heap.insert(pair<int, pair<DELTA_POS, double> > (-1, pair<DELTA_POS, double> (DELTA_POS(0,0,0), min(connIter->max_radius, min(pathInfo.second.second, bestRadii[to])))));
	  }
	  else{
	    DELTA_POS newDisplacement = pathInfo.second.first + (connIter->deltaPos*filter);
	    heap.insert(pair<int, pair<DELTA_POS, double> > (to, pair<DELTA_POS, double> (newDisplacement, 
											  min(dnet->nodes[to].max_radius, min(pathInfo.second.second, connIter->max_radius)))));
	  }
	  connIter++;
	}
      }
    }
    nodeProcessed[nodeIndex] = true;
    nodeIndex++;
  }

  double CUTOFF_RATIO = 0.8; 

  // Identify nodes that satisfy cage criteria
  set<int> cageElements;
  for(unsigned int i = 0; i < vornet->nodes.size(); i++){
    double ratio = bestRadii[i]/(dnet->nodes[i].max_radius);
    if(ratio < CUTOFF_RATIO && bestRadii[i] > probeRadius) {
      cageElements.insert(i);
    }
  }
  
  // Calculate node overlaps
  vector< vector<double> > overlapFractions = calculateNodeOverlap(atmnet, vornet);

  // Cluster nodes in cages
  vector< set<int> > clusters = clusterElements(cageElements, overlapFractions, 0.1);


  //vector<CAGE> cages;
  for(unsigned int i = 0; i < clusters.size(); i++){
    cages.push_back(CAGE(clusters[i], atmnet, vornet));
  }

  if(visualize){
    output << "set num_cages " << clusters.size() << "\n";
    for(unsigned int i = 0; i < cages.size(); i++){
      cages[i].writeToVMD(i, output);
    }
  }
  else{
    vector<double> cageRadii;
    for(unsigned int i = 0; i < cages.size(); i++){
      cageRadii.push_back(cages[i].getRadius());
    }

     sort(cageRadii.begin(), cageRadii.end());
    
    vector<int> counts;
    vector<double> radii;
    double prevRad = -1, tolerance = 0.001;
    int numCageTypes = 0;
    for(unsigned int i = 0; i < cageRadii.size(); i++){
      if(abs(cageRadii[i] - prevRad) < tolerance)
	counts[numCageTypes - 1]++;
      else{
	counts.push_back(1); radii.push_back(cageRadii[i]);
	prevRad = cageRadii[i];	numCageTypes++;
      }
    }

    double totalCageAV = 0;
    output << "Note: Probe radius is subtracted from cage radius" << "\n";
    for(unsigned int i = 0; i < counts.size(); i++){
      output << counts[i] << " cage(s) of radius " << radii[i] - probeRadius << "\n";
      totalCageAV += calcSphereVolume(radii[i] - probeRadius)*counts[i];
    }
    output << "Total Cage AV: " << totalCageAV << " A^3" << "\n";
    double volumeFraction = totalCageAV/calcDeterminant(atmnet->ucVectors);
    double rho_crystal    = calcDensity(atmnet);
    double cageAVPerMass  = volumeFraction/rho_crystal;
    output << "Total Cage AV/Mass: " << cageAVPerMass << " cm^3/g" <<  "\n";
  }
}







/* Determine cage interconnectivity information by determining which cages are directly connected to one another. */
void simplifyCageNetwork(ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet, DIJKSTRA_NETWORK *dnet, double probe_radius){
  fstream output;
  output.open("testing.cage", fstream::out);
  vector<CAGE> cages;
  identifyCages(atmnet, vornet, dnet, probe_radius, true, output, cages);
   
  vector<Point> cageCenters;
  for(unsigned int i = 0; i < cages.size(); i++){
    cageCenters.push_back(cages[i].getCenter());
  }

  vector<int> nodeCageIDs = vector<int>(vornet->nodes.size(), -1);

  for(unsigned int i = 0; i < cages.size(); i++){
    vector<int> ids = cages[i].getNodeIDs();
    for(unsigned int j = 0; j < ids.size(); j++){
      nodeCageIDs[ids[j]] = i;
    }
  }

  set<DELTA_POS, bool(*)(DELTA_POS,DELTA_POS)> emptySet (deltaPosLessThan);
  vector< vector< set<DELTA_POS, bool(*)(DELTA_POS,DELTA_POS)> > > cageConnections = vector< vector< set<DELTA_POS, bool(*)(DELTA_POS,DELTA_POS)> > >(cages.size(), 
											   vector< set<DELTA_POS, bool(*)(DELTA_POS,DELTA_POS)> >(cages.size(), emptySet));
  for(unsigned int i = 0; i < dnet->nodes.size(); i++){
    int nodeID = nodeCageIDs[i];
    if(nodeID != -1){
      vector<CONN>::iterator connIter = dnet->nodes[i].connections.begin();
      while(connIter != dnet->nodes[i].connections.end()){
	int fromID = nodeCageIDs[connIter->from];
	int toID   = nodeCageIDs[connIter->to];
	if((fromID != -1) && (toID != -1)){
	  if(fromID != toID){
	    cageConnections[fromID][toID].insert(connIter->deltaPos);
	  }
	  else{
	    DELTA_POS cageDeltaPos = cages[toID].getNodeOffset(connIter->to) - (cages[fromID].getNodeOffset(i));
	    if(!cageDeltaPos.equals(connIter->deltaPos)){
	      cageConnections[fromID][toID].insert(connIter->deltaPos);
	    }
	  }
	}
	connIter++;
      }
    }
  }

  output << "set vornets(0) {" << "\n"
	 << "{color red} " << "\n";
  for(unsigned int i = 0; i < cageConnections.size(); i++){
    cout << i << "\n";
    for(unsigned int j = 0; j < cageConnections[i].size(); j++){
      if(cageConnections[i][j].size() > 0){
	cout << "\t ->" << j << "\n";
	set<DELTA_POS, bool(*)(DELTA_POS,DELTA_POS)>::iterator iter = cageConnections[i][j].begin();
	while(iter != cageConnections[i][j].end()){
	  DELTA_POS dp = *iter;
	  if((i != j) || (dp.x != 0) || (dp.y != 0) || (dp.z != 0)){
	    cout << "\t\t" <<  dp.x << " " << dp.y << " " << dp.z << "\n";
	    Point other = Point(cageCenters[j][0], cageCenters[j][1], cageCenters[j][2]);
	    atmnet->translatePoint(&other, dp.x, dp.y, dp.z);
	    output << "{line {" << cageCenters[i][0] << " " << cageCenters[i][1] << " " << cageCenters[i][2] << "} " 
		   <<       "{" << other[0] << " " << other[1] << " " << other[2] << "} }" << "\n";
	  }
	  iter++;
	}
      }
    }
  }
  output << "}" << "\n";
  output.close();
}
