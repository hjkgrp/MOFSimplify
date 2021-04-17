//#include "network.h"
#include "channel.h"
#include <cmath>

using namespace std;

/* Create a PORE that does not contain any nodes or connections
 * and spans 0 unit cells.*/
PORE::PORE(){
    nodes = vector<DIJKSTRA_NODE> ();
    connections = vector<CONN> ();
    unitCells = vector<DELTA_POS> ();
    ucNodes = vector< vector<int> > ();
    dimensionality = 0;
    basis[0][0]=0;basis[0][1]=0;basis[0][2]=0;
    basis[1][0]=0;basis[1][1]=0;basis[1][2]=0;
    basis[2][0]=0;basis[2][1]=0;basis[2][2]=0;
}

/* Returns true iff the CHANNEL unit can be depicted within one unit cell*/
bool CHANNEL::isUnicellular(){
    return unitCells.size() == 1;
}

/* Constructs a pore from the provided  nodes in the DIJKSTRA_NETWORK.
 *  Reconstructs the pore by trying to minimize the number of unit cells
 *  required to show a single pore unit.
 */
PORE::PORE(vector<int> nodeIDs, DIJKSTRA_NETWORK *dnet, int dim, int basisVecs[3][3]){
    idMappings = map<int,int> ();
    reverseIDMappings = map<int,int> ();
    nodes = vector<DIJKSTRA_NODE> ();
    connections = vector<CONN> ();
    
    // Reindex the nodes to facilitate storage
    for(unsigned int i = 0; i < nodeIDs.size(); i++){
        idMappings.insert(pair<int,int>(nodeIDs.at(i),i));
        reverseIDMappings.insert(pair<int,int>(i,nodeIDs.at(i)));
    }
    
    // Iterate over all nodes in list
    for(unsigned int i = 0; i < nodeIDs.size(); i++){
        DIJKSTRA_NODE oldNode = dnet->nodes.at(nodeIDs.at(i));
        DIJKSTRA_NODE newNode = DIJKSTRA_NODE(i, oldNode.x, oldNode.y, oldNode.z, oldNode.max_radius);
        
        // Include connection only if end node is in pore
        // Reindex connection start and end ids
        for(unsigned int j = 0; j < oldNode.connections.size(); j++){
            CONN oldConn = oldNode.connections.at(j);
            map<int,int>::iterator id1Iter = idMappings.find(oldConn.from);
            map<int,int>::iterator id2Iter = idMappings.find(oldConn.to);
            if(id2Iter != idMappings.end()){
                CONN newConn = CONN(id1Iter->second,id2Iter->second,oldConn.length,oldConn.max_radius,oldConn.deltaPos);
                newNode.connections.push_back(newConn);
                connections.push_back(newConn);
            }
        }
        nodes.push_back(newNode);
    }
    
    basis[0][0] = basisVecs[0][0];  basis[0][1] = basisVecs[0][1];  basis[0][2] = basisVecs[0][2];
    basis[1][0] = basisVecs[1][0];  basis[1][1] = basisVecs[1][1];  basis[1][2] = basisVecs[1][2];
    basis[2][0] = basisVecs[2][0];  basis[2][1] = basisVecs[2][1];  basis[2][2] = basisVecs[2][2];
    dimensionality = dim;
    
    // Store unit cell vectors
    v_a = dnet->v_a;
    v_b = dnet->v_b;
    v_c = dnet->v_c;
    reconstruct();    // Reconstruct pore
}

/* Provides a vector with IDs of Voronoi nodes correcponding to the current pore */
vector <int> PORE::nodeIds(){
    vector <int> list;
    
    for(unsigned int i=0; i<reverseIDMappings.size(); i++)
    {
        list.push_back(reverseIDMappings.find(i)->second);
    };
    return list;
}





/** Nodes within the provided DIJKSTRA_NETWORK can be classified as either accessible or inaccessible,
 *  where sets of accessible  nodes constitute a CHANNEL and inaccessible POCKET. This function identifies
 *  the POREs that exist for the provided particle diameter and stores them using the provided pointer
 *  to a vector of POREs. CHANNEL and PORE can be distinguished by dimentionality.
 *  In addition, the pointer to the vector of bools is used to a store bool
 *  for each VORONOI_NODE, where infoStorage[i] is true if node #i is accessible. */
/*  WARNING: this function assumes that at some point before VORONOI_NETWORK was pruned
 *  with prune() function. This is to set "active" flag on nodes larger than probe radius */
void PORE::findChannelsAndPockets(DIJKSTRA_NETWORK *dnet, vector<bool> *infoStorage,
                                  vector<PORE> *pores)
{
    //Define three access types
    const int ACCESSIBLE   =  1;
    const int INACCESSIBLE =  0;
    const int UNKNOWN      = -1;
    
    vector<DELTA_POS> directions = vector<DELTA_POS> ();
    
    // Initialize the status of each node as UNKNOWN
    vector<int> accessStatuses = vector<int> (dnet->nodes.size(), UNKNOWN);
    unsigned int nodeIndex = 0;
    
    // Get and Print basic info
    int noAccDijkNodes = 0; // number of active Dijkstra nodes (ones with radius > 0)
    int noAccDijkNodesAssigned = 0; // number of active Dijkstra nodes already assigned to pores.
    for(unsigned int i = 0; i < dnet->nodes.size(); i++)
    {
        if(dnet->nodes.at(i).active == true) noAccDijkNodes++;
        
        //debug line
        //cout << "id,active " << i << "   " << dnet->nodes.at(i).active << endl;
    };
    cout << "\nFinding channels and pockets in Dijkstra network of " << dnet->nodes.size() << " node(s). " << noAccDijkNodes << " are expected to compose pores." << endl;
    
    // Iterate over all nodes
    while(nodeIndex != accessStatuses.size()){
        //Skip node if access status already determined
        if(accessStatuses.at(nodeIndex) != UNKNOWN){
            nodeIndex++;
            continue;
        }
        
        // Start out with 0-dimensional pore w/o any basis vectors
        int dim = 0;
        int basis [3][3] = {{0,0,0},
			{0,0,0},
			{0,0,0}};
        
        // Place starting node on stack with (0,0,0) displacement
        int accessStatus = UNKNOWN;
        DELTA_POS displacement = DELTA_POS(0,0,0);
        map<int,DELTA_POS> visitedNodeDisplacement;
        vector<pair<int,DELTA_POS> > stack;
        stack.push_back(pair<int,DELTA_POS> (nodeIndex,displacement));
        visitedNodeDisplacement.insert(pair<int,DELTA_POS>(nodeIndex, displacement));
        
        while(stack.size() != 0){
            // Remove top-most node
            pair<int,DELTA_POS> nodeInfo = stack.back();
            DIJKSTRA_NODE currentNode = dnet->nodes.at(nodeInfo.first);
            stack.pop_back();
            
            // Follow all edges leading to adjoining nodes
            vector<CONN>::iterator connIter = currentNode.connections.begin();
            while(connIter != currentNode.connections.end()){
                int to = connIter->to;
                DELTA_POS newDisplacement = nodeInfo.second + connIter->deltaPos;
                int localAccessStatus = accessStatuses.at(to);
                
                if(localAccessStatus == UNKNOWN){
                    map<int,DELTA_POS>::iterator visitedNode = visitedNodeDisplacement.find(to);
                    if(visitedNode != visitedNodeDisplacement.end()){
                        if(visitedNode->second.equals(newDisplacement)){
                            // Circling back to previous node
                            // Do nothing
                        }
                        else{
                            // Nodes are ACCESSIBLE b/c same node visited in different unit cell
                            accessStatus = ACCESSIBLE;
                            DELTA_POS direction = newDisplacement - visitedNode->second;
                            
                            
                            /******* Start of Michael's code ******/
                            if(dim == 0){
                                basis[0][0] = direction.x;
                                basis[1][0] = direction.y;
                                basis[2][0] = direction.z;
                                dim++;
                            }
                            else if (dim == 1){
                                double v;
                                
                                if (basis[0][0] != 0) v = 1.0*direction.x/basis[0][0];
                                else if (basis[1][0] != 0) v = 1.0*direction.y/basis[1][0];
                                else if (basis[2][0] != 0) v = 1.0*direction.z/basis[2][0];
                                else {
                                    cerr << "Error: Pore basis vector is zero vector. Exiting..." << "\n";
                                    exit(1);
                                }
                                
                                if(!(basis[0][0]*v == direction.x && basis[1][0]*v == direction.y &&
                                     basis[2][0]*v == direction.z)){
                                    basis[0][1] = direction.x;
                                    basis[1][1] = direction.y;
                                    basis[2][1] = direction.z;
                                    dim++;
                                }
                            }
                            else if (dim == 2){
                                basis[0][2] = direction.x;
                                basis[1][2] = direction.y;
                                basis[2][2] = direction.z;
                                if(abs(calcDeterminant(basis)) < 1e-8){
                                    basis[0][2] = basis[1][2] = basis[2][2] = 0;
                                }
                                else
                                    dim++;
                            }
                            /****** End of Michael's code ******/
                        }
                    }
                    else{
                        // Node being visited for the first time
                        visitedNodeDisplacement.insert(pair<int,DELTA_POS>(to, newDisplacement));
                        stack.push_back(pair<int,DELTA_POS> (to, newDisplacement));
                    }
                }
                else{
                    // The status of a region should be determined
                    // with all of its neighbors. Otherwise, all possible nodes
                    // were not explored.
                    cerr << "Error: Illogical result  when attempting to identify channels/pockets." << "\n";
                    cerr << "Please contact the source code provider with your program input. " << "\n";
                    cerr << "Exiting ..." << "\n";
                    exit(1);
                }
                connIter++;
            }
        }
        
        // All connected nodes are inaccessible if their status is still UNKNOWN
        if((stack.size() == 0) && (accessStatus == UNKNOWN))
            accessStatus = INACCESSIBLE;
        
        // Record the access status for each node in the cycle and create
        // a list of the node ids for channel identification
        map<int,DELTA_POS>::iterator resultIter = visitedNodeDisplacement.begin();
        vector<int> listOfIDs = vector<int> ();
        while(resultIter != visitedNodeDisplacement.end()){
            
            int nodeID = resultIter->first;
            if(accessStatuses.at(nodeID) != UNKNOWN){
                // Each node should only have its accessibility determined once
                cerr << "Error: Accessibility of node was determined more than once." << "\n";
                cerr << "Please contact the source code provider with your program input. " << "\n";
                cerr << "Exiting ..." << "\n";
                exit(1);
            }
            else{
                // Store node's status
                accessStatuses.at(nodeID) = accessStatus;
                listOfIDs.push_back(nodeID);
            }
            resultIter++;
        }
        nodeIndex++;
        
        // Create CHANNEL from ACCESSIBLE nodes
        if(accessStatus == ACCESSIBLE){
            pores->push_back(PORE(listOfIDs, dnet, dim, basis));
            // store pore statistics (number of nodes included in pore)
            noAccDijkNodesAssigned += listOfIDs.size();
            //output << dim << " ";
        }else
        { // OTHERWISE create POCKET from INACCESSIBLE nodes (dim set to 0 indicates pocket)
            // additionally pockets have to be built from active nodes (this is to ensure their
            // size is larger than mix_radius)
            if(dnet->nodes[listOfIDs[0]].active == true)
            {
                pores->push_back(PORE(listOfIDs, dnet, dim, basis));
                // store pore statistics (number of nodes included in pore)
                noAccDijkNodesAssigned += listOfIDs.size();
            };
        };
        
        //more debug
        bool flg;
        for(unsigned int i=0;i<listOfIDs.size();i++)
        {
            if(i==0) flg = dnet->nodes[listOfIDs[0]].active;
            else
            {
                if(flg !=dnet->nodes[listOfIDs[i]].active) cout << "pocket error" << endl;
            };
        };
        //end of debug
    }
    // Store information in provided vector
    infoStorage -> resize(accessStatuses.size());
    for(unsigned int i = 0; i < accessStatuses.size(); i++){
        infoStorage -> at(i) = (accessStatuses.at(i) == ACCESSIBLE);
    }
    
    // Print summary
    int nchannel=0,npocket=0;
    for(unsigned int i = 0; i<pores->size(); i++)
    {
        if(pores->at(i).dimensionality>0) nchannel++; else npocket++;
    };
    cout << "Analyzed and assigned " << nodeIndex << " nodes.";
    cout << "\nIdentified " << nchannel << " channels and " << npocket << " pockets.\n";
    cout << noAccDijkNodesAssigned << " nodes assigned to pores. " << endl;
}

/*  Voronoi nodes within the provided VORONOI_NETWORK can be classified as either
 *  accessible or inaccessible based on a give probe diameter. This is done by constructing
 *  DIJKSTRA_NETWORK and calling the other findChannelsAndPockets function */
void PORE::findChannelsAndPockets(VORONOI_NETWORK *vornet, double minRadius,
                                  vector<bool> *infoStorage, vector<PORE> *pores)
{
    //Remove edges that don't allow the provided particle diameter to
    //pass freely
    //VORONOI_NETWORK newNetwork;
    //pruneVoronoiNetwork(vornet, &newNetwork, minRadius);
    VORONOI_NETWORK newNetwork = vornet->prune(minRadius);
    
    //Build graph data structure
    DIJKSTRA_NETWORK dnet;
    DIJKSTRA_NETWORK::buildDijkstraNetwork(&newNetwork, &dnet);
    findChannelsAndPockets(&dnet, infoStorage, pores);
}


/** Calculates center of mass and internal void radii (distance between the center of mass and its nearest atom)
 *  This function attempts pore reconstruction in cases where pores cross the cell boundaries. 
 */
pair <XYZ, double> PORE::getCenterOfMass(){

/*
 Pore variables
    std::vector<DELTA_POS> unitCells;    // List of displacements of unitcells in single PORE unit
    std::vector< std::vector<int> > ucNodes;  // List of nodes contained in each unit cell
    XYZ v_a, v_b, v_c;              // Unit cell vectors

*/

 vector<XYZ> reconstructedNodeCoords;
 XYZ centerOfMass (0,0,0);
 double dist; // distance from the center of mass to the nearest node

 XYZ tempPt;
 for(int unsigned i=0; i<unitCells.size(); i++)
    {
    for(int unsigned j=0; j<ucNodes[i].size(); j++)
       {
       int nd = ucNodes[i].at(j); // current node

       XYZ deltapos(double(unitCells[i].x), double(unitCells[i].y), double(unitCells[i].z));

       XYZ node(nodes[nd].x, nodes[nd].y, nodes[nd].z);

       node = node + v_a.scale(deltapos.x) + v_b.scale(deltapos.y) + v_c.scale(deltapos.z) ;

       reconstructedNodeCoords.push_back(node);
       centerOfMass = centerOfMass + node;
       };
    };

 centerOfMass = centerOfMass.scale(1.0/double(nodes.size()));

 for(int unsigned i=0; i<reconstructedNodeCoords.size(); i++)
    {
    if(i==0)
      {
      dist = centerOfMass.euclid_dist(reconstructedNodeCoords[i]);
      }
    else
      {
      if(centerOfMass.euclid_dist(reconstructedNodeCoords[i]) < dist) dist = centerOfMass.euclid_dist(reconstructedNodeCoords[i]);
      };
    };

 if(dimensionality>0) cout << "Center of Mass calculation: PORE dimensionality>0; results may be buggy due to PORE reconstruction handling\n";

/* DEBUGGING: Printing out coordinates
 for(int unsigned i=0; i<reconstructedNodeCoords.size(); i++)
    {
    cout << "C    " << reconstructedNodeCoords[i].x << "  " << reconstructedNodeCoords[i].y << "  " << reconstructedNodeCoords[i].z << "\n";
    };
*/

 pair<XYZ, double> p(centerOfMass, dist);
 return p;

} // ends getCenterOfMass() function



/**  This function attempts pore reconstruction in cases where pores cross the cell boundaries. 
 *   It returns a vector of vectors with pore XYZ coordinates
 *   */
 vector< pair <int,XYZ> > PORE::getReconstructedPore(){

/*
 *  Pore variables
 *      std::vector<DELTA_POS> unitCells;    // List of displacements of unitcells in single PORE unit
 *          std::vector< std::vector<int> > ucNodes;  // List of nodes contained in each unit cell
 *              XYZ v_a, v_b, v_c;              // Unit cell vectors
 *
 *              */

 vector< pair <int,XYZ> > reconstructedNodeCoords;

 for(int unsigned i=0; i<unitCells.size(); i++)
    {
    for(int unsigned j=0; j<ucNodes[i].size(); j++)
       {
       int nd = ucNodes[i].at(j); // current node

       XYZ deltapos(double(unitCells[i].x), double(unitCells[i].y), double(unitCells[i].z));

       XYZ node(nodes[nd].x, nodes[nd].y, nodes[nd].z);

       node = node + v_a.scale(deltapos.x) + v_b.scale(deltapos.y) + v_c.scale(deltapos.z) ;

       pair <int,XYZ> p (nd,node);

       reconstructedNodeCoords.push_back(p);
       };
    };

/* DEBUGGING: Check if reconstructured nodes become duplicate atoms

 cout << "DEBUG: Pore/molecule size after reconstruction (XYZ coords) " << reconstructedNodeCoords.size() << endl;  
*/

 if(dimensionality>0) cout << "Calling PORE::getReconstructedPore for a pore with dim>0, it was not intended. DO NOT TRUST\n";

/* DEBUGGING: Printing out coordinates
 *  for(int unsigned i=0; i<reconstructedNodeCoords.size(); i++)
 *      {
 *          cout << "C    " << reconstructedNodeCoords[i].x << "  " << reconstructedNodeCoords[i].y << "  " << reconstructedNodeCoords[i].z << "\n";
 *              };
 *              */

 return reconstructedNodeCoords;

} // ends getReconstructredPore() function




/**  This function attempts pore reconstruction in cases where pores cross the cell boundaries.
 *   In case of pores/molecules that cross boundy, multiple copies are saved each corresponding to one instance in periodic box.
 *   It returns a vector of vectors with pore XYZ coordinates and Ids
 *   */

 vector<  vector< pair <int,XYZ> > >  PORE::getReconstructredPoresWithCrossBoundryCopies(){

/*
 *  Pore variables
 *      std::vector<DELTA_POS> unitCells;    // List of displacements of unitcells in single PORE unit
 *          std::vector< std::vector<int> > ucNodes;  // List of nodes contained in each unit cell
 *              XYZ v_a, v_b, v_c;              // Unit cell vectors
 *
 *              */

 vector< pair <int,XYZ> > reconstructedNodeCoords; /* holds the reconstructure pore/molecule */

 for(int unsigned i=0; i<unitCells.size(); i++)
    {
    for(int unsigned j=0; j<ucNodes[i].size(); j++)
       {
       int nd = ucNodes[i].at(j); // current node

       XYZ deltapos(double(unitCells[i].x), double(unitCells[i].y), double(unitCells[i].z));

       XYZ node(nodes[nd].x, nodes[nd].y, nodes[nd].z);

       node = node + v_a.scale(deltapos.x) + v_b.scale(deltapos.y) + v_c.scale(deltapos.z) ;

       pair <int,XYZ> p (nd,node);

       reconstructedNodeCoords.push_back(p);
       };
    };

 vector<  vector< pair <int,XYZ> > > reconstructedCopies; // This will store copies of all reconstructure pores/molecules


 if(dimensionality>0) cout << "Calling PORE::getReconstructredPoresWithCrossBoundryCopies() for a pore with dim>0, it was not intended. DO NOT TRUST\n";

 for(int unsigned i=0; i<unitCells.size(); i++)
    {
    reconstructedCopies.push_back(reconstructedNodeCoords);

    for(int unsigned j=0; j<reconstructedNodeCoords.size(); j++)
       {

       XYZ deltapos(double(unitCells[i].x), double(unitCells[i].y), double(unitCells[i].z));

       reconstructedCopies[i].at(j).second = reconstructedCopies[i].at(j).second - v_a.scale(deltapos.x) - v_b.scale(deltapos.y) - v_c.scale(deltapos.z) ; 

       };
    };

 return reconstructedCopies;

} // ends getReconstructredPoresWithCrossBoundryCopies() function





/* Prints information about the pore to the provided output stream, including nuber of nodes, the largest included sphere,
 positions and radii of all nodes */
void PORE::printPoreSummary(ostream &out, ATOM_NETWORK *atmNet){
    
    vector <double> poresummary;
    getSimplifiedPocketInfo(atmNet, &poresummary);
    
    out << nodes.size() << "  " << poresummary[0] << "  " << poresummary[1] << "  " << poresummary[2] << "  " << poresummary[3] << "  " << poresummary[4] << "\n";
    
    //  out << nodes.size() << "  " << getIncludedSphereDiameter()  << "\n";
    
    for(unsigned int i = 0; i< nodes.size(); i++){
        //    nodes.at(i).print();
        Point pt = atmNet->xyz_to_abc(nodes.at(i).x, nodes.at(i).y, nodes.at(i).z);
        pt = atmNet->shiftABCInUC(pt);
        out << pt[0] << "  " << pt[1] << "  " << pt[2];
        out << "    " << nodes.at(i).max_radius << "\n";
    };
}

/* Prints information about the CHANNEL to the provided output stream, including
 * the number of nodes, unitcells, and the nodes located in each unit cell. Additional
 * node information is outputted if requested.*/
void CHANNEL::print(ostream &out, bool dispNodeInfo){
    out << "Channel info:" << "\n";
    out << "     # Nodes: " << nodes.size() << "\n";
    
    if(dispNodeInfo){
        out << "     Original Node IDs: ";
        for(unsigned int i = 0; i < nodes.size(); i++){
            out << reverseIDMappings.find(i)->second << " ";
        }
        out << "\n";
        
        out << "     New Node IDs: ";
        for(unsigned int i = 0; i < nodes.size(); i++){
            out << i << " ";
        }
        out << "\n";
        
        out << "  New Node info: " << "\n";
        for(unsigned int i = 0; i< nodes.size(); i++){
            nodes.at(i).print();
        }
    }
    
    out << "     # Unit cells:" << unitCells.size() << "\n";
    
    for(unsigned int i = 0; i < unitCells.size(); i++){
        DELTA_POS position = unitCells.at(i);
        vector<int> ucNode = ucNodes.at(i);
        out << "       Unit cell #: " << i  << "\n"
        << "          Displacement: " << position.x << " " << position.y << " " << position.z << "\n"
        << "          New Node ids: ";
        for(unsigned int j = 0; j < ucNode.size(); j++){
            out << ucNode.at(j) << " ";
        }
        out << "\n";
    }
}


bool(*fn_pt)(DELTA_POS,DELTA_POS) = deltaPosLessThan;




/* Reset the visited positions and set the current position to the origin. */
ReconstructorComparator::ReconstructorComparator(){
    positions = set<DELTA_POS, bool(*)(DELTA_POS,DELTA_POS)> (deltaPosLessThan);
    currentPos = DELTA_POS(0,0,0);
}

/* Set the current position and store the old position. */
void ReconstructorComparator::setPosition(DELTA_POS p){
    positions.insert(p);
    currentPos = p;
}

bool ReconstructorComparator::compare(pair<int,DELTA_POS> p1, pair<int,DELTA_POS> p2) {
    bool p1Current = p1.second.equals(currentPos);
    bool p2Current = p2.second.equals(currentPos);
    if(p1Current && p2Current)
        return false;
    else if (p1Current)
        return false;
    else if (p2Current)
        return true;
    else{
        bool foundP1 = (positions.find(p1.second) != positions.end());
        bool foundP2 = (positions.find(p2.second) != positions.end());
        if(foundP1)
            return false;
        else
            return foundP2;
    }
}



ReconstructorComparator comparer; // Object used to reconstruct CHANNEL
bool compareNodes(pair<int,DELTA_POS> p1, pair<int,DELTA_POS> p2){
    return comparer.compare(p1,p2);
}



/* Reconstructs the PORE by propagating paths until all nodes have been accessed.
 *  Stores each node in the unit cell in which it is encountered.
 *  Attempts to reduce the number of unit cells required for reconstruction by
 *  favoring nodes in already-accessed unit cells over nodes in new unit cells. */
void PORE::reconstruct(){
    vector<bool> haveVisited = vector<bool>(nodes.size(),false);
    vector<DELTA_POS> displacements = vector<DELTA_POS>(nodes.size());
    unsigned int visitCount = 0;
    
    //Search through the PORE, emphasizing nodes located in the current unit
    //cell or in previously visited unit cells. Record the unit cell within which
    //each node is encountered
    comparer = ReconstructorComparator();
    HEAP< pair<int,DELTA_POS> > stack (compareNodes);
    stack.insert(pair<int,DELTA_POS> (0, DELTA_POS(0,0,0))); // Pick the first node as the starting point
    DELTA_POS curPos = DELTA_POS(0,0,0);
    
    // Continue reconstruction until all nodes have been visited
    while(visitCount < nodes.size()){
        if(stack.size() == 0){
            cerr << "Error: Stack empties prior to pore reconstruction completion." << "\n"
            << "Please contact the source code provided with this message." << "\n"
            << "Exiting..." << "\n"
            << "Nnodes = " << nodes.size() << "  visitCount= " << visitCount << "\n";
            
            exit(1);
        }
        
        pair<int,DELTA_POS> best = stack.pop();
        
        if(!haveVisited.at(best.first)){
            visitCount++;
            haveVisited.at(best.first) = true;
            displacements.at(best.first) = best.second;
            
            // If changed unit cell, need to restructure heap
            if(!best.second.equals(curPos)){
                comparer.setPosition(best.second);
                stack.reHeapify();
                curPos = best.second;
            }
            
            // Add all the nodes from the current node to the stack that
            // have not yet been visited
            DIJKSTRA_NODE curNode = nodes.at(best.first);
            for(unsigned int j = 0; j < curNode.connections.size(); j++){
                CONN curConn = curNode.connections.at(j);
                if(!haveVisited.at(curConn.to)){
                    DELTA_POS newPos = curPos + curConn.deltaPos;
                    stack.insert(pair<int,DELTA_POS>(curConn.to,newPos));
                }
            }
        }
    }
    
    map<DELTA_POS, vector<int> ,  bool(*)(DELTA_POS,DELTA_POS)> results (fn_pt);
    
    //Store the list of nodes and unit cells in a usable format
    for(unsigned int i = 0; i < nodes.size(); i++){
        map<DELTA_POS, vector<int> >::iterator iter = results.find(displacements.at(i));
        if(iter == results.end()){
            vector<int> newList; newList.push_back(i);
            results.insert(pair<DELTA_POS, vector<int> > (displacements.at(i),newList));
        }
        else
            iter->second.push_back(i);
    }
    
    // Copy the unit cells visited and their corresponding node ids to intrinsic
    // data structures
    map<DELTA_POS, vector<int> >::iterator rIter = results.begin();
    while(rIter != results.end()){
        unitCells.push_back(rIter->first);
        ucNodes.push_back(rIter->second);
        rIter++;
    }
}

/** Write the commands necessary to draw the CHANNEL in ZeoVis
 *  to the provided output stream. */
void CHANNEL::writeToVMD(int n, fstream &output){
    if(!output.is_open()){
        cerr << "Error: File stream needed to print channel information was not open." << "\n"
        << "Exiting ..." << "\n";
        exit(1);
    }
    else{
        output << "set channels(" << n << ") {" << "\n"
        << "{color $channelColors(" << n << ")}" << "\n";
        
        // Draw the components located in each unit cell
        for(unsigned int i = 0; i < unitCells.size(); i++){
            vector<int> nodeIDs = ucNodes.at(i);
            DELTA_POS disp = unitCells.at(i);
            
            // Iterate over all nodes in the unit cell
            for(unsigned int j = 0; j < nodeIDs.size(); j++){
                DIJKSTRA_NODE curNode = nodes.at(nodeIDs.at(j));
                
                // Find node coordinates
                double xCoord = curNode.x + v_a.x*disp.x + v_b.x*disp.y + v_c.x*disp.z;
                double yCoord = curNode.y + v_a.y*disp.x + v_b.y*disp.y + v_c.y*disp.z;
                double zCoord = curNode.z + v_a.z*disp.x + v_b.z*disp.y + v_c.z*disp.z;
                
                // Command used to draw node
                output << "{sphere {" << xCoord <<  " " << yCoord << " " << zCoord
                << "} radius $nodeRadii(" << nodeIDs.at(j) <<") resolution $sphere_resolution}"
                << "\n";
                
                // Iterate over all connections stemming from the current node
                for(unsigned int k = 0; k < curNode.connections.size(); k++){
                    CONN curConn = curNode.connections.at(k);
                    DIJKSTRA_NODE otherNode = nodes.at(curConn.to);
                    int dx = disp.x + curConn.deltaPos.x;
                    int dy = disp.y + curConn.deltaPos.y;
                    int dz = disp.z + curConn.deltaPos.z;
                    
                    // Find end of edge coordinates
                    double newX = otherNode.x + v_a.x*dx + v_b.x*dy + v_c.x*dz;
                    double newY = otherNode.y + v_a.y*dx + v_b.y*dy + v_c.y*dz;
                    double newZ = otherNode.z + v_a.z*dx + v_b.z*dy + v_c.z*dz;
                    
                    // Command used to draw edge
                    output << "{line {" << xCoord << " " << yCoord << " " << zCoord << "} {"
                    << newX << " " << newY << " " << newZ<< "}}" <<  "\n";
                }
            }
        }
        output << "}" << "\n";
    }
}


/** Write the commands necessary to draw the CHANNEL in ZeoVis
 *  to the provided output stream. Includes a type because features and segments
 *  are drawn using the same command.*/
void CHANNEL::writeToVMD(string type, int n, fstream &output){
    if(!output.is_open()){
        cerr << "Error: File stream needed to print" << type << " information was not open." << "\n"
        << "Exiting ..." << "\n";
        exit(1);
    }
    else{
        output << "set " << type << "s(" << n << ") {" << "\n"
        << "{color $" << type << "Colors(" << n << ")}" << "\n";
        
        // Draw the components located in each unit cell
        for(unsigned int i = 0; i < unitCells.size(); i++){
            vector<int> nodeIDs = ucNodes.at(i);
            DELTA_POS disp = unitCells.at(i);
            
            // Iterate over all nodes in the unit cell
            for(unsigned int j = 0; j < nodeIDs.size(); j++){
                DIJKSTRA_NODE curNode = nodes.at(nodeIDs.at(j));
                
                // Find node coordinates
                double xCoord = curNode.x + v_a.x*disp.x + v_b.x*disp.y + v_c.x*disp.z;
                double yCoord = curNode.y + v_a.y*disp.x + v_b.y*disp.y + v_c.y*disp.z;
                double zCoord = curNode.z + v_a.z*disp.x + v_b.z*disp.y + v_c.z*disp.z;
                
                // Command used to draw node
                output << "{sphere {" << xCoord <<  " " << yCoord << " " << zCoord
                << "} radius $nodeRadii(" << nodeIDs.at(j) <<") resolution $sphere_resolution}" << "\n";
                
                // Iterate over all connections stemming from the current node
                for(unsigned int k = 0; k < curNode.connections.size(); k++){
                    CONN curConn = curNode.connections.at(k);
                    DIJKSTRA_NODE otherNode = nodes.at(curConn.to);
                    int dx = disp.x + curConn.deltaPos.x;
                    int dy = disp.y + curConn.deltaPos.y;
                    int dz = disp.z + curConn.deltaPos.z;
                    
                    // Find end of edge coordinates
                    double newX = otherNode.x + v_a.x*dx + v_b.x*dy + v_c.x*dz;
                    double newY = otherNode.y + v_a.y*dx + v_b.y*dy + v_c.y*dz;
                    double newZ = otherNode.z + v_a.z*dx + v_b.z*dy + v_c.z*dz;
                    
                    // Command used to draw edge
                    output << "{line {" << xCoord << " " << yCoord << " " << zCoord << "} {"
                    << newX << " " << newY << " " << newZ<< "}}" <<  "\n";
                }
            }
        }
        output << "}" << "\n";
    }
}


/* Create a channel from a pore */
CHANNEL::CHANNEL(PORE *p){
    nodes = p->nodes;
    connections = p->connections;
    unitCells = p->unitCells;
    ucNodes = p->ucNodes;
    dimensionality = p->dimensionality;
    basis[0][0]=p->basis[0][0];basis[0][1]=p->basis[0][1];basis[0][2]=p->basis[0][2];
    basis[1][0]=p->basis[1][0];basis[1][1]=p->basis[1][1];basis[1][2]=p->basis[1][2];
    basis[2][0]=p->basis[2][0];basis[2][1]=p->basis[2][1];basis[2][2]=p->basis[2][2];
}


/** Nodes within the provided DIJKSTRA_NETWORK can be classified as either
 *  accessible or inaccessible, where sets of accessible  nodes constitute a
 *  CHANNEL. As a result, this function identifies the CHANNELs that exist for
 *  the provided particle diameter and stores them using the provided pointer
 *  to a vector of channels. In addition, the pointer to the vector of bools is
 *  used to a store bool for each VORONOI_NODE, where infoStorage[i] is true
 * iff node #i is accessible. */
/*  WARNING: this function assumes that at some point before VORONOI_NETWORK was pruned
 *  with prune() function. This is to set "active" flag on nodes larger than probe radius */
void CHANNEL::findChannels(DIJKSTRA_NETWORK *dnet, vector<bool> *infoStorage,
                           vector<CHANNEL> *channels)
{
    
    vector <PORE> pores;
    findChannelsAndPockets(dnet, infoStorage, &pores);
    for(unsigned int i = 0; i<pores.size(); i++)
    {
        if(pores[i].dimensionality>0) channels->push_back(CHANNEL(&pores[i]));
    };
    pores.clear();
}
/* below is an old version before PORE class was introducedd */
/*
 void CHANNEL::findChannels(DIJKSTRA_NETWORK *dnet, vector<bool> *infoStorage,
 vector<CHANNEL> *channels)
 {
 //Define three access types
 const int ACCESSIBLE   =  1;
 const int INACCESSIBLE =  0;
 const int UNKNOWN      = -1;
 
 vector<DELTA_POS> directions = vector<DELTA_POS> ();
 
 // Initialize the status of each node as UNKNOWN
 vector<int> accessStatuses = vector<int> (dnet->nodes.size(), UNKNOWN);
 unsigned int nodeIndex = 0;
 
 // Iterate over all nodes
 while(nodeIndex != accessStatuses.size()){
 //Skip node if access status already determined
 if(accessStatuses.at(nodeIndex) != UNKNOWN){
 nodeIndex++;
 continue;
 }
 
 // Start out with 0-dimensional channel w/o any basis vectors
 int dim = 0;
 int basis [3][3] = {{0,0,0},
 {0,0,0},
 {0,0,0}};
 
 // Place starting node on stack with (0,0,0) displacement
 int accessStatus = UNKNOWN;
 DELTA_POS displacement = DELTA_POS(0,0,0);
 map<int,DELTA_POS> visitedNodeDisplacement;
 vector<pair<int,DELTA_POS> > stack;
 stack.push_back(pair<int,DELTA_POS> (nodeIndex,displacement));
 visitedNodeDisplacement.insert(pair<int,DELTA_POS>(nodeIndex, displacement));
 
 while(stack.size() != 0){
 // Remove top-most node
 pair<int,DELTA_POS> nodeInfo = stack.back();
 DIJKSTRA_NODE currentNode = dnet->nodes.at(nodeInfo.first);
 stack.pop_back();
 
 // Follow all edges leading to adjoining nodes
 vector<CONN>::iterator connIter = currentNode.connections.begin();
 while(connIter != currentNode.connections.end()){
 int to = connIter->to;
 DELTA_POS newDisplacement = nodeInfo.second + connIter->deltaPos;
 int localAccessStatus = accessStatuses.at(to);
 
 if(localAccessStatus == UNKNOWN){
 map<int,DELTA_POS>::iterator visitedNode = visitedNodeDisplacement.find(to);
 if(visitedNode != visitedNodeDisplacement.end()){
 if(visitedNode->second.equals(newDisplacement)){
 // Circling back to previous node
 // Do nothing
 }
 else{
 // Nodes are ACCESSIBLE b/c same node visited in different unit cell
 accessStatus = ACCESSIBLE;
 DELTA_POS direction = newDisplacement - visitedNode->second;
 
 
 // ******* Start of Michael's code ******
 if(dim == 0){
 basis[0][0] = direction.x;
 basis[1][0] = direction.y;
 basis[2][0] = direction.z;
 dim++;
 }
 else if (dim == 1){
 double v;
 
 if (basis[0][0] != 0) v = 1.0*direction.x/basis[0][0];
 else if (basis[1][0] != 0) v = 1.0*direction.y/basis[1][0];
 else if (basis[2][0] != 0) v = 1.0*direction.z/basis[2][0];
 else {
 cerr << "Error: Channel basis vector is zero vector. Exiting..." << "\n";
 exit(1);
 }
 
 if(!(basis[0][0]*v == direction.x && basis[1][0]*v == direction.y &&
 basis[2][0]*v == direction.z)){
 basis[0][1] = direction.x;
 basis[1][1] = direction.y;
 basis[2][1] = direction.z;
 dim++;
 }
 }
 else if (dim == 2){
 basis[0][2] = direction.x;
 basis[1][2] = direction.y;
 basis[2][2] = direction.z;
 if(abs(calcDeterminant(basis)) < 1e-8){
 basis[0][2] = basis[1][2] = basis[2][2] = 0;
 }
 else
 dim++;
 }
 // ****** End of Michael's code ******
 }
 }
 else{
 // Node being visited for the first time
 visitedNodeDisplacement.insert(pair<int,DELTA_POS>(to, newDisplacement));
 stack.push_back(pair<int,DELTA_POS> (to, newDisplacement));
 }
 }
 else{
 // The status of a region should be determined
 // with all of its neighbors. Otherwise, all possible nodes
 // were not explored.
 cerr << "Error: Illogical result  when attempting to identify channels." << "\n";
 cerr << "Please contact the source code provider with your program input. " << "\n";
 cerr << "Exiting ..." << "\n";
 exit(1);
 }
 connIter++;
 }
 }
 
 // All connected nodes are inaccessible if their status is still UNKNOWN
 if((stack.size() == 0) && (accessStatus == UNKNOWN))
 accessStatus = INACCESSIBLE;
 
 // Record the access status for each node in the cycle and create
 // a list of the node ids for channel identification
 map<int,DELTA_POS>::iterator resultIter = visitedNodeDisplacement.begin();
 vector<int> listOfIDs = vector<int> ();
 while(resultIter != visitedNodeDisplacement.end()){
 
 int nodeID = resultIter->first;
 if(accessStatuses.at(nodeID) != UNKNOWN){
 // Each node should only have its accessibility determined once
 cerr << "Error: Accessibility of node was determined more than once." << "\n";
 cerr << "Please contact the source code provider with your program input. " << "\n";
 cerr << "Exiting ..." << "\n";
 exit(1);
 }
 else{
 // Store node's status
 accessStatuses.at(nodeID) = accessStatus;
 listOfIDs.push_back(nodeID);
 }
 resultIter++;
 }
 nodeIndex++;
 
 // Create CHANNEL from ACCESSIBLE nodes
 if(accessStatus == ACCESSIBLE){
 channels->push_back(CHANNEL(listOfIDs, dnet, dim, basis));
 //output << dim << " ";
 }
 }
 // Store information in provided vector
 infoStorage -> resize(accessStatuses.size());
 for(unsigned int i = 0; i < accessStatuses.size(); i++){
 infoStorage -> at(i) = (accessStatuses.at(i) == ACCESSIBLE);
 }
 }
 */

/*  Voronoi nodes within the provided VORONOI_NETWORK can be classified as either
 *  accessible or inaccessible, where sets of accessible  nodes constitute a
 *  CHANNEL. As a result, this function identifies the CHANNELs that exist for
 *  the provided particle diameter and  stores them using the provided pointer
 *  to a vector of channels. In addition, the pointer to the vector of bools is
 *  used to a store bool for each VORONOI_NODE, where infoStorage[i] is true
 *  iff node #i is accessible. */
void CHANNEL::findChannels(VORONOI_NETWORK *vornet, double minRadius,
                           vector<bool> *infoStorage, vector<CHANNEL> *channels)
{
    vector <PORE> pores;
    findChannelsAndPockets(vornet, minRadius, infoStorage, &pores);
    for(unsigned int i = 0; i<pores.size(); i++)
    {
        if(pores[i].dimensionality>0) channels->push_back(CHANNEL(&pores[i]));
    };
    pores.clear();
}

/* below is an earlier version before the PORE class was introduced */
/*
 void CHANNEL::findChannels(VORONOI_NETWORK *vornet, double minRadius,
 vector<bool> *infoStorage, vector<CHANNEL> *channels)
 {
 //Remove edges that don't allow the provided particle diameter to
 //pass freely
 //VORONOI_NETWORK newNetwork;
 //pruneVoronoiNetwork(vornet, &newNetwork, minRadius);
 VORONOI_NETWORK newNetwork = vornet->prune(minRadius);
 
 //Build graph data structure
 DIJKSTRA_NETWORK dnet;
 DIJKSTRA_NETWORK::buildDijkstraNetwork(&newNetwork, &dnet);
 findChannels(&dnet, infoStorage, channels);
 }
 */


/** Stores the ids of all atoms that bound this channel using the provided vector
 *  reference. An atom is considered to bound a channel if a node in the channel
 *  is a member of the atom's Voronoi cell. */
void CHANNEL::findBoundingAtoms(ATOM_NETWORK *atmnet, vector<BASIC_VCELL> &vcells,
                                vector<int> &atomIDs)
{
    atomIDs.clear();
    
    for(unsigned int i = 0; i < vcells.size(); i++){
        BASIC_VCELL cell = vcells[i];
        for(int j = 0; j < cell.getNumNodes(); j++){
            if(idMappings.find(cell.getNodeID(j)) != idMappings.end()){
                atomIDs.push_back(i);
                break;
            }
        }
    }
}

/* Converts the current channel into DIJKSTRA_NETWORK so it can be analyzed (to
 * identify free/max sphere diameters */
/* this function is currently unused and may be removed in the future */
void PORE::buildDijkstraNetwork(DIJKSTRA_NETWORK *dnet){
    //  vector<VOR_NODE> ::iterator niter = vornet->nodes.begin();
    //    int i = 0;
    dnet->nodes.clear();
    
    dnet->nodes=nodes;
    
    // Add copies of all nodes to the network
	//   while(niter != vornet->nodes.end()){
	//       DIJKSTRA_NODE node = DIJKSTRA_NODE(i, niter->x, niter->y, niter->z, niter->rad_stat_sphere);
	//       i++;
	//       niter++;
	//       dnet->nodes.push_back(node);
	//   }
	//
	//  // For each edge, store it in the DIJKSTRA node's list of connections that the connection
	//  // stems from
	//  vector<VOR_EDGE> ::iterator eiter = vornet->edges.begin();
	//  while(eiter != vornet->edges.end()){
	//      DELTA_POS pos = DELTA_POS(eiter->delta_uc_x, eiter->delta_uc_y, eiter->delta_uc_z);
	//      CONN conn = CONN(eiter->from, eiter->to, eiter->length, eiter->rad_moving_sphere,pos);
	//      dnet->nodes.at(conn.from).connections.push_back(conn);
	//      eiter++;
	//  }
	//
    // Copy the unit cell vectors into the DIJKSTRA_NETWORK
    dnet->v_a = v_a;
    dnet->v_b = v_b;
    dnet->v_c = v_c;
}


/* Returns the largest free sphere diameter for the current channel */
/* It is a rather inefficient algorithm that calculates the largest
 * free sphere for each Voronoi node */
pair<double, pair<double,double> > CHANNEL::findFreeIncludedSphereDiameter(){
    
    pair<double, pair<double,double> > maxdidfdif;
    maxdidfdif.first=0;
    maxdidfdif.second.first=0;
    maxdidfdif.second.second=0;
    
    
    for(int i=0;i<nodes.size();i++)
    {
        // loop over all the nodes in the current channel
        if(i==0) {
            maxdidfdif = findFreeIncludedSphereDiameterforNode(i, maxdidfdif);
        }
        else
        {
            pair<double, pair<double,double> > didfdif = findFreeIncludedSphereDiameterforNode(i, maxdidfdif);
            if(didfdif.second.first > maxdidfdif.second.first) maxdidfdif.second = didfdif.second;
            if(didfdif.first > maxdidfdif.first) maxdidfdif.first = didfdif.first;
        };
    };
    
    return maxdidfdif;
}

/* Pointer that connects HEAP class with CHANNEL class */
vector<DIJKSTRA_NODE> *compareConnections_ptr;

/* Function for sorting heap in findFreeSphereDiameter */
bool compareConnections(pair<int,int> C1, pair<int,int> C2){
    
    return (compareConnections_ptr->at(C1.first).connections.at(C1.second).
            max_radius<compareConnections_ptr->at(C2.first).connections.
            at(C2.second).max_radius);
    
}



/* Return the largest free and included sphere diameters (df and dif)  starting
 * from a partiular node  of the current channel
 * the function also accepts point to a vector that flags nodes accessed nodes during
 * execution of the function
 * It uses a pair<pair> object that storage the current 'record' Di/df/dif
 *  This is to speed-up by early discart of paths with more restriction than the current path
 */
pair <double, pair<double,double> > CHANNEL::findFreeIncludedSphereDiameterforNode(int nodeid, pair<double, pair<double,double> > maxdidfdif){
    //Define three access types
    const int ACCESSED   =  1;
    const int UNKNOWN      = -1;
    
    // Make a local copy of nodes vector as it will be modified
    vector<DIJKSTRA_NODE> lnodes = nodes;
    
    // Initialize the unit cell position of each node
    vector<DELTA_POS> directions = vector<DELTA_POS> (lnodes.size(), DELTA_POS(0,0,0));
    
    // Initialize the status of each node as UNKNOWN
    vector<int> accessStatuses = vector<int> (lnodes.size(), UNKNOWN);
    unsigned int nodeIndex = 0;
    
    double maxdi,di,df;
    
    
    // In case of early termination, be ready to return the current node Di
    maxdidfdif.first = 2.0 * lnodes.at(nodeid).max_radius; // .first is use to store the current Di for the current (nodeid) node
    
    // Test if it is worth to run the funtion (e.g. starting node's radius is smaller than the current Df)
    if(lnodes.at(nodeid).max_radius < maxdidfdif.second.first/2.0)
    {
        return maxdidfdif; // if not, just return the current best numbers (maxdidfdif)
    };
    
    // Start analysis for node of id = nodeid
    
    nodeIndex = nodeid;
    
    // accept the initial node
    maxdi = lnodes.at(nodeIndex).max_radius; // (it is not really the max but will be analyzed to get max)
    di=lnodes.at(nodeIndex).max_radius;
    df=lnodes.at(nodeIndex).max_radius;
    
    accessStatuses[nodeIndex] = ACCESSED;
    
    directions[nodeIndex] = DELTA_POS(0,0,0);
    
    compareConnections_ptr = &lnodes;
    HEAP< pair<int,int> > stack(compareConnections);
    
    // loop over all connections and add them to stack
    for(unsigned int i=0; i<lnodes[nodeIndex].connections.size(); i++)
    {
        if(lnodes[nodeIndex].connections.at(i).max_radius > lnodes.at(lnodes[nodeIndex].connections.at(i).to).max_radius)
            lnodes[nodeIndex].connections.at(i).max_radius = lnodes.at(lnodes[nodeIndex].connections.at(i).to).max_radius;
        
        stack.insert(pair<int,int>(nodeIndex,i));
    };
    
    stack.reHeapify();
    
    bool flag=true;
    
    while(stack.size() != 0&&flag){
        //loop with Dijkstra-like propagation
        
        // select the path with largest opening
        pair<int,int> best = stack.pop();
        
        int best_node = lnodes[best.first].connections[best.second].to;
        
        if(accessStatuses[best_node]==UNKNOWN)
        {
            // going to accept the best_node and update di and df
            // update df and di
            if(lnodes[best.first].connections.at(best.second).max_radius<df) df=lnodes[best.first].connections.at(best.second).max_radius;
            
            if(lnodes.at(best_node).max_radius > di) di = lnodes.at(best_node).max_radius;
            
            // termine early if the current Df is smaller than the current record
            
            if(df * 2.0 < maxdidfdif.second.first)
            {
                return maxdidfdif; // if not, just return the current best numbers (maxdidfdif) 
            };
            
            
            // check new connections, add to heap or terminane
            
            accessStatuses[best_node]=ACCESSED; 
            directions[best_node]=directions[best.first];
            directions[best_node]=directions[best_node] + (lnodes[best.first].connections.at(best.second).deltaPos);
            
            // going to add new connections to the heap
            for(unsigned int i=0; i<lnodes[best_node].connections.size(); i++)
            {
                if(lnodes[best_node].connections.at(i).max_radius > lnodes.at(lnodes[best_node].connections.at(i).to).max_radius)
                    lnodes[best_node].connections.at(i).max_radius = lnodes.at(lnodes[best_node].connections.at(i).to).max_radius;
                
                if(lnodes[best_node].connections[i].to!=best.first)
                    stack.insert(pair<int,int>(best_node,i));
            };	
            stack.reHeapify(); 
        }
        
        else if(accessStatuses[best_node]==ACCESSED){
            //selected connection leads to a node that is on stact
            //need to check in which unit cell it is
            
            
            if(!(directions[best_node]).equals(directions[best.first] + lnodes[best.first].connections.at(best.second).deltaPos))
            {
                // if nodes are in different unit cells, a loop is found, the resulting df can be calculated
                if(df>nodes[best.first].connections.at(best.second).max_radius) df=nodes[best.first].connections.at(best.second).max_radius;
                flag=false;
                break;
            };
            // otherwise do nothing
            // need to make sure heap is correct
        }
        
        else{
            //
            // need this section ?
            //
        };
        
    };
    
    // change radii to diameter
    df=df*2.0;
    di=di*2.0;
    maxdi=maxdi*2.0;
    
    pair <double, pair<double,double> > resultdidfdif;
    resultdidfdif.first=maxdi;
    resultdidfdif.second.first=df;
    resultdidfdif.second.second=di;
    
    return resultdidfdif;
}


/* Get the largest included sphere diamter for pocket */

double PORE::getIncludedSphereDiameter(){
    
    double maxdi;
    
    for(unsigned int i=0; i<nodes.size(); i++)
    {
        if(i==0)
        {
            maxdi = nodes.at(i).max_radius;
        }
        else{
            if(nodes.at(i).max_radius > maxdi) maxdi = nodes.at(i).max_radius;
        };
    };
    
    maxdi=maxdi*2.0;
    
    return maxdi;
    
}



/* Return the largest free and included sphere diameters (df and dif) for a path between two nodes within the current
 * pore. Return -1 in nodes are not connected.
 * Node1 and Node2 are indexes provied w.r.t. the original network. Use idMapping inside this function to obtain local indexes 
 * for the considered PORE
 */
 
pair <double, double> PORE::getFreeIncludedSphereDiameterforNodePair(int node1, int node2){
    //Define three access types
    const int ACCESSED   =  1;
    const int UNKNOWN      = -1;

    // flag that is used to terminate the search
    
    bool node2FoundFlag = false;

    // Make a local copy of nodes vector as it will be modified
    vector<DIJKSTRA_NODE> lnodes = nodes;
    
    // Initialize the unit cell position of each node
    vector<DELTA_POS> directions = vector<DELTA_POS> (lnodes.size(), DELTA_POS(0,0,0));
    
    // Initialize the status of each node as UNKNOWN
    vector<int> accessStatuses = vector<int> (lnodes.size(), UNKNOWN);
    unsigned int nodeIndex = 0;
    
    double di,df;

    if(node1 == node2 || node1 < 0 || node2 < 0){
      cerr << "Requesting Df/Dif conculation for nodes of incorrect IDs.\n" << "\n";
      abort();
      };

    // node1 and node2 are provided

    node1 = idMappings.find(node1)->second;
    node2 = idMappings.find(node2)->second;
    
    // Start analysis for node of id = nodeid
    
    nodeIndex = node1; // node1 is w.r.t pore node index (use PORE reverseMapping to retrieve original node
    
    // accept the initial node
    di=lnodes.at(nodeIndex).max_radius;
    df=lnodes.at(nodeIndex).max_radius;
    
    accessStatuses[nodeIndex] = ACCESSED;
    
    directions[nodeIndex] = DELTA_POS(0,0,0);
    
    compareConnections_ptr = &lnodes;
    HEAP< pair<int,int> > stack(compareConnections);
    
    // loop over all connections and add them to stack
    for(unsigned int i=0; i<lnodes[nodeIndex].connections.size(); i++)
    {
        if(lnodes[nodeIndex].connections.at(i).max_radius > lnodes.at(lnodes[nodeIndex].connections.at(i).to).max_radius)
            lnodes[nodeIndex].connections.at(i).max_radius = lnodes.at(lnodes[nodeIndex].connections.at(i).to).max_radius;
        
        stack.insert(pair<int,int>(nodeIndex,i));
    };
    
    stack.reHeapify();
    
    bool flag=true;
    
    while(stack.size() != 0&&flag){
        //loop with Dijkstra-like propagation
        
        // select the path with largest opening
        pair<int,int> best = stack.pop();
        
        int best_node = lnodes[best.first].connections[best.second].to;
        
        if(accessStatuses[best_node]==UNKNOWN)
        {
            // going to accept the best_node and update di and df
            // update df and di
            if(lnodes[best.first].connections.at(best.second).max_radius<df) df=lnodes[best.first].connections.at(best.second).max_radius;
            
            if(lnodes.at(best_node).max_radius > di) di = lnodes.at(best_node).max_radius;

            // in the desired node is reached, and the loop 
            
            if(best_node == node2)
              {
              node2FoundFlag = true; 
              flag = true; // stop flag set
              break;
              };
            
            // otherwise continue

            // check new connections, add to heap or terminane
            
            accessStatuses[best_node]=ACCESSED; 
            directions[best_node]=directions[best.first];
            directions[best_node]=directions[best_node] + (lnodes[best.first].connections.at(best.second).deltaPos);
            
            // going to add new connections to the heap
            for(unsigned int i=0; i<lnodes[best_node].connections.size(); i++)
            {
                if(lnodes[best_node].connections.at(i).max_radius > lnodes.at(lnodes[best_node].connections.at(i).to).max_radius)
                    lnodes[best_node].connections.at(i).max_radius = lnodes.at(lnodes[best_node].connections.at(i).to).max_radius;
                
                if(lnodes[best_node].connections[i].to!=best.first)
                    stack.insert(pair<int,int>(best_node,i));
            };	
            stack.reHeapify(); 
        }
        
        else if(accessStatuses[best_node]==ACCESSED){
            //selected connection leads to a node that is on stact
            //need to check in which unit cell it is
            
            
            if(!(directions[best_node]).equals(directions[best.first] + lnodes[best.first].connections.at(best.second).deltaPos))
            {
                // if nodes are in different unit cells (a loop is found); just update df and continue 
                if(df>nodes[best.first].connections.at(best.second).max_radius) df=nodes[best.first].connections.at(best.second).max_radius;
            };
            // otherwise do nothing
            // need to make sure heap is correct (which should be true as nothing was added)
        }
        
        else{
            //
            // need this section ?
            //
        };
        
    };
    
    // change radii to diameter
    df=df*2.0;
    di=di*2.0;
    
    if(node2FoundFlag == false)
      {
      df = -1;
      di = -1;
      };

    pair<double,double>  resultdfdif;
    resultdfdif.first=df;
    resultdfdif.second=di;
    
    return resultdfdif;
}




/* Get information about the pocket: di, coordinates and encapsulating sphere radius */
void PORE::getSimplifiedPocketInfo(ATOM_NETWORK *atmNet, std::vector <double> *pocketInfo){
    
    double maxdi;
    double maxid = 0;
    pocketInfo->clear();
    
    for(unsigned int i=0; i<nodes.size(); i++)
    {
        if(i==0)
        {
            maxdi = nodes.at(i).max_radius;
        }
        else{
            if(nodes.at(i).max_radius > maxdi) 
            {
                maxdi = nodes.at(i).max_radius;
                maxid = i;
            };
        };
    };
    
    maxdi=maxdi*2.0;
    
    pocketInfo->push_back(maxdi);
    
    Point pt = atmNet->xyz_to_abc(nodes.at(maxid).x, nodes.at(maxid).y, nodes.at(maxid).z);
    pt = atmNet->shiftABCInUC(pt);
    pocketInfo->push_back(pt[0]); pocketInfo->push_back(pt[1]); pocketInfo->push_back(pt[2]);
    
    double pocketR; // radii of pocket (sphere that encapsulates all nodes in the pocket)
    
    pocketR = maxdi*0.5;
    for(unsigned int i=0; i<nodes.size(); i++)
    {
        double dist = calcEuclideanDistance(nodes.at(maxid).x, nodes.at(maxid).y, nodes.at(maxid).z, nodes.at(i).x, nodes.at(i).y, nodes.at(i).z);
        dist = dist + nodes.at(i).max_radius;
        
        if(dist > pocketR) pocketR = dist;
        
    };
    
    pocketInfo->push_back(pocketR);
    
}



/* Functions gets restrictions between provided sements of Voronoi network within the current PORE */
/* It takes the original number of segments (to be specified by user) and a list of original Voronoi nodes and their segment assigments 
   (that is w.r.t. to pre-channel detection) than analyzed segment connectivity within the currnet pore.
*/
//void PORE::getRestrictingDiameters(int nSegments, vector<int> vorNetID, vector< vector<double> > *PLDtable){
void PORE::getRestrictingDiameters(int nSegments, vector<int> vorNetID, vector< vector<double> > *PLDtable, vector< vector< pair<int,int> > > *PLDEdgeTable, 
                                  vector<double> *segmentDi, vector<int> *segmentDiNodeID, vector <double> *segmentDiFinal, vector<int> *segmentDiFinalNodeID){
    //Define three access types
    const int ACCESSED   =  1;
    const int UNKNOWN      = -1;

    // flag that is used to terminate the search
    
    bool node2FoundFlag = false;

    // Make a local copy of nodes vector as it will be modified
    vector<DIJKSTRA_NODE> lnodes = nodes;
    
    // Initialize the unit cell position of each node
    vector<DELTA_POS> directions = vector<DELTA_POS> (lnodes.size(), DELTA_POS(0,0,0));
    
    // Initialize the status of each node as UNKNOWN
    vector<int> accessStatuses = vector<int> (lnodes.size(), UNKNOWN);
    unsigned int nodeIndex = 0;
    
    double di,df;

    if(nSegments < 1){
      cerr << "Number of semgents lower than 1. This function should not be called.\n" << "\n";
      abort();
      };

//debug
//for(unsigned int i=0; i<vorNetID.size();i++) cout << vorNetID[i] << "  ";


    // Need to translate IDs from the original Voronoi network to PORE ids 

    vector <int> PoreNodesSegmentIDs; // will store segment ID for each node in the current PORE

    PoreNodesSegmentIDs.resize(nodes.size(), -1);

    int nSegPresent = 0;
    for(unsigned int i=0; i<nodes.size(); i++)
       {
       int VnodeID = reverseIDMappings.find(i)->second;
       if(vorNetID[VnodeID] >= 0)
         {
         PoreNodesSegmentIDs[i] = vorNetID[VnodeID]; 
         nSegPresent++;
         };

       };

    cout << "Current PORE contains " << nSegPresent << " nodes in predefined segments. Proceeding with Restricting Diameter Calcs\n";

    // Start analysis by accepting all nodes already assigned to segments (segment seed nodes)
    // Store information about largest node pre-assigned to each segment

    for(unsigned int i=0; i<nodes.size(); i++)
       {
       if(PoreNodesSegmentIDs[i] > -1) 
          {
          accessStatuses[i] = ACCESSED;   
          directions[i] = DELTA_POS(0,0,0); // not sure about this line
          
          if(segmentDi->at(PoreNodesSegmentIDs[i]) < nodes.at(i).max_radius*2) 
             {
             segmentDi->at(PoreNodesSegmentIDs[i]) = nodes.at(i).max_radius*2;
             segmentDiNodeID->at(PoreNodesSegmentIDs[i]) = reverseIDMappings.find(i)->second;
             segmentDiFinal->at(PoreNodesSegmentIDs[i]) = segmentDi->at(PoreNodesSegmentIDs[i]);
             };
          };
       };


/* BUG BUG line (1ln) below, it probably clears segmentDiFinal too many times, moving to the above loop */
//    for(int i=0; i<nSegments; i++) segmentDiFinal->at(i) = segmentDi->at(i); // copying current segments into Final
    
/* old stuff to be removed
    // Start analysis for node of id = nodeid
    
    nodeIndex = node1; // node1 is w.r.t pore node index (use PORE reverseMapping to retrieve original node
    
    // accept the initial node
    di=lnodes.at(nodeIndex).max_radius;
    df=lnodes.at(nodeIndex).max_radius;
    
    accessStatuses[nodeIndex] = ACCESSED;
    
    directions[nodeIndex] = DELTA_POS(0,0,0);

 ends old stuff */
    
    compareConnections_ptr = &lnodes;
    HEAP< pair<int,int> > stack(compareConnections);
    
    // loop over all connections and add them to stack
    for(unsigned int n=0;  n<nodes.size(); n++)
      {
      if(accessStatuses[n] == ACCESSED)
        {
        for(unsigned int i=0; i<lnodes[n].connections.size(); i++)
          {
          if(lnodes[n].connections.at(i).max_radius > lnodes.at(lnodes[n].connections.at(i).to).max_radius)
              lnodes[n].connections.at(i).max_radius = lnodes.at(lnodes[n].connections.at(i).to).max_radius;
        
          stack.insert(pair<int,int>(n,i));
          };
        };
      };
    
    stack.reHeapify();

//
// Note: in the following blocks, we probably dont need directions[]
// to be cleaned in the future

    
    while(stack.size() != 0){
        //loop with Dijkstra-like propagation
        
        // select the path with largest opening
        pair<int,int> best = stack.pop();
        
        int best_node = lnodes[best.first].connections[best.second].to;
        
        if(accessStatuses[best_node]==UNKNOWN)
        {
            // if the node (best_node) is not assigned to segment, assign it to the current segment (determined by node best.first)
            // check new connections, add to heap or terminane
            
            accessStatuses[best_node] = ACCESSED; 
            directions[best_node] = directions[best.first];
            directions[best_node] = directions[best_node] + (lnodes[best.first].connections.at(best.second).deltaPos);
            PoreNodesSegmentIDs[best_node] = PoreNodesSegmentIDs[best.first];
           
            if(nodes.at(best_node).max_radius * 2 > segmentDiFinal->at(PoreNodesSegmentIDs[best.first])) 
               {   
               segmentDiFinal->at(PoreNodesSegmentIDs[best.first]) = nodes.at(best_node).max_radius * 2; 
               segmentDiFinalNodeID->at(PoreNodesSegmentIDs[best.first]) = reverseIDMappings.find(best_node)->second;
               };

            // going to add new connections to the heap
            for(unsigned int i=0; i<lnodes[best_node].connections.size(); i++)
            {
                if(lnodes[best_node].connections.at(i).max_radius > lnodes.at(lnodes[best_node].connections.at(i).to).max_radius)
                    lnodes[best_node].connections.at(i).max_radius = lnodes.at(lnodes[best_node].connections.at(i).to).max_radius;
                
                if(lnodes[best_node].connections[i].to!=best.first)
                    stack.insert(pair<int,int>(best_node,i));
            };	
//            stack.reHeapify(); 
//            I think reHeapify is not needed as heap.insert ensures proper sorting
        }
        
        else if(accessStatuses[best_node]==ACCESSED){
            //selected connection leads to a node that is already assigned to a segment

            // if it is of different segment, note the destricting diameter
            if(PoreNodesSegmentIDs[best_node] != PoreNodesSegmentIDs[best.first])
              {
              // record the largest resticting diameter between two segments of the same PORE
              double r = lnodes[best.first].connections.at(best.second).max_radius;

//begin debug
//cout << "PLD table size: " << PLDtable->size() << "   " << PLDtable->at(PoreNodesSegmentIDs[best_node]).size() << "\n";
//cout << "PoreNodesSegmentIDs[best_node]= " << PoreNodesSegmentIDs[best_node] << "   PoreNodesSegmentIDs[best.first]= " << PoreNodesSegmentIDs[best.first] << "\n"; 
//              if(PLDtable->size() > PoreNodesSegmentIDs[best_node] + 1 || PLDtable->at(PoreNodesSegmentIDs[best_node]).size() > PoreNodesSegmentIDs[best.first] + 1)
//    {cout << "alarm...\n";
//};
//end debug
              double stored_d = PLDtable->at(PoreNodesSegmentIDs[best_node]).at(PoreNodesSegmentIDs[best.first]);
              if(2*r> stored_d)
                {
                PLDtable->at(PoreNodesSegmentIDs[best_node]).at(PoreNodesSegmentIDs[best.first]) = 2*r;
                PLDtable->at(PoreNodesSegmentIDs[best.first]).at(PoreNodesSegmentIDs[best_node]) = 2*r;
                pair<int,int> p(reverseIDMappings.find(best_node)->second, reverseIDMappings.find(best.first)->second); 
                PLDEdgeTable->at(PoreNodesSegmentIDs[best_node]).at(PoreNodesSegmentIDs[best.first]) = p;
                PLDEdgeTable->at(PoreNodesSegmentIDs[best.first]).at(PoreNodesSegmentIDs[best_node]) = p;
                }else{
                //cerr << "Trying to store smaller diameter than already stored. This shouldnt happen. Check the code\n";
                };              
              } else
              {
              // both nodes are of the same segment, do nothing

              };
            
           /* old stuff, probably to be removed  
            if(!(directions[best_node]).equals(directions[best.first] + lnodes[best.first].connections.at(best.second).deltaPos))
            {
                // if nodes are in different unit cells (a loop is found); just update df and continue 
                if(df>nodes[best.first].connections.at(best.second).max_radius) df=nodes[best.first].connections.at(best.second).max_radius;
            };
           */
        }
        
        else{
            //
            // need this section ?
            //
        };
        
    };
   
   for(int i=0; i<nSegments; i++)
     {
     if(segmentDi->at(i)!=segmentDiFinal->at(i)) cerr << "Segment Di(" << segmentDi->at(i) << ") is different than Segment Di Final (" << segmentDiFinal->at(i) << ") for segment " << i << ".\n";
     };  

} // ends getRestrictingDiameters()






