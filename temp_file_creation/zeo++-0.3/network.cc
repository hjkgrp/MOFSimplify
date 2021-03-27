// 
//
// Author   : Thomas F. Willems (LBL / UC Berkeley)
// Email    : tfwillems@berkeley.edu
// Date     : July 14 2010
// Updated  : April 29, 2011
// Updated  : many times during 2012 (M. Haranczyk, R. Martin & the rest of the team)

/* Note: Code contained in performVoronoiDecomp() method is largely derived from
 * code originally written by Chris Rycroft, the author of the voro++
 * software package used to perform the Voronoi decomposition.
 */

#include "network.h"
#include "networkinfo.h"
#include "networkanalysis.h"
#include "graphstorage.h"
#include "channel.h"
#include "v_network.h"
#include "material.h"
#include "OMS.h"

using namespace std;
using namespace voro;

/* IMPORTANT - overwriting standard exit function - notifies user that exit was called before exiting */
void exit(int status) {
  printf("NOTICE: calling abort()\n");
  abort();
}

/* converts string to upper case */
std::string toUpperCase(const std::string & s)
{
    std::string ret(s.size(), char());
    for(unsigned int i = 0; i < s.size(); ++i)
        ret[i] = (s[i] <= 'z' && s[i] >= 'a') ? s[i]-('a'-'A') : s[i];
    return ret;
}


/** Decompose the provided network of atoms into a VORONOI_NETWORK that is stored using the provided pointer. If
    the option is specified, information about each VOR_CELL will also be stored using the provied pointer to
    a vector of VOR_CELL instances. The BASIC_VCELL information is stored regardless of the option specified.*/
void* performVoronoiDecomp(bool radial, ATOM_NETWORK *cell, VORONOI_NETWORK *vornet, vector<VOR_CELL> &cells, bool saveVorCells,
			   vector<BASIC_VCELL> &bvcells){
  int i,n;
  double bx,bxy,by,bxz,byz,bz;
  vector<int> atomShifts;
  
  // Set the box dimensions and number of particles
  bx  = cell->v_a.x;
  bxy = cell->v_b.x;
  by  = cell->v_b.y;
  bxz = cell->v_c.x;
  byz = cell->v_c.y;
  bz  = cell->v_c.z;
  n   = cell->numAtoms;
  
  // Print the box dimensions
  printf("Box dimensions:\n"
	 "  va=(%f 0 0)\n"
	 "  vb=(%f %f 0)\n"
	 "  vc=(%f %f %f)\n\n",bx,bxy,by,bxz,byz,bz);
  
  // Check that the input parameters make sense
  if(n<1) {
    char* sentence = new char[300];
    sprintf(sentence, "Error: Invalid number of particles provided for Voronoi decomposition (%d particles were read from file, which is <1)\nExiting ...\n", n);
    fputs(sentence, stderr);
    exit(1);
  }
  if(bx<tolerance||by<tolerance||bz<tolerance){
    fputs("Error: Invalid box dimensions calculated for Voronoi decomposition."
	  " Please check unit cell parameters.\nExiting ...\n",stderr);
    exit(1);
  }
  
  // Compute the internal grid size, aiming to make
  // the grid blocks square with around 6 particles
  // in each
  double ls=pow(n/(9*bx*by*bz),1.0/3.0);
  double nxf=bx*ls+1.1;
  double nyf=by*ls+1.1;
  double nzf=bz*ls+1.1;
  
  // Check the grid is not too huge, using floating point numbers to avoid
  // integer wrap-arounds
  if (nxf*nyf*nzf>max_regions) {
    fprintf(stderr,"voro++: Number of computational blocks exceeds the maximum allowed of %d\n"
		   "Either increase the particle length scale, or recompile with an increased\nmaximum.",max_regions);
    exit(1);
  }
  
  // Now that we are confident that the number of regions is reasonable,
  // create integer versions of them
  int nx=int(nxf);
  int ny=int(nyf);
  int nz=int(nzf);
  printf("Total particles = %d\n\nInternal grid size = (%d %d %d)\n\n",n,nx,ny,nz);

  //a temporary fix to the Voronoi volume check issue:
  //if the user allows, when the volume check fails the structure's atom coordinates can be randomly altered by a very small amount, to try and bypass the error
  int numAttemptsPermitted = 1;
  if(cell->allowAdjustCoordsAndCellFlag) numAttemptsPermitted++;
  for(int attempt=0; attempt<numAttemptsPermitted; attempt++) {
    
    if(radial) {
      puts("Using voro++ with radii for particles.");
      // Create a container with the geometry given above
      container_periodic_poly *rad_con;
      rad_con =  new container_periodic_poly (bx,bxy,by,bxz,byz,bz,nx,ny,nz,memory);
      
      // Read in the particles from the provided ATOM_NETWORK
      vector <ATOM> ::iterator iter = cell->atoms.begin();
      i = 0;
      int da, db, dc;
      while(iter != cell->atoms.end()){ 
        rad_con->put(i,iter->x,iter->y,iter->z,iter->radius, da, db, dc);
        atomShifts.push_back(da); atomShifts.push_back(db); atomShifts.push_back(dc);
        iter++;
        i++;
      }

      // Store the Voronoi Network for later use
      bool volume_correct = storeVoronoiNetwork(*rad_con, cell, vornet, bx, by, bz, bvcells, atomShifts, saveVorCells, cells);
      if(volume_correct) {
        vornet->v_a = cell->v_a;
        vornet->v_b = cell->v_b;
        vornet->v_c = cell->v_c;
        return rad_con;
      } else if(attempt==numAttemptsPermitted-1) {
        printf("Exiting...\n");
        exit(1);
      } else {
        cell->randomlyAdjustCoordsAndCell();
      }
    }
    else {
      puts("Using voro++ without radii for particles.");
      // Create a container with the geometry given above
      container_periodic *no_rad_con;
      no_rad_con = new container_periodic(bx,bxy,by,bxz,byz,bz,nx,ny,nz,memory);

      // Read in the particles from the provided ATOM_NETWORK
      vector <ATOM> ::iterator iter = cell->atoms.begin();
      i = 0;
      int da, db, dc;
      while(iter != cell->atoms.end()){
        no_rad_con->put(i,iter->x,iter->y,iter->z, da, db, dc);
        atomShifts.push_back(da); atomShifts.push_back(db); atomShifts.push_back(dc);
        iter++;
        i++;
      }

      // Store the Voronoi Network for later use
      bool volume_correct = storeVoronoiNetwork(*no_rad_con, cell, vornet, bx, by, bz, bvcells, atomShifts, saveVorCells, cells);
      if(volume_correct) {
        vornet->v_a = cell->v_a;
        vornet->v_b = cell->v_b;
        vornet->v_c = cell->v_c;
        return no_rad_con;
      } else if(attempt==numAttemptsPermitted-1) {
        printf("Exiting...\n");
        exit(1);
      } else {
        cell->randomlyAdjustCoordsAndCell();
      }
    }
  }
}

/** Decompose the provided network of atoms into a VORONOI_NETWORK that is stored using the provided pointer. If
    the option is specified, information about each VOR_CELL will also be stored using the provied pointer to
    a vector of VOR_CELL instances. The BASIC_VCELL information is stored regardless of the option specified. 
    This function is wrapper to above function and has no return type*/
bool performVoronoiDecomp(bool radial, ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet, vector<VOR_CELL> *cells, bool saveVorCells,
			   vector<BASIC_VCELL> *bvcells){

    container_periodic_poly *rad_con = NULL;
    container_periodic *no_rad_con = NULL;
    if (radial)
        rad_con = (container_periodic_poly *)performVoronoiDecomp(radial, atmnet, vornet, *cells, saveVorCells, *bvcells);
     else 
          no_rad_con = (container_periodic *)performVoronoiDecomp (radial, atmnet, vornet, *cells, saveVorCells, *bvcells); 

     delete rad_con;
     delete no_rad_con;
     return true;
}


void createAdvCell(voronoicell &cell, vector<double> coords, int *idMap, VOR_CELL &newCell) {
  int numFaces = cell.number_of_faces();
  vector<int> faceVertices;
  cell.face_vertices(faceVertices);
      
  int index = 0;
  for(int i = 0; i < numFaces; i++){
    vector<Point> faceCoords;
    vector<int>   faceIDs;
    int faceOrder = faceVertices[index];
    index++;
    for(int j = 0; j < faceOrder; j++){
      int verID = faceVertices[index];
      faceCoords.push_back(Point(coords[verID*3], coords[verID*3+1], coords[verID*3+2]));
      faceIDs.push_back(idMap[4*verID]);
      index++;
    }
    newCell.addFace(VOR_FACE(faceCoords, faceIDs));
  }
}

/** Perform the Voronoi decomposition using the voro++
    package and read the results into the provided VORONOI_NETWORK data
    structure. */
template<class c_option>
bool storeVoronoiNetwork(c_option &con, ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet, double bx, double by, double bz,
			 vector<BASIC_VCELL> &basCells, vector<int> &atomShifts, bool storeAdvCells, vector<VOR_CELL> &advCells) {
  voronoi_network vnet (con, VOR_NODE_MERGE_THRESHOLD); // data structure defined in voro++ which is not to be confused with VORONOI_NETWORK
  int id;
  double vvol=0,x,y,z,r;
  voronoicell c(con);

  puts("Performing Voronoi decomposition.");

  basCells.clear(); advCells.clear();
  basCells.resize(atmnet->numAtoms, BASIC_VCELL());
  advCells.resize(atmnet->numAtoms, VOR_CELL());
  vector<int> numNodes;
  int cellIndex = 0;

  vector< vector<double> > vertices;

  // Compute Voronoi cells
  c_loop_all_periodic vl(con);
  vector<int> cellIDs;
  int **cellInfo;
  cellInfo = new int*[atmnet->numAtoms];
    if(vl.start()) {
      do { 
	if(con.compute_cell(c,vl)) {
	  vvol+=c.volume();
	  vl.pos(id,x,y,z,r);
	  	  
	  int *map;
	  vector<double> coords;
	  c.vertices(atmnet->atoms[id].x, atmnet->atoms[id].y, atmnet->atoms[id].z, coords);

	  numNodes.push_back(c.p);
	  cellIDs.push_back(id);
	  vertices.push_back(coords);

	  vnet.add_to_network(c,id,x,y,z,r, map);

	  cellInfo[cellIndex] = map;
	  
	  if(storeAdvCells){
	    VOR_CELL newCell;
	    createAdvCell(c, coords, map, newCell);
	    advCells[id] = newCell;
	  }
	} 
	else {
	  numNodes.push_back(0);
	  cellIDs.push_back(-1);
	  vertices.push_back(vector<double>());
	  
	  cellInfo[cellIndex] = NULL;
	}
	cellIndex++;
      } while(vl.inc());
    }
    else {
      fputs("Error: Unable to begin Voronoi decomposition.\nExiting...\n",stderr);
      exit(1);
    }

    // Carry out the volume check
    printf("Volume check:\n  Total domain volume  = %f\n",bx*by*bz);
    printf("  Total Voronoi volume = %f\n", vvol);
//    if(abs(vvol - bx*by*bz) > 0.00001) {
//    if(abs(vvol - bx*by*bz) > 0.0001) {
//      fputs("Error: Voronoi volume check failed.\nExiting...\n",stderr);
//      exit(1);
//    }
    double box_vol = bx*by*bz;
    double error_percent = 100*abs(vvol - box_vol)/box_vol;
    double error_percent_tolerance = 0.001; // former (before Voro++ fit default = 0.1; 
    if(error_percent > error_percent_tolerance) {
      printf("Error: Voronoi volume check failed (%.3f%% error, > %.3f%% tolerance).\nExiting...\n", error_percent, error_percent_tolerance);
      return false;
//      exit(1);
    }
    
    cout << "Voronoi decomposition finished. Rerouting Voronoi network information." << "\n";
    
    vnet.store_network(vornet->nodes, vornet->edges, false);
    for(int i = 0; i < atmnet->numAtoms; i++){
      if(numNodes[i] == 0){
	continue;
      }
      
      vector<int> nodeIDs;
      vector<Point> nodeLocations; 
      if((int)vertices[i].size() != 3*numNodes[i]){
	cerr << "Error: Improper number of node coordinates in Voronoi decomposition" << "\n"
	     << "Found " << vertices[i].size() << " but expected " << 3*numNodes[i] << "\n"
	     << "Exiting..." << "\n";
	exit(1);
      }

      for(int j = 0; j < numNodes[i]; j++){
	nodeLocations.push_back(Point(vertices[i][3*j], vertices[i][3*j+1], vertices[i][3*j+2])) ;
	nodeIDs.push_back(cellInfo[i][j*4]);
      }
      
      basCells[cellIDs[i]] = BASIC_VCELL(nodeLocations, nodeIDs);
      delete [] cellInfo[i];
    }
    delete [] cellInfo;
    cout << "Finished rerouting information." << "\n";
    return true;
}

/** Extends the provided unit cell in the x, y and z directions using
    the given factors and stores the resulting ATOM_NETWORK using the
    pointer to newCell.For instance, if xfactor = yfactor = zfactor= 1, newCell contains the
    same information as cell. In contrast, if xfactor = yfactor =
    zfactor = 2, newCell is a 2x2x2 replication of cell in 3D space. */
void extendUnitCell(ATOM_NETWORK *cell, ATOM_NETWORK *newCell, int xfactor, int yfactor, int zfactor){
  newCell->a = cell->a*xfactor;
  newCell->b = cell->b*yfactor;
  newCell->c = cell->c*zfactor;
  newCell->alpha = cell->alpha;
  newCell->beta = cell->beta;
  newCell->gamma = cell->gamma;
  newCell->initialize();

  int numAtoms = 0;
  newCell->atoms.clear();
  for(unsigned int i = 0; i<cell->atoms.size(); i++){
    ATOM oldAtom = cell->atoms.at(i);
    for(int j = 0; j < xfactor; j++){
      for(int k = 0; k < yfactor; k++){
	for(int m = 0; m < zfactor; m++){
	  ATOM newAtom;
	  newAtom.specialID = i;
	  newAtom.type = oldAtom.type;
	  newAtom.radius = oldAtom.radius;
	  newAtom.a_coord = oldAtom.a_coord/xfactor + j*1.0/xfactor;
	  newAtom.b_coord = oldAtom.b_coord/yfactor + k*1.0/yfactor;
	  newAtom.c_coord = oldAtom.c_coord/zfactor + m*1.0/zfactor;
	  Point newCoords = newCell->abc_to_xyz(newAtom.a_coord, newAtom.b_coord, newAtom.c_coord);
	  newAtom.x = newCoords[0];
	  newAtom.y = newCoords[1];
	  newAtom.z = newCoords[2];
	  newCell->atoms.push_back(newAtom);
	  numAtoms++;
	}
      }
    }
  }
  newCell->numAtoms = numAtoms;
}



/** Extend the given VORONOI_NETWORK by 10 unit cells along a unit cell vector in the provided
 *  direction. The direction must be either (1,0,0), (0,1,0) or (0,0,1).  The list of atoms 
 *  belonging to each atom are NOT APPROPRIATELY REPLICATED.*/
void extendVorNet(VORONOI_NETWORK *vornet, VORONOI_NETWORK *newNet, DELTA_POS direction, map<int,int> *idAliases, set<int> *sourceNodes){
  idAliases->clear();
  sourceNodes->clear();

  DELTA_POS filter = DELTA_POS(1,1,1) - direction;
  int numIDs = vornet->nodes.size();
  
  for(unsigned int i = 0; i < vornet->edges.size(); i++){
    VOR_EDGE curEdge = vornet->edges.at(i);
    DELTA_POS dirComp = DELTA_POS(curEdge.delta_uc_x, curEdge.delta_uc_y, curEdge.delta_uc_z)*direction;
    if(dirComp.isZero())
      continue;
    else if ((dirComp.x < 0) || (dirComp.y < 0) || (dirComp.z < 0)){
      sourceNodes->insert(curEdge.to);
    }
  }

  newNet->nodes.clear();
  newNet->edges.clear();
  int factor = 10;
  newNet->v_a = vornet->v_a.scale(direction.x * factor);
  newNet->v_b = vornet->v_b.scale(direction.y * factor);
  newNet->v_c = vornet->v_c.scale(direction.z * factor);

  //Repeat copy process for each additional unit cell
  int idCount = 0;
  for(int i = 0; i <= factor; i++){
    
    //Translate all of the Voronoi nodes and assign them new id's
    for(int j = 0; j < numIDs; j++){
      VOR_NODE oldNode = vornet->nodes.at(j);
      VOR_NODE newNode; 
      newNode.x = oldNode.x + i*direction.x*vornet->v_a.x + i*direction.y*vornet->v_b.x + i*direction.z*vornet->v_c.x; 
      newNode.y = oldNode.y + i*direction.x*vornet->v_a.y + i*direction.y*vornet->v_b.y + i*direction.z*vornet->v_c.y; 
      newNode.z = oldNode.z + i*direction.x*vornet->v_a.z + i*direction.y*vornet->v_b.z + i*direction.z*vornet->v_c.z; 
      newNode.rad_stat_sphere = oldNode.rad_stat_sphere;
      newNet->nodes.push_back(newNode);
      if(sourceNodes->find(j) != sourceNodes->end()){
	sourceNodes->insert(idCount);
	idAliases->insert(pair<int,int> (idCount,j));
      }
      idCount++;
    }

    //Translate all of the Voronoi edges, modifying the connectivity information as is appropriate
    for(unsigned int j = 0; j < vornet->edges.size(); j++){
      VOR_EDGE oldEdge = vornet->edges.at(j);
      DELTA_POS oldDirection = DELTA_POS(oldEdge.delta_uc_x, oldEdge.delta_uc_y, oldEdge.delta_uc_z);
      DELTA_POS newDirection = oldDirection*filter;
      DELTA_POS dirComp = oldDirection*direction;

      int changeInTo;
      if(dirComp.isZero())
	changeInTo = 0;
      else if ((dirComp.x < 0) || (dirComp.y < 0) || (dirComp.z < 0))
	changeInTo = -1;
      else
	changeInTo = 1;
     
      int newTo = i + changeInTo;
      if(newTo < 0){
	newDirection = newDirection + (direction*(-1));
	newTo = factor;
      }
      else if (newTo > factor){
	newDirection = newDirection + direction;
	newTo = 0;
      }
 
      VOR_EDGE newEdge;
      newEdge.from = oldEdge.from + i*numIDs; 
      newEdge.to = oldEdge.to + newTo*numIDs; 
      newEdge.rad_moving_sphere = oldEdge.rad_moving_sphere; 
    
      newEdge.delta_uc_x = newDirection.x;
      newEdge.delta_uc_y = newDirection.y; 
      newEdge.delta_uc_z = newDirection.z;
      newEdge.length = oldEdge.length;
      newNet->edges.push_back(newEdge);
    }
  }
}

void calculateFreeSphereParameters(VORONOI_NETWORK *vornet, char *filename, bool extendedPrintout){
  vector<double> freeRadResults;
  vector<double> incRadResults;
  vector<bool> NtoN;

  DELTA_POS directions [3] = {DELTA_POS(1,0,0), DELTA_POS(0,1,0), DELTA_POS(0,0,1)};
  for(unsigned int i = 0; i < 3; i++){
    VORONOI_NETWORK newNet;
    set<int> sourceNodes;
    map<int,int> idAliases;
    extendVorNet(vornet, &newNet, directions[i], &idAliases, &sourceNodes);
   
    DIJKSTRA_NETWORK dnet;
    DIJKSTRA_NETWORK::buildDijkstraNetwork(&newNet,&dnet);

    TRAVERSAL_NETWORK analyzeNet = TRAVERSAL_NETWORK(directions[i].x,directions[i].y,directions[i].z, &dnet);    
    pair<bool,PATH> results = analyzeNet.findMaxFreeSphere(&idAliases, &sourceNodes);
    
    freeRadResults.push_back(2*results.second.max_radius);
    incRadResults.push_back(2*results.second.max_inc_radius);
    NtoN.push_back(results.first);
  }

  fstream output;
  output.setf(ios::fixed,ios::floatfield);  
  output.precision(5);
  output.width(12);
  output.open(filename, fstream::out);
  output << filename << "    " << 2 * findMaxIncludedSphere(vornet) << " ";

  double maxd=0.0; int maxdir=0;
  for(unsigned int i = 0; i < freeRadResults.size(); i++)
     {
     if(i==0) {maxd=freeRadResults[i]; maxdir=i;}
       else
       {
       if(maxd<freeRadResults[i])
         {
         maxd=freeRadResults[i]; 
         maxdir=i;
         }
       else if(maxd==freeRadResults[i])
         {
         if(incRadResults[maxdir]<incRadResults[i]) maxdir=i;
         };
       };
     };

  output << freeRadResults[maxdir] << "  " << incRadResults[maxdir];

  if(extendedPrintout==true){
    output << "  ";
    for(unsigned int i = 0; i < freeRadResults.size(); i++)
       output << freeRadResults[i] << "  ";

    for(unsigned int i = 0; i < incRadResults.size(); i++)
       output << incRadResults[i] << "  ";
    };

 /* 
  for(unsigned int i = 0; i < NtoN.size(); i++)
    output << (NtoN[i] ? "t" : "f") << "  ";
 */
  output << "\n";
}

/* New calculateFreeSphere function that works with MATERIAL class 
 * added by M Haranczyk / Feb 2014 
 */

void NEWcalculateFreeSphereParameters(MATERIAL *Mat){

  if(Mat->doneFlatVoroFlag == false)
    {
    Mat->runVoroFlat();
    };

  vector<double> freeRadResults;
  vector<double> incRadResults;
  vector<bool> NtoN;

  DELTA_POS directions [3] = {DELTA_POS(1,0,0), DELTA_POS(0,1,0), DELTA_POS(0,0,1)};
  for(unsigned int i = 0; i < 3; i++){
    VORONOI_NETWORK newNet;
    set<int> sourceNodes;
    map<int,int> idAliases;
    extendVorNet(&(Mat->vornet), &newNet, directions[i], &idAliases, &sourceNodes);
   
    DIJKSTRA_NETWORK dnet;
    DIJKSTRA_NETWORK::buildDijkstraNetwork(&newNet,&dnet);

    TRAVERSAL_NETWORK analyzeNet = TRAVERSAL_NETWORK(directions[i].x,directions[i].y,directions[i].z, &dnet);    
    pair<bool,PATH> results = analyzeNet.findMaxFreeSphere(&idAliases, &sourceNodes);
    
    freeRadResults.push_back(2*results.second.max_radius);
    incRadResults.push_back(2*results.second.max_inc_radius);
    NtoN.push_back(results.first);
  }

  double maxd=0.0; int maxdir=0;
  for(unsigned int i = 0; i < freeRadResults.size(); i++)
     {
     if(i==0) {maxd=freeRadResults[i]; maxdir=i;}
       else
       {
       if(maxd<freeRadResults[i])
         {
         maxd=freeRadResults[i]; 
         maxdir=i;
         }
       else if(maxd==freeRadResults[i])
         {
         if(incRadResults[maxdir]<incRadResults[i]) maxdir=i;
         };
       };
     };

 Mat->directionDf = freeRadResults; Mat->directionDif = incRadResults;
 Mat->Di = 2 * findMaxIncludedSphere(&(Mat->vornet));
 Mat->Df = freeRadResults[maxdir];
 Mat->Dif = incRadResults[maxdir];

} // ends NEWcalc...



void NEWcalculateFreeSphereParametersPrint(MATERIAL *Mat, char *filename, bool extendedPrintout){
  fstream output;
  output.setf(ios::fixed,ios::floatfield);
  output.precision(5);
  output.width(12);
  output.open(filename, fstream::out);

  output << filename << "    " << Mat->Di << " ";

  output << Mat->Df << "  " << Mat->Dif;

  if(extendedPrintout==true){
    output << "  ";
    for(unsigned int i = 0; i < Mat->directionDf.size(); i++)
       output << Mat->directionDf.at(i) << "  ";

    for(unsigned int i = 0; i < Mat->directionDif.size(); i++)
       output << Mat->directionDif.at(i) << "  ";
    };

 /* 
  for(unsigned int i = 0; i < NtoN.size(); i++)
    output << (NtoN[i] ? "t" : "f") << "  ";
 */

  output << "\n";


}



void viewVoronoiDecomp(ATOM_NETWORK *atmnet, double r_probe, string filename){
  ATOM_NETWORK newAtomNet;
  atmnet->copy(&newAtomNet);
  for(int i = 0; i < newAtomNet.numAtoms; i++){ newAtomNet.atoms[i].radius += r_probe; }
  VORONOI_NETWORK vornet;
  vector<BASIC_VCELL> bvcells;
  vector<VOR_CELL> advCells;
  container_periodic_poly *new_rad_con = (container_periodic_poly *)performVoronoiDecomp(true, &newAtomNet, &vornet, advCells, true, bvcells);
  writeZeoVisFile((char *)filename.data(), &advCells, &newAtomNet, &vornet);
  delete new_rad_con;
}


void loadRadii(ATOM_NETWORK *atmnet){
  vector <ATOM> ::iterator iter = atmnet->atoms.begin();
  while(iter != atmnet->atoms.end()){ 
    iter->radius = lookupRadius(iter->type, true);
    iter++;
  }
}

void loadMass(bool useMassFlag, ATOM_NETWORK *atmnet){
  vector <ATOM> ::iterator iter = atmnet->atoms.begin();
  while(iter != atmnet->atoms.end()){
    if(useMassFlag)
      {
      iter->mass = lookupMass(iter->type);
      }else
      {
      iter->mass = 0.0;
      };
    iter++;
  }
}




/* Print information about topology of the structure 
   It defines a periodic graph representing atom connectivity
   and then reuses findChannel function (required Dijkstra graph for atoms */
void getStructureInformation(char *filename, char *filenameExtendedOutput, ATOM_NETWORK *atmnet, bool extendedOutput){

 // to speed up calculations, identify atoms on the "surface" of the unit cell 

 vector <int> surfaceIDs;
 vector <bool> surfaceFlag;

 double surface_slab_d = 2.5; // cutoff distance from the surface 

 double a_step = surface_slab_d/atmnet->a;
 double b_step = surface_slab_d/atmnet->b;
 double c_step = surface_slab_d/atmnet->c;

 surfaceFlag.resize(atmnet->atoms.size(), false);

 if(a_step>=1.0 || b_step>=1.0 || c_step>=1.0)
   {//no speed up
   for(unsigned int i = 0; i < atmnet->atoms.size(); i++) 
      {
      surfaceIDs.push_back(i);
      surfaceFlag.at(i) = true;
      };
   cout << "Small unit cell. All(" << atmnet->atoms.size() << ") atoms considered to be on the unit cell surface." << endl;
   }else
   {
   for(unsigned int i = 0; i < atmnet->atoms.size(); i++)
      {
      double a = atmnet->atoms.at(i).a_coord; double b = atmnet->atoms.at(i).b_coord; double c = atmnet->atoms.at(i).c_coord;

      if(a < a_step || a > (1-a_step) || b < b_step || b > (1-b_step) || c < c_step || c > (1-c_step) )
        {
        surfaceIDs.push_back(i); surfaceFlag.at(i) = true;
        };
      };
   cout << "Big unit cell. " << surfaceIDs.size() << " out of " << atmnet->atoms.size() << " atoms considered to be on the unit cell surface. Surface definition = " << surface_slab_d << endl;
   };



 // Define connectivity of atoms in terms of Dijkstra graph
 DIJKSTRA_NETWORK DijkstraAtomNetwork;
 DijkstraAtomNetwork.v_a = atmnet->v_a;  DijkstraAtomNetwork.v_b = atmnet->v_b; DijkstraAtomNetwork.v_c = atmnet->v_c;

 for(unsigned int i = 0; i < atmnet->atoms.size(); i++)
    {
    DijkstraAtomNetwork.nodes.push_back(DIJKSTRA_NODE(i, atmnet->atoms.at(i).x, atmnet->atoms.at(i).y, atmnet->atoms.at(i).z,
                                                           lookupCovRadius(atmnet->atoms.at(i).type), true));
    };
 // now need to update CONN connections for each node

 for(unsigned int i = 0; i < DijkstraAtomNetwork.nodes.size(); i++)
    {
    XYZ atom1(DijkstraAtomNetwork.nodes.at(i).x, DijkstraAtomNetwork.nodes.at(i).y, DijkstraAtomNetwork.nodes.at(i).z);
// CONN(int myFrom, int myTo, double len, double maxR, DELTA_POS deltaP);
// DELTA_POS(int myX = 0, int myY = 0, int myZ = 0);

    for(int x=-1; x<2; x++) // loop over all unit cells
    for(int y=-1; y<2; y++)
    for(int z=-1; z<2; z++)
       {
       if(x==0 && y==0 && z==0)
         { // connectivity inside cell
         for(unsigned int j = i+1; j < DijkstraAtomNetwork.nodes.size(); j++)
            {
            XYZ atom2(DijkstraAtomNetwork.nodes.at(j).x, DijkstraAtomNetwork.nodes.at(j).y, DijkstraAtomNetwork.nodes.at(j).z);
            XYZ va = DijkstraAtomNetwork.v_a; XYZ vb = DijkstraAtomNetwork.v_b; XYZ vc = DijkstraAtomNetwork.v_c;
            atom2 = atom2 + va.scale(x) + vb.scale(y) + vc.scale(z);
            double dist = calcEuclideanDistance(atom1.x, atom1.y, atom1.z, atom2.x, atom2.y, atom2.z);
            if(dist < DijkstraAtomNetwork.nodes.at(i).max_radius + DijkstraAtomNetwork.nodes.at(j).max_radius + 0.4)
               { // connected;
               DijkstraAtomNetwork.nodes.at(i).connections.push_back(CONN(i, j, dist, 0.10, DELTA_POS (x,y,z)));
               DijkstraAtomNetwork.nodes.at(j).connections.push_back(CONN(j, i, dist, 0.10, DELTA_POS (x,y,z)));
               };
            };
         }else
         { // connectivity outside cell
         if(surfaceFlag.at(i) == true) // only investigate if surface atom (other will not have a chance to be connected)
           {
           for(unsigned int k = 0; k < surfaceIDs.size(); k++)
              {
              unsigned int j = surfaceIDs.at(k);

              XYZ atom2(DijkstraAtomNetwork.nodes.at(j).x, DijkstraAtomNetwork.nodes.at(j).y, DijkstraAtomNetwork.nodes.at(j).z);
              XYZ va = DijkstraAtomNetwork.v_a; XYZ vb = DijkstraAtomNetwork.v_b; XYZ vc = DijkstraAtomNetwork.v_c; 
              atom2 = atom2 + va.scale(x) + vb.scale(y) + vc.scale(z);

              double dist = calcEuclideanDistance(atom1.x, atom1.y, atom1.z, atom2.x, atom2.y, atom2.z);
              if(dist < DijkstraAtomNetwork.nodes.at(i).max_radius + DijkstraAtomNetwork.nodes.at(j).max_radius + 0.4)
                 { // connected
                 DijkstraAtomNetwork.nodes.at(i).connections.push_back(CONN(i, j, dist, 0.10, DELTA_POS (x,y,z)));
                 };
              }; // ends loop over k
            }; // ends surfaceFlag.at(i) if

        }; // end loop over cells outside the central cell
       }; //ends loop over supercell (x,y,z)
    };


 // Find frameworks
 vector<bool> accessInfo;
 vector<PORE> structureSegments;
 PORE::findChannelsAndPockets(&DijkstraAtomNetwork, &accessInfo, &structureSegments);

 // Save output
 fstream output;
 output.open(filename, fstream::out);
 int nPockets=0, nChan1D=0, nChan2D=0, nChan3D=0;
 for(unsigned int i = 0; i < structureSegments.size(); i++){
    if(structureSegments[i].dimensionality == 0) nPockets++;
    if(structureSegments[i].dimensionality == 1) nChan1D++;
    if(structureSegments[i].dimensionality == 2) nChan2D++;
    if(structureSegments[i].dimensionality == 3) nChan3D++;
    };
 output << filename << "   "<< atmnet->returnChemicalFormula() << "   " << structureSegments.size() << " segments: " << (structureSegments.size() - nPockets) << " framework(s) (1D/2D/3D " 
        << nChan1D << " " << nChan2D << " " << nChan3D
        << " ) and  "  << nPockets << "  molecule(s). Identified dimensionality of framework(s): ";
 for(unsigned int i = 0; i < structureSegments.size(); i++){
    if(structureSegments[i].dimensionality > 0 ) output << structureSegments[i].dimensionality << " ";
    };
 output << "\n";
 output.close();

 // Save extended output: coordinates with framework and molecule ID


// original function
/*
 if(extendedOutput == true)
   {
   fstream output2;
   output2.open(filenameExtendedOutput, fstream::out);

   for(unsigned int i = 0; i < structureSegments.size(); i++)
      {
      for(unsigned int j = 0; j< structureSegments[i].reverseIDMappings.size(); j++)
         {
         int atomID = structureSegments[i].reverseIDMappings.find(j)->second;
         output2 << atmnet->atoms.at(atomID).type << "    ";
//         output2 << atmnet->atoms.at(atomID).x << "   " << atmnet->atoms.at(atomID).y << "   "
//                 << atmnet->atoms.at(atomID).z << "   " ;
         output2 << structureSegments[i].nodes.at(j).x << "   " << structureSegments[i].nodes.at(j).y << "  "
                 << structureSegments[i].nodes.at(j).z << "   ";
         output2 << atmnet->atoms.at(atomID).a_coord << "   "
                 << atmnet->atoms.at(atomID).b_coord << "   " << atmnet->atoms.at(atomID).c_coord << "   "
                 << structureSegments[i].dimensionality << "   " << i << "   ";
         output2 << structureSegments[i].unitCells.size();
         output2 << "\n";
         };
      };

   output2 << "#Center of mass information(XYZ coords, no atoms, ABC coordinates, and radii)\n";
   for(unsigned int i = 0; i < structureSegments.size(); i++)
      {
      pair <XYZ,double> COM = structureSegments[i].getCenterOfMass();
      Point p(COM.first.x, COM.first.y, COM.first.z);
      Point pabc = atmnet->xyz_to_abc(atmnet->shiftXYZInUC(p));
      output2 << "X    " << COM.first.x << "  " << COM.first.y << "  " << COM.first.z << "      ";
      output2 << structureSegments[i].nodes.size() << "    ";
      output2 << pabc[0] << "  " << pabc[1] << "  " << pabc[2] << "     ";
      output2 << i << "  " << 2.0*(COM.second-1.5) << "\n";   // 1.5 stands for average radii of atoms
                                                              // 2.0 multiplyer to convert to diameters
      };

   output2.close();
   };
*/


/* Analysis of cages ONLY 
 * for Ismael  */
/*
 if(extendedOutput == true)
   {
   fstream output2;
   output2.open(filenameExtendedOutput, fstream::out);

   for(unsigned int i = 0; i < structureSegments.size(); i++)
      {
      if(structureSegments[i].dimensionality == 0)
        {
        vector < pair <int,XYZ> > cageReconstructedCoords = structureSegments[i].getReconstructedPore();
        
        for(unsigned int j = 0; j< structureSegments[i].reverseIDMappings.size(); j++)
           {
           int atomID = structureSegments[i].reverseIDMappings.find((cageReconstructedCoords.at(j).first))->second;
           output2 << atmnet->atoms.at(atomID).type << "    ";
           output2 << atmnet->atoms.at(atomID).x << "   " << atmnet->atoms.at(atomID).y << "   "
                   << atmnet->atoms.at(atomID).z << "   " ;
           output2 << atmnet->atoms.at(atomID).radius << "   ";
           output2 << structureSegments[i].nodes.at(j).x << "   " << structureSegments[i].nodes.at(j).y << "  "
                   << structureSegments[i].nodes.at(j).z << "   ";
*/
/* DEBUG
           output2 << atmnet->atoms.at(atomID).a_coord << "   "
                   << atmnet->atoms.at(atomID).b_coord << "   " << atmnet->atoms.at(atomID).c_coord << "   "
                   << structureSegments[i].dimensionality << "   " << i << "   ";
         output2 << structureSegments[i].unitCells.size();
*/
/*
           output2 << cageReconstructedCoords.at(j).second.x << "  " << cageReconstructedCoords.at(j).second.y << "  "
                   << cageReconstructedCoords.at(j).second.z << "  ";
          output2 << "  " << i << "  ";
          output2 << "\n";
            };
         };
      };


   output2.close();
   };
*/

 /* Alternative version of the above, which saves multiple copies of molecules that cross boundry */
//getReconstructredPoresWithCrossBoundryCopies

 if(extendedOutput == true)
   {
   fstream output2;
   output2.open(filenameExtendedOutput, fstream::out);

   for(unsigned int i = 0; i < structureSegments.size(); i++)
      {
      if(structureSegments[i].dimensionality == 0)
        {
        vector< vector < pair <int,XYZ> > > cagesMultipleCopiesReconstructedCoords = structureSegments[i].getReconstructredPoresWithCrossBoundryCopies();

        vector < pair <int,int> > cageBonds;

        for(unsigned int m = 0; m < cagesMultipleCopiesReconstructedCoords[0].size(); m++)
           for(unsigned int n = m + 1; n < cagesMultipleCopiesReconstructedCoords[0].size(); n++)
              {
              int atomID1 = structureSegments[i].reverseIDMappings.find((cagesMultipleCopiesReconstructedCoords[0].at(m).first))->second;
              int atomID2 = structureSegments[i].reverseIDMappings.find((cagesMultipleCopiesReconstructedCoords[0].at(n).first))->second;

              if(cagesMultipleCopiesReconstructedCoords[0].at(m).second.euclid_dist(cagesMultipleCopiesReconstructedCoords[0].at(n).second)
                     < ( lookupCovRadius(atmnet->atoms.at(atomID1).type) + lookupCovRadius(atmnet->atoms.at(atomID2).type) + 0.45 ))
                     {
                     cageBonds.push_back(pair <int,int> (m,n));
                     //if(i==28) cout << m << "   " << n << endl;
                     };
              };
//        cout << "Molecule: " << i << " Atoms: " << cagesMultipleCopiesReconstructedCoords[0].size() << " Bonds: " << cageBonds.size() << endl;

        for(unsigned int k = 0; k< cagesMultipleCopiesReconstructedCoords.size(); k++)
           {
           /* beging loop over all atoms/to be saved in frmaid format */
/*
           for(unsigned int j = 0; j< structureSegments[i].reverseIDMappings.size(); j++)
              {
              int atomID = structureSegments[i].reverseIDMappings.find((cagesMultipleCopiesReconstructedCoords[k].at(j).first))->second;
              output2 << atmnet->atoms.at(atomID).type << "    ";
              output2 << atmnet->atoms.at(atomID).x << "   " << atmnet->atoms.at(atomID).y << "   "
                      << atmnet->atoms.at(atomID).z << "   " ;
              output2 << atmnet->atoms.at(atomID).radius << "   ";
              output2 << structureSegments[i].nodes.at(j).x << "   " << structureSegments[i].nodes.at(j).y << "  "
                      << structureSegments[i].nodes.at(j).z << "   ";
              output2 << cagesMultipleCopiesReconstructedCoords[k].at(j).second.x << "  " << cagesMultipleCopiesReconstructedCoords[k].at(j).second.y << "  "
                   << cagesMultipleCopiesReconstructedCoords[k].at(j).second.z << "  ";
              output2 << "  " << i << "  " << k << " ";
              output2 << "\n";
              };
 */
           /* ends saving file in frameid format */

           /* beging loop over all atoms in SDF format */

           /* SDF header */
           output2 << i << endl;
           output2 << "Zeomolecule" << endl << endl;
           output2 << "  " << structureSegments[i].reverseIDMappings.size() << " " << cageBonds.size() << " 0     0  0  0  0  0  0999 V2000\n"; 
           for(unsigned int j = 0; j< structureSegments[i].reverseIDMappings.size(); j++)
              {
              int atomID = structureSegments[i].reverseIDMappings.find((cagesMultipleCopiesReconstructedCoords[k].at(j).first))->second;
              output2 << "    ";
              output2 << cagesMultipleCopiesReconstructedCoords[k].at(j).second.x << "  " << cagesMultipleCopiesReconstructedCoords[k].at(j).second.y << "  "
                   << cagesMultipleCopiesReconstructedCoords[k].at(j).second.z << "  ";
              output2 << atmnet->atoms.at(atomID).type << "   0  0  0  0  0  0  0  0  0  0  0  0\n";
              };
           for(unsigned int k = 0; k < cageBonds.size(); k++)
              {
              output2 << cageBonds[k].first << "  " << cageBonds[k].second << " 1 0 0 0 0\n" ;
              };
           output2 << "M  END\n$$$$\n";  
           /* ends saving file in SDF format */

           }; // ends for loop over copies of the segment i (each molecule)
         }; // ends if(structureSegments[i].dimensionality == 0)
      };


   output2.close();
   };




/* Save additional atom statistics */
/*
 if(extendedOutput == true)
   {
   fstream output2;
   output2.open(filenameExtendedOutput, fstream::out);


   int nMetals=0, nCs=0, nHs=0;
   std::map<std::string,int> metalmap;

//  mymap['a']="an element";
   for(unsigned int i=0; i < atmnet->atoms.size(); i++)
      {
      if(atmnet->atoms.at(i).type == "C") nCs++;
      if(atmnet->atoms.at(i).type == "H") nHs++;
      if(isMetal(atmnet->atoms.at(i).type) == true) 
        {
        nMetals++;
        if(metalmap[atmnet->atoms.at(i).type] == 0) metalmap[atmnet->atoms.at(i).type] =1;
         else metalmap[atmnet->atoms.at(i).type] = metalmap[atmnet->atoms.at(i).type] +1;
        };
      };

   output2 << filename << "    nCs: " << nCs << "  nHs: " << nHs << "  nMetals: " << nMetals << " nMetalTypes: " << metalmap.size() << "  ";

   std::map<std::string,int>::iterator it;
   for(it=metalmap.begin(); it!=metalmap.end(); ++it)
      {
      output2 << "  " << (*it).first << " " << (*it).second; 
      };


//   for(auto& x: metalmap) {
//      output2 << "  "  << x.first << " " << x.second;
//      };


   output2 << "\n";
   output2.close();
   };
*/

/* quick fix for Greg */
/*
 cout << "Running experimetnal pruning for Greg, uncomment this routine after use !!!!!!!!!!\n";
 if(extendedOutput == true)
   {
   fstream output2;
   output2.open(filenameExtendedOutput, fstream::out);

   output2 << "Processing: " << filenameExtendedOutput << "\n";
   output2 << "Unit_cell: " << atmnet->a <<"  " <<  atmnet->b <<"  " <<  atmnet->c << "  "
           << atmnet->alpha << "  " << atmnet->beta << "   " << atmnet->gamma << "\n";

   for(unsigned int i = 0; i < structureSegments.size(); i++)
      {
      if(structureSegments[i].dimensionality>2) 
      {
      for(unsigned int j = 0; j< structureSegments[i].reverseIDMappings.size(); j++)
         {
         int atomID = structureSegments[i].reverseIDMappings.find(j)->second;
         output2 << atmnet->atoms.at(atomID).type << "    ";
         output2 << atmnet->atoms.at(atomID).a_coord << "   "
                 << atmnet->atoms.at(atomID).b_coord << "   " << atmnet->atoms.at(atomID).c_coord << "\n";
         };
      };

      };

   output2.close();
   };
*/
}



/* Print information about the presence of open metal sites 
   */
   
void getOMSInformation(char *filename, char *filenameExtendedOutput, ATOM_NETWORK *atmnet, bool extendedOutput){

 // General consts
 const double Threshold = PI/16.;

 // Output variables
 int nOMS=0;
 vector< vector<int> > OMS_atomIDs;

 // to speed up calculations, identify atoms on the "surface" of the unit cell 

 vector <int> surfaceIDs;
 vector <bool> surfaceFlag;

 double surface_slab_d = 3.5; // cutoff distance from the surface 

 double a_step = surface_slab_d/atmnet->a;
 double b_step = surface_slab_d/atmnet->b;
 double c_step = surface_slab_d/atmnet->c;

 surfaceFlag.resize(atmnet->atoms.size(), false);

 if(a_step>=1.0 || b_step>=1.0 || c_step>=1.0)
   {//no speed up
   for(unsigned int i = 0; i < atmnet->atoms.size(); i++) 
      {
      surfaceIDs.push_back(i);
      surfaceFlag.at(i) = true;
      };
   cout << "Small unit cell. All(" << atmnet->atoms.size() << ") atoms considered to be on the unit cell surface." << endl;
   }else
   {
   for(unsigned int i = 0; i < atmnet->atoms.size(); i++)
      {
      double a = atmnet->atoms.at(i).a_coord; double b = atmnet->atoms.at(i).b_coord; double c = atmnet->atoms.at(i).c_coord;

      if(a < a_step || a > (1-a_step) || b < b_step || b > (1-b_step) || c < c_step || c > (1-c_step) )
        {
        surfaceIDs.push_back(i); surfaceFlag.at(i) = true;
        };
      };
   cout << "Big unit cell. " << surfaceIDs.size() << " out of " << atmnet->atoms.size() << " atoms considered to be on the unit cell surface. Surface definition = " << surface_slab_d << endl;
   };


 for(unsigned int i = 0; i < atmnet->atoms.size(); i++)
    {
//    DijkstraAtomNetwork.nodes.push_back(DIJKSTRA_NODE(i, atmnet->atoms.at(i).x, atmnet->atoms.at(i).y, atmnet->atoms.at(i).z,
//                                                           lookupCovRadius(atmnet->atoms.at(i).type), true));

    if(isMetal(atmnet->atoms.at(i).type) == true)
     {

     vector< vector <double> > cluster;
     vector <int> cluster_atomIDs;

     XYZ atom1(atmnet->atoms.at(i).x, atmnet->atoms.at(i).y, atmnet->atoms.at(i).z);
     vector <double> origin;
     origin.push_back(atom1.x); origin.push_back(atom1.y); origin.push_back(atom1.z);
     cluster.push_back(origin);
     cluster_atomIDs.push_back(i);

     for(int x=-1; x<2; x++) // loop over all unit cells
     for(int y=-1; y<2; y++)
     for(int z=-1; z<2; z++)
       {
       if(x==0 && y==0 && z==0)
         { // connectivity inside cell
         for(unsigned int j = 0; j < atmnet->atoms.size(); j++)
            {
            if(i!=j)
              {
              XYZ atom2(atmnet->atoms.at(j).x, atmnet->atoms.at(j).y, atmnet->atoms.at(j).z);
              XYZ va = atmnet->v_a; XYZ vb = atmnet->v_b; XYZ vc = atmnet->v_c;
              atom2 = atom2 + va.scale(x) + vb.scale(y) + vc.scale(z);
              double dist = calcEuclideanDistance(atom1.x, atom1.y, atom1.z, atom2.x, atom2.y, atom2.z);
              if(dist < lookupCovRadius(atmnet->atoms.at(i).type) + lookupCovRadius(atmnet->atoms.at(j).type)  + 0.4)
               { // connected;
               vector <double> clusteratom;
               clusteratom.push_back(atom2.x); clusteratom.push_back(atom2.y); clusteratom.push_back(atom2.z);
               cluster.push_back(clusteratom);
               cluster_atomIDs.push_back(j);
               };
              };
            };
         }else
         { // connectivity outside cell
         if(surfaceFlag.at(i) == true) // only investigate if surface atom (other will not have a chance to be connected)
           {
           for(unsigned int k = 0; k < surfaceIDs.size(); k++)
              {
              unsigned int j = surfaceIDs.at(k);

              XYZ atom2(atmnet->atoms.at(j).x, atmnet->atoms.at(j).y, atmnet->atoms.at(j).z);
              XYZ va = atmnet->v_a; XYZ vb = atmnet->v_b; XYZ vc = atmnet->v_c;
              atom2 = atom2 + va.scale(x) + vb.scale(y) + vc.scale(z);

              double dist = calcEuclideanDistance(atom1.x, atom1.y, atom1.z, atom2.x, atom2.y, atom2.z);
              if(dist < lookupCovRadius(atmnet->atoms.at(i).type) + lookupCovRadius(atmnet->atoms.at(j).type)  + 0.4)
               { // connected;
               vector <double> clusteratom;
               clusteratom.push_back(atom2.x); clusteratom.push_back(atom2.y); clusteratom.push_back(atom2.z);
               cluster.push_back(clusteratom);
               cluster_atomIDs.push_back(j);
               };



              }; // ends loop over k
            }; // ends surfaceFlag.at(i) if

        }; // end loop over cells outside the central cell
       }; //ends loop over supercell (x,y,z)

   if(IsExposedMoleculeThreshold(cluster, Threshold)) 
      {
      nOMS++;
      OMS_atomIDs.push_back(cluster_atomIDs);
      };

   }; // ends if isMetal
  }; // finishes loop over all atoms

 // Save output
 fstream output;
 output.open(filename, fstream::out);
 output << filename << " #OMS=  "<< nOMS << "\n";
 output.close();

 // Save extended output: coordinates with framework and molecule ID




/* Save additional atom statistics */

 if(extendedOutput == true)
   {
   fstream output2;
   output2.open(filenameExtendedOutput, fstream::out);

   for(unsigned int i=0; i < OMS_atomIDs.size(); i++)
      {
      output2 << atmnet->atoms.at(OMS_atomIDs[i].at(0)).type << "     ";
      output2 << OMS_atomIDs[i].size()-1 << "  ";
      for(unsigned int j=1; j < OMS_atomIDs[i].size(); j++)
         output2 << atmnet->atoms.at(OMS_atomIDs[i].at(j)).type << " ";
      output2 << "\n";
      };

   output2.close();

   }; // ends extended output


}
