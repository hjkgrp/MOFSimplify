/*
*  PORE SIZE DISTRIBUTION CALCULATION
* Author: Marielle Pinheiro
* Date: Fall 2012
*
* The function modified to fit the new Zeo++ framework; M. Haranczyk, Spring 2014
*
* Info:
*
* The pore size distribution (PSD) function is a further extension of accessible volume function defined in area_and_volume.cc
*
* The PSD function generates sample points across a unit cell; for each point, accessibility is determined.
* Accessible points are stored. Each point is compared to all of the nodes in the Voronoi network. If the point
* is within the radius of one or more node sphere, the point's coordinate and the largest encompassing 
* radius are stored. If the point is not within the radius of any node, the point is stored, and the default
* radius is 0. These stored radii are the lower bounds for our next test.
* Next, a new "ghost" Voronoi cell centered on each sample point is created in relation to existing atoms, and 
* the distance from the sample point to each of the the vertices of the nodes is calculated. The maximum distance
* out of this set is compared to the stored radius value. If distance > node radius, this is the radius of the 
* largest sphere that will encapsulate our sample point; otherwise, the node radius remains the largest sphere.
*
* A histogram with bin size BINSTEP is generated (this is set to 0.1 by default). The output file (PSD histogram) is
* (materialName).psd_histo
*
* A file containing a list of xyz coordinates of accessible points and largest sphere radii is generated for
* visualization in Visit. The output file is (materialName).vpsdpts
*
* A file containing a list of spheres for visualization is also available in one of two forms. 
* With the flag -vpsd, it will produce a list of spheres labeled either 1 (node sphere) or 0 (ghost cell sphere).
* The latter file is named (materialName).vpsdradii
*/

//#include "network.h"
#include <fstream>
#include <iomanip>
#include <cassert>
#include <voro++.hh>

#include "voronoicell.h"
#include "channel.h"
#include "zeo_consts.h"
#include "network.h"
#include "networkaccessibility.h"
#include "material.h"
#include "psd.h"

using namespace std;
using namespace voro;

static const float BINSTEP=0.1;
//#define PI 3.14159265358973
static const float TOL=0.9;
 
 
void calcPoreSizeDistr(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgAtomnet, bool highAccuracy, double r_probe_chan, double r_probe, int numSamples, bool excludePockets, string histogramFile, string pointsFile, string nodeAndRadiiFile, string spheresDistFile, bool visualize, bool visVISIT){

    ofstream output; 
    if (!histogramFile.empty()) output.open(histogramFile.data());
    ofstream outfile;
    if (!pointsFile.empty()) outfile.open(pointsFile.data());
    ofstream nodeAndRadii;
    if (!nodeAndRadiiFile.empty()) nodeAndRadii.open(nodeAndRadiiFile.data());
    ofstream spheresDist;
  // Create an object that handles analysis of accessibility of sampled points
  AccessibilityClass accessAnalysis;
  if(highAccuracy) accessAnalysis.setupAndFindChannels(atmnet, orgAtomnet, highAccuracy, r_probe_chan, r_probe);
    else accessAnalysis.setupAndFindChannels(atmnet, atmnet, highAccuracy, r_probe_chan, r_probe);


  srand(randSeed);

//Vectors to store accessible points, inaccessible points, and accessible points which lie outside Voronoi nodes and edges. 
  vector<Point> axsPoint = vector<Point>();
  vector<Point> inaxsPoint = vector<Point>();

//List of largest sphere diameters for PSD histogram
  vector<double> diameterOfLargestSphere;

//List of sample points and radius of sphere that contains sample point for visualization in Visit
  vector <pair <Point, double> > axsPoint_and_radiiOfNodeSphere;
  vector <pair <Point, double> > axsPoint_and_radiiOfGhostSphere;
//List of node coordinates and radii for visualization
  vector <int> nodeIDs;
  vector <pair <Point, double> > coord_and_radiiOfNodeSphere;
  vector <pair <Point, double> > coord_and_radiiOfGhostSphere;
	vector <double> distanceBetweenNodeSpheres;
	vector <double> distanceBetweenGhostSpheres;

  int count = 0;
  int incount = 0;
  int outcount = 0;

//Generates random points across cell and tests accessibility.
  for(int i = 0; i < numSamples; i++){
    bool overlaps = false;
    bool inside = false;
    double aPoint = (rand()*1.0)/RAND_MAX;
    double bPoint = (rand()*1.0)/RAND_MAX;
    double cPoint = (rand()*1.0)/RAND_MAX;
    Point samplingPoint = atmnet->abc_to_xyz(aPoint, bPoint, cPoint);


    // Calling accessibility object to determine accessibility of the point (this replaced a big chunk of code by Thomas)
    double dist_var;
    pair<bool,bool> answer = (accessAnalysis.isVPointInsideAtomAndNotAccessible(samplingPoint, dist_var));
    inside = answer.first; overlaps = answer.second;
    if(accessAnalysis.needToResample() == true) i--; // the sampled point could not be analyzed in isVPointInsideAtomAndNotAccessible() function, resampling needed

    if(inside == false && excludePockets == false) overlaps = false; // if ignore inacceible pockets, treat the point as accessible (unless inside atom)




// Store sampled points that did not overlap with an atom but were inaccessible
    if(accessAnalysis.needToResample() == false && !inside && overlaps){
	Point abcCoords = Point(aPoint, bPoint, cPoint);
	Point coord = atmnet->abc_to_xyz(abcCoords);
	inaxsPoint.push_back(coord);
    }
// Store accessible points
    if(accessAnalysis.needToResample() == false &&!overlaps) {
	count++;     
	Point abcCoords = Point(aPoint, bPoint, cPoint);
	Point coords = atmnet->abc_to_xyz(abcCoords);
	axsPoint.push_back(coords); 
    }
  } // ends loop over all sampled points


  // Warn user if points were resampled
  int resampleCount =  accessAnalysis.getResampleCount();
  if(resampleCount != 0){
    cerr << "\n" << "\n"
         << "Warning: Resampled " << resampleCount << " points out of " << numSamples
         << " when analyzing " << atmnet->name << "\n"
         << "\n" << "\n";
  }


//Below are functions that go beyond accessible volume calcluation. Each accessible point is analyzed to find the largest encapsulating sphere


//Determine whether sample point is inside/outside node sphere
  for (int a=0; a<axsPoint.size(); a++){
    Point samplePt = axsPoint.at(a);
    Point maxNodeCoord = Point (0, 0, 0);
    double maxRadius = 0;
    for (int b=0; b<accessAnalysis.vornet.nodes.size(); b++)
      {
      bool access = accessAnalysis.accessInfo.at(b);
      if (access){
        Point nodePt = Point(accessAnalysis.vornet.nodes.at(b).x, accessAnalysis.vornet.nodes.at(b).y, accessAnalysis.vornet.nodes.at(b).z);
        double radius = accessAnalysis.vornet.nodes.at(b).rad_stat_sphere + r_probe;
        double sampToNode = accessAnalysis.analyzedAtomNet.calcDistanceXYZ(samplePt[0], samplePt[1], samplePt[2], nodePt[0], nodePt[1], nodePt[2]);
        bool inSphere = (sampToNode <= radius); 
        //if point is inside more than one sphere, chooses largest sphere (i.e. largest radius)
        if (inSphere){
          maxRadius = max(radius, maxRadius);
          //If query zpsd, this stores the largest node sphere and its coordinate for visualization in Visit
          if (visualize){
            if (radius==maxRadius){
               maxNodeCoord = nodePt; 
               }
            }
        } // ends if(inSphere)
      }
     }
    if (maxRadius==0) outcount++; 
    if (maxRadius!=0) incount++;





//Generate new "ghost" Voronoi cell
    double pointParticleRadius = r_probe;
//DEBUG next 4 lines
    voronoicell v(*(accessAnalysis.new_rad_con));
    accessAnalysis.new_rad_con->compute_ghost_cell(v, samplePt[0], samplePt[1], samplePt[2], pointParticleRadius);
    vector<double> vertexCoords;
    v.vertices(vertexCoords);

   //debug
   //cout << "number of ghost vertices = " << vertexCoords.size()/3 << "\n"; cout.flush();
   //for(unsigned int aa=0;aa<vertexCoords.size();aa++) cout << vertexCoords[aa] << " "; cout << "\n";

    vector <Point> pointVertexCoords;
    for (int c=0; c<vertexCoords.size(); c= c+3){
      Point vertexPoint = Point(vertexCoords.at(c), vertexCoords.at(c+1), vertexCoords.at(c+2));
      pointVertexCoords.push_back(vertexPoint);
    }

   //debug
   //cout << "vertex vector constructed \n"; cout.flush();

    double maxDistFromCtoV = 0;
    Point maxGhostCoord = Point (0, 0, 0);
//Calculate distance from center of cell to vertices
   for (int b=0; b<pointVertexCoords.size(); b++){
      Point vertexGlobalCoord = Point((samplePt[0]+pointVertexCoords.at(b)[0]), (samplePt[1]+pointVertexCoords.at(b)[1]), (samplePt[2]+pointVertexCoords.at(b)[2]));
      double calcDistFromCtoV = calcEuclideanDistance(pointVertexCoords.at(b)[0], pointVertexCoords.at(b)[1], pointVertexCoords.at(b)[2], 0.0, 0.0, 0.0) - pointParticleRadius;

      maxDistFromCtoV = max(maxDistFromCtoV, calcDistFromCtoV);
//If query zpsd, this stores the largest ghost sphere and its coordinate for visualization in Visit
      if (visualize){
        if (calcDistFromCtoV==maxDistFromCtoV){
          maxGhostCoord = vertexGlobalCoord;
        } 
      }
    }

  //debug:
  //cout << "getting ready for final comparison and assigment\n"; cout.flush();

//Compares node radius to biggest ghost radius and stores the larger of the two values accordingly
    double origNodeRadius = maxRadius;
    double ghostRadius = maxDistFromCtoV;  
    bool compareRadii = (origNodeRadius >= ghostRadius); 
    if (compareRadii) {
      diameterOfLargestSphere.push_back(2.0*origNodeRadius); 
      if (visualize){
      axsPoint_and_radiiOfNodeSphere.push_back(pair<Point, double>(samplePt, origNodeRadius));
      checkDuplicates(&(accessAnalysis.analyzedAtomNet), &(accessAnalysis.vornet), nodeIDs, coord_and_radiiOfNodeSphere, maxNodeCoord, origNodeRadius);
      }
    }
    if (!compareRadii) {
      diameterOfLargestSphere.push_back(2.0*ghostRadius);
      if (visualize){
        axsPoint_and_radiiOfGhostSphere.push_back(pair<Point, double>(samplePt, ghostRadius));
        coord_and_radiiOfGhostSphere.push_back(pair <Point, double> (maxGhostCoord, maxDistFromCtoV));
      }
    }

  //debug:
  //cout << "finishing pass in axsPoint loop\n"; cout.flush();

  }

  double volumeFraction = count*1.0/numSamples;
  double nodeVolFraction = incount*1.0/numSamples;
  double outFraction = outcount*1.0/numSamples;


  Histogram(output, BINSTEP, MAXBINS, diameterOfLargestSphere, count, nodeVolFraction, outFraction, numSamples);

  printf("Pore size distribution calculated.\n\n");

  if (visualize) {
    printFileCoords_Radii(outfile, axsPoint_and_radiiOfNodeSphere, axsPoint_and_radiiOfGhostSphere, distanceBetweenNodeSpheres, distanceBetweenGhostSpheres, false, -1, -1);
    printf("Points file generated.\n\n");
		if (!visVISIT){
                // here printing ZeoVis stuff
		}
		if (visVISIT) {
		  printFileCoords_Radii(nodeAndRadii, coord_and_radiiOfNodeSphere, coord_and_radiiOfGhostSphere, distanceBetweenNodeSpheres, distanceBetweenGhostSpheres, false, -1, -1);
		  printf("Radii file generated.\n");
		}
  }

  accessAnalysis.deconstruct();

} // end PSD calculation



/* NEW version of PSD function that operates on MATERIAL class, and shares data through the class */ 
void NEWcalcPoreSizeDistr(MATERIAL *Mat, ostream &output){

/* PSD debug section (old code that is not used */
  bool PSDDEBUG = false; // this flag enables execution of old parts of PSD output
                         // that can be used for debuging  

  string pointsFile, nodeAndRadiiFile;
  if(PSDDEBUG == true) { pointsFile = "DEBUG_PSDpointsfile"; nodeAndRadiiFile = "DEBUG_PSDnodeAndRadiiFile"; };
  ofstream outfile;
  if (!pointsFile.empty()) outfile.open(pointsFile.data());
  ofstream nodeAndRadii;
  if (!nodeAndRadiiFile.empty()) nodeAndRadii.open(nodeAndRadiiFile.data());

//The following data structures are used to PSDDEBUG

//List of sample points and radius of sphere that contains sample point for visualization in Visit
  vector <pair <Point, double> > axsPoint_and_radiiOfNodeSphere;
  vector <pair <Point, double> > axsPoint_and_radiiOfGhostSphere;
//List of node coordinates and radii for visualization
  vector <int> nodeIDs;
  vector <pair <Point, double> > coord_and_radiiOfNodeSphere;
  vector <pair <Point, double> > coord_and_radiiOfGhostSphere;
	vector <double> distanceBetweenNodeSpheres;
	vector <double> distanceBetweenGhostSpheres;

/* end PSD debug */


/* PSD codes assumes that AV run was executed and AV data structures in MATERIALS class are filled with data */

  int count = Mat->AVcount; // number of AV samples that turned accessible
  int incount = 0; // count number of points inside Voronoi nodes
  int outcount = 0; 

  cout << "PSD calculation for " << count << " points(" << Mat->AVaxsPoints.size() << ").\n"; 

  double pointParticleRadius; // used to be set to = Mat->AVprobeRadius; 

           // After vivid discussions with Rich, we concluded that the ghost particle radius
           // should be set to the smallest radius in the system
           // this gives the correct voronoi network, but the distance from the vericies of the ghost cell to the particle are not correct
           // representation of the pore diamter
  for(int i=0; i<Mat->accessAnalysis.analyzedAtomNet->atoms.size(); i++)
     {
     if(i == 0) pointParticleRadius = Mat->accessAnalysis.analyzedAtomNet->atoms.at(i).radius;
     else
     if(Mat->accessAnalysis.analyzedAtomNet->atoms.at(i).radius < pointParticleRadius) pointParticleRadius = Mat->accessAnalysis.analyzedAtomNet->atoms.at(i).radius;
     };
  cout << "Radius of ghost cell particle = " << pointParticleRadius << "\n";




//Below are functions that go beyond accessible volume calcluation. Each accessible point is analyzed to find the largest encapsulating sphere

//Determine whether sample point is inside/outside node sphere
  for (int a=0; a<Mat->AVaxsPoints.size(); a++){
    Point samplePt = Mat->accessAnalysis.analyzedAtomNet->abc_to_xyz(Mat->AVaxsPoints.at(a));
    Point maxNodeCoord = Point (0, 0, 0);
    double maxRadius = 0;
    for (int b=0; b<Mat->accessAnalysis.vornet.nodes.size(); b++)
      {
      bool access = Mat->accessAnalysis.accessInfo.at(b);
      if (access){
        Point nodePt = Point(Mat->accessAnalysis.vornet.nodes.at(b).x, Mat->accessAnalysis.vornet.nodes.at(b).y, Mat->accessAnalysis.vornet.nodes.at(b).z);
        double radius = Mat->accessAnalysis.vornet.nodes.at(b).rad_stat_sphere;
         //double radius = Mat->accessAnalysis.vornet.nodes.at(b).rad_stat_sphere + Mat->AVprobeRadius;
                        // commented out in non-inflated version
        double sampToNode = Mat->accessAnalysis.analyzedAtomNet->calcDistanceXYZ(samplePt[0], samplePt[1], samplePt[2], nodePt[0], nodePt[1], nodePt[2]);
        bool inSphere = (sampToNode <= radius); 
        //if point is inside more than one sphere, chooses largest sphere (i.e. largest radius)
        if (inSphere){
          maxRadius = max(radius, maxRadius);
          //If query zpsd, this stores the largest node sphere and its coordinate for visualization in Visit (only for PSD DEBUG)
          if (radius==maxRadius){
             maxNodeCoord = nodePt; 
             }
        } // ends if(inSphere)
      }
     }
    if (maxRadius==0) outcount++; 
    if (maxRadius!=0) incount++;



//Generate new "ghost" Voronoi cell

    // ghost cell particle radius (pointParticleRadius) is now set to the radius of the 
    // smallest atom (before the main PSD loop)

    voronoicell v(*(Mat->accessAnalysis.new_rad_con));
    Mat->accessAnalysis.new_rad_con->compute_ghost_cell(v, samplePt[0], samplePt[1], samplePt[2], pointParticleRadius);
    vector<double> vertexCoords;
    v.vertices(vertexCoords);

    vector <Point> pointVertexCoords;
    for (int a=0; a<vertexCoords.size(); a= a+3){
      Point vertexPoint = Point(vertexCoords.at(a), vertexCoords.at(a+1), vertexCoords.at(a+2));
      pointVertexCoords.push_back(vertexPoint);
    }

    double maxDistFromCtoV = 0;
    Point maxGhostCoord = Point (0, 0, 0);
//Calculate distance from center of cell to vertices
   for (int b=0; b<pointVertexCoords.size(); b++){
      Point vertexGlobalCoord = Point((samplePt[0]+pointVertexCoords.at(b)[0]), (samplePt[1]+pointVertexCoords.at(b)[1]), (samplePt[2]+pointVertexCoords.at(b)[2]));
      double calcDistFromCtoV = calcEuclideanDistance(pointVertexCoords.at(b)[0], pointVertexCoords.at(b)[1], pointVertexCoords.at(b)[2], 0.0, 0.0, 0.0) - pointParticleRadius;
      maxDistFromCtoV = max(maxDistFromCtoV, calcDistFromCtoV);
//If query zpsd, this stores the largest ghost sphere and its coordinate for visualization in Visit (only for PSDDEBUG)
      if (calcDistFromCtoV==maxDistFromCtoV){
         maxGhostCoord = vertexGlobalCoord;
         }
    }

//Compares node radius to biggest ghost radius and stores the larger of the two values accordingly
    double origNodeRadius = maxRadius;
    double ghostRadius = maxDistFromCtoV;  
    bool compareRadii = (origNodeRadius >= ghostRadius); 
    if (compareRadii) {
      Mat->AVaxsPointsPSD.push_back(2.0*origNodeRadius); 
      if (PSDDEBUG){
      axsPoint_and_radiiOfNodeSphere.push_back(pair<Point, double>(samplePt, origNodeRadius));
      checkDuplicates((Mat->accessAnalysis.analyzedAtomNet), &(Mat->accessAnalysis.vornet), nodeIDs, coord_and_radiiOfNodeSphere, maxNodeCoord, origNodeRadius);
      }
    }
    if (!compareRadii) {
      Mat->AVaxsPointsPSD.push_back(2.0*ghostRadius);
      if (PSDDEBUG){
        axsPoint_and_radiiOfGhostSphere.push_back(pair<Point, double>(samplePt, ghostRadius));
        coord_and_radiiOfGhostSphere.push_back(pair <Point, double> (maxGhostCoord, maxDistFromCtoV));
      }
    }
  }

  double volumeFraction = count*1.0/Mat->AVnumSamples;
  double nodeVolFraction = incount*1.0/Mat->AVnumSamples;
  double outFraction = outcount*1.0/Mat->AVnumSamples;


  Histogram(output, BINSTEP, MAXBINS, (Mat->AVaxsPointsPSD), count, nodeVolFraction, outFraction, Mat->AVnumSamples);

  printf("Pore size distribution calculated.\n\n");

  if (PSDDEBUG) {
    printFileCoords_Radii(outfile, axsPoint_and_radiiOfNodeSphere, axsPoint_and_radiiOfGhostSphere, distanceBetweenNodeSpheres, distanceBetweenGhostSpheres, false, -1, -1);
    printf("Points file generated.\n\n");
    printFileCoords_Radii(nodeAndRadii, coord_and_radiiOfNodeSphere, coord_and_radiiOfGhostSphere, distanceBetweenNodeSpheres, distanceBetweenGhostSpheres, false, -1, -1);
    printf("Radii file generated.\n");
    }


} // end PSD calculation






/* backup copy of PSD function

 
void calcPoreSizeDistr(ATOM_NETWORK *atmnet, double r_probe_chan, double r_probe, int numSamples, bool excludePockets, ostream &output, ostream &outfile, ostream &nodeAndRadii, ostream &spheresDist, bool visualize, bool visVISIT){
//Creates temporary copy of atomic network in which each atom's radius has been increased by probe radius.
  ATOM_NETWORK newAtomNet;
  atmnet->copy(&newAtomNet);

  for(int i = 0; i < newAtomNet.numAtoms; i++){ 
    newAtomNet.atoms[i].radius += r_probe; 
  }  

  VORONOI_NETWORK vornet;
  vector<BASIC_VCELL> vorcells;
  vector<VOR_CELL> advCells;

  container_periodic_poly *new_rad_con = (container_periodic_poly *)performVoronoiDecomp(true, &newAtomNet, &vornet, advCells, false, vorcells);

  vector<CHANNEL> channels = vector<CHANNEL>();
  vector<bool> accessInfo = vector<bool> ();

  CHANNEL::findChannels(&vornet, max(0.0, r_probe_chan - r_probe), &accessInfo, &channels);
srand(randSeed);

//Vectors to store accessible points, inaccessible points, and accessible points which lie outside Voronoi nodes and edges. 
  vector<Point> axsPoint = vector<Point>();
  vector<Point> inaxsPoint = vector<Point>();
  vector< pair<int, Point> > resampledInfo = vector< pair<int, Point> > (); 

//List of largest sphere diameters for PSD histogram
  vector<double> diameterOfLargestSphere;

//List of sample points and radius of sphere that contains sample point for visualization in Visit
  vector <pair <Point, double> > axsPoint_and_radiiOfNodeSphere;
  vector <pair <Point, double> > axsPoint_and_radiiOfGhostSphere;
//List of node coordinates and radii for visualization
  vector <int> nodeIDs;
  vector <pair <Point, double> > coord_and_radiiOfNodeSphere;
  vector <pair <Point, double> > coord_and_radiiOfGhostSphere;
	vector <double> distanceBetweenNodeSpheres;
	vector <double> distanceBetweenGhostSpheres;

  int resampleCount = 0;
  int count = 0;
  int incount = 0;
  int outcount = 0;

//Generates random points across cell and tests accessibility.
  for(int i = 0; i < numSamples; i++){
    bool overlaps = false;
    double aPoint = (rand()*1.0)/RAND_MAX;
    double bPoint = (rand()*1.0)/RAND_MAX;
    double cPoint = (rand()*1.0)/RAND_MAX;
    Point samplingPoint = atmnet->abc_to_xyz(aPoint, bPoint, cPoint);

    double newAtomX, newAtomY, newAtomZ;
    int minAtomID;
    bool foundCell = new_rad_con->find_voronoi_cell(samplingPoint[0], samplingPoint[1], samplingPoint[2], newAtomX, newAtomY, newAtomZ, minAtomID);
    if(!foundCell){
	cerr << "Error: Unable to find Voronoi cell for sampled point in AV calculation." << "\n" << "Occurred for structure " << newAtomNet.name << "\n" << "Exiting..." << "\n";
	exit(1);
    }
    ATOM curAtom = atmnet->atoms[minAtomID];

// Adjust sampling point so that it lies within the Voronoi cell of interest constructed previously
    Point newSamplingPoint = (samplingPoint.add(Point(curAtom.x, curAtom.y, curAtom.z).subtract(Point(newAtomX, newAtomY, newAtomZ)))); 

    double minDist = calcEuclideanDistance(newSamplingPoint[0], newSamplingPoint[1], newSamplingPoint[2], curAtom.x, curAtom.y, curAtom.z);
    if(minDist < r_probe + curAtom.radius - 0.00000001) overlaps = true; 
    bool inside = overlaps;
// If necessary, check Voronoi nodes of cell to determine accessibility of point
    if(!overlaps && excludePockets){
      BASIC_VCELL vcell = vorcells[minAtomID]; //stores voronoi cell location at minID point
      Point circCenter = Point(curAtom.x, curAtom.y, curAtom.z);
      double samplingRadius = minDist; 
      Point sampleRay = Point(newSamplingPoint[0]-curAtom.x, newSamplingPoint[1]-curAtom.y, newSamplingPoint[2]-curAtom.z); 
      
// Scan the nodes in the Voronoi cell to find if line can be drawn from the node to the sampling point
      bool foundNode = false;
      if(vcell.getNumNodes() == 0){
	cerr << "Error: Voronoi cell of sampled point does not have any nodes" << "\nPoint: " << newSamplingPoint[0] << " " << newSamplingPoint[1] << " " << newSamplingPoint[2] << "\nVoronoi cell is #" << minAtomID << " in structure " << newAtomNet.name << "\nPlease contact the source code provider." << "\nExiting...\n";
        exit(1);
      }
      for(int k = 0; k < vcell.getNumNodes(); k++){
	Point nodePoint = vcell.getNodeCoord(k); 
	double ptDist= calcEuclideanDistance(nodePoint[0], nodePoint[1], nodePoint[2], circCenter[0], circCenter[1], circCenter[2]);
	bool nodeInsideSphere = (ptDist < samplingRadius);
	if(!nodeInsideSphere){
	  Point otherRay = newSamplingPoint.subtract(nodePoint);
	  double dotProduct = sampleRay.dot_product(otherRay);
	  if(dotProduct > 0) {
	    // Angle is less than 90 degrees and so the line segment intersects twice,
	    // making the path not viable
	  }
	  else {
	    // Angle is at least 90 degrees and so the line segment interesects only once, 
	    // thereby representing a viable path--> overlaps is now FALSE
	    foundNode = true;
	    overlaps = !accessInfo.at(vcell.getNodeID(k));
	    break; 
	  }
	}
      }
// Sampling failed due to lying on Voronoi cell face and numerical inaccurarcy. 
// Record failure, resample and notify user later
      if(!foundNode){
	resampleCount++;
	resampledInfo.push_back(pair<int, Point> (minAtomID, newSamplingPoint));
	i--;
      }
    }

// Store sampled points that did not overlap with an atom but were inaccessible
    if(!inside && overlaps){
	Point abcCoords = Point(aPoint, bPoint, cPoint);
	Point coord = atmnet->abc_to_xyz(abcCoords);
	inaxsPoint.push_back(coord);
    }
// Store accessible points
    if(!overlaps) {
	count++;     
	Point abcCoords = Point(aPoint, bPoint, cPoint);
	Point coords = atmnet->abc_to_xyz(abcCoords);
	axsPoint.push_back(coords); 
    }
  }

  cout<< "There are "<<count<<" accessible samples and "<<resampleCount<<" resamples\n\n";


//Determine whether sample point is inside/outside node sphere
  for (int a=0; a<axsPoint.size(); a++){
    Point samplePt = axsPoint.at(a);
    Point maxNodeCoord = Point (0, 0, 0);
    double maxRadius = 0;
    for (int b=0; b<vornet.nodes.size(); b++)
      {
      bool access = accessInfo.at(b);
      if (access){
        Point nodePt = Point(vornet.nodes.at(b).x, vornet.nodes.at(b).y, vornet.nodes.at(b).z);
        double radius = vornet.nodes.at(b).rad_stat_sphere + r_probe;
        double sampToNode = newAtomNet.calcDistanceXYZ(samplePt[0], samplePt[1], samplePt[2], nodePt[0], nodePt[1], nodePt[2]);
        bool inSphere = (sampToNode <= radius); 
        //if point is inside more than one sphere, chooses largest sphere (i.e. largest radius)
        if (inSphere){
          maxRadius = max(radius, maxRadius);
          //If query zpsd, this stores the largest node sphere and its coordinate for visualization in Visit
          if (visualize){
            if (radius==maxRadius){
               maxNodeCoord = nodePt; 
               }
            }
        } // ends if(inSphere)
      }
     }
    if (maxRadius==0) outcount++; 
    if (maxRadius!=0) incount++;

//Generate new "ghost" Voronoi cell
    double pointParticleRadius = r_probe;
    voronoicell v(*new_rad_con);
    new_rad_con->compute_ghost_cell(v, samplePt[0], samplePt[1], samplePt[2], pointParticleRadius);
    vector<double> vertexCoords;
    v.vertices(vertexCoords);

    vector <Point> pointVertexCoords;
    for (int a=0; a<vertexCoords.size(); a= a+3){
      Point vertexPoint = Point(vertexCoords.at(a), vertexCoords.at(a+1), vertexCoords.at(a+2));
      pointVertexCoords.push_back(vertexPoint);
    }

    double maxDistFromCtoV = 0;
    Point maxGhostCoord = Point (0, 0, 0);
//Calculate distance from center of cell to vertices
   for (int b=0; b<pointVertexCoords.size(); b++){
      Point vertexGlobalCoord = Point((samplePt[0]+pointVertexCoords.at(b)[0]), (samplePt[1]+pointVertexCoords.at(b)[1]), (samplePt[2]+pointVertexCoords.at(b)[2]));
      double calcDistFromCtoV = calcEuclideanDistance(pointVertexCoords.at(b)[0], pointVertexCoords.at(b)[1], pointVertexCoords.at(b)[2], 0.0, 0.0, 0.0) - pointParticleRadius;
      maxDistFromCtoV = max(maxDistFromCtoV, calcDistFromCtoV);
//If query zpsd, this stores the largest ghost sphere and its coordinate for visualization in Visit
      if (visualize){
        if (calcDistFromCtoV==maxDistFromCtoV){
          maxGhostCoord = vertexGlobalCoord;
        } 
      }
    }

//Compares node radius to biggest ghost radius and stores the larger of the two values accordingly
    double origNodeRadius = maxRadius;
    double ghostRadius = maxDistFromCtoV;  
    bool compareRadii = (origNodeRadius >= ghostRadius); 
    if (compareRadii) {
      diameterOfLargestSphere.push_back(2.0*origNodeRadius); 
      if (visualize){
      axsPoint_and_radiiOfNodeSphere.push_back(pair<Point, double>(samplePt, origNodeRadius));
      checkDuplicates(newAtomNet, vornet, nodeIDs, coord_and_radiiOfNodeSphere, maxNodeCoord, origNodeRadius);
      }
    }
    if (!compareRadii) {
      diameterOfLargestSphere.push_back(2.0*ghostRadius);
      if (visualize){
        axsPoint_and_radiiOfGhostSphere.push_back(pair<Point, double>(samplePt, ghostRadius));
        coord_and_radiiOfGhostSphere.push_back(pair <Point, double> (maxGhostCoord, maxDistFromCtoV));
      }
    }
  }

  double volumeFraction = count*1.0/numSamples;
  double nodeVolFraction = incount*1.0/numSamples;
  double outFraction = outcount*1.0/numSamples;


  Histogram(output, BINSTEP, MAXBINS, diameterOfLargestSphere, count, nodeVolFraction, outFraction);

  printf("Pore size distribution calculated.\n\n");

  if (visualize) {
    printFileCoords_Radii(outfile, axsPoint_and_radiiOfNodeSphere, axsPoint_and_radiiOfGhostSphere, distanceBetweenNodeSpheres, distanceBetweenGhostSpheres, false, -1, -1);
    printf("Points file generated.\n\n");
		if (!visVISIT){
                // here printing ZeoVis stuff
		}
		if (visVISIT) {
		  printFileCoords_Radii(nodeAndRadii, coord_and_radiiOfNodeSphere, coord_and_radiiOfGhostSphere, distanceBetweenNodeSpheres, distanceBetweenGhostSpheres, false, -1, -1);
		  printf("Radii file generated.\n");
		}
  }
}



// ends backup copy of PSD
*/


//This function checks for duplicates in the node IDs vector
void checkDuplicates(ATOM_NETWORK *atomnetwork, VORONOI_NETWORK *vornet, vector <int> &listNodeIDs, vector <pair <Point, double> > &nodeAndRadii, Point nodePt, double newRadius){
  int newNodeID = getNodeID(nodePt, atomnetwork, vornet);
	bool unique = true;
  listNodeIDs.push_back(newNodeID);
  nodeAndRadii.push_back(pair <Point, double> (nodePt, newRadius));
  for (int i=0; i<(listNodeIDs.size()-1) && unique; i++){
    if (newNodeID == listNodeIDs.at(i)){
      listNodeIDs.pop_back();
      nodeAndRadii.pop_back();
			unique = false;
    }
  }
}

//This function checks the amount of overlap between spheres. 
//In the circumstance that the volume of intersection overlap is above the threshold, the larger sphere is selected
void calcSphereIntersect(vector<pair <Point, double> > &inputNodeSphereAndRadii, vector<pair <Point, double> > &outputNodeSphereAndRadii){
	int size = inputNodeSphereAndRadii.size(); 
	vector <int> overlap_nodes; 
	for (int a=0; a<inputNodeSphereAndRadii.size(); a++){
		bool noOverlap = true;
		for (int b=a+1; b<inputNodeSphereAndRadii.size() && noOverlap; b++){
			double Radius1= inputNodeSphereAndRadii.at(a). second; double Radius2=inputNodeSphereAndRadii.at(b).second;
			Point coord1 = inputNodeSphereAndRadii.at(a).first; Point coord2 = inputNodeSphereAndRadii.at(b).first;

			double dist = calcEuclideanDistance(coord1[0], coord1[1], coord1[2], coord2[0], coord2[1], coord2[2]);
			double vol1 = (4*PI/3)*(Radius1*Radius1*Radius1);
			double vol2 = (4*PI/3)*(Radius2*Radius2*Radius2);

			bool overlap = (dist < (Radius1 + Radius2));
			if (overlap){
				bool size = (Radius1>Radius2);
				bool inside = (dist <= abs(Radius1-Radius2));
				/*if (inside && size) {
					printf("Node %d encompasses Node %d\n", a, b);
				}
				if (inside && !size){
					printf("Node %d encompasses Node %d\n", b, a);
				}*/
				if (!inside){
				double distsq = dist*dist; 
				double Rad1sq = Radius1*Radius1; double Rad2sq = Radius2*Radius2;
				double volIntersect = (PI/(12*dist))*(Radius1 + Radius2 - dist)*(Radius1 + Radius2 - dist)*(distsq+2*dist*(Radius1 + Radius2) - 3*(Radius1 - Radius2)*(Radius1-Radius2));
				double overlapDist = Radius1 + Radius2 - dist;

				double diff; double volFraction;
				int largerID; double largerRadius; Point largerCoord;
				double smallerRadius; Point smallerCoord;
				if (size){
					diff = 2*Radius2 - overlapDist; 
					volFraction = volIntersect/ vol2; 
				}
				if (!size){
					diff = 2*Radius1 - overlapDist; 
					volFraction = volIntersect/vol1;
				}
				if (volFraction>=TOL) {
					noOverlap = false;
				}
			}
	  }
	}
	if (noOverlap) outputNodeSphereAndRadii.push_back(pair <Point, double> (inputNodeSphereAndRadii.at(a).first, inputNodeSphereAndRadii.at(a).second));
	}
}

//This function calculates the periodic distance between spheres with radii in the range of interest (specified by last 2 arguments in -vpsd)
void calcSpheresDistance(ATOM_NETWORK atomnetwork, vector <pair <Point, double> > &sphereAndRadii, vector <double> &distanceBetweenSpheres){
	for (int a=0; a<sphereAndRadii.size(); a++){
		Point coord1 = sphereAndRadii.at(a).first; double Radius1 = sphereAndRadii.at(a).second;
		double SpheresDistanceCalc = 0;
		double minDistance = 1000000;
		for (int b=0; b<sphereAndRadii.size(); b++){
			Point coord2 = sphereAndRadii.at(b).first; double Radius2 = sphereAndRadii.at(b).second;
			if (Radius1>=2 && Radius1<=2.15){
				if (Radius2>=2 && Radius2<=2.15){
					if (a !=b){
						double checkDistance = atomnetwork.calcDistanceXYZ(coord1[0], coord1[1], coord1[2], coord2[0], coord2[1], coord2[2]);
						minDistance = min(minDistance, checkDistance);
						SpheresDistanceCalc = minDistance; 
					}
				}
			}
		}
		distanceBetweenSpheres.push_back(SpheresDistanceCalc);
	}	
}

//This function will output a histogram file with extension .distr
//R1: bin step R2: Count R3: Cumulative distribution R4: Derivative of Cumulative distribution (PSD)
void Histogram(ostream& output, const double binSize, const int maxBins, vector<double>& diam, int count, double nodefrac, double outfrac, int numSamples){
  assert(binSize > threshold);
  int bins[maxBins]; 
  double cumBins[maxBins];
  double derivBins[maxBins];
  for (int i=0; i<maxBins; i++){
    bins[i] = 0;
    cumBins[i] = 0;
    derivBins[i] = 0;
  }
 
  int bin;
  for (unsigned int i=0; i<diam.size(); i++){
    bin = diam.at(i)/binSize;
    if (bin >= maxBins){
      bin = maxBins - 1;
    }
    bins[bin]++;
    for (int j=0; j<bin+1; j++){
      cumBins[j]++;
    }
  }
//Cumulative distribution: 1 is added to the bin PLUS all preceding bins. 
//Each bin is divided by the value of the first bin.
  double maxC= cumBins[0];
  for (unsigned int b=0; b<maxBins; b++){
    cumBins[b] = (cumBins[b]/maxC);
  }
//Derivative of cumulative distribution: This is pore size distribution. 
  double deriv;
  for (unsigned int n=1; n<maxBins-1; n++){
    double fxh= cumBins[n+1];
    double fx = cumBins[n-1];
    deriv = (fxh-fx)/(2.0*binSize);
    if (deriv!=0) deriv= deriv*-1.0;
    if (deriv >= maxBins){
      deriv = maxBins - 1;
    }
    derivBins[n]=deriv;
  }
 
  double numBin;
  numBin = binSize * 1.0 * maxBins;

  //new format - easier to parse
  output << "Pore size distribution histogram\nBin size (A): "<<binSize<<"\nNumber of bins: "<<maxBins<<"\nFrom: 0\nTo: "<<numBin<<"\nTotal samples: "<<numSamples<<"\nAccessible samples: "<<count<<"\nFraction of sample points in node spheres: "<<nodefrac<<"\nFraction of sample points outside node spheres: "<<outfrac<<"\n\nBin Count Cumulative_dist Derivative_dist\n";
  for (int i=0; i<maxBins; i++){
    output<<binSize*1.0*i<<" "<<bins[i]<<" "<<cumBins[i]<<" "<<derivBins[i]<<"\n";
  }
}



//This function will output a list of coordinates and nodes
//R1: Specifies whether point is in node sphere (1) or ghost sphere (0). R2-R4: xyz coordinates. R5: sphere radius
//*Modification: 0=black (r>2.15) 1= white (r> 2.15) 2=lt blue (r<2) 3= purple (2<=r<=2.15)
void printFileCoords_Radii(ostream& outfile, vector<pair <Point, double> > &pt_node_rad, vector<pair <Point, double> > &pt_ghost_rad, vector <double> &distanceBetweenNodeSpheres, vector <double> &distanceBetweenGhostSpheres, bool overlaps, double startRad, double endRad){
 outfile <<(pt_node_rad.size() + pt_ghost_rad.size())<<"\n\n";
  for (int i=0; i<pt_node_rad.size(); i++){
		if (overlaps){
			if (pt_node_rad.at(i).second>=startRad && pt_node_rad.at(i).second <= endRad) outfile <<"3\t";
			if (pt_node_rad.at(i).second<startRad) outfile <<"2\t";
		  if (pt_node_rad.at(i).second>endRad) outfile <<"1\t"; 
		}
		else outfile <<"1\t";
		outfile << pt_node_rad.at(i).first[0]<< "\t"<< pt_node_rad.at(i).first[1]<< "\t" << pt_node_rad.at(i).first[2]<<"\t"<<pt_node_rad.at(i).second;
		if (overlaps) outfile <<"\t"<< distanceBetweenNodeSpheres.at(i);
		outfile<<"\n";
  }
 for (int i=0; i<pt_ghost_rad.size(); i++){
	if (overlaps){
		if (pt_ghost_rad.at(i).second>=startRad && pt_ghost_rad.at(i).second <= endRad) outfile <<"3\t";
		if (pt_ghost_rad.at(i).second<startRad) outfile <<"2\t";
		if (pt_ghost_rad.at(i).second>endRad) outfile <<"1\t";
	}
	else outfile <<"0\t";
	outfile<< pt_ghost_rad.at(i).first[0]<< "\t"<< pt_ghost_rad.at(i).first[1]<< "\t" << pt_ghost_rad.at(i).first[2]<<"\t"<<pt_ghost_rad.at(i).second;
	if (overlaps) outfile <<"\t"<< distanceBetweenGhostSpheres.at(i);
	outfile<<"\n";
  }
}

