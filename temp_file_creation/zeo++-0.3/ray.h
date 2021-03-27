#ifndef RAY_H
#define RAY_H

//#include "network.h"
#include <iostream>
#include <vector>
#include "networkstorage.h"

class hitdata{
 public:
  bool hit;
  Point hitpoint;
  double dist;
	int id;
  void *object; //information on hit object (should be sphere)
  hitdata()
    {
      hit=false;
      hitpoint = Point(0,0,0);
      dist = 0.0;
      id = -1;
      object = NULL;
    }
};

class ray{
 public:
  Point base;
  Point vector;
  //Methods
};

class Sphere{
 public:  
	Point center;
  double radius;
	//Methods
	void hitSphere(ray r,hitdata& hitsphere);
};

class Plane{
 public:	
	Point point;
	Point normal;
	//Methods
	void hitPlane(ray r,hitdata& hitplane);
	double distToPlane(Point pnt);
};

/** Returns a class hitdata which has important information of the closest object hit. Function calls hitSphere for each sphere in the list **/
hitdata findClosestSphere(std::vector<Sphere>& s,ray& r);

/** Returns a class hitdata which has important information of the closest object hit. Function calls hitPlane for each plane in the list **/
hitdata findClosestPlane(std::vector<Plane>& p,ray& r);

/** Calculates the normals to the faces of the unit cell order is v_a x v_b : v_b x v_c : v_a x v_c as unit vectors and creates a vector with all the planes**/
void calcPlanesToCell(ATOM_NETWORK *cell,std::vector<Plane> &p);

/** Since all the opperations in Zeo++ are done within one unit cell it is important to duplicate the atoms that would overlap with the cell walls on the other side **/
void duplicateSpheresOnFaces(ATOM_NETWORK *cell,std::vector<Sphere>& spheres,std::vector<Plane>& faces);

/** This function converts nodes to spheres it allows for a more flexible implementation of ray_tracing thus making higher abstractions easier we are only concerned with the accessible nodes thus only these will be copied **/
void convertNodeToSphere(VORONOI_NETWORK& vornet,std::vector<Sphere>& nodes,std::vector<bool> *accessInfo);

/** This function converts atoms to spheres it allows for a more flexible implementation of ray_tracing thus making higher abstractions easier **/
void convertAtomToSphere(ATOM_NETWORK& cell,std::vector<Sphere>& atoms);

/** This function finds a nodes radius that encompasses the point raypoint. It will not count the previous node that the raypoint was in by using the reference id**/
bool findSphereOfPoint(Point p,std::vector<Sphere>& spheres,int& id);

/** Calculates the equivalent xyz point from xyz coordinates in next unit_cell assuming that it is only one unit cell away. This covers all cases edges, corners, and faces sadly shiftXYZInUC did not work for all cases **/
Point calcPeriodicPoint(Point xyz_coor,ATOM_NETWORK *cell);

/**Generates a random vector direction and unitizes it. Since we want a random ray orientation in cartesian coordinates, sphereical coordinates are used to create a point **/
Point genRandomVec();

/** This function generates a random point within the unit cell with values 0-1 for a,b, and c **/
Point genRandomPoint();

/** This is a recursive function that will travel within spheres until it hits a point that is no longer within spheres in the direction of the vector. Returns hitdata on the point it hits**/
void rayTraceInsideSphere(ATOM_NETWORK *cell,std::vector<Sphere>& spheres,ray r,hitdata& hitobject);

/** This is a recursive function that will keep on extending a ray unitl it hits an atom or exceeds the maximum distance allowed (MAXRAYDIST)**/
void rayTraceToSphere(ATOM_NETWORK *cell,std::vector<Sphere>& spheres,ray r,std::vector<Plane>& faces,hitdata& hitobject);

/** Returns data from shooting a bunch of rays from accesible regions in the network until they hit an atom. Since the cell is periodic a max ray length can be set so ray does not shoot forever*/
void calcRaysInAV(ATOM_NETWORK *hiaccatmnet, ATOM_NETWORK *orgatmnet, bool highAccuracy, double r_probe_chan,double r_probe,int numSamples, std::ostream &output, bool visualize,std::string option);

void reportRayInfo(std::vector<ray>& axsray);
void reportRays(std::ostream &output, std::vector<ray>& axsray, std::vector<ray>& inaxsray,bool colors);
void reportAtoms(std::ostream &output, std::vector<Sphere>& s);
void reportNodes(std::ostream &output,std::vector<Sphere>& s);
void reportHistogram(std::ostream& output,double binSize,int maxBins,std::vector<ray>& rays);

#endif

