//#include "network.h"
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <voro++.hh>

#include "ray.h"
#include "zeo_consts.h"
#include "channel.h"
#include "area_and_volume.h"
#include "network.h"
#include "networkaccessibility.h"

using namespace std;
using namespace voro;


/** Using |Point_0 - center| - radius^2 = 0 and Point_0 = t(Ray_dir) + Ray_pnt. Solve
    for t which is the distance. **/
void Sphere::hitSphere(ray r,hitdata& hitsphere)
{
  Point newCenter = center-r.base; //map the rand_ray to the origin and move system
  double temp=(r.vector*newCenter)*(r.vector*newCenter)-(newCenter*newCenter)+radius*radius;
  if (temp>0) //check inside of sqrt is not negative
    { 
      if ((r.vector * newCenter) - sqrt(temp) > 0) //If it hits the sphere twice pick the closest
	{ 
	  hitsphere.hit=true;
	  hitsphere.dist=(r.vector*newCenter) - sqrt(temp);
	  hitsphere.hitpoint=r.base+r.vector.scale(hitsphere.dist);
	}
      else if ((r.vector * newCenter) + sqrt(temp) > 0) //Ray was inside sphere (this is possible in some cases)
	{ 
	  hitsphere.hit=true;
	  hitsphere.dist=(r.vector*newCenter) + sqrt(temp);
	  hitsphere.hitpoint=r.base+r.vector.scale(hitsphere.dist);
    	} 
    }
  return;
}

/** Using normal dot (Point_Plane - Point_0| = 0 and Point_Plane = t(Ray_dir) + Ray_pnt. Solve for t which is the distance to plane. **/
void Plane::hitPlane(ray r,hitdata& hitplane)
{
  if (normal*r.vector >= 0){ //ray is either parallel or infront of the plane
    return;
  }	
  
  double dist=((normal*point)-(normal*r.base))/(normal*r.vector);
  if (dist < 0){ //There is no reason that the distance should be negative
    cerr << "Error: Distance = " << dist << " This means ray got outside of unitcell" << endl;
    cerr << "Point: " << r.base << " Vector: " << r.vector << endl;
    abort();
  }
  
  hitplane.hit=true;
  hitplane.dist=dist;
  hitplane.hitpoint=r.base + r.vector.scale(dist); 
  
  return;
}

/** Returns the shortest distance from point to plane **/
double Plane::distToPlane(Point pnt)
{
  return fabs(normal *(pnt - point));
}

/** Returns a class hitdata which has important information of the closest object hit. Function calls hitSphere for each sphere in the list **/
hitdata findClosestSphere(vector<Sphere>& s,ray& r)
{
  hitdata hitsphere;
  hitdata hitobject;
  assert(hitobject.object == NULL);
  for (unsigned int i=0;i<s.size();i++)
    {
      s[i].hitSphere(r,hitsphere);
      if (hitsphere.hit == true)			
	{
	  if(hitobject.hit == false || hitobject.dist > hitsphere.dist) //Conditions for this to be closest sphere found so far
	    {
	      hitobject.hit = true;
	      hitobject.dist = hitsphere.dist;
	      hitobject.hitpoint = hitsphere.hitpoint;
	      hitobject.id = i;
	      hitobject.object = &(s[i]);
	    }
	  hitsphere.hit = false;				
	}
    }
  return hitobject;
} 

/** Returns a class hitdata which has important information of the closest object hit. Function calls hitPlane for each plane in the list **/
hitdata findClosestPlane(vector<Plane>& p,ray& r)
{
  hitdata hitplane;
  hitdata hitobject;

  for (unsigned int i=0;i<p.size();i++)
    {
      p[i].hitPlane(r,hitplane);
      if (hitplane.hit == true)			
	{
	  if(hitobject.hit == false || hitobject.dist > hitplane.dist) //Conditions for this to be closest sphere found so far
	    {
	      hitobject.hit = true;
	      hitobject.dist = hitplane.dist;
	      hitobject.hitpoint = hitplane.hitpoint;
	      hitobject.id = i;
	      hitobject.object = &(p[i]);
	    }
	  hitplane.hit = false;				
	}
    }
  return hitobject;
} 

/** Calculates the normals to the faces of the unit cell order is v_a x v_b : v_b x v_c : v_a x v_c as unit vectors**/
void calcPlanesToCell(ATOM_NETWORK *cell,vector<Plane> &p)
{
  XYZ norm1 = cell->v_a.cross(cell->v_b).unit();
  XYZ norm2 = cell->v_c.cross(cell->v_a).unit();
  XYZ norm3 = cell->v_b.cross(cell->v_c).unit();
  
  Plane face;
  //Each face with abc point 0,0,0
  face.point = cell->abc_to_xyz(Point(0,0,0));
  
  face.normal = Point(norm1.x,norm1.y,norm1.z);
  p.push_back(face);
  
  face.normal = Point(norm2.x,norm2.y,norm2.z);
  p.push_back(face);
  
  face.normal = Point(norm3.x,norm3.y,norm3.z);
  p.push_back(face);
  
  //Each face with point abc point 1,1,1 (Notice that normals are reversed this is so that a vector can easily be determined if
  //it hit within the unit cell
  face.point = cell->abc_to_xyz(Point(1,1,1));
  
  face.normal = Point(-norm1.x,-norm1.y,-norm1.z);
  p.push_back(face);
  
  face.normal = Point(-norm2.x,-norm2.y,-norm2.z);
  p.push_back(face);
  
  face.normal = Point(-norm3.x,-norm3.y,-norm3.z);
  p.push_back(face);
}

/** Since all the opperations in Zeo++ are done within one unit cell it is important to duplicate the atoms that would overlap with the cell walls on the other side **/
void duplicateSpheresOnFace(ATOM_NETWORK *cell,vector<Sphere>& spheres,vector<Plane>& faces)
{
  Sphere temp_sphere;
  bool touch_face[6];
  
  vector<double> permutation_a;
  vector<double> permutation_b;
  vector<double> permutation_c;

  unsigned int size = spheres.size(); //because I will be adding spheres to list (there is a smarter way to do this...)
  for (unsigned int i=0; i<size; i++)
    {
      //Assemble touch_face[] for which faces it touches
      for (unsigned int j=0; j<faces.size(); j++) 
	{
	  if (faces[j].distToPlane(spheres[i].center) < (spheres[i].radius))
	    {
	      touch_face[j] = true;
	    }
	  else
	    {
	      touch_face[j] = false;
	    }
	}    
      
      //Add allowable permutations on the unitcell
      permutation_a.clear(); permutation_a.push_back(0.0);
      permutation_b.clear(); permutation_b.push_back(0.0);
      permutation_c.clear(); permutation_c.push_back(0.0);
      
      if (touch_face[0] == true)
	{
	  permutation_c.push_back(1.0);
	}
      if (touch_face[1] == true)
	{
	  permutation_b.push_back(1.0);
	}
      if (touch_face[2] == true)
	{
	  permutation_a.push_back(1.0);
	}
      if (touch_face[3] == true)
	{
	  permutation_c.push_back(-1.0);
	}
      if (touch_face[4] == true)
	{
	  permutation_b.push_back(-1.0);
	}
      if (touch_face[5] == true)
	{
	  permutation_a.push_back(-1.0);
	}
      
      //Iterate through all the permutations on the sphere
      temp_sphere = spheres[i];
      for (unsigned int a=0; a<permutation_a.size(); a++)
	{
	  for (unsigned int b=0; b<permutation_b.size(); b++)
	    {
	      for (unsigned int c=0; c<permutation_c.size(); c++)
		{
		  if (permutation_a[a] != 0 || permutation_b[b] != 0 || permutation_c[c] != 0)
		    {
		      temp_sphere.center = spheres[i].center;
		      temp_sphere.center = temp_sphere.center + cell->abc_to_xyz(Point(permutation_a[a],permutation_b[b],permutation_c[c]));
		      spheres.push_back(temp_sphere);
		    }
		}
	    }
	}
    }
}

/** This function converts nodes to spheres it allows for a more flexible implementation of ray_tracing thus making higher abstractions easier we are only concerned with the accessible nodes thus only these will be copied **/
void convertNodeToSphere(VORONOI_NETWORK& vornet,vector<Sphere>& nodes, vector<bool>& accessInfo)
{
  Sphere temp_sphere;
  VOR_NODE temp_vornode;
  
  nodes.clear();
#if DEBUG
  cout << "Adding: " << vornet.nodes.size() << " nodes (if all accessible)." << endl;
#endif
  // for each node check if it intersects the unitcell
  for (unsigned int i=0; i<vornet.nodes.size(); i++)
    {
      if (accessInfo[i] ==true)
	{
	  temp_vornode = vornet.nodes[i];
	  temp_sphere.center = Point(temp_vornode.x,temp_vornode.y,temp_vornode.z);
	  temp_sphere.radius = temp_vornode.rad_stat_sphere;
	  nodes.push_back(temp_sphere);
	}
    }
}

/** This function converts atoms to spheres it allows for a more flexible implementation of ray_tracing thus making higher abstractions easier **/
void convertAtomToSphere(ATOM_NETWORK *cell,vector<Sphere>& atoms)
{
  Sphere temp_sphere;
  ATOM temp_atom;
  
  atoms.clear();
#if DEBUG
  cout << "Adding: " << cell->atoms.size() << " atoms." << endl;
#endif
  // for each node check if it intersects the unitcell
  for (unsigned int i=0; i<cell->atoms.size(); i++)
    {
      temp_atom = cell->atoms[i];
      temp_sphere.center = Point(temp_atom.x,temp_atom.y,temp_atom.z);
      temp_sphere.radius = temp_atom.radius;
      atoms.push_back(temp_sphere);
    }
}

/** This function finds a nodes radius that encompasses the point raypoint. It will not count the previous node that the raypoint was in by using the reference id**/
bool findSphereOfPoint(Point p,vector<Sphere>& spheres,int& id)
{
  for (unsigned int i=0; i<spheres.size(); i++){
    if ((calcEuclideanDistance(p,spheres[i].center) < spheres[i].radius) && (int)i != id){ 
      id = i;			
      return true;
    }
  }
  return false;
}

/**Generates a random vector direction and unitizes it. Since we want a random ray orientation in cartesian coordinates, sphereical coordinates are used to create a point **/
Point genRandomVec()
{
  // Randomly sample point on a sphere of radius 1
  double theta = (rand()*1.0/RAND_MAX)*2*PI;
  double cosphi = 1.0 - (rand()*1.0/RAND_MAX)*2.0;
  double phi = acos(cosphi);
  
  // Convert spherical coordinates to xyz coordinates
  double x = sin(phi)*cos(theta); //Don't need abs(sin(phi)) becase phi is from 0 to PI
  double y = sin(phi)*sin(theta);
  double z = cosphi;
  Point temp(x,y,z);
  if (temp.magnitude() == 0){
    temp = genRandomVec();
  }
  return temp.unit();
}

/** This function generates a random point within the unit cell with values 0-1 for a,b, and c **/
Point genRandomPoint()
{
  // Randomly sample ray across the unit cell
  double aPoint = (rand()*1.0)/RAND_MAX;
  double bPoint = (rand()*1.0)/RAND_MAX;
  double cPoint = (rand()*1.0)/RAND_MAX;
  
  return Point(aPoint,bPoint,cPoint);
}

/** This function will travel within a sphere until it hits a point that is no longer within spheres. Returns hitdata on the point it hits**/
void rayTraceInsideSphere(ATOM_NETWORK *cell,vector<Sphere>& spheres,ray r, hitdata& hitobject)
{
  
#if DEBUG
  cout << "Current Ray Length: " << hitobject.dist << endl;
#endif
  
  if (hitobject.dist > MAXRAYDIST){
    return;
  }
  
  if (findSphereOfPoint(r.base,spheres,hitobject.id) == false) //This function also updates the id of the new sphere in hitsphere
    {
      return; //No sphere found that included point.
    }
  hitobject.object = &(spheres[hitobject.id]); //set the object*

  //If inside sphere this it should hit the sphere
  hitdata hitsphere; 
  spheres[hitobject.id].hitSphere(r,hitsphere);
  
  if (hitsphere.hit == false){
    //It is possible that point is an epsilon outside lets check
    if (fabs(calcEuclideanDistance(r.base,spheres[hitobject.id].center) - spheres[hitobject.id].radius) < threshold)
      {
	//It was round-off error so lets continue
	Sphere epsilon_sphere = spheres[hitobject.id];
	if (((r.base - spheres[hitobject.id].center) * r.vector) > 0)
	  { //vector is headed away from sphere
	    epsilon_sphere.radius += threshold;
	  }
	else
	  { // vector is headed into the sphere	    
	    epsilon_sphere.radius += -threshold;
	  }
	epsilon_sphere.hitSphere(r,hitsphere);
      }
    if (hitsphere.hit == false) //Epsilon Case has been covered so thid should not happen
      {
      	cerr << "Error: Ray did not hit a sphere. findSphereOfPoint said that one was within the radius. And it failed to be within an epsilon. Output in vmd style for easy visulization" << endl;
	cout << "draw sphere {" << spheres[hitobject.id].center << "} radius " << spheres[hitobject.id].radius << " resolution 10\n";
	Point p = r.base+r.vector.scale(10);
	cout << "draw line {" << r.base << "} {" << p << "}\n";
	abort();
      }
  }
  
  hitobject.hit = hitsphere.hit;
  hitobject.dist += hitsphere.dist; //Updates the distance traveled in spheres
  r.base = cell->shiftXYZInUC(hitsphere.hitpoint); //Move back to unit cell and set values
  hitobject.hitpoint = hitsphere.hitpoint;
  
  rayTraceInsideSphere(cell,spheres,r,hitobject);
  return;
}

/** This is a recursive function that will keep on extending a ray unitl it hits an atom or exceeds the maximum distance allowed (MAXRAYDIST)**/
void rayTraceToSphere(ATOM_NETWORK *cell,vector<Sphere>& spheres,ray r, vector<Plane>& faces,hitdata& hitobject)
{
  
#if DEBUG
  cout << "Current Ray Length: " << hitobject.dist << endl;
#endif
  
  if (hitobject.dist > MAXRAYDIST){
    return;
  }
  
  

  //Find the closest atom that it will hit
  hitdata hitsphere = findClosestSphere(spheres,r);
  
  //If I hit an Sphere great!
  if (hitsphere.hit == true)
    {
      if (hitsphere.id == hitobject.id)
	{
	  //cerr << "Just hit same atom this is not allowable" << endl;
	}
      hitobject.hit = true;
      hitobject.hitpoint = hitsphere.hitpoint;
      hitobject.dist += hitsphere.dist;
      hitobject.id = hitsphere.id;
      hitobject.object = hitsphere.object;
      return;
    }
  //If I didn't hit a Sphere then I must hit a unitcell face
  hitsphere = findClosestPlane(faces,r);
  
  assert(hitsphere.hit==true); //Thus it hit a plane

  hitobject.hit = true;
  hitobject.hitpoint = hitsphere.hitpoint;
  hitobject.dist += hitsphere.dist;
  hitobject.id = -1;
  hitobject.object = NULL;
  
  //I am moving the point just slightly outside of the unit cell to be shifted back in
  r.base = cell->shiftXYZInUC(hitobject.hitpoint+r.vector.scale(threshold));
  
  //this is a recursive call to get the total length of the ray
  rayTraceToSphere(cell,spheres,r,faces,hitobject);
  return;
}

/** Returns data from shooting a bunch of rays from accesible regions in the network until they hit an atom. Since the cell is periodic a max ray length can be set so ray does not shoot forever*/
void calcRaysInAV(ATOM_NETWORK *hiaccatmnet, ATOM_NETWORK *orgatmnet, bool highAccuracy, double r_probe_chan,double r_probe,int numSamples, ostream &output, bool visualize,string option){

  ATOM_NETWORK *cell;

  if(highAccuracy) cell=orgatmnet; else cell = hiaccatmnet; // this ensures that the loops below loop over original atoms
                                                                // and not "high-acuracy" clusters, which are used only for
                                                                // voronoi decomposition

  // Create an object that handles analysis of accessibility of sampled points
  AccessibilityClass accessAnalysis;
  if(highAccuracy) accessAnalysis.setupAndFindChannels(hiaccatmnet, orgatmnet, highAccuracy, r_probe_chan, r_probe);
    else accessAnalysis.setupAndFindChannels(hiaccatmnet, hiaccatmnet, highAccuracy, r_probe_chan, r_probe);


  srand(randSeed); //randSeed can be set in network.h
  
  //Initialize Ray Trace World
  vector<Sphere> atoms;
  convertAtomToSphere(cell,atoms);
  vector<Sphere> nodes;
  convertNodeToSphere(accessAnalysis.vornet,nodes,accessAnalysis.accessInfo);
  vector<Plane> faces;
  calcPlanesToCell(cell,faces);
  
  //Duplicate the nodes and atoms that intersect with the plane faces	
  duplicateSpheresOnFace(cell,atoms,faces);
  duplicateSpheresOnFace(cell,nodes,faces);	
  
  //Lists of Accessible, Inaccessible, and Resampled Rays
  vector<ray> axsray; //List of successful rays that were shot
  vector<ray> inaxsray;
  vector<ray> resampledray;  // List of resampled rays

#if DEBUG
  //List Atoms, Nodes, and Planes
  cout << "List of Atoms: " << endl;
  for (unsigned int i=0; i<atoms.size(); i++)
    {
      cout << "Center: " << atoms[i].center << " Radius: " << atoms[i].radius << endl;;
    }
  cout << "List of Nodes: " << endl;
  for (unsigned int i=0; i<nodes.size(); i++)
    {
      cout << "Center: " << nodes[i].center << " Radius: " << nodes[i].radius << endl;;
    }
  cout << "List of Cell Wall Faces: " << endl;
  for (unsigned int i=0; i<faces.size(); i++)
    {
      cout << "Point: " << faces[i].point << " Normal: " << faces[i].normal << endl;;
    }
#endif
  
  //Begin Implementaitons of Ray Tracing Algorithm
  cout << "Begin Ray Tracing Analysis: " << endl;
  if (option.compare("atom") == 0)
    {
      cout << "Atom Implementation Chosen: " << endl;
      cout << "Number of Samples: " << numSamples << endl;
      ray rand_ray;
      ray opposite_rand_ray;
      hitdata hitobject;
      for (int i=0; i<numSamples; i++)
	{
	  rand_ray.base =  cell->abc_to_xyz(genRandomPoint());
	  rand_ray.vector = genRandomVec();
	  
	  bool point_accessible = accessAnalysis.isVPointAccessible(rand_ray.base); 
	  
	  if (point_accessible == false)
	    {
              rand_ray.vector = Point(0,0,0);
	      inaxsray.push_back(rand_ray);
	    }
	  else //point is accessible
	    {
	      opposite_rand_ray.base = rand_ray.base;
	      opposite_rand_ray.vector  = rand_ray.vector.scale(-1);
	      
	      rayTraceToSphere(cell,atoms,rand_ray,faces,hitobject);
	      rand_ray.vector = rand_ray.vector.scale(hitobject.dist);
	      
              hitobject.dist = 0.0;
              hitobject.hit = false;
	      hitobject.object = NULL;
	      hitobject.id = -1;
	      hitobject.hitpoint = Point(0,0,0);
	      
	      //Test negative of ray
	      rayTraceToSphere(cell,atoms,opposite_rand_ray,faces,hitobject);
	      opposite_rand_ray.vector = opposite_rand_ray.vector.scale(hitobject.dist);
	      
	      rand_ray.base = rand_ray.base + opposite_rand_ray.vector;
	      rand_ray.vector = rand_ray.vector - opposite_rand_ray.vector;
	      
	      axsray.push_back(rand_ray);
	      
	      hitobject.dist = 0.0;
              hitobject.hit = false;
	      hitobject.object = NULL;
	      hitobject.id = -1;
	      hitobject.hitpoint = Point(0,0,0);
	      
	    }
	}
    }

  if (option.compare("node") == 0)
    {
      cout << "Node Implementation Chosen: " << endl;
      cout << "Number of Samples: " << numSamples << endl;
      int ndx = 0;
      ray rand_ray;
      ray opposite_rand_ray;
      hitdata hitobject;
      for (int i=0; i<numSamples; i++)
	{
	  do {
	    ndx++;
	    if (ndx >= accessAnalysis.accessInfo.size())
	      {
		ndx = 0;
	      }
	  } while(accessAnalysis.accessInfo[ndx] != true);
	  rand_ray.base = nodes[ndx].center;  
	  rand_ray.vector = genRandomVec();
	  
	  opposite_rand_ray.base = rand_ray.base;
	  opposite_rand_ray.vector  = rand_ray.vector.scale(-1);
	  
	  rayTraceToSphere(cell,atoms,rand_ray,faces,hitobject);
	  rand_ray.vector = rand_ray.vector.scale(hitobject.dist);
	  
	  hitobject.dist = 0.0;
	  hitobject.hit = false;
	  hitobject.object = NULL;
	  hitobject.id = -1;
	  hitobject.hitpoint = Point(0,0,0);
	  
	  //Test negative of ray
	  rayTraceToSphere(cell,atoms,opposite_rand_ray,faces,hitobject);
	  opposite_rand_ray.vector = opposite_rand_ray.vector.scale(hitobject.dist);
	  
	  rand_ray.base = rand_ray.base + opposite_rand_ray.vector;
	  rand_ray.vector = rand_ray.vector - opposite_rand_ray.vector;
	  
	  axsray.push_back(rand_ray);
	  
	  hitobject.dist = 0.0;
	  hitobject.hit = false;
	  hitobject.object = NULL;
	  hitobject.id = -1;
	  hitobject.hitpoint = Point(0,0,0);
	}
    }

  
  //Does not require to be found out if point is accessible
  if (option.compare("sphere") == 0)
    {
      cout << "Sphere Implementation Chosen: " << endl;
      cout << "Number of Samples: " << numSamples << endl;
      ray rand_ray;
      ray opposite_rand_ray;
      hitdata hitobject;
      for (int i=0; i<numSamples; i++)
	{
	  rand_ray.base =  cell->abc_to_xyz(genRandomPoint());
	  rand_ray.vector = genRandomVec();
	  opposite_rand_ray.base = rand_ray.base;
	  opposite_rand_ray.vector  = rand_ray.vector.scale(-1);
	      
	  rayTraceInsideSphere(cell,nodes,rand_ray,hitobject);
	  
	  if (hitobject.hit == false)
	    {
	      rand_ray.vector = Point(0,0,0);
	      inaxsray.push_back(rand_ray);
	    }
	  else
	    { //It was traveling within a sphere    
	      rand_ray.vector = rand_ray.vector.scale(hitobject.dist);

              hitobject.dist = 0.0;
              hitobject.hit = false;
	      hitobject.id = -1;
	      hitobject.object = NULL;
	      hitobject.hitpoint = Point(0,0,0);

	      rayTraceInsideSphere(cell,nodes,opposite_rand_ray,hitobject);
	      opposite_rand_ray.vector = opposite_rand_ray.vector.scale(hitobject.dist);

	      rand_ray.base = rand_ray.base + opposite_rand_ray.vector;
	      rand_ray.vector = rand_ray.vector - opposite_rand_ray.vector;

	      axsray.push_back(rand_ray);

	      hitobject.dist = 0.0;
              hitobject.hit = false;
	      hitobject.object = NULL;
	      hitobject.id = -1;
	      hitobject.hitpoint = Point(0,0,0);
	    }
	}
    }
  
  //Does not require to be found out if point is accessible
  if (option.compare("andrew_sphere") == 0)
    {
      cout << "Andrew's Sphere Implementation Chosen: " << endl;
      cout << "Number of Samples: " << numSamples << endl;
      ray rand_ray;
      ray opposite_rand_ray;
      ray init_ray;
      ray opposite_init_ray;
      double distance = 0.0;
      hitdata hitobject;
      bool inside_node = true;
      for (int i=0; i<numSamples; i++)
	{
	  distance = 0.0;
	  
	  rand_ray.base =  cell->abc_to_xyz(genRandomPoint());
	  rand_ray.vector = genRandomVec();
	  
	  //Ray in the opposite direction
	  opposite_rand_ray.base = rand_ray.base;
	  opposite_rand_ray.vector  = rand_ray.vector.scale(-1);
	  
	  //For this implementation I need a reference to the original ray
	  init_ray = rand_ray;
	  opposite_init_ray = opposite_rand_ray;
	  
	  if (findSphereOfPoint(rand_ray.base,nodes,hitobject.id) == false)
	    {
	      inside_node = false; //Point is not within a node
	    }
	  
	  while(distance < MAXRAYDIST) //Iterate Ray untill it reaches max distance
	    {	      
	      if (inside_node == false)
		{
		  rayTraceToSphere(cell,nodes,rand_ray,faces,hitobject);
		  
		  rand_ray.vector = rand_ray.vector.scale(hitobject.dist);
		  rand_ray.base = init_ray.base + init_ray.vector.scale(distance);
		  inaxsray.push_back(rand_ray);
		  
		  //Update Total Info
		  distance += hitobject.dist;
		  rand_ray.vector = init_ray.vector;
		  rand_ray.base = cell->shiftXYZInUC(hitobject.hitpoint+rand_ray.vector.scale(threshold));
		  
		  hitobject.dist = 0.0;
		  hitobject.hit = false;
		  hitobject.hitpoint = Point(0,0,0); 
		  hitobject.id = -1;
		  hitobject.object = NULL;
		  
		  inside_node = true; //Now inside node
		}
	      else
		{ //Shoot Through Nodes    
		  rayTraceInsideSphere(cell,nodes,rand_ray,hitobject);
		  
		  rand_ray.vector = rand_ray.vector.scale(hitobject.dist);
		  rand_ray.base = init_ray.base + init_ray.vector.scale(distance);
		  axsray.push_back(rand_ray);
		  
		  //Update Total Info
		  distance += hitobject.dist;
		  rand_ray.vector = init_ray.vector;
		  rand_ray.base = cell->shiftXYZInUC(hitobject.hitpoint+rand_ray.vector.scale(threshold));
		  
		  hitobject.dist = 0.0;
		  hitobject.hit = false;
		  hitobject.hitpoint = Point(0,0,0);  
		  //Notice that I am not resetting hitobject.id and setting object
		  hitobject.object = &(nodes[hitobject.id]);
		  
		  inside_node = false; //Now outside of nodes
		}
	    }	  
	}
      cout << "Ray Andrew Sphere Implementation Completed:" << endl;
    }
  
  
  //Does not require to be found out if point is accessible
  if (option.compare("andrew_atom") == 0)
    {
      cout << "Andrew's Atom Implementation Chosen: " << endl;
      cout << "Number of Samples: " << numSamples << endl;
      ray rand_ray;
      ray opposite_rand_ray;
      ray init_ray;
      ray opposite_init_ray;
      double distance = 0.0;
      hitdata hitobject;
      bool inside_atom = false;
      bool inside_accessible_region;
      for (int i=0; i<numSamples; i++)
	{
	  distance = 0.0;
	  
	  rand_ray.base =  cell->abc_to_xyz(genRandomPoint());
	  rand_ray.vector = genRandomVec();
	  
	  //Ray in the opposite direction
	  opposite_rand_ray.base = rand_ray.base;
	  opposite_rand_ray.vector  = rand_ray.vector.scale(-1);
	  
	  //For this implementation I need a reference to the original ray
	  init_ray = rand_ray;
	  opposite_init_ray = opposite_rand_ray;
	  
	  if ((findSphereOfPoint(rand_ray.base,atoms,hitobject.id) == true))
	    {
	      inside_atom = true; //Point is within an atom
	    }
	  else //Point is not within an atom
	    {
	      inside_atom = false;
	    }
	  inside_accessible_region = accessAnalysis.isVPointAccessible(rand_ray.base);
	  
	  while(distance < MAXRAYDIST) //Iterate Ray untill it reaches max distance
	    {	      
	      if ((inside_atom == false) && (inside_accessible_region == true))
		{
		  rayTraceToSphere(cell,atoms,rand_ray,faces,hitobject);
		  
		  rand_ray.vector = rand_ray.vector.scale(hitobject.dist);
		  //rand_ray.base = init_ray.base + init_ray.vector.scale(distance);
		  axsray.push_back(rand_ray);
		  
		  //Update Total Info
		  distance += hitobject.dist;
		  rand_ray.vector = init_ray.vector;
		  rand_ray.base = cell->shiftXYZInUC(hitobject.hitpoint+rand_ray.vector.scale(threshold));
		  //Reset all hitobject data
		  hitobject.dist = 0.0;
		  hitobject.hit = false;
		  hitobject.hitpoint = Point(0,0,0); 
		  hitobject.id = -1;
		  hitobject.object = NULL;
		  
		  inside_atom = true; //Now inside atom
		}
	      else if ((inside_atom == false) && (inside_accessible_region == false))
		{ //not within an atom but in inaccessible region
		  rayTraceToSphere(cell,atoms,rand_ray,faces,hitobject);
		  
		  rand_ray.vector = rand_ray.vector.scale(hitobject.dist);
		  //rand_ray.base = init_ray.base + init_ray.vector.scale(distance);
		  inaxsray.push_back(rand_ray);
		  
		  //Update Total Info
		  distance += hitobject.dist;
		  rand_ray.vector = init_ray.vector;
		  rand_ray.base = cell->shiftXYZInUC(hitobject.hitpoint+init_ray.vector.scale(threshold));
		  //Reset all hitobject data
		  hitobject.dist = 0.0;
		  hitobject.hit = false;
		  hitobject.hitpoint = Point(0,0,0); 
		  hitobject.object = NULL;
		  hitobject.id = -1;
		  
		  inside_atom = true; //Now inside atom
		}
	      else // Ray is inside of atoms
		{    
		  rayTraceInsideSphere(cell,atoms,rand_ray,hitobject);
		  		  
		  rand_ray.vector = rand_ray.vector.scale(hitobject.dist);
		  //rand_ray.base = init_ray.base + init_ray.vector.scale(distance);
		  //This is where I would push data on traveling through atoms  
		  //but it is not here yet
   
		  //Update Total Info
		  distance += hitobject.dist;
		  rand_ray.vector = init_ray.vector;
		  rand_ray.base = cell->shiftXYZInUC(hitobject.hitpoint+init_ray.vector.scale(threshold));
		 
		  //Reset all hitobject data besides the last object id hit
		  hitobject.dist = 0.0;
		  hitobject.hit = false;
		  hitobject.hitpoint = Point(0,0,0);  
		  //Notice that I am not resetting hitobject.id and setting object
		  hitobject.object = &(nodes[hitobject.id]);
		  
		  //So at this point I know that I am outside an atom however am
		  //I in an accesible region or not?
		  inside_atom = false; //Now outside of atoms and in accessible region
		  inside_accessible_region = accessAnalysis.isVPointAccessible(rand_ray.base);
		}
	    }	  
	}
      cout << "Ray Andrew Atom Implementation Completed:" << endl;
    }
  
  
  //Post Process Information
  if(visualize)
    {
      reportRays(output, axsray, inaxsray, true);
      //reportAtoms(output,atoms);
      //reportNodes(output,nodes);
    }
  else //no visulization then output a histogram
    {
      //This parameter can be changed in network.h
      reportHistogram(output,BINSIZE,MAXBINS,axsray);
    }
  reportRayInfo(axsray); // A more detailed output of ray info from program

  //FREE MEMORY
  accessAnalysis.deconstruct();
}

/* backup of calcRayinAV

// Returns data from shooting a bunch of rays from accesible regions in the network until they hit an atom. Since the cell is periodic a max ray length can be set so ray does not shoot forever
void calcRaysInAV(ATOM_NETWORK *cell, ATOM_NETWORK *orgcell, bool highAccuracy, double r_probe_chan,double r_probe,int numSamples, ostream &output, bool visualize,string option){
  // Create a temporary copy of the atomic network in which each atom's radius has been increased by the probe radius
  ATOM_NETWORK newAtomNet;
  cell->copy(&newAtomNet);
  for(int i = 0; i < newAtomNet.numAtoms; i++){ newAtomNet.atoms[i].radius += r_probe; }  
  
  // Calculate and store the Voronoi network for this new atomic network
  VORONOI_NETWORK vornet;
  vector<BASIC_VCELL> vorcells;
  vector<VOR_CELL> advCells;
  
  container_periodic_poly *new_rad_con = (container_periodic_poly *)performVoronoiDecomp(true, &newAtomNet, &vornet, advCells, false, vorcells);
  
  vector<CHANNEL> channels = vector<CHANNEL>();
  vector<bool> accessInfo = vector<bool> ();
  CHANNEL::findChannels(&vornet, max(0.0, r_probe_chan - r_probe), &accessInfo, &channels);
  srand(randSeed); //randSeed can be set in network.h
  
  //Initialize Ray Trace World
  vector<Sphere> atoms;
  convertAtomToSphere(&newAtomNet,atoms);
  vector<Sphere> nodes;
  convertNodeToSphere(vornet,nodes,accessInfo);
  vector<Plane> faces;
  calcPlanesToCell(cell,faces);
  
  //Duplicate the nodes and atoms that intersect with the plane faces	
  duplicateSpheresOnFace(cell,atoms,faces);
  duplicateSpheresOnFace(cell,nodes,faces);	
  
  //Lists of Accessible, Inaccessible, and Resampled Rays
  vector<ray> axsray; //List of successful rays that were shot
  vector<ray> inaxsray;
  vector<ray> resampledray;  // List of resampled rays

#if DEBUG
  //List Atoms, Nodes, and Planes
  cout << "List of Atoms: " << endl;
  for (unsigned int i=0; i<atoms.size(); i++)
    {
      cout << "Center: " << atoms[i].center << " Radius: " << atoms[i].radius << endl;;
    }
  cout << "List of Nodes: " << endl;
  for (unsigned int i=0; i<nodes.size(); i++)
    {
      cout << "Center: " << nodes[i].center << " Radius: " << nodes[i].radius << endl;;
    }
  cout << "List of Cell Wall Faces: " << endl;
  for (unsigned int i=0; i<faces.size(); i++)
    {
      cout << "Point: " << faces[i].point << " Normal: " << faces[i].normal << endl;;
    }
#endif
  
  //Begin Implementaitons of Ray Tracing Algorithm
  cout << "Begin Ray Tracing Analysis: " << endl;
  if (option.compare("atom") == 0)
    {
      cout << "Atom Implementation Chosen: " << endl;
      cout << "Number of Samples: " << numSamples << endl;
      ray rand_ray;
      ray opposite_rand_ray;
      hitdata hitobject;
      for (int i=0; i<numSamples; i++)
	{
	  rand_ray.base =  cell->abc_to_xyz(genRandomPoint());
	  rand_ray.vector = genRandomVec();
	  
	  bool point_accessible = accessiblePoint(rand_ray.base,new_rad_con,r_probe,cell,vorcells,accessInfo);				
	  
	  if (point_accessible == false)
	    {
              rand_ray.vector = Point(0,0,0);
	      inaxsray.push_back(rand_ray);
	    }
	  else //point is accessible
	    {
	      opposite_rand_ray.base = rand_ray.base;
	      opposite_rand_ray.vector  = rand_ray.vector.scale(-1);
	      
	      rayTraceToSphere(cell,atoms,rand_ray,faces,hitobject);
	      rand_ray.vector = rand_ray.vector.scale(hitobject.dist);
	      
              hitobject.dist = 0.0;
              hitobject.hit = false;
	      hitobject.object = NULL;
	      hitobject.id = -1;
	      hitobject.hitpoint = Point(0,0,0);
	      
	      //Test negative of ray
	      rayTraceToSphere(cell,atoms,opposite_rand_ray,faces,hitobject);
	      opposite_rand_ray.vector = opposite_rand_ray.vector.scale(hitobject.dist);
	      
	      rand_ray.base = rand_ray.base + opposite_rand_ray.vector;
	      rand_ray.vector = rand_ray.vector - opposite_rand_ray.vector;
	      
	      axsray.push_back(rand_ray);
	      
	      hitobject.dist = 0.0;
              hitobject.hit = false;
	      hitobject.object = NULL;
	      hitobject.id = -1;
	      hitobject.hitpoint = Point(0,0,0);
	      
	    }
	}
    }

  if (option.compare("node") == 0)
    {
      cout << "Node Implementation Chosen: " << endl;
      cout << "Number of Samples: " << numSamples << endl;
      int ndx = 0;
      ray rand_ray;
      ray opposite_rand_ray;
      hitdata hitobject;
      for (int i=0; i<numSamples; i++)
	{
	  do {
	    ndx++;
	    if (ndx >= accessInfo.size())
	      {
		ndx = 0;
	      }
	  } while(accessInfo[ndx] != true);
	  rand_ray.base = nodes[ndx].center;  
	  rand_ray.vector = genRandomVec();
	  
	  opposite_rand_ray.base = rand_ray.base;
	  opposite_rand_ray.vector  = rand_ray.vector.scale(-1);
	  
	  rayTraceToSphere(cell,atoms,rand_ray,faces,hitobject);
	  rand_ray.vector = rand_ray.vector.scale(hitobject.dist);
	  
	  hitobject.dist = 0.0;
	  hitobject.hit = false;
	  hitobject.object = NULL;
	  hitobject.id = -1;
	  hitobject.hitpoint = Point(0,0,0);
	  
	  //Test negative of ray
	  rayTraceToSphere(cell,atoms,opposite_rand_ray,faces,hitobject);
	  opposite_rand_ray.vector = opposite_rand_ray.vector.scale(hitobject.dist);
	  
	  rand_ray.base = rand_ray.base + opposite_rand_ray.vector;
	  rand_ray.vector = rand_ray.vector - opposite_rand_ray.vector;
	  
	  axsray.push_back(rand_ray);
	  
	  hitobject.dist = 0.0;
	  hitobject.hit = false;
	  hitobject.object = NULL;
	  hitobject.id = -1;
	  hitobject.hitpoint = Point(0,0,0);
	}
    }

  
  //Does not require to be found out if point is accessible
  if (option.compare("sphere") == 0)
    {
      cout << "Sphere Implementation Chosen: " << endl;
      cout << "Number of Samples: " << numSamples << endl;
      ray rand_ray;
      ray opposite_rand_ray;
      hitdata hitobject;
      for (int i=0; i<numSamples; i++)
	{
	  rand_ray.base =  cell->abc_to_xyz(genRandomPoint());
	  rand_ray.vector = genRandomVec();
	  opposite_rand_ray.base = rand_ray.base;
	  opposite_rand_ray.vector  = rand_ray.vector.scale(-1);
	      
	  rayTraceInsideSphere(cell,nodes,rand_ray,hitobject);
	  
	  if (hitobject.hit == false)
	    {
	      rand_ray.vector = Point(0,0,0);
	      inaxsray.push_back(rand_ray);
	    }
	  else
	    { //It was traveling within a sphere    
	      rand_ray.vector = rand_ray.vector.scale(hitobject.dist);

              hitobject.dist = 0.0;
              hitobject.hit = false;
	      hitobject.id = -1;
	      hitobject.object = NULL;
	      hitobject.hitpoint = Point(0,0,0);

	      rayTraceInsideSphere(cell,nodes,opposite_rand_ray,hitobject);
	      opposite_rand_ray.vector = opposite_rand_ray.vector.scale(hitobject.dist);

	      rand_ray.base = rand_ray.base + opposite_rand_ray.vector;
	      rand_ray.vector = rand_ray.vector - opposite_rand_ray.vector;

	      axsray.push_back(rand_ray);

	      hitobject.dist = 0.0;
              hitobject.hit = false;
	      hitobject.object = NULL;
	      hitobject.id = -1;
	      hitobject.hitpoint = Point(0,0,0);
	    }
	}
    }
  
  //Does not require to be found out if point is accessible
  if (option.compare("andrew_sphere") == 0)
    {
      cout << "Andrew's Sphere Implementation Chosen: " << endl;
      cout << "Number of Samples: " << numSamples << endl;
      ray rand_ray;
      ray opposite_rand_ray;
      ray init_ray;
      ray opposite_init_ray;
      double distance = 0.0;
      hitdata hitobject;
      bool inside_node = true;
      for (int i=0; i<numSamples; i++)
	{
	  distance = 0.0;
	  
	  rand_ray.base =  cell->abc_to_xyz(genRandomPoint());
	  rand_ray.vector = genRandomVec();
	  
	  //Ray in the opposite direction
	  opposite_rand_ray.base = rand_ray.base;
	  opposite_rand_ray.vector  = rand_ray.vector.scale(-1);
	  
	  //For this implementation I need a reference to the original ray
	  init_ray = rand_ray;
	  opposite_init_ray = opposite_rand_ray;
	  
	  if (findSphereOfPoint(rand_ray.base,nodes,hitobject.id) == false)
	    {
	      inside_node = false; //Point is not within a node
	    }
	  
	  while(distance < MAXRAYDIST) //Iterate Ray untill it reaches max distance
	    {	      
	      if (inside_node == false)
		{
		  rayTraceToSphere(cell,nodes,rand_ray,faces,hitobject);
		  
		  rand_ray.vector = rand_ray.vector.scale(hitobject.dist);
		  rand_ray.base = init_ray.base + init_ray.vector.scale(distance);
		  inaxsray.push_back(rand_ray);
		  
		  //Update Total Info
		  distance += hitobject.dist;
		  rand_ray.vector = init_ray.vector;
		  rand_ray.base = cell->shiftXYZInUC(hitobject.hitpoint+rand_ray.vector.scale(threshold));
		  
		  hitobject.dist = 0.0;
		  hitobject.hit = false;
		  hitobject.hitpoint = Point(0,0,0); 
		  hitobject.id = -1;
		  hitobject.object = NULL;
		  
		  inside_node = true; //Now inside node
		}
	      else
		{ //Shoot Through Nodes    
		  rayTraceInsideSphere(cell,nodes,rand_ray,hitobject);
		  
		  rand_ray.vector = rand_ray.vector.scale(hitobject.dist);
		  rand_ray.base = init_ray.base + init_ray.vector.scale(distance);
		  axsray.push_back(rand_ray);
		  
		  //Update Total Info
		  distance += hitobject.dist;
		  rand_ray.vector = init_ray.vector;
		  rand_ray.base = cell->shiftXYZInUC(hitobject.hitpoint+rand_ray.vector.scale(threshold));
		  
		  hitobject.dist = 0.0;
		  hitobject.hit = false;
		  hitobject.hitpoint = Point(0,0,0);  
		  //Notice that I am not resetting hitobject.id and setting object
		  hitobject.object = &(nodes[hitobject.id]);
		  
		  inside_node = false; //Now outside of nodes
		}
	    }	  
	}
      cout << "Ray Andrew Sphere Implementation Completed:" << endl;
    }
  
  
  //Does not require to be found out if point is accessible
  if (option.compare("andrew_atom") == 0)
    {
      cout << "Andrew's Atom Implementation Chosen: " << endl;
      cout << "Number of Samples: " << numSamples << endl;
      ray rand_ray;
      ray opposite_rand_ray;
      ray init_ray;
      ray opposite_init_ray;
      double distance = 0.0;
      hitdata hitobject;
      bool inside_atom = false;
      bool inside_accessible_region;
      for (int i=0; i<numSamples; i++)
	{
	  distance = 0.0;
	  
	  rand_ray.base =  cell->abc_to_xyz(genRandomPoint());
	  rand_ray.vector = genRandomVec();
	  
	  //Ray in the opposite direction
	  opposite_rand_ray.base = rand_ray.base;
	  opposite_rand_ray.vector  = rand_ray.vector.scale(-1);
	  
	  //For this implementation I need a reference to the original ray
	  init_ray = rand_ray;
	  opposite_init_ray = opposite_rand_ray;
	  
	  if ((findSphereOfPoint(rand_ray.base,atoms,hitobject.id) == true))
	    {
	      inside_atom = true; //Point is within an atom
	    }
	  else //Point is not within an atom
	    {
	      inside_atom = false;
	    }
	  inside_accessible_region = accessiblePoint(rand_ray.base,new_rad_con,r_probe,cell,vorcells,accessInfo); 
	  
	  while(distance < MAXRAYDIST) //Iterate Ray untill it reaches max distance
	    {	      
	      if ((inside_atom == false) && (inside_accessible_region == true))
		{
		  rayTraceToSphere(cell,atoms,rand_ray,faces,hitobject);
		  
		  rand_ray.vector = rand_ray.vector.scale(hitobject.dist);
		  //rand_ray.base = init_ray.base + init_ray.vector.scale(distance);
		  axsray.push_back(rand_ray);
		  
		  //Update Total Info
		  distance += hitobject.dist;
		  rand_ray.vector = init_ray.vector;
		  rand_ray.base = cell->shiftXYZInUC(hitobject.hitpoint+rand_ray.vector.scale(threshold));
		  //Reset all hitobject data
		  hitobject.dist = 0.0;
		  hitobject.hit = false;
		  hitobject.hitpoint = Point(0,0,0); 
		  hitobject.id = -1;
		  hitobject.object = NULL;
		  
		  inside_atom = true; //Now inside atom
		}
	      else if ((inside_atom == false) && (inside_accessible_region == false))
		{ //not within an atom but in inaccessible region
		  rayTraceToSphere(cell,atoms,rand_ray,faces,hitobject);
		  
		  rand_ray.vector = rand_ray.vector.scale(hitobject.dist);
		  //rand_ray.base = init_ray.base + init_ray.vector.scale(distance);
		  inaxsray.push_back(rand_ray);
		  
		  //Update Total Info
		  distance += hitobject.dist;
		  rand_ray.vector = init_ray.vector;
		  rand_ray.base = cell->shiftXYZInUC(hitobject.hitpoint+init_ray.vector.scale(threshold));
		  //Reset all hitobject data
		  hitobject.dist = 0.0;
		  hitobject.hit = false;
		  hitobject.hitpoint = Point(0,0,0); 
		  hitobject.object = NULL;
		  hitobject.id = -1;
		  
		  inside_atom = true; //Now inside atom
		}
	      else // Ray is inside of atoms
		{    
		  rayTraceInsideSphere(cell,atoms,rand_ray,hitobject);
		  		  
		  rand_ray.vector = rand_ray.vector.scale(hitobject.dist);
		  //rand_ray.base = init_ray.base + init_ray.vector.scale(distance);
		  //This is where I would push data on traveling through atoms  
		  //but it is not here yet
   
		  //Update Total Info
		  distance += hitobject.dist;
		  rand_ray.vector = init_ray.vector;
		  rand_ray.base = cell->shiftXYZInUC(hitobject.hitpoint+init_ray.vector.scale(threshold));
		 
		  //Reset all hitobject data besides the last object id hit
		  hitobject.dist = 0.0;
		  hitobject.hit = false;
		  hitobject.hitpoint = Point(0,0,0);  
		  //Notice that I am not resetting hitobject.id and setting object
		  hitobject.object = &(nodes[hitobject.id]);
		  
		  //So at this point I know that I am outside an atom however am
		  //I in an accesible region or not?
		  inside_atom = false; //Now outside of atoms and in accessible region
		  inside_accessible_region = accessiblePoint(rand_ray.base,new_rad_con,r_probe,cell,vorcells,accessInfo); 
		}
	    }	  
	}
      cout << "Ray Andrew Atom Implementation Completed:" << endl;
    }
  
  
  //Post Process Information
  if(visualize)
    {
      reportRays(output, axsray, inaxsray, true);
      //reportAtoms(output,atoms);
      //reportNodes(output,nodes);
    }
  else //no visulization then output a histogram
    {
      //This parameter can be changed in network.h
      reportHistogram(output,BINSIZE,MAXBINS,axsray);
    }
  reportRayInfo(axsray); // A more detailed output of ray info from program

  //FREE MEMORY
  delete new_rad_con;
}


// ends backup of RaysInAV
*/

/** Outputs a data file Ray_Info.txt that gives important information
    about the rays that were found to be accessible such as Point, Vector, and Magnitude **/
void reportRayInfo(vector<ray>& axsray){
  ofstream rayinfo;
  rayinfo.open ("Ray_Info.txt");
  if (rayinfo.good() == true){  
    cout << "Ray_Info.txt: size = " << axsray.size() << endl; 
    rayinfo << "x y z dx dy dz magnitude" << endl;
    for (unsigned int i=0;i<axsray.size();i++){
      ray rand_ray = axsray.at(i);
//      rayinfo << rand_ray.base << " "<< rand_ray.vector << " " << rand_ray.vector.magnitude()  << endl;
      rayinfo << rand_ray.base.vals[0] << " " << rand_ray.base.vals[1] << " "<< rand_ray.base.vals[2] << " "
              << rand_ray.vector.vals[0] << " " << rand_ray.vector.vals[1] << " " << rand_ray.vector.vals[2] << " " 
	      << rand_ray.vector.magnitude()  << endl;
      }
  }
  else {
    cerr << "Ray_Info.txt ran into errors opening" << endl;
    abort();
  }
  rayinfo.close();
  return;
}

/* Print the coordinates contained in the provided vectors in a manner that they can be
 * displayed using ZeoVis and its VMD interface. Accessible points are colored green while
 * inaccessible points are colored red.*/
void reportRays(ostream &output, vector<ray>& axsray, vector<ray>& inaxsray, bool colors){
  if (colors == true)
    {
      //Colors are in order of the rainbow
      output << "{color purple}" << "\n";
      for(unsigned int i = 0; i < axsray.size(); i++)
	{
	  if (axsray[i].vector.magnitude() < 3)
	    {
	      ray rand_ray = axsray.at(i);
	      Point addition = rand_ray.base + rand_ray.vector;
	      output << "{line {" << rand_ray.base << "} {" << addition << "}}" << "\n";
	    }
	}
      output << "{color blue}" << "\n";
      for(unsigned int i = 0; i < axsray.size(); i++)
	{
	  if (axsray[i].vector.magnitude() < 6 && axsray[i].vector.magnitude() >= 3)
	    {
	      ray rand_ray = axsray.at(i);
	      Point addition = rand_ray.base + rand_ray.vector;
	      output << "{line {" << rand_ray.base << "} {" << addition << "}}" << "\n";
	    }
	}
      output << "{color cyan}" << "\n";
      for(unsigned int i = 0; i < axsray.size(); i++)
	{
	  if (axsray[i].vector.magnitude() >= 6 && axsray[i].vector.magnitude() < 9)
	    {
	      ray rand_ray = axsray.at(i);
	      Point addition = rand_ray.base + rand_ray.vector;
	      output << "{line {" << rand_ray.base << "} {" << addition << "}}" << "\n";
	    }
	}
      output << "{color lime}" << "\n";
      for(unsigned int i = 0; i < axsray.size(); i++)
	{
	  if (axsray[i].vector.magnitude() >= 9 && axsray[i].vector.magnitude() < 12)
	    {
	      ray rand_ray = axsray.at(i);
	      Point addition = rand_ray.base + rand_ray.vector;
	      output << "{line {" << rand_ray.base << "} {" << addition << "}}" << "\n";
	    }
	}
      output << "{color orange}" << "\n";
      for(unsigned int i = 0; i < axsray.size(); i++)
	{
	  if (axsray[i].vector.magnitude() >= 12 && axsray[i].vector.magnitude() < 20)
	    {
	      ray rand_ray = axsray.at(i);
	      Point addition = rand_ray.base + rand_ray.vector;
	      output << "{line {" << rand_ray.base << "} {" << addition << "}}" << "\n";
	    }
	}
      output << "{color red}" << "\n";
      for(unsigned int i = 0; i < axsray.size(); i++)
	{
	  if (axsray[i].vector.magnitude() >= 20)
	    {
	      ray rand_ray = axsray.at(i);
	      Point addition = rand_ray.base + rand_ray.vector;
	      output << "{line {" << rand_ray.base << "} {" << addition << "}}" << "\n";
	    }
	}
    }
  else //colors if false
    {
      output << "{color blue}" << "\n";
      for(unsigned int i = 0; i < axsray.size(); i++){
	ray rand_ray = axsray.at(i);
	Point addition = rand_ray.base + rand_ray.vector;
	output << "{line {" << rand_ray.base << "} {" << addition << "}}" << "\n";
      }  
      output << "{color red}" << "\n";
      for(unsigned int i = 0; i < inaxsray.size(); i++){
	ray rand_ray = inaxsray.at(i);
	Point addition = rand_ray.base + rand_ray.vector;
	output << "{line {" <<  rand_ray.base << "} {" << addition << "}}" << "\n";
      }
    }
}

void reportAtoms(ostream &output, vector<Sphere>& s)
{
  output << "{color red}" << "\n";
  for (unsigned int i=0; i < s.size(); i++){
    output << "{sphere {" << s[i].center <<"} radius " << s[i].radius << " resolution 50}\n";
  }
}
void reportNodes(ostream &output, vector<Sphere>& s)
{
  output << "{color green}" << "\n";
  for (unsigned int i=0; i < s.size(); i++){
    output << "{sphere {" << s[i].center <<"} radius " << s[i].radius << " resolution 50}\n";
  }
}
/**
void reportPlane(ostream &output, vector<Plane> p){
  output << "draw color orange" << "\n";
  for (unsigned int i=0; i < p.size(); i++){
    output << "draw traingle {" <<  point << "} radius " << temp.radius << " resolution 10 \n";
  }
}
**/

void reportHistogram(ostream& output,const double binSize,const int maxBins, vector<ray>& rays)
{
  assert(binSize > threshold);
  int bins[maxBins];
  for (int i=0; i<maxBins; i++)
    {
      bins[i] = 0;
    }
  
  //Assemble Hitsogram
  int bin;
  for (unsigned int i=0; i<rays.size(); i++)
    {
      bin = rays[i].vector.magnitude()/binSize;
      if (bin >= maxBins) 
	{
	  bin = maxBins - 1;
	}
      bins[bin]++;
    }
  
  //Output Histogram
  output << "Ray Histogram - Bin Size = " << binSize << " Number of Bins: " << maxBins << "From: 0 To: " << binSize * 1.0 * maxBins<< endl;
  for (int i=0; i<maxBins; i++)
    {
      output << bins[i] << endl;
    }
}
