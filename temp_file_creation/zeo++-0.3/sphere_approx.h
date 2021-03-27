
/* This class is handling replacment of atoms (large spheres) with clusters of small spheres
 * it is done to achive better approximation of S-Voronoi (S-cell) using radial Voronoi decomposition
 * implemented in Voro++ */
/* The class was written by Marielle Pinheiro (March, 2013), modified and merge with Zeo++ by Maciej Haranczyk */

#ifndef SPHERE_APPROX
#define SPHERE_APPROX

#include <iostream>
#include <vector>
#include "networkstorage.h"

class AtomCluster {

  private:
   double smallSphereRadius; //  defines the size of replacement sphere        
   ATOM   orgAtom; // original atom/to be replaced
  public:
	//variables
	ATOM center;
	std::vector <ATOM> atom_vector;
	//constructor
/* old to be delated
	AtomCluster(XYZ point);
	AtomCluster(double x, double y, double z); 
*/
    AtomCluster(ATOM orgatm, double replacementSphereRadii);
	//methods
	//main function
	void replaceAtomByCluster(std::string cluster_type, double large_atom_radius);
    void copyReplacementAtoms(ATOM_NETWORK *atmnet, int atomID, std::vector <ATOM> *newatoms);
  private:
	//translation functions
    ATOM translate_sphere(ATOM sphere, double x_step, double y_step, double z_step, int x_sign, int y_sign, int z_sign);

	void plusMinus_axes(double shift, int i);
	void plusMinus_xy(double x, double y, int i);
	void plusMinus_xz(double x, double z, int i);
	void plusMinus_yz(double y, double z, int i);
	void plusMinus_all(double x, double y, double z, int i);
	void translate_cube(double shift, int i);
	void rotate_xy_cube(double shift, int i);
	void rotate_xz_cube(double shift, int i);
	void rotate_yz_cube(double shift, int i);
	void translate_axes(double shift, int i);
	void translate_dodecahedron(double shift, int i);
	void rotate_icosahedron(double shift, int i);
	void translate_icosidodecahedron(double shift, int i);
	ATOM calc_center(std::vector <ATOM> old_coords, int i_0, int i_1, int i_2, int i_3, int i_4, double shift);
	void calc_centerSpheres(double shift, int i);
	void translate_rhombi(double shift, int i);
    void sphere_spiral(double n, double shift);
	void print_xyz_coords(FILE *output);
};



/* Function that analyzes distribution of atomic radii and execute replacement of large atoms 
   with clusters of small ones */
void setupHighAccuracyAtomNetwork(ATOM_NETWORK *atmnet, std::string AccSetting);


#endif
