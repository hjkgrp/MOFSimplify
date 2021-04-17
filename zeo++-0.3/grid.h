#ifndef GRID_H
#define GRID_H

#include <cstdio>
#include <string>
#include "networkstorage.h"


/* BOV grid functions */
void write_distances(std::FILE *f, double ***grid, int x_grid_steps, int y_grid_steps, int z_grid_steps);
void write_bov(std::FILE *f, char *output_distances_name, int x_grid_steps, int y_grid_steps, int z_grid_steps, double xMin, double yMin, double zMin, double x_box_size, double y_box_size, double z_box_size);
double calculate_distance_function(ATOM_NETWORK *network, int i, int j, int k, double minX, double minY, double minZ, double xGridRes, double yGridRes, double zGridRes, int x_grid_steps, int y_grid_steps, int z_grid_steps, char gridtype);
void generateBOVGrid(ATOM_NETWORK *atmnet, std::string name_f_dist, std::string name_g_dist, std::string name_h_dist, std::string name_f_bov, std::string name_g_bov, std::string name_h_bov);

/* Gaussian cube functions */
void generateGaussianGrid(ATOM_NETWORK *atmnet, std::string cubefilename, bool angstrom_to_bohr, bool useMassFlag);

/*   Gausian cube file with accessibility information  */
void generateGaussianGridWithAccessibilityInfo(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgAtomnet, bool highAccuracy, double probe_radius, 
                                               std::string cubefilename, bool angstrom_to_bohr, bool useMassFlag);

/* Gaussian grid (cube format) for averaging */

void calculateAverageGrid(ATOM_NETWORK *atmnet, std::string inputfilename, std::string cubefilename, bool angstrom_to_bohr, bool useMassFlag);

/* Gaussian grid (cube gormat) with 3D histogram of number of occupancies (1 per frame) */
void calculateAverageGridPerFrame(ATOM_NETWORK *atmnet, std::string inputfilename, std::string cubefilename, bool angstrom_to_bohr, bool useMassFlag);

/* Gaussian grid (cube format) class 
   Multipurpose class written intially to handle averaging tasts */

class GaussianCube {

 // definition of the cube grid

 XYZ o;
 XYZ va,vb,vc; // vectors pointing in a, b, c directons of the unit cell
 XYZ stepabc; // step in abc coordinates

 int na,nb,nc; // number of points in each direction
 int gridsize; //total points

 double ***cube; // grid itself

 public:

/* conctructor */
 GaussianCube()
    {
    gridsize = 0;
    };

  GaussianCube(ATOM_NETWORK *atmnet);

/* allocate */

 void allocate(int size_a, int size_b, int size_c)
  {
  cube = new double**[size_a];
  for(int i=0; i<size_a; i++) {
     cube[i] = new double*[size_b];
     for(int j=0; j<size_b; j++) {
        cube[i][j] = new double[size_c];
        }
     }
  gridsize = na * nb *nc;

  // zeroing the grid
  for(int i=0; i<size_a; i++) 
    for(int j=0; j<size_b; j++)
      for(int k=0; k<size_c; k++)
        cube[i][j][k] = 0.0;

  };

/* deinit */
void deinit()
 {
 delete[] cube;
 };

/* functions in GaussianCube class */

// saves grid to Gaussian cube file
void writeGrid(ATOM_NETWORK *atmnet, std::string cubefilename, bool angstrom_to_bohr, bool useMassFlag);

// calculate distnace grid 
void calculateDistanceGrid(ATOM_NETWORK *atmnet);

// calculate distance grid with accessibility information
void calculateDistanceGridWithAccessibilityInfo(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgatmnet, bool highAccuracy, double probe_radius);

// load a text file with points and project onto a grid to generate 3D histogram
void loadHistogramData(std::string inputfilename);

// load a list of text file with points(frames), and project each frame onto a grid to generate 3D histogram
// the histogram will represent number of frames that had a value >0 for particular grid point
void loadHistogramDataPerFrame(std::string listfilename);

};

#endif
