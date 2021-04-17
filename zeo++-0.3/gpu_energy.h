#ifndef GPU_ENERGY_H
#define GPU_ENERGY_H

#include <utility>
//#include "heap.h"

/* Class used to transform between xyz coordinates and abc coordinates.*/
class UNIT_CELL {
  double va_x;
  double vb_x, vb_y;
  double vc_x, vc_y, vc_z;

  double inv_va_x;
  double inv_vb_x, inv_vb_y;
  double inv_vc_x, inv_vc_y, inv_vc_z;

public:
  /* Construct a unit cell with the provided vector components*/
  UNIT_CELL(double vax, double vbx, double vby, double vcx, double vcy, double vcz);

  /** Transform the coordinates relative to the uc vectors into ones relative to
   *  xyz coordinate system and store them using the provided pointers. */
  void abc_to_xyz(double a, double b, double c, double &x, double &y, double &z);

  /** Transform the coordinates relative to the xyz coordinate system vectors into ones relative to
   *  the unit cell vectors and store them using the provided pointers. */
  void xyz_to_abc(double x, double y, double z, double &a, double &b, double &c);
};

/* Class used to store the displacement and energy of a path when it encounters a particular node.*/
class DISP_INFO {
public:
  bool isReal;
  char a,b,c;      // Displacement components. Chars used to conserve space
  float maxEnergy; // Maximum energy barrier on path so far

  // Constructor that creates a dummy container
  DISP_INFO();

  /** Create a container with the provided displacement components and energy.
   *   Displacements are stored as characters to conserve space. */
  DISP_INFO(int myA, int myB, int myC, float myMaxEnergy);

  // Returns true iff the displacements are equal
  bool equalDisplacement(DISP_INFO other);
};


/* Simple class used to represent a set of three integers. */
class TRIPLET {
public:
  int vals[3];
  int x, y, z;
  
  /* Construct a TRIPLET with the three provided components.*/
  TRIPLET(int myX, int myY, int myZ);

  /* Access one of the TRIPLETs three values.*/
  int& operator[](int index);
  
  /** Add each component of the triplet to that of the provided TRIPLET
   *  and return the result. */
  TRIPLET add(TRIPLET other); 
};


// The number and directions of a grid point's neighbors
const int NUM_DIRECTIONS = 6; 
const TRIPLET DIRECTIONS [6] = {TRIPLET(1,0,0), TRIPLET(-1,0,0), TRIPLET(0,1,0), TRIPLET(0,-1,0), TRIPLET(0,0,1), TRIPLET(0,0,-1)};







/* Return the 1-d index for the provided grid indices*/
int transformIndex(int x, int y, int z, int numX, int numY);

/* Returns the energy of the point located at the provided grid indices. 
 * Assumes the grid indices are within the appropriate range.*/
float getEnergy(int x, int y, int z, int numX, int numY, float *energyGrid);

/* Returns the energy of the point located in the grid referred to by the TRIPLET of indices. 
 * Assumes the grid indices are within the appropriate range.*/
float getEnergy(TRIPLET indices, int numX, int numY, float *energyGrid);

/* Returns true iff the grid point referred to by the provided grid indices
*  is indeed accessible. Assumes the grid indices are within the appropriate range. */
bool isAccessible(int x, int y, int z, int numX, int numY, float *accessGrid);

/* Returns the integer nearest to the provided double.*/
int nearestInt(double num);

/* Translate the coordinate by unit cell increments so that it lies within the 0 to 1 range.*/
double translate_to_original_uc(double x);

/* Adjusts the indices so that they lie within the appropriate range
*  and such that they are adjusted according to periodic boundary conditions.
*  Stores a TRIPLET representing to which unit cell the indices originally referred.*/
void adjustIndices(TRIPLET &gridIndices, TRIPLET &shift, UNIT_CELL &uc);

/* Returns true iff the grid point referred to by the provided triplet of grid indices
*  is indeed accessible. Assumes the grid indices are within the appropriate range. */
bool isAccessible(TRIPLET indices, int numX, int numY, float *accessGrid);

/* Returns true if the first pair has energy >= that of the second pair.*/
bool hasHigherEnergy(std::pair<TRIPLET, DISP_INFO > p1, std::pair<TRIPLET, DISP_INFO > p2);

// Helper function used in calculateMinEnergyBarrier()
bool findMinEnergyBarrier(int startX, int startY, int startZ, double &barrier, 
			  int numX, int numY, int numZ, float *energyGrid, float *accessGrid, UNIT_CELL unit_cell);

/* Calculates the minimum energy barrier of a path that must pass through the point referred to by the provided
 * grid indices. The resulting barrier is stored using the provided pointer. If no path is possible,
 * the barrier is set to DBL_MAX. The function returns true iff a path is found. */
bool calculateMinEnergyBarrier(double &minEnergy, double &barrierEnergy, float box_x, float box_y, float box_z, int numX, int numY, int numZ, float *energyGrid, float *accessGrid,
			       double vax, double vbx, double vby, double vcx, double vcy, double vcz);

#endif
