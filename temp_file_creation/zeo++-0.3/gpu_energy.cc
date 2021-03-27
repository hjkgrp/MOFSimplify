#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <vector>
#include "heap.h"
#include "gpu_energy.h"
//#include "network.h"

using namespace std;

/* Construct a unit cell with the provided vector components*/
UNIT_CELL::UNIT_CELL(double vax, double vbx, double vby, double vcx, double vcy, double vcz){
  va_x = vax; 
  vb_x = vbx; vb_y = vby; 
  vc_x = vcx; vc_y = vcy; vc_z = vcz;
  
  double invDet = 1/(va_x*vb_y*vc_z);
  inv_va_x = 1/va_x;
  inv_vb_x = invDet*-1*vc_z*vb_x; inv_vb_y = 1/vb_y;
  inv_vc_x = invDet*   (vc_y*vb_x-vb_y*vc_x);  inv_vc_y = invDet*-1*vc_y*va_x;  inv_vc_z = 1/vc_z;
}

/** Transform the coordinates relative to the uc vectors into ones relative to
 *  xyz coordinate system and store them using the provided pointers. */
void UNIT_CELL::abc_to_xyz(double a, double b, double c, double &x, double &y, double &z){
  x = a*va_x + b*vb_x + c*vc_x;
  y = b*vb_y + c*vc_y;
  z = c*vc_z;
}

/** Transform the coordinates relative to the xyz coordinate system vectors into ones relative to
 *  the unit cell vectors and store them using the provided pointers. */
void UNIT_CELL::xyz_to_abc(double x, double y, double z, double &a, double &b, double &c){
  a = x*inv_va_x + y*inv_vb_x + z*inv_vc_x;
  b = y*inv_vb_y + z*inv_vc_y;
  c = z*inv_vc_z;
}





// Constructor that creates a dummy container
DISP_INFO::DISP_INFO(){
  isReal = false;
} 

/** Create a container with the provided displacement components and energy.
*   Displacements are stored as characters to conserve space. */
DISP_INFO::DISP_INFO(int myA, int myB, int myC, float myMaxEnergy){
  isReal = true;
  a = (char)myA; 
  b = (char)myB; 
  c = (char)myC;
  maxEnergy = myMaxEnergy;
}

// Returns true iff the displacements are equal
bool DISP_INFO::equalDisplacement(DISP_INFO other){
  return (a == other.a) && (b == other.b) && (c == other.c);
}




/* Construct a TRIPLET with the three provided components.*/
TRIPLET::TRIPLET(int myX, int myY, int myZ){
  vals[0] = x = myX; 
  vals[1] = y = myY; 
  vals[2] = z = myZ;
}

/* Access one of the TRIPLETs three values.*/
int& TRIPLET::operator[](int index){
  if(index < 0 || index > 2){
    cerr << "Error: Invalid index to [] operator for TRIPLET instance" << "\n"
	 << "Exiting..." << "\n";
    exit(1);
  }
  return vals[index];
}

/** Add each component of the triplet to that of the provided TRIPLET
 *  and return the result. */
TRIPLET TRIPLET::add(TRIPLET other){
  return TRIPLET(x + other.x, y + other.y, z + other.z);
}





/* Return the 1-d index for the provided grid indices*/
int transformIndex(int x, int y, int z, int numX, int numY){
  return x + y*numX + z*numX*numY;
}

/* Returns the energy of the point located at the provided grid indices. 
 * Assumes the grid indices are within the appropriate range.*/
float getEnergy(int x, int y, int z, int numX, int numY, float *energyGrid){
  return energyGrid[transformIndex(x, y, z, numX, numY)];
}

/* Returns the energy of the point located in the grid referred to by the TRIPLET of indices. 
 * Assumes the grid indices are within the appropriate range.*/
float getEnergy(TRIPLET indices, int numX, int numY, float *energyGrid){
  return getEnergy(indices[0], indices[1], indices[2], numX, numY, energyGrid);
}

/* Returns true iff the grid point referred to by the provided grid indices
*  is indeed accessible. Assumes the grid indices are within the appropriate range. */
bool isAccessible(int x, int y, int z, int numX, int numY, float *accessGrid){
  cerr << "Error: Accessibility check not implemented" << "\n";
  exit(1);
  return true;
}

/* Returns the integer nearest to the provided double.*/
int nearestInt(double num){
  return (int)(floor(num + 0.5));
}

/* Translate the coordinate by unit cell increments so that it lies within the 0 to 1 range.*/
double translate_to_original_uc(double x){
  double newX =  x - (x < 0.0 ?(-(int)(0.5-x)):((int)(0.5+x)));
  if(newX < 0)
    return newX+1;
  else
    return newX;
} 

/* Adjusts the indices so that they lie within the appropriate range
*  and such that they are adjusted according to periodic boundary conditions.
*  Stores a TRIPLET representing to which unit cell the indices originally referred.*/
void adjustIndices(TRIPLET &gridIndices, TRIPLET &shift, UNIT_CELL &uc){
  double a, b, c;
  uc.xyz_to_abc(gridIndices[0], gridIndices[1], gridIndices[2], a, b, c);
  double newA, newB, newC;
  newA = translate_to_original_uc(a); newB = translate_to_original_uc(b); newC = translate_to_original_uc(c); 
  double x, y, z;
  uc.abc_to_xyz(newA, newB, newC, x, y, z);
  shift = TRIPLET(nearestInt(a - newA), nearestInt(b - newB), nearestInt(c - newC));
  gridIndices = TRIPLET(nearestInt(x), nearestInt(y), nearestInt(z));
}

/* Returns true iff the grid point referred to by the provided triplet of grid indices
*  is indeed accessible. Assumes the grid indices are within the appropriate range. */
bool isAccessible(TRIPLET indices, int numX, int numY, float *accessGrid){
  return isAccessible(indices.x, indices.y, indices.z, numX, numY, accessGrid);
}

/* Returns true if the first pair has energy >= that of the second pair.*/
bool hasHigherEnergy(pair<TRIPLET, DISP_INFO > p1, pair<TRIPLET, DISP_INFO > p2){
  return p1.second.maxEnergy >= p2.second.maxEnergy;
}

// Helper function used in calculateMinEnergyBarrier()
bool findMinEnergyBarrier(int startX, int startY, int startZ, double &barrier, 
			  int numX, int numY, int numZ, float *energyGrid, float *accessGrid, UNIT_CELL unit_cell){
  barrier = DBL_MAX;
  vector<DISP_INFO> visitedNodes = vector<DISP_INFO> (numX*numY*numZ, DISP_INFO());
  
  // Place starting node on heap with (0,0,0) displacement
  HEAP<pair<TRIPLET, DISP_INFO > > heap (hasHigherEnergy); // Lowest energy paths are favored
  heap.insert(pair<TRIPLET, DISP_INFO> (TRIPLET(startX, startY, startZ), DISP_INFO(0, 0, 0, getEnergy(startX, startY, startZ, numX, numY, energyGrid))));

  while(heap.size() != 0){
    // Remove best remaining path
    pair<TRIPLET, DISP_INFO> pathInfo = heap.pop();
    int arrayIndex = transformIndex(pathInfo.first[0], pathInfo.first[1], pathInfo.first[2], numX, numY);
    DISP_INFO nodeInfo = visitedNodes[arrayIndex];

    if(nodeInfo.isReal){
      if(nodeInfo.equalDisplacement(pathInfo.second)){
	// Circling back to previous node so terminate path
	continue;
      }
      else {
	// Found the optimal path
	barrier = min(pathInfo.second.maxEnergy, nodeInfo.maxEnergy);
	return true;
      }
    }
    else{
      // Add all connected neigbors if appropriate
      visitedNodes[arrayIndex] = pathInfo.second;
      
      for(int i = 0; i < NUM_DIRECTIONS; i++){
	TRIPLET dir = DIRECTIONS[i];
	TRIPLET newIndices = dir.add(pathInfo.first);
	TRIPLET change_in_uc(0,0,0);
	adjustIndices(newIndices, change_in_uc, unit_cell);	

	if(isAccessible(newIndices, numX, numY, accessGrid)){
	  int da = change_in_uc[0] + pathInfo.second.a;
	  int db = change_in_uc[1] + pathInfo.second.b;
	  int dc = change_in_uc[2] + pathInfo.second.c;
	  heap.insert(pair<TRIPLET, DISP_INFO> (newIndices, DISP_INFO(da, db, dc, max(pathInfo.second.maxEnergy, getEnergy(newIndices, numX, numY, energyGrid)))));
	}
      }
    }
  } 
  return false;
}

/* Calculates the minimum energy barrier of a path that must pass through the point referred to by the provided
 * grid indices. The resulting barrier is stored using the provided pointer. If no path is possible,
 * the barrier is set to DBL_MAX. The function returns true iff a path is found. */
bool calculateMinEnergyBarrier(double &minEnergy, double &barrierEnergy, float box_x, float box_y, float box_z, int numX, int numY, int numZ, float *energyGrid, float *accessGrid,
			       double vax, double vbx, double vby, double vcx, double vcy, double vcz){
  cout << "Identifying global accessible energy minimum" << "\n";
  // Find global energy minimum and its corresponding grid indices
  minEnergy = DBL_MAX;
  int minZ = -1, minY = -1, minX = -1;
  for(int i = 0; i < numZ; i++){
    for(int j = 0; j < numY; j++){
      for(int k = 0; k < numX; k++){
	if(isAccessible(k, j, i, numX, numY, accessGrid)){
	  float energy = getEnergy(k, j, i, numX, numY, energyGrid);
	  if(energy < minEnergy){
	    minEnergy = energy;
	    minZ = i; minY = j; minX = k;
	  }
	}
      }
    }
  }

  // Ensure that a minimum was found
  if((minX == -1) || (minY == -1) || (minZ == -1)){
    cerr << "Error: Unable to identify global energy minimum when calculating energy barriers." << "\n"
	 << "Aborting energy calculation" << "\n";
    return false;
  }
  cout << "Global accessible energy minimum found" << "\n";
  
  // Construct the unit cell vectors in terms of grid points
  double grid_spacing_x = box_x/numX, grid_spacing_y = box_y/numY, grid_spacing_z = box_z/numZ;
  vax /= grid_spacing_x;
  vbx /= grid_spacing_x; vby /= grid_spacing_y;
  vcx /= grid_spacing_x; vcy /= grid_spacing_y; vcz /= grid_spacing_z;
  UNIT_CELL unit_cell(vax, vbx, vby, vcx, vcy, vcz);

  cout << "Calculating energy barrier for box of dimensions:" << box_x << " by " << box_y << " by " << box_z << "\n"
       << "Grid size: " << numX << " by " << numY << " by " << numZ << "\n"
       << "Grid unitcell vectors:" 
       << "\t v_a:(" << vax << ", 0, 0)" << "\n"
       << "\t v_b:(" << vbx << ", " << vby << ", 0)" << "\n"
       << "\t v_c:{" << vcx << ", " << vcy << ", " << vcz << ")" << "\n" 
       << "\n" << "\n";

  // Calculate the barrier
  bool result = findMinEnergyBarrier(minX, minY, minZ, barrierEnergy, numX, numY, numZ, energyGrid, accessGrid, unit_cell);
  cout << "Energy barrier calculation completed" << "\n" << "\n";
  return result;
}
