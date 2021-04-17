// 
//
// Author   : Thomas F. Willems (LBL / UC Berkeley)
// Email    : tfwillems@berkeley.edu
// Date     : July 14 2010
// Updated  : April 29, 2011

/* Note: Code contained in performVoronoiDecomp() method is largely derived from
 * code originally written by Chris Rycroft, the author of the voro++
 * software package used to perform the Voronoi decomposition.
 */


#ifndef NETWORK_H
#define NETWORK_H


#include <map>
#include <set>
#include <cstring>
#include <string>
#include <vector>



//Options for reading .cif files
#define CIF_RMV_NUMS_FROM_ATOM_TYPE false

//Options for reading .car files
#define CAR_USE_ATOM_TYPE_OVER_NAME true // take atom type from 8th column rather than atom name from 1st


#include "heap.h"
#include <voro++.hh>

#include "networkstorage.h"
#include "voronoicell.h"
#include "graphstorage.h"
#include "string_additions.h"
#include "material.h"

template<class c_option>
bool storeVoronoiNetwork(c_option &con, ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet, double bx, double by, double bz,
			 std::vector<BASIC_VCELL> &basCells, std::vector<int> &atomShifts, bool storeAdvCells, std::vector<VOR_CELL> &advCells);

// A guess for the memory allocation per region
const int memory=16;

// A maximum allowed number of regions, to prevent enormous amounts of memory
// being allocated
const int max_regions=16777216;

// A buffer size
const int bsize=2048;




/* converts string to upper case */
std::string toUpperCase(const std::string & s);



/** Decompose the provided network of atoms into a VORONOI_NETWORK that is stored using the provided pointer. If
    the option is specified, information about each VOR_CELL will also be stored using the provied pointer to
    a vector of VOR_CELL instances. The BASIC_VCELL information is stored regardless of the option specified.*/
void* performVoronoiDecomp(bool radial, ATOM_NETWORK *cell, VORONOI_NETWORK *vornet, std::vector<VOR_CELL> &cells, bool saveVorCells,
			   std::vector<BASIC_VCELL> &bvcells);

/** Decompose the provided network of atoms into a VORONOI_NETWORK that is stored using the provided pointer. If
    the option is specified, information about each VOR_CELL will also be stored using the provied pointer to
    a vector of VOR_CELL instances. The BASIC_VCELL information is stored regardless of the option specified. This function is wrapper to above function, and has no return type*/
bool performVoronoiDecomp(bool radial, ATOM_NETWORK *cell, VORONOI_NETWORK *vornet, std::vector<VOR_CELL> *cells, bool saveVorCells,
			   std::vector<BASIC_VCELL> *bvcells);


void createAdvCell(voro::voronoicell &cell, std::vector<double> coords, int *idMap, VOR_CELL &newCell);

/** Extends the provided unit cell in the x, y and z directions using
    the given factors and stores the resulting ATOM_NETWORK using the
    pointer to newCell.For instance, if xfactor = yfactor = zfactor= 1, newCell contains the
    same information as cell. In contrast, if xfactor = yfactor =
    zfactor = 2, newCell is a 2x2x2 replication of cell in 3D space. */
void extendUnitCell(ATOM_NETWORK *cell, ATOM_NETWORK *newCell, int xfactor, int yfactor, int zfactor);


/** Extend the given VORONOI_NETWORK by 10 unit cells along a unit cell vector in the provided
 *  direction. The direction must be either (1,0,0), (0,1,0) or (0,0,1).  The list of atoms 
 *  belonging to each atom are NOT APPROPRIATELY REPLICATED.*/
void extendVorNet(VORONOI_NETWORK *vornet, VORONOI_NETWORK *newNet, DELTA_POS direction, std::map<int,int> *idAliases, std::set<int> *sourceNodes);

void calculateFreeSphereParameters(VORONOI_NETWORK *vornet, char *filename, bool extendedPrintout);
void NEWcalculateFreeSphereParameters(MATERIAL *Mat);
void NEWcalculateFreeSphereParametersPrint(MATERIAL *Mat, char *filename, bool);

void viewVoronoiDecomp(ATOM_NETWORK *atmnet, double r_probe, std::string filename);

void loadRadii(ATOM_NETWORK *atmnet);

void loadMass(bool useMassFlag, ATOM_NETWORK *atmnet);


/* Print information about topology of the structure */
void getStructureInformation(char *filename, char *filenameExtendedOutput, ATOM_NETWORK *atmnet, bool extendedOutput);

/* Print information about presence of open metal sites  */
void getOMSInformation(char *filename, char *filenameExtendedOutput, ATOM_NETWORK *atmnet, bool extendedOutput);

#endif
