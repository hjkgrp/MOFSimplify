#ifndef AREA_AND_VOLUME_H
#define AREA_AND_VOLUME_H

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <voro++.hh>

#include "networkstorage.h"
#include "geometry.h"
#include "voronoicell.h"
#include "material.h"

//#include <voro++>

/** Returns the density of the provided ATOM_NETWORK, assuming it represents a single
 *  unit cell, in g/cm^3. */
double calcDensity(ATOM_NETWORK *atmnet);

/* Print the coordinates contained in the provided vectors in a manner that they can be
 * displayed using ZeoVis and its VMD interface. Accessible points are colored green while
 * inaccessible points are colored red.*/
void reportPoints(std::ostream &output, std::vector<Point> axsPoints, std::vector<Point> inaxsPoints);
void reportPointsVisIT(std::ostream &output, std::vector<Point> axsPoints, std::vector<Point> inaxsPoints);

/* Print the coordinates contained in the provided vectors in a manner that they can be
 * displayed using ZeoVis and its VMD interface. Useful for debugging errors in the
 * ASA and AV calculations. */
void reportResampledPoints(std::ostream &output, std::vector< std::pair<int, Point> > resampledInfo);

/* Adjust the sampling point to minimize its distance with the central atom. */
void adjustSamplingPoint(Point *samplingPointXYZ, Point newAtomCoordsXYZ, Point origAtomCoordsXYZ, ATOM_NETWORK *atmnet);

// Function defined in network.cc
void* performVoronoiDecomp(bool, ATOM_NETWORK *, VORONOI_NETWORK *, std::vector<VOR_CELL> &, bool,  std::vector<BASIC_VCELL> &);

/* NEW functions that operate on MATERIAL class */
void NEWreportPoints(ostream &output, ATOM_NETWORK *atmnet, vector<Point> *axsPoints, vector<int> *axsPChIDs,
                                                         vector<Point> *inaxsPoints, vector<int> *inaxsPPIDs,
                                                         string type);

void NEWreportPointsValue(ostream &output, ATOM_NETWORK *atmnet, vector<Point> *, vector<int> *, vector<double>*, string);

double NEWcalcAV(MATERIAL *Mat, double r_probe, int numSamples, double low_dist_cutoff, double high_dist_cutoff);
double NEWcalcAV(MATERIAL *Mat, double r_probe, int numSamples);
void NEWcalcAVprint(MATERIAL *Mat, ostream &output, char *filename);

double NEWcalcASA(MATERIAL *Mat, double r_probe, int numSamples);
void NEWcalcASAprint(MATERIAL *Mat, ostream &output, char *filename);

/** Returns the volume accessible to a particle of the provided radius. Accessible volume is defined as any region of space
 *  in which the center of the particle of the provided radius can reach. The first radius is used to determine Voronoi node accessibility,
 *  whereas the second is used to check for overlap. Excludes inaccessible pockets if requested.*/
double calcAV(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgatmnet, bool highAccuracy, double r_probe_chan, double r_probe, int numSamples, bool excludePockets, std::ostream &output, char *filename, bool visualize, bool VisITflag, bool LiverpoolFlag, bool blockingMode, double low_dist_cutoff, double high_dist_cutoff, bool ProbeOccupiableFlag);
double calcAV(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgatmnet, bool highAccuracy, double r_probe, int numSamples, bool excludePockets, std::ostream &output, char *filename, bool visualize, bool VisITflag, bool LiverpoolFlag, bool blockingMode, double low_dist_cutoff, double high_dist_cutoff);
// For python interface
std::string calcAV(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgatmnet, bool highAccuracy, double r_probe_chan, double r_probe, int numSamples, bool excludePockets, double low_dist_cutoff, double high_dist_cutoff);

double calcASA(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgatmnet, bool highAccuracy, double r_probe_chan, double r_probe, double rho_crystal, int numSamples, bool excludePockets, std::ostream &output, char *filename, bool visualize, bool VisITflag, bool LiverpoolFlag, bool ExtendedOutputFlag);
double calcASA(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgatmnet, bool highAccuracy, double r_probe, double rho_crystal, int numSamples, bool excludePockets, std::ostream &output, char *filename, bool visualize, bool VisITflag, bool LiverpoolFlag, bool ExtendedOutputFlag);
// For python interface
std::string calcASA(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgatmnet, bool highAccuracy, double r_probe_chan, double r_probe, int numSamples, bool excludePockets,  bool ExtendedOutputFlag);

void determineAccessibility(ATOM_NETWORK *atmnet, double r_probe_chan, double r_probe, bool excludePockets, std::vector<bool> *isAccessible, VORONOI_NETWORK *pointSet);
void visVoro(char* name, double probeRad, int skel_a, int skel_b, int skel_c, VORONOI_NETWORK* vornet, ATOM_NETWORK* atmnet);

//extract spherical substructures: this is a functionality designed for extracting local substructures of zeolites so that they can be scanned for potential guest molecule binding sites
//this functionality writes out a number of xyz format files containing spherical substructures of the given radius, centred on given probe-accessible Voronoi nodes; if an element_type is given, a simplified Voronoi network is used, based only on atoms of that type
void getLocalSubstructures(char* name, double probeRad, double local_substructure_radius, TRIPLET supercell_steps_local_substructure, std::string element, VORONOI_NETWORK* vornet, ATOM_NETWORK* atmnet, bool radial, bool usingElementOnly);

// Function saving blocking spheres given accessible and nonaccesible points
void blockPockets(ATOM_NETWORK *atmnet, std::ostream &output, std::vector<Point> axsPoints, std::vector<int> axsPChIDs, std::vector<Point> inaxsPoints, std::vector<int> inaxsPPIDs, double probeRad);
// Partner function to the above - identifies the most dense point in a point cloud, to identify where a blocking sphere should be centred
int get_most_dense_index(ATOM_NETWORK *atmnet, std::vector<Point> *points_vector);

/** Returns if a given point is accessible. Accessible points are defined as points that lie inside of the void network that is accessible to a probe of radius r_probe.*/
/* to be removed (backup only)
bool accessiblePoint(Point sample_pnt,voro::container_periodic_poly *new_rad_con, double r_probe, ATOM_NETWORK *atmnet, std::vector<BASIC_VCELL> vorcells,std::vector<bool> accessInfo);
*/


/* NEW ADDITION TEMPORARY IN THIS FILE */
double calcMolecule(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgatmnet, bool highAccuracy, double r_probe_chan, double r_probe, int numSamples, bool excludePockets, std::ostream &output, char *filename, bool visualize, bool VisITflag, bool LiverpoolFlag, bool blockingMode, double low_dist_cutoff, double high_dist_cutoff);



#endif
