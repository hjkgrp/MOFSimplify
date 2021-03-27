#ifndef PSD_H
#define PSD_H

#include <iostream>
#include <vector>
#include <utility>
#include <string>

#include "networkstorage.h"
#include "material.h"

void calcPoreSizeDistr(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgAtomnet, bool highAccuracy, double r_probe_chan, double r_probe, int numSamples, bool excludePockets, std::string histogramFile, std::string pointsFile, std::string nodeAndRadiiFile, std::string spheresDistFile, bool visualize, bool overlapsCheck);

void NEWcalcPoreSizeDistr(MATERIAL *Mat, std::ostream &output); //ATOM_NETWORK *atmnet, ATOM_NETWORK *orgAtomnet, bool highAccuracy, double r_probe_chan, double r_probe, int numSamples, bool excludePockets, std::string histogramFile, std::string pointsFile, std::string nodeAndRadiiFile, std::string spheresDistFile, bool visualize, bool overlapsCheck);

void Histogram(std::ostream& output,const double binSize,const int maxBins, std::vector<double>& radii, int count, double nodefrac, double outfrac, int numSamples);

void Histogram_for_spheres(std::ostream& output, const double binSize, std::vector<double>& nodeDist, std::vector<double> &ghostDist);

void calcSphereIntersect(std::vector<std::pair <Point, double> > &inputNodeSphereAndRadii, std::vector<std::pair <Point, double> > &outputNodeSphereAndRadii);

void calcSpheresDistance(ATOM_NETWORK atomnetwork, std::vector <std::pair <Point, double> > &sphereAndRadii, std::vector <double> &distanceBetweenSpheres);

void printFileCoords_Radii(std::ostream& outfile, std::vector<std::pair <Point, double> > &pt_node_rad, std::vector<std::pair <Point, double> > &pt_ghost_rad, std::vector <double> &distanceBetweenNodeSpheres, std::vector <double> &distanceBetweenGhostSpheres, bool overlaps, double startRad, double endRad);

void checkDuplicates(ATOM_NETWORK *atomnetwork, VORONOI_NETWORK *vornet, std::vector <int> &listNodeIDs, std::vector <std::pair <Point, double> > &nodeAndRadii, Point nodePt, double newRadius);

#endif
