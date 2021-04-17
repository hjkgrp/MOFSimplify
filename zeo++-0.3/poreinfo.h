/** 
  * In poreinfo.* files, there are functions to perform time-series analysis for
  * void space data contain in a set of .poreinfo summary files
  *
  * Initial code by M. Haranczyk, June 2013
  *
  **/

#ifndef POREINFO_H
#define POREINFO_H

#include <cstdio>
#include <string>
#include "general.h"
#include "geometry.h"
#include "networkstorage.h"


/* POREINFO class is used to store and analyze pore information
   from time-series data */

class POREINFO {

  public:
  int dim; // dimentionality of pore
  
  double di; // largest inclluded sphere
  double sa, vol; // surface area and volume 

  Point pos; // pore position, coordinates of Di
  double enc_radius; // radius of encapsulating sphere

  std::vector <NODESPHERE> nodes;

  std::vector <int> prevIDs; // IDs of corresponding pores in previous and next time steps
  std::vector <int> nextIDs;

  POREINFO()
   {
   dim = 0;
   nodes.clear();
   prevIDs.clear();
   nextIDs.clear();
   };

 };


/* Function that performs analysis for a set of .poreinfo files listed in listfilename.
 * Output is written to outputfilename. Atom_network has to be provided to enable distance 
   calculations */

void  analyzePoreInfoFiles(ATOM_NETWORK *atmnet, std::string listfilename, std::string outputfilename);

/* Loads .poreinfo files into POREINFO vector */
void loadPoreInfoFile(std::vector < std::vector<POREINFO> > *, std::string);

/* Function checks if two pores ar connected because of overlap of nodes */
bool arePoresConnected(POREINFO *, POREINFO *);

#endif
