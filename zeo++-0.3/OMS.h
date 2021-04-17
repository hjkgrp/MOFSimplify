// OMS, a library for determining the exposure grade of a molecule
//
// OMS.hh
//
// Routines for single molecule analysis
//
// Author   : Ismael Gómez García
// Email    : 
// Date     : November 27th 2015

#ifndef _OMS_H_
#define _OMS_H_

#include <vector>
#include <cmath>
#include <algorithm>

#include <stdlib.h>

#include <iostream>

#include "Eigen/Dense"

#include "zeo_consts.h"

#define DIM 3

#define PHI_COORD 	0
#define THETA_COORD 	1
#define R_COORD		2

#define X_COORD 	0
#define Y_COORD 	1
#define Z_COORD		2

#define ERR_TOL 0.0000000000000001


using namespace std;
using namespace Eigen;

bool 	IsExposedMolecule 		(vector < vector <double> > MoleculeCoordinates);
bool 	IsExposedMoleculeThreshold	(vector < vector <double> > MoleculeCoordinates, double Threshold );
double 	DegreeOfExposure		(vector < vector <double> > MoleculeCoordinates);



#endif
