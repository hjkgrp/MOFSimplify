#ifndef NETWORKIO_H
#define NETWORKIO_H

//#include "network.h"
#include <iostream>

#include "networkstorage.h"

#define SUPERCELL_SIZE 2

/** Identifies the extension and the prefix present in the provided
    filename and stores them using the two provided character
    pointers. */
void parseFilename(const char * fileName, char *name, char *extension);

/** Ensures that the provided file is of type .cuc, .cssr or .v1 or .cif. 
 *  Otherwise, an error message is displayed and the program is
 *  aborted. */
bool checkInputFile(char * filename);

/** Read the information from the .cif file referred to by filename
    and store it within the provided ATOM_NETWORK **/
//void readCIFFile(char *filename, ATOM_NETWORK *cell, bool radial);
bool readCIFFile(char *filename, ATOM_NETWORK *cell, bool radial);

/** Read the information from the .car file referred to by filename
    and store it within the provided ATOM_NETWORK **/
bool readCARFile(char *filename, ATOM_NETWORK *cell, bool radial);

/** Read the information from the .pdb file referred to by filename
    and store it within the provided ATOM_NETWORK **/
bool readPDBFile(char *filename, ATOM_NETWORK *cell, bool radial);

/** Read the .arc file refererred to by filename and store its
    information within the provided ATOM_NETWORK. */
//void readARCFile(char *filename, ATOM_NETWORK *cell, bool radial);
bool readARCFile(char *filename, ATOM_NETWORK *cell, bool radial);

/** Read the .cuc file refererred to by filename and store its
    information within the provided ATOM_NETWORK. */
//void readCUCFile(char *filename, ATOM_NETWORK *cell, bool radial);
bool readCUCFile(char *filename, ATOM_NETWORK *cell, bool radial);

/** Read the information from the .cssr file referred to by filename
    and store it within the provided ATOM_NETWORK. */
//void readCSSRFile(char *filename, ATOM_NETWORK *cell, bool radial);
bool readCSSRFile(char *filename, ATOM_NETWORK *cell, bool radial);

/** Read the information from the .obcssr file referred to by filename
 *     and store it within the provided ATOM_NETWORK. 
 *     obcssr is Open Babel generated CSSR file */
//void readCSSRFile(char *filename, ATOM_NETWORK *cell, bool radial);
bool readOBCSSRFile(char *filename, ATOM_NETWORK *cell, bool radial);

/** Read the information from the .v1 file referrred to by filename
    and store it within the provided ATOM_NETWORK. */
//void readV1File(char *filename, ATOM_NETWORK *cell, bool radial);
bool readV1File(char *filename, ATOM_NETWORK *cell, bool radial);

/** Read the information from the .dlp file referrred to by filename
    and store it within the provided ATOM_NETWORK. */
bool readDLPFile(char *filename, ATOM_NETWORK *cell, bool radial);


/** Read the VORONOI_NETWORK information located in the provided input stream and 
 *  store it using the provided VORONOI_NETWORK pointer. The input stream must have a file format
 *  corresponding to a .net or .nt2 format. */
void readNet(std::istream *input, VORONOI_NETWORK *vornet);

/* Read the VORONOI_NETWORK located in the provided file and store 
 * it using the provided network pointer. The file must be in the .net/.nt2 
 * file format. */ 
//void readNetFile(char * filename, VORONOI_NETWORK *vornet);
bool readNetFile(char * filename, VORONOI_NETWORK *vornet);

/** Write the information within the provided ATOM_NETWORK in a .cssr
    file format to the provided filename. */
bool writeToCSSR(char *filename, ATOM_NETWORK *cell);

/** Write the information within the provided ATOM_NETWORK in a .cssr
    file format to the provided filename, using labels instead of elements types. */
bool writeToCSSRLabeled(char *filename, ATOM_NETWORK *cell);

/** Write the infomation within the provided ATOM_NETWORK in a .cif
    file format to the provided filename. **/
bool writeToCIF(char *filename, ATOM_NETWORK *cell);

/** Write the information within the provided ATOM_NETWORK in a .v1
    file format to the provided filename. */
bool writeToV1(char * filename, ATOM_NETWORK *cell);

/** Write the information stored within the VORONOI_NETWORK in a .nt2
    file format to the provided filename. Excludes any nodes or nodes with radii
    less than the provided threshold. For the default 0 minRad value, all nodes
    and edges are included. */
bool writeToNt2(char *filename, VORONOI_NETWORK *vornet, double minRad = 0);

/** Write the voronoi node information stored within the VORONOI_NETWORK in XYZ
    file format to the provided filename. Excludes any nodes with radii
    less than the provided threshold. For the default 0 minRad value, all nodes
    are included. */
bool writeToXYZ(char *filename, VORONOI_NETWORK *vornet, double minRad = 0);

/** Write the information stored within the VORONOI_NETWORK in a .nt2
    file format to the provided filename. Includes all nodes and edges. */
// Redundant
//bool writeToNt2(char *filename, VORONOI_NETWORK *vornet);

/** Write the information within the provided ATOM_NETWORK in a .xyz
    file format to the provided filename. */
bool writeToXYZ(char *filename, ATOM_NETWORK *cell, bool is_supercell, bool is_duplicate_perimeter_atoms);

/** Write the information within the provided ATOM_NETWORK in a .vtk
    file format to the provided filename. */
bool writeToVTK(char *filename, ATOM_NETWORK *cell);

/** Write the information within the provided ATOM_NETWORK in a .mop
    file format to the provided filename. */
bool writeToMOPAC(char *filename, ATOM_NETWORK *cell, bool is_supercell);


/** Change the type of the atom to its appropriate form. Used when
 reading .cuc files. */
void changeAtomType(ATOM *atom);

/** Strip atom names from additional strings */
void stripAtomNames(ATOM_NETWORK *cell);

#endif
