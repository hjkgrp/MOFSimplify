/* This file is used to store information about the radii and weights of
 * atoms. The information can be specified in input files or default values can be used
 * for more common types.
 */

#ifndef NETWORKINFO_H
#define NETWORKINFO_H

#include <map>
#include <string>

/* Moved to networkinfo.cc
static std::map <std::string,double> radTable;
static std::map <std::string,double> covRadTable;
static std::map <std::string,double> massTable;
static std::map <std::string,int> atomicNumberTable;
static std::map <std::string,bool> atomicCharacterTable;
*/

/** Initialize atom name stripping */
void initializeStripAtomNameInternalFlag(bool);

/** Fills the radius table with several default values. */
void initializeRadTable();

/** Fills the covalent radius table with several default values. */
void initializeCovRadTable();

/** Fills the mass table with several default values
 ** in units of g/mole. */
void initializeMassTable();

/** Fills metal/nonmetal info for atoms */
void initializeAtomCharacterTable();

/** Fills the atomic number table with atomic number of all elements
 ** */
void initializeAtomicNumberTable();

/** Fills the periodic table table to be used to atom name recognition
 ** */
void initializePT();

/** Reads the radius table from the provided filename. The file must be
    formatted in two columns: atom name and radius. */
void readRadTable(char *filename);

/** Reads the mass table from the provided filename. The file must be
    formatted in two columns: atom name and mass in g/mole. */
void readMassTable(char *filename);

/** Return the radius for the corresponding atom name. If the -nor
    option was specified, returns 0. */
double lookupRadius(std::string atomType, bool radial);

/** Return the covalent radius for the corresponding atom name */
double lookupCovRadius(std::string atomType);

/** Return the mass for the corresponding atom name. */
double lookupMass(std::string atomType); 

/** Return the atomic number for the corresponding symbol */
int lookupAtomicNumber(std::string atomType);

/** Return true is an atom is a metal **/
bool isMetal(std::string atomType);

/** Strip atom name to remove any index strings */
std::string stripAtomName(std::string extAtom);

#endif
