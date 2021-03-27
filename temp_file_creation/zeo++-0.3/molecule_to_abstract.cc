//Richard L Martin 2013/05/16
//C code to read an xyz file and convert it to an abstract linker - it can handle carboxy or Br sites, with the appropriate flag
//Arguments: real_molecule.xyz site_type output_abstract_molecule.xyz

/*------------------------------
       INCLUDES, CONSTS
------------------------------*/

#include <algorithm>
#include <stdlib.h>
#include <vector>
#include <string>
#include <cstring>
#include <cstdio>
using namespace std;
#include "geometry.h"
#include "mindist.h"
#include "networkinfo.h"
#include "rmsd.h"
#include "networkstorage.h"
#include "net.h"
#include "symmetry.h"
#include "string_additions.h"
#include "zeo_consts.h"

const int CARBOXY = 0, BROMINE = 1, DIALCOHOL = 2, A_P_CATECHOL = 3, ALDEHYDE = 4;

const int MAX_CHAR_ARRAY_LENGTH = 500;
const double CARBON_CARBON_BOND_LENGTH = 1.43; //Angstroms
const double MAX_HYDROGEN_BOND_LENGTH = 1.2; //Angstroms
//const double MAX_CARBOXY_BOND_LENGTH = 1.5; //Angstroms
const double MAX_CARBOXY_BOND_LENGTH = 1.6; //Angstroms
const double MAX_BROMINE_BOND_LENGTH = 2.1; //Angstroms
const double MAX_SITE_DISTANCE_DEVIATION = 2.0; //Angstroms
const double MAX_SITE_ANGLE_DEVIATION = 30.0; //degrees
const int NUM_MC = 100000;
const double DEFAULT_RADIUS = 1.35; //Angstroms
const double GRID_RES = 0.25; //Angstroms

/*------------------------------
            CLASSES
------------------------------*/

//int3
class INT3 {
public:
  int x, y, z;
  //constructor
  INT3();
};

/*------------------------------
      METHOD DECLARATIONS
------------------------------*/

/*
void search_for_char(FILE *f, char c);
*/
void read_molecule(string filename, vector<ATOM> *atom_vector);
void read_incon(string filename, vector<ATOM> *atom_vector);
void read_xyz(string filename, vector<ATOM> *atom_vector);
void write_xyz(FILE *outfile, vector<ATOM> *atom_vector, string input_name);
void find_sites(vector<ATOM> *molecule, vector<INT3> *indices, int site_type);
void replace_sites_with_abstract(vector<ATOM> *molecule, vector<INT3> *indices, int site_type, string name);
bool bonded(ATOM a, ATOM b, double d);
double dist(ATOM a, ATOM b);
vector<ATOM> rotate(vector<ATOM> *orig, double mat[3][3]);
double get_rad(string s);
double dist(double x, double y, double z, ATOM b);
XYZ get_vector_from_to(XYZ start, XYZ end);
XYZ midpoint(XYZ a, XYZ b);
XYZ project_onto_line(XYZ initial, XYZ line_start, XYZ line_end);
double radians_to_degrees(double r);
/*
vector<string> split(string line, string delims);
*/
//sorting functions
bool increasing_double (double i, double j) { return (i<j); }

/*------------------------------
             MAIN
------------------------------*/

int main(int argc, char *argv[]) {

//-----hardcoded settings
  bool verbose = false;

//-----check arguments
  int num_args = 4;
  if(argc!=num_args) {
    printf("%d arguments were provided but %d are required:\n", argc, num_args);
    printf("a.out\n");
    printf("real_molecule.xyz\n");
    printf("site_type\n");
    printf("abstract_molecule_output.xyz\n");
    printf("Where site_type is one of:\n\t0 (carboxy)\n\t1 (Br)\n\t2 (dialcohol)\n\t3 (acetenol-protected catechol)\n\t4 (aldehyde)\n");
    printf("Please try again.\n");
    exit(EXIT_FAILURE);
  }

  //by default, we want to populate the table with radii and masses, even if later we will not reference them
  initializeRadTable();
  initializeMassTable();
  initializeCovRadTable();

//-----parse arguments
  string real_mol_name, out_name;
  real_mol_name = argv[1];
  int site_type = atoi(argv[2]);
  out_name = argv[3];
  //read input, saving a vector of atoms found
  vector<ATOM> real_molecule;
  read_molecule(real_mol_name, &real_molecule);
  int num_real_atoms = real_molecule.size();
if(verbose) printf("real molecule contains %d atoms\n", num_real_atoms);
  //find the connection sites
  vector<INT3> real_site_indices;
  find_sites(&real_molecule, &real_site_indices, site_type);
  int num_real_sites = real_site_indices.size();
if(verbose) printf("real molecule contains %d sites\n", num_real_sites);
  //replace the site atoms with an abstract site
  replace_sites_with_abstract(&real_molecule, &real_site_indices, site_type, real_mol_name);

//-----write out the molecule to output
  FILE *outfile = fopen(out_name.c_str(), "w");
  if(outfile==NULL) {
    printf("ERROR: could not open output molecule file with name %s\n", out_name.c_str());
    exit(EXIT_FAILURE);
  }
  write_xyz(outfile, &real_molecule, real_mol_name);
  //printf("Program complete\n\n");
}

/*------------------------------
         OTHER METHODS
------------------------------*/

void replace_sites_with_abstract(vector<ATOM> *molecule, vector<INT3> *indices, int site_type, string name) {
  bool verbose = true;
  int num_sites = indices->size();
  int num_atoms = molecule->size();
  if(site_type==CARBOXY) { //replace the C atom with Q, and remove the O atoms; also remove any H atoms near to the O
    for(int i=0; i<num_sites; i++) {
      INT3 site = indices->at(i);
      molecule->at(site.y).type = "Q"; //replace the C
      molecule->at(site.y).label = "Q";
      molecule->at(site.x).keep = false;
      molecule->at(site.z).keep = false; //remove the Os
      //now find H atoms
      for(int j=0; j<num_atoms; j++) { //find H atoms bonded to this O
        if(molecule->at(j).type=="H" && molecule->at(j).keep) {
          if(bonded(molecule->at(site.x), molecule->at(j), MAX_HYDROGEN_BOND_LENGTH) || bonded(molecule->at(site.z), molecule->at(j), MAX_HYDROGEN_BOND_LENGTH)) {
            molecule->at(j).keep = false;
          }
        }
      }
    }
  } else if(site_type==BROMINE) { //else remove the Br atom and replace it with Q positioned so as to be a carbon bond length away from the adjacent atom
    if(num_sites!=2) {
/*
      printf("ERROR: >2 Br sites (%d) in structure %s so cannot align them - not writing out\n", num_sites, name.c_str());
      exit(EXIT_FAILURE);
*/
      printf("WARNING: >2 Br sites (%d) in structure %s so cannot align them - will write out abstract molecule as it is\n", num_sites, name.c_str());
      for(int i=0; i<num_sites; i++) {
        INT3 site = indices->at(i);
        molecule->at(site.x).type = "Q"; //replace the Br with Q
        molecule->at(site.x).label = "Q";
        molecule->at(site.x).set_xyz(molecule->at(site.y).xyz() + get_vector_from_to(molecule->at(site.y).xyz(), molecule->at(site.x).xyz()).unit().scale(MAX_CARBOXY_BOND_LENGTH)); //shift the Q
      }
    } else { //if dibromo, we can also orient them so that the structure has a chance to be built with nonlinear molecules
      INT3 site_0 = indices->at(0);
      INT3 site_1 = indices->at(1);
      molecule->at(site_0.x).type = "Q"; //replace the Br with Q
      molecule->at(site_0.x).label = "Q";
      molecule->at(site_1.x).type = "Q";
      molecule->at(site_1.x).label = "Q";
      //now, if both Br are anchored to the same atom, we have a problem - we will have to handle it separately:
      if(site_0.y==site_1.y) {
        XYZ site_0_1 = get_vector_from_to(molecule->at(site_0.x).xyz(), molecule->at(site_1.x).xyz()).unit().scale(MAX_CARBOXY_BOND_LENGTH);
        molecule->at(site_1.x).set_xyz(molecule->at(site_1.y).xyz() + site_0_1);
        molecule->at(site_0.x).set_xyz(molecule->at(site_0.y).xyz() + site_0_1.scale(-1));
      } else {
        XYZ anchor_anchor = get_vector_from_to(molecule->at(site_0.y).xyz(), molecule->at(site_1.y).xyz()).unit().scale(MAX_CARBOXY_BOND_LENGTH);
        molecule->at(site_1.x).set_xyz(molecule->at(site_1.y).xyz() + anchor_anchor);
        molecule->at(site_0.x).set_xyz(molecule->at(site_0.y).xyz() + anchor_anchor.scale(-1));
      }
    }
  } else if(site_type==ALDEHYDE) { //replace the O (x) atom with Q only
    for(int i=0; i<num_sites; i++) {
      INT3 site = indices->at(i);
      molecule->at(site.x).type = "Q"; //replace the O
      molecule->at(site.x).label = "Q";
    }
  } else if(site_type==A_P_CATECHOL) { //replace the C (x) atom with Q, and remove the CH3 groups (y/z and surrounding H)
    for(int i=0; i<num_sites; i++) {
      INT3 site = indices->at(i);
      molecule->at(site.x).type = "Q"; //replace the C
      molecule->at(site.x).label = "Q";
      //now find H atoms
      for(int j=0; j<num_atoms; j++) { //find H atoms bonded to this O
        if(molecule->at(j).type=="H" && molecule->at(j).keep) {
          if(bonded(molecule->at(site.y), molecule->at(j), MAX_HYDROGEN_BOND_LENGTH) || bonded(molecule->at(site.z), molecule->at(j), MAX_HYDROGEN_BOND_LENGTH)) {
            molecule->at(j).keep = false;
          }
        }
      }
      //finally, remove the C of CH3
      molecule->at(site.y).keep = false;
      molecule->at(site.z).keep = false;
    }
  } else if(site_type==DIALCOHOL) { //replace the dummy Z (z) atom with Q, and remove the H on y/x
    for(int i=0; i<num_sites; i++) {
      INT3 site = indices->at(i);
      molecule->at(site.z).type = "Q"; //replace the dummy Z
      molecule->at(site.z).label = "Q";
      //now find H atoms
      for(int j=0; j<num_atoms; j++) { //find H atoms bonded to this O
        if(molecule->at(j).type=="H" && molecule->at(j).keep) {
          if(bonded(molecule->at(site.y), molecule->at(j), MAX_HYDROGEN_BOND_LENGTH) || bonded(molecule->at(site.x), molecule->at(j), MAX_HYDROGEN_BOND_LENGTH)) {
            molecule->at(j).keep = false;
          }
        }
      }
    }
  } else {
    printf("ERROR: did not recognize site type %d\n", site_type);
    exit(EXIT_FAILURE);
  }
}

void find_sites(vector<ATOM> *molecule, vector<INT3> *indices, int site_type) {
  bool verbose = true;
  int num_atoms = molecule->size();
  vector<bool> already_checked;
  for(int i=0; i<num_atoms; i++) already_checked.push_back(false);
  if(site_type==CARBOXY) { //CARBOXY SITES
    for(int i=0; i<num_atoms; i++) {
      if(molecule->at(i).type=="O" && !already_checked.at(i)) { //found a new oxygen; potentially a connection site
  if(verbose) printf("found an O atom\n");
        INT3 trial_site;
        trial_site.x = i;
        already_checked.at(trial_site.x) = true;
        //is this a connection site? we need to find if there is one bond only (carbon - ignore hydrogens), and if that carbon is bonded to an oxygen with only one bond (that carbon - ignore hydrogens)
        bool cannot_be_site_first_O = false;
        int num_non_H_bonds_first_O = 0;
        for(int j=0; j<num_atoms && !cannot_be_site_first_O; j++) { //find how many non-H atoms are bonded to this O
          if(molecule->at(j).type!="H" && j!=trial_site.x) {
            if(bonded(molecule->at(trial_site.x), molecule->at(j), MAX_CARBOXY_BOND_LENGTH)) {
  if(verbose) printf("first O is bonded to %s with distance %.3f\n", molecule->at(j).type.c_str(), dist(molecule->at(trial_site.x), molecule->at(j)));
              num_non_H_bonds_first_O++;
              if(num_non_H_bonds_first_O>1 || molecule->at(j).type!="C") cannot_be_site_first_O = true;
              else if(molecule->at(j).type=="C" && !already_checked.at(j)) { //this is a potential C from a carboxy
                trial_site.y = j;
              }
            }
          }
        }
  if(verbose) printf("the first O atom has %d non-H bonds (exactly one required)\n", num_non_H_bonds_first_O);
        //at this point, if we did not find too many bonds, this is still a potential site, and we need to find the other oxygen - otherwise, bypass the rest and move on to test the next oxygen
        bool site_completed = false;
        for(int j=0; j<num_atoms && !cannot_be_site_first_O && !site_completed; j++) {
          if(molecule->at(j).type!="H" && j!=trial_site.y) {
            if(bonded(molecule->at(trial_site.y), molecule->at(j), MAX_CARBOXY_BOND_LENGTH)) {
              if(molecule->at(j).type=="O" && !already_checked.at(j)) { //this is a potential other O from a carboxy
  if(verbose) printf("found a second O atom that is bonded to this C\n");
                trial_site.z = j;
                already_checked.at(trial_site.z) = true;
                bool cannot_be_site_second_O = false;
                int num_non_H_bonds_second_O = 0;
                for(int k=0; k<num_atoms && !cannot_be_site_second_O; k++) { //find how many non-H atoms are bonded to this O
                  if(molecule->at(k).type!="H" && k!=trial_site.z) {
                    if(bonded(molecule->at(trial_site.z), molecule->at(k), MAX_CARBOXY_BOND_LENGTH)) {
  if(verbose) printf("second O is bonded to %s with distance %.3f\n", molecule->at(k).type.c_str(), dist(molecule->at(trial_site.z), molecule->at(k)));
                      num_non_H_bonds_second_O++;
                      if(num_non_H_bonds_second_O>1 || k!=trial_site.y) cannot_be_site_second_O = true;
                    }
                  }
                }
  if(verbose) printf("the second O atom has %d non-H bonds (exactly one required)\n", num_non_H_bonds_second_O);
                if(!cannot_be_site_second_O) {
                  //at this point, if we did not find too many bonds, this is a site!
                  indices->push_back(trial_site);
                  already_checked.at(trial_site.y) = true; //now safe to flag this carbon as having been checked already
                  site_completed = true;
                }
              }
            }
          }
        }
  if(verbose) {
    if(!site_completed) printf("this O atom did not constitute part of a carboxy group\n"); else printf("site found!\n");
  }
      }
    }
  } else if(site_type==BROMINE) { //BROMINE sites
    for(int i=0; i<num_atoms; i++) {
      if(molecule->at(i).type=="Br" && !already_checked.at(i)) { //found a new bromine; potentially a connection site
  if(verbose) printf("found a Br atom\n");
        bool site_completed = false;
        INT3 trial_site;
        trial_site.x = i;
        trial_site.z = -1; //this site is only one atom, and the one it is bonded to (y) - we don't need this third position - this output will be parsed appropriately later
        already_checked.at(trial_site.x) = true;
        //is this a connection site? we need to find if there is one bond only (doesn't matter what element)
        bool cannot_be_site_Br = false;
        int num_bonds_Br = 0;
        for(int j=0; j<num_atoms && !cannot_be_site_Br; j++) { //find how many non-H atoms are bonded to this Br
          if(molecule->at(j).type!="H" && j!=trial_site.x) {
            if(bonded(molecule->at(trial_site.x), molecule->at(j), MAX_BROMINE_BOND_LENGTH)) {
  if(verbose) printf("this Br is bonded to %s with distance %.3f\n", molecule->at(j).type.c_str(), dist(molecule->at(trial_site.x), molecule->at(j)));
              num_bonds_Br++;
              trial_site.y = j;
              if(num_bonds_Br>1) cannot_be_site_Br = true;
            }
          }
        }
  if(verbose) printf("the Br atom has %d bonds (exactly one required)\n", num_bonds_Br);
        if(!cannot_be_site_Br) {
          //at this point, if we did not find too many bonds, this is a site!
          indices->push_back(trial_site);
          site_completed = true;
        }
  if(verbose) {
    if(!site_completed) printf("this Br atom did not constitute a connection site\n"); else printf("site found!\n");
  }
      }
    }
  } else if(site_type==ALDEHYDE) { //ALDEHYDE SITES
    for(int i=0; i<num_atoms; i++) {
      if(molecule->at(i).type=="O" && !already_checked.at(i)) { //found a new oxygen; potentially a connection site
  if(verbose) printf("found an O atom\n");
        INT3 trial_site;
        trial_site.x = i;
        already_checked.at(trial_site.x) = true;
        //is this a connection site? we need to find if there is one bond only (carbon - CHECK FOR hydrogens), and if that carbon is bonded to a hydrogen
        bool cannot_be_site_first_O = false;
        int num_bonds_first_O = 0;
        for(int j=0; j<num_atoms && !cannot_be_site_first_O; j++) { //find how many atoms are bonded to this O
          if(
( molecule->at(j).type!="H" && j!=trial_site.x && bonded(molecule->at(trial_site.x), molecule->at(j), MAX_CARBOXY_BOND_LENGTH) )
||
( molecule->at(j).type=="H" && j!=trial_site.x && bonded(molecule->at(trial_site.x), molecule->at(j), MAX_HYDROGEN_BOND_LENGTH) )
) {
  if(verbose) printf("first O is bonded to %s with distance %.3f\n", molecule->at(j).type.c_str(), dist(molecule->at(trial_site.x), molecule->at(j)));
            num_bonds_first_O++;
            if(num_bonds_first_O!=1) cannot_be_site_first_O = true;
            else if(molecule->at(j).type=="C" && !already_checked.at(j)) { //this is a potential C from an aldehyde
              trial_site.y = j;
            }
          }
        }
  if(verbose) printf("the first O atom has %d bonds (exactly one required)\n", num_bonds_first_O);
        //at this point, if we did not find too many bonds, this is still a potential site, and we need to find the hydrogen AND EXACTLY ONE OTHER
        bool site_completed = false;
        int num_H = 0, num_non_H = 0;
        int id_H = -1;
        for(int j=0; j<num_atoms && !cannot_be_site_first_O && !site_completed; j++) {
          if(molecule->at(j).type=="H" && !already_checked.at(j) && j!=trial_site.x && j!=trial_site.y && bonded(molecule->at(trial_site.y), molecule->at(j), MAX_HYDROGEN_BOND_LENGTH)) {
            //this is a potential H to complete the aldehyde
  if(verbose) printf("found an H atom that is bonded to this C\n");
            num_H++;
            id_H = j;
          } else if(molecule->at(j).type!="H" && j!=trial_site.x && j!=trial_site.y && bonded(molecule->at(trial_site.y), molecule->at(j), MAX_CARBOXY_BOND_LENGTH)) {
            //this is a potential non-H to complete the aldehyde
  if(verbose) printf("found a non-H atom that is bonded to this C\n");
            num_non_H++;
          }
        }
        //it's only an aldehyde if both num_H and num_non_H are 1
        if(num_H==1 && num_non_H==1) {
          trial_site.z = id_H;
          already_checked.at(trial_site.z) = true;
          //at this point, this is a site!
          already_checked.at(trial_site.y) = true; //now safe to flag this carbon as having been checked already
          indices->push_back(trial_site);
          site_completed = true;
        }
  if(verbose) {
    if(!site_completed) printf("this O atom did not constitute part of an aldehyde group\n"); else printf("site found!\n");
  }
      }
    }
  } else if(site_type==A_P_CATECHOL) { //A_P_CATECHOL SITES
    for(int i=0; i<num_atoms; i++) {
      if(molecule->at(i).type=="C" && !already_checked.at(i)) { //found a new carbon; potentially a connection site
  if(verbose) printf("found a C atom\n");
        INT3 trial_site;
        trial_site.x = i;
        already_checked.at(trial_site.x) = true;
        //is this a connection site? we need to find if there are four bonds: two to CH3 and two to OC
        bool cannot_be_site_first_C = false;
        int num_non_H_bonds_first_C = 0;
        vector<int> bonded_ids;
        for(int j=0; j<num_atoms && !cannot_be_site_first_C; j++) { //find how many non-H atoms are bonded to this O
          if(
molecule->at(j).type!="H" && j!=trial_site.x && bonded(molecule->at(trial_site.x), molecule->at(j), MAX_CARBOXY_BOND_LENGTH)
) {
  if(verbose) printf("first C (%d) is bonded to %s (%d) with distance %.3f\n", trial_site.x, molecule->at(j).type.c_str(), j, dist(molecule->at(trial_site.x), molecule->at(j)));
            num_non_H_bonds_first_C++;
            bonded_ids.push_back(j);
          }
        }
        if(num_non_H_bonds_first_C!=4) cannot_be_site_first_C = true;
        else { //find out if those 4 bonds are the 4 we are expecting
          int num_bonds_to_CH3_first_C = 0;
          int num_bonds_to_OC_first_C = 0;
          vector<int> bonded_CH3_ids;
          vector<int> OC_C_ids;
          bool fail = false;
          for(int j=0; j<num_non_H_bonds_first_C; j++) {
            int id = bonded_ids.at(j);
            if(molecule->at(id).type=="O" && (!already_checked.at(id))) { //test for OC
              int num_other_bonds = 0;
              int other_bond_id = -1;
              for(int k=0; k<num_atoms; k++) {
                if(
( k!=id && k!=trial_site.x && molecule->at(k).type=="H" && bonded(molecule->at(id), molecule->at(k), MAX_HYDROGEN_BOND_LENGTH) )
||
( k!=id && k!=trial_site.x && molecule->at(k).type!="H" && bonded(molecule->at(id), molecule->at(k), MAX_CARBOXY_BOND_LENGTH) )
                ) {
  if(verbose) printf("an O on this first C is bonded to %s with distance %.3f\n", molecule->at(k).type.c_str(), dist(molecule->at(id), molecule->at(k)));
                  other_bond_id = k;
                  num_other_bonds++;
                }
              }
              if(num_other_bonds!=1) fail=true; else {OC_C_ids.push_back(other_bond_id); num_bonds_to_OC_first_C++;}
            } else if(molecule->at(id).type=="C") { //test for CH3
              int num_H = 0;
              for(int k=0; k<num_atoms; k++) {
                if(
k!=id && k!=trial_site.x && molecule->at(k).type=="H" && bonded(molecule->at(id), molecule->at(k), MAX_HYDROGEN_BOND_LENGTH)
                ) {
  if(verbose) printf("a C on this first C is bonded to %s with distance %.3f\n", molecule->at(k).type.c_str(), dist(molecule->at(id), molecule->at(k)));
                  num_H++;
                }
              }
              if(num_H!=3) fail=true; else {bonded_CH3_ids.push_back(bonded_ids.at(j)); num_bonds_to_CH3_first_C++;}
            }
          }
          if(num_bonds_to_OC_first_C!=2 || num_bonds_to_CH3_first_C!=2 || fail) {
printf("NOTE: this first C fails here with %d bonds to OC and %d bonds to CH3 groups\n", num_bonds_to_OC_first_C, num_bonds_to_CH3_first_C);
            cannot_be_site_first_C = true;
          } else if(!bonded(molecule->at(OC_C_ids.at(0)), molecule->at(OC_C_ids.at(1)), MAX_CARBOXY_BOND_LENGTH)) cannot_be_site_first_C = true;
          else { //this is a C from an a_p_catechol group
            trial_site.y = bonded_CH3_ids.at(0);
            trial_site.z = bonded_CH3_ids.at(1);
            already_checked.at(trial_site.y) = true;
            already_checked.at(trial_site.z) = true;
            //at this point, this is a site!
            indices->push_back(trial_site);
          }
        }
  if(verbose) {
    if(cannot_be_site_first_C) printf("this C atom did not constitute part of an a_p_catechol group\n"); else printf("site found!\n");
  }
      }
    }
  } else if(site_type==DIALCOHOL) { //DIALCOHOL SITES
    for(int i=0; i<num_atoms; i++) {
      if(molecule->at(i).type=="O" && !already_checked.at(i)) { //found a new oxygen; potentially a connection site
  if(verbose) printf("found an O atom\n");
        INT3 trial_site;
        trial_site.x = i;
        already_checked.at(trial_site.x) = true;
        //is this a connection site? we need to find if it is an OH, bonded to a C, which is bonded to another C-O-H
        bool cannot_be_site_first_O = false;
        int num_C_bonds_first_O = 0;
        int num_H_bonds_first_O = 0;
        int first_C_id = -1;
        for(int j=0; j<num_atoms && !cannot_be_site_first_O; j++) { //find bonds
          if(
(!already_checked.at(j)) && molecule->at(j).type=="C" && j!=trial_site.x && bonded(molecule->at(trial_site.x), molecule->at(j), MAX_CARBOXY_BOND_LENGTH)
) {
  if(verbose) printf("first O is bonded to %s with distance %.3f\n", molecule->at(j).type.c_str(), dist(molecule->at(trial_site.x), molecule->at(j)));
            num_C_bonds_first_O++;
            first_C_id = j;
          } else if(
(!already_checked.at(j)) && molecule->at(j).type=="H" && j!=trial_site.x && bonded(molecule->at(trial_site.x), molecule->at(j), MAX_HYDROGEN_BOND_LENGTH)
) {
  if(verbose) printf("first O is bonded to %s with distance %.3f\n", molecule->at(j).type.c_str(), dist(molecule->at(trial_site.x), molecule->at(j)));
            num_H_bonds_first_O++;
          }
        }
        if(num_C_bonds_first_O!=1 || num_H_bonds_first_O!=1) cannot_be_site_first_O = true;
        else { //we know this is an H-O-C, but does the C connect us to another C-O-H? check bonding environment of first C
          already_checked.at(first_C_id) = true; //we can flag this first C as being considered, since no other C-O-H group will contain it
          vector<int> second_C_id_possibilities;
          for(int j=0; j<num_atoms && !cannot_be_site_first_O; j++) { //find bonds
            if(
(!already_checked.at(j)) && molecule->at(j).type=="C" && j!=first_C_id && bonded(molecule->at(first_C_id), molecule->at(j), MAX_CARBOXY_BOND_LENGTH)
            ) {
if(verbose) printf("first C is bonded to %s with distance %.3f\n", molecule->at(j).type.c_str(), dist(molecule->at(first_C_id), molecule->at(j)));
              second_C_id_possibilities.push_back(j);
            }
          }
          if(second_C_id_possibilities.size()<1) cannot_be_site_first_O = true;
          else { //there are some second C choices to consider
            bool site_found = false; //allow to quit these loops early as soon as a site is found
            for(int j=0; j<second_C_id_possibilities.size() && !site_found; j++) {
              //this second C must be bonded to an OH
              int candidate_second_C_id = second_C_id_possibilities.at(j);
              for(int k=0; k<num_atoms && !cannot_be_site_first_O && !site_found; k++) { //find bonds
                if(
(!already_checked.at(k)) && molecule->at(k).type=="O" && bonded(molecule->at(candidate_second_C_id), molecule->at(k), MAX_CARBOXY_BOND_LENGTH)
              ) {
if(verbose) printf("candidate second C is bonded to %s with distance %.3f\n", molecule->at(k).type.c_str(), dist(molecule->at(candidate_second_C_id), molecule->at(k)));
                  //if this k O is an OH, we are good to go
                  int num_H = 0;
                  for(int l=0; l<num_atoms && !cannot_be_site_first_O && !site_found; l++) { //find bonds
                    if( molecule->at(l).type=="H" && bonded(molecule->at(k), molecule->at(l), MAX_HYDROGEN_BOND_LENGTH) ) {
                      num_H++;
if(verbose) printf("candidate second C's O is bonded to %s with distance %.3f\n", molecule->at(l).type.c_str(), dist(molecule->at(l), molecule->at(k)));
                    }
                  }
                  if(num_H==1) { //satisfied that this is a site
                    trial_site.y = k;
                    //we need to find where the z site is, because it is not an atom!
                    //first find the centroid of site.x, site.y and the two carbons, first_C_id and candidate_second_C_id
                    XYZ site_centroid = ( molecule->at(candidate_second_C_id).xyz() + molecule->at(first_C_id).xyz() + molecule->at(trial_site.y).xyz() + molecule->at(trial_site.x).xyz() ).scale(0.25);
                    double atom_centroid_distance = dist(site_centroid.x,site_centroid.y,site_centroid.z,molecule->at(trial_site.x));
                    XYZ projection_of_centroid = project_onto_line(site_centroid, molecule->at(trial_site.x).xyz(), molecule->at(trial_site.y).xyz());
                    XYZ new_atom_vector = (projection_of_centroid-site_centroid).unit();
                    XYZ new_atom_pos = site_centroid+(new_atom_vector.scale(atom_centroid_distance));
                    //make a new atom at this position
                    ATOM new_atom(new_atom_pos,"Z",0.0);
                    molecule->push_back(new_atom);
                    already_checked.push_back(true);
                    trial_site.z = num_atoms;
                    num_atoms++;
                    already_checked.at(trial_site.y) = true;
                    site_found = true;
                    //at this point, this is a site!
                    indices->push_back(trial_site);
                  }
                }
              }
            }
            //at this stage, this is a dead end
            if(!site_found) cannot_be_site_first_O = true;
          }
        }
  if(verbose) {
    if(cannot_be_site_first_O) printf("this O atom did not constitute part of a dialcohol group\n"); else printf("site found!\n");
  }
      }
    }
  } else {
    printf("ERROR: did not recognize site type %d\n", site_type);
    exit(EXIT_FAILURE);
  }
}

void read_molecule(string filename, vector<ATOM> *atom_vector) {
  //this function inspects the filename, and decides which read function to call - for xyz or incon format
  string extension_start = ".";
  unsigned pos = filename.rfind(extension_start);
  if(pos==-1) {
    printf("ERROR: could not identify extension for the input file named \"%s\"\n", filename.c_str());
    exit(EXIT_FAILURE);
  }
  string xyz_format = "xyz";
  unsigned pos2 = filename.find(xyz_format, pos);
  if(pos2!=-1) {
    //we found xyz in the file extension - treat it as an xyz format file
    read_xyz(filename, atom_vector);
  } else {
    //not xyz .. is it incon?
    string incon_format = "incon";
    unsigned pos3 = filename.find(incon_format, pos);
    if(pos3!=-1) {
      //we found incon in the file extension - treat it as an incon format file
      read_incon(filename, atom_vector);
    } else {
      printf("ERROR: could not parse either \"%s\" or \"%s\" in the extension for the input file named \"%s\"\n", xyz_format.c_str(), incon_format.c_str(), filename.c_str());
      exit(EXIT_FAILURE);
    }
  }
}

void read_incon(string filename, vector<ATOM> *atom_vector) {
  bool verbose = false;
  int status = 0;
  //open input
  FILE *input;
	input = fopen(filename.c_str(), "r");
	if(input==NULL) {
		printf("ERROR: could not open file with name %s. Please try again.\n", filename.c_str());
		exit(EXIT_FAILURE);
	}
  //read input
  int num_atoms_read = 0;
  char ch1[MAX_CHAR_ARRAY_LENGTH] = "";
  bool done_reading_atoms = false;
  bool started_reading_atoms = false;
//  char *atom_line_second_column = new char[MAX_CHAR_ARRAY_LENGTH];
//  atom_line_second_column = "core";
  char atom_line_second_column[MAX_CHAR_ARRAY_LENGTH] = "core";
  while(fgets(ch1, MAX_CHAR_ARRAY_LENGTH, input)!=NULL && !done_reading_atoms) {

//    vector<string> token;
//    token = split(ch1," \r\t");
    double x=0,y=0,z=0;
    char *name = new char[MAX_CHAR_ARRAY_LENGTH];
    char *dummy = new char[MAX_CHAR_ARRAY_LENGTH];
    status = sscanf(ch1, "%s %s %lf %lf %lf", name, dummy, &x, &y, &z);
//    if(token.size()==5) { //a line comprising an atom name, a string, and three coords
    if(strcmp(atom_line_second_column,dummy)==0) { //line contains an atom in incon format!
      started_reading_atoms = true;
      ATOM a;
//      char *name = new char[MAX_CHAR_ARRAY_LENGTH];
//      char *dummy = new char[MAX_CHAR_ARRAY_LENGTH];
//      status = sscanf(ch1, "%s %s %lf %lf %lf", name, dummy, &a.x, &a.y, &a.z);
      a.x = x; a.y = y; a.z = z;

      string name_string(name);
      a.label = name_string;
      //parse out any digits from the name of the element
      int name_length = name_string.length();
      int digit_pos = -1;
      for(int i=0; i<name_length && digit_pos==-1; i++) {
        if(isdigit(name_string[i])) digit_pos = i;
      }
      if(digit_pos==0) { //first character was a digit - this is no good
        if(started_reading_atoms) done_reading_atoms = true; //if we have already started parsing atom data, and encounter a line without a valid atom name, then we should stop
      } else {
        if(digit_pos>0) {
          a.type = name_string.substr(0,digit_pos);
        } else a.type = name_string;
        //good to go - store in vector
        a.radius = get_rad(a.type);
if(verbose) printf("atom %s assigned radius %.3f\n", a.type.c_str(), a.radius);
        delete[] name;
        delete[] dummy;    
        atom_vector->push_back(a);
      }
    } else {
//printf("DEBUG: line did not have 5 tokens (%d tokens): %s\n", num_tokens, ch1);
//printf("DEBUG: line did not have %s as the second component: %s\n", atom_line_second_column, ch1);
      if(started_reading_atoms) done_reading_atoms = true; //if we have already started parsing atom data, and encounter a line without atom data, then we should stop
    }
  }
  if(!done_reading_atoms) {
    printf("ERROR: finished reading data but did not finishing parsing atom info correctly - managed to read %d atoms\n", (int)(atom_vector->size()));
    exit(EXIT_FAILURE);
  }  
}

void read_xyz(string filename, vector<ATOM> *atom_vector) {
  bool verbose = false;
  int status = 0;
  //open input
  FILE *input;
	input = fopen(filename.c_str(), "r");
	if(input==NULL) {
		printf("ERROR: could not open file with name %s. Please try again.\n", filename.c_str());
		exit(EXIT_FAILURE);
	}
  //read input
  char ch1[MAX_CHAR_ARRAY_LENGTH] = "";
  int num_input_atoms = -1;
  if(fgets(ch1, MAX_CHAR_ARRAY_LENGTH, input)!=NULL) {
    string str = string(ch1);
    int pos=0;
    char c = str[pos];
    while(c<=0) { //for some reason, xyz files from marvinsketch have three characters at the beginning with ascii codes <0 - very strange, but we need to remove these to succesfully parse any data
      pos++;
      c = str[pos];
    }
    char ch2[MAX_CHAR_ARRAY_LENGTH] = "";
    str.copy(ch2, str.size()-pos, pos);
    status = sscanf(ch2, "%d", &num_input_atoms);
    if(num_input_atoms<=0) {
      printf("ERROR: number of atoms to be read from file is %d<=0\n", num_input_atoms);
      exit(EXIT_FAILURE);
    }
  } else {
    printf("ERROR: could not read string\n");
    exit(EXIT_FAILURE);
  }
  search_for_char(input, '\n'); //now we are on a line containing an atom
  for(int i=0; i<num_input_atoms; i++) {
    ATOM a;
    char ch[MAX_CHAR_ARRAY_LENGTH];
    if(fgets(ch, MAX_CHAR_ARRAY_LENGTH, input)!=NULL) {
      string str = string(ch);
      int pos=0;
      char c = str[pos];
      while(c<=0) { //for some reason, xyz files from marvinsketch have three characters at the beginning with ascii codes <0 - very strange, but we need to remove these to succesfully parse any data
        pos++;
        c = str[pos];
      }
      char ch2[MAX_CHAR_ARRAY_LENGTH] = "";
      str.copy(ch2, str.size()-pos, pos);
      char *name = new char[MAX_CHAR_ARRAY_LENGTH];      
      status = sscanf(ch2, "%s %lf %lf %lf", name, &a.x, &a.y, &a.z);
      string name_string = string(name);
      a.label = name_string;
      //parse out any digits from the name of the element
      int name_length = name_string.length();
      int digit_pos = -1;
      for(int i=0; i<name_length && digit_pos==-1; i++) {
        if(isdigit(name_string[i])) digit_pos = i;
      }
      if(digit_pos==0) { //first character was a digit - this is no good
        printf("ERROR: could not parse label from atom beginning with a digit in read_xyz: %s: %s\n", filename.c_str(), name_string.c_str());
        exit(EXIT_FAILURE);
      } else {
        if(digit_pos>0) {
          a.type = name_string.substr(0,digit_pos);
        } else a.type = name_string;
//printf("DEBUG: with name_string=%s, digit_pos=%d, therefore a.type=%s\n", name_string.c_str(), digit_pos, a.type.c_str());
      }
      a.radius = get_rad(a.type);
if(verbose) printf("atom %s assigned radius %.3f\n", a.type.c_str(), a.radius);
      delete[] name;
    } else {
      printf("ERROR: could not read expected atom coord string from %s - %d out of %d atom coords were read\n", filename.c_str(), i, num_input_atoms);
      exit(EXIT_FAILURE);
    }
    char match = 0;
    atom_vector->push_back(a);
  }
  fclose(input);
}

void write_xyz(FILE *outfile, vector<ATOM> *atom_vector, string input_name) {
  int num_atoms_total = atom_vector->size();
  int num_atoms_to_print = 0;
  for(int i=0; i<num_atoms_total; i++) {
    if(atom_vector->at(i).keep) num_atoms_to_print++;
  }
  fprintf(outfile, "%d\nabstract molecule created from input %s\n", num_atoms_to_print, input_name.c_str());
  for(int i=0; i<num_atoms_total; i++) {
    if(atom_vector->at(i).keep) {
      ATOM a = atom_vector->at(i);
      fprintf(outfile, "%s %.6f %.6f %.6f\n", a.label.c_str(), a.x, a.y, a.z);
//      fprintf(outfile, "%s %.6f %.6f %.6f\n", a.type.c_str(), a.x, a.y, a.z);
    }
  }
  fclose(outfile);
}
/*
void search_for_char(FILE *f, char target) {
	char c = getc(f);
	while(c!=target && c!=EOF) {
		c = getc(f);
	}
	if(c==EOF) {
		printf("ERROR: The required character (%c, with int value %d) was not found in this file.\n", target, (int)target);
		exit(EXIT_FAILURE);
	}
}
*/
bool bonded(ATOM a, ATOM b, double d) {
  if(dist(a,b)<d) return true; else return false;
}

double dist(ATOM a, ATOM b) {
  XYZ a_xyz(a.x,a.y,a.z);
  XYZ b_xyz(b.x,b.y,b.z);
  return (a_xyz-b_xyz).magnitude();
/*
  double dx = a.xyz.x-b.xyz.x;
  double dy = a.xyz.y-b.xyz.y;
  double dz = a.xyz.z-b.xyz.z;
  return sqrt((dx*dx)+(dy*dy)+(dz*dz));
*/
}

double dist(double x, double y, double z, ATOM b) {
  XYZ a_xyz(x,y,z);
  XYZ b_xyz(b.x,b.y,b.z);
  return (a_xyz-b_xyz).magnitude();
/*
  double dx = x-b.xyz.x;
  double dy = y-b.xyz.y;
  double dz = z-b.xyz.z;
  return sqrt((dx*dx)+(dy*dy)+(dz*dz));
*/
}

double get_rad(string s) { //CCDC radii
  if(s=="H") return 1.09;
  else if(s=="He") return 1.4;
  else if(s=="Li") return 1.82;
  else if(s=="Be") return 2;
  else if(s=="B") return 2;
  else if(s=="C") return 1.7;
  else if(s=="N") return 1.55;
  else if(s=="O") return 1.52;
  else if(s=="F") return 1.47;
  else if(s=="Ne") return 1.54;
  else if(s=="Na") return 2.27;
  else if(s=="Mg") return 1.73;
  else if(s=="Al") return 2;
  else if(s=="Si") return 2.1;
  else if(s=="P") return 1.8;
  else if(s=="S") return 1.8;
  else if(s=="Cl") return 1.75;
  else if(s=="Ar") return 1.88;
  else if(s=="K") return 2.75;
  else if(s=="Ca") return 2;
  else if(s=="Sc") return 2;
  else if(s=="Ti") return 2;
  else if(s=="V") return 2;
  else if(s=="Cr") return 2;
  else if(s=="Mn") return 2;
  else if(s=="Fe") return 2;
  else if(s=="Co") return 2;
  else if(s=="Ni") return 1.63;
  else if(s=="Cu") return 1.4;
  else if(s=="Zn") return 1.39;
  else if(s=="Ga") return 1.87;
  else if(s=="Ge") return 2;
  else if(s=="As") return 1.85;
  else if(s=="Se") return 1.9;
  else if(s=="Br") return 1.85;
  else if(s=="Kr") return 2.02;
  else if(s=="Rb") return 2;
  else if(s=="Sr") return 2;
  else if(s=="Y") return 2;
  else if(s=="Zr") return 2;
  else if(s=="Nb") return 2;
  else if(s=="Mo") return 2;
  else if(s=="Tc") return 2;
  else if(s=="Ru") return 2;
  else if(s=="Rh") return 2;
  else if(s=="Pd") return 1.63;
  else if(s=="Ag") return 1.72;
  else if(s=="Cd") return 1.58;
  else if(s=="In") return 1.93;
  else if(s=="Sn") return 2.17;
  else if(s=="Sb") return 2;
  else if(s=="Te") return 2.06;
  else if(s=="I") return 1.98;
  else if(s=="Xe") return 2.16;
  else if(s=="Cs") return 2;
  else if(s=="Ba") return 2;
  else if(s=="La") return 2;
  else if(s=="Ce") return 2;
  else if(s=="Pr") return 2;
  else if(s=="Nd") return 2;
  else if(s=="Pm") return 2;
  else if(s=="Sm") return 2;
  else if(s=="Eu") return 2;
  else if(s=="Gd") return 2;
  else if(s=="Tb") return 2;
  else if(s=="Dy") return 2;
  else if(s=="Ho") return 2;
  else if(s=="Er") return 2;
  else if(s=="Tm") return 2;
  else if(s=="Yb") return 2;
  else if(s=="Lu") return 2;
  else if(s=="Hf") return 2;
  else if(s=="Ta") return 2;
  else if(s=="W") return 2;
  else if(s=="Re") return 2;
  else if(s=="Os") return 2;
  else if(s=="Ir") return 2;
  else if(s=="Pt") return 1.72;
  else if(s=="Au") return 1.66;
  else if(s=="Hg") return 1.55;
  else if(s=="Tl") return 1.96;
  else if(s=="Pb") return 2.02;
  else if(s=="Bi") return 2;
  else if(s=="Po") return 2;
  else if(s=="At") return 2;
  else if(s=="Rn") return 2;
  else if(s=="Fr") return 2;
  else if(s=="Ra") return 2;
  else if(s=="Ac") return 2;
  else if(s=="Th") return 2;
  else if(s=="Pa") return 2;
  else if(s=="U") return 1.86;
  else if(s=="Np") return 2;
  else if(s=="Pu") return 2;
  else if(s=="Am") return 2;
  else if(s=="Cm") return 2;
  else if(s=="Bk") return 2;
  else if(s=="Cf") return 2;
  else if(s=="Es") return 2;
  else if(s=="Fm") return 2;
  else if(s=="Md") return 2;
  else if(s=="No") return 2;
  else if(s=="Lr") return 2;
  else if(s=="Rf") return 2;
  else if(s=="Db") return 2;
  else if(s=="Sg") return 2;
  else if(s=="Bh") return 2;
  else if(s=="Hs") return 2;
  else if(s=="Mt") return 2;
  else if(s=="Ds") return 2;
  else {
    printf("WARNING: could not find radius for element type %s - using default of %.3f\n", s.c_str(), DEFAULT_RADIUS);
    return DEFAULT_RADIUS;
  }
}

INT3::INT3() {
  x = 0; y = 0; z = 0;
}

double radians_to_degrees(double r) {return 360.0*r/(2.0*PI);}
/*
vector<string> split(string line, string delims) {
  vector<string> token;
  string temp=line;
  int ndx;
  while(!temp.empty()) {
    ndx = temp.find_first_of(delims);
    if(ndx>0) {
      token.push_back(temp.substr(0,ndx));
    }
    else if(ndx ==-1) {
      token.push_back(temp);
      return token;
    }
    temp=temp.substr(ndx+1);
  }
  return token;
}
*/
//OLD CODE BELOW HERE

/*
ATOM::ATOM(double x0, double y0, double z0, string s0, double r0) {
  XYZ xyz1(x0,y0,z0);
  xyz = xyz1;
  name = s0;
  rad = r0;
  keep = true;
}

ATOM::ATOM(XYZ xyz0, string s0, double r0) {
  xyz = xyz0;
  name = s0;
  rad = r0;
  keep = true;
}

ATOM::ATOM() {
  name = "";
  rad = 0;
  keep = true;
}
*/

/*
bool increasing_rmsd (FIT i, FIT j) { return (i.rmsd<j.rmsd); }
*/

//ATOM
/*
class ATOM {
public:
  bool keep;
  XYZ xyz; //the Cartesian position of this ATOM
  double rad;
  string name;
  //constructor
  ATOM(double x0, double y0, double z0, string s0, double r0); //separate Cartesian values constructor
  ATOM(XYZ xyz0, string s0, double r0); //XYZ constructor
  ATOM();
};
*/

//FIT
/*
class FIT{
public:
  vector<ATOM> molecule;
  double rmsd;
  int perm_ID;
};
*/

/*
//XYZ
class XYZ {
public:
  double x, y, z;
  //constructor
  XYZ(double x0, double y0, double z0);
  XYZ();
  //methods
  double angle_between(XYZ other);
  double dot_product(XYZ other);
  XYZ cross_product(XYZ other);
  double magnitude();
  XYZ scale(double factor);
  XYZ unit();
  XYZ operator-(XYZ other); 
  XYZ operator+(XYZ other);
};
*/

/*
XYZ::XYZ(double x0, double y0, double z0) {
  x = x0; y = y0; z = z0;
}

XYZ::XYZ() {
  x = 0; y = 0; z = 0;
}

double XYZ::angle_between(XYZ other) {
  double cos_angle = dot_product(other)/(magnitude()*other.magnitude());
  if(cos_angle>1) {
    cos_angle = 1;
  } else if(cos_angle<-1) {
    cos_angle = -1; //these steps are necessary for rounding issues when cos_angle should be exactly 1 or -1, but falls very slightly on the wrong side
  }
  double output = acos(cos_angle);
  if(isnan(output)) {
    return 0;
  } else return output; //this step is necessary in the case of identity, which produces NAN
}

double XYZ::magnitude(){
  return sqrt(x*x+y*y+z*z);
}

XYZ XYZ::cross_product(XYZ other){
  return XYZ(y*other.z-z*other.y,z*other.x-x*other.z,x*other.y-y*other.x);
}

double XYZ::dot_product(XYZ other){
  return x*other.x + y*other.y + z*other.z; 
}

XYZ XYZ::scale(double factor){
  return XYZ(factor*x, factor*y, factor*z);
}

XYZ XYZ::unit(){
  return XYZ(x/magnitude(),y/magnitude(),z/magnitude());
}

XYZ XYZ::operator-(XYZ other){
   return XYZ(x-other.x, y-other.y, z-other.z); 
}

XYZ XYZ::operator+(XYZ other){
   return XYZ(x + other.x, y + other.y, z + other.z); 
}

XYZ get_vector_from_to(XYZ start, XYZ end) {return end-start;}

XYZ midpoint(XYZ a, XYZ b) {return (a+b).scale(0.5);}

XYZ project_onto_line(XYZ initial, XYZ line_start, XYZ line_end) {
  XYZ projection;
  XYZ initial_to_line_start_vector = get_vector_from_to(initial, line_start);
  XYZ line_vector = get_vector_from_to(line_start, line_end);
  double line_length = line_vector.magnitude();
  double movement_along_line = (-1*initial_to_line_start_vector.dot_product(line_vector))/(line_length*line_length);
  projection.x = line_start.x + (movement_along_line*line_vector.x);
  projection.y = line_start.y + (movement_along_line*line_vector.y);
  projection.z = line_start.z + (movement_along_line*line_vector.z);
  return projection;
}
*/

