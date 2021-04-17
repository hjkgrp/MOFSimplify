//-----Richard Luis Martin, Lawrence Berkeley National Laboratory, 2013
/* this file contains classes and functions required for handling topology (net) information */

#ifndef NET_H
#define NET_H

#include <string>
#include <vector>
#include <sstream>
using namespace std;
#include "geometry.h"

/* defines and consts */
#define MAX_DUPE_CON_ANGLE_DIFF PI/36.0
#define MAX_DUPE_CON_MAG_DIFF 5
#define MIN_UNIT_ATOM_NETWORK_SIDE_LENGTH 0.5
#define NEARLY_PARALLEL_ANGLE_DEGREES 10.0
#define GENERAL_BOND_DISTANCE 1.9 //if two atoms are within this distance, they are considered to be bonded
#define IS_PART_OF_SITE_DISTANCE 2.5 //if an atom and a site marker are within this distance, the atom is considered to be part of the site - if two atoms are both part of some site, then they will not be counted as collisions (because they may be forming a legitimate bond)
#define TOLERANCE_PROPORTION 1.05 //5% tolerance in distance measurements etc.
#define SAME_POSITION_DISTANCE 0.1 //how far apart sites should be in distinct aligned molecules in order to the alignments to be considered unique (note: setting the tolerance too high can cause the correct alignment to be discarded ...)
#define DISTANCE_TOLERANCE 0.01 //used for distances etc. when an absolute difference between things is calculated
#define SITE_DEVIATION_WARNING 1.5 //1.0 //if sites are further apart than this in the best framework, issue a warning

/* class to store vertex position and connectivity information */
class VERTEX {
public:
  XYZ abc; //the fractional position of this VERTEX
  int expected_edges;
  std::vector<XYZ> edges; //data on the edges for this vertex
  std::vector<XYZ> dummy_edges; //a single dummy edge is used only for 2-c vertices, as a means of providing an orthogonal direction to use when fitting 2-c molecules (which also have a corresponding dummy site) - it is stored as a vector to enable easy checking of whether it is declared
  std::string name;
  /* constructors */
  VERTEX(double a0, double b0, double c0); //fractional constructor
  VERTEX();
};

/* class to store vertex connection information */
class CONNECTION{
public:
  int v1, v2;
  int e1, e2;
  int a, b, c; //the periodicity of this connection

  //constructor
  CONNECTION(int v1_, int v2_, int e1_, int a_, int b_, int c_);
};

/* class to store a non-periodic collection of atoms, i.e. a molecule, including the position of connection sites to attach to other molecules, and the order of sites used w.r.t. the underlying net */
class MOLECULE {
public:
  vector<XYZ> atoms_xyz;
  vector<string> atoms_type;
  vector<string> atoms_label;
  vector<int> sites;
  vector<int> dummy_sites;
  vector<int> permutation; //stores the order in which the sites correspond to the edges of a vertex in the net, i.e. edge[i] corresponds to site permutation[i]
  XYZ com; //centre of mass, or just centroid, of site positions
};

/* class to store possible orientations of a molecule to some reference points */
class FIT{
public:
  MOLECULE fit;
  double rmsd;
  int perm_ID;
};

/* method declarations */
XYZ get_mol_site_CoM(MOLECULE *mol);
bool bonded_xyz(XYZ xyz1, XYZ xyz2);
bool bonded_xyz(XYZ xyz1, XYZ xyz2, double threshold);
bool matches(CONNECTION ci, CONNECTION cj);
void swap(vector<int> *vec, int a, int b);
void recursive_visit_vertices(int this_v, vector<CONNECTION> *two_way_connections, vector<bool> *visited_vertex);
bool recursive_find_loops(int this_v, int a, int b, int c, XYZ xyz, vector<CONNECTION> *two_way_connections, vector<MOLECULE> *assembled_molecules, vector<bool> *edge_traversed, vector<int> *vertex_visit_a, vector<int> *vertex_visit_b, vector<int> *vertex_visit_c, vector<XYZ> *vertex_visit_xyz, vector<bool> *vertex_visited, vector<int> *loops_a, vector<int> *loops_b, vector<int> *loops_c, vector<XYZ> *loops_xyz, vector<int> *unit_cell_vector_IDs, vector<XYZ> *unit_cell_vectors, int dimensionality);
bool try_unit_cell_vector_assign(int a, int b, int c, XYZ xyz, vector<int> *loops_a, vector<int> *loops_b, vector<int> *loops_c, vector<XYZ> *loops_xyz, vector<int> *unit_cell_vector_IDs, vector<XYZ> *unit_cell_vectors, int num_uc_vectors_required);
bool loop_is_unique(int a, int b, int c, vector<int> *loops_a, vector<int> *loops_b, vector<int> *loops_c);
MOLECULE translate(MOLECULE m, XYZ t);
MOLECULE rotate(MOLECULE orig, double mat[3][3]);
bool molecule_alignment_chemistry_is_unique(MOLECULE *mol, vector<MOLECULE> *existing);
bool overlaps_xyz(XYZ xyz1, XYZ xyz2);
void permute(vector<int> vec, int pos, vector< vector<int> > *out_vec);
int recursive_test_dimensionality(int this_v, int a, int b, int c, vector<CONNECTION> *two_way_connections, vector<bool> *edge_traversed, vector<int> *vertex_visit_a, vector<int> *vertex_visit_b, vector<int> *vertex_visit_c, vector<bool> *vertex_visited, vector<bool> *is_periodic);
void get_com(MOLECULE *m);
int assign_dummy_site(MOLECULE *m);

/* TO_FIX from file_and_string_ops.h - need to sort out any possible redundancies here */
void write_connections(FILE *f, vector<CONNECTION> *two_way_connections, vector<MOLECULE> *assembled_molecules, string name);
void write_molecule(FILE *rotmol, MOLECULE *rotated, string name, int vertex_ID, int sym_op, bool write_out_sites);
void write_framework(FILE *f, MOLECULE *framework, string name);
void search_for_char(FILE *f, char c);
string convertToString(int number);
string toLowerCase(string s);

/* TO_FIX - REDUNDANT RESOLVE THIS */
vector<XYZ> get_periodic_images_of_uc_abc_position(XYZ e_abc);

#endif

