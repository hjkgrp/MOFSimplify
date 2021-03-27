//-----Richard Luis Martin, Lawrence Berkeley National Laboratory, 2013
//-----a code to construct frameworks based on a specified topology

//TO DO:
// - some net input files are highly irregular, and fail (address this somehow; perhaps "node" and "atom" cannot be handled in the ways specified in the cgd input routine?); example fails: oft, xaa

/* includes etc. */
#include <vector>
#include <string>
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <math.h>
#include <float.h>
#include <queue>
#include <sys/time.h>
using namespace std;
#include "geometry.h"
#include "mindist.h"
#include "networkinfo.h"
#include "rmsd.h"
#include "networkstorage.h"
#include "net.h"
#include "symmetry.h"
#include "networkio.h"
#include "zeo_consts.h"

//global, const, struct
string global_net; 
const bool write_old_ml_file = false; //just for debugging, writes out ratio of different building blocks in the structure in an old format
const int MIN_DIM = 2;
typedef struct {int seconds; int microseconds;} seconds_microseconds;

//should make these variables
bool allow_poor_quality_framework = false;
bool bypass_collision_detection = false;

//methods
void construct_framework_from_alignment(vector< pair<int,int> > this_alignment, int index_of_this_framework, vector<int> *count_each_molecule, double *previous_lowest_deviation, ATOM_NETWORK *framework, int *assembly_method, int num_full_vertices, ATOM_NETWORK full_cell, vector< vector<int> > mol_index_at_each_basic_vertex, vector< vector< vector<MOLECULE> > > basic_oriented_molecules, string prefix, string net_name, vector<CONNECTION> two_way_connections, double unit_edge_length, double layer_separation_A, bool layer_override, int dimensionality, bool a_periodicity, bool b_periodicity, bool c_periodicity, bool disconnected_net);
bool iterate_alignment_code(vector<int> *example_alignment_code, vector< vector< pair<int,int> > > *alignment_options_each_basic_vertex);
seconds_microseconds get_time_difference(struct timeval start, struct timeval end);

/*------------------------------
        MAIN BEGINS HERE
------------------------------*/

int main(int argc, char *argv[]) {

  bool verbose = false; //set to true for debugging purposes, to write out individually aligned building block files, write extra info to terminal, etc.

  int num_args = 4;
  if(argc<num_args) {
    printf("%d arguments were provided but at least %d are required:\n", argc, num_args);
    printf("./this_code\n");
    printf("net.cgd (a topology definition file - e.g. from the RCSR http://rcsr.anu.edu.au/)\n");
    printf("augment_2c (0 or 1, representing whether to insert 2-c vertices along each edge, e.g. for pcu net MOF-5)\n");
    printf("output_file_prefix\n");
    printf("OPTIONAL: at least one mol.xyz, if you want to build a framework - otherwise only the net will be returned\n");
    printf("OPTIONAL: for 2D (layer) structures, a separation distance in Angstroms, i.e., the 3rd unit cell side length\n");
    printf("Please try again.\n");
    exit(EXIT_FAILURE);
  }

  //start the clock
  struct timeval time_at_start, time_at_end;
  gettimeofday(&time_at_start, NULL);

  //set up
  bool poor_quality_framework = false;

  //by default, we want to populate the table with radii and masses, even if later we will not reference them
  initializeRadTable();
  initializeMassTable();
  initializeCovRadTable();

//let's have a new line to begin the terminal output
printf("\n");

  //parse arguments
  string net_name(argv[1]);
  bool augment_2c = !(atoi(argv[2])==0);
  string prefix(argv[3]);
  double layer_separation_A = -1;
  bool layer_override = false;
  if(argc>num_args) { //attempt to read layer separation override
    double l = atof(argv[argc-1]);
    if(l<TOLERANCE) { //i.e. if too small to be valid, or if string was not a number (which gives 0 as result)
      //do nothing
    } else {
      layer_separation_A = l;
      layer_override = true;
    }
  }
  if(verbose) {
    if(layer_override) printf("DEBUG: layer separation override of %.3f was requested\n", layer_separation_A);
    else printf("DEBUG: no layer separation override was provided\n");
  }
  vector<string> mol_names;
  for(int i=num_args; i<argc; i++) {
    if(layer_override && i==argc-1) {
      //do nothing - we have determined that the final arg is a layer override and not a file name
    } else {
      string m(argv[i]);
      mol_names.push_back(m);
    }
  }
  
/*------------------------------
         NET HANDLING
------------------------------*/

  //make cell
  FILE *cgd = fopen(net_name.c_str(),"r");
  if(cgd==NULL) {
    printf("ERROR: could not open cgd file with name %s\n", net_name.c_str());
    exit(EXIT_FAILURE);
  }
  ATOM_NETWORK basic_cell;
  bool edges_provided_exhaustively = read_cgd(cgd, &basic_cell, &global_net); //we need to know whether the edges were provided exhaustively in the input file or not, so that we can reconstruct the missing edges if required!
  if(!edges_provided_exhaustively) add_missing_edges(&basic_cell);
  //if augmenting with 2c vertices, add the vertices and normal edges here - the dummy edges will be handled later - but only do this if the 2c vertices were not already provided!
  bool augmented_2c_already = false;
  int num_basic_vertices = basic_cell.vertices.size();
  for(int i=0; i<num_basic_vertices && !augmented_2c_already; i++) {
    if(basic_cell.vertices.at(i).edges.size()==2) augmented_2c_already = true; //are there already 2c vertices? if so, they were provided in the input net file, and so we don't need to augment
  }
  if(augment_2c) {
    if(!augmented_2c_already) {
      add_2c_vertices_and_normal_edges(&basic_cell);
      printf("augmented net with additional 2-c vertices as requested\n");
    } else printf("did not need to augment net with additional 2-c vertices as requested, since net file explicitly provided them\n");
  }
  fclose(cgd);

  //how much cgd data was read
  num_basic_vertices = basic_cell.vertices.size(); //update this one incase we added 2c vertices
  int num_basic_edges = 0;
  for(int i=0; i<num_basic_vertices; i++) {
    num_basic_edges+=basic_cell.vertices.at(i).edges.size();
  }

  //what is the unit edge length? we need this to scale the final unit cell
  double unit_edge_length = get_unit_edge_length(&basic_cell);

  //for each basic vertex, find its symmetric positions and push them to a new cell
  ATOM_NETWORK full_cell = basic_cell;
  full_cell.vertex_basic_indices.clear();
  full_cell.vertex_symmetry_operators.clear();
  full_cell.vertices.clear();
  for(int i=0; i<num_basic_vertices; i++) {
    XYZ test = basic_cell.vertices.at(i).abc;
    vector<XYZ> equiv = GetEquivalentPositions(basic_cell.sym_ID,&test);
    int num_equiv = equiv.size();
    for(int j=0; j<num_equiv; j++) {
      XYZ here = equiv.at(j);
      //for this symmetry position, we want to know if it is unique before we commit it
      int num_v = full_cell.vertices.size();
      bool unique = true;
      for(int k=0; k<num_v && unique; k++) {
        unique = IsUniqueVertex(&here, full_cell);
      }
      if(unique) {
        VERTEX v(here.x, here.y, here.z);
        v.expected_edges = basic_cell.vertices.at(i).expected_edges;
        //this vertex is unique - that means that we also want to get its edges by symmetry!
        int num_e = basic_cell.vertices.at(i).edges.size();
        for(int k=0; k<num_e; k++) {
          XYZ test_e = basic_cell.vertices.at(i).edges.at(k);
          vector<XYZ> equiv_e = GetEquivalentPositions(basic_cell.sym_ID,&test_e);
          v.edges.push_back(equiv_e.at(j)); //this is nice - we know that we are working on symmetric vertex image with index j - so we can just keep the corresponding edge image - shifting it to unit cell by the same amount that we used for the vertex!
        }
        int num_de = basic_cell.vertices.at(i).dummy_edges.size();
        for(int k=0; k<num_de; k++) {
          XYZ test_de = basic_cell.vertices.at(i).dummy_edges.at(k);
          vector<XYZ> equiv_de = GetEquivalentPositions(basic_cell.sym_ID,&test_de);
          v.dummy_edges.push_back(equiv_de.at(j)); //this is nice - we know that we are working on symmetric vertex image with index j - so we can just keep the corresponding edge image - shifting it to unit cell by the same amount that we used for the vertex!
        }
        full_cell.vertices.push_back(v);
        full_cell.vertex_symmetry_operators.push_back(j);
        full_cell.vertex_basic_indices.push_back(i);
      }
    }
  }
  int num_full_vertices = full_cell.vertices.size();
  int num_full_edges = 0;
  if(verbose) printf("by symmetry there are %d non-overlapping vertex positions; in fractional coords:\n", num_full_vertices);
  vector<bool> basic_indices_still_basic; //keep track of which basic vertices are basic, and which were in fact equivalent by symmetry (should only happen for augmented 2c vertices, but theoretically could also happen for a sloppy topology file)
  for(int i=0; i<num_basic_vertices; i++) basic_indices_still_basic.push_back(false);
  for(int i=0; i<num_full_vertices; i++) {
    XYZ here = full_cell.vertices.at(i).abc;
    int basic_ID = full_cell.vertex_basic_indices.at(i);
    if(verbose) printf("\t%.3f %.3f %.3f (basic index %d following symmetry operation %d)\n", here.x, here.y, here.z, basic_ID, full_cell.vertex_symmetry_operators.at(i));
    basic_indices_still_basic.at(basic_ID) = true;
    num_full_edges+=full_cell.vertices.at(i).edges.size();
  }

  //now we can remove 'basic' vertices that we now know are redundant by symmetry (probably only happens for 2-c augmented vertices, but run it anyway in the general case, just in case
  vector<VERTEX> temp_vertices;
  vector<int> new_IDs_for_old_basic_vertex_IDs;
  for(int i=0; i<num_basic_vertices; i++) {
    if(basic_indices_still_basic.at(i)) {
      new_IDs_for_old_basic_vertex_IDs.push_back(temp_vertices.size());
      temp_vertices.push_back(basic_cell.vertices.at(i));
    } else new_IDs_for_old_basic_vertex_IDs.push_back(-1); //this means the basic vertex was deleted, and we should not refer to it in future - if we do, we would be referring to a different data point, so putting -1 here makes sure that doesn't happen (we would get seg fault instead of looking at the wrong data point (which might have been seg fault anyway))
  }
  //now update the stored data
  basic_cell.vertices = temp_vertices;
  num_basic_vertices = basic_cell.vertices.size();
  for(int i=0; i<num_full_vertices; i++) {
    int old_basic_ID = full_cell.vertex_basic_indices.at(i);
    full_cell.vertex_basic_indices.at(i) = new_IDs_for_old_basic_vertex_IDs.at(old_basic_ID);
  }
  
  //identify the two-way connections in the full_cell
  vector<CONNECTION> two_way_connections;
  bool disconnected_net = find_two_way_connections(&full_cell, &two_way_connections);
  if(disconnected_net) printf("NOTICE: the net file comprises disconnected nets - the connection-based assembly method can not be performed\n");
  int num_two_way_connections = two_way_connections.size();

  //identify whether this is a 2D (layer) structure
  vector<bool> edge_traversed; for(int i=0; i<num_two_way_connections; i++) edge_traversed.push_back(false);
  vector<int> vertex_visit_a, vertex_visit_b, vertex_visit_c; vector<bool> vertex_visited; for(int i=0; i<num_full_vertices; i++) {vertex_visit_a.push_back(0); vertex_visit_b.push_back(0); vertex_visit_c.push_back(0); vertex_visited.push_back(false);}
  vector<bool> is_periodic; for(int i=0; i<3; i++) is_periodic.push_back(false);
  int dimensionality = recursive_test_dimensionality(0, 0, 0, 0, &two_way_connections, &edge_traversed, &vertex_visit_a, &vertex_visit_b, &vertex_visit_c, &vertex_visited, &is_periodic);
  bool a_periodicity = false, b_periodicity = false, c_periodicity = false;
  if(is_periodic.at(0)) a_periodicity = true;
  if(is_periodic.at(1)) b_periodicity = true;
  if(is_periodic.at(2)) c_periodicity = true;
  printf("from net information, detected that this is a %d-dimensional net\n", dimensionality);
  if(dimensionality<MIN_DIM) {
    printf("ERROR: net dimensionality of %d is less than the minimum, %d - cannot handle this kind of net\n", dimensionality, MIN_DIM);
    exit(EXIT_FAILURE);
  }

  //use two-way connections to augment both the basic and full nets, if necessary - but remember, not if we read a special net file which already provides this data
  if(augment_2c && !augmented_2c_already) add_2c_dummy_edges(&basic_cell, &full_cell, &two_way_connections);

  //now move everything to the unit cell
  for(int i=0; i<num_basic_vertices; i++) {
    XYZ here = basic_cell.vertices.at(i).abc;
    XYZ here_uc = trans_to_origuc(here);
    XYZ shift = here_uc - here;
    //we know the shift, so apply it to vertex, edge and dummy edge positions
    if(verbose) printf("DEBUG: basic vertex %d was at %.3f %.3f %.3f\n", i, basic_cell.vertices.at(i).abc.x, basic_cell.vertices.at(i).abc.y, basic_cell.vertices.at(i).abc.z);
    basic_cell.vertices.at(i).abc = here_uc;
    if(verbose) printf("DEBUG: basic vertex %d now at %.3f %.3f %.3f\n", i, basic_cell.vertices.at(i).abc.x, basic_cell.vertices.at(i).abc.y, basic_cell.vertices.at(i).abc.z);
    int num_e = basic_cell.vertices.at(i).edges.size();
    for(int j=0; j<num_e; j++) {
      XYZ here_e = basic_cell.vertices.at(i).edges.at(j);
      if(verbose) printf("DEBUG: basic vertex %d edge %d was at %.3f %.3f %.3f\n", i, j, basic_cell.vertices.at(i).edges.at(j).x, basic_cell.vertices.at(i).edges.at(j).y, basic_cell.vertices.at(i).edges.at(j).z);
      basic_cell.vertices.at(i).edges.at(j) = here_e+shift;
      if(verbose) printf("DEBUG: basic vertex %d edge %d now at %.3f %.3f %.3f\n", i, j, basic_cell.vertices.at(i).edges.at(j).x, basic_cell.vertices.at(i).edges.at(j).y, basic_cell.vertices.at(i).edges.at(j).z);
    }
    int num_dummy_e = basic_cell.vertices.at(i).dummy_edges.size();
//printf("%d dummy edges\n", num_dummy_e);
    for(int j=0; j<num_dummy_e; j++) {
      XYZ here_e = basic_cell.vertices.at(i).dummy_edges.at(j);
      if(verbose) printf("DEBUG: basic vertex %d dummy_edge %d was at %.3f %.3f %.3f\n", i, j, basic_cell.vertices.at(i).dummy_edges.at(j).x, basic_cell.vertices.at(i).dummy_edges.at(j).y, basic_cell.vertices.at(i).dummy_edges.at(j).z);
      basic_cell.vertices.at(i).dummy_edges.at(j) = here_e+shift;
      if(verbose) printf("DEBUG: basic vertex %d dummy_edge %d now at %.3f %.3f %.3f\n", i, j, basic_cell.vertices.at(i).dummy_edges.at(j).x, basic_cell.vertices.at(i).dummy_edges.at(j).y, basic_cell.vertices.at(i).dummy_edges.at(j).z);
    }
  }
  for(int i=0; i<num_full_vertices; i++) {
    XYZ here = full_cell.vertices.at(i).abc;
    XYZ here_uc = trans_to_origuc(here);
    XYZ shift = here_uc - here;
    //we know the shift, so apply it to vertex and edge positions
    full_cell.vertices.at(i).abc = here_uc;
    int num_e = full_cell.vertices.at(i).edges.size();
    for(int j=0; j<num_e; j++) {
      XYZ here_e = full_cell.vertices.at(i).edges.at(j);
      full_cell.vertices.at(i).edges.at(j) = here_e+shift;
    }
    int num_dummy_e = full_cell.vertices.at(i).dummy_edges.size();
//printf("%d full dummy edges\n", num_dummy_e);
    for(int j=0; j<num_dummy_e; j++) {
      XYZ here_e = full_cell.vertices.at(i).dummy_edges.at(j);
      full_cell.vertices.at(i).dummy_edges.at(j) = here_e+shift;
    }
  }

  //replace "V" with a number in the xyz output files? this will color-code vertices based on their connectivity
  bool rename_vertices_by_connectivity = true;  

  //write unit cell to vtk file
  string basic_uc_name = prefix+"_net_unit_cell.vtk";
  FILE *basic_uc = fopen(basic_uc_name.c_str(),"w");
  if(basic_uc==NULL) {
    printf("ERROR: could not open output unit cell file with name %s\n", basic_uc_name.c_str());
    exit(EXIT_FAILURE);
  }
  write_unit_cell(basic_uc, &basic_cell, net_name, true); //is_net? true
  fclose(basic_uc);
  printf("\tunit cell written to %s\n", basic_uc_name.c_str());

  //write full edges to vtk file
  string full_e_name = prefix+"_net_full_edges.vtk";
  FILE *full_e = fopen(full_e_name.c_str(),"w");
  if(full_e==NULL) {
    printf("ERROR: could not open output net edge file with name %s\n", full_e_name.c_str());
    exit(EXIT_FAILURE);
  }
  write_edges(full_e, &full_cell, net_name);
  fclose(full_e);
  printf("\t%d edges (after symmetry operations) written to %s\n", num_full_edges, full_e_name.c_str());

  //write full vertices to xyz file
  string full_v_name = prefix+"_net_full_vertices.xyz";
  FILE *full_v = fopen(full_v_name.c_str(),"w");
  if(full_v==NULL) {
    printf("ERROR: could not open output net vertex file with name %s\n", full_v_name.c_str());
    exit(EXIT_FAILURE);
  }
  write_vertices(full_v, &full_cell, net_name, rename_vertices_by_connectivity, true); //is_net? true
  fclose(full_v);
  printf("\t%d vertices (after symmetry operations) written to %s\n", num_full_vertices, full_v_name.c_str());

  if(verbose) { //extra output files for debugging etc.
    //write basic edges to vtk file
    string basic_e_name = prefix+"_net_basic_edges.vtk";
    FILE *basic_e = fopen(basic_e_name.c_str(),"w");
    if(basic_e==NULL) {
      printf("ERROR: could not open output net edge file with name %s\n", basic_e_name.c_str());
      exit(EXIT_FAILURE);
    }
    write_edges(basic_e, &basic_cell, net_name);
    fclose(basic_e);
    printf("\t%d edges (before symmetry operations) written to %s\n", num_basic_edges, basic_e_name.c_str());

    //write basic vertices to xyz file
    string basic_v_name = prefix+"_net_basic_vertices.xyz";
    FILE *basic_v = fopen(basic_v_name.c_str(),"w");
    if(basic_v==NULL) {
      printf("ERROR: could not open output net vertex file with name %s\n", basic_v_name.c_str());
      exit(EXIT_FAILURE);
    }
    write_vertices(basic_v, &basic_cell, net_name, rename_vertices_by_connectivity, true); //is_net? true
    fclose(basic_v);
    printf("\t%d vertices (before symmetry operations) written to %s\n", num_basic_vertices, basic_v_name.c_str());

    //write full vertices and edges to cif (edges are shortened, so as to present a 'valid' abstract cif containing non-overlapping dummy atoms)
    string abstract_cif_name = prefix+"_net_abstract.cif";
    FILE *abstract_cif = fopen(abstract_cif_name.c_str(),"w");
    if(abstract_cif==NULL) {
      printf("ERROR: could not open output abstract cif file with name %s\n", abstract_cif_name.c_str());
      exit(EXIT_FAILURE);
    }
    write_abstract_cif(abstract_cif, &full_cell, net_name);
    fclose(abstract_cif);
    printf("\tabstract cif written to %s\n", abstract_cif_name.c_str());
  }

/*------------------------------
       FRAMEWORK HANDLING (optional)
------------------------------*/

  int num_mols = mol_names.size();
  if(num_mols>0) {
    //handle multiple molecule input
    printf("about to parse multiple molecules from command line\n");
    vector<MOLECULE> original_molecules;
    vector<int> count_each_molecule;
    int max_num_sites = 12;
    vector< vector<int> > indices_of_n_site_molecules;
    for(int i=0; i<=max_num_sites; i++) {
      vector<int> v;
      indices_of_n_site_molecules.push_back(v);
    }
    for(int i=0; i<num_mols; i++) {
      //parse molecule, getting sites
      FILE *mol_file = fopen(mol_names.at(i).c_str(),"r");
      if(mol_file==NULL) {
        printf("ERROR: could not open molecule file with name %s\n", mol_names.at(i).c_str());
        exit(EXIT_FAILURE);
      }
      MOLECULE mol;
      read_xyz(mol_file, &mol, mol_names.at(i).c_str());
      fclose(mol_file);
      vector<int> sites, dummy_sites;
      int num_dummy_sites = 0; //if a 2-c molecule contains a 'J' atom, this is used to specify the position of the dummy site, i.e. the 'width' of the molecule which is aligned along the dummy edge
      for(int j=0; j<mol.atoms_type.size(); j++) {
        if(mol.atoms_type.at(j)[0]=='Q') sites.push_back(j); //allow anything beginning with 'Q' to be a site, so that distinct sites can be handled (matched to each other, e.g. Q could refer to a carboxy connection, while Q1 refers to a nitrogen connection)
        else if(mol.atoms_type.at(j)[0]=='J') { //anything beginning with 'J' is a dummy site
          dummy_sites.push_back(j);
          num_dummy_sites++;
        }
      }
      int num_sites = sites.size();
      printf("\tread molecule ID %d from file %s containing %d sites (and %d dummy sites)\n", i, mol_names.at(i).c_str(), num_sites, num_dummy_sites);
      if(num_sites<=0) {
        printf("ERROR: molecule read from %s contains no connection sites\nNOTE: currently, sites are handled as follows: any element symbol beginning with 'Q' represents a site, and sites in adjacent molecules are then attempted to be perfectly overlaid; 'J' is uesd in 2-c molecules if you wish to specify the third, dummy direction used to align the 2-c molecule to the net\n", mol_names.at(i).c_str());
        exit(EXIT_FAILURE);
      } else if(num_sites!=2 && num_dummy_sites!=0) {
        printf("ERROR: %d dummy sites were read, but these are only permitted if number of sites is 2 (%d sites were read)\n", num_dummy_sites, num_sites);
        exit(EXIT_FAILURE);
      } else if(num_dummy_sites!=0 && num_dummy_sites!=1) {
        printf("ERROR: %d dummy sites were read, but only 0 or 1 are permitted\n", num_dummy_sites);
        exit(EXIT_FAILURE);
      } else if(num_sites>max_num_sites) {
        printf("ERROR: %d sites were read, but maximum number of sites is %d\n", num_sites, max_num_sites);
        exit(EXIT_FAILURE);
      }
      indices_of_n_site_molecules.at(num_sites).push_back(i);
      mol.sites = sites;
      //get centre of mass
      get_com(&mol);
      //if the molecule is 2-c, we need at least 3 sites to do rmsd alignment - insert an orthogonal dummy site (if it was not explicitly provided)!
      if(num_sites==2 && num_dummy_sites==0) {
        dummy_sites.push_back(assign_dummy_site(&mol));
      }
      mol.dummy_sites = dummy_sites;
      //finished with this molecule now
      original_molecules.push_back(mol);
      count_each_molecule.push_back(0); //this keeps track of how many of this molecule we used - so far none, until we do alignment ...
    }

    //handle bipartite and higher order segmented frameworks - if there are more molecules with a particular connectivity than there are distinct basic vertices, we will need to split up the basic vertices
    vector< vector<int> > indices_of_n_site_basic_vertices, indices_of_n_site_full_vertices;
    for(int i=0; i<=max_num_sites; i++) {
      vector<int> v;
      indices_of_n_site_basic_vertices.push_back(v);
      indices_of_n_site_full_vertices.push_back(v);
    }
    for(int i=0; i<num_basic_vertices; i++) {
      int num_sites = basic_cell.vertices.at(i).edges.size();
      if(num_sites<=0) {
        printf("ERROR: basic cell's vertex read from %s contains no connection sites (edges)\n", net_name.c_str());
        exit(EXIT_FAILURE);
      }
      if(num_sites>max_num_sites) {
        printf("ERROR: maximum number of sites is %d\n", max_num_sites);
        exit(EXIT_FAILURE);
      }
      indices_of_n_site_basic_vertices.at(num_sites).push_back(i);
    }
    for(int i=0; i<num_full_vertices; i++) {
      int num_sites = full_cell.vertices.at(i).edges.size();
      if(num_sites<=0) {
        printf("ERROR: full cell's vertex read from %s contains no connection sites (edges)\n", net_name.c_str());
        exit(EXIT_FAILURE);
      }
      if(num_sites>max_num_sites) {
        printf("ERROR: maximum number of sites is %d\n", max_num_sites);
        exit(EXIT_FAILURE);
      }
      indices_of_n_site_full_vertices.at(num_sites).push_back(i);
    }
    for(int i=0; i<max_num_sites; i++) {
      int num_mols_here = indices_of_n_site_molecules.at(i).size();
      int num_basic_vertices_here = indices_of_n_site_basic_vertices.at(i).size();
      int num_full_vertices_here = indices_of_n_site_full_vertices.at(i).size();
      if(num_mols_here>0 || num_basic_vertices_here>0) {
if(verbose) {
  printf("DEBUG: there are %d molecules and %d basic vertices with %d sites\n", num_mols_here, num_basic_vertices_here, i);
  printf("\tmolecule IDs:");
  for(int j=0; j<num_mols_here; j++) {
    printf(" %d", indices_of_n_site_molecules.at(i).at(j));
  }
  printf("\n\tbasic vertex IDs:");
  for(int j=0; j<num_basic_vertices_here; j++) {
    printf(" %d", indices_of_n_site_basic_vertices.at(i).at(j));
  }
  printf("\n");
}
        if(num_mols_here<=0) {
          printf("ERROR: there are no molecules with %d sites to fit to basic vertices\n", i);
          exit(EXIT_FAILURE);
        } else if(num_basic_vertices_here<=0) {
          printf("ERROR: there are no basic vertices with %d sites for molecules to fit to\n", i);
          exit(EXIT_FAILURE);
        } else if(num_basic_vertices_here<num_mols_here) {
          printf("there are fewer vertices in the net that are unique by symmetry (%d) than molecules (%d) with %d sites - a bipartite (or higher order) segmentation of the net is necessary\n", num_basic_vertices_here, num_mols_here, i);
          //which kinds of net segmentation can we deal with, or are valid?
          if(num_full_vertices_here<num_mols_here) { //no matter what we do, there are just not enough full vertices for the number of molecules provided
            printf("ERROR: too many molecules (%d) were provided to fit to the %d vertices (after symmetry operations) with %d sites; there is no segmentation of the net that will allow a framework to be constructed\nNOTE: it may be possible to construct this model by providing a supercell of the net instead\n", num_mols_here, num_full_vertices_here, i);
            exit(EXIT_FAILURE);          
          } else if(num_basic_vertices_here!=1) { //e.g., if there are two n-c vertices and three n-c molecules, it is not clear how to do a segmentation
            printf("ERROR: currently net segmentation is not handled for anything other than 1 (%d) symmetrically unique vertices - this may be introduced in a later release\n", num_basic_vertices_here);
            exit(EXIT_FAILURE);
          } else if(num_mols_here!=2) { //else, we know we have just one basic vertex - but if we have to assign anything other than two molecules to it, it is not clear how to do it - but we can potentially do it recursively
            printf("ERROR: currently net segmentation is not handled for anything other than 2 (%d) molecules - this may be introduced in a later release\n", num_mols_here);
            exit(EXIT_FAILURE);
          } else { //else, we know that we have enough full vertices, and exactly one basic vertex and two molecules
            //add basic vertices to the net until we have one per full vertex
            int basic_vertex_ID = indices_of_n_site_basic_vertices.at(i).at(0); //copy a vertex - safe to take ID 0 since there is exactly one entry
            while(num_basic_vertices_here<num_full_vertices_here) {
              VERTEX copy_vertex = basic_cell.vertices.at(basic_vertex_ID);
              int copy_vertex_ID = basic_cell.vertices.size();
              //push to basic cell, and update local variables to reference this new vertex
              basic_cell.vertices.push_back(copy_vertex);
              num_basic_vertices_here++;
              num_basic_vertices++;
              indices_of_n_site_basic_vertices.at(i).push_back(copy_vertex_ID);
            }
            //now, these basic vertices are not yet assigned to any full vertices in the framework - need to update the full vertices to reference them
if(verbose) printf("DEBUG: after adding copies of symmetrically unique vertices, there are now %d to which the %d molecules can be aligned; need to reference these symmetrically unique vertices in the %d vertices after symmetry which have %d sites\n", num_basic_vertices_here, num_mols_here, num_full_vertices_here, i);
            for(int j=0; j<num_full_vertices_here; j++) {
              int this_full_vertex_ID = indices_of_n_site_full_vertices.at(i).at(j);
if(verbose) printf("\tvertex %d is originally based on symmetrically unique vertex %d and symmetry operation %d\n", this_full_vertex_ID, full_cell.vertex_basic_indices.at(this_full_vertex_ID), full_cell.vertex_symmetry_operators.at(this_full_vertex_ID));
              int new_basic_vertex_ID = indices_of_n_site_basic_vertices.at(i).at(j);
              full_cell.vertex_basic_indices.at(this_full_vertex_ID) = new_basic_vertex_ID;
if(verbose) printf("\t\tnow based on symmetrically unique vertex %d and symmetry operation %d\n", full_cell.vertex_basic_indices.at(this_full_vertex_ID), full_cell.vertex_symmetry_operators.at(this_full_vertex_ID));
            }
          }
        }
      }
    }

/*------------------------------
    BUILDING BLOCK ALIGNMENT (part of optional framework handling)
------------------------------*/

    if(disconnected_net) printf("NOTICE: the net file comprises disconnected nets - the connection-based assembly method can not be performed\n");
    //we want to get fitted molecules for each of the basic vertices
if(verbose) printf("DEBUG: attempting to create oriented molecules for each of %d basic vertices\n", num_basic_vertices);
//    vector< vector<MOLECULE> > basic_oriented_molecules;
    vector< vector< vector<MOLECULE> > > basic_oriented_molecules;
    vector< vector<int> > mol_index_at_each_basic_vertex;
    printf("there are %d symmetrically unique vertices, which will be fit to molecules as follows:\n", num_basic_vertices);
    for(int i=0; i<num_basic_vertices; i++) {
//      vector<MOLECULE> basic_oriented_molecules_here;
      vector< vector<MOLECULE> > basic_oriented_molecules_here;
      vector<int> mol_index_at_each_basic_vertex_here;
      //find which molecule to orient to this vertex
      int num_sites = basic_cell.vertices.at(i).edges.size();
      int num_mols_here = indices_of_n_site_molecules.at(num_sites).size();
if(verbose) printf("DEBUG: symmetrically unique vertex ID %d has %d sites - there are %d molecules satisfying this number of sites\n", i, num_sites, num_mols_here);
      if(num_mols_here==-1) {
        printf("ERROR: could not find a molecule with %d sites to fit to this vertex\n", num_sites);
        exit(EXIT_FAILURE);
      }
      for(int j=0; j<num_mols_here; j++) {
        int mol_index = indices_of_n_site_molecules.at(num_sites).at(j);
        vector<MOLECULE> all_fits = get_multiple_best_RMSD_fits(original_molecules.at(indices_of_n_site_molecules.at(num_sites).at(j)), &basic_cell, i, global_net, prefix, mol_index); //get all the permissible fits of our molecule with this vertex
        printf("\tsymmetrically unique vertex %d has %d possible fits with molecule ID %d\n", i, (int)(all_fits.size()), indices_of_n_site_molecules.at(num_sites).at(j));
        basic_oriented_molecules_here.push_back(all_fits);
        mol_index_at_each_basic_vertex_here.push_back(mol_index);
      }
      basic_oriented_molecules.push_back(basic_oriented_molecules_here);
      mol_index_at_each_basic_vertex.push_back(mol_index_at_each_basic_vertex_here);
    }

//WE HAVE TO DO IT THIS WAY - ENUMERATE ALL THE POSSIBILITIES, AND THEN RUN EACH OF THEM
    //scheme for building all possible frameworks based on distinct molecular alignments begins here
    vector< vector< pair<int,int> > > alignment_options_each_basic_vertex;
    for(int bv=0; bv<num_basic_vertices; bv++) {
if(verbose) printf("DEBUG: basic vertex %d has these alignment options:\n", bv);
      vector< pair<int,int> > alignment_options_this_basic_vertex; //an 'alignment option' is a pair consisting of the molecule index and the alignment index of that molecule
      int num_molecule_choices_this_basic_vertex = basic_oriented_molecules.at(bv).size();
      for(int mc=0; mc<num_molecule_choices_this_basic_vertex; mc++) {
        int num_alignment_choices_this_molecule = basic_oriented_molecules.at(bv).at(mc).size();
        for(int ac=0; ac<num_alignment_choices_this_molecule; ac++) {
          pair<int,int> p;
          p.first = mc;
          p.second = ac;
if(verbose) printf("\tmolecule choice %d alignment choice %d\n", mc, ac);
          alignment_options_this_basic_vertex.push_back(p);
        }
      }
if(verbose) printf("DEBUG: i.e. basic vertex %d has %d alignment options\n", bv, (int)(alignment_options_this_basic_vertex.size()));
if(alignment_options_this_basic_vertex.size()<=0) {
  printf("ERROR: basic vertex %d has %d alignment options; must be at least one\n", bv, (int)(alignment_options_this_basic_vertex.size()));
  exit(EXIT_FAILURE);
}
      alignment_options_each_basic_vertex.push_back(alignment_options_this_basic_vertex);
    }

    //at this point, we should have a sensible list of the candidate frameworks to build - let's iterate over them and find the best
    double previous_lowest_deviation = -1;
    ATOM_NETWORK framework;
    int count = 0;
    int assembly_method = -1;
    vector<int> example_alignment_code;
    for(int bv=0; bv<num_basic_vertices; bv++) example_alignment_code.push_back(0);
    bool continuing = true;
    while(continuing) {
      vector< pair<int,int> > this_alignment;
      for(int bv=0; bv<num_basic_vertices; bv++) { //fill in the example_alignment based on the example_alignment_code
        this_alignment.push_back( alignment_options_each_basic_vertex.at(bv).at(example_alignment_code.at(bv)) );
      }
      //does this alignment combination use all the provided molecules?
      bool all_mols_used = true;
      for(int m=0; m<num_mols && all_mols_used; m++) {
        bool this_mol_used = false;
        for(int bv=0; bv<num_basic_vertices && !this_mol_used; bv++) {
          if(m==indices_of_n_site_molecules.at(basic_cell.vertices.at(bv).edges.size()).at(this_alignment.at(bv).first)) this_mol_used = true;
        }
        all_mols_used = this_mol_used;
      }
      if(all_mols_used) { //this alignment uses all molecules - but are the connection sites correct (e.g., do not want to connect a Q1 to a Q2, representing say, a carboxy connection site and a nitrogen-terminated linker)
        //here we will check whether constructing with this alignment is possible (based on the different kinds of connection sites that may be specified)
        bool possible = true;
        for(int i=0; i<num_two_way_connections && possible; i++) {
          CONNECTION c = two_way_connections.at(i);
          int basic_ID_m1 = full_cell.vertex_basic_indices.at(c.v1);
          MOLECULE m1 = basic_oriented_molecules.at(basic_ID_m1).at(this_alignment.at(basic_ID_m1).first).at(this_alignment.at(basic_ID_m1).second);
          int basic_ID_m2 = full_cell.vertex_basic_indices.at(c.v2);
          MOLECULE m2 = basic_oriented_molecules.at(basic_ID_m2).at(this_alignment.at(basic_ID_m2).first).at(this_alignment.at(basic_ID_m2).second);
          int m1_site_index = m1.permutation.at(c.e1);
          int m2_site_index = m2.permutation.at(c.e2);
          int m1_atom_index = m1.sites.at(m1_site_index);
          int m2_atom_index = m2.sites.at(m2_site_index);
          if(m1.atoms_type.at(m1_atom_index)!=m2.atoms_type.at(m2_atom_index)) possible=false;
        }
if(verbose) {
  if(possible) printf("NOTICE (decider): this alignment should be possible\n");
  else printf("NOTICE (decider): this alignment should NOT be possible\n");
}

        if(possible) {
          //VALID ALIGNMENT: DO THIS ALIGNMENT
  if(verbose) {
    printf("DEBUG: working on alignment combination %d...\n", count);
    printf("DEBUG: \there we will perform the alignment given by:\n");
    for(int bv=0; bv<num_basic_vertices; bv++) {
      printf("\t\tbv %d = molecule choice %d alignment choice %d\n", bv, this_alignment.at(bv).first, this_alignment.at(bv).second);
    }
  }
          construct_framework_from_alignment(this_alignment, count, &count_each_molecule, &previous_lowest_deviation, &framework, &assembly_method, num_full_vertices, full_cell, mol_index_at_each_basic_vertex, basic_oriented_molecules, prefix, net_name, two_way_connections, unit_edge_length, layer_separation_A, layer_override, dimensionality, a_periodicity, b_periodicity, c_periodicity, disconnected_net);
          count++;
        }
      }
      //now, iterate
      continuing = iterate_alignment_code(&example_alignment_code, &alignment_options_each_basic_vertex);
    }
    if(count==0) printf("it was not possible to construct any framework models which utilized all provided molecules\n");
    else printf("in total, %d molecule to basic vertex alignment combinations were used to construct frameworks\n", count);

//let's have a new line to separate the details from the results
printf("\n");

    //if we achieved a valid framework, write it out
    if(previous_lowest_deviation>SITE_DEVIATION_WARNING) poor_quality_framework = true;
//    if(framework.vertices.size()>0) {
    if(framework.numAtoms>0 && (!poor_quality_framework || allow_poor_quality_framework)) {
      printf("the best framework exhibits a maximum site deviation of %.3f, using the", previous_lowest_deviation);
      if(assembly_method==0) printf(" net-"); else printf(" connection-");
      printf("based assembly method\n");

      //write out number of each component 
      if(write_old_ml_file) {
        FILE *ml_ratio_verbose = fopen((net_name+".ml").c_str(), "w"); //for debugging etc.
        if(ml_ratio_verbose==NULL) {
          printf("ERROR: could not open output ratio file with name %s\n", (net_name+".ml").c_str());
          exit(EXIT_FAILURE);
        }
        for(int i=0; i<num_mols; i++) {
          fprintf(ml_ratio_verbose, "%d\n", count_each_molecule.at(i));
        }
        fclose(ml_ratio_verbose);
      }
      string ratio_file_name = prefix+"_ratio.txt";
      FILE *ratio = fopen(ratio_file_name.c_str(), "w"); //nicer version for general use
      if(ratio==NULL) {
        printf("ERROR: could not open output ratio file with name %s\n", ratio_file_name.c_str());
        exit(EXIT_FAILURE);
      }
      for(int i=0; i<num_mols; i++) {
        printf("cell contains %d instances of molecule ID %d (%s)\n", count_each_molecule.at(i), i, mol_names.at(i).c_str());
        fprintf(ratio, "%d %s\n", count_each_molecule.at(i), mol_names.at(i).c_str());
      }
      fclose(ratio);
      printf("\tratio of building blocks written to %s\n", ratio_file_name.c_str());

/*
      //replace "V" with a number in the xyz output files? this will color-code vertices based on their connectivity
      rename_vertices_by_connectivity = false; //we don't want this for molecular output
      //print cssr ...
      string frame_file_name = prefix+"_framework.cssr"; //write framework molecule to cssr file
      FILE *frame_file = fopen(frame_file_name.c_str(),"w");
      if(frame_file==NULL) {
        printf("ERROR: could not open output cssr file with name %s\n", frame_file_name.c_str());
        exit(EXIT_FAILURE);
      }
      write_cssr(frame_file, &framework, net_name, frame_file_name);
      fclose(frame_file);
      printf("\tcrystallographic structure file written to %s\n", frame_file_name.c_str());
      //...print vtk ...
      string vtk_file_name = prefix+"_framework.vtk"; //write framework's unit cell
      FILE *vtk_file = fopen(vtk_file_name.c_str(),"w");
      if(vtk_file==NULL) {
        printf("ERROR: could not open output vtk file with name %s\n", vtk_file_name.c_str());
        exit(EXIT_FAILURE);
      }
      write_unit_cell(vtk_file, &framework, frame_file_name, false); //is_net? false
      fclose(vtk_file);
      printf("\tunit cell written to %s\n", vtk_file_name.c_str());
      //...and print xyz
      string xyz_file_name = prefix+"_framework.xyz"; //write framework's atoms
      FILE *xyz_file = fopen(xyz_file_name.c_str(),"w");
      if(xyz_file==NULL) {
        printf("ERROR: could not open output xyz file with name %s\n", xyz_file_name.c_str());
        exit(EXIT_FAILURE);
      }
      write_vertices(xyz_file, &framework, frame_file_name, rename_vertices_by_connectivity, false); //is_net? false
      fclose(xyz_file);
      printf("\tatom coordinates written to %s\n", xyz_file_name.c_str());
*/

      //print cssr ...
      string frame_file_name = prefix+"_framework.cssr"; //write framework molecule to cssr file
      writeToCSSR(const_cast<char *>(frame_file_name.c_str()), &framework);
      printf("\tcrystallographic structure file written to %s\n", frame_file_name.c_str());
      //...print cssr ...
      string frame_file_name_labeled = prefix+"_framework_labeled.cssr"; //write framework molecule to cssr file with labels!
      writeToCSSRLabeled(const_cast<char *>(frame_file_name_labeled.c_str()), &framework);
      printf("\tcrystallographic structure file with labels (if they were provided in input files) written to %s\n", frame_file_name_labeled.c_str());
      //...print vtk ...
      string vtk_file_name = prefix+"_framework.vtk"; //write framework's unit cell
      writeToVTK(const_cast<char *>(vtk_file_name.c_str()), &framework);
      printf("\tunit cell written to %s\n", vtk_file_name.c_str());
      //...and print xyz
      string xyz_file_name = prefix+"_framework.xyz"; //write framework's atoms
      writeToXYZ(const_cast<char *>(xyz_file_name.c_str()), &framework, false, false); //falses indicate 'is_supercell' and 'is_duplicate_perimeter_atoms'
      printf("\tatom coordinates written to %s\n", xyz_file_name.c_str());

    } else { //notify that no framework could be achieved
      if(poor_quality_framework && !allow_poor_quality_framework) {
        printf("#####\nREJECTED: this framework has a site deviation larger than the tolerance (%.3f > %.3f) - framework could be incorrect, checking the output is recommended\n#####\n", previous_lowest_deviation, SITE_DEVIATION_WARNING);
      }
      printf("no framework could be achieved\n");
    }
  }

/*------------------------------
           COMPLETE!
------------------------------*/

  //end
  gettimeofday(&time_at_end, NULL);
  seconds_microseconds time_so_far = get_time_difference(time_at_start, time_at_end);
  int milliseconds = (int)round(((double)(time_so_far.microseconds))/1000);
  printf("%d.%03d seconds elapsed\n", time_so_far.seconds, milliseconds);
  printf("Program complete\n\n");
}

/*------------------------------
       FRAMEWORK ASSEMBLY (part of optional framework handling)
------------------------------*/

void construct_framework_from_alignment(vector< pair<int,int> > this_alignment, int index_of_this_framework, vector<int> *count_each_molecule, double *previous_lowest_deviation, ATOM_NETWORK *framework, int *assembly_method, int num_full_vertices, ATOM_NETWORK full_cell, vector< vector<int> > mol_index_at_each_basic_vertex, vector< vector< vector<MOLECULE> > > basic_oriented_molecules, string prefix, string net_name, vector<CONNECTION> two_way_connections, double unit_edge_length, double layer_separation_A, bool layer_override, int dimensionality, bool a_periodicity, bool b_periodicity, bool c_periodicity, bool disconnected_net) {
  bool verbose = false;
  //we want to use symmetry information on the selected rotated molecules to get molecules oriented to each of the vertices in the full cell
if(verbose) printf("DEBUG: attempting to use symmetry on these oriented molecules for each of %d full vertices\n", num_full_vertices);
  vector<MOLECULE> full_oriented_molecules;
  int num_mols = count_each_molecule->size();
  vector<int> local_count_each_molecule;
  for(int i=0; i<num_mols; i++) local_count_each_molecule.push_back(0);
  for(int i=0; i<num_full_vertices; i++) {
    int sym_op = full_cell.vertex_symmetry_operators.at(i);
    int basic_ID = full_cell.vertex_basic_indices.at(i);
    int mol_choice = this_alignment.at(basic_ID).first;
    int mol_index = mol_index_at_each_basic_vertex.at(basic_ID).at(mol_choice);
    int alignment_choice = this_alignment.at(basic_ID).second;
if(verbose) {
  printf("DEBUG: full vertex ID %d\n", i);
  printf("DEBUG: ... has sym_op %d\n", sym_op);
  printf("DEBUG: ... has basic_ID %d\n", basic_ID);
  printf("DEBUG: ... ... and this alignment indicates using molecule choice %d, which gives us molecule ID %d\n", mol_choice, mol_index);
  printf("DEBUG: ... ... with alignment choice %d\n", alignment_choice);
}
    MOLECULE rotated = apply_symmetry_operator(basic_oriented_molecules.at(basic_ID).at(mol_choice).at(alignment_choice), sym_op, full_cell.sym_ID, &full_cell); //apply the appropriate symmetry to the appropriate basic oriented molecule - use the underlying cell info to apply the sym op in the underlying net fractional space!
    full_oriented_molecules.push_back(rotated);
    local_count_each_molecule.at(mol_index)++;
if(verbose) { //write out each rotated molecule for debugging
  string rotmol_name = prefix+"_rotated_molecule_framework_ID_"+convertToString(index_of_this_framework)+"_basic_"+convertToString(basic_ID)+"_sym_"+convertToString(sym_op)+".xyz"; //write molecule to xyz file
  FILE *rotmol = fopen(rotmol_name.c_str(),"w");
  if(rotmol==NULL) {
    printf("ERROR: could not open output rotated molecule file with name %s\n", rotmol_name.c_str());
    exit(EXIT_FAILURE);
  }
  write_molecule(rotmol, &rotated, net_name, basic_ID, sym_op, true); //oriented to basic vertex basic_ID with symmetry sym_op - this will be noted in the file - true indicates 'write out sites'
  fclose(rotmol);
  printf("\tmolecule oriented to full vertex %d (basic vertex %d with symmetry operation %d) written to %s\n", i, basic_ID, sym_op, rotmol_name.c_str());
}
  }
if(verbose) printf("DEBUG: %d molecules oriented to full vertices were achieved\n", (int)(full_oriented_molecules.size()));

  //now, we have all the molecules we need in the correct orientation to construct the framework - we just need to position them relative to each other to achieve a framework
  //new construction strategy 2013-10-18 - build many frameworks, using each possible molecule-vertex alignment, and each of the two construction approaches; output whichever method gives us the lowest site deviation!
  vector<ATOM_NETWORK> all_frameworks;
  vector<int> all_framework_assembly_methods;
  vector<double> all_max_site_deviations;


/* - this was moved up to be outside of this function, and the variable names were changed accordingly - keeping this as a backup for now
  //here we will check whether constructing with this alignment is possible (based on the different kinds of connection sites that may be specified)
  bool possible = true;
  int num_cons = two_way_connections.size();
  for(int i=0; i<num_cons && possible; i++) {
    CONNECTION c = two_way_connections.at(i);
    MOLECULE m1 = full_oriented_molecules.at(c.v1);
    MOLECULE m2 = full_oriented_molecules.at(c.v2);
    int m1_site_index = m1.permutation.at(c.e1);
    int m2_site_index = m2.permutation.at(c.e2);
    int m1_atom_index = m1.sites.at(m1_site_index);
    int m2_atom_index = m2.sites.at(m2_site_index);
    if(m1.atoms_type.at(m1_atom_index)!=m2.atoms_type.at(m2_atom_index)) possible=false;
  }
  if(possible) printf("NOTICE (double check): this alignment should be possible\n");
  else printf("NOTICE (double check): this alignment should NOT be possible\n");

//  for(int method=0; method<2 && possible; method++) { //for each alignment option, we have two methods of construction to try
  for(int method=0; method<1 && possible; method++) { //DEBUG: net-based only
//  for(int method=1; method<2 && possible; method++) { //DEBUG: connection-based only
*/


  for(int method=0; method<2; method++) { //for each alignment option, we have two methods of construction to try
//  for(int method=0; method<1; method++) { //DEBUG: net-based only
//  for(int method=1; method<2; method++) { //DEBUG: connection-based only
    if(disconnected_net && method==1) {
      //do nothing - we cannot assemble disconnected nets with the connection-based method
    } else {
      ATOM_NETWORK this_framework;
      vector<MOLECULE> assembled_molecules;
      if(method==0) { //NEW METHOD: framework based on the preserving the underlying unit cell and positioning all vertices accordingly
if(verbose) printf("\ttrying net-based assembly method: based on preserving the underlying net unit cell\n");
        //method first needs to know the 'mean_edge_length' of the net, which, as we are assuming edge-transitivity, is essentially just the 'edge_length'
        vector<MOLECULE> assemble_two_molecules;
        double edge_length = determine_edge_length(&full_oriented_molecules, &full_cell, &assemble_two_molecules, &two_way_connections);
        //go ahead and build, if we can (a valid edge length was identified)
        if(edge_length>0) this_framework = build_framework(&full_oriented_molecules, &full_cell, &assembled_molecules, edge_length, unit_edge_length);
      } else { //OLD METHOD: framework based on 'snapping' vertices together - may not preserve the unit cell
if(verbose)printf("\ttrying connection-based assembly method: based on connecting adjacent molecules and then inferring a unit cell\n");
        //method produces a 'mean_edge_length'
        double mean_edge_length = -1;
        //go ahead and build
        this_framework = connect_molecules(&full_oriented_molecules, &full_cell, &assembled_molecules, &two_way_connections, &mean_edge_length, dimensionality, a_periodicity, b_periodicity, c_periodicity, unit_edge_length);
      }
      //it's possible we returned an empty/incomplete framework in an error case - need to handle this safely; only proceed if valid framework was provided
      if(this_framework.valid) {
        //at this point, whichever assembly method we used, we have a valid framework
        //given that we have a framework, we need to do/check some more things
        //A) if this is a 2D (layer) structure, and an override for layer separation distance was provided, use it
        if(dimensionality==2 && layer_override) {
  if(verbose) printf("DEBUG: using provided layer separation distance of %.3f A to override the layer separation provided in the net file\n", layer_separation_A);
          int num_override = 0;
          if(!a_periodicity) {
            this_framework.a = layer_separation_A;
            num_override++;
          }
          if(!b_periodicity) {
            this_framework.b = layer_separation_A;
            num_override++;
          } 
          if(!c_periodicity) {
            this_framework.c = layer_separation_A;
            num_override++;
          }
          if(num_override!=1) {
            printf("ERROR: was expecting to override exactly 1 cell side length value, but %d were overwritten - this is a bug\n", num_override);
            exit(EXIT_FAILURE);
          }
          //now need to re-initialize the cell to update matrices etc!
          this_framework.initialize();
          //this involves updating the fractional coords of all the atoms too ...
          this_framework.update_atom_fractional_coords();
        }
        //B) check for collision in the framework (i.e. the final unit cell and list of atoms within it) based on list of two way connections and xyz coords of molecule atoms
        bool collision = true; //assume not a good structure - only if the matrices are valid do we check for collisions
        if(this_framework.valid && !bypass_collision_detection) { //check again - layer change may have caused a problem
          collision = check_for_collision(&this_framework, &assembled_molecules, &two_way_connections);
        }
        if(!collision || bypass_collision_detection) { //if no collisions were incurred in this structure, we can determine the 'goodness', or rather the 'badness' of the achieved framework by looking at the maximum distance (respecting the periodic boundary) between sites that should be connected
  if(verbose) printf("DEBUG: no collision was detected in the framework - this is a candidate for writing out - determining site deviation\n");
          int num_two_way_connections = two_way_connections.size();
          double max_site_deviation = -1;
          for(int i=0; i<num_two_way_connections; i++) {
            CONNECTION c = two_way_connections.at(i);
            XYZ xyz1 = assembled_molecules.at(c.v1).atoms_xyz.at(assembled_molecules.at(c.v1).sites.at(assembled_molecules.at(c.v1).permutation.at(c.e1)));
            XYZ xyz2 = assembled_molecules.at(c.v2).atoms_xyz.at(assembled_molecules.at(c.v2).sites.at(assembled_molecules.at(c.v2).permutation.at(c.e2)));
            XYZ abc1 = this_framework.xyz_to_abc_returning_XYZ(xyz1);
            XYZ abc2 = this_framework.xyz_to_abc_returning_XYZ(xyz2);
  //          double d = calcDistanceABC(abc1, abc2, &this_framework);
            double d = this_framework.calcDistanceABC(abc1, abc2);
            if(d>max_site_deviation || max_site_deviation<0) max_site_deviation = d;
          }
  if(verbose) printf("DEBUG: maximum site deviation for alignment option %d and assembly method %d (0 is net-based, 1 is connection-based) is %.3f\n", index_of_this_framework/*alignment_option*/, method, max_site_deviation);
          all_max_site_deviations.push_back(max_site_deviation);
          //before pushing back the framework, make sure it is centred
          int lowest_c = 0;
          int chosen_full_vertex = -1;
          for(int i=0; i<num_full_vertices; i++) if(full_cell.vertices.at(i).edges.size()<lowest_c || lowest_c==0) {chosen_full_vertex = i; lowest_c = full_cell.vertices.at(i).edges.size();}
          MOLECULE anchor = assembled_molecules.at(chosen_full_vertex);
          XYZ anchor_CoM(0,0,0);
          for(int s=0; s<anchor.sites.size(); s++) anchor_CoM = anchor_CoM + anchor.atoms_xyz.at(anchor.sites.at(s));
          anchor_CoM = anchor_CoM.scale(1.0/anchor.sites.size());
          XYZ anchor_CoM_abc = this_framework.xyz_to_abc_returning_XYZ(anchor_CoM);
          XYZ target_abc = full_cell.vertices.at(chosen_full_vertex).abc;
          XYZ abc_shift = target_abc - anchor_CoM_abc;
  //        for(int a=0; a<this_framework.vertices.size(); a++) this_framework.vertices.at(a).abc = trans_to_origuc(this_framework.vertices.at(a).abc+abc_shift);
          for(int a=0; a<this_framework.atoms.size(); a++) {
            XYZ orig_atom_abc(this_framework.atoms.at(a).a_coord, this_framework.atoms.at(a).b_coord, this_framework.atoms.at(a).c_coord);
            XYZ new_atom_abc = trans_to_origuc(orig_atom_abc+abc_shift);
            XYZ new_atom_xyz = this_framework.abc_to_xyz_returning_XYZ(new_atom_abc);
            this_framework.atoms.at(a).a_coord = new_atom_abc.x;
            this_framework.atoms.at(a).b_coord = new_atom_abc.y;
            this_framework.atoms.at(a).c_coord = new_atom_abc.z;
            this_framework.atoms.at(a).set_xyz(new_atom_xyz);
          }
          //push finished framework
          all_frameworks.push_back(this_framework);
          all_framework_assembly_methods.push_back(method);
        } else if(verbose) printf("\t\ta collision was detected during this framework alignment option - this is not a candidate for writing out\n");
      } else if(verbose) printf("\t\tan incomplete or empty framework was returned - ignoring this alignment combination\n");
    }
  }

  //now, we have evaluated all possible frameworks, and have a list of candidates - find the one with the lowest max_site_deviation and write it out
  int num_valid_frameworks = all_frameworks.size();
  if(num_valid_frameworks==0) {
    if(verbose) printf("NOTICE: no valid frameworks could be achieved for alignment option %d using either assembly method\n", index_of_this_framework);
  } else if(num_valid_frameworks!=all_max_site_deviations.size()) {
    printf("ERROR: %d valid frameworks were produced but an unequal number %d of maximum site deviations were calculated\n", num_valid_frameworks, (int)(all_max_site_deviations.size()));
    exit(EXIT_FAILURE);
  } else {
    int index_best_framework = -1;
    double best_deviation = -1;
    for(int i=0; i<num_valid_frameworks; i++) {
      double dev = all_max_site_deviations.at(i);
if(verbose) {
  if(all_framework_assembly_methods.at(i)==0) printf("\tnet-"); else printf("\tconnection-");
  printf("based assembly method results in a maximum site deviation of %.3f\n", dev);
}
      if(dev<best_deviation || best_deviation<0) {best_deviation = dev; index_best_framework = i;}
    }
    //now, only if the best framework here is better that anything we had before, do we 'return' it using the pointer
    if(best_deviation<(*previous_lowest_deviation) || (*previous_lowest_deviation)<0) {
if(verbose) {
  printf("\tthe new best framework so far has site deviation %.3f, using the", best_deviation);
  if(all_framework_assembly_methods.at(index_best_framework)==0) printf(" net-"); else printf(" connection-");
  printf("based assembly method\n");
}
      (*framework) = all_frameworks.at(index_best_framework);
      (*count_each_molecule) = local_count_each_molecule;
      (*previous_lowest_deviation) = best_deviation;
      (*assembly_method) = all_framework_assembly_methods.at(index_best_framework);
    }
  }
}

/*------------------------------
     OTHER MINOR FUNCTIONS
------------------------------*/

//this function attempts to iterate on a provided alignment code, returning false if there are no more to explore
bool iterate_alignment_code(vector<int> *example_alignment_code, vector< vector< pair<int,int> > > *alignment_options_each_basic_vertex) {
  bool done = false;
  bool success = false;
  //iterating the code involves iteratively incrementing the first entry, and if it is not valid, setting it to zero and iterating the next
  int iteration_position = 0;
  while(!done) {
    example_alignment_code->at(iteration_position)++;
    if(example_alignment_code->at(iteration_position)>=alignment_options_each_basic_vertex->at(iteration_position).size()) { //overstepped - set to zero and try the next one
      example_alignment_code->at(iteration_position) = 0;
      iteration_position++;
      //if the new iteration_position is invalid, we have exhausted our search
      if(iteration_position>=example_alignment_code->size()) done = true;
    } else { //didn't overstep - this is a valid new code
      success = true;
      done = true;
    }
  }
  return success;
}

//compute the difference between two timestamps
seconds_microseconds get_time_difference(struct timeval start, struct timeval end) {
  seconds_microseconds difference;
  if(start.tv_sec==end.tv_sec) {
    difference.seconds = 0;
    difference.microseconds = end.tv_usec-start.tv_usec;
  } else {
    difference.microseconds = 1000000-start.tv_usec;
    difference.seconds = end.tv_sec-(start.tv_sec+1);
    difference.microseconds += end.tv_usec;
    if(difference.microseconds >= 1000000) {
      difference.microseconds -= 1000000;
      difference.seconds += 1;
    }
  }
  return difference;
}

