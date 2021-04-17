//-----Richard Luis Martin, Lawrence Berkeley National Laboratory, 2013
/* this file contains classes and functions required for handling topology (net) information */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "net.h"
#include "networkstorage.h"
#include "symmetry.h"

/* FIT sorting comparator */
bool increasing_rmsd (FIT i, FIT j) { return (i.rmsd<j.rmsd); }

/* VERTEX constructors */
VERTEX::VERTEX(){
}
VERTEX::VERTEX(double a0, double b0, double c0){
  expected_edges = 0;
  XYZ abc1(a0,b0,c0);
  abc = abc1;
  name = "V";
}

/* CONNECTION constructor */
CONNECTION::CONNECTION(int v1_, int v2_, int e1_, int a_, int b_, int c_) {
  v1 = v1_;
  v2 = v2_;
  e1 = e1_;
  e2 = -1; //IMPORTANT! when the CONNECTION is constructed, we don't yet know the edge ID of the other side
  a = a_;
  b = b_;
  c = c_;
}

/* method bodies */

/* take a molecule, and return a new MOLECULE rotated by the given matrix */
MOLECULE rotate(MOLECULE orig, double mat[3][3]) {
  MOLECULE rot = orig;
  int num_atom = orig.atoms_xyz.size();
  for(int i=0; i<num_atom+1; i++) {
    XYZ o;
    if(i<num_atom) o = orig.atoms_xyz.at(i); else o = orig.com; //handle the centre of mass too!
    XYZ n(
mat[0][0]*o.x + mat[0][1]*o.y + mat[0][2]*o.z,
mat[1][0]*o.x + mat[1][1]*o.y + mat[1][2]*o.z,
mat[2][0]*o.x + mat[2][1]*o.y + mat[2][2]*o.z);
    if(i<num_atom) rot.atoms_xyz.at(i) = n; else rot.com = n; //handle the centre of mass too!
  }
  return rot;
}

/* returns center of mass of a molecule's connection sites */
XYZ get_mol_site_CoM(MOLECULE *mol) {
  int num_mol_sites_no_dummy = mol->sites.size();
  XYZ mol_CoM(0,0,0);
  for(int i=0; i<num_mol_sites_no_dummy; i++) {
    mol_CoM = mol_CoM+mol->atoms_xyz.at(mol->sites.at(i));
  }
  mol_CoM = mol_CoM.scale(1.0/num_mol_sites_no_dummy);
  return mol_CoM;
}

/* returns whether two connections are equivalent */
bool matches(CONNECTION ci, CONNECTION cj) {
  bool match = false;
  if(ci.v1==cj.v2 && ci.v2==cj.v1) {
    if(ci.a+cj.a==0 && ci.b+cj.b==0 && ci.c+cj.c==0) match=true;
  }
  return match;
}

/* returns whether two Cartesian points are bonded */
bool bonded_xyz(XYZ xyz1, XYZ xyz2) {
  return bonded_xyz(xyz1, xyz2, GENERAL_BOND_DISTANCE);
}
bool bonded_xyz(XYZ xyz1, XYZ xyz2, double threshold) {
  return (xyz1-xyz2).magnitude()<threshold;
}
/* returns whether two Cartesian points overlap */
bool overlaps_xyz(XYZ xyz1, XYZ xyz2) {
  return (xyz1-xyz2).magnitude()<DISTANCE_TOLERANCE;
}

/* move a molecule by a given shift */
MOLECULE translate(MOLECULE m, XYZ t) {
  MOLECULE new_m = m;
  int num_atoms = m.atoms_xyz.size();
  for(int i=0; i<num_atoms; i++) {
    XYZ new_pos = m.atoms_xyz.at(i)+t;
    new_m.atoms_xyz.at(i) = new_pos;
  }
  XYZ new_com = m.com+t;
  new_m.com = new_com;
  return new_m;
}

/* calculate permutations of a set of integers recursively */
void permute(vector<int> vec, int pos, vector< vector<int> > *out_vec) {
  if (pos == vec.size()) out_vec->push_back(vec);
  else {
    for (int j = pos; j < vec.size(); j++) {
      swap(&vec, pos, j);
      permute(vec, pos+1, out_vec);
      swap(&vec, pos, j); //backtrack
    }
  }
}

/* simple vector swap method for integers */
void swap(vector<int> *vec, int a, int b) {
  int temp = vec->at(a);
  vec->at(a) = vec->at(b);
  vec->at(b) = temp;
}

/* places dummy sites in a molecule where necessary to enable alignment */
int assign_dummy_site(MOLECULE *m) {
  int num_sites = m->sites.size();
  int num_atoms = m->atoms_xyz.size();
  if(num_sites!=2) {
    printf("ERROR: should not call assign_dummy_site() with num_sites!=2\n");
    exit(EXIT_FAILURE);
  }
  //ok, we assign the dummy site to create a 'regular T' shape with respect to the 2 actual sites, i.e., find the midpoint of the sites and the half-width, and move in an orthogonal direction by that same distance
  XYZ site1 = m->atoms_xyz.at(m->sites.at(0));
  XYZ site2 = m->atoms_xyz.at(m->sites.at(1));
  XYZ midpoint = (site1+site2).scale(0.5);
  double half_width = get_vector_from_to(midpoint, site1).magnitude();
  //find the orthogonal direction by reference to the atoms of this molecule
  int furthest_atom_index = -1;
  double furthest_dist = -1;
  XYZ furthest_atom_xyz(0,0,0);
  XYZ furthest_projection(0,0,0);
  for(int i=0; i<num_atoms; i++) {
    XYZ this_atom_xyz = m->atoms_xyz.at(i);
    XYZ projection = project_onto_line(this_atom_xyz, site1, site2);
    double distance = get_vector_from_to(projection, this_atom_xyz).magnitude();
    if(furthest_dist<0 || distance>furthest_dist) {
      if(distance>DISTANCE_TOLERANCE) { //if the molecule is essentially a rod, we will not be able to find an appropriate orthogonal direction - using a tolerance, we can identify this scenario and later just choose an arbitrary orthogonal direction (since it won't matter for a rod anyway)
        furthest_dist = distance;
        furthest_atom_index = i;
        furthest_atom_xyz = this_atom_xyz;
        furthest_projection = projection;
      }
    }
  }
  if(furthest_atom_index==-1) { //could not find orthogonal atom, perhaps because of rod-shape molecule
    XYZ rand1(10,0,0);
    XYZ rand2(0,10,0);
    XYZ rand3(0,0,10); //three possible orthogonal positions - they are not in a line, so one should be valid!
    vector<XYZ> rands;
    rands.push_back(rand1);
    rands.push_back(rand2);
    rands.push_back(rand3);
    int num_rands = rands.size();
//    bool found = false;
    for(int i=0; i<num_rands/* && !found*/; i++) {
      XYZ this_rand = rands.at(i);
      XYZ projection = project_onto_line(this_rand, site1, site2);
      double distance = get_vector_from_to(projection, this_rand).magnitude();
      if(distance>DISTANCE_TOLERANCE) {
        XYZ projection_vector_norm = get_vector_from_to(this_rand, projection).unit();
        XYZ orthogonal_site_xyz = midpoint+projection_vector_norm.scale(half_width);
        m->atoms_xyz.push_back(orthogonal_site_xyz);
        m->atoms_type.push_back("J");
        m->atoms_label.push_back("J");
        return num_atoms;
//        m->dummy_sites.push_back(orthogonal_site_xyz); //done
//        found = true;
      }
    }
  } else { //else we found a good atom, so we can use it
    XYZ projection_vector_norm = get_vector_from_to(furthest_atom_xyz, furthest_projection).unit();
    XYZ orthogonal_site_xyz = midpoint+projection_vector_norm.scale(half_width);
    m->atoms_xyz.push_back(orthogonal_site_xyz);
    m->atoms_type.push_back("J");
    m->atoms_label.push_back("J");
    return num_atoms;
//    m->dummy_sites.push_back(orthogonal_site_xyz); //done
  }
}

/* defines a MOLECULE's center of mass (site centroid, excluding dummy sites) */
void get_com(MOLECULE *m) {
  bool verbose = false;
  int num_sites = m->sites.size();
  XYZ com(0,0,0);
  for(int i=0; i<num_sites; i++) {
    int site_index = m->sites.at(i);
    //update centre of mass accordingly
    com = com + m->atoms_xyz.at(site_index);
  }
  //finalise centre of mass by normalising it
  com = com.scale(1.0/((double)num_sites));
  m->com = com;
  if(verbose) printf("DEBUG: centre of mass of this molecule is at %.3f %.3f %.3f\n", com.x, com.y, com.z);
}

/* a recursive algorithm to test whether a provided set of two-way connections constitutes a disconnected net */
void recursive_visit_vertices(int this_v, vector<CONNECTION> *two_way_connections, vector<bool> *visited_vertex) {
  visited_vertex->at(this_v) = true;
  int num_c = two_way_connections->size();
  for(int i=0; i<num_c; i++) {
    CONNECTION c = two_way_connections->at(i);
    if(c.v1==this_v && !visited_vertex->at(c.v2)) recursive_visit_vertices(c.v2, two_way_connections, visited_vertex);
    if(c.v2==this_v && !visited_vertex->at(c.v1)) recursive_visit_vertices(c.v1, two_way_connections, visited_vertex);
  }
}

/* new, hopefully more efficient recursive loop finder */
bool recursive_find_loops(int this_v, int a, int b, int c, XYZ xyz, vector<CONNECTION> *two_way_connections, vector<MOLECULE> *assembled_molecules, vector<bool> *edge_traversed, vector<int> *vertex_visit_a, vector<int> *vertex_visit_b, vector<int> *vertex_visit_c, vector<XYZ> *vertex_visit_xyz, vector<bool> *vertex_visited, vector<int> *loops_a, vector<int> *loops_b, vector<int> *loops_c, vector<XYZ> *loops_xyz, vector<int> *unit_cell_vector_IDs, vector<XYZ> *unit_cell_vectors, int dimensionality) {
  bool verbose = false;
  //visit this vertex
  if(vertex_visited->at(this_v)) {
    //previously visited this vertex - have we made a loop?
    int old_a = vertex_visit_a->at(this_v);
    int old_b = vertex_visit_b->at(this_v);
    int old_c = vertex_visit_c->at(this_v);
    XYZ old_xyz = vertex_visit_xyz->at(this_v);
    int period_a = a-old_a;
    int period_b = b-old_b;
    int period_c = c-old_c;
    XYZ period_xyz = xyz-old_xyz;
    if(period_a!=0 || period_b!=0 || period_c!=0) { //it's a loop!
if(verbose) printf("DEBUG: a loop was found with periodicity (%d %d %d) and shift %.3f %.3f %.3f\n", period_a, period_b, period_c, period_xyz.x, period_xyz.y, period_xyz.z);
      bool found_all_uc_vectors = false;
      //is this new loop unique?
      if(loop_is_unique(period_a, period_b, period_c, loops_a, loops_b, loops_c)) {
if(verbose) printf("DEBUG: new loop is unique, testing for unit cell vector assignment ... \n");
        //does this new loop allow us to define a unit cell vector?
        found_all_uc_vectors = try_unit_cell_vector_assign(period_a, period_b, period_c, period_xyz, loops_a, loops_b, loops_c, loops_xyz, unit_cell_vector_IDs, unit_cell_vectors, dimensionality);
        //push loop
if(verbose) printf("DEBUG: ... pushing the new loop\n");
        loops_a->push_back(period_a);
        loops_b->push_back(period_b);
        loops_c->push_back(period_c);
        loops_xyz->push_back(period_xyz);
      } else {
if(verbose) printf("DEBUG: new loop is a duplicate\n");
      }
      return found_all_uc_vectors; //quit the recursive method here - we have a loop, time to try and find some other loops
    } else {
      //else do nothing - it's not a loop, so just proceed with journey from this vertex
    }
  } else {
    //haven't previously visited this vertex; visit it
    vertex_visited->at(this_v) = true;
    vertex_visit_a->at(this_v) = a;
    vertex_visit_b->at(this_v) = b;
    vertex_visit_c->at(this_v) = c;
    vertex_visit_xyz->at(this_v) = xyz;
  }
  //at this point, we have either completed a loop, or visited a vertex that did not complete a loop; only if a loop was found do we quit out - else proceed to recursively traverse other connections
  bool found_all_uc_vectors = false; //if we find them all during this run, terminate early
  int num_c = two_way_connections->size();
  for(int i=0; i<num_c && !found_all_uc_vectors; i++) {
    if(!edge_traversed->at(i)) {
      //only traverse new edges
      CONNECTION con = two_way_connections->at(i);
      if(con.v1==this_v) {
        //traverse this edge
        edge_traversed->at(i) = true;
        //don't visit the vertex - that will happen in next call
        int new_a = a+con.a;
        int new_b = b+con.b;
        int new_c = c+con.c;
        int a_edge_index = con.e1;
        int b_edge_index = con.e2;
        int a_vertex_index = con.v1;
        int b_vertex_index = con.v2;
        MOLECULE a_mol = assembled_molecules->at(a_vertex_index);
        int a_site_index = a_mol.permutation.at(a_edge_index);
        int a_atom_index = a_mol.sites.at(a_site_index);
        XYZ a_site = a_mol.atoms_xyz.at(a_atom_index);
        MOLECULE b_mol = assembled_molecules->at(b_vertex_index);
        int b_site_index = b_mol.permutation.at(b_edge_index);
        int b_atom_index = b_mol.sites.at(b_site_index);
        XYZ b_site = b_mol.atoms_xyz.at(b_atom_index);
        XYZ shift = b_site-a_site;
        XYZ new_xyz = xyz-shift;
        found_all_uc_vectors = recursive_find_loops(con.v2, new_a, new_b, new_c, new_xyz, two_way_connections, assembled_molecules, edge_traversed, vertex_visit_a, vertex_visit_b, vertex_visit_c, vertex_visit_xyz, vertex_visited, loops_a, loops_b, loops_c, loops_xyz, unit_cell_vector_IDs, unit_cell_vectors, dimensionality);
      } else if(con.v2==this_v) {
        //traverse this edge, pointing the other way
        edge_traversed->at(i) = true;
        //don't visit the vertex - that will happen in next call
        int new_a = a-con.a;
        int new_b = b-con.b;
        int new_c = c-con.c;
        int a_edge_index = con.e1;
        int b_edge_index = con.e2;
        int a_vertex_index = con.v1;
        int b_vertex_index = con.v2;
        MOLECULE a_mol = assembled_molecules->at(a_vertex_index);
        int a_site_index = a_mol.permutation.at(a_edge_index);
        int a_atom_index = a_mol.sites.at(a_site_index);
        XYZ a_site = a_mol.atoms_xyz.at(a_atom_index);
        MOLECULE b_mol = assembled_molecules->at(b_vertex_index);
        int b_site_index = b_mol.permutation.at(b_edge_index);
        int b_atom_index = b_mol.sites.at(b_site_index);
        XYZ b_site = b_mol.atoms_xyz.at(b_atom_index);
        XYZ shift = b_site-a_site;
        XYZ new_xyz = xyz+shift;
        found_all_uc_vectors = recursive_find_loops(con.v1, new_a, new_b, new_c, new_xyz, two_way_connections, assembled_molecules, edge_traversed, vertex_visit_a, vertex_visit_b, vertex_visit_c, vertex_visit_xyz, vertex_visited, loops_a, loops_b, loops_c, loops_xyz, unit_cell_vector_IDs, unit_cell_vectors, dimensionality);
      }
    }
  }
  return found_all_uc_vectors;
}

/* determines whether a loop is unique - considers equivalance, reverse, and other multiplicity */
bool loop_is_unique(int a, int b, int c, vector<int> *loops_a, vector<int> *loops_b, vector<int> *loops_c) {
  bool verbose = false;
  bool unique = true; //unique until proven otherwise
  int num_loops = loops_a->size();
if(verbose) printf("DEBUG: testing new loop uniqueness versus %d existing loops\n", num_loops);
  for(int i=0; i<num_loops && unique; i++) {
if(verbose) printf("DEBUG: checking against existing loop %d: (%d %d %d)\n", i, loops_a->at(i), loops_b->at(i), loops_c->at(i));
    //first, test for plain equivalence
    if(a==loops_a->at(i) && b==loops_b->at(i) && c==loops_c->at(i)) {
      unique = false;
if(verbose) printf("DEBUG: trivial equivalence\n");
    } else {
      //second, we are definitely unique if there are any zero-nonzero pairings
      vector<double> ratios;
      bool zero_nonzero_pairing = false;
      if((a==0 && loops_a->at(i)!=0) || (a!=0 && loops_a->at(i)==0)) zero_nonzero_pairing = true; else if(a!=0 && loops_a->at(i)!=0) ratios.push_back(((double)(a))/loops_a->at(i));
      if((b==0 && loops_b->at(i)!=0) || (b!=0 && loops_b->at(i)==0)) zero_nonzero_pairing = true; else if(b!=0 && loops_b->at(i)!=0) ratios.push_back(((double)(b))/loops_b->at(i));
      if((c==0 && loops_c->at(i)!=0) || (c!=0 && loops_c->at(i)==0)) zero_nonzero_pairing = true; else if(c!=0 && loops_c->at(i)!=0) ratios.push_back(((double)(c))/loops_c->at(i));
      if(!zero_nonzero_pairing) {
        //third, given that there was no zero-nonzero pairing, check the ratios to see if one is a multiple of the other
        int num_ratios = ratios.size();
        if(num_ratios==0) {
          printf("ERROR: passed trivial equivalence filter but both loops are completely zero (no ratios between elements could be found)\n");
          exit(EXIT_FAILURE);
        }
        double ratio = ratios.at(0); //there should always be at least one ratio, otherwise we would have been caught by the trivial equivalence filter
        unique = false; //here we assume not unique, and it becomes unique if any ratios are different
        for(int j=1; j<num_ratios && !unique; j++) {
          if(ratio!=ratios.at(j)) unique = true;
        }
if(verbose && !unique) printf("DEBUG: equivalence by ratio (%.3f)\n", ratio);
      }
    }
  }
  return unique;
}

/* new, hopefully more efficient recursive loop finder */
int recursive_test_dimensionality(int this_v, int a, int b, int c, vector<CONNECTION> *two_way_connections, vector<bool> *edge_traversed, vector<int> *vertex_visit_a, vector<int> *vertex_visit_b, vector<int> *vertex_visit_c, vector<bool> *vertex_visited, vector<bool> *is_periodic) {
  bool verbose = false;
  //visit this vertex
  if(vertex_visited->at(this_v)) {
    //previously visited this vertex - have we made a loop?
    int old_a = vertex_visit_a->at(this_v);
    int old_b = vertex_visit_b->at(this_v);
    int old_c = vertex_visit_c->at(this_v);
    int period_a = a-old_a;
    int period_b = b-old_b;
    int period_c = c-old_c;
    if(period_a!=0 || period_b!=0 || period_c!=0) { //it's a loop!
if(verbose) printf("DEBUG: a loop was found with periodicity (%d %d %d)\n", period_a, period_b, period_c);
      //update periodicity
      if(period_a!=0) is_periodic->at(0) = true;
      if(period_b!=0) is_periodic->at(1) = true;
      if(period_c!=0) is_periodic->at(2) = true;
      int dimensionality = 0;
      for(int i=0; i<3; i++) {
        if(is_periodic->at(i)) dimensionality++;
      }
      return dimensionality; //quit the recursive method here - we have a loop, time to try and find some other loops
    } else {
      //else do nothing - it's not a loop, so just proceed with journey from this vertex
    }
  } else {
    //haven't previously visited this vertex; visit it
    vertex_visited->at(this_v) = true;
    vertex_visit_a->at(this_v) = a;
    vertex_visit_b->at(this_v) = b;
    vertex_visit_c->at(this_v) = c;
  }
  //at this point, we have either completed a loop, or visited a vertex that did not complete a loop; only if a loop was found do we quit out - else proceed to recursively traverse other connections - if we find three-dimensionality during this run, terminate early
  int num_c = two_way_connections->size();
  int dimensionality = 0;
  for(int i=0; i<3; i++) {
    if(is_periodic->at(i)) dimensionality++;
  }
  for(int i=0; i<num_c && dimensionality<3; i++) {
    if(!edge_traversed->at(i)) {
      //only traverse new edges
      CONNECTION con = two_way_connections->at(i);
      if(con.v1==this_v) {
        //traverse this edge
        edge_traversed->at(i) = true;
        //don't visit the vertex - that will happen in next call
        int new_a = a+con.a;
        int new_b = b+con.b;
        int new_c = c+con.c;
        dimensionality = recursive_test_dimensionality(con.v2, new_a, new_b, new_c, two_way_connections, edge_traversed, vertex_visit_a, vertex_visit_b, vertex_visit_c, vertex_visited, is_periodic);
      } else if(con.v2==this_v) {
        //traverse this edge, pointing the other way
        edge_traversed->at(i) = true;
        //don't visit the vertex - that will happen in next call
        int new_a = a-con.a;
        int new_b = b-con.b;
        int new_c = c-con.c;
        dimensionality = recursive_test_dimensionality(con.v1, new_a, new_b, new_c, two_way_connections, edge_traversed, vertex_visit_a, vertex_visit_b, vertex_visit_c, vertex_visited, is_periodic);
      }
    }
  }
  return dimensionality;
}

/* tests whether the positions and types of atoms in an aligned molecule are unique compared to previously generated alignments */
bool molecule_alignment_chemistry_is_unique(MOLECULE *mol, vector<MOLECULE> *existing) {
  bool verbose = false;
  bool unique = true;
  int num_existing = existing->size();
  int num_atoms = mol->atoms_xyz.size();
  int num_sites_no_dummy = mol->sites.size();
  int num_sites_dummy = mol->dummy_sites.size();
  int num_sites = num_sites_no_dummy + num_sites_dummy;
  //1) we first check if the sites are aligned
  vector<XYZ> mol_site_positions;
  for(int i=0; i<num_sites; i++) {
    if(i<num_sites_no_dummy)
      mol_site_positions.push_back(mol->atoms_xyz.at(mol->sites.at(i)));
    else
      mol_site_positions.push_back(mol->dummy_sites.at(i-num_sites_no_dummy));
  }
if(verbose) {
  printf("DEBUG: mol sites are positioned as follows:\n");
  for(int i=0; i<num_sites; i++) printf("\t(%.3f,%.3f,%.3f)\n", mol_site_positions.at(i).x, mol_site_positions.at(i).y, mol_site_positions.at(i).z);
}
  for(int i=0; i<num_existing && unique; i++) {
    vector<XYZ> existing_site_positions;
    for(int j=0; j<num_sites; j++) {
      if(j<num_sites_no_dummy) {
        existing_site_positions.push_back(existing->at(i).atoms_xyz.at(existing->at(i).sites.at(j)));
      } else {
        existing_site_positions.push_back(existing->at(i).dummy_sites.at(j-num_sites_no_dummy));
      }
    }
if(verbose) {
  printf("DEBUG: comparing against stored mol sites as follows:\n");
  for(int d=0; d<num_sites; d++) printf("\t(%.3f,%.3f,%.3f)\n", existing_site_positions.at(d).x, existing_site_positions.at(d).y, existing_site_positions.at(d).z);
}
    vector<bool> existing_site_partner_found;
    for(int j=0; j<num_sites; j++) existing_site_partner_found.push_back(false);
    double furthest_site_match = -1;
    for(int m=0; m<num_sites; m++) {
      //for each site in the molecule, find the closest site in the existing one, that has not already been assigned
      double closest_distance = -1;
      int closest_id = -1;
      for(int e=0; e<num_sites; e++) {
        if(!existing_site_partner_found.at(e)) {
          double distance = get_vector_from_to(mol_site_positions.at(m), existing_site_positions.at(e)).magnitude();
          if(distance<closest_distance || closest_distance<0) {
            closest_distance = distance;
            closest_id = e;
          }
        }
      }
      existing_site_partner_found.at(closest_id) = true;
      if(closest_distance>furthest_site_match) furthest_site_match = closest_distance;
    }
    //now we have matched all sites in the molecule to the existing one - if the furthest distance we stored is shorter than the threshold, these are equivalent
if(verbose) printf("DEBUG: furthest site match distance was %.3f\n", furthest_site_match);
    if(furthest_site_match<SAME_POSITION_DISTANCE) unique = false;
    //1b) now for the chemistry bit - if the sites are not uniquely oriented between these two specific molecules, we might still have unique chemistry (because of functional groups or other asymmetry)
    if(!unique) {
      unique = true; //we again assume that the molecule is unique until proven to be otherwise
      vector<bool> existing_atom_partner_found;
      for(int j=0; j<num_atoms; j++) existing_atom_partner_found.push_back(false);
      double furthest_atom_match = -1;
      for(int m=0; m<num_atoms; m++) {
        //for each atom in the molecule, find the closest atom of the same type in the existing one, that has not already been assigned
        double closest_distance = -1;
        int closest_id = -1;
        for(int e=0; e<num_atoms; e++) {
          if(!existing_atom_partner_found.at(e) && existing->at(i).atoms_type.at(e)==mol->atoms_type.at(m)) {
            double distance = get_vector_from_to(mol->atoms_xyz.at(m), existing->at(i).atoms_xyz.at(e)).magnitude();
            if(distance<closest_distance || closest_distance<0) {
              closest_distance = distance;
              closest_id = e;
            }
          }
        }
        existing_atom_partner_found.at(closest_id) = true;
        if(closest_distance>furthest_site_match) furthest_atom_match = closest_distance;
      }
      //now we have matched all atoms in the molecule to the existing one - if the furthest distance we stored is shorter than the threshold, these are equivalent
if(verbose) printf("DEBUG: furthest atom match distance was %.3f\n", furthest_atom_match);
      if(furthest_atom_match<SAME_POSITION_DISTANCE) unique = false;
    }
  }
if(verbose) {
  if(unique) printf("DEBUG: returning unique\n"); else printf("DEBUG: returning NOT unique\n");
}
  return unique;
}

/* determine whether a new loop allows us to assign any unit cell vectors */
bool try_unit_cell_vector_assign(int a, int b, int c, XYZ xyz, vector<int> *loops_a, vector<int> *loops_b, vector<int> *loops_c, vector<XYZ> *loops_xyz, vector<int> *unit_cell_vector_IDs, vector<XYZ> *unit_cell_vectors, int num_uc_vectors_required) {
  bool verbose = false;
  int num_uc_vectors_assigned = unit_cell_vectors->size();
  for(int uc=0; uc<3 && num_uc_vectors_assigned<num_uc_vectors_required; uc++) {
    //is this unit cell vector assigned already?
    bool this_uc_assigned = false;
    for(int i=0; i<num_uc_vectors_assigned && !this_uc_assigned; i++) {
      if(unit_cell_vector_IDs->at(i)==uc) this_uc_assigned = true;
    }
    if(!this_uc_assigned) {
      //if not assigned already, can we assign it based on this loop alone?
      if(
(uc==0 && a!=0 && b==0 && c==0) ||
(uc==1 && a==0 && b!=0 && c==0) ||
(uc==2 && a==0 && b==0 && c!=0)
      ) {
        XYZ vec(0,0,0);
        if(uc==0) vec=xyz.scale(1.0/a); else if(uc==1) vec=xyz.scale(1.0/b); else vec=xyz.scale(1.0/c);
if(verbose) printf("DEBUG: we can assign unit cell vector %c as %.3f %.3f %.3f\n", 'a'+uc, vec.x, vec.y, vec.z);
        unit_cell_vector_IDs->push_back(uc);
        unit_cell_vectors->push_back(vec);
        this_uc_assigned = true;
        num_uc_vectors_assigned++;
      } else {
        //we could not assign it based on this loop alone - but we might be able to based on other loops ...
        int num_loops = loops_xyz->size();
        for(int i=0; i<num_loops && !this_uc_assigned; i++) {
          XYZ other_xyz = loops_xyz->at(i);
          int other_a = loops_a->at(i);
          int other_b = loops_b->at(i);
          int other_c = loops_c->at(i);
if(verbose) printf("DEBUG: comparing to previous loop (%d %d %d) %.3f %.3f %.3f to try and isolate vector %c\n", other_a, other_b, other_c, other_xyz.x, other_xyz.y, other_xyz.z, 'a'+uc);
          //can we isolate this unit cell vector with a linear combination of these two loops?
          //bool can_isolate = true;
          if(
(uc!=0 && a!=0 && other_a==0) ||
(uc!=1 && b!=0 && other_b==0) ||
(uc!=2 && c!=0 && other_c==0)
          ) {
if(verbose) printf("DEBUG: cannot isolate vector %c because some other unit cell vector cannot be reduced to zero\n", 'a'+uc);
          } else if(
(uc==0 && a==0 && other_a==0) ||
(uc==1 && b==0 && other_b==0) ||
(uc==2 && c==0 && other_c==0)
          ) {
if(verbose) printf("DEBUG: cannot isolate vector %c because both coefficients of this vector are zero\n", 'a'+uc);
          } else { //there is another way it might not be possible - if the ratios of the two other unit cell vector coefficients are not equal
if(verbose) printf("DEBUG: may be able to isolate vector %c\n", 'a'+uc);
            vector<double> amounts_of_other_loop_to_subtract;
            vector<double> amount_of_other_loop_that_would_reduce_uc_to_zero;
            if(!(a==0 && other_a==0)) {double r = ((double)(a))/other_a; if(uc!=0) amounts_of_other_loop_to_subtract.push_back(r); else amount_of_other_loop_that_would_reduce_uc_to_zero.push_back(r);}
            if(!(b==0 && other_b==0)) {double r = ((double)(b))/other_b; if(uc!=1) amounts_of_other_loop_to_subtract.push_back(r); else amount_of_other_loop_that_would_reduce_uc_to_zero.push_back(r);}
            if(!(c==0 && other_c==0)) {double r = ((double)(c))/other_c; if(uc!=2) amounts_of_other_loop_to_subtract.push_back(r); else amount_of_other_loop_that_would_reduce_uc_to_zero.push_back(r);}
            int num_amounts = amounts_of_other_loop_to_subtract.size();
            if(num_amounts==0) {
              printf("ERROR: detected that unit cell vector has already been determined in determination routine\n");
              exit(EXIT_FAILURE);
            } else if(num_amounts>2) {
              printf("ERROR: pushed too many solutions to unit cell vector determination routine\n");
              exit(EXIT_FAILURE);
            } else {
              bool can_isolate = true;
              if(num_amounts==2) {
                if(fabs(amounts_of_other_loop_to_subtract.at(0)-amounts_of_other_loop_to_subtract.at(1))>DISTANCE_TOLERANCE) {
                  can_isolate = false;
if(verbose) printf("DEBUG: cannot isolate vector %c because no linear combination will isolate this vector\n", 'a'+uc);
                }
              }
              if(can_isolate) { //we should be able to isolate it - but are we reducing to zero?
                double ratio = amounts_of_other_loop_to_subtract.at(0);
                bool can_isolate_nonzero = true;
                if(amount_of_other_loop_that_would_reduce_uc_to_zero.size()!=0) {
                  if(fabs(ratio-amount_of_other_loop_that_would_reduce_uc_to_zero.at(0))<DISTANCE_TOLERANCE) {
                    can_isolate_nonzero = false;
if(verbose) printf("DEBUG: cannot isolate vector %c because this linear combination reveals that the two vectors are multiples of each other\n", 'a'+uc);
                  }
                }
                if(can_isolate_nonzero) {
                  XYZ pre_scale_vector = xyz-(other_xyz.scale(ratio));
if(verbose) printf("DEBUG: this loop %.3f %.3f %.3f, minus %.3f times other loop %.3f %.3f %.3f, gives isolated %c loop %.3f %.3f %.3f\n", xyz.x, xyz.y, xyz.z, ratio, other_xyz.x, other_xyz.y, other_xyz.z, 'a'+uc, pre_scale_vector.x, pre_scale_vector.y, pre_scale_vector.z);
                  XYZ vec(0,0,0);
                  if(uc==0) vec=pre_scale_vector.scale(1.0/(a-(other_a*ratio))); else if(uc==1) vec=pre_scale_vector.scale(1.0/(b-(other_b*ratio))); else vec=pre_scale_vector.scale(1.0/(c-(other_c*ratio)));
if(verbose) printf("DEBUG: we can assign unit cell vector %c as %.3f %.3f %.3f based on a linear combination with a previous loop\n", 'a'+uc, vec.x, vec.y, vec.z);
                  unit_cell_vector_IDs->push_back(uc);
                  unit_cell_vectors->push_back(vec);
                  this_uc_assigned = true;
                  num_uc_vectors_assigned++;
                  //also push this new loop so we can refer to it for subsequent linear combinations - only if it is unique!
                  int new_a=0, new_b=0, new_c=0;
                  if(uc==0) new_a=1; else if(uc==1) new_b=1; else new_c=1;
                  if(loop_is_unique(new_a, new_b, new_c, loops_a, loops_b, loops_c)) {
if(verbose) printf("DEBUG: pushing the newly defined, unique loop\n");
                    loops_a->push_back(new_a);
                    loops_b->push_back(new_b);
                    loops_c->push_back(new_c);
                    loops_xyz->push_back(vec);
                    num_loops++;
                  }
                }
              }
            }
          } 
        }
      }
    }
  }
  return (num_uc_vectors_assigned==num_uc_vectors_required);
}

/* TO_FIX - a reduced version of the 'images' routine */
vector<XYZ> get_periodic_images_of_uc_abc_position(XYZ e_abc) {
  vector<XYZ> periodic_images;
  for(int i=-2; i<=2; i++) {
    for(int j=-2; j<=2; j++) {
      for(int k=-2; k<=2; k++) {
        XYZ image(e_abc.x+i, e_abc.y+j, e_abc.z+k);
        periodic_images.push_back(image);
      }
    }
  }
  return periodic_images;
}

/* TO_FIX - should this I/O function go in a different file? */
void write_molecule(FILE *rotmol, MOLECULE *rotated, string name, int vertex_ID, int sym_op, bool write_out_sites) {
  bool verbose = false;
  int num_atom = rotated->atoms_xyz.size();
  int num_atom_to_write = num_atom;
  int num_sites = rotated->sites.size();
  if(!write_out_sites) num_atom_to_write-=num_sites;
  fprintf(rotmol, "%d\n%s", num_atom_to_write, name.c_str());
  if(vertex_ID==-1 && sym_op==-1) { //obviously not a symmetry operation ...
    fprintf(rotmol,"\n");
  } else {
    fprintf(rotmol," - molecule rotated to align with basic vertex ID %d and symmetry operator %d\n", vertex_ID, sym_op);
  }
  for(int i=0; i<num_atom; i++) {
    bool write_this_one = true;
    if(!write_out_sites) {
      for(int j=0; j<num_sites && write_this_one; j++) {
        if(i==rotated->sites.at(j)) write_this_one = false;
      }
    }
    if(write_this_one) {
      fprintf(rotmol, "%s %.3f %.3f %.3f\n", rotated->atoms_type.at(i).c_str(), rotated->atoms_xyz.at(i).x, rotated->atoms_xyz.at(i).y, rotated->atoms_xyz.at(i).z);
    } else {
if(verbose) printf("DEBUG: not writing atom ID %d (%s) because it is a site\n", i, rotated->atoms_type.at(i).c_str());
    }
  }
}

/* TO_FIX - general purpose int to string conversion */
string convertToString(int number) {
  stringstream ss;
  ss << number;
  return ss.str();
}

/* TO_FIX - general purpose file parser based on individual characters */
void search_for_char(FILE *f, char target) {
	char c = getc(f);
//printf("DEBUG: found %c with int value %d\n", c, (int)c);
	while(c!=target && c!=EOF) {
		c = getc(f);
//printf("DEBUG: found %c with int value %d\n", c, (int)c);
	}
	if(c==EOF) {
		printf("ERROR: The required character (%c, with int value %d) was not found in this file.\n", target, (int)target);
		exit(EXIT_FAILURE);
	}
}


