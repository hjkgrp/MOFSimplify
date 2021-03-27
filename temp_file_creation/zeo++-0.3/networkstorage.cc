/* Stores information about the basic data types that constitute a unit cell,
 * such as atoms (ATOM) and networks of atoms (ATOM_NETWORK). Also stores
 * information about the underlying Voronoi network constituents such as nodes
 * (VOR_NODE), edges (VOR_EDGE), and the network itself (VORONOI_NETWORK).
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <map>
#include <algorithm>
#include <float.h>
#include "networkstorage.h"
#include "networkinfo.h"
#include "geometry.h"
#include "zeo_consts.h"
#include "rmsd.h"
#include "symmetry.h"
#include "string_additions.h"
#include "ray.h"
using namespace std;

ATOM_NETWORK::ATOM_NETWORK() {highAccuracyFlag = false; valid = false; allowAdjustCoordsAndCellFlag = false;} //the flag 'valid', from the framework builder, describes whether the unit cell makes sense - it has to be false initially, and becomes true if the matrix could be inverted; allowAdjustCoordsAndCellFlag is for the temporary Voronoi volume error fix

ATOM::ATOM() {
  type = "";
  label = "";
  radius = 0;
  charge = 0; //because this is not always provided, a default value of zero is set
  keep = true;
}

ATOM::ATOM(XYZ xyz, string s, double r) {
  x = xyz.x;
  y = xyz.y;
  z = xyz.z;
  type = s;
  label = s;
  radius = r;
  charge = 0;
  keep = true;
}

ATOM::ATOM(XYZ xyz, string s, string l, double r) {
  x = xyz.x;
  y = xyz.y;
  z = xyz.z;
  type = s;
  label = l;
  radius = r;
  charge = 0;
  keep = true;
}

/** Print the information about this atom to the provided output stream. */
void ATOM::print(ostream &out){
    out << "   type:" << type << "   x:" << x << "   y:" << y << "   z:"  << z
    << "   a:" << a_coord << "   b:" << b_coord << "   c:" << c_coord
    << "   radius:" << radius << "\n";
}

XYZ ATOM::xyz() {
  XYZ xyz(x,y,z);
  return xyz;
}
void ATOM::set_xyz(XYZ new_xyz) {
  x = new_xyz.x;
  y = new_xyz.y;
  z = new_xyz.z;
}

/* updates the fractional coords of all atoms in the cell (for example, if you modify the unit cell parameters and want to preserve the xyz coordinates of the atoms) */
void ATOM_NETWORK::update_atom_fractional_coords() {
  for(int i=0; i<numAtoms; i++) {
    XYZ xyz(atoms.at(i).x, atoms.at(i).y, atoms.at(i).z);
    XYZ abc = trans_to_origuc(xyz_to_abc_returning_XYZ(xyz));
    atoms.at(i).a_coord = abc.x;
    atoms.at(i).b_coord = abc.y;
    atoms.at(i).c_coord = abc.z;
  }
}

/* randomly changes the atom coordinates and unit cell definition by a small degree */
void ATOM_NETWORK::randomlyAdjustCoordsAndCell() {
  double shiftMagnitude = TOLERANCE;
  printf("NOTICE: attempting random vector shift of all atom coordinates by %e (and unit cell parameters by up to this amount) to overcome Voronoi volume check failure (this option can be disabled by not using the -allowAdjustCoordsAndCell flag)\n", shiftMagnitude);
  printf("NOTICE: original cell dimensions and angles: %e %e %e; %e %e %e\n", a, b, c, alpha, beta, gamma);
  make(
(((rand()*1.0/RAND_MAX)-0.5)*2.0*shiftMagnitude)+a,
(((rand()*1.0/RAND_MAX)-0.5)*2.0*shiftMagnitude)+b,
(((rand()*1.0/RAND_MAX)-0.5)*2.0*shiftMagnitude)+c,
(((rand()*1.0/RAND_MAX)-0.5)*2.0*shiftMagnitude)+alpha,
(((rand()*1.0/RAND_MAX)-0.5)*2.0*shiftMagnitude)+beta,
(((rand()*1.0/RAND_MAX)-0.5)*2.0*shiftMagnitude)+gamma
);
  printf("NOTICE: randomly adjusted cell dimensions and angles: %e %e %e; %e %e %e\n", a, b, c, alpha, beta, gamma);
  for(int at=0; at<numAtoms; at++) {
    //shift each Cartesian coord in some direction, by distance TOLERANCE
    Point smallShift = genRandomVec().scale(shiftMagnitude);
    Point resultOfShiftInABCNoPBC = xyz_to_abc(Point(
      atoms.at(at).x+smallShift[0],
      atoms.at(at).y+smallShift[1],
      atoms.at(at).z+smallShift[2]
    ));
//    Point resultOfShiftInABC = Point(trans_to_origuc(resultOfShiftInABCNoPBC[0]),trans_to_origuc(resultOfShiftInABCNoPBC[1]),trans_to_origuc(resultOfShiftInABCNoPBC[2]));
    Point resultOfShiftInABC = resultOfShiftInABCNoPBC;
    Point resultOfShiftInXYZ = abc_to_xyz(resultOfShiftInABC);
    //now write the new coords to the atmnet
//printf("DEBUG: atom %d was at %e %e %e, now at %e %e %e\n", at, atoms.at(at).x, atoms.at(at).y, atoms.at(at).z, resultOfShiftInXYZ[0], resultOfShiftInXYZ[1], resultOfShiftInXYZ[2]);
    atoms.at(at).x = resultOfShiftInXYZ[0];
    atoms.at(at).y = resultOfShiftInXYZ[1];
    atoms.at(at).z = resultOfShiftInXYZ[2];
    atoms.at(at).a_coord = resultOfShiftInABC[0];
    atoms.at(at).b_coord = resultOfShiftInABC[1];
    atoms.at(at).c_coord = resultOfShiftInABC[2];
  }
}

/* ATOM_NETWORK pseudo-constructor - assigns values and calls initialization method */
void ATOM_NETWORK::make(double a0, double b0, double c0, double alpha0, double beta0, double gamma0){
  a = a0;
  b = b0;
  c = c0;
  alpha = alpha0;
  beta = beta0;
  gamma = gamma0;
  initialize();
}

/* an alternative way of making a cell - if we know the unit cell vectors */
void ATOM_NETWORK::make(XYZ va, XYZ vb, XYZ vc) {
  v_a = va;
  v_b = vb;
  v_c = vc;
  a = v_a.magnitude();
  b = v_b.magnitude();
  c = v_c.magnitude();
  alpha = 360.0*v_b.angle_between(v_c)/(2.0*PI);
  beta = 360.0*v_a.angle_between(v_c)/(2.0*PI);
  gamma = 360.0*v_a.angle_between(v_b)/(2.0*PI);
  initMatrices();
  distanceCalculator = MIN_PER_DISTANCE(v_a.x, v_b.x, v_b.y, v_c.x, v_c.y, v_c.z);
}

/** Print the information about this network of atoms to the
 *  provided output stream, including the information about each
 *  atom in the network. */
void ATOM_NETWORK::print(ostream &out){
    out << "Name: " << name << "\n"
    << "A: " << a <<"     B: " << b << "     C: " << c << "\n"
    << "Alpha: " << alpha <<"     Beta: " << beta << "     Gamma: " << gamma
    << "\n";
    out << "v_a: "; v_a.print();
    out << "v_b: "; v_b.print();
    out << "v_c: "; v_c.print();
    out << "Number of atoms: " << numAtoms << "\n";
    for(int i = 0; i<numAtoms; i++)
        atoms.at(i).print();
}

/** Copy the data contained in this ATOM_NETWORK to a new network using
 the provided pointer.
 */
void ATOM_NETWORK::copy(ATOM_NETWORK *newNet){
    newNet->a = a; newNet->b = b; newNet->c = c;
    newNet->alpha = alpha; newNet->beta = beta; newNet->gamma = gamma;
    newNet->v_a = v_a; newNet->v_b = v_b; newNet->v_c = v_c;
    newNet->numAtoms = numAtoms; newNet->name = name;
    newNet->atoms.clear(); newNet->atoms = atoms;
    newNet->IDmapping.clear(); newNet->IDmapping = IDmapping;

    newNet->vertices.clear(); newNet->vertices = vertices;
    newNet->orphan_edge_starts.clear(); newNet->orphan_edge_starts = orphan_edge_starts;
    newNet->orphan_edge_ends.clear(); newNet->orphan_edge_ends = orphan_edge_ends;
    newNet->vertex_basic_indices.clear(); newNet->vertex_basic_indices = vertex_basic_indices;
    newNet->vertex_symmetry_operators.clear(); newNet->vertex_symmetry_operators = vertex_symmetry_operators;
    newNet->sym_ID = sym_ID;
    newNet->sym_name = sym_name;
    newNet->valid = valid;

    newNet->initialize();
}

/** Calculate the unit cell vectors based on the provided values
 of its side lengths and angles. v_a corresponds to the cartesian
 x-axis.
 */
void ATOM_NETWORK::initialize(){
  double tempd, talpha, tbeta, tgamma;
  talpha = 2*PI/360.0*alpha;
  tbeta  = 2*PI/360.0*beta;
  tgamma = 2*PI/360.0*gamma;
  double cosb = cos(tbeta);
  double cosg = cos(tgamma);
  double sing = sin(tgamma);
  tempd=(cos(talpha)-cosg*cosb)/sing;
  v_a.x=a;
  v_a.y=0.0;
  v_a.z=0.0;  
  v_b.x=b*cosg; 
  if(fabs(v_b.x) < TOLERANCE)
    v_b.x=0;  
  v_b.y=b*sing;
  v_b.z=0.0;  
  v_c.x=c*cosb;
  if(fabs(v_c.x) < TOLERANCE)
    v_c.x = 0;
  v_c.y=c*tempd;
  if(fabs(v_c.y) < TOLERANCE)
    v_c.y = 0;
  v_c.z=c*sqrt(1.0-(cosb*cosb)-(tempd*tempd));
  initMatrices();
  distanceCalculator = MIN_PER_DISTANCE(v_a.x, v_b.x, v_b.y, v_c.x, v_c.y, v_c.z);
}

/* Store the initialized unit cell vectors in matrix form. */
void ATOM_NETWORK::initMatrices(){
    ucVectors[0][0] = v_a.x; ucVectors[1][0] = v_a.y; ucVectors[2][0] = v_a.z;
    ucVectors[0][1] = v_b.x; ucVectors[1][1] = v_b.y; ucVectors[2][1] = v_b.z;
    ucVectors[0][2] = v_c.x; ucVectors[1][2] = v_c.y; ucVectors[2][2] = v_c.z;
    valid = tryInvertMatrix(ucVectors,invUCVectors);
}

MIN_PER_DISTANCE ATOM_NETWORK::getDistCalc() const {
    return distanceCalculator;
}

/** Determine the smallest supercell dimensions such that a sphere of a given
 *  diameter does not overlap with itself across the periodic boundary */
TRIPLET ATOM_NETWORK::getSmallestSupercell(double diam) {
    //1) set lower bound based on lengths of cell axes
    int na = 1+ (diam/a);
    int nb = 1+ (diam/b);
    int nc = 1+ (diam/c);
    
    int fewest_cells = -1; //no supercell yet successful
    TRIPLET smallest_supercell(-1,-1,-1);
    
    //2) search all possible supercell sizes to find the smallest one satisfying
    // the minimum image convention for this radius
    TRIPLET lb(na, nb, nc);
    vector<TRIPLET> supercells;
    supercells.push_back(lb);
    while(supercells.size()>0) {
        //3) take the oldest candidate on the vector
        TRIPLET s = supercells.at(0); //oldest supercell candidate, s
        for(int i=0; i<supercells.size()-1; i++)
            supercells.at(i) = supercells.at(i+1); //shift vector up
        supercells.pop_back(); //delete last, which is now duplicated at last-1
        int num_cells = s.x*s.y*s.z;
        //4) is s a potential new smallest supercell?
        if(num_cells<fewest_cells || fewest_cells<0) {
            //5) does s satisfy the minimum image convention?
            int status = check_sphere_overlap(s.x, s.y, s.z, diam, this); //For time being
            if(status==-1) {
                printf("WARNING: bad unit cell angles!\n");
                return smallest_supercell;
            } else if(status==1) { //acceptable!
                fewest_cells = num_cells;
                smallest_supercell = s;
                //printf("smallest satisfactory supercell so far: (%d %d %d) =
                //%d cells\n", s.a, s.b, s.c, num_cells);
            } else { //unacceptable - try larger supercells in each direction
                TRIPLET s2(s.x+1, s.y, s.z);
                TRIPLET s3(s.x, s.y+1, s.z);
                TRIPLET s4(s.x, s.y, s.z+1);
                supercells.push_back(s2);
                supercells.push_back(s3);
                supercells.push_back(s4);
            }
        }
    }
    return smallest_supercell;
}

// abc_to_xyz and xyz_to_abc functions which need ATOM_NETWORK are moved from geometry
/** Convert coordinates relative to the unit cell (a,b,c) to absolute cartesian
 coordinates (x,y,z). */
const Point ATOM_NETWORK::abc_to_xyz(double a, double b, double c) const {
    //Point xyzCoords;
    /*
     xyzCoords.x = a*cell->v_a.x+b*cell->v_b.x+c*cell->v_c.x;
     xyzCoords.y = a*cell->v_a.y+b*cell->v_b.y+c*cell->v_c.y;
     xyzCoords.z = a*cell->v_a.z+b*cell->v_b.z+c*cell->v_c.z;
     if(fabs(xyzCoords.x) < 0.00001) xyzCoords.x = 0;
     if(fabs(xyzCoords.y) < 0.00001) xyzCoords.y = 0;
     if(fabs(xyzCoords.z) < 0.00001) xyzCoords.z = 0;
     */
    
    // Use only non-zero elements in computation
    double xt = a*v_a.x+b*v_b.x+c*v_c.x;
    double yt = b*v_b.y+c*v_c.y;
    double zt = c*v_c.z;
    return Point(xt, yt, zt);
}
const XYZ ATOM_NETWORK::abc_to_xyz_returning_XYZ(double a, double b, double c) const {
  Point p = abc_to_xyz(a, b, c);
  return XYZ(p[0], p[1], p[2]);
}

/** Convert coordinates relative to the unit cell (a,b,c) to absolute cartesian
 coordinates (x,y,z). */
const Point ATOM_NETWORK::abc_to_xyz(Point abcPt) const {
    return abc_to_xyz(abcPt[0], abcPt[1], abcPt[2]);
}

const Point ATOM_NETWORK::abc_to_xyz (const XYZ& temp) const {
    return abc_to_xyz(temp.x,temp.y,temp.z);
}
const XYZ ATOM_NETWORK::abc_to_xyz_returning_XYZ(const XYZ& temp) const {
  return abc_to_xyz_returning_XYZ(temp.x,temp.y,temp.z);
}

/** Convert coordinates relative to the cartesian axes to those
 relative to the unit cell vectors. Returns an Point instance of the
 form (a_coord, b_coord, c_coord). */
const Point ATOM_NETWORK::xyz_to_abc(double xi, double yi, double zi) const{
    //Point abcCoords;
    /*
     abcCoords.x = x*cell->invUCVectors[0][0]+y*cell->invUCVectors[0][1]+z*cell->invUCVectors[0][2];
     abcCoords.y = x*cell->invUCVectors[1][0]+y*cell->invUCVectors[1][1]+z*cell->invUCVectors[1][2];
     abcCoords.z = x*cell->invUCVectors[2][0]+y*cell->invUCVectors[2][1]+z*cell->invUCVectors[2][2];
     if(fabs(abcCoords.x) < 0.0001) abcCoords.x = 0;
     if(fabs(abcCoords.y) < 0.0001) abcCoords.y = 0;
     if(fabs(abcCoords.z) < 0.0001) abcCoords.z = 0;
     */
    
    // Use only non-zero elements in computation
    double xt = xi*invUCVectors[0][0]+yi*invUCVectors[0][1]+zi*invUCVectors[0][2];
    double yt = yi*invUCVectors[1][1]+zi*invUCVectors[1][2];
    double zt = zi*invUCVectors[2][2];
    return Point(xt, yt, zt);
}
const XYZ ATOM_NETWORK::xyz_to_abc_returning_XYZ(double a, double b, double c) const {
  Point p = xyz_to_abc(a, b, c);
  return XYZ(p[0], p[1], p[2]);
}

/* calculates the unit edge length in a net */
double get_unit_edge_length(ATOM_NETWORK *c) {
  int num_v = c->vertices.size();
  double unit_edge_length = -1;
  for(int i=0; i<num_v; i++) {
    VERTEX v = c->vertices.at(i);
    int num_e = v.edges.size();
    XYZ v_xyz = c->abc_to_xyz_returning_XYZ(v.abc);
    for(int j=0; j<num_e; j++) {
      XYZ e_xyz = c->abc_to_xyz_returning_XYZ(v.edges.at(j));
      double edge_length = (e_xyz-v_xyz).magnitude();
      if(unit_edge_length<0) unit_edge_length = edge_length;
      else if(fabs(unit_edge_length-edge_length)>DISTANCE_TOLERANCE) {
        printf("ERROR: found a basic edge length of %.3f which is sufficiently different to the previous length of %.3f; at the moment, nets with more than one edge length are not handled\n", edge_length, unit_edge_length);
        exit(EXIT_FAILURE);
      }
    }
  }
  return unit_edge_length;
}

/** Convert coordinates relative to the cartesian axes to those
 relative to the unit cell vectors. Returns an Point instance of the
 form (a_coord, b_coord, c_coord). */
const Point ATOM_NETWORK::xyz_to_abc (Point xyzPt) const {
    return xyz_to_abc(xyzPt[0], xyzPt[1], xyzPt[2]);
}

const Point ATOM_NETWORK::xyz_to_abc(const XYZ& temp) const {
    return xyz_to_abc(temp.x,temp.y,temp.z);
}
const XYZ ATOM_NETWORK::xyz_to_abc_returning_XYZ(const XYZ& temp) const {
  return xyz_to_abc_returning_XYZ(temp.x,temp.y,temp.z);
}

/** Calculates the minimum distance between the two points whose coordinates are relative to the unit cell vectors.    */
double ATOM_NETWORK::calcDistanceABC(double a1, double b1, double c1, double a2, double b2, double c2) const{
    return getDistCalc().minimum_periodic_distance(a1, b1, c1, a2, b2, c2);
}
double ATOM_NETWORK::calcDistanceABC(XYZ a, XYZ b) {
  return getDistCalc().minimum_periodic_distance(a.x, a.y, a.z, b.x, b.y, b.z);
}

/** Calculates the minimum distance between the point (x,y,z) and the other point whose coordinates
 *  are relative to the unit cell vectors. */
double ATOM_NETWORK::calcDistanceXYZABC(double x1, double y1, double z1, double a2, double b2, double c2) const {
    Point abcCoord = xyz_to_abc(x1, y1, z1);
    return calcDistanceABC(abcCoord[0], abcCoord[1], abcCoord[2], a2, b2, c2);
}

/** Calculates the minimum distance between the point (x,y,z) and the provided atom.  */
double ATOM_NETWORK::calcDistance(double x, double y, double z, ATOM *atm) const {
    return calcDistanceXYZABC(x, y, z, atm->a_coord, atm->b_coord, atm->c_coord);
}

/** Calculates the minimum distance between the two provided atoms.  */
double ATOM_NETWORK::calcDistance(const ATOM& atm1, const ATOM& atm2) const {
    return calcDistanceABC(atm1.a_coord, atm1.b_coord, atm1.c_coord, atm2.a_coord, atm2.b_coord, atm2.c_coord);
}

/** Calculates the minimum distance between the points (x1,y1,z1) and (x2,y2,z2).  */
double ATOM_NETWORK::calcDistanceXYZ(double x1, double y1, double z1, double x2, double y2, double z2) const {
    Point abcCoord = xyz_to_abc(x1, y1, z1);
    //cout << abcCoord.vals[0] <<  abcCoord.vals[1] << abcCoord.vals[2] << endl;
    //cout << x2 <<  y2 << z2 << endl;
    //cout << calcDistanceXYZABC(x2, y2, z2, abcCoord[0], abcCoord[1], abcCoord[2]);
    return calcDistanceXYZABC(x2, y2, z2, abcCoord[0], abcCoord[1], abcCoord[2]);
}

/** Rich edit: for static point (x1,y1,z1), returns the closest periodic image of point (x2,y2,z2).*/
const XYZ ATOM_NETWORK::getClosestPointInABC(double x1, double y1, double z1, double x2, double y2, double z2){
    Point abcCoordStatic = xyz_to_abc(x1, y1, z1);
    Point abcCoordMobile = xyz_to_abc(x2, y2, z2);
    XYZ answer;
    getDistCalc().closest_periodic_image(abcCoordStatic[0], abcCoordStatic[1], abcCoordStatic[2],
                                         abcCoordMobile[0], abcCoordMobile[1], abcCoordMobile[2],
                                         answer.x, answer.y, answer.z);
    return answer;
}

/** Modify the provided (x,y,z) Point so that its coordinates reflect unit cell translations
 *  by the provided amounts along each unit cell axis. */
void ATOM_NETWORK::translatePoint(Point *origPoint, double da, double db, double dc){
    (*origPoint)[0] = (*origPoint)[0] +  da*v_a.x + db*v_b.x + dc*v_c.x;
    (*origPoint)[1] = (*origPoint)[1] +  da*v_a.y + db*v_b.y + dc*v_c.y;
    (*origPoint)[2] = (*origPoint)[2] +  da*v_a.z + db*v_b.z + dc*v_c.z;
}

/** Shifts the provided Point whose coordinates are relative to the unit cell vectors
 *  such that it lies within the unitcell.  */
const Point ATOM_NETWORK::shiftABCInUC(Point abcCoords) {
    return Point(trans_to_origuc(abcCoords[0]), trans_to_origuc(abcCoords[1]), trans_to_origuc(abcCoords[2]));
}

/** Shifts the provided Point whose coordinates are relative to the x,y,z vectors
 *  such that it lies within the unitcell   */
const Point ATOM_NETWORK::shiftXYZInUC(Point xyzCoords) {
    Point abcCoords = shiftABCInUC(xyz_to_abc(xyzCoords));
    return abc_to_xyz(abcCoords[0], abcCoords[1], abcCoords[2]);
}

/** Shift the coordinates of the provided Point using the unit cell vectors until
 *  the Euclidean distance between the Point and (x,y,z) is minimal. Returns the resulting point.
 */
const Point ATOM_NETWORK::minimizePointDistance(Point origPoint, double dx, double dy, double dz){
    Point abc_one  = xyz_to_abc(origPoint);
    Point abc_two  = xyz_to_abc(dx, dy, dz);
    
    double minDa = DBL_MAX, minDb = DBL_MAX, minDc = DBL_MAX, best_a = DBL_MAX, best_b = DBL_MAX, best_c = DBL_MAX;
    getDistCalc().closest_periodic_image(abc_two[0], abc_two[1], abc_two[2], abc_one[0], abc_one[1], abc_one[2],
                                         minDa, minDb, minDc, best_a, best_b, best_c);
    return abc_to_xyz(best_a, best_b, best_c);
}

/* Identifies tetrahedra of the given atom type, calculates their tetrahedrality and returns as a vector of doubles */
vector<double> ATOM_NETWORK::find_tetrahedra(string element) {
    vector<double> tetras;
    double d_ij = -1, d_ik = -1, d_il = -1, d_jk = -1, d_jl = -1, d_kl = -1; //stores distances between each pair of atoms in a potential tetrahedral arrangement
    double maxDist = 5, minDist = 0.1; //four atoms are considered to be in a tetrahedral arrangement if they are each within this distance range of each other - if this is true, we then calculate the index of tetrahedral distortion
    for(int i=0; i<numAtoms; i++) {
        if(atoms.at(i).type.compare(element)==0) {
            for(int j=i+1; j<numAtoms; j++) {
                if(atoms.at(j).type.compare(element)==0) {
                    d_ij = getDistCalc().minimum_periodic_distance(atoms.at(i).a_coord, atoms.at(i).b_coord, atoms.at(i).c_coord, atoms.at(j).a_coord, atoms.at(j).b_coord, atoms.at(j).c_coord);
                    if(d_ij>minDist && d_ij<maxDist) {
                        for(int k=j+1; k<numAtoms; k++) {
                            if(atoms.at(k).type.compare(element)==0) {
                                d_ik = getDistCalc().minimum_periodic_distance(atoms.at(i).a_coord, atoms.at(i).b_coord, atoms.at(i).c_coord, atoms.at(k).a_coord, atoms.at(k).b_coord, atoms.at(k).c_coord);
                                if(d_ik>minDist && d_ik<maxDist) {
                                    d_jk = getDistCalc().minimum_periodic_distance(atoms.at(j).a_coord, atoms.at(j).b_coord, atoms.at(j).c_coord, atoms.at(k).a_coord, atoms.at(k).b_coord, atoms.at(k).c_coord);
                                    if(d_jk>minDist && d_jk<maxDist) {
                                        for(int l=k+1; l<numAtoms; l++) {
                                            if(atoms.at(l).type.compare(element)==0) {
                                                d_il = getDistCalc().minimum_periodic_distance(atoms.at(i).a_coord, atoms.at(i).b_coord, atoms.at(i).c_coord, atoms.at(l).a_coord, atoms.at(l).b_coord, atoms.at(l).c_coord);
                                                if(d_il>minDist && d_il<maxDist) {
                                                    d_jl = getDistCalc().minimum_periodic_distance(atoms.at(j).a_coord, atoms.at(j).b_coord, atoms.at(j).c_coord, atoms.at(l).a_coord, atoms.at(l).b_coord, atoms.at(l).c_coord);
                                                    if(d_jl>minDist && d_jl<maxDist) {
                                                        d_kl = getDistCalc().minimum_periodic_distance(atoms.at(k).a_coord, atoms.at(k).b_coord, atoms.at(k).c_coord, atoms.at(l).a_coord, atoms.at(l).b_coord, atoms.at(l).c_coord);
                                                        if(d_kl>minDist && d_kl<maxDist) {
                                                            //we have found four atoms within the tolerance
                                                            double t = CalculateTetrahedrality4Atoms(atoms.at(i), atoms.at(j), atoms.at(k), atoms.at(l));
                                                            tetras.push_back(t);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    //sort before returning
    sort(tetras.begin(), tetras.end());
    return tetras;
}

/* Returns Tetrahedrality index for a tetrahedron defined by four atoms */
double ATOM_NETWORK::CalculateTetrahedrality4Atoms(const ATOM& atm1, const ATOM& atm2, const ATOM& atm3, const ATOM& atm4) const {
    
    vector <double> edges; // tetrahedra edges' lengths
    double d_mean=0.0; // mean value
    double Tindex=0.0;
    
    edges.push_back(calcDistance(atm1,atm2));
    edges.push_back(calcDistance(atm1,atm3));
    edges.push_back(calcDistance(atm1,atm4));
    edges.push_back(calcDistance(atm2,atm3));
    edges.push_back(calcDistance(atm2,atm4));
    edges.push_back(calcDistance(atm3,atm4));
    
    for(int i=0;i<6;i++) d_mean+=edges[i];
    
    d_mean=d_mean/6.0;
    
    for(int i=0;i<5;i++) {
        for(int j=i+1;j<6;j++) {
            Tindex+=((edges[i]-edges[j])*(edges[i]-edges[j]))/(15.0*d_mean*d_mean);
        }
    }
    
    return Tindex;
}


/* Return chemical formula */
std::string ATOM_NETWORK::returnChemicalFormula(){

 std::string formula;

 for(int i=0; i<MAX_ATOMIC_NUMBER; i++) atomic_composition[i] = 0;

 for(int i=0; i<numAtoms; i++) {
    atomic_composition[lookupAtomicNumber(atoms.at(i).type)-1]+=1;
    };

 for(int i=0; i<numAtoms; i++) {

    if(atomic_composition[lookupAtomicNumber(atoms.at(i).type)-1]>0)
      {
      std::ostringstream os ;
      os << atomic_composition[lookupAtomicNumber(atoms.at(i).type)-1];
      formula += atoms.at(i).type + os.str();
      atomic_composition[lookupAtomicNumber(atoms.at(i).type)-1]=0;
      };

    };

 return formula;

} // end returnChemicalFormula

/* Determine the smallest supercell dimensions such that a sphere of a given
 * diameter does not overlap with itself across the periodic boundary */
/* Marked for deletion
 TRIPLET getSmallestSupercell(double diam, ATOM_NETWORK *atmnet) {
 //1) set lower bound based on lengths of cell axes
 int na = 1+ (diam/atmnet->a);
 int nb = 1+ (diam/atmnet->b);
 int nc = 1+ (diam/atmnet->c);
 
 int fewest_cells = -1; //no supercell yet successful
 TRIPLET smallest_supercell(-1,-1,-1);
 
 //2) search all possible supercell sizes to find the smallest one satisfying
 // the minimum image convention for this radius
 TRIPLET lb(na, nb, nc);
 vector<TRIPLET> supercells;
 supercells.push_back(lb);
 while(supercells.size()>0) {
 //3) take the oldest candidate on the vector
 TRIPLET s = supercells.at(0); //oldest supercell candidate, s
 for(int i=0; i<supercells.size()-1; i++)
 supercells.at(i) = supercells.at(i+1); //shift vector up
 supercells.pop_back(); //delete last, which is now duplicated at last-1
 int num_cells = s.x*s.y*s.z;
 //4) is s a potential new smallest supercell?
 if(num_cells<fewest_cells || fewest_cells<0) {
 //5) does s satisfy the minimum image convention?
 int status = check_sphere_overlap(s.x, s.y, s.z, diam, atmnet);
 if(status==-1) {
 printf("WARNING: bad unit cell angles!\n");
 return smallest_supercell;
 } else if(status==1) { //acceptable!
 fewest_cells = num_cells;
 smallest_supercell = s;
 //printf("smallest satisfactory supercell so far: (%d %d %d) =
 //%d cells\n", s.a, s.b, s.c, num_cells);
 } else { //unacceptable - try larger supercells in each direction
 TRIPLET s2(s.x+1, s.y, s.z);
 TRIPLET s3(s.x, s.y+1, s.z);
 TRIPLET s4(s.x, s.y, s.z+1);
 supercells.push_back(s2);
 supercells.push_back(s3);
 supercells.push_back(s4);
 }
 }
 }
 return smallest_supercell;
 }
 * End: Marked for deletion
 */

/// Default constructor
VOR_EDGE::VOR_EDGE(){}

VOR_EDGE::VOR_EDGE(int myFrom, int myTo, double rad,
                   int dx, int dy, int dz, double len):
    from(myFrom), to(myTo), rad_moving_sphere(rad),
    delta_uc_x(dx), delta_uc_y(dy), delta_uc_z(dz),
    length(len) {}

/// Copy constructor
VOR_EDGE::VOR_EDGE(const VOR_EDGE& orig):
    from(orig.from), to(orig.to), rad_moving_sphere(orig.rad_moving_sphere),
    delta_uc_x(orig.delta_uc_x), delta_uc_y(orig.delta_uc_y),
    delta_uc_z(orig.delta_uc_z) {}

//Default constructor for VOR_NODE
VOR_NODE::VOR_NODE(){}

VOR_NODE::VOR_NODE(double myX, double myY, double myZ,
                   double rad, vector<int> ids){
    x = myX; y = myY; z = myZ;
    rad_stat_sphere = rad;
    atomIDs = ids;
}

//Copy constructor for VOR_NODE
/*VOR_NODE::VOR_NODE(const VOR_NODE& orig):
    x(orig.x), y(orig.y), z(orig.z),
    rad_stat_sphere(orig.rad_stat_sphere), atomIDs(orig.atomIDs) {}
    */

/// Voronoi network constructor
VORONOI_NETWORK::VORONOI_NETWORK (){}
VORONOI_NETWORK::VORONOI_NETWORK(const XYZ& inp_va, const XYZ& inp_vb, const XYZ& inp_vc,
                                 const vector<VOR_NODE>& inp_nodes, 
                                 const vector<VOR_EDGE>& inp_edges): 
    edges(inp_edges), nodes(inp_nodes), v_a(inp_va), v_b(inp_vb), v_c(inp_vc){}

/// Copy constructor for VORONOI_NETWORK
VORONOI_NETWORK::VORONOI_NETWORK (const VORONOI_NETWORK& net):
    v_a(net.v_a), v_b(net.v_b), v_c(net.v_c), edges(net.edges), nodes(net.nodes) {}

/** Copy the data contained in this VORONOI_NETWORK to a new network using
 the provided pointer. Deprecated, use copy constructor.
 */
void VORONOI_NETWORK::copy(VORONOI_NETWORK *newNet){
    newNet->v_a = v_a; newNet->v_b = v_b; newNet->v_c = v_c;
    newNet->edges.clear(); newNet->edges = edges;
    newNet->nodes.clear(); newNet->nodes = nodes;
}


/* Marked for deletion
 void VORONOI_NETWORK::filterVornetEdges(vector<int> nodeIDs,
 VORONOI_NETWORK *oldNet,
 VORONOI_NETWORK *newNet)
 {
 vector<bool> includeNodes = vector<bool>(oldNet->nodes.size(), false);
 for(unsigned int i = 0; i < nodeIDs.size(); i++){
 includeNodes[nodeIDs[i]] = true;
 }
 
 vector<VOR_NODE> newNodes = vector<VOR_NODE>();
 for(unsigned int i = 0; i < oldNet->nodes.size(); i++){
 newNodes.push_back(oldNet->nodes[i]);
 }
 
 vector<VOR_EDGE> newEdges = vector<VOR_EDGE>();
 
 for(unsigned int i = 0; i < oldNet->edges.size(); i++){
 VOR_EDGE edge = oldNet->edges[i];
 if(includeNodes[edge.from] && includeNodes[edge.to]){
 newEdges.push_back(edge);
 }
 }
 
 newNet->nodes = newNodes;
 newNet->edges = newEdges;
 newNet->v_a = oldNet->v_a;
 newNet->v_b = oldNet->v_b;
 newNet->v_c = oldNet->v_c;
 }
 * End: marked for deletion
 */

/* Removes all the edges between nodes that are not both contained in nodeIDs.
 However, these nodes remain in the Voronoi network
 */
const VORONOI_NETWORK VORONOI_NETWORK::filterEdges(vector<int> nodeIDs)
{
    vector<bool> includeNodes = vector<bool>(nodes.size(), false);
    for(unsigned int i = 0; i < nodeIDs.size(); i++){
        includeNodes[nodeIDs[i]] = true;
    }
    
    vector<VOR_NODE> newNodes = vector<VOR_NODE>();
    for(unsigned int i = 0; i < nodes.size(); i++){
        newNodes.push_back(nodes[i]);
    }
    
    vector<VOR_EDGE> newEdges = vector<VOR_EDGE>();
    for(unsigned int i = 0; i < edges.size(); i++){
        VOR_EDGE edge = edges[i];
        if(includeNodes[edge.from] && includeNodes[edge.to]){
            newEdges.push_back(edge);
        }
    }
    
    return VORONOI_NETWORK(v_a, v_b, v_c, newNodes, newEdges);
}

VORONOI_NETWORK filterVoronoiNetwork(const VORONOI_NETWORK* vornet, double minRadius)
{
    // Add all nodes whose radius is greater than the provided minimum to a list
    map<int,int> idMappings;
    vector<VOR_NODE> newNodes;
    
    int i = 0;
    int idCount = 0;
    vector<VOR_NODE>::const_iterator nodeIter = vornet->nodes.begin();
    for(;nodeIter != vornet->nodes.end(); i++, nodeIter++)
        if(nodeIter->rad_stat_sphere > minRadius){
            newNodes.push_back(*nodeIter);
            idMappings.insert(pair<int,int> (i, idCount));
            idCount++;
        }
    
    // Add all edges whose radius is greater than the provided minimum to a list
    vector<VOR_EDGE> newEdges;
    
    vector<VOR_EDGE>::const_iterator edgeIter = vornet->edges.begin();
    for(;edgeIter != vornet->edges.end(); edgeIter++)
        if((edgeIter->rad_moving_sphere > minRadius) && 
                    (idMappings.find(edgeIter->from) != idMappings.end()) && 
                    (idMappings.find(edgeIter->to) != idMappings.end())){
            int from = idMappings.find(edgeIter->from)->second;
            int to = idMappings.find(edgeIter->to)->second;
            newEdges.push_back(
                               VOR_EDGE(from, to, edgeIter->rad_moving_sphere,
                                        edgeIter->delta_uc_x, edgeIter->delta_uc_y,
                                        edgeIter->delta_uc_z, edgeIter->length)
                              );
        }
    
    return VORONOI_NETWORK(vornet->v_a, vornet->v_b, vornet->v_c, newNodes, newEdges);
    
}

/** Copies all edges and nodes within the provided VORONOI_NETWORK
 to a new network iff a sphere with the specified radius can pass.*/
/* Deprecated. Use the function above. */
void filterVoronoiNetwork(VORONOI_NETWORK *vornet, VORONOI_NETWORK *newVornet, double minRadius){
    map<int,int> idMappings;
    vector<VOR_NODE>::iterator nodeIter = vornet->nodes.begin();
    vector<VOR_NODE> newNodes;
    int i = 0;
    int idCount = 0;
    
    // Add all nodes whose radius is greater than the provided minimum to a list
    while(nodeIter != vornet->nodes.end()){
        if(nodeIter->rad_stat_sphere > minRadius){
            newNodes.push_back(*nodeIter);
            idMappings.insert(pair<int,int> (i, idCount));
            idCount++;
        }
        i++;
        nodeIter++;
    }
    
    // Copy nodes that met requirement
    newVornet->nodes = newNodes;
    
    vector<VOR_EDGE>::iterator edgeIter = vornet->edges.begin();
    vector<VOR_EDGE> newEdges;
    
    // Add all edges whose radius is greater than the provided minimum to a list
    while(edgeIter != vornet->edges.end()){
        if((edgeIter->rad_moving_sphere > minRadius) && (
                                                         idMappings.find(edgeIter->from) != idMappings.end()) && (
                                                                                                                  idMappings.find(edgeIter->to) != idMappings.end())){
            VOR_EDGE newEdge;
            newEdge.from = idMappings.find(edgeIter->from)->second;
            newEdge.to   = idMappings.find(edgeIter->to)->second;
            newEdge.rad_moving_sphere = edgeIter->rad_moving_sphere;
            newEdge.delta_uc_x = edgeIter->delta_uc_x;
            newEdge.delta_uc_y = edgeIter->delta_uc_y;
            newEdge.delta_uc_z = edgeIter->delta_uc_z;
            newEdge.length = edgeIter->length;
            newEdges.push_back(newEdge);
        }
        edgeIter++;
    }
    
    // Copy edges that met requirement
    newVornet->edges = newEdges;
    
    // Copy unitcell vectors to new network
    newVornet->v_a = vornet->v_a;
    newVornet->v_b = vornet->v_b;
    newVornet->v_c = vornet->v_c;
}


/** Returns a copy of the VORNOI_NETWORK instance
 but removes the edges that do not allow a sphere
 with the provided radius to pass. */
const VORONOI_NETWORK VORONOI_NETWORK::prune(const double& minRadius)
{
    vector<VOR_EDGE> newEdges;
    // Add edges whose radius is greater than the input minimum to a list
    for(vector<VOR_EDGE>::iterator edgeIter = edges.begin();
        edgeIter != edges.end(); edgeIter++) {
        if(edgeIter->rad_moving_sphere > minRadius){
            //VOR_EDGE newEdge = *edgeIter;
            
            // further check: keep edges only if they connect accessible nodes
            if( nodes[edgeIter->from].rad_stat_sphere > minRadius && nodes[edgeIter->to].rad_stat_sphere > minRadius )
                newEdges.push_back(*edgeIter);
            
        }
    };
    
    vector<VOR_NODE> newNodes = nodes;
    for(unsigned int i = 0; i < nodes.size(); i++)
    {
        if(nodes[i].rad_stat_sphere > minRadius) newNodes[i].active = true; else newNodes[i].active = false;
    };
    
    return VORONOI_NETWORK(v_a, v_b, v_c, newNodes, newEdges);
}


/** Stores a copy of the original VORNOI_NETWORK into the other provided
 VORONOI_NETWORK but removes the edges that do not allow a sphere
 with the provided radius to pass. */
/* Marked for deletion
 void pruneVoronoiNetwork(VORONOI_NETWORK *vornet,
 VORONOI_NETWORK *newVornet,
 double minRadius){
 newVornet->nodes = vornet->nodes;
 
 vector<VOR_EDGE>::iterator edgeIter = vornet->edges.begin();
 vector<VOR_EDGE> newEdges;
 
 // Add all edges whose radius is greater than the provided minimum to a list
 while(edgeIter != vornet->edges.end()){
 if(edgeIter->rad_moving_sphere > minRadius){
 VOR_EDGE newEdge;
 newEdge.from = edgeIter->from;
 newEdge.to   = edgeIter->to;
 newEdge.rad_moving_sphere = edgeIter->rad_moving_sphere;
 newEdge.delta_uc_x = edgeIter->delta_uc_x;
 newEdge.delta_uc_y = edgeIter->delta_uc_y;
 newEdge.delta_uc_z = edgeIter->delta_uc_z;
 newEdge.length = edgeIter->length;
 newEdges.push_back(newEdge);
 }
 edgeIter++;
 }
 
 // Copy edges that met requirement
 newVornet->edges = newEdges;
 
 // Copy unitcell vectors to new network
 newVornet->v_a = vornet->v_a;
 newVornet->v_b = vornet->v_b;
 newVornet->v_c = vornet->v_c;
 }
 * End: Marked for deletion
 */

/** Stores a copy of the original VORNOI_NETWORK into the other provided
 *  VORONOI_NETWORK but removes the edges that are connected to specified
 *  nodes (specified by ID list)
 */
void pruneVoronoiNetworkfromEdgeList(VORONOI_NETWORK *vornet,
                                     VORONOI_NETWORK *newVornet,
                                     vector <int> ids){
    
    newVornet->nodes = vornet->nodes;
    
    vector<VOR_EDGE>::iterator edgeIter = vornet->edges.begin();
    vector<VOR_EDGE> newEdges;
    
    // Add all edges whose are not connected to specified nodes
    while(edgeIter != vornet->edges.end()){
        int flag=0;
        for(unsigned int i=0; i<ids.size(); i++){
            // if edge connects to node of specified id, flag it
            if(edgeIter->from==ids[i]||edgeIter->to==ids[i]) flag++;
        };
        
        if(flag==0){ // only keep unflagged edges
            VOR_EDGE newEdge;
            newEdge.from = edgeIter->from;
            newEdge.to   = edgeIter->to;
            newEdge.rad_moving_sphere = edgeIter->rad_moving_sphere;
            newEdge.delta_uc_x = edgeIter->delta_uc_x;
            newEdge.delta_uc_y = edgeIter->delta_uc_y;
            newEdge.delta_uc_z = edgeIter->delta_uc_z;
            newEdge.length = edgeIter->length;
            newEdges.push_back(newEdge);
        }
        edgeIter++;
    }
    
    
    newVornet->edges = newEdges;
    newVornet->v_a = vornet->v_a;
    newVornet->v_b = vornet->v_b;
    newVornet->v_c = vornet->v_c;
}


/* Attempt to substitute every other Si atom with an Al atom. ATOM_NETWORK may
 * only consist of Si and O atoms, where each Si atom must be bonded to exactly
 * 4 oxygen atoms and each oxygen atom must be bonded to exactly 2 Si atoms.
 * Returns true iff the substitution was successful and stores the number of
 * substitutions using the provided reference. The provided boolean specifies
 * whether the seeded silicon atom is substituted or not.
 * Since only 2 configurations are possible if the structure is consistent,
 * changing this parameter enables generation of all configurations. */
bool substituteAtoms(ATOM_NETWORK *origNet, ATOM_NETWORK *newNet,
                     bool substituteSeed, int *numSubstitutions, bool radial){
    int numAtoms = origNet->numAtoms;
    double max_bond_length = 1.95;
    
    vector< vector<int> > bonds = vector< vector<int> >(numAtoms, vector<int>());
    for(int i = 0; i < numAtoms; i++){
        ATOM atom_one = origNet->atoms[i];
        for(int j = i + 1; j < numAtoms; j++){
            ATOM atom_two = origNet->atoms[j];
            if(origNet->calcDistance(atom_one, atom_two) < max_bond_length){
                if(atom_one.type.compare(atom_two.type) == 0){
                    cerr << "Atomic network substitution aborted because atoms of same type "
                    "are bonded to one another" << "\n" << "Occurred for type "
                    << atom_one.type << " between atoms " << i << " and " << j << "\n";
                    return false;
                }
                bonds[i].push_back(j);
                bonds[j].push_back(i);
            }
        }
    }
    
    for(int i = 0; i < numAtoms; i++){
        if(origNet->atoms[i].type.compare("Si") == 0){
            if(bonds[i].size() != 4){
                cerr << "Atomic network substitution aborted because Si atom bonded to "
                << bonds[i].size() << " other atoms instead of 4" << "\n"
                << "Occurred for atom " << i << "\n";
                return false;
            }
        }
        else if(origNet->atoms[i].type.compare("O") == 0){
            if(bonds[i].size() != 2){
                cerr << "Atomic network substitution aborted because O atom bonded to "
                << bonds[i].size() << " other atoms instead of 2" << "\n"
                << "Occurred for atom " << i << "\n";
                return false;
            }
        }
        else{
            cerr << "Atomic network substitution aborted because atom type other than "
            "Si or O detected" << "\n" << "Occurred for atom " << i << "\n";
            return false;
        }
    }
    
    int firstSiID = -1;
    for(int i = 0; i < numAtoms; i++){
        if(origNet->atoms[i].type.compare("Si") == 0){
            firstSiID = i;
            break;
        }
    }
    
    if(firstSiID == -1){
        cerr << "Error: Atom substitution failed because structure does not "
        "contain any Si atoms \n";
        return false;
    }
    
    vector<bool> atomProc = vector<bool>(numAtoms, false); // Whether atom has been processed
    vector<bool> atomSubs = vector<bool>(numAtoms, false); // Whether atom has been substituted
    int numProc = 0; // Number of processed atoms
    
    vector< pair<int, bool> > atomsToProcess; // Stack of atom/substituion info
    vector<int> fromIDs; // Stack of where each atom/subst came from
    
    pair<int, bool> seed(firstSiID, substituteSeed);
    atomsToProcess.push_back(seed);  // Seed with the first Si atom and not substituting it
    fromIDs.push_back(-1); // Starting node
    
    while(atomsToProcess.size() != 0){
        pair<int, bool> atomInfo = atomsToProcess.back(); atomsToProcess.pop_back();
        int id = atomInfo.first; bool subst = atomInfo.second;
        int fromID = fromIDs.back(); fromIDs.pop_back();
        
        // Atom previously processed. Ensure that atom substitution agrees with previous results for Si atoms
        if(atomProc[id] ){
            if((origNet->atoms[id].type.compare("Si") == 0) && (atomSubs[id] != subst)){
                cerr << "Atomic network substitution aborted because every other Si atom "
                "substitution criteria failed" << "\n";
                return false;
            }
        }
        // First time atom is processed
        else {
            numProc++;
            // Mark atom as processed
            atomProc[id] = true;
            
            // If oxygen, should not be substituted but reverse substituion state for next set of Si atoms
            if(origNet->atoms[id].type.compare("O") == 0){
                atomSubs[id] = false;
                subst = !subst;
            }
            // If silicon, pass on substitution state to oxygen atoms, who will then reverse it
            // Record substitution status
            else if(origNet->atoms[id].type.compare("Si") == 0){
                atomSubs[id] = subst;
            }
            // Invalid type
            else {
                cerr << "Atomic network substitution aborted because atom type other than Si or O detected" << "\n"
                << "Occurred for atom " << id << "\n";
                return false;
            }
            
            // Added bonded atoms to stack with appropriate subsitution state
            for(unsigned int j = 0; j < bonds[id].size(); j++){
                if(bonds[id][j] != fromID){
                    atomsToProcess.push_back(pair<int, bool>(bonds[id][j], subst));
                    fromIDs.push_back(id);
                }
            }
        }
    }
    
    if(numProc != numAtoms){
        cerr << "Atom network substituion failed because not all atoms are interconnected" << "\n"
        << "Visited " << numProc << " out of " << numAtoms << " atoms " << "\n";
        return false;
    }
    
    // Perform final consistency check
    for(int i = 0; i < numAtoms; i++){
        // Just check pair of atoms bonded to each oxygen
        if(origNet->atoms[i].type.compare("O") == 0){
            if(atomSubs[bonds[i][0]] && atomSubs[bonds[i][1]]){
                cerr << "Atom network substitution failed final consistency check." << "\n"
                << "Oxygen atom #" << i << " is bound to atoms " <<  bonds[i][0] << " and " << bonds[i][1] << " which are of same type. " << "\n"
                << "Algorithmic issue likely." << "\n"
                << "Please contact the source code provider" << "\n";
                return false;
            }
        }
    }
    
    // Copy original information
    origNet->copy(newNet);
    
    // Change atoms that should be subsituted
    *numSubstitutions = 0;
    for(int i = 0; i < numAtoms; i++){
        if(atomSubs[i]){
            newNet->atoms[i].type   = "Al";
            newNet->atoms[i].radius = lookupRadius("Al", radial);
            (*numSubstitutions)++;
        }
    }
    
    // Successful substitution
    return true;
}


ptrdiff_t myrandom (ptrdiff_t i) { return rand();} // Random number generator function for fractional substitution
ptrdiff_t (*p_myrandom)(ptrdiff_t) = myrandom;     // Pointer to it

bool comparePairs(pair<int,int> p_one, pair<int,int> p_two){
    return p_one.first < p_two.first;
}



/* Attempt to substitute the specified fraction of Si atom with an Al atom. ATOM_NETWORK may only consist of
 *  Si and O atoms, where each Si atom must be bonded to exactly 4 oxygen atoms and each oxygen atom must
 *  be bonded to exactly 2 Si atoms. Returns true iff the substitution was successful and stores the number of substitutions using the
 *  provided reference. The provided boolean specifies whether the seeded silicon atom is substituted or not. The function works
 *  by first substituting every other Si atom and then reverting some of the substituted atoms back to Si. The provided
 *  random number generator seed is used to choose which atoms to switch back.*/
/* bool radial is always set to true in the new version (from 0.3 on) */
bool fracSubstituteAtoms(ATOM_NETWORK *origNet, ATOM_NETWORK *newNet,
                         bool substituteSeed, double frac, int randSeed,
                         int *numSubstitutions, double *fracSub, bool radial){
    int numAtoms = origNet->numAtoms;
    double max_bond_length = 1.95;
    
    if(frac > 0.5){
        cerr << "Fractional atomic network substitution failed because the fraction can not exceed 0.5" << "\n";
        return false;
    }
    if(frac < 0.0){
        cerr << "Fractional atomic network substitution failed because of invalid negative fraction " << frac << "\n";
        return false;
    }
    
    vector< vector<int> > bonds = vector< vector<int> >(numAtoms, vector<int>());
    for(int i = 0; i < numAtoms; i++){
        ATOM atom_one = origNet->atoms[i];
        for(int j = i + 1; j < numAtoms; j++){
            ATOM atom_two = origNet->atoms[j];
            if(origNet->calcDistance(atom_one, atom_two) < max_bond_length){
                if(atom_one.type.compare(atom_two.type) == 0){
                    cerr << "Fractional atomic network substitution aborted because atoms of same type are bonded to one another" << "\n"
                    << "Occurred for type " << atom_one.type << " between atoms " << i << " and " << j << "\n";
                    return false;
                }
                bonds[i].push_back(j);
                bonds[j].push_back(i);
            }
        }
    }
    
    int numSi = 0;
    
    for(int i = 0; i < numAtoms; i++){
        if(origNet->atoms[i].type.compare("Si") == 0){
            if(bonds[i].size() != 4){
                cerr << "Fractional atomic network substitution aborted because Si atom "
                "bonded to " << bonds[i].size() << " other atoms instead of 4" << "\n"
                << "Occurred for atom " << i << "\n";
                return false;
            }
            numSi++;
        }
        else if(origNet->atoms[i].type.compare("O") == 0){
            if(bonds[i].size() != 2){
                cerr << "Fractional atomic network substitution aborted because O atom "
                "bonded to " << bonds[i].size() << " other atoms instead of 2" << "\n"
                << "Occurred for atom " << i << "\n";
                return false;
            }
        }
        else{
            cerr << "Fractional atomic network substitution aborted because atom type "
            "other than Si or O detected" << "\n"
            << "Occurred for atom " << i << "\n";
            return false;
        }
    }
    
    int firstSiID = -1;
    for(int i = 0; i < numAtoms; i++){
        if(origNet->atoms[i].type.compare("Si") == 0){
            firstSiID = i;
            break;
        }
    }
    
    if(firstSiID == -1){
        cerr << "Error: Fractional atomic network substitution failed because "
        "structure does not contain any Si atoms \n";
        return false;
    }
    
    vector<bool> atomProc = vector<bool>(numAtoms, false); // Whether atom has been processed
    vector<bool> atomSubs = vector<bool>(numAtoms, false); // Whether atom has been substituted
    int numProc = 0; // Number of processed atoms
    
    vector< pair<int, bool> > atomsToProcess; // Stack of atom/substituion info
    vector<int> fromIDs; // Stack of where each atom/subst came from
    
    pair<int, bool> seed(firstSiID, substituteSeed);
    atomsToProcess.push_back(seed);  // Seed with the first Si atom and not substituting it
    fromIDs.push_back(-1); // Starting node
    
    while(atomsToProcess.size() != 0){
        pair<int, bool> atomInfo = atomsToProcess.back(); atomsToProcess.pop_back();
        int id = atomInfo.first; bool subst = atomInfo.second;
        int fromID = fromIDs.back(); fromIDs.pop_back();
        
        // Atom previously processed.
        // Ensure that atom substitution agrees with previous results for Si atoms
        if(atomProc[id] ){
            if((origNet->atoms[id].type.compare("Si") == 0) && (atomSubs[id] != subst)){
                cerr << "Fractional atomic network substitution aborted because every other "
                "Si atom substitution criteria failed" << "\n";
                return false;
            }
        }
        // First time atom is processed
        else {
            numProc++;
            // Mark atom as processed
            atomProc[id] = true;
            
            // If oxygen, should not be substituted but reverse substituion state for next set of Si atoms
            if(origNet->atoms[id].type.compare("O") == 0){
                atomSubs[id] = false;
                subst = !subst;
            }
            // If silicon, pass on substitution state to oxygen atoms, who will then reverse it
            // Record substitution status
            else if(origNet->atoms[id].type.compare("Si") == 0){
                atomSubs[id] = subst;
            }
            // Invalid type
            else {
                cerr << "Fractional atomic network substitution aborted because atom type "
                "other than Si or O detected" << "\n"
                << "Occurred for atom " << id << "\n";
                return false;
            }
            
            // Added bonded atoms to stack with appropriate subsitution state
            for(unsigned int j = 0; j < bonds[id].size(); j++){
                if(bonds[id][j] != fromID){
                    atomsToProcess.push_back(pair<int, bool>(bonds[id][j], subst));
                    fromIDs.push_back(id);
                }
            }
        }
    }
    
    if(numProc != numAtoms){
        cerr << "Fractional atom network substituion failed because not all atoms "
        "are interconnected" << "\n" << "Visited " << numProc << " out of "
        << numAtoms << " atoms " << "\n";
        return false;
    }
    
    // Perform final consistency check
    for(int i = 0; i < numAtoms; i++){
        // Just check pair of atoms bonded to each oxygen
        if(origNet->atoms[i].type.compare("O") == 0){
            if(atomSubs[bonds[i][0]] && atomSubs[bonds[i][1]]){
                cerr << "Fractional atom network substitution failed final consistency check."
                << "\n" << "Oxygen atom #" << i << " is bound to atoms " <<  bonds[i][0]
                << " and " << bonds[i][1] << " which are of same type. " << "\n"
                << "Algorithmic issue likely." << "\n"
                << "Please contact the source code provider" << "\n";
                return false;
            }
        }
    }
    
    // Copy original information
    origNet->copy(newNet);
    
    int totalSubs = nearestInt(frac*numSi); // How many substituions to actually perform
    *fracSub = (1.0*totalSubs)/numSi;        // Store the actual substitution fraction
    
    srand(randSeed); // Seed random number generator
    
    // Record all originally substituted atoms
    vector<pair<int, int> > subIDs;
    for(int i = 0; i < numAtoms; i++){
        if(atomSubs[i]){
            subIDs.push_back(pair<int,int> (rand(), i));
        }
    }
    
    sort(subIDs.begin(), subIDs.end(), comparePairs); // Sort the pairs by their random values
    
    // Change atoms that should be substituted
    for(*numSubstitutions = 0; *numSubstitutions < totalSubs; (*numSubstitutions)++){
        int id = subIDs[*numSubstitutions].second;
        newNet->atoms[id].type   = "Al";
        newNet->atoms[id].radius = lookupRadius("Al", radial);
    }
    
    
    // Change Oxygen type if O atom connected to Al
    for(int i=0; i < numAtoms; i++){
        if(origNet->atoms[i].type.compare("O") == 0)
        {
            if(newNet->atoms[(bonds[i].at(0))].type.compare("Al") == 0 || newNet->atoms[(
                                                                                         bonds[i].at(1))].type.compare("Al") == 0)
            {
                newNet->atoms[i].type   = "O_Al";
            };
        };
    };
    
    
    
    // Successful substitution
    return true;
}


/* Maciek's substitution function (based on Thomas -fsub function) */
/* Attempt to substitute the specified fraction of Si atom with an Al atom.
 * ATOM_NETWORK may only consist of Si and O atoms, where each Si atom must be
 * bonded to exactly 4 oxygen atoms and each oxygen atom must be bonded to
 * exactly 2 Si atoms. Returns true if the substitution was successful and
 * stores the number of substitutions using the provided reference. The
 * provided boolean specifies whether the seeded silicon atom is substituted
 * or not. This function does not require initial 50/50 substitution as in
 * Thomas version, and therefore is suitable for systems in which not every
 * other atom is a Si atom */
/* bool radial is always set to true in the new version (from 0.3 on) */
bool fracSubstituteAtoms_Maciek(ATOM_NETWORK &origNet, ATOM_NETWORK &newNet,
                                bool substituteSeed, double frac, int randSeed,
                                int &numSubstitutions, double &fracSub, bool radial){
    int numAtoms = origNet.numAtoms;
    double max_bond_length = 1.95;
    
    if(frac > 0.50){
        cerr << "Fractional atomic network substitution failed because the fraction "
        "can not exceed 0.5" << "\n";
        return false;
    }
    if(frac < 0.0){
        cerr << "Fractional atomic network substitution failed because of invalid "
        "negative fraction " << frac << "\n";
        return false;
    }
    
    vector< vector<int> > bonds = vector< vector<int> >(numAtoms, vector<int>());
    for(int i = 0; i < numAtoms; i++){
        ATOM atom_one = origNet.atoms[i];
        for(int j = i + 1; j < numAtoms; j++){
            ATOM atom_two = origNet.atoms[j];
            if(origNet.calcDistance(atom_one, atom_two) < max_bond_length){
                if(atom_one.type.compare(atom_two.type) == 0){
                    cerr << "Fractional atomic network substitution aborted because atoms of "
                    "same type are bonded to one another" << "\n" << "Occurred for type "
                    << atom_one.type << " between atoms " << i << " and " << j << "\n";
                    return false;
                }
                bonds[i].push_back(j);
                bonds[j].push_back(i);
            }
        }
    }
    
    int numSi = 0;
    
    for(int i = 0; i < numAtoms; i++){
        if(origNet.atoms[i].type.compare("Si") == 0){
            if(bonds[i].size() != 4){
                cerr << "Fractional atomic network substitution aborted because Si atom "
                "bonded to " << bonds[i].size() << " other atoms instead of 4" << "\n"
                << "Occurred for atom " << i << "\n";
                return false;
            }
            numSi++;
        }
        else if(origNet.atoms[i].type.compare("O") == 0){
            if(bonds[i].size() != 2){
                cerr << "Fractional atomic network substitution aborted because O atom "
                "bonded to " << bonds[i].size() << " other atoms instead of 2" << "\n"
                << "Occurred for atom " << i << "\n";
                return false;
            }
        }
        else{
            cerr << "Fractional atomic network substitution aborted because atom type "
            "other than Si or O detected" << "\n"
            << "Occurred for atom " << i << "\n";
            return false;
        }
    } // analysis of the input structure completed
    
    // save IDs of Si atoms in SiatomsIDs vector
    vector <int> SiatomsIDs;
    for(int i = 0; i < numAtoms; i++){
        if(origNet.atoms[i].type.compare("Si") == 0){
            SiatomsIDs.push_back(i);
        }
    }
    
    
    int totalSubs = nearestInt(frac*numSi); // How many substituions to actually perform
    fracSub = (1.0*totalSubs)/numSi;        // Store the actual substitution fraction
    
    srand(randSeed); // Seed random number generator
    
    // Record all Si atoms (each of them can be substituted)
    vector<pair<int, int> > subIDs;
    for(int i = 0; i < numSi; i++){
        subIDs.push_back(pair<int,int> (rand(), i));
    }
    
    sort(subIDs.begin(), subIDs.end(), comparePairs); // Sort the pairs by their random values
    
    vector<bool> atomSubs = vector<bool>(numAtoms, false); // Whether atom has been substituted
    
    // Change atoms that should be substituted
    for(numSubstitutions = 0; numSubstitutions < totalSubs; numSubstitutions++){
        atomSubs[SiatomsIDs[subIDs[numSubstitutions].second]]=true;
    } // now atomsSubs points all Si atoms substituted in the intial step
    
    
    // Verify if the random distribution of Al atoms satisfy the condition not to have two neighboring Al atoms (Lowenstein's rule)
    // Fix if neccessary
    int counter=0;
    int problematicAl=-1;
    
    int max_try=10000; // number of swaps when trying to find a configuration satisfying Lowenstein's rule
    
    do{
        problematicAl=-1;
        // Step 1: identify atom with problems
        for(int i = 0; i < numAtoms; i++){
            if(origNet.atoms[i].type.compare("Si") == 0 &&atomSubs[i]==true)
            {
                // find list of 4 neighboring Si/Al atoms
                vector <int> nlist;
                for(int j=0;j<4;j++)
                {
                    if(bonds[bonds[i][j]][0]==i) {nlist.push_back(bonds[bonds[i][j]][1]);} else {nlist.push_back(bonds[bonds[i][j]][0]);};
                }
                
                if(atomSubs[nlist[0]]||atomSubs[nlist[1]]||atomSubs[nlist[2]]||atomSubs[nlist[3]]) // if any of neighbors is Al, move Al from position i
                {
                    problematicAl=i; // Al atom i needs to be this shifted;
                    break;
                }
            }
        }
        
        // Step 2: identify a suitable place for Al atom
        int suitableSi=-1;
        if(problematicAl>-1)
        {
            
            for(int i = 0; i < numAtoms; i++){
                if(origNet.atoms[i].type.compare("Si") == 0 && atomSubs[i]==false)
                {
                    // find list of 4 neighboring Si/Al atoms
                    vector <int> nlist;
                    for(int j=0;j<4;j++)
                    {
                        if(bonds[bonds[i][j]][0]==i) {
                            nlist.push_back(bonds[bonds[i][j]][1]);
                        } else {
                            nlist.push_back(bonds[bonds[i][j]][0]);
                        };
                    }
                    
                    if(!atomSubs[nlist[0]] && !atomSubs[nlist[1]] && 
                       !atomSubs[nlist[2]]&&!atomSubs[nlist[3]]) // if all of neighbors are Si, position i can be used
                    {
                        suitableSi=i; // Si atom i can be replaced with Al;
                        break;
                    }
                }
            }
        } // end of Step 2
        
        // Step 3: Swap 1 & 2
        if(problematicAl>-1&&suitableSi>-1)
        {
            atomSubs[problematicAl]=false;
            atomSubs[suitableSi]=true;
        }else{
            
            if(problematicAl>-1) 
            {
                cerr << "Cannot fix a random Al distribution by swapping Al atoms" << "\n";
                return false;
            }
        }
        // Step 4: increase counter
        counter++;
    }while(problematicAl>-1&&counter<max_try*numSi);
    
    
    if(counter==(max_try*numSi))
    {
        cerr << "Could not fix the initial Al distribution in " << max_try*numSi << "steps\n";
        return false;
    }
    
    
    
    // Now update the original structure according to the generate distribution of Al atoms
    
    // Copy original information
    origNet.copy(&newNet);
    
    
    // Change atoms that should be substituted
    for(int i = 0; i < numAtoms; i++){
        if(atomSubs[i]){
            newNet.atoms[i].type   = "Al";
            newNet.atoms[i].radius = lookupRadius("Al", radial);
        }
    }
    
    
    
    // Change Oxygen type if O atom connected to Al
    for(int i=0; i < numAtoms; i++){
        if(origNet.atoms[i].type.compare("O") == 0)
        {
            if(newNet.atoms[(bonds[i].at(0))].type.compare("Al") == 0 || newNet.atoms[(bonds[i].at(1))].type.compare("Al") == 0)
            {
                newNet.atoms[i].type   = "O_Al";
            };
        };  
    };
    
    
    // Successful substitution
    return true;
}





/* Returns the integer nearest to the provided double.*/
int nearestInt(double num){
    return (int)(floor(num + 0.5));
}

/* Determines whether a specific supercell size satisfies the non-overlapping sphere requirement */
int check_sphere_overlap(int num_a, int num_b, int num_c, double diam, ATOM_NETWORK *atmnet) {
    //check each neighbouring cell in three dimensions and find shortest distance to image
    //it is sufficient to check only those images (a b c) where the leftmost non-zero term is one, e.g. (0 1 -1) but we can skip (0 -1 1), since it is equivalent in magnitude
    double min_d = -1.0; //set min distance to invalid quantity - if we cannot find any result we can return the error case
    bool overlaps = false;
    for(int a=0; a<=1 && !overlaps; a++) {
        for(int b=-1; b<=1 && !overlaps; b++) {
            for(int c=-1; c<=1 && !overlaps; c++) {
                if (
                    (a==0 && b==0 && c==1)
                    ||
                    (a==0 && b==1)
                    ||
                    (a==1)
                    )
                {
                    Point image_abc(num_a*a, num_b*b, num_c*c);
                    Point image_xyz = atmnet->abc_to_xyz(image_abc);
                    double d=image_xyz.magnitude();
                    if(d<min_d || min_d<0) {
                        min_d = d;
                        if(min_d<diam+0.001) overlaps = true;
                    }
                }
            }
        }
    }
    if(min_d<0) return -1; //bad_angles!
    else if(overlaps) return 0; //wrong
    else return 1; //correct
}

/* Returns the id of the VOR_NODE in the provided VORONOI_NETWORK whose coordinates
 * match those of the provided Point.
 */
int getNodeID(Point pt, ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet){
    double minDist = DBL_MAX;
    int minID = -1;
    for(unsigned int i = 0; i < vornet->nodes.size(); i++){
        VOR_NODE curNode = vornet->nodes.at(i);
        double dist = atmnet->calcDistanceXYZ(pt[0], pt[1], pt[2], curNode.x, curNode.y, curNode.z);
        if(dist < threshold)
            return i;
        else {
            if(dist < minDist){
                minDist = dist;
                minID = i;
            }
        }
    }
    
    cerr << "Warning : When identifying Voronoi node, the distance exceeded the threshold of " << threshold << "\n"
    << "Occurred during analysis of " << atmnet->name << "\n"
    << "Closest node was within " << minDist << "\n"
    << "Proceeding with analysis" << "\n";
    return minID;
}

/* returns whether two fractional points are bonded periodically */
bool bonded_abc(XYZ abc1, XYZ abc2, ATOM_NETWORK *c) {
  return c->getDistCalc().minimum_periodic_distance(abc1.x, abc1.y, abc1.z, abc2.x, abc2.y, abc2.z)<GENERAL_BOND_DISTANCE;
}

/* returns whether two fractional points overlap periodically */
bool overlaps_abc(XYZ abc1, XYZ abc2, ATOM_NETWORK *c) {
  return c->getDistCalc().minimum_periodic_distance(abc1.x, abc1.y, abc1.z, abc2.x, abc2.y, abc2.z)<DISTANCE_TOLERANCE;
}

/* returns whether two fractional points (one being an atom, the other a site marker) are close enough for the atom to be considered part of the site for collision detection purposes */
bool is_part_of_site_abc(XYZ abc1, XYZ abc2, ATOM_NETWORK *c) {
  return c->getDistCalc().minimum_periodic_distance(abc1.x, abc1.y, abc1.z, abc2.x, abc2.y, abc2.z)<IS_PART_OF_SITE_DISTANCE;
}

/* put a molecule in its fixed position according to the net - returns the total number of molecules placed so far */
int place_molecule(vector<MOLECULE> *molecules, vector<CONNECTION> *two_way_connections, ATOM_NETWORK *cell, vector<bool> *placed, vector<bool> *connected, vector<MOLECULE> *framework, int num_placed, double *sum_edge_lengths) {
  bool verbose = false;
  //1) if no molecules have been placed yet, just fix the first one in position with its first edge at the origin
  if(num_placed==0) {
if(verbose) printf("DEBUG: detected that no molecules have yet been placed - positioning first molecule with first edge at the origin\n");
    int v_index = two_way_connections->at(0).v1;
    int e_index = two_way_connections->at(0).e1;
    MOLECULE placing = molecules->at(v_index);
    int site_index = placing.permutation.at(e_index);
    int atom_index = placing.sites.at(site_index);
    XYZ site = placing.atoms_xyz.at(atom_index);
    XYZ shift = site.scale(-1);
if(verbose) printf("DEBUG: v=%d, e=%d; site=%.3f %.3f %.3f; shift=%.3f %.3f %.3f\n", v_index, e_index, site.x, site.y, site.z, shift.x, shift.y, shift.z);
    framework->at(v_index) = (translate(placing, shift));
    placed->at(v_index) = true;
    return num_placed+1;
  }
  //2) otherwise, need to find a new molecule to place
if(verbose) printf("DEBUG: detected that %d molecules have already been placed\n", num_placed);
  int num_molecules = molecules->size();
  int num_connections = two_way_connections->size();
  int molecule_index = -1;
  for(int i=0; i<num_molecules && molecule_index==-1; i++) {
    if(!placed->at(i)) { //if this molecule has not yet been placed - it is a candidate for placement, but we must have an adjacent molecule already placed which we can attach it to
      int connection_index = -1;
      for(int j=0; j<num_connections && connection_index==-1; j++) { //is there any two way connection which connects this new molecule and any previously placed one?
        if(!connected->at(j)) { //if this connection was already made, we don't need to check it again
          CONNECTION c = two_way_connections->at(j);
          int placing_edge_index = -1, fixed_edge_index = -1;
          int placing_vertex_index = -1, fixed_vertex_index = -1;
          if(c.v1==i) { //vertex 1 in this connection is the molecule we are placing
            if(placed->at(c.v2)) { //vertex 2 in this connection is a molecule which has already been placed
              placing_edge_index = c.e1;
              fixed_edge_index = c.e2;
              placing_vertex_index = c.v1;
              fixed_vertex_index = c.v2;
            }
          } else if(c.v2==i) { //vertex 2 in this connection is the molecule we are placing
            if(placed->at(c.v1)) { //vertex 1 in this connection is a molecule which has already been placed
              placing_edge_index = c.e2;
              fixed_edge_index = c.e1;
              placing_vertex_index = c.v2;
              fixed_vertex_index = c.v1;
            }
          }
          if(fixed_edge_index!=-1 && placing_edge_index!=-1 && fixed_vertex_index!=-1 && placing_vertex_index!=-1) { //this connection connects this molecule and a fixed molecule, and we have stored the relevant edge and vertex indices - we can proceed to place this molecule!
if(verbose) printf("DEBUG: found that new vertex %d edge %d is connected to fixed vertex %d edge %d\n", placing_vertex_index, placing_edge_index, fixed_vertex_index, fixed_edge_index);
            MOLECULE placing = molecules->at(placing_vertex_index);
            int placing_site_index = placing.permutation.at(placing_edge_index);
            int placing_atom_index = placing.sites.at(placing_site_index);
            MOLECULE fixed = framework->at(fixed_vertex_index);
            int fixed_site_index = fixed.permutation.at(fixed_edge_index);
            int fixed_atom_index = fixed.sites.at(fixed_site_index);
if(verbose) printf("DEBUG: %s in placing connecting to %s in fixed\n", placing.atoms_type.at(placing_atom_index).c_str(), fixed.atoms_type.at(fixed_atom_index).c_str());
            if(placing.atoms_type.at(placing_atom_index)!=fixed.atoms_type.at(fixed_atom_index)) {//only make the connection if the sites are the same type!
if(verbose) printf("DEBUG: trying to connect conflicting connection sites - aborting construction attempt\n");
              return num_placed;
            }
//if(1) exit(1);
            XYZ placing_site = placing.atoms_xyz.at(placing_atom_index);
            XYZ fixed_site = fixed.atoms_xyz.at(fixed_atom_index);
            XYZ shift = fixed_site-placing_site;
if(verbose) printf("DEBUG: site index in fixed = %d, site_index in placing = %d\n", fixed_site_index, placing_site_index);
if(verbose) printf("DEBUG: fixed_site=%.3f %.3f %.3f; placing_site=%.3f %.3f %.3f; shift=%.3f %.3f %.3f\n", fixed_site.x, fixed_site.y, fixed_site.z, placing_site.x, placing_site.y, placing_site.z, shift.x, shift.y, shift.z);
            framework->at(placing_vertex_index) = (translate(placing, shift));
            placed->at(placing_vertex_index) = true; //record this new vertex as now being fixed in place
            connected->at(j) = true; //record that this connection was used
            double com_dist = (fixed.com - framework->at(placing_vertex_index).com).magnitude();
if(verbose) printf("DEBUG: made a connection indicating a topological edge length of %.3f\n", com_dist);
            XYZ start_abc = cell->vertices.at(fixed_vertex_index).abc;
            XYZ end_abc = cell->vertices.at(fixed_vertex_index).edges.at(fixed_edge_index);
if(verbose) printf("\tDEBUG: this edge in the underlying net goes from fractional %.3f %.3f %.3f to fractional %.3f %.3f %.3f\n", start_abc.x, start_abc.y, start_abc.z, end_abc.x, end_abc.y, end_abc.z);
            XYZ shift_abc = end_abc-start_abc;
            XYZ shift_xyz = framework->at(placing_vertex_index).com-fixed.com;
if(verbose) printf("\tDEBUG: this edge in the underlying net had a fractional shift %.3f %.3f %.3f; here we observe a cartesian shift %.3f %.3f %.3f\n", shift_abc.x, shift_abc.y, shift_abc.z, shift_xyz.x, shift_xyz.y, shift_xyz.z);
            XYZ unit_cell_parameter_suggestion(shift_xyz.x/shift_abc.x, shift_xyz.y/shift_abc.y, shift_xyz.z/shift_abc.z);
if(verbose) printf("\tDEBUG: suggests unit cell dimensions %.3f %.3f %.3f\n", unit_cell_parameter_suggestion.x, unit_cell_parameter_suggestion.y, unit_cell_parameter_suggestion.z);
            (*sum_edge_lengths)+=com_dist;
            return num_placed+1;
          }
        }
      }
      if(connection_index==-1) {
if(verbose) printf("DEBUG: could not find any connection between vertex %d and any fixed vertex\n", i);
      }
    }
  }
  printf("WARNING: could not find a molecule which has not been placed and which is adjacent to an already placed molecule - this may happen if the net is interpenetrated (not currently handled with the connection-based alignment method)\n");
  return num_placed;
}

/* creates a unit cell based on supplied vectors */
void create_unit_cell_from_vectors(vector<XYZ> *vecs, ATOM_NETWORK *cell) {
  bool verbose = false;
  //DEBUG: complex version - attempts selection of appropriate vectors for each axis
  int num_vecs = vecs->size();
  if(num_vecs!=3) {
    printf("ERROR: create_unit_cell_from_vectors() called with !=3 (%d) vectors\n", num_vecs);
    exit(EXIT_FAILURE);
  }
  vector<int> assignments, scales;
  vector<bool> is_assigned;
  for(int i=0; i<num_vecs; i++) { 
    if(i<2) {
      assignments.push_back(-1);
      scales.push_back(1);
    }
    is_assigned.push_back(false);
  }
  //1) find which vector is closest to the x axis, and which way around it is; repeat for y
  for(int ax=0; ax<2; ax++) {
    XYZ axis(0,0,0);
    if(ax==0) axis.x=1;
    if(ax==1) axis.y=1;
    int index_closest_to_axis = -1;
    int scale = 1;
    double best_angle = -1;
    for(int i=0; i<num_vecs; i++) {
      if(!is_assigned.at(i)) { //candidate
        XYZ cand = vecs->at(i);
        XYZ reverse_cand = cand.scale(-1);
        double angle = axis.angle_between(cand);
        if(angle<best_angle || best_angle<0) {
          best_angle = angle;
          index_closest_to_axis = i;
          scale = 1;
        }
        angle = axis.angle_between(reverse_cand);
        if(angle<best_angle || best_angle<0) {
          best_angle = angle;
          index_closest_to_axis = i;
          scale = -1;
        }
      }
    }
    assignments.at(ax) = index_closest_to_axis;
    scales.at(ax) = scale;
    is_assigned.at(index_closest_to_axis) = true;
    string axis_name("ERROR");
    if(ax==0) axis_name = "x";
    if(ax==1) axis_name = "y";
    if(verbose) {
      if(scale==-1) printf("DEBUG: reversed vector %d is closest to the %s axis\n", index_closest_to_axis, axis_name.c_str());
      else printf("DEBUG: vector %d is closest to the %s axis\n", index_closest_to_axis, axis_name.c_str());
    }
  }
  //2) z axis is different - we know x and y, and want to find z to satisfy the right hand rule!
  XYZ x = vecs->at(assignments.at(0)).scale(scales.at(0));
  XYZ y = vecs->at(assignments.at(1)).scale(scales.at(1));
  XYZ z(0,0,1);
  bool found_z = false;
  for(int i=0; i<num_vecs; i++) {
    if(!is_assigned.at(i)) { //found z
      if(found_z) {
        printf("ERROR: after setting x and y vectors, more than one vector remains to be assigned to z\n");
        exit(EXIT_FAILURE);
      }
      z = vecs->at(i);
      found_z = true;
    }
  }
  XYZ cross = x.cross(y);
  if(verbose) printf("DEBUG: x = %.3f %.3f %.3f, y = %.3f %.3f %.3f, cross = %.3f %.3f %.3f; need to know which orientation of z = %.3f %.3f %.3f satisfies right hand rule\n", x.x, x.y, x.z, y.x, y.y, y.z, cross.x, cross.y, cross.z, z.x, z.y, z.z);
  int z_scale = 1;
  if(cross.angle_between(z)>cross.angle_between(z.scale(-1))) z_scale = -1;
  z = z.scale(z_scale);
  if(verbose) printf("DEBUG: selected z to be %.3f %.3f %.3f\n", z.x, z.y, z.z);
  //3) build the cell
  cell->make(x,y,z);
}

/* populate the vector of two-way connections - this vector can be inspected to find which vertices in the net are connected to each other - returns true if net is disconnected */
bool find_two_way_connections(ATOM_NETWORK *cell, vector<CONNECTION> *two_way_connections) {
  bool verbose = false;
  //1) for each VERTEX in the ATOM_NETWORK, look at the edge positions and find which VERTEX it overlaps with periodically - it might be itself!
  vector<CONNECTION> connections;
  int num_v = cell->vertices.size();
  if(verbose) printf("DEBUG: num_v = %d\n", num_v);
  for(int i=0; i<num_v; i++) {
    VERTEX vert = cell->vertices.at(i);
    int num_e = vert.edges.size();
    if(verbose) printf("DEBUG: num_e = %d\n", num_e);
    for(int j=0; j<num_e; j++) {
      XYZ e = vert.edges.at(j);
      int overlap_v_index = -1;
      for(int k=0; k<num_v && overlap_v_index==-1; k++) {
        XYZ v = cell->vertices.at(k).abc;
if(verbose) printf("DEBUG: checking for overlap between edge end position %.3f %.3f %.3f and other vertex position %.3f %.3f %.3f\n", e.x, e.y, e.z, v.x, v.y, v.z);
        if(overlaps_abc(e, v, cell))
          overlap_v_index = k;
      }
      if(overlap_v_index==-1) {
        printf("ERROR: could not construct net: no vertex could be found which overlaps periodically with vertex %d edge %d - check cgd file (did you forget the \":H\" in the symmetry group?)\n", i, j);
        exit(EXIT_FAILURE);
      } else {
        XYZ v = cell->vertices.at(overlap_v_index).abc;
        XYZ shift = e-v;
        int a_shift = (int)(round(shift.x));
        int b_shift = (int)(round(shift.y));
        int c_shift = (int)(round(shift.z));
        if(verbose) printf("DEBUG: vertex %d edge %d at %.3f %.3f %.3f overlaps with vertex %d at %.3f %.3f %.3f with periodicity (%d %d %d)\n", i, j, e.x, e.y, e.z, overlap_v_index, cell->vertices.at(overlap_v_index).abc.x, cell->vertices.at(overlap_v_index).abc.y, cell->vertices.at(overlap_v_index).abc.z, a_shift, b_shift, c_shift);
        CONNECTION c(i, overlap_v_index, j, a_shift, b_shift, c_shift);
        connections.push_back(c);
      }
    }
  }
  int num_connections = connections.size();
  if(verbose) printf("DEBUG: there are %d one-way connections\n", num_connections);
  if(num_connections%2!=0) {
    printf("ERROR: there are an odd number (%d) of one-way connections between vertices - this should not be the case because connections are expressed redundantly\n", num_connections);
    exit(EXIT_FAILURE);
  }
  //2) connections only go one way right now - we cannot yet position molecules because we need to know which edge is used in each VERTEX!
  vector<bool> considered;
  for(int i=0; i<num_connections; i++) considered.push_back(false);
  for(int i=0; i<num_connections; i++) {
    if(!considered.at(i)) { //we are looking at one-way connections - if we have not yet found a match between this connection and another, try and find one
      bool matched = false;
      CONNECTION ci = connections.at(i);
      for(int j=0; j<num_connections; j++) { //we continue searching for matches even after finding one - because we want to handle the error case where edges are declared ambiguously (would be caused by invalid cgd net)
        if(!considered.at(j)) {
          CONNECTION cj = connections.at(j);
          if(matches(ci,cj)) {
            if(matched) {
              printf("ERROR: determined that vertex %d edge %d overlaps with more than one vertex!\n", ci.v1, ci.e1);
              exit(EXIT_FAILURE);
            }
            CONNECTION cij = ci;
            cij.e2 = cj.e1;
            if(verbose) printf("DEBUG: vertex %d edge %d overlaps with vertex %d edge %d with periodicity (%d %d %d)\n", cij.v1, cij.e1, cij.v2, cij.e2, cij.a, cij.b, cij.c);
            two_way_connections->push_back(cij);
            considered.at(i) = true;
            considered.at(j) = true;
            matched = true;
          }
        }
      }
      if(!matched) {
        printf("ERROR: could not find corresponding connection for vertex %d edge %d overlapping with vertex %d\n", ci.v1, ci.e1, ci.v2);
        exit(EXIT_FAILURE);
      }
    }
  }
  int num_two_way_connections = two_way_connections->size();
  if(verbose) printf("DEBUG: there are %d two-way connections\n", num_two_way_connections);
  //3) now we know all the two-way connections, we can recursively traverse until we establish whether the net is disconnected - we need to know this in order to avoid trying the connection-based assembly method later, only to find that it fails
  vector<bool> visited_vertex;
  for(int i=0; i<num_v; i++) visited_vertex.push_back(false);
  recursive_visit_vertices(0, two_way_connections, &visited_vertex);
  bool is_disconnected = false;
  for(int i=0; i<num_v && !is_disconnected; i++) {
    if(!visited_vertex.at(i)) is_disconnected = true;
  }
  return is_disconnected;
}

/* take a MOLECULE, and return a list of new MOLECULEs aligned well with a specified VERTEX in some ATOM_NETWORK */
vector<MOLECULE> get_multiple_best_RMSD_fits(MOLECULE mol, ATOM_NETWORK *cell, int vertex_ID, string net, string prefix, int mol_ID) {
  bool verbose = false; //for printing data to terminal
  //1) can we fit this MOLECULE to this VERTEX?
  int num_fit_sites_no_dummy = cell->vertices.at(vertex_ID).edges.size();
  int num_fit_sites_dummy = cell->vertices.at(vertex_ID).dummy_edges.size();
  int num_fit_sites = num_fit_sites_no_dummy + num_fit_sites_dummy;
  int num_mol_sites_no_dummy = mol.sites.size();
  int num_mol_sites_dummy = mol.dummy_sites.size();
  int num_mol_sites = num_mol_sites_no_dummy + num_mol_sites_dummy;
  if(num_fit_sites_dummy!=num_mol_sites_dummy || num_fit_sites_no_dummy!=num_mol_sites_no_dummy) { //we can only fit the molecule to the given vertex is the number of sites is the same in each
    printf("ERROR: cannot fit molecule with %d sites and %d dummy sites to a vertex with %d sites and %d dummy sites!\n", num_mol_sites_no_dummy, num_mol_sites_dummy, num_fit_sites_no_dummy, num_fit_sites_dummy);
    exit(EXIT_FAILURE);
  }

  //2) how many ways can we do it?
  vector<int> permute_base_no_dummy, permute_base_dummy;
  for(int i=0; i<num_fit_sites_no_dummy; i++) permute_base_no_dummy.push_back(i);
  for(int i=0; i<num_fit_sites_dummy; i++) permute_base_dummy.push_back(i+num_fit_sites_no_dummy);
  vector< vector<int> > permutations_no_dummy, permutations_dummy;

  permute(permute_base_no_dummy, 0, &permutations_no_dummy); //fills in a vector with all possible permutations by which we can try and orient our two objects, e.g. (0123), (0132), ...
  permute(permute_base_dummy, 0, &permutations_dummy);

  //2b) now collapse permutations for non-dummy and dummy together to get acceptable permutations
  vector< vector<int> > permutations;
  int num_perm_no_dummy = permutations_no_dummy.size();
  int num_perm_dummy = permutations_dummy.size();
  for(int i=0; i<num_perm_no_dummy; i++) {
    for(int j=0; j<num_perm_dummy; j++) {
      vector<int> new_perm = permutations_no_dummy.at(i);
      for(int k=0; k<num_fit_sites_dummy; k++) new_perm.push_back(permutations_dummy.at(j).at(k));
      permutations.push_back(new_perm);
    }
  }

  int num_perm = permutations.size();
  if(verbose) {
    printf("there are %d permutations of length %d:\n", num_perm, num_fit_sites);
    for(int i=0; i<num_perm; i++) {
      printf("\t");
      for(int j=0; j<permutations.at(i).size(); j++) {
        printf("%d ", permutations.at(i).at(j));
      }
      printf("\n");
    }
  }
  //3) which way gives us the lowest rmsd?
  double best_rmsd = -1; //dummy value
  int best_rmsd_ID = -1;
  vector<FIT> vector_of_fits;
  vector<MOLECULE> vector_of_centred_molecules;
  for(int p=0; p<num_perm; p++) { //try each permutation!
    vector<int> perm = permutations.at(p);
    double fixed[num_fit_sites][3], moving[num_fit_sites][3];
    XYZ fixed_CoM = origin, moving_CoM = origin;
    for(int i=0; i<num_fit_sites; i++) {
      XYZ fix(0,0,0);
      if(i<num_fit_sites_no_dummy)
        fix = cell->abc_to_xyz_returning_XYZ(cell->vertices.at(vertex_ID).edges.at(i)); //the edges are always in the same order
      else
        fix = cell->abc_to_xyz_returning_XYZ(cell->vertices.at(vertex_ID).dummy_edges.at(i-num_fit_sites_no_dummy));
      fixed[i][0] = fix.x;
      fixed[i][1] = fix.y;
      fixed[i][2] = fix.z;
      fixed_CoM = fixed_CoM+fix;
      XYZ move(0,0,0);
      if(perm.at(i)<num_fit_sites_no_dummy)
        move = mol.atoms_xyz.at(mol.sites.at(perm.at(i))); //permutation encodes in what order the sites should go to correspond to the edges
      else
//        move = mol.dummy_sites.at(perm.at(i)-num_fit_sites_no_dummy);
        move = mol.atoms_xyz.at(mol.dummy_sites.at(perm.at(i)-num_fit_sites_no_dummy));
      moving[i][0] = move.x;
      moving[i][1] = move.y;
      moving[i][2] = move.z;
      moving_CoM = moving_CoM+move;
    }

    if(verbose) {
      printf("PERMUTATION %d:\n\t", p);
      for(int i=0; i<num_fit_sites; i++) {
        printf("%d ", perm.at(i));
      }
      printf("\n");
      printf("\tfits the following positions:\n");
      for(int i=0; i<num_fit_sites; i++) {
        printf("\t\t%.3f %.3f %.3f to %.3f %.3f %.3f\n", moving[i][0], moving[i][1], moving[i][2], fixed[i][0], fixed[i][1], fixed[i][2]);
      }
    }
    //3a) evaluate this permutation
    double mov_com[3], mov_to_ref[3]; //leave these blank, they are filled in by the rmsd method
    double rotation_matrix[3][3];
    for(int i=0; i<3; i++) for(int j=0; j<3; j++) rotation_matrix[i][j] = 0; //dummy values - filled in by the rmsd method
    double rmsd = 0; //dummy value - filled in by the rmsd method
    calculate_rotation_rmsd(fixed, moving, num_fit_sites, mov_com, mov_to_ref, rotation_matrix, &rmsd); //calling rmsd code - writes the corresponding rotation matrix
    //3b) some permutations are not possible, and produce NAN
    bool valid = true;
    if(isnan(rmsd)) valid=false;
    for(int i=0; i<3 && valid; i++) for(int j=0; j<3 && valid; j++) if(isnan(rotation_matrix[i][j])) valid=false;
    if(valid) {
      FIT f;
      f.fit = rotate(mol, rotation_matrix);
      MOLECULE centred_aligned_mol = translate(f.fit, get_mol_site_CoM(&f.fit).scale(-1));
      f.fit.permutation = permutations.at(p);
      f.perm_ID = p;
      f.rmsd = rmsd;
//try and select alignments appropriately ...
      if(best_rmsd<0 || rmsd<best_rmsd) {
        best_rmsd = rmsd;
//        best_rmsd_ID = p;
        best_rmsd_ID = vector_of_fits.size();
      }
      vector_of_fits.push_back(f);
      vector_of_centred_molecules.push_back(centred_aligned_mol);
    } else if(verbose) printf("\tTHIS PERMUTATION COULD NOT PRODUCE A VALID ROTATION MATRIX\n");
  }

  //now return all molecules which are within a tolerance of the best RMSD
  int num_fits = vector_of_fits.size();
if(verbose) printf("DEBUG: num valid rotation matrices = %d; best_rmsd_ID = %d\n", num_fits, best_rmsd_ID);
  vector<MOLECULE> vec, vector_of_unique_centred_molecules;
  vector_of_unique_centred_molecules.push_back(vector_of_centred_molecules.at(best_rmsd_ID));
  FIT f0 = vector_of_fits.at(best_rmsd_ID);
  vec.push_back(f0.fit);
if(verbose) {
  string rotmol_name = prefix+"_molecule_ID_"+convertToString(mol_ID)+"_basic_vertex_ID_"+convertToString(vertex_ID)+"_permutation_"+convertToString(f0.perm_ID)+".xyz"; //write molecule to xyz file
  FILE *rotmol = fopen(rotmol_name.c_str(),"w");
  if(rotmol==NULL) {
    printf("ERROR: could not open output rotated molecule file with name %s\n", rotmol_name.c_str());
    exit(EXIT_FAILURE);
  }
  MOLECULE rotated = translate(f0.fit, get_mol_site_CoM(&f0.fit).scale(-1));
  write_molecule(rotmol, &rotated, net, vertex_ID, -1, true); //oriented to basic vertex basic_ID with symmetry "-1" (not applicable here) - this will be noted in the file - true indicates 'write out sites'
  fclose(rotmol);
  printf("\tmolecule %d oriented to basic vertex %d written to %s\n", mol_ID, vertex_ID, rotmol_name.c_str());
}
  for(int i=0; i<num_fits; i++) {
    if(i!=best_rmsd_ID) { //don't check the one we pushed already
      if(vector_of_fits.at(i).rmsd<best_rmsd*TOLERANCE_PROPORTION) { //acceptable RMSD?
        if(molecule_alignment_chemistry_is_unique(&vector_of_centred_molecules.at(i), &vector_of_unique_centred_molecules)) { //unique molecule?
          vector_of_unique_centred_molecules.push_back(vector_of_centred_molecules.at(i));
          FIT f = vector_of_fits.at(i);
          vec.push_back(f.fit);
if(verbose) {
  string rotmol_name = prefix+"_molecule_ID_"+convertToString(mol_ID)+"_basic_vertex_ID_"+convertToString(vertex_ID)+"_permutation_"+convertToString(f.perm_ID)+".xyz"; //write molecule to xyz file
  FILE *rotmol = fopen(rotmol_name.c_str(),"w");
  if(rotmol==NULL) {
    printf("ERROR: could not open output rotated molecule file with name %s\n", rotmol_name.c_str());
    exit(EXIT_FAILURE);
  }
  MOLECULE rotated = translate(f.fit, get_mol_site_CoM(&f.fit).scale(-1));
  write_molecule(rotmol, &rotated, net, vertex_ID, -1, true); //oriented to basic vertex basic_ID with symmetry "-1" (not applicable here) - this will be noted in the file - true indicates 'write out sites'
  fclose(rotmol);
  printf("\tmolecule %d oriented to basic vertex %d written to %s\n", mol_ID, vertex_ID, rotmol_name.c_str());
}
        }
      }
    }
  }
if(verbose) printf("DEBUG: of %d permutations, %d unique alignments were generated\n", num_perm, (int)(vec.size()));
  return vec;
}

void put_atoms_in_atom_network(ATOM_NETWORK *cell, vector<MOLECULE> *assembled_molecules) {
  cell->numAtoms = 0;
  int num_assembled_molecules = assembled_molecules->size();
  for(int i=0; i<num_assembled_molecules; i++) {
    int num_at = assembled_molecules->at(i).atoms_xyz.size();
    for(int j=0; j<num_at; j++) {
      bool is_site = false, is_dummy_site = false;
      int num_sites = assembled_molecules->at(i).sites.size();
      int num_dummy_sites = assembled_molecules->at(i).dummy_sites.size();
      for(int k=0; k<num_sites && !is_site; k++) {
        is_site = (j == assembled_molecules->at(i).sites.at(k));
      }
      for(int k=0; k<num_dummy_sites && !is_dummy_site && !is_site; k++) {
        is_dummy_site = (j == assembled_molecules->at(i).dummy_sites.at(k));
      }
      if(!is_site && !is_dummy_site) { //don't save sites! they are dummy atoms and we don't need them in the ATOM_NETWORK (they will still be accessible in the vector<MOLECULE> assembled_molecules)
        XYZ xyz = assembled_molecules->at(i).atoms_xyz.at(j);
        XYZ abc = trans_to_origuc(cell->xyz_to_abc_returning_XYZ(xyz));
        string name = assembled_molecules->at(i).atoms_type.at(j);
        string label = assembled_molecules->at(i).atoms_label.at(j);
        ATOM a(xyz, name, label, lookupRadius(name, true));
        a.a_coord = abc.x;
        a.b_coord = abc.y;
        a.c_coord = abc.z;
        cell->atoms.push_back(a);
        cell->numAtoms++;
      }
    }
  }
}

/* net-based alignment method, one of the two ways of creating framework models from aligned MOLECULEs - in this method, we arrange a vector of MOLECULEs in 3D space preserving the unit cell of the underlying ATOM_NETWORK (may cause MOLECULEs to be slightly distorted in order to fit to the net */
ATOM_NETWORK build_framework(vector<MOLECULE> *molecules, ATOM_NETWORK *cell, vector<MOLECULE> *assembled_molecules, double mean_edge_length, double unit_edge_length) {
  //1) set up the unit cell based on the underlying cell and the calculated edge length as a scaling factor
  ATOM_NETWORK alt_framework;
  alt_framework.make(mean_edge_length*cell->a/unit_edge_length, mean_edge_length*cell->b/unit_edge_length, mean_edge_length*cell->c/unit_edge_length, cell->alpha, cell->beta, cell->gamma);
  //2) for each vertex fractional position in the underlying cell, put the appropriate oriented molecule at the corresponding fractional position in the new cell, using the centre of mass of the molecule
  int num_mols = molecules->size();
  int num_vertices = cell->vertices.size();
  if(num_mols!=num_vertices) {
    printf("ERROR: the number of oriented molecules (%d) is not equal to the number of vertices in the underlying cell (%d)\n", num_mols, num_vertices);
  }
  for(int i=0; i<num_mols; i++) {
    XYZ abc = cell->vertices.at(i).abc;
    XYZ new_xyz = alt_framework.abc_to_xyz_returning_XYZ(abc);
    MOLECULE old_m = molecules->at(i);
    XYZ shift = new_xyz - old_m.com;
    assembled_molecules->push_back(translate(old_m, shift));
  }
  //3) put atoms into the cell, converting to fractional coords
  put_atoms_in_atom_network(&alt_framework, assembled_molecules);
  return alt_framework;
}

/* a cut down version of the 'connect_molecules' method which just connects two adjacent molecules, returning the resulting edge length */
double determine_edge_length(vector<MOLECULE> *molecules, ATOM_NETWORK *cell, vector<MOLECULE> *assembled_molecules, vector<CONNECTION> *two_way_connections) {
  bool verbose = false;
  int num_placed = 0;
  bool placed_ok = true;
  double sum_edge_lengths = 0;
  vector<bool> placed;
  int num_mols = molecules->size();
  for(int i=0; i<num_mols; i++) placed.push_back(false);
  vector<bool> connected;
  int num_cons = two_way_connections->size();
  for(int i=0; i<num_cons; i++) connected.push_back(false);
  for(int i=0; i<num_mols; i++) assembled_molecules->push_back(molecules->at(i)); //this vector starts out the same as the molecules vector, but the molecules are then each shifted into their relative positions to form the assembled_molecules
  while(num_placed<2 && placed_ok) {
    if(verbose) printf("DEBUG: num_placed = %d; placing...\n", num_placed);
    int new_num_placed = place_molecule(molecules, two_way_connections, cell, &placed, &connected, assembled_molecules, num_placed, &sum_edge_lengths);
    if(verbose) printf("DEBUG: after calling place_molecule, %d molecules have been placed\n", new_num_placed);
    if(new_num_placed==num_placed+1) {
      num_placed = new_num_placed;
    } else placed_ok = false;
  }
  if(!placed_ok || num_placed!=2) {
//    printf("ERROR: could not find a molecule which has not been placed and which is adjacent to an already placed molecule in order to connected two building blocks together to determine edge length - is the net valid?\n");
//    exit(EXIT_FAILURE);
if(verbose) printf("NOTICE: could not place any molecules - returning invalid distance\n");
    return -1;
  }
  return sum_edge_lengths;
}

/* connection-based alignment method, one of the two ways of creating framework models from aligned MOLECULEs - in this method, the MOLECULEs are connected one by one, and the unit cell is then assigned afterwards (may be different angles to the net, etc.) */
ATOM_NETWORK connect_molecules(vector<MOLECULE> *molecules, ATOM_NETWORK *cell, vector<MOLECULE> *assembled_molecules, vector<CONNECTION> *two_way_connections, double *mean_edge_length, int dimensionality, bool a_periodicity, bool b_periodicity, bool c_periodicity, double unit_edge_length) {
  bool verbose = false;
  ATOM_NETWORK unit_cell_any_matrix;
  //1) keep track of which molecules have been fixed in place already, and which connections were made
  vector<bool> placed;
  int num_mols = molecules->size();
  for(int i=0; i<num_mols; i++) placed.push_back(false);
  vector<bool> connected;
  int num_cons = two_way_connections->size();
  for(int i=0; i<num_cons; i++) connected.push_back(false);
  //2) fix each molecule in place
  for(int i=0; i<num_mols; i++) assembled_molecules->push_back(molecules->at(i)); //this vector starts out the same as the molecules vector, but the molecules are then each shifted into their relative positions to form the assembled_molecules
  int num_placed = 0;
  bool placed_ok = true;
  double sum_edge_lengths = 0;
  if(verbose) printf("DEBUG: about to place molecules\n");
  while(num_placed<num_mols && placed_ok) {
    int new_num_placed = place_molecule(molecules, two_way_connections, cell, &placed, &connected, assembled_molecules, num_placed, &sum_edge_lengths);
    if(verbose) printf("DEBUG: after calling place_molecule, %d molecules have been placed\n", new_num_placed);
    if(new_num_placed==num_placed+1) {
      num_placed = new_num_placed;
    } else placed_ok = false;
  }
  if(num_placed>1) (*mean_edge_length) = ((double)(sum_edge_lengths))/((double)(num_placed-1));
  if(num_placed!=num_mols) {
//    printf("ERROR: only %d of %d molecules were fixed in place during connection-based assembly method - this method should have been automatically bypassed if a disconnected net was provided\n", num_placed, num_mols);
//    exit(EXIT_FAILURE);
if(verbose) printf("NOTICE: only %d of %d molecules were fixed in place during connection-based assembly method - returning empty framework\n", num_placed, num_mols);
    return unit_cell_any_matrix;
  } else {
    if(verbose) printf("DEBUG: all %d of %d molecules were fixed in place successfully!\n", num_placed, num_mols);
  }

  //3) new new idea - look for loops in the two-way connections, and calculate their shifts, using this information to determine the net?
  vector<bool> edge_traversed;
  for(int i=0; i<num_cons; i++) edge_traversed.push_back(false);
  vector<int> vertices_visited;
  vector<int> loops_a, loops_b, loops_c;
  vector<XYZ> loops_xyz;
  vector<int> unit_cell_vector_IDs;
  vector<XYZ> unit_cell_vectors;
  vector<int> vertex_visit_a, vertex_visit_b, vertex_visit_c;
  vector<XYZ> vertex_visit_xyz;
  vector<bool> vertex_visited;
  int num_v = cell->vertices.size();
  for(int i=0; i<num_v; i++) {vertex_visited.push_back(false); vertex_visit_a.push_back(0); vertex_visit_b.push_back(0); vertex_visit_c.push_back(0); vertex_visit_xyz.push_back(origin);}
  bool found_all_uc_vectors = recursive_find_loops(0, 0, 0, 0, origin, two_way_connections, assembled_molecules, &edge_traversed, &vertex_visit_a, &vertex_visit_b, &vertex_visit_c, &vertex_visit_xyz, &vertex_visited, &loops_a, &loops_b, &loops_c, &loops_xyz, &unit_cell_vector_IDs, &unit_cell_vectors, dimensionality);
if(verbose) {
  int num_loops = loops_xyz.size();
  for(int i=0; i<num_loops; i++) {
    XYZ l = loops_xyz.at(i);
    printf("DEBUG: unique loop was found with periodicity (%d %d %d) and Cartesian shift %.3f %.3f %.3f\n", loops_a.at(i), loops_b.at(i), loops_c.at(i), l.x, l.y, l.z);
  }
}
  if(found_all_uc_vectors) {
if(verbose) printf("DEBUG: found all uc vectors from loop analysis\n");
  } else {
    printf("ERROR: did not find all uc vectors from loop analysis\n");
    exit(EXIT_FAILURE);
  }

  //if we are dealing with a 2D net (a layer), then we need to introduce a periodic connection running in that direction
  if(dimensionality==2) {
if(verbose) printf("DEBUG: imposing a three dimensional structure on this two-dimensional framework in order to define a unit cell\n");
    int num_override = 0;
    if(!a_periodicity) {
      XYZ vec((*mean_edge_length)*cell->a/unit_edge_length,0,0);
      unit_cell_vectors.push_back(vec);
      unit_cell_vector_IDs.push_back(0);
      num_override++;
    }
    if(!b_periodicity) {
      XYZ vec(0,(*mean_edge_length)*cell->b/unit_edge_length,0);
      unit_cell_vectors.push_back(vec);
      unit_cell_vector_IDs.push_back(1);
      num_override++;
    }
    if(!c_periodicity) {
      XYZ vec(0,0,(*mean_edge_length)*cell->c/unit_edge_length);
      unit_cell_vectors.push_back(vec);
      unit_cell_vector_IDs.push_back(2);
      num_override++;
    }
    if(num_override!=1) {
      printf("ERROR: was expecting to override exactly 1 cell side length value, but %d were overwritten - this is a bug\n", num_override);
      exit(EXIT_FAILURE);
    }
  }
  vector<XYZ> sorted_unit_cell_vectors;
  for(int i=0; i<3; i++) {
    bool pushed = false;
    for(int j=0; j<3 && !pushed; j++) {
      if(unit_cell_vector_IDs.at(j)==i) {
        sorted_unit_cell_vectors.push_back(unit_cell_vectors.at(j));
        pushed = true;
      }
    }
  }
  create_unit_cell_from_vectors(&sorted_unit_cell_vectors, &unit_cell_any_matrix);

  //6) put atoms into the cell, converting to fractional coords
  put_atoms_in_atom_network(&unit_cell_any_matrix, assembled_molecules);
  //7) that ATOM_NETWORK was constructed with matrices based on periodic connections, that may not conform to the standard axes relations - build a new cell with a normal matrix, based on this one, and copy the atom fractional coords into this system
  ATOM_NETWORK unit_cell;
  unit_cell_any_matrix.copy(&unit_cell);
  unit_cell.make(unit_cell_any_matrix.a, unit_cell_any_matrix.b, unit_cell_any_matrix.c, unit_cell_any_matrix.alpha, unit_cell_any_matrix.beta, unit_cell_any_matrix.gamma);
  return unit_cell;
}

/* collision detection routine */
bool check_for_collision(ATOM_NETWORK *framework, vector<MOLECULE> *assembled_molecules, vector<CONNECTION> *two_way_connections) {
  bool verbose = false;
if(verbose) printf("DEBUG: checking for collisions...\n");
  XYZ centre_abc(0.5,0.5,0.5);
  XYZ centre_xyz = framework->abc_to_xyz_returning_XYZ(centre_abc);
  bool collision = false;
  int num_assembled_molecules = assembled_molecules->size();
  int num_c = two_way_connections->size();
  //is it possible to have more than one periodic bond between atoms? depends on unit cell - we need to check all x,y,z in the range -1 to 1, except for redundant combinations, handled below
  bool multiple_periodic_connections_possible = false;
  for(int z=-1; z<2 && !multiple_periodic_connections_possible; z++) {
    for(int y=-1; y<2 && !multiple_periodic_connections_possible; y++) {
      for(int x=0; x<2 && !multiple_periodic_connections_possible; x++) { //x cannot be -1
        if((x==0 && y==0 && z==1) || (x==0 && y==1) || (x==1)) {
          XYZ image_abc(x,y,z);
          double self_distance = framework->abc_to_xyz_returning_XYZ(image_abc).magnitude();
          if(self_distance<=GENERAL_BOND_DISTANCE*2.0) multiple_periodic_connections_possible = true;
        }
      }
    }
  }

  //new and much faster system - build a set of vectors storing the info we need and just query it according to the basic collision rules - and don't include sites or H atoms
  vector <int> atom_molecule_id_vector;
  vector <int> atom_position_id_vector;
  vector <XYZ> atom_abc_vector;
  vector <bool> atom_is_close_to_site_vector;
/*
int debug1=198, debug2=730;
*/
  for(int i=0; i<num_assembled_molecules; i++) { //for each molecule
    int num_at_i = assembled_molecules->at(i).atoms_xyz.size();
    int num_sites_i = assembled_molecules->at(i).sites.size();
    int num_dummy_sites_i = assembled_molecules->at(i).dummy_sites.size();
    vector <XYZ> local_atom_abc_vector;
    vector <XYZ> local_site_abc_vector;
    for(int j=0; j<num_at_i; j++) { //for each atom
      if(assembled_molecules->at(i).atoms_type.at(j)!="H") { //if not H
        XYZ xyz_j = assembled_molecules->at(i).atoms_xyz.at(j);
        XYZ abc_j = framework->xyz_to_abc_returning_XYZ(xyz_j);
        bool is_site_j = false;
        bool is_dummy_site_j = false;
        for(int s=0; s<num_sites_i && !is_site_j; s++) { //find if site
          is_site_j = (j == assembled_molecules->at(i).sites.at(s));
        }
        for(int s=0; s<num_dummy_sites_i && !is_dummy_site_j && !is_site_j; s++) { //find if dummy site
          is_dummy_site_j = (j == assembled_molecules->at(i).dummy_sites.at(s));
        }
        if(is_site_j) {
          local_site_abc_vector.push_back(abc_j);
        } else if(is_dummy_site_j) { //do nothing - don't want to write this one out, and don't want to bypass collision detection near to dummy sites
        } else {
          local_atom_abc_vector.push_back(abc_j);
          atom_abc_vector.push_back(abc_j);
          atom_molecule_id_vector.push_back(i);
          atom_position_id_vector.push_back(j);
        }
      }
    }
    int num_actual_at_i = local_atom_abc_vector.size();
    for(int j=0; j<num_actual_at_i; j++) {
      bool is_close_to_site = false;
      for(int k=0; k<num_sites_i && !is_close_to_site; k++) {
        is_close_to_site = is_part_of_site_abc(local_atom_abc_vector.at(j), local_site_abc_vector.at(k), framework);
/*
if(atom_is_close_to_site_vector.size()==debug1 || atom_is_close_to_site_vector.size()==debug2) {
  if(is_close_to_site) printf("DEBYG: overall atom index %d (fractional %.3f/%.3f/%.3f) is close to site index %d (fractional %.3f/%.3f/%.3f) in molecule index %d\n", (int)(atom_is_close_to_site_vector.size()), local_atom_abc_vector.at(j).x, local_atom_abc_vector.at(j).y, local_atom_abc_vector.at(j).z, k, local_site_abc_vector.at(k).x, local_site_abc_vector.at(k).y, local_site_abc_vector.at(k).z, i);
  else printf("DEBYG: overall atom index %d (fractional %.3f/%.3f/%.3f) is NOT close to site index %d (fractional %.3f/%.3f/%.3f) in molecule index %d (distance = %.3f)\n", (int)(atom_is_close_to_site_vector.size()), local_atom_abc_vector.at(j).x, local_atom_abc_vector.at(j).y, local_atom_abc_vector.at(j).z, k, local_site_abc_vector.at(k).x, local_site_abc_vector.at(k).y, local_site_abc_vector.at(k).z, i, framework->getDistCalc().minimum_periodic_distance(local_atom_abc_vector.at(j).x, local_atom_abc_vector.at(j).y, local_atom_abc_vector.at(j).z, local_site_abc_vector.at(k).x, local_site_abc_vector.at(k).y, local_site_abc_vector.at(k).z));
}
*/
      }
      atom_is_close_to_site_vector.push_back(is_close_to_site);
    }
  }
  int num_total_atoms = atom_abc_vector.size();
  for(int i=0; i<num_total_atoms && !collision; i++) {
    for(int j=i+1; j<num_total_atoms && !collision; j++) {
      if(atom_molecule_id_vector.at(i)==atom_molecule_id_vector.at(j)) { //atoms are in the same molecule
        if(!bonded_xyz(assembled_molecules->at(atom_molecule_id_vector.at(i)).atoms_xyz.at(atom_position_id_vector.at(i)), assembled_molecules->at(atom_molecule_id_vector.at(j)).atoms_xyz.at(atom_position_id_vector.at(j)))) { //if they are not supposed to be bonded
          if(bonded_abc(atom_abc_vector.at(i), atom_abc_vector.at(j), framework)) {
            if(verbose) printf("NOTICE: non-bonded atoms in building block %d are bonded periodically, constituting a collision\n", atom_molecule_id_vector.at(i));
            collision = true; //collision if they are bonded
          }
        } else if(multiple_periodic_connections_possible) { //else they are supposed to be bonded - so we need to make sure that they are not bonded more than once periodically, if that is a possibility given the unit cell parameters
          XYZ shift_i_to_centre = centre_abc - trans_to_origuc(atom_abc_vector.at(i));
          XYZ uc_j = trans_to_origuc(atom_abc_vector.at(j)+shift_i_to_centre);
          int num_periodic_bonds = 0;
          for(int da=-1; da<2 && num_periodic_bonds<2; da++) {
            for(int db=-1; db<2 && num_periodic_bonds<2; db++) {
              for(int dc=-1; dc<2 && num_periodic_bonds<2; dc++) {
                XYZ image(da,db,dc);
                XYZ image_j_xyz = framework->abc_to_xyz_returning_XYZ(uc_j+image);
                if(bonded_xyz(image_j_xyz, centre_xyz)) num_periodic_bonds++;
              }
            }
          }
          if(num_periodic_bonds>1) {
            if(verbose) printf("NOTICE: more than one (%d) bonds were detected between atoms in building block %d, constituting a collision\n", num_periodic_bonds, atom_molecule_id_vector.at(i));
            collision = true;
          }
        }
      } else if(bonded_abc(atom_abc_vector.at(i), atom_abc_vector.at(j), framework)) { //atoms are in different molecules, and they are bonded
        if(!atom_is_close_to_site_vector.at(i) || !atom_is_close_to_site_vector.at(j)) {
          if(verbose) printf("NOTICE: a bond was detected between atoms that are not both part of a site: %s (fractional %.3f/%.3f/%.3f) and %s (fractional %.3f/%.3f/%.3f) in building blocks %d and %d respectively (overall indices %d and %d), constituting a collision (distance = %.3f)\n", assembled_molecules->at(atom_molecule_id_vector.at(i)).atoms_type.at(atom_position_id_vector.at(i)).c_str(), atom_abc_vector.at(i).x, atom_abc_vector.at(i).y, atom_abc_vector.at(i).z, assembled_molecules->at(atom_molecule_id_vector.at(j)).atoms_type.at(atom_position_id_vector.at(j)).c_str(), atom_abc_vector.at(j).x, atom_abc_vector.at(j).y, atom_abc_vector.at(j).z, atom_molecule_id_vector.at(i), atom_molecule_id_vector.at(j), i, j, framework->getDistCalc().minimum_periodic_distance(atom_abc_vector.at(i).x, atom_abc_vector.at(i).y, atom_abc_vector.at(i).z, atom_abc_vector.at(j).x, atom_abc_vector.at(j).y, atom_abc_vector.at(j).z));

          collision = true; //if either atom is not part of a site, this is a collision
        } else { //else they are both part of sites - so we need to check whether these molecules are supposed to be connected
          bool mols_connected = false;
          for(int k=0; k<num_c && !mols_connected; k++) {
            CONNECTION c = two_way_connections->at(k);
            int v1 = c.v1;
            int v2 = c.v2;
            if((v1==atom_molecule_id_vector.at(i) && v2==atom_molecule_id_vector.at(j)) || (v2==atom_molecule_id_vector.at(i) && v1==atom_molecule_id_vector.at(j))) {
              mols_connected = true;
            }
          }
          if(!mols_connected) {
            if(verbose) printf("NOTICE: a bond was detected between disconnected building blocks %d and %d, constituting a collision\n", atom_molecule_id_vector.at(i), atom_molecule_id_vector.at(j));
            collision = true;
          }
        }
      }
    }
  }
  return collision;
}

/* take a molecule, and return a new MOLECULE after applying a specified symmetry operation */
MOLECULE apply_symmetry_operator(MOLECULE orig, int sym_op, int sym_ID, ATOM_NETWORK *underlying_net) {
  MOLECULE rot = orig;
  int num_atom = orig.atoms_xyz.size();
  for(int i=0; i<num_atom+1; i++) {
    XYZ pos;
    if(i<num_atom) pos = orig.atoms_xyz.at(i); else pos = orig.com; //handle the centre of mass too!
    XYZ abc_pos = underlying_net->xyz_to_abc_returning_XYZ(pos);
    XYZ new_pos = underlying_net->abc_to_xyz_returning_XYZ(GetEquivalentPositions(sym_ID, &abc_pos).at(sym_op));
    if(i<num_atom) rot.atoms_xyz.at(i) = new_pos; else rot.com = new_pos; //handle the centre of mass too!
  }
  return rot;
}

/* parse a cgd (topology, Systre format) file */
bool read_cgd(FILE *cgd, ATOM_NETWORK *cell, string *net) {
  bool verbose = false; //write parsing progress to terminal
  if(verbose) printf(" --- called read_cgd() --- \n");
  char *line = new char[maxline];
  vector<string> token;
  int line_num = 0;
  bool read_ok = false;
  int atom_index = 0;
  bool uses_atom = false;
  bool uses_node = false; //two different kinds of file format use these different flags, and order/format the data differently ...
  //flags for when line breaks separate the data to be read
  bool anticipating_name = false, anticipating_cell = false, anticipating_group = false, anticipating_atom = false, anticipating_node = false, anticipating_edge = false;
  //reading loop begins here
  while(fgets(line, maxline, cgd)!=NULL && !read_ok) {
    if(verbose) printf("read line \"%s\"\n", line);
    line_num++;
    token = split(string(line)," ()\r\t\n");
    int num_tokens = token.size();
    if(num_tokens==0) {
      //do nothing on this empty line, just skip to next line
    }
    else if(token.at(0)=="name" || token.at(0)=="NAME" || token.at(0)=="id" || token.at(0)=="ID") {
      if(num_tokens>1) { //read from this line if there are fields, else read from next line...
        (*net) = token.at(1); //returning
      } else anticipating_name = true;
      if(verbose) printf("parsed topology name %s\n", net->c_str());
    }
    else if(anticipating_name) {
      (*net) = token.at(0); //returning
      anticipating_name = false;
    }
    else if(token.at(0)=="cell" || token.at(0)=="CELL") {
      if(num_tokens>6) { //read from this line if there are fields, else read from next line...
        parse_cell(&token, 1, cell);
      } else anticipating_cell = true;
    }
    else if(anticipating_cell && num_tokens>5) { //right now potential line breaks between cell definition and data fields are only handled when all the data fields are on the same line
      parse_cell(&token, 0, cell);
      anticipating_cell = false;
    }
    else if(token.at(0)=="group" || token.at(0)=="GROUP") {
      if(num_tokens>1) { //read from this line if there are fields, else read from next line...
        parse_group(&token, 1, cell);
      } else anticipating_group = true;
    }
    else if(anticipating_group) {
      parse_group(&token, 0, cell);
      anticipating_group = false;
    }
    else if(token.at(0)=="atom" || token.at(0)=="ATOM") {
      uses_atom = true;
      if(num_tokens>5) { //read from this line if there are fields, else read from next line...
        parse_atom(&token, 1, cell, &atom_index, &line_num, line, cgd);
      } else anticipating_atom = true;
    }
    else if(anticipating_atom && num_tokens>4) { //right now potential line breaks between atom definition and data fields are only handled when all the data fields are on the same line
      parse_atom(&token, 0, cell, &atom_index, &line_num, line, cgd);
      anticipating_atom = false;
    }
    else if(token.at(0)=="node" || token.at(0)=="NODE") { //a different kind of atom, does not list nodes directly
      uses_node = true;
      if(num_tokens>5) { //read from this line if there are fields, else read from next line...
        parse_node(&token, 1, cell, &atom_index);
      } else anticipating_node = true;
    }
    else if(anticipating_node && num_tokens>4) { //right now potential line breaks between node definition and data fields are only handled when all the data fields are on the same line
      parse_node(&token, 0, cell, &atom_index);
      anticipating_node = false;
    }
    else if(token.at(0)=="edge" || token.at(0)=="EDGE") { //should only occur when the different kind of atom (node) is used - handled differently to edges above
      if(num_tokens>6) { //read from this line if there are fields, else read from next line...
        parse_edge(&token, 1, cell, atom_index, uses_node);
      } else anticipating_edge = true;
    }
    else if(anticipating_edge && num_tokens>5) {
      parse_edge(&token, 0, cell, atom_index, uses_node);
      anticipating_edge = false;
    }
    else if(token.at(0)=="#") { //comment - ignore
      //do nothing
    }
    else if(token.at(0)=="end" || token.at(0)=="END") {
      if(verbose) printf(" --- successfully finished parsing net file --- \n");
      read_ok = true;
    }
    if(uses_node && uses_atom) {
      printf("NET ERROR: detected that both atom and node flags are used in the input net file - this is currently assumed to indicate an invalid input file\n");
      exit(EXIT_FAILURE);
    }
  }
  if(!read_ok) {
    printf("NET WARNING: net file parsing ended before \"end\" line was read\n");
  } else printf("net file parsed correctly\n");
  int num_v = cell->vertices.size();
  bool edges_provided_exhaustively = true;
  for(int i=0; i<num_v && edges_provided_exhaustively; i++) {
    if(cell->vertices.at(i).expected_edges!=cell->vertices.at(i).edges.size()) edges_provided_exhaustively = false;
  }
  delete[] line;
  return edges_provided_exhaustively;
}

void parse_cell(vector<string> *token, int first_index, ATOM_NETWORK *cell) {
  bool verbose = false;
  cell->make(convertToDouble(token->at(first_index+0)), convertToDouble(token->at(first_index+1)), convertToDouble(token->at(first_index+2)), convertToDouble(token->at(first_index+3)), convertToDouble(token->at(first_index+4)), convertToDouble(token->at(first_index+5)));
  if(verbose) printf("parsed cell params %.3f %.3f %.3f with angles %.3f %.3f %.3f\n", cell->a, cell->b, cell->c, cell->alpha, cell->beta, cell->gamma);
}

void parse_atom(vector<string> *token, int first_index, ATOM_NETWORK *cell, int *atom_index, int *line_num, char *line, FILE *cgd) {
  bool verbose = false;
  VERTEX v(convertToDouble(token->at(first_index+2)), convertToDouble(token->at(first_index+3)), convertToDouble(token->at(first_index+4)));
  if(verbose) printf("parsed atom labelled %d at %.3f %.3f %.3f\n", convertToInt(token->at(first_index+0)), v.abc.x, v.abc.y, v.abc.z);
  v.expected_edges = convertToInt(token->at(first_index+1)); //now parse the edges
  for(int i=0; i<v.expected_edges; i++) {
    if(fgets(line, maxline, cgd)!=NULL) {
      (*line_num)++;
      vector<string> local_token = split(string(line)," ()\r\t");
      if(local_token.at(0)=="edge" || local_token.at(0)=="EDGE") {
        XYZ e(convertToDouble(local_token.at(2)), convertToDouble(local_token.at(3)), convertToDouble(local_token.at(4)));
        v.edges.push_back(e);
        if(verbose) printf("parsed edge labelled %d at %.3f %.3f %.3f\n", convertToInt(local_token.at(1)), e.x, e.y, e.z);
      }
      else {
        printf("NET ERROR: atom with %d edges was declared but the string \"%s\" was read instead of the data for edge ID %d\n", v.expected_edges, local_token.at(0).c_str(), i);
        exit(EXIT_FAILURE);
      }
    } else {
      printf("NET WARNING: %d edges were expected but file ended after reading %d edges\n", v.expected_edges, i+1);
      //exit(EXIT_FAILURE);
    }
  }
  //here we can now explicitly read in dummy edges, if we read a 2c vertex
  if(v.expected_edges==2) { //i.e., if 2c
    if(fgets(line, maxline, cgd)!=NULL) {
      (*line_num)++;
      vector<string> local_token = split(string(line)," ()\r\t");
      if(local_token.at(0)=="dummy_edge" || local_token.at(0)=="DUMMY_EDGE") {
        XYZ e(convertToDouble(local_token.at(2)), convertToDouble(local_token.at(3)), convertToDouble(local_token.at(4)));
        v.dummy_edges.push_back(e);
        if(verbose) printf("parsed dummy edge labelled %d at %.3f %.3f %.3f\n", convertToInt(local_token.at(1)), e.x, e.y, e.z);
      }
      else {
        printf("NET ERROR: dummy edge for 2c atom with index %d was declared but the string \"%s\" was read instead of the data for the dummy edge\n", (*atom_index), local_token.at(0).c_str());
        exit(EXIT_FAILURE);
      }
    } else {
      printf("NET WARNING: dummy edge was expected for 2c atom with index %d, but file ended instead\n", (*atom_index));
      //exit(EXIT_FAILURE);
    }    
  }
  cell->vertices.push_back(v);
  cell->vertex_symmetry_operators.push_back(0); //the vertices read from file are the basic ones, before any symmetry operator (or you could say after operator 0, which does nothing)
  cell->vertex_basic_indices.push_back((*atom_index));
  (*atom_index)++;
}

void parse_node(vector<string> *token, int first_index, ATOM_NETWORK *cell, int *atom_index) {
  bool verbose = false;
  VERTEX v(convertToDouble(token->at(first_index+2)), convertToDouble(token->at(first_index+3)), convertToDouble(token->at(first_index+4)));
  if(verbose) printf("parsed node labelled %d at %.3f %.3f %.3f\n", convertToInt(token->at(first_index+0)), v.abc.x, v.abc.y, v.abc.z);
  v.expected_edges = convertToInt(token->at(first_index+1)); //this many edges are expected, but note that with this file format, the edges only come once all nodes are read and, in all likelihood, are not given exhaustively
  cell->vertices.push_back(v);
  cell->vertex_symmetry_operators.push_back(0); //the vertices read from file are the basic ones, before any symmetry operator (or you could say after operator 0, which does nothing)
  cell->vertex_basic_indices.push_back((*atom_index));
  (*atom_index)++;
}

void parse_group(vector<string> *token, int first_index, ATOM_NETWORK *cell) {
  bool verbose = false;
  string sym_name = token->at(first_index+0);
  if(verbose) printf("parsed group name %s\n", sym_name.c_str());
  int sym_ID = get_sym_ID(sym_name);
  if(verbose) printf("\ti.e. group ID %d\n", sym_ID);
  cell->sym_ID = sym_ID;
  cell->sym_name = sym_name;
}

void parse_edge(vector<string> *token, int first_index, ATOM_NETWORK *cell, int atom_index, bool uses_node) {
  bool verbose = false;
  if(!uses_node) {
    printf("NET ERROR: parsed unexpected \"edge\" field, when edges where expected to be provided with atoms\n");
    exit(EXIT_FAILURE);
  }
  XYZ start(convertToDouble(token->at(first_index+0)), convertToDouble(token->at(first_index+1)), convertToDouble(token->at(first_index+2)));
  XYZ start_xyz = cell->abc_to_xyz_returning_XYZ(start);
  XYZ end(convertToDouble(token->at(first_index+3)), convertToDouble(token->at(first_index+4)), convertToDouble(token->at(first_index+5)));
  XYZ end_xyz = cell->abc_to_xyz_returning_XYZ(end);
  //1) detect which vertex it starts at
  int v_index = -1;
  for(int i=0; i<atom_index && v_index==-1; i++) {
    if((start_xyz-cell->abc_to_xyz_returning_XYZ(cell->vertices.at(i).abc)).magnitude()<DISTANCE_TOLERANCE) { //match to parsed vertex
      v_index = i;
    }
  }
  if(v_index==-1) {
    printf("NET ERROR: could not match this edge start position to a previously parsed vertex\n");
    exit(EXIT_FAILURE);
  }
  //2) update the cell accordingly
  cell->vertices.at(v_index).edges.push_back(end);
  if(verbose) printf("parsed edge assigned to vertex %d at %.3f %.3f %.3f\n", v_index+1, end.x, end.y, end.z);
  //3) detect which vertex it ends at
  v_index = -1;
  for(int i=0; i<atom_index && v_index==-1; i++) {
    if((end_xyz-cell->abc_to_xyz_returning_XYZ(cell->vertices.at(i).abc)).magnitude()<DISTANCE_TOLERANCE) { //match to parsed vertex
      v_index = i;
    }
  }
  if(v_index==-1) {
/*
    printf("ERROR: could not match this edge end position to a previously parsed vertex\n");
    exit(EXIT_FAILURE);
*/
    if(verbose) printf("DEBUG: could not match this edge end position to a previously parsed vertex - creating an orphan edge to temporarily store this start/end until its true position by symmetry is known\n");
    //this is a specific scenario - the cgd file contains a definition of an edge, for which the end point is not a basic vertex - normally this is the default case, but for this cgd format the edges are not given exhaustively and so we need to know of each edge in both directions in order to construct the complete net including symmetry operations upon these edges
    //accordingly, we need to store this orphan edge, by saving its start and end in the cell as an orphan - this vertex is the symmetry image of some basic vertex, and we will assign this orphaned edge to that vertex when we know it
    cell->orphan_edge_starts.push_back(end);
    cell->orphan_edge_ends.push_back(start);
  } else {
    //4) update the cell accordingly
    cell->vertices.at(v_index).edges.push_back(start);
    if(verbose) printf("parsed edge assigned to vertex %d at %.3f %.3f %.3f\n", v_index+1, start.x, start.y, start.z);
  }
}

void add_missing_edges(ATOM_NETWORK *cell) { //this function is called when the vertices in a cell do not have the correct number of edges, due to edges not being provided exhaustively in the input file - we need to add them in by looking at symmetry!
  bool verbose = false;
  int num_vertices = cell->vertices.size();
  for(int i=0; i<num_vertices; i++) {
    VERTEX v = cell->vertices.at(i);
    XYZ vertex_xyz = cell->abc_to_xyz_returning_XYZ(v.abc);
    int num_edges = (int)(v.edges.size());
    int num_expect_edges = v.expected_edges;
if(verbose) printf("\n");
if(verbose) printf("DEBUG: vertex ID %d at %.3f %.3f %.3f [%.3f %.3f %.3f] has %d edges, %d were expected\n", i, v.abc.x, v.abc.y, v.abc.z, vertex_xyz.x, vertex_xyz.y, vertex_xyz.z, num_edges, num_expect_edges);
    vector<double> read_edge_lengths;
    if(num_edges<num_expect_edges) {
if(verbose) printf("DEBUG: BASIC VERTEX ID %d has fewer edges than expected - filling in missing edges by symmetry operations\n", i);
if(verbose) printf("\tDEBUG: edges read were as follows:\n");
      vector<XYZ> new_edges_abc, new_edges_xyz; //a new vector of edges - fill it up and hopefully it will be the right length, once duplicates are removed
      for(int j=0; j<num_edges; j++) {
        XYZ this_read_edge = v.edges.at(j); //the jth edge that was read - add it to the new vector if it is not a duplicate of one already added
        XYZ edge_xyz = cell->abc_to_xyz_returning_XYZ(this_read_edge);
        double read_edge_length = (edge_xyz-vertex_xyz).magnitude();
if(verbose) printf("\tDEBUG: %.3f %.3f %.3f [%.3f %.3f %.3f], length = %.3f\n", this_read_edge.x, this_read_edge.y, this_read_edge.z, edge_xyz.x, edge_xyz.y, edge_xyz.z, read_edge_length);
        int num_read_edge_lengths = read_edge_lengths.size();
        bool may_be_new_edge_length = true;
        for(int k=0; k<num_read_edge_lengths && may_be_new_edge_length; k++) {
          if(fabs(read_edge_length-read_edge_lengths.at(k))>DISTANCE_TOLERANCE) may_be_new_edge_length = false;
        }
        if(may_be_new_edge_length) read_edge_lengths.push_back(read_edge_length);
        bool this_read_edge_is_dupe = false;
        int num_new_edges = new_edges_abc.size();
        for(int k=0; k<num_new_edges && !this_read_edge_is_dupe; k++) {
          if(overlaps_xyz(edge_xyz, new_edges_xyz.at(k))) this_read_edge_is_dupe = true; //we should not reject edges because their end points overlap periodically - think pcu! they are valid unless they overlap in xyz!
        }
        if(!this_read_edge_is_dupe) { //the jth read edge is not a duplicate - add it to the vector
          new_edges_abc.push_back(this_read_edge);
          new_edges_xyz.push_back(edge_xyz);
        }
        //now, identify all the images of the jth edge by symmetry, and periodic images thereof (!), and 1) check that the length is correct, and 2) perform the same duplication check - note that even if j was a dupe, we are checking its images here - it should be redundant but just in case...
        vector<XYZ> equiv = GetEquivalentPositions(cell->sym_ID,&this_read_edge);
        int num_sym_images = equiv.size();
if(verbose) printf("\tDEBUG: num_sym_images of this edge = %d\n", num_sym_images);
        for(int k=0; k<num_sym_images; k++) {
          XYZ this_base_sym_edge = trans_to_origuc(equiv.at(k)); //the kth symmetry image of the jth edge that was read, translated to the unit cell - look at each periodic image, and add them to the new vector if they are not a duplicate of one already added
          vector<XYZ> periodic_images = get_periodic_images_of_uc_abc_position(this_base_sym_edge);
          int num_periodic_images = periodic_images.size();
          for(int im=0; im<num_periodic_images; im++) {
            XYZ this_sym_edge = periodic_images.at(im);
            XYZ sym_edge_xyz = cell->abc_to_xyz_returning_XYZ(this_sym_edge);
            double sym_edge_length = (sym_edge_xyz-vertex_xyz).magnitude();
            bool length_ok = false;
            if(fabs(sym_edge_length-read_edge_length)<DISTANCE_TOLERANCE) length_ok = true;
            if(length_ok) {
              bool this_sym_edge_is_dupe = false;
              num_new_edges = new_edges_abc.size();
if(verbose) printf("\t\tDEBUG: edge at %.3f %.3f %.3f [%.3f %.3f %.3f] length is ok - checking for overlap against the existing %d edges...\n", this_sym_edge.x, this_sym_edge.y, this_sym_edge.z, sym_edge_xyz.x, sym_edge_xyz.y, sym_edge_xyz.z, num_new_edges);
              for(int l=0; l<num_new_edges && !this_sym_edge_is_dupe; l++) {
                //if(overlaps_abc(this_sym_edge, new_edges.at(l), cell)) this_sym_edge_is_dupe = true;
if(verbose) printf("\t\t\tDEBUG: checking against existing edge %d at %.3f %.3f %.3f [%.3f %.3f %.3f]...\n", l, new_edges_abc.at(l).x, new_edges_abc.at(l).y, new_edges_abc.at(l).z, new_edges_xyz.at(l).x, new_edges_xyz.at(l).y, new_edges_xyz.at(l).z);
                if(overlaps_xyz(sym_edge_xyz, new_edges_xyz.at(l))) this_sym_edge_is_dupe = true; //we should not reject edges because their end points overlap periodically - think pcu! they are valid unless they overlap in xyz!
              }
              if(!this_sym_edge_is_dupe) { //the kth symmetry image of the jth read edge is not a duplicate - add it to the vector
                new_edges_abc.push_back(this_sym_edge);
                new_edges_xyz.push_back(sym_edge_xyz);
if(verbose) printf("\t\t\tDEBUG: this sym edge is unique\n");
              } else {
if(verbose) printf("\t\t\tDEBUG: this sym edge is a duplicate!\n");
              }
            }
          }
        }
      }
      int num_new_edges = new_edges_abc.size();
if(verbose) printf("\tDEBUG: after examining symmetry, the following %d edges were identified:\n", num_new_edges);
      for(int j=0; j<num_new_edges; j++) {
        XYZ this_new_edge_abc = new_edges_abc.at(j);
        XYZ this_new_edge_xyz = new_edges_xyz.at(j);
if(verbose) printf("\t\tDEBUG: %.3f %.3f %.3f [%.3f %.3f %.3f]\n", this_new_edge_abc.x, this_new_edge_abc.y, this_new_edge_abc.z, this_new_edge_xyz.x, this_new_edge_xyz.y, this_new_edge_xyz.z);
      }
      if(num_new_edges==num_expect_edges) {
if(verbose) printf("DEBUG: the correct number of edges, %d, have been identified from symmetry and periodic operations\n", num_new_edges);
        v.edges = new_edges_abc;
        cell->vertices.at(i) = v;
      } else {
if(verbose) printf("DEBUG: %d edges were expected but after examining edges by symmetry, %d were identified - can not complete the net with just this information - need to consider any orphan vertices\n", num_expect_edges, num_new_edges);
        //this is where we need the orphan edges - they may start at symmetry images of this vertex, and if they do, we can use the orphan edge here
        int num_orphan_edges = cell->orphan_edge_starts.size();
        int num_orphan_edges_check = cell->orphan_edge_ends.size();
        if(num_orphan_edges_check!=num_orphan_edges) {
          printf("ERROR: attempting to interpret orphan edges in order to complete the basic net, but the number of start points (%d) is not equal to the number of end points (%d)\n", num_orphan_edges, num_orphan_edges_check);
          exit(EXIT_FAILURE);
        }
if(verbose) printf("DEBUG: there are %d orphan edges to consider!\n", num_orphan_edges);
        for(int j=0; j<num_orphan_edges; j++) {
          XYZ orphan_vertex = cell->orphan_edge_starts.at(j); //the start is treated as a vertex here
          vector<XYZ> equiv_orphan_vertices = GetEquivalentPositions(cell->sym_ID,&orphan_vertex);
          int num_equiv_orphan_vertices = equiv_orphan_vertices.size();
          for(int k=0; k<num_equiv_orphan_vertices; k++) { //this orphan vertex does not overlap with this basic vertex - however, it might be equivalent by symmetry, and if it is, we need to check the edge end points by symmetry too, in order to complete the basic cell
            XYZ this_sym_orphan_vertex = equiv_orphan_vertices.at(k);
            if(overlaps_abc(this_sym_orphan_vertex, v.abc, cell)) { //this image of the orphan vertex overlaps with the basic vertex - it is a candidate for a missing edge
if(verbose) printf("DEBUG: orphan vertex ID %d at %.3f %.3f %.3f has a symmetry image at %.3f %.3f %.3f which overlaps with this basic vertex at %.3f %.3f %.3f\n", j, orphan_vertex.x, orphan_vertex.y, orphan_vertex.z, this_sym_orphan_vertex.x, this_sym_orphan_vertex.y, this_sym_orphan_vertex.z, v.abc.x, v.abc.y, v.abc.z);
              XYZ orphan_edge = cell->orphan_edge_ends.at(j);
              vector<XYZ> equiv_orphan_edges = GetEquivalentPositions(cell->sym_ID,&orphan_edge); //get all the symmetry positions of this orphan edge, but we only need to look at the particular symmetry operation in question, hopefully
              XYZ this_sym_orphan_edge = equiv_orphan_edges.at(k);
if(verbose) printf("\tDEBUG: the corresponding edge end symmetry image is at %.3f %.3f %.3f\n", this_sym_orphan_edge.x, this_sym_orphan_edge.y, this_sym_orphan_edge.z);
              XYZ sym_orphan_edge_xyz = cell->abc_to_xyz_returning_XYZ(this_sym_orphan_edge);
              double sym_orphan_edge_length = (sym_orphan_edge_xyz-vertex_xyz).magnitude();
              bool length_ok = false;
              int num_read_edge_lengths = read_edge_lengths.size();
              for(int l=0; l<num_read_edge_lengths && !length_ok; l++) {
                if(fabs(sym_orphan_edge_length-read_edge_lengths.at(l))<DISTANCE_TOLERANCE) {
                  length_ok = true;
if(verbose) printf("\t\tDEBUG: length is ok - %.3f is similar to desired %.3f\n", sym_orphan_edge_length, read_edge_lengths.at(l));
                }
              }
              if(length_ok) {
                bool this_sym_orphan_edge_is_dupe = false;
                int num_new_edges = new_edges_abc.size();
                for(int l=0; l<num_new_edges && !this_sym_orphan_edge_is_dupe; l++) {
                  if(overlaps_xyz(sym_orphan_edge_xyz, new_edges_xyz.at(l))) this_sym_orphan_edge_is_dupe = true; //we should not reject edges because their end points overlap periodically - think pcu! they are valid unless they overlap in xyz!
                }
                if(!this_sym_orphan_edge_is_dupe) { //the kth symmetry image of the jth orphan edge is not a duplicate - add it to the vector
if(verbose) printf("\t\t\tDEBUG: this orphan edge is unique\n");
                  new_edges_abc.push_back(this_sym_orphan_edge);
                  new_edges_xyz.push_back(sym_orphan_edge_xyz);
                } else {
if(verbose) printf("\t\t\tDEBUG: this orphan edge is a duplicate!\n");
                }
              }
            }
          }
        }
        //now we have exhausted the orphan edges too - if the cell is not complete now, then there has to be a problem
        num_new_edges = new_edges_abc.size();
if(verbose) printf("DEBUG: after examining symmetry AND orphan edges, the following %d edges were identified:\n", num_new_edges);
        for(int j=0; j<num_new_edges; j++) {
          XYZ this_new_edge_abc = new_edges_abc.at(j);
          XYZ this_new_edge_xyz = new_edges_xyz.at(j);
if(verbose) printf("\tDEBUG: %.3f %.3f %.3f [%.3f %.3f %.3f]\n", this_new_edge_abc.x, this_new_edge_abc.y, this_new_edge_abc.z, this_new_edge_xyz.x, this_new_edge_xyz.y, this_new_edge_xyz.z);
        }
        if(num_new_edges==num_expect_edges) {
if(verbose) printf("DEBUG: the correct number of edges, %d, have been identified from symmetry and periodic operations, and orphan vertices\n", num_new_edges);
          v.edges = new_edges_abc;
          cell->vertices.at(i) = v;
        } else {
          printf("ERROR: %d edges were expected but after examining edges by symmetry AND orphan edges, %d were identified - can not complete the net\n", num_expect_edges, num_new_edges);
          exit(EXIT_FAILURE);
        }
      }
    } else if(num_edges==num_expect_edges) {
if(verbose) printf("DEBUG: the correct number of edges, %d, have been identified - no action required for this vertex\n", num_edges);
    } else {
      printf("ERROR: more edges than expected - error case in file parsing, probably a bug!\n");
      exit(EXIT_FAILURE);
    }
  }
}

void add_2c_vertices_and_normal_edges(ATOM_NETWORK *cell) { //just add the 2-c vertices and normal edges - dummy edges are handled later, once we know all the vertices by symmetry, and can make a better guess of the orthogonal direction of the 2-c vertices
  bool verbose = false;
  if(verbose) printf("augmenting underlying net with 2-c vertices\n");
  int num_v = cell->vertices.size();
  vector<XYZ> candidates;
  vector<XYZ> from_abc_vector, to_abc_vector; //, orthogonal_abc_vector;
  for(int i=0; i<num_v; i++) {
    XYZ v_abc = cell->vertices.at(i).abc;
    XYZ v_xyz = cell->abc_to_xyz_returning_XYZ(v_abc);
    int num_e = cell->vertices.at(i).edges.size();
    for(int j=0; j<num_e; j++) {
      XYZ e_abc = cell->vertices.at(i).edges.at(j);
      XYZ e_xyz = cell->abc_to_xyz_returning_XYZ(e_abc);
      XYZ this_edge_vector = get_vector_from_to(v_xyz, e_xyz);
      double this_edge_length = this_edge_vector.magnitude();
      XYZ this_candidate_abc = v_abc+(e_abc-v_abc).scale(0.5);
      XYZ this_candidate_xyz = cell->abc_to_xyz_returning_XYZ(this_candidate_abc);
      candidates.push_back(this_candidate_abc);
      from_abc_vector.push_back(v_abc);
      to_abc_vector.push_back(e_abc);
      //forget dummy edges for now - they are handled later
    }
  }
  int num_candidates = candidates.size();
  //find and remove any duplicates
  vector<bool> keep_vector;
  int num_new_vertices = 0;
  for(int i=0; i<num_candidates; i++) {
    //check against all previous candidates to see if we should reject this one
    bool keep = true;
    for(int j=0; j<i && keep; j++) {
      if(overlaps_abc(candidates.at(i),candidates.at(j),cell)) keep = false;
    }
    keep_vector.push_back(keep);
    if(keep) num_new_vertices++;
  }
  if(num_new_vertices>0) {
    //halve the edge lengths of all previous vertices to match these new short edges!
    for(int i=0; i<num_v; i++) {
      XYZ v_abc = cell->vertices.at(i).abc;
      int num_e = cell->vertices.at(i).edges.size();
      for(int j=0; j<num_e; j++) {
        XYZ e_abc = cell->vertices.at(i).edges.at(j);
        XYZ new_e_abc = v_abc+(e_abc-v_abc).scale(0.5);
        cell->vertices.at(i).edges.at(j) = new_e_abc;
      }
    }
    //set up final vertices
    for(int i=0; i<num_candidates; i++) {
      if(keep_vector.at(i)) {
        if(verbose) printf("accepted candidate ID %d: %.3f %.3f %.3f\n", i, candidates.at(i).x, candidates.at(i).y, candidates.at(i).z);
        VERTEX v(candidates.at(i).x, candidates.at(i).y, candidates.at(i).z);
        v.expected_edges = 2; //it is a 2-c vertex
        v.edges.push_back(from_abc_vector.at(i));
        v.edges.push_back(to_abc_vector.at(i));
        cell->vertices.push_back(v);
      }
    }
  }
}

void add_2c_dummy_edges(ATOM_NETWORK *basic_cell, ATOM_NETWORK *full_cell, vector<CONNECTION> *two_way_connections) { //the vertices are already in place - now we must add the dummy edges
  bool verbose = false;
  if(verbose) printf("DEBUG: adding dummy edges to 2-c vertices to represent orthogonal point of reference for molecular alignment\n");
  vector<int> basic_IDs_examined;
  int num_two_way_connections = two_way_connections->size();
  int num_full_v = full_cell->vertices.size();
  for(int i=0; i<num_full_v; i++) {
    VERTEX v = full_cell->vertices.at(i);
    XYZ v_xyz = full_cell->abc_to_xyz_returning_XYZ(v.abc);
    int num_e = v.edges.size();
    if(num_e==2) { //i.e. if a 2-c vertex, we need to add a dummy orthogonal edge
      int basic_ID = full_cell->vertex_basic_indices.at(i);
      bool examined_already = false;
      int num_examined_already = basic_IDs_examined.size();
      int sym_op = full_cell->vertex_symmetry_operators.at(i);
      for(int ex=0; ex<num_examined_already && !examined_already; ex++) {
        if(basic_IDs_examined.at(ex)==basic_ID) examined_already = true;
      }
      if(!examined_already) { //only need to find the direction if this vertex is not a symmetry image of another
        basic_IDs_examined.push_back(basic_ID);
        if(sym_op!=0) { //double check that the symmetry operator is zero (no sym operator), as it should be
          printf("ERROR: found a 2-c vertex ID %d based on underlying 2-c vertex ID %d, which has not already been examined, but the symmetry operator is non-zero (%d)\n", i, basic_ID, sym_op);
          exit(EXIT_FAILURE);
        }
        //we need to know which two vertices (in the full cell) this 2-c vertex is connected to, so that we can look at their other edges and decide a direction
        int v1 = -1, v2 = -1;
        int e1 = -1, e2 = -1;
        for(int j=0; j<num_two_way_connections && (v1==-1 || v2==-1); j++) {
          CONNECTION c = two_way_connections->at(j);
          if(i==c.v1) {
            //this connection links i with something
            if(v1==-1) {v1 = c.v2; e1 = c.e2;} else if(v2==-1) {v2 = c.v2; e2 = c.e2;}
          } else if(i==c.v2) {
            //this connection links something with i
            if(v1==-1) {v1 = c.v1; e1 = c.e1;} else if(v2==-1) {v2 = c.v1; e2 = c.e1;}
          }
        }
        if(v1==-1 || v2==-1) { //failed to identify the two adjacent full cell vertices
          printf("ERROR: could not identify the two full cell vertices connected to 2-c vertex %d\n", i);
          exit(EXIT_FAILURE);
        }
        XYZ v1_xyz = full_cell->abc_to_xyz_returning_XYZ(full_cell->vertices.at(v1).abc);
        XYZ v2_xyz = full_cell->abc_to_xyz_returning_XYZ(full_cell->vertices.at(v2).abc);
if(verbose) printf("DEBUG: the two full vertices which define this 2c vertex are v1, id %d at %.3f %.3f %.3f, and v2, id %d at %.3f %.3f %.3f\n", v1, v1_xyz.x, v1_xyz.y, v1_xyz.z, v2, v2_xyz.x, v2_xyz.y, v2_xyz.z);
if(verbose) printf("DEBUG: v_xyz (2c vertex position) = %.3f %.3f %.3f\n", v_xyz.x, v_xyz.y, v_xyz.z);
        XYZ this_edge_vector = get_vector_from_to(v_xyz, v1_xyz);
        double this_edge_length = this_edge_vector.magnitude(); //already half of the original net's edge length
        //now, look at the edges of these two adjacent vertices until we find a pair that are neither parallel to the edges of this vertex, nor mutually parallel, and use them to identify an appropriate orthogonal direction
        bool non_parallel_direction_found = false;
        bool parallel_direction_found = false;
        XYZ dummy_edge_abc; dummy_edge_abc.x = 1001; dummy_edge_abc.y = 1001; dummy_edge_abc.z = 1001; //dummy data
        int num_edges_v1 = full_cell->vertices.at(v1).edges.size();
        int num_edges_v2 = full_cell->vertices.at(v2).edges.size();
        for(int j=0; j<num_edges_v1 && !non_parallel_direction_found; j++) {
          if(j!=e1) { //don't select the edge that connects to the 2-c vertex in question ...
            XYZ other_edge_1_end_abc = full_cell->vertices.at(v1).edges.at(j);
            XYZ other_edge_1_end_xyz = full_cell->abc_to_xyz_returning_XYZ(other_edge_1_end_abc);
            XYZ other_edge_1_vector = get_vector_from_to(v1_xyz, other_edge_1_end_xyz);
            //determine whether this edge points in a different direction, if it does, we can project the end point and use that direction as the reference
            double angle_between_degrees = 360.0*this_edge_vector.angle_between(other_edge_1_vector)/(2.0*PI);
            double acute_angle_between_degrees = min(angle_between_degrees, 180.0-angle_between_degrees);
if(verbose) printf("DEBUG: considering edge id %d of v1 which points to position %.3f %.3f %.3f (%.3f degrees angular difference)\n", j, other_edge_1_end_xyz.x, other_edge_1_end_xyz.y, other_edge_1_end_xyz.z, acute_angle_between_degrees);
            if(acute_angle_between_degrees>15.0) { //acceptable angle deviation
              XYZ other_edge_1_projection_xyz = project_onto_line(other_edge_1_end_xyz, v_xyz, v1_xyz);
              XYZ projection_1_vector_norm = get_vector_from_to(other_edge_1_end_xyz, other_edge_1_projection_xyz).unit();
if(verbose) printf("DEBUG: projection vector is from %.3f %.3f %.3f and norm is %.3f %.3f %.3f\n", other_edge_1_projection_xyz.x, other_edge_1_projection_xyz.y, other_edge_1_projection_xyz.z, projection_1_vector_norm.x, projection_1_vector_norm.y, projection_1_vector_norm.z);
              for(int k=0; k<num_edges_v2 && !non_parallel_direction_found; k++) {
                if(k!=e2) { //don't select the edge that connects to the 2-c vertex in question ...
                  XYZ other_edge_2_end_abc = full_cell->vertices.at(v2).edges.at(k);
                  XYZ other_edge_2_end_xyz = full_cell->abc_to_xyz_returning_XYZ(other_edge_2_end_abc);
                  XYZ other_edge_2_vector = get_vector_from_to(v2_xyz, other_edge_2_end_xyz);
                  //determine whether this edge points in a different direction, if it does, we can project the end point and use that direction as the reference
                  double angle_between_degrees = 360.0*this_edge_vector.angle_between(other_edge_2_vector)/(2.0*PI);
                  double acute_angle_between_degrees = min(angle_between_degrees, 180.0-angle_between_degrees);
if(verbose) printf("DEBUG: considering edge id %d of v2 which points to position %.3f %.3f %.3f (%.3f degrees angular difference)\n", k, other_edge_2_end_xyz.x, other_edge_2_end_xyz.y, other_edge_2_end_xyz.z, acute_angle_between_degrees);
                  if(acute_angle_between_degrees>15.0) { //acceptable angle deviation
                    //NOTE: need to set new v_xyz, and project using this and v2_xyz, otherwise we project onto the wrong line (this line may, by periodic boundary conditions, differ)
                    XYZ new_v_xyz = full_cell->abc_to_xyz_returning_XYZ(full_cell->vertices.at(v2).edges.at(e2));
if(verbose) printf("DEBUG: new_v_xyz (image of 'basic' 2c vertex) = %.3f %.3f %.3f\n", new_v_xyz.x, new_v_xyz.y, new_v_xyz.z);
                    XYZ other_edge_2_projection_xyz = project_onto_line(other_edge_2_end_xyz, new_v_xyz, v2_xyz);
                    XYZ projection_2_vector_norm = get_vector_from_to(other_edge_2_end_xyz, other_edge_2_projection_xyz).unit();
if(verbose) printf("DEBUG: projection vector is from %.3f %.3f %.3f and norm is %.3f %.3f %.3f\n", other_edge_2_projection_xyz.x, other_edge_2_projection_xyz.y, other_edge_2_projection_xyz.z, projection_2_vector_norm.x, projection_2_vector_norm.y, projection_2_vector_norm.z);
                    //at this point, we now have two normalised directions that are not parallel with the 2-c vertex's edges - but are they parallel with each other?
                    angle_between_degrees = 360.0*projection_1_vector_norm.angle_between(projection_2_vector_norm)/(2.0*PI);
                    acute_angle_between_degrees = min(angle_between_degrees, 180.0-angle_between_degrees);
                    if(acute_angle_between_degrees>15.0) { //acceptable angle deviation
                      //not parallel! easy case - find the midpoint vector, and reverse it
                      XYZ dummy_vector_norm = (projection_1_vector_norm+projection_2_vector_norm).unit();
                      XYZ dummy_edge_xyz = v_xyz+dummy_vector_norm.scale(this_edge_length);
                      dummy_edge_abc = full_cell->xyz_to_abc_returning_XYZ(dummy_edge_xyz);
if(verbose) printf("DEBUG: non-parallel combo found: %.3f %.3f %.3f and %.3f %.3f %.3f leading to %.3f %.3f %.3f - with edge length of %.3f and a vector position %.3f %.3f %.3f, the edge ends at %.3f %.3f %.3f (in abc, %.3f %.3f %.3f)\n", projection_1_vector_norm.x, projection_1_vector_norm.y, projection_1_vector_norm.z, projection_2_vector_norm.x, projection_2_vector_norm.y, projection_2_vector_norm.z, dummy_vector_norm.x, dummy_vector_norm.y, dummy_vector_norm.z, this_edge_length, v_xyz.x, v_xyz.y, v_xyz.z, dummy_edge_xyz.x, dummy_edge_xyz.y, dummy_edge_xyz.z, dummy_edge_abc.x, dummy_edge_abc.y, dummy_edge_abc.z);
                      non_parallel_direction_found = true;
                    } else if(!parallel_direction_found) {
                      //else parallel; harder case - take the cross product of either, say 1; save one of these only, in case we do not find a parallel one
                      XYZ dummy_vector_norm = (this_edge_vector.unit().cross(projection_1_vector_norm)).unit();
                      XYZ dummy_edge_xyz = v_xyz+dummy_vector_norm.scale(this_edge_length);
                      dummy_edge_abc = full_cell->xyz_to_abc_returning_XYZ(dummy_edge_xyz);
if(verbose) printf("DEBUG: parallel combo found: %.3f %.3f %.3f and %.3f %.3f %.3f leading to %.3f %.3f %.3f - with edge length of %.3f and a vector position %.3f %.3f %.3f, the edge ends at %.3f %.3f %.3f (in abc, %.3f %.3f %.3f)\n", projection_1_vector_norm.x, projection_1_vector_norm.y, projection_1_vector_norm.z, projection_2_vector_norm.x, projection_2_vector_norm.y, projection_2_vector_norm.z, dummy_vector_norm.x, dummy_vector_norm.y, dummy_vector_norm.z, this_edge_length, v_xyz.x, v_xyz.y, v_xyz.z, dummy_edge_xyz.x, dummy_edge_xyz.y, dummy_edge_xyz.z, dummy_edge_abc.x, dummy_edge_abc.y, dummy_edge_abc.z);
                      parallel_direction_found = true;
                    }
                  }
                }
              }
            }
          }
        }
        //at this point we know the selected direction, and can implement it
        if(dummy_edge_abc.x>1000 && dummy_edge_abc.y>1000 && dummy_edge_abc.z>1000) {
          printf("ERROR: could not find a candidate orthogonal direction with which to orient full cell 2-c vertex %d's dummy edge\n", i);
          exit(EXIT_FAILURE);
        } else if(verbose) {
          printf("DEBUG: orthogonal position abc chosen to be %.3f %.3f %.3f\n", dummy_edge_abc.x, dummy_edge_abc.y, dummy_edge_abc.z);
        }
        //so update the full cell vertex with the new data, and store the basic vertex info locally only
        full_cell->vertices.at(i).dummy_edges.push_back(dummy_edge_abc);
        basic_cell->vertices.at(basic_ID).dummy_edges.push_back(dummy_edge_abc);
      } else { //else this vertex is a symmetry image of another - in that case we can use this info to assign the orientation based on that of the underlying vertex
        XYZ test_dummy_e = basic_cell->vertices.at(basic_ID).dummy_edges.at(0);
        vector<XYZ> equiv_dummy_e = GetEquivalentPositions(basic_cell->sym_ID,&test_dummy_e);
/*
        XYZ here_uc = trans_to_origuc(v.abc);
        XYZ shift = here_uc - v.abc;
        full_cell->vertices.at(i).dummy_edges.push_back(equiv_dummy_e.at(sym_op)+shift); //this is nice - we know that we are working on symmetric vertex image with index j - so we can just keep the corresponding edge image - shifting it to unit cell by the same amount that we used for the vertex!
*/
        full_cell->vertices.at(i).dummy_edges.push_back(equiv_dummy_e.at(sym_op)); //this is nice - we know that we are working on symmetric vertex image with index j - so we can just keep the corresponding edge image - shifting it to unit cell by the same amount that we used for the vertex!
      }
    }
  }
}

void write_abstract_cif(FILE *cif, ATOM_NETWORK *cell, string name) {
  //this function writes out a cif file with abstracted cell information
  //the basic vertices and edges are written, i.e. before symmetry operations, and the symmetry group is written to the file
  //edges are written as vertices - they are shortened to 1/4 length, so that they do not overlap with real vertices, and so that their direction can be preserved
  fprintf(cif, "#******************************************\n#\n# CIF file created by Richard L. Martin, Lawrence Berkeley National Laboratory, 2013\n#\n#*******************************************\n\n");
  fprintf(cif, "_cell_length_a\t\t%.3f(0)\n", cell->a);
  fprintf(cif, "_cell_length_b\t\t%.3f(0)\n", cell->b);
  fprintf(cif, "_cell_length_c\t\t%.3f(0)\n", cell->c);
  fprintf(cif, "_cell_angle_alpha\t\t%.3f(0)\n", cell->alpha);
  fprintf(cif, "_cell_angle_beta\t\t%.3f(0)\n", cell->beta);
  fprintf(cif, "_cell_angle_gamma\t\t%.3f(0)\n", cell->gamma);
//  fprintf(cif, "_symmetry_space_group_name_H-M\t\t%s\n", cell->sym_name.c_str());
fprintf(cif, "_symmetry_space_group_name_H-M\t\tP-1\n");
//  fprintf(cif, "_symmetry_Int_Tables_number\t\t%d\n", cell->sym_ID);
  fprintf(cif, "_symmetry_Int_Tables_number\t\t1\n");
  fprintf(cif, "_symmetry_cell_setting\t\t");
  
  //Determine the Crystal System
  if (cell->alpha == 90 && cell->beta == 90 && cell->gamma == 90){
    if (cell->a == cell->b || cell->b == cell->c || cell->a == cell->c){
      if (cell->a == cell->b && cell->b == cell->c){
        fprintf(cif, "Isometric\n\n");
      }
      else {
        fprintf(cif, "Tetragonal\n\n");
      }
    }
    else{
    	fprintf(cif, "Orthorhombic\n\n");
    }
  }
  else if(cell->alpha == cell->beta || cell->beta == cell->gamma || cell->alpha == cell->gamma){
    fprintf(cif, "Monoclinic\n\n");
  }
  else{
    fprintf(cif, "Triclinic\n\n");
  }
  
  fprintf(cif, "loop_\n");
  fprintf(cif, "_symmetry_equiv_pos_as_xyz\n");
fprintf(cif, "'+x,+y,+z'\n\n"); //can we do something smarter here?
  fprintf(cif, "loop_\n");
  fprintf(cif, "_atom_site_label\n");
  fprintf(cif, "_atom_site_type_symbol\n");
  fprintf(cif, "_atom_site_fract_x\n");
  fprintf(cif, "_atom_site_fract_y\n");
  fprintf(cif, "_atom_site_fract_z\n");
  int num_vertices = cell->vertices.size();
  int count = 1;
  for (int i=0; i<num_vertices; i++){
    VERTEX v = cell->vertices.at(i);
    int num_edges = v.edges.size();
    if(num_edges>2) { //don't write out augmented, 2-c vertices...
//      fprintf(cif, "%d\t%d\t%.6f\t%.6f\t%.6f\n", count, count, v.abc.x, v.abc.y, v.abc.z);
      fprintf(cif, "%s\t%s\t%.6f\t%.6f\t%.6f\n", "C", "C", v.abc.x, v.abc.y, v.abc.z);
      count++;
      for (int j=0; j<num_edges; j++) { //... or dummy edges...
        XYZ edge = v.edges.at(j);
        XYZ quarter_edge = ((edge-v.abc).scale(1.0/4.0))+v.abc;
//        fprintf(cif, "%d\t%d\t%.6f\t%.6f\t%.6f\n", count, count, quarter_edge.x, quarter_edge.y, quarter_edge.z);
        fprintf(cif, "%s\t%s\t%.6f\t%.6f\t%.6f\n", "H", "H", quarter_edge.x, quarter_edge.y, quarter_edge.z);
        count++;        
      }
    }
  }
}

void write_unit_cell(FILE *orig_uc, ATOM_NETWORK *cell, string name, bool is_net) {
  //printf("DEBUG: called write_unit_cell()\n");
  //For each corner in abc, get the coords in xyz
  vector<XYZ> corners; XYZ p;
  p = cell->abc_to_xyz_returning_XYZ(0, 0, 0); corners.push_back(p);
  p = cell->abc_to_xyz_returning_XYZ(0, 0, 1); corners.push_back(p);
  p = cell->abc_to_xyz_returning_XYZ(0, 1, 0); corners.push_back(p);
  p = cell->abc_to_xyz_returning_XYZ(0, 1, 1); corners.push_back(p);
  p = cell->abc_to_xyz_returning_XYZ(1, 0, 0); corners.push_back(p);
  p = cell->abc_to_xyz_returning_XYZ(1, 0, 1); corners.push_back(p);
  p = cell->abc_to_xyz_returning_XYZ(1, 1, 0); corners.push_back(p);
  p = cell->abc_to_xyz_returning_XYZ(1, 1, 1); corners.push_back(p);
  // Write header, information about the cell, and footer
  int num_c = corners.size();
  fprintf(orig_uc, "# vtk DataFile Version 2.0\n%s", name.c_str());
  if(is_net) fprintf(orig_uc, " - vtk format representation of original unit cell boundary\n"); else fprintf(orig_uc, " - vtk format representation of structure's unit cell boundary\n");
  fprintf(orig_uc, "ASCII\nDATASET POLYDATA\nPOINTS %d double\n", num_c);
  for(int i=0; i<num_c; i++) fprintf(orig_uc, "%.3f %.3f %.3f\n", corners.at(i)[0], corners.at(i)[1], corners.at(i)[2]);
  fprintf(orig_uc, "LINES 12 36\n2 0 1\n2 0 2\n2 1 3\n2 2 3\n2 4 5\n2 4 6\n2 5 7\n2 6 7\n2 0 4\n2 1 5\n2 2 6\n2 3 7\n");
}

void write_vertices(FILE *orig_v, ATOM_NETWORK *cell, string name, bool rename_vertices_by_connectivity, bool is_net) {
  //printf("DEBUG: called write_vertices()\n");
  //For each vertex in abc, get the coords in xyz
  vector<XYZ> vertices; XYZ p;
  //assign vertices names based on their connectivity, if desired
  if(rename_vertices_by_connectivity) {
    for(int i=0; i<cell->vertices.size(); i++) {
      char* temp_name = new char[100];
      sprintf(temp_name, "%d", cell->vertices.at(i).expected_edges);
//printf("renaming vertex: %s to %s\n", cell->vertices.at(i).name.c_str(), temp_name);
      cell->vertices.at(i).name = temp_name;
      delete[] temp_name;
    }
  }
  for(int i=0; i<cell->vertices.size(); i++) {
    p = cell->abc_to_xyz_returning_XYZ(cell->vertices.at(i).abc);
    vertices.push_back(p);
  }
  // Write header, and vertex coords
  int num_v = vertices.size();
  double radius = 1.0; //give particles a dummy radius
  fprintf(orig_v, "%d\n%s", num_v, name.c_str());
  if(is_net) fprintf(orig_v, " - xyz format representation of original vertices\n"); else fprintf(orig_v, " - xyz format representation of structure\n");
  for(int i=0; i<num_v; i++) fprintf(orig_v, "%s %.3f %.3f %.3f %.3f\n", cell->vertices.at(i).name.c_str(), vertices.at(i)[0], vertices.at(i)[1], vertices.at(i)[2], radius);
}

void write_edges(FILE *orig_e, ATOM_NETWORK *cell, string name) {
  bool verbose = false;
  //printf("DEBUG: called write_edges_dir()\n");
  //For each vertex in abc, get the coords in xyz, and find all edge end points
  vector<XYZ> vertices;
  XYZ v;
  vector< vector<XYZ> > edges;
  int total_num_e = 0;
  for(int i=0; i<cell->vertices.size(); i++) {
    XYZ vabc = cell->vertices.at(i).abc;
    v = cell->abc_to_xyz_returning_XYZ(vabc);
    vertices.push_back(v);
    int num_e = cell->vertices.at(i).edges.size();
    int num_dummy_e = cell->vertices.at(i).dummy_edges.size();
if(verbose) printf("DEBUG: this vertex has %d edges and %d dummy_edges\n", num_e, num_dummy_e);
    vector<XYZ> these_edges;
    for(int j=0; j<num_e; j++) {
      XYZ e;
      XYZ eabc = cell->vertices.at(i).edges.at(j);
      e = cell->abc_to_xyz_returning_XYZ(eabc); //treat the edge coords as positions
      these_edges.push_back(e);
      total_num_e++;
    }
    for(int j=0; j<num_dummy_e; j++) {
      XYZ e;
      XYZ eabc = cell->vertices.at(i).dummy_edges.at(j);
      e = cell->abc_to_xyz_returning_XYZ(eabc); //treat the edge coords as positions
      these_edges.push_back(e);
      total_num_e++;
    }
    edges.push_back(these_edges);
  }
  // Write header, information about the cell, and footer
  int num_v = vertices.size();
  fprintf(orig_e, "# vtk DataFile Version 2.0\n%s - vtk format representation of original edges\nASCII\nDATASET POLYDATA\nPOINTS %d double\n", name.c_str(), total_num_e+num_v);
  for(int i=0; i<num_v; i++) { //for each vertex ...
    fprintf(orig_e, "%.3f %.3f %.3f\n", vertices.at(i)[0], vertices.at(i)[1], vertices.at(i)[2]); //...print the vertex coords ...
    int num_e = edges.at(i).size();
    for(int j=0; j<num_e; j++) fprintf(orig_e, "%.3f %.3f %.3f\n", edges.at(i).at(j)[0], edges.at(i).at(j)[1], edges.at(i).at(j)[2]); //...and print the edge coords
  }
  //now this is trickier - print the edge connectivity with respect to vtk points
  fprintf(orig_e, "LINES %d %d\n", total_num_e, total_num_e*3);
  int v_index = 0;
  for(int i=0; i<num_v; i++) { //for each vertex ...
    int num_e = edges.at(i).size();
    for(int j=0; j<num_e; j++) { //...for each edge of this vertex ...
      fprintf(orig_e, "2 %d %d\n", v_index, v_index+j+1); //...print that this vertex is connected to each of its edges
    }
    v_index+=num_e+1;
  }
}

void read_xyz(FILE *input, MOLECULE *mol, const char *filename) {
  int length = 100;
  char ch1[length];
  int status = 0;
  int num_atoms = 0;
  if(fgets(ch1, length, input)!=NULL) {
    string str = string(ch1);
//    printf("DEBUG: read char array as \"%s\", parsed to string as \"%s\"\n", ch, str.c_str());
    int pos=0;
    char c = str[pos];

//@@@ WARNING - there is a memory leak here and I cannot tell why

    while(c<=0) { //for some reason, xyz files from marvinsketch have three characters at the beginning with ascii codes <0 - very strange, but we need to remove these to succesfully parse any data
      pos++;
      c = str[pos];
    }
    char ch2[length];
    str.copy(ch2, str.size()-pos, pos);
    status = sscanf(ch2, "%d", &num_atoms);
    //printf("DEBUG: parsed num_atoms = %d from string %s (status = %d)\n", num_atoms, ch2, status);
  } else {
    printf("ERROR: could not read string\n");
    exit(EXIT_FAILURE);
  }
  //printf("(status = %d) DEBUG: read that there will be %d atoms\n", status, num_atoms);
  search_for_char(input, '\n'); //now we are on a line containing an atom
  for(int i=0; i<num_atoms; i++) {
    XYZ a;
    char e[length];
    char ch[length];
    if(fgets(ch, length, input)!=NULL) {
      string str = string(ch);
      //printf("DEBUG: read char array as \"%s\", parsed to string as \"%s\"\n", ch, str.c_str());
      int pos=0;
      char c = str[pos];
      while(c<=0) { //for some reason, xyz files from marvinsketch have three characters at the beginning with ascii codes <0 - very strange, but we need to remove these to succesfully parse any data
        pos++;
        c = str[pos];
      }
      char ch2[length];
      str.copy(ch2, str.size()-pos, pos);
      status = sscanf(ch2, "%s %lf %lf %lf", e, &a.x, &a.y, &a.z);
    } else {
      printf("ERROR: could not read expected atom coord string from %s - %d out of %d atom coords were read\n", filename, i, num_atoms);
      exit(EXIT_FAILURE);
    }

    mol->atoms_xyz.push_back(a);
    string name_string(e);
    mol->atoms_label.push_back(name_string);
    //parse out any digits from the name of the element
    int name_length = name_string.length();
    int digit_pos = -1;
    for(int i=0; i<name_length && digit_pos==-1; i++) {
      if(isdigit(name_string[i])) digit_pos = i;
    }
    if(digit_pos==0) { //first character was a digit - this is no good
      printf("ERROR: could not parse label from atom beginning with a digit in read_xyz: %s: %s\n", filename, name_string.c_str());
      exit(EXIT_FAILURE);
    } else {
      string element(e);
      if(digit_pos>0) {
        element = name_string.substr(0,digit_pos);
      }
      mol->atoms_type.push_back(element);
    }
  }
}


