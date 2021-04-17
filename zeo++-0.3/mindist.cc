#include <iostream>
#include <cfloat>
#include <voro++.hh>

#include "mindist.h"
#include "geometry.h"
//#include "network.h"

using namespace std;
using namespace voro;

/* This file contains the functions necessary to calculate the minimum periodic distance between two points.
 * Although the calculation is trivial when the unit cell vectors are orthogonal, unexpected complications arise
 * when the unit cell vectors are non-orthogonal. As a result, this file contains a class which performs
 * the calculation. Each ATOM_NETWORK ultimately stores an instance of this class. //which is then used in geomety.cc
 */


/* Compute in which unit cells the Voronoi cells that constitute the unit cell lie. Store the resulting information using the provided
 * vector references.*/
void compute_images(double in_bx, double in_bxy, double in_by, double in_bxz, double in_byz, double in_bz,  
		    vector<int> &i_vals, vector<int> &j_vals, vector<int> &k_vals){

  unitcell uc(in_bx, in_bxy, in_by, in_bxz, in_byz, in_bz);
  vector<int> ijk_vals; vector<double> vd;
  uc.images(ijk_vals, vd);
  unsigned int i, j;
  i_vals.clear(); j_vals.clear(); k_vals.clear();
  for(i = j = 0; i < vd.size(); i++, j+=3){
    i_vals.push_back(ijk_vals[j]);
    j_vals.push_back(ijk_vals[j+1]);
    k_vals.push_back(ijk_vals[j+2]);
  }
  
}

/* Convert coordinates relative to the unit cell vectors into xyz coordinates and
*  store the result using the provided references. */
void abc_to_xyz(double a, double b, double c, double &x, double &y, double &z){
    //cout << "bx: " << bx << "     bxy: " << bxy << "    bxz: " << bxz << endl;
    //cout << "by: " << by << "     byz: " << byz << endl;
    //cout << "bz: " << bz << endl;
  x = a*bx+b*bxy+c*bxz;
  y = b*by+c*byz;
  z = c*bz;
}


void MIN_PER_DISTANCE::abc_to_xyz(double a, double b, double c, double &x, double &y, double &z){
    //cout << "bx: " << bx << "     bxy: " << bxy << "    bxz: " << bxz << endl;
    //cout << "by: " << by << "     byz: " << byz << endl;
    //cout << "bz: " << bz << endl;
  x = a*bx+b*bxy+c*bxz;
  y = b*by+c*byz;
  z = c*bz;
}
   
MIN_PER_DISTANCE::MIN_PER_DISTANCE(){}

MIN_PER_DISTANCE::MIN_PER_DISTANCE(double va_x, double vb_x, double vb_y, double vc_x, double vc_y, double vc_z) {
  bx  = va_x;  bxy = vb_x;  by = vb_y; bxz = vc_x; byz = vc_y;  bz  = vc_z;
  ivals = vector<int>(); jvals = vector<int>(); kvals = vector<int>();
  compute_images(bx, bxy, by, bxz, byz, bz, ivals, jvals, kvals);
}
  
/* Returns the minimum periodic distance between the two points whose coordinates are relative to the
 *  unit cell vectors. Also stores the minimum displacment vector between the two points using the provided references.*/
double MIN_PER_DISTANCE::minimum_periodic_distance(double a1, double b1, double c1, double a2, double b2, double c2, double &minDa, double &minDb, double &minDc){
  double va_a = a1, va_b = b1, va_c = c1;
  double vb_a = a2, vb_b = b2, vb_c = c2;
  double vc_a = 0.5, vc_b = 0.5, vc_c = 0.5;
  
  double vd_a = vb_a - va_a + vc_a;
  double vd_b = vb_b - va_b + vc_b;
  double vd_c = vb_c - va_c + vc_c;
  double ve_a = trans_to_origuc(vd_a);
  double ve_b = trans_to_origuc(vd_b);
  double ve_c = trans_to_origuc(vd_c);
  
  double minDist = DBL_MAX;
  
  for(unsigned int index = 0; index < ivals.size(); index++){
    double vi_a = ivals[index], vi_b = jvals[index], vi_c = kvals[index];
    double da = vc_a - (ve_a+vi_a);
    double db = vc_b - (ve_b+vi_b);
    double dc = vc_c - (ve_c+vi_c);
    double dx, dy, dz;
    abc_to_xyz(da, db, dc, dx, dy, dz);
    //cout << "a: " << da << "   b: " << db << "    c: " << dc << endl;
    //cout << "x: " << dx << "   y: " << dy << "    z: " << dz << endl;
    double dist = sqrt(dx*dx + dy*dy + dz*dz);
    if(dist < minDist){
      minDist = dist;
      minDa = -da; minDb = -db; minDc = -dc;
    }
  }
  //cout << "minDist here: " << minDist << endl;
  return minDist;
}

/* Rich edit: returns the closest periodic image of point 2 to point 1, whose coordinates are relative to the
 *  unit cell vectors. Also stores the minimum displacment vector between the two points using the provided references.*/
void MIN_PER_DISTANCE::closest_periodic_image(double a1, double b1, double c1, double a2, double b2, double c2, 
					      double &minDa, double &minDb, double &minDc, double &return_a, double &return_b, double &return_c){
  double va_a = a1, va_b = b1, va_c = c1;
  double vb_a = a2, vb_b = b2, vb_c = c2;
  double vc_a = 0.5, vc_b = 0.5, vc_c = 0.5;
  
  double vd_a = vb_a - va_a + vc_a;
  double vd_b = vb_b - va_b + vc_b;
  double vd_c = vb_c - va_c + vc_c;
  double ve_a = trans_to_origuc(vd_a);
  double ve_b = trans_to_origuc(vd_b);
  double ve_c = trans_to_origuc(vd_c);
  
  double minDist = DBL_MAX;
  
  for(unsigned int index = 0; index < ivals.size(); index++){
    double vi_a = ivals[index], vi_b = jvals[index], vi_c = kvals[index];
    double da = vc_a - (ve_a+vi_a);
    double db = vc_b - (ve_b+vi_b);
    double dc = vc_c - (ve_c+vi_c);
    double dx, dy, dz;
    abc_to_xyz(da, db, dc, dx, dy, dz);
    double dist = sqrt(dx*dx + dy*dy + dz*dz);
    if(dist < minDist){
      minDist = dist;
      minDa = -da; minDb = -db; minDc = -dc;
/*      return_a = ve_a+vi_a;
      return_b = ve_b+vi_b;
      return_c = ve_c+vi_c;*/ //incorrect - returns closest point with respect to the centre, c, not with respect to the stationary point, a
      return_a = a1+minDa;
      return_b = b1+minDb;
      return_c = c1+minDc;
    }
  }
}

/* Returns the minimum periodic distance between the two points whose coordinates are relative to the
 *  unit cell vectors. */
double MIN_PER_DISTANCE::minimum_periodic_distance(double a1, double b1, double c1, double a2, double b2, double c2){
  double minDa, minDb, minDc;
  //cout << "a1: " << a1 << "     b1: " << b1 << "    c1: " << c1 << endl;
  //cout << "a2: " << a2 << "     b2: " << b2 << "    c2: " << c2 << endl;
  //cout << "mindA: " << minDa << "     mindB: " << minDb << "    mindC: " << minDc << endl;
  //cout << minimum_periodic_distance(a1, b1, c1, a2, b2, c2, minDa, minDb, minDc) << endl;
  //cout << "mindA: " << minDa << "     mindB: " << minDb << "    mindC: " << minDc << endl;
  return minimum_periodic_distance(a1, b1, c1, a2, b2, c2, minDa, minDb, minDc);
}

/* Rich edit: returns the closest periodic image of point 2 to point 1, whose coordinates are relative to the
 *  unit cell vectors. */
void MIN_PER_DISTANCE::closest_periodic_image(double a1, double b1, double c1, double a2, double b2, double c2, double &return_a, double &return_b, double &return_c){
  double minDa, minDb, minDc;
  closest_periodic_image(a1, b1, c1, a2, b2, c2, minDa, minDb, minDc, return_a, return_b, return_c);
}

/* Print the periodic images checked by the distance calculator when determining the minimum periodic distance.*/
void MIN_PER_DISTANCE::print_images(){
  cout << "Printing images: " << "\n";
  for(unsigned int i = 0; i < ivals.size(); i++){
    cout << ivals[i] << " " << jvals[i] << " " << kvals[i] << "\n";
  }
  cout << "Images printed" << "\n" << "\n";
}
