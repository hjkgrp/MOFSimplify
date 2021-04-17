#ifndef MINDIST_H
#define MINDIST_H

#include <vector>

/* This file contains the functions necessary to calculate the minimum periodic distance between two points.
 * Although the calculation is trivial when the unit cell vectors are orthogonal, unexpected complications arise
 * when the unit cell vectors are non-orthogonal. As a result, this file contains a class which performs
 * the calculation. Each ATOM_NETWORK ultimately stores an instance of this class which is then used in geomety.cc
 */

/* Begin: Deprecated. */
static double bx,bxy,by,bxz,byz,bz;
/* End: Deprecated. */

/* Compute in which unit cells the Voronoi cells that constitute the unit cell lie. Store the resulting information using the provided
 * vector references.*/
void compute_images(double in_bx, double in_bxy, double in_by, double in_bxz, double in_byz, double in_bz,  
		    std::vector<int> &i_vals, std::vector<int> &j_vals, std::vector<int> &k_vals);

/* Begin: Deprecated. The function replicated in MIN_PER_DISTANCE class.*/
/* Convert coordinates relative to the unit cell vectors into xyz coordinates and
*  store the result using the provided references. */
void abc_to_xyz(double a, double b, double c, double &x, double &y, double &z);
/* End: Deprecated. */

/* Class used to calculate the minimum periodic distance between two points. Vital for accurate calculations
 * when using non-orthogonal unit cells. */
class MIN_PER_DISTANCE {
public:
    std::vector<int> ivals, jvals, kvals;
    double bx,bxy,by,bxz,byz,bz;
   MIN_PER_DISTANCE();
   MIN_PER_DISTANCE(double va_x, double vb_x, double vb_y, double vc_x, double vc_y, double vc_z);
  
    /* Convert coordinates relative to the unit cell vectors into xyz coordinates and
    *  store the result using the provided references. */
    void abc_to_xyz(double a, double b, double c, double &x, double &y, double &z);
  /* Returns the minimum periodic distance between the two points whose coordinates are relative to the
  *  unit cell vectors. Also stores the minimum displacment vector between the two points using the provided references.*/
   double minimum_periodic_distance(double a1, double b1, double c1, double a2, double b2, double c2, double &minDa, double &minDb, double &minDc);

  /* Rich edit: returns the closest periodic image of point 2 to point 1, whose coordinates are relative to the
  *  unit cell vectors. Also stores the minimum displacment vector between the two points using the provided references.*/
  void closest_periodic_image(double a1, double b1, double c1, double a2, double b2, double c2, 
			      double &minDa, double &minDb, double &minDc, double &return_a, double &return_b, double &return_c);

  /* Returns the minimum periodic distance between the two points whose coordinates are relative to the
  *  unit cell vectors. */
  double minimum_periodic_distance(double a1, double b1, double c1, double a2, double b2, double c2);

  /* Rich edit: returns the closest periodic image of point 2 to point 1, whose coordinates are relative to the
  *  unit cell vectors. */
  void closest_periodic_image(double a1, double b1, double c1, double a2, double b2, double c2, double &return_a, double &return_b, double &return_c);

  /* Print the periodic images checked by the distance calculator when determining the minimum periodic distance.*/
  void print_images();
};

#endif
