/** This file contains all of the necessary geometric object classes and
 *  procedures that are used throughout the package.
 */

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <iostream>
#include <vector>
#include <utility>

/** Returns the determinant of the provided 3x3 matrix of doubles. */
double calcDeterminant(double matrix[3][3]);

/** Returns the determinant of the provided 3x3 matrix of integers. */
int calcDeterminant(int matrix[3][3]);

/** Store the result of inverting the provided matrix into the new
    matrix. */
void invertMatrix(double matrix [3][3], double newMatrix[3][3]);

/* tries to invert a matrix, returning whether or not it succeeded */
bool tryInvertMatrix(double matrix [3][3], double newMatrix[3][3]);

/** Simple class used to store a set of three double precision numbers and perform 
 *  some trivial arithmetic operations. */
class Point{
public:
  double vals[3];

  Point (double myX = 0, double myY = 0, double myZ = 0);
  const Point add(Point other) const ;
  const bool equals(Point other) const;
  const Point subtract(Point other) const;
  const double magnitude() const;
  const Point scale(const double factor)const ;
  const double dot_product(Point other) const;
  double& operator[](const int index);
  const Point operator-(Point other) const;
  const Point operator+(Point other) const;
  const double operator*(Point other) const;
  const Point unit() const;
  const Point cross(Point other) const ;
  void print(std::ostream &out) const;
};

/**overrides the ostream opperator so that you can easily print a Point object**/
std::ostream &operator<<(std::ostream &out, Point &obj);

bool pointIsLess(Point p1, Point p2);


// Class TRIPLET and XYZ are moved from networkstorage */
/* Simple class used to represent a set of three integers. */
class TRIPLET {
public:
  int vals[3];
  int x, y, z;
  
  /* Construct a TRIPLET with the three provided components.*/
  TRIPLET(int myX, int myY, int myZ);

  /* Access one of the TRIPLETs three values.*/
  int& operator[](int index);
  
  /** Add each component of the triplet to that of the provided TRIPLET
   *  and return the result. */
  const TRIPLET operator+(const TRIPLET& other); 
};

/** Data structure used to store a set of three points */
class XYZ{
public:
  double x, y ,z;

  XYZ(double myX = 0.0, double myY = 0.0, double myZ = 0.0);
  
  void print(std::ostream& out = std::cout) const;

  const double magnitude() const;
  /** Return a instance whose three components are equal to this XYZ's components
   *  multiplied by the provided factor. */
  void set(double newX=0.0, double newY=0.0, double newZ=0.0);
  const XYZ scale(const double factor) const;
  void scale(const double factor, XYZ* result);
  const double dot_product(const XYZ& other) const;
  double& operator[](const int index);
  const XYZ operator-(const XYZ& other) const;
  const XYZ operator+(const XYZ& other) const;
  const double operator*(const XYZ& other) const;
  const XYZ unit() const;
  const XYZ cross(const XYZ& other) const; 
  double angle_between(const XYZ& other) const;
  double euclid_dist(const XYZ& other) const;
};

/* the origin */
const XYZ origin(0,0,0);

/* Translate the coordinate by unit cell increments so that it lies within the 0 to 1 range.*/
double trans_to_origuc(double x);
XYZ trans_to_origuc(XYZ x);

//functions inherited from framework builder
XYZ get_vector_from_to(XYZ start, XYZ end);
XYZ midpoint(XYZ a, XYZ b);
XYZ project_onto_line(XYZ initial, XYZ line_start, XYZ line_end);
XYZ project_onto_plane(XYZ initial, XYZ line_start, XYZ line_end, XYZ third);
XYZ RotatePointAboutLine(XYZ p, double theta, XYZ p1, XYZ p2);

/* Calculate the volume of a sphere of the provided radius.*/
double calcSphereVolume(double radius);

/** Calculates the Euclidean distance between (x1,y1,z1) and (x2,y2,z2). */
double calcEuclideanDistance(double x1, double y1, double z1, double x2, double y2, double z2);

/** Calculates the Euclidean distance between the provided Points. */
double calcEuclideanDistance(Point p1, Point p2);

/** Calculates the number of intersections between the line that travels through Points p1 and p2
 *  and the sphere with the provided characterstics. Returns a pair made up of the number of intersections
 *  and a vector containing the intersection points. Refer to http://local.wasp.uwa.edu.au/~pbourke/geometry/sphereline
 *  for more information.
 */
std::pair<int, std::vector<Point> > findLineIntersections(Point p1, Point p2, Point circCenter, double rad);

std::pair<int, std::vector<Point> > findLineSegIntersections(Point p1, Point p2, Point circCenter, double rad);

/* Returns the shortest distance along a sphere of the provided circle radius 
 * between the two provided points. Refer to http://en.wikipedia.org/wiki/Great-circle_distance. */ 
double calcSphereDistance(Point p1, Point p2, double circRad);

/* Returns the projection of the provided Point onto the circle with the specified center and radius. */
Point projectPointOnSphere(Point p, double circRad, Point circCenter);

/* Returns a pair of doubles representing the longitude and latitude of the given (x,y,z) point.
 * Refer to http://en.wikipedia.org/wiki/Geodetic_system#geodetic_to.2Ffrom_ECEF_coordinates */
std::pair<double,double> findLongAndLat(Point pt);

/** Returns the shortest distance from point to plane **/
double distToPlane(Point pnt,Point p,Point normal);

#endif
