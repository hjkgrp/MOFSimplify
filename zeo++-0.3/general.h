
/* Created by Maciej Haranczyk on 4/25/2011 */
/* this file stores general datatypes used
 * in various sections of the code */

#ifndef GENERAL_H
#define GENERAL_H


//strucutre that stores edges connecting segments
//typedef struct {
/*
class SEGCONN {
public:
int from, to;  // nodes involved in connection
int from_seg,to_seg; // segments that are connected
double max_radius;
int merged; // connection to be merged
}
*/
//} SEGCONN;

typedef struct {
     double x,y,z,r;
} SPHERE;


class NODESPHERE {
  public:
  double a, b, c, r;
  NODESPHERE(){};
 };

#endif
