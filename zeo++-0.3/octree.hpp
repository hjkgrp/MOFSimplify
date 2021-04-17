#ifndef OCTREE_H
#define OCTREE_H

#include "network.h"

#define MAXOBJECTSPERNODE 8
#define MAXOCTREEDEPTH 10 
#define MINBOXSIZE 0.1

typedef Point Vector

class Octree{
private:
  Vector m_ABC[3];
  Vector m_Normals_ABC[6];
  Node* pm_Head;
public:
  //Constructors
  Octree();
  //Destructors
  ~Octree();
  void Print();
};

class Node{
private:
  Node* p_Parent;
  Node** pm_Children;
  vector<Point>* pm_Objects;
  Node** pm_Neighbors;
  Point m_ABC;
  double m_Length;
  unsigned int m_Depth;
  //Private Methods
  void AddChildren();
  void OctreeConstruct(vector<Point>& objects);
  bool Threshold(vector<Point>& objects);
  bool Subset(Sphere object);
  bool Subset(Point object);
public:
  //Constructors
  Node();
  //Destructors
  ~Node();
  //Public Methods
  void Print();
  Node* GetOctant(Point abc); //Travels from Top to Bottom to find node
  Node* GetNextNode(unsigned int planeHit, Point abc); //Returns pointer to next node
};

#endif OCTREE_H
