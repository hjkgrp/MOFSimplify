
/* Functions to handle symmetry operatons and equivalent atom detection

   Spacegroups taken from http://it.iucr.org/A/

   added by M. Haranczyk, 11/11/11

*/

#include <cstdlib>
#include <cstdio>

#include "networkstorage.h"
#include "geometry.h"
#include "symmetry.h"
#include "zeo_consts.h"

/* Checks whether a given XYZ position is unique to any existing VERTEX in a ATOM_NETWORK */
bool IsUniqueVertex(XYZ *p, ATOM_NETWORK &cell) {
  double dist;
  for(int v=0; v<cell.vertices.size(); v++) {
//    dist=calcDistanceABC(p->x,p->y,p->z,cell.vertices.at(v).abc.x,cell.vertices.at(v).abc.y,cell.vertices.at(v).abc.z,&cell);
    dist=cell.calcDistanceABC(p->x,p->y,p->z,cell.vertices.at(v).abc.x,cell.vertices.at(v).abc.y,cell.vertices.at(v).abc.z);
    if(dist<DISTANCE_TOLERANCE) return false;
  }
  return true;  
}

/* Checks if a given point p is equivalent by symmetry with any of other points */
bool IsEquivalent(XYZ *p,vector <XYZ> *EqPos,ATOM_NETWORK &atomNet)
{
//cout << "EqPos vector side " << EqPos->size() << endl;
//cout << endl << endl;
  double dist;
  for(unsigned int pos=0; pos<EqPos->size(); pos++)
    {
    dist=atomNet.calcDistanceABC(p->x,p->y,p->z,EqPos->at(pos).x,EqPos->at(pos).y,EqPos->at(pos).z);
//cout << dist << " ";
    if(dist<TOLERANCE) return true;
    };
return false;
}

/* Identifies VERTICES that are equivalent by symmetry */
/*
vector < vector<int> >  IdentifyEquivalentVertices(ATOM_NETWORK &cell, int spacegroup) {
  vector <int> VertexGroup; VertexGroup.resize(cell.vertices.size(),-1); // for each atom stores information if it waswas assigned to a group
  vector < vector<int> > Groups;
  vector < vector<XYZ> > EqPosGroup;

  // First atom sets a new group of atoms equivalent by position

  VertexGroup[0]=0;
  vector <int> v(1,0);
  Groups.push_back(v);
  XYZ at(cell.vertices[0].abc.x, cell.vertices[0].abc.y, cell.vertices[0].abc.z);
  EqPosGroup.push_back(GetEquivalentPositions(spacegroup,&at));

  // Second and other atoms
  for(int i=1; i<cell.vertices.size(); i++) {
    double dist;

    for(unsigned int gr=0; gr<EqPosGroup.size(); gr++) {
      for(unsigned int pos=0; pos<EqPosGroup[gr].size(); pos++) {
//         dist=calcDistanceABC(cell.vertices[i].abc.x,cell.vertices[i].abc.y,cell.vertices[i].abc.z,EqPosGroup[gr].at(pos).x,EqPosGroup[gr].at(pos).y,EqPosGroup[gr].at(pos).z,&cell);
         dist=cell.calcDistanceABC(cell.vertices[i].abc.x,cell.vertices[i].abc.y,cell.vertices[i].abc.z,EqPosGroup[gr].at(pos).x,EqPosGroup[gr].at(pos).y,EqPosGroup[gr].at(pos).z);
         if(dist<TOLERANCE) {
            if(cell.vertices[i].name!=cell.vertices[Groups[gr].at(0)].name)
              {
              cerr << "Vertices of different names occupy positions equivalent by symmetry.\n";
              exit(EXIT_FAILURE); 
              };
            VertexGroup[i]=gr;
            Groups[gr].push_back(i);
            break;
            };
         };
         if(dist<TOLERANCE) break; // leave the outer loop if equivalent position was found
      }; // outer loop 

    if(dist>=TOLERANCE) // if dist larger than TOLERANCE, equivalent position has not been found and a new group of equivalent position is to be created
      {
      VertexGroup[i]=EqPosGroup.size();
      XYZ at(cell.vertices[i].abc.x, cell.vertices[i].abc.y, cell.vertices[i].abc.z);
      EqPosGroup.push_back(GetEquivalentPositions(spacegroup,&at));
      vector <int> v(1,i);
      Groups.push_back(v);
      };

    }; // loop other all atoms


 // Print Symmetry group summary 
 for(unsigned int i=0; i<Groups.size(); i++)
   {
   cout << "Group " << i << "  size= " << Groups[i].size() << "\n";
   };

  return Groups; 

}
*/

/* this function identifies atoms that are equivalent by symmetry */
vector < vector<int> >  IdentifyEquivalentAtoms(ATOM_NETWORK &atomNet, int spacegroup)
{
  vector <int> AtomGroup; AtomGroup.resize(atomNet.numAtoms,-1); // for each atom stores information if it waswas assigned to a group
  vector < vector<int> > Groups;
  vector < vector<XYZ> > EqPosGroup;

  // First atom sets a new group of atoms equivalent by position

  AtomGroup[0]=0;
  vector <int> v(1,0);
  Groups.push_back(v);
  XYZ at(atomNet.atoms[0].a_coord, atomNet.atoms[0].b_coord, atomNet.atoms[0].c_coord);
  EqPosGroup.push_back(GetEquivalentPositions(spacegroup,&at));

  // Second and other atoms
  for(int i=1; i<atomNet.numAtoms; i++)
    {
    double dist;

    for(unsigned int gr=0; gr<EqPosGroup.size(); gr++) 
      {
      for(unsigned int pos=0; pos<EqPosGroup[gr].size(); pos++)
         {
         dist=atomNet.calcDistanceXYZABC(atomNet.atoms[i].x,atomNet.atoms[i].y,atomNet.atoms[i].z,
                                    EqPosGroup[gr].at(pos).x,EqPosGroup[gr].at(pos).y,EqPosGroup[gr].at(pos).z);
         if(dist<TOLERANCE) 
            {
            if(atomNet.atoms[i].type!=atomNet.atoms[Groups[gr].at(0)].type)
              {
              cerr << "Atoms of different types occupy positions equivalent by symmetry.\n";
              exit(1); 
              };
            AtomGroup[i]=gr;
            Groups[gr].push_back(i);
            break;
            };
         };
         if(dist<TOLERANCE) break; // leave the outer loop if equivalent position was found
      }; // outer loop 

    if(dist>=TOLERANCE) // if dist larger than TOLERANCE, equivalent position has not been found and a new group of equivalent position is to be created
      {
      AtomGroup[i]=EqPosGroup.size();
      XYZ at(atomNet.atoms[i].a_coord, atomNet.atoms[i].b_coord, atomNet.atoms[i].c_coord);
      EqPosGroup.push_back(GetEquivalentPositions(spacegroup,&at));
      vector <int> v(1,i);
      Groups.push_back(v);
      };

    }; // loop other all atoms


 // Print Symmetry group summary 
 for(unsigned int i=0; i<Groups.size(); i++)
   {
   cout << "Group " << i << "  size= " << Groups[i].size() << "\n";
   };

  return Groups; 

}

/* Returns a vector of positions equivalent by symmetry to an input position
Spacegroups were taken from http://it.iucr.org/Ab/ch7o1v0001/    
For groups with multiple definitions, we took first axis choice, second origin choice
and rhombo-cell tather than hexagonal one. */
vector <XYZ> GetEquivalentPositions(int spacegroup, XYZ *pt) {
  vector <XYZ> EqPosVec; // vector with equivalent positions for a given point - this is returned
  XYZ p;
  double x,y,z;
  vector <XYZ> Sets;
  XYZ s;
  x=pt->x; y=pt->y; z=pt->z;
  switch(spacegroup)
  {
case 1 : // P 1
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 2: // P -1
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 3: // P 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 4: // P 2 1
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 5: // C 2
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 6: // P m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 7: // P c
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 8: // C m
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 9: // C c
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 10: // P 2 m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 11: // P 2 1 m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 12: // C 2 m
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 13: // P 2 c
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 14: // P 2 1 c
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 15: // C 2 c
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 16: // P 2 2 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 17: // P 2 2 2 1
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 18: // P 2 1 2 1 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 19: // P 2 1 2 1 2 1
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 20: // C 2 2 2 1
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 21: // C 2 2 2
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 22: // F 2 2 2
   s.set(0,0,0); Sets.push_back(s);
   s.set(0,1.0/2,1.0/2); Sets.push_back(s);
   s.set(1.0/2,0,1.0/2); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 23: // I 2 2 2
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 24: // I 2 1 2 1 2 1
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 25: // P m m 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 26: // P m c 2 1
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 27: // P c c 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 28: // P m a 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 29: // P c a 2 1
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 30: // P n c 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 31: // P m n 2 1
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 32: // P b a 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 33: // P n a 2 1
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 34: // P n n 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 35: // C m m 2
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 36: // C m c 2 1
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 37: // C c c 2
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 38: // A m m 2
   s.set(0,0,0); Sets.push_back(s);
   s.set(0,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 39: // A e m 2
   s.set(0,0,0); Sets.push_back(s);
   s.set(0,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 40: // A m a 2
   s.set(0,0,0); Sets.push_back(s);
   s.set(0,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 41: // A e a 2
   s.set(0,0,0); Sets.push_back(s);
   s.set(0,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 42: // F m m 2
   s.set(0,0,0); Sets.push_back(s);
   s.set(0,1.0/2,1.0/2); Sets.push_back(s);
   s.set(1.0/2,0,1.0/2); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 43: // F d d 2
   s.set(0,0,0); Sets.push_back(s);
   s.set(0,1.0/2,1.0/2); Sets.push_back(s);
   s.set(1.0/2,0,1.0/2); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/4,-y+1.0/4,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/4,y+1.0/4,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 44: // I m m 2
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 45: // I b a 2
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 46: // I m a 2
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 47: // P m m m 
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 48: // P n n n
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 49: // P c c m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 50: // P b a n
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 51: // P m m a
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 52: // P n n a
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 53: // P m n a
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 54: // P c c a 
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 55: // P b a m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 56: // P c c n
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 57: // P b c m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 58: // P n n m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 59: // P m m n
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 60: // P b c n
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 61: // P b c a
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 62: // P n m a 
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 63: // C m c m
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 64: // C m c e
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 65: // C m m m
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 66: // C c c m
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 67: // C m m e
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 68: // C c c e
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 69: // F m m m
   s.set(0,0,0); Sets.push_back(s);
   s.set(0,1.0/2,1.0/2); Sets.push_back(s);
   s.set(1.0/2,0,1.0/2); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 70: // F d d d
   s.set(0,0,0); Sets.push_back(s);
   s.set(0,1.0/2,1.0/2); Sets.push_back(s);
   s.set(1.0/2,0,1.0/2); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+3.0/4,-y+3.0/4,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+3.0/4,y,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+3.0/4,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/4,y+1.0/4,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/4,-y,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/4,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 71: // I m m m
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 72: // I b a m
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 73: // I b c a
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 74: // I m m a
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 75: // P 4
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 76: // P 41
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 77: // P 42
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 78: // P 43
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 79: // I 4
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 80: // I 41
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x+1.0/2,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 81: // P -4
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 82: // I -4
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 83: // P 4/m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 84: // P 42/m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 85: // P 4/n
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 86: // P 42/n
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 87: // I 4/m
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 88: // I 41/a
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+3.0/4,x+1.0/4,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+3.0/4,-x+3.0/4,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/4,-x+3.0/4,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/4,x+1.0/4,-z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 89: // P 4 2 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 90: // P 4 21 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 91: // P 41 2 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 92: // P 41 21 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x+1.0/2,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x+1.0/2,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,-z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 93: // P 42 2 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 94: // P 42 21 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 95: // P 43 2 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 96: // P 43 21 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x+1.0/2,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x+1.0/2,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 97: // I 4 2 2
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 98: // I 41 2 2
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x+1.0/2,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,-z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 99: // P 4 m m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 100: // P 4 b m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 101: // P 42 c m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 102: // P 42 n m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 103: // P 4 c c
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 104: // P 4 n c
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 105: // P 42 m c
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 106: // P 42 b c
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 107: // I 4 m m
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 108: // I 4 c m
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 109: // I 41 md
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x+1.0/2,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x+1.0/2,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 110: // I 41 c d
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x+1.0/2,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x+1.0/2,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 111: // P -4 2 m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 112: // P -4 2 c
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 113: // P -4 21 m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 114: // P -4 21 c
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 115: // P -4 m 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 116: // P -4 c 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 117: // P -4 b 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 118: // P -4 n 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 119: // I -4 m 2
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 120: // I -4 c 2
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 121: // I -4 2 m
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 122: // I -4 2 d
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 123: // P 4/m m m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 124: // P 4/m c c
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 125: // P 4/n b m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 126: // P 4/n n c 
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 127: // P 4/m b m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 128: // P 4/m n c
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 129: // P 4/n m m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 130: // P 4/n c c
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 131: // P 42/m m c
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 132: // P 42/m c m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 133: // P 42/n b c
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 134: // P 42/n n m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 135: // P 42/m b c
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 136: // P 42/m n m 
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 137: // P 42/n m c
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 138: // P 42/n c m 
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 139: // I 4/m m m
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 140: // I 4/m c m
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 141: // I 41/a m d
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/4,x+3.0/4,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/4,-x+1.0/4,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/4,x+3.0/4,-z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/4,-x+1.0/4,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+3.0/4,-x+1.0/4,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+3.0/4,x+3.0/4,-z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+3.0/4,-x+1.0/4,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+3.0/4,x+3.0/4,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 142: // I 41/a c d
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/4,x+3.0/4,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/4,-x+1.0/4,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/4,x+3.0/4,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/4,-x+1.0/4,-z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+3.0/4,-x+1.0/4,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+3.0/4,x+3.0/4,-z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+3.0/4,-x+1.0/4,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+3.0/4,x+3.0/4,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 143: // P 3
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 144: // P 31
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 145: // P 32
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 146: // R 3
   s.set(0,0,0); Sets.push_back(s);
   s.set(2.0/3,1.0/3,1.0/3); Sets.push_back(s);
   s.set(1.0/3,2.0/3,2.0/3); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 147: // P -3
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 148: // R -3
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 149: // P 3 1 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 150: // P 3 2 1
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 151: // P 31 1 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,-z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 152: // P 31 2 1
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,-z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,-z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 153: // P 32 1 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,-z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 154: // P 32 2 1
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,-z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,-z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 155: // R 3 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 156: // P 3 m 1
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 157: // P 3 1 m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 158: // P 3 c 1
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 159: // P 3 1 c
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 160: // R 3 m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 161: // R 3 c
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,y+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,z+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 162: // P -3 1 m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 163: // P -3 1 c
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 164: // P -3 m 1
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 165: // P -3 c 1
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 166: // R -3 m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 167: // R -3 c
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,-y+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-z+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,y+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,z+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 168: // P 6
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 169: // P 61
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,z+5.0/6); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,z+1.0/6); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 170: // P 65
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,z+1.0/6); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,z+5.0/6); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 171: // P 62
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 172: // P 64
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 173: // P 63
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 174: // P -6
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 175: // P 6/m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 176: // P 63/m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 177: // P 6 2 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 178: // P 61 2 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,z+5.0/6); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,z+1.0/6); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,-z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z+5.0/6); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,-z+1.0/6); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 179: // P 65 2 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,z+1.0/6); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,z+5.0/6); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,-z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z+1.0/6); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,-z+5.0/6); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 180: // P 62 2 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,-z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,-z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 181: // P 64 2 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,-z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z+1.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,-z+2.0/3); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 182: // P 63 2 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 183: // P 6 m m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 184: // P 6 c c
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 185: // P 63 c m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 186: // P 63 m c
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 187: // P -6 m 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 188: // P -6 c 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 189: // P -6 2 m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 190: // P -6 2 c
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 191: // P 6/m m m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 192: // P 6/m c c
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 193: // P 63/m c m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 194: // P 63/m m c
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 195: // P 2 3
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 196: // F 2 3
   s.set(0,0,0); Sets.push_back(s);
   s.set(0,1.0/2,1.0/2); Sets.push_back(s);
   s.set(1.0/2,0,1.0/2); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 197: // I 2 3
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 198: // P 21 3
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,-x+1.0/2,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,-x,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-z+1.0/2,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-z,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 199: // I 21 3
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,-x+1.0/2,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,-x,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-z+1.0/2,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-z,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 200: // P m -3
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 201: // P n -3
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,-x+1.0/2,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,x,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,z,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-z+1.0/2,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,x+1.0/2,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,-x,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-z,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,z+1.0/2,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 202: // F m -3
   s.set(0,0,0); Sets.push_back(s);
   s.set(0,1.0/2,1.0/2); Sets.push_back(s);
   s.set(1.0/2,0,1.0/2); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 203: // F d -3
   s.set(0,0,0); Sets.push_back(s);
   s.set(0,1.0/2,1.0/2); Sets.push_back(s);
   s.set(1.0/2,0,1.0/2); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+3.0/4,-y+3.0/4,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+3.0/4,y,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+3.0/4,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x+3.0/4,-y+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+3.0/4,-x+3.0/4,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+3.0/4,x,-y+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+3.0/4,z,-x+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z+3.0/4,-x+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+3.0/4,-z+3.0/4,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/4,y+1.0/4,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/4,-y,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/4,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x+1.0/4,y+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/4,x+1.0/4,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/4,-x,y+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/4,-z,x+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z+1.0/4,x+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/4,z+1.0/4,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 204: // I m -3
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 205: // P a -3
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,-x+1.0/2,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,-x,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-z+1.0/2,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-z,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,x+1.0/2,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,x,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,z+1.0/2,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,z,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 206: // I a -3
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,-x+1.0/2,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,-x,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-z+1.0/2,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-z,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,x+1.0/2,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,x,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,z+1.0/2,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,z,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 207: // P 4 3 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 208: // P 42 3 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,z+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,z+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-z+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-z+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,y+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,-y+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,y+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,-y+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 209: // F 4 3 2
   s.set(0,0,0); Sets.push_back(s);
   s.set(0,1.0/2,1.0/2); Sets.push_back(s);
   s.set(1.0/2,0,1.0/2); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 210: // F 41 3 2
   s.set(0,0,0); Sets.push_back(s);
   s.set(0,1.0/2,1.0/2); Sets.push_back(s);
   s.set(1.0/2,0,1.0/2); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,-x,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,x+1.0/2,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,z+1.0/2,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-z,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+3.0/4,x+1.0/4,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/4,-x+1.0/4,-z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/4,-x+3.0/4,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+3.0/4,x+3.0/4,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+3.0/4,z+1.0/4,-y+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+3.0/4,z+3.0/4,y+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/4,-z+1.0/4,-y+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/4,-z+3.0/4,y+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+3.0/4,y+1.0/4,-x+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/4,-y+3.0/4,x+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+3.0/4,y+3.0/4,x+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/4,-y+1.0/4,-x+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 211: // I 4 3 2
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 212: // P 43 3 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,-x+1.0/2,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,-x,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-z+1.0/2,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-z,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/4,x+3.0/4,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/4,-x+1.0/4,-z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+3.0/4,-x+3.0/4,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+3.0/4,x+1.0/4,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/4,z+3.0/4,-y+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+3.0/4,z+1.0/4,y+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/4,-z+1.0/4,-y+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+3.0/4,-z+3.0/4,y+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/4,y+3.0/4,-x+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+3.0/4,-y+3.0/4,x+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+3.0/4,y+1.0/4,x+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/4,-y+1.0/4,-x+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 213: // P 41 3 2
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,-x+1.0/2,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,-x,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-z+1.0/2,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-z,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+3.0/4,x+1.0/4,-z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+3.0/4,-x+3.0/4,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/4,-x+1.0/4,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/4,x+3.0/4,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+3.0/4,z+1.0/4,-y+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/4,z+3.0/4,y+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+3.0/4,-z+3.0/4,-y+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/4,-z+1.0/4,y+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+3.0/4,y+1.0/4,-x+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/4,-y+1.0/4,x+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/4,y+3.0/4,x+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+3.0/4,-y+3.0/4,-x+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 214: // I 41 3 2
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,-x+1.0/2,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,-x,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-z+1.0/2,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-z,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+3.0/4,x+1.0/4,-z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+3.0/4,-x+3.0/4,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/4,-x+1.0/4,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/4,x+3.0/4,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+3.0/4,z+1.0/4,-y+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/4,z+3.0/4,y+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+3.0/4,-z+3.0/4,-y+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/4,-z+1.0/4,y+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+3.0/4,y+1.0/4,-x+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/4,-y+1.0/4,x+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/4,y+3.0/4,x+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+3.0/4,-y+3.0/4,-x+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 215: // P -4 3 m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 216: // F -4 3 m
   s.set(0,0,0); Sets.push_back(s);
   s.set(0,1.0/2,1.0/2); Sets.push_back(s);
   s.set(1.0/2,0,1.0/2); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 217: // I -4 3 m
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 218: // P -4 3 n
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,z+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,z+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-z+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-z+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,y+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,-y+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,y+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,-y+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 219: // F -4 3 c
   s.set(0,0,0); Sets.push_back(s);
   s.set(0,1.0/2,1.0/2); Sets.push_back(s);
   s.set(1.0/2,0,1.0/2); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,z+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,z+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-z+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-z+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,y+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,-y+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,y+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,-y+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 220: // I -4 3 d
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,-x+1.0/2,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,-x,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-z+1.0/2,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-z,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/4,x+1.0/4,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/4,-x+3.0/4,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+3.0/4,-x+1.0/4,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+3.0/4,x+3.0/4,-z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/4,z+1.0/4,y+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+3.0/4,z+3.0/4,-y+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/4,-z+3.0/4,y+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+3.0/4,-z+1.0/4,-y+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/4,y+1.0/4,x+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+3.0/4,-y+1.0/4,-x+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+3.0/4,y+3.0/4,-x+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/4,-y+3.0/4,x+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 221: // P m -3 m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 222: // P n -3 n
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,-x+1.0/2,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,x,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,z,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-z+1.0/2,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,z,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-z+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-z+1.0/2,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,y,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-y+1.0/2,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,-y+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,x+1.0/2,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,-x,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-z,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,z+1.0/2,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-z,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,z+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,z+1.0/2,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-y,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,y+1.0/2,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,-y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,y+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 223: // P m -3 n
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,z+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,z+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-z+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-z+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,y+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,-y+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,y+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,-y+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-z+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-z+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,z+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,z+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,-y+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,y+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,-y+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,y+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 224: // P n -3 m
   s.set(0,0,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,-x+1.0/2,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,x,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,z,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-z+1.0/2,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,z+1.0/2,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,z+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-z,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,y+1.0/2,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,-y,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,y+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,x+1.0/2,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,-x,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-z,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,z+1.0/2,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-z+1.0/2,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-z+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,z,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,-y+1.0/2,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,y,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-y+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 225: // F m -3 m
   s.set(0,0,0); Sets.push_back(s);
   s.set(0,1.0/2,1.0/2); Sets.push_back(s);
   s.set(1.0/2,0,1.0/2); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 226: // F m -3 c
   s.set(0,0,0); Sets.push_back(s);
   s.set(0,1.0/2,1.0/2); Sets.push_back(s);
   s.set(1.0/2,0,1.0/2); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,z+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,z+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-z+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-z+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,y+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,-y+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,y+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,-y+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-z+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-z+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,z+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,z+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,-y+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,y+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,-y+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,y+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 227: // F d -3 m
   s.set(0,0,0); Sets.push_back(s);
   s.set(0,1.0/2,1.0/2); Sets.push_back(s);
   s.set(1.0/2,0,1.0/2); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+3.0/4,-y+1.0/4,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/4,y+1.0/2,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+3.0/4,-z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,-x+3.0/4,-y+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+3.0/4,-x+1.0/4,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/4,x+1.0/2,-y+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/4,z+1.0/2,-x+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-z+3.0/4,-x+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+3.0/4,-z+1.0/4,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+3.0/4,x+1.0/4,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/4,-x+1.0/2,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,x+3.0/4,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+3.0/4,z+1.0/4,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,z+3.0/4,y+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/4,-z+1.0/2,y+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+3.0/4,y+1.0/4,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/4,-y+1.0/2,x+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,y+3.0/4,x+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/4,y+3.0/4,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+3.0/4,-y+1.0/2,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/4,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,x+1.0/4,y+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/4,x+3.0/4,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+3.0/4,-x+1.0/2,y+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+3.0/4,-z+1.0/2,x+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,z+1.0/4,x+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/4,z+3.0/4,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/4,-x+3.0/4,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+3.0/4,x+1.0/2,-z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-x+1.0/4,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/4,-z+3.0/4,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-z+1.0/4,-y+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+3.0/4,z+1.0/2,-y+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/4,-y+3.0/4,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+3.0/4,y+1.0/2,-x+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,-y+1.0/4,-x+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 228: // F d -3 c
   s.set(0,0,0); Sets.push_back(s);
   s.set(0,1.0/2,1.0/2); Sets.push_back(s);
   s.set(1.0/2,0,1.0/2); Sets.push_back(s);
   s.set(1.0/2,1.0/2,0); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/4,-y+3.0/4,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+3.0/4,y+1.0/2,-z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/4,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,-x+1.0/4,-y+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/4,-x+3.0/4,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+3.0/4,x+1.0/2,-y+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+3.0/4,z+1.0/2,-x+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-z+1.0/4,-x+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/4,-z+3.0/4,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+3.0/4,x+1.0/4,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-x+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/4,-x,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x+3.0/4,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+3.0/4,z+1.0/4,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,z+3.0/4,y+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-z+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/4,-z,y+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+3.0/4,y+1.0/4,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/4,-y,x+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,y+3.0/4,x+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,-y+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+3.0/4,y+1.0/4,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/4,-y+1.0/2,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+3.0/4,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,x+3.0/4,y+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+3.0/4,x+1.0/4,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/4,-x+1.0/2,y+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/4,-z+1.0/2,x+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,z+3.0/4,x+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+3.0/4,z+1.0/4,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/4,-x+3.0/4,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,x+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+3.0/4,x,-z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+1.0/4,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/4,-z+3.0/4,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-z+1.0/4,-y+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,z+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+3.0/4,z,-y+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/4,-y+3.0/4,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+3.0/4,y,-x+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-y+1.0/4,-x+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,y+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 229: // I m -3 m
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,z,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,z,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-y,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,y,x); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

case 230: // I a -3 d
   s.set(0,0,0); Sets.push_back(s);
   s.set(1.0/2,1.0/2,1.0/2); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,-y,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,y+1.0/2,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,-y+1.0/2,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,x,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,-x+1.0/2,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,-x,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,x+1.0/2,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,z,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,z+1.0/2,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,-z+1.0/2,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,-z,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+3.0/4,x+1.0/4,-z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+3.0/4,-x+3.0/4,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/4,-x+1.0/4,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/4,x+3.0/4,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+3.0/4,z+1.0/4,-y+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/4,z+3.0/4,y+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+3.0/4,-z+3.0/4,-y+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/4,-z+1.0/4,y+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+3.0/4,y+1.0/4,-x+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/4,-y+1.0/4,x+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/4,y+3.0/4,x+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+3.0/4,-y+3.0/4,-x+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/2,y,-z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,-y+1.0/2,z+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/2,y+1.0/2,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z,-x,-y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/2,x+1.0/2,y); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/2,x,-y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z,-x+1.0/2,y+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-z,-x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-z+1.0/2,x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/2,z+1.0/2,x); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/2,z,-x+1.0/2); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+1.0/4,-x+3.0/4,z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+1.0/4,x+1.0/4,z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y+3.0/4,x+3.0/4,-z+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y+3.0/4,-x+1.0/4,-z+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+1.0/4,-z+3.0/4,y+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+3.0/4,-z+1.0/4,-y+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x+1.0/4,z+1.0/4,y+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+3.0/4,z+3.0/4,-y+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+1.0/4,-y+3.0/4,x+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-z+3.0/4,y+3.0/4,-x+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+3.0/4,-y+1.0/4,-x+3.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(z+1.0/4,y+1.0/4,x+1.0/4); p=p+(Sets[i]); EqPosVec.push_back(p);

     };
    break; // end of 230

//MANUAL ENTRY

case 1000: // R -3 m :H
   s.set(0,0,0); Sets.push_back(s);
   s.set(2.0/3,1.0/3,1.0/3); Sets.push_back(s);
   s.set(1.0/3,2.0/3,2.0/3); Sets.push_back(s);
   for(unsigned int i=0; i<Sets.size(); i++)
     {

p.set(x,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-x+y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x,-y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(y,-x+y,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x-y,x,-z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-y,-x,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(-x+y,y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
p.set(x,x-y,z); p=p+(Sets[i]); EqPosVec.push_back(p);
     };
   break;

  default: // default switch to handle unknown groups 
    cerr << "Unknown space group: << " << spacegroup << ". Check symmetry.c for definitions\n";
    exit(EXIT_FAILURE);
  }; // end of switch

return EqPosVec;
}

/* Return the symmetry group ID number based on matching a string */
int get_sym_ID(string s) {
  printf("symmetry group = %s\n", s.c_str());
/*

  //this block edits symmetry group strings - the removel of ":H" is wrong because it does not represent the correct symmetry, but this kind of string operation might be needed later

  //remove ":H" from the end
  int length = s.length();
  printf("DEBUG: before truncating string = %s\n", s.c_str());
  for(int i=length-2; i>=0; i--) {
    if(s[i]==':' && s[i+1]=='H') {
      s.erase(i,2);
      i--; //additional i-- moves cursor again, to prevent checking the end string character
    }
  }
  printf("after truncating, symmetry group = %s\n", s.c_str());
*/
  //scan list of space group strings
  if(s=="P1") return 1;
  if(s=="P-1") return 2;
  if(s=="P2") return 3;
  if(s=="P21") return 4;
  if(s=="C2") return 5;
  if(s=="Pm") return 6;
  if(s=="Pc") return 7;
  if(s=="Cm") return 8;
  if(s=="Cc") return 9;
  if(s=="P2m") return 10;
  if(s=="P21m") return 11;
  if(s=="C2m") return 12;
  if(s=="P2c") return 13;
  if(s=="P21c") return 14;
  if(s=="C2c") return 15;
  if(s=="P222") return 16;
  if(s=="P2221") return 17;
  if(s=="P21212") return 18;
  if(s=="P212121") return 19;
  if(s=="C2221") return 20;
  if(s=="C222") return 21;
  if(s=="F222") return 22;
  if(s=="I222") return 23;
  if(s=="I212121") return 24;
  if(s=="Pmm2") return 25;
  if(s=="Pmc21") return 26;
  if(s=="Pcc2") return 27;
  if(s=="Pma2") return 28;
  if(s=="Pca21") return 29;
  if(s=="Pnc2") return 30;
  if(s=="Pmn21") return 31;
  if(s=="Pba2") return 32;
  if(s=="Pna21") return 33;
  if(s=="Pnn2") return 34;
  if(s=="Cmm2") return 35;
  if(s=="Cmc21") return 36;
  if(s=="Ccc2") return 37;
  if(s=="Amm2") return 38;
  if(s=="Aem2") return 39;
  if(s=="Ama2") return 40;
  if(s=="Aea2") return 41;
  if(s=="Fmm2") return 42;
  if(s=="Fdd2") return 43;
  if(s=="Imm2") return 44;
  if(s=="Iba2") return 45;
  if(s=="Ima2") return 46;
  if(s=="Pmmm") return 47;
  if(s=="Pnnn") return 48;
  if(s=="Pccm") return 49;
  if(s=="Pban") return 50;
  if(s=="Pmma") return 51;
  if(s=="Pnna") return 52;
  if(s=="Pmna") return 53;
  if(s=="Pcca") return 54;
  if(s=="Pbam") return 55;
  if(s=="Pccn") return 56;
  if(s=="Pbcm") return 57;
  if(s=="Pnnm") return 58;
  if(s=="Pmmn") return 59;
  if(s=="Pbcn") return 60;
  if(s=="Pbca") return 61;
  if(s=="Pnma") return 62;
  if(s=="Cmcm") return 63;
  if(s=="Cmce") return 64;
  if(s=="Cmmm") return 65;
  if(s=="Cccm") return 66;
  if(s=="Cmme") return 67;
  if(s=="Ccce") return 68;
  if(s=="Fmmm") return 69;
  if(s=="Fddd") return 70;
  if(s=="Immm") return 71;
  if(s=="Ibam") return 72;
  if(s=="Ibca") return 73;
  if(s=="Imma") return 74;
  if(s=="P4") return 75;
  if(s=="P41") return 76;
  if(s=="P42") return 77;
  if(s=="P43") return 78;
  if(s=="I4") return 79;
  if(s=="I41") return 80;
  if(s=="P-4") return 81;
  if(s=="I-4") return 82;
  if(s=="P4/m") return 83;
  if(s=="P42/m") return 84;
  if(s=="P4/n") return 85;
  if(s=="P42/n") return 86;
  if(s=="I4/m") return 87;
  if(s=="I41/a") return 88;
  if(s=="P422") return 89;
  if(s=="P4212") return 90;
  if(s=="P4122") return 91;
  if(s=="P41212") return 92;
  if(s=="P4222") return 93;
  if(s=="P42212") return 94;
  if(s=="P4322") return 95;
  if(s=="P43212") return 96;
  if(s=="I422") return 97;
  if(s=="I4122") return 98;
  if(s=="P4mm") return 99;
  if(s=="P4bm") return 100;
  if(s=="P42cm") return 101;
  if(s=="P42nm") return 102;
  if(s=="P4cc") return 103;
  if(s=="P4nc") return 104;
  if(s=="P42mc") return 105;
  if(s=="P42bc") return 106;
  if(s=="I4mm") return 107;
  if(s=="I4cm") return 108;
  if(s=="I41md") return 109;
  if(s=="I41cd") return 110;
  if(s=="P-42m") return 111;
  if(s=="P-42c") return 112;
  if(s=="P-421m") return 113;
  if(s=="P-421c") return 114;
  if(s=="P-4m2") return 115;
  if(s=="P-4c2") return 116;
  if(s=="P-4b2") return 117;
  if(s=="P-4n2") return 118;
  if(s=="I-4m2") return 119;
  if(s=="I-4c2") return 120;
  if(s=="I-42m") return 121;
  if(s=="I-42d") return 122;
  if(s=="P4/mmm") return 123;
  if(s=="P4/mcc") return 124;
  if(s=="P4/nbm") return 125;
  if(s=="P4/nnc") return 126;
  if(s=="P4/mbm") return 127;
  if(s=="P4/mnc") return 128;
  if(s=="P4/nmm") return 129;
  if(s=="P4/ncc") return 130;
  if(s=="P42/mmc") return 131;
  if(s=="P42/mcm") return 132;
  if(s=="P42/nbc") return 133;
  if(s=="P42/nnm") return 134;
  if(s=="P42/mbc") return 135;
  if(s=="P42/mnm") return 136;
  if(s=="P42/nmc") return 137;
  if(s=="P42/ncm") return 138;
  if(s=="I4/mmm") return 139;
  if(s=="I4/mcm") return 140;
  if(s=="I41/amd") return 141;
  if(s=="I41/acd") return 142;
  if(s=="P3") return 143;
  if(s=="P31") return 144;
  if(s=="P32") return 145;
  if(s=="R3") return 146;
  if(s=="P-3") return 147;
  if(s=="R-3") return 148;
  if(s=="P312") return 149;
  if(s=="P321") return 150;
  if(s=="P3112") return 151;
  if(s=="P3121") return 152;
  if(s=="P3212") return 153;
  if(s=="P3221") return 154;
  if(s=="R32") return 155;
  if(s=="P3m1") return 156;
  if(s=="P31m") return 157;
  if(s=="P3c1") return 158;
  if(s=="P31c") return 159;
  if(s=="R3m") return 160;
  if(s=="R3c") return 161;
  if(s=="P-31m") return 162;
  if(s=="P-31c") return 163;
  if(s=="P-3m1") return 164;
  if(s=="P-3c1") return 165;
  if(s=="R-3m") return 166;
  if(s=="R-3c") return 167;
  if(s=="P6") return 168;
  if(s=="P61") return 169;
  if(s=="P65") return 170;
  if(s=="P62") return 171;
  if(s=="P64") return 172;
  if(s=="P63") return 173;
  if(s=="P-6") return 174;
  if(s=="P6/m") return 175;
  if(s=="P63/m") return 176;
  if(s=="P622") return 177;
  if(s=="P6122") return 178;
  if(s=="P6522") return 179;
  if(s=="P6222") return 180;
  if(s=="P6422") return 181;
  if(s=="P6322") return 182;
  if(s=="P6mm") return 183;
  if(s=="P6cc") return 184;
  if(s=="P63cm") return 185;
  if(s=="P63mc") return 186;
  if(s=="P-6m2") return 187;
  if(s=="P-6c2") return 188;
  if(s=="P-62m") return 189;
  if(s=="P-62c") return 190;
  if(s=="P6/mmm") return 191;
  if(s=="P6/mcc") return 192;
  if(s=="P63/mcm") return 193;
  if(s=="P63/mmc") return 194;
  if(s=="P23") return 195;
  if(s=="F23") return 196;
  if(s=="I23") return 197;
  if(s=="P213") return 198;
  if(s=="I213") return 199;
  if(s=="Pm-3") return 200;
  if(s=="Pn-3") return 201;
  if(s=="Fm-3") return 202;
  if(s=="Fd-3") return 203;
  if(s=="Im-3") return 204;
  if(s=="Pa-3") return 205;
  if(s=="Ia-3") return 206;
  if(s=="P432") return 207;
  if(s=="P4232") return 208;
  if(s=="F432") return 209;
  if(s=="F4132") return 210;
  if(s=="I432") return 211;
  if(s=="P4332") return 212;
  if(s=="P4132") return 213;
  if(s=="I4132") return 214;
  if(s=="P-43m") return 215;
  if(s=="F-43m") return 216;
  if(s=="I-43m") return 217;
  if(s=="P-43n") return 218;
  if(s=="F-43c") return 219;
  if(s=="I-43d") return 220;
  if(s=="Pm-3m") return 221;
  if(s=="Pn-3n") return 222;
  if(s=="Pm-3n") return 223;
  if(s=="Pn-3m") return 224;
  if(s=="Fm-3m") return 225;
  if(s=="Fm-3c") return 226;
  if(s=="Fd-3m") return 227;
  if(s=="Fd-3c") return 228;
  if(s=="Im-3m") return 229;
  if(s=="Ia-3d") return 230;
//MANUAL ENTRY
  if(s=="R-3m:H") return 1000;
  printf("WARNING: could not parse symmetry group string \"%s\" to find the corresponding ID number\n", s.c_str());
  return -1; //could not parse string!
}

