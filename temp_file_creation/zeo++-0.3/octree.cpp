#include "network.h"

Octree::Octree()
{
  m_ABC[0] = Vector(1,0,0);
  m_ABC[1] = Vector(0,1,0);
  m_ABC[2] = Vector(0,0,1);
    
  pm_Head = new Node();

  InitializeOctree();
}

Node::Node()
{
  p_Parent = NULL;
  pm_Children = NULL;
  pm_Objects = NULL;
  pm_Neighbors = NULL;
  m_Length = 0.0;
  m_Depth = 0;
}

Node::~Node()
{
  if (pm_Children != NULL)
    {
      for (unsigned int i=0;i<8; i++)
	{
	  delete *pm_Children[i];
	}
    }
  delete Node;
}

void Node::Print()
{
  if (pm_Children != NULL)
    {
      for (unsigned int i=0;i<8; i++)
	{
	  pm_Children[i]->Print();
	}
    }
  else
    {
      std::cout << "\nDepth " << m_Depth << std::endl;
      for (unsigned int i=0; i<pm_objects->size(); i++)
	{
	  std::cout << "Element " << i << ": " << pm_objects->at(i) << std::endl;
	}
    }
}

bool Node::Threshold(vector<Point>& objects)
{
  if ((object.size() < MAXOBJECTSPERNODE) || m_Depth > MAXOCTREEDEPTH || m_Length > MINBOXSIZE)
    {
      return true;
    }
  return false;
}

Node* Node::OctreeConstruct(vector<Point>& objects)
{
  if (Threshold(objects) == true)
    {
      pm_Children = NULL;
      m_Objects = new vector<Point>(objects.begin(),objects.end());
    }
  else //Threshold condition is not satisfied
    {
      InitializeChildren();

      for (unsigned int i=0; i<8; i++) //For each octant reinvoke octree
	{  
	  vector<Point> subset_Objects;
	  for (unsigned int j=0; j<objects.size(); j++)
	    {
	      if (Subset(objects[j]) == true)
		{
		  subset_Objects.push_back(objects[j]);
		}
	    }

	  pm_Children[i]->OctreeConstruct(subset_Objects);
	}
    }
}

Node* Node::InitializeChildren()
{
  pm_Children = new Node* [8];
  pm_Neighbors = new Node* [26];
  for (unsigned int i=0; i<8; i++)
    {
      pm_Children[i] = new Node();
      pm_Children[i]->p_Parent = this;
      pm_Children[i]->m_Depth = m_Depth + 1;
      pm_Children[i]->m_Length =  0.5 * m_Length;
    }

  //Now I will add all node specific information Eg. Location, and pointers which will significantly speed up the code
  double delta_h = 0.5 * m_Length;
  for (unsinged int i=0; i<8; i++)
    {
      switch(i)
	{
	case 0:
	  pm_Children[i]->m_ABC = m_ABC + Point(delta_h,-delta_h,delta_h);
	  pm_Children[i]->pm_Neighbors[0] = pm_Neighbors[0]; //U
	  pm_Children[i]->pm_Neighbors[1] = pm_Children[4]; //D
	  pm_Children[i]->pm_Neighbors[2] = pm_Neighbors[2]; //L
	  pm_Children[i]->pm_Neighbors[3] = pm_Children[2]; //R
	  pm_Children[i]->pm_Neighbors[4] = pm_Neighbors[4]; //F
	  pm_Children[i]->pm_Neighbors[5] = pm_Children[1]; //B
	  pm_Children[i]->pm_Neighbors[6] = pm_Neighbors[6]; //UL
	  pm_Children[i]->pm_Neighbors[7] = pm_Neighbors[0]; //UR
	  pm_Children[i]->pm_Neighbors[8] = pm_Neighbors[8]; //UF
	  pm_Children[i]->pm_Neighbors[9] = pm_Neighbors[0]; //UB
	  pm_Children[i]->pm_Neighbors[10] = pm_Neighbors[2]; //DL
	  pm_Children[i]->pm_Neighbors[11] = pm_Children[6]; //DR
	  pm_Children[i]->pm_Neighbors[12] = pm_Neighbors[4]; //DF
	  pm_Children[i]->pm_Neighbors[13] = pm_Children[5]; //DB
	  pm_Children[i]->pm_Neighbors[14] = pm_Neighbors[14]; //LF
	  pm_Children[i]->pm_Neighbors[15] = pm_Neighbors[2]; //LB
	  pm_Children[i]->pm_Neighbors[16] = pm_Neighbors[4]; //RF
	  pm_Children[i]->pm_Neighbors[17] = pm_Children[3]; //RB
	  pm_Children[i]->pm_Neighbors[18] = pm_Neighbors[18]; //ULF
	  pm_Children[i]->pm_Neighbors[19] = pm_Neighbors[6]; //ULB
	  pm_Children[i]->pm_Neighbors[20] = pm_Neighbors[8]; //URF
	  pm_Children[i]->pm_Neighbors[21] = pm_Neighbors[0]; //URB
	  pm_Children[i]->pm_Neighbors[22] = pm_Neighbors[14]; //DLF
	  pm_Children[i]->pm_Neighbors[23] = pm_Neighbors[2]; //DLB
	  pm_Children[i]->pm_Neighbors[24] = pm_Neighbors[4]; //DRF
	  pm_Children[i]->pm_Neighbors[25] = pm_Children[7]; //DRB
	  break;
	case 1:
	  pm_Children[i]->m_ABC = m_ABC + Point-(delta_h,-delta_h,delta_h);
	  pm_Children[i]->pm_Neighbors[0] = pm_Neighbors[0]; //U
	  pm_Children[i]->pm_Neighbors[1] = pm_Children[5]; //D
	  pm_Children[i]->pm_Neighbors[2] = pm_Neighbors[2]; //L
	  pm_Children[i]->pm_Neighbors[3] = pm_Children[3]; //R
	  pm_Children[i]->pm_Neighbors[4] = pm_Children[0]; //F
	  pm_Children[i]->pm_Neighbors[5] = pm_Neighbors[5]; //B
	  pm_Children[i]->pm_Neighbors[6] = pm_Neighbors[6]; //UL
	  pm_Children[i]->pm_Neighbors[7] = pm_Neighbors[0]; //UR
	  pm_Children[i]->pm_Neighbors[8] = pm_Neighbors[0]; //UF
	  pm_Children[i]->pm_Neighbors[9] = pm_Neighbors[9]; //UB
	  pm_Children[i]->pm_Neighbors[10] = pm_Neighbors[2]; //DL
	  pm_Children[i]->pm_Neighbors[11] = pm_Children[7]; //DR
	  pm_Children[i]->pm_Neighbors[12] = pm_Children[4]; //DF
	  pm_Children[i]->pm_Neighbors[13] = pm_Neighbors[5]; //DB
	  pm_Children[i]->pm_Neighbors[14] = pm_Neighbors[2]; //LF
	  pm_Children[i]->pm_Neighbors[15] = pm_Neighbors[15]; //LB
	  pm_Children[i]->pm_Neighbors[16] = pm_Children[2]; //RF
	  pm_Children[i]->pm_Neighbors[17] = pm_Neighbors[5]; //RB
	  pm_Children[i]->pm_Neighbors[18] = pm_Neighbors[6]; //ULF
	  pm_Children[i]->pm_Neighbors[19] = pm_Neighbors[19]; //ULB
	  pm_Children[i]->pm_Neighbors[20] = pm_Neighbors[0]; //URF
	  pm_Children[i]->pm_Neighbors[21] = pm_Neighbors[9]; //URB
	  pm_Children[i]->pm_Neighbors[22] = pm_Neighbors[2]; //DLF
	  pm_Children[i]->pm_Neighbors[23] = pm_Neighbors[15]; //DLB
	  pm_Children[i]->pm_Neighbors[24] = pm_Children[6]; //DRF
	  pm_Children[i]->pm_Neighbors[25] = pm_Neighbors[5]; //DRB
	  break;
	case 2:
	  pm_Children[i]->m_ABC = m_ABC + Point(delta_h,delta_h,delta_h);
	  pm_Children[i]->pm_Neighbors[0] = pm_Neighbors[0]; //U
	  pm_Children[i]->pm_Neighbors[1] = pm_Children[6]; //D
	  pm_Children[i]->pm_Neighbors[2] = pm_Children[0]; //L
	  pm_Children[i]->pm_Neighbors[3] = pm_Neighbors[3]; //R
	  pm_Children[i]->pm_Neighbors[4] = pm_Neighbors[4]; //F
	  pm_Children[i]->pm_Neighbors[5] = pm_Children[3]; //B
	  pm_Children[i]->pm_Neighbors[6] = pm_Neighbors[0]; //UL
	  pm_Children[i]->pm_Neighbors[7] = pm_Neighbors[7]; //UR
	  pm_Children[i]->pm_Neighbors[8] = pm_Neighbors[8]; //UF
	  pm_Children[i]->pm_Neighbors[9] = pm_Neighbors[0]; //UB
	  pm_Children[i]->pm_Neighbors[10] = pm_Children[4]; //DL
	  pm_Children[i]->pm_Neighbors[11] = pm_Neighbors[3]; //DR
	  pm_Children[i]->pm_Neighbors[12] = pm_Neighbors[4]; //DF
	  pm_Children[i]->pm_Neighbors[13] = pm_Children[7]; //DB
	  pm_Children[i]->pm_Neighbors[14] = pm_Neighbors[4]; //LF
	  pm_Children[i]->pm_Neighbors[15] = pm_Children[1]; //LB
	  pm_Children[i]->pm_Neighbors[16] = pm_Neighbors[16]; //RF
	  pm_Children[i]->pm_Neighbors[17] = pm_Neighbors[3]; //RB
	  pm_Children[i]->pm_Neighbors[18] = pm_Neighbors[8]; //ULF
	  pm_Children[i]->pm_Neighbors[19] = pm_Neighbors[0]; //ULB
	  pm_Children[i]->pm_Neighbors[20] = pm_Neighbors[20]; //URF
	  pm_Children[i]->pm_Neighbors[21] = pm_Neighbors[7]; //URB
	  pm_Children[i]->pm_Neighbors[22] = pm_Neighbors[4]; //DLF
	  pm_Children[i]->pm_Neighbors[23] =pm_Children[5]; //DLB
	  pm_Children[i]->pm_Neighbors[24] = pm_Neighbors[16]; //DRF
	  pm_Children[i]->pm_Neighbors[25] = pm_Neighbors[3]; //DRB
	  break;
	case 3:
	  pm_Children[i]->m_ABC = m_ABC + Point(-delta_h,delta_h,delta_h);
	  pm_Children[i]->pm_Neighbors[0] = pm_Neighbors[0]; //U
	  pm_Children[i]->pm_Neighbors[1] = pm_Children[7]; //D
	  pm_Children[i]->pm_Neighbors[2] = pm_Children[1]; //L
	  pm_Children[i]->pm_Neighbors[3] = pm_Neighbors[3]; //R
	  pm_Children[i]->pm_Neighbors[4] = pm_Children[2]; //F
	  pm_Children[i]->pm_Neighbors[5] = pm_Neighbors[5]; //B
	  pm_Children[i]->pm_Neighbors[6] = pm_Neighbors[0]; //UL
	  pm_Children[i]->pm_Neighbors[7] = pm_Neighbors[7]; //UR
	  pm_Children[i]->pm_Neighbors[8] = pm_Neighbors[0]; //UF
	  pm_Children[i]->pm_Neighbors[9] = pm_Neighbors[9]; //UB
	  pm_Children[i]->pm_Neighbors[10] = pm_Children[5]; //DL
	  pm_Children[i]->pm_Neighbors[11] = pm_Neighbors[3]; //DR
	  pm_Children[i]->pm_Neighbors[12] = pm_Children[6]; //DF
	  pm_Children[i]->pm_Neighbors[13] = pm_Neighbors[5]; //DB
	  pm_Children[i]->pm_Neighbors[14] = pm_Children[0]; //LF
	  pm_Children[i]->pm_Neighbors[15] = pm_Neighbors[5]; //LB
	  pm_Children[i]->pm_Neighbors[16] = pm_Neighbors[3]; //RF
	  pm_Children[i]->pm_Neighbors[17] = pm_Neighbors[17]; //RB
	  pm_Children[i]->pm_Neighbors[18] = pm_Neighbors[0]; //ULF
	  pm_Children[i]->pm_Neighbors[19] = pm_Neighbors[9]; //ULB
	  pm_Children[i]->pm_Neighbors[20] = pm_Neighbors[7]; //URF
	  pm_Children[i]->pm_Neighbors[21] = pm_Neighbors[21]; //URB
	  pm_Children[i]->pm_Neighbors[22] = pm_Children[4]; //DLF
	  pm_Children[i]->pm_Neighbors[23] = pm_Neighbors[5]; //DLB
	  pm_Children[i]->pm_Neighbors[24] = pm_Neighbors[3]; //DRF
	  pm_Children[i]->pm_Neighbors[25] = pm_Neighbors[17]; //DRB
	  break;
	case 4:
	  pm_Children[i]->m_ABC = m_ABC + Point(delta_h,-delta_h,delta_h);
	  pm_Children[i]->pm_Neighbors[0] = pm_Children[0]; //U
	  pm_Children[i]->pm_Neighbors[1] = pm_Neighbors[1]; //D
	  pm_Children[i]->pm_Neighbors[2] = pm_Neighbors[2]; //L
	  pm_Children[i]->pm_Neighbors[3] = pm_Children[6]; //R
	  pm_Children[i]->pm_Neighbors[4] = pm_Neighbors[4]; //F
	  pm_Children[i]->pm_Neighbors[5] = pm_Children[5]; //B
	  pm_Children[i]->pm_Neighbors[6] = pm_Neighbors[2]; //UL
	  pm_Children[i]->pm_Neighbors[7] = pm_Children[2]; //UR
	  pm_Children[i]->pm_Neighbors[8] = pm_Neighbors[4]; //UF
	  pm_Children[i]->pm_Neighbors[9] = pm_Children[1]; //UB
	  pm_Children[i]->pm_Neighbors[10] = pm_Neighbors[10]; //DL
	  pm_Children[i]->pm_Neighbors[11] = pm_Neighbors[1]; //DR
	  pm_Children[i]->pm_Neighbors[12] = pm_Neighbors[12]; //DF
	  pm_Children[i]->pm_Neighbors[13] = pm_Neighbors[1]; //DB
	  pm_Children[i]->pm_Neighbors[14] = pm_Neighbors[14]; //LF
	  pm_Children[i]->pm_Neighbors[15] = pm_Neighbors[2]; //LB
	  pm_Children[i]->pm_Neighbors[16] = pm_Neighbors[4]; //RF
	  pm_Children[i]->pm_Neighbors[17] = pm_Children[7]; //RB
	  pm_Children[i]->pm_Neighbors[18] = pm_Neighbors[14]; //ULF
	  pm_Children[i]->pm_Neighbors[19] = pm_Neighbors[2]; //ULB
	  pm_Children[i]->pm_Neighbors[20] = pm_Neighbors[4]; //URF
	  pm_Children[i]->pm_Neighbors[21] = pm_Children[3]; //URB
	  pm_Children[i]->pm_Neighbors[22] = pm_Neighbors[22]; //DLF
	  pm_Children[i]->pm_Neighbors[23] = pm_Neighbors[10]; //DLB
	  pm_Children[i]->pm_Neighbors[24] = pm_Neighbors[12]; //DRF
	  pm_Children[i]->pm_Neighbors[25] = pm_Neighbors[1]; //DRB
	  break;
	case 5:
	  pm_Children[i]->m_ABC = m_ABC + Point(-delta_h,-delta_h,delta_h);
	  pm_Children[i]->pm_Neighbors[0] = pm_Children[1]; //U
	  pm_Children[i]->pm_Neighbors[1] = pm_Neighbors[1]; //D
	  pm_Children[i]->pm_Neighbors[2] = pm_Neighbors[2]; //L
	  pm_Children[i]->pm_Neighbors[3] = pm_Children[7]; //R
	  pm_Children[i]->pm_Neighbors[4] = pm_Children[4]; //F
	  pm_Children[i]->pm_Neighbors[5] = pm_Neighbors[5]; //B
	  pm_Children[i]->pm_Neighbors[6] = pm_Neighbors[2]; //UL
	  pm_Children[i]->pm_Neighbors[7] = pm_Children[3]; //UR
	  pm_Children[i]->pm_Neighbors[8] = pm_Children[0]; //UF
	  pm_Children[i]->pm_Neighbors[9] = pm_Neighbors[5]; //UB
	  pm_Children[i]->pm_Neighbors[10] = pm_Neighbors[10]; //DL
	  pm_Children[i]->pm_Neighbors[11] = pm_Neighbors[1]; //DR
	  pm_Children[i]->pm_Neighbors[12] = pm_Neighbors[1]; //DF
	  pm_Children[i]->pm_Neighbors[13] = pm_Neighbors[13]; //DB
	  pm_Children[i]->pm_Neighbors[14] = pm_Neighbors[2]; //LF
	  pm_Children[i]->pm_Neighbors[15] = pm_Neighbors[15]; //LB
	  pm_Children[i]->pm_Neighbors[16] = pm_Children[6]; //RF
	  pm_Children[i]->pm_Neighbors[17] = pm_Neighbors[5]; //RB
	  pm_Children[i]->pm_Neighbors[18] = pm_Neighbors[2]; //ULF
	  pm_Children[i]->pm_Neighbors[19] = pm_Neighbors[15]; //ULB
	  pm_Children[i]->pm_Neighbors[20] = pm_Children[2]; //URF
	  pm_Children[i]->pm_Neighbors[21] = pm_Neighbors[5]; //URB
	  pm_Children[i]->pm_Neighbors[22] = pm_Neighbors[10]; //DLF
	  pm_Children[i]->pm_Neighbors[23] = pm_Neighbors[23]; //DLB
	  pm_Children[i]->pm_Neighbors[24] = pm_Neighbors[1]; //DRF
	  pm_Children[i]->pm_Neighbors[25] = pm_Neighbors[13]; //DRB
	  break;
	case 6:
	  pm_Children[i]->m_ABC = m_ABC + Point(delta_h,delta_h,delta_h);
	  pm_Children[i]->pm_Neighbors[0] = pm_Children[2]; //U
	  pm_Children[i]->pm_Neighbors[1] = pm_Neighbors[1]; //D
	  pm_Children[i]->pm_Neighbors[2] = pm_Children[4]; //L
	  pm_Children[i]->pm_Neighbors[3] = pm_Neighbors[3]; //R
	  pm_Children[i]->pm_Neighbors[4] = pm_Neighbors[4]; //F
	  pm_Children[i]->pm_Neighbors[5] = pm_Children[7]; //B
	  pm_Children[i]->pm_Neighbors[6] = pm_Children[0]; //UL
	  pm_Children[i]->pm_Neighbors[7] = pm_Neighbors[3]; //UR
	  pm_Children[i]->pm_Neighbors[8] = pm_Neighbors[4]; //UF
	  pm_Children[i]->pm_Neighbors[9] = pm_Children[3]; //UB
	  pm_Children[i]->pm_Neighbors[10] = pm_Neighbors[1]; //DL
	  pm_Children[i]->pm_Neighbors[11] = pm_Neighbors[11]; //DR
	  pm_Children[i]->pm_Neighbors[12] = pm_Neighbors[12]; //DF
	  pm_Children[i]->pm_Neighbors[13] = pm_Neighbors[1]; //DB
	  pm_Children[i]->pm_Neighbors[14] = pm_Neighbors[4]; //LF
	  pm_Children[i]->pm_Neighbors[15] = pm_Children[5]; //LB
	  pm_Children[i]->pm_Neighbors[16] = pm_Neighbors[16]; //RF
	  pm_Children[i]->pm_Neighbors[17] = pm_Neighbors[3]; //RB
	  pm_Children[i]->pm_Neighbors[18] = pm_Neighbors[4]; //ULF
	  pm_Children[i]->pm_Neighbors[19] = pm_Children[1]; //ULB
	  pm_Children[i]->pm_Neighbors[20] = pm_Neighbors[16]; //URF
	  pm_Children[i]->pm_Neighbors[21] = pm_Neighbors[3]; //URB
	  pm_Children[i]->pm_Neighbors[22] = pm_Neighbors[12]; //DLF
	  pm_Children[i]->pm_Neighbors[23] = pm_Neighbors[1]; //DLB
	  pm_Children[i]->pm_Neighbors[24] = pm_Neighbors[24]; //DRF
	  pm_Children[i]->pm_Neighbors[25] = pm_Neighbors[11]; //DRB
	  break;
	case 7:
	  pm_Children[i]->m_ABC = m_ABC + Point(-delta_h,delta_h,delta_h);
	  pm_Children[i]->pm_Neighbors[0] = pm_Children[3]; //U
	  pm_Children[i]->pm_Neighbors[1] = pm_Neighbors[1]; //D
	  pm_Children[i]->pm_Neighbors[2] = pm_Children[5]; //L
	  pm_Children[i]->pm_Neighbors[3] = pm_Neighbors[3]; //R
	  pm_Children[i]->pm_Neighbors[4] = pm_Children[6]; //F
	  pm_Children[i]->pm_Neighbors[5] = pm_Neighbors[5]; //B
	  pm_Children[i]->pm_Neighbors[6] = pm_Children[1]; //UL
	  pm_Children[i]->pm_Neighbors[7] = pm_Neighbors[3]; //UR
	  pm_Children[i]->pm_Neighbors[8] = pm_Children[2]; //UF
	  pm_Children[i]->pm_Neighbors[9] = pm_Neighbors[5]; //UB
	  pm_Children[i]->pm_Neighbors[10] = pm_Neighbors[1]; //DL
	  pm_Children[i]->pm_Neighbors[11] = pm_Neighbors[11]; //DR
	  pm_Children[i]->pm_Neighbors[12] = pm_Neighbors[1]; //DF
	  pm_Children[i]->pm_Neighbors[13] = pm_Neighbors[13]; //DB
	  pm_Children[i]->pm_Neighbors[14] = pm_Children[4]; //LF
	  pm_Children[i]->pm_Neighbors[15] = pm_Neighbors[5]; //LB
	  pm_Children[i]->pm_Neighbors[16] = pm_Neighbors[3]; //RF
	  pm_Children[i]->pm_Neighbors[17] = pm_Neighbors[17]; //RB
	  pm_Children[i]->pm_Neighbors[18] = pm_Children[0]; //ULF
	  pm_Children[i]->pm_Neighbors[19] = pm_Neighbors[5]; //ULB
	  pm_Children[i]->pm_Neighbors[20] = pm_Neighbors[3]; //URF
	  pm_Children[i]->pm_Neighbors[21] = pm_Neighbors[17]; //URB
	  pm_Children[i]->pm_Neighbors[22] = pm_Neighbors[1]; //DLF
	  pm_Children[i]->pm_Neighbors[23] = pm_Neighbors[13]; //DLB
	  pm_Children[i]->pm_Neighbors[24] = pm_Neighbors[11]; //DRF
	  pm_Children[i]->pm_Neighbors[25] = pm_Neighbors[25]; //DRB
	  break;
	}
    }
}

bool Node::Subset(Point object)
{
  //If the center lies within the cube then it is a member
   if ( m_ABC[0]+m_Length > object[0] && m_ABC[0]-m_Length < object[0])
    {
      if ( m_ABC[1]+m_Length > object[1] && m_ABC[1]-m_Length < object[1])
	{
	  if ( m_ABC[2]+m_Length > object[2] && m_ABC[2]-m_Length < object[2])
	    {
	      return true;
	    }
	}
    } 
   return false;
}

bool Node::Subset(Sphere object)
{
  //If the center lies within the cube then it is a member (Of course this needs to be Extended
   if ( m_ABC[0]+m_Length > object.center[0] && m_ABC[0]-m_Length < object.center[0])
    {
      if ( m_ABC[1]+m_Length > object.center[1] && m_ABC[1]-m_Length < object.center[1])
	{
	  if ( m_ABC[2]+m_Length > object.center[2] && m_ABC[2]-m_Length < object.center[2])
	    {
	      return true;
	    }
	}
    } 
   return false;
}

