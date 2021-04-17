
/* 
  Here are functions that determine accessibility of points 

*/

#include "network.h"
#include "channel.h"
#include "networkaccessibility.h"

using namespace std;
using namespace voro;

//
//
///
//
//  PLEASE NOTE THAT BELOW ARE NON-INFLATED VERSIONS OF THESE FUNCTIONS
//  WHICH WILL BECOME STANDARD IN FUTURE RELEASES
//
//
//
//



/* Sets up accessibility class, which is used to determine if a point is accessible or not
   The class need two atomic networks (one is original and another "high accuracy" one
   which has large atoms replaced by small atoms.
   The analysis is done on "inflated" atoms by r_probe_sample
   Detection of channels and accessible pockets is done using r_probe_chan */
void AccessibilityClass::setupAndFindChannels(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgatmnet, bool highAccuracy, double r_probe_chan, double r_probe_sampl){

 r_probe = r_probe_sampl;

 highAccuracyFlag = highAccuracy;

 // Create a temporary copy of the atomic network in which each atom's radius has been increased by the probe radius
 if(highAccuracy)
   {
   atmnet->copy(&analyzedAtomNet);
   orgatmnet->copy(&orgAtomNet);
   }else{
   orgatmnet->copy(&analyzedAtomNet);
   orgatmnet->copy(&orgAtomNet);
   };
 for(unsigned int i = 0; i < orgAtomNet.atoms.size(); i++){ orgAtomNet.atoms[i].radius += r_probe; }
 for(unsigned int i = 0; i < analyzedAtomNet.atoms.size(); i++){ analyzedAtomNet.atoms[i].radius += r_probe; }

 // Calculate and store the Voronoi network for this new atomic network

 new_rad_con = (container_periodic_poly *)performVoronoiDecomp(true, &analyzedAtomNet, &vornet, advCells, false, vorcells);

//  Print vornet summary
 int noAccVorNodes = 0;
 for(unsigned int i = 0; i < vornet.nodes.size(); i++)
    {
    if(vornet.nodes.at(i).rad_stat_sphere > 0) noAccVorNodes++;
    };
 cout << "Voronoi network with " << vornet.nodes.size() << " nodes. " << noAccVorNodes << " of them are accessible. " << endl;

// CHANNEL::findChannels(&vornet, max(0.0, r_probe_chan - r_probe), &accessInfo, &channels);

 PORE::findChannelsAndPockets(&vornet, max(0.0, r_probe_chan - r_probe), &accessInfo, &pores);

 channelMapping.resize(accessInfo.size(),-1);
 pocketMapping.resize(accessInfo.size(),-1);
 n_channels = 0; n_pockets = 0; // number of channels and pockets
 for(unsigned int i = 0; i < pores.size(); i++)
   {
   if(pores[i].dimensionality>0)
    { // Channels
     for(unsigned int j = 0; j < pores[i].nodes.size(); j++)
       {
       channelMapping[pores[i].reverseIDMappings.find(j)->second] = n_channels;
       };
     n_channels++;
    }else
     { // Pockets
     for(unsigned int j = 0; j < pores[i].nodes.size(); j++)
       {
       pocketMapping[pores[i].reverseIDMappings.find(j)->second] = n_pockets;
       };
     n_pockets++;
    };

   };

/* DEBUG

for(unsigned int i = 0; i < accessInfo.size(); i++)
  {
  cout << "DEBUG: NodeID, radius, channelID, pocketID = " << i << "  ,   " << vornet.nodes.at(i).rad_stat_sphere << " ,  " << channelMapping[i] << "  ,  " << pocketMapping[i] << endl;
  };

END DEBUG */

// cout << "Accessibility setup: no channels = " << n_channels << " and no pockets = " << n_pockets << "/n";

} 


/* Remove nodes that will not be used for analysis */
void AccessibilityClass::removeOverlappedNodes(){
  // Remove all nodes from Voronoi cells that lie within the sampling sphere
  for(unsigned int i = 0; i < vorcells.size(); i++){
    vorcells[i].removeOverlappedNodes(i, &analyzedAtomNet, 0); // 0 is probe radious, 0 because atoms are already inflated
    };
}




/* checks if the provided point is in accessible volume */
pair <bool,bool> AccessibilityClass::isVPointInsideAtomAndNotAccessible(Point samplingPoint,double& mindist_domod){
return isPointInsideAtomAndNotAccessible(samplingPoint, mindist_domod, -1);
}

/* checks if the provided point is in accessible volume */
/* return true if a point is accessible */
bool AccessibilityClass::isVPointAccessible(Point samplingPoint){
double mindist_domod;
pair <bool,bool> answer = isPointInsideAtomAndNotAccessible(samplingPoint, mindist_domod, -1);

if(answer.first == false && answer.second == false) return true; else return false;

}

/* checks if thr provided point on a surface of atomID is accessible */ 
pair <bool,bool> AccessibilityClass::isSPointInsideAtomAndNotAccessible(Point samplingPoint, int atomID){
double mindist_domod;
return isPointInsideAtomAndNotAccessible(samplingPoint, mindist_domod, atomID);
}

/* Checks is provide point is inside any atom and if it is accessible */
/* atomID is an atom in the original network from where the sampling point is
   this is to check if the point is inside atom in SA calculaton */
pair <bool,bool> AccessibilityClass::isPointInsideAtomAndNotAccessible(Point samplingPoint, double& mindist_domod, int atomID){

 bool inside = false, overlaps = false; // flags to state if a point is inside atom, inaccessible

 needToResampleFlag = false; // since new point is tested, we clear this variable

 Point smplPoint; // temporary sampling point

 double newAtomX, newAtomY, newAtomZ;
 int minAtomID;
 bool foundCell = new_rad_con->find_voronoi_cell(samplingPoint[0], samplingPoint[1], samplingPoint[2], newAtomX, newAtomY, newAtomZ, minAtomID);
 if(!foundCell){
     cerr << "Error: Unable to find Voronoi cell for sampled point." << "\n"
          << "Exiting..." << "\n";
     exit(1);
    };

 tempMinDistAtomID = minAtomID; // store this information so it is not needed to be recomputed
 tempPoint = samplingPoint;

 // The routine first checks if the provided point is inside atoms


 // if in SA routine, check if the sampled point is within other atom
 if(atomID >= 0){
  if(highAccuracyFlag == false)
    {
    // old check from Thomas code (pre-high accuracy)
 
    // If point in Voronoi cell of different atom, probe-atom overlap occurs because of d^2-r^2 criterion.
    if(minAtomID != atomID)
      overlaps = true;

    }else{
    // new check if high accuracy is requested
    if(analyzedAtomNet.IDmapping[minAtomID] != atomID)
      overlaps = true;
    };
  }; // finishing check for SA


 ATOM curAtom = analyzedAtomNet.atoms[minAtomID];

 // Adjust sampling point so that it lies within the Voronoi cell of interest constructed previously
 smplPoint = (samplingPoint.add(Point(curAtom.x, curAtom.y, curAtom.z).subtract(Point(newAtomX, newAtomY, newAtomZ))));


 double minDist = calcEuclideanDistance(smplPoint[0], smplPoint[1], smplPoint[2], curAtom.x, curAtom.y, curAtom.z);
 if(minDist < curAtom.radius - 0.00000001) 
   overlaps = true;


 if(highAccuracyFlag)  // additional check for high accuracy calculations to check if a point is within the original atom
   {
   curAtom = orgAtomNet.atoms[analyzedAtomNet.IDmapping[minAtomID]];
   minDist = orgAtomNet.calcDistance(smplPoint[0], smplPoint[1], smplPoint[2], &curAtom);
   if(minDist < curAtom.radius - 0.00000001) 
     overlaps = true;
   };

 tempMinDist = minDist; // store temporary (to be used in AV within range function

 inside = overlaps;

 mindist_domod = minDist-curAtom.radius; // DO mod distance for probe-occupiable space

 if(inside == true) return pair<bool,bool> (inside,overlaps); // if the point is inside an atom
                                                              // terminate and return the answer



 // If the point is outside of atoms
 // The routine then checks if the point is in accessible or inaccessible volume/surface 

 curAtom = analyzedAtomNet.atoms[minAtomID];  // making sure we look at the correct atom net (hiAcc or regular)

 // Adjust sampling point so that it lies within the Voronoi cell of interest constructed previously
 samplingPoint = (samplingPoint.add(Point(curAtom.x, curAtom.y, curAtom.z).subtract(Point(newAtomX, newAtomY, newAtomZ))));

 minDist = calcEuclideanDistance(samplingPoint[0], samplingPoint[1], samplingPoint[2], curAtom.x, curAtom.y, curAtom.z);

 // If necessary, check Voronoi nodes of cell to determine accessibility of point
 if(!overlaps){
   BASIC_VCELL vcell = vorcells[minAtomID];
   Point circCenter = Point(curAtom.x, curAtom.y, curAtom.z);
   double samplingRadius = minDist;
   Point sampleRay = Point(samplingPoint[0]-curAtom.x, samplingPoint[1]-curAtom.y, samplingPoint[2]-curAtom.z);

   // Scan the nodes in the Voronoi cell to find if line can be drawn from the node to the sampling point
   bool foundNode = false;
   needToResampleFlag = false;
   if(vcell.getNumNodes() == 0){
     cerr << "Error: Voronoi cell of sampled point does not have any nodes" << "\n"
          << "Point: " << samplingPoint[0] << " " << samplingPoint[1] << " " << samplingPoint[2] << "\n"
          << "Voronoi cell is #" << minAtomID << "\n"
          << "Please contact the source code provider." << "\n"
          << "Exiting..." << "\n";
     exit(1);
   }
   for(int k = 0; k < vcell.getNumNodes(); k++){
           Point nodePoint = vcell.getNodeCoord(k);
           bool nodeInsideSphere = (calcEuclideanDistance(nodePoint[0], nodePoint[1], nodePoint[2], circCenter[0], circCenter[1], circCenter[2]) < samplingRadius);
           bool nodeInsideOtherAtom = (vornet.nodes[vcell.getNodeID(k)].rad_stat_sphere < 0.0 );
           if(!nodeInsideSphere&&!nodeInsideOtherAtom){
             Point otherRay = samplingPoint.subtract(nodePoint);
             double dotProduct = sampleRay.dot_product(otherRay);
             if(dotProduct > 0) {
               // Angle is less than 90 degrees and so the line segment intersects twice,
               // making the path not viable
             }
             else {
               // Angle is at least 90 degrees and so the line segment interesects only once, 
               // thereby representing a viable path
               foundNode = true;
               tempNodeID = vcell.getNodeID(k);
               overlaps = !accessInfo.at(vcell.getNodeID(k));
               break;
             }
           }
    }

    // Sampling failed due to lying on Voronoi cell face and numerical inaccurarcy. 
    // Record failure, resample and notify user later
    if(!foundNode){
       resampleCount++;
       resampledInfo.push_back(pair<int, Point> (minAtomID, samplingPoint));
       needToResampleFlag = true;
       };

 };

return pair<bool,bool> (inside,overlaps);

}

// ///////////////////////////////////////////////////////
// 
//
//
/* BELOW ARE NON-INFLATED VERSIONS OF FUNCTIONS */
//
//
//
// ////////////////////////////////////////////////////////

void AccessibilityClassNINF::AccessibilityClassSetup(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgatmnet, bool highAccuracy,
                                                    voro::container_periodic_poly *p,VORONOI_NETWORK *v, std::vector<BASIC_VCELL> *bvc, std::vector<VOR_CELL> *vc){
 highAccuracyFlag = highAccuracy;
 if(highAccuracy)
   {
   analyzedAtomNet = atmnet;
   orgAtomNet = orgatmnet;
   }else{
   analyzedAtomNet = orgatmnet;
   orgAtomNet = orgatmnet;
   };

 /* Making a local copy of the Voronoi Network (because it will be modified during the channel detection */
 new_rad_con = p;
 vornet = *v;
 advCells = *vc;
 vorcells = *bvc;
 };


void AccessibilityClassNINF::FindChannels(double r_pr){

 if(alreadySegmentedFlag == true)
   {
   if(r_pr != r_probe)
     {
     cerr << "Trying to segment the net that has been already semented with different r\n"
          << "(use new accessibility class)\n";
     abort();
     };
   } else
   { // everything seems OK; proceed with segmentation

   r_probe = r_pr;

   //  Print vornet summary
   int noAccVorNodes = 0;
   for(unsigned int i = 0; i < vornet.nodes.size(); i++)
      {
      if(vornet.nodes.at(i).rad_stat_sphere > r_probe) noAccVorNodes++;
      };
   cout << "Voronoi network with " << vornet.nodes.size() << " nodes. " << noAccVorNodes << " of them are accessible. " << endl;

   PORE::findChannelsAndPockets(&vornet, r_probe, &accessInfo, &pores);

   channelMapping.resize(accessInfo.size(),-1);
   pocketMapping.resize(accessInfo.size(),-1);
   n_channels = 0; n_pockets = 0; // number of channels and pockets
   for(unsigned int i = 0; i < pores.size(); i++)
      {
      if(pores[i].dimensionality>0)
        { // Channels
        for(unsigned int j = 0; j < pores[i].nodes.size(); j++)
          {
          channelMapping[pores[i].reverseIDMappings.find(j)->second] = n_channels;
          };
        n_channels++;
        }else
        { // Pockets
        for(unsigned int j = 0; j < pores[i].nodes.size(); j++)
          {
          pocketMapping[pores[i].reverseIDMappings.find(j)->second] = n_pockets;
          };
        n_pockets++;
        };
     };
   alreadySegmentedFlag = true;
   }; // ends else{

/* DEBUG

for(unsigned int i = 0; i < accessInfo.size(); i++)
  {
  cout << "DEBUG: NodeID, radius, channelID, pocketID = " << i << "  ,   " << vornet.nodes.at(i).rad_stat_sphere << " ,  " << channelMapping[i] << "  ,  " << pocketMapping[i] << endl;
  };

END DEBUG */

// cout << "Accessibility setup: no channels = " << n_channels << " and no pockets = " << n_pockets << "/n";

} 


/* Remove nodes that will not be used for analysis */
void AccessibilityClassNINF::removeOverlappedNodes(){
  // Remove all nodes from Voronoi cells that lie within the sampling sphere
  for(unsigned int i = 0; i < vorcells.size(); i++){
    vorcells[i].removeOverlappedNodes(i, analyzedAtomNet, r_probe); // 0 is probe radious, 0 because atoms are already inflated
                                                                    // changes to r_probe as we moved away from inflating
    };
}




/* checks if the provided point is in accessible volume */
pair <bool,bool> AccessibilityClassNINF::isVPointInsideAtomAndNotAccessible(Point samplingPoint){
return isPointInsideAtomAndNotAccessible(samplingPoint, -1);
}

/* checks if the provided point is in accessible volume */
/* return true if a point is accessible */
bool AccessibilityClassNINF::isVPointAccessible(Point samplingPoint){
pair <bool,bool> answer = isPointInsideAtomAndNotAccessible(samplingPoint, -1);

if(answer.first == false && answer.second == false) return true; else return false;

}

/* checks if thr provided point on a surface of atomID is accessible */ 
pair <bool,bool> AccessibilityClassNINF::isSPointInsideAtomAndNotAccessible(Point samplingPoint, int atomID){
return isPointInsideAtomAndNotAccessible(samplingPoint, atomID);
}

/* Checks if the provided point is inside any atom and if it is accessible */
/* atomID is an atom in the original network from where the sampling point is
   this is to check if the point is inside atom in SA calculaton */
pair <bool,bool> AccessibilityClassNINF::isPointInsideAtomAndNotAccessible(Point samplingPoint, int atomID){

 bool inside = false, overlaps = false; // flags to state if a point is inside atom, inaccessible

 needToResampleFlag = false; // since new point is tested, we clear this variable

 Point smplPoint; // temporary sampling point

 double newAtomX, newAtomY, newAtomZ;
 int minAtomID;
 bool foundCell = new_rad_con->find_voronoi_cell(samplingPoint[0], samplingPoint[1], samplingPoint[2], newAtomX, newAtomY, newAtomZ, minAtomID);
 if(!foundCell){
     cerr << "Error: Unable to find Voronoi cell for sampled point." << "\n"
          << "Exiting..." << "\n";
     exit(1);
    };

 tempMinDistAtomID = minAtomID; // store this information so it is not needed to be recomputed
 tempPoint = samplingPoint;

 // The routine first checks if the provided point is inside atoms


 // if in SA routine, check if the sampled point is within other atom
 if(atomID >= 0){
  if(highAccuracyFlag == false)
    {
    // old check from Thomas code (pre-high accuracy)
 
    // If point in Voronoi cell of different atom, probe-atom overlap occurs because of d^2-r^2 criterion.
    if(minAtomID != atomID)
      overlaps = true;

    }else{
    // new check if high accuracy is requested
    if(analyzedAtomNet->IDmapping[minAtomID] != atomID)
      overlaps = true;
    };
  }; // finishing check for SA


 ATOM curAtom = analyzedAtomNet->atoms[minAtomID];

 // Adjust sampling point so that it lies within the Voronoi cell of interest constructed previously
 smplPoint = (samplingPoint.add(Point(curAtom.x, curAtom.y, curAtom.z).subtract(Point(newAtomX, newAtomY, newAtomZ))));


 double minDist = calcEuclideanDistance(smplPoint[0], smplPoint[1], smplPoint[2], curAtom.x, curAtom.y, curAtom.z);
 if(minDist < curAtom.radius + r_probe - 0.00000001) 
   overlaps = true;


 if(highAccuracyFlag)  // additional check for high accuracy calculations to check if a point is within the original atom
   {
   curAtom = orgAtomNet->atoms[analyzedAtomNet->IDmapping[minAtomID]];
   minDist = orgAtomNet->calcDistance(smplPoint[0], smplPoint[1], smplPoint[2], &curAtom);
   if(minDist < curAtom.radius + r_probe - 0.00000001) 
     overlaps = true;
   };

 tempMinDist = minDist; // store temporary (to be used in AV within range function

 inside = overlaps;

 if(inside == true) return pair<bool,bool> (inside,overlaps); // if the point is inside an atom
                                                              // terminate and return the answer



 // If the point is outside of atoms
 // The routine then checks if the point is in accessible or inaccessible volume/surface 

 curAtom = analyzedAtomNet->atoms[minAtomID];  // making sure we look at the correct atom net (hiAcc or regular)

 // Adjust sampling point so that it lies within the Voronoi cell of interest constructed previously
 samplingPoint = (samplingPoint.add(Point(curAtom.x, curAtom.y, curAtom.z).subtract(Point(newAtomX, newAtomY, newAtomZ))));

 minDist = calcEuclideanDistance(samplingPoint[0], samplingPoint[1], samplingPoint[2], curAtom.x, curAtom.y, curAtom.z);

 // If necessary, check Voronoi nodes of cell to determine accessibility of point
 if(!overlaps){
   BASIC_VCELL vcell = vorcells[minAtomID];
   Point circCenter = Point(curAtom.x, curAtom.y, curAtom.z);
   double samplingRadius = minDist;
   Point sampleRay = Point(samplingPoint[0]-curAtom.x, samplingPoint[1]-curAtom.y, samplingPoint[2]-curAtom.z);

   // Scan the nodes in the Voronoi cell to find if line can be drawn from the node to the sampling point
   bool foundNode = false;
   needToResampleFlag = false;
   if(vcell.getNumNodes() == 0){
     cerr << "Error: Voronoi cell of sampled point does not have any nodes" << "\n"
          << "Point: " << samplingPoint[0] << " " << samplingPoint[1] << " " << samplingPoint[2] << "\n"
          << "Voronoi cell is #" << minAtomID << "\n"
          << "Please contact the source code provider." << "\n"
          << "Exiting..." << "\n";
     exit(1);
   }
   for(int k = 0; k < vcell.getNumNodes(); k++){
           Point nodePoint = vcell.getNodeCoord(k);
           bool nodeInsideSphere = (calcEuclideanDistance(nodePoint[0], nodePoint[1], nodePoint[2], circCenter[0], circCenter[1], circCenter[2]) < samplingRadius);
           bool nodeInsideOtherAtom = (vornet.nodes[vcell.getNodeID(k)].rad_stat_sphere < r_probe );
           if(!nodeInsideSphere&&!nodeInsideOtherAtom){
             Point otherRay = samplingPoint.subtract(nodePoint);
             double dotProduct = sampleRay.dot_product(otherRay);
             if(dotProduct > 0) {
               // Angle is less than 90 degrees and so the line segment intersects twice,
               // making the path not viable
             }
             else {
               // Angle is at least 90 degrees and so the line segment interesects only once, 
               // thereby representing a viable path
               foundNode = true;
               tempNodeID = vcell.getNodeID(k);
               overlaps = !accessInfo.at(vcell.getNodeID(k));
               break;
             }
           }
    }

    // Sampling failed due to lying on Voronoi cell face and numerical inaccurarcy. 
    // Record failure, resample and notify user later
    if(!foundNode){
       resampleCount++;
       resampledInfo.push_back(pair<int, Point> (minAtomID, samplingPoint));
       needToResampleFlag = true;
       };

 };

return pair<bool,bool> (inside,overlaps);

}


/* perform segmentation of the Vornet */
void AccessibilityClassNINF::segmentPoresBasedOnRadius(double seg_r){

 if(alreadySegmentedFlag == false || seg_r<=r_probe)
   {
   cerr << "This function requires initial segmnetation (use new accessibility class(NINF))\n"
        << "the source has most likely bugs. Contact the source code provider.\n"
        << "Also segmenting r needs to be larger than one used for inital segmentation\n";
   abort();
   } else
   { // everything seems OK; proceed with segmentation

   VORONOI_NETWORK localVorNet = vornet;
   std::vector<bool> localAccessInfo;

   std::vector<PORE> segments;

   PORE::findChannelsAndPockets(&localVorNet, seg_r, &localAccessInfo, &segments);

   segmentMapping.resize(localAccessInfo.size(),-1);
   n_segments = 0; // number of segments
   for(unsigned int i = 0; i < segments.size(); i++)
      {
      for(unsigned int j = 0; j < segments[i].nodes.size(); j++)
        {
        segmentMapping[segments[i].reverseIDMappings.find(j)->second] = n_segments;
        };
      n_segments++;
     };

   segments.clear();

   }; // ends else{

 cout << "Additional segmentation: n_segments = " << n_segments << "\n";

} 


/* Performs segmentation of the Vornet based on segments read from a provided file
   These segments are defined by a set of spheres that overlap with node centers
   Format: fractional coordinate, ID, diamter */
bool AccessibilityClassNINF::segmentPoresBasedOnFile(string filename){

  fstream input;
  input.open(filename.c_str());

  int linecount=0;
  int maxid=-1;

  if(!input.is_open()){
    cout << "\n" << "Failed to open segment file  " << filename << "\n";
    cout << "Exiting function ..." << "\n";
    return false;
  }
  else{
    cout << "Reading segment file " << filename << "\n";

    segmentMapping.resize(vornet.nodes.size(), -1);

    // Read and store information about sphere defining segments 
    while(!input.eof()){
      double a,b,c, r;
      int id;

      input >> a >> b >> c >> id >> r;

      if(input.eof()) 
        break;

      r = r*0.5;
      if(id>maxid) maxid = id;

      for(unsigned int i=0; i< vornet.nodes.size(); i++)
         {
         if(vornet.nodes.at(i).active == true)
           {
           if(orgAtomNet->calcDistanceXYZABC(vornet.nodes.at(i).x, vornet.nodes.at(i).y ,vornet.nodes.at(i).z, a,b,c) <= r) 
             {
             // node is within segment sphere
             if(segmentMapping[i] != -1 && segmentMapping[i] != id)
               {
               cerr << "Voronoi node has been assigned to more than one segment. This is wrong.\n";
               abort();

               }else{
               segmentMapping[i] = id;
               };
             };
           };
         };

      linecount++;

    }
    input.close();
  }
  cout << "Segment file: " << linecount << " lines read.";
  cout << "Max segment ID = " << maxid << "\n";

  n_segments = maxid+1; // this assumes that the segmentation file has correct IDs (from 0 to maxid)

  return true;
}



/* Performs segmentation of the Vornet based on molecules present in the unit cell 
   These segments are defined by the volume bounding the molecules as defined in Ismael's code
    */
bool AccessibilityClassNINF::segmentPoresBasedOnMoleculesPresent(string filename){

/* EDIT THIS FUNCTION
  fstream input;
  input.open(filename.c_str());

  int linecount=0;
  int maxid=-1;

  if(!input.is_open()){
    cout << "\n" << "Failed to open segment file  " << filename << "\n";
    cout << "Exiting function ..." << "\n";
    return false;
  }
  else{
    cout << "Reading segment file " << filename << "\n";

    segmentMapping.resize(vornet.nodes.size(), -1);

    // Read and store information about sphere defining segments 
    while(!input.eof()){
      double a,b,c, r;
      int id;

      input >> a >> b >> c >> id >> r;

      if(input.eof())
        break;

      r = r*0.5;
      if(id>maxid) maxid = id;

      for(unsigned int i=0; i< vornet.nodes.size(); i++)
         {
         if(vornet.nodes.at(i).active == true)
           {
           if(orgAtomNet->calcDistanceXYZABC(vornet.nodes.at(i).x, vornet.nodes.at(i).y ,vornet.nodes.at(i).z, a,b,c) <= r)
             {
             // node is within segment sphere
             if(segmentMapping[i] != -1 && segmentMapping[i] != id)
               {
               cerr << "Voronoi node has been assigned to more than one segment. This is wrong.\n";
               abort();

               }else{
               segmentMapping[i] = id;
               };
             };
           };
         };

      linecount++;

    }
    input.close();
  }
  cout << "Segment file: " << linecount << " lines read.";
  cout << "Max segment ID = " << maxid << "\n";

  n_segments = maxid+1; // this assumes that the segmentation file has correct IDs (from 0 to maxid)
*/
  return true;
} // ends segmentPoresBasedOnMoleculesPresent()



/* analyzed connectivity between segments to get PLDs 
 * segments are defined by voids larger than specified r */
void AccessibilityClassNINF::calculatePLDbasedOnRadius(double seg_r){

 segmentPoresBasedOnRadius(seg_r);

 calculatePLD();

} // ends calculatePLDbasedOnRadius()

/* analyzed connectivity between segments to get PLDs 
 * segments are defined by a list of spheres in a file  */
void AccessibilityClassNINF::calculatePLDbasedOnFile(string segment_filename){

 segmentPoresBasedOnFile(segment_filename);

 calculatePLD();

} // ends calculatePLDbasedOnRadius()


/* analyzes connectivity between segments  to get PLDs 
 * segments are defined by the segments of the void space defined by molecules present in the system  */
void AccessibilityClassNINF::calculatePLDbasedOnMoleculesPresent(string structure_filename){

 segmentPoresBasedOnMoleculesPresent(structure_filename);

 calculatePLD();

} // ends calculatePLDbasedOnRadius()

/* analyze connectivity between segments to get PLD */
void AccessibilityClassNINF::calculatePLD(){

 if(n_segments < 2) 
   { // exit is there are no pre-defined segments
   cerr << "Number of segments used as seed for flood fill algorithm is lower than 2\n";
   return;
   } else
   {
   // setup results tables
   PLDtable.clear();
   PLDEdgeTable.clear();

   segmentDi.resize(n_segments, -1);
   segmentDiFinal.resize(n_segments, -1);

   segmentDiNodeID.resize(n_segments, -1);
   segmentDiFinalNodeID.resize(n_segments, -1);

   vector<double> restrDiam;
   restrDiam.resize(n_segments, -1);
   pair<int,int> dummy_pair (-1,-1);
   vector< pair<int,int> > restrEdge;
   restrEdge.resize(n_segments, dummy_pair);

   for(int i=0; i < n_segments; i++) 
     {
     PLDtable.push_back(restrDiam);
     PLDEdgeTable.push_back(restrEdge);
     };

   cout << "There are " << n_segments << " in PLD calculation.\n";

   for(unsigned int i = 0; i < pores.size(); i++)
     {
     cout << "Analyzing pore " << i << " of size " << pores[i].nodes.size() << " nodes.\n";
     pores[i].getRestrictingDiameters(n_segments, segmentMapping, &PLDtable, &PLDEdgeTable, &segmentDi, &segmentDiNodeID, &segmentDiFinal, &segmentDiFinalNodeID);
     };
  }; // ends if(n_segments < 2) ...
}

/* prints PLD table into the provided output */
void AccessibilityClassNINF::reportPLD(ostream &output){ 

 output <<  n_segments << " segments\n";

 output << "segmentDi: ";
 for(int i=0; i < n_segments; i++) output << segmentDiFinal[i] << " ";

 output << "\n====Pairwise===PLDs===between===segments==============\n";
 output.setf(ios::fixed);
 output.precision(3);
 for(int i=0; i < n_segments; i++)
    {
    for(int j=0; j < n_segments; j++)
      {
      double pld;
      if(PLDtable[i].at(j)==-1) pld = 0.000; else pld = PLDtable[i].at(j);
      output <<  pld << "   ";
      };
    output << "\n";
    };

}


/* uploads data required for visualization to the provided vector structures */
/* string reqvisdata indicates which data points are requested */
void AccessibilityClassNINF::getPLDvisData(vector<Point> *NodeFracCoord, vector<int> *NodeSegmentIDs, vector<double> *NodeSegmentValue, string reqvisdata){

 // preparing data structures
 NodeFracCoord->clear();
 NodeSegmentIDs->clear();
 NodeSegmentValue->clear();

 if(reqvisdata == "INITSEGMAP"){
   // initial segment mapping (prints nodes initially assigned to semgents)
   for(unsigned int i=0; i<vornet.nodes.size(); i++)
     {
     if(segmentMapping[i]>=0) 
       {
       NodeFracCoord->push_back(orgAtomNet->xyz_to_abc(vornet.nodes.at(i).x, vornet.nodes.at(i).y, vornet.nodes.at(i).z));
       NodeSegmentIDs->push_back(segmentMapping[i]);
       NodeSegmentValue->push_back(segmentDi[segmentMapping[i]]);
       };
     };
   };

 if(reqvisdata == "FINALSEGMAP"){
   // initial segment mapping (prints nodes initially assigned to semgents)
   for(unsigned int i=0; i<vornet.nodes.size(); i++)
     {
     if(segmentMappingFinal[i]>=0)
       {
       NodeFracCoord->push_back(orgAtomNet->xyz_to_abc(vornet.nodes.at(i).x, vornet.nodes.at(i).y, vornet.nodes.at(i).z));
       NodeSegmentIDs->push_back(segmentMappingFinal[i]);
       NodeSegmentValue->push_back(segmentDiFinal[segmentMappingFinal[i]]);
       };
     };
   };

 if(reqvisdata == "INITSEGDINODE"){
   // initial segment mapping (prints nodes initially assigned to semgents)
   for(unsigned int j=0; j<segmentDiNodeID.size(); j++)
     {
     int i=segmentDiNodeID[j];
     NodeFracCoord->push_back(orgAtomNet->xyz_to_abc(vornet.nodes.at(i).x, vornet.nodes.at(i).y, vornet.nodes.at(i).z));
     NodeSegmentIDs->push_back(j);
     NodeSegmentValue->push_back(segmentDi[j]);
     };
   };

 if(reqvisdata == "PLDNODES"){
   // initial segment mapping (prints nodes initially assigned to semgents)
   for(int i=0; i<n_segments; i++)
   for(int j=i+1; j<n_segments; j++)
     {
     if(PLDtable[i].at(j)> 0.0)
       { // segments i and j are connected
       int node1 = PLDEdgeTable.at(i).at(j).first;
       int node2 = PLDEdgeTable.at(i).at(j).second;  
       for(unsigned int k=0; k<vornet.edges.size(); k++){
       int count=0;
       if((vornet.edges.at(k).from==node1&&vornet.edges.at(k).to==node2)||(vornet.edges.at(k).from==node2&&vornet.edges.at(k).to==node1))
         { // edge found
         if(count > 0) cerr << "Two nodes connected by more than one edge. Visualization may contain artefacts.\n";
         XYZ nd1(vornet.nodes.at(vornet.edges.at(k).from).x, vornet.nodes.at(vornet.edges.at(k).from).y, vornet.nodes.at(vornet.edges.at(k).from).z);
         XYZ nd2(vornet.nodes.at(vornet.edges.at(k).to).x, vornet.nodes.at(vornet.edges.at(k).to).y, vornet.nodes.at(vornet.edges.at(k).to).z);

         nd2 = nd2 + vornet.v_a.scale(vornet.edges.at(k).delta_uc_x) + vornet.v_b.scale(vornet.edges.at(k).delta_uc_y) + vornet.v_c.scale(vornet.edges.at(k).delta_uc_z);

         XYZ nd = midpoint(nd1, nd2);
         
         NodeFracCoord->push_back(orgAtomNet->xyz_to_abc(nd.x, nd.y, nd.z));
         NodeSegmentIDs->push_back(0); // no ID for edge
         NodeSegmentValue->push_back(PLDtable[i].at(j));
         count++;
         };
       };
     };
     };
   };

 if(reqvisdata == "DFSPHERES"){
   // display spheres in between connected cages
   // this is a temp function without PBC 
   for(int i=0; i<n_segments; i++)
   for(int j=i+1; j<n_segments; j++)
     {
     int dinode1 =segmentDiNodeID[i];
     int dinode2 =segmentDiNodeID[j];

     if(PLDtable[i].at(j)> 0.0 && calcEuclideanDistance(vornet.nodes.at(dinode1).x, vornet.nodes.at(dinode1).y, vornet.nodes.at(dinode1).z,
                               vornet.nodes.at(dinode2).x, vornet.nodes.at(dinode2).y, vornet.nodes.at(dinode2).z) < 15.0 )
       { // segments i and j are connected
       int count=0;
         if(count > 0) cerr << "Two nodes connected by more than one edge. Visualization may contain artefacts.\n";

         XYZ nd1(vornet.nodes.at(dinode1).x, vornet.nodes.at(dinode1).y, vornet.nodes.at(dinode1).z);
         XYZ nd2(vornet.nodes.at(dinode2).x, vornet.nodes.at(dinode2).y, vornet.nodes.at(dinode2).z);

         XYZ nd = midpoint(nd1, nd2);

         NodeFracCoord->push_back(orgAtomNet->xyz_to_abc(nd.x, nd.y, nd.z));
         NodeSegmentIDs->push_back(0); // no ID for edge
         NodeSegmentValue->push_back(PLDtable[i].at(j));
         count++;
       };
     };
   };


}

