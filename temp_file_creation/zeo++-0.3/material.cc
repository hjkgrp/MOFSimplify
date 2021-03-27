/* Functions of the MATERIAL class
 * (some fragments taken from main()
 * Added by M. Haranczyk, Feb 2014 
 *
 */

#include <vector>
#include <voro++.hh>

#include "network.h"
#include "networkio.h"
#include "networkinfo.h"
#include "networkanalysis.h"
#include "grid.h"
#include "channel.h"
#include "feature.h"
#include "holograms.h"
#include "instructions.h"
#include "ray.h"
#include "psd.h"
#include "area_and_volume.h"
#include "sphere_approx.h"
#include "poreinfo.h"
#include "zeojobs.h"
#include "material.h"

using namespace std;
using namespace voro;


void MATERIAL::runVoroFlat()
 {
 if(doneFlatVoroFlag == false) {
   // Perform the relevant Voronoi decomposition
   cout << "Starting Voronoi decomposition" << "\n";
   if(radial)
     rad_con = (container_periodic_poly *)performVoronoiDecomp(true, &atmnet, &vornet, cells, outputZvis, bvcells);
   else
     no_rad_con = (container_periodic *)performVoronoiDecomp (false, &atmnet, &vornet, cells, outputZvis, bvcells);

   doneFlatVoroFlag = true;
   cout << "Finished Voronoi decomposition" << "\n";
   };

 };


/* ACCESSIBLE VOLUME FUNCTIONS (incl. blocking&PSD) */

void MATERIAL::AVcalc(double r, int sampleDensity,ostream &output, char *filename){
  runVoroFlat();
  accessAnalysis.AccessibilityClassSetup(&atmnet, &orgAtomnet, highAccuracy, rad_con, &vornet, &bvcells, &cells);
  accessAnalysis.FindChannels(r);
  MATERIAL *pMaterial = this;
  if(!AVdoneFlag) NEWcalcAV(pMaterial, r, sampleDensity);
  AVdoneFlag = true;

  NEWcalcAVprint(pMaterial, output,filename); 
  };

void MATERIAL::AVblockPockets(ostream &output){
  if(AVdoneFlag == false)
    {
    cerr << "Cannot execute blocking before AV run.\n";
    } else
    {
    if(AVblockDoneFlag == false) blockPockets(&atmnet, output, AVaxsPoints, AVaxsPointsChannelIDs, AVinaxsPoints, AVinaxsPointsPocketIDs, AVprobeRadius);
    AVblockDoneFlag = true;
    };
  };


void MATERIAL::AVcalcPoreSizeDistr(ostream &output){
  if(AVdoneFlag == false)
    {
    cerr << "Cannot execute PSD before AV run.\n";
    } else
    {
    MATERIAL *pMaterial = this;
    if(AVPSDdoneFlag == false) NEWcalcPoreSizeDistr(pMaterial, output);
    AVPSDdoneFlag = true;
    };
  };

void MATERIAL::AVreportPoints(ostream &output){
  NEWreportPoints(output, &atmnet, &AVaxsPoints, &AVaxsPointsChannelIDs, &AVinaxsPoints, &AVinaxsPointsPocketIDs, VisSetting);
  };

void MATERIAL::AVreportPSDPoints(ostream &output){
  NEWreportPointsValue(output, &atmnet, &AVaxsPoints, &AVaxsPointsChannelIDs, &AVaxsPointsPSD, VisSetting);
  };

/* SURFACE AREA FUNCTIONS */

void MATERIAL::ASAcalc(double r, int sampleDensity, ostream &output, char *filename){
  runVoroFlat();
  accessAnalysis.AccessibilityClassSetup(&atmnet, &orgAtomnet, highAccuracy, rad_con, &vornet, &bvcells, &cells);
  accessAnalysis.FindChannels(r);
  MATERIAL *pMaterial = this;
  if(!ASAdoneFlag) NEWcalcASA(pMaterial, r, sampleDensity);
  ASAdoneFlag = true;

  NEWcalcASAprint(pMaterial, output,filename);
  };

void MATERIAL::ASAreportPoints(ostream &output){
  NEWreportPoints(output, &atmnet, &ASAaxsPoints, &ASAaxsPointsChannelIDs, &ASAinaxsPoints, &ASAinaxsPointsPocketIDs, VisSetting);
  };


/* visVORO FUNCTION */
void MATERIAL::visualizeVoroNet(char* name, double r, int skel_a, int skel_b, int skel_c){
  visVoro(name, r, skel_a, skel_b, skel_c, &vornet, &atmnet);
  }; // ends visVoro


/* PORE LIMITING ANALYSIS FUNCTIONS */

// PLDcalc uses seg_r to segment the void space into POREs, then takes the network available to r and analysis PLD between the segments identified for seg_r
void MATERIAL::PLDcalc(double r, double seg_r, string seg_file, ostream &output, char *filename){
  runVoroFlat();
  accessAnalysis.AccessibilityClassSetup(&atmnet, &orgAtomnet, highAccuracy, rad_con, &vornet, &bvcells, &cells);
  accessAnalysis.FindChannels(r);

  if(seg_r > 0)
    { // segment the network based on provide radii of void
    accessAnalysis.calculatePLDbasedOnRadius(seg_r);
    accessAnalysis.reportPLD(output);
    } else
    { // segment the network based on provided segment definitions from a file
    accessAnalysis.calculatePLDbasedOnFile(seg_file);
    accessAnalysis.reportPLD(output);    
    };
/*
  MATERIAL *pMaterial = this;
  if(!ASAdoneFlag) NEWcalcASA(pMaterial, r, sampleDensity);
  ASAdoneFlag = true;
  NEWcalcASAprint(pMaterial, output,filename);
*/
  };

// PLDcalc uses molecules present in the periodic box to segment the void space into POREs, then takes the network available to r and analysie PLD between the segments 
void MATERIAL::PLDcalcFromMolecules(double r, ostream &output, char *filename){
  runVoroFlat();
  accessAnalysis.AccessibilityClassSetup(&atmnet, &orgAtomnet, highAccuracy, rad_con, &vornet, &bvcells, &cells);
  accessAnalysis.FindChannels(r);

  accessAnalysis.calculatePLDbasedOnMoleculesPresent(filename);
  accessAnalysis.reportPLD(output);
  };

// Function to save visualization of the calculated PLD
void MATERIAL::PLDvisualize(string basefilename, string visformat){

  fstream output;

  // Three data structures storing information about 
  vector<Point> NodesFracCoord;
  vector<int> NodesSegmentIDs;
  vector<double> NodesSize;

  string filename_vis;
  if(visformat == "ZEOVIS") filename_vis = basefilename + ".zpld_segments";
  if(visformat == "VISIT") filename_vis = basefilename + ".vpld_segments";
  if(visformat == "LIVERPOOL") filename_vis = basefilename + ".lpld_segments";
  output.open(filename_vis.data(), fstream::out);

  accessAnalysis.getPLDvisData(&NodesFracCoord, &NodesSegmentIDs, &NodesSize, "INITSEGMAP");
  NEWreportPointsValue(output, &atmnet, &NodesFracCoord, &NodesSegmentIDs, &NodesSize, visformat);

  output.close();

  if(visformat == "ZEOVIS") filename_vis = basefilename + ".zpld_segmentdi";
  if(visformat == "VISIT") filename_vis = basefilename + ".vpld_segmentdi";
  if(visformat == "LIVERPOOL") filename_vis = basefilename + ".lpld_segmentdi";
  output.open(filename_vis.data(), fstream::out);

  accessAnalysis.getPLDvisData(&NodesFracCoord, &NodesSegmentIDs, &NodesSize, "INITSEGDINODE");
  NEWreportPointsValue(output, &atmnet, &NodesFracCoord, &NodesSegmentIDs, &NodesSize, visformat);

  output.close();

  if(visformat == "ZEOVIS") filename_vis = basefilename + ".zpld_segmentpld";
  if(visformat == "VISIT") filename_vis = basefilename + ".vpld_segmentpld";
  if(visformat == "LIVERPOOL") filename_vis = basefilename + ".lpld_segmentpld";
  output.open(filename_vis.data(), fstream::out);

  accessAnalysis.getPLDvisData(&NodesFracCoord, &NodesSegmentIDs, &NodesSize, "PLDNODES");
  NEWreportPointsValue(output, &atmnet, &NodesFracCoord, &NodesSegmentIDs, &NodesSize, visformat);

  output.close();

  if(visformat == "ZEOVIS") filename_vis = basefilename + ".zpld_segmentdf";
  if(visformat == "VISIT") filename_vis = basefilename + ".vpld_segmentdf";
  if(visformat == "LIVERPOOL") filename_vis = basefilename + ".lpld_segmentdf";
  output.open(filename_vis.data(), fstream::out);

  accessAnalysis.getPLDvisData(&NodesFracCoord, &NodesSegmentIDs, &NodesSize, "DFSPHERES");
  NEWreportPointsValue(output, &atmnet, &NodesFracCoord, &NodesSegmentIDs, &NodesSize, visformat);

  output.close();
}


