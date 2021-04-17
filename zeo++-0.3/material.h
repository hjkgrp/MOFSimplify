/** 
  * In material.* files, there are functions to control execution of Material class
  * This is the main class to store information about a structure and its porosity
  *
  * Initial code by M. Haranczyk, Feb 2014
  *
  **/

#ifndef MATERIAL_H
#define MATERIAL_H

#include <cstdio>
#include <string>
#include "general.h"
#include "geometry.h"
#include "networkstorage.h"
#include "networkaccessibility.h"
//#include "area_and_volume.h"

using namespace voro;

/* MATERIAL class is used to store the structure and perform all its analysis
*/

class MATERIAL {

  public:

    /* definitions */
    ATOM_NETWORK atmnet; ATOM_NETWORK orgAtomnet;
    VORONOI_NETWORK vornet;
    vector< VOR_CELL> cells;
    vector< BASIC_VCELL> bvcells;
    bool radial, outputZvis, buildBasicCells, useMass, allowAdjustCoordsAndCell;

    // High accuracy setting
    bool highAccuracy;
    string AccSetting;

    AccessibilityClassNINF accessAnalysis; 

    // Visualization setting
    bool VisFlag;
    string VisSetting;

    // Declase object req. for Voronoi decomposition

    container_periodic_poly *rad_con;
    container_periodic   *no_rad_con;

    /* General status flags */

    bool doneFlatVoroFlag;
    bool doneInflatedVoroFlag;
    
//    bool doneAccessibilityAnalysis;

    
    /* Data structures storing the results */

    // free sphere stuff
    double Di, Df, Dif;
    vector <double> directionDf; vector <double> directionDif; // Df and Dif along 3 crystallographic directions

    // accessible volume ()

    double AVprobeRadius; // radius of the probe used for sampling (needed for blocking pockets)
    int AVnumSamples; // total number of samples
    int AVcount; // number of samples in accessible space
    int AVcount_inaxs; // number of samples in inaccessible space
    int AVcount_within_range; // number of samples within selected distance range
    int AVcount_outside_range;
    bool AVwithin_rangeFlag;

    vector<Point> AVaxsPoints; // stores accessible points (fractional coordinates)
    vector<int> AVaxsPointsChannelIDs; // stores the corresponding channel IDs 
    vector<Point> AVinaxsPoints; // inaccessible points (fractonal coordinates)
    vector<int> AVinaxsPointsPocketIDs; // stores the corresponding pocket IDs
    vector<int> AVcount_inChannel;  // stores the number of points assigned to particular channel
    vector<int> AVcount_inPocket; // stores the nubmer of points assigned to particular pocket

    vector<double> AVaxsPointsPSD; // stores diameter of pore associated with each accessible point

    bool AVrequestBlockingPockets; // flag that indicates that blocking spheres need to be saved
    bool AVrequestPSD; // flag that indicates that PSD calculation need to be executed

    /* AV related flags */
    bool AVdoneFlag;  
    bool AVPSDdoneFlag;
    bool AVblockDoneFlag;

    // accessible surface area ();
    double ASAprobeRadius; // radius of the probe used for sampling (for inflated atoms)
    int ASAnumSamples;  // total number of samples (sum over all atoms)
    //int ASAcount;
    //int ASAcount_inaxs;
    double ASAtotal; // accessible surface 
    double ASAtotal_inaxs; // inaccessible ASA  

    vector<Point> ASAaxsPoints; // stores accessible points (fractional coordinates)
    vector<int> ASAaxsPointsChannelIDs; // stores the corresponding channel IDs 
    vector<Point> ASAinaxsPoints; // inaccessible points (fractonal coordinates)
    vector<int> ASAinaxsPointsPocketIDs; // stores the corresponding pocket IDs
    vector<double> ASA_inChannel;  // stores the SA assigned to particular channel
    vector<double> ASA_inPocket; // stores the SA assigned to particular pocket

    vector< pair <double,double> > ASAhistogramSAperAtom; // stores SA contributions per atom (accessible/nonaccessible)

    /* ASA relared flags */
    bool ASAdoneFlag;



  /* Init. of the class */
  MATERIAL()
   {
   radial = true; outputZvis = false; buildBasicCells = false; useMass = true; allowAdjustCoordsAndCell = false;

   AVrequestBlockingPockets = false; // by default dont save blocking spheres in AV run
   AVrequestPSD = false; // by default do not run PSD after AV run
 
   highAccuracy = true; AccSetting = "DEF";
   cout << "-ha setting is default now (DEF setting). Please use \"-noha\" to override.\n";

   VisFlag = false;

   doneFlatVoroFlag = false;
   doneInflatedVoroFlag = false;
//   doneAccessibilityAnalysis = false;
   AVdoneFlag = false;
   AVPSDdoneFlag = false;
   AVblockDoneFlag = false;

   ASAdoneFlag = false;

   // Declase object req. for Voronoi decomposition

   container_periodic_poly *rad_con = NULL;
   container_periodic   *no_rad_con = NULL;


   };

  /* Deinit. of the class */
  void deinit()
   {
   accessAnalysis.deconstruct();
   //FREE MEMORY
   delete rad_con;
   delete no_rad_con;
   };

  void runVoroFlat();

  void test(double r){
  runVoroFlat();
  accessAnalysis.AccessibilityClassSetup(&atmnet, &orgAtomnet, highAccuracy, rad_con, &vornet, &bvcells, &cells);
  accessAnalysis.FindChannels(r);
  };

  void AVcalc(double r, int sampleDensity,ostream &output, char *filename);

  void AVblockPockets(ostream &output);
  void AVcalcPoreSizeDistr(ostream &output);

  void AVreportPoints(ostream &output);
  void AVreportPSDPoints(ostream &output);

  void ASAcalc(double r, int sampleDensity,ostream &output, char *filename);
  void ASAreportPoints(ostream &output); 

  void visualizeVoroNet(char *, double, int, int, int);

  void PLDcalc(double r, double seg_r, string seg_file, ostream &output, char *filename);
  void PLDcalcFromMolecules(double r, ostream &output, char *filename);
  void PLDvisualize(string basefilename, string visformat);

 };



#endif
