#ifndef ACCESSIBILITY_H
#define ACCESSIBILITY_H

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <voro++.hh>

#include "networkstorage.h"
#include "geometry.h"
#include "voronoicell.h"
#include "channel.h"

/* Accessibility class handels determination of accessibility of the void space and MC-sampled points */

/* This class also handles high-accuracy calcluations where large atoms are replaced by clusters of small atoms
 when settign up this class with setupAndFindChannels() function, the first atom net provided is the high accuracy
 one based on clusters, the second is the original atom network. If the high accuracy flag is not set, both
 atom network pointers should point to the original net */

/* This version used inflated atoms; WILL BECOME OBSOLETE */

class AccessibilityClass{
    
public:
    
    ATOM_NETWORK orgAtomNet;
    ATOM_NETWORK analyzedAtomNet; // this to store the network to be analyzed (either original or the high accuracy one)
    
    bool highAccuracyFlag;
    
    VORONOI_NETWORK vornet;
    std::vector<BASIC_VCELL> vorcells;
    std::vector<VOR_CELL> advCells;
    
    std::vector<PORE> pores;
    //  vector<CHANNEL> channels;
    //  vector<POCKET> pockets;
    int n_channels, n_pockets;
    std::vector<bool> accessInfo; // flags stating if nodes are accessible
    std::vector<int> channelMapping; // maps node IDs to channel IDs
    std::vector<int> pocketMapping;  // maps node IDs to pocket IDs
    
    double r_probe;
    
    voro::container_periodic_poly *new_rad_con;
    
    double tempMinDist; // temporary variable to store min. dis. in accessibility calcluations
    Point tempPoint; // this array stores coordiantes of the last sampled point
    int tempMinDistAtomID; // ID of the closest atom to the last investigated point
    int tempNodeID;  // ID of the voronoi node used to determine accessibility
    
    std::vector< std::pair<int, Point> > resampledInfo; // List of resampled points and the id of the Voronoi cell to which they belong
    int resampleCount;
    bool needToResampleFlag; // this flag is needed after accessibility functions have been moved to a separate class
    
public:
    
    
    AccessibilityClass(){needToResampleFlag = false; resampleCount = 0;};
    
    void setupAndFindChannels(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgatmnet, bool highAccuracy, double r_probe_chan, double r_probe_sampl);
    
    /* checks if the provided point is in accessible volume */
    std::pair <bool,bool> isVPointInsideAtomAndNotAccessible(Point samplingPoint, double& mindist_domod);
    
    /* checks if thr provided point on a surface of atomID is accessible */
    std::pair <bool,bool> isSPointInsideAtomAndNotAccessible(Point samplingPoint, int atomID);
    
    bool isVPointAccessible(Point samplingPoint);
    
    /* Checks is provide point is inside any atom and if it is accessible */
    /* the functions returns value = true if atoms overlaps/is not accessible) */
    /* atomID is an atom in the original network from where the sampling point is
     this is to check if the point is inside atom in SA calculaton */
    /* DO mod also stores distance to the nearest atom for further evaluation
     * in probe-occupiable volume */
    std::pair <bool,bool> isPointInsideAtomAndNotAccessible(Point samplingPoint, double& mindist_domod, int atomID);
    
    /* Returns the last calculated minDist (calculated in isPointInsideAtomAndNotAccessible() )  */
    double lastMinDist(){return tempMinDist;};
    /* Returns channel or pocket ID for the last calcluated point */
    std::pair <int,int> lastChannelOrPocket()
    {
        if(channelMapping[tempNodeID]<0 && pocketMapping[tempNodeID]<0){
            std::cerr << "CoP_error: cannot determine point accessiblity. Consider running at higher accuracy (-ha flag).(minDist=" << tempMinDist << ")." << std::endl;
            std::cerr << "NodeID= " << tempNodeID << "  minDistAtomID= " << tempMinDistAtomID << " resampleFlag= " << needToResampleFlag << std::endl;
            abort();
        }
        return std::pair <int,int> (channelMapping[tempNodeID],pocketMapping[tempNodeID]);
    };
    
    /* Returns information if a point need to be resampled due to node not found in isPointInsideAtomAndNotAccessible() */
    bool needToResample()
    {
/*
        if(needToResampleFlag) {
            needToResampleFlag = false; // change the value after first used of the function
            return true;
        }
        else return false;
*/

        if(needToResampleFlag==true) std::cout << "Resample flag is raised. Resample count = " << resampleCount << std::endl;

        return needToResampleFlag;
    }
    
    /* Return number of resampled ponits */
    int getResampleCount(){return resampleCount;};
    
    /* remove Voronoi nodes that are not used in analysis */
    void removeOverlappedNodes();
    
    /* deconstruct class */
    void deconstruct(){delete new_rad_con;};
};

// ////////////////////////////////////////////////////////////
//
/* NEW VERSION OF THE CLASS THAT DOES NOT USE INFLATED ATOMS */
/* WILL TOTALLY REPLACE THE OLD VERSION WITH INFLATED ATOMS */
//
// ////////////////////////////////////////////////////////////
class AccessibilityClassNINF{
    
public:
    
    ATOM_NETWORK *orgAtomNet;
    ATOM_NETWORK *analyzedAtomNet; // this to store the network to be analyzed (either original or the high accuracy one)
    
    bool highAccuracyFlag;
    
    VORONOI_NETWORK vornet;
    std::vector<BASIC_VCELL> vorcells;
    std::vector<VOR_CELL> advCells;
    
    std::vector<PORE> pores;
    //  vector<CHANNEL> channels;
    //  vector<POCKET> pockets;
    int n_channels, n_pockets;
    std::vector<bool> accessInfo; // flags stating if nodes are accessible
    std::vector<int> channelMapping; // maps node IDs to channel IDs
    std::vector<int> pocketMapping;  // maps node IDs to pocket IDs
    
    double r_probe;
    
    voro::container_periodic_poly *new_rad_con;
    
    double tempMinDist; // temporary variable to store min. dis. in accessibility calcluations
    Point tempPoint; // this array stores coordiantes of the last sampled point
    int tempMinDistAtomID; // ID of the closest atom to the last investigated point
    int tempNodeID;  // ID of the voronoi node used to determine accessibility
    
    std::vector< std::pair<int, Point> > resampledInfo; // List of resampled points and the id of the Voronoi cell to which they belong
    int resampleCount;
    bool needToResampleFlag; // this flag is needed after accessibility functions have been moved to a separate class

    /* Status Flags */
    bool alreadySegmentedFlag;


    /* Functions and variables allowing additional segmentation of the Voronoi network
     * which may be required for PLD-type of calculations */
    int n_segments; 
    std::vector<int>segmentMapping; // maps node ID to segments (inital), this mapping is used to seed flood-fill-like algorithm
                                    // that performs segmentation of the network (the latter typically after inital accessibility analysis)
    std::vector<int>segmentMappingFinal; // maps node ID to segments (final, after PLD analysis)

    std::vector<double>segmentDi; // stores Di for each segment (at inital step)
    std::vector<int>segmentDiNodeID; // stores node ID corresponding to Di for the corresponding segment
    std::vector<double>segmentDiFinal; // stores Di for each segment (after flood-fill)
    std::vector<int>segmentDiFinalNodeID; // stores node ID corresponding to the Final Di for the corresponding segment

    std::vector< vector<double> > PLDtable; // stores PLD (restricting diameters) between segments
    std::vector< vector< pair<int,int> > > PLDEdgeTable; // stores pair of ints defining segment-connecting edges (node id1, node2)
    
public:
    
    
    AccessibilityClassNINF(){needToResampleFlag = false; resampleCount = 0; alreadySegmentedFlag = false;};

    void AccessibilityClassSetup(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgatmnet, bool highAccuracy, voro::container_periodic_poly *,VORONOI_NETWORK *, std::vector<BASIC_VCELL> *, std::vector<VOR_CELL> *);
    
    void FindChannels(double r_pr);
    
    /* checks if the provided point is in accessible volume */
    std::pair <bool,bool> isVPointInsideAtomAndNotAccessible(Point samplingPoint);
    
    /* checks if thr provided point on a surface of atomID is accessible */
    std::pair <bool,bool> isSPointInsideAtomAndNotAccessible(Point samplingPoint, int atomID);
    
    bool isVPointAccessible(Point samplingPoint);
    
    /* Checks is provide point is inside any atom and if it is accessible */
    /* the functions returns value = true if atoms overlaps/is not accessible) */
    /* atomID is an atom in the original network from where the sampling point is
     this is to check if the point is inside atom in SA calculaton */
    std::pair <bool,bool> isPointInsideAtomAndNotAccessible(Point samplingPoint, int atomID);
    
    /* Returns the last calculated minDist (calculated in isPointInsideAtomAndNotAccessible() )  */
    double lastMinDist(){return tempMinDist;};
    /* Returns channel or pocket ID for the last calcluated point */
    std::pair <int,int> lastChannelOrPocket()
    {
        if(channelMapping[tempNodeID]<0 && pocketMapping[tempNodeID]<0){
            std::cerr << "CoP_error: cannot determine point accessiblity. Consider running at higher accuracy (-ha flag).(minDist=" << tempMinDist << ")." << std::endl;
            std::cerr << "NodeID= " << tempNodeID << "  minDistAtomID= " << tempMinDistAtomID << " resampleFlag= " << needToResampleFlag << std::endl;
            abort();
        }
        return std::pair <int,int> (channelMapping[tempNodeID],pocketMapping[tempNodeID]);
    };
    
    /* Returns information if a point need to be resampled due to node not found in isPointInsideAtomAndNotAccessible() */
    bool needToResample()
    {
/*
        if(needToResampleFlag) {
            needToResampleFlag = false; // change the value after first used of the function
            return true;
        }
        else return false;
*/

        if(needToResampleFlag==true) std::cout << "Resample flag is raised. Resample count = " << resampleCount << std::endl;

        return needToResampleFlag;
    }
    
    /* Return number of resampled ponits */
    int getResampleCount(){return resampleCount;};
    
    /* remove Voronoi nodes that are not used in analysis */
    void removeOverlappedNodes();
    
    /* perform segmentation of the Vornet */
    void segmentPoresBasedOnRadius(double seg_r);
    bool segmentPoresBasedOnFile(string filename);
    bool segmentPoresBasedOnMoleculesPresent(string);

    void calculatePLDbasedOnRadius(double seq_r);
    void calculatePLDbasedOnFile(string);
    void calculatePLDbasedOnMoleculesPresent(string);
    void calculatePLD(); // core function that performs calculation

    void reportPLD(ostream &output); // prints the calculated PLD table

    void getPLDvisData(vector<Point> *, vector<int> *, vector<double> *, string);

    /* deconstruct class */
    void deconstruct(){};
};



#endif

