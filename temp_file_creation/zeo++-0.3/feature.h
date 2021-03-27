#ifndef FEATURE_H
#define FEATURE_H

#include <iostream>
#include <fstream>
#include <vector>

#include "networkstorage.h"
#include "graphstorage.h"
#include "channel.h"
// was 1.2
#define NODE_R_EDGE_R_RATIO 1.1

// for segment merging
// May 20, was 0.95
#define FEAT_NODE_R_EDGE_R_RATIO 0.8

#define N_MC_SAMPLE 100000

#define N_MC_SA_SAMPLE 500 // num sampling points for each sphere
#define NBINS_SA_HOLO 1
#define SA_HOLO_STEP 1.50    // size of bins in SA hologram

// scaling on node radius for segmentation
#define SPHERE_SCALE 0.7

const unsigned int randSeed2 = 994879221;

class FEATURE : public CHANNEL {
  ATOM_NETWORK   *UnitCellPtr;     // pointer to the original unit cell, to be used for distance calcs with periodic boundary conditions 
  int nsegments;                   // Number of segments
  std::vector <int> segments;           // segment ID for each node
//moved_to_public  vector <double> segment_max_r;   // largest radius in a segment
  std::vector <SEGCONN> segment_conn; // List of connections between segments using floodfill approach
                                 // (this is usually smaller than the true number of connections)
  std::vector <SEGCONN> segment_conn_fulllist; // lists all connections (edges) briging different segments
                                          // this will be used to determine topology of segments
  std::vector < std::vector<int> > segment_faces_conn;   // list of connections of 'faces' for each segment, faces are clusters of nodes that connect to another segment
  std::vector < std::vector<double> > segment_faces_max_r; // max r for each segment's face  
  int nfeatures;                   // Number of features after merging
  std::vector <int> features;           // feature ID for each segment after merging
  std::vector <int> feature_node_mapping; // feature id for each node
// move to public  vector <double> feature_max_r;   // largest radius in a feature after merging
  std::vector <SEGCONN> feature_conn; // List of all connections between features after merging
  std::vector <int> segment_merge_list; // for each segment list features it is assigned to
  std::vector <double> feature_volume; // volume of each feature

public:
  std::vector <double> segment_max_r;   // largest radius in a segment
  std::vector <double> feature_max_r;   // largest radius in a feature after merging

  FEATURE(std::vector<int> nodeIDs, DIJKSTRA_NETWORK *dnet, int dim, int basisVecs[3][3]) : CHANNEL(nodeIDs, dnet, dim, basisVecs){
  }

  void setUnitCellPointer(ATOM_NETWORK *p){
    UnitCellPtr = p;
  };

    int createFeatures(ATOM_NETWORK *uc,VORONOI_NETWORK *origVNet, DIJKSTRA_NETWORK *origDNet, std::fstream &output, int initIndex,char *name);


    int createSegments(ATOM_NETWORK *uc,VORONOI_NETWORK *origVNet, DIJKSTRA_NETWORK *origDNet, std::fstream &output, int initIndex);

  void segmentChannel(ATOM_NETWORK *);
  int segment_findMaxNode();
  void segment_growSegment_start(ATOM_NETWORK *,int);
  void segment_growSegment_cont(ATOM_NETWORK *,int,int);
  void merge_segments();
  void segment_saveVis(std::fstream &);

  void segment_distBasedSegmentation(ATOM_NETWORK *);
  int validateSegmentation();
  int validateSegment(int);

  void segment_identify_connections();
  void segment_identify_connections_recurr(std::vector <int> *,std::vector <int> *,int *,int ,int , std::vector <SEGCONN> *,int);

  void merge_newSegmentGrow(int);
  int merge_findNotAssignedFragment();
  void merge_nothing();
  void merge_singlenode_segments();
  int number_nodes_in_segment(int);

  void update_node_segmentinfo(std::vector <int> *,int *,int);

  double calculateVolume();
  void calculateVolumeForAllFeatures(VORONOI_NETWORK *,DIJKSTRA_NETWORK *);

  double calculateSurfaceArea(ATOM_NETWORK *,std::vector<DIJKSTRA_NODE> *,std::vector <int> *,int,char *,int);
  double calculateSurfaceArea_OLD(ATOM_NETWORK *,std::vector<DIJKSTRA_NODE> *,std::vector <int> *,int,char *,int);
};

#endif
