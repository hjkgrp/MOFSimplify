/* routines for segmentation of Voronai channel systems */
/* added my Maciej Haranczyk, 4/20/2012                 */

#ifndef SEGMENT_H
#define SEGMENT_H

#include <iostream>
#include <vector>

#include "networkstorage.h"
#include "feature.h"
#define SCALE_SPHERE 0.7

//void segmentChannels(VORONOI_NETWORK *vornet, float probeRad, vector<bool> *accessInfo, char *name);

/*
class SEGMENT_CONNECTION{
  int seg_from, seg_to;
  int node_from, node_to;
};


class SEGMENT{
  int id;                            // SEGMENT id
  CHANNEL *parent;                   // Pointer to parent 
  map<int,int> idMappings;           // (node id, new index) pairs
  map<int,int> reverseIDMappings;    // (new index, node id) pairs
  vector<DIJKSTRA_NODE> nodes;       // List of nodes in SEGMENT indexed by new index
  vector<CONN> connections;          // List of connections between nodes contained in same SEGMENT

  SEGMENT(vector<int> nodeIDs, DIJKSTRA_NETWORK *dnet, CHANNEL *parent){
    for(int i = 0; i < nodeIDs.size(); i++){
      DIJKSTRA_NODE oldNode = dnet->nodes[nodeIDs[i]];
      idMappings.insert(pair<int, int> (oldNode.id, i));
      reverseIDMappings.insert(pair<int,int>(i, oldNode.id));
    }
   
  }

};
*/

/*
void createSegments(CHANNEL chan);
*/

/** This is a function that segments channels based on feature (bad description...) **/
void segmentChannels(ATOM_NETWORK *uc,std::vector<FEATURE>  *channel_list,std::ostream &output);

/* This function runs segmentation of channels and provides two vectors: seg_ids - with segment IDs for each node in voronoi network
 *   and segmaxr with max radius assigned to each segment (if featflag set to 0) */
/* the same function is used to do this for features (featflag=1) */
void segmentChannels_forHolograms(ATOM_NETWORK *uc,std::vector<FEATURE>  *channel_list,std::vector <int> *seg_ids,std::vector <double> *segmaxr,int featflag) ;

/* This function save active segment files to run simluations for particular segments */
void segmentChannels_saveSegments(ATOM_NETWORK *uc, VORONOI_NETWORK *vornet,std::vector<int> *accessInfo, std::vector<double> *segment_radii_vector,int num_segments, char *name,const char *feattype);

#endif
