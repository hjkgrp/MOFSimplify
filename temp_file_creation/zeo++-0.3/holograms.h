// short header file for holograms.cc defining the default variables, bin upper bound values, and the functions used

#ifndef HOLOGRAMS_H
#define HOLOGRAMS_H

#include <vector>

#include "networkstorage.h"
//for all
#define NUM_BINS 16 //NOTE: we now just define NUM_BINS 16
#define IS_COUNT 1 //use this to output feature counts instead of feature proportions
/*
//for by_atoms
#define IS_INTER_ATOMIC_DISTANCE 1 //measure distance between representative atoms, rather than from some arbitrary centre voronoi node
#define IS_COMBINED_SEGMENTS 1 //treat each segment as the same - allowing construction of atom cages for entire accessible structure
#define ROUNDING_ERROR 0.001
#define ATOM_GRAB_DISTANCE 0.01
*/
//for pre_segment
#define IS_FULL_VORONOI 0 //use this to look at all nodes and edges regardless of their accessibility by their own size or by being cut off by a thin channel
int get_bin(float measure, float *upper_bounds);

//CONSTANTS
//defaults come from profiling known and random hypothetical zeolites - please see JCIM paper
const float default_edge_length_bins[16] = {0.059978, 0.138895, 0.210915, 0.308332, 0.406551, 0.517595, 0.616633, 0.732605, 0.859930, 1.017810, 1.196410, 1.436500, 1.757160, 2.231210, 3.080860, 1000.0};
const float default_node_radii_bins[16] = {1.72506, 1.82695, 1.9322, 2.0294, 2.10543, 2.19328, 2.30978, 2.42251, 2.51981, 2.64958, 2.78254, 2.97481, 3.17345, 3.40024, 3.83535, 1000.0};

void analyze_accessible_voronoi_pre_segment(VORONOI_NETWORK *vornet, float probeRad, std::vector<bool> *accessInfo, char *name, const char *bin_directory = NULL);

#endif
