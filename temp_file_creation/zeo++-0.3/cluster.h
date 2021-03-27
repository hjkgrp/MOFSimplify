/* Clusters of voronoi nodes are analyzed and the compacted to a single 
 * voronoi node for each cluster. 
 * Useful in material science where the high symmetry voronoi nodes
 * are taken as interstitial sites
 */
#include <vector>

#include "networkstorage.h"

class Vector_XYZ{
public:
    std::vector<XYZ> nodes;
    Vector_XYZ(){}
};

/* Function to compact voronoi network
 * Args:
 *  Input: 
 *      inp_vor: Pointer to input voronoi network
 *      cluster_rad: Maximum Radius of cluster
 *                   Default is 0.5 Angstrom
 *  Output:
 *      Returns compacted voronoi network
 */
VORONOI_NETWORK cluster_reduce(const VORONOI_NETWORK* inp_vor, const float cluster_rad=0.5);

void simplify_ha_vornet(ATOM_NETWORK*);
/*
class CLUSTER {
public:
    std::vector<DIJ
    */

void high_accuracy_vornodes_reduction(ATOM_NETWORK*, Vector_XYZ*);

void high_accuracy_vornodes_reduction(ATOM_NETWORK*, std::vector<XYZ>*);

/* Function to prune high accuracy voronoi network.
 * Removes the voronoi nodes within the higher radius atoms
 * Args:
 *  Input: 
 *      ha_vor: Pointer to high accuracy voronoi network
 *      low_atm_net: Pointer to original atom network
 *      ha_atm_net: Pointer to high accuracy atom network
 *  Output:
 *      Returns compacted voronoi network
 */
void prune_high_accuracy_voronoi_network(VORONOI_NETWORK* ha_vor, ATOM_NETWORK* low_atm_net, 
                                        ATOM_NETWORK* high_atm_net, double delta, bool print=false);

/* Function to identify the nearest high accuracy voronoi node that has a high radius.
 * Removes the voronoi nodes within the higher radius atoms
 * Args:
 *  Input: 
 *      ha_vornet: Pointer to high accuracy voronoi network
 *      vornet: Pointer to original voronoi network
 *      atm_net: Pointer to atom network
 *      cutoff: Radius of cutoff sphere from each voronoi ode
 *  Output:
 *      red_vornet: Pointer to reduced voronoi network
 */
void nearest_largest_diameter_ha_vornet(VORONOI_NETWORK* ha_vornet, 
        VORONOI_NETWORK* vornet, ATOM_NETWORK* atmnet, 
        VORONOI_NETWORK* red_vornet, float cutoff=0.25);

/* Function to simplify the high accuracy voronoi network such that 
 * the small spheres pertaining to vornodes originate from different
 * atoms in the original atom network.  
 * Args:
 *  Input: 
 *      ha_vornet: Pointer to high accuracy voronoi network
 *      ha_atment: Pointer to high accuracy atom network
 *  Output:
 *      red_vornet: Pointer to reduced voronoi network
 */
void simplify_high_accuracy_vornet(VORONOI_NETWORK* ha_vornet, 
        ATOM_NETWORK* ha_atmnet, VORONOI_NETWORK* red_vornet);

/* Function to prune the voronoi network based on geometry
 * such that within a 0.1Ang^2 grid only one voronoi node is retained. 
 * Implemented that such that no two nodes are with 0.1 Ang distance 
 * Args:
 *  Input: 
 *      vornet: Pointer to input voronoi network
 *      atment: Pointer to atom network
 *      dist: Pruning threshold (float)
 *  Output:
 *      red_vornet: Pointer to reduced voronoi network
 */
void geometry_pruning(VORONOI_NETWORK* ha_vornet, ATOM_NETWORK* atmnet, 
                      float cutoff, VORONOI_NETWORK* red_vornet);

/* Function to prune the high accuracy voronoi network such that nodes within 
 * the original atoms are pruned. Nodes within 0.1Ang^2 from surface are 
 * retained. 
 * Args:
 *  Input: 
 *      ha_vornet: Pointer to high accuracy voronoi network
 *      atment: Pointer to original atom network
 *  Output:
 *      red_vornet: Pointer to reduced voronoi network
 */
void ha_prune_within_atom(VORONOI_NETWORK* ha_vornet, ATOM_NETWORK* atmnet, 
                      float cutoff, VORONOI_NETWORK* red_vornet);
