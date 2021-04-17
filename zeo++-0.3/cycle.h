/* Cycles for each vertex in voronoi graph are computed
 * Useful in material science where the voronoi nodes
 * are taken as tetrahedral interstitial sites. 
 * An octahedral interstitial site is computed as centroid of vertices
 * present in cycle of length 4 is needed.
 */
#include <vector>

#include "networkstorage.h"
#include "graphstorage.h"
#include "geometry.h"

class CYCLE{
public:
    double length;
    std::vector<DIJKSTRA_NODE> nodes;
    CYCLE();
};

/* Function to compute girth of voronoi network
 * Args:
 *  Input: 
 *      inp_vor: Pointer to input voronoi network
 *  Output:
 *      Returns length of smallest cycle
 */
double girth(const VORONOI_NETWORK* inp_vor);

/* Function to compute girth of voronoi network with integer weights
 * Args:
 *  Input: 
 *      inp_vor: Pointer to input voronoi network
 *      weight_flag: Denotes wether the weights have value other than 1
 *      range: Max. value of weights. weights = {1,2,...,range}
 *  Output:
 *      Returns length of smallest cycle
 */
//int girth(const VORONOI_NETWORK* inp_vor, bool weight_flag=false, int range=1);

/* Function to compute a cycle of given length for each voronoi node
 * Args:
 *  Input: 
 *      vornet: Pointer to input voronoi network
 *      cycle_length: Length of cycle to be computed
 *                    If 0, minimum length cycle is computed
 *      cycles: Pointer to cycle array
 *      weight_flag: Denotes wether the weights have value other than 1
 *      range: Max. value of weights. weights = {1,2,...,range}
 *  Output:
 *      Returns false if the cycle of given length cannot be found
 */
//bool compute_cycle(const VORONOI_NETWORK* vornet, int cycle_legnth, vector<CYCLE>* cycles, bool weight_flag=false, int range=1); 

/* Function to compute cycles of length 4 for each voronoi node
 * Args:
 *  Input: 
 *      vornet: Pointer to input voronoi network
 *      cycle_length: Length of cycle to be computed
 *                    If 0, minimum length cycle is computed
 *                    Use default(0) now
 *      cycles: Pointer to cycle array
 *      weight_flag: Denotes wether the weights have value other than 1
 *                  For future
 *      weight_range: Max. value of weights. weights = {1,2,...,range}
 *                  If weight_flag=false, weight_range = 1
 *                  For future
 *  Output:
 *      Returns false if the cycle of given length cannot be found
 */
bool compute_4cycle(VORONOI_NETWORK* vornet, vector<CYCLE>* cycles, bool weight_flag=false, int weigth_range=1); 

/* Computes the centroid of the dijskstra nodes in the cycle */
void centroid(const CYCLE* const, XYZ*, vector<int>*);

/* Computes the centroid of the nodes in the VORO_FACE */
void face_center(ATOM_NETWORK*, vector<XYZ>*);


