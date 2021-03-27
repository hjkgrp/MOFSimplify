//-----Richard Luis Martin 2012/12/15
//-----This code generates so-called Voronoi holograms, which are 3D vector
//-----representations of the accessible void space within a material.
//-----With respect to a certain probe size, the accessible voronoi network
//-----is calculated and each edge of the network is characterized in terms of
//-----length, start and end radii. The hologram stores the count of each edge
//-----type with respect to these three properties.
//-----The code outputs: 1) "*.stats", which describes the contents of the hologram; 2) "*_holo.txt" which contains only the non-zero entries; 3) "*_complete_holo.txt" which contains the complete holograms written on one line.

//#include "network.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "holograms.h"

using namespace std;


/*
void analyze_accessible_voronoi_by_atoms(ATOM_NETWORK *atmnet, VORONOI_NETWORK *vornet, float probeRad, vector<int> *accessInfo, vector<double> *segment_radii_vector, int num_segments, char *name) { //accessInfo now contains for each node, segment number (1+) or 0 if inaccessible

	cout << "NUM_BINS = " << NUM_BINS << "\n";

	unsigned int i,j,k,l,m;

	//make forced assignment by IS_COMBINED_SEGMENTS
	vector<int> *backup_accessInfo;
	int backup_num_segments;
	if(IS_COMBINED_SEGMENTS==1) {
		if(num_segments>1) {
			backup_accessInfo = new vector<int>;
			backup_num_segments = num_segments;
			for(i=0; i<accessInfo->size(); i++) {
				backup_accessInfo->push_back(accessInfo->at(i));
				if(accessInfo->at(i)>1) accessInfo->at(i) = 1;
			}
			num_segments = 1;
		}
	}

	//make a list of atoms which correspond to each segment
	char **atom_in_this_segment; //segments*atoms
	char ***atom_in_this_segment_per_node; //segments*nodes*atoms
	int *num_atoms_in_this_segment, *num_atoms_in_this_segment_per_node, *centre_node_access_index, *num_nodes_in_this_segment;
	atom_in_this_segment = new char*[num_segments];
	atom_in_this_segment_per_node = new char**[num_segments];
	num_atoms_in_this_segment = new int[num_segments];
	num_atoms_in_this_segment_per_node = new int[num_segments];
	centre_node_access_index = new int[num_segments];
	num_nodes_in_this_segment = new int[num_segments];
	unsigned int node_count = 0;

	vector<VOR_NODE> accessible_voronoi_nodes = vector<VOR_NODE>();
	for(j=0; j<num_segments; j++) { //work on each segment in order
		atom_in_this_segment[j] = new char[atmnet->numAtoms];
		atom_in_this_segment_per_node[j] = new char*[accessInfo->size()];
		for(i=0; i<atmnet->numAtoms; i++) atom_in_this_segment[j][i]=0; //easiest way to keep track of which atoms are in this segment
		num_atoms_in_this_segment[j] = 0;
		num_atoms_in_this_segment_per_node[j] = 0;
		num_nodes_in_this_segment[j] = 0;
		float max_node_radius = -1; //used to find which node is the largest in this segment - need this as the reference point for atom distance calculation
		for(i=0; i<accessInfo->size(); i++) {
			atom_in_this_segment_per_node[j][i] = new char[atmnet->numAtoms];
			for(k=0; k<atmnet->numAtoms; k++) atom_in_this_segment_per_node[j][i][k] = 0;
			if(accessInfo->at(i)==j+1) { //for each node in this segment...
				if(vornet->nodes.at(i).rad_stat_sphere>max_node_radius || max_node_radius<0) {
					max_node_radius = vornet->nodes.at(i).rad_stat_sphere;
					centre_node_access_index[j] = i;
				}
				for(k=0; k<atmnet->numAtoms; k++) { //...find the corresponding atoms
					if(calcDistance(vornet->nodes.at(i).x, vornet->nodes.at(i).y, vornet->nodes.at(i).z, &(atmnet->atoms.at(k)), atmnet)<vornet->nodes.at(i).rad_stat_sphere+ATOM_GRAB_DISTANCE) {
						atom_in_this_segment_per_node[j][i][k]=1;
						num_atoms_in_this_segment_per_node[j]++;
						if(atom_in_this_segment[j][k]==0) {
							atom_in_this_segment[j][k]=1;
							num_atoms_in_this_segment[j]++;
						}
					}
				}
				num_nodes_in_this_segment[j]++;
				node_count++;
			}
		}
		cout << "Segment index " << j << " has " << num_atoms_in_this_segment[j] << " corresponding atoms, or " << num_atoms_in_this_segment_per_node[j] << " pairs of atoms-nodes.\n";
		cout << "Segment index " << j << " has max_node_radius " << max_node_radius << "\n";
	}

	//now bin each pair of atoms in a segment by their edge lengths following a pre-set binning system
	int int_edge_length_bins[NUM_BINS];
	float float_edge_length_bins[NUM_BINS];
	float *float_edge_length_bin_upper_bounds;
	float_edge_length_bin_upper_bounds = new float[NUM_BINS];
	//a 1D binning system unlike for original (3D) holograms
	int int_edge_stats_bins[NUM_BINS];
	float float_edge_stats_bins[NUM_BINS];

	//build bin data 'manually' for this task
	for(i=0; i<NUM_BINS-1; i++) {
		if(IS_INTER_ATOMIC_DISTANCE==0)
			float_edge_length_bin_upper_bounds[i] = 3.0 + 0.25*((float)i);
		else
			float_edge_length_bin_upper_bounds[i] = 3.0 + 0.5*((float)i); //idea for atom triplets (pharmacophores)
	}	
	float_edge_length_bin_upper_bounds[NUM_BINS-1] = 1000; //i.e. 'infinity'

	//DEBUG print bins
	printf("'edge' length bins:\n");
	for(i=0; i<NUM_BINS; i++) {
		printf("\t%.3f\n", float_edge_length_bin_upper_bounds[i]);
	}

	//instantiate arrays
	for(i=0; i<NUM_BINS; i++) {
		int_edge_stats_bins[i] = 0;
		float_edge_stats_bins[i] = 0;
	}

	//prepare to write to seg_atom_holo output file
	char outputFile2 [256];
	strcpy(outputFile2,name);
	strcat(outputFile2,".seg_atom_holo");
	fstream output2; output2.open(outputFile2, fstream::out);
	output2 << name << "\n";
	output2 << num_segments << "\n";

	//prepare to write complete_atom_holo
	char outputFile4 [256];
	strcpy(outputFile4,name);
	strcat(outputFile4,".complete_atom_holo");
	fstream output4; output4.open(outputFile4, fstream::out);

	//only do anything if there are any nodes!
	if(node_count>0) {
	//for VisIt
	float x_offset = -1*(min(0.0,atmnet->v_b.x) + min(0.0,atmnet->v_c.x));
	float y_offset = -1*min(0.0,atmnet->v_c.y);
	printf("DEBUG: x_offset = %.3f, y_offset = %.3f (for VisIt)\n", x_offset, y_offset);
		//now handle the 3D edge-based descriptors on a segment by segment basis
		for(l=0; l<(unsigned int)num_segments; l++) {
			XYZ *best_coords, *best_coords_nodes;
			best_coords = new XYZ[num_atoms_in_this_segment[l]];
			best_coords_nodes = new XYZ[num_nodes_in_this_segment[l]];
			int count_temp = 0;
			for(i=0; i<accessInfo->size(); i++) {
				if(accessInfo->at(i)==l+1) {
					best_coords_nodes[count_temp] = getClosestPointInABC(vornet->nodes.at(centre_node_access_index[l]).x, vornet->nodes.at(centre_node_access_index[l]).y, vornet->nodes.at(centre_node_access_index[l]).z, vornet->nodes.at(i).x, vornet->nodes.at(i).y, vornet->nodes.at(i).z, atmnet);
					count_temp++;
				}
			}
			//prepare to print out the atoms for VisIt
			vector<XYZ> final_atom_coords;
			//seg_atom_holo
			output2 << l << "\n";
			float radius = segment_radii_vector->at(l);
			output2 << radius << "\n";
			//complete_atom_holo
			output4 << name << "_" << l; //to avoid problem with columns
			//new atom-based hologram generation: find the distance from each corresponding atom to the centre point
			int count = 0;
int num_cage_atoms_before_checking_images = 0;
			for(i=0; i<atmnet->numAtoms; i++) {
				//addition to this function - keep track of the positions pushed
				vector<XYZ> pushed_images;
				if(atom_in_this_segment[l][i]==1) {
					float distance = calcDistance(vornet->nodes.at(centre_node_access_index[l]).x, vornet->nodes.at(centre_node_access_index[l]).y, vornet->nodes.at(centre_node_access_index[l]).z, &(atmnet->atoms.at(i)), atmnet);
					best_coords[count] = getClosestPointInABC(vornet->nodes.at(centre_node_access_index[l]).x, vornet->nodes.at(centre_node_access_index[l]).y, vornet->nodes.at(centre_node_access_index[l]).z, atmnet->atoms.at(i).x, atmnet->atoms.at(i).y, atmnet->atoms.at(i).z, atmnet);
					Point best_xyz = abc_to_xyz(best_coords[count].x, best_coords[count].y, best_coords[count].z, atmnet);
					XYZ base_atom;
					base_atom.x = best_xyz.x;
					base_atom.y = best_xyz.y;
					base_atom.z = best_xyz.z;
					final_atom_coords.push_back(base_atom);
					num_cage_atoms_before_checking_images++;
					if(IS_INTER_ATOMIC_DISTANCE==0) {
						int bin = get_bin(distance, float_edge_length_bin_upper_bounds);
						if(IS_COUNT==0) int_edge_stats_bins[bin]+=100;
						else int_edge_stats_bins[bin]++;
					}
					//addition to this function, check for other valid positions
					pushed_images.push_back(best_coords[count]); //in abc
					int i_diff, j_diff, k_diff;
					for(j=0; j<accessInfo->size(); j++) {
						//for each node, if related to this atom, get distance and check images for the right distance
						if(atom_in_this_segment_per_node[l][j][i]==1) {
							//firstly move this node to be in the same region as the central node!
							int node_find = 0;
							for(k=0; k<j; k++) {
								if(accessInfo->at(k)==l+1) node_find++; //basically the number of segment nodes before this one is the index
							}
							Point new_node_xyz = abc_to_xyz(best_coords_nodes[node_find].x, best_coords_nodes[node_find].y, best_coords_nodes[node_find].z, atmnet);
							float pushed_distance = sqrt( ((best_xyz.x-new_node_xyz.x)*(best_xyz.x-new_node_xyz.x)) + ((best_xyz.y-new_node_xyz.y)*(best_xyz.y-new_node_xyz.y)) + ((best_xyz.z-new_node_xyz.z)*(best_xyz.z-new_node_xyz.z)) );
							for(i_diff=-1; i_diff<=1; i_diff++) {
								for(j_diff=-1; j_diff<=1; j_diff++) {
									for(k_diff=-1; k_diff<=1; k_diff++) {
										if(!(i_diff==0 && j_diff==0 && k_diff==0)) {
											//so we are looking at periodic images of this 'best' abc
											XYZ candidate_image;
											candidate_image.x = best_coords[count].x+i_diff;
											candidate_image.y = best_coords[count].y+j_diff;
											candidate_image.z = best_coords[count].z+k_diff;
											char already_pushed = 0;
											for(k=0; k<pushed_images.size() && already_pushed==0; k++) {
												if(candidate_image.x==pushed_images.at(k).x && candidate_image.y==pushed_images.at(k).y && candidate_image.z==pushed_images.at(k).z) already_pushed = 1; 
											}
											if(already_pushed==0) {
												Point image_xyz = abc_to_xyz(candidate_image.x, candidate_image.y, candidate_image.z, atmnet);
												float image_distance = sqrt( ((image_xyz.x-new_node_xyz.x)*(image_xyz.x-new_node_xyz.x)) + ((image_xyz.y-new_node_xyz.y)*(image_xyz.y-new_node_xyz.y)) + ((image_xyz.z-new_node_xyz.z)*(image_xyz.z-new_node_xyz.z)) );
												if(image_distance<vornet->nodes.at(j).rad_stat_sphere+ATOM_GRAB_DISTANCE) {
													pushed_images.push_back(candidate_image);
													XYZ image_atom;
													image_atom.x = image_xyz.x;
													image_atom.y = image_xyz.y;
													image_atom.z = image_xyz.z;
													final_atom_coords.push_back(image_atom);
												}
											}
										}
									}
								}
							}
						}
					}
					count++;
				}

			}

			//output atoms for VisIt
			char output_visit_atoms [256], int_string[10];
			sprintf(int_string, "%d", l);
			strcpy(output_visit_atoms,name);
			strcat(output_visit_atoms,"_");
			strcat(output_visit_atoms,int_string);
			strcat(output_visit_atoms,".xyz");
			fstream output3; output3.open(output_visit_atoms, fstream::out);
			output3 << final_atom_coords.size() << "\nCartesian coords of segment's representative atoms\n";
			char output_tensor_atoms [256];
			strcpy(output_tensor_atoms,name);
			strcat(output_tensor_atoms,"_tensor_");
			strcat(output_tensor_atoms,int_string);
			strcat(output_tensor_atoms,".xyz");
			fstream output5; output5.open(output_tensor_atoms, fstream::out);
			output5 << final_atom_coords.size() << "\n";
			for(i=0; i<final_atom_coords.size(); i++) {
				if(i<num_cage_atoms_before_checking_images) output3 << "N " << final_atom_coords.at(i).x + x_offset << " " << final_atom_coords.at(i).y + y_offset << " " << final_atom_coords.at(i).z << " " << 0.1 << "\n";
				else output3 << "O " << final_atom_coords.at(i).x + x_offset << " " << final_atom_coords.at(i).y + y_offset << " " << final_atom_coords.at(i).z << " " << 0.1 << "\n";
				output5 << "1\t" << final_atom_coords.at(i).x + x_offset << "\t" << final_atom_coords.at(i).y + y_offset << "\t" << final_atom_coords.at(i).z << "\n";
			}
			output3.close();
			output5.close();

			//now calculate pharmacophores by looping over the xyz coords in final_atom_coords
			printf("DEBUG: pharmacophore data for segment %d\n", l);
			int pharmacophores[NUM_BINS][NUM_BINS][NUM_BINS];
			for(i=0; i<NUM_BINS; i++) {
				for(j=0; j<NUM_BINS; j++) {
					for(k=0; k<NUM_BINS; k++) {
						pharmacophores[i][j][k] = 0;
					}
				}
			}
			float shortest = -1, longest = -1;
			int num_minimal_equilaterals = 0, num_right_angles = 0;
			int num_obtuse = 0, num_acute = 0, num_right = 0, num_over_135 = 0;
			int num_triangles = 0;
			float shortest_max = -1, longest_min = -1;
			for(i=0; i<final_atom_coords.size(); i++) {
				for(j=i+1; j<final_atom_coords.size(); j++) {
					for(k=j+1; k<final_atom_coords.size(); k++) {
						//note the loop system ensures that i<j<k
						float i_j_dist = sqrt( ((final_atom_coords.at(i).x-final_atom_coords.at(j).x)*(final_atom_coords.at(i).x-final_atom_coords.at(j).x)) + ((final_atom_coords.at(i).y-final_atom_coords.at(j).y)*(final_atom_coords.at(i).y-final_atom_coords.at(j).y)) + ((final_atom_coords.at(i).z-final_atom_coords.at(j).z)*(final_atom_coords.at(i).z-final_atom_coords.at(j).z)) );
						float i_k_dist = sqrt( ((final_atom_coords.at(i).x-final_atom_coords.at(k).x)*(final_atom_coords.at(i).x-final_atom_coords.at(k).x)) + ((final_atom_coords.at(i).y-final_atom_coords.at(k).y)*(final_atom_coords.at(i).y-final_atom_coords.at(k).y)) + ((final_atom_coords.at(i).z-final_atom_coords.at(k).z)*(final_atom_coords.at(i).z-final_atom_coords.at(k).z)) );
						float min_dist = min(i_j_dist, i_k_dist);
						float max_dist = max(i_j_dist, i_k_dist);
						float j_k_dist = sqrt( ((final_atom_coords.at(j).x-final_atom_coords.at(k).x)*(final_atom_coords.at(j).x-final_atom_coords.at(k).x)) + ((final_atom_coords.at(j).y-final_atom_coords.at(k).y)*(final_atom_coords.at(j).y-final_atom_coords.at(k).y)) + ((final_atom_coords.at(j).z-final_atom_coords.at(k).z)*(final_atom_coords.at(j).z-final_atom_coords.at(k).z)) );
						//now this third dist is either shorter than both, longer than both, or in between - assume in between then modify if needed
						float mid_dist = j_k_dist;
						if(j_k_dist>max_dist) {
							mid_dist = max_dist;
							max_dist = j_k_dist;
						} else if(j_k_dist<min_dist) {
							mid_dist = min_dist;
							min_dist = j_k_dist;
						}
						//NOTE: it is possible, given these three distances, to calculate the three angles in the triangle - we could use this information as a shape descriptor
						//NOTE: also an angle near 180 degrees, i.e. min+mid is slightly greater than max, might not be interesting, just like a minimal equilateral triangle might be very common and so not that interesting
						//printf("\t%.3f %.3f %.3f\n", min_dist, mid_dist, max_dist);
						if(min_dist<shortest || shortest<0) shortest = min_dist;
						if(max_dist>longest || longest<0) longest = max_dist;
						if(max_dist<shortest_max || shortest_max<0) shortest_max = max_dist;
						if(min_dist>longest_min || longest_min<0) longest_min = min_dist;
						//debug: get stats
						if(max_dist-2.7<ROUNDING_ERROR) num_minimal_equilaterals++;
						if(abs((min_dist*min_dist)+(mid_dist*mid_dist)-(max_dist*max_dist))<ROUNDING_ERROR) num_right_angles++;
						num_triangles++;
						//now get the angle opposite the longest side
						float cos_gamma = ((max_dist*max_dist)-((min_dist*min_dist)+(mid_dist*mid_dist)))/(-2.0*min_dist*mid_dist);
						float gamma = (acos(cos_gamma)*360.0)/(2.0*PI);
						if(gamma>135+ROUNDING_ERROR) num_over_135++;
						if(gamma>90+ROUNDING_ERROR) num_obtuse++;
						else if(gamma<90-ROUNDING_ERROR) num_acute++;
						else num_right++;
						//now bin them
						int bin_1 = get_bin(min_dist, float_edge_length_bin_upper_bounds);
						int bin_2 = get_bin(mid_dist, float_edge_length_bin_upper_bounds);
						int bin_3 = get_bin(max_dist, float_edge_length_bin_upper_bounds);
						pharmacophores[bin_1][bin_2][bin_3]++;
					}
				}
			}

			int iter = 0;
			printf("DEBUG: pharmacophore non-zero entries:\n");
			int iter2 = 0;
			for(i=0; i<NUM_BINS; i++) {
				for(j=i+1; j<NUM_BINS; j++) {
					for(k=i+1; k<NUM_BINS; k++) {
						if(pharmacophores[i][j][k] > 0) {
							output2 << "# " << i << " " << j << " " << k << " " << pharmacophores[i][j][k] << "\n"; //seg_atom_holo
							iter++;
						}
						output4 << " " << pharmacophores[i][j][k]; //complete_atom_holo
						iter2++;
					}
				}
			}
			printf("\t%d pharmacophore non-zero entries (out of %d) for segment %d\n", iter, iter2, l);
			output4 << "\n";

		} //end for each segment l
	} //end if any nodes

	//finish writing to seg_atom_holo file
	output2 << "\n";
	output2.close();
	//complete_atom_holo
	output4.close();

	//undo forced assignment by IS_COMBINED_SEGMENTS
	if(IS_COMBINED_SEGMENTS==1) {
		if(backup_num_segments>1) {
			for(i=0; i<backup_accessInfo->size(); i++) {
				accessInfo->at(i) = backup_accessInfo->at(i);
			}
			num_segments = backup_num_segments;
		}
	}
}

void analyze_accessible_voronoi_with_segments(VORONOI_NETWORK *vornet, float probeRad, vector<int> *accessInfo, vector<double> *segment_radii_vector, int num_segments, char *name) { //accessInfo now contains for each node, segment number (1+) or 0 if inaccessible
	
	cout << "NUM_BINS = " << NUM_BINS << "\n";

	//find info and count accessible nodes
	unsigned int i,j,k,l;
	unsigned int node_count = 0, edge_count = 0;
	float min_node_radius = -1, max_node_radius = -1, average_node_radius = 0, min_edge_radius = -1, max_edge_radius = -1, average_edge_radius = 0, min_edge_length = -1, max_edge_length = -1, average_edge_length = 0;
	vector<VOR_NODE> accessible_voronoi_nodes = vector<VOR_NODE>();
	for(i=0; i<accessInfo->size(); i++) {
		if(accessInfo->at(i)>0) {
			average_node_radius+=vornet->nodes.at(i).rad_stat_sphere;
			accessible_voronoi_nodes.push_back(vornet->nodes.at(i));
			if(vornet->nodes.at(i).rad_stat_sphere<min_node_radius || min_node_radius<0) min_node_radius = vornet->nodes.at(i).rad_stat_sphere;
			if(vornet->nodes.at(i).rad_stat_sphere>max_node_radius || max_node_radius<0) max_node_radius = vornet->nodes.at(i).rad_stat_sphere;
			node_count ++;
		}
	}
	cout << node_count << " accessible nodes in total." << "\n";
	if(node_count>0) {
		average_node_radius/=node_count;
	}

	//based on accessible nodes, determine accessible edges
	vector<VOR_EDGE> accessible_voronoi_edges = vector<VOR_EDGE>();
	vector<unsigned int> accessible_voronoi_edge_indices = vector<unsigned int>();
	vector<int> accessible_voronoi_edge_segment_indices = vector<int>();
	vector <int> segment_num_edges(num_segments);
	for(i=0; i<(unsigned int)num_segments; i++) {
		segment_num_edges[i] = 0;
	}
	for(i=0; i<vornet->edges.size(); i++) {
		if(accessInfo->at(vornet->edges.at(i).from)>0 && accessInfo->at(vornet->edges.at(i).to)>0 && vornet->edges.at(i).rad_moving_sphere>probeRad) { //i.e. if both nodes for this edge are accessible, and if the radius is larger than the probe radius
			average_edge_radius+=vornet->edges.at(i).rad_moving_sphere;
			average_edge_length+=vornet->edges.at(i).length;
			if(vornet->edges.at(i).rad_moving_sphere<min_edge_radius || min_edge_radius<0) min_edge_radius = vornet->edges.at(i).rad_moving_sphere;
			if(vornet->edges.at(i).rad_moving_sphere>max_edge_radius || max_edge_radius<0) max_edge_radius = vornet->edges.at(i).rad_moving_sphere;
			if(vornet->edges.at(i).length<min_edge_length || min_edge_length<0) min_edge_length = vornet->edges.at(i).length;
			if(vornet->edges.at(i).length>max_edge_length || max_edge_length<0) max_edge_length = vornet->edges.at(i).length;
			accessible_voronoi_edges.push_back(vornet->edges.at(i));
			accessible_voronoi_edge_indices.push_back(i);
			if(accessInfo->at(vornet->edges.at(i).from) == accessInfo->at(vornet->edges.at(i).to)) { //if to and from have same segment
				int segment_index = accessInfo->at(vornet->edges.at(i).from)-1; //convert to index and push back - this edge exists within a segment
				accessible_voronoi_edge_segment_indices.push_back(segment_index);
				segment_num_edges[segment_index]++;
			} else accessible_voronoi_edge_segment_indices.push_back(-1); //push back invalid number - this edge connects two segments
			edge_count ++;
		}
	}
	cout << edge_count << " accessible edges in total." << "\n";
	if(edge_count>0) {
		average_edge_radius/=edge_count;
		average_edge_length/=edge_count;
	}

	//write stats to terminal
	cout << "min_node_radius = " << min_node_radius << " max_node_radius = " << max_node_radius << " average_node_radius = " << average_node_radius << "\n";
//	cout << "min_edge_radius = " << min_edge_radius << " max_edge_radius = " << max_edge_radius << " average_edge_radius = " << average_edge_radius << "\n";
	cout << "min_edge_length = " << min_edge_length << " max_edge_length = " << max_edge_length << " average_edge_length = " << average_edge_length << "\n";

	//now bin the nodes and edges following a pre-set binning system
	int int_node_radii_bins[NUM_BINS], int_edge_radii_bins[NUM_BINS], int_edge_length_bins[NUM_BINS];
	float float_node_radii_bins[NUM_BINS], float_edge_radii_bins[NUM_BINS], float_edge_length_bins[NUM_BINS];
	float *float_node_radii_bin_upper_bounds, *float_edge_radii_bin_upper_bounds, *float_edge_length_bin_upper_bounds;
	float_node_radii_bin_upper_bounds = new float[NUM_BINS], float_edge_radii_bin_upper_bounds = new float[NUM_BINS], float_edge_length_bin_upper_bounds = new float[NUM_BINS];
	//a 3D binning system considering for each edge its 1) length, 2/3) radii of largest then smallest attached node
	int int_edge_stats_bins[NUM_BINS][NUM_BINS][NUM_BINS];
	float float_edge_stats_bins[NUM_BINS][NUM_BINS][NUM_BINS];

	//BIN READING
	//this process does one of two things:
	//1) no bin_directory is specified - the default bins are used
	//2) a bin_directory is given - grab the files from there
  if(!bin_directory) { //case 1)
		printf("Bin directory not specified - using defaults.\n");
		for(i=0; i<NUM_BINS-1; i++) {
			float_edge_length_bin_upper_bounds[i] = default_edge_length_bins[i];
			float_node_radii_bin_upper_bounds[i] = default_node_radii_bins[i];
    }
	} else { //case 2)
  	FILE *edge_length_bin_file, *node_radii_bin_file;
    char *edge_length_bin_string = new char[100];
    strcpy(edge_length_bin_string, bin_directory);
    strcat(edge_length_bin_string, "edge_length_bins");
    char *node_radii_bin_string = new char[100];
    strcpy(node_radii_bin_string, bin_directory);
    strcat(node_radii_bin_string, "node_radii_bins");
		printf("Bin directory specified - trying to open files %s and %s.\n", edge_length_bin_string, node_radii_bin_string);
	  edge_length_bin_file = fopen(edge_length_bin_string, "r");
	  node_radii_bin_file = fopen(node_radii_bin_string, "r");
    char error = 0;
    if(edge_length_bin_file==NULL) {
      printf("ERROR: could not open bins file %s\n", edge_length_bin_string);
      error = 1;
    }
    if(node_radii_bin_file==NULL) {
      printf("ERROR: could not open bins file %s\n", node_radii_bin_string);
      error = 1;
    }
    if(error==1) exit(EXIT_FAILURE);
		printf("Files for bin bounds opened successfully.\n");
		//read in the values for each bin from the files
		for(i=0; i<NUM_BINS-1; i++) {
			float temp1, temp2;
			int status;
			status = fscanf(edge_length_bin_file, "%f", &temp1);
			float_edge_length_bin_upper_bounds[i] = temp1;
      if(status==-1) {
        printf("ERROR: could not read edge length bin bound number %d from file %s\n", i+1, edge_length_bin_string);
        exit(EXIT_FAILURE);
      }
			status = fscanf(node_radii_bin_file, "%f", &temp2);
			float_node_radii_bin_upper_bounds[i] = temp2;
      if(status==-1) {
        printf("ERROR: could not read node radii bin bound number %d from file %s\n", i+1, node_radii_bin_string);
        exit(EXIT_FAILURE);
      }
		}
		fclose(edge_length_bin_file);
		fclose(node_radii_bin_file);
		//we've read the bins in so can make the holograms
	}
  //now check that the bins make sense
  int error = 0;
  for(i=1; i<NUM_BINS-1 && error==0; i++) {
    if(float_edge_length_bin_upper_bounds[i]<float_edge_length_bin_upper_bounds[i-1]) error = 1;
    if(error==1) {
      printf("ERROR: edge length bins do not constitute a non-decreasing series (entry index %d = %f; entry index %d = %f)\n", i-1, float_edge_length_bin_upper_bounds[i-1], i, float_edge_length_bin_upper_bounds[i]);
      if(!bin_directory) printf("\t\tPlease contact the developers since this problem occurred for the default bins\n");
      exit(EXIT_FAILURE);
    }
    if(float_node_radii_bin_upper_bounds[i]<float_node_radii_bin_upper_bounds[i-1]) error = 1;
    if(error==1) {
      printf("ERROR: node radii bins do not constitute a non-decreasing series (entry index %d = %f; entry index %d = %f)\n", i-1, float_node_radii_bin_upper_bounds[i-1], i, float_node_radii_bin_upper_bounds[i]);
      if(!bin_directory) printf("\t\tPlease contact the developers since this problem occurred for the default bins\n");
      exit(EXIT_FAILURE);
    }
  }
  //now safe to move on
	float_node_radii_bin_upper_bounds[NUM_BINS-1] = 1000; //i.e. 'infinity'
	float_edge_radii_bin_upper_bounds[NUM_BINS-1] = 1000; //i.e. 'infinity'
	float_edge_length_bin_upper_bounds[NUM_BINS-1] = 1000; //i.e. 'infinity'

	//instantiate arrays
	for(i=0; i<NUM_BINS; i++) {
		int_node_radii_bins[i] = 0;
		int_edge_radii_bins[i] = 0;
		int_edge_length_bins[i] = 0;
		float_node_radii_bins[i] = 0;
		float_edge_radii_bins[i] = 0;
		float_edge_length_bins[i] = 0;
		for(j=0; j<NUM_BINS; j++) {
			for(k=0; k<NUM_BINS; k++) {
				int_edge_stats_bins[i][j][k] = 0;
				float_edge_stats_bins[i][j][k] = 0;
			}
		}
	}

	//prepare to write to seg_holo output file
	char outputFile2 [256];
	strcpy(outputFile2,name);
	strcat(outputFile2,".seg_holo");
	fstream output2; output2.open(outputFile2, fstream::out);
	output2 << name << "\n";
	output2 << num_segments << "\n";

	//only do anything if there are any nodes!
	if(node_count>0) {
		//do 1D binning for general stats
		for(i=0; i<accessible_voronoi_nodes.size(); i++) {
			float radius = accessible_voronoi_nodes.at(i).rad_stat_sphere;
			int bin_index = -1;
			bin_index = get_bin(radius, float_node_radii_bin_upper_bounds);
			if(IS_COUNT==0) int_node_radii_bins[bin_index]+=100; //100 so we can get percentages
			else int_node_radii_bins[bin_index]++;
		}
		for(i=0; i<accessible_voronoi_edges.size(); i++) {
			float radius = accessible_voronoi_edges.at(i).rad_moving_sphere, length = accessible_voronoi_edges.at(i).length;
			int bin_index = -1;
			bin_index = get_bin(radius, float_edge_radii_bin_upper_bounds);
			if(IS_COUNT==0) int_edge_radii_bins[bin_index]+=100; //100 so we can get percentages
			else int_edge_radii_bins[bin_index]++;
			bin_index = -1;
			bin_index = get_bin(length, float_edge_length_bin_upper_bounds);
			if(IS_COUNT==0) int_edge_length_bins[bin_index]+=100; //100 so we can get percentages
			else int_edge_length_bins[bin_index]++;
		}
		//now convert binning sums into binning percentages
		for(i=0; i<NUM_BINS; i++) {
			float_node_radii_bins[i] = ((float)int_node_radii_bins[i])/accessible_voronoi_nodes.size();
			float_edge_radii_bins[i] = ((float)int_edge_radii_bins[i])/accessible_voronoi_edges.size();
			float_edge_length_bins[i] = ((float)int_edge_length_bins[i])/accessible_voronoi_edges.size();
		}
		//now handle the 3D edge-based descriptors on a segment by segment basis
		for(l=0; l<(unsigned int)num_segments; l++) {
			output2 << l << "\n";
			float radius = segment_radii_vector->at(l);
			output2 << radius << "\n";
			for(i=0; i<accessible_voronoi_edges.size(); i++) {
				if(accessible_voronoi_edge_segment_indices.at(i)==(int)l) { //check that this edge is in the segment being considered
					float from_radii = vornet->nodes.at(vornet->edges.at(accessible_voronoi_edge_indices.at(i)).from).rad_stat_sphere;
					float to_radii = vornet->nodes.at(vornet->edges.at(accessible_voronoi_edge_indices.at(i)).to).rad_stat_sphere;
					float min_radii = min(from_radii,to_radii), max_radii = max(from_radii,to_radii);
					float length = accessible_voronoi_edges.at(i).length;
					if(IS_COUNT==0) int_edge_stats_bins[get_bin(length, float_edge_length_bin_upper_bounds)][get_bin(max_radii, float_node_radii_bin_upper_bounds)][get_bin(min_radii, float_node_radii_bin_upper_bounds)]+=100; //100 so we can get percentages
					else int_edge_stats_bins[get_bin(length, float_edge_length_bin_upper_bounds)][get_bin(max_radii, float_node_radii_bin_upper_bounds)][get_bin(min_radii, float_node_radii_bin_upper_bounds)]++; //100 so we can get percentages
				}
			}
			//now convert binning sums into binning percentages
			for(i=0; i<NUM_BINS; i++) {
				for(j=0; j<NUM_BINS; j++) {
					for(k=0; k<=j; k++) { //because k is never larger than j
						float_edge_stats_bins[i][j][k] = ((float)int_edge_stats_bins[i][j][k])/segment_num_edges[l];
						if(float_edge_stats_bins[i][j][k] > 0) {
							//write this point to seg_holo file
							if(IS_COUNT==0) output2 << "# " << i << " " << j << " " << k << " " << float_edge_stats_bins[i][j][k] << "\n";
							else output2 << "# " << i << " " << j << " " << k << " " << int_edge_stats_bins[i][j][k] << "\n";
						}
					}
				}
				//now reset the segment-specific counts
				for(j=0; j<NUM_BINS; j++) {
					for(k=0; k<NUM_BINS; k++) {
						int_edge_stats_bins[i][j][k] = 0;
						float_edge_stats_bins[i][j][k] = 0;
					}
				}
			}
		}
	} //end if any nodes

	//finish writing to seg_holo file
	output2 << "\n";
	output2.close();

	//write stats to terminal - at this stage if no nodes, all percentages will be zero
	cout << "edge_length bin percentages:";
	for(i=0; i<NUM_BINS; i++) {
		if(IS_COUNT==0) cout << " " << float_edge_length_bins[i];
		else cout << " " << int_edge_length_bins[i];
	}
	cout << "\n";
	cout << "node_radius bin percentages:";
	for(i=0; i<NUM_BINS; i++) {
		if(IS_COUNT==0) cout << " " << float_node_radii_bins[i];
		else cout << " " << int_node_radii_bins[i];
	}
	cout << "\n";
//	cout << "edge_radius bin percentages:";
//	for(i=0; i<NUM_BINS; i++) {
//		if(IS_COUNT==0) cout << " " << float_edge_radii_bins[i];
//		else cout << " " << int_edge_radii_bins[i];
//	}
//	cout << "\n";

	//now write to stats output file
	char outputFile [256];
	strcpy(outputFile,name);
	strcat(outputFile,".stats");
	fstream output; output.open(outputFile, fstream::out);
	output << name << ": " << node_count << " " << min_node_radius << " " << max_node_radius << " " << average_node_radius << " " << edge_count << " " << min_edge_radius << " " << max_edge_radius << " " << average_edge_radius << " " << min_edge_length << " " << max_edge_length << " " << average_edge_length << " ";
	for(i=0; i<NUM_BINS; i++) {
		if(IS_COUNT==0) output << " " << float_node_radii_bins[i];
		else output << " " << int_node_radii_bins[i];
	}
	output << " ";
	for(i=0; i<NUM_BINS; i++) {
		if(IS_COUNT==0) output << " " << float_edge_radii_bins[i];
		else output << " " << int_edge_radii_bins[i];
	}
	output << " ";
	for(i=0; i<NUM_BINS; i++) {
		if(IS_COUNT==0) output << " " << float_edge_length_bins[i];
		else output << " " << int_edge_length_bins[i];
	}
	output << "\n";
	output.close();
}
*/

//the main code for analyzing the accessible network and categorizing the edges, producing hologram representations
void analyze_accessible_voronoi_pre_segment(VORONOI_NETWORK *vornet, float probeRad, vector<bool> *accessInfo, char *name, const char *bin_directory) {

//1) inspect the accessibility of the network - specifically, we need to find the accessible edges so they can be categorized

	cout << "NUM_BINS = " << NUM_BINS << "\n";

	//find info and count accessible nodes
	unsigned int i,j,k;
	unsigned int node_count = 0, edge_count = 0;
	float min_node_radius = -1, max_node_radius = -1, average_node_radius = 0, min_edge_radius = -1, max_edge_radius = -1, average_edge_radius = 0, min_edge_length = -1, max_edge_length = -1, average_edge_length = 0;
	vector<VOR_NODE> accessible_voronoi_nodes = vector<VOR_NODE>();
	for(i=0; i<accessInfo->size(); i++) {
		if(accessInfo->at(i) || IS_FULL_VORONOI==1) {
			average_node_radius+=vornet->nodes.at(i).rad_stat_sphere;
			accessible_voronoi_nodes.push_back(vornet->nodes.at(i));
			if(vornet->nodes.at(i).rad_stat_sphere<min_node_radius || min_node_radius<0) min_node_radius = vornet->nodes.at(i).rad_stat_sphere;
			if(vornet->nodes.at(i).rad_stat_sphere>max_node_radius || max_node_radius<0) max_node_radius = vornet->nodes.at(i).rad_stat_sphere;
			node_count ++;
		}
	}
	cout << node_count << " accessible nodes in total." << "\n";
	if(node_count>0) {
		average_node_radius/=node_count;
	}
	//based on accessible nodes, determine accessible edges
	vector<VOR_EDGE> accessible_voronoi_edges = vector<VOR_EDGE>();
	vector<unsigned int> accessible_voronoi_edge_indices = vector<unsigned int>();
	for(i=0; i<vornet->edges.size(); i++) {
		if((accessInfo->at(vornet->edges.at(i).from) && accessInfo->at(vornet->edges.at(i).to) && vornet->edges.at(i).rad_moving_sphere>probeRad) || IS_FULL_VORONOI==1) { //i.e. if both nodes for this edge are accessible, and if the radius is larger than the probe radius
			average_edge_radius+=vornet->edges.at(i).rad_moving_sphere;
			average_edge_length+=vornet->edges.at(i).length;
			if(vornet->edges.at(i).rad_moving_sphere<min_edge_radius || min_edge_radius<0) min_edge_radius = vornet->edges.at(i).rad_moving_sphere;
			if(vornet->edges.at(i).rad_moving_sphere>max_edge_radius || max_edge_radius<0) max_edge_radius = vornet->edges.at(i).rad_moving_sphere;
			if(vornet->edges.at(i).length<min_edge_length || min_edge_length<0) min_edge_length = vornet->edges.at(i).length;
			if(vornet->edges.at(i).length>max_edge_length || max_edge_length<0) max_edge_length = vornet->edges.at(i).length;
			accessible_voronoi_edges.push_back(vornet->edges.at(i));
			accessible_voronoi_edge_indices.push_back(i);
			edge_count ++;
		}
	}
	cout << edge_count << " accessible edges in total." << "\n";
	if(edge_count>0) {
		average_edge_radius/=edge_count;
		average_edge_length/=edge_count;
	}

	//write stats to terminal
	cout << "min_node_radius = " << min_node_radius << " max_node_radius = " << max_node_radius << " average_node_radius = " << average_node_radius << "\n";
//	cout << "min_edge_radius = " << min_edge_radius << " max_edge_radius = " << max_edge_radius << " average_edge_radius = " << average_edge_radius << "\n";
	cout << "min_edge_length = " << min_edge_length << " max_edge_length = " << max_edge_length << " average_edge_length = " << average_edge_length << "\n";

//2) now we have the accessible edges, categorize them in terms of length, start and end radii, using three separate histograms

	//now bin the nodes and edges following a pre-set binning system
	int int_node_radii_bins[NUM_BINS], int_edge_radii_bins[NUM_BINS], int_edge_length_bins[NUM_BINS];
	float float_node_radii_bins[NUM_BINS], float_edge_radii_bins[NUM_BINS], float_edge_length_bins[NUM_BINS];
	//geometric progression custom binning for edge lengths
	float *float_node_radii_bin_upper_bounds, *float_edge_radii_bin_upper_bounds, *float_edge_length_bin_upper_bounds;
	float_node_radii_bin_upper_bounds = new float[NUM_BINS], float_edge_radii_bin_upper_bounds = new float[NUM_BINS], float_edge_length_bin_upper_bounds = new float[NUM_BINS];
	//a 3D binning system considering for each edge its 1) length, 2/3) radii of largest/smallest attached node
	int int_edge_stats_bins[NUM_BINS][NUM_BINS][NUM_BINS];
	float float_edge_stats_bins[NUM_BINS][NUM_BINS][NUM_BINS];

	//BIN READING
	//this process does one of two things:
	//1) no bin_directory is specified - the default bins are used
	//2) a bin_directory is given - grab the files from there
  if(!bin_directory) { //case 1)
		printf("Bin directory not specified - using defaults.\n");
		for(i=0; i<NUM_BINS-1; i++) {
			float_edge_length_bin_upper_bounds[i] = default_edge_length_bins[i];
			float_node_radii_bin_upper_bounds[i] = default_node_radii_bins[i];
    }
	} else { //case 2)
  	FILE *edge_length_bin_file, *node_radii_bin_file;
    char *edge_length_bin_string = new char[100];
    strcpy(edge_length_bin_string, bin_directory);
    strcat(edge_length_bin_string, "edge_length_bins");
    char *node_radii_bin_string = new char[100];
    strcpy(node_radii_bin_string, bin_directory);
    strcat(node_radii_bin_string, "node_radii_bins");
		printf("Bin directory specified - trying to open files %s and %s.\n", edge_length_bin_string, node_radii_bin_string);
	  edge_length_bin_file = fopen(edge_length_bin_string, "r");
	  node_radii_bin_file = fopen(node_radii_bin_string, "r");
    char error = 0;
    if(edge_length_bin_file==NULL) {
      printf("ERROR: could not open bins file %s\n", edge_length_bin_string);
      error = 1;
    }
    if(node_radii_bin_file==NULL) {
      printf("ERROR: could not open bins file %s\n", node_radii_bin_string);
      error = 1;
    }
    if(error==1) exit(EXIT_FAILURE);
		printf("Files for bin bounds opened successfully.\n");
		//read in the values for each bin from the files
		for(i=0; i<NUM_BINS-1; i++) {
			float temp1, temp2;
			int status;
			status = fscanf(edge_length_bin_file, "%f", &temp1);
			float_edge_length_bin_upper_bounds[i] = temp1;
      if(status==-1) {
        printf("ERROR: could not read edge length bin bound number %d from file %s\n", i+1, edge_length_bin_string);
        exit(EXIT_FAILURE);
      }
			status = fscanf(node_radii_bin_file, "%f", &temp2);
			float_node_radii_bin_upper_bounds[i] = temp2;
      if(status==-1) {
        printf("ERROR: could not read node radii bin bound number %d from file %s\n", i+1, node_radii_bin_string);
        exit(EXIT_FAILURE);
      }
		}
		fclose(edge_length_bin_file);
		fclose(node_radii_bin_file);
		//we've read the bins in so can make the holograms
	}
  //now check that the bins make sense
  int error = 0;
  for(i=1; i<NUM_BINS-1 && error==0; i++) {
    if(float_edge_length_bin_upper_bounds[i]<float_edge_length_bin_upper_bounds[i-1]) error = 1;
    if(error==1) {
      printf("ERROR: edge length bins do not constitute a non-decreasing series (entry index %d = %f; entry index %d = %f)\n", i-1, float_edge_length_bin_upper_bounds[i-1], i, float_edge_length_bin_upper_bounds[i]);
      if(!bin_directory) printf("\t\tPlease contact the developers since this problem occurred for the default bins\n");
      exit(EXIT_FAILURE);
    }
    if(float_node_radii_bin_upper_bounds[i]<float_node_radii_bin_upper_bounds[i-1]) error = 1;
    if(error==1) {
      printf("ERROR: node radii bins do not constitute a non-decreasing series (entry index %d = %f; entry index %d = %f)\n", i-1, float_node_radii_bin_upper_bounds[i-1], i, float_node_radii_bin_upper_bounds[i]);
      if(!bin_directory) printf("\t\tPlease contact the developers since this problem occurred for the default bins\n");
      exit(EXIT_FAILURE);
    }
  }

//3) prepare the 3D vectors that will store the hologram data, and fill them
  //now safe to move on
	float_node_radii_bin_upper_bounds[NUM_BINS-1] = 1000; //i.e. 'infinity'
	float_edge_length_bin_upper_bounds[NUM_BINS-1] = 1000; //i.e. 'infinity'
	for(i=0; i<NUM_BINS; i++) {
		int_node_radii_bins[i] = 0;
		int_edge_length_bins[i] = 0;
		float_node_radii_bins[i] = 0;
		float_edge_length_bins[i] = 0;
		for(j=0; j<NUM_BINS; j++) {
			for(k=0; k<NUM_BINS; k++) {
				int_edge_stats_bins[i][j][k] = 0;
				float_edge_stats_bins[i][j][k] = 0;
			}
		}
	}
	//only do anything if there are any nodes!
	if(node_count>0) {
		for(i=0; i<accessible_voronoi_nodes.size(); i++) {
			float radius = accessible_voronoi_nodes.at(i).rad_stat_sphere;
			int bin_index = -1;
			bin_index = get_bin(radius, float_node_radii_bin_upper_bounds);
			if(IS_COUNT==1) int_node_radii_bins[bin_index]++; //1 for count
			else int_node_radii_bins[bin_index]+=100; //100 so we can get percentages
		}
		for(i=0; i<accessible_voronoi_edges.size(); i++) {
			float radius = accessible_voronoi_edges.at(i).rad_moving_sphere, length = accessible_voronoi_edges.at(i).length;
			int bin_index = -1;
			bin_index = get_bin(radius, float_edge_radii_bin_upper_bounds);
			if(IS_COUNT==1) int_edge_radii_bins[bin_index]++; //1 for count
			else int_edge_radii_bins[bin_index]+=100; //100 so we can get percentages
			bin_index = -1;
			bin_index = get_bin(length, float_edge_length_bin_upper_bounds);
			if(IS_COUNT==1) int_edge_length_bins[bin_index]++; //1 for count
			else int_edge_length_bins[bin_index]+=100; //100 so we can get percentages
			//now handle the 3D edge-based descriptors
			float from_radii = vornet->nodes.at(vornet->edges.at(accessible_voronoi_edge_indices.at(i)).from).rad_stat_sphere;
			float to_radii = vornet->nodes.at(vornet->edges.at(accessible_voronoi_edge_indices.at(i)).to).rad_stat_sphere;
			float min_radii = min(from_radii,to_radii), max_radii = max(from_radii,to_radii);

			if(IS_COUNT==1) int_edge_stats_bins[get_bin(length, float_edge_length_bin_upper_bounds)][get_bin(max_radii, float_node_radii_bin_upper_bounds)][get_bin(min_radii, float_node_radii_bin_upper_bounds)]++; //1 for count
			else int_edge_stats_bins[get_bin(length, float_edge_length_bin_upper_bounds)][get_bin(max_radii, float_node_radii_bin_upper_bounds)][get_bin(min_radii, float_node_radii_bin_upper_bounds)]+=100; //100 so we can get percentages
		}
		//now convert binning sums into binning proportions
		for(i=0; i<NUM_BINS; i++) {
			float_node_radii_bins[i] = ((float)int_node_radii_bins[i])/accessible_voronoi_nodes.size();
			float_edge_radii_bins[i] = ((float)int_edge_radii_bins[i])/accessible_voronoi_edges.size();
			float_edge_length_bins[i] = ((float)int_edge_length_bins[i])/accessible_voronoi_edges.size();
			for(j=0; j<NUM_BINS; j++) {
				for(k=0; k<=j; k++) { //because k is never larger than j
					float_edge_stats_bins[i][j][k] = ((float)int_edge_stats_bins[i][j][k])/accessible_voronoi_edges.size();
				}
			}
		}
	} //end if any nodes

	//write stats to terminal - at this stage if no nodes, all percentages will be zero
	if(IS_COUNT==1) {
		cout << "edge_length bin counts:";
		for(i=0; i<NUM_BINS; i++) {
		        cout << " " << int_edge_length_bins[i];
		}
		cout << "\n";
		cout << "node_radius bin counts:";
		for(i=0; i<NUM_BINS; i++) {
		        cout << " " << int_node_radii_bins[i];
		}
		cout << "\n";
//		cout << "edge_radius bin counts:";
//		for(i=0; i<NUM_BINS; i++) {
//		        cout << " " << int_edge_radii_bins[i];
//		}
//		cout << "\n";
	} else {
		cout << "node_radius bin percentages:";
		for(i=0; i<NUM_BINS; i++) {
			cout << " " << float_node_radii_bins[i];
		}
		cout << "\n";
//		cout << "edge_radius bin percentages:";
//		for(i=0; i<NUM_BINS; i++) {
//			cout << " " << float_edge_radii_bins[i];
//		}
//		cout << "\n";
		cout << "edge_length bin percentages:";
		for(i=0; i<NUM_BINS; i++) {
			cout << " " << float_edge_length_bins[i];
		}
		cout << "\n";
	}

//4) file output section; writes three files: "*.stats", which describes the contents of the hologram; "*.holo.txt", which is a concise representation
//storing only the non-zero entries, and also in a format which can be visualized in the VisIt software package; "*complete_holo.txt", which is a one-line
//print out of the complete 3D hologram, including zero-valued entries.

	//now write to output file
	char outputFile [256];
	strcpy(outputFile,name);
	strcat(outputFile,".stats");
  FILE *output;
	output = fopen(outputFile, "w");
  if(output==NULL) {
    printf("ERROR: cannot open file %s to write Voronoi hologram statistics\n", outputFile);
    exit(EXIT_FAILURE);
  }
	fprintf(output, "%s statistics:\n %d Nodes:\n  min_radius %.6f\n  max_radius %.6f\n  average_radius %.6f\n %d Edges:\n  min_radius %.6f\n  max_radius %.6f\n  average_radius %.6f\n  min_length %.6f\n  max_length %.6f\n  average_length %.6f\n", name, node_count, min_node_radius, max_node_radius, average_node_radius, edge_count, min_edge_radius, max_edge_radius, average_edge_radius, min_edge_length, max_edge_length, average_edge_length);
	if(IS_COUNT==1) {
    fprintf(output, " Node_radii_bin_values:\n ");
		for(i=0; i<NUM_BINS; i++) {
		  fprintf(output, " %d", int_node_radii_bins[i]);
		}
    fprintf(output, "\n");
		fprintf(output, " Edge_length_bin_values:\n ");
		for(i=0; i<NUM_BINS; i++) {
		  fprintf(output, " %d", int_edge_length_bins[i]);
		}
    fprintf(output, "\n");
		fclose(output);
	} else {
    fprintf(output, " Node_radii_bin_values:\n ");
		for(i=0; i<NUM_BINS; i++) {
		  fprintf(output, " %.3f", float_node_radii_bins[i]);
		}
    fprintf(output, "\n");
		fprintf(output, " Edge_length_bin_values:\n ");
		for(i=0; i<NUM_BINS; i++) {
		  fprintf(output, " %.3f", float_edge_length_bins[i]);
		}
    fprintf(output, "\n");
		fclose(output);
	}
/*
	fstream output; output.open(outputFile, fstream::out);
	output << name << ": " << node_count << " " << min_node_radius << " " << max_node_radius << " " << average_node_radius << " " << edge_count << " " << min_edge_radius << " " << max_edge_radius << " " << average_edge_radius << " " << min_edge_length << " " << max_edge_length << " " << average_edge_length << " ";
	if(IS_COUNT==1) {
		for(i=0; i<NUM_BINS; i++) {
		        output << " " << int_node_radii_bins[i];
		}
		output << " ";
		for(i=0; i<NUM_BINS; i++) {
		        output << " " << int_edge_radii_bins[i];
		}
		output << " ";
		for(i=0; i<NUM_BINS; i++) {
		        output << " " << int_edge_length_bins[i];
		}
		output << "\n";
		output.close();
	} else {
		for(i=0; i<NUM_BINS; i++) {
			output << " " << float_node_radii_bins[i];
		}
		output << " ";
		for(i=0; i<NUM_BINS; i++) {
			output << " " << float_edge_radii_bins[i];
		}
		output << " ";
		for(i=0; i<NUM_BINS; i++) {
			output << " " << float_edge_length_bins[i];
		}
		output << "\n";
		output.close();
	}
*/

	//now write to second output file!
	char outputFile2 [256];
	strcpy(outputFile2,name);
	strcat(outputFile2,"_holo.txt");
/*  int num_bins_set = 0;
	for(i=0; i<NUM_BINS; i++) {
		for(j=0; j<NUM_BINS; j++) {
			for(k=0; k<=j; k++) { //because k is never larger than j
				if(int_edge_stats_bins[i][j][k] > 0) num_bins_set++;
      }
    }
  }*/
  FILE *output2;
	output2 = fopen(outputFile2, "w");
  if(output2==NULL) {
    printf("ERROR: cannot open file %s to write Voronoi hologram xyz format file\n", outputFile2);
    exit(EXIT_FAILURE);
  }
//	fstream output2; output2.open(outputFile2, fstream::out);
//	output2 << num_bins_set << "\n" << name << " Voronoi hologram as xyz format\n";
	fprintf(output2, "x y z frequency\n");
	if(IS_COUNT==1) {
		for(i=0; i<NUM_BINS; i++) {
			for(j=0; j<NUM_BINS; j++) {
				for(k=0; k<=j; k++) { //because k is never larger than j
					if(int_edge_stats_bins[i][j][k] > 0) {
//            output2 << "# " << i << " " << j << " " << k << " " << int_edge_stats_bins[i][j][k] << "\n";
            fprintf(output2, "%d %d %d %d\n", i, j, k, int_edge_stats_bins[i][j][k]);
					}
				}
			}
		}
		fprintf(output2, "\n");
		fclose(output2);
	} else {
		for(i=0; i<NUM_BINS; i++) {
			for(j=0; j<NUM_BINS; j++) {
				for(k=0; k<=j; k++) { //because k is never larger than j
					if(int_edge_stats_bins[i][j][k] > 0) {
//						output2 << "# " << i << " " << j << " " << k << " " << float_edge_stats_bins[i][j][k] << "\n";
            fprintf(output2, "%d %d %d %.3f\n", i, j, k, float_edge_stats_bins[i][j][k]);
					}
				}
			}
		}
		fprintf(output2, "\n");
		fclose(output2);
	}

	//now write to third output file!
	char outputFile3 [256];
	strcpy(outputFile3,name);
	strcat(outputFile3,"_complete_holo.txt");
  FILE *output3;
	output3 = fopen(outputFile3, "w");
  if(output3==NULL) {
    printf("ERROR: cannot open file %s to write Voronoi hologram complete data file\n", outputFile3);
    exit(EXIT_FAILURE);
  }
//	fstream output3; output3.open(outputFile3, fstream::out);
	fprintf(output3, "%s:", name);
	if(IS_COUNT==1) {
		for(i=0; i<NUM_BINS; i++) {
			for(j=0; j<NUM_BINS; j++) {
				for(k=0; k<=j; k++) { //because k is never larger than j
					fprintf(output3, " %d", int_edge_stats_bins[i][j][k]);
				}
			}
		}
		fprintf(output3, "\n");
		fclose(output3);
	} else {
		for(i=0; i<NUM_BINS; i++) {
			for(j=0; j<NUM_BINS; j++) {
				for(k=0; k<=j; k++) { //because k is never larger than j
					fprintf(output3, " %.3f", float_edge_stats_bins[i][j][k]);
				}
			}
		}
		fprintf(output3, "\n");
		fclose(output3);
	}

}

//for a given float value, determine which bin it falls into by inspecting the upper bounds on each bin
int get_bin(float measure, float *upper_bounds) {
	int i;
	for(i=0; i<NUM_BINS-1; i++) {
		if(measure<upper_bounds[i]) {
			return i;
		}
	}
	return NUM_BINS-1; //i.e. the final bin will contain everything too large for the others
}

