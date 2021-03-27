//-----Richard Luis Martin 2012/12/15
//-----Program to read a list of hologram files and output n diverse structures, using the first structure as the 'seed'

#include "holo_similarity.h"

typedef struct {int candidate; int comparator;} tuple;

float **similarities, **tiebreak_similarities;
char **file_locations;
char *is_valid;

main(int argc, char *argv[]) {

	int i,j;

	//read command line
	int needed_args = 8;
	int optional_args = 2;
	if(argc!=needed_args && argc!=needed_args+optional_args) {
		printf("Number of command line arguments should be %d. Please try again, e.g.:\n", needed_args);
		printf("./a.out list_of_all_holos_files num_structures diversity_output_filename num_desired max_similarity DIFFUSION_RANGE SIMILARITY_MEASURE\n\t(optional args: run_index num_runs, for distributed seeding)\n");
		printf("Program will identify diverse structures until either num_desired are chosen, or max_similarity is exceeded; give either criterion a negative value to ignore it\n");
		print_similarity_descriptions();
		exit(EXIT_FAILURE);
	}

	//read file and num_structures
	FILE *in, *out, *holo;
	in = fopen(argv[1], "r");
	if(in==NULL) {
		printf("ERROR: could not open file with name %s. Please try again.\n", argv[1]);
		exit(EXIT_FAILURE);
	}
	int num_structures = atoi(argv[2]);
	int num_desired = atoi(argv[4]);
	float max_similarity = atof(argv[5]);
	printf("Program will identify diverse structures until either %d are chosen, or %.6f similarity is exceeded (if either value is <0, this condition will not be considered)\n", num_desired, max_similarity);
	const int NUM_BINS = 16; //finalised
	const int DIFFUSION_RANGE = atoi(argv[6]);
	const int SIMILARITY_MEASURE = atoi(argv[7]);
	if(SIMILARITY_MEASURE>MAX_ALLOWED || SIMILARITY_MEASURE<0) {
		printf("ERROR: invalid input for SIMILARITY_MEASURE\n");
		exit(EXIT_FAILURE);
	}
  printf("SIMILARITY_MEASURE = %s\n", similarity_coefficient_names[SIMILARITY_MEASURE-1]);

	int run_index = 0, num_runs = 0, seed = 0;
	if(argc==needed_args+optional_args) {
		run_index = atoi(argv[8]);
		num_runs = atoi(argv[9]);
    if(num_runs<=0) {
		  printf("ERROR: num_runs should be positive, but has value %d. Please try again.\n", num_runs);
		  exit(EXIT_FAILURE);
    }
    if(run_index>=num_runs || run_index<0) {
		  printf("ERROR: run_index should be in range 0 to num_runs-1, i.e. 0 to %d, but has value %d. Please try again.\n", num_runs-1, run_index);
		  exit(EXIT_FAILURE);
    }
		seed = rintf(((float)(num_structures*run_index))/((float)(num_runs)));
    printf("Running %d times with seeds evenly sampled from file list\n", num_runs);
	} else printf("Using first structure as seed\n");

	//set up specific structure stats
	hologram *hologram_array;
	hologram_array = new hologram[num_structures];
	file_locations = new char*[num_structures];
	is_valid = new char[num_structures];
	similarities = new float*[num_structures];
	//initialise arrays
	for(i=0; i<num_structures; i++) {
		file_locations[i] = new char[100];
		is_valid[i] = 0;
		hologram_array[i].num_entries = 0;
	}
	
	//loop over lines in file and update stats
  int status = 0;
	for(i=0; i<num_structures; i++) {
 		status = fscanf(in, "%s", file_locations[i]);
	  if(status==-1) {
		  printf("ERROR: could not read file name for structure index %d. Please try again.\n", i);
		  exit(EXIT_FAILURE);
	  }
    holo = fopen(file_locations[i], "r");
	  if(holo==NULL) {
		  printf("ERROR: could not open file with name %s. Please try again.\n", file_locations[i]);
		  exit(EXIT_FAILURE);
	  }
	  coords_int grid_point;
	  int quantity;
	  search_for_char(holo, '\n'); //this puts us at the end of the header line for this structure
    status = 0;
	  while(status!=-1) {
		  status = fscanf(holo, "%d %d %d %d", &grid_point.x, &grid_point.y, &grid_point.z, &quantity);
      if(status!=-1) {
		    hologram_array[i].grid_points.push_back(grid_point);
		    hologram_array[i].quantities.push_back(quantity);
		    hologram_array[i].num_entries++;
      }
	  }
		if(hologram_array[i].num_entries>0) is_valid[i] = 1;
    fclose(holo);
	}
	fclose(in);

	//need a big array for each of the two structures being compared at any moment; fill with values incorporating diffusion based on nearby entries
	unsigned int k=0,l=0,m=0;
	float ***array_i, ***array_j;
	array_i = new float**[NUM_BINS];
	array_j = new float**[NUM_BINS];
	for(k=0; k<NUM_BINS; k++) {
		array_i[k] = new float*[NUM_BINS];
		array_j[k] = new float*[NUM_BINS];
		for(l=0; l<NUM_BINS; l++) {
			array_i[k][l] = new float[NUM_BINS];
			array_j[k][l] = new float[NUM_BINS];
		}
	}
	//now we have all the holograms read in, begin the loop
	vector <tuple> indices_of_selected;
	vector <float> lowest_max_sims;
	while(!is_valid[seed]) {
		seed++;
		if(seed>=num_structures) seed = 0;
	}
	for(k=0; k<NUM_BINS; k++) {
		for(l=0; l<NUM_BINS; l++) {
			for(m=0; m<NUM_BINS; m++) {
				array_i[k][l][m] = 0;
			}
		}
	}
	build_array(array_i, hologram_array[seed], NUM_BINS, DIFFUSION_RANGE);
	//keep track of num_selected
	int num_selected = 0;
	tuple t = {seed, -1};
	indices_of_selected.push_back(t); //i.e. 'seed' the diversity selection with either the first valid structure in the list (index 0), or the distributed seed index
	lowest_max_sims.push_back(-1);
	num_selected++;
	printf("Seeding the diversity selection with structure index %d, %s\n", seed, file_locations[seed]);
  printf("Number of structures selected so far = %d\n", num_selected);
	char exceeded_max_similarity = 0;
	while((num_selected<num_desired || num_desired<0) && (exceeded_max_similarity==0 || max_similarity<0)) {
		similarities[num_selected-1] = new float[num_structures];
		for(j=0; j<num_structures; j++) {
			similarities[num_selected-1][j] = -1;
		}
		//1) calculate the similarities between the most-recently selected structure and all others (keeping track of previously calculated similarities)
		for(j=0; j<num_structures; j++) {
			if(is_valid[j]) {
				for(k=0; k<NUM_BINS; k++) {
					for(l=0; l<NUM_BINS; l++) {
						for(m=0; m<NUM_BINS; m++) {
							array_j[k][l][m] = 0;
						}
					}
				}
				build_array(array_j, hologram_array[j], NUM_BINS, DIFFUSION_RANGE);
    		float similarity = get_similarity(array_i, array_j, NUM_BINS, DIFFUSION_RANGE, SIMILARITY_MEASURE);
				similarities[num_selected-1][j] = similarity;
			} else similarities[num_selected-1][j] = -1;
		}
		//2) find the structure which has the lowest maximum similarity w.r.t. previously selected structures
		float current_lowest_max_sim = -1;
		vector<tuple> candidates;
		int num_candidates=0;
		for(j=0; j<num_structures; j++) {
			float max_sim = -1;
			int temp_comparator = -1;
			for(i=0; i<num_selected; i++) {
				if((similarities[i][j]>max_sim || max_sim<0) && similarities[i][j]>=0) {
					max_sim = similarities[i][j];
					temp_comparator = i;
				}
			}
			if((max_sim<current_lowest_max_sim || current_lowest_max_sim<0) && max_sim>=0) {
				current_lowest_max_sim = max_sim;
				num_candidates = 1;
				candidates.clear();
				tuple t = {j,temp_comparator};
				candidates.push_back(t);
			} else if(max_sim<=current_lowest_max_sim && num_candidates>0 && max_sim>=0) { //tie-break scenario, where at least one structure has already been selected as valid
				num_candidates++;
				tuple t = {j,temp_comparator};
				candidates.push_back(t);
			}
		}
		if(current_lowest_max_sim>max_similarity && !(max_similarity<0)) {
			exceeded_max_similarity = 1;
			printf("Not proceeding because max_similarity has been exceeded (%.6f>%.6f)\n", current_lowest_max_sim, max_similarity);
		} else {
			//2b) in the case of a tie, we don't modify the similarities but instead work with a different array of similarities
			int current_selection, current_comparator;
			int num_tiebreaks = 0;
			while(num_candidates>1) {
				tiebreak_similarities = new float*[num_selected];
				for(i=0; i<num_selected; i++) {
					tiebreak_similarities[i] = new float[num_structures];
					for(j=0; j<num_structures; j++) {
						tiebreak_similarities[i][j] = -1;
					}
				}
				num_tiebreaks++;
				//at this point rather than just the earliest, lowest value, we have a vector of all structures with the same lowest value (tie-break by increasing diffusion)
				if(DIFFUSION_RANGE+num_tiebreaks>NUM_BINS/2) {
					//nothing changes from here onwards - if the structures are still tied they are likely similar - it is easier for now to just run to this point and choose the alphabetically first structure
					num_candidates = 1;
					printf("%d structures all have max similarity of %.6f to any selected structure; force-ending tiebreak because increasing diffusion will have no effect\n", num_candidates, current_lowest_max_sim);
				} else {
					printf("%d structures all have max similarity of %.6f to any selected structure; tiebreak based on increasing diffusion to %d and recalculating similarities:\n", num_candidates, current_lowest_max_sim, DIFFUSION_RANGE+num_tiebreaks);

					for(i=0; i<num_selected; i++) { //for each selected structure, remake the array based in more diffusion
						for(k=0; k<NUM_BINS; k++) {
							for(l=0; l<NUM_BINS; l++) {
								for(m=0; m<NUM_BINS; m++) {
									array_i[k][l][m] = 0;
								}
							}
						}
						build_array(array_i, hologram_array[i], NUM_BINS, DIFFUSION_RANGE+num_tiebreaks);
						int n;
						for(n=0; n<num_candidates; n++) {
							int j = candidates.at(n).candidate;
							if(is_valid[j]) {
								for(k=0; k<NUM_BINS; k++) {
									for(l=0; l<NUM_BINS; l++) {
										for(m=0; m<NUM_BINS; m++) {
											array_j[k][l][m] = 0;
										}
									}
								}
								build_array(array_j, hologram_array[j], NUM_BINS, DIFFUSION_RANGE+num_tiebreaks);
            		float similarity = get_similarity(array_i, array_j, NUM_BINS, DIFFUSION_RANGE, SIMILARITY_MEASURE);
								tiebreak_similarities[i][j] = similarity; //write to [i][j] so we keep track of j
							} else tiebreak_similarities[i][j] = -1;
						}
					} //end of foreach num_selected
					//now calculate the smallest maximum similarity
					float tiebreak_current_lowest_max_sim = -1;
					vector<tuple> tiebreak_candidates;
					int tiebreak_num_candidates=0;
					int n;
					for(n=0; n<num_candidates; n++) {
						int j = candidates.at(n).candidate;
						int temp_comparator = candidates.at(n).comparator;
						float max_sim = -1;
						for(i=0; i<num_selected; i++) {
							if((tiebreak_similarities[i][j]>max_sim || max_sim<0) && tiebreak_similarities[i][j]>=0) max_sim = tiebreak_similarities[i][j];
						}
						printf("\tstructure %d (%s) has tiebreak similarity of %.6f\n", j, file_locations[j], max_sim);
						if((max_sim<tiebreak_current_lowest_max_sim || tiebreak_current_lowest_max_sim<0) && max_sim>=0) {
							tiebreak_current_lowest_max_sim = max_sim;
							tiebreak_num_candidates = 1;
							tiebreak_candidates.clear();
							tuple t = {j,temp_comparator};
							tiebreak_candidates.push_back(t);
						} else if(max_sim<=tiebreak_current_lowest_max_sim && tiebreak_num_candidates>0 && max_sim>=0) { //tie-break scenario, where at least one structure has already been selected as valid
							tiebreak_num_candidates++;
							tuple t = {j,temp_comparator};
							tiebreak_candidates.push_back(t);
						}
					}
					//now see if we need to repeat, or can select one structure now
					num_candidates = tiebreak_num_candidates;
					candidates = tiebreak_candidates;
					for(i=0; i<num_selected; i++) {
						delete tiebreak_similarities[i];
					}
				}
			} //end while(>1)
			//now this loop is over, num_candidates is 1 or lower
			if(num_candidates<=0) {
				printf("ERROR: could not select the structure with the lowest maximum similarity (probably because only invalid structures remain)\n");
				exit(EXIT_FAILURE);
			} else { //else == 1
				current_selection = candidates.at(0).candidate;
				current_comparator = candidates.at(0).comparator;
			}
			//3) select this structure
			tuple t = {current_selection, current_comparator};
			indices_of_selected.push_back(t);
			lowest_max_sims.push_back(current_lowest_max_sim);
			printf("Adding structure index %d, %s, similarity %.6f\n", current_selection, file_locations[current_selection], current_lowest_max_sim);
			num_selected++;
      printf("Number of structures selected so far = %d\n", num_selected);
			//4) update array_i accordingly
			for(k=0; k<NUM_BINS; k++) {
				for(l=0; l<NUM_BINS; l++) {
					for(m=0; m<NUM_BINS; m++) {
						array_i[k][l][m] = 0;
					}
				}
			}
			build_array(array_i, hologram_array[current_selection], NUM_BINS, DIFFUSION_RANGE);
		} //end else (haven't exceeded max_similarity)

	}
	if(num_desired!=num_selected && exceeded_max_similarity==0) {
		printf("ERROR: finished loop without exceeding %.6f similarity, having selected %d structures, where %d were desired\n", max_similarity, num_selected, num_desired);
		exit(EXIT_FAILURE);
	}
	printf("%d structures selected.\n", num_selected);

	//write output
	out = fopen(argv[3], "w");
	if(out==NULL) {
		printf("ERROR: could not open file with name %s. Please try again.\n", argv[3]);
		exit(EXIT_FAILURE);
	}
	for(i=0; i<num_selected; i++) {
		if(indices_of_selected.at(i).comparator == -1) {
			fprintf(out, "%s\n", file_locations[indices_of_selected.at(i).candidate]);
		} else {
			fprintf(out, "%s %.6f %s\n", file_locations[indices_of_selected.at(i).candidate], lowest_max_sims[i], file_locations[indices_of_selected.at(indices_of_selected.at(i).comparator).candidate]);
		}
	}
	fclose(out);

	//end program
	printf("Program complete.\n\n");
}

