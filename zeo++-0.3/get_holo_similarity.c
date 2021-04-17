//-----Richard Luis Martin 2012/12/15
//-----Program to simply read two _holo.txt files and output their similarity

#include "holo_similarity.h"

char *is_valid;

main(int argc, char *argv[]) {

	int i,j;

	//read command line
	int needed_args = 5;
	if(argc!=needed_args) {
		printf("Number of command line arguments should be %d. Please try again, e.g.:\n", needed_args);
		printf("./a.out 1_holo.txt 2_holo.txt DIFFUSION_RANGE SIMILARITY_MEASURE\n");
		print_similarity_descriptions();
		exit(EXIT_FAILURE);
	}

	//read file and num_structures
	FILE *in1, *in2;
	in1 = fopen(argv[1], "r");
	if(in1==NULL) {
		printf("ERROR: could not open file with name %s. Please try again.\n", argv[1]);
		exit(EXIT_FAILURE);
	}
	in2 = fopen(argv[2], "r");
	if(in2==NULL) {
		printf("ERROR: could not open file with name %s. Please try again.\n", argv[2]);
		exit(EXIT_FAILURE);
	}
	const int NUM_BINS = 16; //finalised value
	const int DIFFUSION_RANGE = atoi(argv[3]);
	const int SIMILARITY_MEASURE = atoi(argv[4]);
	if(SIMILARITY_MEASURE>MAX_ALLOWED || SIMILARITY_MEASURE<1) {
		printf("ERROR: invalid input for SIMILARITY_MEASURE\n");
		exit(EXIT_FAILURE);
	}
  printf("SIMILARITY_MEASURE = %s\n", similarity_coefficient_names[SIMILARITY_MEASURE-1]);

	//set up specific structure stats
	int num_structures = 2;
	hologram hologram_array[num_structures];
	is_valid = new char[num_structures];
	//initialise arrays
	for(i=0; i<num_structures; i++) {
		is_valid[i] = 0;
		hologram_array[i].num_entries = 0;
	}

	coords_int grid_point;
	int quantity;
  int status = 0;
	int num_holos_i = 1, num_holos_j = 1; //each structure has one holo only
	search_for_char(in1, '\n'); //this puts us at the end of the header line for this structure
  status = 0;
	while(status!=-1) {
		status = fscanf(in1, "%d %d %d %d", &grid_point.x, &grid_point.y, &grid_point.z, &quantity);
    if(status!=-1) {
		  hologram_array[0].grid_points.push_back(grid_point);
		  hologram_array[0].quantities.push_back(quantity);
		  hologram_array[0].num_entries++;
    }
	}
	printf("Read in file location %s\n", argv[1]);
	printf("\tthis structure has %d total hologram points\n", hologram_array[0].num_entries);
	if(hologram_array[0].num_entries>0) is_valid[0] = 1;
	fclose(in1);
	search_for_char(in2, '\n'); //this puts us at the end of the header line for this structure
	status = 0;
	while(status!=-1) {
		status = fscanf(in2, "%d %d %d %d", &grid_point.x, &grid_point.y, &grid_point.z, &quantity);
    if(status!=-1) {
	    hologram_array[1].grid_points.push_back(grid_point);
	    hologram_array[1].quantities.push_back(quantity);
	    hologram_array[1].num_entries++;
    }
	}
	printf("Read in file location %s\n", argv[2]);
	printf("\tthis structure has %d total hologram points\n", hologram_array[1].num_entries);
	if(hologram_array[1].num_entries>0) is_valid[1] = 1;
	fclose(in2);

	//now only proceed if both structures are valid
	if(is_valid[0] && is_valid[1]) {
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
		build_array(array_i, hologram_array[0], NUM_BINS, DIFFUSION_RANGE);
		build_array(array_j, hologram_array[1], NUM_BINS, DIFFUSION_RANGE);
		float similarity = get_similarity(array_i, array_j, NUM_BINS, DIFFUSION_RANGE, SIMILARITY_MEASURE);
		printf("Similarity: %.6f\n", similarity);
	} else printf("At least one structure is invalid\n");
}

