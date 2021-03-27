//-----Richard Luis Martin 2012/12/15
//-----Header file for regularly-used hologram similarity functions

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
//using namespace std;

typedef struct {int x; int y; int z;} coords_int; //standard triplet of ints
typedef struct {int num_entries; std::vector <coords_int> grid_points; std::vector <int> quantities;} hologram;

void search_for_char(FILE *f, char c);
void print_similarity_descriptions();
void build_array(float ***array, hologram holo, const int NUM_BINS, const int DIFFUSION_RANGE);
float get_similarity(float ***array_i, float ***array_j, const int NUM_BINS, const int DIFFUSION_RANGE, const int SIMILARITY_MEASURE);

const int BIN_TAN = 1, MTW = 2, MTU = 3, MT_80_20 = 4, BARONI_URBANI = 5, CONT_TAN = 6, MTW_COMBO = 7, MTU_COMBO = 8, MISMATCH = 9, BRAY_CURTIS = 10, MAX_ALLOWED = 10;
const char* similarity_coefficient_names[10] = {
  "binary Tanimoto",
  "binary modified Tanimoto weighted (MTW)",
  "binary modified Tanimoto unweighted (MTU)",
  "binary modified Tanimoto 80:20 weighted presence:absence)",
  "binary Baroni-Urbani",
  "continuous Tanimoto",
  "continuous MTW",
  "continuous MTU",
  "continuous Tanimoto of mismatch",
  "continuous BrayCurtis"
};

void search_for_char(FILE *f, char target) {
	char c = getc(f);
	while(c!=target && c!=EOF) {
		c = getc(f);
	}
	if(c==EOF) {
		printf("ERROR: The required character (%c, with int value %d) was not found in this file.\n", target, (int)target);
		exit(EXIT_FAILURE);
	}
}

void print_similarity_descriptions() {
	printf("Where SIMILARITY_MEASURE is:\n");
  for(int i=0; i<MAX_ALLOWED; i++) {
    printf("\t%d = %s\n", i+1, similarity_coefficient_names[i]);
  }
}

void build_array(float ***array, hologram holo, const int NUM_BINS, const int DIFFUSION_RANGE) {
	int k,l,m;
	for(k=0; k<NUM_BINS; k++) {
		for(l=0; l<NUM_BINS; l++) {
			for(m=0; m<NUM_BINS; m++) {
				array[k][l][m] = 0;
			}
		}
	}
	int i, hi, hj, hk;
	for(i=0; i<holo.num_entries; i++) {
		hi = holo.grid_points.at(i).x;
		hj = holo.grid_points.at(i).y;
		hk = holo.grid_points.at(i).z;
		float value = (float)(holo.quantities.at(i));
		//loop over nearby points
		int i_diff, j_diff, k_diff;
		int new_hi, new_hj, new_hk;
		for(i_diff=-1*DIFFUSION_RANGE; i_diff<=DIFFUSION_RANGE; i_diff++) {
			new_hi = hi+i_diff;
			if(!(new_hi<0 || new_hi>=NUM_BINS)) {
				for(j_diff=-1*DIFFUSION_RANGE; j_diff<=DIFFUSION_RANGE; j_diff++) {
					new_hj = hj+j_diff;
					if(!(new_hj<0 || new_hj>=NUM_BINS)) {
						for(k_diff=-1*DIFFUSION_RANGE; k_diff<=DIFFUSION_RANGE; k_diff++) {
							new_hk = hk+k_diff;
							if(!(new_hk<0 || new_hk>=NUM_BINS)) {
								if(!(new_hk>new_hj)) { //i.e. don't write to the side of the grid which wasn't written to before!
									float distance = sqrt( (i_diff)*(i_diff) + (j_diff)*(j_diff) + (k_diff)*(k_diff) );
									distance+=1; //avoid dividing by zero
									array[new_hi][new_hj][new_hk] += value/distance; //increment this neighbouring point by the inverse distance times the value of this hologram point
								}
							}
						}
					}
				}
			}
		}
	}
}

float get_similarity(float ***array_i, float ***array_j, const int NUM_BINS, const int DIFFUSION_RANGE, const int SIMILARITY_MEASURE) {
	unsigned int k,l,m;
	float similarity;
	float sum_ab = 0.0, sum_aa = 0.0, sum_bb = 0.0;
	float sum_diff = 0.0;
	float presence_tan, absence_tan;
	int n11 = 0, n00 = 0, n = 0, ni = 0, nj = 0;
	float p_hat;
	float sum_a_plus_b = 0.0, sum_diff_a_b = 0.0;
	int sum_d = 0, sum_a = 0, sum_b = 0, sum_c = 0;
	//now loop over both arrays and update the sums for cont_tan
	for(k=0; k<NUM_BINS; k++) {
		for(l=0; l<NUM_BINS; l++) {
			for(m=0; m<NUM_BINS; m++) {
				if(!(m>l)) { //i.e. don't read from the empty side of the grid
					sum_ab += array_i[k][l][m] * array_j[k][l][m];
					sum_aa += array_i[k][l][m] * array_i[k][l][m];
					sum_bb += array_j[k][l][m] * array_j[k][l][m];
					sum_diff += (array_i[k][l][m] - array_j[k][l][m])*(array_i[k][l][m] - array_j[k][l][m]);
					n11 += (array_i[k][l][m]>0 && array_j[k][l][m]>0);
					n00 += ((!(array_i[k][l][m]>0)) && (!(array_j[k][l][m]>0)));
					ni += (array_i[k][l][m]>0);
					nj += (array_j[k][l][m]>0);
					n++; //num_features
					sum_b += (array_i[k][l][m]>0);
					sum_c += (array_j[k][l][m]>0);
					sum_a += (array_i[k][l][m]>0)*(array_j[k][l][m]>0);
					sum_d += ((array_i[k][l][m]>0)+(array_j[k][l][m]>0)==0);
				}
			}
		}
	}
	if(SIMILARITY_MEASURE==CONT_TAN) {
		float denominator = sum_aa+sum_bb-sum_ab;
		if(denominator<=0.0) similarity = -1;
		else similarity = sum_ab/denominator;
	} else if(SIMILARITY_MEASURE==BIN_TAN) {
		float denominator = (float)(ni+nj-n11);
		if(denominator<=0.0) similarity = -1;
		else similarity = ((float)(n11))/denominator;
	} else if(SIMILARITY_MEASURE==MISMATCH) {
		float denominator = sum_aa+sum_bb;
		if(denominator<=0.0) similarity = -1;
		else similarity = 1 - (sum_diff/denominator);
	} else if(SIMILARITY_MEASURE==MTU || SIMILARITY_MEASURE==MTW) {
		if(SIMILARITY_MEASURE==MTW) p_hat = ((float)(n11+(n-n00)))/((float)(2*n)); //weighted version
		else p_hat = 0.5; //unweighted version
		if(n-n00 == 0) presence_tan = 1; //as instructed in paper
		else presence_tan = ((float)(n11))/((float)(n-n00));
		if(n-n11 == 0) absence_tan = 1; //as instructed in paper
		else absence_tan = ((float)(n00))/((float)(n-n11));
		similarity = (((2.0-p_hat)/3.0)*presence_tan)+(((1.0+p_hat)/3.0)*absence_tan);
	} else if(SIMILARITY_MEASURE==MTU_COMBO || SIMILARITY_MEASURE==MTW_COMBO) {
		if(SIMILARITY_MEASURE==MTW_COMBO) p_hat = ((float)(n11+(n-n00)))/((float)(2*n)); //weighted version
		else p_hat = 0.5; //unweighted version
		float denominator = sum_aa+sum_bb-sum_ab;
		if(denominator <= 0) presence_tan = 1; //as instructed in paper
		else presence_tan = sum_ab/denominator;
		if(n-n11 == 0) absence_tan = 1; //as instructed in paper
		else absence_tan = ((float)(n00))/((float)(n-n11));
		similarity = (((2.0-p_hat)/3.0)*presence_tan)+(((1.0+p_hat)/3.0)*absence_tan);
	} else if(SIMILARITY_MEASURE==MT_80_20) {
		if(n-n00 == 0) presence_tan = 1; //as instructed in paper
		else presence_tan = ((float)(n11))/((float)(n-n00));
		if(n-n11 == 0) absence_tan = 1; //as instructed in paper
		else absence_tan = ((float)(n00))/((float)(n-n11));
		similarity = (0.8*presence_tan)+(0.2*absence_tan);
	} else if(SIMILARITY_MEASURE==BRAY_CURTIS) {
		if(sum_a_plus_b<=0.0) similarity = -1;
		else similarity = 1 - (sum_diff_a_b/sum_a_plus_b); //because this is a distance function
	} else if(SIMILARITY_MEASURE==BARONI_URBANI) {
		float denominator = sqrt(float(sum_a*sum_d))+sum_a+sum_b+sum_c;
		if(denominator<=0.0) similarity = -1;
		else similarity = sqrt(float((sum_a*sum_d)+sum_a))/denominator;
	}
	return similarity;
}

