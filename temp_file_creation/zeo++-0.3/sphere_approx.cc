//Sphere clusters 
//Author: Marielle Pinheiro
//March 28, 2013
//Create class AtomCluster that will replace a large atom with a cluster
//	of smaller atoms in a specified configuration

//#include "network.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <fstream>

#include "zeo_consts.h"
#include "networkstorage.h"
#include "sphere_approx.h"


// Bharat changed the SPHERE_APPROX_MIN_R from 0.6 to 0.5
#define SPHERE_APPROX_MIN_R 0.5
#define SPHERE_APPROX_MAX_R 2.8

#define LARGE_RATIO 1.3  // ratio between radii with triggers high accuract setting (replacement of an atom with a large cluster)

using namespace std;


//Constructor sets parameters for the small atom that is the center for the cluster
//This will be the origin for all translations
AtomCluster::AtomCluster(ATOM orgatm, double replacementSphereRadii) {

 smallSphereRadius = replacementSphereRadii;

 orgAtom = orgatm;

 atom_vector.clear();

 center.x = orgAtom.x;
 center.y = orgAtom.y;
 center.z = orgAtom.z;
 center.radius = smallSphereRadius;

}

/* old constructures to be removed
AtomCluster::AtomCluster(XYZ point) {
	center.x = point.x;
	center.y = point.y;
	center.z = point.z;
	center.radius = smallSphereRadius;
}

AtomCluster::AtomCluster(double x, double y, double z) {
	center.x = x;
	center.y = y;
	center.z = z;
	center.radius = smallSphereRadius;
}
*/


/* function that updates information about replaced atoms (small particle information is updated:
   mass update, moving to UC, calculation of fractional coordiantes) and copied to a new ATOM 
   vector provided by a pointer */
void AtomCluster::copyReplacementAtoms(ATOM_NETWORK *atmnet, int atomID, std::vector <ATOM> *newatoms){

 double newMass = orgAtom.mass/((double) atom_vector.size());
 double newCharge = orgAtom.charge/((double) atom_vector.size());

 for(unsigned int i = 0; i < atom_vector.size(); i++)
    {
    atom_vector[i].mass = newMass;
    atom_vector[i].charge = newCharge;
    atom_vector[i].radius = smallSphereRadius;
    atom_vector[i].type = orgAtom.type;
    atom_vector[i].specialID = orgAtom.specialID;

    Point pt(atom_vector[i].x, atom_vector[i].y, atom_vector[i].z);

    atmnet->shiftXYZInUC(pt);
    atom_vector[i].x = pt[0];
    atom_vector[i].y = pt[1];
    atom_vector[i].z = pt[2];

    pt = atmnet->xyz_to_abc(pt);
    atom_vector[i].a_coord = pt[0];
    atom_vector[i].b_coord = pt[1];
    atom_vector[i].c_coord = pt[2];

    newatoms->push_back(atom_vector[i]);
    atmnet->IDmapping.push_back(atomID);
    };


} // ends copyReplacementAtoms()



//Based on specified 3-letter code, generate vector of atoms for sphere cluster approximation
//This is the main function that calls all successive functions
void AtomCluster::replaceAtomByCluster(string cluster_type, double large_atom_radius){
	double shift = large_atom_radius - smallSphereRadius;

//Type: OCC: Octahedron + cube (15 spheres)
        if (cluster_type == "OCC"){
                for (int i = 0; i<15; i++){
                        atom_vector.push_back(center);
                }
                translate_cube(shift, 0);
                plusMinus_axes(shift, 8);
        }

//Type: FCC = 4 cubes (33 spheres) 
	if (cluster_type == "FCC" || cluster_type == "ACC"){

		//fills vector with spheres positioned at origin
		for (int i = 0; i<32; i++){
			atom_vector.push_back(center);
		}

		//Translate cube vertices for all 4 cubes
		translate_cube(shift, 0);
		rotate_xy_cube(shift, 8);
		rotate_xz_cube(shift, 16);
		rotate_yz_cube(shift, 24);

//Type: ACC = 4 cubes w/ axis spheres (39 spheres)
		if (cluster_type == "ACC"){
			vector <ATOM> axes;
			for (int i = 0; i<6; i++){
				atom_vector.push_back(center);
			}
			plusMinus_axes(shift, 32);
		}
		atom_vector.push_back(center);
	}

//Type: AQC = Axes, quadrants, cube (27 spheres)
//	Spheres along +- axes, 45 degrees from xy, xz, yz axes,
//	and cube with edges aligned to axes 
	else if (cluster_type == "AQC"){
		for (int i = 0; i<27; i++){
			atom_vector.push_back(center);
		}
		translate_axes(shift, 0);
		translate_cube(shift, 18);
	}

//Type: DDH = Dodecahedron (21 spheres); pentagonal faces
	else if (cluster_type == "DDH" || cluster_type == "TIH"){
		for (int i = 0; i<20; i++){
			atom_vector.push_back(center);
		}
		translate_dodecahedron(shift, 0);
//Type: TIH = Truncated icosahedron (33 spheres)
//	dodecahedron w/ center of faces filled w/ spheres
		if (cluster_type == "TIH"){
			vector <ATOM> ico;
			for (int i = 0; i<12; i++){
				atom_vector.push_back(center);
			}
			rotate_icosahedron(shift, 20);
		}
		atom_vector.push_back(center);
	}

//Type: ICH = Icosidodecahedron (31 spheres)
// Quasiregular polyhedron. pentagonal/triangular faces
	else if (cluster_type == "ICH" || cluster_type == "ICC"){
		for (int i = 0; i<30; i++){
			atom_vector.push_back(center);
		}
		translate_icosidodecahedron(shift, 0);
//Type: ICC = Icosidodecahedron w/ center face spheres (43 spheres)
		if (cluster_type == "ICC"){
	//		vector <ATOM> center_spheres;
		//	for (int i = 0; i<atom_vector.size(); i++){
		//		center_spheres.push_back(atom_vector.at(i));
		//	}
			calc_centerSpheres(shift, 0);
		}
		atom_vector.push_back(center);
	}

//Type: RIH = Rhombicosidodecahedron (61 spheres)
// Semiregular polyhedron. Square, pentagonal, triangular faces
	else if (cluster_type == "RIH"){
		for (int i = 0; i<61; i++){
			atom_vector.push_back(center);
		}

		translate_rhombi(shift, 0);
	}
//Type: S4 = 4 spheres in a spiral
        else if (cluster_type == "S4"){
                sphere_spiral(4, shift);
        }
//Type: S10 = 10 spheres in a spiral
        else if (cluster_type == "S10"){
                sphere_spiral(10, shift);
        }
//Type: S20 = 20 spheres in a spiral
        else if (cluster_type == "S20"){
                sphere_spiral(20, shift);
        }
//Type: S30 = 30 spheres in a spiral
        else if (cluster_type == "S30"){
                sphere_spiral(30, shift);
        }
//Type: S40 = 40 spheres in a spiral
        else if (cluster_type == "S40"){
                sphere_spiral(40, shift);
        }
//Type: S50 = 50 spheres in a spiral
        else if (cluster_type == "S50"){
                sphere_spiral(50, shift);
        }
//Type: S100 = 100 spheres in a spiral
        else if (cluster_type == "S100"){
                sphere_spiral(100, shift);
        }
//Type: S500 = 500 spheres in a spiral
        else if (cluster_type == "S500"){
                sphere_spiral(500, shift);
        }
//Type: S1000 = 1000 spheres in a spiral
        else if (cluster_type == "S1000"){
                sphere_spiral(1000, shift);
        }
//Type: S1000 = 1000 spheres in a spiral
        else if (cluster_type == "S10000"){
                sphere_spiral(10000, shift);
        }
}

//Shifts atom in designated direction

ATOM AtomCluster::translate_sphere(ATOM sphere, double x_step, double y_step, double z_step, int x_sign, int y_sign, int z_sign){
	ATOM new_sphere_pos; new_sphere_pos.radius = sphere.radius;
//sign: -1 = negative, 0 = no change, 1 = positive
	if (x_sign == -1) new_sphere_pos.x = sphere.x - x_step;
	else if (x_sign == 1) new_sphere_pos.x = sphere.x + x_step;
	else if (x_sign == 0) new_sphere_pos.x = sphere.x;

	if (y_sign == -1) new_sphere_pos.y = sphere.y - y_step;
	else if (y_sign == 1) new_sphere_pos.y = sphere.y + y_step;
	else if (y_sign == 0) new_sphere_pos.y = sphere.y;

	if (z_sign == -1) new_sphere_pos.z = sphere.z - z_step;
	else if (z_sign == 1) new_sphere_pos.z = sphere.z + z_step;
	else if (z_sign == 0) new_sphere_pos.z = sphere.z;

	return new_sphere_pos;
}

//translates along all +- axes
void AtomCluster::plusMinus_axes(double shift, int i){
	atom_vector.at(i) = translate_sphere(atom_vector.at(i), shift, 0, 0, 1, 0, 0);
	atom_vector.at(i+1) = translate_sphere(atom_vector.at(i+1), shift, 0, 0, -1, 0, 0);
	atom_vector.at(i+2) = translate_sphere(atom_vector.at(i+2), 0, shift, 0, 0, 1, 0);
	atom_vector.at(i+3) = translate_sphere(atom_vector.at(i+3), 0, shift, 0, 0, -1, 0);
	atom_vector.at(i+4) = translate_sphere(atom_vector.at(i+4), 0, 0, shift, 0, 0, 1);
	atom_vector.at(i+5) = translate_sphere(atom_vector.at(i+5), 0, 0, shift, 0, 0, -1);
}

//For specified vector shift (+-x, +-y, 0)
void AtomCluster::plusMinus_xy(double x, double y, int i){
	atom_vector.at(i) = translate_sphere(atom_vector.at(i), x, y, 0, 1, 1, 0);
	atom_vector.at(i+1) = translate_sphere(atom_vector.at(i+1), x, y, 0, -1, 1, 0);
	atom_vector.at(i+2) = translate_sphere(atom_vector.at(i+2), x, y, 0, -1, -1, 0);
	atom_vector.at(i+3) = translate_sphere(atom_vector.at(i+3), x, y, 0, 1, -1, 0);
}

//For specified vector shift (+-x, 0, +-z)
void AtomCluster::plusMinus_xz(double x, double z, int i){
	atom_vector.at(i) = translate_sphere(atom_vector.at(i), x, 0, z, 1, 0, 1);
	atom_vector.at(i+1) = translate_sphere(atom_vector.at(i+1), x, 0, z, -1, 0, 1);
	atom_vector.at(i+2) = translate_sphere(atom_vector.at(i+2), x, 0, z, -1, 0, -1);
	atom_vector.at(i+3) = translate_sphere(atom_vector.at(i+3), x, 0, z, 1, 0, -1);
}

//For specified vector shift (0, +-y, +-z)
void AtomCluster::plusMinus_yz(double y, double z, int i){
	atom_vector.at(i) = translate_sphere(atom_vector.at(i), 0, y, z, 0, 1, 1);
	atom_vector.at(i+1) = translate_sphere(atom_vector.at(i+1), 0, y, z, 0, -1, 1);
	atom_vector.at(i+2) = translate_sphere(atom_vector.at(i+2), 0, y, z, 0, -1, -1);
	atom_vector.at(i+3) = translate_sphere(atom_vector.at(i+3), 0, y, z, 0, 1, -1);
}

//For specified vector shift (+-x, +-y, +-z)
void AtomCluster::plusMinus_all(double x, double y, double z, int i){
	atom_vector.at(i) = translate_sphere(atom_vector.at(i), x, y, z, 1, 1, 1); 
	atom_vector.at(i+1) = translate_sphere(atom_vector.at(i+1), x, y, z, -1, 1, 1);
	atom_vector.at(i+2) = translate_sphere(atom_vector.at(i+2), x, y, z, 1, -1, 1);
	atom_vector.at(i+3) = translate_sphere(atom_vector.at(i+3), x, y, z, 1, 1, -1);
	atom_vector.at(i+4) = translate_sphere(atom_vector.at(i+4), x, y, z, -1, -1, 1);
	atom_vector.at(i+5) = translate_sphere(atom_vector.at(i+5), x, y, z, -1, 1, -1);
	atom_vector.at(i+6) = translate_sphere(atom_vector.at(i+6), x, y, z, 1, -1, -1);
	atom_vector.at(i+7) = translate_sphere(atom_vector.at(i+7), x, y, z, -1, -1, -1);
}

//cube with corners in xyz quadrants
void AtomCluster::translate_cube(double shift, int i){
	double corner_vector = 1.0/(sqrt(3));
	double cube_shift = corner_vector * shift;
	plusMinus_all(cube_shift, cube_shift, cube_shift, i);
}

//cube with corners along xz and yz
void AtomCluster::rotate_xy_cube(double shift, int i){

	double axis_shift = shift * sqrt(2.0/3.0);
	double height_shift = shift * (1.0/sqrt(3));

	plusMinus_yz(axis_shift, height_shift, i);
	plusMinus_xz(axis_shift, height_shift, i+4);

}

//cube with corners along xy and yz
void AtomCluster::rotate_xz_cube(double shift, int i){

	double axis_shift = shift * sqrt(2.0/3.0);
	double height_shift = shift * (1.0/sqrt(3));

	plusMinus_xy(axis_shift, height_shift, i);
	plusMinus_yz(height_shift, axis_shift, i+4);

}

//cube with corners along yz and xz
void AtomCluster::rotate_yz_cube(double shift, int i){
	double axis_shift = shift * sqrt(2.0/3.0);
	double height_shift = shift * (1.0/sqrt(3));

	plusMinus_xy(height_shift, axis_shift, i);
	plusMinus_xz(height_shift, axis_shift, i+4);

}

//18-sphere axes/quadrant translation (combines w/cube)
void AtomCluster::translate_axes(double shift, int i){
	double axis_shift = shift;
	double quad_shift = 0.5 * sqrt(2) * shift;

	plusMinus_axes(axis_shift, i);
	plusMinus_xy(quad_shift, quad_shift, i+6);
	plusMinus_xz(quad_shift, quad_shift, i+10);
	plusMinus_yz(quad_shift, quad_shift, i+14);
}

//Dodecahedron (20 vertices)
void AtomCluster::translate_dodecahedron(double shift, int i){
	double scale = 1/(sqrt(3));
	double golden_ratio = (1+sqrt(5))/2.0;
	double inverse_ratio = 1.0/golden_ratio;
	double scaled_shift = scale * shift;
	double scaled_ratio = golden_ratio * scaled_shift;
	double scaled_inverse = inverse_ratio * scaled_shift;

	double half_ratio = 0.5*shift * (1+sqrt(5))/2.0;
	double half_shift = 0.5*shift;

	plusMinus_all(scaled_shift, scaled_shift, scaled_shift, i);
	plusMinus_xy(scaled_inverse, scaled_ratio, i+8);
	plusMinus_yz(scaled_inverse, scaled_ratio, i+12);
	plusMinus_xz(scaled_ratio, scaled_inverse, i+16);
}

//Icosahedron (12 vertices) that fills faces of dodecahedron
void AtomCluster::rotate_icosahedron(double shift, int i){
	double half_ratio = 0.5*shift * (1+sqrt(5))/2.0;
	double half_shift = 0.5*shift;

	plusMinus_xy(half_ratio, half_shift, i);
	plusMinus_yz(half_ratio, half_shift, i+4);
	plusMinus_xz(half_shift, half_ratio, i+8);
}

//Icosidodecahedron (Triangular/pentagonal faces, 30 vertices)
void AtomCluster::translate_icosidodecahedron(double shift, int i){
	double golden_ratio = (1+sqrt(5))/2.0;
	double half_ratio = golden_ratio * 0.5;
	double oth_ratio = (1+golden_ratio)*0.5;

	double scale = 1.0/golden_ratio;
	double scaled_shift = shift * scale;
	double scaled_ratio = golden_ratio * scaled_shift;
	double scaled_half = half_ratio * scaled_shift;
	double scaled_oth = oth_ratio * scaled_shift;

	plusMinus_axes(scaled_ratio, i);
	plusMinus_all(0.5 * scaled_shift, scaled_half, scaled_oth, i+6);
	plusMinus_all(scaled_half, scaled_oth, 0.5 * scaled_shift, i+14);
	plusMinus_all(scaled_oth, 0.5 * scaled_shift, scaled_half, i+22);
}

//calculates center of pentagon (to fill faces of icosidodecahedron)
ATOM AtomCluster::calc_center(vector <ATOM> old_coords, int i_0, int i_1, int i_2, int i_3, int i_4, double shift){
	ATOM origin = center;

	vector <ATOM> coords;

	coords.push_back(old_coords.at(i_0));
	coords.push_back(old_coords.at(i_1));
	coords.push_back(old_coords.at(i_2));
	coords.push_back(old_coords.at(i_3));
	coords.push_back(old_coords.at(i_4));

	double x_sum = 0;
	double y_sum = 0;
	double z_sum = 0;
	for (int i = 0; i<coords.size(); i++){
		x_sum = x_sum + coords.at(i).x;
		y_sum = y_sum + coords.at(i).y;
		z_sum = z_sum + coords.at(i).z;
	}

	double n = coords.size();

	ATOM center_sphere;

//translation vector
	double x = (x_sum/n)-origin.x;
	double y = (y_sum/n)-origin.y;
	double z = (z_sum/n)-origin.z;
	
	double norm = sqrt(x*x + y*y + z*z);
	double scaled_shift = shift/norm;

	center_sphere.x = x*scaled_shift + origin.x;
	center_sphere.y = y*scaled_shift + origin.y;
	center_sphere.z = z*scaled_shift + origin.z;
	center_sphere.radius = smallSphereRadius;

	return center_sphere;
}

//Place spheres in centers of icosidodecahedron pentagonal faces
void AtomCluster::calc_centerSpheres(double shift, int i){

	ATOM new_center;
//Top and bottom 4 (w/ z vertices) 
//Uses existing coordinates from icosidodecahedron 
//Pentagon 1: (0, 0, T), (1/2, +- T/2, (1+T)/2), ((1+T)/2, +-1/2, T/2) 
	new_center = calc_center(atom_vector, i+4, i+6, i+8, i+22, i+24, shift);
	atom_vector.push_back(new_center);
//Pentagon 2: (0, 0, T), (-1/2, +- T/2, (1+T)/2), (-(1+T)/2, +-1/2, T/2) 
	new_center = calc_center(atom_vector, i+4, i+7, i+10, i+23, i+26, shift);
	atom_vector.push_back(new_center);
//Pentagon 3: (0, 0, -T), (1/2, +- T/2, -(1+T)/2), ((1+T)/2, +-1/2, -T/2) 
	new_center = calc_center(atom_vector, i+5, i+9, i+12, i+25, i+28, shift);
	atom_vector.push_back(new_center);
//Pentagon 4: (0, 0, -T), (-1/2, +- T/2, -(1+T)/2), (-(1+T)/2, +-1/2, -T/2) 
	new_center = calc_center(atom_vector, i+5, i+11, i+13, i+27, i+29, shift);
	atom_vector.push_back(new_center);

//mid-top 4 (w/ y vertices)
//Pentagon 5: (+-1/2, T/2, (1+T)/2), (+-T/2, (1+T)/2, 1/2), (0, T, 0)
	new_center = calc_center(atom_vector, i+2, i+6, i+7, i+14, i+15, shift);
	atom_vector.push_back(new_center);
//Pentagon 6: (+-1/2, T/2, -(1+T)/2), (+-T/2, (1+T)/2, -1/2), (0, T, 0)
	new_center = calc_center(atom_vector, i+2, i+9, i+11, i+17, i+19, shift);
	atom_vector.push_back(new_center);
//Pentagon 7: (+-1/2, -T/2, (1+T)/2), (+-T/2, -(1+T)/2, 1/2), (0, -T, 0)
	new_center = calc_center(atom_vector, i+3, i+8, i+10, i+16, i+18, shift);
	atom_vector.push_back(new_center);
//Pentagon 8: (+-1/2, -T/2, -(1+T)/2), (+-T/2, -(1+T)/2, -1/2), (0, -T, 0)
	new_center = calc_center(atom_vector, i+3, i+12, i+13, i+20, i+21, shift);
	atom_vector.push_back(new_center);

//Middle 4 (w/ x vertices) 
//Pentagon 9: (T, 0, 0), (T/2, (1+T)/2, +-1/2), ((1+T)/2, 1/2, +-T/2)
	new_center = calc_center(atom_vector, i, i+14, i+17, i+22, i+25, shift);
	atom_vector.push_back(new_center);
//Pentagon 10: (T, 0, 0), (T/2, -(1+T)/2, +-1/2), ((1+T)/2, -1/2, +-T/2)
	new_center = calc_center(atom_vector, i, i+16, i+20, i+24, i+28, shift);
	atom_vector.push_back(new_center);
//Pentagon 11: (-T, 0, 0), (-T/2, (1+T)/2, +-1/2), (-(1+T)/2, 1/2, +-T/2)
	new_center = calc_center(atom_vector, i+1, i+15, i+19, i+23, i+27, shift);
	atom_vector.push_back(new_center);
//Pentagon 12: (-T, 0, 0), (-T/2, -(1+T)/2, +-1/2), (-(1+T)/2, -1/2, +-T/2)
	new_center = calc_center(atom_vector, i+1, i+18, i+21, i+26, i+29, shift);
	atom_vector.push_back(new_center);
	
}

//Rhombicosidodecahedron, 60 vertices
void AtomCluster::translate_rhombi(double shift, int i){
	double golden_ratio = (1+sqrt(5))/2.0;
	double ratio_square = golden_ratio * golden_ratio;
	double ratio_cube = ratio_square * golden_ratio;
	double length = sqrt(2+(ratio_cube*ratio_cube));
	double scale = 1.0/length;

	double scaled_shift = shift * scale;
	double scaled_ratio = scaled_shift * golden_ratio;
	double scaled_cube = scaled_shift * ratio_cube;
	double scaled_square = scaled_shift * ratio_square;
	double scaled_plus = scaled_shift * (2.0 + golden_ratio);

	plusMinus_all(scaled_shift, scaled_shift, scaled_cube, i);
	plusMinus_all(scaled_cube, scaled_shift, scaled_shift, i+8);
	plusMinus_all(scaled_shift, scaled_cube, scaled_shift, i+16);
	plusMinus_all(scaled_square, scaled_ratio, 2.0*scaled_ratio, i+24);
	plusMinus_all(2.0*scaled_ratio, scaled_square, scaled_ratio, i+32);
	plusMinus_all(scaled_ratio, 2.0*scaled_ratio, scaled_square, i+40);
	plusMinus_xz(scaled_plus, scaled_square, i+48);
	plusMinus_xy(scaled_square, scaled_plus, i+52);
	plusMinus_yz(scaled_square, scaled_plus, i+56);
}

//Generates spheres in a spiral 
void AtomCluster::sphere_spiral(double n, double shift){
        double golden_angle = PI * (3.0-sqrt(5));

        for (int i = 0; i<n; i++){
                double theta = golden_angle * i;
                double z_shift = (1-(1/n))*(1-((2*i)/(n-1)));
                double z_sq = z_shift * z_shift;
                double radius = sqrt(1-z_sq);
                double x_shift = radius * cos(theta);
                double y_shift = radius * sin(theta);

                ATOM new_sphere = translate_sphere(center, shift * x_shift, shift * y_shift, shift * z_shift, 1, 1, 1);
                atom_vector.push_back(new_sphere);
        }
        atom_vector.push_back(center);
}


//outputs an xyz format file with the approximation sphere coordinates
void AtomCluster::print_xyz_coords(FILE *output){
	int size = atom_vector.size();
	fprintf(output, "%d\n\n", size);
	for (int a = 0; a<atom_vector.size(); a++){
		ATOM atom = atom_vector.at(a);
		fprintf(output, "H %f %f %f %f\n", atom.x, atom.y, atom.z, atom.radius);
	}
	fclose(output);
}












/* Function that analyzes distribution of atomic radii and execute replacement of large atoms 
   with clusters of small ones */
void setupHighAccuracyAtomNetwork(ATOM_NETWORK *atmnet, std::string AccSetting){

 double minr,maxr;

/* analyze the original atom network */
 for(unsigned int i = 0; i < atmnet->atoms.size(); i++)
   {
   if(i == 0) {
     minr = atmnet->atoms.at(i).radius; maxr = minr;
     } else{
     if( atmnet->atoms.at(i).radius < minr) minr = atmnet->atoms.at(i).radius;
     if( atmnet->atoms.at(i).radius > maxr) maxr = atmnet->atoms.at(i).radius;
    };
   };
 
 cout << "Radii analysis: the smallest atom r = " << minr << " while the largest atoms r = " << maxr << ".\n";


 if(minr < SPHERE_APPROX_MIN_R || maxr > SPHERE_APPROX_MAX_R)
   {
   cerr << "HIGH ACCURACY CANNOT BE APPLIED!\n"
        << "The current version of high accuracy routine is using parameters derived from analysis of \n" 
        << "CCDC database. The radii defined for the current system are outside of predefined range.\n"
        << "If you are a pro-user, change the source code (sphere_approx.cc) or contact the authors.\n"
        << "Exiting the -ha routines without any changes..." << "\n";
   }else
 {

/* start replacing the atoms */

 std::vector <ATOM> newAtomList;

 for(unsigned int i = 0; i < atmnet->atoms.size(); i++)
   {

   if(atmnet->atoms.at(i).radius == minr)
     {
     // in case of the smallest atoms, copy them without modifications
     newAtomList.push_back(atmnet->atoms[i]);
     atmnet->IDmapping.push_back(i);
     }else
     {
     // in case of larger atoms, run atom modification routines

     if(AccSetting == "OCC" ||
        AccSetting == "FCC" ||
        AccSetting == "ACC" ||
        AccSetting == "AQC" ||
        AccSetting == "DDH" ||
        AccSetting == "TIH" ||
        AccSetting == "ICH" ||
        AccSetting == "ICC" ||
        AccSetting == "RIH" ||
        AccSetting == "S4" ||
        AccSetting == "S10" ||
        AccSetting == "S20" ||
        AccSetting == "S30" ||
        AccSetting == "S40" ||
        AccSetting == "S50" ||
        AccSetting == "S100" ||
        AccSetting == "S500" ||
        AccSetting == "S1000" ||
        AccSetting == "S10000")
        {
        AtomCluster cluster(atmnet->atoms[i], minr);
        cluster.replaceAtomByCluster(AccSetting, atmnet->atoms.at(i).radius);
        cluster.copyReplacementAtoms(atmnet, i, &newAtomList);
        }
     else if(AccSetting == "HI")
        {
        AccSetting = "S50";
        AtomCluster cluster(atmnet->atoms[i], minr);
        cluster.replaceAtomByCluster(AccSetting, atmnet->atoms.at(i).radius);
        cluster.copyReplacementAtoms(atmnet, i, &newAtomList);
        }
     else if(AccSetting == "MED")
        {
        AccSetting = "S30";
        AtomCluster cluster(atmnet->atoms[i], minr);
        cluster.replaceAtomByCluster(AccSetting, atmnet->atoms.at(i).radius);
        cluster.copyReplacementAtoms(atmnet, i, &newAtomList);
        }
     else if(AccSetting == "LOW")
        {
        AccSetting = "S10";
        AtomCluster cluster(atmnet->atoms[i], minr);
        cluster.replaceAtomByCluster(AccSetting, atmnet->atoms.at(i).radius);
        cluster.copyReplacementAtoms(atmnet, i, &newAtomList);
        }
     else // equal to (AccSetting == "DEF")
        {
        if(atmnet->atoms.at(i).radius/minr < LARGE_RATIO)
          {
          AccSetting = "S30";
          } else{ 
          AccSetting = "S50";
          };
        AtomCluster cluster(atmnet->atoms[i], minr);
        cluster.replaceAtomByCluster(AccSetting, atmnet->atoms.at(i).radius);
        cluster.copyReplacementAtoms(atmnet, i, &newAtomList);
        };

     };
   }; // ends loop over all atoms in original atom list


/* finishing with replacement of old atom list with the new one */

 atmnet->atoms = newAtomList;
 atmnet->numAtoms = newAtomList.size();

 }; // ends if statement checking if minr and maxr are within the accaptable range

} // ends setupHighAccuracyAtomNetwork()











