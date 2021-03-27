/* Distance grid functions 
*  Added by Marielle Pinheiro, Richard Martin and Maciej Haranczyk
*  Fall 2012/Winter 2013
*
*  Distance grids are used for analysis and visualization
*  Grids are written in two formats: BOV (for Visit) and Gaussian cube (for VMD)
*  
*  BOV distance grids are written in a number of ways depending on distance function definition
*  the latter is done to test Voronoi approximation and it commented out in the official release
*
*
*
*  BOV Grid example, general scenario
*  The variables in this scenario assume that we will have information
*  found in cssr files: Unit cell dimensions a, b, c, and unit cell
*  angles alpha, beta, gamma. 
*  F: for location p, sphere center s--> d(s,p) < d(s,x_i) - r_i
*  G: for location p, sphere center s--> [d(s,p)]^2 < [d(s,x_i)]^2 - (r_i)^2
*  H: for location p, sphere center s--> sqrt(G+r_bar^2)-r_bar, where r_bar is some standard atom radius (below 1.35)
*
*
*
*  Gaussian cube file
*  only standard distance is considered
*/

//#include "network.h"
#include <cmath>
#include <cstring>
#include <fstream>

#include "grid.h"
#include "geometry.h"
#include "networkinfo.h"
#include "networkaccessibility.h"
#define GRIDRES 0.15  // standard low accuracy, fast to calculate grid resolution
//#define GRIDRES 0.05  // increase grid resolution for better contours but at significantly greater expense (cubic increase in cost)

using namespace std;


/* BOV grid function was originally written by Richard Martin for another project (Int. J. High Perf Comput Appl 26 347-357 2012)
*  therefore it uses slightly different convention for box definition. Therefore you notice conversion of variables in the following 
*  lines of code */
void generateBOVGrid(ATOM_NETWORK *atmnet, string name_f_dist, string name_g_dist, string name_h_dist, string name_f_bov, string name_g_bov, string name_h_bov) {

//-----declare large 3D array as a 3D pointer (i.e., ***) - we have to declare it as pointers otherwise we cannot declare an arbitrarily large grid - but if the grid is too large, you may not have enough RAM to store it, so watch the system monitor!
  double ***grid_distance;

//-----set the extent of the grid. If we know coordinate H (i.e. the coordinate in the abc plane that has all non-0 coordinates) we will know the maximum extent. 
//For my purposes, I defined a cube in terms of corners. A, B, and C correspond to defined unit cell vectors.
//D is the bottom corner between the B and C axes. E is the top corner between the A and C axes.
//F is the top corner between the B and C axes. G is the outermost corner. The innermost corner is simply (0, 0, 0).

  double origin_x = 0; 
  double origin_y = 0; 
  double origin_z = 0; 
  double Ax  = atmnet->v_a.x;
  double Bx = atmnet->v_b.x;
  double By  = atmnet->v_b.y;
  double Cx = atmnet->v_c.x;
  double Cy = atmnet->v_c.y;
  double Cz  = atmnet->v_c.z;

  double Dx = Ax + Bx; double Dy = By; double Dz = 0;
  double Ex = Ax + Cx; double Ey = Cy; double Ez = Cz;
  double Fx = Bx + Cx; double Fy = Cy; double Fz = Cz;
  double Gx = Ax + Bx + Cx; double Gy = By + Cy; double Gz = Cz;

  vector <double> xCoords;
  xCoords.push_back(origin_x); xCoords.push_back(Ax); xCoords.push_back(Bx); xCoords.push_back(Cx); xCoords.push_back(Dx); xCoords.push_back(Ex); xCoords.push_back(Fx); xCoords.push_back(Gx);
  vector <double> yCoords;
  yCoords.push_back(origin_y);yCoords.push_back(By); yCoords.push_back(Cy); yCoords.push_back(Gy);
	vector <double> zCoords;
	zCoords.push_back(origin_y); zCoords.push_back(Cz);

  double xMin = 1000; double xMax = 0;
  for (int a=0; a<xCoords.size(); a++){
    xMin = min(xMin, xCoords.at(a));
    xMax = max(xMax, xCoords.at(a));
  }

  double yMin = 1000; double yMax = 0;
  for (int b=0; b<yCoords.size(); b++){
    yMin = min(yMin, yCoords.at(b));
    yMax = max(yMax, yCoords.at(b));
  }
  double zMin = 1000; double zMax = 0;
  for (int c=0; c<zCoords.size(); c++){
    zMin = min(zMin, zCoords.at(c));
    zMax = max(zMax, zCoords.at(c));
  }

  double XGridDist = xMax - xMin; double YGridDist = yMax - yMin; double ZGridDist = zMax - zMin;


//-----set size that grid will take, i.e. the number of grid points in each dimension
  
  double tempXgrid = XGridDist/GRIDRES; double tempYgrid = YGridDist/GRIDRES; double tempZgrid = ZGridDist/GRIDRES;

  int x_grid_steps = ceil(tempXgrid);
  int y_grid_steps = ceil(tempYgrid); 
  int z_grid_steps = ceil(tempZgrid);

  double shiftInX = x_grid_steps - tempXgrid;
  double shiftInY = y_grid_steps - tempYgrid;
  double shiftInZ = z_grid_steps - tempZgrid;

  double xGridRes = XGridDist/x_grid_steps; double yGridRes = YGridDist/y_grid_steps; double zGridRes = ZGridDist/z_grid_steps;

  x_grid_steps++;
  y_grid_steps++;
  z_grid_steps++;

//-----create large 3D array
  printf("Declaring 3D array.\n\n");
	grid_distance = new double**[x_grid_steps];
	for(int i=0; i<x_grid_steps; i++) {
		grid_distance[i] = new double*[y_grid_steps];
		for(int j=0; j<y_grid_steps; j++) {
			grid_distance[i][j] = new double[z_grid_steps];
		}
	}

//-----fill it with values for the F function
  printf("Filling 3D array F with values.\n\n");
	for(int i=0; i<x_grid_steps; i++) {
		for(int j=0; j<y_grid_steps; j++) {
			for(int k=0; k<z_grid_steps; k++) {
				grid_distance[i][j][k] = calculate_distance_function(atmnet, i,j,k, xMin, yMin, 0, xGridRes, yGridRes, zGridRes, x_grid_steps, y_grid_steps, z_grid_steps, 'f'); 
			}
		}
	}

//-----Print structure_f.distances and structure_f.bov
  printf("Printing F grid.\n\n");
  FILE *file_distances_f, *file_bov_f; 
  char* distance_filename_f = new char [100];
  strcpy(distance_filename_f, name_f_dist.c_str()); 
  file_distances_f = fopen(distance_filename_f,"w");
  write_distances(file_distances_f, grid_distance, x_grid_steps, y_grid_steps, z_grid_steps);
  fclose(file_distances_f); 
  char* bov_filename_f = new char [100];
  strcpy(bov_filename_f, name_f_bov.c_str()); 
  file_bov_f = fopen(bov_filename_f,"w");
  write_bov(file_bov_f, distance_filename_f, x_grid_steps, y_grid_steps, z_grid_steps, xMin, yMin, zMin, XGridDist, YGridDist, ZGridDist);
  fclose(file_bov_f);

//-----fill it with values for the G function
/*
  printf("Filling 3D array G with values.\n\n");
	for(int i=0; i<x_grid_steps; i++) {
		for(int j=0; j<y_grid_steps; j++) {
			for(int k=0; k<z_grid_steps; k++) {
				grid_distance[i][j][k] = calculate_distance_function(atmnet, i,j,k, xMin, yMin, 0, xGridRes, yGridRes, zGridRes, x_grid_steps, y_grid_steps, z_grid_steps, 'g'); 
			}
		}
	}

//-----Print structure_g.distances and structure_g.bov
  printf("Printing G grid.\n\n");
  FILE *file_distances_g, *file_bov_g; 
  char* distance_filename_g = new char [100];
  strcpy(distance_filename_g, name_g_dist.c_str()); 
  file_distances_g = fopen(distance_filename_g,"w");
  write_distances(file_distances_g, grid_distance, x_grid_steps, y_grid_steps, z_grid_steps);
  fclose(file_distances_g); 
  char* bov_filename_g = new char [100];
  strcpy(bov_filename_g, name_g_bov.c_str()); 
  file_bov_g = fopen(bov_filename_g,"w");
  write_bov(file_bov_g, distance_filename_g, x_grid_steps, y_grid_steps, z_grid_steps, xMin, yMin, zMin, XGridDist, YGridDist, ZGridDist);
  fclose(file_bov_g);
*/


//-----fill it with values for the H function
/*
  printf("Filling 3D array H with values.\n\n");
	for(int i=0; i<x_grid_steps; i++) {
		for(int j=0; j<y_grid_steps; j++) {
			for(int k=0; k<z_grid_steps; k++) {
				grid_distance[i][j][k] = calculate_distance_function(atmnet, i,j,k, xMin, yMin, 0, xGridRes, yGridRes, zGridRes, x_grid_steps, y_grid_steps, z_grid_steps, 'h'); 
			}
		}
	}


//-----Print structure_h.distances and structure_h.bov

  printf("Printing H grid.\n\n");
  FILE *file_distances_h, *file_bov_h; 
  char* distance_filename_h = new char [100];
  strcpy(distance_filename_h, name_h_dist.c_str()); 
  file_distances_h = fopen(distance_filename_h,"w");
  write_distances(file_distances_h, grid_distance, x_grid_steps, y_grid_steps, z_grid_steps);
  fclose(file_distances_h); 
  char* bov_filename_h = new char [100];
  strcpy(bov_filename_h, name_h_bov.c_str()); 
  file_bov_h = fopen(bov_filename_h,"w");
  write_bov(file_bov_h, distance_filename_h, x_grid_steps, y_grid_steps, z_grid_steps, xMin, yMin, zMin, XGridDist, YGridDist, ZGridDist);
  fclose(file_bov_h);
*/


//-----free memory!
  printf("Freeing memory.\n\n");
	for(int i=0; i<x_grid_steps; i++) {
		for(int j=0; j<y_grid_steps; j++) {
			delete[] grid_distance[i][j];
		}
    delete[] grid_distance[i];
	}
  delete[] grid_distance; //important to have one "delete[]" per "new"
  delete[] distance_filename_f;
  //delete[] distance_filename_g;
  //delete[] distance_filename_h;
  delete[] bov_filename_f;
  //delete[] bov_filename_g;
  //delete[] bov_filename_h;

//---end successfully
  printf("Program complete.\n\n");


}

/* writes grid file */
void write_distances(FILE *f, double ***grid, int x_grid_steps, int y_grid_steps, int z_grid_steps) {
	vector <double> visit_friendly_array;
	for(int k=0; k<z_grid_steps; k++) {
		for(int j=0; j<y_grid_steps; j++) {
			for(int i=0; i<x_grid_steps; i++) { //for some reason, the convention is to write out the data in the z,y,x order, i.e., when you read the data, x iterates the fastest, and z the slowest
				visit_friendly_array.push_back(grid[i][j][k]);
			}
		}
	}
	fwrite(&visit_friendly_array[0], sizeof(double), x_grid_steps*y_grid_steps*z_grid_steps, f);
}

/* writes BOV grid info file */
void write_bov(FILE *f, char *output_distances_name, int x_grid_steps, int y_grid_steps, int z_grid_steps, double xMin, double yMin, double zMin, double x_box_size, double y_box_size, double z_box_size) {
	fprintf(f, "TIME: 99\nDATA_FILE: %s\nDATA_SIZE: %d %d %d\nDATA_FORMAT: DOUBLE\nVARIABLE: time\nDATA_ENDIAN: LITTLE\nCENTERING: nodal\nBRICK_ORIGIN: %f. %f. %f.\nBRICK_SIZE: %f. %f. %f.", output_distances_name, x_grid_steps, y_grid_steps, z_grid_steps, xMin, yMin, zMin, x_box_size, y_box_size, z_box_size);
}



/* calculates distance function for BOV grid definition
   there are 3 distance definitions possible but two of them are commented out in the original function */
double calculate_distance_function(ATOM_NETWORK *network, int i, int j, int k, double minX, double minY, double minZ, double xGridRes, double yGridRes, double zGridRes, int x_grid_steps, int y_grid_steps, int z_grid_steps, char gridtype) {

  double xPosition = minX + xGridRes*(i);
  double yPosition = minY + yGridRes*(j);
  double zPosition = minZ + zGridRes*(k);
  double dist;

  Point fractionalCoord = network->xyz_to_abc(xPosition, yPosition, zPosition);

  bool inStructure = (fractionalCoord[0]>=-0.01 && fractionalCoord[0]<=1.01 && fractionalCoord[1]>=-0.01 && fractionalCoord[1]<=1.01 && fractionalCoord[2]>=-0.01 && fractionalCoord[2]<=1.01);
  if (!inStructure) dist = 0;
  else {
    double minAtomDist = 1000;
    for (int i=0; i<network->numAtoms; i++){
      double atmDist = network->calcDistanceXYZ(xPosition, yPosition, zPosition, network->atoms[i].x, network->atoms[i].y, network->atoms[i].z);
      double radius = network->atoms[i].radius;
      if (gridtype=='f') atmDist -= radius;
      else{
        atmDist = atmDist*atmDist - radius*radius ;
        if (gridtype=='h') atmDist = sqrt(atmDist + 1.35*1.35)-1.35;
      }
      minAtomDist = min(minAtomDist, atmDist);
    }
    dist = minAtomDist;
  }
  return dist;
}






/*   Gausian cube file section */

void generateGaussianGrid(ATOM_NETWORK *atmnet, string cubefilename, bool angstrom_to_bohr, bool useMassFlag){

 GaussianCube cube(atmnet);

 cube.calculateDistanceGrid(atmnet);
 cube.writeGrid(atmnet,cubefilename, angstrom_to_bohr, useMassFlag);
 cube.deinit();

}


/*   Gausian cube file with accessibility information  */
void generateGaussianGridWithAccessibilityInfo(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgAtomnet, bool highAccuracy, double probe_radius, string cubefilename, bool angstrom_to_bohr, bool useMassFlag){

 GaussianCube cube(atmnet);

 cube.calculateDistanceGridWithAccessibilityInfo(atmnet, orgAtomnet, highAccuracy, probe_radius);

 if(highAccuracy==true) cube.writeGrid(orgAtomnet,cubefilename, angstrom_to_bohr, useMassFlag);
                  else  cube.writeGrid(atmnet,cubefilename, angstrom_to_bohr, useMassFlag);
 cube.deinit();

}



/* calculates a grid with 3D histogram of points defined by fractional coordinate (specified by inputfile) */
void calculateAverageGrid(ATOM_NETWORK *atmnet, std::string inputfilename, std::string cubefilename, bool angstrom_to_bohr, bool useMassFlag){

 GaussianCube cube(atmnet);

 cube.loadHistogramData(inputfilename);
 cube.writeGrid(atmnet,cubefilename, angstrom_to_bohr, useMassFlag);
 cube.deinit();

} //ends calculateAverageGrid()



/* calculates a grid with 3D histogram of points defined by fractional coordinate (specified by inputfile) */
void calculateAverageGridPerFrame(ATOM_NETWORK *atmnet, std::string inputfilename, std::string cubefilename, bool angstrom_to_bohr, bool useMassFlag){

 GaussianCube cube(atmnet);

 cube.loadHistogramDataPerFrame(inputfilename);
 cube.writeGrid(atmnet,cubefilename, angstrom_to_bohr, useMassFlag);
 cube.deinit();

} //ends calculateAverageGrid()




/*  Below are functions for GaussianCube class */

GaussianCube::GaussianCube(ATOM_NETWORK *atmnet){

 // calculate desired cube size
 na=(int)ceil(atmnet->a/GRIDRES)+1;
 nb=(int)ceil(atmnet->b/GRIDRES)+1;
 nc=(int)ceil(atmnet->c/GRIDRES)+1;

 // define the grid and allocate memory
 gridsize = na*nb*nc;

 if(gridsize<2)
    {
    cout << "Grid size is 1 or less. Aborting..." << "\n";
//    abort();
    } else
    {
    cout << "Gaussian cube grid - " << na << " x " << nb << " x " << nc << " = " << gridsize << " points.\n";
    };

 allocate(na, nb, nc);

 o = XYZ (0.0,0.0,0.0);
 va=atmnet->v_a;
 vb=atmnet->v_b;
 vc=atmnet->v_c;

 va=va.scale(1.0/(na-1.0)); // substracting 1.0 from na in denominator to make sure va reaches 1.0 after na steps (starting from 0)
 vb=vb.scale(1.0/(nb-1.0));
 vc=vc.scale(1.0/(nc-1.0));

}

/* saves the grid into Gaussian cube file to be displayed in VisIT or VMD */
void GaussianCube::writeGrid(ATOM_NETWORK *atmnet, std::string cubefilename, bool angstrom_to_bohr, bool useMassFlag){

 FILE *fp1;
 fp1=fopen(cubefilename.c_str(),"w");

 double toB=1; //scale factor - update it if we want to scale the output
 if(angstrom_to_bohr) toB=1.8903592 ; // scale factor to convert from A to Bohr

 fprintf(fp1,"\nThis is distance grid\n");

 fprintf(fp1,"%d % 13.6lf % 13.6lf % 13.6lf\n",atmnet->numAtoms,0.0,0.0,0.0);

 fprintf(fp1,"%d % 13.6lf % 13.6lf % 13.6lf\n",na,va.x*toB,va.y*toB,va.z*toB);
 fprintf(fp1,"%d % 13.6lf % 13.6lf % 13.6lf\n",nb,vb.x*toB,vb.y*toB,vb.z*toB);
 fprintf(fp1,"%d % 13.6lf % 13.6lf % 13.6lf\n",nc,vc.x*toB,vc.y*toB,vc.z*toB);

 for(int i=0;i<atmnet->numAtoms;i++)
   {
   if(useMassFlag == true)
     {
     fprintf(fp1,"%d % 13.6lf % 13.6lf % 13.6lf % 13.6lf\n",lookupAtomicNumber(atmnet->atoms[i].type),lookupMass(atmnet->atoms[i].type),
                                                              atmnet->atoms[i].x*toB, atmnet->atoms[i].y*toB, atmnet->atoms[i].z*toB);
     }else{
     fprintf(fp1,"%d % 13.6lf % 13.6lf % 13.6lf % 13.6lf\n", 1, 1.0,
                                                              atmnet->atoms[i].x*toB, atmnet->atoms[i].y*toB, atmnet->atoms[i].z*toB);
     };
   };

 fprintf(fp1," 1    1\n");

 int count2=0; // for formatting only
 for(int x=0;x<na;x++)
   for(int y=0;y<nb;y++)
     for(int z=0;z<nc;z++)
       {
       fprintf(fp1," % 13.6E ",cube[x][y][z]);

       // count2 is for formatting only (max 6 columns with data)
       count2++;
       if(z==(nc-1)) {fprintf(fp1,"\n");count2=0;};
       if(count2==6) {fprintf(fp1,"\n");count2=0;};
       };

} // ends writeGrid()





/* calculates distnace grid based on positions of atoms */
void GaussianCube::calculateDistanceGrid(ATOM_NETWORK *atmnet){

 for(int x=0; x<na; x++)
   for(int y=0; y<nb; y++)
     for(int z=0; z<nc; z++)
       {
       XYZ point(o.x+x*va.x+y*vb.x+z*vc.x,
                 o.y+x*va.y+y*vb.y+z*vc.y,
                 o.z+x*va.z+y*vb.z+z*vc.z);
       double dist=10000;

       for(int i=0;i<atmnet->numAtoms; i++){
           double atmDist = atmnet->calcDistanceXYZ(point.x, point.y, point.z, atmnet->atoms[i].x, atmnet->atoms[i].y, atmnet->atoms[i].z);
           double radius = atmnet->atoms[i].radius;
           atmDist -= radius;
           if(atmDist<dist) dist=atmDist;
           };


       cube[x][y][z]=dist;
       };

} // ends calculateDistanceGrid()




/* load a text file with points and project onto a grid to generate 3D histogram */
void GaussianCube::loadHistogramData(std::string inputfilename){

 fstream input;
 input.open(inputfilename.c_str());
 if(input.is_open()==false){
 cerr << "Error: CSSR failed to open " << inputfilename << endl;
 }else{

 // reading in datapoints (expect points in "Liverpool format" (fractional coordinates)

 int n = 0; // number of points 
 double a,b,c;
 int code;
 int pocketID;
 string ptStatus; 
 while(!input.eof())
    {
    input >> a;

    if(input.eof())
      {
      n--;
      break;
      };

    input >> b >> c >> code >> ptStatus >> pocketID;

    n++;

    // analysis of the read-in data
    a = trans_to_origuc(a);
    b = trans_to_origuc(b);
    c = trans_to_origuc(c);    

    // projecting into a grid

    int ca, cb, cc;

    ca = floor(a * na);
    cb = floor(b * nb);
    cc = floor(c * nc);

    cube[ca][cb][cc] += 1.0; // increase count at grid point ca, cb, cc

    };

 cout << n << " lines read." << "\n";

 input.close();
 };



} // ends loadHistogramData()


// load a list of text file with points(frames), and project each frame onto a grid to generate 3D histogram
// the histogram will represent number of frames that had a value >0 for particular grid point
void GaussianCube::loadHistogramDataPerFrame(std::string listfilename){

 fstream flist;

 flist.open(listfilename.c_str());
 if(flist.is_open()==false){
   cerr << "Error: A file with frames (" << listfilename << ") failed to open. \n" ;
   }else{
   // file with a list of frames is open. proceed.

   int nfiles=0;
   cout << "Loading filenames from " << listfilename << endl;
   while(!flist.eof())
    {
    string inputfilename;

    flist >> inputfilename;
    if(flist.eof()){
      nfiles--;
      break;
      };

    nfiles++;
 
    fstream input;
    input.open(inputfilename.c_str());
    if(input.is_open()==false){
      cerr << "Error: CSSR failed to open " << inputfilename << endl;
    }else{

    // reading in datapoints (expect points in "Liverpool format" (fractional coordinates)

    int n = 0; // number of points 
    double a,b,c;
    int code;
    string ptStatus; 
    int pocketID;
    while(!input.eof())
       {
       input >> a;

       if(input.eof())
         {
         n--;
         break;
         };

       input >> b >> c >> code >> ptStatus >> pocketID;

       n++;

       // analysis of the read-in data
       a = trans_to_origuc(a);
       b = trans_to_origuc(b);
       c = trans_to_origuc(c);    

       // projecting into a grid

       int ca, cb, cc;

       ca = floor(a * na);
       cb = floor(b * nb);
       cc = floor(c * nc);

//       cube[ca][cb][cc] += 1.0; // increase count at grid point ca, cb, cc

       if(cube[ca][cb][cc]-floor(cube[ca][cb][cc]) == 0.0) cube[ca][cb][cc] += 0.5; // adding 0.5 serves as a flag

       };

    cout << "File " << inputfilename << ":  " << n << " lines read." << "\n";

    input.close();
    };

   for(int i=0; i<na; i++)
      for(int j=0; j< nb; j++)
        for(int k=0; k<nc; k++)
           cube[i][j][k] = ceil(cube[i][j][k]);

   }; // ends while loop over a file list

  cout << nfiles << " frames loaded.\n";

  };

} // ends loadHistogramDataPerFrame()





/* calculates distnace grid based on positions of atoms 
   but also supplies with information if space is accessible to a probe
   nonaccessible volume is highlighted by negative distance
*/
void GaussianCube::calculateDistanceGridWithAccessibilityInfo(ATOM_NETWORK *atmnet, ATOM_NETWORK *orgatmnet, bool highAccuracy, double probe_radius){

 AccessibilityClass accessAnalysis;
 if(highAccuracy) accessAnalysis.setupAndFindChannels(atmnet, orgatmnet, highAccuracy, probe_radius, probe_radius);
   else accessAnalysis.setupAndFindChannels(atmnet, atmnet, highAccuracy, probe_radius, probe_radius);

 for(int x=0; x<na; x++)
   for(int y=0; y<nb; y++)
     for(int z=0; z<nc; z++)
       {
       bool overlaps = false;
       bool inside = false;
       Point point(o.x+x*va.x+y*vb.x+z*vc.x,
                 o.y+x*va.y+y*vb.y+z*vc.y,
                 o.z+x*va.z+y*vb.z+z*vc.z);
       double dist_var;
       pair<bool,bool> answer = (accessAnalysis.isVPointInsideAtomAndNotAccessible(point, dist_var));
       inside = answer.first; overlaps = answer.second;

       if(accessAnalysis.needToResample() == true)
          {
          cout << "Need to resample in grid calc. Abort." << "Contact the author" << endl;
          abort();
          };

       double dist;

//       dist = accessAnalysis.lastMinDist();
//       if(inside == false && overlaps == true) dist = dist * -1.0;
  
       if(inside == true)
         {
         dist = 0.0;
         } else
         {
         dist = accessAnalysis.lastMinDist()-probe_radius;
         if(overlaps == true) dist = dist * -1.0;
         };

       cube[x][y][z]=dist;
       };

} // ends calculateDistanceGrid()




