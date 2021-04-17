
/* This function prints program instructions when Zeo++ is executed without parameters.
 */

//#include "network.h"
#include <iostream>
#include "instructions.h"

using namespace std;

void PrintInstructions() {
cout << "Network commandline invocation syntax:" << "\n" << "\n"
     << "./network [-cssr [outputfile_cssr]] " << "\n"
     <<	"          [-cif  [outputfile_cif]] " << "\n"
     <<	"          [-v1   [outputfile_v1]] " << "\n"
     << "          [-xyz  [outputfile_xyz]] " << "\n"
     << "          [-superxyz  [outputfile_xyz]] " << "\n"
     << "          [-vtk  [outputfile_vtk]] " << "\n"
     << "          [-vis  [outputfile_vtk/xyz]] " << "\n"
     << "          [-nt2  [outputfile_nt2]] " << "\n"
     << "          [-mopac  [outputfile_mop]] " << "\n"
     << "          [-supermopac  [outputfile_mop]] " << "\n"
     << "          [-res  [outputfile_res]] " << "\n"
     << "          [-zvis [outputfile_zvis]] " << "\n"
     <<	"          [-axs probe_radius [outputfile_axs]] " << "\n"
     <<	"          [-sa chan_radius probe_radius num_samples [outputfile_sa]] " << "\n"
     <<	"          [-vol chan_radius probe_radius num_samples [outputfile_vol]] " << "\n"
     << "          [-volpo chan_radius probe_radius num_samples [outputfile_vol]] " << "\n"
     << "          [-psd chan_radius probe_radius num_samples [outputfile_sa]] " << "\n"
     << "          [-ray chan_radius probe_radius num_samples [outputfile_ray]] " << "\n"
     <<	"          [-chan probe_radius [outputfile_chan]] " << "\n"
     <<	"          [-visVoro probe_radius [unit_cell_shifts: a b c]] " << "\n"
     <<	"          [-sphericalSubstructures probe_radius sphere_radius [element_type]]" << "\n"
     << "          [-findTetrahedra element_type]" << "\n" 
     << "          [-cellmulti sphere_radius]" << "\n" 
     <<	"          [-r [radfile]] " << "\n"
     << "          [-nor ] " << "\n"
     <<	"          INPUTFILE " << "\n" << "\n" << "\n";

 cout << " Note: For all of the program features that have an optional outputfile, the output is written to the filename when it is provided. Otherwise, " << "\n"
      << " the output is  written to ($file_prefix).($feature_name). For instance, if the provided file is ABW.cuc and -cssr is specified w/o any output file, " << "\n"
      << " the resulting input will be written by default to ABW.cssr " << "\n" << "\n"
      << " -cssr [outputfile_cssr] : Output a .cssr representation of the input file. " << "\n" << "\n"
      << " -cif  [outputfile_cif]  : Output a .cif representation of the input file. " << "\n" << "\n"
      << " -v1   [outputfile_v1]   : Output a .v1 representation of the input file. " << "\n" << "\n"
      << " -xyz  [outputfile_xyz]  : Output a .xyz representation of the input file. " << "\n" << "\n"
      << " -superxyz  [outputfile_xyz]  : Output a .xyz representation of a supercell of the input file. " << "\n" << "\n"
      << " -vtk  [outputfile_vtk]  : Output a .vtk representation of the unit cell of the input file. " << "\n" << "\n"
      << " -vis  [outputfile_vtk/xyz]  : For visualization: Output .xyz representations of the input file and a supercell of the input file, with duplicated perimeter atoms, and a .vtk representation of the unit cell of the input file. " << "\n" << "\n"
      << " -nt2  [outputfile_nt2]  : Output a .nt2 representation of the Voronoi network for the given unit cell. " << "\n" << "\n"
      << " -mopac  [outputfile_mop]  : Output a .mop representation of the input file. " << "\n" << "\n"
      << " -supermopac  [outputfile_mop]  : Output a .mop representation of a supercell of the input file. " << "\n" << "\n"
      << " -r    [radfile] : Use radii while performing the Vornoi decomposition (default behavior). If no .rad file is provided, " << "\n"
      << "                   default radii contained within the program are used, if possible. If a .rad " << "\n"
      << "                   file is provided, the default radii are overriden. The .rad file consists of" << "\n"
      << "                   two columns, containing the atomic types and radii respectively. For an example " << "\n"
      << "                   file, please refer to example.rad. " << "\n"
      << " -nor             : Do not use radii while performing the Vornoi decomposition. " << "\n" << "\n"
      << " -mass massfile :  Replace the built-in mass table with that provided in a .mass file. The .mass file consists" << "\n"
      << "                   of two columns, containing the atomic types and masses respectively. For an example file, " << "\n"
      << "                   please refer to example.mass " << "\n" << "\n";

 cout << " -res  [outputfile_res] : Output the maximum free sphere results to the provided file. " << "\n"
      << "                          Format is: maxIncDiam maxDiam incDiam " << "\n"    
      << "                          where  maxIncDiam is the maximum included sphere diamter, maxDiam is the maximum free sphere diameter" << "\n"
      << "	                    and incDiam is the maximum included diameters for the maximum free sphere path. "  << "\n" << "\n"
      << " -resex  [outputfile_res] : Extended -res printout. Additionally prints free and included sphere along each of a,b, and c axis " << "\n" << "\n"
      << " -zvis [outputfile_zvis] : Output a file used by ZeoVis to depict important features about the atomic network and Voronoi network. " << "\n" << "\n"
      << " -axs probe_radius [outputfile_axs] : Outputs a file with N lines, where N = # of Voronoi nodes. " << "\n" 
      << "	  		              Each line is true or false, representing the accessibility of the ith Voronoi node. " << "\n" << "\n";

 cout << " -sa chan_radius probe_radius num_samples_per_atom  [outputfile_sa]  : Determine the amount of accessible surface area. The accessibility of nodes is determined " << "\n"
      << "                                   	                               using $chan_radius, and probe-atom overlap is determiend using $probe_radius. A total of " << "\n"
      << "	                                                               $num_samples_per_atom points are sampled at a distance of $probe_radius from the atoms and Monte " << "\n"
      << "	                                                               Carlo integration is used to determine the final result. " << "\n" << "\n";

 cout << " -vol chan_radius probe_radius num_samples_total [outputfile_vol]  : Determine the amount of accessible volume. The accessibility of nodes is determined " << "\n"
      << "                                                                    using $chan_radius, and probe-atom overlap is determined using $probe_radius. A total of " << "\n"
      << "	                                                            $num_samples_total points are sampled across the unit cell and Monte Carlo integration is used " << "\n"
      << "	                                                            to determine the final result." <<  "\n" << "\n";

 cout << " -volpo chan_radius probe_radius num_samples_total [outputfile_vol]  : Similar to -vol command but return probe-occupiable values. \n\n";

 cout << " -psd chan_radius probe_radius num_samples_total [outputfile_psd]  : Calculate pore size distribution histogram. The accessibility of nodes is determined " << "\n"
      << "                                                                    using $chan_radius, and probe-atom overlap is determined using $probe_radius. A total of " << "\n"
      << "                                                                  $num_samples_total points are sampled across the unit cell and Monte Carlo integration is used " << "\n"
      << "                                                                  to determine the final result." <<  "\n" << "\n";

 cout << " -block probe_radius num_samples_total                       : Calculate blocking spheres to be used in MC simluations. The accessibility of nodes is determined " << "\n"
      << "                                                                   using $probe_radius. A total of $num_samples_total points are sampled across " << "\n"
      << "                                                                   the unit cell and Monte Carlo integration is used to determine the final result." <<  "\n" << "\n";

 cout << " -ray_atom or -zray_atom chan_radius probe_radius num_samples_total [outputfile_ray] : \n Shoot lines through accesible volume of cell until hits an atom and returns outputfile.r_atom with histogram of ray lengths (can be manipulated in source) also full information of ray lengths can be output in a file Ray_Info.txt in source" << "\n\n";
 cout << " -ray_sphere or -zray_sphere chan_radius probe_radius num_samples_total [outputfile_ray] : \n Shoot lines through voronoi nodes until it reaches a point outside of vornode and returns outputfile.r_atom with hitsogram of ray lengths (can be manipulated in source) also full information of ray lengths can be output in a file Ray_Info.txt in source" << "\n\n";
 cout << " -ray_andrew_sphere or -zray_andrew_sphere chan_radius probe_radius num_samples_total [outputfile_ray] : \n Randomly picks a point within the unitcell and shoots a ray until it reaches MAXRAYDIST and returns outputfile.r_atom with hitsogram of ray segments within vornodes (can be manipulated in source) also full information of ray lengths can be output in a file Ray_Info.txt in source" << "\n\n";
 cout << " -ray_andrew_atom or -zray_andrew_atom chan_radius probe_radius num_samples_total [outputfile_ray] : \n Randomly picks a point within the unitcell and shoots a ray until it reaches MAXRAYDIST and returns outputfile.r_atom with hitsogram of ray sements within accessible volume (can be manipulated in source) also full information of ray lengths can be output in a file Ray_Info.txt in source" << "\n\n";
 cout << " -sub [outputfile_sub] : Attempt to replace every other Si atom within the network with an Al atom consistently. " << "\n"
      << "                         If successful, write the resulting network to a .cssr file. No output is written if the substitution fails." << "\n" << "\n"
      << " -fsub frac [outputfile_fsub] : Attempt to randomly replace $frac of Si atoms within the network with an Al atom consistently. The specified fraction " << "\n"
      << "                                must not exceed 0.5. If successful, write the resulting network to a .cssr file. No output is written if the substitution fails." << "\n" << "\n"
      << " -fsubM frac seed [outputfile_fsub] : Attempt to randomly replace $frac of Si atoms within the network with Al atoms. The specified fraction " << "\n"
      << "                                must not exceed 0.5. Random seed specified by seed (use 0 fo default). If successful, write the resulting network to a .cssr file." << "\n" 
      << "                                This function uses different algorithms than -sub and -fsub, and may have problems reaching high $frac (email authors for more details). No output is written if the substitution fails." << "\n" << "\n"
      << " -gridBOV : Write distance grid in BOV format (inputfilename_f.bov). " << "\n"
      << " -gridG  : Write distance grid in Gaussian cube format (inputfilename.cube). " << "\n"
      << " -gridGBohr  : Write distance grid in Gaussian cube format, with distances converted from Angstrom to Bohr (inputfilename.cube). " << "\n"
      << " -chan probe_radius [outputfile_chan] : Determine the channels available to a probe of size $probe_radius. " << "\n" 
      << " -visVoro probe_radius [unit_cell_shifts: a b c] : For visualization: draw the Voronoi network as xyz and vtk formats, including only the accessible portion, which can be drawn in an adjacent unit cell with the provided shifts. " << "\n" 
      << " -sphericalSubstructures probe_radius sphere_radius [element_type] : For zeolites: Write out a number of xyz format files containing spherical substructures of the given radius, centred on given probe-accessible Voronoi nodes; if an element_type is given, a simplified Voronoi network is used, based only on atoms of that type. " << "\n" 
      << " -findTetrahedra element_type: Write to terminal the tetrahedrality (distortion) index for each tetrahedral arrangement of the specified atom type in the material" << "\n" 
      << " -cellmulti sphere_radius: Write to terminal the number of cells required in each crystallographic axis in order to construct a supercell, where a sphere with the provided radius will not overlap with itself periodically" << "\n" 
      << "		                        Outputs a file used by ZeoVis to draw the resulting channels. " << "\n" << "\n"
      << " INPUTFILE    : Must be in the .cssr, .cif, .cuc, .car, .arc or .v1 file format. " << "\n" << "\n" << "\n";
}
