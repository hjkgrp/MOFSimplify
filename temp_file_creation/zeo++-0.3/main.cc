#include <vector>
#include <voro++.hh>

#include "network.h"
#include "networkio.h"
#include "networkinfo.h"
#include "networkanalysis.h"
#include "grid.h"
#include "channel.h"
#include "feature.h"
#include "holograms.h"
#include "instructions.h"
#include "ray.h"
#include "psd.h"
#include "area_and_volume.h"
#include "sphere_approx.h"
#include "poreinfo.h"
#include "zeojobs.h"
#include "arguments.h" 

using namespace std;
using namespace voro;

/** Processes the command line parameters and invokes the appropriate
    functions accordingly. */
int main(int argc, char * argv[]){
  bool error = false;
  if(argc == 1){
    PrintInstructions();
//    exit(0);
    error = true;
  }

  if(!error) {
    char *filename = argv[argc-1];  
    if(!checkInputFile(filename)) error=true;
    if(!error) {
      char *name = new char [256];
      char *extension = new char [256];
      parseFilename(filename, name, extension);
     
      string inputFile = argv[argc - 1];
      vector< vector<string> > commands = vector< vector<string> > ();
      vector<string> currentCommand = vector<string> ();

      int argcount = 1;
      while(argcount < argc - 1){
        if(argv[argcount][0] == '-'){
          if(currentCommand.size() != 0)
	    commands.push_back(currentCommand);
          currentCommand = vector<string> ();
        }
        
        currentCommand.push_back(argv[argcount]);
        argcount++;
      }
      if(currentCommand.size() != 0)
        commands.push_back(currentCommand);

      int numCommands = commands.size();
      ATOM_NETWORK atmnet; ATOM_NETWORK orgAtomnet;
      VORONOI_NETWORK vornet;
      vector< VOR_CELL> cells;
      vector< BASIC_VCELL> bvcells;
      bool radial = true, outputZvis = false, buildBasicCells = false, useMass = true, highAccuracy = false, allowAdjustCoordsAndCell = false, stripAtomNamesFlag = false;

      string AccSetting;

      //by default, we want to populate the table with radii and masses, even if later we will not reference them
      initializeRadTable();
      initializeMassTable();
      initializeCovRadTable();
      initializePT();

      //detect -e, which will direct the entire Zeo++ execution to a new route (within jobs/tasks classes)
      for(int i = 0; i < numCommands; i++){
        vector<string> command = commands[i];
        if(command[0].compare("-e") == 0){
          cout << "Entering the new execution route\n";
          zeoJob job;
          bool setup_error = job.setup(commands, name, extension); // setup run parameters for a job
          if(setup_error == false) job.execute();
          job.deinit();
        return 0;
        };
      }



      // Search for commands that affect Voronoi decomposition beforehand
      for(int i = 0; i < numCommands; i++){
        vector<string> command = commands[i];
        if(command[0].compare("-r") == 0){
          processRadialParameters(command);
          radial = true;
        }
        else if(command[0].compare("-nor") == 0){
          radial = false;
        }
        else if(command[0].compare("-mass") == 0){
          processMassParameters(command);
          useMass = true;
        }
        else if(command[0].compare("-nomass") == 0){
          useMass = false;
        }
        else if(commands[i][0].compare("-zvis") == 0){
          outputZvis = true;
        }
        else if(commands[i][0].compare("-sa") == 0 || commands[i][0].compare("-vol") == 0){
          buildBasicCells = true;
        }
        else if(commands[i][0].compare("-ha") == 0){
          AccSetting = processAccuracyParameters(command);
          highAccuracy = true;
        }
        else if(commands[i][0].compare("-stripatomnames") == 0){
          stripAtomNamesFlag = true;
          initializeStripAtomNameInternalFlag(true);
          cout << "Striping atom names expected during input file read-in\n";
        }
        else if(commands[i][0].compare("-allowAdjustCoordsAndCell") == 0){
          allowAdjustCoordsAndCell = true;
          atmnet.allowAdjustCoordsAndCellFlag = allowAdjustCoordsAndCell;
          orgAtomnet.allowAdjustCoordsAndCellFlag = allowAdjustCoordsAndCell;
        }
      }

      //Read input file information into the ATOM_NETWORK
      if(strcmp(extension, "cuc") == 0){
        if(!readCUCFile(filename, &atmnet, radial)) error=true;
      }
      else if(strcmp(extension, "arc") == 0){
        if(!readARCFile(filename, &atmnet, radial)) error=true;
      }
      else if(strcmp(extension, "cssr") == 0){
        if(!readCSSRFile(filename, &atmnet, radial)) error=true;
      }
      else if(strcmp(extension, "obcssr") == 0){
        if(!readOBCSSRFile(filename, &atmnet, radial)) error=true;
      }
      else if(strcmp(extension, "cif") == 0){
        if(!readCIFFile(filename, &atmnet, radial)) error=true;
      }
      else if(strcmp(extension, "car") == 0){
        if(!readCARFile(filename, &atmnet, radial)) error=true;
      }
      else if (strcmp(extension, "v1") == 0){
        if(!readV1File(filename, &atmnet, radial)) error=true;
      }
      else if (strcmp(extension, "dlp") == 0){
        if(!readDLPFile(filename, &atmnet, radial)) error=true;
      }
      else if (strcmp(extension, "pdb") == 0){
        if(!readPDBFile(filename, &atmnet, radial)) error=true;
      }
      //Should never be reached because of checkInputFile()
      else{
        cerr << "Invalid input file extension " << extension << "\n"
	     << "Exiting ..." << "\n";
  //      exit(1);
        error = true;
      }


      if(!error) {

        //  Strip atom names with non-atomic-name indexes
        if(stripAtomNamesFlag)
         {
         stripAtomNames(&atmnet);  
         cout << "Striping atom names\n";
         };
        // Read radii after reading atoms if -r  requested
        if(radial)
          loadRadii(&atmnet);

        initializeAtomicNumberTable();
        initializeAtomCharacterTable();

        loadMass(useMass, &atmnet); // assigns masses to atoms (assigns 0s in useMass flag = false)
        // High accuracy options - if the system contains particles of different radii, the employed radical
        // Voronoi decomposition may not be accurate. it can be, however, improved if large particles are 
        // replaced by clusters of smaller one (idea taken from Carolyn Philiphs/Chris Rycroft). 
        atmnet.copy(&orgAtomnet); // keep a copy of original atomnet (if high accuracy is not set, this is the same as analyzed network)
        if(radial && highAccuracy)
          {
          // atmnet.copy(&orgAtomnet); // keep a copy of original atomnet
          // calling the following function will modify atmnet - replace large atoms with clusters of small ones
          setupHighAccuracyAtomNetwork(&atmnet, AccSetting);
          };
        // please note that orgAtomnet is empty unless highAccuracy flag is set !! 
 

        // Declase object req. for Voronoi decomposition

        container_periodic_poly *rad_con = NULL;
        container_periodic   *no_rad_con = NULL;

        // Perform the relevant Voronoi decomposition
        cout << "Starting Voronoi decomposition" << "\n";
        if(radial) 
          rad_con = (container_periodic_poly *)performVoronoiDecomp(true, &atmnet, &vornet, cells, outputZvis, bvcells);
        else 
          no_rad_con = (container_periodic *)performVoronoiDecomp (false, &atmnet, &vornet, cells, outputZvis, bvcells); 
        cout << "Finished Voronoi decomposition" << "\n";


        // Execute each command
        for(int i = 0; i < numCommands && !error; i++){
          vector<string> command = commands[i];  cout << "Command " << i << "  " << command[0] << "\n";
          if(command[0].compare("-cssr") == 0){
            if(strcmp(extension, ".v1") == 0){
	      cerr << "Can not output a .cssr file when a .v1 file is used as input." << "\n"
	           << "Exiting ..." << "\n";
      //	exit(1);
              error = true; break;
            }
            string filename = processFilename(command, name, ".cssr", 0, 1);
            if(filename.empty()) {error=true; break;}
            if(!writeToCSSR((char *)filename.data(), &atmnet)) {error=true; break;}
          }
          else if(command[0].compare("-cif") == 0){
            string filename = processFilename(command, name, ".cif", 0, 1);
            if(filename.empty()) {error=true; break;}
            //if(!writeToCIF((char *)filename.data(), name, &atmnet)) {error=true; break;}
            if(!writeToCIF((char *)filename.data(), &atmnet)) {error=true; break;}
          }
          else if(command[0].compare("-xyz") == 0){
            string filename = processFilename(command, name, ".xyz", 0, 1);
            if(filename.empty()) {error=true; break;}
            bool is_supercell = false; bool is_duplicate_perimeter_atoms = false;
            if(!writeToXYZ((char *)filename.data(), &atmnet, is_supercell, 
                        is_duplicate_perimeter_atoms)) {
                error=true; break;
            }
          }
          else if(command[0].compare("-superxyz") == 0){
            string filename = processFilename(command, name, "_supercell.xyz", 0, 1);
            if(filename.empty()) {error=true; break;}
            bool is_supercell = true; bool is_duplicate_perimeter_atoms = false;
            if(!writeToXYZ((char *)filename.data(), &atmnet, is_supercell, 
                        is_duplicate_perimeter_atoms)) {
                error=true; break;
            }
          }
          else if(command[0].compare("-vis") == 0){
            string filename1 = processFilename(command, name, "_vis.xyz", 0, 1);
            if(filename1.empty()) {error=true; break;}
            bool is_supercell = false; bool is_duplicate_perimeter_atoms = true;
            if(!writeToXYZ((char *)filename1.data(), &atmnet, is_supercell, 
                        is_duplicate_perimeter_atoms)) {
                error=true; break;
            }
            string filename2 = processFilename(command, name, "_supercell_vis.xyz", 0, 1);
            if(filename2.empty()) {error=true; break;}
            is_supercell = true;
            if(!writeToXYZ((char *)filename2.data(), &atmnet, is_supercell, 
                        is_duplicate_perimeter_atoms)) {
                error=true; break;
            }
            string filename3 = processFilename(command, name, "_vis.vtk", 0, 1);
            if(filename3.empty()) {error=true; break;}
            if(!writeToVTK((char *)filename3.data(), &atmnet)) {error=true; break;}
          }
          else if(command[0].compare("-vtk") == 0){
            string filename = processFilename(command, name, ".vtk", 0, 1);
            if(filename.empty()) {error=true; break;}
            if(!writeToVTK((char *)filename.data(), &atmnet)) {error=true; break;}
          }
          else if(command[0].compare("-v1") == 0){
            string filename =  processFilename(command, name, ".v1", 0, 1);
            if(filename.empty()) {error=true; break;}
            if(!writeToV1((char *)filename.data(), &atmnet)) {error=true; break;}
          }
          else if(command[0].compare("-nt2") == 0){
            string filename =  processFilename(command, name, ".nt2", 0, 1);
            if(filename.empty()) {error=true; break;}
            if(!writeToNt2((char *)filename.data(), &vornet)) {error=true; break;}
          }
          else if(command[0].compare("-mopac") == 0){
            string filename =  processFilename(command, name, ".mop", 0, 1);
            if(filename.empty()) {error=true; break;}
            bool is_supercell = false;
            if(!writeToMOPAC((char *)filename.data(), &atmnet, is_supercell)) {error=true; break;}
          }
          else if(command[0].compare("-supermopac") == 0){
            string filename = processFilename(command, name, "_supercell.mop", 0, 1);
            if(filename.empty()) {error=true; break;}
            bool is_supercell = true;
            if(!writeToMOPAC((char *)filename.data(), &atmnet, is_supercell)) {error=true; break;}
          }
          else if(command[0].compare("-res") == 0){ // free sphere calcluation
            string filename = processFilename(command, name, ".res", 0, 1);
            if(filename.empty()) {error=true; break;}
            calculateFreeSphereParameters(&vornet, (char *)filename.data(), false);
          }
          else if(command[0].compare("-resex") == 0){ // free sphere calcluaton with extended output
            string filename = processFilename(command, name, ".res", 0, 1);
            if(filename.empty()) {error=true; break;}
            calculateFreeSphereParameters(&vornet, (char *)filename.data(), true);
          }
          else if(command[0].compare("-strinfo") == 0 || command[0].compare("-strinfoex") == 0){ 
                                               // print information about the structure/topology
                                               // extended option saves framework/molecules with their ids next to each atom
            string filename = processFilename(command, name, ".strinfo", 0, 1);
            string filenameExtOutput = processFilename(command, name, ".frameid", 0, 1); // filename to save ext. output
            if(filename.empty()) {error=true; break;}
            if(command[0].compare("-strinfo") == 0)
              getStructureInformation((char *)filename.data(), (char *)filenameExtOutput.data(), &atmnet, false); // last argument 4 ext. output
              else
              getStructureInformation((char *)filename.data(), (char *)filenameExtOutput.data(), &atmnet, true);
          }
          else if(command[0].compare("-oms") == 0 || command[0].compare("-omsex") == 0){
                                               // print information about open metal site present 
                                               // extended option saves detailed information for each atom 
            string filename = processFilename(command, name, ".oms", 0, 1);
            string filenameExtOutput = processFilename(command, name, ".omsex", 0, 1); // filename to save ext. output
            if(filename.empty()) {error=true; break;}
            if(command[0].compare("-oms") == 0)
              getOMSInformation((char *)filename.data(), (char *)filenameExtOutput.data(), &atmnet, false); // last argument 4 ext. output
              else
              getOMSInformation((char *)filename.data(), (char *)filenameExtOutput.data(), &atmnet, true);
          }
          else if(command[0].compare("-zvis") == 0){
            string filename = processFilename(command, name, ".zvis", 0, 1);
            if(filename.empty()) {error=true; break;}
            writeZeoVisFile((char *)filename.data(), &cells, &atmnet, &vornet);
          }
          else if(command[0].compare("-zvor") == 0){
            double r_probe  = strtod(command[1].data(), NULL);
            string output   = processFilename(command, name, ".zvor", 1, 2);
            if(output.empty()) {error=true; break;}
            viewVoronoiDecomp(&atmnet, r_probe, output);
          }
          else if(command[0].compare("-axs") == 0){
            string filename =  processFilename(command, name, ".axs", 1, 2);
            if(filename.empty()) {error=true; break;}
            vector<bool> accessInfo;
            vector<CHANNEL> channels;
            double probeRad = strtod(command[1].data(), NULL);
            CHANNEL::findChannels(&vornet, probeRad, &accessInfo, &channels);
            
      //      analyze_accessible_voronoi_pre_segment(&vornet, probeRad, &accessInfo, name);
            fstream output; output.open(filename.c_str(), fstream::out);
            for(unsigned int i = 0; i < accessInfo.size(); i++)
	      output << (accessInfo.at(i) ? "true" : "false") << "\n";
            output.close();
          }

      //Rich: holograms as a separate flag, -holo
          else if(command[0].compare("-holo") == 0){
            vector<bool> accessInfo;
            vector<CHANNEL> channels;
            if(command.size()<2 || command.size()>3) {
              printf("Error: -holo option accepts 1 (probe radius) or 2 (probe radius and then bin directory) arguments but %d arguments were supplied.\n", (int)(command.size() - 1));
              printf("Exiting...\n");
  //            exit(1);
              error=true; break;
            }
            double probeRad = strtod(command[1].data(), NULL);
            CHANNEL::findChannels(&vornet, probeRad, &accessInfo, &channels);
            if(command.size()==3) { //i.e. a directory was specified
              char *bin_directory = (char *)command[2].data();
              if(bin_directory[strlen(bin_directory)-1]!='/') {
                strcat(bin_directory, "/");
              }
              analyze_accessible_voronoi_pre_segment(&vornet, probeRad, &accessInfo, name, bin_directory); //this line generates and outputs Voronoi holograms with bins from a specified directory
            }
            analyze_accessible_voronoi_pre_segment(&vornet, probeRad, &accessInfo, name); //this line generates and outputs Voronoi holograms with default bins 
          }

      /* Rich: uses segments
          else if(command[0].compare("-simpl") == 0){    // This feature is added by Maciek
            string filename =  processFilename(command, name, ".seglist", 1, 2);
            if(filename.empty()) {error=true; break;}
            //  double chan_radius  = strtod(command[1].data(), NULL);
            double probe_radius = strtod(command[1].data(), NULL);
            vector<bool> accessInfo;
            vector<CHANNEL> channels;

            fstream output;
            output.open(filename.data(), fstream::out);
            output << calcDeterminant(atmnet.ucVectors) << "\n"
	           << calcDensity(&atmnet) << "\n";
            CHANNEL::findChannels(&vornet, probe_radius, &accessInfo, &channels);
            segmentChannels(&atmnet,&channels,output);
            output.close();

            // Added by Rich just to get a simple density and volume file
            fstream output2;
            string densvol = processFilename(command, name, ".densvol", 1, 2);
            if(densvol.empty()) {error=true; break;}
            output2.open(densvol.data(), fstream::out);
            output2 << name << " " << calcDeterminant(atmnet.ucVectors) << " " << calcDensity(&atmnet) << "\n";
            output2.close();
            }
      */
      /* Rich: uses segments
            else if(command[0].compare("-featholo") == 0){    // This feature is added by Maciek&Rich
                                                            // this should be called segmentholo
            string filename =  processFilename(command, name, ".featholo", 1, 2);
            if(filename.empty()) {error=true; break;}
            double probe_radius = strtod(command[1].data(), NULL);
            vector<bool> accessInfo;
            vector<CHANNEL> channels;
            fstream output;
      //    output.open(filename.data(), fstream::out);
            CHANNEL::findChannels(&vornet, probe_radius, &accessInfo, &channels);

            vector <int> vornode_segment_id_info;
            vector <double> segment_max_r;
            vornode_segment_id_info.resize(vornet.nodes.size(),0);
            segmentChannels_forHolograms(&atmnet,&channels,&vornode_segment_id_info,&segment_max_r,0); // 0 is for segments
            segmentChannels_saveSegments(&atmnet,&vornet, &vornode_segment_id_info,&segment_max_r,segment_max_r.size(),name,(char *)"segment");
      //Rich: if you want segment holograms use this function
	      //analyze_accessible_voronoi_with_segments(&vornet, probe_radius, &vornode_segment_id_info,&segment_max_r,segment_max_r.size(),name);
      //Rich: if you want atom cages, use this function
	      //analyze_accessible_voronoi_by_atoms(&atmnet, &vornet, probe_radius, &vornode_segment_id_info,&segment_max_r,segment_max_r.size(),name);
      //      output.close();
          }
      */
      /* Rich: uses segments
          else if(command[0].compare("-realfeatholo") == 0){    // This feature is added by Maciek
                                                                // it does the same as -featholo but for features not segments
            string filename =  processFilename(command, name, ".realfeatholo", 1, 2);
            if(filename.empty()) {error=true; break;}
            double probe_radius = strtod(command[1].data(), NULL);
            vector<bool> accessInfo;
            vector<CHANNEL> channels;
            fstream output;
      //    output.open(filename.data(), fstream::out);
            CHANNEL::findChannels(&vornet, probe_radius, &accessInfo, &channels);

            vector <int> vornode_feature_id_info;
            vector <double> feature_max_r;
            vornode_feature_id_info.resize(vornet.nodes.size(),0);
            segmentChannels_forHolograms(&atmnet,&channels,&vornode_feature_id_info,&feature_max_r,1); // 1 is for features
            segmentChannels_saveSegments(&atmnet,&vornet, &vornode_feature_id_info,&feature_max_r,feature_max_r.size(),name,(char *)"feature");
      //Rich: if you want segment holograms use this function
	      //analyze_accessible_voronoi_with_segments(&vornet, probe_radius, &vornode_segment_id_info,&segment_max_r,segment_max_r.size(),name);
      //Rich: if you want atom cages, use this function
	      //analyze_accessible_voronoi_by_atoms(&atmnet, &vornet, probe_radius, &vornode_segment_id_info,&segment_max_r,segment_max_r.size(),name);
      //      output.close();
      }
      */
          else if(command[0].compare("-chan") == 0){
            string filename =  processFilename(command, name, ".chan", 1, 2);
            if(filename.empty()) {error=true; break;}
            double probe_radius = strtod(command[1].data(), NULL);
            vector<bool> accessInfo;
            vector<CHANNEL> channels;
            CHANNEL::findChannels(&vornet, probe_radius, &accessInfo, &channels);

            fstream output;
            output.open(filename.data(), fstream::out);
            output << filename.data() << "   " << channels.size() << " channels identified of dimensionality ";
            for(unsigned int i = 0; i < channels.size(); i++){
	            output << channels[i].dimensionality << " ";
            }
	          output << "\n"; 
	          // Per channel analysis

            double di=0, df=0, dif=0;
	          for(unsigned int i=0; i < channels.size(); i++){
	            pair <double, pair<double,double> > maxdidfdif = channels[i].findFreeIncludedSphereDiameter();
              output << "Channel  " << i << "  "  << maxdidfdif.first << "  " << maxdidfdif.second.first << "  " << maxdidfdif.second.second << "\n";
              if(maxdidfdif.first>di) di = maxdidfdif.first;
              if(maxdidfdif.second.first>df) df = maxdidfdif.second.first;
              if(maxdidfdif.second.second>dif) dif = maxdidfdif.second.second;
	          }
            output << filename << " summary(Max_of_columns_above)   " << di << " " << df << "  " << dif << "  probe_rad: " << probe_radius << "  probe_diam: " << 2*probe_radius << "\n";

            output.close();
          }

          else if(command[0].compare("-zchan") == 0){
            string filename =  processFilename(command, name, ".zchan", 1, 2);
            if(filename.empty()) {error=true; break;}
            double probe_radius = strtod(command[1].data(), NULL);
            
            fstream output; 
            output.open(filename.data(), fstream::out);
            vector<bool> accessInfo;
            vector<CHANNEL> channels;
            CHANNEL::findChannels(&vornet, probe_radius, &accessInfo, &channels);
            for(unsigned int i = 0; i < channels.size(); i++){
	      channels.at(i).writeToVMD(i, output);
            }
            output << "set num_channels " << channels.size() << "\n";
            output.close();
          }
          // Extended printout on pores (channels and pockets)
          // This include di and positions of nodes
          else if(command[0].compare("-poreinfo") == 0){
            int numSamplesAV = 100;
            int numSamplesASA = 2500;
            string filename =  processFilename(command, name, ".poreinfo", 1, 2);
            if(filename.empty()) {error=true; break;}
            double probe_radius = strtod(command[1].data(), NULL);
            vector<bool> accessInfo;
            vector<PORE> pores;
            PORE::findChannelsAndPockets(&vornet, probe_radius, &accessInfo, &pores);

            int  n_channels = 0; int n_pockets = 0; // number of channels and pockets
            for(unsigned int i = 0; i < pores.size(); i++)
              {
              if(pores[i].dimensionality>0)
                { // Channels
                n_channels++;
                }else
                { // Pockets
                n_pockets++;
                };
              };

            fstream output;
            output.open(filename.data(), fstream::out);
            output << filename.data() << " Pores:  " << pores.size() << " Channels:  " << n_channels << " Pockets:  " << n_pockets << "\n";

            double volume = calcDeterminant(atmnet.ucVectors); // Unit cell volume/Units of A^3
            int numSamples = (int)(volume*numSamplesAV);
            calcAV(&atmnet, &orgAtomnet, highAccuracy, probe_radius, probe_radius, numSamples, true, output, (char *)filename.c_str(), 0, 0, 0, 0, -1,-1, 0);
            calcASA(&atmnet, &orgAtomnet, highAccuracy, probe_radius, probe_radius, calcDensity(&atmnet), numSamplesASA, true, output, (char *)filename.c_str(), 0, 0, 0, 0);

            for(unsigned int i = 0; i < pores.size(); i++)
              if(pores[i].dimensionality>0) pores[i].printPoreSummary(output, &atmnet); // Channels
           
            for(unsigned int i = 0; i < pores.size(); i++)
              if(pores[i].dimensionality==0) pores[i].printPoreSummary(output, &atmnet); // Pockets


/*
            for(unsigned int i = 0; i < channels.size(); i++){
              output << channels[i].dimensionality << " ";
              };
            output << "\n";
            // Per channel analysis
            for(unsigned int i=0; i < channels.size(); i++){
              pair <double, pair<double,double> > maxdidfdif = channels[i].findFreeIncludedSphereDiameter();
              output << "Channel  " << i << "  "  << maxdidfdif.first << "  " << maxdidfdif.second.first << "  " << maxdidfdif.second.second << "\n";
              };
*/
            output.close();
          }
          else if(command[0].compare("-poreinfoSummary")==0){  // analyzes a series of .poreinfo files  
            string filename_output =  processFilename(command, name, ".poreinfoSummary", 1, 2);
            if(filename_output.empty()) {error=true; break;}
            string filename_InputData = string(command[1]);
//            calculateAverageGrid(&atmnet, filename_InputData, filename_Gaussian_cube, angstrom_to_bohr, useMass);
            analyzePoreInfoFiles(&atmnet, filename_InputData, filename_output);
          }
          else if((command[0].compare("-sa") == 0) || (command[0].compare("-saex") == 0) || (command[0].compare("-zsa") == 0) || (command[0].compare("-vsa") == 0) || (command[0].compare("-lsa")==0)){
            bool visualize = (command[0].compare("-zsa") == 0 || command[0].compare("-vsa") == 0 || command[0].compare("-lsa")==0);
            bool visVisITflag = (command[0].compare("-vsa")==0 || command[0].compare("-lsa")==0);
            bool LiverpoolFlag = (command[0].compare("-lsa")==0);
            bool ExtendedOutputFlag = (command[0].compare("-saex") == 0);
            string suffix;
            if(visualize)
	      {
	      if(visVisITflag) {
                if(LiverpoolFlag) suffix = ".lsa"; else suffix = ".vsa";
                } else {suffix = ".zsa";};
	      }
            else
	      suffix = ".sa";
            string filename = processFilename(command, name, suffix, 3, 4);
            if(filename.empty()) {error=true; break;}
            double chan_radius  = strtod(command[1].data(), NULL);
            double probe_radius = strtod(command[2].data(), NULL);
            int numSamples      = int(strtod(command[3].data(), NULL));
          
            fstream output; 
            output.open(filename.data(), fstream::out);
            calcASA(&atmnet, &orgAtomnet, highAccuracy, chan_radius, probe_radius, calcDensity(&atmnet), numSamples, true, output, (char *)filename.data(), visualize, visVisITflag, LiverpoolFlag, ExtendedOutputFlag); 
            output.close();
          }
          else if( (command[0].compare("-vol") == 0) || (command[0].compare("-zvol") == 0) || (command[0].compare("-vvol") == 0) || (command[0].compare("-lvol") == 0) || (command[0].compare("-volpo") == 0)) {
            bool visualize = (command[0].compare("-zvol") == 0 || command[0].compare("-vvol") == 0 || command[0].compare("-lvol") == 0);
            bool VisITflag = (command[0].compare("-vvol") == 0 || command[0].compare("-lvol") == 0);
            bool LiverpoolFlag = (command[0].compare("-lvol") == 0);
            bool ProbeOccupiableFlag = (command[0].compare("-volpo") == 0);
            string suffix;
            if(visualize){
                if(VisITflag == true)
                  {
                  if(LiverpoolFlag == true) suffix = ".lvol"; else suffix = ".vvol";
                  } else
                  {
            	  suffix = ".zvol";
                  };
                }
            else
            	suffix = ".vol";
            if(ProbeOccupiableFlag == true) suffix +="po";
            if(command.size()<4 || command.size()>7) { //running with additional distance arguments tell you the fraction of accessible volume which is within this distance range away from any atom CENTRE
              printf("Error: -vol option accepts between 3 and 6 arguments 1 but %d arguments were supplied.\n", (int)(command.size() - 1));
              printf("Note:\t-vol chan_radius probe_radius num_samples [outputfile_vol]\n");
              printf("or:\t-vol chan_radius probe_radius num_samples low_distance_range high_distance_range [outputfile_vol]\n");
              printf("Exiting...\n");
  //            exit(1);
              error=true; break;
            }
            if(command.size()>5 && command.size()<8 && ProbeOccupiableFlag == true) {
              printf("Error: -volpo does not run with up to 4 arguments  but %d arguments were supplied.\n", (int)(command.size() - 1));
              printf("Exiting...\n");
              error=true; break;
            }
            double chan_radius  = strtod(command[1].data(), NULL);
            double probe_radius = strtod(command[2].data(), NULL);
            int numSamples      = int(strtod(command[3].data(), NULL));
            if(numSamples<0) {
              printf("ERROR: cannot call -vol flag with negative numSamples (arg was: %d)\n", numSamples);
  //            exit(EXIT_FAILURE);
              error=true; break;
            }
            double low_dist_range = -1;
            double high_dist_range = -1; //negative values indicate if we are not sampling only within a specific distance of the surface
            bool within_range =false;
            if(command.size()>5) { //specifying ranges
              within_range = true;
              suffix+="wr";
              low_dist_range  = strtod(command[4].data(), NULL);
              high_dist_range = strtod(command[5].data(), NULL);
              if(low_dist_range<0 || high_dist_range<0) {
                printf("ERROR: cannot call -vol flag with negative distance ranges (args were: %.3f %.3f)\n", low_dist_range, high_dist_range);
  //              exit(EXIT_FAILURE);
              error=true; break;
              }
            }
            string filename;
            if(command.size()==5 || command.size()==7) //specifying output filename
              filename = string(command[command.size()-1]);
            else filename = string(name).append(suffix);
            fstream output; 
            output.open(filename.data(), fstream::out);
            bool blockingMode = false;
            calcAV(&atmnet, &orgAtomnet, highAccuracy, chan_radius, probe_radius, numSamples, true, output, (char *)filename.data(), visualize, VisITflag, LiverpoolFlag, blockingMode, low_dist_range, high_dist_range, ProbeOccupiableFlag);
            output.close();
          }
          else if( (command[0].compare("-block") == 0) ) {
            string suffix;
            suffix = ".block";
            double probe_radius  = strtod(command[1].data(), NULL);
            int numSamples      = int(strtod(command[2].data(), NULL));
            if(numSamples<0) {
              printf("ERROR: cannot call -block flag with negative numSamples (arg was: %d)\n", numSamples);
  //            exit(EXIT_FAILURE);
              error=true; break;
            }
            double low_dist_range = -1;
            double high_dist_range = -1; //negative values indicate if we are not sampling only within a specific distance of the surface
            string filename;
            if(command.size()==4) //specifying output filename
              filename = string(command[command.size()-1]);
            else filename = string(name).append(suffix);

            double volume = calcDeterminant(atmnet.ucVectors); // Unit cell volume/Units of A^3
            numSamples = (int)(volume*numSamples);
            fstream output; 
            output.open(filename.data(), fstream::out);
            calcAV(&atmnet, &orgAtomnet, highAccuracy, probe_radius, probe_radius, numSamples, true, output, (char *)filename.data(), 0, 0, 0, 1, low_dist_range, high_dist_range, 0);
            output.close();
          }

          else if(command[0].compare("-gridGAI")==0 || command[0].compare("-gridGAIBohr")==0){  // saves distance grid in Gausian cube format
                                                                                                // supplementing it with accessibility info
            double probe_radius  = strtod(command[1].data(), NULL);
            string filename_Gaussian_cube =  processFilename(command, name, ".cube", 1, 2);
            bool angstrom_to_bohr = command[0].compare("-gridGAIBohr")==0;
            if(filename_Gaussian_cube.empty()) {error=true; break;}
            generateGaussianGridWithAccessibilityInfo(&atmnet, &orgAtomnet, highAccuracy, probe_radius, filename_Gaussian_cube, angstrom_to_bohr, useMass);
          }

          else if(command[0].compare("-gridBOV")==0){  // saves distance grid in BOV format
            string filename_f_dist =  processFilename(command, name, "_f.distances", 0, 1);
            string filename_f_bov =  processFilename(command, name, "_f.bov", 0, 0);
            string filename_g_dist =  processFilename(command, name, "_g.distances", 0, 1);  // g and h grids are for testing only
            string filename_g_bov =  processFilename(command, name, "_g.bov", 0, 1);         // the code is commented out
            string filename_h_dist =  processFilename(command, name, "_h.distances", 0, 1);
            string filename_h_bov =  processFilename(command, name, "_h.bov", 0, 1);
            if(filename_f_dist.empty() || filename_f_bov.empty() || filename_g_dist.empty() || filename_g_bov.empty() || filename_h_dist.empty() || filename_h_bov.empty()) {error=true; break;}
            generateBOVGrid(&atmnet, filename_f_dist, filename_g_dist, filename_h_dist, filename_f_bov, filename_g_bov, filename_h_bov);
          }
          else if(command[0].compare("-gridG")==0 || command[0].compare("-gridGBohr")==0){  // saves distance grid in Gausian cube format
            string filename_Gaussian_cube =  processFilename(command, name, ".cube", 0, 1);
            bool angstrom_to_bohr = command[0].compare("-gridGBohr")==0;
            if(filename_Gaussian_cube.empty()) {error=true; break;}
            generateGaussianGrid(&atmnet, filename_Gaussian_cube, angstrom_to_bohr, useMass);
          }
          else if(command[0].compare("-gridGprojdata")==0 || command[0].compare("-gridGprojdataBohr")==0){  // projects a file with datapoints (frac. coordinates) onto a gaussian grid (3d histogram) 
            string filename_Gaussian_cube =  processFilename(command, name, ".cube", 1, 2);
            bool angstrom_to_bohr = command[0].compare("-gridGprojdataBohr")==0;
            if(filename_Gaussian_cube.empty()) {error=true; break;}
            string filename_InputData = string(command[1]);
            // cout << "projdata " << filename_InputData << "  " << filename_Gaussian_cube << "\n";
            calculateAverageGrid(&atmnet, filename_InputData, filename_Gaussian_cube, angstrom_to_bohr, useMass);
          }
          else if(command[0].compare("-gridGprojdataperframe")==0 || command[0].compare("-gridGprojdataperframeBohr")==0){  // projects a file with datapoints (frac. coordinates) onto a gaussian grid (3d histogram)
                                                                     // it adds 1 per each frame if a point is occupied 
            string filename_Gaussian_cube =  processFilename(command, name, ".cube", 1, 2);
            bool angstrom_to_bohr = command[0].compare("-gridGprojdataperframeBohr")==0;
            if(filename_Gaussian_cube.empty()) {error=true; break;}
            string filename_InputData = string(command[1]);
            // cout << "projdata " << filename_InputData << "  " << filename_Gaussian_cube << "\n";
            calculateAverageGridPerFrame(&atmnet, filename_InputData, filename_Gaussian_cube, angstrom_to_bohr, useMass);
          }
          else if((command[0].compare("-psd") == 0) || (command[0].compare("-zpsd") == 0) || command[0].compare("-vpsd") == 0 ) {
            string suffix_histogram = ".psd_histo";
            string filename_histogram = processFilename(command, name, suffix_histogram, 3, 4);
            if(filename_histogram.empty()) {error=true; break;}
            double chan_radius  = strtod(command[1].data(), NULL);
            double probe_radius = strtod(command[2].data(), NULL);
            int numSamples = int(strtod(command[3].data(), NULL));
            //fstream output;
            //fstream outfile;
            //fstream nodesRadii;
            //fstream spheresDist;
            //output.open(filename_histogram.data(), fstream::out);
            bool visualize = (command[0].compare("-zpsd") == 0 || command[0].compare("-vpsd")==0);
                              bool visVisITflag = (command[0].compare("-vpsd")==0);
            if(visualize){  // -vpsd or -zpsd, visualization run
              string filename_points;
              string filename_radii;

                            if (visVisITflag){
                               string suffix_points = ".vpsdpts";
                               string suffix_radii = ".vpsdradii";
                               filename_points = processFilename(command, name, suffix_points, 3, 4);
                               filename_radii = processFilename(command, name, suffix_radii, 3, 4);
                                if(filename_points.empty() || filename_radii.empty()) {error=true; break;}
                               //outfile.open(filename_points.data(), fstream::out);
                               //nodesRadii.open(filename_radii.data(), fstream::out);
                               calcPoreSizeDistr(&atmnet, &orgAtomnet, highAccuracy, chan_radius, probe_radius, numSamples, true, filename_histogram, filename_points, filename_radii, "", visualize, visVisITflag);
                               }
                            if (!visVisITflag){
                               // here goes ZeoVis
                               }
            }
            else {  // -psd: regular PSD run, no visualization
              calcPoreSizeDistr(&atmnet, &orgAtomnet, highAccuracy,chan_radius, probe_radius, numSamples, true, filename_histogram, "", "", "", false, false);
            }
            //output.close();
            //outfile.close();
            //nodesRadii.close();
            //spheresDist.close();
            } // ends PSD
           else if((command[0].compare("-ray_atom") == 0) || (command[0].compare("-zray_atom") == 0)) {
            bool visualize = (command[0].compare("-zray_atom") == 0);
            string suffix;
            if(visualize)
	      suffix = ".zray_atom";
            else
	      suffix = ".ray_atom";
            string filename =  processFilename(command, name, suffix, 3, 4);
            if(filename.empty()) {error=true; break;}
            double chan_radius  = strtod(command[1].data(), NULL);
            double probe_radius = strtod(command[2].data(), NULL);
            int numSamples      = int(strtod(command[3].data(), NULL));
            string atom="atom";
            fstream output; 
            output.open(filename.data(), fstream::out);
            calcRaysInAV(&atmnet, &orgAtomnet, highAccuracy, chan_radius, probe_radius, numSamples, output, visualize,atom); 
            output.close();
            }

          else if((command[0].compare("-ray_sphere") == 0) || (command[0].compare("-zray_sphere") == 0)) {
            bool visualize = (command[0].compare("-zray_sphere") == 0);
            string suffix;
            if(visualize)
	      suffix = ".zray_sphere";
            else
	      suffix = ".ray_sphere";
            string filename =  processFilename(command, name, suffix, 3, 4);
            if(filename.empty()) {error=true; break;}
            double chan_radius  = strtod(command[1].data(), NULL);
            double probe_radius = strtod(command[2].data(), NULL);
            int numSamples      = int(strtod(command[3].data(), NULL));
            string sphere="sphere";
            fstream output; 
            output.open(filename.data(), fstream::out);
            calcRaysInAV(&atmnet, &orgAtomnet, highAccuracy, chan_radius, probe_radius, numSamples, output, visualize,sphere); 
            output.close();
            }

          else if((command[0].compare("-ray_node") == 0) || (command[0].compare("-zray_node") == 0)) {
            bool visualize = (command[0].compare("-zray_node") == 0);
            string suffix;
            if(visualize)
	      suffix = ".zray_node";
            else
	      suffix = ".ray_node";
            string filename =  processFilename(command, name, suffix, 3, 4);
            if(filename.empty()) {error=true; break;}
            double chan_radius  = strtod(command[1].data(), NULL);
            double probe_radius = strtod(command[2].data(), NULL);
            int numSamples      = int(strtod(command[3].data(), NULL));
            string node="node";
            fstream output; 
            output.open(filename.data(), fstream::out);
            calcRaysInAV(&atmnet, &orgAtomnet, highAccuracy, chan_radius, probe_radius, numSamples, output, visualize,node); 
            output.close();
            }
          else if((command[0].compare("-ray_andrew_sphere") == 0) || (command[0].compare("-zray_andrew_sphere") == 0)) {
            bool visualize = (command[0].compare("-zray_andrew_sphere") == 0);
            string suffix;
            if(visualize)
	      suffix = ".zray_andrew_sphere";
            else
	      suffix = ".ray_andrew_sphere";
            string filename =  processFilename(command, name, suffix, 3, 4);
            if(filename.empty()) {error=true; break;}
            double chan_radius  = strtod(command[1].data(), NULL);
            double probe_radius = strtod(command[2].data(), NULL);
            int numSamples      = int(strtod(command[3].data(), NULL));
            string option="andrew_sphere";
            fstream output; 
            output.open(filename.data(), fstream::out);
            calcRaysInAV(&atmnet, &orgAtomnet, highAccuracy, chan_radius, probe_radius, numSamples, output, visualize,option); 
            output.close();
          }
          else if((command[0].compare("-ray_andrew_atom") == 0) || (command[0].compare("-zray_andrew_atom") == 0)) {
            bool visualize = (command[0].compare("-zray_andrew_atom") == 0);
            string suffix;
            if(visualize)
	      suffix = ".zray_andrew_atom";
            else
	      suffix = ".ray_andrew_atom";
            string filename =  processFilename(command, name, suffix, 3, 4);
            if(filename.empty()) {error=true; break;}
            double chan_radius  = strtod(command[1].data(), NULL);
            double probe_radius = strtod(command[2].data(), NULL);
            int numSamples      = int(strtod(command[3].data(), NULL));
            string option="andrew_atom";
            fstream output; 
            output.open(filename.data(), fstream::out);
            calcRaysInAV(&atmnet, &orgAtomnet, highAccuracy, chan_radius, probe_radius, numSamples, output, visualize,option); 
            output.close();
          }
          else if(command[0].compare("-r") == 0 || command[0].compare("-nor") == 0){
            // Preprocessed radii
          }
          else if(command[0].compare("-ha") == 0 ){
            // Preprocessed high accuracy settings
          }
          else if(command[0].compare("-mass") == 0 || command[0].compare("-nomass") == 0){
            // Preprocessed atomic masses
          }
          else if(command[0].compare("-stripatomnames") == 0){
            // Preprocessed atomic masses
            //
          }
          else if((command[0].compare("-cage") == 0) || (command[0].compare("-zcage") == 0)){
            bool visualize = (command[0].compare("-zcage") == 0);
            string suffix;
            if(visualize)
	      suffix = ".zcage";
            else
	      suffix = ".cage";

            string filename =  processFilename(command, name, suffix, 1, 2);
            if(filename.empty()) {error=true; break;}
            double probe_radius = strtod(command[1].data(), NULL);
            
            DIJKSTRA_NETWORK dnet;
            DIJKSTRA_NETWORK::buildDijkstraNetwork(&vornet,&dnet);

            fstream output; vector<CAGE> cages;
            output.open(filename.data(), fstream::out);
            identifyCages(&atmnet, &vornet, &dnet, probe_radius, visualize, output, cages);
            output.close();
          }
          else if(command[0].compare("-zseg") == 0){
            string filename     = processFilename(command, name, ".zseg", 1, 2);
            if(filename.empty()) {error=true; break;}
            double probe_radius = strtod(command[1].data(), NULL);
            vector<bool> accessInfo;     vector<CHANNEL> channels;
            //VORONOI_NETWORK newNetwork;  pruneVoronoiNetwork(&vornet, &newNetwork, probe_radius);
            VORONOI_NETWORK newNetwork = vornet.prune(probe_radius);
            DIJKSTRA_NETWORK dnet;       DIJKSTRA_NETWORK::buildDijkstraNetwork(&newNetwork, &dnet);
            
            CHANNEL::findChannels(&dnet, &accessInfo, &channels);

            fstream output;  output.open(filename.data(), fstream::out);
            int initIndex = 0;
            vector<int> segCounts;
            for(unsigned int i = 0; i < channels.size(); i++){ 
	      vector<int> nodeIDs = vector<int>();
	      map<int,int>::iterator iter = channels[i].idMappings.begin();
	      while(iter != channels[i].idMappings.end()){
	        nodeIDs.push_back(iter->first);
	        iter++;
	      }

	      FEATURE seg = FEATURE(nodeIDs, &dnet, channels[i].dimensionality, channels[i].basis);
	      int numSegs = seg.createSegments(&atmnet,&newNetwork, &dnet, output, initIndex);
	      segCounts.push_back(numSegs);
	      initIndex += numSegs;
            }
            output.close();

            //fstream output;  
            output.open("seginfo.data", fstream::out);
            output << "set seg_counts {";
            for(unsigned int i = 0; i < segCounts.size(); i++)
	      output << segCounts[i] << " ";
            output << "}" << "\n";
            output.close();
            }
          else if(command[0].compare("-zfeat") == 0){
            string filename     = processFilename(command, name, ".zfeat", 1, 2);
            if(filename.empty()) {error=true; break;}
            double probe_radius = strtod(command[1].data(), NULL);
            vector<bool> accessInfo;     vector<CHANNEL> channels;
            //VORONOI_NETWORK newNetwork;  pruneVoronoiNetwork(&vornet, &newNetwork, probe_radius);
            VORONOI_NETWORK newNetwork = vornet.prune(probe_radius);
            DIJKSTRA_NETWORK dnet;       DIJKSTRA_NETWORK::buildDijkstraNetwork(&newNetwork, &dnet);

            CHANNEL::findChannels(&dnet, &accessInfo, &channels);
            
            fstream output;  output.open(filename.data(), fstream::out);
            int initIndex = 0;
            vector<int> featCounts;
            for(unsigned int i = 0; i < channels.size(); i++){
	      vector<int> nodeIDs = vector<int>();
	      map<int,int>::iterator iter = channels[i].idMappings.begin();
	      while(iter != channels[i].idMappings.end()){
	        nodeIDs.push_back(iter->first);
	        iter++;
	      }

	      FEATURE feat = FEATURE(nodeIDs, &dnet, channels[i].dimensionality, channels[i].basis);
	      int numFeats = feat.createFeatures(&atmnet,&newNetwork, &dnet, output, initIndex, name);
	      featCounts.push_back(numFeats);
	      initIndex += numFeats;
            }
            output.close();
            
            output.open("featinfo.data", fstream::out);
            output << "set feat_counts {";
            for(unsigned int i = 0; i < featCounts.size(); i++)
	      output << featCounts[i] << " ";
            output << "}" << "\n";
            output.close();
            }
          else if(command[0].compare("-sub") == 0){
            string filename = processFilename(command, name, "_sub.cssr", 0, 1);
            if(filename.empty()) {error=true; break;}
            ATOM_NETWORK newNetwork;
            newNetwork.allowAdjustCoordsAndCellFlag = allowAdjustCoordsAndCell;
            int numSubs = 0;
            bool subSeed = false; // Change to true to generate other possible configuration
            bool success = substituteAtoms(&atmnet, &newNetwork, subSeed, &numSubs, radial);
            if(!success){
	      cerr << "Error: Atom substitution unsucessful." << "\n"
	           << "Output not written to file." << "\n";
            }
            else {
	      writeToCSSR((char *)filename.data(), &newNetwork);
	      cout << "Substituted " << numSubs << " atoms when analyzing " << name << "\n";
            }
          }
          else if(command[0].compare("-fsub") == 0){
            string filename = processFilename(command, name, "_fsub.cssr", 1, 2);
            if(filename.empty()) {error=true; break;}
            double frac     = strtod(command[1].data(), NULL);
            ATOM_NETWORK newNetwork;
            newNetwork.allowAdjustCoordsAndCellFlag = allowAdjustCoordsAndCell;
            int numSubs = 0;
            double subFrac = 0.0;
            
            // Change either parameter to generate other possible configurations
            int  randSeed = 105543297; // Random number generator seed. 
            bool subSeed = false;      // Substitute first Si when creating initial 50/50 configuration 
            bool success = fracSubstituteAtoms(&atmnet, &newNetwork, subSeed, frac, randSeed, &numSubs, &subFrac, radial);
            if(!success){
	      cerr << "Error: Atom substitution unsucessful." << "\n"
	           << "Output not written to file." << "\n";
            }
            else {
	      writeToCSSR((char *)filename.data(), &newNetwork);
	      cout << "Substituted " << numSubs << " Si atoms for a fraction of " << subFrac << " when analyzing " << name << "\n";
            }
          }
          /* Maciek's version of substitution */
          /* This function has two modifications w.r.t -fsub. It allows chaging random seed to generate different configurations
             and removes topology check which ensured that every other atom is Si */ 
          /* usage: -fsubM fraction random_seed(0 for default) outname */
          else if(command[0].compare("-fsubM") == 0){
            string filename = processFilename(command, name, "_fsub.cssr", 2, 3);
            if(filename.empty()) {error=true; break;}
            double frac     = strtod(command[1].data(), NULL);
            int rsd         = int(strtod(command[2].data(), NULL));
            ATOM_NETWORK newNetwork;
            newNetwork.allowAdjustCoordsAndCellFlag = allowAdjustCoordsAndCell;
            int numSubs = 0;
            double subFrac = 0.0;

            // Change either parameter to generate other possible configurations
            int  randSeed = 105543297; // Random number generator seed. 
         
            if(rsd!=0) 
                {
                randSeed = rsd; // this is a quick fix to renerate random structures
                cout << "Overriding default random seed with " << rsd << endl;
                };
            bool subSeed = false;      // Substitute first Si when creating initial 50/50 configuration 
            bool success = fracSubstituteAtoms_Maciek(atmnet, newNetwork, subSeed, frac, randSeed, numSubs, subFrac, radial);
            if(!success){
              cerr << "Error: Atom substitution unsucessful." << "\n"
                   << "Output not written to file." << "\n";
            }
            else {
              writeToCSSR((char *)filename.data(), &newNetwork);
              cout << "Substituted " << numSubs << " Si atoms for a fraction of " << subFrac << " when analyzing " << name << "\n";
            }
          }
          //create skeleton diagram (prints Voronoi network to file, and accessible Voronoi network to a separate file, in an adjacent cell if specified with shift - xyz and vtk output files generated)
          else if(command[0].compare("-visVoro") == 0){
            if(command.size()!=2 && command.size()!=5) {
              printf("Error: -visVoro option accepts 1 (probe radius) or 4 (probe radius and then a, b and c shifts for illustrating accessible part of network) arguments but %d arguments were supplied.\n", (int)(command.size() - 1));
              printf("Exiting...\n");
  //            exit(1);
              error=true; break;
            }
            double probeRad = strtod(command[1].data(), NULL);
            int skel_a = 0, skel_b = 0, skel_c = 0;
            if(command.size()==5) {
              //if shift was provided
              skel_a = strtod(command[2].data(), NULL), skel_b = strtod(command[3].data(), NULL), skel_c = strtod(command[4].data(), NULL);
            }
            visVoro(name, probeRad, skel_a, skel_b, skel_c, &vornet, &atmnet);
          }
          //extract spherical substructures: this is a functionality designed for extracting local substructures of zeolites so that they can be scanned for potential guest molecule binding sites
          //this functionality writes out a number of xyz format files containing spherical substructures of the given radius, centred on given probe-accessible Voronoi nodes; if an element_type is given, a simplified Voronoi network is used, based only on atoms of that type
          else if(command[0].compare("-sphericalSubstructures") == 0) {
            if(command.size()!=3 && command.size()!=4) {
              printf("Error: -sphericalSubstructures option accepts 2 or 3 (probe_radius, sphere_radius, [element_type]) argument but %d arguments were supplied.\n", (int)(command.size() - 1));
              printf("Exiting...\n");
  //            exit(1);
              error=true; break;
            }
            bool usingElementOnly = false;
            string element("NONE");
            if(command.size()==4) {
              usingElementOnly = true;
              element = command[3];
      printf("DEBUG: extracting spheres from Voronoi network of %s atoms only\n", element.c_str());
            }
            double probeRad = strtod(command[1].data(), NULL);
            double sphereRad = strtod(command[2].data(), NULL);
            //utilise new cellmulti function; the smallest supercell satisfying sphereRad/2.0 is the same as the number of additional cells in each direction required to extract sphereRad around any point in the unit cell
            TRIPLET supercell = atmnet.getSmallestSupercell(sphereRad/2.0);
            getLocalSubstructures(name, probeRad, sphereRad, supercell, element, &vornet, &atmnet, radial, usingElementOnly);
          }

          //flag to print to terminal the tetrahedrality of, e.g. Zn4 groups inside a material
          else if(command[0].compare("-findTetrahedra") == 0) {
            if(command.size()!=2) {
              printf("Error: -findTetrahedra option accepts 1 (element type) argument but %d arguments were supplied.\n", (int)(command.size() - 1));
              printf("Exiting...\n");
  //            exit(1);
              error=true; break;
            }
            string element(command[1]);
            vector<double> tetras = atmnet.find_tetrahedra(element);
            //report
            printf("%s %d %s tetrahedra", name, (int)(tetras.size()), element.c_str());
            for(int i=0; i<(int)(tetras.size()); i++) printf(" %.3f", tetras.at(i));
            printf("\n");
          }

          //flag to write to terminal the size of supercell required in order to satisfy 'minimum image convention' (correct terminology?) - how many cells in each axis required so that no sphere of a given radius overlaps with itself periodically
          else if(command[0].compare("-cellmulti") == 0) {
            if(command.size()!=2) {
              printf("Error: -cellmulti option accepts 1 (sphere radius) argument but %d arguments were supplied.\n", (int)(command.size() - 1));
              printf("Exiting...\n");
  //            exit(1);
              error=true; break;
            }
            double rad = strtod(command[1].data(), NULL);
            double diam=2.0*rad;
            TRIPLET smallest_supercell = atmnet.getSmallestSupercell(diam);
            printf("%s supercell: %d %d %d = %d cells\n", name, smallest_supercell.x, smallest_supercell.y, smallest_supercell.z, smallest_supercell.x*smallest_supercell.y*smallest_supercell.z);
          }


else if(command[0].compare("-defining") == 0) {
  int num_v = vornet.nodes.size();
  for(int i=0; i<num_v; i++) {
    VOR_NODE v = vornet.nodes.at(i);
    int num_a = v.atomIDs.size();
    printf("node ID %d is defined by %d atoms:", i, num_a);
    bool hasMetal = false, onlyMetal = true;
    vector<string> elements;
    for(int j=0; j<num_a; j++) {
//      printf(" %d", v.atomIDs.at(j));
      string s = atmnet.atoms.at(v.atomIDs.at(j)).type;
      if(isMetal(s)) hasMetal = true; else onlyMetal = false;
      elements.push_back(s);
    }
    sort(elements.begin(), elements.end());
    int num_e = elements.size();
    for(int j=0; j<num_e; j++) {
      printf(" %s", elements.at(j).c_str());
    }
    if(hasMetal && onlyMetal) printf("\tMETAL"); else if(hasMetal) printf("\tMETAL-ORGANIC"); else printf("\tORGANIC");
    printf("\n");
  }
}


          else{
            cerr << "Error: Invalid option " << command[0] << " specified." << "\n"
	         << "Exiting..." << "\n"; 
      //      exit(EXIT_FAILURE);
            error=true; break; //memory leaks avoided by allowing main to terminate, rather than calling exit()
          }
        }
        //FREE MEMORY
        delete rad_con;
        delete no_rad_con;
      }
      //FREE MEMORY
      delete[] name;
      delete[] extension;
      return error ? error+1 : 0;
    }
  }
}


