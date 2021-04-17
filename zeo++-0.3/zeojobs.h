/** 
  * In zeojob.* files, there are functions to control execution of Zeo++ jobs
  *
  * Initial code by M. Haranczyk, Feb 2014
  *
  **/

#ifndef ZEOJOBS_H
#define ZEOJOBS_H

#include <cstdio>
#include <string>
#include "general.h"
#include "zeo_consts.h"
#include "geometry.h"
#include "networkstorage.h"
#include "channel.h"
#include "arguments.h"
#include "material.h"


/* zeoJob class is used to control execution of jobs in Zeo++
   it handles task management, structure data and output

*/

class zeoJob {
  /* Status flags */
  bool setupErrorFlag;

  /* Argument list  */
  std::vector< std::vector<std::string> > commands;
  char name[256];
  char extension[256];

  /* Run setting */
  double probeRadius;  // this stores a probe radius used in all Zeo++ analysis
  bool probeRadiusFlag;// set to true is "-probe" is set as argument

  bool volumeRunFlag; // sets to true if "-vol" is in argument list

  double nMCsamplesScaleFactor; // scaling factor for the default number of samples
  bool nMCsamplesScaleFactorFlag; // set to true if "-sqf" is used
  int nMCAVsampleDensity; // density of samples used in AV
  int nMCASAsampleDensity;//                and ASA cals

  bool extendedOutputFlag; // set to true if extended output flag is set

  public:

  MATERIAL Material;

  /* Constructor */
  zeoJob()
   {
   setupErrorFlag = false;
   probeRadiusFlag = false;
   volumeRunFlag = false;
   extendedOutputFlag = false;
   nMCsamplesScaleFactorFlag = false;
   nMCsamplesScaleFactor = 1.0;
   nMCAVsampleDensity = MCAVppA3; // use defauls from zeo_consts.h
   nMCASAsampleDensity = MCASAppA2; 
   };

  /* Deconstructor */
  void deinit()
   {
   Material.deinit();
   };


  bool setup(std::vector< std::vector<std::string> > commands_string, const char *name_string, const char*extension_string)
   {
   commands = commands_string;
   strcpy(name, name_string);
   strcpy(extension, extension_string);

   /* Initial setup */
   bool error = false;
   int numCommands = commands.size();

   /* Display welcome message */
   cout << "Entering new Zeo++ execution route. Input parameters are:\n";
   for(unsigned int i=0;i<commands.size();i++)
     for(unsigned int j=0;j<commands[i].size();j++)
       cout << commands[i][j] << " ";
   cout << name << "   " << extension << "\n";
   /* Ends display message */


   /* Pre-process the argument list:
    * Search for commands that affect Voronoi decomposition beforehand */  
   for(int i = 0; i < numCommands; i++){
     vector<string> command = commands[i];
     if(command[0].compare("-r") == 0){
       processRadialParameters(command);
       Material.radial = true;
     }
     else if(command[0].compare("-nor") == 0){
       //Material.radial = false;
       cerr << "-nor is disabled. Set radii to zero instead.\n";
     }
     else if(command[0].compare("-mass") == 0){
       processMassParameters(command);
       Material.useMass = true;
     }
     else if(command[0].compare("-nomass") == 0){
       Material.useMass = false;
     }
     else if(commands[i][0].compare("-zvis") == 0){
       Material.outputZvis = true;
     }
     else if(commands[i][0].compare("-sa") == 0 || commands[i][0].compare("-vol") == 0){
       Material.buildBasicCells = true;
       if(commands[i][0].compare("-vol") == 0) volumeRunFlag = true;
     }
     else if(commands[i][0].compare("-noha") == 0){
       Material.highAccuracy = false;
     }
     else if(commands[i][0].compare("-ha") == 0){
       Material.AccSetting = processAccuracyParameters(command);
       Material.highAccuracy = true;
     }
     else if(commands[i][0].compare("-ext") == 0 || commands[i][0].compare("-eo") == 0 || commands[i][0].compare("-extendedOutput") == 0){
       extendedOutputFlag = true;
     }
     else if(commands[i][0].compare("-vo") == 0 || commands[i][0].compare("-visual") == 0 || commands[i][0].compare("-visout") == 0){
       Material.VisSetting = processVisualizationParameters(command);
       Material.VisFlag = true;
     }
     else if(commands[i][0].compare("-block") == 0){
       Material.AVrequestBlockingPockets = true; 
     }
     else if(commands[i][0].compare("-psd") == 0){
       Material.AVrequestPSD = true;
     }
     else if(commands[i][0].compare("-allowAdjustCoordsAndCell") == 0){
       Material.allowAdjustCoordsAndCell = true;
       Material.atmnet.allowAdjustCoordsAndCellFlag = Material.allowAdjustCoordsAndCell;
       Material.orgAtomnet.allowAdjustCoordsAndCellFlag = Material.allowAdjustCoordsAndCell;
     }
     else if(commands[i][0].compare("-probe") == 0){
       probeRadius = strtod(commands[i][1].data(), NULL);
       probeRadiusFlag = true;
     }
     else if(commands[i][0].compare("-sqf") == 0){
       nMCsamplesScaleFactor = strtod(commands[i][1].data(), NULL);
       nMCsamplesScaleFactorFlag = true;
       nMCAVsampleDensity = (int)(nMCAVsampleDensity * nMCsamplesScaleFactor); 
       nMCASAsampleDensity = (int)(nMCASAsampleDensity * nMCsamplesScaleFactor);
       if(nMCAVsampleDensity < 1 || nMCASAsampleDensity < 1)
         {
         error = true;
         cerr << "Density of MC samples below 1. Please use larger value with -sqf argument. Current density for AV = " 
              << nMCAVsampleDensity << " and ASA = " << nMCASAsampleDensity <<".\n";
         };
     }
   }; // ends loop over all commands 

   /* Pre-process the argument list to make sure all arguemtns are specified */
   if(probeRadiusFlag == false){
     for(int i = 0; i < numCommands; i++){
       vector<string> command = commands[i];
       if(command[0].compare("-zvor") == 0 ||
          command[0].compare("-axs") == 0 ||
          command[0].compare("-holo") == 0 ||
          command[0].compare("-chan") == 0 ||
          command[0].compare("-zchan") == 0 ||
          command[0].compare("-poreinfo") == 0 ||
          command[0].compare("-test") == 0 ||
          command[0].compare("-vol") == 0 ||
          command[0].compare("-sa") == 0 ||
          command[0].compare("-visVoro") == 0 || 
          command[0].compare("-pld") == 0 ||
          command[0].compare("-pldf") == 0||
          command[0].compare("-pldmol") == 0)
       {
       error = true;
       cerr << command[0] << " option requires probe radius. Please use \"-probe R\" (where R is radius)\n";
       }
     }
   }

   /* Pre-process the arugment list to make sure all run dependiences are fulfilled */
   if(volumeRunFlag == false){
     for(int i = 0; i < numCommands; i++){
       vector<string> command = commands[i];
       if(command[0].compare("-block") == 0 || command[0].compare("-psd") == 0)
         {
         error = true;
         cerr << command[0] << " option depends on accessible volume calculation. Please add \"-vol nS\" (where nS is the number of MC samples per cubic Ang).\n";
         };
     };
   }   

   /* Load structure */
   //Read input file information into the ATOM_NETWORK
   if(error == false){
   char filename[256];
   strcpy(filename,name); strcat(filename,"."); strcat(filename,extension);
   cout << filename << endl;
   if(strcmp(extension, "cuc") == 0){
     if(!readCUCFile(filename, &(Material.atmnet), Material.radial)) error=true;
   }
   else if(strcmp(extension, "arc") == 0){
     if(!readARCFile(filename, &(Material.atmnet), Material.radial)) error=true;
   }
   else if(strcmp(extension, "cssr") == 0){
     if(!readCSSRFile(filename, &(Material.atmnet), Material.radial)) error=true;
   }
   else if(strcmp(extension, "cif") == 0){
     if(!readCIFFile(filename, &(Material.atmnet), Material.radial)) error=true;
   }
   else if(strcmp(extension, "car") == 0){
     if(!readCARFile(filename, &(Material.atmnet), Material.radial)) error=true;
   }
   else if (strcmp(extension, "v1") == 0){
     if(!readV1File(filename, &(Material.atmnet), Material.radial)) error=true;
   }
   else if(strcmp(extension, "pdb") == 0){
     if(!readPDBFile(filename, &(Material.atmnet), Material.radial)) error=true;
   }
   else if (strcmp(extension, "dlp") == 0){
     if(!readDLPFile(filename, &(Material.atmnet), Material.radial)) error=true;
   }
   //Should never be reached because of checkInputFile()
   else{
     cerr << "Invalid input file extension " << extension << "\n"
          << "Exiting ..." << "\n";
     error = true;
   };
   };

   /* Initialize radii, mass and accuracy setting */
   if(!error) {

     // Read radii after reading atoms if -r  requested
     if(Material.radial)
       loadRadii(&(Material.atmnet));

     loadMass(Material.useMass, &(Material.atmnet)); // assigns masses to atoms (assigns 0s in useMass flag = false)

     // High accuracy options - if the system contains particles of different radii, the employed radical
     // Voronoi decomposition may not be accurate. it can be, however, improved if large particles are 
     // replaced by clusters of smaller one (idea taken from Carolyn Philiphs/Chris Rycroft). 
     Material.atmnet.copy(&(Material.orgAtomnet)); // keep a copy of original atomnet 
                                                   //(if high accuracy is not set, this is the same as analyzed network)
     if(Material.radial && Material.highAccuracy)
       {
       // atmnet.copy(&orgAtomnet); // keep a copy of original atomnet
       // calling the following function will modify atmnet - replace large atoms with clusters of small ones
       setupHighAccuracyAtomNetwork(&(Material.atmnet), Material.AccSetting);
       };
     // please note that orgAtomnet is empty unless highAccuracy flag is set !! 

   }; // ends if(!error)

   setupErrorFlag = error;
   return error;

   }; // ends setup()


   void execute(){
    /* Status variables */
    bool error = false;
    int numCommands = commands.size();

    if(setupErrorFlag == true)
      {
      cerr << "ZEOJOB.execute() errors due to setupErrorFlag\n";
      error = true;
      };


    /* *************************************************************
     *    NOTES for execution()
     *
     *
     *  -ha option becomes default, as a consequence Material.atmnet 
     *             becomes modified (sphere approx). Please use Material.orgAtomnet
     *             to access original atom network (as read from the input file)
     *
     *  -Voronoi decomposition is now run only per request, please make sure to 
     *             execute Material.runVoroFlat()
     *
     */


    // Execute each command
    for(int i = 0; i < numCommands && !error; i++){
      vector<string> command = commands[i];  cout << "Command " << i << "  " << command[0] << "\n";
      if(command[0].compare("-cssr") == 0){
        if(strcmp(extension, ".v1") == 0){
          cerr << "Can not output a .cssr file when a .v1 file is used as input." << "\n"
               << "Exiting ..." << "\n";
  //        exit(1);
          error = true; break;
        }
        string filename = processFilename(command, name, ".cssr", 0, 1);
        if(filename.empty()) {error=true; break;}
        if(!writeToCSSR((char *)filename.data(), &(Material.orgAtomnet))) {error=true; break;}
      }
      else if(command[0].compare("-cif") == 0){
        string filename = processFilename(command, name, ".cif", 0, 1);
        if(filename.empty()) {error=true; break;}
        //if(!writeToCIF((char *)filename.data(), name, &atmnet)) {error=true; break;}
        if(!writeToCIF((char *)filename.data(), &(Material.orgAtomnet))) {error=true; break;}
      }
      else if(command[0].compare("-xyz") == 0){
        string filename = processFilename(command, name, ".xyz", 0, 1);
        if(filename.empty()) {error=true; break;}
        bool is_supercell = false; bool is_duplicate_perimeter_atoms = false;
        if(!writeToXYZ((char *)filename.data(), &(Material.orgAtomnet), is_supercell,
                    is_duplicate_perimeter_atoms)) {
            error=true; break;
        }
      }
      else if(command[0].compare("-superxyz") == 0){
        string filename = processFilename(command, name, "_supercell.xyz", 0, 1);
        if(filename.empty()) {error=true; break;}
        bool is_supercell = true; bool is_duplicate_perimeter_atoms = false;
        if(!writeToXYZ((char *)filename.data(), &(Material.orgAtomnet), is_supercell,
                    is_duplicate_perimeter_atoms)) {
            error=true; break;
        }
      }
      else if(command[0].compare("-vis") == 0){
        string filename1 = processFilename(command, name, "_vis.xyz", 0, 1);
        if(filename1.empty()) {error=true; break;}
        bool is_supercell = false; bool is_duplicate_perimeter_atoms = true;
        if(!writeToXYZ((char *)filename1.data(), &(Material.orgAtomnet), is_supercell,
                    is_duplicate_perimeter_atoms)) {
            error=true; break;
        }
        string filename2 = processFilename(command, name, "_supercell_vis.xyz", 0, 1);
        if(filename2.empty()) {error=true; break;}
        is_supercell = true;
        if(!writeToXYZ((char *)filename2.data(), &(Material.orgAtomnet), is_supercell,
                    is_duplicate_perimeter_atoms)) {
            error=true; break;
        }
        string filename3 = processFilename(command, name, "_vis.vtk", 0, 1);
        if(filename3.empty()) {error=true; break;}
        if(!writeToVTK((char *)filename3.data(), &(Material.orgAtomnet))) {error=true; break;}
      }
      else if(command[0].compare("-vtk") == 0){
        string filename = processFilename(command, name, ".vtk", 0, 1);
        if(filename.empty()) {error=true; break;}
        if(!writeToVTK((char *)filename.data(), &(Material.orgAtomnet))) {error=true; break;}
      }
      else if(command[0].compare("-v1") == 0){
        string filename =  processFilename(command, name, ".v1", 0, 1);
        if(filename.empty()) {error=true; break;}
        if(!writeToV1((char *)filename.data(), &(Material.orgAtomnet))) {error=true; break;}
      }
      else if(command[0].compare("-nt2") == 0){
        string filename =  processFilename(command, name, ".nt2", 0, 1);
        if(filename.empty()) {error=true; break;}
        Material.runVoroFlat();
        if(!writeToNt2((char *)filename.data(), &(Material.vornet))) {error=true; break;}
      }
      else if(command[0].compare("-mopac") == 0){
        string filename =  processFilename(command, name, ".mop", 0, 1);
        if(filename.empty()) {error=true; break;}
        bool is_supercell = false;
        if(!writeToMOPAC((char *)filename.data(), &(Material.orgAtomnet), is_supercell)) {error=true; break;}
      }
      else if(command[0].compare("-supermopac") == 0){
        string filename = processFilename(command, name, "_supercell.mop", 0, 1);
        if(filename.empty()) {error=true; break;}
        bool is_supercell = true;
        if(!writeToMOPAC((char *)filename.data(), &(Material.orgAtomnet), is_supercell)) {error=true; break;}
      }
      else if(command[0].compare("-res") == 0){ // free sphere calcluation
        string filename = processFilename(command, name, ".res", 0, 1);
        if(filename.empty()) {error=true; break;}
        Material.runVoroFlat();
        NEWcalculateFreeSphereParameters(&Material);
        NEWcalculateFreeSphereParametersPrint(&Material, (char *)filename.data(), extendedOutputFlag);
      }
      else if(command[0].compare("-strinfo") == 0 || command[0].compare("-strinfoex") == 0){
                                           // print information about the structure/topology
                                           // extended option saves framework/molecules with their ids next to each atom
        string filename = processFilename(command, name, ".strinfo", 0, 1);
        string filenameExtOutput = processFilename(command, name, ".frameid", 0, 1); // filename to save ext. output
        if(filename.empty()) {error=true; break;}
        if(command[0].compare("-strinfo") == 0)
          getStructureInformation((char *)filename.data(), (char *)filenameExtOutput.data(), &(Material.orgAtomnet), false); // last argument 4 ext. output
          else
          getStructureInformation((char *)filename.data(), (char *)filenameExtOutput.data(), &(Material.orgAtomnet), true);
      }
      else if(command[0].compare("-zvis") == 0){
        string filename = processFilename(command, name, ".zvis", 0, 1);
        if(filename.empty()) {error=true; break;}
        Material.runVoroFlat();
        writeZeoVisFile((char *)filename.data(), &(Material.cells), &(Material.atmnet), &(Material.vornet));
      }
      else if(command[0].compare("-zvor") == 0){  // THIS NEEDS CLEARINING
                                                  // viewVoronoiDecomp runs inflated Voronoi decomp
        string output   = processFilename(command, name, ".zvor", 0, 1);
        if(output.empty()) {error=true; break;}
        viewVoronoiDecomp(&(Material.atmnet), probeRadius, output);
      }
      else if(command[0].compare("-axs") == 0){
        string filename =  processFilename(command, name, ".axs", 0, 1);
        if(filename.empty()) {error=true; break;}
        vector<bool> accessInfo;
        vector<CHANNEL> channels;
        Material.runVoroFlat();
        CHANNEL::findChannels(&(Material.vornet), probeRadius, &accessInfo, &channels);

  //      analyze_accessible_voronoi_pre_segment(&vornet, probeRad, &accessInfo, name);
        fstream output; output.open(filename.c_str(), fstream::out);
        for(unsigned int i = 0; i < accessInfo.size(); i++)
          output << (accessInfo.at(i) ? "true" : "false") << "\n";
        output.close();
      }
      else if(command[0].compare("-visVoro") == 0){ // for visualization
        if(command.size()!=1 && command.size()!=4) {
          printf("Error: -visVoro option accepts 0 or 3 (a, b and c shifts for illustrating accessible part of network) arguments but %d arguments were supplied.\n", (int)(command.size() - 1));
          printf("Exiting...\n");
          error=true; break;
        }
        int skel_a = 0, skel_b = 0, skel_c = 0;
        if(command.size()==4) {
          //if shift was provided
          skel_a = strtod(command[1].data(), NULL), skel_b = strtod(command[2].data(), NULL), skel_c = strtod(command[3].data(), NULL);
        };
        Material.runVoroFlat();
       //        visVoro(name, probeRad, skel_a, skel_b, skel_c, &vornet, &atmnet);
        Material.visualizeVoroNet(name, probeRadius, skel_a, skel_b, skel_c);
      }

  //Rich: holograms as a separate flag, -holo
      else if(command[0].compare("-holo") == 0){
        vector<bool> accessInfo;
        vector<CHANNEL> channels;
        if(command.size()<1 || command.size()>2) {
          printf("Error: -holo option accepts no or 1 (bin directory) argument but %d arguments were supplied.\n", (int)(command.size() - 1));
          printf("Exiting...\n");
 //            exit(1);
          error=true; break;
        }
        Material.runVoroFlat();
        CHANNEL::findChannels(&(Material.vornet), probeRadius, &accessInfo, &channels);
        if(command.size()==2) { //i.e. a directory was specified
          char *bin_directory = (char *)command[1].data();
          if(bin_directory[strlen(bin_directory)-1]!='/') {
            strcat(bin_directory, "/");
          }
          analyze_accessible_voronoi_pre_segment(&(Material.vornet), probeRadius, &accessInfo, name, bin_directory); //this line generates and outputs Voronoi holograms with bins from a specified directory
        }
        analyze_accessible_voronoi_pre_segment(&(Material.vornet), probeRadius, &accessInfo, name); //this line generates and outputs Voronoi holograms with default bins 
      }
      else if(command[0].compare("-chan") == 0 || command[0].compare("-zchan") == 0){
        bool visualize = false;
        string filename =  processFilename(command, name, ".chan", 0, 1);
        if(command[0].compare("-zchan") == 0)
          {
          visualize = true;
          filename =  processFilename(command, name, ".zchan", 0, 1);
          };

        if(filename.empty()) {error=true; break;}
        vector<bool> accessInfo;
        vector<CHANNEL> channels;
        Material.runVoroFlat();
        CHANNEL::findChannels(&(Material.vornet), probeRadius, &accessInfo, &channels);

        fstream output;
        output.open(filename.data(), fstream::out);
        if(visualize == false)
          {
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
          output << filename << " summary(Max_of_columns_above)   " << di << " " << df << "  " << dif << 
                                "  probe_rad: " << probeRadius << "  probe_diam: " << 2*probeRadius << "\n";
          }else{
          for(unsigned int i = 0; i < channels.size(); i++){
            channels.at(i).writeToVMD(i, output);
            }
          output << "set num_channels " << channels.size() << "\n";

          };
        output.close();
      } // ends -chan
      else if(command[0].compare("-poreinfo") == 0){
        int numSamplesAV = 100;
        int numSamplesASA = 2500;
        string filename =  processFilename(command, name, ".poreinfo", 0, 1);
        if(filename.empty()) {error=true; break;}
        vector<bool> accessInfo;
        vector<PORE> pores;
        Material.runVoroFlat();
        PORE::findChannelsAndPockets(&(Material.vornet), probeRadius, &accessInfo, &pores);

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

        double volume = calcDeterminant(Material.atmnet.ucVectors); // Unit cell volume/Units of A^3
        int numSamples = (int)(volume*numSamplesAV);
        calcAV(&(Material.atmnet), &(Material.orgAtomnet), Material.highAccuracy, probeRadius, probeRadius, numSamples, true, output, (char *)filename.c_str(), 0, 0, 0, 0, -1,-1, 0);
        calcASA(&(Material.atmnet), &(Material.orgAtomnet), Material.highAccuracy, probeRadius, probeRadius, calcDensity(&(Material.atmnet)), numSamplesASA, true, output, (char *)filename.c_str(), 0, 0, 0, 0);

        for(unsigned int i = 0; i < pores.size(); i++)
          if(pores[i].dimensionality>0) pores[i].printPoreSummary(output, &(Material.atmnet)); // Channels

        for(unsigned int i = 0; i < pores.size(); i++)
          if(pores[i].dimensionality==0) pores[i].printPoreSummary(output, &(Material.atmnet)); // Pockets


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
      } // ends -poreinfo
      else if(command[0].compare("-poreinfoSummary")==0){  // analyzes a series of .poreinfo files  
        string filename_output =  processFilename(command, name, ".poreinfoSummary", 1, 2);
        if(filename_output.empty()) {error=true; break;}
        string filename_InputData = string(command[1]);
        analyzePoreInfoFiles(&(Material.atmnet), filename_InputData, filename_output);
      }

      else if(command[0].compare("-test")==0){  // test new functionality   
        //string filename_output =  processFilename(command, name, ".poreinfoSummary", 1, 2);
        //if(filename_output.empty()) {error=true; break;}
        //string filename_InputData = string(command[1]);
        Material.test(probeRadius);
      }
      else if(command[0].compare("-vol")==0){
        string basefilename = processFilename(command, name, "", 0, 1);
        string filename_output = basefilename + ".vol";
        fstream output;
        output.open(filename_output.data(), fstream::out);
        Material.AVcalc(probeRadius, nMCAVsampleDensity, output, name);
        output.close();
        if(Material.VisFlag == true)
           {
           string filename_vis;
           if(Material.VisSetting == "ZEOVIS") filename_vis = basefilename + ".zvol";
           if(Material.VisSetting == "VISIT") filename_vis = basefilename + ".vvol";
           if(Material.VisSetting == "LIVERPOOL") filename_vis = basefilename + ".lvol";
           output.open(filename_vis.data(), fstream::out);
           Material.AVreportPoints(output);
           output.close();
           };
        // additional calculations for accessible/inaccessible points from AV calc
        if(Material.AVrequestBlockingPockets == true) 
           {
           string filename_block = basefilename + ".block";
           output.open(filename_block.data(), fstream::out); 
           Material.AVblockPockets(output);
           output.close();
           };
        if(Material.AVrequestPSD == true)
           {
           string filename_block = basefilename + ".psd";
           output.open(filename_block.data(), fstream::out);
           Material.AVcalcPoreSizeDistr(output);
           output.close();
           if(Material.VisFlag == true)
              {
              string filename_vis;
              if(Material.VisSetting == "VISIT") filename_vis = basefilename + ".vpsd";
              if(Material.VisSetting == "LIVERPOOL") filename_vis = basefilename + ".lpsd";
              output.open(filename_vis.data(), fstream::out);
              Material.AVreportPSDPoints(output);
              output.close();
              };
           };
      } // ends -vol
      else if(command[0].compare("-sa")==0){
        string basefilename = processFilename(command, name, "", 0, 1);
        string filename_output = basefilename + ".sa";
        fstream output;
        output.open(filename_output.data(), fstream::out);
        Material.ASAcalc(probeRadius, nMCASAsampleDensity, output, name);
        output.close();
        if(Material.VisFlag == true)
           {
           string filename_vis;
           if(Material.VisSetting == "ZEOVIS") filename_vis = basefilename + ".zsa";
           if(Material.VisSetting == "VISIT") filename_vis = basefilename + ".vsa";
           if(Material.VisSetting == "LIVERPOOL") filename_vis = basefilename + ".lsa";
           output.open(filename_vis.data(), fstream::out);
           Material.ASAreportPoints(output);
           output.close();
           };
      } // ends -sa

    /* Functionality related with PLD calculations */
      else if(command[0].compare("-pld")==0 || command[0].compare("-pldf")==0 || command[0].compare("-pldmol")==0){
        string basefilename = processFilename(command, name, "", 1, 2);
        string filename_output = basefilename + ".pld";
        double segmentingRadius  = -1; // radius used for inital segmentation of strucutre
        string segmentingFile; // file with definition of segments
        if(command[0].compare("-pld")==0)
          {
          segmentingRadius  = strtod(command[1].data(), NULL);
          } else
          {
          segmentingFile = command[1].data();
          };

        fstream output;
        output.open(filename_output.data(), fstream::out);
        Material.PLDcalc(probeRadius, segmentingRadius, segmentingFile, output, name);
        output.close();
        if(Material.VisFlag == true)
           {
           Material.PLDvisualize(basefilename, Material.VisSetting);
           };
      } // ends -pld


    /* Zeolite specific functionality */

      /* Silica substitution functions */
      else if(command[0].compare("-fsub") == 0){
        /* This Si substitution function changes every second Si into Al
           then convert them back to Si according to the request fraction (frac argument)
           This function fails for system in which topology does not allow such initial substitution
           and therefore the next block of function presents an alternative */
        string filename = processFilename(command, name, "_fsub.cssr", 1, 2);
        if(filename.empty()) {error=true; break;}
        double frac     = strtod(command[1].data(), NULL);
        ATOM_NETWORK newNetwork;
//        newNetwork.allowAdjustCoordsAndCellFlag = allowAdjustCoordsAndCell;
        bool success = false;
        int numSubs = 0;
        double subFrac = 0.0;
        // Change either parameter to generate other possible configurations
        int  randSeed = 105543297; // Random number generator seed. 
        bool subSeed = false;      // Substitute first Si when creating initial 50/50 configuration 
        if(frac < 0 || frac > 0.5){
           cerr << "Error: Frac argument needs to be a number from the range of (0,0.5>\n" << "/n";
           } else
          success = fracSubstituteAtoms(&(Material.orgAtomnet), &newNetwork, subSeed, frac, randSeed, &numSubs, &subFrac, true);
        if(!success){
          cerr << "Error: Atom substitution unsucessful." << "\n"
               << "Output not written to file." << "\n";
        }
        else {
          writeToCSSR((char *)filename.data(), &newNetwork);
          cout << "Substituted " << numSubs << " Si atoms for a fraction of " << subFrac << " when analyzing " << name << "\n";
        }
      } // ends -fsub
      else if(command[0].compare("-fsubM") == 0){
      /* Maciek's version of substitution */
      /* This function has two modifications w.r.t -fsub. It allows chaging random seed to generate different configurations
         and removes topology check which ensured that every other atom is Si (in some cases achieving the requested Al fraction 
         will not be possible. */
      /* usage: -fsubM fraction random_seed(0 for default) outname */
        string filename = processFilename(command, name, "_fsub.cssr", 2, 3);
        if(filename.empty()) {error=true; break;}
        double frac     = strtod(command[1].data(), NULL);
        int rsd         = int(strtod(command[2].data(), NULL));
        ATOM_NETWORK newNetwork;
//        newNetwork.allowAdjustCoordsAndCellFlag = allowAdjustCoordsAndCell;
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
        bool success = fracSubstituteAtoms_Maciek((Material.orgAtomnet), newNetwork, subSeed, frac, randSeed, numSubs, subFrac, true);
        if(!success){
          cerr << "Error: Atom substitution unsucessful." << "\n"
               << "Output not written to file." << "\n";
        }
        else {
          writeToCSSR((char *)filename.data(), &newNetwork);
          cout << "Substituted " << numSubs << " Si atoms for a fraction of " << subFrac << " when analyzing " << name << "\n";
        }
      } // ends -fsubM


   /* SOME More-or-less Random functions */
      //flag to print to terminal the tetrahedrality of, e.g. Zn4 groups inside a material
      else if(command[0].compare("-findTetrahedra") == 0) {
        if(command.size()!=2) {
          printf("Error: -findTetrahedra option accepts 1 (element type) argument but %d arguments were supplied.\n", (int)(command.size() - 1));
          printf("Exiting...\n");
          error=true; break;
        }
        string element(command[1]);
        vector<double> tetras = Material.orgAtomnet.find_tetrahedra(element);
        //report
        printf("%s %d %s tetrahedra", name, (int)(tetras.size()), element.c_str());
        for(int i=0; i<(int)(tetras.size()); i++) printf(" %.3f", tetras.at(i));
        printf("\n");
      } // ends -findTetrahedra

      //flag to write to terminal the size of supercell required in order to satisfy 'minimum image convention' (correct terminology?) 
      // - how many cells in each axis required so that no sphere of a given radius overlaps with itself periodically
      else if(command[0].compare("-cellmulti") == 0) {
        if(command.size()!=2) {
          printf("Error: -cellmulti option accepts 1 (sphere radius) argument but %d arguments were supplied.\n", (int)(command.size() - 1));
          printf("Exiting...\n");
          error=true; break;
        }
        double rad = strtod(command[1].data(), NULL);
        double diam=2.0*rad;
        TRIPLET smallest_supercell = Material.orgAtomnet.getSmallestSupercell(diam);
        printf("%s supercell: %d %d %d = %d cells\n", name, smallest_supercell.x, smallest_supercell.y, smallest_supercell.z, smallest_supercell.x*smallest_supercell.y*smallest_supercell.z);
      }

    }; // end for loop over all commands

   }; // ends execute()



 };



#endif
