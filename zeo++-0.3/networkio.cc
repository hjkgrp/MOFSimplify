#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <fstream>
#include <string>
#include <sstream>

#include "networkio.h"
#include "string_additions.h"
#include "networkinfo.h"
#include "geometry.h"
#include "symbcalc.h"
#include "zeo_consts.h"
#include "network.h"        // Try to eliminate cyclic dependency in future

using namespace std;

/** Identifies the extension and the prefix present in the provided
    filename and stores them using the two provided character
    pointers. */
void parseFilename(const char * fileName, char *name, char *extension){
  string s(fileName);
  size_t index = s.find_last_of(".");
  if(index == string::npos){
    cerr << "Improper input filename " << fileName << "\n";
    cerr << "No . extension found. Exiting ..." << "\n";
    exit(1);
  }
  else{
    string prefix = s.substr(0, index);
    string suffix = s.substr(index + 1);
    strncpy(name, prefix.data(), prefix.size()); name[prefix.size()] = '\0';
    strncpy(extension, suffix.data(), suffix.size()); extension[suffix.size()] = '\0';
  }
}

/** Ensures that the provided file is of type .arc, .car, .cuc, .cssr or .v1 or .cif. 
 *  Otherwise, an error message is displayed and the program is
 *  aborted. */
bool checkInputFile(char * filename){
  string file(filename);
  string fileTypes [] = {".cuc", ".arc", ".cssr", ".obcssr", ".v1",".cif",".car",".dlp",".pdb"};
  int numTypes = 8;

  for(int i = 0; i < numTypes; i++){
    if(file.find(fileTypes[i]) != string::npos)
      return true;
  }

  cerr << "Invalid input filename " << filename << "\n" << "Exiting ..." << "\n";
//  exit(1);
  return false;
}

/** Read the information from the .cif file referred to by filename
    and store it within the provided ATOM_NETWORK. */
//void readCIFFile(char *filename, ATOM_NETWORK *cell, bool radial){
bool readCIFFile(char *filename, ATOM_NETWORK *cell, bool radial){
  string line;
  // Known markers in CIF File if you change the order it will affect the code
  string descriptor[] = {"data_",
			 "_cell_length_a", //case 1
			 "_cell_length_b", //2
			 "_cell_length_c", //3
			 "_cell_angle_alpha", //4
			 "_cell_angle_beta", //5
			 "_cell_angle_gamma", //6
			 "loop_", //7
			 "_symmetry_equiv_pos_as_xyz", //8
			 "_space_group_symop_operation_xyz", //9
			 "_atom_site_label", //10
			 "_atom_site_type_symbol", //11
			 "_atom_site_fract_x", //12
			 "_atom_site_fract_y", //13
			 "_atom_site_fract_z", //14
			 "_atom_site_charge", //15
       "_symmetry_Int_Tables_number", //16
			 "_atom_site_Cartn_x", //17
			 "_atom_site_Cartn_y", //18
			 "_atom_site_Cartn_z", //19
                         "_symmetry_equiv_pos_site_id", // 20 
			 "NULL"};  
  int ndx;
  vector<string> list = strAry2StrVec(descriptor);
  vector<string> token;

  vector<string> sym_x;
  vector<string> sym_y;
  vector<string> sym_z;
  vector<string> atom_label;
  vector<string> atom_type;
  vector<double> atom_x;
  vector<double> atom_y;
  vector<double> atom_z;
  vector<double> atom_charge;

  int symmetry_Int_Table_number = -1; //set to dummy value so we can check if anything was read
  bool symmetry_equiv_pos_site_id_Flag = false; // this is to read symmetry lines from lines that have sym. id.
 
  // Try opening the file if it opens proceed with processing
  ifstream ciffile;
  cout << "Opening File: " << filename << endl;
  ciffile.open(filename);
  bool read_a = false, read_b = false, read_c = false, read_alpha = false, read_beta = false, read_gamma = false, initialized_cell = false; //keep track of when all cell params are parsed, so the cell can be created
  if(ciffile.is_open()) {
    while(!ciffile.eof()) {
      if(read_a && read_b && read_c && read_alpha && read_beta && read_gamma && !initialized_cell) {cell->initialize(); initialized_cell=true;}
      getline(ciffile,line);
//printf("DEBUG: read line %s\n", line.c_str());
      token = split(line," ()\r\t");
      exception: //I needed an easy way to jump out of the _loop command if an unknown command was found
      if (token.size() > 0) {
	      //Where all non-loop commands should be added
	      if(token[0].substr(0,5).compare(list[0]) == 0){ //name of unit cell
	        cell->name=token[0].substr(5);
	      }
	      else if (token[0].compare(list[1]) == 0){ //length a
	        cell->a=convertToDouble(token[1]);
          read_a = true;
	      }
	      else if (token[0].compare(list[2]) == 0){ //length b
	        cell->b=convertToDouble(token[1]);
          read_b = true;
	      }
	      else if (token[0].compare(list[3]) == 0){ //length c
	        cell->c=convertToDouble(token[1]);
          read_c = true;
	      }
	      else if (token[0].compare(list[4]) == 0){ //alpha
	        cell->alpha=convertToDouble(token[1]);
          read_alpha = true;
	      }
	      else if (token[0].compare(list[5]) == 0){ //beta
	        cell->beta=convertToDouble(token[1]);
          read_beta = true;
	      }
	      else if (token[0].compare(list[6]) == 0){ //gamma
	        cell->gamma=convertToDouble(token[1]);
          read_gamma = true;
	      }
        else if (token[0].compare(list[16]) == 0){ //_symmetry_Int_Tables_number
	        symmetry_Int_Table_number = convertToInt(token[1]);
	      }
	      else if (token[0].compare(list[7]) == 0){ //loop_
//printf("DEBUG: inside a \"loop_\" section\n");
	        vector<string> column_labels;
	        getline(ciffile,line);
	        token = split(line," \r\t");
          bool tokenized = false, in_loop = true;
          if(token.size()>0) tokenized=true;
//	        while (token[0].at(0)=='_') { //collect all of the collumn labels
	        while (tokenized && in_loop) { //collect all of the collumn labels
            if(token[0].at(0)=='_') {
	            column_labels.push_back(token[0]);
//printf("DEBUG: within loop, parsed a column header \"%s\"\n", token[0].c_str());
	            getline(ciffile,line);
//printf("DEBUG: read line \"%s\" - tokenizing ...\n", line.c_str());
	            token = split(line," \r\t");
              if(token.size()<=0) tokenized=false;
            } else in_loop = false;
	        }
          if(!tokenized) {
printf("\n#####\n##### WARNING: parsed a loop in cif file, but data was not present\n#####\n\n");
          } else {
//printf("DEBUG: exited loop successfully, having read data line \"%s\"\n", line.c_str());
            //collecting the data associated with each column and put
            // in correct place
            token = split(line," ,'\r\t"); //This is needed to split data
            bool need_to_convert_to_fractional = false;
            while(token.size()>0){
              if (token[0].at(0) =='_' || token[0].compare("loop_")==0 || token[0].at(0) == '#'){
                goto exception; // unexpected input
              }
              if (token.size()!=column_labels.size() && column_labels[0].compare(list[8]) != 0 && column_labels[0].compare(list[9]) != 0 && column_labels[0].compare(list[20]) != 0){
                goto exception; //unexpected input
              }
              for (unsigned int i=0; i<column_labels.size(); i++){
                switch (strCmpList(list,column_labels[i])){
  //printf("SYM DEBUG: checking for symmetry section ...\n");
	              //Where all loop commands should be added
                case 8: //_symmetry_equiv_pos_as_xyz and
                case 9: //_space_group_symop_operation_xyz have the same meaning
  //printf("SYM DEBUG: symmetry section found!\n");
	              if (!((token.size()==3&&symmetry_equiv_pos_site_id_Flag==false)||(token.size()==4&&symmetry_equiv_pos_site_id_Flag==true))){
	                cerr << "Error: Expected 3 strings for _symmetry_equiv_pos_as_xyz (or 4 if _symmetry_equiv_pos_site_id is present)" << endl;
  //	              abort();
                  ciffile.close();
                  return false;
	              }
                      if(token.size()==3) {
	                sym_x.push_back(token[0]);
	                sym_y.push_back(token[1]);
	                sym_z.push_back(token[2]);}
                        else {
                        sym_x.push_back(token[1]);
                        sym_y.push_back(token[2]);
                        sym_z.push_back(token[3]);
                        };
  //printf("SYM DEBUG: pushing back %s %s %s\n", token[0].c_str(), token[1].c_str(), token[2].c_str());
	              break;
                case 10:
	              atom_label.push_back(token[i]);
	              break;
                case 11:
	              atom_type.push_back(token[i]);
	              break;
                case 12:
	              atom_x.push_back(trans_to_origuc(convertToDouble(token[i])));
	              break;
                case 13:
	              atom_y.push_back(trans_to_origuc(convertToDouble(token[i])));
	              break;
                case 14:
	              atom_z.push_back(trans_to_origuc(convertToDouble(token[i])));
	              break;
                case 15:
	              atom_charge.push_back(convertToDouble(token[i]));
	              break;
                case 17:
	              atom_x.push_back(convertToDouble(token[i]));
                need_to_convert_to_fractional = true;
	              break;
                case 18:
	              atom_y.push_back(convertToDouble(token[i]));
                need_to_convert_to_fractional = true;
	              break;
                case 19:
	              atom_z.push_back(convertToDouble(token[i]));
                need_to_convert_to_fractional = true;
	              break;
                case 20:
                      symmetry_equiv_pos_site_id_Flag = true;
                      break;
                }	 
              }  
              //now we've finished parsing this line, we might need to convert Cartesian to fractional coords
              if(need_to_convert_to_fractional) {
                Point tempcoor;
                int this_ID = atom_x.size()-1;
        	      tempcoor = cell->xyz_to_abc(atom_x.at(this_ID),atom_y.at(this_ID),atom_z.at(this_ID));
//printf("DEBUG: atom at Cartesian %.3f %.3f %.3f written to fractional %.3f %.3f %.3f\n", atom_x.at(this_ID),atom_y.at(this_ID),atom_z.at(this_ID), tempcoor[0], tempcoor[1], tempcoor[2]);
                atom_x.at(this_ID) = tempcoor[0];
                atom_y.at(this_ID) = tempcoor[1];
                atom_z.at(this_ID) = tempcoor[2];
              }
              getline(ciffile,line);
              token = split(line," ',\r\t");
            }
            column_labels.clear();
          }
        }
        token.clear();
      }
    }	  
    ciffile.close();
    
    //If no symmetry info was provided, we can assume that none is required (i.e., structure has 'P 1' symmetry), ONLY IF we did not read a symmetry_Int_Tables_Number!=1
    if (sym_x.size()==0 || sym_y.size()==0 || sym_z.size()==0){
      //no symmetry provided
      if (symmetry_Int_Table_number>=0 && symmetry_Int_Table_number!=1) {
        //read a symmetry table number, but it was not 1
        printf("ERROR:\n\tcif file provided no symmetry operations; however, a symmetry_Int_Tables_Number of %d was provided,\n\tindicating symmetry which is not 'P 1' (no symmetry operations);\n\tcannot proceed with insufficient symmetry information\nExiting ...\n", symmetry_Int_Table_number);
//        exit(EXIT_FAILURE);
        return false;
      } else {
        sym_x.push_back("x");
        sym_y.push_back("y");
        sym_z.push_back("z");
      }
    }

    //Now determine whether it is correct to use atom_label or atom_type (atom_label used only if atom_type has no data)
    bool CIF_USE_LABEL = atom_type.size()==0;

    //Now determine whether charges are used - only if they were provided
    bool CIF_CHARGES = atom_charge.size()>0;

    //Parse out the numbers from atom labels and types, if specified (defined in network.h)
    if (CIF_RMV_NUMS_FROM_ATOM_TYPE){
      if (!CIF_USE_LABEL){
	      for (unsigned int i=0; i<atom_type.size(); i++){
	        atom_type[i] = split(atom_type[i],"0123456789").at(0);
	      }
      }
      else{
	      for (unsigned int i=0; i<atom_label.size(); i++){
	        atom_label[i] = split(atom_label[i],"0123456789").at(0);
	      }
      }
    }
    
    //Error checking
    if (sym_x.size()==0 || sym_y.size()==0 || sym_z.size()==0 ){
      cerr << "Error: No .cif symmetry information given" << endl;
      cerr << "DEBUG SHOULDN'T HAPPEN" << endl;
      return false;
//      abort();
    }
    if (atom_label.size()==0 && CIF_USE_LABEL == true){
      cerr << "Error: No ''_atom_site_label'' or ''_atom_site_type_symbol'' information " << endl;
      cerr << "This structure appears to be invalid: there are no atom identifying element types or labels - if this is not the case, then this is a bug; please contact the developers." << endl;
      return false;
    }
    if (atom_type.size()==0 && CIF_USE_LABEL == false){
      cerr << "Error: No ''_atom_site_type_symbol'' information" << endl;
      cerr << "DEBUG SHOULDN'T HAPPEN" << endl;
      return false;
//      abort();
    }
    if (atom_x.size()!=atom_y.size() || atom_y.size()!=atom_z.size() || atom_x.size()==0){
      cerr << "Error: Atom coordinates not read properly (" << atom_x.size() << " x coords, " << atom_y.size() << " y coords, " << atom_z.size() << " z coords)" << endl;
//      abort();
      return false;
    }

    //If no name for the unit cell is given the filename becomes its name
    if (cell->name.size() == 0){
      cell->name = *filename;
    }
    
    //Now to fully assemble the ATOM_NETWORK initialize constructor
//    cell->initialize(); //this should now happen much earlier, during parsing
    
    ATOM tempatom;
    Point tempcoor;
    for (unsigned int i=0; i<atom_x.size(); i++){ //atoms
      if (CIF_USE_LABEL){
	      tempatom.type = atom_label[i];
      }
      else {
	      tempatom.type = atom_type[i];
      }
      tempatom.radius = lookupRadius(tempatom.type, radial);
      if(CIF_CHARGES) {
        tempatom.charge = atom_charge[i];
      } else tempatom.charge = 0;
      for (unsigned int j=0;j<sym_x.size(); j++){ //symmetries
	      tempatom.a_coord = trans_to_origuc(symbCalc(sym_x[j],atom_x[i],atom_y[i],atom_z[i]));
	      tempatom.b_coord = trans_to_origuc(symbCalc(sym_y[j],atom_x[i],atom_y[i],atom_z[i]));
	      tempatom.c_coord = trans_to_origuc(symbCalc(sym_z[j],atom_x[i],atom_y[i],atom_z[i]));
	      tempcoor = cell->abc_to_xyz (tempatom.a_coord,tempatom.b_coord,tempatom.c_coord);
	      tempatom.x = tempcoor[0]; 
	      tempatom.y = tempcoor[1]; 
	      tempatom.z = tempcoor[2];
	      tempatom.specialID = (i*sym_x.size())+j;
	
	      //make sure that no duplicate atoms are writen
	      int match=1;
	      for (unsigned int k=0;k<cell->atoms.size(); k++){
                if(cell->calcDistance(cell->atoms[k],tempatom)<thresholdLarge) { match=0;
/*
	        if (tempatom.a_coord-thresholdLarge < cell->atoms[k].a_coord && cell->atoms[k].a_coord < tempatom.a_coord+thresholdLarge)
          if (tempatom.b_coord-thresholdLarge < cell->atoms[k].b_coord && cell->atoms[k].b_coord < tempatom.b_coord+thresholdLarge )
          if (tempatom.c_coord-thresholdLarge < cell->atoms[k].c_coord && cell->atoms[k].c_coord < tempatom.c_coord+thresholdLarge ){
            match=0;
*/
          }
        }
        if (match == 1){
          cell->atoms.push_back(tempatom);
        }
      }
    }
    cell->numAtoms = cell->atoms.size();
  }
  else{
    cout << "Failed to open: " << filename << endl;
    ciffile.close();
//    exit(1);
    return false;
  }
  return true;
}

/** Read the .arc file refererred to by filename and store its
    information within the provided ATOM_NETWORK. */
//void readARCFile(char *filename, ATOM_NETWORK *cell, bool radial){
bool readARCFile(char *filename, ATOM_NETWORK *cell, bool radial){
  FILE * input;
  input = fopen(filename, "r");
  int numAtoms=0;

  if(input==NULL){
    cout << "\n" << "Failed to open .arc input file " << filename << "\n";
    cout << "Exiting ..." << "\n";
//    exit(1);
    return false;
  }
  else{
    cout << "Reading input file " << filename << "\n";

    // Parse down to the FINAL GEOMETRY OBTAINED line
    char this_line[500];
    char found_final = 0;
    while(found_final == 0) {
      if(fgets(this_line, 500, input)!=NULL) {
        char str1[100], str2[100], str3[100];
        int status = sscanf(this_line, "%s %s %s", str1, str2, str3);
        if(status!=-1) {
          if(strcmp(str1,"FINAL")==0 && strcmp(str2,"GEOMETRY")==0 && strcmp(str3,"OBTAINED")==0) {
            found_final = 1;
          }
        }
      } else {
        printf("ERROR: finished parsing ARC file before finding geometry section\n");
//        exit(EXIT_FAILURE);
        fclose(input);
        return false;
      }
    }

    // Now parse following lines, trying to extract atom info
    int found_atoms = 0;
    double x, y, z, charge;
    char element[100], str1[100], str2[100], str3[100];
    while(found_atoms==0) {
      if(fgets(this_line, 500, input)!=NULL) {
        int status = sscanf(this_line, "%s %lf %s %lf %s %lf %s %lf", element, &x, str1, &y, str2, &z, str3, &charge);
        if(status==8) { //i.e. exactly 8 fields are required to be read
          found_atoms = 1;
        }
      } else {
        printf("ERROR: finished parsing ARC file before finding individual atom information\n");
//        exit(EXIT_FAILURE);
        fclose(input);
        return false;
      }
    }

    // At this point, we have the first atom info in memory
    ATOM newAtom;
    while(found_atoms==1) {
      //save previously discovered atom data - just xyz data for now
      newAtom.x = x;
      newAtom.y = y;
      newAtom.z = z;
      newAtom.type = string(element);
      newAtom.radius = lookupRadius(newAtom.type, radial);
      newAtom.charge = charge;
      cell->atoms.push_back(newAtom);
      numAtoms++;
      //try to find next atom data
      if(fgets(this_line, 500, input)!=NULL) {
        int status = sscanf(this_line, "%s %lf %s %lf %s %lf %s %lf", element, &x, str1, &y, str2, &z, str3, &charge);
        if(status!=8) found_atoms = 0; //if we can't read all 8 fields, we have, presumably, reached the unit cell params
      } else {
        printf("ERROR: finished parsing ARC file before finding unit cell info\n");
//        exit(EXIT_FAILURE);
        fclose(input);
        return false;
      }
    }

    // Now we have read all the atoms, we read the unit cell vectors
    XYZ v_a, v_b, v_c;
    for(int i=0; i<3; i++) {
      if(i==0) {
        v_a.x = x; v_a.y = y; v_a.z = z;
      } else if(i==1) {
        v_b.x = x; v_b.y = y; v_b.z = z;
      } else if(i==2) {
        v_c.x = x; v_c.y = y; v_c.z = z;
      }
      if(i!=2) {
        if(fgets(this_line, 500, input)!=NULL) {
          int status = sscanf(this_line, "%s %lf %s %lf %s %lf %s", element, &x, str1, &y, str2, &z, str3);
          if(status!=7) { //didn't read exactly 7 fields
            printf("ERROR: could not read exactly three unit cell vectors\n");
//            exit(EXIT_FAILURE);
            fclose(input);
            return false;
          }
        }
      }
    }
    cell->numAtoms = numAtoms;
    fclose(input);

    // Set up cell
    cell->v_a = v_a;
    cell->v_b = v_b;
    cell->v_c = v_c;
    double alpha = v_b.angle_between(v_c);
    double beta = v_a.angle_between(v_c);
    double gamma = v_a.angle_between(v_b);
    cell->alpha = alpha*360.0/(2.0*PI);
    cell->beta = beta*360.0/(2.0*PI);
    cell->gamma = gamma*360.0/(2.0*PI); //required since it expects to store degrees
    cell->a = v_a.magnitude();
    cell->b = v_b.magnitude();
    cell->c = v_c.magnitude();
    cell->initMatrices();
    cell->name = filename;
    cell->name.erase(cell->name.end()-4, cell->name.end());
    
//   cout << "number of atoms read " << numAtoms << "\n";

    // Convert atom coords to abc and update the xyz values to be within the UC
    for(int i=0; i<numAtoms; i++) {
      Point newCoords = cell->xyz_to_abc(cell->atoms.at(i).x,cell->atoms.at(i).y,cell->atoms.at(i).z);
      cell->atoms.at(i).a_coord = trans_to_origuc(newCoords[0]); cell->atoms.at(i).b_coord = trans_to_origuc(newCoords[1]); cell->atoms.at(i).c_coord = trans_to_origuc(newCoords[2]);
      newCoords = cell->abc_to_xyz(cell->atoms.at(i).a_coord,cell->atoms.at(i).b_coord,cell->atoms.at(i).c_coord);
      cell->atoms.at(i).x = newCoords[0]; cell->atoms.at(i).y = newCoords[1]; cell->atoms.at(i).z = newCoords[2];
      //cout << i << "  " << cell->atoms.at(i).type << "  " << newCoords[0] << "  " << newCoords[1] << "  " << newCoords[2] <<"\n";
    }
  }
  return true;
}


/** Read the .cuc file refererred to by filename and store its
    information within the provided ATOM_NETWORK. */
//void readCUCFile(char *filename, ATOM_NETWORK *cell, bool radial){
bool readCUCFile(char *filename, ATOM_NETWORK *cell, bool radial){
  fstream input;
  input.open(filename);
  char garbage[256];
  int numAtoms;

  if(!input.is_open()){
    cout << "\n" << "Failed to open .cuc input file " << filename << "\n";
    cout << "Exiting ..." << "\n";
//    exit(1);
    return false;
  }
  else{
    cout << "Reading input file " << filename << "\n";

    // Read and store information about the unit cell
    cell->name = filename;
    cell->name.erase(cell->name.end()-4, cell->name.end());
    input.getline(garbage, 256);
    input >> garbage;
    input >> cell->a >> cell->b >> cell->c;
    input >> cell->alpha >> cell->beta >> cell->gamma;
    cell->initialize(); // Initializes the unit cell vectors using the
			// angles and side lengths

    // Read and store information about each atom
    numAtoms = 0;
    while(!input.eof()){
      ATOM newAtom;
      input >> newAtom.type;
      if(newAtom.type.empty())
	break;
      
      changeAtomType(&newAtom); // Converts to atom type
      input >> newAtom.a_coord >> newAtom.b_coord >> newAtom.c_coord;
      newAtom.a_coord = trans_to_origuc(newAtom.a_coord);
      newAtom.b_coord = trans_to_origuc(newAtom.b_coord);
      newAtom.c_coord = trans_to_origuc(newAtom.c_coord);
      Point newCoords = cell->abc_to_xyz (newAtom.a_coord,newAtom.b_coord,newAtom.c_coord);
      newAtom.x = newCoords[0]; newAtom.y = newCoords[1]; newAtom.z = newCoords[2];
      newAtom.radius = lookupRadius(newAtom.type, radial);
      newAtom.label=newAtom.type;
      cell->atoms.push_back(newAtom);
      numAtoms++;
    }   
    cell->numAtoms= numAtoms; 
    input.close();
  }
  return true;
}

/** Read the information from the .cssr file referred to by filename
    and store it within the provided ATOM_NETWORK. */
//void readCSSRFile(char *filename, ATOM_NETWORK *cell, bool radial){
bool readCSSRFile(char *filename, ATOM_NETWORK *cell, bool radial){
  string garbage; //this is a garbage string to read line
  int i,j;

  fstream input;
  input.open(filename);
  if(input.is_open()==true){
    // Read the header information
    cout << "Reading input file: " << filename << endl;;
    input >> cell->a >> cell->b >> cell->c;
    input >> cell->alpha >> cell->beta >> cell->gamma;
    getline(input,garbage);

    string numStr;
    bool longCSSR=false;  // this flag wiill enable switching to a different read routine
                          // to handle files with number of atoms larger than 10000 
    bool CartCoords=false; // this flag enables reading Cartesian coordinates

    input >> numStr >> CartCoords;
    getline(input,garbage);
//    input >> cell->numAtoms;
    if(numStr.compare("****") == 0) longCSSR=true;

//    input >> i >> cell->name;  // some files have '0' preceeding name
//    getline(input,garbage);
    getline(input,cell->name);
    cell->initialize();

    if(longCSSR == false){
      cell->numAtoms = atoi(numStr.c_str());
    
      // Read and store information about each atom
      for(j=0;j<cell->numAtoms;j++) {
        ATOM newAtom;

        if(CartCoords==false) {
           input >> newAtom.specialID >> newAtom.type >> newAtom.a_coord >> newAtom.b_coord >>newAtom.c_coord;
           newAtom.a_coord = trans_to_origuc(newAtom.a_coord);
           newAtom.b_coord = trans_to_origuc(newAtom.b_coord);
           newAtom.c_coord = trans_to_origuc(newAtom.c_coord);
           Point newCoords=cell->abc_to_xyz(newAtom.a_coord,newAtom.b_coord,newAtom.c_coord);
           newAtom.x=newCoords[0]; newAtom.y=newCoords[1]; newAtom.z=newCoords[2];}
         else{ // read-in cartesian
           input >> newAtom.specialID >> newAtom.type >> newAtom.x >> newAtom.y >> newAtom.z;
           Point newCoords=cell->xyz_to_abc(newAtom.x,newAtom.y,newAtom.z);
           newAtom.a_coord=newCoords[0]; newAtom.b_coord=newCoords[1]; newAtom.c_coord=newCoords[2];
           newAtom.a_coord = trans_to_origuc(newAtom.a_coord);
           newAtom.b_coord = trans_to_origuc(newAtom.b_coord);
           newAtom.c_coord = trans_to_origuc(newAtom.c_coord);
           newCoords=cell->abc_to_xyz(newAtom.a_coord,newAtom.b_coord,newAtom.c_coord);
           newAtom.x=newCoords[0]; newAtom.y=newCoords[1]; newAtom.z=newCoords[2];
           };
         
        newAtom.radius=lookupRadius(newAtom.type, radial);
        cell->atoms.push_back(newAtom);
//        getline(input,garbage);
        //read these numbers after the coords, since the final field contains charge information
        int empty_int = 0;
        for(int k=0; k<8; k++) input >> empty_int;
        input >> newAtom.charge;
        }
      }
    else{ // start longCSSR=true
      cout << "Long CSSR file. Switching to another reading routine.\n" ;
      int na=1;
      while(!input.eof())
        {
        ATOM newAtom;
        newAtom.specialID = na;
        input >> garbage;
        if(input.eof()) 
          {
          na--;
          break;
          };
        if(CartCoords==false) {
          input >> newAtom.type >> newAtom.a_coord >> newAtom.b_coord >>newAtom.c_coord;
          newAtom.a_coord = trans_to_origuc(newAtom.a_coord);
          newAtom.b_coord = trans_to_origuc(newAtom.b_coord);
          newAtom.c_coord = trans_to_origuc(newAtom.c_coord);
          Point newCoords=cell->abc_to_xyz(newAtom.a_coord,newAtom.b_coord,newAtom.c_coord);
          newAtom.x=newCoords[0]; newAtom.y=newCoords[1]; newAtom.z=newCoords[2];}
        else{  // input in cartesian
          input >> newAtom.type >> newAtom.x >> newAtom.y >> newAtom.z;
          Point newCoords=cell->xyz_to_abc(newAtom.x,newAtom.y,newAtom.z);
          newAtom.a_coord=newCoords[0]; newAtom.b_coord=newCoords[1]; newAtom.c_coord=newCoords[2];
          newAtom.a_coord = trans_to_origuc(newAtom.a_coord);
          newAtom.b_coord = trans_to_origuc(newAtom.b_coord);
          newAtom.c_coord = trans_to_origuc(newAtom.c_coord);
          newCoords=cell->abc_to_xyz(newAtom.a_coord,newAtom.b_coord,newAtom.c_coord);
          newAtom.x=newCoords[0]; newAtom.y=newCoords[1]; newAtom.z=newCoords[2];
          };

        newAtom.radius=lookupRadius(newAtom.type, radial);
        //read these numbers after the coords, since the final field contains charge information
        int empty_int = 0;
        for(int k=0; k<8; k++) input >> empty_int;
        input >> newAtom.charge;
        cell->atoms.push_back(newAtom);
//        cout << na << "  " << newAtom.type << " " << newAtom.a_coord << "  " << newAtom.b_coord << "  " << newAtom.c_coord << endl;
        na++;
        };
      cell->numAtoms = na;
      cout << na << " atoms read." << endl;
    };// end longCSSR=true

    input.close();
  }
  else{
    cerr << "Error: CSSR failed to open " << filename << endl;
//    abort();
    return false;
  }
  return true;
}




/** Read the information from the .obcssr file referred to by filename
    and store it within the provided ATOM_NETWORK. 
    obcssr are Open Babel generated cssr files  */
bool readOBCSSRFile(char *filename, ATOM_NETWORK *cell, bool radial){
  string garbage; //this is a garbage string to read line
  int i,j;

  fstream input;
  input.open(filename);
  if(input.is_open()==true){
    // Read the header information
    cout << "Reading input file: " << filename << endl;;
    for(int i=0; i<6; i++) input >> garbage; 
    input >> cell->a >> cell->b >> cell->c;
    getline(input,garbage);
    input >> garbage >> garbage;
    input >> cell->alpha >> cell->beta >> cell->gamma;
    getline(input,garbage);

    string numStr;
    bool longCSSR=false;  // this flag wiill enable switching to a different read routine
                          // to handle files with number of atoms larger than 10000 
    bool CartCoords=false; // this flag enables reading Cartesian coordinates

    cout << "Attempt to read OpenBabel CSSR file. Atom connectivity and charge columns will be omitted" << endl; 

    input >> numStr >> CartCoords;
    getline(input,garbage);
//    input >> cell->numAtoms;
    if(numStr.compare("****") == 0) longCSSR=true;

//    input >> i >> cell->name;  // some files have '0' preceeding name
//    getline(input,garbage);
    getline(input,cell->name);
    cell->initialize();

    if(longCSSR == false){
      cell->numAtoms = atoi(numStr.c_str());
    
      // Read and store information about each atom
      for(j=0;j<cell->numAtoms;j++) {
        ATOM newAtom;

        if(CartCoords==false) {
           input >> newAtom.specialID >> newAtom.type >> newAtom.a_coord >> newAtom.b_coord >>newAtom.c_coord;
           newAtom.a_coord = trans_to_origuc(newAtom.a_coord);
           newAtom.b_coord = trans_to_origuc(newAtom.b_coord);
           newAtom.c_coord = trans_to_origuc(newAtom.c_coord);
           Point newCoords=cell->abc_to_xyz(newAtom.a_coord,newAtom.b_coord,newAtom.c_coord);
           newAtom.x=newCoords[0]; newAtom.y=newCoords[1]; newAtom.z=newCoords[2];}
         else{ // read-in cartesian
           input >> newAtom.specialID >> newAtom.type >> newAtom.x >> newAtom.y >> newAtom.z;
           Point newCoords=cell->xyz_to_abc(newAtom.x,newAtom.y,newAtom.z);
           newAtom.a_coord=newCoords[0]; newAtom.b_coord=newCoords[1]; newAtom.c_coord=newCoords[2];
           newAtom.a_coord = trans_to_origuc(newAtom.a_coord);
           newAtom.b_coord = trans_to_origuc(newAtom.b_coord);
           newAtom.c_coord = trans_to_origuc(newAtom.c_coord);
           newCoords=cell->abc_to_xyz(newAtom.a_coord,newAtom.b_coord,newAtom.c_coord);
           newAtom.x=newCoords[0]; newAtom.y=newCoords[1]; newAtom.z=newCoords[2];
           };
         
        newAtom.radius=lookupRadius(newAtom.type, radial);
        cell->atoms.push_back(newAtom);
//uncommented the following line to reads Marcos files
        getline(input,garbage);
        //read these numbers after the coords, since the final field contains charge information
//        int empty_int = 0;
//        for(int k=0; k<8; k++) input >> empty_int;
//        input >> newAtom.charge;
        //  TEMP addition to read Marco's cssrs
//        input >> empty_int;





        }
      }
    else{ // start longCSSR=true
      cout << "Long CSSR file. Switching to another reading routine.\n" ;
      int na=1;
      while(!input.eof())
        {
        ATOM newAtom;
        newAtom.specialID = na;
        input >> garbage;
        if(input.eof()) 
          {
          na--;
          break;
          };
        if(CartCoords==false) {
          input >> newAtom.type >> newAtom.a_coord >> newAtom.b_coord >>newAtom.c_coord;
          newAtom.a_coord = trans_to_origuc(newAtom.a_coord);
          newAtom.b_coord = trans_to_origuc(newAtom.b_coord);
          newAtom.c_coord = trans_to_origuc(newAtom.c_coord);
          Point newCoords=cell->abc_to_xyz(newAtom.a_coord,newAtom.b_coord,newAtom.c_coord);
          newAtom.x=newCoords[0]; newAtom.y=newCoords[1]; newAtom.z=newCoords[2];}
        else{  // input in cartesian
          input >> newAtom.type >> newAtom.x >> newAtom.y >> newAtom.z;
          Point newCoords=cell->xyz_to_abc(newAtom.x,newAtom.y,newAtom.z);
          newAtom.a_coord=newCoords[0]; newAtom.b_coord=newCoords[1]; newAtom.c_coord=newCoords[2];
          newAtom.a_coord = trans_to_origuc(newAtom.a_coord);
          newAtom.b_coord = trans_to_origuc(newAtom.b_coord);
          newAtom.c_coord = trans_to_origuc(newAtom.c_coord);
          newCoords=cell->abc_to_xyz(newAtom.a_coord,newAtom.b_coord,newAtom.c_coord);
          newAtom.x=newCoords[0]; newAtom.y=newCoords[1]; newAtom.z=newCoords[2];
          };

        newAtom.radius=lookupRadius(newAtom.type, radial);
        //read these numbers after the coords, since the final field contains charge information
        int empty_int = 0;
        for(int k=0; k<8; k++) input >> empty_int;
        input >> newAtom.charge;
        cell->atoms.push_back(newAtom);
//  TEMP addition to read Marco's cssrs
input >> empty_int;
//        cout << na << "  " << newAtom.type << " " << newAtom.a_coord << "  " << newAtom.b_coord << "  " << newAtom.c_coord << endl;
        na++;
        };
      cell->numAtoms = na;
      cout << na << " atoms read." << endl;
    };// end longCSSR=true

    input.close();
  }
  else{
    cerr << "Error: CSSR failed to open " << filename << endl;
//    abort();
    return false;
  }
  return true;
}









/** Read the information from the .car file referred to by filename
    and store it within the provided ATOM_NETWORK. */
bool readCARFile(char *filename, ATOM_NETWORK *cell, bool radial){
  string garbage; //this is a garbage string to read line
  int i,j;

  fstream input;
  input.open(filename);
  if(input.is_open()==true){
    // Read the header information
    cout << "Reading input file: " << filename << endl;
    getline(input, garbage); // skip the first comment line

    string PBCline;

    input >> PBCline; getline(input, garbage);
    if(PBCline.compare("PBC=ON") != 0)
      {
      cerr << "This .car file does not have a periodic structure. Exiting...\n";
      return false;
      };

    getline(input, garbage); // skipping two next lines
    getline(input, garbage);

    input >> garbage; 
    input >> cell->a >> cell->b >> cell->c;
    input >> cell->alpha >> cell->beta >> cell->gamma;
    string SYMMline;
    input >> SYMMline;
    getline(input,garbage);

    if(SYMMline.compare("(P1)") !=0)
      {
      cerr << "The current .car reader does only work for (P1) symmetry.\n";
      return false;
      };

    cell->name = filename;
    cell->initialize();

    bool end=false;
    int na=0;
    while(!end)
      {
      string str1, str2, str3, str4;

      input >> str1;

      if(str1.compare("end") == 0 || str1.compare("END") == 0)
        {
        end = true;  
        } else
        {
        ATOM newAtom;
        input >> newAtom.x >> newAtom.y >> newAtom.z;
        input >> str2 >> str3 >> str4;
        input >> newAtom.type >> newAtom.charge;

        if(!CAR_USE_ATOM_TYPE_OVER_NAME) newAtom.type = str1;

        Point newCoords = cell->xyz_to_abc(newAtom.x, newAtom.y, newAtom.z);
        newAtom.a_coord = newCoords[0]; newAtom.b_coord = newCoords[1]; newAtom.c_coord = newCoords[2];
        newAtom.radius = lookupRadius(newAtom.type, radial);
        cell->atoms.push_back(newAtom);
        na++;
        };

      };

    cell->numAtoms = na;
    cout << na << " atoms read." << endl;

    input.close();
  }
  else{
    cerr << "Error: CAR failed to open " << filename << endl;
//    abort();
    return false;
  }
  return true;
}

/** Read the information from the .pdb file referred to by filename
    and store it within the provided ATOM_NETWORK. 
    This function is written to handle example .pdb files from RASPA
 */
bool readPDBFile(char *filename, ATOM_NETWORK *cell, bool radial){
  string garbage; //this is a garbage string to read line
  int i,j;

  fstream input;
  input.open(filename);
  if(input.is_open()==true){
    // Read the header information
    cout << "Reading input file: " << filename << endl;
    getline(input, garbage); // skip the first comment line

    string PBCline;

    input >> PBCline; 
    if(PBCline.compare("CRYST1") != 0)
      {
      cerr << "This .pdb files does not contain CRYST1 in the second line. File format not compatible. Exiting...\n";
      return false;
      };

    input >> cell->a >> cell->b >> cell->c;
    input >> cell->alpha >> cell->beta >> cell->gamma;
    getline(input,garbage); // reading rest of line

    cell->name = filename;
    cell->initialize();

    bool end=false;
    int na=0;
    while(!end)
      {
      string str1, str2, str3, str4;

      input >> str1; // reads ATOM or ENDMDL keyword

      if(str1.compare("ENDMDL") == 0)
        {
        end = true;  
        } else
        {
        ATOM newAtom;
        input >> str2 ; // reading ID#
        input >> newAtom.type;
        input >> str4; // reading MDL
        input >> newAtom.x >> newAtom.y >> newAtom.z;
        input >> str2 >> str3 >> str4; // ignore 2 numbers and a string with label

//        if(!CAR_USE_ATOM_TYPE_OVER_NAME) newAtom.type = str1;

        Point newCoords = cell->xyz_to_abc(newAtom.x, newAtom.y, newAtom.z);
        newAtom.a_coord = newCoords[0]; newAtom.b_coord = newCoords[1]; newAtom.c_coord = newCoords[2];
        newAtom.radius = lookupRadius(newAtom.type, radial);
        cell->atoms.push_back(newAtom);
        na++;
        };

      };

    cell->numAtoms = na;
    cout << na << " atoms read." << endl;

    input.close();
  }
  else{
    cerr << "Error: PDB failed to open " << filename << endl;
//    abort();
    return false;
  }
  return true;
}


/** Read the information from the .v1 file referrred to by filename
    and store it within the provided ATOM_NETWORK. */
//void readV1File(char *filename, ATOM_NETWORK *cell, bool radial){
bool readV1File(char *filename, ATOM_NETWORK *cell, bool radial){
  fstream input;
  char garbage[256];
  int i;
  input.open(filename);
  if(!input.is_open()){
    cout<< "Failed to open .v1 file " << filename << "\n";
    cout<<"Exiting ..." << "\n";
//    exit(1);
    return false;
  }
  else{
    cout  << "Reading input file " << filename << "\n";
    input.getline(garbage, 256);

    // Read and store information about the unit cell. While the
    // vector components are known, the side lengths and angles remain unknown
    input >> garbage >> cell->v_a.x >> cell->v_a.y >> cell->v_a.z;
    input >> garbage >> cell->v_b.x >> cell->v_b.y >> cell->v_b.z;
    input >> garbage >> cell->v_c.x >> cell->v_c.y >> cell->v_c.z; 
    input >> cell->numAtoms;
    cell->initMatrices();

    //Recover information about unit cell angles and side lengths from
    //vector components. Essentially reverses the steps performed in
    //the initialize() method for ATOM_NETWORK instances as contained
    //in the file networkstorage.cc 
    cell->a = cell->v_a.x;
    cell->b = sqrt(cell->v_b.x*cell->v_b.x+cell->v_b.y*cell->v_b.y);
    cell->c = sqrt(cell->v_c.x*cell->v_c.x+cell->v_c.y*cell->v_c.y+cell->v_c.z*cell->v_c.z);
    cell->beta = acos(cell->v_c.x/cell->c)*360.0/(2.0*PI);
    cell->gamma = acos(cell->v_b.x/cell->b)*360.0/(2.0*PI);
    cell->alpha = 360.0/(2*PI)*acos((cell->v_c.y/cell->c*sin(2.0*PI*cell->gamma/360.0)) 
				    +cos(2.0*PI/360.0*cell->gamma)*cos(2.0*PI/360.0*cell->beta));

    // Read and store information about each atom. The coordinates
    // relative to the unit cell vectors remain unknown
    for(i=0; i<cell->numAtoms; i++){
      ATOM newAtom;
      input >> newAtom.type >> newAtom.x >> newAtom.y >> newAtom.z;
      Point abcCoords = cell->xyz_to_abc(newAtom.x, newAtom.y, newAtom.z);
      newAtom.a_coord = trans_to_origuc(abcCoords[0]); newAtom.b_coord = trans_to_origuc(abcCoords[1]); newAtom.c_coord = trans_to_origuc(abcCoords[2]);
      newAtom.radius = lookupRadius(newAtom.type, radial);
      cell->atoms.push_back(newAtom);
    }    
    input.close();
  }
  return true;
}


/** Read the information from the .pld file referrred to by filename
    and store it within the provided ATOM_NETWORK. 
    .dlp is a frame from DL_poly HISTORY file */
bool readDLPFile(char *filename, ATOM_NETWORK *cell, bool radial){
  fstream input;
  char garbage[256];
  int i;
  input.open(filename);
  if(!input.is_open()){
    cout<< "Failed to open .dlp file " << filename << "\n";
    cout<<"Exiting ..." << "\n";
//    exit(1);
    return false;
  }
  else{
    cout  << "Reading input file " << filename << "\n";
    input.getline(garbage, 256);

    // Read and store information about the unit cell. While the
    // vector components are known, the side lengths and angles remain unknown
    input >> cell->v_a.x >> cell->v_a.y >> cell->v_a.z;
    input >> cell->v_b.x >> cell->v_b.y >> cell->v_b.z;
    input >> cell->v_c.x >> cell->v_c.y >> cell->v_c.z;
     
//    input >> cell->numAtoms;
    cell->initMatrices();

    //Recover information about unit cell angles and side lengths from
    //vector components. Essentially reverses the steps performed in
    //the initialize() method for ATOM_NETWORK instances as contained
    //in the file networkstorage.cc 
    cell->a = cell->v_a.x;
    cell->b = sqrt(cell->v_b.x*cell->v_b.x+cell->v_b.y*cell->v_b.y);
    cell->c = sqrt(cell->v_c.x*cell->v_c.x+cell->v_c.y*cell->v_c.y+cell->v_c.z*cell->v_c.z);
    cell->beta = acos(cell->v_c.x/cell->c)*360.0/(2.0*PI);
    cell->gamma = acos(cell->v_b.x/cell->b)*360.0/(2.0*PI);
    cell->alpha = 360.0/(2*PI)*acos((cell->v_c.y/cell->c*sin(2.0*PI*cell->gamma/360.0)) 
				    +cos(2.0*PI/360.0*cell->gamma)*cos(2.0*PI/360.0*cell->beta));

    // Read and store information about each atom
    int numAtoms = 0;
    while(!input.eof()){
      ATOM newAtom;
      input >> newAtom.type;
      if(newAtom.type.empty())
	break;
	      
      input.getline(garbage,256); // reading remaining data from a line

      input >> newAtom.x >> newAtom.y >> newAtom.z;
      input.getline(garbage,256); // read end of line

      Point newCoords = cell->xyz_to_abc (newAtom.x,newAtom.y,newAtom.z);
      newAtom.a_coord = newCoords[0]; newAtom.b_coord = newCoords[1]; newAtom.c_coord = newCoords[2];
      newAtom.a_coord = trans_to_origuc(newAtom.a_coord);
      newAtom.b_coord = trans_to_origuc(newAtom.b_coord);
      newAtom.c_coord = trans_to_origuc(newAtom.c_coord);
      newAtom.radius = lookupRadius(newAtom.type, radial);
      cell->atoms.push_back(newAtom);
      numAtoms++;
    }   
    cell->numAtoms= numAtoms; 

    input.close();
  }
  return true;
}



/** Read the VORONOI_NETWORK information located in the provided input stream and 
 *  store it using the provided VORONOI_NETWORK pointer. The input stream must have a file format
 *  corresponding to a .net or .nt2 format. */
void readNet(istream *input, VORONOI_NETWORK *vornet){
  char buff [256];
  input->getline(buff,256); // Read line "Vertex table:"
  VOR_NODE node; 
  string garbage;
  
  // Read information about each Voronoi node
  int i = 0;
  while(true){
    (*input) >> garbage;
    if(strcmp(garbage.data(), "Edge") == 0)
      break;
    
    (*input) >> node.x >> node.y >> node.z >> node.rad_stat_sphere;
    
    // Read node connectivity information
    char *connectBuff = new char [256];
    char *origBuff = connectBuff;
    input->getline(connectBuff,256);
    connectBuff+= 1; //Skip space character
    char *currentChar = connectBuff;
    vector<int> nearestAtoms;
    
    while(true){
      if(*currentChar == ' ' || *currentChar == '\0'){
	char nextID[256];
	strncpy(nextID,connectBuff,currentChar-connectBuff);
	nextID[currentChar-connectBuff] = '\0';
	nearestAtoms.push_back(atoi(nextID));
	connectBuff = currentChar + 1;
      }
      if(*currentChar == '\0')
	break;
      currentChar++;
    } 

    delete [] origBuff;
    node.atomIDs = nearestAtoms;
    vornet->nodes.push_back(node);
    i++;
  }
  
  input->getline(buff,256); // Reads remainder of line "Edge table:"
  
  //Read information about each Voronoi edge
  VOR_EDGE edge;
  while(!input->eof()){
    (*input) >> edge.from >> garbage >> edge.to >> edge.rad_moving_sphere  
	     >> edge.delta_uc_x >> edge.delta_uc_y >> edge.delta_uc_z >> edge.length;
    vornet->edges.push_back(edge);
  }
  vornet->edges.pop_back();
}

/* Read the VORONOI_NETWORK located in the provided file and store 
 * it using the provided network pointer. The file must be in the .net/.nt2 
 * file format. */ 
//void readNetFile(char * filename, VORONOI_NETWORK *vornet){
bool readNetFile(char * filename, VORONOI_NETWORK *vornet){
  fstream input;
  input.open(filename);
  if(!input.is_open()){
    cout<< "Failed to open .nt2 file " << filename << "\n";
    cout<<"Exiting ..." << "\n";
//    exit(1);
    return false;
  }
  else{
    readNet(&input, vornet);
  }
  return true;
}

/** Write the information within the provided ATOM_NETWORK in a .cssr
    file format to the provided filename. */
bool writeToCSSR(char *filename, ATOM_NETWORK *cell){
  fstream output;
  output.open(filename, fstream::out);
  if(!output.is_open()){
    cerr << "Error: Failed to open .cssr output file " << filename << endl;
	//cerr << "Exiting ..." << "\n";
    //exit(1);
    return false;
  }
  else{
    cout   << "Writing atom network information to " << filename << "\n";
    
    // Write information about the unit cell
    output << "\t\t\t\t" << cell->a << "  " << cell->b << "  " << cell->c << "\n";
    output << "\t\t" << cell->alpha<<"  "<< cell->beta <<"  " << cell->gamma <<"  SPGR =  1 P 1\t\t OPT = 1" << "\n";
    output << cell->numAtoms << "   0 " << "\n";
    output << "0 " << cell->name << "\t" << ": " << cell->name << "\n";
    output.setf(ios::fixed, ios::floatfield);
    int i;
    ATOM atm;
    
    // Write information about each atom
    for(i = 0; i<cell->numAtoms; i++){
      atm = cell->atoms.at(i);
      output << " " << i+1 << " " << cell->atoms.at(i).type << " " << atm.a_coord << " "
	     << atm.b_coord << " " << atm.c_coord << "  0  0  0  0  0  0  0  0  " << atm.charge <<  "\n";
    }
    output.close();
    return true;
  }
}

/** Write the information within the provided ATOM_NETWORK in a .cssr
    file format to the provided filename, using labels instead of element types. */
bool writeToCSSRLabeled(char *filename, ATOM_NETWORK *cell){
  fstream output;
  output.open(filename, fstream::out);
  if(!output.is_open()){
    cerr << "Error: Failed to open .cssr output file " << filename << endl;
        //cerr << "Exiting ..." << "\n";
    //exit(1);
    return false;
  }
  else{
    cout   << "Writing atom network information to " << filename << "\n";

    // Write information about the unit cell
    output << "\t\t\t\t" << cell->a << "  " << cell->b << "  " << cell->c << "\n";
    output << "\t\t" << cell->alpha<<"  "<< cell->beta <<"  " << cell->gamma <<"  SPGR =  1 P 1\t\t OPT = 1" << "\n";
    output << cell->numAtoms << "   0 " << "\n";
    output << "0 " << cell->name << "\t" << ": " << cell->name << "\n";
    output.setf(ios::fixed, ios::floatfield);
    int i;
    ATOM atm;

    // Write information about each atom
    for(i = 0; i<cell->numAtoms; i++){
      atm = cell->atoms.at(i);
      output << " " << i+1 << " " << cell->atoms.at(i).label << " " << atm.a_coord << " "
             << atm.b_coord << " " << atm.c_coord << "  0  0  0  0  0  0  0  0  " << atm.charge <<  "\n";
    }
    output.close();
    return true;
  }
}

/* Computes the formula of the input ATOM_NETWORK and returns it as a c++ string */
string get_formula(ATOM_NETWORK const *const cell)
{
    vector<string> atomtypes;
    map<string,int>  typecount;

    for (vector<ATOM>::const_iterator it=cell->atoms.begin(); it!=cell->atoms.end(); ++it){
        if (find(atomtypes.begin(), atomtypes.end(), it->type) != atomtypes.end()){
            ++(typecount[it->type]);
        }
        else{
            atomtypes.push_back(it->type);
            typecount[it->type] = 1;
        }
    }
    string formula;
    for (map<string,int>::iterator it=typecount.begin(); it!=typecount.end(); ++it){
        formula.append(it->first);
        stringstream ss;
        ss << it->second;
        formula.append(ss.str());
    }
    return formula;
}

/* Generates time stamp */
string get_timestamp()
{
    char buff[80];
    time_t rawtime;
    struct tm* timeinfo;

    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buff, 80, "%F_%T", timeinfo);

    return string(buff);
}


/** Write the infomation within the provided ATOM_NETWORK in a .cif
    file format to the provided filename. **/
bool writeToCIF(char *filename,  ATOM_NETWORK *cell){
  fstream output;
  output.open(filename, fstream::out);
  if(!output.is_open()){
    cerr << "Error: Failed to open .cif output file " << filename << endl;
	//cerr << "Exiting ..." << "\n";
    //exit(1);
    return false;
  }
  else{
    cout   << "Writing atom network information to " << filename << "\n";
    //output << "data_" << origname << endl; //we need this data_(name) line otherwise the cif file will not be considered valid in some software packages
    //Instead of supplying argument to write data_..., using timestamp with number of atoms in the structure
    string formula = get_formula(cell);
    string time_stamp = get_timestamp();
    output << "data_" << formula << "_" << time_stamp << endl;
    output << "#******************************************" << endl;
    output << "#" << endl;
    output << "# CIF file created by Zeo++" << endl;
    output << "# Zeo++ is an open source package to" << endl;
    output << "# analyze microporous materials" << endl;
    output << "#" << endl;
    output << "#*******************************************" << "\n\n";
    output << "_cell_length_a\t\t" << cell->a << "   " << endl;     //  removed (0) 
    output << "_cell_length_b\t\t" << cell->b << "   " << endl;
    output << "_cell_length_c\t\t" << cell->c << "   " << endl;
    output << "_cell_angle_alpha\t\t" << cell->alpha << "   " << endl;
    output << "_cell_angle_beta\t\t" << cell->beta << "   " << endl;
    output << "_cell_angle_gamma\t\t" << cell->gamma << "   \n\n";
    output << "_symmetry_space_group_name_H-M\t\t" << "'P1'" << endl;
    output << "_symmetry_Int_Tables_number\t\t" << "1" << endl;
    output << "_symmetry_cell_setting\t\t";
    
    //Determine the Crystal System
    if (cell->alpha == 90 && cell->beta == 90 && cell->gamma == 90){
      if (cell->a == cell->b || cell->b == cell->c || cell->a == cell->c){
	      if (cell->a == cell->b && cell->b == cell->c){
	        output << "Isometric\n" << endl;
	      }
	      else {
	        output << "Tetragonal\n" << endl;
        }
      }
      else{
      	output << "Orthorhombic\n" << endl;
      }
    }
    else if(cell->alpha == cell->beta || cell->beta == cell->gamma || cell->alpha == cell->gamma){
      output << "Monoclinic\n" << endl;
    }
    else{
      output << "Triclinic\n" << endl;
    }
    
    output << "loop_" << endl;
    output << "_symmetry_equiv_pos_as_xyz" << endl;
    output << "'+x,+y,+z'\n" << endl;
    output << "loop_" << endl;
    output << "_atom_site_label" << endl;
    output << "_atom_site_type_symbol" << endl;
    output << "_atom_site_fract_x" << endl;
    output << "_atom_site_fract_y" << endl;
    output << "_atom_site_fract_z" << endl;
    for (unsigned int i=0; i<cell->atoms.size(); i++){
      ATOM *temp=&(cell->atoms.at(i));
      output << temp->specialID << "\t" << temp->type << "\t" << trans_to_origuc(temp->a_coord) << "\t" << trans_to_origuc(temp->b_coord) << "\t" << trans_to_origuc(temp->c_coord) << endl;
    }
    output.close();
    return true;
  }
}


/** Write the information within the provided ATOM_NETWORK in a .v1
    file format to the provided filename. */
bool writeToV1(char * filename, ATOM_NETWORK *cell){
  fstream output;
  output.open(filename, fstream::out);
  if(!output.is_open()){
    cerr << "Error: Failed to open .v1 output file " << filename << endl;
    //cerr << "Exiting ..." << "\n";
    //exit(1);
    return false;
  }
  else{
    cout << "Writing atom network information to " << filename << "\n";

    // Write information about the unit cell
    output << "Unit cell vectors:" << "\n";
    output.precision(8);
    output << "va= " << cell->v_a.x << " " << cell->v_a.y << " " << cell->v_a.z << "\n";
    output << "vb= " << cell->v_b.x << " " << cell->v_b.y << " " << cell->v_b.z << "\n";
    output << "vc= " << cell->v_c.x << " " << cell->v_c.y << " " << cell->v_c.z << "\n";
    output << cell->numAtoms << "\n";
    
    // Write information about each atom
    vector <ATOM> ::iterator iter = cell->atoms.begin();
    while(iter != cell->atoms.end()){
      output << iter->type << " " << iter->x << " " << iter->y << " " << iter->z << "\n";
      iter++;
    }
  }
  output.close();
  return true;
}

/** Write the information stored within the VORONOI_NETWORK in a .nt2
    file format to the provided filename. Excludes any nodes or nodes with radii
    less than the provided threshold. For the default 0, all nodes and endges are included*/
bool writeToNt2(char *filename, VORONOI_NETWORK *vornet, double minRad){
  fstream output;
  output.open(filename, fstream::out);
  if(!output.is_open()){
    cerr << "Error: Failed to open .net2 output file " << filename << "\n";
    //cerr << "Exiting ..." << "\n";
    //exit(1);
    return false;
  }
  else{
    cout << "Writing Voronoi network information to " << filename << "\n";
  
    // Write Voronoi node information
    output << "Vertex table:" << "\n";
    vector<VOR_NODE> ::iterator niter = vornet->nodes.begin();
    int i = 0;
    while(niter != vornet->nodes.end()){
      if(niter->rad_stat_sphere > minRad){
	output << i << " " << niter-> x << " " << niter-> y << " " 
	       << niter-> z << " " << niter->rad_stat_sphere;
      
	//Write Voronoi node/atom pairing information in i j k....z format
	output << " ";
	for(unsigned int j = 0; j < niter->atomIDs.size(); j++){
	  output << niter->atomIDs.at(j);
	  if(j < niter->atomIDs.size()-1)
	    output << " ";
	}
	output <<  "\n";
      }
      i++;
      niter++;
    }

    // Write Voronoi edge information
    output << "\n" << "Edge table:" << "\n";
    vector<VOR_EDGE> ::iterator eiter = vornet->edges.begin();
    while(eiter != vornet->edges.end()){
      if(eiter->rad_moving_sphere > minRad){
	output << eiter->from << " -> " << eiter->to << " " << eiter->rad_moving_sphere
	       << " " << eiter->delta_uc_x << " " << eiter->delta_uc_y << " "
	       << eiter->delta_uc_z << " " << eiter->length << "\n";
      }
      eiter++;
    }
  }
  output.close();
  return true;
}

/** Write the voronoi noide information within the VORONOI_NETWORK in .xyz 
    file format to the provided filename. Excludes any nodes with radii
    less than the provided threshold. For the default 0, all nodes are included*/
bool writeToXYZ(char *filename, VORONOI_NETWORK *vornet, double minRad){
    fstream output;
    output.open(filename, fstream::out);
    if(!output.is_open()){
        cerr << "Error: Failed to open .net2 output file " << filename << "\n";
        return false;
    }
    else{
        cout << "Writing Voronoi network information to " << filename << "\n";
        
        // Write Voronoi node information
        //Compute the # of nodes to be written
        int i = 0;
        for (vector<VOR_NODE>::const_iterator iter = vornet->nodes.begin(); 
                iter != vornet->nodes.end(); iter++)
            if (iter->rad_stat_sphere > minRad){
                i++;
            }
        output << i << "\n\n";
        for (vector<VOR_NODE>::const_iterator iter = vornet->nodes.begin(); 
                iter != vornet->nodes.end(); iter++)
            if (iter->rad_stat_sphere > minRad)
                output << "X  " << iter->x << " " << iter->y << " "
                       << iter->z << " " << iter->rad_stat_sphere << "\n";
    }
    output.close();
    return true;
}

/** Write the information stored within the VORONOI_NETWORK in a .nt2
    file format to the provided filename. Includes all nodes and edges. */
/* redudant
 bool writeToNt2(char *filename, VORONOI_NETWORK *vornet){
  return nwriteToNt2(filename, vornet, 0);
}   */

/** Write the information within the provided ATOM_NETWORK in a .xyz
    file format to the provided filename. */
//updated this function to permit supercell (2x2x2) and duplication of atoms on the unit cell perimeter, if desired
bool writeToXYZ(char *filename, ATOM_NETWORK *cell, bool is_supercell, bool is_duplicate_perimeter_atoms){
  int num_cells = SUPERCELL_SIZE; //defined in networkio.h
  if(!is_supercell) num_cells = 1;
  fstream output;
  output.open(filename, fstream::out);
  if(!output.is_open()){
    cerr << "Error: Failed to open .xyz output file " << filename << endl;
    //cerr << "Exiting ..." << "\n";
    //exit(1);
    return false;
  }
  else{
    cout << "Writing atom network information to " << filename << "\n";

    //1) construct a new vector of atoms, which contains any supercell projections and duplicates
    vector<ATOM> new_atoms;
    for(int i=0; i<cell->numAtoms; i++) {
      ATOM atom = cell->atoms.at(i);
      Point orig_abc(atom.a_coord, atom.b_coord, atom.c_coord);
      Point uc_abc = cell->shiftABCInUC(orig_abc);
      //first, determine all the supercell images of this atom as required
      for(int a=0; a<num_cells; a++) {
        for(int b=0; b<num_cells; b++) {
          for(int c=0; c<num_cells; c++) {
            atom.a_coord=uc_abc[0]+a;
            atom.b_coord=uc_abc[1]+b;
            atom.c_coord=uc_abc[2]+c;
            // Do xyz conversion later, when needed for output
            new_atoms.push_back(atom);
            //now we have this atom, be it the original or a supercell projection, find the duplicates
            if(is_duplicate_perimeter_atoms) {
              double dupe_thresh = 0.001;
              if(atom.a_coord<dupe_thresh) {
                ATOM dupe = atom;
                dupe.a_coord = atom.a_coord+num_cells;
                new_atoms.push_back(dupe);
              }
              if(atom.b_coord<dupe_thresh) {
                ATOM dupe = atom;
                dupe.b_coord = atom.b_coord+num_cells;
                new_atoms.push_back(dupe);
              }
              if(atom.c_coord<dupe_thresh) {
                ATOM dupe = atom;
                dupe.c_coord = atom.c_coord+num_cells;
                new_atoms.push_back(dupe);
              }
              if(atom.a_coord<dupe_thresh && atom.b_coord<dupe_thresh) {
                ATOM dupe = atom;
                dupe.a_coord = atom.a_coord+num_cells;
                dupe.b_coord = atom.b_coord+num_cells;
                new_atoms.push_back(dupe);
              }
              if(atom.a_coord<dupe_thresh && atom.c_coord<dupe_thresh) {
                ATOM dupe = atom;
                dupe.a_coord = atom.a_coord+num_cells;
                dupe.c_coord = atom.c_coord+num_cells;
                new_atoms.push_back(dupe);
              }
              if(atom.b_coord<dupe_thresh && atom.c_coord<dupe_thresh) {
                ATOM dupe = atom;
                dupe.b_coord = atom.b_coord+num_cells;
                dupe.c_coord = atom.c_coord+num_cells;
                new_atoms.push_back(dupe);
              }
              if(atom.a_coord<dupe_thresh && atom.b_coord<dupe_thresh && atom.c_coord<dupe_thresh) {
                ATOM dupe = atom;
                dupe.a_coord = atom.a_coord+num_cells;
                dupe.b_coord = atom.b_coord+num_cells;
                dupe.c_coord = atom.c_coord+num_cells;
                new_atoms.push_back(dupe);
              }
              if(atom.a_coord>((double)(num_cells))-dupe_thresh) {
                ATOM dupe = atom;
                dupe.a_coord = atom.a_coord-num_cells;
                new_atoms.push_back(dupe);
              }
              if(atom.b_coord>((double)(num_cells))-dupe_thresh) {
                ATOM dupe = atom;
                dupe.b_coord = atom.b_coord-num_cells;
                new_atoms.push_back(dupe);
              }
              if(atom.c_coord>((double)(num_cells))-dupe_thresh) {
                ATOM dupe = atom;
                dupe.c_coord = atom.c_coord-num_cells;
                new_atoms.push_back(dupe);
              }
              if(atom.a_coord>((double)(num_cells))-dupe_thresh && atom.b_coord>((double)(num_cells))-dupe_thresh) {
                ATOM dupe = atom;
                dupe.a_coord = atom.a_coord-num_cells;
                dupe.b_coord = atom.b_coord-num_cells;
                new_atoms.push_back(dupe);
              }
              if(atom.a_coord>((double)(num_cells))-dupe_thresh && atom.c_coord>((double)(num_cells))-dupe_thresh) {
                ATOM dupe = atom;
                dupe.a_coord = atom.a_coord-num_cells;
                dupe.c_coord = atom.c_coord-num_cells;
                new_atoms.push_back(dupe);
              }
              if(atom.b_coord>((double)(num_cells))-dupe_thresh && atom.c_coord>((double)(num_cells))-dupe_thresh) {
                ATOM dupe = atom;
                dupe.b_coord = atom.b_coord-num_cells;
                dupe.c_coord = atom.c_coord-num_cells;
                new_atoms.push_back(dupe);
              }
              if(atom.a_coord>((double)(num_cells))-dupe_thresh && atom.b_coord>((double)(num_cells))-dupe_thresh && atom.c_coord>((double)(num_cells))-dupe_thresh) {
                ATOM dupe = atom;
                dupe.a_coord = atom.a_coord-num_cells;
                dupe.b_coord = atom.b_coord-num_cells;
                dupe.c_coord = atom.c_coord-num_cells;
                new_atoms.push_back(dupe);
              }
            }
          }
        }
      }
    }

    //2) write atom data
    output << new_atoms.size() << "\n" << "\n";
    for(int i=0; i<new_atoms.size(); i++) {
      Point p = cell->abc_to_xyz(new_atoms.at(i).a_coord, new_atoms.at(i).b_coord, new_atoms.at(i).c_coord);
      output << new_atoms.at(i).type << " " << p[0] << " " << p[1] << " " << p[2] << "\n";
    }
  }
  output.close();
  return true;
}

/** Write the boundary of the unit cell expressed within the provided ATOM_NETWORK in a .vtk
    file format to the provided filename. */
bool writeToVTK(char *filename, ATOM_NETWORK *cell){
  fstream output;
  output.open(filename, fstream::out);
  if(!output.is_open()){
    cerr << "Error: Failed to open .vtk output file " << filename << endl;
    //cerr << "Exiting ..." << "\n";
    //exit(1);
    return false;
  }
  else{
    cout << "Writing unit cell information to " << filename << "\n";

    //For each corner in abc, get the coords in xyz
    vector<Point> corners; Point p;
    p = cell->abc_to_xyz(0, 0, 0); corners.push_back(p);
    p = cell->abc_to_xyz(0, 0, 1); corners.push_back(p);
    p = cell->abc_to_xyz(0, 1, 0); corners.push_back(p);
    p = cell->abc_to_xyz(0, 1, 1); corners.push_back(p);
    p = cell->abc_to_xyz(1, 0, 0); corners.push_back(p);
    p = cell->abc_to_xyz(1, 0, 1); corners.push_back(p);
    p = cell->abc_to_xyz(1, 1, 0); corners.push_back(p);
    p = cell->abc_to_xyz(1, 1, 1); corners.push_back(p);
    // Write header, information about the cell, and footer
    output << "# vtk DataFile Version 2.0\nvtk format representation of unit cell boundary\nASCII\nDATASET POLYDATA\nPOINTS 8 double\n";
    for(int i=0; i<8; i++) output << corners.at(i)[0] << " " << corners.at(i)[1] << " " << corners.at(i)[2] << "\n";
    output << "LINES 12 36\n2 0 1\n2 0 2\n2 1 3\n2 2 3\n2 4 5\n2 4 6\n2 5 7\n2 6 7\n2 0 4\n2 1 5\n2 2 6\n2 3 7\n";
  }
  output.close();
  return true;
}


/** Write the information within the provided ATOM_NETWORK in a .mop
    file format to the provided filename. */
/** In .mop file, unit cell vectors are listed as Tv and the shape of the unit cell is kept fixed **/
bool writeToMOPAC(char *filename, ATOM_NETWORK *cell, bool is_supercell){
  int num_cells = SUPERCELL_SIZE; //defined in networkio.h
  if(!is_supercell) num_cells = 1;
  fstream output;
  output.open(filename, fstream::out);
  if(!output.is_open()){
    cout << "Error: Failed to open .mop output file " << filename << endl;
    //cout << "Exiting ..." << "\n";
    //exit(1);
    return false;
  }
  else{
    cout << "Writing atom network information to " << filename << "\n";

    //Write two empty lines by default
    output << "\n" << "\n";

/*
    // Write information about the atoms
    vector<ATOM>::iterator atomIter = cell->atoms.begin();
    while(atomIter != cell->atoms.end()){
      output << atomIter->type << "  " << atomIter->x << " +1 "
             << atomIter->y << " +1 " << atomIter->z << " +1\n";
      atomIter++;
    }
*/
    for(int i=0; i<cell->numAtoms; i++) {
      //first, determine all the supercell images of this atom as required
      for(int a=0; a<num_cells; a++) {
        for(int b=0; b<num_cells; b++) {
          for(int c=0; c<num_cells; c++) {
            ATOM atom = cell->atoms.at(i);
            atom.a_coord = trans_to_origuc(atom.a_coord)+a;
            atom.b_coord = trans_to_origuc(atom.b_coord)+b;
            atom.c_coord = trans_to_origuc(atom.c_coord)+c;
            Point p = cell->abc_to_xyz(atom.a_coord, atom.b_coord, atom.c_coord);
            output << atom.type << "  " << p[0] << " +1 "
             << p[1] << " +1 " << p[2] << " +1\n";
          }
        }
      }
    }      

    //Write unit cell box information
    output << "Tv " << cell->v_a.x*num_cells << " +1 ";
    if(cell->v_a.y==0.0) { output << " 0.0 0 "; } else { output << cell->v_a.y*num_cells << " +1 ";};
    if(cell->v_a.z==0.0) { output << " 0.0 0 \n"; } else { output << cell->v_a.z*num_cells << " +1 \n";};

    output << "Tv ";
    if(cell->v_b.x==0.0) { output << " 0.0 0 "; } else { output << cell->v_b.x*num_cells << " +1 ";};
    output << cell->v_b.y*num_cells << " +1 ";    
    if(cell->v_b.z==0.0) { output << " 0.0 0 \n"; } else { output << cell->v_b.z*num_cells << " +1 \n";};

    output << "Tv ";
    if(cell->v_c.x==0.0) { output << " 0.0 0 "; } else { output << cell->v_c.x*num_cells << " +1 ";};
    if(cell->v_c.y==0.0) { output << " 0.0 0 "; } else { output << cell->v_c.y*num_cells << " +1 ";};
    output << cell->v_c.z*num_cells << " +1 \n\n";

  }
  output.close();
  return true;
}


/** Change the type of the atom to its appropriate form. Used when
 reading .cuc files. */
void changeAtomType(ATOM *atom){
    switch(atom->type[0]){
        case 'O': case 'o':
            atom->type = "O";
            break;
        case 'T': case 't':
            atom->type = "Si";
            break;
        case 'A': case 'a':
            atom->type = "Si";
            break;
        case 'H': case 'h':
            atom->type = "H";
            break;
        case 'S': case 's':
            if(tolower(atom->type[1]) == 'i')
                atom->type = "Si";
            else
                atom->type = "S";
            break;
        default:
            cerr<< "Error: Atom name not recognized " << atom->type << "\n";
            //	<< "Exiting..." << "\n";
            //    exit(1);
            break;
    }
}



/** Strip atom names from additional indexes following atom string */
void stripAtomNames(ATOM_NETWORK *cell)
{
 for(unsigned int na = 0; na < cell->atoms.size(); na++)
  cell->atoms[na].type = stripAtomName(cell->atoms[na].type);
// for (vector<ATOM>::const_iterator it=cell->atoms.begin(); it!=cell->atoms.end(); ++it){
//   it->type = stripAtomName(it->type);
//   cout << stripAtomName(it->type) << "\n";
//    }
}
