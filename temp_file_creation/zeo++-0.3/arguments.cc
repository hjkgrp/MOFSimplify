/* Functions that process input arguments */

#include "network.h"
#include "networkinfo.h"
#include "arguments.h"

//using namespace std;



string processFilename(vector<string> command, string name, string file_type, unsigned int num_args_1, unsigned int num_args_2){
    if(command.size() == (num_args_1 + 1)){
        return string(name).append(file_type);
    }
    else if (command.size() == (num_args_2 + 1)){
        return string(command[num_args_2]);
    }
    else {
        cerr << "Error: " << command[0] << " option accepts " << num_args_1 << " or " << num_args_2
        << " arguments but " << command.size() - 1 << " arguments were supplied. " << "\n"
        << "Exiting..." << "\n";
        //    exit(1);
        return "";
    }
}


/** Assuming that -r was part of the command line parameters, this
 function examines the list of arguments to determine whether
 default radii or radii from a file will be used. If a file was
 provided, it reads the file and stores them within the radii table
 for future use. */
void processRadialParameters(vector<string> args){
    cout << "Loading radii. "<< "\n";
    if(args.size() == 1){
        //do nothing - radTable will have already been initialised by default
    }
    else if(args.size() == 2){
        if(args[1].find(".rad") == string::npos){
            cerr << "Error: If an argument is provided to -r, it must be a file of type .rad containing the types and radii of each atom" << "\n"
            << "Exiting..." << "\n";
            exit(1);
        }
        else{
            readRadTable((char *)args[1].c_str());
        }
    }
    else{
        cerr << "Error: -r option accepts at most 1 argument but " << args.size() - 1 << " were provided." << "\n"
        << "Exiting..." << "\n";
        exit(1);
    }
}

/** Assuming that -mass was part of the command line parameters, this
 function examines the list of arguments to determine whether
 default atomic masses or masses from a file will be used. If a file was
 provided, it reads the file and stores them within the mass table
 for future use. */
void processMassParameters(std::vector<std::string> args){
    if(args.size() == 1){
        //do nothing - massTable will have already been initialised by default
    }
    else if(args.size() == 2){
        if(args[1].find(".mass") == string::npos){
            cerr << "Error: If an argument is provided to -mass, it must be a file of type .mass containing the types and masses of each atom" << "\n"
            << "Exiting..." << "\n";
            exit(1);
        }
        else{
            readMassTable((char *)args[1].c_str());
        }
    }
    else{
        cerr << "Error: -mass option requires 1 argument but " << (args.size()-1) << " were supplied" << "\n"
        << "Exiting..." << "\n";
        exit(1);
    };
    
} // ends processMassParameters()



/** Assuming that -ha was part of the command line parameters, this
 function examines the list of arguments to determine the requested
 accuracy setting. If no parameters were specified, default settings
 are assumed. */
string processAccuracyParameters(vector<string> args){
    string accSetting;
    cout << "High accuracy requested (DEF (default) settings). "<< "\n";
    if(args.size() == 1){
        accSetting = "DEF";
    }
    else if(args.size() == 2){
        
        string tempSetting = toUpperCase(args[1]);
        
        if(tempSetting == "OCC" ||
           tempSetting == "FCC" ||
           tempSetting == "ACC" ||
           tempSetting == "AQC" ||
           tempSetting == "DDH" ||
           tempSetting == "TIH" ||
           tempSetting == "ICH" ||
           tempSetting == "ICC" ||
           tempSetting == "RIH" ||
           tempSetting == "S4" ||
           tempSetting == "S10" ||
           tempSetting == "S20" ||
           tempSetting == "S30" ||
           tempSetting == "S40" ||
           tempSetting == "S50" ||
           tempSetting == "S100" ||
           tempSetting == "S500" ||
           tempSetting == "S1000" ||
           tempSetting == "S10000" ||
           tempSetting == "DEF" ||
           tempSetting == "HI" ||
           tempSetting == "MED" ||
           tempSetting == "LOW") { cout << "Selected accuracy setting: " << tempSetting << "\n"; accSetting = tempSetting;}
        else{
            cerr << "Error: -ha option of " << tempSetting << " is not recognized.\n"
            << "Available options: FCC ACC AQC DDH TIH ICH ICC RIH S4 S10 S20 S30 S40 S50 S100 S500 S1000 S10000 HI MED LOW DEF\n"
            << "Exiting..." << "\n";
            exit(1);
            
        };
        
    }
    else{
        cerr << "Error: -ha option accepts at most 1 argument but " << args.size() - 1 << " were provided." << "\n"
        << "Exiting..." << "\n";
        exit(1);
    }
    
    return accSetting;
    
}

/* Analyze visualization settings for writing output files */
string processVisualizationParameters(vector<string> args){
    string visSetting;
    cout << "Additional files for visualization requested."<< "\n";
    if(args.size() == 1){
        visSetting = "VISIT";
    }
    else if(args.size() == 2){

        string tempSetting = toUpperCase(args[1]);

        if(tempSetting == "VISIT" ||
           tempSetting == "LIVERPOOL" ||
           tempSetting == "ZEOVIS" ||
           tempSetting == "FRAC" ||
           tempSetting == "CART" ||
           tempSetting == "C" ||
           tempSetting == "F" ||
           tempSetting == "L" ||
           tempSetting == "LIV" ) { cout << "Selected visualization setting: " << tempSetting << "\n"; visSetting = tempSetting;}
        else{
            cerr << "Error: -vo (-visual) option of " << tempSetting << " is not recognized.\n"
            << "Available options: VISIT, CART, C (Caart coord.) / FRAC, L, LIV, LIVERPOOL (frac. coord) / ZEOVIS \n"
            << "Exiting..." << "\n";
            exit(1);

        };
    if(tempSetting == "CART" ||tempSetting == "C") visSetting = "VISIT";
    if(tempSetting == "L" ||tempSetting == "F" ||tempSetting == "FRAC" ||tempSetting == "LIV" ) visSetting = "LIVERPOOL";
    }
    else{
        cerr << "Error: -vo (-visual) option accepts at most 1 argument but " << args.size() - 1 << " were provided." << "\n"
        << "Exiting..." << "\n";
        exit(1);
    }

    return visSetting;

}
