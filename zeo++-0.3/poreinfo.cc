/** 
  * In poreinfo.* files, there are functions to perform time-series analysis for
  * void space data contain in a set of .poreinfo summary files
  *
  * Initial code by M. Haranczyk, June 2013
  *
  **/


//#include "network.h"
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>

#include "grid.h"
#include "geometry.h"
#include "networkinfo.h"
#include "poreinfo.h"

using namespace std;


/* Function that performs analysis for a set of .poreinfo files listed in listfilename.
 * Output is written to outputfilename. Atom_network has to be provided to enable distance 
   calculations */

void  analyzePoreInfoFiles(ATOM_NETWORK *atmnet, std::string listfilename, std::string outputfilename){

 vector< vector<POREINFO> > PoreInfoFrames;

 fstream flist;
 int nfiles = 1;

 flist.open(listfilename.c_str());
 if(flist.is_open()==false){
   cerr << "Error: A file with .poreinfo frames (" << listfilename << ") failed to open. \n" ;
   }else{
   // file with a list of frames is open. proceed.
   while(!flist.eof())
    {
    string inputfilename;

    flist >> inputfilename;
    if(flist.eof()){
      nfiles--;
      break;
      };

    nfiles++;

    loadPoreInfoFile(&PoreInfoFrames, inputfilename);

    cout << "File " << inputfilename  << " read." << "\n";

    };

   flist.close();

   }; // ends while loop over a file list

 cout << nfiles << " frames loaded.\n";


} // ends analyzePoreInfoFiles()



/* Function loads .poreinfo file into a vector of poreinfo frames 

 */
void loadPoreInfoFile(vector < vector<POREINFO> > *poreinfoframes, std::string poreinfofilename){

 vector <POREINFO> PoreInfoFrame;
 int nPores, nChannels, nPockets;
 int n;
 string temp,temp0,temp1,temp2;


 fstream input;
 input.open(poreinfofilename.c_str());
 if(input.is_open()==false){
   cerr << "Error: .poreinfo failed to open: " << poreinfofilename << endl;
   }else{

   input >> temp >> temp0 >> nPores >> temp1 >> nChannels >> temp2 >> nPockets; getline(input,temp);

   PoreInfoFrame.resize(nPores);

   getline(input,temp); // skipping line with AV data

   input >> temp >> n;
   if(n != nChannels)
      {
      cerr << "Number of AV channels if different than -chan channels, consider running with -ha; " << n << "  " << nChannels << "\n";
      };

   if(nChannels == 0 )
     {
     } else
     {
     input >> temp2;
     for(int i = 0; i < nChannels; i++)
       {
       input >> PoreInfoFrame.at(i).vol;
       PoreInfoFrame.at(i).dim = 1;
       };
     };
   getline(input,temp); // reading end of line character

   input >> temp >> n;
   if(n != nPockets)
      {
      cerr << "Number of AV pockets if different than -chan pockets, consider running with -ha; " << n << "  " << nPockets << "\n";
      };

   if(nPockets == 0 )
     {
     } else
     {
     input >> temp2;
     for(int i = 0; i < nPockets; i++)
       {
       input >> PoreInfoFrame.at(i+nChannels).vol;
       PoreInfoFrame.at(i+nChannels).dim = 0;
       };
     };
   getline(input,temp); // reading end of line character

   getline(input, temp); // skipping line with ASA data

   input >> temp >> n;
   if(n != nChannels)
      {
      cerr << "Number of ASA channels if different than -chan channels, consider running with -ha;" << n << "  " << nChannels << "\n";
      };

   if(nChannels == 0 )
     {
     } else
     {
     input >> temp2;
     for(int i = 0; i < nChannels; i++)
       {
       input >> PoreInfoFrame.at(i).sa;
       };
     };
   getline(input,temp); // reading end of line character

   input >> temp >> n;
   if(n != nPockets)
      {
      cerr << "Number of ASA pockets if different than -chan pockets, consider running with -ha; " << n << "  " << nPockets << "\n";
      };

   if(nPockets == 0 )
     {
     } else
     {
     input >> temp2;
     for(int i = 0; i < nPockets; i++)
       {
       input >> PoreInfoFrame.at(i+nChannels).sa;
       };
     };
   getline(input,temp); // reading end of line character

   // channel and pocket AV and ASA data read
   // proceed to reading node information

   for(int i = 0; i < nChannels; i++)
     {
     int nNodes; NODESPHERE ns;
     input >> nNodes >> PoreInfoFrame.at(i).di >> PoreInfoFrame.at(i).pos[0] >> PoreInfoFrame.at(i).pos[1] 
           >> PoreInfoFrame.at(i).pos[2] >> PoreInfoFrame.at(i).enc_radius;
     for(int j = 0; j < nNodes; j++)
        {
        input >> ns.a >> ns.b >> ns.c >> ns.r;
        };
     PoreInfoFrame.at(i).nodes.push_back(ns);
     };

   for(int i = 0; i < nPockets; i++)
     {
     int nNodes; NODESPHERE ns;
     input >> nNodes >> PoreInfoFrame.at(i+nChannels).di >> PoreInfoFrame.at(i+nChannels).pos[0] >> PoreInfoFrame.at(i+nChannels).pos[1]
           >> PoreInfoFrame.at(i+nChannels).pos[2] >> PoreInfoFrame.at(i+nChannels).enc_radius;
     for(int j = 0; j < nNodes; j++)
        {
        input >> ns.a >> ns.b >> ns.c >> ns.r;
        };
     PoreInfoFrame.at(i+nChannels).nodes.push_back(ns);
     };

   }; // ends reading in the file

 input.close();

 poreinfoframes->push_back(PoreInfoFrame);


} // ends loadPoreInfoFile


