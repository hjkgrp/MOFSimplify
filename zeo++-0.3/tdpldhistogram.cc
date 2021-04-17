/* 
 *  Time-dependent PLD and Histogram tool
 *
 *  added by M. Haranczyk
 *           May 2014
 *
 *  It reads a list of .pld files from Zeo++ and perform analysis
 *
 *
 */


#include <stdlib.h>
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

double bin_size = 0.1;

int NBINS = 201;

int n_segments;

vector<int> DiType; // stores a vector with pore type ID
vector<double> DiTable;  // vector storing Di for each segment
vector< vector<double> > PLDTable; // vector storing PLD table (Df/PLD for each segment pair)

vector< vector<double> > TDDiTable;  // vector storing frames of Di information (Di for each segment)
vector< vector< vector<double> > > TDPLDTable; // vector storing frames of PLD tables (Df/PLD for each segment pair)

vector< vector<int> > PLDActiveTable; // highlights connectivity between segments to be monitored


/* Calculates histogram based on vector */
vector<int> histogram(vector <double> valuelist){

 vector<int> histogram;
 histogram.resize(NBINS,0);

 for(unsigned int i=0; i<valuelist.size(); i++)
    {

    histogram[ int(floor(valuelist[i] / bin_size))]++;

    };

 /* histogram sanity check */

 int nhist=0;
 for(unsigned int k=0;k<NBINS;k++) nhist=nhist+histogram[k];
 if(nhist!=valuelist.size())
   {
   cerr << "Histogram sanity check failed, exit\n";
   abort();
   }; 

 return histogram;
}

/* Calculates Mean and Std. Dev. based on distribution */
pair<double,double> normal_distr_parameters(vector <double> valuelist){

 pair<double,double> mean_stddev;

 mean_stddev.first = 0.0;
 for(unsigned int i=0; i<valuelist.size(); i++)
    {
    mean_stddev.first += valuelist[i];
    };
 mean_stddev.first = mean_stddev.first / (double) valuelist.size();
 
 mean_stddev.second = 0.0;

 for(unsigned int i=0; i<valuelist.size(); i++)
    {
    mean_stddev.second += ( (valuelist[i] - mean_stddev.first)*(valuelist[i] - mean_stddev.first) );
    }; 
 mean_stddev.second = mean_stddev.second /(double) valuelist.size();

 mean_stddev.second = sqrt(mean_stddev.second);

 return mean_stddev;

}

int main(int argc, char * argv[]){

 if(argc != 4){

 cout << "This tool requires 2 arguments: filename1 filename2 filename3\n" 
      << "    filename1 is .pld file that defines topology used in analysis\n"
      << "    filename2 is a file with a list of .pld files to be analyzed\n"
      << "    filename3 is one specifying types of cages (5th column with type id)"
      << "    \n";

// cout << "FORGET ABOUT BINSIZE IN THIS VERIOSN\n"; // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ override

 return 0;

 } else
 { // assume correct input parameters have been provided
   // proceed with analysis
 string summaryfilename, listfilename, poreidfile;

 summaryfilename = argv[1];
 listfilename = argv[2];
 poreidfile = argv[3];

// bin_size = atof(argv[3]); 
 bin_size = 0.1;           // @@@@@@@@@@@@@@@@@2  overridding to match NBINS

 cout << "Proceeding with arguments: " << summaryfilename << " " << listfilename << " " << bin_size << "\n";

 // First open summaryfilename to setup required parameters

 fstream sumfile;
 sumfile.open(summaryfilename.c_str());
 if(sumfile.is_open()==false){
   cerr << "Error: A file with summary/topology  (" << summaryfilename << ") failed to open. \n" ;
   return 1;
   }else
   {
   string buffer;
   sumfile >> n_segments;
   getline(sumfile, buffer);
   cout << "Expecting " << n_segments << " segments.\n";

   DiTable.resize(n_segments, -1);
   vector<int> vinttemp;
   vinttemp.resize(n_segments, -1);
   for(int i=0; i<n_segments; i++) 
     {
     PLDTable.push_back(DiTable);
     PLDActiveTable.push_back(vinttemp);
     };

   /* Read in pore id */

   fstream poreidf;
   poreidf.open(poreidfile.c_str());
   if(poreidf.is_open()==false)
     {
     cerr << "poreid file error\n";
     return 1;
     }
     else{
     int pid;
     for(int i=0; i<n_segments; i++)
        {
        poreidf >> pid;
        DiType.push_back(pid);
        };
     };

   getline(sumfile, buffer); // skipping sgment Di line
   getline(sumfile, buffer); // skipping divider line

   for(int i=0; i<n_segments; i++)
     {
     for(int j=0;j<n_segments; j++)
       {
       double d;
       sumfile >> d;
       if(d>1.0) { PLDActiveTable.at(i).at(j) = 1;} else { PLDActiveTable.at(i).at(j) = 0;};
       };
     getline(sumfile, buffer); // reading end of line
     };   
  
   };

 /* double check connectivity - pore type connections */

/*
   for(int i=0; i<n_segments; i++)
     {
     for(int j=i+1;j<n_segments; j++)
       {
       if(PLDActiveTable.at(i).at(j) == 1)
         {
         if(DiType[i]==DiType[j])
           {
           cerr << "Segments " << i << "  " << j << " of the same type are connected.\nEXIT\n";
           return 1;
           };
         }
       };
     };

  */


 /* loading data */

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
//      nfiles--;
      break;
      };

    nfiles++;
 
    fstream input;
    input.open(inputfilename.c_str());
    if(input.is_open()==false){
      cerr << "Error: .pld failed to open " << inputfilename << endl;
      }else{

      // reading frame in .pld format
      string buffer;
      int ns;
      input >> ns;
      if(ns !=n_segments)
        {
        cerr << "Number of segments in this frame is wrong\n";
        return 1;
        };
     
      getline(input, buffer);

      input >> buffer; // skipping a string
      for(int i=0; i<n_segments; i++)
        {
        input >> DiTable[i];
        };
 
      TDDiTable.push_back(DiTable);

      getline(input, buffer); // skipping sgment Di line
      getline(input, buffer); // skipping divider line

      for(int i=0; i<n_segments; i++)
        {
        for(int j=0;j<n_segments; j++)
          {
          input >> PLDTable.at(i).at(j);
          };
        getline(input, buffer); // reading end of line
        };

      TDPLDTable.push_back(PLDTable);      

      input.close();
      }; // ends input.is_open() for reading-in the current frame
   }; // ends while loop over a file list

  cout << nfiles << " frames loaded.\n";









  // ANALYSIS BLOCK COMES HERE
  vector <double> di1values;
//  vector <double> di2values;
  vector <double> dfvalues;

  vector <double> dfnonconvalues; // stores unusual connections

  vector <double> trackdi1, trackdi2, trackdf, trackdfall;

  for(unsigned int k=0; k<TDDiTable.size(); k++)
     {
     for(int i=0; i<n_segments; i++)
        {
//        cout << "k and i= " << k << "  " << i << "\n";

        //tracking
//        if(DiType[i]==204) di1values.push_back(TDDiTable[k].at(i));
//        if(DiType[i]==168) di2values.push_back(TDDiTable[k].at(i));

        di1values.push_back(TDDiTable[k].at(i));

        for(int j=i+1; j<n_segments; j++)
           {
           if(PLDActiveTable.at(i).at(j)==1) 
             {
             dfvalues.push_back(TDPLDTable[k].at(i).at(j));

//             if(i==0&&j==8) trackdf.push_back(TDPLDTable[k].at(i).at(j));
//             if(i==0) trackdfall.push_back(TDPLDTable[k].at(i).at(j));

             } else
             {
             if(TDPLDTable[k].at(i).at(j)>0)
                {
//                cerr << "Unusual segment connection detected " << i << "  " << j << " with diameter " << TDPLDTable[k].at(i).at(j) << "\n"; 
                dfnonconvalues.push_back(TDPLDTable[k].at(i).at(j));
                };
             };
           
           };
        };

     };

  cout << "Sizes of vectors: divalue, dfvalues and dfnonconvalues " << di1values.size() 
       << "     " << dfvalues.size() << "   " << dfnonconvalues.size() << "\n";

  vector<int> di1hist = histogram(di1values);
//  vector<int> di2hist = histogram(di2values);

  vector<int> dfhist = histogram(dfvalues);

  vector<int> dfnonconhist = histogram(dfnonconvalues);

  for(unsigned int i=0;i<di1hist.size();i++)
    {
    if(i==0)
      {
      di1hist[i]=di1values.size();
//      di2hist[i]=di2values.size();
      dfhist[i]=dfvalues.size();
      dfnonconhist[i]=dfnonconvalues.size();
      } else
      {
      di1hist[0]-=di1hist[i];
//      di2hist[0]-=di2hist[i];
      dfhist[0]-=dfhist[i];
      dfnonconhist[0]-=dfnonconhist[i];

      };
    };

  cout << "\n\nDi1 and  Dfi histogram, and normalized Di1, Df and Dfnoncon histogram\n";

  for(unsigned int i=0;i<di1hist.size();i++)
    {
    if(i<81)
    cout << i*bin_size << "   " << di1hist[i]  << "    "<< dfhist[i] << "       ";
    cout << double(di1hist[i])/double(di1values.size()) << "    " 
         << double(dfhist[i])/double(dfvalues.size()) << "    "
         << double(dfnonconhist[i])/double(dfvalues.size()) << "\n";
    };


  cout << "\n\n Gaussian fitting\n";

  pair<double,double> DiG = normal_distr_parameters(di1values);
  pair<double,double> DfG = normal_distr_parameters(dfvalues);

  cout << "Di: mean = " << DiG.first << "    std.dev.= " << DiG.second << "\n";
  cout << "Df: mean = " << DfG.first << "    std.dev.= " << DfG.second << "\n";

/*

  cout << "\n\n\nTD data\nFrame Di1 Di2  Df_D1-D2    Df1 Df2 Df3 Df4 (for cage D1)\n\n";

  for(unsigned int i=0; i<trackdi1.size(); i++)
     {
     cout << i <<  " "<< trackdi1[i] << "  " << trackdi2[i] << "  " <<trackdf[i] << "  " <<trackdfall[4*i] << "  "<<trackdfall[4*i+1] << "  "<<trackdfall[4*i+2] << "  "<<trackdfall[4*i+3] << "\n ";

     };

  cout << "\n\nHistogram for Di1 and Df_D1-D2\n\n";

  dihist = histogram(trackdi1);
  dfhist = histogram(trackdf);
  for(unsigned int i=0;i<dihist.size();i++)
    {
    if(i<81)
    cout << i*bin_size << "   " << dihist[i] << "    " << dfhist[i];
    cout << "      " << double(dihist[i])/double(trackdi1.size()) << "    " << double(dfhist[i])/double(trackdf.size()) << "\n";
    };

*/

  }; // ends block corresponding to successful opening of file with frame list

 }; // ends if(argc!=4

}









