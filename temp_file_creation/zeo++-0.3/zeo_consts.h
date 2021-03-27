//Define constants for Zeo++
 
#ifndef ZEOCONSTS_H
#define ZEOCONSTS_H


const double PI=3.14159265358973;
const double AVOGRADOS_NUMBER = 6.0221415e23;
const double TOLERANCE=0.0001;
const double threshold=0.0000001;
const double thresholdLarge=0.001;  // Large threshold used to detect special positions when deadling with CIF files
const double VOR_NODE_MERGE_THRESHOLD=0.02; // Threshold used to merge nearby nodes in Voronoi network
                                            // this greatly simplifies the network
const int maxline = 500; //for char arrays

// MC Sampling parameters
const int MCAVppA3 = 200; // by default use 200 samples per cubic A for volume calcs
const int MCASAppA2 = 50; // by default use 50 samples per square A for ASA calcs
                          // our tests show that these give convergance within 0.2 %

inline double SQR(double x) {
     return x*x;
}
inline int SQR(int x) {
     return x*x;
}

//Options for ray.cc
#define DEBUG  false
const int MAXRAYDIST=100;
const float BINSIZE=0.1;
const int MAXBINS=1000;

 //#define SQR(x) ((x)*(x))
 //
#endif
