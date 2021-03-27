//#include "network.h"
#include <iostream>
#include <sstream>
#include <cstdlib>

#include "string_additions.h"

using namespace std;

/** Compares string with a vector of strings. Returns
	the index if it finds a match, otherwise -1 **/
int strCmpList(vector<string> list,string str){
  for (unsigned int ndx=0;ndx<list.size();ndx++){
    if(list[ndx].compare(str) == 0){
      return ndx;
    }
  }
  return -1;
}

/** Function that will split a string based on the 
	delimiters used if there are two delimeters in a 
	row it will skip them **/
vector<string> split(string line, string delimeter){
  vector<string> token;
  string temp=line;
  int ndx;

  while (!temp.empty()){
    ndx = temp.find_first_of(delimeter);
    if (ndx>0){
      token.push_back(temp.substr(0,ndx));
    }
    else if (ndx ==-1){
      token.push_back(temp);
      return token;
    }
    temp=temp.substr(ndx+1);
  }

  return token;
}

/** While this is messy it was messy to convert an array
	strings to a vector of strings **/
vector<string> strAry2StrVec(string list[]){
  vector<string> veclist;
  int ndx=0;
  while (list[ndx]!="NULL"){
    veclist.push_back(list[ndx]);
    ndx++;
  }
  return veclist;
}

/** Function takes in a string and returns a double **/
double convertToDouble(string const& str){
  istringstream i(str);
  double x;
  if (!(i >> x)){
    cout << "Bad string to double conversion" << endl;
    exit(0);
  }
  return x;
}

/** Function takes in a string and returns an int **/
int convertToInt(string const& str){
  istringstream i(str);
  int x;
  if (!(i >> x)){
    cout << "Bad string to int conversion" << endl;
    exit(0);
  }
  return x;
}

/** Function takes in a double and returns a string **/
string doubleToString(float const& dbl){
  ostringstream buffer;
  if (!(buffer << dbl)){
    cout << "Bad double to string conversion" << endl;
    exit(0);
  }
  return buffer.str();
}

//convert int to string
string intAsString(int number) {
   std::ostringstream sin;
   sin << number;
   std::string val = sin.str();
   return val;
}

