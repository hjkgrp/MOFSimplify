//#include "network.h"
#include <vector>

#include "symbcalc.h"
#include "string_additions.h"

using namespace std;

/** This function is used to evaluate a string expresion
It understands /,+,-,*,x,y,z,and 0-9,. the calculator returns
the numerical result as a double**/
float symbCalc(string expresion,float x,float y,float z){
  float val=0;
  vector<string> numbers = split(expresion,"+-/*");
  
  //This will replace all x,y,z with there "value"
  for (unsigned int i=0;i<numbers.size();i++){
    if (numbers[i].compare("x")==0){
      numbers[i] = doubleToString(x);
    }
    else if (numbers[i].compare("y")==0){
      numbers[i] = doubleToString(y);
    }
    else if (numbers[i].compare("z")==0){
      numbers[i] = doubleToString(z);
    }
  }

  vector<double> nums;
  for (unsigned int i=0;i<numbers.size();i++){
    nums.push_back(convertToDouble(numbers[i]));
  }
  
  vector<string> opperator = split(expresion,"1234567890.xyz");
  
  //To get rid of the +- before the expresion
  if (opperator.size() == numbers.size()){
    if (opperator[0].compare("-")==0){
      nums[0]=-nums[0]; //negate the element
    }
    opperator.erase(opperator.begin());
  }
  
  //deal with - signs
  for (unsigned int i=0;i<opperator.size();i++){
    if (opperator[i].compare("-")==0){
      nums[i+1]=-nums[i+1]; //negate the element
    }
  }
  
  //deal with * signs
  for (unsigned int i=0;i<opperator.size();i++){
    if (opperator[i].compare("*")==0){
      nums[i+1]=nums[i]*nums[i+1];
      opperator.erase(opperator.begin()+i);
      nums.erase(nums.begin()+i); //I chose i+1 becuase i stored value in i
    }
  }

  //deal with / signs
  for (unsigned int i=0;i<opperator.size();i++){
    if (opperator[i].compare("/")==0){
      nums[i]=nums[i]/nums[i+1];
      opperator.erase(opperator.begin()+i);
      nums.erase(nums.begin()+i+1); //I chose i+1 becuase i stored value in i  
    }
  }

  //since only + is left lets add up all the values!!
  for (unsigned int i=0;i<nums.size();i++){
    val = val + nums[i];
  }
  
  return val;
}

