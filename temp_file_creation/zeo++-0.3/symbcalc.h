#ifndef SYMBCALC_H
#define SYMBCALC_H

#include <string>

/** This function is used to evaluate a string expresion
It understands /,+,-,*,x,y,z,and 0-9 the calculator returns
the numerical result as a double**/
float symbCalc(std::string expresion,float x,float y,float z);

#endif
