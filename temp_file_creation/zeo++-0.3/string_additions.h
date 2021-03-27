#ifndef STRING_ADDITIONS_H
#define STRING_ADDITIONS_H

#include <string>
#include <vector>


/** Compares string with a vector of strings. Returns
	the index if it finds a match, otherwise -1 **/
int strCmpList(std::vector<std::string> list,std::string str);

/** Function that will split a string based on the 
	delimiters used if there are two delimeters in a 
	row it will skip them **/
std::vector<std::string> split(std::string line, std::string delimeter);

/** While this is messy it was messy to convert an array
	strings to a vector of strings **/
std::vector<std::string> strAry2StrVec(std::string list[]);

/** Function takes in a string and returns a double **/
double convertToDouble(std::string const& str);

/** Function takes in a string and returns an int **/
int convertToInt(std::string const& str);

/** Function takes in a double and returns a string **/
std::string doubleToString(float const& dbl);

//convert int to string
std::string intAsString(int number);

#endif
