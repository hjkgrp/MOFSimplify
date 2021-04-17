/* Functions that help process input arguments (moved from main) */

#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#include <vector>
#include <string>

using namespace std;

string processFilename(vector<string> command, string name, string file_type, unsigned int num_args_1, unsigned int num_args_2);
void processRadialParameters(vector<string> args);
void processMassParameters(std::vector<std::string> args);
string processAccuracyParameters(vector<string> args);
string processVisualizationParameters(vector<string> args);

#endif
