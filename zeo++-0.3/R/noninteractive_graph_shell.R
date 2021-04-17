#Pore size distribution graphs: main source code, shell version
#V.1.0 March 6, 2013
#Author: Marielle Pinheiro

#	To run this file: R --slave --args (directory path, input file name, 
#	output file name, custom parameters (Y/N))

#Load functions
source("noninteractive_graph_functions.R")

##Initial plot parameters 
## These get overwritten by program to optimized limits
## For customized parameters, change desired variables in
## parameters file
start_val = 1
stop_val = 3
main = "PSD Distribution"
xmin = 0
ymin = 0
xmax = 16
ymax = 3
xlab = "Diameter (A)"
ylab = "PSD"

#1. Initialize with command line arguments: directory path, file path/ pattern, 
#		output file/ path, custom parameters (Y/N), single or batch (S/B)
args <- commandArgs(TRUE)
if (length(args) != 4){
	cat("Wrong number of arguments. Arguments should be: Directory path, input file, output file, customize (Y/N).\n")
} else {
	dir_path = args[1]
	name_path = args[2]
	output_name = args[3]
	custom_par = args[4]

	if (custom_par == "Y" || custom_par == "y"){
		source("parameters.R")
	}

	#2. Set working directory to where your files are located
	setwd(dir_path)
	#3. Generate plot (if valid, otherwise returns error)
	gen_single(name_path, output_name, custom_par)
}
