#Pore size distribution graphs, default functions
#V.1.0 March 6, 2013
#Author: Marielle Pinheiro

#This function loads the file into R
read_files <-function(file, name)
{
	## skips first 5 lines of text
	read = read.table(file, skip = 5)
	if (is.nan(read[1,3]) == TRUE) {
		cat("Cannot plot function; distribution is 0.\n")
		return(FALSE)
	} else {	##Assigns variable names for your distributions (AFT, GIS, etc)
		assign(name, read, envir = .GlobalEnv)
		return(TRUE)
	}
}

#This function finds the highest diameter for the plot (plus a couple extra bins)
get_length <-function(variable){
	new_distr = which(variable$V3!=0)
	newlength = length(new_distr)+3
	xmax = newlength/10.0
	return(xmax)
}

#This function finds the maximum peak height for the plot
get_height <-function(variable){
	ymax = max(variable$V4) + 0.2
	return(ymax)
}

#This function plots a single file
plot_distr <-function(variable, parameters, output_name)
{
	X11()
	plot(variable$V1, variable$V4, type = "l", main = parameters[1], xlim = as.numeric(c(parameters[2], parameters[3])), ylim = as.numeric(c(parameters[4], parameters[5])), xlab = parameters[6], ylab = parameters[7])
	dev.copy2eps(file = output_name)
	dev.off()
}

#This function is the entire single plot function
gen_single <- function(name_path, output_name, custom_par)
{
	#1. Load file
	check_true <- TRUE
	var_name = substr(name_path, start= start_val, stop= stop_val)
	check_true = read_files(name_path, var_name)
	if (check_true == TRUE){
		variable = get(var_name)
		if (custom_par == "N" || custom_par == "n"){ 
			main = paste(var_name, "PSD Distribution")
			xmax = get_length(variable)
			ymax = get_height(variable)
		}
		#2. Generate single plot
		plot_parameters = c(main, xmin, xmax, ymin, ymax, xlab, ylab)
		plot_distr(variable, plot_parameters, output_name)
	}
}

