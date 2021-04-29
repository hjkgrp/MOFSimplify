This folder contains all of the files necessary to reproduce the thermal stability regressor temperatures. 

The file "TGA_info_log.xlsx" from the previous folder contains information about the compound name within the manuscript, in addition to where the TGA trace is within the manuscript. After obtaining the TGA figure, we used the webplotdigitizer tool to obtain the thermal collapse temperature, as determined by manual inspection. 

Each MOF (as identified by a refcode) that was mapped to a TGA trace, contains extracted data. This data is reported inside of the "digitized_data/" folder. For each MOF (labeled by its refcode), four points are stored: two points at the start of MOF decomposition, and two points at the end of MOF decomposition. Thus, the x-axis is temperature (in degrees Celsius), and the y-axis is mass lost (either absolute or relative). The first two points make one line, and the second two points make a second line. These two lines cross each other at a given temperature, which is the onset temperature. We identify this temperature by fitting the two sets of two points to two distinct lines and finding the x-coordinate of their intersection.

Thus, each CSV in the digitized_data/ folder serves as the raw data used to extract the collapse temperatures.
