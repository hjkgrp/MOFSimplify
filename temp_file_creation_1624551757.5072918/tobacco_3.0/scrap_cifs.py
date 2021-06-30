# Code written by Gianmarco to get the names of all the cif files in a given database folder
# for use in populating the dropdowns in the building block section of mofSimplify
	# just copy and paste the terminal output into the appropriate section of index.html

import glob
import os


os.chdir('edges_database')
edges_array = []
for filepath in sorted(glob.iglob('*.cif')): # sorted gives alphabetical order
    edges_array.append(filepath[:-4]) # don't want the .cif part

os.chdir('../nodes_database')
nodes_array = []
for filepath in sorted(glob.iglob('*.cif')):
    nodes_array.append(filepath[:-4])

os.chdir('../template_database')
nets_array = []
for filepath in sorted(glob.iglob('*.cif')):
    nets_array.append(filepath[:-4])

print('edges')
print(edges_array)

print('\n')
print('nodes')
print(nodes_array)

print('\n')
print('nets')
print(nets_array)
