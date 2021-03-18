from molSimplify.Informatics.MOF.MOF_descriptors import *
get_primitive('temp_cif.cif', 'temp_cif_primitive.cif')
full_names, full_descriptors = get_MOF_descriptors('temp_cif_primitive.cif',3,path='/', xyzpath='temp_cif.xyz')
