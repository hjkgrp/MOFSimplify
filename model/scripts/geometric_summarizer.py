import sys
import os
import subprocess
import pandas as pd
import numpy as np


''' The geometric descriptors are largest included sphere (Di), 
largest free sphere (Df), largest included sphere along free path (Dif),
crystal density (rho), volumetric surface area (VSA), gravimetric surface (GSA), 
volumetric pore volume (VPOV) and gravimetric pore volume (GPOV). 
Also, we include cell volume as a descriptor.

All Zeo++ calculations use a pore radius of 1.86 angstrom, and zeo++ is called by subprocess.
'''

dict_list = []
base_dir = sys.argv[1] #base_dir must be an absolute path
if base_dir[-1] != '/':
    base_dir+='/'
for cif_file in os.listdir(base_dir+'/primitive/'):
    print('---- now on ----, '+cif_file)
    if '.cif' not in cif_file:
        continue
    basename = cif_file.strip('.cif')
    largest_included_sphere, largest_free_sphere, largest_included_sphere_along_free_sphere_path  = np.nan, np.nan, np.nan
    unit_cell_volume, crystal_density, VSA, GSA  = np.nan, np.nan, np.nan, np.nan
    VPOV, GPOV = np.nan, np.nan
    POAV, PONAV, GPOAV, GPONAV, POAV_volume_fraction, PONAV_volume_fraction = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    if (os.path.exists(base_dir+'geometric/'+str(basename)+'_pd.txt') & os.path.exists(base_dir+'geometric/'+str(basename)+'_sa.txt') &
        os.path.exists(base_dir+'geometric/'+str(basename)+'_av.txt') & os.path.exists(base_dir+'geometric/'+str(basename)+'_pov.txt')
        ):
        with open(base_dir+'geometric/'+str(basename)+'_pd.txt') as f:
            pore_diameter_data = f.readlines()
            for row in pore_diameter_data:
                largest_included_sphere = float(row.split()[1]) # largest included sphere
                largest_free_sphere = float(row.split()[2]) # largest free sphere
                largest_included_sphere_along_free_sphere_path = float(row.split()[3]) # largest included sphere along free sphere path
        with open(base_dir+'geometric/'+str(basename)+'_sa.txt') as f:
            surface_area_data = f.readlines()
            for i, row in enumerate(surface_area_data):
                if i == 0:
                    unit_cell_volume = float(row.split('Unitcell_volume:')[1].split()[0]) # unit cell volume
                    crystal_density = float(row.split('Unitcell_volume:')[1].split()[0]) # crystal density
                    VSA = float(row.split('ASA_m^2/cm^3:')[1].split()[0]) # volumetric surface area
                    GSA = float(row.split('ASA_m^2/g:')[1].split()[0]) # gravimetric surface area
        with open(base_dir+'geometric/'+str(basename)+'_pov.txt') as f:
            pore_volume_data = f.readlines()
            for i, row in enumerate(pore_volume_data):
                if i == 0:
                    density = float(row.split('Density:')[1].split()[0])
                    POAV = float(row.split('POAV_A^3:')[1].split()[0]) # Probe accessible pore volume
                    PONAV = float(row.split('PONAV_A^3:')[1].split()[0]) # Probe non-accessible probe volume
                    GPOAV = float(row.split('POAV_cm^3/g:')[1].split()[0])
                    GPONAV = float(row.split('PONAV_cm^3/g:')[1].split()[0])
                    POAV_volume_fraction = float(row.split('POAV_Volume_fraction:')[1].split()[0]) # probe accessible volume fraction
                    PONAV_volume_fraction = float(row.split('PONAV_Volume_fraction:')[1].split()[0]) # probe non accessible volume fraction
                    VPOV = POAV_volume_fraction+PONAV_volume_fraction
                    GPOV = VPOV/density
    else:
        print(basename, 'not all 4 files exist!', 'sa: ',os.path.exists(base_dir+'geometric/'+str(basename)+'_sa.txt'), 
              'pd: ',os.path.exists(base_dir+'geometric/'+str(basename)+'_pd.txt'), 'av: ', os.path.exists(base_dir+'geometric/'+str(basename)+'_av.txt'),
              'pov: ', os.path.exists(base_dir+'geometric/'+str(basename)+'_pov.txt'))
    geo_dict = {'name':basename, 'cif_file':cif_file, 'Di':largest_included_sphere, 'Df': largest_free_sphere, 'Dif': largest_included_sphere_along_free_sphere_path,
                'rho': crystal_density, 'VSA':VSA, 'GSA': GSA, 'VPOV': VPOV, 'GPOV':GPOV, 'POAV_vol_frac':POAV_volume_fraction, 
                'PONAV_vol_frac':PONAV_volume_fraction, 'GPOAV':GPOAV,'GPONAV':GPONAV,'POAV':POAV,'PONAV':PONAV}
    dict_list.append(geo_dict)
geo_df = pd.DataFrame(dict_list)
geo_df.to_csv(base_dir+'/geometric_parameters.csv',index=False)



