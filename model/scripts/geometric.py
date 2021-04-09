import sys
import os
import subprocess


''' The geometric descriptors are largest included sphere (Di), 
largest free sphere (Df), largest included sphere along free path (Dif),
crystal density (rho), volumetric surface area (VSA), gravimetric surface (GSA), 
volumetric pore volume (VPOV) and gravimetric pore volume (GPOV). 
Also, we include cell volume as a descriptor.

All Zeo++ calculations use a pore radius of 1.86 angstrom, and zeo++ is called by subprocess.
'''
base_dir = sys.argv[1] #base_dir must be an absolute path
if base_dir[-1] != '/':
    base_dir+='/'
if not os.path.exists(base_dir+'/primitive/'):
    os.mkdir(base_dir+'/geometric/')
for cif_file in os.listdir(base_dir+'/primitive/'):
    print('---- now on ----, '+cif_file)
    if '.cif' not in cif_file:
        continue
    basename = cif_file.strip('.cif')
    cmd1 = '/Users/adityanandy/zeo++-0.3/network -ha -res '+base_dir+'geometric/'+str(basename)+'_pd.txt '+base_dir+'/primitive/'+cif_file
    # print(cmd1)
    # sard
    cmd2 = '/Users/adityanandy/zeo++-0.3/network -sa 1.86 1.86 10000 '+base_dir+'geometric/'+str(basename)+'_sa.txt '+base_dir+'/primitive/'+cif_file
    cmd3 = '/Users/adityanandy/zeo++-0.3/network -ha -vol 1.86 1.86 10000 '+base_dir+'geometric/'+str(basename)+'_av.txt '+base_dir+'/primitive/'+cif_file
    cmd4 = '/Users/adityanandy/zeo++-0.3/network -volpo 1.86 1.86 10000 '+base_dir+'geometric/'+str(basename)+'_pov.txt '+base_dir+'/primitive/'+cif_file
    process1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=None, shell=True)
    process2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, stderr=None, shell=True)
    process3 = subprocess.Popen(cmd3, stdout=subprocess.PIPE, stderr=None, shell=True)
    process4 = subprocess.Popen(cmd4, stdout=subprocess.PIPE, stderr=None, shell=True)
    output1 = process1.communicate()[0]
    output2 = process2.communicate()[0]
    output3 = process3.communicate()[0]
    output4 = process4.communicate()[0]
