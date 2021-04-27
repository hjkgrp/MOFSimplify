#!/usr/bin/env python
# -*- coding: utf-8 -*-
import flask
import pandas as pd
import tensorflow.keras
import numpy as np
import time
from tensorflow.keras import backend as K
import json
import string
import os
import shutil
import subprocess
import molSimplify
import pickle
import molSimplify.Classes.mol3D as ms_mol3D
import molSimplify.Informatics.RACassemble as ms_RAC
import molSimplify.python_nn.tf_ANN as ms_ANN
import pathlib 
import sys
from molSimplify.Scripts.generator import startgen_pythonic
from molSimplify.Scripts.molSimplify_io import getlicores
from bokeh.plotting import figure
from bokeh.resources import CDN
from bokeh.models import ColumnDataSource, SingleIntervalTicker, LinearAxis
from bokeh.embed import file_html
from bokeh.models import Span
from bokeh.models import ColorBar, LinearColorMapper, LogColorMapper, HoverTool
from bokeh.models.markers import Circle
from bokeh.palettes import Inferno256#,RdBu11#,Viridis256
from flask import jsonify
from molSimplify.Informatics.MOF.MOF_descriptors import get_primitive, get_MOF_descriptors;

cmap_bokeh = Inferno256
# cmap_bokeh = RdBu11 # optional select other colormaps

app = flask.Flask(__name__)

MOFSIMPLIFY_PATH = os.path.abspath('.')

@app.route('/demo')
def serve_demo():
    # Go to localhost:8000/demo to see this.
    return 'you are on the demo page'

@app.route('/logo.png')
def serve_logo():
    return flask.send_from_directory('.', 'logo.png')

@app.route('/truncated_logo.png')
def serve_truncated_logo():
    return flask.send_from_directory('.', 'truncated_logo.png')

@app.route('/')
def serve_homepage():
    # Serves homepage
    return flask.send_from_directory('.', 'index.html')

@app.route('/about.html')
def serve_about():
    # Serves homepage
    return flask.send_from_directory('.', 'about.html')

@app.route('/mof_examples/<path:path>')
def serve_example(path):
    # Serves example
    return flask.send_from_directory('mof_examples', path)

@app.route('/temp_file_creation/tobacco_3.0/output_cifs/<path:path>')
def serve_bb(path):
    # Serves bb-generated MOF
    return flask.send_from_directory('temp_file_creation/tobacco_3.0/output_cifs', path)

@app.route('/ris_files/<path:path>')
def serve_ris(path):
    # Serves homepage
    return flask.send_from_directory('ris_files', path)

@app.route('/how_to_cite.html')
def serve_cite():
    # Serves homepage
    return flask.send_from_directory('.', 'how_to_cite.html')

@app.route('/libraries/<path:path>')
def serve_library_files(path):
    # Serves libraries
    return flask.send_from_directory('libraries', path)

@app.route('/bbcif/<path:path>')
def serve_bbcif(path):
    # Serves the building block generated MOF
    return flask.send_from_directory('temp_file_creation/tobacco_3.0/output_cifs', path);

def listdir_nohidden(path): # used for bb_generate. Ignores hidden files
    myList = os.listdir(path);
    for i in myList:
        if i.startswith('.'):
            myList.remove(i)
    return myList

@app.route('/get_bb_generated_MOF', methods=['POST']) # Gianmarco Terrones addition
def bb_generate():
    # generates a MOF using the building blocks and net specified by the user
    # uses ToBaCCo code, version 3.0
    # returns the constructed MOF's name to the front end

    # To begin, always go to main directory
    os.chdir(MOFSIMPLIFY_PATH)

    # Grab data
    mydata = json.loads(flask.request.get_data())

    linker = mydata['linker']
    sbu = mydata['sbu']
    net = mydata['net']

    os.chdir("temp_file_creation/tobacco_3.0")

    # clear the edges, nodes, templates, and output cifs folders to start fresh
        # when running python tobacco.py, it looks in these folders

    shutil.rmtree('edges')
    os.mkdir('edges')
    shutil.rmtree('nodes')
    os.mkdir('nodes')
    shutil.rmtree('templates')
    os.mkdir('templates')
    shutil.rmtree('output_cifs')
    os.mkdir('output_cifs')

    # copy over the linker, sbu, and net specified by the user in the edges, nodes, and templates folders

    shutil.copy('edges_database/' + linker + '.cif', 'edges/' + linker + '.cif')
    shutil.copy('nodes_database/' + sbu + '.cif', 'nodes/' + sbu + '.cif')
    shutil.copy('template_database/' + net + '.cif', 'templates/' + net + '.cif')


    # run the command to construct the MOF
    subprocess.run(['python', 'tobacco.py'])

    # if successful, there will be an output cif in the folder output_cifs
    if listdir_nohidden('output_cifs') == []: # no files in folder
        print('Construction failed.')
        return 'FAILED'
    
    constructed_MOF = listdir_nohidden('output_cifs')
    constructed_MOF = constructed_MOF[0] # getting the first, and only, element out of the list

    dictionary = {};

    dictionary['mof_name'] = constructed_MOF

    # getting the primitive cell using molSimplify
    get_primitive('output_cifs/' + constructed_MOF, 'output_cifs/primitive_' + constructed_MOF);

    # mof_file = open('output_cifs/primitive_' + constructed_MOF, 'r'); # reading the primitive file for information about the MOF
    # mof_info = mof_file.read();
    # mof_file.close();

    # dictionary['mof_info'] = mof_info;

    json_object = json.dumps(dictionary, indent = 4); # TODO dict not necessary now

    return json_object

    
# Note: the h5 model for the solvent stability prediction and the thermal stability prediction should be trained on the same version of Terachem (here, 1.14)
# the two h5 models show up in solvent_ANN.py and thermal_ANN.py, respectively
@app.route('/predict_solvent_stability', methods=['POST']) # Gianmarco Terrones addition
def ss_predict():
    # Generates solvent stability prediction
    # To do this, need to generate RAC featurization and Zeo++ geometry information for the MOF
    # Then, apply Aditya's model to make prediction

    print('TIME CHECK 1')

    # To begin, always go to main directory
    os.chdir(MOFSIMPLIFY_PATH)

    # Grab data
    mydata = json.loads(flask.request.get_data())

    os.chdir("temp_file_creation") # changing directory

    # Write the data back to a cif file
    cif_file = open('temp_cif.cif', 'w')
    cif_file.write(mydata)
    cif_file.close()


    # delete the RACs folder, then remake it (to start fresh for this prediction)
    shutil.rmtree('RACs')
    os.mkdir('RACs')

    # doing the same with the Zeo++ folder
    shutil.rmtree('zeo++')
    os.mkdir('zeo++')

    os.chdir("RACs") # move to RACs folder

    print('TIME CHECK 2')

    # Next, running MOF featurization
    try:
        get_primitive('../temp_cif.cif', 'temp_cif_primitive.cif');
    except ValueError:
        return 'FAILED'

    try:
        full_names, full_descriptors = get_MOF_descriptors('temp_cif_primitive.cif',3,path= str(pathlib.Path().absolute()), xyzpath= 'temp_cif.xyz');
    except ValueError:
        return 'FAILED'
    except NotImplementedError:
        return 'FAILED'
    except AssertionError:
        return 'FAILED'

    if (len(full_names) <= 1) and (len(full_descriptors) <= 1): # this is a featurization check from MOF_descriptors.py
        return 'FAILED'

    print('TIME CHECK 3')

    # At this point, have the RAC featurization. Need geometry information next.

    # Run Zeo++
    os.chdir('..')
    os.chdir('zeo++')

    import time # debugging
    timeStarted = time.time() # save start time (debugging)

    cmd1 = '../../zeo++-0.3/network -ha -res temp_cif_pd.txt ' + '../RACs/temp_cif_primitive.cif'
    cmd2 = '../../zeo++-0.3/network -sa 1.86 1.86 10000 temp_cif_sa.txt ' + '../RACs/temp_cif_primitive.cif'
    cmd3 = '../../zeo++-0.3/network -ha -vol 1.86 1.86 10000 temp_cif_av.txt ' + '../RACs/temp_cif_primitive.cif'
    cmd4 = '../../zeo++-0.3/network -volpo 1.86 1.86 10000 temp_cif_pov.txt '+ '../RACs/temp_cif_primitive.cif'
    # four parallelized Zeo++ commands
    process1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=None, shell=True)
    process2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, stderr=None, shell=True)
    process3 = subprocess.Popen(cmd3, stdout=subprocess.PIPE, stderr=None, shell=True)
    process4 = subprocess.Popen(cmd4, stdout=subprocess.PIPE, stderr=None, shell=True)
    output1 = process1.communicate()[0]
    output2 = process2.communicate()[0]
    output3 = process3.communicate()[0]
    output4 = process4.communicate()[0]

    timeDelta = time.time() - timeStarted # get execution time
    print('Finished process in ' + str(timeDelta) + ' seconds')

    # Have written output of Zeo++ commands to files. Now, code below extracts information from those files

    ''' The geometric descriptors are largest included sphere (Di), 
    largest free sphere (Df), largest included sphere along free path (Dif),
    crystal density (rho), volumetric surface area (VSA), gravimetric surface (GSA), 
    volumetric pore volume (VPOV) and gravimetric pore volume (GPOV). 
    Also, we include cell volume as a descriptor.

    All Zeo++ calculations use a pore radius of 1.86 angstrom, and zeo++ is called by subprocess.
    '''

    dict_list = []
    # base_dir = sys.argv[1] #base_dir must be an absolute path
    # if base_dir[-1] != '/':
    #     base_dir+='/'
    # for cif_file in os.listdir(base_dir+'/primitive/'):
    #     print('---- now on ----, '+cif_file)
    #     if '.cif' not in cif_file:
    #         continue
    #     basename = cif_file.strip('.cif')
    cif_file = 'temp.cif' # techincally, calculations were with the primitive, but I'll just call it temp
    basename = cif_file.strip('.cif')
    largest_included_sphere, largest_free_sphere, largest_included_sphere_along_free_sphere_path  = np.nan, np.nan, np.nan
    unit_cell_volume, crystal_density, VSA, GSA  = np.nan, np.nan, np.nan, np.nan
    VPOV, GPOV = np.nan, np.nan
    POAV, PONAV, GPOAV, GPONAV, POAV_volume_fraction, PONAV_volume_fraction = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    if (os.path.exists('temp_cif_pd.txt') & os.path.exists('temp_cif_sa.txt') &
        os.path.exists('temp_cif_av.txt') & os.path.exists('temp_cif_pov.txt')
        ):
        with open('temp_cif_pd.txt') as f:
            pore_diameter_data = f.readlines()
            for row in pore_diameter_data:
                largest_included_sphere = float(row.split()[1]) # largest included sphere
                largest_free_sphere = float(row.split()[2]) # largest free sphere
                largest_included_sphere_along_free_sphere_path = float(row.split()[3]) # largest included sphere along free sphere path
        with open('temp_cif_sa.txt') as f:
            surface_area_data = f.readlines()
            for i, row in enumerate(surface_area_data):
                if i == 0:
                    unit_cell_volume = float(row.split('Unitcell_volume:')[1].split()[0]) # unit cell volume
                    crystal_density = float(row.split('Unitcell_volume:')[1].split()[0]) # crystal density
                    VSA = float(row.split('ASA_m^2/cm^3:')[1].split()[0]) # volumetric surface area
                    GSA = float(row.split('ASA_m^2/g:')[1].split()[0]) # gravimetric surface area
        with open('temp_cif_pov.txt') as f:
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
        print('Not all 4 files exist, so at least one Zeo++ call failed!', 'sa: ',os.path.exists('temp_cif_sa.txt'), 
              '; pd: ',os.path.exists('temp_cif_pd.txt'), '; av: ', os.path.exists('temp_cif_av.txt'),
              '; pov: ', os.path.exists('temp_cif_pov.txt'))
        return 'FAILED'
    geo_dict = {'name':basename, 'cif_file':cif_file, 'Di':largest_included_sphere, 'Df': largest_free_sphere, 'Dif': largest_included_sphere_along_free_sphere_path,
                'rho': crystal_density, 'VSA':VSA, 'GSA': GSA, 'VPOV': VPOV, 'GPOV':GPOV, 'POAV_vol_frac':POAV_volume_fraction, 
                'PONAV_vol_frac':PONAV_volume_fraction, 'GPOAV':GPOAV,'GPONAV':GPONAV,'POAV':POAV,'PONAV':PONAV}
    dict_list.append(geo_dict)
    geo_df = pd.DataFrame(dict_list)
    geo_df.to_csv('geometric_parameters.csv',index=False)


    print('TIME CHECK 4')

    # Applying the model next

    # Merging geometric information with get_MOF_descriptors files (lc_descriptors.csv, sbu_descriptors.csv, linker_descriptors.csv)

    # from IPython.display import display # debugging

    lc_df = pd.read_csv("../../temp_file_creation/RACs/lc_descriptors.csv") 
    sbu_df = pd.read_csv("../../temp_file_creation/RACs/sbu_descriptors.csv")
    linker_df = pd.read_csv("../../temp_file_creation/RACs/linker_descriptors.csv")

    lc_df = lc_df.mean().to_frame().transpose() # averaging over all rows. Convert resulting Series into a Dataframe, then transpose
    sbu_df = sbu_df.mean().to_frame().transpose()
    linker_df = linker_df.mean().to_frame().transpose()

    # print('check U')
    # display(geo_df)
    # print(geo_df.columns)

    # print('check V')
    # display(lc_df) # debugging
    # print(lc_df.columns)

    merged_df = pd.concat([geo_df, lc_df, sbu_df, linker_df], axis=1)

    # print('check W')
    # display(merged_df) # debugging

    # merged_df = merged_df.merge(sbu_df, how='outer') 
    # merged_df = merged_df.merge(linker_df, how='outer') 

    # print('check X')
    # display(merged_df) # debugging

    merged_df.to_csv('../merged_descriptors.csv',index=False) # written in /temp_file_creation

    os.chdir('..')
    os.chdir('..')
    os.chdir('model/solvent/ANN')

    print('TIME CHECK 5')

    os.system('python solvent_ANN.py > solvent_prediction.txt')
    # import for_GT
    # prediction = for_GT.main()

    f = open("solvent_prediction.txt", "r")
    line = f.read()
    line = line.split('[')
    line = line[2]
    line = line.split(']')
    prediction = line[0] # isolating just the prediction, since the model spits out the prediction like [[PREDICTION]], as in, in hard brackets
    f.close()

    print('TIME CHECK 6')

    return prediction

@app.route('/predict_thermal_stability', methods=['POST']) # Gianmarco Terrones addition
def ts_predict():
    # Generates thermal stability prediction 
    # To do this, need to generate RAC featurization and Zeo++ geometry information for the MOF
    # Then, apply Aditya's model to make prediction

    print('TIME CHECK 1')

    # To begin, always go to main directory 
    os.chdir(MOFSIMPLIFY_PATH)

    # Grab data
    mydata = json.loads(flask.request.get_data())

    os.chdir("temp_file_creation") # changing directory

    # Write the data back to a cif file
    cif_file = open('temp_cif.cif', 'w')
    cif_file.write(mydata)
    cif_file.close()

    # delete the RACs folder, then remake it (to start fresh for this prediction)
    shutil.rmtree('RACs')
    os.mkdir('RACs')

    os.chdir("RACs") # move to RACs folder

    print('TIME CHECK 2')

    # Next, running MOF featurization
    try:
        get_primitive('../temp_cif.cif', 'temp_cif_primitive.cif');
    except ValueError:
        return 'FAILED'

    try:
        full_names, full_descriptors = get_MOF_descriptors('temp_cif_primitive.cif',3,path= str(pathlib.Path().absolute()), xyzpath= 'temp_cif.xyz');
    except ValueError:
        return 'FAILED'
    except NotImplementedError:
        return 'FAILED'
    except AssertionError:
        return 'FAILED'

    if (len(full_names) <= 1) and (len(full_descriptors) <= 1): # this is a featurization check from MOF_descriptors.py
        return 'FAILED'

    print('TIME CHECK 3')

    # At this point, have the RAC featurization. Need geometry information next.

    # Run Zeo++
    os.chdir('..')
    os.chdir('zeo++')

    import time # debugging
    timeStarted = time.time() # save start time (debugging)

    cmd1 = '../../zeo++-0.3/network -ha -res temp_cif_pd.txt ' + '../RACs/temp_cif_primitive.cif'
    cmd2 = '../../zeo++-0.3/network -sa 1.86 1.86 10000 temp_cif_sa.txt ' + '../RACs/temp_cif_primitive.cif'
    cmd3 = '../../zeo++-0.3/network -ha -vol 1.86 1.86 10000 temp_cif_av.txt ' + '../RACs/temp_cif_primitive.cif'
    cmd4 = '../../zeo++-0.3/network -volpo 1.86 1.86 10000 temp_cif_pov.txt '+ '../RACs/temp_cif_primitive.cif'
    # four parallelized Zeo++ commands
    process1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=None, shell=True)
    process2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, stderr=None, shell=True)
    process3 = subprocess.Popen(cmd3, stdout=subprocess.PIPE, stderr=None, shell=True)
    process4 = subprocess.Popen(cmd4, stdout=subprocess.PIPE, stderr=None, shell=True)
    output1 = process1.communicate()[0]
    output2 = process2.communicate()[0]
    output3 = process3.communicate()[0]
    output4 = process4.communicate()[0]

    timeDelta = time.time() - timeStarted # get execution time
    print('Finished process in ' + str(timeDelta) + ' seconds')

    # Have written output of Zeo++ commands to files. Now, code below extracts information from those files

    ''' The geometric descriptors are largest included sphere (Di), 
    largest free sphere (Df), largest included sphere along free path (Dif),
    crystal density (rho), volumetric surface area (VSA), gravimetric surface (GSA), 
    volumetric pore volume (VPOV) and gravimetric pore volume (GPOV). 
    Also, we include cell volume as a descriptor.

    All Zeo++ calculations use a pore radius of 1.86 angstrom, and zeo++ is called by subprocess.
    '''

    dict_list = []
    # base_dir = sys.argv[1] #base_dir must be an absolute path
    # if base_dir[-1] != '/':
    #     base_dir+='/'
    # for cif_file in os.listdir(base_dir+'/primitive/'):
    #     print('---- now on ----, '+cif_file)
    #     if '.cif' not in cif_file:
    #         continue
    #     basename = cif_file.strip('.cif')
    cif_file = 'temp.cif' # techincally, calculations were with the primitive, but I'll just call it temp
    basename = cif_file.strip('.cif')
    largest_included_sphere, largest_free_sphere, largest_included_sphere_along_free_sphere_path  = np.nan, np.nan, np.nan
    unit_cell_volume, crystal_density, VSA, GSA  = np.nan, np.nan, np.nan, np.nan
    VPOV, GPOV = np.nan, np.nan
    POAV, PONAV, GPOAV, GPONAV, POAV_volume_fraction, PONAV_volume_fraction = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    if (os.path.exists('temp_cif_pd.txt') & os.path.exists('temp_cif_sa.txt') &
        os.path.exists('temp_cif_av.txt') & os.path.exists('temp_cif_pov.txt')
        ):
        with open('temp_cif_pd.txt') as f:
            pore_diameter_data = f.readlines()
            for row in pore_diameter_data:
                largest_included_sphere = float(row.split()[1]) # largest included sphere
                largest_free_sphere = float(row.split()[2]) # largest free sphere
                largest_included_sphere_along_free_sphere_path = float(row.split()[3]) # largest included sphere along free sphere path
        with open('temp_cif_sa.txt') as f:
            surface_area_data = f.readlines()
            for i, row in enumerate(surface_area_data):
                if i == 0:
                    unit_cell_volume = float(row.split('Unitcell_volume:')[1].split()[0]) # unit cell volume
                    crystal_density = float(row.split('Unitcell_volume:')[1].split()[0]) # crystal density
                    VSA = float(row.split('ASA_m^2/cm^3:')[1].split()[0]) # volumetric surface area
                    GSA = float(row.split('ASA_m^2/g:')[1].split()[0]) # gravimetric surface area
        with open('temp_cif_pov.txt') as f:
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
        print('Not all 4 files exist, so at least one Zeo++ call failed!', 'sa: ',os.path.exists('temp_cif_sa.txt'), 
              '; pd: ',os.path.exists('temp_cif_pd.txt'), '; av: ', os.path.exists('temp_cif_av.txt'),
              '; pov: ', os.path.exists('temp_cif_pov.txt'))
        return 'FAILED'
    geo_dict = {'name':basename, 'cif_file':cif_file, 'Di':largest_included_sphere, 'Df': largest_free_sphere, 'Dif': largest_included_sphere_along_free_sphere_path,
                'rho': crystal_density, 'VSA':VSA, 'GSA': GSA, 'VPOV': VPOV, 'GPOV':GPOV, 'POAV_vol_frac':POAV_volume_fraction, 
                'PONAV_vol_frac':PONAV_volume_fraction, 'GPOAV':GPOAV,'GPONAV':GPONAV,'POAV':POAV,'PONAV':PONAV}
    dict_list.append(geo_dict)
    geo_df = pd.DataFrame(dict_list)
    geo_df.to_csv('geometric_parameters.csv',index=False)


    print('TIME CHECK 4')

    # Applying the model next

    # Merging geometric information with get_MOF_descriptors files (lc_descriptors.csv, sbu_descriptors.csv, linker_descriptors.csv)

    lc_df = pd.read_csv("../../temp_file_creation/RACs/lc_descriptors.csv") 
    sbu_df = pd.read_csv("../../temp_file_creation/RACs/sbu_descriptors.csv")
    linker_df = pd.read_csv("../../temp_file_creation/RACs/linker_descriptors.csv")

    lc_df = lc_df.mean().to_frame().transpose() # averaging over all rows. Convert resulting Series into a Dataframe, then transpose
    sbu_df = sbu_df.mean().to_frame().transpose()
    linker_df = linker_df.mean().to_frame().transpose()

    merged_df = pd.concat([geo_df, lc_df, sbu_df, linker_df], axis=1)
    merged_df.to_csv('../merged_descriptors.csv',index=False) # written in /temp_file_creation

    os.chdir('..')
    os.chdir('..')
    os.chdir('model/thermal/ANN')

    print('TIME CHECK 5')

    os.system('python thermal_ANN.py > thermal_prediction.txt')

    f = open("thermal_prediction.txt", "r")
    prediction = f.read()
    f.close()

    print('TIME CHECK 6')

    return prediction

@app.route('/plot_thermal_stability', methods=['POST']) # Gianmarco Terrones addition
def plot_thermal_stability():
    # returns a plot of the distribution of thermal breakdown temperatures of the MOFs our ANN was trained on
    # additionally, displays the position of the current MOF's thermal breakdown temperature

    # To begin, always go to main directory 
    os.chdir(MOFSIMPLIFY_PATH)

    # debugging
    print('Check A')
    print(os.getcwd())

    # Grab data
    mydata = json.loads(flask.request.get_data()) # this is the current MOF's predicted thermal breakdown temperature
    mydata = mydata[:-3] # getting rid of the celsius symbol, left with just the number
    mydata = float(mydata)
    print(mydata)

    # Getting the temperature data
    temps_df = pd.read_csv("model/thermal/ANN/adjusted_TSD_df_all.csv")

    print('Check A2')

    import matplotlib
    matplotlib.use('Agg') # noninteractive backend
    import matplotlib.pyplot as plt
    plt.close("all")
    import scipy.stats as stats

    # in training data, smallest T breakdown is 35, and largest T breakdown is 654

    # use stats.gaussian_kde to estimate the probability density function from the histogram
    density = stats.gaussian_kde(temps_df['T'])
    x = np.arange(30,661,1) # in training data, smallest T breakdown is 35, and largest T breakdown is 654
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.plot(x, density(x))
    print('Check A3')
    print(mydata)
    print(density(mydata))
    print('sanity check')
    print(density(0))
    print(density(700))
    plt.plot(mydata, density(mydata), "or") # the current MOF's predicted thermal breakdown temperature

    print('Check B')

    print('Check C2')

    # labelling axes
    # myHist.set_xlabel('Breakdown temperature (°C)')
    # myHist.set_ylabel('Number of MOFs in training data')
    # myHist.set_title('Current MOF\'s breakdown temperature relative to others')
    ax.set_xlabel('Breakdown temperature (°C)')
    ax.set_ylabel('Frequency in the training data')
    ax.set_title('Current MOF\'s breakdown temperature relative to others')

    import mpld3

    return mpld3.fig_to_html(fig)


@app.route('/get_components', methods=['POST']) # Gianmarco Terrones addition
def get_components():
    # Uses Aditya's MOF code to get linkers and sbus
    # Returns a dictionary with the linker and sbu xyz files's text, along with information about the number of linkers and sbus

    # To begin, always go to main directory 
    os.chdir(MOFSIMPLIFY_PATH);

    # Grab data
    mydata = json.loads(flask.request.get_data());

    os.chdir("temp_file_creation"); # changing directory

    # Write the data back to a cif file
    cif_file = open('temp_cif.cif', 'w');
    cif_file.write(mydata);
    cif_file.close();


    # delete the RACs folder, then remake it (to start fresh for this prediction)
    shutil.rmtree('RACs');
    os.mkdir('RACs');

    os.chdir("RACs"); # move to RACs folder

    # Next, running MOF featurization
    try:
        get_primitive('../temp_cif.cif', 'temp_cif_primitive.cif');
    except ValueError:
        return 'FAILED'

    try:
        full_names, full_descriptors = get_MOF_descriptors('temp_cif_primitive.cif',3,path= str(pathlib.Path().absolute()), xyzpath= 'temp_cif.xyz');
    except ValueError:
        return 'FAILED'
    except NotImplementedError:
        return 'FAILED'
    except AssertionError:
        return 'FAILED'

    if (len(full_names) <= 1) and (len(full_descriptors) <= 1): # this is a featurization check from MOF_descriptors.py
        return 'FAILED'

    # At this point, have the RAC featurization. 

    # will return a json object
    # the fields are string representations of the linkers and sbus, however many there are

    dictionary = {};
    
    linker_num = 0;
    while True:
        if not os.path.exists('linkers/temp_cif_primitive_linker_' + str(linker_num) + '.xyz'):
            break
        else:
            linker_file = open('linkers/temp_cif_primitive_linker_' + str(linker_num) + '.xyz', 'r');
            linker_info = linker_file.read();
            linker_file.close();

            dictionary['linker_' + str(linker_num)] = linker_info;

            linker_num = linker_num + 1;


    sbu_num = 0;
    while True:
        if not os.path.exists('sbus/temp_cif_primitive_sbu_' + str(sbu_num) + '.xyz'):
            break
        else:
            sbu_file = open('sbus/temp_cif_primitive_sbu_' + str(sbu_num) + '.xyz', 'r');
            sbu_info = sbu_file.read();
            sbu_file.close();

            dictionary['sbu_' + str(sbu_num)] = sbu_info;

            sbu_num = sbu_num + 1;


    dictionary['total_linkers'] = linker_num;
    dictionary['total_sbus'] = sbu_num;


    os.chdir("linkers"); # move to linkers folder

    # Identifying which linkers and sbus have different connectivities

    # Code below uses molecular graph determinants 
    # two ligands that are the same (via connectivity) will have the same molecular graph determinant. 
    # Molecular graph determinant fails to distinguish between isomers, but so would RCM (Reverse Cuthill McKee)

    # This simple script is meant to be run within the linkers directory, and it will give a bunch of numbers. 
    # If those numbers are the same, the linker is the same, if not, the linkers are different, etc


    from molSimplify.Classes.mol3D import mol3D
    import glob
    MOF_of_interest = 'temp_cif_primitive'
    XYZs = sorted(glob.glob('*'+MOF_of_interest+'*xyz'))
    det_list = []

    for xyz in XYZs:
        net = xyz.replace('xyz', 'net') # substring replacement; getting the appropriate .net file
        linker_mol = mol3D()
        linker_mol.readfromxyz(xyz)
        with open(net,'r') as f:
            data = f.readlines()
        for i, line in enumerate(data):
            if i==0:
                # Skip first line
                continue
            elif i == 1:
                graph = np.array([[int(val) for val in line.strip('\n').split(',')]])
            else:
                graph = np.append(graph, [np.array([int(val) for val in line.strip('\n').split(',')])],axis=0)
        linker_mol.graph = graph
        safedet = linker_mol.get_mol_graph_det(oct=False)
        det_list.append(safedet)
    #### linkers with the same molecular graph determinant are the same
    #### molecular graph determinant does not catch isomers
    linker_det_list = det_list

    unique_linker_det = set(linker_det_list) # getting the unique determinants
    unique_linker_indices = []
    for item in unique_linker_det:
        unique_linker_indices.append(linker_det_list.index(item)) # indices of the unique linkers in the list of linkers



    os.chdir("../sbus"); # move to sbus folder

    XYZs = sorted(glob.glob('*'+MOF_of_interest+'*xyz'))
    det_list = []
    for xyz in XYZs:
        net = xyz.replace('xyz', 'net') # substring replacement; getting the appropriate .net file
        linker_mol = mol3D()
        linker_mol.readfromxyz(xyz)
        with open(net,'r') as f:
            data = f.readlines()
        for i, line in enumerate(data):
            if i==0:
                # Skip first line
                continue
            elif i == 1:
                graph = np.array([[int(val) for val in line.strip('\n').split(',')]])
            else:
                graph = np.append(graph, [np.array([int(val) for val in line.strip('\n').split(',')])],axis=0)
        linker_mol.graph = graph
        safedet = linker_mol.get_mol_graph_det(oct=False)
        det_list.append(safedet)

    sbu_det_list = det_list

    unique_sbu_det = set(sbu_det_list) # getting the unique determinants
    unique_sbu_indices = []
    for item in unique_sbu_det:
        unique_sbu_indices.append(sbu_det_list.index(item)) # indices of the unique sbus in the list of linkers

    print('det lists')
    print(linker_det_list)
    print(sbu_det_list)


    # adding the unique indices to the dictionary
    dictionary['unique_linker_indices'] = unique_linker_indices;
    dictionary['unique_sbu_indices'] = unique_sbu_indices;

    json_object = json.dumps(dictionary, indent = 4);

    return json_object


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8000)
