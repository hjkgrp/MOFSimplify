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

@app.route('/neighbor/<path:path>')
def serve_neighbor(path):
    # Serves the neighbor CoRE MOF
    return flask.send_from_directory('CoRE2019', path);

@app.route('/neighbor_info/<path:path>')
def serve_neighbor_txt(path):
    # Serves the neighbor CoRE MOF information in txt file format
    return flask.send_from_directory('temp_file_creation/latent_neighbor', path);

def listdir_nohidden(path): # used for bb_generate. Ignores hidden files
    myList = os.listdir(path);
    for i in myList:
        if i.startswith('.'):
            myList.remove(i)
    return myList

@app.route('/get_bb_generated_MOF', methods=['POST']) 
def bb_generate():
    # generates a MOF using the building blocks and net specified by the user
    # uses ToBaCCo code, version 3.0
    # returns the constructed MOF's name to the front end

    # To begin, always go to main directory
    os.chdir(MOFSIMPLIFY_PATH)

    # Grab data
    my_data = json.loads(flask.request.get_data())

    linker = my_data['linker']
    sbu = my_data['sbu']
    net = my_data['net']

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

    json_object = json.dumps(dictionary, indent = 4);

    return json_object

    
# Note: the h5 model for the solvent stability prediction and the thermal stability prediction should be trained on the same version of TensorFlow (here, 1.14)
# the two h5 models show up in solvent_ANN.py and thermal_ANN.py, respectively
@app.route('/predict_solvent_stability', methods=['POST']) 
def ss_predict():
    # Generates solvent stability prediction.
    # To do this, need to generate RAC featurization and Zeo++ geometry information for the MOF.
    # Then, apply Aditya's model to make prediction.

    print('TIME CHECK 1')

    # To begin, always go to main directory.
    os.chdir(MOFSIMPLIFY_PATH)

    # Grab data
    my_data = json.loads(flask.request.get_data())

    os.chdir("temp_file_creation") # changing directory

    # Write the data back to a cif file.
    cif_file = open('temp_cif.cif', 'w')
    cif_file.write(my_data)
    cif_file.close()


    # Delete the RACs folder, then remake it (to start fresh for this prediction).
    shutil.rmtree('RACs')
    os.mkdir('RACs')

    # Doing the same with the Zeo++ folder.
    shutil.rmtree('zeo++')
    os.mkdir('zeo++')

    os.chdir("RACs") # move to RACs folder

    print('TIME CHECK 2')
    import time # debugging
    timeStarted = time.time() # save start time (debugging)

    # Next, running MOF featurization
    try:
        get_primitive('../temp_cif.cif', 'temp_cif_primitive.cif');
    except ValueError:
        return 'FAILED'

    timeDelta = time.time() - timeStarted # get execution time
    print('Finished process in ' + str(timeDelta) + ' seconds')

    timeStarted = time.time() # save start time (debugging)

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

    timeDelta = time.time() - timeStarted # get execution time
    print('Finished process in ' + str(timeDelta) + ' seconds')

    print('TIME CHECK 3')

    # At this point, have the RAC featurization. Need geometry information next.

    # Run Zeo++
    os.chdir('..')
    os.chdir('zeo++')

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

    # Have written output of Zeo++ commands to files. Now, code below extracts information from those files.

    ''' The geometric descriptors are largest included sphere (Di), 
    largest free sphere (Df), largest included sphere along free path (Dif),
    crystal density (rho), volumetric surface area (VSA), gravimetric surface (GSA), 
    volumetric pore volume (VPOV) and gravimetric pore volume (GPOV). 
    Also, we include cell volume as a descriptor.

    All Zeo++ calculations use a pore radius of 1.86 angstrom, and zeo++ is called by subprocess.
    '''

    timeStarted = time.time() # save start time (debugging)

    dict_list = []
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


    timeDelta = time.time() - timeStarted # get execution time
    print('Finished process in ' + str(timeDelta) + ' seconds')

    print('TIME CHECK 4')

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
    os.chdir('model/solvent/ANN')


    # Here, I do a check to see if the current MOF is in the training data.
    # If it is, then I return the known truth for the MOF, rather than make a prediction.

    # Will iterate through the rows of the train pandas dataframe
    train_df = pd.read_csv('dropped_connectivity_dupes/train.csv')

    in_train = False

    # print('In train check')

    for index, row in train_df.iterrows(): # iterate through rows

        row_match = True # gets set to false if any values don't match 

        for col in merged_df.columns: # iterate through columns
            if col == 'name' or col == 'cif_file':
                continue # skip these
            # print('Column: ')
            # print(col)
            # print(row[col])
            # print(merged_df.iloc[0][col])
            if np.absolute(row[col] - merged_df.iloc[0][col]) > 0.05 * np.absolute(merged_df.iloc[0][col]): # row[col] != merged_df.iloc[0][col] was leading to some same values being idenfitied as different b/c of some floating 10^-15 values 
                row_match = False

                # debugging
                # print('TAGBUZ check')
                # if row['CoRE_name'] == 'TAGBUZ_clean':
                #     print(col)
                #     print(row[col])
                #     print(merged_df.iloc[0][col])

                break
        
        if row_match: # all columns for row match!
            in_train = True
            match_truth = row['flag'] # the flag for the MOF that matches the current MOF
            print(row['CoRE_name'])
            break

    if in_train:
        myDict = {'in_train': True, 'truth': match_truth}
        return myDict

    print('TIME CHECK 5')


    # Applying the model next

    timeStarted = time.time() # save start time (debugging)

    os.system('python solvent_ANN.py > solvent_prediction.txt')
    # import for_GT
    # prediction = for_GT.main()

    timeDelta = time.time() - timeStarted # get execution time
    print('Finished process in ' + str(timeDelta) + ' seconds')

    f = open("solvent_prediction.txt", "r")
    line = f.readline()
    line = line.split('[')
    line = line[2]
    line = line.split(']')
    prediction = line[0] # isolating just the prediction, since the model spits out the prediction like [[PREDICTION]], as in, in hard brackets
    f.readline() # skip a line
    neighbor_names = f.readline()
    neighbor_distances = f.readline()
    f.close()

    # Next, some hacky stuff to convert strings back into lists
    neighbor_names = neighbor_names.split('\', \'')
    neighbor_names[0] = neighbor_names[0][2:]
    neighbor_names[-1] = neighbor_names[-1][:-3]

    neighbor_distances = neighbor_distances.split(', ')
    neighbor_distances[0] = neighbor_distances[0][2:]
    neighbor_distances[-1] = neighbor_distances[-1][:-2]    

    print('check check')
    print(neighbor_names) # debugging
    print(neighbor_distances) # debugging
    print(type(neighbor_names)) # debugging
    print(type(neighbor_distances)) # debugging    

    results = {'prediction': prediction,
        'neighbor_names': neighbor_names,
        'neighbor_distances': neighbor_distances,
        'in_train': False} # a prediction was made. Requested MOF was not in the training data.

    print('TIME CHECK 6')

    return results

@app.route('/predict_thermal_stability', methods=['POST']) 
def ts_predict():
    # Generates thermal stability prediction.
    # To do this, need to generate RAC featurization and Zeo++ geometry information for the MOF.
    # Then, apply Aditya's model to make prediction.

    print('TIME CHECK 1')

    # To begin, always go to main directory.
    os.chdir(MOFSIMPLIFY_PATH)

    # Grab data
    my_data = json.loads(flask.request.get_data())

    os.chdir("temp_file_creation") # changing directory

    # Write the data back to a cif file
    cif_file = open('temp_cif.cif', 'w')
    cif_file.write(my_data)
    cif_file.close()

    # Delete the RACs folder, then remake it (to start fresh for this prediction).
    shutil.rmtree('RACs')
    os.mkdir('RACs')

    os.chdir("RACs") # move to RACs folder

    print('TIME CHECK 2')

    # Next, running MOF featurization.
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

    # Have written output of Zeo++ commands to files. Now, code below extracts information from those files.

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

    # Applying the model next.

    # Merging geometric information with get_MOF_descriptors files (lc_descriptors.csv, sbu_descriptors.csv, linker_descriptors.csv).

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



    # Here, I do a check to see if the current MOF is in the training data.
    # If it is, then I return the known truth for the MOF, rather than make a prediction.

    # Will iterate through the rows of the train pandas dataframe
    train_df = pd.read_csv('train.csv')

    in_train = False

    # print('In train check')

    for index, row in train_df.iterrows(): # iterate through rows

        row_match = True # gets set to false if any values don't match 

        for col in merged_df.columns: # iterate through columns
            if col == 'name' or col == 'cif_file':
                continue # skip these
            # print('Column: ')
            # print(col)
            # print(row[col])
            # print(merged_df.iloc[0][col])
            if np.absolute(row[col] - merged_df.iloc[0][col]) > 0.05 * np.absolute(merged_df.iloc[0][col]): # row[col] != merged_df.iloc[0][col] was leading to some identical values being identified as different b/c of some small differences 
                row_match = False

                # debugging
                print('BETDAH check')
                if row['CoRE_name'] == 'BETDAH_clean':
                    print(col)
                    print(row[col])
                    print(merged_df.iloc[0][col])

                break
        
        if row_match: # all columns for row match!
            in_train = True
            match_truth = row['T'] # the flag for the MOF that matches the current MOF

            match_truth = np.round(match_truth,1) # round to 1 decimal

            # adding units
            degree_sign= u'\N{DEGREE SIGN}'
            match_truth = str(match_truth) + degree_sign + 'C' # degrees Celsius

#            print(row['CoRE_name'])
            break

    if in_train:
        myDict = {'in_train': True, 'truth': match_truth}
        return myDict



    print('TIME CHECK 5')

    timeStarted = time.time() # save start time (debugging)

    os.system('python thermal_ANN.py > thermal_prediction.txt')

    timeDelta = time.time() - timeStarted # get execution time
    print('Finished process in ' + str(timeDelta) + ' seconds')

    f = open("thermal_prediction.txt", "r")
    prediction = f.readline()
    f.readline() # skip a line
    neighbor_names = f.readline()
    neighbor_distances = f.readline()
    f.close()

    # Next, some hacky stuff to convert strings back into lists.
    neighbor_names = neighbor_names.split('\', \'')
    neighbor_names[0] = neighbor_names[0][2:]
    neighbor_names[-1] = neighbor_names[-1][:-3]

    neighbor_distances = neighbor_distances.split(', ')
    neighbor_distances[0] = neighbor_distances[0][2:]
    neighbor_distances[-1] = neighbor_distances[-1][:-2]    

    print('check check')
    print(neighbor_names) # debugging
    print(neighbor_distances) # debugging
    print(type(neighbor_names)) # debugging
    print(type(neighbor_distances)) # debugging    

    results = {'prediction': prediction,
        'neighbor_names': neighbor_names,
        'neighbor_distances': neighbor_distances,
        'in_train': False} # a prediction was made. Requested MOF was not in the training data.

    print('TIME CHECK 6')

    return results

@app.route('/plot_thermal_stability', methods=['POST']) 
def plot_thermal_stability():
    # Returns a plot of the distribution of thermal breakdown temperatures of the MOFs our ANN was trained on.
    # Additionally, displays the position of the current MOF's thermal breakdown temperature.

    # To begin, always go to main directory 
    os.chdir(MOFSIMPLIFY_PATH)

    # Grab data
    info = json.loads(flask.request.get_data()) 
    my_data = info['temperature'] # this is the current MOF's predicted thermal breakdown temperature
    my_data = my_data[:-3] # getting rid of the celsius symbol, left with just the number
    my_data = float(my_data)
    print('checkerino')
    print(my_data)

    # Getting the temperature data
    temps_df = pd.read_csv("model/thermal/ANN/adjusted_TSD_df_all.csv")

    import matplotlib
    matplotlib.use('Agg') # noninteractive backend
    import matplotlib.pyplot as plt
    plt.close("all")
    import scipy.stats as stats

    # In training data, smallest T breakdown is 35, and largest T breakdown is 654.

    # use stats.gaussian_kde to estimate the probability density function from the histogram
    density = stats.gaussian_kde(temps_df['T'])
    x = np.arange(30,661,1) # in training data, smallest T breakdown is 35, and largest T breakdown is 654
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.plot(x, density(x))
    plt.plot(my_data, density(my_data), "or") # the current MOF's predicted thermal breakdown temperature

    ax.set_xlabel('Breakdown temperature (°C)')
    ax.set_ylabel('Frequency in the training data')
    if info['prediction']: # MOF wasn't in training data, and its ANN predicted breakdown temperature is used
        ax.set_title('Current MOF\'s predicted breakdown temperature relative to others')
    else: # MOF was in the training data, and its reported breakdown temperature is used
        ax.set_title('Current MOF\'s breakdown temperature relative to others')

    import mpld3

    return mpld3.fig_to_html(fig)

@app.route('/thermal_stability_percentile', methods=['POST']) 
def thermal_stability_percentile():
    # Returns what percentile the thermal breakdown temperature of the selected MOF lies in
    # with respect to the MOFs used to train the ANN for thermal stability predictions.

    # To begin, always go to main directory. 
    os.chdir(MOFSIMPLIFY_PATH)

    # Grab data.
    my_data = json.loads(flask.request.get_data()) # this is the current MOF's predicted thermal breakdown temperature
    my_data = my_data[:-3] # getting rid of the celsius symbol, left with just the number
    my_data = float(my_data)
    print(my_data)

    # Getting the temperature data.
    temps_df = pd.read_csv("model/thermal/ANN/adjusted_TSD_df_all.csv")

    # Will find what percentile our prediction belongs to, by checking the 100 percentiles and seeing which is closest to our prediction.
    difference = np.Infinity

    breakdown_Ts = temps_df['T']
    for i in np.arange(0,100.1,1): # 0,1, ..., 99, 100
        current_percentile = np.percentile(breakdown_Ts, i) # ith percentile
        current_difference = np.absolute(my_data - current_percentile) # absolute difference
        if current_difference < difference:
            difference = current_difference
            our_percentile = i

    our_percentile = int(our_percentile) # no decimal points

    return str(our_percentile)

# I use https://stackoverflow.com/questions/3252194/numpy-and-line-intersections to get code for the intersection of the two TGA lines
# for use in /TGA_plot call. perp and seg_intersect are from that stackoverflow page

#
# line segment intersection using vectors
# see Computer Graphics by F.S. Hill
#
def perp( a ) :
    b = np.empty_like(a)
    b[0] = -a[1]
    b[1] = a[0]
    return b

# line segment a given by endpoints a1, a2
# line segment b given by endpoints b1, b2
# return 
def seg_intersect(a1,a2, b1,b2) :
    da = a2-a1
    db = b2-b1
    dp = a1-b1
    dap = perp(da)
    denom = np.dot( dap, db)
    num = np.dot( dap, dp )
    return (num / denom.astype(float))*db + b1

# Makes the TGA plot for the current thermal ANN nearest neighbor.
@app.route('/TGA_plot', methods=['POST'])
def TGA_plot():

    # To begin, always go to main directory 
    os.chdir(MOFSIMPLIFY_PATH)

    # Grab data
    my_data = json.loads(flask.request.get_data()); # This is the neighbor complex
    my_data = my_data[:6] # only want the first six letters

    # Grab data 
    slopes_df = pd.read_csv("TGA/raw_TGA_digitization_data/digitized_csv/" + my_data + ".csv")

    x_values = []
    y_values = []
    for i in range(4): # 0, 1, 2, 3
        x_values.append(slopes_df.iloc[[i]]['T (degrees C)'][i])
        y_values.append(slopes_df.iloc[[i]]['mass (arbitrary units)'][i])

    # Making the four points.
    p1 = np.array( [x_values[0], y_values[0]] )
    p2 = np.array( [x_values[1], y_values[1]] )

    p3 = np.array( [x_values[2], y_values[2]] )
    p4 = np.array( [x_values[3], y_values[3]] )

    intersection_point = seg_intersect(p1, p2, p3, p4)

    # Instantiating the figure object. 
    graph = figure(title = "Simplified literature TGA plot of selected thermal ANN neighbor")  
         
    # The points to be plotted.
    xs = [[x_values[0], x_values[1],intersection_point[0]], [x_values[2], x_values[3],intersection_point[0]]] 
    ys = [[y_values[0], y_values[1],intersection_point[1]], [y_values[2], y_values[3],intersection_point[1]]] 
        
    # Plotting the graph.
    graph.multi_line(xs, ys) 
    graph.circle([intersection_point[0]], [intersection_point[1]], size=20, color="navy", alpha=0.5)
    graph.xaxis.axis_label = 'Temperature (°C)'
    graph.yaxis.axis_label = 'Percentage mass remaining or Mass'    

    return file_html(graph,CDN,'my plot')

@app.route('/get_components', methods=['POST']) 
def get_components():
    # Uses Aditya's MOF code to get linkers and sbus
    # Returns a dictionary with the linker and sbu xyz files's text, along with information about the number of linkers and sbus
    # Also in the dictionary: SMILES string for each of the linkers and sbus

    # To begin, always go to main directory.
    os.chdir(MOFSIMPLIFY_PATH);

    # Grab data
    my_data = json.loads(flask.request.get_data());

    os.chdir("temp_file_creation"); # changing directory

    # Write the data back to a cif file.
    cif_file = open('temp_cif.cif', 'w');
    cif_file.write(my_data);
    cif_file.close();

    # Delete the RACs folder, then remake it (to start fresh for this prediction).
    shutil.rmtree('RACs');
    os.mkdir('RACs');

    os.chdir("RACs"); # move to RACs folder

    # Next, running MOF featurization.
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
    # The fields are string representations of the linkers and sbus, however many there are.

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

    # Identifying which linkers and sbus have different connectivities.

    # Code below uses molecular graph determinants.
    # two ligands that are the same (via connectivity) will have the same molecular graph determinant. 
    # Molecular graph determinant fails to distinguish between isomers, but so would RCM (Reverse Cuthill McKee).

    # This simple script is meant to be run within the linkers directory, and it will give a bunch of numbers. 
    # If those numbers are the same, the linker is the same, if not, the linkers are different, etc.

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

    #### Linkers with the same molecular graph determinant are the same.
    #### Molecular graph determinant does not catch isomers.
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

    # adding the unique indices to the dictionary
    dictionary['unique_linker_indices'] = unique_linker_indices
    dictionary['unique_sbu_indices'] = unique_sbu_indices


    # In this next section, getting the SMILES strings for all of the linkers and sbus using pybel.
    # Write the smiles strings to file, then read the file.
    import pybel

    os.chdir("../linkers"); # move to linkers folder

    for i in range(dictionary['total_linkers']): # 0, 1, 2, ..., numberoflinkersminus1
        smilesFile = pybel.Outputfile('smi', 'temp_cif_primitive_linker_' + str(i) + '.txt') # smi refers to SMILES
        smilesFile.write(next(pybel.readfile('xyz', 'temp_cif_primitive_linker_' + str(i) + '.xyz'))) # writes SMILES string to the text file

        # Next, get the SMILES string from the text file.
        f = open('temp_cif_primitive_linker_' + str(i) + '.txt', 'r')
        line = f.readline()
        line = line.split('\t') # split at tabs
        smiles_ID = line[0]
        dictionary['linker_' + str(i) + '_SMILES'] = smiles_ID
        f.close()

    os.chdir("../sbus"); # move to sbus folder

    for i in range(dictionary['total_sbus']):
        smilesFile = pybel.Outputfile('smi', 'temp_cif_primitive_sbu_' + str(i) + '.txt') # smi refers to SMILES
        smilesFile.write(next(pybel.readfile('xyz', 'temp_cif_primitive_sbu_' + str(i) + '.xyz'))) # writes SMILES string to the text file

        # Next, get the SMILES string from the text file.
        f = open('temp_cif_primitive_sbu_' + str(i) + '.txt', 'r')
        line = f.readline()
        line = line.split('\t') # split at tabs
        smiles_ID = line[0]
        dictionary['sbu_' + str(i) + '_SMILES'] = smiles_ID
        f.close()


    json_object = json.dumps(dictionary, indent = 4);

    return json_object

@app.route('/solvent_neighbor_flag', methods=['POST']) 
def is_stable():
    # Returns the flag (whether or not stable upon solvent removal) and DOI of the neighbor sent over from the front end.

    # Grab data.
    my_data = json.loads(flask.request.get_data()); # This is the neighbor complex

    print(my_data)
    
    # To begin, always go to main directory. 
    os.chdir(MOFSIMPLIFY_PATH);

    os.chdir('model/solvent/ANN/dropped_connectivity_dupes')
    solvent_flags_df = pd.read_csv('train.csv')

    this_neighbor = solvent_flags_df[solvent_flags_df['CoRE_name'] == my_data] # getting the row with the MOF of interest

    this_neighbor_flag = this_neighbor['flag'] # getting the flag value

    this_neighbor_flag =  str(this_neighbor_flag.iloc[0]) # extract the flag value and return it as a string

    os.chdir(MOFSIMPLIFY_PATH);
    # next, getting DOI
    os.chdir('model/solvent/ANN/dropped_connectivity_dupes')
    solvent_flags_df = pd.read_csv('train.csv')
    this_neighbor = solvent_flags_df[solvent_flags_df['CoRE_name'] == my_data] # getting the row with the MOF of interest
    this_neighbor_doi = this_neighbor['doi'] # getting the doi
    this_neighbor_doi = this_neighbor_doi.iloc[0]

    myDict = {'flag': this_neighbor_flag, 'doi': this_neighbor_doi}

    return myDict

@app.route('/thermal_neighbor_T', methods=['POST']) 
def breakdown_T():
    # Returns the thermal breakdown temperature and DOI of the neighbor sent over from the front end.

    # Grab data
    my_data = json.loads(flask.request.get_data()); # This is the neighbor complex

    print(my_data)
    
    # To begin, always go to main directory 
    os.chdir(MOFSIMPLIFY_PATH);

    os.chdir('model/thermal/ANN')
    breakdown_T_df = pd.read_csv('train.csv')
    this_neighbor = breakdown_T_df[breakdown_T_df['CoRE_name'] == my_data] # getting the row with the MOF of interest
    this_neighbor_T = this_neighbor['T'] # getting the breakdown temperature value
    this_neighbor_T =  str(round(this_neighbor_T.iloc[0], 1)) # extract the breakdown temperature value and return it as a string. Want just one decimal place


    os.chdir(MOFSIMPLIFY_PATH);
    os.chdir('TGA')
    TGA_df = pd.read_excel('TGA_info_log.xlsx')
    this_neighbor = TGA_df[TGA_df['CoRE_name'] == my_data] # getting the row with the MOF of interest
    this_neighbor_doi = this_neighbor['doi'] # getting the doi
    this_neighbor_doi = this_neighbor_doi.iloc[0]

    myDict = {'T': this_neighbor_T, 'doi': this_neighbor_doi}

    return myDict

@app.route('/neighbor_writer', methods=['POST']) 
def neighbor_writer():
    # Writes information to a txt file about the selected latent space nearest neighbor.

    # Grab data
    my_data = json.loads(flask.request.get_data()); # This is a dictionary with information about the neighbor and the MOF that was analyzed by an ANN

    # To begin, always go to main directory. 
    os.chdir(MOFSIMPLIFY_PATH);

    print('neighbor writer check')
    print(my_data)

    prediction_type = my_data['prediction_type']
    current_MOF = my_data['current_MOF']
    selected_neighbor = my_data['selected_neighbor']
    latent_space_distance = my_data['latent_space_distance']

    # next, getting neighbor_truth
    if prediction_type == 'solvent':
        os.chdir('model/solvent/ANN/dropped_connectivity_dupes')
        solvent_flags_df = pd.read_csv('train.csv')
        this_neighbor = solvent_flags_df[solvent_flags_df['CoRE_name'] == selected_neighbor] # getting the row with the MOF of interest
        this_neighbor_flag = this_neighbor['flag'] # getting the flag value
        neighbor_truth = str(this_neighbor_flag.iloc[0]) # extract the flag value and convert it into a string
        if neighbor_truth == '1':
            neighbor_truth = 'Stable upon solvent removal'
        else:
            neighbor_truth = 'Unstable upon solvent removal'

    else: # thermal stability
        os.chdir('model/thermal/ANN')
        breakdown_T_df = pd.read_csv('train.csv')
        this_neighbor = breakdown_T_df[breakdown_T_df['CoRE_name'] == selected_neighbor] # getting the row with the MOF of interest
        this_neighbor_T = this_neighbor['T'] # getting the breakdown temperature value
        neighbor_truth = str(round(this_neighbor_T.iloc[0], 1)) # extract the flag value and convert it into a string. Want just one decimal place


    os.chdir(MOFSIMPLIFY_PATH);
    # next, getting DOI
    if prediction_type == 'solvent':
        os.chdir('model/solvent/ANN/dropped_connectivity_dupes')
        solvent_flags_df = pd.read_csv('train.csv')
        this_neighbor = solvent_flags_df[solvent_flags_df['CoRE_name'] == selected_neighbor] # getting the row with the MOF of interest
        this_neighbor_doi = this_neighbor['doi'] # getting the doi
        this_neighbor_doi = this_neighbor_doi.iloc[0]
    else: # thermal
        os.chdir('TGA')
        TGA_df = pd.read_excel('TGA_info_log.xlsx')
        this_neighbor = TGA_df[TGA_df['CoRE_name'] == selected_neighbor] # getting the row with the MOF of interest
        this_neighbor_doi = this_neighbor['doi'] # getting the doi
        this_neighbor_doi = this_neighbor_doi.iloc[0]

    # Now, writing all of this information to a string
    os.chdir(MOFSIMPLIFY_PATH);

    # Delete and remake the latent_neighbor folder each time, so only one file is ever in it
    shutil.rmtree('temp_file_creation/latent_neighbor')
    os.mkdir('temp_file_creation/latent_neighbor')

    with open('temp_file_creation/latent_neighbor/' + prediction_type + '__' + current_MOF + '__' + selected_neighbor + '.out','w') as f:
        f.write('Prediction type: ' + prediction_type + '\n')
        f.write('Current MOF: ' + current_MOF + '\n')
        f.write('Selected CoRE nearest neighbor in latent space: ' + selected_neighbor + '\n')
        f.write('Latent space distance: ' + latent_space_distance + '\n')
        f.write('Property for nearest neighbor: ' + neighbor_truth + '\n')
        f.write('Neighbor DOI: ' + this_neighbor_doi + '\n')

    print('neighbor writer check 2')

    return 'Success!'


@app.route('/TGA_maker', methods=['POST']) 
def TGA_maker():
    # Making the TGA plot and saving it, for it to be downloaded
    from bokeh.io import export_png

    # To begin, always go to main directory. 
    os.chdir(MOFSIMPLIFY_PATH);

    # Grab data
    my_data = json.loads(flask.request.get_data()); # This is the neighbor complex
    my_data = my_data[:6] # only want the first six letters

    print('neighbor writer check 3')

    # Grab data 
    slopes_df = pd.read_csv("TGA/raw_TGA_digitization_data/digitized_csv/" + my_data + ".csv")

    x_values = []
    y_values = []
    for i in range(4): # 0, 1, 2, 3
        x_values.append(slopes_df.iloc[[i]]['T (degrees C)'][i])
        y_values.append(slopes_df.iloc[[i]]['mass (arbitrary units)'][i])

    # Making the four points.
    p1 = np.array( [x_values[0], y_values[0]] )
    p2 = np.array( [x_values[1], y_values[1]] )

    p3 = np.array( [x_values[2], y_values[2]] )
    p4 = np.array( [x_values[3], y_values[3]] )

    intersection_point = seg_intersect(p1, p2, p3, p4)

    # Instantiating the figure object. 
    graph = figure(title = "Simplified literature TGA plot of selected thermal ANN neighbor")  
         
    # The points to be plotted.
    xs = [[x_values[0], x_values[1],intersection_point[0]], [x_values[2], x_values[3],intersection_point[0]]] 
    ys = [[y_values[0], y_values[1],intersection_point[1]], [y_values[2], y_values[3],intersection_point[1]]] 
        
    # Plotting the graph.
    graph.multi_line(xs, ys) 
    graph.circle([intersection_point[0]], [intersection_point[1]], size=20, color="navy", alpha=0.5)
    graph.xaxis.axis_label = 'Temperature (°C)'
    graph.yaxis.axis_label = 'Percentage mass remaining or Mass' 

    print('neighbor writer check 4')

    export_png(graph, filename='temp_file_creation/latent_neighbor/' + my_data + "_simplified_TGA.png")
    
    return 'Success!'

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8000)
