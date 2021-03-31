#!/usr/bin/env python
# -*- coding: utf-8 -*-
import flask
import pandas as pd
import keras
import numpy as np
import time
from keras import backend as K
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
# app.config['SECRET_KEY'] = "secretkeythatisverysecret,averysecretsecretkey"
# next, session code added by Gianmarco
# SESSION_TYPE = 'redis'
# app.config.from_object(__name__)
#sess = Session(app)

# Display text for neural network predictions
results_string='''<b>Results</b>
<br> ΔE<sub>H-L</sub> = %.1f kcal/mol predicted for the provided structure.
<br> This corresponds to S=%s being the ground state.
<br>
<br> <b>Value</b>
<br> When |ΔE<sub>H-L</sub>| < 5 kcal/mol, the complex is a potential spin-crossover complex (SCO), which can have important implications in catalytic reaction cycles.
<br> Our machine learning (ML) models [<a href="https://doi.org/10.1039/C7SC01247K">1</a>] took %.3f seconds to predict the geometry, spin splitting energy, and more for this complex!
<br> By comparison, density functional theory (DFT) optimizations on the two structures to get these properties would take about an hour using one Nvidia 970 GPU.
<br> molSimplify can help accelerate discovery of new materials, such as SCOs! [<a href="https://doi.org/10.1021/acs.jpclett.8b00170">2</a>]
<br>
<br> <b>Details</b>
<br> The ML models are trained on data from density functional theory (DFT) with the B3LYP functional, the LANL2DZ effective core potential for metal atoms, and the 6-31G* basis set for all other atoms.
<br> Warning: the training data do not include 4d transition metals, so results for these are extrapolations and may not be as accurate.
'''

licores = getlicores()

# Pre-load neural networks to save time later
valid_predictors = ['ls_ii', 'split']
predictor_vars = {}
# Load time is 0.5 seconds per predictor (only happens once, when server is started)
for predictor_name in valid_predictors:
    curr_dict = {}
    curr_dict['train_vars'] = ms_ANN.load_ANN_variables(predictor_name)
    curr_dict['norm_data'] = ms_ANN.load_normalization_data(predictor_name)
    curr_dict['my_ANN'] = ms_ANN.load_keras_ann(predictor_name)
    curr_dict['my_ANN']._make_predict_function()
    curr_dict['train_x'] = ms_ANN.load_training_data(predictor_name)
    predictor_vars[predictor_name] = curr_dict

# Pre-load PCA for plotting
pca_plot_df = pd.read_csv('data/pca_for_bokeh.csv')
with open('data/rac_names.txt','r') as file1:
    lines = file1.readlines()
pca_variables = [x.strip('\n') for x in lines]
# pca_model = pickle.load(open('data/PCA_model.pkl','rb')) # removed
with open('data/PCA_model.pkl','rb') as f: # next 4 lines were added by Gianmarco Terrones, for Python 3
    u = pickle._Unpickler(f)
    u.encoding = 'latin1'
    pca_model = u.load()

# KDE data for plotting
datadf = pd.read_csv('data/kde_data.csv')

def perform_ANN_prediction(RAC_dataframe, predictor_name, RAC_column='RACs'):
    # Function that performs ANN predictions (written quite inefficiently)
    assert type(RAC_dataframe) == pd.DataFrame
    start_time = time.time()

    curr_dict = predictor_vars[predictor_name]
    train_vars = curr_dict['train_vars']
    train_mean_x, train_mean_y, train_var_x, train_var_y = curr_dict['norm_data']
    my_ANN = curr_dict['my_ANN']
    train_x = curr_dict['train_x']
    #train_vars = load_ANN_variables(predictor_name)
    #train_mean_x, train_mean_y, train_var_x, train_var_y = load_normalization_data(predictor_name)
    #my_ANN = load_keras_ann(predictor_name)

    # Check if any RAC elements are missing from the provided dataframe
    missing_labels = [i for i in train_vars if i not in RAC_dataframe.columns]

    if len(missing_labels) > 0:
        # Try checking if there is anything in the column `RAC_column`. If so, deserialize it and re-run.
        if RAC_column in RAC_dataframe.columns:
            deserialized_RACs = pd.DataFrame.from_records(RAC_dataframe[RAC_column].values, index=RAC_dataframe.index)
            deserialized_RACs = deserialized_RACs.astype(float)
            RAC_dataframe = RAC_dataframe.join(deserialized_RACs)
            return perform_ANN_prediction(RAC_dataframe, predictor_name, RAC_column='RACs')
        else:
            raise ValueError('Please supply missing variables in your RAC dataframe: %s' % missing_labels)
    if 'alpha' in train_vars:
        if any(RAC_dataframe.alpha > 1):
            raise ValueError('Alpha is too large - should be between 0 and 1.')
    RAC_subset_for_ANN = RAC_dataframe.loc[:,train_vars].astype(float)
    normalized_input = ms_ANN.data_normalize(RAC_subset_for_ANN, train_mean_x, train_var_x)
    ANN_prediction = my_ANN.predict(normalized_input)
    rescaled_output = ms_ANN.data_rescale(ANN_prediction, train_mean_y, train_var_y)
    # Get latent vectors for training data and queried data
    #train_x = load_training_data(predictor_name)
    #train_x = pd.DataFrame(train_x, columns=train_vars).astype(float)
    #get_outputs = K.function([my_ANN.layers[0].input, K.learning_phase()],
    #                         [my_ANN.layers[len(my_ANN.layers) - 2].output])
    #normalized_train = ms_ANN.data_normalize(train_x, train_mean_x, train_var_x)
    #training_latent = get_outputs([normalized_train, 0])[0]
    #query_latent = get_outputs([normalized_input, 0])[0]


    # Append all results to dataframe
    results_allocation = [None for i in range(len(RAC_dataframe))]
    for i in range(len(RAC_dataframe)):
        results_dict = {}
        #min_latent_distance = min(np.linalg.norm(training_latent - query_latent[i][:], axis=1))
        #results_dict['%s_latent_vector' % predictor_name] = query_latent[i]
        #results_dict['%s_min_latent_distance' % predictor_name] = min_latent_distance
        output_value = rescaled_output[i]
        if len(output_value) == 1: # squash array of length 1 to the value it contains
            output_value = output_value[0]
        results_dict['%s_prediction' % predictor_name] = output_value
        results_allocation[i] = results_dict
    results_df = pd.DataFrame(results_allocation, index=RAC_dataframe.index)
    RAC_dataframe_with_results = RAC_dataframe.join(results_df)
    return RAC_dataframe_with_results

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

# @app.route('/get_example_mol')
# def serve_example_mol():
#     #print('testing testing') # debugging
#     #subprocess.run('pwd') # debugging
#     #subprocess.run('ls') # debugging

#     # TODO
#     return 'example placeholder'


    
@app.route('/predict_solvent_stability', methods=['POST']) # Gianmarco Terrones addition
def ss_predict():
    # Generates solvent stability prediction TODO
    # To do this, need to generate RAC featurization and Zeo++ geometry information for the MOF
    # Then, apply Aditya's model to make prediction

    # To begin, always go to main directory (this directory will vary depending on computer)
    os.chdir("/Users/gianmarcoterrones/Research/mofSimplify/")

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

    # Next, running MOF featurization
    #try:
    get_primitive('../temp_cif.cif', 'temp_cif_primitive.cif');
    # except ValueError:
    #     return 'FAILED'

    # try:
    full_names, full_descriptors = get_MOF_descriptors('temp_cif_primitive.cif',3,path= str(pathlib.Path().absolute()), xyzpath= 'temp_cif.xyz')
    # except ValueError:
    #     return 'FAILED'
    # except NotImplementedError:
    #     return 'FAILED'
    # except AssertionError:
    #     return 'FAILED'

    # if (len(full_names) <= 1) and (len(full_descriptors) <= 1): # this is a featurization check from MOF_descriptors.py
    #     return 'FAILED'

    # At this point, have the RAC featurization. Need geometry information next.

    # Run Zeo++
    os.chdir('..')
    os.chdir('zeo++-0.3')
    zeo_output = subprocess.check_output(['./network', '-ha', '-res', '../RACs/temp_cif_primitive.cif']) # Zeo++ command (see http://www.zeoplusplus.org/examples.html)
    py_file = open('txt_file_bin/geometry.txt', 'w')
    py_file.write(zeo_output.decode("utf-8")) # convert from bytes to string
    py_file.close()

    # TODO apply model next


    # will return a json object
    # one field will be the prediction
    # the other fields are string representations of the linkers and sbus, however many there are

    dictionary = {
        'prediction': 'test ss_predict' 
    }

    os.chdir('..')
    os.chdir('RACs')
    
    linker_num = 0;
    while True:
        if not os.path.exists('linkers/temp_cif_primitive_linker_' + str(linker_num) + '.xyz'):
            break
        else:
            linker_file = open('linkers/temp_cif_primitive_linker_' + str(linker_num) + '.xyz', 'r')
            linker_info = linker_file.read()
            linker_file.close()

            dictionary['linker_' + str(linker_num)] = linker_info

            linker_num = linker_num + 1;


    sbu_num = 0;
    while True:
        if not os.path.exists('sbus/temp_cif_primitive_sbu_' + str(sbu_num) + '.xyz'):
            break
        else:
            sbu_file = open('sbus/temp_cif_primitive_sbu_' + str(sbu_num) + '.xyz', 'r')
            sbu_info = sbu_file.read()
            sbu_file.close()

            dictionary['sbu_' + str(sbu_num)] = sbu_info

            sbu_num = sbu_num + 1


    dictionary['total_linkers'] = linker_num
    dictionary['total_sbus'] = sbu_num

    json_object = json.dumps(dictionary, indent = 4)

    return json_object

@app.route('/predict_thermal_stability', methods=['POST']) # Gianmarco Terrones addition
def ts_predict():
    # Generates thermal stability prediction TODO
    # To do this, need to generate RAC featurization and Zeo++ geometry information for the MOF
    # Then, apply Aditya's model to make prediction

    return 'test ts_predict'

@app.route('/generate', methods=['POST'])
def generate_mol():
    # Grab data
    mydata = json.loads(flask.request.get_data())

    good_request = True

    # Rename ligands appropriately
    ligand_rename_dict = {'pyridine': 'pyr', 
        'methyl isocyanide': 'misc',
        'phenyl isocyanide': 'pisc',
        'acetylacetonate': 'acac',
        'bipyridine': 'bipy'}

    for lig_type in ['ax1', 'ax2', 'eq']:
        if mydata[lig_type] in ligand_rename_dict:
            mydata[lig_type] = ligand_rename_dict[mydata[lig_type]]

    ligoccs = []
    ligs = []

    if mydata['eq'] in licores: # Check for higher denticity
        if isinstance(licores[mydata['eq']][2],list):
            n = int(4/len(licores[mydata['eq']][2]))
            for item in range(n):
                ligoccs.append(1)
                ligs.append(mydata['eq'])
        else: # Monodentate
            for item in range(4):
                ligoccs.append(1)
                ligs.append(mydata['eq'])
    else: # Force monodentate for now
        for item in range(4):
            ligoccs.append(1)
            ligs.append(mydata['eq'])

    ligs.append(mydata['ax1'])
    ligs.append(mydata['ax2'])
    ligoccs.append(1)
    ligoccs.append(1)

    ligoccs = ','.join([str(x) for x in ligoccs])
    if any([True for x in ligs if len(x) == 0]): # Catch cases where empty ligand passed.
        good_request = False

    ligs = ','.join([str(x) for x in ligs])
    print(ligs,ligoccs)
    # Generates an xyz file from parameters provided in POST request
    rundir = os.path.join(os.getcwd(), 'geos/')
    jobname = 'run_generate_query'
    mytext = string.Template('''-rundir $rundir
-ffoption no
-skipANN True
-distort 0
-name $jobname
-coord 6
-core $metal
-ligocc $ligoccs
-geometry oct
-spin $spin
-ligloc True
-oxstate $ox
-lig $ligs
-ff uff''')
    mytext = mytext.substitute(rundir=rundir, metal=mydata['metal'], spin=mydata['spin'], ox=mydata['ox'], ligoccs=ligoccs, ligs=ligs, jobname=jobname)
    inputstring_to_dict = lambda s: {i[0]:' '.join(i[1:]) for i in [i.split(' ') for i in s.split('\n')]}
    mytext_dict = inputstring_to_dict(mytext)
    try:
        strfiles, emsg, this_diag = startgen_pythonic(mytext_dict)
        print("printing xyz")
        # print(this_diag.mol.writexyz('', writestring=True))
        if good_request: 
            outstring = this_diag.mol.writexyz('', writestring=True)
        else:
            outstring = mydata['geometry']
    except UnboundLocalError:
        print('CAUGHT UNBOUND LOCAL ERROR')
        good_request = False
        outstring = mydata['geometry']

    http_response = {}
    http_response['geometry'] = outstring
    if good_request:
        http_response['message'] = "Geometry Successfully Generated!"
    else:
        http_response['message'] = "Please either specify a valid monodentate SMILES string or select a pre-populated ligand."

    return jsonify(http_response)
    



@app.route('/nn_predict', methods=['POST'])
def serve_nn_prediction():
    # Run an NN prediction for the xyz file, ox state, and HFX contained in the POST request
    try:
        mydata = json.loads(flask.request.get_data())
        xyzstr = mydata['geometry']
        mymol = ms_mol3D.mol3D()
        mymol.readfromstring(xyzstr)
        rac_names, rac_vals = ms_RAC.get_descriptor_vector(mymol)
        mydict = {rac_name:rac_val for rac_name, rac_val in zip(rac_names, rac_vals)}
        mydict['alpha'] = float(mydata['hfx'])
        mydict['ox'] = int(mydata['ox'])
        mytable = pd.DataFrame([pd.Series(mydict)])
	    # Local ANN prediction
        start_time = time.time()
        mytable = perform_ANN_prediction(mytable, 'split')
        SSE_prediction = mytable.split_prediction.iloc[0]

	    # Sketchy math to find ground state spin
        quantum_spins = {1: '0', 2: '1/2', 3: '1', 4: '3/2', 5: '2', 6: '5/2'}
        mult = int(mydata['spin'])
        ground_mult = 6 - mult if SSE_prediction < 0 else mult
        quantum_spin_str = quantum_spins[ground_mult]
        time_taken = time.time() - start_time
        myresult = results_string % (SSE_prediction, quantum_spin_str, time_taken)
        http_response = {}
        http_response['resulttext'] = myresult
        http_response['result'] = str(SSE_prediction) 
        return jsonify(http_response)
    except:
        http_response = {}
        http_response['resulttext'] = 'Error'
        http_response['result'] = str(-24.5) 
        return jsonify(http_response)

# Version 1 - KDE
def make_plot(x1):
    source = ColumnDataSource(datadf)
    print(datadf.columns.values)
    x_min = min(x1, min(datadf.sse.values))
    x_max = max(x1, max(datadf.sse.values))
    x_range = x_max - x_min
    x_min, x_max = (x_min - 0.05*x_range, x_max + 0.05*x_range)
    p = figure(title = "Your Complex Spin Splitting Energy Compared to Training Data", 
               sizing_mode="scale_both", plot_width=400, plot_height=200, x_axis_type=None, x_range=(x_min, x_max))
    p.yaxis.axis_label = 'KDE of Training Data'
    p.line('sse', 'kde', source=source, line_width=2)
    vertical = Span(location=x1,dimension='height',
                    line_color='green',line_dash='dashed',line_width=3)
    p.yaxis.major_tick_line_color = None  # turn off y-axis major ticks
    p.yaxis.minor_tick_line_color = None  # turn off y-axis minor ticks
    p.yaxis.major_label_text_font_size = '0pt'  # turn off y-axis tick labels
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.add_layout(vertical)
    ticker = SingleIntervalTicker(interval=10, num_minor_ticks=2)
    xaxis = LinearAxis(ticker=ticker)
    p.add_layout(xaxis, 'below')
    p.xaxis.axis_label = 'Spin Splitting Energy (kcal/mol)'
    return p

# Version 2 - PCA
def make_plotpca(sse,geometry,ox,alpha):
    print(pca_plot_df.columns.values)
    ##### Apply Data
    my_mol = ms_mol3D.mol3D()
    my_mol.readfromstring(geometry)
    RACs = my_mol.get_features()
    RACs['ox'] = ox
    RACs['alpha'] = alpha
    new_vect = np.array([RACs[x] for x in pca_variables])
    plot_vals = pca_model.transform(np.array(new_vect).reshape(1,-1))
    x1 = plot_vals[0][0] # New X
    y1 = plot_vals[0][1] # New Y
    ########
    x_min = min(x1, min(pca_plot_df.PC1.values))
    x_max = max(x1, max(pca_plot_df.PC1.values))
    y_min = min(y1, min(pca_plot_df.PC2.values))
    y_max = max(y1, max(pca_plot_df.PC2.values))
    x_range = x_max - x_min
    y_range = y_max - y_min
    x_min, x_max = (x_min - 0.05*x_range, x_max + 0.05*x_range)
    y_min, y_max = (y_min - 0.05*y_range, y_max + 0.05*y_range)
    #### Colormap
    colormax = max(np.abs(sse),np.abs(pca_plot_df['Spin Splitting Energy (kcal/mol)'].values).max())
    colormax = 50
    mapper = LinearColorMapper(palette=cmap_bokeh, low=-colormax, high=colormax)
    colors = {'field':np.array(list(pca_plot_df['Spin Splitting Energy (kcal/mol)'].values) + [sse]),
                'transform':mapper}
    cticker = SingleIntervalTicker(interval=10, num_minor_ticks=0)
    color_bar = ColorBar(color_mapper=mapper, ticker=cticker, label_standoff=5, 
                         border_line_color=None, location=(0,0))#, title='SSE (kcal/mol)')
    ##### Plotting
    source = ColumnDataSource({
        'PC1':pca_plot_df.PC1.values,
        'PC2':pca_plot_df.PC2.values,
        'chemname':pca_plot_df['Chemical Name'].values,
        'color': pca_plot_df['Spin Splitting Energy (kcal/mol)'].values})
    TOOLS="crosshair,pan,wheel_zoom,zoom_in,zoom_out,box_zoom,reset,tap,"
    p = figure(title = "Your Complex (green) Compared to Training Data in Feature Space", 
               sizing_mode="scale_both", plot_width=400, plot_height=400,
                x_range=(x_min, x_max),y_range=(y_min, y_max),tools=TOOLS)
    # p.line('sse', 'kde', source=source, line_width=2)
    source2 = ColumnDataSource({
        'PC1':[x1],
        'PC2':[y1],
        'chemname':['Your Complex'],
        'color': [sse]})
    p.circle('PC1','PC2',source=source,fill_color={'field':'color','transform':mapper},line_color={'field':'color','transform':mapper},
             size=10)
    p.circle('PC1','PC2',source=source2,fill_color={'field':'color','transform':mapper},line_color='mediumseagreen',
             size=20,line_width=3)
    # Datalabels with hover animation.
    p.add_tools(HoverTool(
        tooltips = [('Chemical Name', '@chemname{%s}'),
                     ('Spin Splitting Energy (kcal/mol)', '@color{%0.2d}'),
                     ],
        formatters={ 'chemname':'printf',
            'color': 'printf'
        }
    ))
    p.yaxis.major_tick_line_color = None  # turn off y-axis major ticks
    p.yaxis.minor_tick_line_color = None  # turn off y-axis minor ticks
    p.xaxis.major_tick_line_color = None  # turn off x-axis major ticks
    p.xaxis.minor_tick_line_color = None  # turn off x-axis minor ticks
    p.yaxis.major_label_text_font_size = '0pt'  # turn off y-axis tick labels
    p.xaxis.major_label_text_font_size = '0pt' # turn off x-axis tick labels
    p.xgrid.grid_line_color = None # Remove grid
    p.ygrid.grid_line_color = None # Remove grid
    p.outline_line_width = 2
    p.outline_line_alpha = 1
    p.outline_line_color = "black"
    p.xaxis.axis_label = 'PC1'
    p.yaxis.axis_label = 'PC2'
    p.add_layout(color_bar,'right')
    return p

@app.route('/plot', methods=['POST'])
def plot():
    mydata = json.loads(flask.request.get_data())
    plt = make_plot(float(mydata['sseplot'])) # V1 KDE
    return file_html(plt,CDN,'my plot')

@app.route('/plotpca', methods=['POST'])
def plotpca():
    mydata = json.loads(flask.request.get_data())
    plt = make_plotpca(float(mydata['sseplot']),mydata['geometry'],
                    mydata['ox'],mydata['hfx'])
    return file_html(plt,CDN,'my plot pca')

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8000)
