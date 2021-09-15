#!/usr/bin/env python
# -*- coding: utf-8 -*-
import flask
import tensorflow as tf
import pandas as pd
import numpy as np
import time
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
import stat
import keras
import keras.backend as K
import sklearn
import sklearn.preprocessing
from sklearn.metrics import pairwise_distances
from molSimplify.Scripts.generator import startgen_pythonic
from molSimplify.Scripts.molSimplify_io import getlicores
from bokeh.plotting import figure
from bokeh.resources import CDN
from bokeh.models import ColumnDataSource, SingleIntervalTicker, LinearAxis
from bokeh.embed import file_html
from bokeh.models import Span
from bokeh.models import ColorBar, LinearColorMapper, LogColorMapper, HoverTool
from bokeh.models.markers import Circle
from bokeh.palettes import Inferno256
from flask import jsonify, render_template, redirect, request, url_for, session
from functools import partial
from keras.callbacks import EarlyStopping
import flask_login
from flask_login import LoginManager, UserMixin, login_required, current_user
from molSimplify.Informatics.MOF.MOF_descriptors import get_primitive, get_MOF_descriptors
from flask_cors import CORS
import bson
from datetime import datetime
from pymongo import MongoClient
from werkzeug.utils import secure_filename

cmap_bokeh = Inferno256

MOFSIMPLIFY_PATH = os.path.abspath('.') # the main directory
MOFSIMPLIFY_PATH += '/'
USE_SPLASH_PAGE = False

app = flask.Flask(__name__)
app.secret_key = 'hjk_secret_key_mofsimplify_2021' # secret key
app.config['MAX_CONTENT_LENGTH'] = 1024 * 1024 * 20 # Upload max 20 megabytes
app.config['UPLOAD_EXTENSIONS'] = ['.jpg', '.jpeg', '.png', '.pdf', '.tiff', '.tif', '.eps']
cors = CORS(app)

### splash page management: https://stackoverflow.com/questions/37275262/anonym-password-protect-pages-without-username-with-flask
login_manager = LoginManager()
login_manager.init_app(app)
users = {'user1':{'password':'MOFSimplify!Beta2021'}}

def precision(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
    precision = true_positives / (predicted_positives + K.epsilon())
    return precision
def recall(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    recall = true_positives / (possible_positives + K.epsilon())
    return recall
def f1(y_true, y_pred):
    p = precision(y_true, y_pred)
    r = recall(y_true, y_pred)
    return 2 * ((p * r) / (p + r + K.epsilon()))

# global variables for the ANN models, loaded once when server first started
solvent_ANN_path = MOFSIMPLIFY_PATH + 'model/solvent/ANN/' 
thermal_ANN_path = MOFSIMPLIFY_PATH + 'model/thermal/ANN/'
dependencies = {'precision':precision,'recall':recall,'f1':f1}
# loading the ANN models to save time later
tf_session = tf.Session()
from tensorflow import keras as tf_keras
tf_keras.backend.set_session(tf_session)
solvent_model = keras.models.load_model(solvent_ANN_path + 'final_model_flag_few_epochs.h5',custom_objects=dependencies)
thermal_model = keras.models.load_model(thermal_ANN_path + 'final_model_T_few_epochs.h5',custom_objects=dependencies)

class User(UserMixin):
  pass

@login_manager.user_loader
def user_loader(username):
  if username not in users:
    return
  user = User()
  user.id = username
  return user

@login_manager.request_loader
def request_loader(request):
  username = request.form.get('username')
  if username not in users:
    return
  user = User()
  user.id = username

  user.is_authenticated = request.form['password'] == users[username]['password']

  return user

### End of splash page management

@app.route('/new_user', methods=['GET'])
def set_ID():
    # sets the session user ID. This is used to generate unique folders, so that multiple users can use the website at a time 
    # specifically, copies of the temp_file_creation folder
    session['ID'] = time.time() # a unique ID for this session
    # print('MY ID CHECK 1')
    # print(session['ID'])

    # make a version of the temp_file_creation folder for this user
    new_folder = MOFSIMPLIFY_PATH + '/temp_file_creation_' + str(session['ID'])
    shutil.copytree('temp_file_creation', new_folder)
    os.remove(new_folder + '/temp_cif.cif') # remove this, for sole purpose of updating time stamp on the new folder (copytree doesn't)

    # print('walk check')
    # delete all temp_file_creation clone folders that haven't been used for a while, to prevent folder accumulation
    for root, dirs, files in os.walk(MOFSIMPLIFY_PATH):
        for dir in dirs:
            # print(dir)
            target_str = 'temp_file_creation'
            if len(dir) > len(target_str) and target_str in dir and file_age_in_seconds(dir) > 7200: # 7200s is two hours
                # target_str in dir since all copies start with temp_file_creation
                # len(dir) > len(target_str) to prevent deleting the original temp_file_creation folder
                shutil.rmtree(dir)


    return str(session['ID']) # return a string

@app.route('/get_ID', methods=['GET'])
def get_ID():
    # gets the session user ID. This is used for getting building block generated MOFs
    return str(session['ID']) # return a string


@app.route('/demo')
def serve_demo():
    # Go to localhost:8000/demo to see this.
    return 'you are on the demo page'

@app.route('/logo.png')
def serve_logo():
    return flask.send_from_directory('images', 'logo.png')

@app.route('/truncated_logo.png')
def serve_truncated_logo():
    return flask.send_from_directory('images', 'truncated_logo.png')

@app.route('/TGA_graphic.png')
def serve_TGA_graphic():
    return flask.send_from_directory('images', 'TGA_graphic.png')

@app.route('/MOF5_background.png')
def serve_MOF5():
    return flask.send_from_directory('images', 'MOF5_background.png')

@app.route('/banner_light')
def serve_banner_light():
    return flask.send_from_directory('images', 'MOF_light.webp') # google's webp format. It is optimized for websites and loads fast

@app.route('/banner_dark')
def serve_banner_dark():
    return flask.send_from_directory('images', 'MOF_dark.webp') # google's webp format. It is optimized for websites and loads fast

@app.route('/MOF_logo.png')
def serve_MOFSimplify_logo():
    return flask.send_from_directory('images', 'MOF_logo.png')

@app.route('/images/bg.jpg')
def serve_bg():
    # Hack to show the background image on the success/failure screens.
    return flask.send_from_directory('./splash_page/images', 'bg.jpg')

## Handle feedback
@app.route('/process_feedback', methods=['POST'])
def process_feedback():
    client = MongoClient('18.18.63.68',27017) # connect to mongodb
    db = client.feedback
    collection = db.MOFSimplify
    fields = ['feedback_form_name', 'rating', 'email', 'reason', 'comments', 'cif_file_name', 'structure']
    #$meta_fields = ['IP', 'datetime', 'cif_file', 'MOF_name']
    final_dict = {}
    for field in fields:
        final_dict[field] = request.form.get(field)

    # Populate special fields
    uploaded_file = request.files['file']
    final_dict['filetype'] = uploaded_file.content_type
    filename = secure_filename(uploaded_file.filename)
    final_dict['filename'] = filename
    final_dict['file'] = uploaded_file.read()
    file_ext = os.path.splitext(filename)[1].lower()
    if file_ext not in app.config['UPLOAD_EXTENSIONS']:
        return flask.send_from_directory('./splash_page/', 'error.html')
    final_dict['ip'] = request.remote_addr
    final_dict['timestamp'] = datetime.now().isoformat()
    
    print(final_dict)
    collection.insert(final_dict) # insert the dictionary into the mongodb collection
    with open('sample.bson', 'wb') as outfile:
        outfile.write(bson.encode(final_dict))
    return flask.send_from_directory('./splash_page/', 'success.html')


## Splash page management
@app.route('/', methods=['GET', 'POST'])
@app.route('/<path:path>', methods=['GET', 'POST'])
def index(path='index.html'):
  if request.method == 'POST':
    username = 'user1'
    if request.form.get('password') == users[username]['password']:
      user = User()
      user.id = 'user1'
      flask_login.login_user(user)
  # print('is user authenticated?')
  # print(current_user.is_authenticated)
  # print('input check')
  # print(request.form.get('password'))
  if current_user.is_authenticated:
    return flask.send_from_directory('.', 'index.html')
  elif request.form.get('password') == None:
    return flask.send_from_directory('./splash_page/', path)
  else:
    return flask.send_from_directory('./splash_page/', 'index_wrong_password.html')

@app.route('/about.html')
def serve_about():
    # Serves homepage
    return flask.send_from_directory('.', 'about.html')

@app.route('/mof_examples/<path:path>') # needed for fetch
def serve_example(path):
    # Serves example
    return flask.send_from_directory('mof_examples', path)

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

    path_parts = path.split('~')
    cif_name = path_parts[0]
    user_ID = path_parts[1]
    # print('ID check')
    # print(session)
    # print(session['ID'])
    return flask.send_from_directory('temp_file_creation_' + user_ID + '/tobacco_3.0/output_cifs', cif_name);

@app.route('/neighbor/<path:path>') # needed for fetch
def serve_neighbor(path):
    # Serves the neighbor CoRE MOF
    return flask.send_from_directory('CoRE2019', path);

# @app.route('/neighbor_info/<path:path>')
# def serve_neighbor_txt(path):
#     # Serves the neighbor CoRE MOF information in txt file format
#     return flask.send_from_directory('temp_file_creation_' + str(session['ID']) + '/latent_neighbor', path);

@app.route('/merged_descriptors/<path:path>')
def serve_descriptors(path):
    # Serves the csv of descriptors for the selected MOF
    return flask.send_from_directory('temp_file_creation_' + str(session['ID']) + '/merged_descriptors', path);

def listdir_nohidden(path): # used for bb_generate. Ignores hidden files
    myList = os.listdir(path);
    for i in myList:
        if i.startswith('.'):
            myList.remove(i)
    return myList

def file_age_in_seconds(pathname): 
    # print('age_check')
    # print(time.time() - os.stat(pathname)[stat.ST_MTIME]) # time since last modification
    return time.time() - os.stat(pathname)[stat.ST_MTIME]

@app.route('/curr_users', methods=['GET'])
def curr_num_users():
    # Returns the current number of users on MOFSimplify. 

    sum = 0

    for root, dirs, files in os.walk(MOFSIMPLIFY_PATH):
        for dir in dirs:
            target_str = 'temp_file_creation'
            if len(dir) > len(target_str) and target_str in dir: 
                # target_str in dir since all copies start with temp_file_creation
                # len(dir) > len(target_str) to prevent counting the original temp_file_creation folder
                sum += 1

    return str(sum+1)

@app.route('/get_bb_generated_MOF', methods=['POST']) 
def bb_generate():
    # generates a MOF using the building blocks and net specified by the user
    # uses ToBaCCo code, version 3.0
    # returns the constructed MOF's name to the front end

    tobacco_folder = MOFSIMPLIFY_PATH + "temp_file_creation_" + str(session['ID']) + "/tobacco_3.0/"

    # Grab data
    my_data = json.loads(flask.request.get_data())

    linker = my_data['linker']
    sbu = my_data['sbu']
    net = my_data['net']

    # clear the edges, nodes, templates, and output cifs folders to start fresh
        # when running python tobacco.py, it looks in these folders

    shutil.rmtree(tobacco_folder + 'edges')
    os.mkdir(tobacco_folder + 'edges')
    shutil.rmtree(tobacco_folder + 'nodes')
    os.mkdir(tobacco_folder + 'nodes')
    shutil.rmtree(tobacco_folder + 'templates')
    os.mkdir(tobacco_folder + 'templates')
    shutil.rmtree(tobacco_folder + 'output_cifs')
    os.mkdir(tobacco_folder + 'output_cifs')

    # copy over the linker, sbu, and net specified by the user in the edges, nodes, and templates folders

    shutil.copy(tobacco_folder + 'edges_database/' + linker + '.cif', tobacco_folder + 'edges/' + linker + '.cif')
    shutil.copy(tobacco_folder + 'nodes_database/' + sbu + '.cif', tobacco_folder + 'nodes/' + sbu + '.cif')
    shutil.copy(tobacco_folder + 'template_database/' + net + '.cif', tobacco_folder + 'templates/' + net + '.cif')


    # run the command to construct the MOF
    os.chdir(tobacco_folder) 
    # note: os.chdir here could cause issues if multiple users are using the website and try to make a building block generated MOF at the same time, since MOFSimplify server might chdir when someone else is in the middle of an operation
    # luckily, it is a quick operation, so this is unlikely
    subprocess.run(['python', 'tobacco.py']) 
    os.chdir(MOFSIMPLIFY_PATH)

    # if successful, there will be an output cif in the folder output_cifs
    if listdir_nohidden(tobacco_folder + 'output_cifs') == []: # no files in folder
        print('Construction failed.')
        return 'FAILED'
    
    constructed_MOF = listdir_nohidden(tobacco_folder + 'output_cifs')
    constructed_MOF = constructed_MOF[0] # getting the first, and only, element out of the list

    dictionary = {};

    dictionary['mof_name'] = constructed_MOF

    # getting the primitive cell using molSimplify
    get_primitive(tobacco_folder + 'output_cifs/' + constructed_MOF, tobacco_folder + 'output_cifs/primitive_' + constructed_MOF);

    json_object = json.dumps(dictionary, indent = 4);

    return json_object


### Next, the prediction functions ### 

# This function takes in two dataframes df_train and df_newMOF, one for the training data (many rows) and one for the new MOF (one row) for which a prediction is to be generated.
# This function also takes in fnames (the feature names) and lname (the target property name).
# This function normalizes the X values from the pandas dataframes and returns them as X_train and X_newMOF.
# It also "normalizes" y_train, which are the solvent stability flags in the training data dataframe, and returns x_scaler (which scaled X_train).
def normalize_data_solvent(df_train, df_newMOF, fnames, lname, unit_trans=1, debug=False):
    _df_train = df_train.copy().dropna(subset=fnames+lname)
    _df_newMOF = df_newMOF.copy().dropna(subset=fnames) 
    X_train, X_newMOF = _df_train[fnames].values, _df_newMOF[fnames].values # takes care of ensuring ordering is same for both X
    y_train = _df_train[lname].values
    if debug:
        print("training data reduced from %d -> %d because of nan." % (len(df_train), y_train.shape[0]))
    x_scaler = sklearn.preprocessing.StandardScaler()
    x_scaler.fit(X_train)
    X_train = x_scaler.transform(X_train)
    X_newMOF = x_scaler.transform(X_newMOF)
    y_train = np.array([1 if x == 1 else 0 for x in y_train.reshape(-1, )])
    return X_train, X_newMOF, y_train, x_scaler

def standard_labels(df, key="flag"):
    flags = [1 if row[key] == 1 else 0 for _, row in df.iterrows()]
    df[key] = flags
    return df

def run_solvent_ANN(user_id, path, MOF_name, solvent_ANN):

    RACs = ['D_func-I-0-all','D_func-I-1-all','D_func-I-2-all','D_func-I-3-all',
     'D_func-S-0-all', 'D_func-S-1-all', 'D_func-S-2-all', 'D_func-S-3-all',
     'D_func-T-0-all', 'D_func-T-1-all', 'D_func-T-2-all', 'D_func-T-3-all',
     'D_func-Z-0-all', 'D_func-Z-1-all', 'D_func-Z-2-all', 'D_func-Z-3-all',
     'D_func-chi-0-all', 'D_func-chi-1-all', 'D_func-chi-2-all',
     'D_func-chi-3-all', 'D_lc-I-0-all', 'D_lc-I-1-all', 'D_lc-I-2-all',
     'D_lc-I-3-all', 'D_lc-S-0-all', 'D_lc-S-1-all', 'D_lc-S-2-all',
     'D_lc-S-3-all', 'D_lc-T-0-all', 'D_lc-T-1-all', 'D_lc-T-2-all',
     'D_lc-T-3-all', 'D_lc-Z-0-all', 'D_lc-Z-1-all', 'D_lc-Z-2-all',
     'D_lc-Z-3-all', 'D_lc-chi-0-all', 'D_lc-chi-1-all', 'D_lc-chi-2-all',
     'D_lc-chi-3-all', 'D_mc-I-0-all', 'D_mc-I-1-all', 'D_mc-I-2-all',
     'D_mc-I-3-all', 'D_mc-S-0-all', 'D_mc-S-1-all', 'D_mc-S-2-all',
     'D_mc-S-3-all', 'D_mc-T-0-all', 'D_mc-T-1-all', 'D_mc-T-2-all',
     'D_mc-T-3-all', 'D_mc-Z-0-all', 'D_mc-Z-1-all', 'D_mc-Z-2-all',
     'D_mc-Z-3-all', 'D_mc-chi-0-all', 'D_mc-chi-1-all', 'D_mc-chi-2-all',
     'D_mc-chi-3-all', 'f-I-0-all', 'f-I-1-all', 'f-I-2-all', 'f-I-3-all',
     'f-S-0-all', 'f-S-1-all', 'f-S-2-all', 'f-S-3-all', 'f-T-0-all', 'f-T-1-all',
     'f-T-2-all', 'f-T-3-all', 'f-Z-0-all', 'f-Z-1-all', 'f-Z-2-all', 'f-Z-3-all',
     'f-chi-0-all', 'f-chi-1-all', 'f-chi-2-all', 'f-chi-3-all', 'f-lig-I-0',
     'f-lig-I-1', 'f-lig-I-2', 'f-lig-I-3', 'f-lig-S-0', 'f-lig-S-1', 'f-lig-S-2',
     'f-lig-S-3', 'f-lig-T-0', 'f-lig-T-1', 'f-lig-T-2', 'f-lig-T-3', 'f-lig-Z-0',
     'f-lig-Z-1', 'f-lig-Z-2', 'f-lig-Z-3', 'f-lig-chi-0', 'f-lig-chi-1',
     'f-lig-chi-2', 'f-lig-chi-3', 'func-I-0-all', 'func-I-1-all',
     'func-I-2-all', 'func-I-3-all', 'func-S-0-all', 'func-S-1-all',
     'func-S-2-all', 'func-S-3-all', 'func-T-0-all', 'func-T-1-all',
     'func-T-2-all', 'func-T-3-all', 'func-Z-0-all', 'func-Z-1-all',
     'func-Z-2-all', 'func-Z-3-all', 'func-chi-0-all', 'func-chi-1-all',
     'func-chi-2-all', 'func-chi-3-all', 'lc-I-0-all', 'lc-I-1-all', 'lc-I-2-all',
     'lc-I-3-all', 'lc-S-0-all', 'lc-S-1-all', 'lc-S-2-all', 'lc-S-3-all',
     'lc-T-0-all', 'lc-T-1-all', 'lc-T-2-all', 'lc-T-3-all', 'lc-Z-0-all',
     'lc-Z-1-all', 'lc-Z-2-all', 'lc-Z-3-all', 'lc-chi-0-all', 'lc-chi-1-all',
     'lc-chi-2-all', 'lc-chi-3-all', 'mc-I-0-all', 'mc-I-1-all', 'mc-I-2-all',
     'mc-I-3-all', 'mc-S-0-all', 'mc-S-1-all', 'mc-S-2-all', 'mc-S-3-all',
     'mc-T-0-all', 'mc-T-1-all', 'mc-T-2-all', 'mc-T-3-all', 'mc-Z-0-all',
     'mc-Z-1-all', 'mc-Z-2-all', 'mc-Z-3-all', 'mc-chi-0-all', 'mc-chi-1-all',
     'mc-chi-2-all', 'mc-chi-3-all']
    geo = ['Df','Di', 'Dif','GPOAV','GPONAV','GPOV','GSA','POAV','POAV_vol_frac',
      'PONAV','PONAV_vol_frac','VPOV','VSA','rho']
     
    other = ['cif_file','name','filename']

    ANN_path = path + 'model/solvent/ANN/'
    temp_file_path = path + 'temp_file_creation_' + user_id + '/'
    df_train = pd.read_csv(ANN_path+'dropped_connectivity_dupes/train.csv')
    df_train = df_train.loc[:, (df_train != df_train.iloc[0]).any()]
    df_newMOF = pd.read_csv(temp_file_path + 'merged_descriptors/' + MOF_name + '_descriptors.csv') # assumes that temp_file_creation/ is in parent folder
    features = [val for val in df_train.columns.values if val in RACs+geo]

    df_train = standard_labels(df_train, key="flag")

    ### The normalize_data_solvent function is expecting a dataframe with each MOF in a separate row, and features in columns
    ### At this location, use get_MOF_descriptors to get features
        # Look at the files that are generated: lc_descriptors.csv, sbu_descriptors.csv, linker_descriptors.csv
    ### Then store those features in a usable form (df)
    ### Need to merge with geometry features from Zeo++
        # done in app.py

    ### Utilize the function below to normalize the RACs + geos of the new MOF
    # newMOF refers to the MOF that has been uploaded to mofSimplify, for which a prediction will be generated
    X_train, X_newMOF, y_train, x_scaler = normalize_data_solvent(df_train, df_newMOF, features, ["flag"], unit_trans=1, debug=False)
    # Order of values in X_newMOF matters, but this is taken care of in normalize_data_solvent.
    X_train.shape, y_train.reshape(-1, ).shape
    model = solvent_ANN

    from tensorflow.python.keras.backend import set_session
    with tf_session.as_default(): # session stuff is needed because the model was loaded from h5 a while ago
        with tf_session.graph.as_default():
            ### new_MOF_pred will be a decimal value between 0 and 1, below 0.5 is unstable, above 0.5 is stable
            new_MOF_pred = np.round(model.predict(X_newMOF),2) # round to 2 decimals
            # print('success!')

            # Define the function for the latent space. This will depend on the model. We want the layer before the last, in this case this was the 12th one.
            get_latent = K.function([model.layers[0].input],
                                    [model.layers[12].output]) # Last layer before dense-last

            # Get the latent vectors for the training data first, then the latent vectors for the test data.
            training_latent = get_latent([X_train, 0])[0]
            design_latent = get_latent([X_newMOF, 0])[0]

            # print(training_latent.shape,design_latent.shape)

    # Compute the pairwise distances between the test latent vectors and the train latent vectors to get latent distances
    d1 = pairwise_distances(design_latent,training_latent,n_jobs=30)
    df1 = pd.DataFrame(data=d1, columns=df_train['CoRE_name'].tolist())
    df1.to_csv(temp_file_path + 'solvent_test_latent_dists.csv')

    # Want to find the closest points (let's say the closest 5 points); so, smallest values in df1
    neighbors = 5 # number of closest points

    # will make arrays of length neighbors, where each entry is the next closest neighbor (will do this for both names and distances)
    neighbors_names = []
    neighbors_distances = []

    df_reformat = df1.min(axis='index')

    for i in range(neighbors):
        name = df_reformat.idxmin() # name of next closest complex in the traiing data
        distance = df_reformat.min() # distance of the next closest complex in the training data to the new MOF
        df_reformat = df_reformat.drop(name) # dropping the next closest complex, in order to find the next-next closest complex

        neighbors_names.append(name)
        neighbors_distances.append(str(distance))

    return str(new_MOF_pred[0][0]), neighbors_names, neighbors_distances

# This function takes in two dataframes df_train and df_newMOF, one for the training data (many rows) and one for the new MOF (one row) for which a prediction is to be generated.
# This function also takes in fnames (the feature names) and lname (the target property name).
# This function normalizes the X values from the pandas dataframes and returns them as X_train and X_newMOF.
# It also normalizes y_train, which are the thermal breakdown temperatures in the training data dataframe, and returns x_scaler (which scaled X_train) and y_scaler (which scaled y_train).
def normalize_data_thermal(df_train, df_newMOF, fnames, lname, unit_trans=1, debug=False): # assumes gets Pandas dataframes with MOFs as rows and features as columns
    _df_train = df_train.copy().dropna(subset=fnames+lname)
    _df_newMOF = df_newMOF.copy().dropna(subset=fnames) 
    X_train, X_newMOF = _df_train[fnames].values, _df_newMOF[fnames].values # takes care of ensuring ordering is same for both X
    y_train = _df_train[lname].values
    if debug:
        print("training data reduced from %d -> %d because of nan." % (len(df_train), y_train.shape[0]))
    x_scaler = sklearn.preprocessing.StandardScaler()
    x_scaler.fit(X_train)
    X_train = x_scaler.transform(X_train)
    X_newMOF = x_scaler.transform(X_newMOF)
    y_scaler = sklearn.preprocessing.StandardScaler()
    y_scaler.fit(y_train)
    y_train = y_scaler.transform(y_train)
    return X_train, X_newMOF, y_train, x_scaler, y_scaler

def run_thermal_ANN(user_id, path, MOF_name, thermal_ANN):

    RACs = ['D_func-I-0-all','D_func-I-1-all','D_func-I-2-all','D_func-I-3-all',
     'D_func-S-0-all', 'D_func-S-1-all', 'D_func-S-2-all', 'D_func-S-3-all',
     'D_func-T-0-all', 'D_func-T-1-all', 'D_func-T-2-all', 'D_func-T-3-all',
     'D_func-Z-0-all', 'D_func-Z-1-all', 'D_func-Z-2-all', 'D_func-Z-3-all',
     'D_func-chi-0-all', 'D_func-chi-1-all', 'D_func-chi-2-all',
     'D_func-chi-3-all', 'D_lc-I-0-all', 'D_lc-I-1-all', 'D_lc-I-2-all',
     'D_lc-I-3-all', 'D_lc-S-0-all', 'D_lc-S-1-all', 'D_lc-S-2-all',
     'D_lc-S-3-all', 'D_lc-T-0-all', 'D_lc-T-1-all', 'D_lc-T-2-all',
     'D_lc-T-3-all', 'D_lc-Z-0-all', 'D_lc-Z-1-all', 'D_lc-Z-2-all',
     'D_lc-Z-3-all', 'D_lc-chi-0-all', 'D_lc-chi-1-all', 'D_lc-chi-2-all',
     'D_lc-chi-3-all', 'D_mc-I-0-all', 'D_mc-I-1-all', 'D_mc-I-2-all',
     'D_mc-I-3-all', 'D_mc-S-0-all', 'D_mc-S-1-all', 'D_mc-S-2-all',
     'D_mc-S-3-all', 'D_mc-T-0-all', 'D_mc-T-1-all', 'D_mc-T-2-all',
     'D_mc-T-3-all', 'D_mc-Z-0-all', 'D_mc-Z-1-all', 'D_mc-Z-2-all',
     'D_mc-Z-3-all', 'D_mc-chi-0-all', 'D_mc-chi-1-all', 'D_mc-chi-2-all',
     'D_mc-chi-3-all', 'f-I-0-all', 'f-I-1-all', 'f-I-2-all', 'f-I-3-all',
     'f-S-0-all', 'f-S-1-all', 'f-S-2-all', 'f-S-3-all', 'f-T-0-all', 'f-T-1-all',
     'f-T-2-all', 'f-T-3-all', 'f-Z-0-all', 'f-Z-1-all', 'f-Z-2-all', 'f-Z-3-all',
     'f-chi-0-all', 'f-chi-1-all', 'f-chi-2-all', 'f-chi-3-all', 'f-lig-I-0',
     'f-lig-I-1', 'f-lig-I-2', 'f-lig-I-3', 'f-lig-S-0', 'f-lig-S-1', 'f-lig-S-2',
     'f-lig-S-3', 'f-lig-T-0', 'f-lig-T-1', 'f-lig-T-2', 'f-lig-T-3', 'f-lig-Z-0',
     'f-lig-Z-1', 'f-lig-Z-2', 'f-lig-Z-3', 'f-lig-chi-0', 'f-lig-chi-1',
     'f-lig-chi-2', 'f-lig-chi-3', 'func-I-0-all', 'func-I-1-all',
     'func-I-2-all', 'func-I-3-all', 'func-S-0-all', 'func-S-1-all',
     'func-S-2-all', 'func-S-3-all', 'func-T-0-all', 'func-T-1-all',
     'func-T-2-all', 'func-T-3-all', 'func-Z-0-all', 'func-Z-1-all',
     'func-Z-2-all', 'func-Z-3-all', 'func-chi-0-all', 'func-chi-1-all',
     'func-chi-2-all', 'func-chi-3-all', 'lc-I-0-all', 'lc-I-1-all', 'lc-I-2-all',
     'lc-I-3-all', 'lc-S-0-all', 'lc-S-1-all', 'lc-S-2-all', 'lc-S-3-all',
     'lc-T-0-all', 'lc-T-1-all', 'lc-T-2-all', 'lc-T-3-all', 'lc-Z-0-all',
     'lc-Z-1-all', 'lc-Z-2-all', 'lc-Z-3-all', 'lc-chi-0-all', 'lc-chi-1-all',
     'lc-chi-2-all', 'lc-chi-3-all', 'mc-I-0-all', 'mc-I-1-all', 'mc-I-2-all',
     'mc-I-3-all', 'mc-S-0-all', 'mc-S-1-all', 'mc-S-2-all', 'mc-S-3-all',
     'mc-T-0-all', 'mc-T-1-all', 'mc-T-2-all', 'mc-T-3-all', 'mc-Z-0-all',
     'mc-Z-1-all', 'mc-Z-2-all', 'mc-Z-3-all', 'mc-chi-0-all', 'mc-chi-1-all',
     'mc-chi-2-all', 'mc-chi-3-all']
    geo = ['Df','Di', 'Dif','GPOAV','GPONAV','GPOV','GSA','POAV','POAV_vol_frac',
      'PONAV','PONAV_vol_frac','VPOV','VSA','rho']
     
    other = ['cif_file','name','filename']

    ANN_path = path + 'model/thermal/ANN/'
    temp_file_path = path + 'temp_file_creation_' + user_id + '/'
    df_train_all = pd.read_csv(ANN_path+"train.csv").append(pd.read_csv(ANN_path+"val.csv"))
    df_train = pd.read_csv(ANN_path+"train.csv")
    df_train = df_train.loc[:, (df_train != df_train.iloc[0]).any()]
    df_newMOF = pd.read_csv(temp_file_path + 'merged_descriptors/' + MOF_name + '_descriptors.csv') # Assume temp_file_creation/ in parent directory
    features = [val for val in df_train.columns.values if val in RACs+geo]

    X_train, X_newMOF, y_train, x_scaler, y_scaler = normalize_data_thermal(df_train, df_newMOF, features, ["T"], unit_trans=1, debug=False)
    X_train.shape, y_train.reshape(-1, ).shape 

    model = thermal_ANN

    from tensorflow.python.keras.backend import set_session
    with tf_session.as_default():
        with tf_session.graph.as_default():
            new_MOF_pred = y_scaler.inverse_transform(model.predict(X_newMOF))
            new_MOF_pred = np.round(new_MOF_pred,1) # round to 1 decimal

            # isolating just the prediction, since the model spits out the prediction like [[PREDICTION]], as in, in hard brackets
            new_MOF_pred = new_MOF_pred[0][0]
            new_MOF_pred = str(new_MOF_pred)

            # adding units
            degree_sign= u'\N{DEGREE SIGN}'
            new_MOF_pred = new_MOF_pred + degree_sign + 'C' # degrees Celsius

            # Define the function for the latent space. This will depend on the model. We want the layer before the last, in this case this was the 8th one.
            get_latent = K.function([model.layers[0].input],
                                    [model.layers[8].output]) # Last layer before dense-last

            # Get the latent vectors for the training data first, then the latent vectors for the test data.
            training_latent = get_latent([X_train, 0])[0]
            design_latent = get_latent([X_newMOF, 0])[0]

            print(training_latent.shape,design_latent.shape)

    # Compute the pairwise distances between the test latent vectors and the train latent vectors to get latent distances
    d1 = pairwise_distances(design_latent,training_latent,n_jobs=30)
    df1 = pd.DataFrame(data=d1, columns=df_train['CoRE_name'].tolist())
    df1.to_csv(temp_file_path + 'solvent_test_latent_dists.csv')

    # Want to find the closest points (let's say the closest 5 points); so, smallest values in df1
    neighbors = 5 # number of closest points

    # will make arrays of length neighbors, where each entry is the next closest neighbor (will do this for both names and distances)
    neighbors_names = []
    neighbors_distances = []

    df_reformat = df1.min(axis='index')

    for i in range(neighbors):
        name = df_reformat.idxmin() # name of next closest complex in the traiing data
        distance = df_reformat.min() # distance of the next closest complex in the training data to the new MOF
        df_reformat = df_reformat.drop(name) # dropping the next closest complex, in order to find the next-next closest complex

        neighbors_names.append(name)
        neighbors_distances.append(str(distance))

    return new_MOF_pred, neighbors_names, neighbors_distances

# This function is used by both ss_predict() and ts_predict() to generate RACs and Zeo++ descriptors
# These descriptors are subsequently used in ss_predict() and ts_predict() for the ANN models
# Inputs are the name of the MOF and the structure (cif file text) of the MOF for which descriptors are to be generated
# The third input indicates the type of prediction (solvent removal or thermal)
def descriptor_generator(name, structure, prediction_type):
    print('TIME CHECK 2')
    import time # debugging
    timeStarted = time.time() # save start time (debugging)

    temp_file_folder = MOFSIMPLIFY_PATH + "temp_file_creation_" + str(session['ID']) + '/'
    cif_folder = temp_file_folder + 'cifs/'

    # Write the data back to a cif file.
    try:
        cif_file = open(cif_folder + name + '.cif', 'w')
    except FileNotFoundError:
        return 'FAILED'
    cif_file.write(structure)
    cif_file.close()

    RACs_folder = temp_file_folder +  prediction_type + '_RACs/'
    zeo_folder = temp_file_folder + prediction_type + '_zeo++/'

    # Delete the RACs folder, then remake it (to start fresh for this prediction).
    shutil.rmtree(RACs_folder)
    os.mkdir(RACs_folder)

    # Doing the same with the Zeo++ folder.
    shutil.rmtree(zeo_folder)
    os.mkdir(zeo_folder)

    # Next, running MOF featurization
    try:
        get_primitive(cif_folder + name + '.cif', cif_folder + name + '_primitive.cif');
    except ValueError:
        return 'FAILED'

    # try:
    #     full_names, full_descriptors = get_MOF_descriptors(cif_folder + name + '_primitive.cif',3,path= RACs_folder, xyzpath= RACs_folder + name + '.xyz');
    #         # makes the linkers and sbus folders
    # except ValueError:
    #     return 'FAILED'
    # except NotImplementedError:
    #     return 'FAILED'
    # except AssertionError:
    #     return 'FAILED'

    # if (len(full_names) <= 1) and (len(full_descriptors) <= 1): # this is a featurization check from MOF_descriptors.py
    #     return 'FAILED'

    timeDelta = time.time() - timeStarted # get execution time
    print('Finished process in ' + str(timeDelta) + ' seconds')

    print('TIME CHECK 3')

    # At this point, have the RAC featurization. Need geometry information next.

    # Run Zeo++

    timeStarted = time.time() # save start time (debugging)

    cmd1 = MOFSIMPLIFY_PATH + 'zeo++-0.3/network -ha -res ' + zeo_folder + name + '_pd.txt ' + cif_folder + name + '_primitive.cif'
    cmd2 = MOFSIMPLIFY_PATH + 'zeo++-0.3/network -sa 1.86 1.86 10000 ' + zeo_folder + name + '_sa.txt ' + cif_folder + name + '_primitive.cif'
    cmd3 = MOFSIMPLIFY_PATH + 'zeo++-0.3/network -volpo 1.86 1.86 10000 ' + zeo_folder + name + '_pov.txt '+ cif_folder + name + '_primitive.cif'
    cmd4 = 'python ' + MOFSIMPLIFY_PATH + 'model/RAC_getter.py %s %s %s' %(cif_folder, name, RACs_folder)

    # four parallelized Zeo++ commands
    process1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=None, shell=True)
    process2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, stderr=None, shell=True)
    process3 = subprocess.Popen(cmd3, stdout=subprocess.PIPE, stderr=None, shell=True)
    process4 = subprocess.Popen(cmd4, stdout=subprocess.PIPE, stderr=None, shell=True)

    output1 = process1.communicate()[0]
    output2 = process2.communicate()[0]
    output3 = process3.communicate()[0]
    output4 = process4.communicate()[0]

    # Have written output of Zeo++ commands to files. Now, code below extracts information from those files.

    ''' The geometric descriptors are largest included sphere (Di), 
    largest free sphere (Df), largest included sphere along free path (Dif),
    crystal density (rho), volumetric surface area (VSA), gravimetric surface (GSA), 
    volumetric pore volume (VPOV) and gravimetric pore volume (GPOV). 
    Also, we include cell volume as a descriptor.

    All Zeo++ calculations use a pore radius of 1.86 angstrom, and zeo++ is called by subprocess.
    '''

    dict_list = []
    cif_file = name + '_primitive.cif' 
    basename = cif_file.strip('.cif')
    largest_included_sphere, largest_free_sphere, largest_included_sphere_along_free_sphere_path  = np.nan, np.nan, np.nan
    unit_cell_volume, crystal_density, VSA, GSA  = np.nan, np.nan, np.nan, np.nan
    VPOV, GPOV = np.nan, np.nan
    POAV, PONAV, GPOAV, GPONAV, POAV_volume_fraction, PONAV_volume_fraction = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

    if (os.path.exists(zeo_folder + name + '_pd.txt') & os.path.exists(zeo_folder + name + '_sa.txt') &
        os.path.exists(zeo_folder + name + '_pov.txt')):
        with open(zeo_folder + name + '_pd.txt') as f:
            pore_diameter_data = f.readlines()
            for row in pore_diameter_data:
                largest_included_sphere = float(row.split()[1]) # largest included sphere
                largest_free_sphere = float(row.split()[2]) # largest free sphere
                largest_included_sphere_along_free_sphere_path = float(row.split()[3]) # largest included sphere along free sphere path
        with open(zeo_folder + name + '_sa.txt') as f:
            surface_area_data = f.readlines()
            for i, row in enumerate(surface_area_data):
                if i == 0:
                    unit_cell_volume = float(row.split('Unitcell_volume:')[1].split()[0]) # unit cell volume
                    crystal_density = float(row.split('Unitcell_volume:')[1].split()[0]) # crystal density
                    VSA = float(row.split('ASA_m^2/cm^3:')[1].split()[0]) # volumetric surface area
                    GSA = float(row.split('ASA_m^2/g:')[1].split()[0]) # gravimetric surface area
        with open(zeo_folder + name + '_pov.txt') as f:
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
        print('Not all 3 files exist, so at least one Zeo++ call failed!', 'sa: ',os.path.exists(zeo_folder + name + '_sa.txt'), 
              '; pd: ',os.path.exists(zeo_folder + name + '_pd.txt'), '; pov: ', os.path.exists(zeo_folder + name + '_pov.txt'))
        return 'FAILED'
    geo_dict = {'name':basename, 'cif_file':cif_file, 'Di':largest_included_sphere, 'Df': largest_free_sphere, 'Dif': largest_included_sphere_along_free_sphere_path,
                'rho': crystal_density, 'VSA':VSA, 'GSA': GSA, 'VPOV': VPOV, 'GPOV':GPOV, 'POAV_vol_frac':POAV_volume_fraction, 
                'PONAV_vol_frac':PONAV_volume_fraction, 'GPOAV':GPOAV,'GPONAV':GPONAV,'POAV':POAV,'PONAV':PONAV}
    dict_list.append(geo_dict)
    geo_df = pd.DataFrame(dict_list)
    geo_df.to_csv(zeo_folder + 'geometric_parameters.csv',index=False)

    # error handling for cmd4
    with open(RACs_folder + 'RAC_getter_log.txt', 'r') as f:
        if f.readline() == 'FAILED':
            print('RAC generation failed.')
            return 'FAILED'


    timeDelta = time.time() - timeStarted # get execution time
    print('Finished process in ' + str(timeDelta) + ' seconds')

    print('TIME CHECK 4')

    timeStarted = time.time() # save start time

    # Merging geometric information with get_MOF_descriptors files (lc_descriptors.csv, sbu_descriptors.csv, linker_descriptors.csv)
    try:
        lc_df = pd.read_csv(RACs_folder + "lc_descriptors.csv") 
        sbu_df = pd.read_csv(RACs_folder + "sbu_descriptors.csv")
        linker_df = pd.read_csv(RACs_folder + "linker_descriptors.csv")
    except Exception: # csv files have been deleted
        return 'FAILED' 

    lc_df = lc_df.mean().to_frame().transpose() # averaging over all rows. Convert resulting Series into a Dataframe, then transpose
    sbu_df = sbu_df.mean().to_frame().transpose()
    linker_df = linker_df.mean().to_frame().transpose()

    merged_df = pd.concat([geo_df, lc_df, sbu_df, linker_df], axis=1)

    merged_df.to_csv(temp_file_folder + '/merged_descriptors/' + name + '_descriptors.csv',index=False) # written in /temp_file_creation_SESSIONID

    if prediction_type == 'solvent':

        ANN_folder = MOFSIMPLIFY_PATH + 'model/solvent/ANN/'
        train_df = pd.read_csv(ANN_folder + 'dropped_connectivity_dupes/train.csv')

    if prediction_type == 'thermal':
        ANN_folder = MOFSIMPLIFY_PATH + 'model/thermal/ANN/'
        train_df = pd.read_csv(ANN_folder + 'train.csv')

    # Here, I do a check to see if the current MOF is in the training data.
    # If it is, then I return the known truth for the MOF, rather than make a prediction.

    # Will iterate through the rows of the train pandas dataframe


    in_train = False

    for index, row in train_df.iterrows(): # iterate through rows of the training data MOFs

        row_match = True # gets set to false if any values don't match 
        matching_MOF = None # gets assigned a value if the current MOF is in the training data

        for col in merged_df.columns: # iterate through columns of the single new MOF we are predicting on (merged_df is just one row)
            if col == 'name' or col == 'cif_file' or col == 'Dif':
                continue # skip these
                # Dif was sometimes differing between new Zeo++ call and training data value, for the same MOF
            # print('Column: ')
            # print(col)
            # print(row[col])
            # print(merged_df.iloc[0][col])

            # If for any property a training MOF and the new MOF we are predicting on differ too much, we know they are not the same MOF
            # So row_match is set to false for this training MOF
            if np.absolute(row[col] - merged_df.iloc[0][col]) > 0.05 * np.absolute(merged_df.iloc[0][col]): # row[col] != merged_df.iloc[0][col] was leading to some same values being idenfitied as different b/c of some floating 10^-15 values 
                row_match = False

                # print(row['CoRE_name'])
                # # debugging
                # if row['CoRE_name'] == 'HISJAW_clean': (check with AFUKIX for thermal)
                #     print('HISJAW check')
                #     print(col)
                #     print(row[col])
                #     print(merged_df.iloc[0][col])

                break
        
        if row_match and prediction_type == 'solvent': # all columns for the row match! Some training MOF is the same as the new MOF
            in_train = True
            match_truth = row['flag'] # the flag for the MOF that matches the current MOF
            matching_MOF = row['CoRE_name']
            print(row['CoRE_name'])
            break
        if row_match and prediction_type =='thermal': # all columns for the row match! Some training MOF is the same as the new MOF
            in_train = True
            match_truth = row['T'] # the flag for the MOF that matches the current MOF
            match_truth = np.round(match_truth,1) # round to 1 decimal
            matching_MOF = row['CoRE_name']
            # adding units
            degree_sign= u'\N{DEGREE SIGN}'
            match_truth = str(match_truth) + degree_sign + 'C' # degrees Celsius
            break

    if in_train:
        myDict = {'in_train': True, 'truth': match_truth, 'match': matching_MOF}
        return myDict

    timeDelta = time.time() - timeStarted # get execution time
    print('Finished process in ' + str(timeDelta) + ' seconds')

    print('TIME CHECK 5')

    return [temp_file_folder, ANN_folder]


    
# Note: the h5 model for the solvent stability prediction and the thermal stability prediction should be trained on the same version of TensorFlow (here, 1.14)
# the two h5 models show up in solvent_ANN.py and thermal_ANN.py, respectively
@app.route('/predict_solvent_stability', methods=['POST']) 
def ss_predict():
    # Generates solvent stability prediction.
    # To do this, need to generate RAC featurization and Zeo++ geometry information for the MOF.
    # Then, apply Aditya's model to make prediction.

    print('TIME CHECK 1')
    import time
    timeStarted = time.time() # save start time (debugging)

    # Grab data
    my_data = json.loads(flask.request.get_data())
    structure = my_data['structure']
    name = my_data['name']
    if name == 'Example MOF':
        name = 'HKUST-1' # spacing in name was causing issues down the line
    if name[-4:] == '.cif':
        name = name[:-4]

    timeDelta = time.time() - timeStarted # get execution time
    print('Finished process in ' + str(timeDelta) + ' seconds')

    output = descriptor_generator(name, structure, 'solvent') # generate descriptors
    if output == 'FAILED': # model failure
        return 'FAILED'
    elif isinstance(output, dict): # MOF was in the training data
        return output
    else: # grabbing some variables. We will make a prediction
        temp_file_folder = output[0]
        ANN_folder = output[1]

    # Applying the model next

    timeStarted = time.time() # save start time (debugging)
    prediction, neighbor_names, neighbor_distances = run_solvent_ANN(str(session['ID']), MOFSIMPLIFY_PATH, name, solvent_model)

    # print('check check')
    # print(neighbor_names) # debugging
    # print(neighbor_distances) # debugging
    # print(type(neighbor_names)) # debugging
    # print(type(neighbor_distances)) # debugging
    # print(prediction) # debugging
    # print(type(prediction))  # debugging  

    results = {'prediction': prediction,
        'neighbor_names': neighbor_names,
        'neighbor_distances': neighbor_distances,
        'in_train': False} # a prediction was made. Requested MOF was not in the training data.

    # print(results) # debugging
    # print(type(results)) # debugging

    timeDelta = time.time() - timeStarted # get execution time
    print('Finished process in ' + str(timeDelta) + ' seconds') 

    print('TIME CHECK 6')

    return results

@app.route('/predict_thermal_stability', methods=['POST']) 
def ts_predict():
    # Generates thermal stability prediction.
    # To do this, need to generate RAC featurization and Zeo++ geometry information for the MOF.
    # Then, apply Aditya's model to make prediction.


    # Grab data
    my_data = json.loads(flask.request.get_data())
    structure = my_data['structure']
    name = my_data['name']
    if name == 'Example MOF':
        name = 'HKUST-1' # spacing in name was causing issues down the line
    if name[-4:] == '.cif':
        name = name[:-4]

    output = descriptor_generator(name, structure, 'thermal') # generate descriptors
    if output == 'FAILED': # model failure
        return 'FAILED'
    elif isinstance(output, dict): # MOF was in the training data
        return output
    else: # grabbing some variables. We will make a prediction
        temp_file_folder = output[0]
        ANN_folder = output[1]

    # Applying the model next

    timeStarted = time.time() # save start time (debugging)

    prediction, neighbor_names, neighbor_distances = run_thermal_ANN(str(session['ID']), MOFSIMPLIFY_PATH, name, thermal_model)

    timeDelta = time.time() - timeStarted # get execution time
    print('Finished process in ' + str(timeDelta) + ' seconds')

    print('check check')
    print(neighbor_names) # debugging
    print(neighbor_distances) # debugging
    print(type(neighbor_names)) # debugging
    print(type(neighbor_distances)) # debugging    

    results = {'prediction': prediction,
        'neighbor_names': neighbor_names,
        'neighbor_distances': neighbor_distances,
        'in_train': False} # a prediction was made. Requested MOF was not in the training data.

    print('TIME CHECK 6 ' + str(session['ID']))

    return results

@app.route('/plot_thermal_stability', methods=['POST']) 
def plot_thermal_stability():
    # Returns a plot of the distribution of thermal breakdown temperatures of the MOFs our ANN was trained on.
    # Additionally, displays the position of the current MOF's thermal breakdown temperature.

    # Grab data
    info = json.loads(flask.request.get_data()) 
    my_data = info['temperature'] # this is the current MOF's predicted thermal breakdown temperature
    my_data = my_data[:-3] # getting rid of the celsius symbol, left with just the number
    my_data = float(my_data)
    # print('checkerino')
    # print(my_data)

    # Getting the temperature data
    temps_df = pd.read_csv(MOFSIMPLIFY_PATH + "model/thermal/ANN/adjusted_TSD_df_all.csv")

    import matplotlib
    matplotlib.use('Agg') # noninteractive backend
    import matplotlib.pyplot as plt
    plt.close("all")
    plt.rcParams.update({'font.size': 16}) # large font
    import scipy.stats as stats

    # In training data, smallest T breakdown is 35, and largest T breakdown is 654.

    # use stats.gaussian_kde to estimate the probability density function from the histogram
    density = stats.gaussian_kde(temps_df['T'])
    x = np.arange(30,661,1) # in training data, smallest T breakdown is 35, and largest T breakdown is 654
    fig = plt.figure(figsize=(10, 6), dpi=80) # set figure size
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

    # Grab data.
    my_data = json.loads(flask.request.get_data()) # this is the current MOF's predicted thermal breakdown temperature
    my_data = my_data[:-3] # getting rid of the celsius symbol, left with just the number
    my_data = float(my_data)
    print(my_data)

    # Getting the temperature data.
    temps_df = pd.read_csv(MOFSIMPLIFY_PATH + "model/thermal/ANN/adjusted_TSD_df_all.csv")

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

    # Grab data
    my_data = json.loads(flask.request.get_data()); # This is the neighbor complex

    # cut off these endings, in order to access the TGA file correctly:
    # _clean, _ion_b, _neutral_b, _SL, _charged, _clean_h, _manual, _auto, _charged, etc
    cut_index = my_data.find('_') # gets index of the first underscore
    my_data = my_data[:cut_index] 

    # Grab data 
    slopes_df = pd.read_csv(MOFSIMPLIFY_PATH + "TGA/raw_TGA_digitization_data/digitized_csv/" + my_data + ".csv")

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

    # Want lines to extend beyond intersection point. I will make it 20%
    x_extension = [None] * 2 # will hold the x values of the two extension points
    y_extension = [None] * 2 # will hold the y values of the two extension points
    slope = [None] * 2 # will hold the slopes of the two lines

    x_extension[0] = intersection_point[0] + (intersection_point[0] - x_values[0])*0.2 # x_values[0] is the smallest value, x_values[3] is the largest value
    x_extension[1] = intersection_point[0] - (x_values[3] - intersection_point[0])*0.2
    slope[0] = (y_values[1] - y_values[0])/(x_values[1] - x_values[0])
    slope[1] = (y_values[3] - y_values[2])/(x_values[3] - x_values[2])
    y_extension[0] = intersection_point[1] + slope[0] * (x_extension[0] - intersection_point[0]) # y2 = y1 + m(x2-x1)
    y_extension[1] = intersection_point[1] + slope[1] * (x_extension[1] - intersection_point[0]) # y2 = y1 + m(x2-x1)

    # Instantiating the figure object. 
    graph = figure(title = "Simplified literature TGA plot of selected thermal ANN neighbor")  
         
    # The points to be plotted.
    xs = [[x_values[0], x_values[1],intersection_point[0], x_extension[0]], [x_values[2], x_values[3],intersection_point[0], x_extension[1]]] 
    ys = [[y_values[0], y_values[1],intersection_point[1], y_extension[0]], [y_values[2], y_values[3],intersection_point[1], y_extension[1]]] 
        
    # Plotting the graph.
    graph.multi_line(xs, ys, line_dash='dashed') 
    graph.circle([intersection_point[0]], [intersection_point[1]], size=20, fill_color="red", line_color='black')
    graph.star(x=x_values, y=y_values, size=20, fill_color="yellow", line_color='black') # four stars
    graph.xaxis.axis_label = 'Temperature (°C)'
    graph.yaxis.axis_label = 'Percentage mass remaining or Mass'    

    return file_html(graph,CDN,'my plot')

@app.route('/get_components', methods=['POST']) 
def get_components():
    # Uses Aditya's MOF code to get linkers and sbus
    # Returns a dictionary with the linker and sbu xyz files's text, along with information about the number of linkers and sbus
    # Also in the dictionary: SMILES string for each of the linkers and sbus

    # Grab data
    my_data = json.loads(flask.request.get_data());

    structure = my_data['structure']
    name = my_data['name']
    if name == 'Example MOF':
        name = 'HKUST-1' # spacing in name was causing issues down the line
    if name[-4:] == '.cif':
        name = name[:-4]

    temp_file_folder = MOFSIMPLIFY_PATH + "temp_file_creation_" + str(session['ID']) + '/'
    cif_folder = temp_file_folder + 'cifs/'

    # Write the data back to a cif file.
    try:
        cif_file = open(cif_folder + name + '.cif', 'w')
    except FileNotFoundError:
        return 'FAILED'
    cif_file.write(structure)
    cif_file.close()

    RACs_folder = temp_file_folder +  'components_RACs/'

    # Delete the RACs folder, then remake it (to start fresh for this prediction).
    shutil.rmtree(RACs_folder)
    os.mkdir(RACs_folder)

    # Next, running MOF featurization
    try:
        get_primitive(cif_folder + name + '.cif', cif_folder + name + '_primitive.cif');
    except ValueError:
        return 'FAILED'

    try:
        full_names, full_descriptors = get_MOF_descriptors(cif_folder + name + '_primitive.cif',3,path= RACs_folder, xyzpath= RACs_folder + name + '.xyz');
            # makes the linkers and sbus folders
    except ValueError:
        return 'FAILED'
    except NotImplementedError:
        return 'FAILED'
    except AssertionError:
        return 'FAILED'

    if (len(full_names) <= 1) and (len(full_descriptors) <= 1): # this is a featurization check from MOF_descriptors.py
        return 'FAILED'

    # At this point, have the RAC featurization. 

    # will return a json object (dictionary)
    # The fields are string representations of the linkers and sbus, however many there are.

    dictionary = {};
    
    linker_num = 0;
    while True:
        if not os.path.exists(RACs_folder + 'linkers/' + name + '_primitive_linker_' + str(linker_num) + '.xyz'):
            break
        else:
            linker_file = open(RACs_folder + 'linkers/' + name + '_primitive_linker_' + str(linker_num) + '.xyz', 'r');
            linker_info = linker_file.read();
            linker_file.close();

            dictionary['linker_' + str(linker_num)] = linker_info;

            linker_num = linker_num + 1;


    sbu_num = 0;
    while True:
        if not os.path.exists(RACs_folder + 'sbus/' + name + '_primitive_sbu_' + str(sbu_num) + '.xyz'):
            break
        else:
            sbu_file = open(RACs_folder + 'sbus/' + name + '_primitive_sbu_' + str(sbu_num) + '.xyz', 'r');
            sbu_info = sbu_file.read();
            sbu_file.close();

            dictionary['sbu_' + str(sbu_num)] = sbu_info;

            sbu_num = sbu_num + 1;


    dictionary['total_linkers'] = linker_num;
    dictionary['total_sbus'] = sbu_num;


    # Identifying which linkers and sbus have different connectivities.

    # Code below uses molecular graph determinants.
    # two ligands that are the same (via connectivity) will have the same molecular graph determinant. 
    # Molecular graph determinant fails to distinguish between isomers, but so would RCM (Reverse Cuthill McKee).

    # This simple script is meant to be run within the linkers directory, and it will give a bunch of numbers. 
    # If those numbers are the same, the linker is the same, if not, the linkers are different, etc.

    from molSimplify.Classes.mol3D import mol3D
    import glob
    MOF_of_interest = name + '_primitive'
    XYZs = sorted(glob.glob(RACs_folder + 'linkers/*'+MOF_of_interest+'*xyz'))
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

    XYZs = sorted(glob.glob(RACs_folder + 'sbus/*'+MOF_of_interest+'*xyz'))
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

    linkers_folder = RACs_folder + 'linkers/'

    for i in range(dictionary['total_linkers']): # 0, 1, 2, ..., numberoflinkersminus1
        smilesFile = pybel.Outputfile('smi', linkers_folder + name + '_primitive_linker_' + str(i) + '.txt') # smi refers to SMILES
        smilesFile.write(next(pybel.readfile('xyz', linkers_folder + name + '_primitive_linker_' + str(i) + '.xyz'))) # writes SMILES string to the text file

        # Next, get the SMILES string from the text file.
        f = open(linkers_folder + name + '_primitive_linker_' + str(i) + '.txt', 'r')
        line = f.readline()
        line = line.split('\t') # split at tabs
        smiles_ID = line[0]
        dictionary['linker_' + str(i) + '_SMILES'] = smiles_ID
        f.close()

    sbus_folder = RACs_folder + 'sbus/'

    for i in range(dictionary['total_sbus']):
        smilesFile = pybel.Outputfile('smi', sbus_folder + name + '_primitive_sbu_' + str(i) + '.txt') # smi refers to SMILES
        smilesFile.write(next(pybel.readfile('xyz', sbus_folder + name + '_primitive_sbu_' + str(i) + '.xyz'))) # writes SMILES string to the text file

        # Next, get the SMILES string from the text file.
        f = open(sbus_folder + name + '_primitive_sbu_' + str(i) + '.txt', 'r')
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
    
    my_folder = MOFSIMPLIFY_PATH + 'model/solvent/ANN/dropped_connectivity_dupes/'
    solvent_flags_df = pd.read_csv(my_folder + 'train.csv')

    this_neighbor = solvent_flags_df[solvent_flags_df['CoRE_name'] == my_data] # getting the row with the MOF of interest

    this_neighbor_flag = this_neighbor['flag'] # getting the flag value
    this_neighbor_flag =  str(this_neighbor_flag.iloc[0]) # extract the flag value and return it as a string

    # next, getting DOI
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
    ANN_folder = MOFSIMPLIFY_PATH + 'model/thermal/ANN/'

    breakdown_T_df = pd.read_csv(ANN_folder + 'train.csv')
    this_neighbor = breakdown_T_df[breakdown_T_df['CoRE_name'] == my_data] # getting the row with the MOF of interest
    this_neighbor_T = this_neighbor['T'] # getting the breakdown temperature value
    this_neighbor_T =  str(round(this_neighbor_T.iloc[0], 1)) # extract the breakdown temperature value and return it as a string. Want just one decimal place


    TGA_folder = MOFSIMPLIFY_PATH + 'TGA/'
    TGA_df = pd.read_excel(TGA_folder + 'TGA_info_log.xlsx')
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

    print('neighbor writer check')
    print(my_data)

    prediction_type = my_data['prediction_type']
    current_MOF = my_data['current_MOF']
    selected_neighbor = my_data['selected_neighbor']
    latent_space_distance = my_data['latent_space_distance']

    # next, getting neighbor_truth
    if prediction_type == 'solvent':
        my_folder = MOFSIMPLIFY_PATH + 'model/solvent/ANN/dropped_connectivity_dupes/'
        solvent_flags_df = pd.read_csv(my_folder + 'train.csv')
        this_neighbor = solvent_flags_df[solvent_flags_df['CoRE_name'] == selected_neighbor] # getting the row with the MOF of interest
        this_neighbor_flag = this_neighbor['flag'] # getting the flag value
        neighbor_truth = str(this_neighbor_flag.iloc[0]) # extract the flag value and convert it into a string
        if neighbor_truth == '1':
            neighbor_truth = 'Stable upon solvent removal'
        else:
            neighbor_truth = 'Unstable upon solvent removal'

    else: # thermal stability
        my_folder = MOFSIMPLIFY_PATH + 'model/thermal/ANN/'
        breakdown_T_df = pd.read_csv(my_folder + 'train.csv')
        this_neighbor = breakdown_T_df[breakdown_T_df['CoRE_name'] == selected_neighbor] # getting the row with the MOF of interest
        this_neighbor_T = this_neighbor['T'] # getting the breakdown temperature value
        neighbor_truth = str(round(this_neighbor_T.iloc[0], 1)) # extract the flag value and convert it into a string. Want just one decimal place


    # next, getting DOI
    if prediction_type == 'solvent':
        # os.chdir('model/solvent/ANN/dropped_connectivity_dupes')
        # solvent_flags_df = pd.read_csv('train.csv')
        this_neighbor = solvent_flags_df[solvent_flags_df['CoRE_name'] == selected_neighbor] # getting the row with the MOF of interest
        this_neighbor_doi = this_neighbor['doi'] # getting the doi
        this_neighbor_doi = this_neighbor_doi.iloc[0]
    else: # thermal
        TGA_df = pd.read_excel(MOFSIMPLIFY_PATH + 'TGA/TGA_info_log.xlsx')
        this_neighbor = TGA_df[TGA_df['CoRE_name'] == selected_neighbor] # getting the row with the MOF of interest
        this_neighbor_doi = this_neighbor['doi'] # getting the doi
        this_neighbor_doi = this_neighbor_doi.iloc[0]

    # Now, writing all of this information to a string

    # Delete and remake the latent_neighbor folder each time, so only one file is ever in it
    shutil.rmtree(MOFSIMPLIFY_PATH + 'temp_file_creation_' + str(session['ID']) + '/latent_neighbor')
    os.mkdir(MOFSIMPLIFY_PATH + 'temp_file_creation_' + str(session['ID']) + '/latent_neighbor')

    with open(MOFSIMPLIFY_PATH + 'temp_file_creation_' + str(session['ID']) + '/latent_neighbor/' + prediction_type + '-' + current_MOF + '-' + selected_neighbor + '.txt','w') as f:
        f.write('Prediction type: ' + prediction_type + '\n')
        f.write('Current MOF: ' + current_MOF + '\n')
        f.write('Selected CoRE nearest neighbor in latent space: ' + selected_neighbor + '\n')
        f.write('Latent space distance: ' + latent_space_distance + '\n')
        if prediction_type == 'solvent':
            f.write('Property for nearest neighbor: ' + neighbor_truth + '\n')
        else:
            degree_sign= u'\N{DEGREE SIGN}'
            f.write('Property for nearest neighbor: ' + neighbor_truth + ' ' + degree_sign + 'C\n')
        f.write('Neighbor DOI: ' + this_neighbor_doi + '\n')

    with open(MOFSIMPLIFY_PATH + 'temp_file_creation_' + str(session['ID']) + '/latent_neighbor/' + prediction_type + '-' + current_MOF + '-' + selected_neighbor + '.txt','r') as f:
        contents = f.read()

    print('neighbor writer check 2')

    my_dict = {'contents': contents, 'file_name': prediction_type + '-' + current_MOF + '-' + selected_neighbor + '.txt'}

    return my_dict


# @app.route('/TGA_maker', methods=['POST']) 
# def TGA_maker():
#     # Making the TGA plot and saving it, for it to be downloaded
#     from bokeh.io import export_png

#     # Grab data
#     my_data = json.loads(flask.request.get_data()); # This is the neighbor complex
#     my_data = my_data[:6] # only want the first six letters

#     print('neighbor writer check 3')

#     # Grab data 
#     slopes_df = pd.read_csv(MOFSIMPLIFY_PATH + "TGA/raw_TGA_digitization_data/digitized_csv/" + my_data + ".csv")

#     x_values = []
#     y_values = []
#     for i in range(4): # 0, 1, 2, 3
#         x_values.append(slopes_df.iloc[[i]]['T (degrees C)'][i])
#         y_values.append(slopes_df.iloc[[i]]['mass (arbitrary units)'][i])

#     # Making the four points.
#     p1 = np.array( [x_values[0], y_values[0]] )
#     p2 = np.array( [x_values[1], y_values[1]] )

#     p3 = np.array( [x_values[2], y_values[2]] )
#     p4 = np.array( [x_values[3], y_values[3]] )

#     intersection_point = seg_intersect(p1, p2, p3, p4)

#     # Instantiating the figure object. 
#     graph = figure(title = "Simplified literature TGA plot of selected thermal ANN neighbor")  
         
#     # The points to be plotted.
#     xs = [[x_values[0], x_values[1],intersection_point[0]], [x_values[2], x_values[3],intersection_point[0]]] 
#     ys = [[y_values[0], y_values[1],intersection_point[1]], [y_values[2], y_values[3],intersection_point[1]]] 
        
#     # Plotting the graph.
#     graph.multi_line(xs, ys) 
#     graph.circle([intersection_point[0]], [intersection_point[1]], size=20, color="navy", alpha=0.5)
#     graph.xaxis.axis_label = 'Temperature (°C)'
#     graph.yaxis.axis_label = 'Percentage mass remaining or Mass' 

#     print('neighbor writer check 4')

#     export_png(graph, filename='temp_file_creation_' + str(session['ID']) + '/latent_neighbor/' + my_data + "_simplified_TGA.png")
    
#     return 'Success!'

@app.route('/get_descriptors', methods=['POST']) 
def descriptor_getter():
    # Grab data
    name = json.loads(flask.request.get_data()); # This is the selected MOF
    if name == 'Example MOF':
        name = 'HKUST-1' # spacing in name was causing issues down the line
    if name[-4:] == '.cif':
        name = name[:-4]

    temp_file_folder = MOFSIMPLIFY_PATH + "temp_file_creation_" + str(session['ID']) + '/'
    descriptors_folder = temp_file_folder + "merged_descriptors/"

    with open(descriptors_folder + name + '_descriptors.csv', 'r') as f:
        contents = f.read()

    return contents

@app.route('/get_TGA', methods=['POST']) 
def TGA_getter():
    # Grab data
    name = json.loads(flask.request.get_data()); # This is the thermal ANN latent space nearest neighbor
    
    # cut off these endings, in order to access the TGA file correctly:
    # _clean, _ion_b, _neutral_b, _SL, _charged, _clean_h, _manual, _auto, _charged, etc
    cut_index = name.find('_') # gets index of the first underscore
    name = name[:cut_index] 

    tga_folder = MOFSIMPLIFY_PATH + "TGA/raw_TGA_digitization_data/digitized_csv/"

    with open(tga_folder + name + '.csv', 'r') as f:
        contents = f.read()

    return contents

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8000, debug=True)
    #app.run(host='0.0.0.0', port=8000, threaded=False, processes=10)
