from numpy.random import seed
seed(1)
import pandas as pd
import numpy as np
import pickle
import sys, os
from keras import backend as K
import tensorflow as tf
from scipy.spatial import distance_matrix
from sklearn.metrics import pairwise_distances
import json
import sklearn
tf.compat.v1.disable_eager_execution()

# Load in whatever model is of interest here
from keras.models import load_model
model = load_model('final_model_T_few_epochs.h5')

# Standard function to normalize data
def normalize_data(df_train, df_test, fnames, lname, unit_trans=1, debug=False):
    _df_train = df_train.copy().dropna(subset=fnames+lname)
    _df_test = df_test.copy().dropna(subset=fnames+lname)
    X_train, X_test = _df_train[fnames].values, _df_test[fnames].values
    y_train, y_test = _df_train[lname].values, _df_test[lname].values
    if debug:
        print("training data reduced from %d -> %d because of nan." % (len(df_train), y_train.shape[0]))
        print("test data reduced from %d -> %d because of nan." % (len(df_test), y_test.shape[0]))
    x_scaler = sklearn.preprocessing.StandardScaler()
    x_scaler.fit(X_train)
    X_train = x_scaler.transform(X_train)
    X_test = x_scaler.transform(X_test)
    y_scaler = sklearn.preprocessing.StandardScaler()
    y_scaler.fit(y_train)
    y_train = y_scaler.transform(y_train)
    y_test = y_scaler.transform(y_test)
    return X_train, X_test, y_train, y_test, x_scaler, y_scaler

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

# Load in the data, keep relevant features
path = os.getcwd()+'/'
df_train_all = pd.read_csv(path+"/train.csv").append(pd.read_csv(path+"/val.csv"))
df_train = pd.read_csv(path+"/train.csv")
df_train = df_train.loc[:, (df_train != df_train.iloc[0]).any()]
df_val = pd.read_csv(path+"/val.csv")
df_test = pd.read_csv(path+"/test.csv")
features = [val for val in df_train.columns.values if val in RACs+geo]

# Scale the data
X_train, X_test, y_train, y_test, x_scaler, y_scaler = normalize_data(df_train, df_test, features, ["T"], unit_trans=1, debug=False)
_, X_val, _, y_val, _, _ = normalize_data(df_train, df_val, features, ["T"], unit_trans=1, debug=False)

# Define the function for the latent space. This will depend on the model. We want the layer before the last, in this case this was the 8th one.
get_latent = K.function([model.layers[0].input],
                        [model.layers[8].output])

# Get the latent vectors for the training data first, then the latent vectors for the test data.
training_latent = get_latent([X_train, 0])[0]
design_latent = get_latent([X_test, 0])[0]

print(training_latent.shape,design_latent.shape)

# Compute the pairwise distances between the test latent vectors and the train latent vectors to get latent distances
d1 = pairwise_distances(design_latent,training_latent,n_jobs=30)
df1 = pd.DataFrame(data=d1, index=df_test['CoRE_name'].tolist(), columns=df_train['CoRE_name'].tolist())
df1.to_csv('thermal_test_latent_dists.csv')

# Get train train latent dists.
d2 = pairwise_distances(training_latent,training_latent,n_jobs=30)
df2 = pd.DataFrame(data=d2, index=df_train['CoRE_name'].tolist(), columns=df_train['CoRE_name'].tolist())
df2.to_csv('thermal_train_latent_dists.csv')
