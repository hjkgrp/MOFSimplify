# Purpose of this file: check the number of layers in our models
# Not used in mofSimplify, just used in development
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

from keras.models import load_model
dependencies = {'precision':precision,'recall':recall,'f1':f1}
solvent_model = load_model('../model/solvent/ANN/final_model_flag_few_epochs.h5', custom_objects=dependencies)
thermal_model = load_model('../model/thermal/ANN/final_model_T_few_epochs.h5',custom_objects=dependencies)
print()
print('solvent model summary')
print(solvent_model.summary())
print()
print('thermal model summary')
print(thermal_model.summary())
print()
print('Test test A')
print(solvent_model.layers[0])
print(thermal_model.layers[0])
print('Test test B')
print(solvent_model.layers[3])
print(thermal_model.layers[3])
print('Test test C')
print(solvent_model.layers[10])
print(thermal_model.layers[10])
print('Test test D')
# print(solvent_model.layers[16]) # 16 is out of range
# print(thermal_model.layers[11]) # 11 is out of range
print(len(solvent_model.layers))
print(len(thermal_model.layers))