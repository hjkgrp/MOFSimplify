#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 12:22:00 2023

@author: mattrivera
"""

### Import necessary starting packages and initial datasets
import sys
import pandas as pd
import tensorflow as tf
import keras


### Set working directory ###
training_dir=sys.argv[1]



### Import training data ###
x_train=pd.read_csv(training_dir+'x_train.csv')
y_train=pd.read_csv(training_dir+'y_train.csv')




### Hyperparameter dict ###
hyperspace = {'hidden_nodes':1,
              'activation': 'relu',
              'layers': 1,
              'learning_rate': 0.01,
              'dropout': 0.05
    }





 ### Model builder function ###       
def final_model(hyperspace):
    model = keras.Sequential()
    dropout_rate=hyperspace['dropout']
    model.add(keras.layers.Dropout(rate=dropout_rate, input_shape=(input_shape,)))

    # Tune the number of units in the first Dense layer
    # Choose an optimal value between 25-1000
    activation_function='relu'
    nodes=hyperspace['hidden_nodes']
   
    for i in range(0,hyperspace['layers']):
        model.add(keras.layers.Dense(units=nodes, activation=activation_function))
        model.add(keras.layers.Dropout(rate=dropout_rate))
    model.add(keras.layers.Dense(units=1))
    # Tune the learning rate for the optimizer
    hp_learning_rate = hyperspace['learning_rate']
    
    model.compile(optimizer=keras.optimizers.Adam(learning_rate=hp_learning_rate),
                  loss=keras.losses.MeanAbsoluteError(),
                  metrics=['mae'])    

    return model






### Make model object with optimum hyperparameters ###
input_shape=x_train.shape[1]
optimum_model=final_model(hyperspace)

### Train the model ###
optimum_model.fit(x=x_train,y=y_train,epochs=5000)


### Save the model object ###
optimum_model.save(training_dir+'ethane_model.h5')


    
    
    
