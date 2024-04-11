#!/usr/bin/env python

"""
  NAME:
    final_project_complete.py
    
  PURPOSE:
    Run Python functions for training and running the DCGAN

    NOTE: BE SURE TO ACTIVATE THE pytf ENVIRONMENT BEFORE 
        RUNNING

  SYNTAX:
    python final_project_complete.py

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2022/02/27
      Written (modification of the Boundless colab tutorial
        from Tensorflow)
"""

import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import tensorflow as tf
from tensorflow.keras import layers
from tensorflow.keras.layers import *
from tensorflow.keras import Input
from tensorflow.keras import Model
import time
import sys
import random
import h5py
from sklearn.model_selection import train_test_split


#tf.debugging.set_log_device_placement(True)

print("GPU available:", tf.config.list_physical_devices('GPU'))
#print("GPU available:", tftest.is_gpu_available())
print("Built with CUDA:", tf.test.is_built_with_cuda())

inp1 = np.array([i-1 for i in range(3000)], dtype = float)
inp2 = np.array([i-1 for i in range(3000)], dtype = float)
inp3 = np.array([i-1 for i in range(3000)], dtype = float)
#sumd = np.array([(inptt[0] + inptt[1]) \
sumd = np.array([(inptt[0] + inptt[1] + inptt[2]) \
                 for inptt in zip(inp1, inp2, inp3)], dtype = float)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Prep the neural network
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


l_load_model = True
if(l_load_model):
    model = tf.keras.models.load_model('test_model2.keras')
else:

    """
    input_shape = 2
    model = tf.keras.Sequential([
            tf.keras.layers.Dense(units = 4, activation = 'relu', 
                                  input_shape = (input_shape,)),
            tf.keras.layers.Dense(units = 1, activation = 'linear')
    
        ])
    """


    x1 = Input(shape = (1,))
    x2 = Input(shape = (1,))
    x3 = Input(shape = (1,))

    input_layer = concatenate([x1, x2, x3], name = 'input')
    hidden_layer1 = Dense(units = 4, activation = 'relu', name = 'hidden1')(input_layer)
    hidden_layer2 = Dense(units = 4, activation = 'relu', name = 'hidden2')(hidden_layer1)
    prediction = Dense(1, activation = 'linear')(hidden_layer2)

    model = Model(inputs = [x1, x2, x3], outputs = prediction)
    #model = Model(inputs = [x1, x2, x3], outputs = prediction)





    model.compile(loss = 'mean_squared_error', optimizer = 'adam', metrics = ['mae'])
    model.summary()
    model.fit([inp1, inp2, inp3], sumd, epochs=300, batch_size = 32, verbose = 2)
    #model.fit([inp1, inp2, inp3], sumd, epochs=300, batch_size = 32, verbose = 2)

    """
    model = tf.keras.Sequential([
            tf.keras.layers.Dense(units = 8, activation = 'relu', 
                                  input_shape = (input_shape,)),
            tf.keras.layers.Dense(units = 16, activation = 'relu'), 
            tf.keras.layers.Dense(units = 8, activation = 'relu'), 
            tf.keras.layers.Dense(units = 1, activation = 'linear')
            #tf.keras.layers.Dense(units = 32, activation = 'relu'), 
    
        ])
    
    model.compile(optimizer = 'adam',\
                  loss = 'mae')
                  #loss = 'mae')

    losses = model.fit(x_train, y_train, \
                       validation_data = (x_val, y_val), \
                       batch_size = 128, \
                       epochs = 100,
            )
    """
    
    
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    # Train the neural network
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Test the neural network
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#test_out = model([np.expand_dims(np.array(4), 0), \
#                  np.expand_dims(np.array(5), 0)]).numpy()
#test_out = model([np.expand_dims(np.array(4), 0), \
#                  np.expand_dims(np.array(5), 0), \
#                  np.expand_dims(np.array(6), 0)]).numpy()
test_out = model.predict([np.expand_dims(np.array(4), 0), \
                          np.expand_dims(np.array(5), 0), \
                          np.expand_dims(np.array(6), 0)])
print(test_out)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Test saving the neural network
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

model.save('test_model2.keras')



print("SUCCESS")
