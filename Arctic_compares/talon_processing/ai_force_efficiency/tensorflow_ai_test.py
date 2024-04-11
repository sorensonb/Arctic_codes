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


# input variables
# - OMI AI
# - OMI SZA
# - NSIDC ICE
# - MODIS COD
# - MODIS cloud top temp?

# output variable:
# - CERES SWF

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Prep the data
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Scale all input variables to 0 - 1

data_path = '/home/blake.sorenson/OMI/arctic_comp/comp_data/'

fdates = [
        '201807040005',\
        '201807040144',\
        '201807040322',\
        '201807041633',\
        '201807041812',\
        '201807041951',\
        '201807042130',\
        '201807042309',\
        '201807050048',\
        '201807050227',\
        '201807051538',\
        '201807051717',\
        '201807051856',\
        '201807052034',\
        '201807052213',\
        '201807052352'
    ]

files = [data_path + 'colocated_subset_' + fdd + '.hdf5' for fdd in fdates]

#files = ['/home/blake.sorenson/OMI/arctic_comp/comp_data/colocated_subset_201807052213.hdf5']


# NOTE: added the "==254" section so that it temporarily removes land data.
#       Want to test the system over only ocean and/or ice to start with
def select_data(data, min_ai, max_ai, minlat, min_swf, max_swf, max_cod):
    local_data = np.ma.masked_invalid(data['omi_uvai_raw'])
    return np.ma.masked_where( (local_data < min_ai) | \
                               (local_data > max_ai) | \
                               (data['omi_lat'][:,:] < minlat) | \
                               ( (data['ceres_swf'][:,:] < min_swf) | \
                               (data['ceres_swf'][:,:] > max_swf) ) | \
                               (data['modis_cod'][:,:] == -999) | \
                               (data['modis_cod'][:,:] > max_cod) | \
                               (np.isnan(data['modis_cod'][:,:]) == True) | \
                               (data['modis_cld_top_pres'][:,:] == -999) | \
                               (data['modis_cld_top_pres'][:,:] > 1020) | \
                               (np.isnan(data['modis_cld_top_pres'][:,:]) == True) | \
                               (data['nsidc_ice'][:,:] == -999) | \
                               (data['nsidc_ice'][:,:] == 251) | \
                               #(data['nsidc_ice'][:,:] == 254) | \
                               (data['nsidc_ice'][:,:] == 253), \
                               local_data\
    )
    

# Figure out the total size to insert the data
# ---------------------------------------------
#minlat = 70.
min_ai = -2.0
max_ai = 1.0
minlat = 65.    # 2024/01/10: changed to 65.
min_swf = 0.
max_swf = 3000.
max_cod = 70.
min_ice = 0.
max_ice = 500.

total_size = 0
for ff in files:
    data = h5py.File(ff,'r')

    local_data = select_data(data, min_ai, max_ai, minlat, min_swf, \
        max_swf, max_cod)

    ##!#local_data = np.ma.masked_invalid(data['omi_uvai_raw'])
    ##!#local_data = np.ma.masked_where( (abs(local_data) > 12) | \
    ##!#                                 (data['omi_lat'][:,:] < minlat) | \
    ##!#                                 ( (data['ceres_swf'][:,:] < 0) | \
    ##!#                                 (data['ceres_swf'][:,:] > 3000) ) | \
    ##!#                                 (data['modis_cod'][:,:] == -999) | \
    ##!#                                 (data['modis_cod'][:,:] > 70) | \
    ##!#                                 (np.isnan(data['modis_cod'][:,:]) == True) | \
    ##!#                                 (data['nsidc_ice'][:,:] == -999) | \
    ##!#                                 (data['nsidc_ice'][:,:] == 251) | \
    ##!#                                 (data['nsidc_ice'][:,:] == 253), \
    ##!#                                 local_data\
    ##!#)

    local_size = local_data.compressed().shape[0]
    total_size += local_size
    print(ff, local_size, total_size)

    data.close()

combined_data = {}
combined_data['omi_uvai_pert'] = np.full(total_size, np.nan)
#combined_data['omi_uvai_raw']  = np.full(total_size, np.nan)
#combined_data['modis_cld']     = np.full(total_size, np.nan)
combined_data['modis_cld_top_pres']     = np.full(total_size, np.nan)
combined_data['modis_cod']     = np.full(total_size, np.nan)
combined_data['ceres_swf']     = np.full(total_size, np.nan)
#combined_data['modis_ch7']     = np.full(total_size, np.nan)
combined_data['omi_sza']       = np.full(total_size, np.nan)
combined_data['omi_lat']       = np.full(total_size, np.nan)
combined_data['nsidc_ice']     = np.full(total_size, np.nan)

print("Loading data")

# Loop back over the files and insert the data into the structure
# ---------------------------------------------------------------
total_size = 0
beg_idx = 0
end_idx = 0
for ff in files:

    data = h5py.File(ff,'r')


    local_data = select_data(data, min_ai, max_ai, minlat, min_swf, \
        max_swf, max_cod)


    # NOTE: Changed the omi variable here from "pert" to "raw" on 20230623.
    #       This move should allow for coloc data to be read after 2020
    ##!#local_data = np.ma.masked_invalid(data['omi_uvai_raw'])

    ##!#local_data = np.ma.masked_where( (abs(local_data) > 12) | \
    ##!#                                 (data['omi_lat'][:,:] < minlat) | \
    ##!#                                 ( (data['ceres_swf'][:,:] < 0) | \
    ##!#                                 (data['ceres_swf'][:,:] > 3000) ) | \
    ##!#                                 (data['modis_cod'][:,:] == -999) | \
    ##!#                                 (data['modis_cod'][:,:] > 70) | \
    ##!#                                 (np.isnan(data['modis_cod'][:,:]) == True) | \
    ##!#                                 (data['nsidc_ice'][:,:] == -999) | \
    ##!#                                 (data['nsidc_ice'][:,:] == 251) | \
    ##!#                                 (data['nsidc_ice'][:,:] == 253), 
    ##!#                                 local_data\
    ##!#)

    local_size = local_data.compressed().shape[0]

    beg_idx = end_idx
    end_idx = beg_idx + local_size

    for tkey in combined_data.keys():
        combined_data[tkey][beg_idx:end_idx] = \
            data[tkey][~local_data.mask]

    print(local_size)
    total_size += local_size

    data.close()

print(combined_data['omi_uvai_pert'].shape)

# Fix the cloud top pressure so that '0' values (which I assume mean 'clear-sky'?)
# are at 1025 rather than 0. This will hopefully give continuity from low clouds
# to no clouds
combined_data['modis_cld_top_pres'] = \
    np.where(combined_data['modis_cld_top_pres'] == 0., 1025., \
    combined_data['modis_cld_top_pres'])


print(np.min(combined_data['omi_uvai_pert']), np.max(combined_data['omi_uvai_pert']))
print(np.min(combined_data['ceres_swf']), np.max(combined_data['ceres_swf']))
print(np.min(combined_data['modis_cod']), np.max(combined_data['modis_cod']))
print(np.min(combined_data['modis_cld_top_pres']), np.max(combined_data['modis_cld_top_pres']))
print(np.min(combined_data['omi_sza']), np.max(combined_data['omi_sza']))
print(np.min(combined_data['nsidc_ice']), np.max(combined_data['nsidc_ice']))

combined_data['nsidc_ice'][:] = \
    np.where(combined_data['nsidc_ice'][:] == 254., 101., combined_data['nsidc_ice'][:])

min_max_dict = {}

key_variables = ['omi_uvai_pert', 'omi_sza', 'modis_cod', 'modis_cld_top_pres', 'nsidc_ice', 'ceres_swf']

for key in key_variables:
    min_max_dict[key] = {}
    min_max_dict[key]['min'] = np.min(combined_data[key])
    min_max_dict[key]['max'] = np.max(combined_data[key])

    drange = min_max_dict[key]['max'] - min_max_dict[key]['min']
    combined_data[key] = ((combined_data[key] - min_max_dict[key]['min']) / drange) * 100.
    #combined_data[key] = ((combined_data[key] - min_max_dict[key]['min']) / drange) * 100.

print(np.min(combined_data['omi_uvai_pert']), np.max(combined_data['omi_uvai_pert']))
print(np.min(combined_data['ceres_swf']), np.max(combined_data['ceres_swf']))
print(np.min(combined_data['modis_cod']), np.max(combined_data['modis_cod']))
print(np.min(combined_data['modis_cld_top_pres']), np.max(combined_data['modis_cld_top_pres']))
print(np.min(combined_data['omi_sza']), np.max(combined_data['omi_sza']))
print(np.min(combined_data['nsidc_ice']), np.max(combined_data['nsidc_ice']))

pcnt_test = 0.25
num_test = int(combined_data['omi_uvai_pert'].shape[0] * pcnt_test)
num_train = combined_data['omi_uvai_pert'].shape[0] - num_test
ranges = np.arange(0, combined_data['omi_uvai_pert'].shape[0])

train_idxs, test_idxs = train_test_split(ranges, test_size = num_test)

print(train_idxs.shape, test_idxs.shape)

# Input format: OMI SZA, NSIDC ICE, MODIS COD
x_train = np.array([combined_data['omi_sza'][train_idxs], \
                    combined_data['nsidc_ice'][train_idxs], \
                    combined_data['modis_cod'][train_idxs], \
                    combined_data['modis_cld_top_pres'][train_idxs]\
                  ]).T
y_train = combined_data['ceres_swf'][train_idxs]

x_val   = np.array([combined_data['omi_sza'][test_idxs], \
                    combined_data['nsidc_ice'][test_idxs], \
                    combined_data['modis_cod'][test_idxs], \
                    combined_data['modis_cld_top_pres'][test_idxs]\
                  ]).T
y_val   = combined_data['ceres_swf'][test_idxs]

print(np.max(y_train), np.max(y_val))


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Prep the neural network
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Inputs
# - OMI SZA
# - MODIS COD
# - NSIDC ICE
# Future inputs
# - MODIS CLD TOP TEMP/PRESS
# - OMI ROW NUMBER

l_load_model = True
if(l_load_model):
    model = tf.keras.models.load_model('test_model_noland2.keras')
    model.summary()
else:

    x1 = Input(shape = (1,))
    x2 = Input(shape = (1,))
    x3 = Input(shape = (1,))
    x4 = Input(shape = (1,))

    input_layer = concatenate([x1, x2, x3, x4], name = 'input')
    hidden_layer1 = Dense(units = 8, activation = 'relu', name = 'hidden1')(input_layer)
    hidden_layer2 = Dense(units = 8, activation = 'relu', name = 'hidden2')(hidden_layer1)
    hidden_layer4 = Dense(units = 8, activation = 'relu', name = 'hidden4')(hidden_layer2)
    hidden_layer3 = Dense(units = 4, activation = 'relu', name = 'hidden3')(hidden_layer4)
    prediction = Dense(1, activation = 'linear')(hidden_layer4)

    model = Model(inputs = [x1, x2, x3, x4], outputs = prediction)
    #model = Model(inputs = [x1, x2, x3], outputs = prediction)





    model.compile(loss = 'mean_squared_error', optimizer = 'adam', metrics = ['mae'])
    model.summary()
    losses = model.fit([combined_data['omi_sza'][train_idxs], \
                        combined_data['nsidc_ice'][train_idxs], \
                        combined_data['modis_cod'][train_idxs], \
                        combined_data['modis_cld_top_pres'][train_idxs]], \
                        combined_data['ceres_swf'][train_idxs], epochs=100, batch_size = 32, verbose = 2)
    #model.fit([inp1, inp2, inp3], sumd, epochs=300, batch_size = 32, verbose = 2)


    """
    input_shape = 4
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
    
    model.summary()
    
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    #
    # Train the neural network
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    losses = model.fit(x_train, y_train, \
                       validation_data = (x_val, y_val), \
                       batch_size = 128, \
                       epochs = 100,
            )
    """

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Test the neural network
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#calc = model.predict(x_train[0:10,:])
test_out = np.array([model.predict([np.expand_dims(np.array(combined_data['omi_sza'][test_idxs][xx]), 0), \
                          np.expand_dims(np.array(combined_data['nsidc_ice'][test_idxs][xx]), 0), \
                          np.expand_dims(np.array(combined_data['modis_cod'][test_idxs][xx]), 0), \
                          np.expand_dims(np.array(combined_data['modis_cld_top_pres'][test_idxs][xx]), 0)]) \
                    for xx in range(10)]).squeeze()

print(  ((combined_data['ceres_swf'][test_idxs][0:10] / 100) * drange) + min_max_dict['ceres_swf']['min'])
print(  ((test_out / 100) * drange) + min_max_dict['ceres_swf']['min'])

##!#print(x_train[0:10,:])
##!#print(    ((y_val[0:10].T / 100) * drange) + min_max_dict['ceres_swf']['min'])
##!#print(    ((calc / 100) * drange) + min_max_dict['ceres_swf']['min'])

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Test saving the neural network
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

#model.save('test_model_allsfc.keras')
model.save('test_model_noland2.keras')

# HERE: Read in a single one of the OMI colocation files, scale the input data
#       based on the scaling from above, and then use the model to calculate
#       the predicted SWF values at each gridpoint in the OMI file. Want
#       to then plot the results and see if there is any spatial information
#       that makes it through the model
#       
#       ALSO: Must remember that this is for ice/ocean values ONLY. Land data
#       are not included right now.
test_time = '201807052213'
infile = '/home/blake.sorenson/OMI/arctic_comp/comp_data/colocated_subset_' + \
    test_time + '.hdf5'

# Set up input data
data = h5py.File(infile, 'r')

# Set up output data
outname = 'test_calc_out.hdf5'
dset = h5py.File(outname, 'w')

calc_swf = np.full(data['ceres_swf'].shape, np.nan)

for ii in range(calc_swf.shape[0]):
    for jj in range(calc_swf.shape[1]):
        # Check the value here
        # --------------------  
        if( (data['omi_lat'][ii,jj] >= minlat) & \
            (data['modis_cod'][ii,jj] != -999) & \
            (data['modis_cod'][ii,jj] < max_cod) & \
            (np.isnan(data['modis_cod'][ii,jj]) == False) & \
            (data['modis_cld_top_pres'][ii,jj] != -999) & \
            (data['modis_cld_top_pres'][ii,jj] <= 1020) & \
            (np.isnan(data['modis_cld_top_pres'][ii,jj]) == False) & \
            (data['nsidc_ice'][ii,jj] != -999) & \
            (data['nsidc_ice'][ii,jj] <= 100)):


            # Scale the input values here
            # ---------------------------
            new_sza = ( ( data['omi_sza'][ii,jj] - min_max_dict['omi_sza']['min']) / \
                        ( min_max_dict['omi_sza']['max'] - min_max_dict['omi_sza']['min'] )) * 100.
            new_ice = ( ( data['nsidc_ice'][ii,jj] - min_max_dict['nsidc_ice']['min']) / \
                        ( min_max_dict['nsidc_ice']['max'] - min_max_dict['nsidc_ice']['min'] )) * 100.
            new_cod = ( ( data['modis_cod'][ii,jj] - min_max_dict['modis_cod']['min']) / \
                        ( min_max_dict['modis_cod']['max'] - min_max_dict['modis_cod']['min'] )) * 100.
            new_cpr = ( ( data['modis_cld_top_pres'][ii,jj] - min_max_dict['modis_cld_top_pres']['min']) / \
                        ( min_max_dict['modis_cld_top_pres']['max'] - min_max_dict['modis_cld_top_pres']['min'] )) * 100.

            # Use the model to estimate the SWF
            # ---------------------------------
            calc_swf[ii,jj] = np.array([model.predict([np.expand_dims(np.array(new_sza), 0), \
                                      np.expand_dims(np.array(new_ice), 0), \
                                      np.expand_dims(np.array(new_cod), 0), \
                                      np.expand_dims(np.array(new_cpr), 0)]) \
                                ]).squeeze()
            calc_swf[ii,jj] = ((calc_swf[ii,jj] / 100) * \
                    (min_max_dict['ceres_swf']['max'] - min_max_dict['ceres_swf']['min'])) + \
                    min_max_dict['ceres_swf']['min']

            print(ii, jj, data['ceres_swf'][ii,jj], calc_swf[ii,jj])

# Add the calculated values to the output file
# --------------------------------------------
dset.create_dataset('calc_swf', data = calc_swf[:,:])
dset.close()

print("SUCCESS")
