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
import time
import sys
import importlib, final_project_lib
from sklearn.model_selection import train_test_split

from final_project_lib import *
import objects

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Main code
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


##!#build_trop_training_dataset('2018')
##!#
##!#sys.exit()

####automate_TROP_input_process('201804', minlat = 55, remove_trop_files = True, \
####                            begin_date = '20180401', end_date = '20180430', \
####                            remove_data_date_dir = True)
# Read training images
# --------------------
#dat_name = 'default'
dat_name = 'OMI'
month = 7
l_only_highAI = True
l_only_lowAI = False
l_stretch_all = False
l_include_tropomi = True
num_epochs = 20000
batch_size = 32
low_res = False

l_load_model = True

objects.dat_name = dat_name
#read_train_dataset(dat_name = dat_name, low_res = low_res)
read_train_dataset(dat_name = dat_name, month = month, low_res = low_res, \
    include_tropomi = l_include_tropomi, l_remove_lower_passes = True, \
    only_highAI = l_only_highAI, only_lowAI = l_only_lowAI, \
    stretch_all = l_stretch_all)

# Extract testing images
# ----------------------
num_test = 100
objects.num_test_image = num_test
ranges = np.arange(0, objects.train_images.shape[0])
#test_idxs = random.sample(range(0,objects.train_images.shape[0]) , num_test)
#tf_no_idxs = [num not in test_idxs for num in ranges]
#train_idxs = ranges[tf_no_idxs]

train_idxs, test_idxs = train_test_split(ranges, test_size = num_test)

print(train_idxs.shape, test_idxs.shape)

objects.test_images  = objects.train_images[test_idxs,:,:]
objects.train_images = objects.train_images[train_idxs,:,:]
#objects.test_images  = objects.train_images[test_idxs,:,:,:]
#objects.train_images = objects.train_images[train_idxs,:,:,:]
if(dat_name == 'OMI'):
    objects.min_loss = 60
    objects.test_lat     = objects.train_lat[test_idxs,:,:]
    objects.test_lon     = objects.train_lon[test_idxs,:,:]
    objects.test_time    = objects.train_time[test_idxs]
    objects.train_lat    = objects.train_lat[train_idxs,:,:]
    objects.train_lon    = objects.train_lon[train_idxs,:,:]
    objects.train_time   = objects.train_time[train_idxs]

# Set constants from the training dataset size
# --------------------------------------------
objects.buffer_size = objects.train_images.shape[0]
#objects.batch_size = 256
#objects.batch_size = 128
#objects.batch_size = 64
objects.batch_size = batch_size
objects.image_size = objects.train_images.shape[1]
#objects.image_size = 28

print(objects.buffer_size, objects.batch_size, objects.image_size)

# Train the model
# ----------------
if(l_load_model):
    objects.generator     = tf.keras.models.load_model('test_generator_20240430_20000_highAI_july.keras')
    objects.discriminator = tf.keras.models.load_model('test_discriminator_20240430_20000_highAI_july.keras')
    #objects.generator     = tf.keras.models.load_model('test_generator_20240430_20000_lowAI_july.keras')
    #objects.discriminator = tf.keras.models.load_model('test_discriminator_20240430_20000_lowAI_july.keras')

    #objects.generator     = tf.keras.models.load_model('test_generator_20240423_10000_lowAI.keras')
    #objects.discriminator = tf.keras.models.load_model('test_discriminator_20240423_10000_lowAI.keras')
else:
    l_load_checkpoints = False
    train_model(EPOCHS = num_epochs, l_load_checkpoints = l_load_checkpoints)

print("PRINTING GENERATOR STUFF")
objects.generator.summary()
print("PRINTING DISCRIMINATOR STUFF")
objects.discriminator.summary()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Image completion stuff
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Set up image mask
# -----------------
np_mask = np.ones((objects.image_size,objects.image_size))
np_mask[:,16:22] = 0.0
objects.np_mask = np_mask

# Generate a figure showing the original, masked, and completed OMI images
# ------------------------------------------------------------------------
#idxs = np.array([15,25,35,45])
#plot_validate_multiple(idxs)

# Generate a figure showing the error statistics
# ----------------------------------------------
#valid_idxs = np.arange(num_test)
valid_idxs = np.array([20,30,40,50])
plot_validate_multiple(valid_idxs, complete = True, test = True, save = True)

#plot_validate_stats_multiple(50, complete = True, noise = True, \
#    test = True, save = True)
