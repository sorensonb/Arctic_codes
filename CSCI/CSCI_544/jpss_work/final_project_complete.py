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
import random
import importlib, final_project_lib

from final_project_lib import *
import objects

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Main code
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

####automate_TROP_input_process('201804', minlat = 55, remove_trop_files = True, \
####                            begin_date = '20180401', end_date = '20180430', \
####                            remove_data_date_dir = True)
# Read training images
# --------------------
#dat_name = 'default'
dat_name = 'OMI'
low_res = False
objects.dat_name = dat_name
read_train_dataset(dat_name = dat_name, low_res = low_res)
#read_train_dataset(dat_name = dat_name, month = 7, low_res = low_res)


# Extract testing images
# ----------------------
num_test = 100
objects.num_test_image = num_test
ranges = np.arange(0, objects.train_images.shape[0])
test_idxs = random.sample(range(0,objects.train_images.shape[0]) , num_test)
tf_no_idxs = [num not in test_idxs for num in ranges]
train_idxs = ranges[tf_no_idxs]

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
objects.batch_size = 128
#objects.batch_size = 64
#objects.batch_size = 256
objects.image_size = objects.train_images.shape[1]
#objects.image_size = 28

print(objects.buffer_size, objects.batch_size, objects.image_size)

# Train the model
# ---------------
train_model(EPOCHS = 1000)

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
valid_idxs = np.arange(num_test)
plot_validate_multiple(valid_idxs, complete = True, test = True, save = True)
