#!/usr/bin/env python
"""
  NAME:

  PURPOSE:

  PYTHON VERSION:

  MODULES:
    
"""

import numpy as np
import sys
from glob import glob
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt
#import matplotlib.colors as cm
#import matplotlib.colors as mc
#import matplotlib.cm as cm
#import matplotlib.dates as mdates
import cartopy.crs as ccrs
import cartopy.feature as cfeature
#from matplotlib.patches import Polygon
#import matplotlib.path as mpath
#import matplotlib.patches as mpatches
#import matplotlib.colors as color
#from matplotlib.colors import rgb2hex,Normalize
#from matplotlib.cm import ScalarMappable
#from matplotlib.colorbar import ColorbarBase
import matplotlib.gridspec as gridspec
from scipy import stats
import h5py
import subprocess
from scipy.stats import pearsonr,spearmanr
#from sklearn.linear_model import HuberRegressor
#from sklearn.preprocessing import StandardScaler
#import pandas as pd
#from scipy.signal import argrelextrema, find_peaks
import time
import tensorflow as tf
from tensorflow.keras import layers

sys.path.append('/home/bsorenson')
from python_lib import circle, plot_subplot_label, plot_lat_circles
sys.path.append('/home/bsorenson/Research/OMI/')
from OMILib import *

# contrast_stretch performs contrast stretching on any data passed to the
# function.
def contrast_stretch(data):
    y1 = 0.
    y2 = 255.
    x1 = -11. # or -3?
    x2 = 12.
    #x1 = np.min(data)
    #x2 = np.max(data)
    stretch_data = y1 + ((data - x1)/(x2 - x1)) * (y2 - y1)
    return stretch_data

# contrast_stretch performs contrast stretching on any data passed to the
# function.
def contrast_unstretch(data):
    y1 = -11. # or -3.?
    y2 = 12.
    x1 = 0.
    x2 = 255.
    #x1 = np.min(data)
    #x2 = np.max(data)
    unstretch_data = y1 + ((data - x1)/(x2 - x1)) * (y2 - y1)
    return unstretch_data

def select_OMI_section(OMI_data, begin_idx):
    x_max = begin_idx
    x_min = x_max - 300
    localLAT = OMI_data['LAT'][x_min:x_max, :]
    localLON = OMI_data['LON'][x_min:x_max, :]
    localAI  = OMI_data['UVAI'][x_min:x_max, :]

    # Extract every 10th latitude and longitude
    # -----------------------------------------
    new_LAT = localLAT[::5, :]
    new_LON = localLON[::5, :]

    # Average the OMI data along each row down to 60 x 60
    # ---------------------------------------------------
    new_AI = localAI.reshape(-1, 5, localAI.shape[1]).mean(axis = 1)
   
    return new_LAT, new_LON, new_AI
 

# Set up to work with the output from readOMI_swath_hdf
def resample_OMI_swath(OMI_data):
    
    # Pull only the top 600 values from the OMI data swath
    # ----------------------------------------------------
    #if(OMI_data['UVAI'].shape[0] < 1600):
    x_max = np.where(OMI_data['UVAI'][:,59].mask == False)[0][-1]
    #x_max = OMI_data['UVAI'].shape[0]

    LAT1, LON1, AI1 = select_OMI_section(OMI_data, x_max)
    x_max = x_max - 300
    LAT2, LON2, AI2 = select_OMI_section(OMI_data, x_max)

    new_OMI = {}
    new_OMI['scene_1'] = {}
    new_OMI['scene_2'] = {}
    new_OMI['scene_1']['LAT']  = LAT1
    new_OMI['scene_1']['LON']  = LON1
    new_OMI['scene_1']['UVAI'] = AI1
    new_OMI['scene_2']['LAT']  = LAT2
    new_OMI['scene_2']['LON']  = LON2
    new_OMI['scene_2']['UVAI'] = AI2

    return new_OMI

def build_training_dataset(date_str):

    # Figure out if date string is a year, month, or day
    # --------------------------------------------------
    if(len(date_str) == 4):
        str_fmt = '%Y'
        out_fmt = '%Ym'
    elif(len(date_str) == 6):
        str_fmt = '%Y%m'
        out_fmt = '%Ym%m'
    elif(len(date_str) == 8):
        str_fmt = '%Y%m%d'
        out_fmt = '%Ym%m%d'
    else:
        print("INVALID DATE STRING. RETURNING")
        return

    dt_date_str = datetime.strptime(date_str, str_fmt)
         
    # For whatever date period is desired, use glob
    # to find the number of associated files
    # ---------------------------------------------
    files = glob(dt_date_str.strftime(\
        '/home/bsorenson/data/OMI/H5_files/OMI-A*_' + out_fmt + '*.he5'))

    # The length of the returned list (times two) will
    # be the size of array needed for the netCDF file
    # ------------------------------------------------
    n_files = len(files)
    n_images = n_files * 2
    n_x = n_y = 60
    print(n_files)

    # ----------------------
    # Set up the netCDF file
    # ----------------------
    file_name = 'train_dataset_' + date_str + '.nc'
    nc = Dataset(file_name,'w',format='NETCDF4')
    
    # Use the dimension size variables to actually create dimensions in 
    # the file.
    # ----------------------------------------------------------------- 
    n_image_dim = nc.createDimension('n_img',n_images)
    n_lat_dim   = nc.createDimension('n_x',n_x)
    n_lon_dim   = nc.createDimension('n_y',n_y)

    # Create variables for the three dimensions
    # -----------------------------------------
    images = nc.createVariable('images','i2',('n_img'))
    images.description = 'Number of OMI image subsets in the dataset'
    ##x_dim = nc.createVariable('x_dim','i2',('n_x'))
    ##x_dim.description = 'Indices in along-track direction'
    ##x_dim.units = 'None'
    ##y_dim = nc.createVariable('y_dim','i2',('n_y'))
    ##y_dim.description = 'Indices in across-track direction'
    ##y_dim.units = 'None'

    # Create a variable for the AI images, dimensioned using all three
    # dimensions.
    # ----------------------------------------------------------------  
    AI = nc.createVariable('AI','f4',('n_img','n_x','n_y'))
    AI.description = 'UV Aerosol Index'
    AI.units = 'None'
    LAT = nc.createVariable('LAT','f4',('n_img','n_x','n_y'))
    LAT.description = 'Latitude'
    LAT.units = 'degrees N'
    LON = nc.createVariable('LON','f4',('n_img','n_x','n_y'))
    LON.description = 'Longitude'
    LON.units = 'degrees E'
    TIME = nc.createVariable('TIME','f8',('n_img'))
    TIME.description = 'Seconds since 00:00 UTC 01 January 2005'
    TIME.units = 'seconds'

    base_date = datetime(year = 2005, month = 1, day = 1, hour = 0, minute = 0)

    # Loop over the glob list
    # -----------------------
    for ii, fname in enumerate(files):
        dt_name = datetime.strptime(fname.strip().split('/')[-1][20:34],'%Ym%m%dt%H%M')

        time_diff = (dt_name - base_date).total_seconds()

        # Read and resample each file
        # ---------------------------
        OMI_data = readOMI_swath_hdf(dt_name.strftime('%Y%m%d%H%M'), 'control', latmin = 23, \
            skiprows = None)

        OMI_resample = resample_OMI_swath(OMI_data)

        # Insert the data for the scene(s) from the current file into the
        # netCDF file
        # ---------------------------------------------------------------
        TIME[ii*2] = time_diff
        TIME[ii*2+1] = time_diff
        AI[(ii*2),:,:]    = OMI_resample['scene_1']['UVAI'] 
        LAT[(ii*2),:,:]   = OMI_resample['scene_1']['LAT'] 
        LON[(ii*2),:,:]   = OMI_resample['scene_1']['LON'] 
        AI[(ii*2)+1,:,:]  = OMI_resample['scene_2']['UVAI'] 
        LAT[(ii*2)+1,:,:] = OMI_resample['scene_2']['LAT'] 
        LON[(ii*2)+1,:,:] = OMI_resample['scene_2']['LON'] 

    # Close the netCDF file, which actually saves it.
    # -----------------------------------------------
    nc.close()

    print('netCDF file saved')

# Set up to work with the netCDF file object from build_training_dataset
def plot_train_data(OMI_data, idx, map = False, save = False):
    base_date = datetime(year = 2005, month = 1, day = 1, hour = 0, minute = 0)
   
    local_time = base_date + timedelta(seconds = OMI_data['TIME'][idx]/1) 
   
    plt.close()
    fig = plt.figure()
    if(map):

        avg_lat = np.nanmean(OMI_data['LAT'][idx,:,:])

        if(avg_lat > 70):
            mapcrs = ccrs.NorthPolarStereo()
        else:
            mapcrs = ccrs.LambertConformal(central_latitude = avg_lat,\
                central_longitude = \
                np.nanmean(OMI_data['LON'][idx,:,:]))
        ax = fig.add_subplot(1,1,1, projection = mapcrs)
        
        ax.coastlines()
        ax.add_feature(cfeature.BORDERS)
        mesh = ax.pcolormesh(OMI_data['LON'][idx,:,:], \
            OMI_data['LAT'][idx,:,:], OMI_data['AI'][idx,:,:], 
            transform = ccrs.PlateCarree())
        
    else:
        ax = fig.add_subplot(1,1,1)
        mesh = ax.imshow(OMI_data['AI'][idx])
    
    cbar = plt.colorbar(mesh, ax = ax, label = 'AI')
    ax.set_title(local_time.strftime('%Y-%m-%d %H:%M'))
    plt.show()
 
# Plot the data in just a pcolormesh
# ----------------------------------
def plot_resample_OMI_nomap(date_str, save = False):

    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')

    # Read the original AI data
    # -------------------------
    OMI_data = readOMI_swath_hdf(date_str, 'control', latmin = 23, \
        skiprows = None)

    print("OMI_data.shape", OMI_data['UVAI'].shape)
    
    # Resample the OMI data
    # ---------------------
    OMI_resample = resample_OMI_swath(OMI_data)

    # Set up the figure
    # -----------------
    plt.close('all')
    fig = plt.figure(figsize = (12,6))

    # Plot the original and resampled data
    ax0 = fig.add_subplot(1,3,1)
    ax1 = fig.add_subplot(1,3,2)
    ax2 = fig.add_subplot(1,3,3)

    # Plot the data
    # -------------
    ax0.pcolormesh(OMI_data['UVAI'])
    ax1.pcolormesh(OMI_resample['scene_1']['UVAI'])
    ax2.pcolormesh(OMI_resample['scene_2']['UVAI'])

    plt.suptitle(date_str)

    plt.show()

# Plot the data on a map
# ----------------------
def plot_resample_OMI(date_str, save = False):

    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')

    # Read the original AI data
    # -------------------------
    OMI_data = readOMI_swath_hdf(date_str, 'control', latmin = 50, \
        skiprows = None)

    print("OMI_data.shape", OMI_data['UVAI'].shape)
    
    # Resample the OMI data
    # ---------------------
    OMI_resample = resample_OMI_swath(OMI_data)

    # Set up the figure
    # -----------------
    plt.close('all')
    fig = plt.figure(figsize = (8,6))

    # Plot the original and resampled data
    mapcrs  = ccrs.NorthPolarStereo()
    datacrs = ccrs.PlateCarree()
    ax0 = fig.add_subplot(1,2,1, projection = mapcrs)
    ax1 = fig.add_subplot(1,2,2, projection = mapcrs)

    # Plot the data
    # -------------
    ax0.pcolormesh(OMI_data['LON'], OMI_data['LAT'], OMI_data['UVAI'], \
        transform = datacrs, cmap = 'jet')
    ax0.set_extent([-180, 180, 65, 90], datacrs)
    ax0.coastlines()

    ax1.pcolormesh(OMI_resample['LON'], OMI_resample['LAT'], \
        OMI_resample['UVAI'], transform = datacrs, cmap = 'jet')
    ax1.set_extent([-180, 180, 65, 90], datacrs)
    ax1.coastlines()

    plt.suptitle(date_str)

    plt.show()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =  
#
# DCGAN stuff
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =  

# This method returns a helper function to compute cross entropy loss
cross_entropy = tf.keras.losses.BinaryCrossentropy(from_logits=True)
BUFFER_SIZE = 60000
BATCH_SIZE = 64
#BATCH_SIZE = 256
noise_dim = 100
num_examples_to_generate = 16
x_dim = 15

##!#def make_generator_model():
##!#    model = tf.keras.Sequential()
##!#    model.add(layers.Dense(7*7*256, use_bias=False, input_shape=(100,)))
##!#    model.add(layers.BatchNormalization())
##!#    model.add(layers.LeakyReLU())
##!#
##!#    model.add(layers.Reshape((7, 7, 256)))
##!#    assert model.output_shape == (None, 7, 7, 256)  # Note: None is the batch size
##!#
##!#    model.add(layers.Conv2DTranspose(128, (5, 5), strides=(1, 1), padding='same', use_bias=False))
##!#    assert model.output_shape == (None, 7, 7, 128)
##!#    model.add(layers.BatchNormalization())
##!#    model.add(layers.LeakyReLU())
##!#
##!#    model.add(layers.Conv2DTranspose(64, (5, 5), strides=(2, 2), padding='same', use_bias=False))
##!#    assert model.output_shape == (None, 14, 14, 64)
##!#    model.add(layers.BatchNormalization())
##!#    model.add(layers.LeakyReLU())
##!#
##!#    model.add(layers.Conv2DTranspose(1, (5, 5), strides=(2, 2), padding='same', use_bias=False, activation='tanh'))
##!#    assert model.output_shape == (None, 28, 28, 1)
##!#
##!#    return model

##!#def make_discriminator_model():
##!#    model = tf.keras.Sequential()
##!#    model.add(layers.Conv2D(64, (5, 5), strides=(2, 2), padding='same',
##!#                                     input_shape=[28, 28, 1]))
##!#    model.add(layers.LeakyReLU())
##!#    model.add(layers.Dropout(0.3))
##!#
##!#    model.add(layers.Conv2D(128, (5, 5), strides=(2, 2), padding='same'))
##!#    model.add(layers.LeakyReLU())
##!#    model.add(layers.Dropout(0.3))
##!#
##!#    model.add(layers.Flatten())
##!#    model.add(layers.Dense(1))
##!#
##!#    return model

def make_generator_model():
    model = tf.keras.Sequential()
    model.add(layers.Dense(15*15*256, use_bias=False, input_shape=(100,)))
    #model.add(layers.Dense(15*15*64, use_bias=False, input_shape=(100,)))
    model.add(layers.BatchNormalization())
    model.add(layers.LeakyReLU())

    model.add(layers.Reshape((15, 15, 256)))
    assert model.output_shape == (None, 15, 15, 256)  # Note: None is the batch size
    #model.add(layers.Reshape((15, 15, 64)))
    #assert model.output_shape == (None, 15, 15, 64)  # Note: None is the batch size

    model.add(layers.Conv2DTranspose(128, (5, 5), strides=(1, 1), padding='same', use_bias=False))
    assert model.output_shape == (None, 15, 15, 128)
    #model.add(layers.Conv2DTranspose(32, (5, 5), strides=(1, 1), padding='same', use_bias=False))
    #assert model.output_shape == (None, 15, 15, 32)
    model.add(layers.BatchNormalization())
    model.add(layers.LeakyReLU())

    model.add(layers.Conv2DTranspose(64, (5, 5), strides=(2, 2), padding='same', use_bias=False))
    assert model.output_shape == (None, 30, 30, 64)
    #model.add(layers.Conv2DTranspose(16, (5, 5), strides=(2, 2), padding='same', use_bias=False))
    #assert model.output_shape == (None, 30, 30, 16)
    model.add(layers.BatchNormalization())
    model.add(layers.LeakyReLU())

    model.add(layers.Conv2DTranspose(1, (5, 5), strides=(2, 2), padding='same', use_bias=False, activation='tanh'))
    assert model.output_shape == (None, 60, 60, 1)

    return model

def make_discriminator_model():
    model = tf.keras.Sequential()
    model.add(layers.Conv2D(64, (5, 5), strides=(2, 2), padding='same',
                                     input_shape=[60, 60, 1]))
    #model.add(layers.Conv2D(16, (5, 5), strides=(2, 2), padding='same',
    #                                 input_shape=[60, 60, 1]))
    model.add(layers.LeakyReLU())
    model.add(layers.Dropout(0.3))

    model.add(layers.Conv2D(128, (5, 5), strides=(2, 2), padding='same'))
    #model.add(layers.Conv2D(32, (5, 5), strides=(2, 2), padding='same'))
    model.add(layers.LeakyReLU())
    model.add(layers.Dropout(0.3))

    model.add(layers.Flatten())
    model.add(layers.Dense(1))

    return model

def generator_loss(fake_output):
    return cross_entropy(tf.ones_like(fake_output), fake_output)

def discriminator_loss(real_output, fake_output):
    real_loss = cross_entropy(tf.ones_like(real_output), real_output)
    fake_loss = cross_entropy(tf.zeros_like(fake_output), fake_output)
    total_loss = real_loss + fake_loss
    return total_loss


##!## Notice the use of `tf.function`
##!## This annotation causes the function to be "compiled".
##!#@tf.function
##!#def train_step(images):
##!#    noise = tf.random.normal([BATCH_SIZE, noise_dim])
##!#
##!#    with tf.GradientTape() as gen_tape, tf.GradientTape() as disc_tape:
##!#        generated_images = generator(noise, training=True)
##!#
##!#        real_output = discriminator(images, training=True)
##!#        fake_output = discriminator(generated_images, training=True)
##!#
##!#        gen_loss = generator_loss(fake_output)
##!#        disc_loss = discriminator_loss(real_output, fake_output)
##!#
##!#    gradients_of_generator = gen_tape.gradient(gen_loss, generator.trainable_variables)
##!#    gradients_of_discriminator = disc_tape.gradient(disc_loss, discriminator.trainable_variables)
##!#
##!#    generator_optimizer.apply_gradients(zip(gradients_of_generator, generator.trainable_variables))
##!#    discriminator_optimizer.apply_gradients(zip(gradients_of_discriminator, discriminator.trainable_variables))
##!#
##!#def train(dataset, epochs):
##!#    for epoch in range(epochs):
##!#        start = time.time()
##!#
##!#        for image_batch in dataset:
##!#            train_step(image_batch)
##!#
##!#        # Produce images for the GIF as you go
##!#        display.clear_output(wait=True)
##!#        generate_and_save_images(generator,
##!#                                 epoch + 1,
##!#                                 seed)
##!#
##!#        # Save the model every 15 epochs
##!#        if (epoch + 1) % 15 == 0:
##!#            checkpoint.save(file_prefix = checkpoint_prefix)
##!#
##!#        print ('Time for epoch {} is {} sec'.format(epoch + 1, time.time()-start))
##!#
##!#    # Generate after the final epoch
##!#    display.clear_output(wait=True)
##!#    generate_and_save_images(generator,
##!#                             epochs,
##!#                             seed)

def generate_and_save_images(model, epoch, test_input):
  # Notice `training` is set to False.
  # This is so all layers run in inference mode (batchnorm).
  predictions = model(test_input, training=False)

  fig = plt.figure(figsize=(4, 4))

  for i in range(predictions.shape[0]):
      plt.subplot(4, 4, i+1)
      plt.imshow(predictions[i, :, :, 0] * 127.5 + 127.5, cmap='gray')
      plt.axis('off')

  plt.savefig('image_at_epoch_{:04d}.png'.format(epoch))
  plt.show()



