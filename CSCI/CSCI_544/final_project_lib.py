#!/usr/bin/env python
"""
  NAME:
    final_project_lib.py
    
  PURPOSE:
    Contains functions necessary for implementing an image-completion
    Deep Convolutional Generative Adversarial Network (DCGAN) for
    the CSCI 543 machine learning course.

    NOTE: the original examples on wh may be found at;
        https://www.tensorflow.org/hub/tutorials/boundless

        https://www.tensorflow.org/tutorials/generative/dcgan


  PYTHON VERSION:
    3.10.2

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>         - 2022/05/05:
        Written 
    
"""

import numpy as np
import sys
from glob import glob
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.gridspec as gridspec
import scipy
from scipy.interpolate import interp2d
from scipy import stats
import h5py
import subprocess
from scipy.stats import pearsonr,spearmanr
import os
from netCDF4 import Dataset
import time
import tensorflow as tf
from tensorflow.keras import layers

# Custom modules
home_dir = os.environ['HOME']
import objects
sys.path.append(home_dir)
from python_lib import circle, plot_subplot_label, plot_lat_circles
sys.path.append(home_dir + '/Research/TROPOMI')
from TROPOMI_Lib import *
sys.path.append(home_dir + '/Research/OMI')
from OMILib import readOMI_swath_hdf, write_swath_to_HDF5, download_OMI_all_HDF

sigma = [0.5, 0.5]
EPOCHS = 10
buffer_size = 60000   # size of training dataset
batch_size = 256

base_date = datetime(year = 2005, month = 1, day = 1, hour = 0, minute = 0)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Dataset generation functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# contrast_stretch performs contrast stretching on any data passed to the
# function.
def contrast_stretch(data):
    y1 = 0.
    y2 = 255.
    x1 = objects.min_AI
    x2 = objects.max_AI
    stretch_data = y1 + ((data - x1)/(x2 - x1)) * (y2 - y1)
    return stretch_data

# contrast_unstretch undoes the contrast stretching to get from 
# DCGAN-space to AI space 
def contrast_unstretch(data):
    y1 = objects.min_AI
    y2 = objects.max_AI
    x1 = 0.
    x2 = 255.
    unstretch_data = y1 + ((data - x1)/(x2 - x1)) * (y2 - y1)
    return unstretch_data

# Extracts a 300-pixel section of an OMI swath and
# averages every 5 pixels
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


def automate_TROP_input_process(date_str, minlat = 55, \
        save_dir = home_dir + '/CSCI/CSCI_544/final_project/prep_data/', \
        copy_to_talon = False, remove_trop_files = False, \
        remove_data_date_dir = False, \
        begin_date = None, end_date = None):
    
    # Figure out if date string is a year, month, or day
    # --------------------------------------------------
    if(len(date_str) == 4):
        str_fmt = '%Y'
        out_fmt = '%Ym'
        time_adder = 'year'
    elif(len(date_str) == 6):
        str_fmt = '%Y%m'
        out_fmt = '%Ym%m'
        time_adder = 'month'
    elif(len(date_str) == 8):
        str_fmt = '%Y%m%d'
        out_fmt = '%Ym%m%d'
        time_adder = 'day'
    else:
        print("INVALID DATE STRING. RETURNING")
        return


    dt_date_str = datetime.strptime(date_str, str_fmt)
    if(begin_date is None):
        begin_date = dt_date_str
    else:
        # Assume that the user wants a day as the beginner
        begin_date = datetime.strptime(begin_date, '%Y%m%d')

    if(end_date is None):
        if(time_adder == 'year'):
            end_date = dt_date_str + relativedelta(years = 1)
        elif(time_adder == 'month'):
            end_date = dt_date_str + relativedelta(months = 1)
            #end_date = dt_date_str + timedelta(month = 1)
        else:
            end_date = dt_date_str + timedelta(days = 1)
    else:
        end_date = datetime.strptime(end_date, '%Y%m%d')  + \
            timedelta(days = 1)

    local_date = begin_date

    while(local_date < end_date):

        print(local_date)
        # This function takes a single day str as input 
        download_OMI_all_HDF(local_date.strftime('%Y%m%d'))

        # For this date, grab the file names
        # ---------------------------------------------
        files = glob(local_date.strftime(\
            home_dir + '/data/OMI/H5_files/OMI-A*_%Ym%m%dt*.he5'))

        # The length of the returned list (times two) will
        # be the size of array needed for the netCDF file
        # ------------------------------------------------
        n_files = len(files)
        n_images = n_files * 2
        n_x = n_y = 60
        print(n_files)

        # Loop over all of the OMI files and process accordingly
        # ------------------------------------------------------
        for ii, fname in enumerate(files):
            dt_name = datetime.strptime(fname.strip().split('/')[-1][20:34],\
                '%Ym%m%dt%H%M')
            local_str = dt_name.strftime('%Y%m%d%H%M')
            print(local_str)

            flag_val = generate_paired_TROP_input(local_str, minlat = minlat, \
                                       remove_trop_files = remove_trop_files, 
                                       remove_data_date_dir = remove_data_date_dir)

        local_date = local_date + timedelta(days = 1)



# remove_data_date_dir: gets rid of the 201807052213 directory
# that contains the trop_prep and omi_shawn files, since these
# are both contained within the gzipped tar files.
def generate_paired_TROP_input(date_str, minlat = 55, \
        save_dir = home_dir + '/CSCI/CSCI_544/final_project/prep_data/', \
        copy_to_talon = False, remove_trop_files = False, \
        remove_data_date_dir = False):

    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')

    # Step 1: Read in the OMI data from the H5 file
    # ---------------------------------------------
    OMI_data = readOMI_swath_hdf(date_str, 'control', \
        latmin = minlat, skiprows = None)

    # Step 2: Call the resample function to extract just the data in the 
    #         top 600 valid swaths/rows. Need the following variables
    #   - LAT
    #   - LON
    #   - LATCRNR
    #   - LONCRNR
    #   - SZA
    #   - VZA
    #   - AZM
    #   - GPQF?
    # ------------------------------------------------------------------
    x_max = np.where(OMI_data['UVAI'][:,59].mask == False)[0][-1]
    x_min = x_max - 300
    #x_min = x_max - 6

    # Step 4: Read in the TROPOMI file and subset the data according
    #         to the latitudes given in the OMI subset grabbed above
    # ---------------------------------------------------------------
    trop_file_name = \
        generate_TROPOMI_prep_data(date_str, copy_to_raindrop = False, \
        minlat = minlat, remove_empty_scans = False, \
        trop_time = None, remove_large_files = remove_trop_files, \
        data_dir = save_dir, 
        generate_gzip_output = False, return_trop_name = True)

    if(isinstance(trop_file_name, str) is False):
        print("Error identified in generate_TROPOMI_prep_data")
        return -1

    # Extract the file path here
    # ---------------------------
    final_path = '/'.join(trop_file_name.strip().split('/')[:-1]) + '/'
    print(final_path)

    # Step 3: Write these extracted parameters to an HDF5 file, using
    #         the original variable names from the H5 file
    # ---------------------------------------------------------------
    write_swath_to_HDF5(OMI_data, 'control', save_path = final_path, minlat = minlat, \
        shawn_path = home_dir + '/data/OMI/shawn_files/', \
        scan_min = x_min, scan_max = x_max, \
        remove_empty_scans = False)

    # Step 5: Write the subsetted TROPOMI values to an output file
    # ------------------------------------------------------------

    # Step N: GZIP the subsetted OMI file and subsetted TROPOMI file
    # --------------------------------------------------------------
    os.chdir(save_dir)
    cmnd = dt_date_str.strftime('tar -cvzf ' + \
        'combined_prep_%Y%m%d%H%M.tar.gz ' + date_str + '/')
    print(cmnd)
    os.system(cmnd)

    if(copy_to_talon):
        # Secure copy the gzipped file to Raindrop
        # ----------------------------------------
        cmnd = dt_date_str.strftime('scp ' + \
            'combined_prep_%Y%m%d%H%M.tar.gz ' + \
            'blake.sorenson@134.129.128.241:' + \
            '~/CSCI_544/final_project/prep_data/')
        #cmnd = dt_date_str.strftime('scp ' + \
        #    'combined_subsets_%Y%m%d%H%M.tar.gz ' + \
        #    'bsorenson@raindrop.atmos.und.edu:' + \
        #    '/home/bsorenson/OMI/arctic_comp/comp_data/')
        print(cmnd)
        os.system(cmnd)

    # Delete the data directory if desired, since this data is
    # housed within the gzipped tar file
    # --------------------------------------------------------
    if(remove_data_date_dir):
        cmnd = 'rm -r ' + date_str + '/'
        print(cmnd)
        os.system(cmnd)

    return 0

# Reads associated OMI files for a given year and generates the OMI
# dataset.
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

    latmin = 5

    # Loop over the glob list
    # -----------------------
    for ii, fname in enumerate(files):
        dt_name = datetime.strptime(fname.strip().split('/')[-1][20:34],\
            '%Ym%m%dt%H%M')

        time_diff = (dt_name - base_date).total_seconds()

        # Read and resample each file
        # ---------------------------
        OMI_data = readOMI_swath_hdf(dt_name.strftime('%Y%m%d%H%M'), 'control', \
            latmin = latmin, \
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

# Reads associated OMI files for a given year and generates the OMI
# dataset.
def build_trop_training_dataset(date_str):

    # Figure out if date string is a year, month, or day
    # --------------------------------------------------
    if(len(date_str) == 4):
        str_fmt = '%Y'
        out_fmt = '%Y'
    elif(len(date_str) == 6):
        str_fmt = '%Y%m'
        out_fmt = '%Y%m'
    elif(len(date_str) == 8):
        str_fmt = '%Y%m%d'
        out_fmt = '%Y%m%d'
    else:
        print("INVALID DATE STRING. RETURNING")
        return

    dt_date_str = datetime.strptime(date_str, str_fmt)

    # For whatever date period is desired, use glob
    # to find the number of associated files
    # ---------------------------------------------
    files = glob(dt_date_str.strftime(\
        'coloc_data/colocated_tropomi_' + out_fmt + '*.hdf5'))

    # First, loop over all the data and figure out how many
    # swaths are acceptable
    # -----------------------------------------------------
    n_files = 0
    num_acceptable = 0

    print(files)

    # Loop over the glob list
    # -----------------------
    for ii, fname in enumerate(files):
        # Read data from each subsetted and colocated TROP file
        # -----------------------------------------------------
        data = h5py.File(fname)

        # Fix the missing regridded TROP pixels
        # -------------------------------------
        mask_trop = np.ma.masked_where(data['trop_ai'][:,:] == -999., data['trop_ai'][:,:])
        new_trop = mask_trop.reshape(-1, 5, mask_trop.shape[1]).mean(axis = 1)

        # Figure out how many rows contain missing values
        # ----------------------------------------------
        if(not np.ma.is_masked(new_trop)):
            num_acceptable += 1
        else:
            #rows_with_masked = np.array([True in new_trop.mask[jj,:] \
            #    for jj in range(new_trop.shape[0])])
    
            # For these rows, figure out how many values per row are missing
            # --------------------------------------------------------------
            num_masked_per_row = np.array([\
                len(np.where(new_trop.mask[jj,:] == True)[0]) \
                for jj in range(new_trop.shape[0])])

            print(fname, np.max(num_masked_per_row))

            if(np.max(num_masked_per_row) < 4): 
                num_acceptable += 1

        data.close()

    print('Total files:', len(files))
    print('Total acceptable:', num_acceptable)
    
    # The length of the returned list (times two) will
    # be the size of array needed for the netCDF file
    # ------------------------------------------------
    n_files = num_acceptable
    n_images = n_files * 2
    n_x = n_y = 60
    print('FILE DIMENSION:', n_files)

    # ----------------------
    # Set up the netCDF file
    # ----------------------
    file_name = 'train_dataset_trop_' + date_str + '.nc'
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
    images.description = 'Number of TROPOMI image subsets in the dataset'

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
    SZA = nc.createVariable('SZA','f4',('n_img','n_x','n_y'))
    SZA.description = 'Solar Zenith Angle'
    SZA.units = 'degrees'
    TIME = nc.createVariable('TIME','f8',('n_img'))
    TIME.description = 'Seconds since 00:00 UTC 01 January 2005'
    TIME.units = 'seconds'

    good_idx = 0

    # Loop over the glob list
    # -----------------------
    for ii, fname in enumerate(files):
        dt_name = datetime.strptime(fname.strip().split('/')[-1][18:30],\
            '%Y%m%d%H%M')

        time_diff = (dt_name - base_date).total_seconds()

        # Read data from each subsetted and colocated TROP file
        # -----------------------------------------------------
        data = h5py.File(fname)

        # Fix the missing regridded TROP pixels
        # -------------------------------------
        mask_trop = np.ma.masked_where(data['trop_ai'][:,:] == -999., data['trop_ai'][:,:])
        new_trop = mask_trop.reshape(-1, 5, mask_trop.shape[1]).mean(axis = 1)

        # Extract every 10th latitude and longitude
        # -----------------------------------------
        new_lat = data['omi_lat'][:,:]
        new_lat = new_lat.reshape(-1, 5, new_lat.shape[1]).mean(axis = 1)
        new_lon = data['omi_lon'][:,:]
        new_lon = new_lon.reshape(-1, 5, new_lon.shape[1]).mean(axis = 1)
        new_sza = data['omi_sza'][:,:]
        new_sza = new_sza.reshape(-1, 5, new_sza.shape[1]).mean(axis = 1)
       
        l_save_data = False
 
        # Figure out how many rows contain missing values
        # ----------------------------------------------
        if(not np.ma.is_masked(new_trop)):
            l_save_data = True
        else:
    
            # For these rows, figure out how many values per row are missing
            # --------------------------------------------------------------
            num_masked_per_row = np.array([\
                len(np.where(new_trop.mask[jj,:] == True)[0]) \
                for jj in range(new_trop.shape[0])])

            if(np.max(num_masked_per_row) < 4): 
                l_save_data = True

                for jj in range(new_trop.shape[0]):
                    where_mask = np.where(new_trop[jj,:].mask == True)[0]
                   
                    # See if this current line is masked
                    # ----------------------------------
                    if(len(where_mask) != 0):
 
                        # Check if the mask is on only one or both sides of
                        # the image
                        # -------------------------------------------------
                        min_mask = np.min(where_mask)
                        max_mask = np.max(where_mask)
    
                        if((min_mask < 10) & (max_mask < 10)):
                            # Mask is on only one side of the image. 
                            new_trop[jj,:][np.where(new_trop[jj,:].mask == True)] = \
                                new_trop[jj,:][np.where((new_trop[jj,:].mask == False))[0][0]]
                            
                        elif((min_mask > 50) & (max_mask > 50)):
                            # Mask is on the other side of the image. 
                            new_trop[jj,:][np.where(new_trop[jj,:].mask == True)] = \
                                new_trop[jj,:][np.where((new_trop[jj,:].mask == False))[0][-1]]
                        else:
                            # Mask is split


                            # Address the first side of the image first
                            # -----------------------------------------
                            new_trop[jj,:30][np.where(new_trop[jj,:30].mask == True)] = \
                                new_trop[jj,:30][np.where((new_trop[jj,:30].mask == False))[0][0]]

                            # Now, address the other side of the image
                            # ----------------------------------------
                            new_trop[jj,30:][np.where(new_trop[jj,30:].mask == True)] = \
                                new_trop[jj,30:][np.where((new_trop[jj,30:].mask == False))[0][-1]]


                    ##!#new_trop[jj,:][np.where(new_trop[jj,:].mask == True)] = \
                    ##!#    new_trop[jj,:][np.where((new_trop[jj,:].mask == False))[0][-1]]


        if(l_save_data):
            print(fname, good_idx)
    
            # Insert the data for the scene(s) from the current file into the
            # netCDF file
            # ---------------------------------------------------------------
            TIME[good_idx*2] = time_diff
            TIME[good_idx*2+1] = time_diff
            AI[(good_idx*2),:,:]    = new_trop
            LAT[(good_idx*2),:,:]   = new_lat[:,:]
            LON[(good_idx*2),:,:]   = new_lon[:,:]
            SZA[(good_idx*2),:,:]   = new_sza[:,:]
            AI[(good_idx*2)+1,:,:]  = new_trop[::-1,:]
            LAT[(good_idx*2)+1,:,:] = new_lat[::-1,:]
            LON[(good_idx*2)+1,:,:] = new_lon[::-1,:]
            SZA[(good_idx*2)+1,:,:] = new_sza[::-1,:]

            good_idx += 1

        data.close()


    # Close the netCDF file, which actually saves it.
    # -----------------------------------------------
    nc.close()

    print('netCDF file saved')






# Reads the three OMI datasets for 2005, 2006, and 2007, combines the
# datasets, removes images with bad data, and selects only data from
# a single month if desired. 
def read_train_dataset(dat_name = 'default', month = None, low_res = False):
    res_add = ''
    if(dat_name == 'OMI'):
        # Open all 3 training datasets
        # ----------------------------
        data1 = Dataset('./train_dataset_2005.nc','r')
        data2 = Dataset('./train_dataset_2006.nc','r')
        data3 = Dataset('./train_dataset_2007.nc','r')

        begin_size = data1['AI'].shape[0] + data2['AI'].shape[0] + \
            data3['AI'].shape[0] 
        print("Before thinning: ", begin_size)

        # Remove any swaths with any missing pixels
        # -----------------------------------------
        masks1 = np.array([(True in np.array(data1['AI'][ii,:,:].mask)) \
            for ii in range(np.array(data1['AI']).shape[0])])
        masks2 = np.array([(True in np.array(data2['AI'][ii,:,:].mask)) \
            for ii in range(np.array(data2['AI']).shape[0])])
        masks3 = np.array([(True in np.array(data3['AI'][ii,:,:].mask)) \
            for ii in range(np.array(data3['AI']).shape[0])])
      
        testdata1 = data1['AI'][np.where(masks1 == False)]
        testdata2 = data2['AI'][np.where(masks2 == False)]
        testdata3 = data3['AI'][np.where(masks3 == False)]
        testLAT1  = data1['LAT'][np.where(masks1 == False)]
        testLAT2  = data2['LAT'][np.where(masks2 == False)]
        testLAT3  = data3['LAT'][np.where(masks3 == False)]
        testLON1  = data1['LON'][np.where(masks1 == False)]
        testLON2  = data2['LON'][np.where(masks2 == False)]
        testLON3  = data3['LON'][np.where(masks3 == False)]
        testTIME1 = data1['TIME'][np.where(masks1 == False)]
        testTIME2 = data2['TIME'][np.where(masks2 == False)]
        testTIME3 = data3['TIME'][np.where(masks3 == False)]

        # Close the netCDF objects.
        # -------------------------
        data1.close()
        data2.close()
        data3.close()
   
        total_size = testdata1.shape[0] + testdata2.shape[0] + \
            testdata3.shape[0]

        print("After removing unknown bad rows: ", total_size)

        # Insert the masked data from the three individual datasets into
        # single arrays.
        # --------------------------------------------------------------
        total_data = np.full((total_size, 60, 60), np.nan)
        total_LAT  = np.full((total_size, 60, 60), np.nan)
        total_LON  = np.full((total_size, 60, 60), np.nan)
        total_TIME = np.full((total_size), np.nan)
 
        total_data[:testdata1.shape[0], :, :] = testdata1
        total_data[testdata1.shape[0]:testdata1.shape[0] + \
            testdata2.shape[0], :, :] = testdata2
        total_data[testdata1.shape[0] + testdata2.shape[0] : \
            testdata1.shape[0] + testdata2.shape[0] + \
            testdata3.shape[0], :, :] = testdata3

        total_LAT[:testdata1.shape[0], :, :] = testLAT1
        total_LAT[testdata1.shape[0]:testdata1.shape[0] + \
            testdata2.shape[0], :, :] = testLAT2
        total_LAT[testdata1.shape[0] + testdata2.shape[0] : \
            testdata1.shape[0] + testdata2.shape[0] + \
            testdata3.shape[0], :, :] = testLAT3

        total_LON[:testdata1.shape[0], :, :] = testLON1
        total_LON[testdata1.shape[0]:testdata1.shape[0] + \
            testdata2.shape[0], :, :] = testLON2
        total_LON[testdata1.shape[0] + testdata2.shape[0] : \
            testdata1.shape[0] + testdata2.shape[0] + \
            testdata3.shape[0], :, :] = testLON3

        total_TIME[:testdata1.shape[0]] = testTIME1
        total_TIME[testdata1.shape[0]:testdata1.shape[0] + \
            testdata2.shape[0]] = testTIME2
        total_TIME[testdata1.shape[0] + testdata2.shape[0] : \
            testdata1.shape[0] + testdata2.shape[0] + \
            testdata3.shape[0]] = testTIME3

        # Only use one month?
        # -------------------
        if(month is not None):
            print("months")

            # First, remove all 2nd swath sections
            # ------------------------------------
            total_data = total_data[::2]
            total_LAT  = total_LAT[::2]
            total_LON  = total_LON[::2]
            total_TIME = total_TIME[::2]

            # Then, use the times to remove any first swaths that
            # are not from the desired month
            # ---------------------------------------------------
            local_times = np.array([base_date + \
                timedelta(seconds = ttime) for ttime in total_TIME[:]])        
            months = np.array([ltime.month for ltime in local_times]) 
            if(isinstance(month, list)):    
                if(len(month) == 2):
                    indices = np.where((months == month[0]) | \
                        (months == month[1]))
                elif(len(month) == 3):
                    indices = np.where((months == month[0]) | \
                                       (months == month[1]) | \
                                       (months == month[2]))
            else:
                indices = np.where(months == month)

            total_data = total_data[indices]
            total_LAT  = total_LAT[indices]
            total_LON  = total_LON[indices]
            total_TIME = total_TIME[indices]

            print("After month screening: ", total_data.shape[0])

        # Clip the data to remove the extra low AI data?
        # ----------------------------------------------
        total_data = np.clip(total_data, -6.0, 10.0)
   
        if(low_res):
            print("Scaling back the data")
            res_add = '_lowres'
            temp_data = np.full((total_data.shape[0], 28, 28), np.nan)
            temp_lat  = np.full((total_data.shape[0], 28, 28), np.nan)
            temp_lon  = np.full((total_data.shape[0], 28, 28), np.nan)

            x_dim = 60 / 28
            x_range = np.arange(60)
            interp_range = np.arange(0, 60, x_dim)
    
            # Loop over each variable, interpolating 
            for ii in range(total_data.shape[0]):
                interp_data = interp2d(x_range, x_range, total_data[ii,:,:]) 
                interp_lat  = interp2d(x_range, x_range, total_LAT[ii,:,:]) 
                interp_lon  = interp2d(x_range, x_range, total_LON[ii,:,:]) 

                temp_data[ii,:,:] = interp_data(interp_range, interp_range)
                temp_lat[ii,:,:]  = interp_lat(interp_range, interp_range)
                temp_lon[ii,:,:]  = interp_lon(interp_range, interp_range)

            total_data = temp_data
            total_LAT = temp_lat
            total_LON = temp_lon
 
        # Scale the AI data from -1 to 1
        # ------------------------------
        objects.max_AI = np.max(total_data)
        objects.min_AI = np.min(total_data)
        total_data = contrast_stretch(total_data[:,:,:])
        total_data = total_data.reshape(total_data.shape[0], \
            total_data.shape[1], total_data.shape[2], 1).astype('float32')
        objects.train_images = (total_data - 127.5) / 127.5

        objects.train_lat  = total_LAT
        objects.train_lon  = total_LON
        objects.train_time = total_TIME

    else:
        print('reading training images')
        
        (train_images, train_labels), (_, _) = tf.keras.datasets.mnist.load_data()
        
        train_images = train_images.reshape(train_images.shape[0], 28, 28, 1).astype('float32')
        objects.train_images = (train_images - 127.5) / 127.5  # Normalize the images to [-1, 1]
        
   
    # Make figure for the data distribution
    # -------------------------------------
    plt.close('all')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.hist(objects.train_images.flatten(), bins = 100)
    ax.set_yscale('log')
    ax.set_ylabel('Counts')
    outname = 'train_data_dist_'+objects.dat_name+res_add+'.png'
    fig.savefig(outname)
    print("Saved image", outname)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Plotting functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Plot some OMI AI data from the training dataset
# Set up to work with the netCDF file object from build_training_dataset
def plot_train_data(OMI_data, idx, plot_map = False, save = False):
   
    local_time = base_date + timedelta(seconds = OMI_data['TIME'][idx]/1) 
   
    plt.close()
    fig = plt.figure()
    if(plot_map):

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

def generate_and_save_images(model, epoch, test_input):
  # Notice `training` is set to False.
  # This is so all layers run in inference mode (batchnorm).
  predictions = model(test_input, training=False)

  fig = plt.figure(figsize=(4, 4))

  for i in range(predictions.shape[0]):
      plt.subplot(4, 4, i+1)
      plt.imshow(predictions[i, :, :, 0] * 127.5 + 127.5, cmap='gray')
      plt.axis('off')

  plt.savefig('image_at_epoch_{:04d}_original.png'.format(epoch))

# Saves a generated image (not complete, just generator output)
def save_gen_image(model, epoch, test_input):

  # Notice `training` is set to False.
  # This is so all layers run in inference mode (batchnorm).
  predictions = model(test_input, training=False)

  fig = plt.figure(figsize=(4, 4))
  ax1 = fig.add_subplot(1,1,1)
  ax1.imshow(predictions[0,:,:,0] * 127.5 + 127.5, cmap = 'gray')
  ax1.axis('off')

  outname = 'gen_image_epoch_{:04d}_original.png'.format(epoch)
  plt.savefig(outname)
  print("Saved image", outname)

# Saves an image from the training dataset
def save_training_image(dataset, idx, test = False, mask = None):

  if(mask is not None):
    plot_data = (dataset[idx,:,:,0] * 127.5 + 127.5)
    plot_data = np.ma.masked_where(mask == 0, plot_data)
  else:
    plot_data = (dataset[idx,:,:,0] * 127.5 + 127.5)

  fig = plt.figure(figsize=(4, 4))
  ax1 = fig.add_subplot(1,1,1)
  ax1.imshow(plot_data, cmap = 'gray')
  ax1.axis('off')

  if(test):
    adder = 'test'
  else:
    adder = 'train'
  if(mask is not None):
    outname = adder + '_image_{:05d}_mask.png'.format(idx)
    plt.savefig(outname)
    print("Saved image", outname)
  else:
    outname = adder + '_image_{:05d}.png'.format(idx)
    plt.savefig(outname)
    print("Saved image", outname)

# Generates a figure for a DCGAN-completed image. Contains both the 
# completed image and the generator output
def save_complete_image(iter_num = 1, idx = 3, smooth = False, save = False):

    # Make completed image
    # --------------------
    complete_img = objects.complete_img
    if(len(complete_img.shape) > 2):
        complete_img = complete_img[idx,:,:]
        gen_img = objects.gen_img[idx,:,:,0]
        out_add = '_im'+str(idx)   
    else:
        out_add = ''

    complete_img = complete_img[:,:]
    gen_img = objects.gen_img[0,:,:,0]

    smooth_add = ''
    if(objects.dat_name == 'OMI' and smooth):
        smooth_add = '_smooth'
        gen_img = \
            scipy.ndimage.filters.gaussian_filter(\
            gen_img, sigma, mode='constant')

    fig = plt.figure(figsize = (8,4))
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
    ax1.imshow(complete_img, cmap = 'gray')
    ax2.imshow(gen_img, cmap = 'gray', vmin = np.min(complete_img), \
        vmax = np.max(complete_img))
    ax1.axis('off')
    ax2.axis('off')

    if(save):
        
        outname = 'complete_image_test_v{:04d}'.format(iter_num) + \
            out_add + smooth_add + '.png'
        fig.savefig(outname)
        print("Saved image", outname)

# Generates a figure showing the completion losses while completing
# an image
def save_complete_losses(idx, losses):

    if(losses.shape[0] > 1):
        losses = losses[0,:]

    plt.close('all')
    fig = plt.figure(figsize = (4,4))
    ax1 = fig.add_subplot(1,1,1)
    ax1.plot(losses)
    ax1.set_title('Completion of image ' + str(idx))
    ax1.set_xlabel('Iterations')
    ax1.set_ylabel('Complete loss')
    fig.tight_layout()
    outname = 'complete_losses_'+str(idx) + '.png'
    fig.savefig(outname)
    print("Saved image", outname)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# TensorFlow Model Functions.
#
# Modified from the TensorFlow DCGAN examples
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# This method returns a helper function to compute cross entropy loss
cross_entropy = tf.keras.losses.BinaryCrossentropy(from_logits=True)

def make_generator_model():
    model = tf.keras.Sequential()
    model.add(layers.Dense(int(objects.image_size / 4) * \
        int(objects.image_size / 4) * int(objects.batch_size), \
        use_bias=False, input_shape=(100,)))
    #model.add(layers.Dense(7*7*256, use_bias=False, input_shape=(100,)))
    model.add(layers.BatchNormalization())
    model.add(layers.LeakyReLU())

    #model.add(layers.Reshape((7, 7, 256)))
    #assert model.output_shape == (None, 7, 7, 256)  # Note: None is the batch size
    model.add(layers.Reshape((int(objects.image_size / 4), \
        int(objects.image_size / 4), int(objects.batch_size))))
    assert model.output_shape == (None, int(objects.image_size / 4), \
        int(objects.image_size / 4), int(objects.batch_size))  # Note: None is the batch size

    #model.add(layers.Conv2DTranspose(128, (5, 5), strides=(1, 1), padding='same', use_bias=False))
    #assert model.output_shape == (None, 7, 7, 128)
    model.add(layers.Conv2DTranspose(int(objects.batch_size / 2), (5, 5), \
        strides=(1, 1), padding='same', use_bias=False))
    assert model.output_shape == (None, int(objects.image_size / 4), \
        int(objects.image_size / 4), int(objects.batch_size / 2))
    model.add(layers.BatchNormalization())
    model.add(layers.LeakyReLU())

    #model.add(layers.Conv2DTranspose(64, (5, 5), strides=(2, 2), padding='same', use_bias=False))
    #assert model.output_shape == (None, 14, 14, 64)
    model.add(layers.Conv2DTranspose(int(objects.batch_size / 4), (5, 5), \
        strides=(2, 2), padding='same', use_bias=False))
    assert model.output_shape == (None, int(objects.image_size / 2), \
        int(objects.image_size / 2), int(objects.batch_size / 4))
    model.add(layers.BatchNormalization())
    model.add(layers.LeakyReLU())

    model.add(layers.Conv2DTranspose(1, (5, 5), strides=(2, 2), \
        padding='same', use_bias=False, activation='tanh'))
    #assert model.output_shape == (None, 28, 28, 1)
    assert model.output_shape == (None, int(objects.image_size), \
        int(objects.image_size), 1)

    return model

def make_discriminator_model():
    model = tf.keras.Sequential()
    #model.add(layers.Conv2D(64, (5, 5), strides=(2, 2), padding='same',
    #                                 input_shape=[28, 28, 1]))
    model.add(layers.Conv2D(int(objects.batch_size / 4), (5, 5), \
        strides=(2, 2), padding='same',
        input_shape=[int(objects.image_size), \
        int(objects.image_size), 1]))
    model.add(layers.LeakyReLU())
    model.add(layers.Dropout(0.3))

    #model.add(layers.Conv2D(128, (5, 5), strides=(2, 2), padding='same'))
    model.add(layers.Conv2D(int(objects.batch_size / 2), \
        (5, 5), strides=(2, 2), padding='same'))
    model.add(layers.LeakyReLU())
    model.add(layers.Dropout(0.3))

    model.add(layers.Flatten())
    model.add(layers.Dense(1))

    return model

def discriminator_loss(real_output, fake_output):
    real_loss = cross_entropy(tf.ones_like(real_output), real_output)
    fake_loss = cross_entropy(tf.zeros_like(fake_output), fake_output)
    total_loss = real_loss + fake_loss
    return total_loss

def generator_loss(fake_output):
    return cross_entropy(tf.ones_like(fake_output), fake_output)

# Notice the use of `tf.function`
# This annotation causes the function to be "compiled".
@tf.function
def train_step(images):
    noise = tf.random.normal([objects.batch_size, objects.noise_dim])

    with tf.GradientTape() as gen_tape, tf.GradientTape() as disc_tape:
        
        generated_images = objects.generator(noise, training=True)

        real_output = objects.discriminator(images, training=True)
        fake_output = objects.discriminator(generated_images, training=True)

        gen_loss = generator_loss(fake_output)
        disc_loss = discriminator_loss(real_output, fake_output)

    gradients_of_generator = gen_tape.gradient(gen_loss, \
        objects.generator.trainable_variables)
    gradients_of_discriminator = disc_tape.gradient(disc_loss, \
        objects.discriminator.trainable_variables)

    objects.generator_optimizer.apply_gradients(zip(\
        gradients_of_generator, objects.generator.trainable_variables))
    objects.discriminator_optimizer.apply_gradients(zip(\
        gradients_of_discriminator, objects.discriminator.trainable_variables))

def train(dataset, epochs):
    for epoch in range(epochs):
        start = time.time()

        for image_batch in dataset:
          train_step(image_batch)
      
        # Produce images for the GIF as you go
        plt.close('all')
        generate_and_save_images(objects.generator,
                                 epoch + 1,
                                 objects.seed)
      
        # Save the model every 15 epochs
        if (epoch + 1) % 15 == 0:
          objects.checkpoint.save(file_prefix = objects.checkpoint_prefix)
      
        print ('Time for epoch {} is {} sec'.format(epoch + 1, \
            time.time()-start))
    
    # Generate after the final epoch
    generate_and_save_images(objects.generator,
                             epochs,
                             objects.seed)

    print('done training')

# Trains the model for a certain number of epochs
def train_model(EPOCHS = 10):

    objects.epochs = EPOCHS

    print('making generator')

    objects.generator = make_generator_model()
    objects.discriminator = make_discriminator_model()
    
    objects.generator_optimizer = \
        tf.keras.optimizers.Adam(objects.learn_rate)
    objects.discriminator_optimizer = \
        tf.keras.optimizers.Adam(objects.learn_rate)
    
    print('setting checkpoints')
    checkpoint_dir = './training_checkpoints'
    objects.checkpoint_prefix = os.path.join(checkpoint_dir, "ckpt")
    objects.checkpoint = tf.train.Checkpoint(\
            generator_optimizer=objects.generator_optimizer,
            discriminator_optimizer=objects.discriminator_optimizer,
            generator=objects.generator,
            discriminator=objects.discriminator)
    
    # Batch and shuffle the data
    train_dataset = tf.data.Dataset.from_tensor_slices(\
        objects.train_images).shuffle(objects.buffer_size).batch(\
        objects.batch_size)

    # You will reuse this seed overtime (so it's easier)
    # to visualize progress in the animated GIF)
    objects.seed = tf.random.normal([objects.num_examples_to_generate, \
        objects.noise_dim])
    
    train(train_dataset, EPOCHS)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# TensorFlow image completion functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def run_complete(plot_losses = False, smooth = False):
    lam = tf.constant(0.1)

    mask = tf.constant(objects.np_mask, dtype = tf.float32)

    # Set up starting seeds 
    # ---------------------
    if(len(objects.process_img.shape) > 3):
        zhats = tf.Variable(np.random.uniform(-1, 1, \
            size=(objects.process_img.shape[0], objects.noise_dim)))
    else:
        zhats = tf.Variable(np.random.uniform(-1, 1, \
            size=(1, objects.noise_dim)))
  
    print('zhats = ', zhats.shape)
 
    iter_nums = 2000
    v_num = 0.0
    losses = np.full((objects.process_img.shape[0], iter_nums), np.nan)
    num_loss = 1e5
    ii = 0
    for ii in range(iter_nums):
        
        with tf.GradientTape() as gen_tape:
            gen_tape.watch(zhats)
        
            # Set up loss calculations
            # ------------------------
            complete_loss = tf.reduce_sum(tf.compat.v1.layers.flatten(\
                tf.abs(tf.multiply(mask[np.newaxis,...,np.newaxis], \
                objects.generator(zhats, training = False)) - tf.multiply(\
                mask[np.newaxis,...,np.newaxis], \
                tf.Variable(objects.process_img))))) + \
                lam * generator_loss(\
                objects.discriminator(objects.generator(\
                    zhats, training = False), training = False))
   
        # Calculate the loss gradients with respect to the seed
        # ----------------------------------------------------- 
        grad_complete_loss = gen_tape.gradient(complete_loss, zhats)
        num_loss =  complete_loss.numpy()
        losses[:, ii] = num_loss
   
        # Modify the seeds
        # ---------------- 
        v_prev = np.copy(v_num)
        if(zhats.shape[0] > 1):
            v_num = objects.beta1 * v_prev - \
                objects.complete_learn_rate * grad_complete_loss
            zhats = zhats.assign_add(\
                -objects.beta1 * v_prev + (1.0 + objects.beta1) * v_num)
        else:
            v_num = objects.beta1 * v_prev - objects.complete_learn_rate * \
                grad_complete_loss[0]
            zhats = zhats.assign_add(\
                tf.expand_dims(-objects.beta1 * v_prev + \
                    (1.0 + objects.beta1) * v_num, axis = 0))

        zhats = tf.Variable(np.clip(zhats, -1, 1))
   
        # Save a completed image every 100 iterations
        # ------------------------------------------- 
        if(ii % 100 == 0):
            zhat_gen_img = objects.generator(zhats, training = False) 
            objects.gen_img = zhat_gen_img
            objects.complete_img = np.squeeze(objects.process_img) * \
                objects.np_mask + (1. - objects.np_mask) * \
                np.squeeze(zhat_gen_img)
            plt.close('all')
            save_complete_image(iter_num = ii, idx = 1, smooth = smooth, \
                save = True)
        ii += 1

    # Put the generated image in objects.
    # -----------------------------------
    zhat_gen_img = objects.generator(zhats, training = False) 
    objects.gen_img = zhat_gen_img

    # Put the completed image in objects
    # ----------------------------------
    objects.complete_img = np.squeeze(objects.process_img) * \
        objects.np_mask + (1. - objects.np_mask) * np.squeeze(zhat_gen_img)

    # Save a figure for the losses
    # ----------------------------
    objects.losses = losses
    if(plot_losses):
        save_complete_losses(idx, losses)

def complete_image(idx, test = False, plot_losses = True, smooth = False):

    # Prep the indices. Use the index/indices to extract the
    # images for completion, whether from testing dataset or training
    # ---------------------------------------------------------------
    if(not isinstance(idx, np.ndarray)): 
        idx = np.array([idx])
   
    if(test):
        objects.process_img = objects.test_images[idx,:,:,:]
        save_training_image(objects.test_images, idx[0], test = test, mask = None)
    else:     
        objects.process_img = objects.train_images[idx,:,:,:]
        save_training_image(objects.train_images, idx[0], test = test, mask = None)
    
    print(objects.process_img.shape)

    run_complete(plot_losses = plot_losses, smooth = smooth)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Validation functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Calculates errors between the true AI and DCGAN AI, as well as between
# true AI and noise AI, if desired.
def calc_validate_stats(idx, noise = False):

    # Extract the mask from the np_mask
    # ---------------------------------
    mask_mask = np.ma.masked_where(objects.np_mask == 0, objects.np_mask)

    # Mask the original images
    # ------------------------
    mask_test = np.array([ob_tst[mask_mask.mask] for ob_tst in \
        objects.test_images[idx,:,:,0]])

    if(noise):
        # Make normally-distributed noise in the mask
        # -------------------------------------------
        mask_complete = np.array([np.random.normal(np.mean(mt), \
            np.std(mt), mt.shape[0]) for mt in mask_test])
    else:
        # Extract the completed AI within the mask region 
        # -----------------------------------------------
        mask_complete = np.array([ob_im[mask_mask.mask] for ob_im in \
            objects.complete_img])

    # Calculate errors between the true AI and the DCGAN/noise 
    # --------------------------------------------------------
    objects.errors = np.nanmean(((mask_complete - mask_test) / mask_test), \
        axis = 1) * 100.
    

# Complete - whether or not to process the completion code. Set to False if images
# have already been generated.
def plot_validate_stats(idx, complete = True, noise = False, test = True, \
        save = True):
    
    # Add mask ranges here?
    # ---------------------
    beg_masks = np.array([14, 13, 13, 12, 12, 11])
    end_masks = np.array([15, 15, 16, 16, 17, 17])
    
    if(noise):
        noise_add = '_noise'
        means   = np.full((2, len(beg_masks)), np.nan)
        medians = np.full((2, len(beg_masks)), np.nan)
        stds    = np.full((2, len(beg_masks)), np.nan)
    else:
        noise_add = ''
        means   = np.full((1, len(beg_masks)), np.nan)
        medians = np.full((1, len(beg_masks)), np.nan)
        stds    = np.full((1, len(beg_masks)), np.nan)

    range_num = 50

    objects.combined_errors = np.full((len(idx), len(beg_masks)), np.nan)

    plt.close('all')
    fig = plt.figure(figsize = (9,7))
    for ii in range(len(beg_masks)):
        print('Iteration',ii)

        # Set up panel
        # ------------
        ax = fig.add_subplot(2,3,ii+1)

        # Set up current mask
        # -------------------
        np_mask = np.ones((objects.image_size,objects.image_size))
        np_mask[:, beg_masks[ii]:end_masks[ii]] = 0.0
        objects.np_mask = np_mask

        # Calculate statistics for generated images
        # -----------------------------------------
        if(complete):
            if(not noise):
                # Use complete_image to generate the completed images
                # ---------------------------------------------------
                print("Completing images")
                complete_image(idx, test = test, plot_losses = False)

        calc_validate_stats(idx, noise = False)
        objects.combined_errors[:, ii] = objects.errors

        # Insert DCGAN error stats into arrays
        # ------------------------------------
        means[0,ii]   = np.mean(objects.errors)
        medians[0,ii] = np.median(objects.errors)
        stds[0,ii]    = np.std(objects.errors)

        print('  Generator')
        print('    mean error', means[0,ii])
        print('    med  error', medians[0,ii])
        print('    std  error', stds[0,ii])
        

        # Plot the current generated statistics in panel
        # ----------------------------------------------
        ax.hist(objects.combined_errors[:,ii], bins = 50, \
            range = (-range_num, range_num), alpha = 0.5, \
            label = 'DCGAN') 

        # Calculate statistics for random noise, if desired
        # -------------------------------------------------
        if(noise):       
            calc_validate_stats(idx, noise = noise)

            # Insert noise error stats into arrays
            # ------------------------------------
            means[1,ii]   = np.mean(objects.errors)
            medians[1,ii] = np.median(objects.errors)
            stds[1,ii]    = np.std(objects.errors)

            print('  Noise')
            print('    mean error', means[1,ii])
            print('    med  error', medians[1,ii]) 
            print('    std  error', stds[1,ii])

            # Plot the current generated statistics in panel
            # ----------------------------------------------
            ax.hist(objects.errors, bins = 50, \
                range = (-range_num, range_num), alpha = 0.5, \
                label = 'Noise') 

        ax.set_xlim(-range_num, range_num)
        ax.axvline(0, color = 'black', linestyle = ':')
        ax.set_title('Mask size = '+str(int(end_masks[ii] - \
            beg_masks[ii])))
        ax.set_xlabel('Average percent error within mask')
        ax.set_ylabel('Counts')
        ax.legend()

    # Save the error statistics to a CSV file
    # ---------------------------------------
    np.savetxt('./mean_error_stats' + objects.dat_name + '_epoch' + \
        str(int(objects.epochs)) + noise_add + '.csv', means, \
        delimiter = ',')
    np.savetxt('./std_ error_stats' + objects.dat_name + '_epoch' + \
        str(int(objects.epochs)) + noise_add + '.csv', stds, \
        delimiter = ',')
    np.savetxt('./median_error_stats' + objects.dat_name + '_epoch' + \
        str(int(objects.epochs)) + noise_add + '.csv', medians, \
        delimiter = ',')
    print("Saved error statistics")

    fig.tight_layout()
    if(save):
        outname = 'combined_stats_' + objects.dat_name + '_epoch' + \
            str(int(objects.epochs)) + noise_add + '.png'
        fig.savefig(outname)
        print("Saved image", outname) 
    else:
        plt.close('all')
    
def validate_image(idx, test = True, plot_map = False, save = False):

    # Make the generated image here
    # -----------------------------
    complete_image(idx, test = test, plot_losses = True)

    # Unscale the data to get it back into AI space
    # ---------------------------------------------
    print("UNSCALE HERE")

    # Calculate percent error within the masked region between the original
    # and generated images
    # ---------------------------------------------------------------------
    if(test):
        working_img = objects.test_images[idx,:,:,0]
        adder = 'test'
    else:
        working_img = objects.train_images[idx,:,:,0]
        adder = 'train'

    mask_img = np.ma.masked_where(objects.np_mask == 0.0, working_img)
    base_img  = working_img * (1. - objects.np_mask)
    gen_img   = objects.complete_img * (1. - objects.np_mask)
    avg_error = np.nanmean((gen_img - base_img) / base_img) * 100.

    print(idx, avg_error)

    if(save):
        # Set up figure. Four panels:
        #   - Original panel
        #   - Original panel with mask
        #   - Generated image
        #   - Loss time series?
        # ----------------------------
        fig = plt.figure()
        if(plot_map):
            if(test):
                lons = objects.test_lon[idx,:,:]
                lats = objects.test_lat[idx,:,:]
            else:
                lons = objects.train_lon[idx,:,:]
                lats = objects.train_lat[idx,:,:]
            mapcrs  = ccrs.NorthPolarStereo()
            datacrs = ccrs.PlateCarree()
 
            avg_lat = np.nanmean(objects.train_lat[idx,:,:])

            if(avg_lat > 70):
                mapcrs = ccrs.NorthPolarStereo()
            else:
                mapcrs = ccrs.LambertConformal(central_latitude = avg_lat,\
                    central_longitude = \
                    np.nanmean(objects.train_lon[idx,:,:]))

            ax0 = fig.add_subplot(1,3,1, projection = mapcrs)
            ax1 = fig.add_subplot(1,3,2, projection = mapcrs)
            ax2 = fig.add_subplot(1,3,3, projection = mapcrs)
            
            ax0.coastlines()
            ax1.coastlines()
            ax2.coastlines()
            ax0.add_feature(cfeature.BORDERS)
            ax1.add_feature(cfeature.BORDERS)
            ax2.add_feature(cfeature.BORDERS)

            mesh0 = ax0.pcolormesh(lons,lats, working_img, \
                transform = ccrs.PlateCarree())
            mesh1 = ax1.pcolormesh(lons,lats, mask_img,  
                transform = ccrs.PlateCarree())
            mesh2 = ax2.pcolormesh(lons,lats, objects.complete_img,  \
                transform = ccrs.PlateCarree())
        else:
            ax0 = fig.add_subplot(1,3,1)
            ax1 = fig.add_subplot(1,3,2)
            ax2 = fig.add_subplot(1,3,3)

            ax0.imshow(working_img, cmap = 'gray')
            ax1.imshow(mask_img, cmap = 'gray')
            ax2.imshow(objects.complete_img, cmap = 'gray')
        
            ax0.axis('off')
            ax1.axis('off')
            ax2.axis('off')

        ax0.set_title('Original')
        ax1.set_title('Masked')
        ax2.set_title('Generated')
        
        fig.tight_layout()  
       
        if(objects.dat_name == 'OMI'): 
            adder = +'_min'+\
            str(int(objects.min_AI)) + '_max'+str(int(objects.max_AI))
        else:
            adder = '_numbers'
        outname = 'validate_' + str(idx) + '_e'+str(objects.epochs) + \
            adder + '.png'
        fig.savefig(outname)
        print("Saved image", outname)

def plot_validate_multiple(idx, complete = True, test = True, smooth = False, \
        save = True):

    if(len(idx) > 5):
        test_idxs = random.sample(range(0,len(idx)) , 5)

    if(complete):
        # Make the generated image here
        # -----------------------------
        complete_image(idx, test = test, smooth = smooth, plot_losses = True)
  
    # Set up figure
    # -------------
    plt.close('all')
    fig = plt.figure(figsize = (9,9)) 
    axs = fig.subplots(len(idx), 4)
   
    for ii in range(len(idx)):
        mask_img = np.ma.masked_where(objects.np_mask == 0.0, \
            objects.test_images[idx[ii],:,:,0])
        # panel 1 - original
        # panel 2 - masked
        # panel 3 - generated
        # panel 4 - losses
        smooth_add = ''
        if(objects.dat_name == 'OMI' and smooth):
            smooth_add = '_smooth'
            gen_img = \
                scipy.ndimage.filters.gaussian_filter(\
                objects.complete_img[ii,:,:], sigma, mode='constant')
        else:
            gen_img = objects.complete_img[ii,:,:]

        axs[ii,0].imshow(objects.test_images[idx[ii],:,:,0], cmap = 'gray')
        axs[ii,1].imshow(mask_img, cmap = 'gray')
        axs[ii,2].imshow(gen_img, cmap = 'gray')
        
        axs[ii,0].set_title('Original ' + str(idx[ii]))
        axs[ii,1].set_title('Masked ' + str(idx[ii]))
        axs[ii,2].set_title('Completed ' + str(idx[ii]))

        axs[ii,0].axis('off')
        axs[ii,1].axis('off')
        axs[ii,2].axis('off')

        axs[ii,3].plot(objects.losses[ii])
        axs[ii,3].set_xlabel('Iterations')
        axs[ii,3].set_ylabel('Complete loss')

    fig.subplots_adjust(wspace = None, hspace = None)
    fig.tight_layout()

    if(save):
        outname = 'combined_validate_' + objects.dat_name + '_epoch' + \
            str(int(objects.epochs)) + '.png'
        fig.savefig(outname)
        print("Saved image", outname)
    else:
        print('not saving')
       
