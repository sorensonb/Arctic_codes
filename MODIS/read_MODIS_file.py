#!/usr/bin/env python

"""
  NAME:
    read_MODIS_file.py

  PURPOSE:
    Show examples of different ways to read MODIS files using either
    pyhdf or netCDF4
  
  SYNTAX:
    $ python read_MODIS_file.py    

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2023/07/18:
      Written

"""

import numpy as np
from netCDF4 import Dataset
from pyhdf import SD
import matplotlib.pyplot as plt

test_filename = 'MYD04_L2.A2018001.0000.061.2018003204016.hdf'

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Method 1: Read the file using netCDF4     *RECOMMENDED*
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

##!## Open the file using netCDF4
##!## ---------------------------
##!#data1 = Dataset(test_filename)
##!#
##!## If desired, you can print the names of all the variables in the file using
##!## the following line:
##!##print(data1.variables.keys())
##!#
##!## Access the latitude and longitude data
##!## --------------------------------------
##!#lat = data1['Latitude'][:,:]
##!#lon = data1['Longitude'][:,:]
##!#
##!## Access one of the AOD parameters
##!## --------------------------------
##!#aod_db = data1['AOD_550_Dark_Target_Deep_Blue_Combined'][:,:]
##!#
##!## NOTE: netCDF4 automatically accounts for the scale factor and offset values
##!##       when reading an attribute from a file. Thus, no further action
##!##       is needed to use the AOD data.
##!#
##!## Close the file
##!## --------------
##!#data1.close()

##!## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
##!##
##!## Method 2: Read the file using pyhdf
##!##
##!## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
##!#
##!#data2 = SD.SD(test_filename)
##!#
##!## Access the latitude and longitude data
##!## --------------------------------------
##!#lat = data2.select('Latitude').get()
##!#lon = data2.select('Longitude').get()
##!#
##!## Access one of the AOD parameters
##!## --------------------------------
##!#aod_db_hdf = data2.select('AOD_550_Dark_Target_Deep_Blue_Combined').get()
##!#
##!## Remove any missing values. Missing values are defined using the 'fill value'
##!## parameter.
##!## ----------------------------------------------------------------------------
##!#fill_value = data2.select(\
##!#    'AOD_550_Dark_Target_Deep_Blue_Combined').attributes().get('_FillValue')
##!#
##!#aod_db_hdf = np.ma.masked_where(aod_db_hdf == fill_value, aod_db_hdf)
##!#
##!## Use the scale factor to correct the data. Sometimes, NASA data files
##!## will also contain 'add_offset' values that are also used to convert
##!## the stored integer values to real data. However, it does not appear
##!## that any of the variables in the MYD04 data file use offsets,
##!## so just use the scale factor to get the data.
##!## --------------------------------------------------------------------
##!#scale_factor = data2.select(\
##!#    'AOD_550_Dark_Target_Deep_Blue_Combined').attributes().get('scale_factor')
##!#
##!#aod_db_hdf = aod_db_hdf * scale_factor
##!#
##!#data2.end()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# BEGIN USER INPUT 
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# For this example, base_path is the path to the raw data files
# -------------------------------------------------------------
#base_path = home_dir + '/data/OMI/H5_files/'

# total_list is the list of files that are being averaged into the a 
# grid.
# ------------------------------------------------------------------
#total_list = subprocess.check_output('ls '+base_path+\
#        'OMI-Aura_L2-OMAERUV_'+year+'m'+date+'t'+time+'*.he5',\
#        shell=True).decode('utf-8').strip().split('\n')

total_list = ['MYD04_L2.A2018001.0000.061.2018003204016.hdf']

# lat_bounds and lon_bounds give the minimum and maximum latitude and 
# longitude values of the overall grid
# ---------------------------------------------------------------------
lat_bounds = [35,60]
lon_bounds = [-180,180]

# resolution gives the lat/lon resolution of the grid. 1.0 creates a 
# 1x1 degree grid, 0.5 gives a half degree grid, 0.25 gives a quarter
# degree grid, 2.0 gives a 2 degree grid...
# ---------------------------------------------------------------------
resolution = 1.0

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# END USER INPUT
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# Set up values for gridding the AI data
# These are all automatically generated from the values set above
# ---------------------------------------------------------------
lat_gridder = lat_bounds[0] * (1. / resolution)
max_idx = int(360 / resolution)
multiplier = int(1 / resolution)
adder = int(180 / resolution)

# Use the user-defined lat/lon bounds to define the grid 
# ------------------------------------------------------
lat_ranges = np.arange(lat_bounds[0],lat_bounds[1], resolution)
lon_ranges = np.arange(lon_bounds[0],lon_bounds[1], resolution)

# Set up grid arrays to hold the data and the counts
# --------------------------------------------------
g_UVAI_354      = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
count_UVAI_354  = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))

# Loop over each OMI file
# -----------------------
for fileI in range(len(total_list)):

    # Open the current data file. NOTE: for files that are not HDF5
    # files, this will have to be changed (for example, if opening a netCDF 
    # file)
    # ---------------------------------------------------------------------
    #data = h5py.File(total_list[fileI],'r')
    data1 = Dataset(test_filename)
    print(total_list[fileI])

    lat = data1['Latitude'][:,:]
    lon = data1['Longitude'][:,:]
    
    # Access one of the AOD parameters
    # --------------------------------
    aod_db = data1['AOD_550_Dark_Target_Deep_Blue_Combined'][:,:]

    ## Extract the latitude, longitude, and data values from the file object
    ## ---------------------------------------------------------------------
    #LAT     = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude'][:,:]
    #LON     = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude'][:,:]
    #UVAI    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex'][:,:]

    # Loop over the current OMI file
    # ------------------------------
    for i in range(aod_db.shape[0]):
        for j in range(aod_db.shape[1]):
            # Make sure that the current file sits within the desired lat/lon range
            if((lat[i,j] > lat_bounds[0]) & (lat[i,j] < lat_bounds[1]) & \
                (lon[i,j] > lon_bounds[0]) & (lon[i,j] < lon_bounds[1])):
                # Make sure that the current data values is not missing or bad
                # ------------------------------------------------------------
                #if(UVAI[i,j]>-2e5):
                # Use the pre-defined grid information to figure out 
                # where the current pixel should fit within the grid
                # --------------------------------------------------
                index1 = int(np.floor(lat[i,j]*multiplier - lat_gridder))
                index2 = int(np.floor(lon[i,j]*multiplier + adder))
                    
                # Add checks to make sure the calculated grid indices
                # are not outside the bounds of the grid array.
                # --------------------------------------------------- 
                if(index1 < 0): index1 = 0
                if(index1 > len(lat_ranges)-1): index1 = len(lat_ranges) - 1
                if(index2 < 0): index2 = 0                                                                                            
                if(index2 > (max_idx - 1)): index2 = max_idx - 1

                print(lat[i,j], lon[i,j], lat_ranges[index1], lon_ranges[index2])

                # Add the current OMI value to the running total
                # ----------------------------------------------
                g_UVAI_354[index2, index1] += aod_db[i,j]

                # Increment the running total of ob counts
                # ----------------------------------------
                count_UVAI_354[index2,index1] +=  1

    # Close the current file object
    # -----------------------------
    data1.close()

# Mask any grid boxes that have no data. These empty boxes are denoted
# by boxes with 0 ob counts.
# --------------------------------------------------------------------
g_UVAI_354     = np.ma.masked_where(count_UVAI_354 == 0, g_UVAI_354)
count_UVAI_354 = np.ma.masked_where(count_UVAI_354 == 0, count_UVAI_354)

# Divide the summed total observation values in each grid box by 
# the number of counts in each box. This calculates the binned average
# --------------------------------------------------------------------
g_UVAI_354 /= count_UVAI_354



