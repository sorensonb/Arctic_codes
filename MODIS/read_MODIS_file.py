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

# Open the file using netCDF4
# ---------------------------
data1 = Dataset(test_filename)

# If desired, you can print the names of all the variables in the file using
# the following line:
#print(data1.variables.keys())

# Access the latitude and longitude data
# --------------------------------------
lat = data1['Latitude'][:,:]
lon = data1['Longitude'][:,:]

# Access one of the AOD parameters
# --------------------------------
aod_db = data1['AOD_550_Dark_Target_Deep_Blue_Combined'][:,:]

# NOTE: netCDF4 automatically accounts for the scale factor and offset values
#       when reading an attribute from a file. Thus, no further action
#       is needed to use the AOD data.

# Close the file
# --------------
data1.close()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Method 2: Read the file using pyhdf
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

data2 = SD.SD(test_filename)

# Access the latitude and longitude data
# --------------------------------------
lat = data2.select('Latitude').get()
lon = data2.select('Longitude').get()

# Access one of the AOD parameters
# --------------------------------
aod_db_hdf = data2.select('AOD_550_Dark_Target_Deep_Blue_Combined').get()

# Remove any missing values. Missing values are defined using the 'fill value'
# parameter.
# ----------------------------------------------------------------------------
fill_value = data2.select(\
    'AOD_550_Dark_Target_Deep_Blue_Combined').attributes().get('_FillValue')

aod_db_hdf = np.ma.masked_where(aod_db_hdf == fill_value, aod_db_hdf)

# Use the scale factor to correct the data. Sometimes, NASA data files
# will also contain 'add_offset' values that are also used to convert
# the stored integer values to real data. However, it does not appear
# that any of the variables in the MYD04 data file use offsets,
# so just use the scale factor to get the data.
# --------------------------------------------------------------------
scale_factor = data2.select(\
    'AOD_550_Dark_Target_Deep_Blue_Combined').attributes().get('scale_factor')

aod_db_hdf = aod_db_hdf * scale_factor

data2.end()
