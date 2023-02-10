#!/usr/bin/env python
"""
  NAME:
    plot_single_TROPOMI.py

  PURPOSE:
    Plot data from a raw OMI HDF5 file using cartopy. See the 'path_dict'
    dictionary keys for accepted variables.

    NOTE: If running on a computer that is not Blake Sorenson's office computer,
    change the 'base_path' variable to point to the directory containing the
    HDF5 OMI files.

  SYNTAX:
    $ python plot_single_TROPOMI.py date variable row_max

       - date: for a single-swath image, use a format of YYYYMMDDHH
               for a daily average, use a format of YYYYMMDD
       - variable: variable name from the HDF5 file. This code is only
               set up to handle certain variables right now, so to see
               which variables are currently allowed, see either the syntax
               message or the keys to any of the dictionaries found below
       - row_max:  the maximum row to plot data through. For example, if you
               only want to plot rows 1 to 23, use a row_max value of 23.
               Otherwise, use a value of 60 to plot all rows.

  EXAMPLE:
    $ python plot_single_TROPOMI.py 201807260741 UVAerosolIndex 60 

  MODIFICATIONS:
  Jacob Zanker- conversion from OMI- 2021/11/19
    Blake Sorenson <blake.sorenson@und.edu>     - 2021/04/01:
      Written (modification of AerosolIndexPlotting_Procedure_August112014t0046.pro by rcontreras) 

"""

import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as color
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import subprocess
import h5py
import netCDF4 as nc4
from netCDF4 import Dataset

# Set up dictionary to hold file naming variables
name_dict = {
    'SurfaceAlbedo':           'alb',\
    'NormRadiance':            'nrad',\
    'Reflectivity':            'refl',\
    'CloudFraction':           'cldfrc',\
    'UVAerosolIndex':          'uvai',\
    'PixelQualityFlags':       'pxlqf',\
    'MeasurementQualityFlags': 'msmtqf',\
    'FinalAlgorithmFlags':     'falgf',\
    'ViewingZenithAngle':      'vza',\
    'SolarZenithAngle':        'sza',\
    'RelativeZenithAngle':     'rza',\
    'GroundPixelQualityFlags': 'grndqf'
}

var_dict = {
    'SurfaceAlbedo':           {'min': 0.0,  'max': 1.0},\
    'Reflectivity':            {'min': 0.0,  'max': 1.0},\
    'NormRadiance':            {'min': 0.0,  'max': 0.2},\
    'CloudFraction':           {'min': None, 'max': None},\
    'UVAerosolIndex':          {'min': -2.0, 'max': 10.0 },\
    'PixelQualityFlags':       {'min': None, 'max': None},\
    'MeasurementQualityFlags': {'min': None, 'max': None},\
    'FinalAlgorithmFlags':     {'min':    0, 'max':    8},\
    'ViewingZenithAngle':      {'min': 0.0,  'max': 180.},\
    'SolarZenithAngle':        {'min': None, 'max': None},\
    'RelativeZenithAngle':     {'min': None, 'max': None},\
    'GroundPixelQualityFlags': {'min': None, 'max': None},\
}

if(len(sys.argv)<4):
    print("SYNTAX: python plot_single_TROPOMI_backup.py date variable row_max")
    print("\n        date: YYYYMMDDHHMM for single scan")
    print("              YYYYMMDD for a daily average")
    print("\n        Accepted variable names:")
    for name in name_dict.keys():
        print("        - "+name)
    print("\n        row_max: the highest row from which to plot data")
    print("                 to plot all data, use 60")
    sys.exit()

plot_time = sys.argv[1]
variable = sys.argv[2]
row_max=int(sys.argv[3])
channel_idx = 0
str_wave = ''

# This is the path that points to the TROPOMI files. This must be changed
# if running on a new system.
base_path = '/home/jzanker/'

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
# Select the needed data files
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Extract date information to use when finding the data files
year = plot_time[:4]
date = plot_time[4:8]
if(len(plot_time)==14):
    time = plot_time[8:]
else:
    time = ''

total_list = subprocess.check_output('ls '+base_path+'S5P_OFFL_L2__AER_AI_'+\
    year+date+'T'+time+'*.nc',shell=True).decode('utf-8').strip().split('\n')

# Set up values for gridding the AI data
latmin = 60 #25 for CA
lat_gridder = latmin * 4.

lat_ranges = np.arange(latmin,90.1,0.25)
lon_ranges = np.arange(-180,180.1,0.25)

print(total_list[0])
data = Dataset(total_list[0], mode='r') 
LAT   = data.groups['PRODUCT'].variables['latitude'][:][0,:,:]
LON   = data.groups['PRODUCT'].variables['longitude'][:][0,:,:]
UVAI  = data.groups['PRODUCT'].variables['aerosol_index_354_388'][0,:,:]
if(len(UVAI.shape) == 3):
    # If 3 dimensions, assume that the smallest dimension is the wavelength
    # dimension. Find that dimension and grab the first index of that
    # dimension, which is the 354 nm channel.
    min_dim = np.min(UVAI.shape)
    if(UVAI.shape[2] == min_dim):
        UVAI = UVAI[:,:,channel_idx]
    elif(UVAI.shape[1] == min_dim):
        UVAI = UVAI[:,channel_idx,:]
    else:
        UVAI = UVAI[channel_idx,:,:]
XTRACK = data.groups['PRODUCT'].variables['qa_value'][:][0,:,:]
mask_UVAI = np.ma.masked_where((XTRACK < -2e5) | (UVAI < -2e5) | (LAT < latmin), UVAI)
plot_lon = LON
plot_lat = LAT

data.close()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
# Plot the gridded data
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
print("Time to plot")
mapcrs = ccrs.NorthPolarStereo()
datacrs = ccrs.PlateCarree()
colormap = plt.cm.jet

# Set up the polar stereographic projection map
fig1 = plt.figure(figsize=(8,8))
if(latmin<45):
    ax = plt.axes(projection = ccrs.Miller())
else:
    ax = plt.axes(projection = ccrs.NorthPolarStereo(central_longitude = 0.))
ax.gridlines()
ax.coastlines(resolution = '50m')


plt.title('TROPOMI ' + variable + str_wave + ' '+plot_time)
mesh = ax.pcolormesh(plot_lon, plot_lat,mask_UVAI,transform = datacrs,cmap = colormap,\
        vmin = var_dict[variable]['min'], vmax = var_dict[variable]['max'],\
        shading='auto')

# Center the figure over the Arctic
ax.set_extent([-180,180,latmin,90],ccrs.PlateCarree())
#For western USA
#ax.set_extent([-130, -110, 25, 50], datacrs)
#ax.set_extent([-121.5, -119.5, 39.5, 42.5], datacrs) # Values for 2021/08/05 plume case
#ax.set_extent([-118.0, -114.0, 36.0, 39.0], datacrs) # Values for 2021/08/06 plume case
#ax.set_extent([-122.0, -119.5, 39.5, 42.0], datacrs) # Values for 2021/08/05 plume case
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.STATES)
# Depending on the desired variable, set the appropriate colorbar ticks
if(variable == 'NormRadiance'):
    tickvals = np.arange(0.0,1.010,0.10)
elif(variable == 'Reflectivity'):
    tickvals = np.arange(0.0,1.1,0.10)
elif(variable == 'FinalAlgorithmFlags'):
    tickvals = np.arange(0,9)
elif(variable == 'SolarZenithAngle'):
    tickvals = np.arange(0,90,10)
else:
    tickvals = np.arange(-2.0,11.1,1.0)

cbar = plt.colorbar(mesh,ticks = tickvals,orientation='horizontal',pad=0,\
    aspect=50,shrink = 0.850)
cbar.ax.tick_params(labelsize=14)
cbar.set_label(variable,fontsize=16,weight='bold')

save = True
if(save == True):
    out_name = 'tropomi_single_pass_'+name_dict[variable] + str_wave + '_'+\
        plot_time+'_rows_0to'+str(row_max)+'.png'
    plt.savefig(out_name)
    print('Saved image '+out_name)
else:
    plt.show()
total_list.clear()
