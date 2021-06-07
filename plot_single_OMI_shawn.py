#!/usr/bin/env python
"""
  NAME:
    plot_single_OMI.py

  PURPOSE:
    Plot data from a raw OMI HDF5 file using cartopy. See the 'path_dict'
    dictionary keys for accepted variables.

    NOTE: If running on a computer that is not Blake Sorenson's office computer,
    change the 'base_path' variable to point to the directory containing the
    HDF5 OMI files.

  SYNTAX:
    $ python plot_single_OMI.py date variable row_max

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
    $ python plot_single_OMI.py 201807260741 UVAerosolIndex 60 

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2021/04/01:
      Written (modification of AerosolIndexPlotting_Procedure_August112014t0046.pro by rcontreras) 

"""

import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as color
import cartopy.crs as ccrs
import subprocess
import h5py

var_dict = {
    'SurfaceAlbedo':           {'min': 0.0,  'max': 1.0},\
    'Reflectivity':            {'min': 0.0,  'max': 1.0},\
    'NormRadiance':            {'min': 0.0,  'max': 0.1},\
    'CloudFraction':           {'min': None, 'max': None},\
    'UVAerosolIndex':          {'min': -2.0, 'max': 3.0 },\
    'PixelQualityFlags':       {'min': None, 'max': None},\
    'MeasurementQualityFlags': {'min': None, 'max': None},\
    'ViewingZenithAngle':      {'min': 0.0,  'max': 180.},\
    'SolarZenithAngle':        {'min': None, 'max': None},\
    'RelativeZenithAngle':     {'min': None, 'max': None},\
    'GroundPixelQualityFlags': {'min': None, 'max': None},\
}

if(len(sys.argv)<2):
    print("SYNTAX: python plot_single_OMI_shawn.py date")
    print("\n        date: YYYYMMDDHHMM for single scan")
    print("\n        date: YYYYMMDD     for single day")
    sys.exit()

plot_time = sys.argv[1]

# This is the path that points to the HDF5 OMI files. This must be changed
# if running on a new system.
base_path = '/home/bsorenson/Research/OMI/shawn_files/'

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
# Select the needed data files
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Extract date information to use when finding the data files
year = plot_time[:4]
date = plot_time[4:8]
if(len(plot_time)==12):
    time = plot_time[8:12]
else:
    time = ''

total_list = subprocess.check_output('ls '+base_path+year+date+time+\
    '*',shell=True).decode('utf-8').strip().split('\n')

# Set up values for gridding the AI data
latmin = 60 
lat_gridder = latmin * 4.

lat_ranges = np.arange(latmin,90.1,0.25)
lon_ranges = np.arange(-180,180.1,0.25)

# Set up blank grid arrays to hold the counts and the data
UVAI = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
count = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
# Grid the data
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
for fileI in range(len(total_list)):
    # read in data directly from HDF5 files
    print(total_list[fileI])
    with open(total_list[fileI]) as infile:
        for line in infile:
            templine = line.strip().split()
            lat     = float(templine[0])
            lon     = float(templine[1])
            cleanAI = float(templine[4])
            if((cleanAI>-2e5)):
                #if((j != 52) & (PDATA[i,j]>-2e5)):
                # Only plot if XTrack flag is met
                if(lat > latmin):

                    index1 = int(np.floor(lat*4 - lat_gridder))
                    index2 = int(np.floor(lon*4 + 720.))
                      
                    if(index1 < 0): index1 = 0
                    if(index1 > len(lat_ranges)-1): index1 = len(lat_ranges) - 1
                    if(index2 < 0): index2 = 0                                                                                            
                    if(index2 > 1439): index2 = 1439

                    UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + cleanAI)/(count[index2,index1]+1)
                    count[index2, index1] = count[index2,index1] + 1

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
fig1, axs = plt.subplots(2, 1, figsize=(8,8))
if(latmin<45):
    axs[0] = plt.axes(projection = ccrs.Miller())
else:
    axs[0] = plt.axes(projection = ccrs.NorthPolarStereo(central_longitude = 0.))
axs[0].gridlines()
axs[0].coastlines(resolution = '50m')


# Use meshgrid to convert the 1-d lat/lon arrays into 2-d, which is needed
# for pcolormesh.
plot_lat, plot_lon = np.meshgrid(lat_ranges,lon_ranges)
mask_UVAI = np.ma.masked_where(count == 0, UVAI)

plt.title('OMI ' + ' '+plot_time)
mesh = axs[0].pcolormesh(plot_lon, plot_lat,mask_UVAI,transform = datacrs,cmap = colormap,\
        vmin = -1.0, vmax = 3.0)

# Center the figure over the Arctic
axs[0].set_extent([-180,180,latmin,90],ccrs.PlateCarree())

# Depending on the desired variable, set the appropriate colorbar ticks
tickvals = np.arange(-2.0,4.1,0.5)

cbar = plt.colorbar(mesh,ticks = tickvals,orientation='horizontal',pad=0,\
    aspect=50,shrink = 0.850)
cbar.ax.tick_params(labelsize=14)
cbar.set_label('UV Aerosol Index',fontsize=16,weight='bold')

save = False
if(save == True):
    out_name = 'omi_single_pass_ai'+ '_'+\
        plot_time+'_shawn.png'
    plt.savefig(out_name)
    print('Saved image '+out_name)
else:
    plt.show()
