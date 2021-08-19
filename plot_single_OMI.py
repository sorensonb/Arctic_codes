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

# Set up variable dictionary. This is used in pulling the variables
# from the HDF5 structure.
path_dict = {
    'SurfaceAlbedo':           'Data Fields/',\
    'NormRadiance':            'Data Fields/',\
    'Reflectivity':            'Data Fields/',\
    'CloudFraction':           'Data Fields/',\
    'UVAerosolIndex':          'Data Fields/',\
    'PixelQualityFlags':       'Data Fields/',\
    'MeasurementQualityFlags': 'Data Fields/',\
    'ViewingZenithAngle':      'Geolocation Fields/',\
    'SolarZenithAngle':        'Geolocation Fields/',\
    'RelativeZenithAngle':     'Geolocation Fields/',\
    'GroundPixelQualityFlags': 'Geolocation Fields/'
}

# Set up dictionary to hold file naming variables
name_dict = {
    'SurfaceAlbedo':           'alb',\
    'NormRadiance':            'nrad',\
    'Reflectivity':            'refl',\
    'CloudFraction':           'cldfrc',\
    'UVAerosolIndex':          'uvai',\
    'PixelQualityFlags':       'pxlqf',\
    'MeasurementQualityFlags': 'msmtqf',\
    'ViewingZenithAngle':      'vza',\
    'SolarZenithAngle':        'sza',\
    'RelativeZenithAngle':     'rza',\
    'GroundPixelQualityFlags': 'grndqf'
}

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

if(len(sys.argv)<4):
    print("SYNTAX: python plot_single_OMI.py date variable row_max")
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

# Check which channels are being used, which only matters for the radiances,
# reflectivities, and albedoes.
if((variable[:12] == 'NormRadiance')  | \
    (variable[:12] == 'Reflectivity') | \
    (variable[:13] == 'SurfaceAlbedo')):
    if(variable[-3:] == '354'):
        variable = variable[:-3]
        str_wave = '354nm'
        channel_idx = 0
    elif(variable[-3:] == '388'):
        variable = variable[:-3]
        str_wave = '388nm'
        channel_idx = 1
    elif(variable[-3:] == '500'):
        variable = variable[:-3]
        str_wave = '500nm'
        channel_idx = 2
    else:
        print("Channel not selected. Defaulting to 354 nm")
        str_wave = '354nm'
        channel_idx = 0

# This is the path that points to the HDF5 OMI files. This must be changed
# if running on a new system.
base_path = '/home/bsorenson/data/OMI/H5_files/'

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
# Select the needed data files
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Extract date information to use when finding the data files
year = plot_time[:4]
date = plot_time[4:8]
if(len(plot_time)==13):
    time = plot_time[9:]
elif(len(plot_time)==12):
    time = plot_time[8:]
else:
    time = ''

total_list = subprocess.check_output('ls '+base_path+'OMI-Aura_L2-OMAERUV_'+\
    year+'m'+date+'t'+time+'*.he5',shell=True).decode('utf-8').strip().split('\n')

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
    data = h5py.File(total_list[fileI],'r')

    LAT   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude']
    LON   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude']
    PDATA = data['HDFEOS/SWATHS/Aerosol NearUV Swath/' + path_dict[variable] + variable]
    if(len(PDATA.shape) == 3):
        # If 3 dimensions, assume that the smallest dimension is the wavelength
        # dimension. Find that dimension and grab the first index of that
        # dimension, which is the 354 nm channel.
        min_dim = np.min(PDATA.shape)
        if(PDATA.shape[2] == min_dim):
            PDATA = PDATA[:,:,channel_idx]
        elif(PDATA.shape[1] == min_dim):
            PDATA = PDATA[:,channel_idx,:]
        else:
            PDATA = PDATA[channel_idx,:,:]
    XTRACK = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags']

    # Loop over the values and rows 
    for i in range(PDATA.shape[0]):
        for j in range(0,row_max):
            if((PDATA[i,j]>-2e5)):
            #if((j != 52) & (PDATA[i,j]>-2e5)):
            # Only plot if XTrack flag is met
                if((XTRACK[i,j] == 0) | (XTRACK[i,j] == 4)):
                    # Print values to text file
                    if(LAT[i,j] > latmin):

                        index1 = int(np.floor(LAT[i,j]*4 - lat_gridder))
                        index2 = int(np.floor(LON[i,j]*4 + 720.))
                        
                        if(index1 < 0): index1 = 0
                        if(index1 > len(lat_ranges)-1): index1 = len(lat_ranges) - 1
                        if(index2 < 0): index2 = 0                                                                                            
                        if(index2 > 1439): index2 = 1439

                        UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + PDATA[i,j])/(count[index2,index1]+1)
                        count[index2, index1] = count[index2,index1] + 1
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

# Determine the percentage of grid boxes that are actually filled
# with values.
total_boxes = mask_UVAI.size
total_good = mask_UVAI.count()
pcnt_good = (total_good / total_boxes) * 100.
print("Total_boxes = ",total_boxes,"Total good = ",total_good)
print("Percent good = ",pcnt_good)

plt.title('OMI ' + variable + str_wave + ' '+plot_time)
mesh = axs[0].pcolormesh(plot_lon, plot_lat,mask_UVAI,transform = datacrs,cmap = colormap,\
        vmin = var_dict[variable]['min'], vmax = var_dict[variable]['max'])

# Center the figure over the Arctic
axs[0].set_extent([-180,180,latmin,90],ccrs.PlateCarree())

# Depending on the desired variable, set the appropriate colorbar ticks
if(variable == 'NormRadiance'):
    tickvals = np.arange(0.0,0.101,0.01)
elif(variable == 'Reflectivity'):
    tickvals = np.arange(0.0,1.1,0.10)
else:
    tickvals = np.arange(-2.0,4.1,0.5)

cbar = plt.colorbar(mesh,ticks = tickvals,orientation='horizontal',pad=0,\
    aspect=50,shrink = 0.850)
cbar.ax.tick_params(labelsize=14)
cbar.set_label('UV Aerosol Index',fontsize=16,weight='bold')

save = True 
if(save == True):
    out_name = 'omi_single_pass_'+name_dict[variable] + str_wave + '_'+\
        plot_time+'_rows_0to'+str(row_max)+'.png'
    plt.savefig(out_name)
    print('Saved image '+out_name)
else:
    plt.show()
