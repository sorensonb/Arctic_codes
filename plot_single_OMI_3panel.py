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
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import subprocess
import h5py

def get_ice_flags(value):
    return int(format(value,"016b")[-15:-8],2)

var_dict = {
    'UVAerosolIndex':          {'min': -2.0, 'max': 3.0 },\
}

if(len(sys.argv)<2):
    print("SYNTAX: python plot_single_OMI.py date")
    print("\n        date: YYYYMMDDHHMM for single scan")
    print("              YYYYMMDD for a daily average")
    sys.exit()

plot_time = sys.argv[1]
variable = 'UVAerosolIndex'
row_max = 60
channel_idx = 0
str_wave = ''

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
# Select the needed data files
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Read in data files for JZ and control files
# - - - - - - - - - - - - - - - - - - - - - -

# Extract date information to use when finding the data files
year = plot_time[:4]
date = plot_time[4:8]
if(len(plot_time)==12):
    time = plot_time[8:12]
else:
    time = ''

# This is the path that points to the HDF5 OMI files. This must be changed
# if running on a new system.
base_path = '/home/bsorenson/data/OMI/H5_files/'

# Pull in file names for control and JZ data files
total_h5_list = subprocess.check_output('ls '+base_path+'OMI-Aura_L2-OMAERUV_'+\
    year+'m'+date+'t'+time+'*.he5',shell=True).decode('utf-8').strip().split('\n')

# Read in data files for JZ and control files
# - - - - - - - - - - - - - - - - - - - - - -

# This is the path that points to Shawn's clean OMI files. This must be changed
# if running on a new system.
base_path = '/home/bsorenson/data/OMI/shawn_files/'

total_sh_list = subprocess.check_output('ls '+base_path+year+date+time+\
    '*',shell=True).decode('utf-8').strip().split('\n')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
# Grid the control and JZ data
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Set up values for gridding the AI data
latmin = 60 
lat_gridder = latmin * 4.

lat_ranges = np.arange(latmin,90.1,0.25)
lon_ranges = np.arange(-180,180.1,0.25)

# Set up blank grid arrays to hold the counts and the data
UVAI_c   = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
UVAI_JZ  = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
count_c  = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
count_JZ = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))

for fileI in range(len(total_h5_list)):
    # read in data directly from HDF5 files
    print(total_h5_list[fileI])
    data = h5py.File(total_h5_list[fileI],'r')

    LAT    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude']
    LON    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude']
    UVAI   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex']
    GPQF   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/GroundPixelQualityFlags']
    XTRACK = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags']
    AZM    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/RelativeAzimuthAngle']

    # Loop over the values and rows 
    for i in range(UVAI.shape[0]):
        for j in range(0,row_max):
            # Only plot if XTrack flag is met
            if((LAT[i,j] > latmin) & \
               (UVAI[i,j]>-2e5) & \
               ((XTRACK[i,j] == 0) | (XTRACK[i,j] == 4))):

                # If these first conditions are met, grid these data into the
                # control averages.
                index1 = int(np.floor(LAT[i,j]*4 - lat_gridder))
                index2 = int(np.floor(LON[i,j]*4 + 720.))
                        
                if(index1 < 0): index1 = 0
                if(index1 > len(lat_ranges)-1): index1 = len(lat_ranges) - 1
                if(index2 < 0): index2 = 0                                                                                            
                if(index2 > 1439): index2 = 1439

                UVAI_c[index2, index1] = (UVAI_c[index2,index1]*count_c[index2,index1] + UVAI[i,j])/(count_c[index2,index1]+1)
                count_c[index2, index1] = count_c[index2,index1] + 1

                pixel_val = get_ice_flags(GPQF[i,j])
                # After inserting the control values, check the JZ criteria
                # before inserting into the JZ arrays
                if((((pixel_val >= 0) & (pixel_val <= 101)) | (pixel_val == 104)) & \
                    (AZM[i,j] > 100) ):
                    UVAI_JZ[index2, index1] = (UVAI_JZ[index2,index1]*count_JZ[index2,index1] + UVAI[i,j])/(count_JZ[index2,index1]+1)
                    count_JZ[index2, index1] = count_JZ[index2,index1] + 1
    data.close()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
# Grid the Shawn data
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Set up blank grid arrays to hold the counts and the data
UVAI_sh  = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
count_sh = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))

for fileI in range(len(total_sh_list)):
    # read in data directly from HDF5 files
    print(total_sh_list[fileI])
    with open(total_sh_list[fileI]) as infile:
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

                    UVAI_sh[index2, index1] = (UVAI_sh[index2,index1]*count_sh[index2,index1] + cleanAI)/(count_sh[index2,index1]+1)
                    count_sh[index2, index1] = count_sh[index2,index1] + 1


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
# Plot the gridded data
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
print("Time to plot")
mapcrs = ccrs.NorthPolarStereo(central_longitude = 0)
datacrs = ccrs.PlateCarree()
colormap = plt.cm.jet

# Plot the control data in the first subplot
# - - - - - - - - - - - - - - - - - - - - - -

# Set up the polar stereographic projection map
#fig1, axs = plt.subplots(1, 3,subplot_kw={'projection':mapcrs},figsize=(11,8.5))
fig1 = plt.figure(1,figsize=(12,5))
gs = gridspec.GridSpec(nrows=1, ncols=3, hspace = 0.00, wspace=0.08)
plt.suptitle(plot_time)
#if(latmin<45):
#    axs[0] = plt.axes(projection = ccrs.Miller())
#else:
#axs[0] = plt.axes(projection = mapcrs)
ax0 = plt.subplot(gs[0,0],projection=mapcrs)
ax0.gridlines()
ax0.coastlines(resolution = '50m')

# Use meshgrid to convert the 1-d lat/lon arrays into 2-d, which is needed
# for pcolormesh.
plot_lat, plot_lon = np.meshgrid(lat_ranges,lon_ranges)
mask_UVAI = np.ma.masked_where(count_c == 0, UVAI_c)

ax0.set_title('Control',fontsize=12)
mesh0 = ax0.pcolormesh(plot_lon, plot_lat,mask_UVAI,transform = datacrs,cmap = colormap,\
        vmin = var_dict[variable]['min'], vmax = var_dict[variable]['max'])

# Center the figure over the Arctic
ax0.set_extent([-180,180,latmin,90],ccrs.PlateCarree())

# Depending on the desired variable, set the appropriate colorbar ticks
tickvals = np.arange(-2.0,4.1,1.0)

cbar0 = plt.colorbar(mesh0,ticks = tickvals,orientation='horizontal',pad=0,\
    aspect=50)
#cbar0.ax.tick_params(labelsize=14)
cbar0.set_label('UV Aerosol Index',weight='bold')

# Plot the JZ data in the second subplot
# - - - - - - - - - - - - - - - - - - - - - -

#axs[1] = plt.axes(projection = mapcrs)
ax1 = plt.subplot(gs[0,1],projection=mapcrs)
ax1.gridlines()
ax1.coastlines(resolution = '50m')

# Use meshgrid to convert the 1-d lat/lon arrays into 2-d, which is needed
# for pcolormesh.
plot_lat, plot_lon = np.meshgrid(lat_ranges,lon_ranges)
mask_UVAI_JZ = np.ma.masked_where(count_JZ == 0, UVAI_JZ)

ax1.set_title('Screened',fontsize=12)
mesh1 = ax1.pcolormesh(plot_lon, plot_lat,mask_UVAI_JZ,transform = datacrs,cmap = colormap,\
        vmin = var_dict[variable]['min'], vmax = var_dict[variable]['max'])

# Center the figure over the Arctic
ax1.set_extent([-180,180,latmin,90],ccrs.PlateCarree())

# Depending on the desired variable, set the appropriate colorbar ticks
tickvals = np.arange(-2.0,4.1,1.0)

cbar1 = plt.colorbar(mesh1,ticks = tickvals,orientation='horizontal',pad=0,\
    aspect=50)#,shrink = 0.850)
#cbar1.ax.tick_params(labelsize=14)
cbar1.set_label('UV Aerosol Index',weight='bold')

# Plot the Shawn data in the second subplot
# - - - - - - - - - - - - - - - - - - - - - -

#axs[2] = plt.axes(projection = mapcrs)
ax2 = plt.subplot(gs[0,2],projection=mapcrs)
ax2.gridlines()
ax2.coastlines(resolution = '50m')

# Use meshgrid to convert the 1-d lat/lon arrays into 2-d, which is needed
# for pcolormesh.
plot_lat, plot_lon = np.meshgrid(lat_ranges,lon_ranges)
mask_UVAI_sh = np.ma.masked_where(count_sh == 0, UVAI_sh)

ax2.set_title('Perturbation',fontsize=12)
mesh2 = ax2.pcolormesh(plot_lon, plot_lat,mask_UVAI_sh,transform = datacrs,cmap = colormap,\
        vmin = var_dict[variable]['min'], vmax = var_dict[variable]['max'])

# Center the figure over the Arctic
ax2.set_extent([-180,180,latmin,90],ccrs.PlateCarree())

# Depending on the desired variable, set the appropriate colorbar ticks
tickvals = np.arange(-2.0,4.1,1.0)

cbar2 = plt.colorbar(mesh2,ticks = tickvals,orientation='horizontal',pad=0,\
    aspect=50)#,shrink = 0.850)
#cbar2.ax.tick_params(labelsize=14)
cbar2.set_label('UV Aerosol Index Perturbation',weight='bold')

save = True 
if(save == True):
    out_name = 'omi_single_pass_uvai_' +\
        plot_time+'_3panel.png'
    plt.savefig(out_name,dpi=300)
    print('Saved image '+out_name)
else:
    plt.show()
