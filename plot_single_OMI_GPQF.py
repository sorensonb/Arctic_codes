#!/usr/bin/env python
"""
  NAME:

  PURPOSE:

  SYNTAX:

  EXAMPLE:

  MODIFICATIONS:

"""

import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as color
import cartopy.crs as ccrs
##import matplotlib.colors as colors
import subprocess
import h5py

def get_ice_flags(value):
    return int(format(value,"016b")[-15:-8],2)


if(len(sys.argv)<3):
    print("SYNTAX: python plot_single_OMI_GPQF.py date row_max")
    sys.exit()

plot_time = sys.argv[1]
row_max=int(sys.argv[2])

#new_time = '20130611t1322_new'
base_path = '/home/bsorenson/data/OMI/H5_files/'

n_p = 1440
nl = 720

#manually enter the longitude limits
lonmin = -180
lonmax = 180

latmax =  90

#date = sys.argv[1]

# Set up the dictionary
latmin = 60 

# Set up values for gridding the AI data
lat_gridder = latmin * 4.

lat_ranges = np.arange(latmin,90.1,0.25)
lon_ranges = np.arange(-180,180.1,0.25)

# Set up cartopy variables
mapcrs = ccrs.NorthPolarStereo()
datacrs = ccrs.PlateCarree()
colormap = plt.cm.jet

# Set up the polar stereographic projection map
fig1, axs = plt.subplots(2, 1, figsize=(8,8))
#fig1.set_size_inches(8,8)
#fig1 = plt.figure(figsize=(8,8))
if(latmin<45):
    ax = plt.axes(projection = ccrs.Miller())
    #m = Basemap(projection='mill',lon_0=0,resolution='l')
else:
    axs[0] = plt.axes(projection = ccrs.NorthPolarStereo())
    #ax = plt.axes(projection = ccrs.NorthPolarStereo(central_longitude = 45.))

#fig = plt.gcf()
axs[0].gridlines()
axs[0].coastlines(resolution = '50m')


year = plot_time[:4]
date = plot_time[4:8]
if(len(plot_time)==13):
    time = plot_time[9:]
elif(len(plot_time)==12):
    time = plot_time[8:]
else:
    time = ''
total_list = subprocess.check_output('ls '+base_path+'OMI-Aura_L2-OMAERUV_'+year+'m'+date+'t'+time+'*.he5',\
          shell=True).decode('utf-8').strip().split('\n')


g_GPQF = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
count_GPQF = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))

for fileI in range(len(total_list)):
    # read in data directly from HDF5 files
    data = h5py.File(total_list[fileI],'r')
    print(total_list[fileI])

    LAT_flat  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude'][:,0:row_max].flatten()
    LON_flat  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude'][:,0:row_max].flatten()
    GPQF_flat = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/GroundPixelQualityFlags'][:,0:row_max].flatten()

    GPQF_decode = np.array([get_ice_flags(val) for val in GPQF_flat])

    for i in range(GPQF_flat.shape[0]):
        if(LAT_flat[i] > latmin):

            index1 = int(np.floor(LAT_flat[i]*4 - lat_gridder))
            index2 = int(np.floor(LON_flat[i]*4 + 720.))
        
            if(index1 < 0): index1 = 0
            if(index1 > len(lat_ranges)-1): index1 = len(lat_ranges) - 1
            if(index2 < 0): index2 = 0                                                                                            
            if(index2 > 1439): index2 = 1439

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            # Plot ice as one color
            if(GPQF_decode[i] == 0):                         g_GPQF[index2,index1] = 8
            elif((GPQF_decode[i] > 0) & (GPQF_decode[i] < 101)): g_GPQF[index2,index1] = 9
            elif(GPQF_decode[i] == 101):                       g_GPQF[index2,index1] = 10
            elif(GPQF_decode[i] == 103):                       g_GPQF[index2,index1] = 11
            elif(GPQF_decode[i] == 104):                       g_GPQF[index2,index1] = 12
            elif(GPQF_decode[i] == 124):                       g_GPQF[index2,index1] = 13
            else:                                            g_GPQF[index2,index1] = -9

            count_GPQF[index2,index1] = count_GPQF[index2,index1] + 1
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

            ##!## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            ##!## Plot gradients in ice concentration
            ##!#if((GPQF_flat[i] > 0) & (GPQF_flat[i] < 101)):
            ##!#    g_GPQF[index2, index1] = (g_GPQF[index2,index1]* \
            ##!#        count_GPQF[index2,index1] + GPQF_decode[i]) / \
            ##!#        (count_GPQF[index2,index1] + 1)
            ##!#else:
            ##!#    g_GPQF[index2, index1] = GPQF_decode[i]
            ##!#count_GPQF[index2,index1] = count_GPQF[index2,index1] + 1
            ##!## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    data.close()
    
print("Time to plot")

#create an generic lat and lon
#define the lat which will be gridded into 0.25 x 0.25 bins
LATalt = lat_ranges
LONalt = lon_ranges
#LATalt = np.arange(-90,90,0.25)
#LONalt = np.arange(-180,180,0.25)


cbar_labels = ['Shallow ocean','Land','Shallow inland water',\
                'Ocean coastline / lake shoreline',\
                'ephemeral (intermittent) water','Deep inland water',\
                'Continental shelf ocean','Deep ocean','Snow-free land',\
                'Sea ice','Permanent ice\n(Greenland, Antarctica)','Dry snow',\
                'Ocean','Mixed pixels\nat coastline','Suspect ice value','Corners']

#cmap = plt.cm.jet

# Set up plot

plot_lat, plot_lon = np.meshgrid(LATalt,LONalt)

mask_UVAI = np.ma.masked_where(count_GPQF[:] == 0, g_GPQF)
mask_UVAI = np.ma.masked_where(mask_UVAI[:] == -9, mask_UVAI)

print(np.max(mask_UVAI), np.min(mask_UVAI))
cmap = plt.get_cmap('jet', np.max(mask_UVAI) - np.min(mask_UVAI)+1)

plt.title('OMI GroundPixelQualityFlags\n' + plot_time)
#plt.title('NRAD500 - NRAD354')
#plt.title('OMI Reflectivity - Surface Albedo '+plot_time)
mesh = axs[0].pcolormesh(plot_lon, plot_lat,mask_UVAI,transform = datacrs, \
    cmap = cmap, vmin = np.min(mask_UVAI) - 0.5, vmax = np.max(mask_UVAI) + 0.5)
axs[0].set_extent([-180,180,60,90],ccrs.PlateCarree())
#axs[0].set_xlim(-0630748.535086173,2230748.438879491)
#axs[0].set_ylim(-2513488.8763307533,0343353.899053069)
#ax.set_xlim(-4170748.535086173,4167222.438879491)
#ax.set_ylim(-2913488.8763307533,2943353.899053069)
#cbar = plt.colorbar(mesh,ticks = np.arange(-2.0,4.1,0.01),orientation='horizontal',pad=0,\
cbar = plt.colorbar(mesh,orientation='horizontal',pad=0,\
    aspect=50,shrink = 0.905, ticks = np.arange(np.min(mask_UVAI), np.max(mask_UVAI) + 1))
cbar.ax.set_xticklabels(cbar_labels[int(np.min(mask_UVAI)):int(np.max(mask_UVAI))+1],\
    fontsize=8,rotation=35)

save = True 
if(save == True):
    out_name = 'omi_single_pass_gpqf_'+plot_time+'_rows_0to'+str(row_max)+'.png'       
    plt.savefig(out_name)
    print('Saved image '+out_name)

plt.show()
