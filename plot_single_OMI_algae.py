#!/usr/bin/env python
"""
  NAME:
    aerosol_trend.py

  PURPOSE:
    Read in OMI data and plot AI values over the northern hemisphere from
    48 degrees North to the North Pole. Plots one day of data.

  SYNTAX:
    ./aerosol_trend.py date
    
    - date is in the format of CCYYMMDD

  EXAMPLE:
    ./aerosol_reader.py 20140811

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2018/09/07:
      Written (modification of AerosolIndexPlotting_Procedure_August112014t0046.pro by rcontreras) 

"""

import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as color
import matplotlib.patches as mpatches
import cartopy.crs as ccrs
import subprocess
import h5py

# Set up variable dictionary
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

if(len(sys.argv)<3):
    print("SYNTAX: python plot_single_OMI.py date row_max")
    sys.exit()

plot_time = sys.argv[1]
row_max=int(sys.argv[2])
channel_idx = 0
str_wave = ''

#new_time = '20130611t1322_new'
base_path = '/home/bsorenson/data/OMI/H5_files/'
#new_base_path = '/home/shared/OMAERUV/OMAERUV_Parameters/new_files/'

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

# Calculate average and standard deviations
#newch1=np.zeros(shape=(len(lat_ranges),len(lon_ranges)))

##UVAI = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
newch1=np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
#newch1=np.zeros(shape=(1440,720))

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

##total_list = ['/home/shared/OMAERUV/H5_files_20190212_download/OMI-Aura_L2-OMAERUV_2008m0427t0052-o20122_v003-2017m0721t120210.he5',\
##              '/home/shared/OMAERUV/H5_files_20190212_download/OMI-Aura_L2-OMAERUV_2008m0427t0231-o20123_v003-2017m0721t120217.he5',\
##              '/home/shared/OMAERUV/H5_files_20190212_download/OMI-Aura_L2-OMAERUV_2008m0427t0410-o20124_v003-2017m0721t120234.he5',\
##              '/home/shared/OMAERUV/H5_files_20190212_download/OMI-Aura_L2-OMAERUV_2008m0427t0549-o20125_v003-2017m0721t121042.he5']

min_diff = 30
max_diff = -30

g_NRAD_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
g_NRAD_388 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
g_NRAD_500 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
g_REFL_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
g_REFL_388 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
g_SALB_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
g_UVAI_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
count_NRAD_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
count_NRAD_388 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
count_NRAD_500 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
count_REFL_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
count_REFL_388 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
count_SALB_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
count_UVAI_354 = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))

for fileI in range(len(total_list)):
    # read in data directly from HDF5 files
    data = h5py.File(total_list[fileI],'r')
    print(total_list[fileI])

    LAT     = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude']
    LON     = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude']
    NRAD    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/NormRadiance']
    REFL    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/Reflectivity']
    SALB    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/SurfaceAlbedo']
    UVAI    = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex']
    XTRACK  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags']


    ###ALBEDO = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/SurfaceAlbedo']
    ###REFLECTANCE = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/Reflectivity']
    ###CLD = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/CloudFraction']
    ###AI  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex']
    ###PIXEL= data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/PixelQualityFlags']
    ###MSMNT= data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/MeasurementQualityFlags']
    ###VZA  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/ViewingZenithAngle']
    ###SZA  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/SolarZenithAngle']
    ###RAZ  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/RelativeAzimuthAngle']
    ###GRND = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/GroundPixelQualityFlags']

    #albedo = ALBEDO[:,:,0]   
    #reflectance = REFLECTANCE[:,:,0]   
    counter = 0
    #AI = AI[:,:,0]   
    # Loop over the values and rows 
    #for i in range(0,int(CBA2)):
    #for i in range(albedo.shape[0]):
    for i in range(NRAD.shape[0]):
        for j in range(0,row_max):
            #if((albedo[i,j]>-20) & (reflectance[i,j]>-20)):
            if(NRAD[i,j,0]>-2e5):
            #if(plotAI[i,j]>-20):
                # Only plot if XTrack flag is met
                if((XTRACK[i,j] == 0) | (XTRACK[i,j] == 4)):
                    # Print values to text file
                    if(LAT[i,j] > latmin):
                        counter+=1
                         ##fout.write("{0:.6f} {1:.6f} {2:.6f} {3:.6f} {4:.6f} {5:.6f} {6:.6f} {7:.6f} {8:.6f} {9:.6f} {10:.6f} {11:.6f}\n".format(\
                         ##    LAT[i,j],LON[i,j],AI[i,j],0.5,SZA[i,j],VZA[i,j],RAZ[i,j], \
                         ##    ALBEDO[i,j,0],ALBEDO[i,j,1],REFLECTANCE[i,j,0],\
                         ##    REFLECTANCE[i,j,1],CLD[i,j]))


                    #if((plotXTrack[i,j] == 0) | ((plotXTrack[i,j] & 4 == 4))):
                        index1 = int(np.floor(LAT[i,j]*4 - lat_gridder))
                        index2 = int(np.floor(LON[i,j]*4 + 720.))
                        #index1 = int(np.floor(plotLAT[i,j]*4 + 360.))
                        #index2 = int(np.floor(plotLON[i,j]*4 + 720.))
                    
                        if(index1 < 0): index1 = 0
                        if(index1 > len(lat_ranges)-1): index1 = len(lat_ranges) - 1
                        if(index2 < 0): index2 = 0                                                                                            
                        if(index2 > 1439): index2 = 1439
               
                        #diff = reflectance[i,j] - albedo[i,j]
                        #if(diff<min_diff):
                        #    min_diff = diff
                        #if(diff>max_diff):
                        #    max_diff = diff
               
                        #if(diff<0.2): 
                        #    UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + AI[i,j])/(count[index2,index1]+1)
                        g_NRAD_354[index2, index1] = (g_NRAD_354[index2,index1]*count_NRAD_354[index2,index1] + NRAD[i,j,0])/(count_NRAD_354[index2,index1]+1)
                        g_NRAD_388[index2, index1] = (g_NRAD_388[index2,index1]*count_NRAD_388[index2,index1] + NRAD[i,j,1])/(count_NRAD_388[index2,index1]+1)
                        g_NRAD_500[index2, index1] = (g_NRAD_500[index2,index1]*count_NRAD_500[index2,index1] + NRAD[i,j,2])/(count_NRAD_500[index2,index1]+1)
                        g_REFL_354[index2, index1] = (g_REFL_354[index2,index1]*count_REFL_354[index2,index1] + REFL[i,j,0])/(count_REFL_354[index2,index1]+1)
                        g_REFL_388[index2, index1] = (g_REFL_388[index2,index1]*count_REFL_388[index2,index1] + REFL[i,j,1])/(count_REFL_388[index2,index1]+1)
                        g_SALB_354[index2, index1] = (g_SALB_354[index2,index1]*count_SALB_354[index2,index1] + SALB[i,j,0])/(count_SALB_354[index2,index1]+1)
                        g_UVAI_354[index2, index1] = (g_UVAI_354[index2,index1]*count_UVAI_354[index2,index1] + UVAI[i,j])/(count_UVAI_354[index2,index1]+1)
                        #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + diff)/(count[index2,index1]+1)
                            #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + plotAI[i,j])/(count[index2,index1]+1)
                            #if((ii==1309) and (ii2==59)): print UVAI[index1,index2]
                        #count[index2, index1] = count[index2,index1] + 1
                        count_NRAD_354[index2,index1] = count_NRAD_354[index2,index1] + 1
                        count_NRAD_388[index2,index1] = count_NRAD_388[index2,index1] + 1
                        count_NRAD_500[index2,index1] = count_NRAD_500[index2,index1] + 1
                        count_REFL_354[index2,index1] = count_REFL_354[index2,index1] + 1
                        count_REFL_388[index2,index1] = count_REFL_388[index2,index1] + 1
                        count_SALB_354[index2,index1] = count_SALB_354[index2,index1] + 1
                        count_UVAI_354[index2,index1] = count_UVAI_354[index2,index1] + 1
    
    # Calculate the row-average AI for the secondary plot
#mask_n354_avgs = np.nanmean(np.ma.masked_where(g_NRAD_354[:,:] < -2e5, g_NRAD_354[:,:]),axis=0)
#mask_n388_avgs = np.nanmean(np.ma.masked_where(g_NRAD_388[:,:] < -2e5, g_NRAD_388[:,:]),axis=0)
#mask_n500_avgs = np.nanmean(np.ma.masked_where(g_NRAD_500[:,:] < -2e5, g_NRAD_500[:,:]),axis=0)
#mask_r354_avgs = np.nanmean(np.ma.masked_where(g_REFL_354[:,:] < -2e5, g_REFL_354[:,:]),axis=0)
#mask_r388_avgs = np.nanmean(np.ma.masked_where(g_REFL_388[:,:] < -2e5, g_REFL_388[:,:]),axis=0)

# CALCULATIONS
plot_calc = 1.0 - (1.0 / (g_NRAD_500 - g_NRAD_354))
#plot_calc = 1.0 - (g_NRAD_388 / g_NRAD_354)

#UVAI = plot_calc
#count = count_NRAD_500

#UVAI = g_REFL_354
#count = count_REFL_354

#fout.close()
#print("Min diff = ",min_diff)
#print("Max diff = ",max_diff)
print("Time to plot")
#create an generic lat and lon
#define the lat which will be gridded into 0.25 x 0.25 bins
LATalt = lat_ranges
LONalt = lon_ranges
#LATalt = np.arange(-90,90,0.25)
#LONalt = np.arange(-180,180,0.25)

cmap = plt.cm.jet
#v_min = np.min(UVAI)
#v_max = np.max(UVAI)
#v_min = -1.000  # AI
#v_max = 3.000
#v_min = 0.0  # Albedo
#v_max = 1.0
#v_min = -0.4  # Reflectance - Albedo
#v_max = 1.2
#v_min = 0.0  # VZA
#v_max = 180.0

# Center the colorbar on zero
#norm = MidpointNormalize(midpoint=mid_val,vmin=v_min,vmax=v_max)
###norm = Normalize(vmin=v_min,vmax=v_max)
###mapper = ScalarMappable(norm=norm,cmap=cmap)
###nodatacolor="black"


# Find the values at the comparison points in Greenland
# For file 200804271224
x_index_70 = 1285
y_index_N40 = 9
x_index_80 = 1351
y_index_N50 = 30

#data.close()

plot_lat, plot_lon = np.meshgrid(LATalt,LONalt)
mask_UVAI = np.ma.masked_where(((count_NRAD_354 == 0)), plot_calc)
##!#mask_rad500 = np.ma.masked_where(((count_UVAI_354 == 0)), g_UVAI_354)
##!#mask_rad500 = np.ma.masked_where(((g_SALB_354 > 0.09)), mask_rad500)
##!#mask_rad500 = np.ma.masked_where(((g_REFL_354 > 0.18)), mask_rad500)
##!#mask_rad500   = np.ma.masked_where(((g_REFL_354 < 0.09)), mask_rad500)
##!#mask_UVAI = np.ma.masked_where(((g_UVAI_354 < 0.6)), mask_rad500)

plt.title('OMI Algae ' + plot_time)
#plt.title('NRAD500 - NRAD354')
#plt.title('OMI Reflectivity - Surface Albedo '+plot_time)
mesh = axs[0].pcolormesh(plot_lon, plot_lat,mask_UVAI,transform = datacrs,cmap = colormap,\
        #vmin = -2.0, vmax = 3.0)
        vmin = 0, vmax = 200)
axs[0].set_extent([-180,180,65,90],ccrs.PlateCarree())
#axs[0].set_extent([10,55,65,80],ccrs.PlateCarree())
#axs[0].set_xlim(-0630748.535086173,2230748.438879491)
#axs[0].set_ylim(-2513488.8763307533,0343353.899053069)
#ax.set_xlim(-4170748.535086173,4167222.438879491)
#ax.set_ylim(-2913488.8763307533,2943353.899053069)
cbar = plt.colorbar(mesh,orientation='horizontal',pad=0,\
#cbar = plt.colorbar(mesh,ticks = np.arange(-2.0,4.1,0.5),orientation='horizontal',pad=0,\
#cbar = plt.colorbar(mesh,orientation='horizontal',pad=0,\
    aspect=50,shrink = 0.905,label='Aerosol Index')
    #aspect=50,shrink = 0.905,label='Normalized 500 nm Radiance')

### Draw scatter region box on map
##lonmin,lonmax,latmin,latmax = 16.5,54.5,67.5,76.5
##
##axs[0].plot([lonmin, lonmin, lonmax, lonmax, lonmin],\
##            [latmax, latmin, latmin, latmax, latmax],\
##            transform = ccrs.NorthPolarStereo(central_longitude = 35.))

save = False 
if(save == True):
    out_name = 'omi_single_pass_algae_ai_'+plot_time+'_rows_0to'+str(row_max)+'_zoom.png'       
    plt.savefig(out_name)
    print('Saved image '+out_name)

plt.show()
