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
import cartopy.crs as ccrs
##import matplotlib.colors as colors
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

if(len(sys.argv)<4):
    print("SYNTAX: python plot_single_OMI.py date variable row_max")
    sys.exit()

plot_time = sys.argv[1]
variable = sys.argv[2]
row_max=int(sys.argv[3])
channel_idx = 0
str_wave = ''

# Check which channels are being used
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

#new_time = '20130611t1322_new'
base_path = '/home/bsorenson/data/OMI/H5_files/'
#new_base_path = '/home/shared/OMAERUV/OMAERUV_Parameters/new_files/'

#define the directory to the aerosol index files
##pathlength = '/home/shared/OMAERUV/OMAERUV_Parameters/Aerosol_Index/'
##pathlength1 = '/home/shared/OMAERUV/OMAERUV_Parameters/Latitude/'
##pathlength2 = '/home/shared/OMAERUV/OMAERUV_Parameters/Longitude/'
##pathlength3 = '/home/shared/OMAERUV/OMAERUV_Parameters/X_Track_Quality_Flags/'

###outfilename = 'omi_single_pass_'+plot_time+'_values_80.txt'
###fout = open(outfilename,'w')

n_p = 1440
nl = 720

#infile = sys.argv[1]

#manually enter the longitude limits
lonmin = -180
lonmax = 180

latmax =  90

#date = sys.argv[1]

# Set up the dictionary
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# The dictionary is set up so that each lat and lon pair are used as the
# keys. For example, key '48x98' contains the data for latitude ranges
# 48.0-49.0 and longitude ranges of 98.0-99.0
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

latmin = 60 

# Set up values for gridding the AI data
lat_gridder = latmin * 4.

lat_ranges = np.arange(latmin,90.1,0.25)
lon_ranges = np.arange(-180,180.1,0.25)

### Look through the directories for data with the desired date
##ailist      = glob.glob(pathlength+'UVAerosolIndex_'+date+'*.bin')
##latlist     = glob.glob(pathlength1+'Latitude_'+ date+'*.bin')
##lonlist     = glob.glob(pathlength2+'Longitude_'+ date+'*.bin')
##xtracklist  = glob.glob(pathlength3+'XTrackQualityFlags_'+date+'*.bin')
##num_files   = len(ailist)

mapcrs = ccrs.NorthPolarStereo()
datacrs = ccrs.PlateCarree()
colormap = plt.cm.jet

NorthAmerica = True
# Set up the polar stereographic projection map
fig1, axs = plt.subplots(2, 1, figsize=(8,8))
#fig1.set_size_inches(8,8)
#fig1 = plt.figure(figsize=(8,8))
if(latmin<45):
    ax = plt.axes(projection = ccrs.Miller())
    #m = Basemap(projection='mill',lon_0=0,resolution='l')
else:
    axs[0] = plt.axes(projection = ccrs.NorthPolarStereo(central_longitude = 0.))
    #ax = plt.axes(projection = ccrs.NorthPolarStereo(central_longitude = 45.))
    #m = Basemap(projection='npstere',boundinglat=latmin,lon_0=0,resolution='l')
#if(NorthAmerica is True):
#    m = Basemap(projection='lcc',width=12000000,height=9000000,\
#                rsphere=(6378137.00,6356752.3142),\
#                resolution='l',area_thresh=1000.,\
#                lat_1=45.,lat_2=55.,lat_0=50.,lon_0=-107.)
#fig = plt.gcf()
axs[0].gridlines()
axs[0].coastlines(resolution = '50m')
#m.drawcoastlines()
#m.drawparallels(np.arange(-80.,81.,20.))
#m.drawmeridians(np.arange(-180.,181.,20.))

#ax = plt.axes(projection=ccrs.Miller())
#ax.coastlines()
##ax.add_feature(cfeature.BORDERS)



RED   = np.array([  0., 30.,135.,  0.,  0.,255.,255.,200.,255.])
GREEN = np.array([  0.,144.,206.,255.,255.,255.,  0.,  0.,255.])
BLUE  = np.array([255.,255.,250.,255.,  0.,  0.,  0.,  0.,255.])

RED   = RED/256.
GREEN = GREEN/256.
BLUE  = BLUE/256.

ai_range=np.array([-1.0,-0.5,0.0,0.5,1.0,2.0,2.5,3.0,4.0])
val = ['-1.0','-0.5','0.0','0.5','1.0','2.0','2.5','3.0','4.0']
color_num=9
# Calculate average and standard deviations
#newch1=np.zeros(shape=(len(lat_ranges),len(lon_ranges)))

UVAI = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
count = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
newch1=np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
#newch1=np.zeros(shape=(1440,720))

#ai_list = subprocess.check_output('ls /home/shared/OMAERUV/'+             \
#          'OMAERUV_Parameters_20190212_download/Aerosol_Index/UVAerosolIndex_'+plot_time+'*.bin',\
#          shell=True).decode('utf-8').strip().split('\n')
#lat_list = subprocess.check_output('ls /home/shared/OMAERUV/'+             \
#          'OMAERUV_Parameters_20190212_download/Latitude/Latitude_'+lat_time+'*.bin',\
#          shell=True).decode('utf-8').strip().split('\n')
#lon_list = subprocess.check_output('ls /home/shared/OMAERUV/'+             \
#          'OMAERUV_Parameters_20190212_download/Longitude/Longitude_'+plot_time+'*.bin',\
#          shell=True).decode('utf-8').strip().split('\n')
#x_list = subprocess.check_output('ls /home/shared/OMAERUV/'+             \
#          'OMAERUV_Parameters_20190212_download/X_Track_Quality_Flags/XTrackQualityFlags_'+plot_time+'*.bin',\
#          shell=True).decode('utf-8').strip().split('\n')
###ai_list = subprocess.check_output('ls /home/shared/OMAERUV/'+             \
###          'OMAERUV_Parameters_20190212_download/jzhang_blake/UVAerosolIndex_'+plot_time+'*.bin',\
###          shell=True).decode('utf-8').strip().split('\n')
###albedo_list = subprocess.check_output('ls /home/shared/OMAERUV/'+             \
###          'OMAERUV_Parameters_20190212_download/jzhang_blake/UVAerosolIndex_'+plot_time+'*.bin',\
###          shell=True).decode('utf-8').strip().split('\n')
###lat_list = subprocess.check_output('ls /home/shared/OMAERUV/'+             \
###          'OMAERUV_Parameters_20190212_download/jzhang_blake/Latitude_'+plot_time+'*.bin',\
###          shell=True).decode('utf-8').strip().split('\n')
###lon_list = subprocess.check_output('ls /home/shared/OMAERUV/'+             \
###          'OMAERUV_Parameters_20190212_download/jzhang_blake/Longitude_'+plot_time+'*.bin',\
###          shell=True).decode('utf-8').strip().split('\n')
###x_list = subprocess.check_output('ls /home/shared/OMAERUV/'+             \
###          'OMAERUV_Parameters_20190212_download/jzhang_blake/XTrackQualityFlags_'+plot_time+'*.bin',\
###          shell=True).decode('utf-8').strip().split('\n')


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
    for i in range(PDATA.shape[0]):
        for j in range(0,row_max):
            #if((albedo[i,j]>-20) & (reflectance[i,j]>-20)):
            if(PDATA[i,j]>-2e5):
            #if(plotAI[i,j]>-20):
                # Only plot if XTrack flag is met
                if((XTRACK[i,j] == 0) | (XTRACK[i,j] == 4)):
                    # Print values to text file
                    if(LAT[i,j] > latmin):
                        counter+=1
#                        fout.write("{0:.6f} {1:.6f} {2:.6f} {3:.6f} {4:.6f} {5:.6f} {6:.6f} {7:.6f} {8:.6f} {9:.6f} {10:.6f} {11:.6f}\n".format(\
#                            LAT[i,j],LON[i,j],AI[i,j],0.5,SZA[i,j],VZA[i,j],RAZ[i,j], \
#                            ALBEDO[i,j,0],ALBEDO[i,j,1],REFLECTANCE[i,j,0],\
#                            REFLECTANCE[i,j,1],CLD[i,j]))


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
                        UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + PDATA[i,j])/(count[index2,index1]+1)
                        #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + diff)/(count[index2,index1]+1)
                            #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + plotAI[i,j])/(count[index2,index1]+1)
                            #if((ii==1309) and (ii2==59)): print UVAI[index1,index2]
                        count[index2, index1] = count[index2,index1] + 1

# Calculate the row-average AI for the secondary plot
mask_avgs = np.nanmean(np.ma.masked_where(PDATA[:,:] < -2e5, PDATA[:,:]),axis=0)

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
v_min = -1.000  # AI
v_max = 3.000
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

##print("For point 70 N, -40 E")
##print("LAT   = ",LAT[x_index_70,y_index_N40])
##print("LON   = ",LON[x_index_70,y_index_N40])
##print("AI    = ",AI[x_index_70,y_index_N40])
##print("CLD   = ",CLD[x_index_70,y_index_N40])
##print("ALB   = ",ALBEDO[x_index_70,y_index_N40,0])
##print("REFL  = ",REFLECTANCE[x_index_70,y_index_N40,0])
##print("VZA   = ",VZA[x_index_70,y_index_N40])
##print("SZA   = ",SZA[x_index_70,y_index_N40])
##print("XTRCK = ",XTRACK[x_index_70,y_index_N40])
##print("PIXEL = ",PIXEL[x_index_70,y_index_N40,0])
##print("MSMNT = ",MSMNT[x_index_70])
##print("GRND  = ",bin(GRND[x_index_70,y_index_N40]))
##print(" ")
##print("For point 80 N, -50 E")
##print("LAT   = ",LAT[x_index_80,y_index_N50])
##print("LON   = ",LON[x_index_80,y_index_N50])
##print("AI    = ",AI[x_index_80,y_index_N50])
##print("CLD   = ",CLD[x_index_80,y_index_N50])
##print("ALB   = ",ALBEDO[x_index_80,y_index_N50,0])
##print("REFL  = ",REFLECTANCE[x_index_80,y_index_N50,0])
##print("VZA   = ",VZA[x_index_80,y_index_N50])
##print("SZA   = ",SZA[x_index_80,y_index_N50])
##print("XTRCK = ",XTRACK[x_index_80,y_index_N50])
##print("PIXEL = ",PIXEL[x_index_80,y_index_N50,0])
##print("MSMNT = ",MSMNT[x_index_80])
##print("GRND  = ",bin(GRND[x_index_80,y_index_N50]))
##print("counter = ",counter)
data.close()

## Convert integer to binary
#bin_val = bin(GRND[x_index_80,y_index_N50])
## Convet binary value back to integer
#int_val = int(bin_val,2)
 
#data.close()

plot_lat, plot_lon = np.meshgrid(LATalt,LONalt)
mask_UVAI = np.ma.masked_where(count == 0, UVAI)

plt.title('OMI ' + variable + str_wave + ' '+plot_time)
#plt.title('OMI Reflectivity - Surface Albedo '+plot_time)
mesh = axs[0].pcolormesh(plot_lon, plot_lat,mask_UVAI,transform = datacrs,cmap = colormap,\
        vmin = var_dict[variable]['min'], vmax = var_dict[variable]['max'])
axs[0].set_extent([-180,180,latmin,90],ccrs.PlateCarree())
axs[0].set_xlim(-3430748.535086173,3430748.438879491)
axs[0].set_ylim(-3413488.8763307533,3443353.899053069)
#ax.set_xlim(-4170748.535086173,4167222.438879491)
#ax.set_ylim(-2913488.8763307533,2943353.899053069)
cbar = plt.colorbar(mesh,ticks = np.arange(-2.0,4.1,0.5),orientation='horizontal',pad=0,\
    aspect=50,shrink = 0.905,label=variable)

save = False
if(save == True):
    out_name = 'omi_single_pass_'+name_dict[variable] + str_wave + '_'+plot_time+'_rows_0to'+str(row_max)+'.png'       
    plt.savefig(out_name)
    print('Saved image '+out_name)

plt.show()
