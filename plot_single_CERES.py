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
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
import matplotlib.colors as color
##import matplotlib.colors as colors
from matplotlib.colors import rgb2hex,Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.colorbar import ColorbarBase
import subprocess
from netCDF4 import Dataset

if(len(sys.argv)<3):
    print("SYNTAX: python plot_single_OMI.py date row_max")
    sys.exit()

plot_time = sys.argv[1]
#new_time = '20130611t1322_new'
base_path = '/home/bsorenson/data/CERES/SSFLevel2/'

n_p = 1440
nl = 720

#infile = sys.argv[1]

#manually enter the longitude limits
lonmin = -180
lonmax = 180
latmax =  90

lat_ranges = np.arange(-90,90,1.0)
lon_ranges = np.arange(-180,180,1.0)

latmin = -89
NorthAmerica = True
# Set up the polar stereographic projection map
fig1 = plt.figure(figsize=(8,8))
if(latmin<45):
    m = Basemap(projection='mill',lon_0=0,resolution='l')
else:
    m = Basemap(projection='npstere',boundinglat=latmin,lon_0=0,resolution='l')
if(NorthAmerica is True):
    m = Basemap(projection='lcc',width=12000000,height=9000000,\
                rsphere=(6378137.00,6356752.3142),\
                resolution='l',area_thresh=1000.,\
                lat_1=45.,lat_2=55.,lat_0=50.,lon_0=-107.)
fig = plt.gcf()
m.drawcoastlines()
m.drawparallels(np.arange(-80.,81.,20.))
m.drawmeridians(np.arange(-180.,181.,20.))

#ax = plt.axes(projection=ccrs.Miller())
#ax.coastlines()
##ax.add_feature(cfeature.BORDERS)


year = plot_time[:4]
date = plot_time[4:8]
if(len(plot_time)==13):
    time = plot_time[9:]
elif(len(plot_time)==12):
    time = plot_time[8:]
else:
    time = ''
total_list = [base_path+'CERES_SSF_Terra-XTRK_Edition4A_Subset_2015062900-2015063011.nc']
#total_list = subprocess.check_output('ls '+base_path+'OMI-Aura_L2-OMAERUV_'+year+'m'+date+'t'+time+'*.he5',\
#          shell=True).decode('utf-8').strip().split('\n')

##total_list = ['/home/shared/OMAERUV/H5_files_20190212_download/OMI-Aura_L2-OMAERUV_2008m0427t0052-o20122_v003-2017m0721t120210.he5',\
##              '/home/shared/OMAERUV/H5_files_20190212_download/OMI-Aura_L2-OMAERUV_2008m0427t0231-o20123_v003-2017m0721t120217.he5',\
##              '/home/shared/OMAERUV/H5_files_20190212_download/OMI-Aura_L2-OMAERUV_2008m0427t0410-o20124_v003-2017m0721t120234.he5',\
##              '/home/shared/OMAERUV/H5_files_20190212_download/OMI-Aura_L2-OMAERUV_2008m0427t0549-o20125_v003-2017m0721t121042.he5']

min_diff = 30
max_diff = -30

swf_grid = np.zeros(shape=(n_p,nl))
count = np.zeros(shape=(n_p,nl))
newch1=np.zeros(shape=(1440,720))
lat_checkmin = 35.
lat_checkmax = 70.
lon_checkmin = -140.
lon_checkmax = -80.
for fileI in range(len(total_list)):
    # read in data directly from HDF5 files
    print(total_list[fileI])
    data = Dataset(total_list[fileI],'r')
    lat   = 90.-data.variables['Colatitude_of_CERES_FOV_at_surface'][:]
    lon   = data.variables['Longitude_of_CERES_FOV_at_surface'][:]
    lon[lon>179.99] = -360.+lon[lon>179.99]
    swf   = data.variables['CERES_SW_TOA_flux___upwards'][:]

    counter = 0
    # Loop over the values and rows 
    for i in range(swf.shape[0]):
        #if((albedo[i,j]>-20) & (reflectance[i,j]>-20)):
        if(swf[i]<5000):
            if(((lat[i]>lat_checkmin) & (lat[i]<lat_checkmax))  & \
               ((lon[i]>lon_checkmin) & (lon[i]<lon_checkmax))):
                index1 = int(np.floor(lat[i]*4 + 360.))
                index2 = int(np.floor(lon[i]*4 + 720.))
                
                if(index1 < 0): index1 = 0
                if(index1 > 719): index1 = 719
                if(index2 < 0): index2 = 0                                                                                            
                if(index2 > 1439): index2 = 1439
          
                swf_grid[index2, index1] = (swf_grid[index2,index1]*count[index2,index1] + swf[i])/(count[index2,index1]+1)
                #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + diff)/(count[index2,index1]+1)
                    #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + plotAI[i,j])/(count[index2,index1]+1)
                    #if((ii==1309) and (ii2==59)): print UVAI[index1,index2]
                count[index2, index1] = count[index2,index1] + 1
    data.close()
print("Time to plot")
#create an generic lat and lon
#define the lat which will be gridded into 0.25 x 0.25 bins
LATalt = np.arange(-90,90,0.25)
LONalt = np.arange(-180,180,0.25)

cmap = plt.get_cmap('jet')
v_min = 0.000  # AI
v_max = 600.000
#v_min = 0.0  # Albedo
#v_max = 1.0
#v_min = -0.4  # Reflectance - Albedo
#v_max = 1.2
#v_min = 0.0  # VZA
#v_max = 180.0

# Center the colorbar on zero
#norm = MidpointNormalize(midpoint=mid_val,vmin=v_min,vmax=v_max)
norm = Normalize(vmin=v_min,vmax=v_max)
mapper = ScalarMappable(norm=norm,cmap=cmap)
nodatacolor="black"

for ii in range(0,n_p-1):
    for jj in range(0,nl-1):
        if(count[ii,jj]>0):
            colors = mapper.to_rgba(swf_grid[ii,jj])
            
            lon0 = LONalt[ii]
            lat0 = LATalt[jj]
            lat1 = LATalt[jj+1]
            lon1 = LONalt[ii+1]
       
            if(lat1>=latmin):
     
                y = [lat0,lat1,lat1,lat0]
                x = [lon0,lon0,lon1,lon1]
                mx, my = m(x,y)
                mxy = zip(mx,my)
                pair1 = (mx[0],my[0])
                pair2 = (mx[1],my[1])
                pair3 = (mx[2],my[2])
                pair4 = (mx[3],my[3])
                #mxy = zip(mx,my)
                # Plot the box on the map using color
                #colors = [RED[index_val],GREEN[index_val],BLUE[index_val]]
                color2 = rgb2hex(colors)
                poly = Polygon([pair1,pair2,pair3,pair4],facecolor=color2,edgecolor=color2)
                plt.gca().add_patch(poly)
                #if((i==1309) and (j==59)): print x, y
                #print x, y
                #if(fi == 0) then begin 
        

print("Generating Plot")
plt.title('CERES SWF '+plot_time)
#plt.title('OMI Reflectivity - Surface Albedo '+plot_time)
cax = fig.add_axes([0.16,0.075,0.7,0.025])
cb = ColorbarBase(cax,cmap=cmap,norm=norm,orientation='horizontal')
cb.ax.set_xlabel('SWF')
#cb.ax.set_xlabel('Reflectivity - Surface Albedo')
#out_name = 'omi_single_pass_ai_200804270052_to_0549_composite_rows_0to'+str(row_max)+'.png'       
out_name = 'ceres_single_pass_swf_'+plot_time+'.png'       
#out_name = 'omi_single_pass_refl_albedo_diff_'+plot_time+'_rows_0to'+str(row_max)+'.png'       
plt.savefig(out_name,dpi=300)
print('Saved image '+out_name)
plt.show()
