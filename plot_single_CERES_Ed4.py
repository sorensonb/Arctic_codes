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
##import matplotlib.colors as colors
import cartopy.crs as ccrs
import subprocess
from datetime import datetime
from dateutil.relativedelta import relativedelta
from netCDF4 import Dataset

var_dict = {
    'SWF': 'CERES_SW_TOA_flux___upwards',
    'LWF': 'CERES_LW_TOA_flux___upwards'
}

max_dict = {
    'SWF': 600.,
    'LWF': 280.
}

min_dict = {
    'SWF': 0.,
    'LWF': 100.
}

if(len(sys.argv)<3):
    print("SYNTAX: python plot_single_CERES_Ed4.py date variable")
    print("        date: YYYYMMDDHH")
    sys.exit()

plot_time = sys.argv[1]
invar     = sys.argv[2]

# Set up datetime objects to use when selecting the correct data
day = False
day_adder = ''
if(len(plot_time) == 10):
    str_fmt = "%Y%m%d%H"
    hour_adder = 1
elif(len(plot_time) == 8):
    print("Daily average")
    day = True
    str_fmt = "%Y%m%d"
    hour_adder = 23
    day_adder = 'daily_'
else:
    print("INVALID PLOT TIME")
    sys.exit()

start_date = datetime.strptime(plot_time,str_fmt) - relativedelta(hours = hour_adder)
end_date   = start_date + relativedelta(hours = hour_adder)

base_path = '/home/bsorenson/data/CERES/SSFLevel2/'

n_p = 1440
nl = 720

#infile = sys.argv[1]

#manually enter the longitude limits
lonmin = -180
lonmax = 180
latmax =  90
latmin = 60 

# Set up values for gridding the AI data
lat_gridder = latmin * 4.

lat_ranges = np.arange(latmin,90.1,0.25)
lon_ranges = np.arange(-180,180.1,0.25)

mapcrs = ccrs.NorthPolarStereo()
datacrs = ccrs.PlateCarree()
colormap = plt.cm.jet
# Set up the polar stereographic projection map
fig = plt.figure(figsize=(8,8))
if(latmin<45):
    ax = plt.axes(projection = ccrs.Miller())
    #m = Basemap(projection='mill',lon_0=0,resolution='l')
else:
    ax = plt.axes(projection = ccrs.NorthPolarStereo(central_longitude = 0.))

ax.gridlines()
ax.coastlines(resolution = '50m')


year = plot_time[:4]
date = plot_time[4:8]
if(len(plot_time)==13):
    time = plot_time[9:]
elif(len(plot_time)==12):
    time = plot_time[8:]
else:
    time = ''
total_list = [base_path+'CERES_SSF_Aqua-XTRK_Edition4A_Subset_2008042200-2008042223.nc']
total_list = subprocess.check_output('ls '+base_path+'CERES_SSF_Aqua-XTRK_Edition4A_Subset_'+year+date+'*.nc',\
          shell=True).decode('utf-8').strip().split('\n')

##total_list = ['/home/shared/OMAERUV/H5_files_20190212_download/OMI-Aura_L2-OMAERUV_2008m0427t0052-o20122_v003-2017m0721t120210.he5',\
##              '/home/shared/OMAERUV/H5_files_20190212_download/OMI-Aura_L2-OMAERUV_2008m0427t0231-o20123_v003-2017m0721t120217.he5',\
##              '/home/shared/OMAERUV/H5_files_20190212_download/OMI-Aura_L2-OMAERUV_2008m0427t0410-o20124_v003-2017m0721t120234.he5',\
##              '/home/shared/OMAERUV/H5_files_20190212_download/OMI-Aura_L2-OMAERUV_2008m0427t0549-o20125_v003-2017m0721t121042.he5']

swf_grid = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
count = np.zeros(shape=(len(lon_ranges),len(lat_ranges)))
base_date = datetime(year=1970,month=1,day=1)

for fileI in range(len(total_list)):
    # read in data directly from HDF5 files
    print(total_list[fileI])
    data = Dataset(total_list[fileI],'r')
    lat   = 90. - data.variables['Colatitude_of_CERES_FOV_at_surface'][:]
    lon   = data.variables['Longitude_of_CERES_FOV_at_surface'][:]
    lon[lon>179.99] = -360.+lon[lon>179.99]
    flux  = data.variables[var_dict[invar]][:]
    time  = data.variables['time'][:]

    counter = 0
    # Loop over the values and rows 
    for i in range(flux.shape[0]):
        #if((albedo[i,j]>-20) & (reflectance[i,j]>-20)):
        local_time = base_date + relativedelta(days = time[i])
        if((local_time >= start_date) & (local_time < end_date)):
            if((flux[i] < 5000) and (flux[i] > 0)):
                index1 = int(np.floor(lat[i]*4 - lat_gridder))
                index2 = int(np.floor(lon[i]*4 + 720.))
                
                if(index1 < 0): index1 = 0
                if(index1 > 719): index1 = 719
                if(index2 < 0): index2 = 0                                                                                            
                if(index2 > 1439): index2 = 1439
          
                swf_grid[index2, index1] = (swf_grid[index2,index1]*count[index2,index1] + flux[i])/(count[index2,index1]+1)
                #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + diff)/(count[index2,index1]+1)
                #UVAI[index2, index1] = (UVAI[index2,index1]*count[index2,index1] + plotAI[i,j])/(count[index2,index1]+1)
                #if((ii==1309) and (ii2==59)): print UVAI[index1,index2]
                count[index2, index1] = count[index2,index1] + 1
    data.close()
print("Time to plot")
#create an generic lat and lon
#define the lat which will be gridded into 0.25 x 0.25 bins
LATalt = lat_ranges
LONalt = lon_ranges

cmap = plt.cm.jet
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
nodatacolor="black"

plot_lat, plot_lon = np.meshgrid(LATalt,LONalt)
mask_flux = np.ma.masked_where(count == 0, swf_grid)

print("Generating Plot")
plt.title('CERES ' + invar + ' ' +plot_time)
#plt.title('OMI Reflectivity - Surface Albedo '+plot_time)
mesh = ax.pcolormesh(plot_lon, plot_lat,mask_flux,transform = datacrs,cmap = colormap,\
        vmin = min_dict[invar], vmax = max_dict[invar])
ax.set_extent([-180,180,latmin,90],ccrs.PlateCarree())
        #vmin = var_dict[variable]['min'], vmax = var_dict[variable]['max'])
cbar = plt.colorbar(mesh,orientation='horizontal',pad=0,\
    aspect=50,shrink = 0.905,label=invar)

save = True 
if(save == True):
    out_name = 'ceres_single_pass_' + day_adder + invar.lower() +'_'+plot_time+'.png'       
    plt.savefig(out_name,dpi=300)
    print('Saved image '+out_name)
else:
    plt.show()
