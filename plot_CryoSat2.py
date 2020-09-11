#!/usr/bin/env python
"""
  NAME:

  PURPOSE:

  SYNTAX:

  EXAMPLE:

  MODIFICATIONS:
    2020/09/10

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import cartopy.feature as cfeature
import cartopy.crs as ccrs
from scipy.ndimage import gaussian_filter


if(len(sys.argv) < 3):
    print("SYNTAX: python plot_CryoSat2.py file variable")
    sys.exit()

min_dict = {
    'sea_ice_thickness': 0.0,
    'ice_con': 1.
}

max_dict = {
    'sea_ice_thickness': 5.0,
    'ice_con': 100.
}

tick_dict = {
    'sea_ice_thickness': [1,2,3,4,5],
    'ice_con': [1,20,40,60,80,100]
}

tick_label_dict = {
    'sea_ice_thickness': ['1','2','3','4','5'],
    'ice_con': ['1','20','40','60','80','100']
}

infile = sys.argv[1]
variable = sys.argv[2]

#new_time = '20130611t1322_new'
data_path = '/home/bsorenson/data/CryoSat2/'
#new_base_path = '/home/shared/OMAERUV/OMAERUV_Parameters/new_files/'

# Set up projections
colormap = plt.cm.ocean
datacrs = ccrs.PlateCarree()
mapcrs = ccrs.NorthPolarStereo(central_longitude=45.)

# Read in the data
data = Dataset(infile,'r')

# Pull out lat and lon
latitude  = data['lat'][:,:]
longitude = data['lon'][:,:]
ice_thick = data['sea_ice_thickness'][:,:]
ice_conc  = data['ice_con'][:,:]
var_name = data[variable].long_name
data.close()

#ice_thick[ice_thick < -999.] = None
#ice_conc[ice_conc < -999.] = None

# mask bad data
#mask_thick = np.ma.masked_where(ice_thick == None, ice_thick)
#mask_conc  = np.ma.masked_where(ice_conc == None, ice_conc)
mask_thick = np.ma.masked_where(ice_thick < -999., ice_thick)
mask_conc  = np.ma.masked_where(ice_conc < -999., ice_conc)
smooth_thick = gaussian_filter(mask_thick,sigma=1)

#smooth_thick[smooth_thick < 0.0] = None

# Make the figure
plt.close()
ax = plt.axes(projection = mapcrs)
ax.gridlines()
ax.coastlines()
ax.set_extent([-180,180,60,90],datacrs)
if(variable == 'sea_ice_thickness'):
    #mesh = ax.pcolormesh(longitude,latitude,mask_thick,transform = datacrs,\
    mesh = ax.pcolormesh(longitude,latitude,smooth_thick,transform = datacrs,\
            cmap=colormap,vmin=min_dict[variable],vmax=max_dict[variable])
elif(variable == 'ice_con'):
    mesh = ax.pcolormesh(longitude,latitude,mask_conc,transform = datacrs,\
            cmap=colormap,vmin=min_dict[variable],vmax=max_dict[variable])
#CS = ax.contour(longitude,latitude,smooth_thick,[0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0],transform = datacrs)
CS = ax.contour(longitude,latitude,smooth_thick,[0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0],transform = datacrs)

# Adjust and make it look good
ax.add_feature(cfeature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
cbar = plt.colorbar(mesh,ticks = tick_dict[variable],orientation='horizontal',pad=0,aspect=50,label=var_name)
cbar.ax.set_xticklabels(tick_label_dict[variable])
ax.set_xlim(-4170748.535086173,4167222.438879491)
ax.set_ylim(-2913488.8763307533,2943353.899053069)
ax.set_title(infile)
plt.show()

