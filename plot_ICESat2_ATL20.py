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
from matplotlib.widgets import Slider
import h5py
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime

##!#mapcrs = ccrs.NorthPolarStereo(central_longitude = 45.)
##!#datacrs = ccrs.PlateCarree()
##!#colormap = plt.cm.jet
##!#
##!#def main(infile):
##!#    # Read in the ICESat2 data
##!#    data = h5py.File(infile,'r')
##!#    # Assume daily data are being analyzed
##!#    ex = Explorer(data)
##!#    ex.show() 
##!#
##!#class Explorer(object):
##!#    def __init__(self,data):
##!#        nt = len(data['daily'].keys())
##!#
##!#        all_data = np.zeros((nt,448,304))
##!#
##!#        # Insert ice freeboard data into array
##!#        for xi in range(31):
##!#            in_key = 'day'+str(xi+1).zfill(2)
##!#            all_data[xi,:,:] = data['daily/'+in_key+'/mean_fb'][:,:]
##!#        
##!#        mask_data = np.ma.masked_where(all_data > 1000, all_data)
##!#        self.data = mask_data
##!#
##!#        # Pull out lat and lon
##!#        latitude  = data['grid_lat'][:,:]
##!#        longitude = data['grid_lon'][:,:]
##!#
##!#        # Flip the longitudes
##!#        longitude = longitude + 180.
##!#        longitude[longitude > 180.] = longitude[longitude > 180.] - 360.
##!#   
##!#        self.latitude = latitude
##!#        self.longitude = longitude
##!# 
##!#        self.fig = plt.figure()
##!#        self.ax = plt.subplot(projection = mapcrs)
##!#        self.ax.gridlines()
##!#        self.ax.coastlines()
##!#        self.ax.set_extent([-180,180,60,90])
##!#
##!#        self.ax.set_title('Daily ICESat2 Sea Ice Freeboard')
##!#
##!#        self.sliderax = self.fig.add_axes([0.2,0.02,0.65,0.04])
##!#        self.slider = Slider(self.sliderax, 'Day',1,nt,valinit=1,valstep = 1)
##!#        self.slider.on_changed(self.update)
##!#
##!#        self.im = self.ax.pcolormesh(longitude,latitude,mask_data[0,:,:],\
##!#            transform = datacrs, cmap = colormap)
##!#
##!#        cbar = plt.colorbar(self.im,ticks = np.arange(0,1.01,0.1),orientation='horizontal',pad=0,\
##!#            aspect=100,shrink = 1.0,label='Sea Ice Freeboard [m]')
##!#        # Adjust and make it look good
##!#        #self.ax.add_feature(cfeature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
##!#        self.ax.set_xlim(-4170748.535086173,4167222.438879491)
##!#        self.ax.set_ylim(-2913488.8763307533,2943353.899053069)
##!#
##!#    def update(self,i):
##!#        self.im.set_array(self.data[int(i-1),:,:].ravel())
##!#        #self.im = self.ax.pcolormesh(self.longitude,self.latitude,self.data[int(i-1),:,:],\
##!#        #    transform = datacrs, cmap = colormap)
##!#        self.fig.canvas.draw_idle()
##!#
##!#    def show(self):
##!#        plt.show()

# NOTE: Data located in /home/bsorenson/data/ICESat2/ATL20/
#       Downloaded from nsidc.org/data/atl20
if(len(sys.argv)<3):
    print("SYNTAX: python plot_ICESat2.py file day")
    print("       day: 1,2,3, etc.")
    print("       if set to 'month', will plot monthly data instead")
    sys.exit()


##!#var_dict = {
##!#    'freeboard': 'gt2r/sea_ice_segments/heights/height_segment_height',
##!#}
##!#
##!#max_dict = {
##!#    'sea_ice_height': 5.0,
##!#    'sea_ice_conc': 102.
##!#}

# Grab command line arguments
infile = sys.argv[1]
date_key = sys.argv[2]

# Extract date info from title
str_year  = infile.strip().split('_')[1][:4]
str_month = infile.strip().split('_')[1][4:6]
##!#main(infile)

#new_time = '20130611t1322_new'
data_path = '/home/bsorenson/data/ICESat2/'

# Read in the data
data = h5py.File(infile,'r')

# Pull out lat and lon
latitude  = data['grid_lat'][:,:]
longitude = data['grid_lon'][:,:]

# Flip the longitudes
longitude = longitude + 180.
longitude[longitude > 180.] = longitude[longitude > 180.] - 360.

# Get the data key from the command line
if(date_key == 'month'):
    total_key = 'monthly/mean_fb'
    dtm_date = datetime(int(str_year),int(str_month),1)
    plot_title = 'Monthly Averaged Sea Ice Freeboard\n'+dtm_date.strftime('%b %Y')
    file_saver = '_month'
else:
    total_key = 'daily/day'+str(int(date_key)+1).zfill(2)+'/mean_fb'
    dtm_date = datetime(int(str_year),int(str_month),int(date_key)+1)
    plot_title = 'Daily Sea Ice Freeboard\n'+dtm_date.strftime("%d %b %Y")
    file_saver = '_day'

plot_data = data[total_key][:,:]
#plot_data = data['monthly/mean_fb'][:,:]
#plot_data = data[var_dict[variable]][:]
data.close()

# mask bad data
mask_data = np.ma.masked_where(plot_data > 1000, plot_data)

mapcrs = ccrs.NorthPolarStereo(central_longitude = 0.)
datacrs = ccrs.PlateCarree()
colormap = plt.cm.jet

plt.close()
ax = plt.axes(projection = mapcrs)
ax.gridlines()
ax.coastlines()
ax.set_extent([-180,180,60,90])
mesh = ax.pcolormesh(longitude,latitude,mask_data,\
        transform = datacrs, cmap = colormap\
       , vmin = 0.0, vmax = 1.0)
#cbar = plt.colorbar(mesh,ticks = tick_dict[variable],orientation='horizontal',pad=0,\
cbar = plt.colorbar(mesh,ticks = np.arange(0,1.01,0.1),orientation='horizontal',pad=0,\
    aspect=50,shrink = 0.905,label='Sea Ice Freeboard [m]')
#CS = ax.contour(longitude,latitude,smooth_thick,[0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0],transform = datacrs)

# Adjust and make it look good
ax.add_feature(cfeature.LAND,zorder=100,edgecolor='darkgrey',facecolor='darkgrey')
ax.set_xlim(-4170748.535086173,4167222.438879491)
ax.set_ylim(-2913488.8763307533,2943353.899053069)

# Set the plot title
ax.set_title(plot_title)

# Save image, if desired
save = True 
if(save == True):
    outname = 'icesat2_freeboard'+file_saver+'_'+dtm_date.strftime('%Y%m%d')+'.png'
    plt.savefig(outname,dpi=300)
    print("Saved image ",outname)
else:
    plt.show()

###
###for i,avg_arr in enumerate(split_data):
###    if(avg_arr.count() > 0):
###        avg_split[i] = np.nanmean(avg_arr)
###
####[np.nanmean(avg_arr) for avg_arr in np.array_split(test_mask,len(test_mask)/10)]
###
###
#### Make the figure
###plt.close()
###plt.figure()
####plt.scatter(np.arange(len(plot_data)),plot_data,s=4)
###plt.scatter(np.arange(len(avg_split)),avg_split,s=4)
###plt.ylim(np.min(plot_data)*0.98,max_dict[variable])
###plt.title(variable)
###plt.grid()
###plt.show()
###
###data.close()
