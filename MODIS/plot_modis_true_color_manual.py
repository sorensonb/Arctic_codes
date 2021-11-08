#!/usr/bin/env python

"""
  NAME:
    plot_modis_true_color_manual.py

  PURPOSE:
    Plot true color imagery from a MODIS L1B data file containing 1-km 
    reflectance data, with the colors calculate manually using the RGB
    channel data in the MODIS file.
  
    NOTE: The MODIS channel and true color functions are designed to work with
    HDF MODIS files retriefed from 
    the data ordering website at this address:
    https://ladsweb.modaps.eosdis.nasa.gov/search/order/1/MODIS:Aqua
    
  MODULES:
    The only non-standard cartopy modules used in this program are:
    - netCDF4
    - cartopy
    - satpy  

  SYNTAX:
    ./plot_modis_true_color_manual.py modis_filename

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>     - 2021/11/08:
      Written (modification of lab_4.py by J. Marquis)

"""
import numpy as np
import sys
from netCDF4 import Dataset
from datetime import datetime,timedelta
from dateutil.relativedelta import relativedelta
from pyhdf import SD
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from glob import glob

if(len(sys.argv) != 2):
    print("SYNTAX: ./plot_modis_true_color_manual.py modis_filename")
    sys.exit()

mapcrs = ccrs.LambertConformal()
datacrs = ccrs.PlateCarree()

# Open the MODIS file
# -------------------
filename = sys.argv[1]
modis = SD.SD(filename)

# Determine whether this is Aqua MODIS or Terra MODIS data
# --------------------------------------------------------
if(filename.strip().split('/')[-1].find('MOD') >= 0):
    sat_name = "Terra"
else:
    sat_name = "Aqua"

# Extract date and time information from the filename
# ---------------------------------------------------
name_split = filename.strip().split('/')[-1].split('.')
year = name_split[1][1:5]
jday = name_split[1][5:]
time = name_split[2]

# Use the extracted information to create a datetime object to
# hold the MODIS overpass time
# ------------------------------------------------------------
base_date = datetime(year=int(year),month=1,day=1)
modis_date = base_date + timedelta(days = int(jday) - 1, \
    hours = int(time[:2]), minutes = int(time[2:]))


# Extract the equator crossing date and time from the file metadata
# -----------------------------------------------------------------
dat = modis.attributes().get('CoreMetadata.0').split()
indx = dat.index('EQUATORCROSSINGDATE')+9
cross_date = dat[indx][1:len(dat[indx])-1]
cross_time = filename.strip().split('/')[-1].split('.')[2]

# Extract the MODIS latitude and longitude grids
# ----------------------------------------------
lat5 = modis.select('Latitude').get()
lon5 = modis.select('Longitude').get()

# Pull out the red, green, and blue reflectance data
# --------------------------------------------------
red   = modis.select('EV_250_Aggr1km_RefSB').get()[0]
green = modis.select('EV_500_Aggr1km_RefSB').get()[1]
blue  = modis.select('EV_500_Aggr1km_RefSB').get()[0]

# The RGB reflectances are on a much higher resolution than the lats and
# lons, so use only every 5th value
# ----------------------------------------------------------------------
red   = red[::5,::5]
green = green[::5,::5]
blue  = blue[::5,::5]

# Extract the scales and offsets for each channel from the file. These are 
# used to convert these reflectance values into usable data
# ------------------------------------------------------------------------
red_scale    = modis.select('EV_250_Aggr1km_RefSB').attributes().get('reflectance_scales')[0]
red_offset   = modis.select('EV_250_Aggr1km_RefSB').attributes().get('reflectance_offsets')[0]
green_scale  = modis.select('EV_500_Aggr1km_RefSB').attributes().get('reflectance_scales')[1]
green_offset = modis.select('EV_500_Aggr1km_RefSB').attributes().get('reflectance_offsets')[1]
blue_scale   = modis.select('EV_500_Aggr1km_RefSB').attributes().get('reflectance_scales')[0]
blue_offset  = modis.select('EV_500_Aggr1km_RefSB').attributes().get('reflectance_offsets')[0]

# Close the satellite file object
# -------------------------------
modis.end()

# Use the scales and offset calibration values to convert from counts to
# reflectance
# ----------------------------------------------------------------------
red   = (red - red_offset) * red_scale
green = (green - green_offset) * green_scale
blue  = (blue - blue_offset) * blue_scale

red   = red*(255./1.1) 
green = green*(255./1.1) 
blue  = blue*(255./1.1) 

# Create color scales for each RGB channel
# ----------------------------------------
old_val = np.array([0,30,60,120,190,255])
ehn_val = np.array([0,110,160,210,240,255])
red     = np.interp(red,old_val,ehn_val) / 255.
old_val = np.array([0,30,60,120,190,255])
ehn_val = np.array([0,110,160,200,230,240])
green   = np.interp(green,old_val,ehn_val) / 255.
old_val = np.array([0,30,60,120,190,255])
ehn_val = np.array([0,100,150,210,240,255])
blue    = np.interp(blue,old_val,ehn_val) / 255.

# Combine the three RGB channels into 1 3-d array
# -----------------------------------------------
image = np.zeros((red.shape[0],red.shape[1],3))
image[:,:,0] = red
image[:,:,1] = green
image[:,:,2] = blue

# Convert the color values into a format usable for plotting
# ----------------------------------------------------------
colortuple = tuple(np.array([image[:,:,0].flatten(), \
    image[:,:,1].flatten(), image[:,:,2].flatten()]).transpose().tolist())

# Set up the figure
# -----------------
plt.close('all') 
fig1 = plt.figure()
ax = plt.axes(projection = mapcrs)

# Mask any missing values in the data
# -----------------------------------
image = np.ma.masked_where(np.isnan(image),image)

# Plot the data on the image using the calculated colors
# ------------------------------------------------------
ax.pcolormesh(lon5,lat5,image[:,:,0],color= colortuple, \
    shading='auto', transform = ccrs.PlateCarree()) 

# Add state boundaries, nation boundaries, and coastlines
# -------------------------------------------------------
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.STATES)
ax.coastlines()

zoom = False
if(zoom):
    ax.set_extent([lon_lims[0], lon_lims[1], lat_lims[0], lat_lims[1]], \
        datacrs)
    zoom_add = '_zoom'
else:
    zoom_add = ''

# Add a title to the figure
# -------------------------
ax.set_title(sat_name + ' MODIS ' + modis_date.strftime('%Y-%m-%d %H:%M'))

# Save the figure
# ---------------
outname = 'modis_' + sat_name + '_true_color_' + \
    modis_date.strftime('%Y%m%d%H%M') + zoom_add + '_manual.png'
plt.savefig(outname,dpi=300)
print("Saved image",outname)
plt.show()
