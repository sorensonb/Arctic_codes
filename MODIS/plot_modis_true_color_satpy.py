#!/usr/bin/env python

"""
  NAME:
    plot_modis_true_color_satpy.py

  PURPOSE:
    Plot true color imagery from a MODIS L1B data file containing 1-km 
    reflectance data, with the colors and mapping information calculated
    automatically by the satpy module
 
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
    ./python plot_modis_true_color_satpy.py modis_filename

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>   - 2021/11/08:
      Written

"""
import numpy as np
import sys
from netCDF4 import Dataset
from datetime import datetime,timedelta
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from satpy.scene import Scene
from satpy.writers import get_enhanced_image
from glob import glob

mapcrs = ccrs.LambertConformal()
datacrs = ccrs.PlateCarree()

if(len(sys.argv) != 2):
    print("SYNTAX: ./plot_modis_true_color_satpy.py modis_filename")
    sys.exit()

# Determine the correct MODIS file associated with the date
# ---------------------------------------------------------
filename = sys.argv[1]

print(filename)
day_filenames = glob(filename)
cmpst_add = ''

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

# Use satpy (Scene) to open the file
# ----------------------------------
scn = Scene(reader = 'modis_l1b', filenames = day_filenames)

# Load true-color data 
# --------------------
scn.load(['true_color'])

# Save the unedited, uncorrected, raw swath image
# -----------------------------------------------
scn.save_dataset('true_color','test_image_true3.png')

# Set up latitude and longitude bounds for the image
# --------------------------------------------------
lat_lims = [39.0, 42.5]
lon_lims = [-123., -119.]

# Set the map projection and center the data
# ------------------------------------------
my_area = scn['true_color'].attrs['area'].compute_optimal_bb_area({\
    'proj':'lcc', 'lon_0': lon_lims[0], 'lat_0': lat_lims[0], \
    'lat_1': lat_lims[0], 'lat_2': lat_lims[0]})
new_scn = scn.resample(my_area)

# Enhance the image for plotting
# ------------------------------
var = get_enhanced_image(new_scn['true_color']).data
var = var.transpose('y','x','bands')

# Extract the map projection from the data for plotting
# -----------------------------------------------------
crs = new_scn['true_color'].attrs['area'].to_cartopy_crs()

# Plot the true-color data
# ------------------------
plt.close('all')
ax = plt.axes(projection=crs)
ax.imshow(var.data, transform = crs, extent=(var.x[0], var.x[-1], var.y[-1], \
    var.y[0]), origin='upper')

# Zoom in the figure if desired
# -----------------------------
zoom = True 
if(zoom):
    ax.set_extent([lon_lims[0],lon_lims[1],lat_lims[0],lat_lims[1]],\
                   crs = ccrs.PlateCarree())
    zoom_add = '_zoom'
else:
    zoom_add = ''

# Add a title to the figure
# -------------------------
ax.set_title(sat_name + ' MODIS ' + modis_date.strftime('%Y-%m-%d %H:%M'))

# Save the figure
# ---------------
outname = 'modis_' + sat_name + '_true_color_' + \
    modis_date.strftime('%Y%m%d%H%M') + zoom_add + '_satpy.png'
#plt.savefig(outname,dpi=300)
print("Saved image",outname)
plt.show()
