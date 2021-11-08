#!/usr/bin/env python

"""
  NAME:

  PURPOSE:
  
    NOTE: The MODIS channel and true color functions are designed to work with
    HDF MODIS files retriefed from 
    the data ordering website at this address:
    https://ladsweb.modaps.eosdis.nasa.gov/search/order/1/MODIS:Aqua


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
#dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
#filename = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis']
filename = sys.argv[1]

print(filename)
##!#if(composite):
##!#    day_filenames = glob(filename[:50]+'*')
##!#    cmpst_add = '_composite'
##!#else:
day_filenames = glob(filename)
cmpst_add = ''

# Extract the modis true-color plot limits
# ----------------------------------------
##lat_lims = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis_Lat']
##lon_lims = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis_Lon']

# Use satpy (Scene) to open the file
# ----------------------------------
scn = Scene(reader = 'modis_l1b', filenames = day_filenames)

# Load true-color data 
scn.load(['true_color'])

scn.save_dataset('true_color','test_image_true2.png')
scn.show('true_color')

plt.show()
##!#if(save):
##!#    outname = 'modis_true_color_' + date_str + zoom_add + cmpst_add + '_satpy.png'
##!#    plt.savefig(outname,dpi=300)
##!#    print("Saved image",outname)
##!#else:
##!#    plt.show()
