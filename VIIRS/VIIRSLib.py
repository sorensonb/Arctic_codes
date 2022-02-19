"""
  NAME:
    VIIRSLib

  PURPOSE:

"""

import numpy as np
import numpy.ma as ma
import sys
from netCDF4 import Dataset
from datetime import datetime,timedelta
from dateutil.relativedelta import relativedelta
from scipy.stats import pearsonr,spearmanr
import subprocess
from scipy import stats
from pyhdf import SD
import pandas as pd
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as cm
from matplotlib.dates import DateFormatter
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from satpy import find_files_and_readers
from satpy.scene import Scene
from satpy.writers import get_enhanced_image
from glob import glob

sys.path.append('/home/bsorenson/')
from python_lib import plot_trend_line, plot_subplot_label, plot_figure_text, \
    nearest_gridpoint, aerosol_event_dict, init_proj

mapcrs = ccrs.LambertConformal(central_longitude = -110, central_latitude = 40)

def plot_VIIRS_granule_DNB(filename, ax, band = 'M15'):
    
    # Read the filename to determine the file type
    # --------------------------------------------
    total_split = filename.split('/')
    name_split = total_split[-1].split('.')
    dataset_name = name_split[0]
    if(dataset_name == 'VNP46A1'):
        cmap = 'cividis'

        viirs_ds = h5py.File(filename, 'r')
        lat_max = viirs_ds['HDFEOS/GRIDS/VNP_Grid_DNB/'].attrs['NorthBoundingCoord']
        lat_min = viirs_ds['HDFEOS/GRIDS/VNP_Grid_DNB/'].attrs['SouthBoundingCoord']
        lon_max = viirs_ds['HDFEOS/GRIDS/VNP_Grid_DNB/'].attrs['EastBoundingCoord']
        lon_min = viirs_ds['HDFEOS/GRIDS/VNP_Grid_DNB/'].attrs['WestBoundingCoord']
        dnb = viirs_ds['HDFEOS/GRIDS/VNP_Grid_DNB/Data Fields/DNB_At_Sensor_Radiance_500m']

        # NOTE: Need to convert these values to actual radiance using the 
        # scale factor
        dnb_scale  = viirs_ds['HDFEOS/GRIDS/VNP_Grid_DNB/Data Fields/DNB_At_Sensor_Radiance_500m'].attrs['scale_factor']
        dnb_offset = viirs_ds['HDFEOS/GRIDS/VNP_Grid_DNB/Data Fields/DNB_At_Sensor_Radiance_500m'].attrs['add_offset']

        dnb = (dnb - dnb_offset) * dnb_scale

        lats = np.linspace(lat_max[0], lat_min[0], dnb.shape[1])
        lons = np.linspace(lon_min[0], lon_max[0], dnb.shape[0])
        x, y = np.meshgrid(lons, lats)

    elif(dataset_name[:5] == 'VNP02'):

        # Extract whether DNB or MOD
        # --------------------------
        data_type = dataset_name[5:]

        # Determine if geolocation data is present in same path
        # -----------------------------------------------------
        print("Searching for geolocation data for",'/'.join(total_split[:-1]) \
            + '/VNP03' + data_type + '.' + '.'.join(name_split[1:4]))
        geoloc_list = glob('/'.join(total_split[:-1]) + '/VNP03' + data_type + \
            '.' + '.'.join(name_split[1:4]) + '*') 
        if(len(geoloc_list) == 0):
            print("ERROR: no geolocation found for file",filename)
            return
        else:
            viirs_ds   = h5py.File(filename, 'r')
            viirs_ds_g = h5py.File(geoloc_list[0], 'r')

            # Extract the DNB values
            # ----------------------
            if(data_type == 'DNB'):
                cmap = 'Greys_r'
                dnb = viirs_ds['observation_data/DNB_observations'][:]

                dnb = np.ma.masked_where(dnb < 0, dnb)
                dnb = dnb * 1e9
                #dnb = np.log10(dnb)
                #vmax = np.log10(38.2)
                vmax = 10.0

            else:
                cmap = 'Greys_r'
                dnb = viirs_ds['observation_data/' + band][:]
                dnb_scale = viirs_ds['observation_data/' + band].attrs['scale_factor']
                dnb_offset = viirs_ds['observation_data/' + band].attrs['add_offset']

                dnb = np.ma.masked_where(dnb < 0, dnb)
                dnb = (dnb - dnb_offset) * dnb_scale
                print(np.min(dnb), np.max(dnb))
                vmax = None

            # Extract the lat and lons
            # ------------------------
            x = viirs_ds_g['geolocation_data']['longitude'][:]
            y = viirs_ds_g['geolocation_data']['latitude'][:]

            viirs_ds.close()
            viirs_ds_g.close() 

    ax.pcolormesh(x,y, dnb, cmap = cmap, vmin = None, vmax = vmax, \
        transform = ccrs.PlateCarree(), shading = 'auto')
    viirs_ds.close()

def plot_VIIRS_DNB(filename, band = 'M12', zoom = True):

    mapcrs = init_proj('202107232155')
    plt.close('all')
    fig = plt.figure(figsize = (6,6))
    ax = fig.add_subplot(1,1,1, projection = mapcrs)
    if(isinstance(filename, list)):

        for ff in filename:
            plot_VIIRS_granule_DNB(ff,ax, band = band)

    elif(isinstance(filename, str)):
        plot_VIIRS_granule_DNB(filename,ax, band = band)

    if(zoom):
        ax.set_extent([-122.0,-119.5,39.5,42.0],\
                       ccrs.PlateCarree())
        
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.STATES)
    plt.show()

