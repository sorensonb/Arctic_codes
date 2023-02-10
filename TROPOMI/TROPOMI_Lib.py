#!/usr/bin/env python
"""
  NAME:
    TROPOMI_Lib.py

  PURPOSE:
    House functions related to reading variables from TROPOMI data files

  MODIFICATIONS:
    Blake Sorenson <blake.sorenson@und.edu>         - 2022/12/20:
        Written (modification of plot_single_TROPOMI.py by Jacob Zanker)

"""

import sys
from glob import glob
import numpy as np
import subprocess
import netCDF4 as nc4
from netCDF4 import Dataset
from datetime import datetime, timedelta
import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

mapcrs = ccrs.NorthPolarStereo()
datacrs = ccrs.PlateCarree()

home_dir = os.environ['HOME']
sys.path.append(home_dir)
from python_lib import *

data_dir = home_dir + '/data/TROPOMI/'

def read_TROPOMI(date_str, minlat = 60.):

    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')

    # Find the matching file
    # ----------------------
    files = glob(dt_date_str.strftime(data_dir + 'TROPOMI-Sentinel*_%Ym%m%dt%H%M*.nc'))
    if(len(files) == 0):
        print("ERROR: Unable to find TROPOMI file for ",date_str)
        print("RETURNING")
        return -1

    infile = files[0]

    # Open the file
    # -------------
    data = Dataset(infile, 'r')

    # Mask the data outside of the minlat requirement
    # -----------------------------------------------
    uvai  = data['SCIDATA/UVAerosolIndex'][:,:]
    lat   = data['GEODATA/latitude'][:,:]
    lon   = data['GEODATA/longitude'][:,:]
    time  = data['GEODATA/delta_time'][:]

    # Pull out the base file time
    # ---------------------------
    base_date = datetime.strptime(data['GEODATA/delta_time'].units, \
        'milliseconds since %Y-%m-%d %H:%M:%S')

    # Calculate the individual times
    # ------------------------------
    abs_time = np.array([base_date + timedelta(milliseconds = int(tval)) \
        for tval in data['GEODATA/delta_time'][:]])

    data.close()

    uvai = np.ma.masked_where(lat < minlat, uvai)

    # Prep output dictionary
    # ----------------------
    TROP_data = {}
    TROP_data['date_str'] = date_str
    TROP_data['lat']  = lat
    TROP_data['lon']  = lon
    TROP_data['uvai'] = uvai
    TROP_data['times'] = abs_time
    
    return TROP_data

def plot_TROPOMI(TROP_data, ax = None, minlat = 60, vmin = -2, vmax = 3, \
        circle_bound = True, zoom = True, save = False):

    if(isinstance(TROP_data, str)):
        dt_date_str = datetime.strptime(TROP_data, '%Y%m%d%H%M')
        TROP_data = read_TROPOMI(TROP_data)
    else:
        dt_date_str = datetime.strptime(TROP_data['date_str'], '%Y%m%d%H%M')

    in_ax = True
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, projection = mapcrs)

    plot_data = np.ma.masked_where(TROP_data['lat'] < minlat, TROP_data['uvai'])

    labelsize = 10
    mesh = ax.pcolormesh(TROP_data['lon'], TROP_data['lat'], plot_data, \
        cmap = 'jet', vmin = vmin, vmax = vmax, \
        transform = ccrs.PlateCarree(), shading = 'auto')
    cbar = plt.colorbar(mesh, ax = ax, orientation='vertical',\
        pad=0.04, fraction = 0.040, extend = 'both')
    cbar.set_label('UV Aerosol Index', size = labelsize, weight = 'bold')
    #cbar.ax.tick_params(labelsize = labelticksize)
    
    if(circle_bound):
        ax.set_boundary(circle, transform = ax.transAxes)

    if(zoom):
        ax.set_extent([-180.0,180.0,minlat,90.0],\
                       ccrs.PlateCarree())

    if(not in_ax):
        plt.show()

def plot_TROPOMI_figure(date_str, minlat = 65., vmin = None, vmax = None, \
        circle_bound = True, ptitle = '', zoom = True, \
        save = False):

    plt.close('all')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection = mapcrs)

    # Read the data for this granule
    # ------------------------------
    TROP_data = read_TROPOMI(date_str, minlat = minlat)

    if(TROP_data == -1):
        print("ERROR with TROPOMI read. Quitting")
        return

    # Plot the data for this granule
    # ------------------------------
    plot_TROPOMI(TROP_data, ax = ax, minlat = minlat, vmin = vmin, vmax = vmax, \
        circle_bound = True, zoom = True, save = False)

    ax.set_title('TROPOMI UVAI\n' + TROP_data['date_str'])
 
    ax.coastlines()

    fig.tight_layout()

    if(save):
        outname = 'tropomi_ai_' + date_str + '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()
##!## This function takes in a TROPOMI filename and returns a dictionary containing
##!## the data from the file. Assume that each TROPOMI file contains a single 
##!## variable (AI, NO2, etc.)
##!## ----------------------------------------------------------------------------- 
##!#def read_TROPOMI(filename):
##!#
##!#    # Open the netCDF file
##!#    # --------------------
##!#    data = Dataset(filename, 'r')  
##!# 
##!#    # Extract the metadata
##!#    # --------------------
##!#    lat       = data['PRODUCT/latitude'][0,:,:]    # Latitude
##!#    lon       = data['PRODUCT/longitude'][0,:,:]   # Longitude
##!#    time      = data['PRODUCT/time'][0]            # Base time of file. Seconds since 2010-01-01 00:00:00
##!#    base_tstr = ' '.join(data['PRODUCT/time'].units.split()[2:])
##!#    dtime     = data['PRODUCT/delta_time'][0,:]    # Milliseconds since 
##!#
##!#    # Depending on the filename, extract the matching data variable
##!#    # -------------------------------------------------------------
##!#    uvai = data['PRODUCT/aerosol_index_354_388'][0,:,:]
##!#    #uvai = data['PRODUCTS/aerosol_index_340_380'][0,:,:]
##!#
##!#    # Calculate the correct times
##!#    # ---------------------------
##!#    base_date = datetime.strptime(base_tstr, '%Y-%m-%d %H:%M:%S')
##!#    file_time = base_date + timedelta(seconds = int(time))
##!#
##!#    times = np.array([file_time + timedelta(milliseconds = int(ttime)) for ttime in dtime.data])
##!#
##!#    # Close the netCDF file
##!#    # ---------------------
##!#    data.close()
##!#
##!#    # Insert the variables into the dictionary
##!#    # ----------------------------------------
##!#    TROPOMI_dict = {}
##!#
##!#    TROPOMI_dict['lat'] = lat
##!#    TROPOMI_dict['lon'] = lon
##!#    TROPOMI_dict['time'] = times
##!#    TROPOMI_dict['AI'] = uvai
##!#
##!#    return TROPOMI_dict


