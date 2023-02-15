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

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Downloading functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# This downloads the TROPOMI l1b HDF5 file that is closest to the passed
# date string from the LAADS DAAC archive. 
# --------------------------------------------------------------------
def download_TROPOMI_file(date_str, dest_dir = data_dir):

    base_url = 'https://measures.gesdisc.eosdis.nasa.gov/data/AER/TROPOMAER.1'

    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')

    # For each desired channel, figure out the closest file time
    # to the input date
    # ----------------------------------------------------------
    try:
        files = listFD(dt_date_str.strftime(base_url + '/%Y/%j'), ext = '.nc')
    except subprocess.CalledProcessError:
        print("ERROR: No TROPOMI files for the input DTG",date_str)
        return -2

    if(len(files) == 0):
        print("ERROR: No TROPOMI files returned from the request. Exiting")
        return -1
    
    # Remove the timestamps from the file strings
    # -------------------------------------------
    total_files = files[::2]
    files_only = [tfile.strip().split('/')[-1] for tfile in total_files]

    # Use the actual file times to get timestamps
    # -------------------------------------------
    file_dates = [datetime.strptime(tfile[33:49],'%Ym%m%dt%H%M%S') for tfile in files_only]

    # Figure out where the closest files to the desired time are
    # ----------------------------------------------------------
    time_diffs = np.array([abs((dt_date_str - ddate).total_seconds()) \
        for ddate in file_dates])

    # Extract the index of the matching TROPOMI file
    # --------------------------------------------
    file_idx = np.argmin(time_diffs)
    local_file = files_only[file_idx]
    found_file = total_files[file_idx]
 
    # Check if the file is already downloaded
    # ---------------------------------------
    if(os.path.exists(dest_dir + local_file)):
        print(found_file + ' already exists. Not downloading')
    else: 
        # Download the file
        # -----------------
        base_cmnd = "wget --load-cookies ~/urs_cookies --save-cookies "+\
            "~/.urs_cookies --keep-session-cookies --content-disposition "
        cmnd = dt_date_str.strftime(base_cmnd + found_file)
        print(cmnd)
        os.system(cmnd)

        # Move the file to the destination folder
        # ---------------------------------------
        cmnd = "mv " + local_file + " " + dest_dir
        print(cmnd) 
        os.system(cmnd)

    return found_file

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Reading functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def convert_TROPOMI_to_HDF5(date_str, save_path = home_dir + '/Research/TROPOMI/'):

    # Read the TROPOMI data
    TROP_data = read_TROPOMI(date_str, minlat = 65.)
    
    # Prep the data for output for the temporary pre-colocating file
    # --------------------------------------------------------------
    print(TROP_data['uvai'].shape)

    # Create a new netCDF dataset to write to the file
    # ------------------------------------------------
    outfile = save_path + 'tropomi_uvai_coloc_prep_' + date_str + '.hdf5'
    dset = h5py.File(outfile,'w')
 
    dset.create_dataset('latitude',  data = TROP_data['lat'][~TROP_data['uvai'].mask])
    dset.create_dataset('longitude', data = TROP_data['lon'][~TROP_data['uvai'].mask])
    dset.create_dataset('uvai',      data = TROP_data['uvai'].compressed())

    # Save, write, and close the HDF5 file
    # --------------------------------------
    dset.close()

    print("Saved file ",outfile)

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

def read_TROPOMI_coloc(date_str, minlat = 60.):

    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')

    # Find the matching file
    # ----------------------
    files = glob(dt_date_str.strftime(home_dir + '/Research/TROPOMI/colocated_tropomi_%Y%m%d%H%M.hdf5'))
    if(len(files) == 0):
        print("ERROR: Unable to find TROPOMI colocation file for ",date_str)
        print("RETURNING")
        return -1

    infile = files[0]

    # Open the file
    # -------------
    data = h5py.File(infile, 'r')


    # Prep output dictionary
    # ----------------------
    TROP_coloc = {}
    TROP_coloc['date_str'] = date_str
    for key in data.keys():
        TROP_coloc[key] = data[key][:,:]

    data.close()
    
    return TROP_coloc

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Plotting functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
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
    ax.coastlines()
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

def plot_TROPOMI_coloc(TROP_coloc, var, ax = None, minlat = 60, vmin = -2, vmax = 3, \
        circle_bound = True, zoom = True, save = False):

    if(isinstance(TROP_coloc, str)):
        dt_date_str = datetime.strptime(TROP_coloc, '%Y%m%d%H%M')
        TROP_coloc = read_TROPOMI_coloc(TROP_coloc)
    else:
        dt_date_str = datetime.strptime(TROP_coloc['date_str'], '%Y%m%d%H%M')

    in_ax = True
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, projection = mapcrs)

    plot_data = np.ma.masked_where(TROP_coloc['trop_lat'] < minlat, TROP_coloc[var])

    labelsize = 10
    mesh = ax.pcolormesh(TROP_coloc['trop_lon'], TROP_coloc['trop_lat'], plot_data, \
        cmap = 'jet', vmin = vmin, vmax = vmax, \
        transform = ccrs.PlateCarree(), shading = 'auto')
    cbar = plt.colorbar(mesh, ax = ax, orientation='vertical',\
        pad=0.04, fraction = 0.040, extend = 'both')
    cbar.set_label(var, size = labelsize, weight = 'bold')
    ax.coastlines()
    #cbar.ax.tick_params(labelsize = labelticksize)
    
    if(circle_bound):
        ax.set_boundary(circle, transform = ax.transAxes)

    if(zoom):
        ax.set_extent([-180.0,180.0,minlat,90.0],\
                       ccrs.PlateCarree())

    if(not in_ax):
        plt.show()

def plot_TROPOMI_coloc_figure(date_str, var, minlat = 65., vmin = None, vmax = None, \
        circle_bound = True, ptitle = '', zoom = True, \
        save = False):

    plt.close('all')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection = mapcrs)

    # Read the data for this granule
    # ------------------------------
    TROP_coloc = read_TROPOMI_coloc(date_str, minlat = minlat)

    if(TROP_coloc == -1):
        print("ERROR with TROPOMI read. Quitting")
        return

    # Plot the data for this granule
    # ------------------------------
    plot_TROPOMI_coloc(TROP_coloc, var, ax = ax, minlat = minlat, vmin = vmin, vmax = vmax, \
        circle_bound = True, zoom = True, save = False)

    ax.set_title('TROPOMI ' + var + '\n' + TROP_coloc['date_str'])
 
    ax.coastlines()

    fig.tight_layout()

    if(save):
        outname = 'tropomi_coloc_' + var + '_' + date_str + '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

def plot_compare_OMI_TROPOMI(date_str, minlat = 65., save = False):

    # Read high-res TROPOMI data
    # --------------------------
    TROP_full = read_TROPOMI(date_str)

    # Read colocated TROPOMI data (contains OMI data)
    # -----------------------------------------------
    TROP_coloc = read_TROPOMI_coloc(date_str)

    if(TROP_coloc == -1):
        print("NO PLOTTING")
        return

    # Set up figure
    # -------------
    plt.close('all')
    fig = plt.figure(figsize = (9,9))
    ax1 = fig.add_subplot(2,2,1, projection = mapcrs)
    ax2 = fig.add_subplot(2,2,2, projection = mapcrs)
    ax3 = fig.add_subplot(2,2,3, projection = mapcrs)
    ax4 = fig.add_subplot(2,2,4)

    # Plot spatial OMI and TROPOMI data
    # ---------------------------------
    labelsize = 10
    mesh = ax1.pcolormesh(TROP_coloc['trop_lon'], TROP_coloc['trop_lat'], \
        TROP_coloc['omi_uvai_pert'], transform = datacrs, shading = 'auto', \
        cmap = 'jet', vmin = -2, vmax = 3)
    cbar = plt.colorbar(mesh, ax = ax1, orientation='vertical',\
        pad=0.04, fraction = 0.040, extend = 'both')
    cbar.set_label('UV Aerosol Index', size = labelsize, weight = 'bold')
    ax1.coastlines()
    ax1.set_extent([-180, 180, minlat, 90], datacrs)
    ax1.set_boundary(circle, transform = ax1.transAxes)

    plot_TROPOMI(TROP_full, ax = ax2, minlat = minlat, vmin = -2, vmax = 3, \
        circle_bound = True, zoom = True, save = False)

    plot_TROPOMI_coloc(TROP_coloc, 'trop_ai', ax = ax3, minlat = minlat, \
        vmin = -2, vmax = 3, circle_bound = True, zoom = True, save = False)
 
    # Plot scatter of OMI vs TROPOMI
    # ------------------------------
    
    fig.tight_layout()
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


