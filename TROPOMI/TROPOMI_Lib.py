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
from scipy import stats
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


trop_to_omi_dict = {
    '201807051819': '201807051856',
    '201807052142': '201807052213',
    '201908110044': '201908110033', 
    #'202108010207': '202108010117',
    '202108010207': '202108010256',
    '202108010349': '202108010435',
    '202108010530': '202108010614',
    '202108011539': '202108011607',
    '202108011902': '202108011925',
    '202108012044': '202108012103',
    '202108012225': '202108012242',
    # TROP time      OMI time

}

json_time_database = home_dir + '/Research/Arctic_compares/json_comp_times.txt'
json_file_database = home_dir + '/Research/Arctic_compares/json_comp_files.txt'

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Automation
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# This function automates everything for a single date string:
# the data downloading/matching, the first-look comparison-making, the 
# hdf5 file-making, and the data copying to raindrop.
# Flags:
#   - include_tropomi: adds the colocated TROPOMI data to the colocation
#   - file package
def automate_TROPOMI_preprocess(date_str, download = True, images = True, \
        process = True, omi_dtype = 'ltc3', include_tropomi = True,\
        remove_bad_OMI = False):

    if(home_dir + '/Research/OMI/' not in sys.path):
        sys.path.append(home_dir + '/Research/OMI/')

    from OMILib import download_OMI_all_HDF, identify_OMI_HDF_swaths

    if(isinstance(date_str, str)):
        date_str = [date_str]

    ## Reload the json file here
    #with open(json_file_database, 'r') as fin:
    #    aerosol_event_dict = json.load(fin)

    # date_str should consist of an array of YYYYMMDD strings
    for dstr in date_str:

        dt_date_str = datetime.strptime(dstr, '%Y%m%d')

        print(dstr)

        # Determine if TROPOMI data are already downloaded for this date
        # --------------------------------------------------------------
        #trop_files = glob(dt_date_str.strftime(data_dir + \
        #    'TROPOMI*_%Ym%m%dt*-o*.nc'))
        #print(trop_files)

        #download_TROPOMI_match_OMI(date_str, \
        #    save_path = home_dir + '/Research/TROPOMI'):

        ## If no TROPOMI data
        #if(len(trop_files) == 0):
 
        # Look for OMI HDF5 files for this date
        # -------------------------------------
        ##!#omi_files = glob(dt_date_str.strftime(\
        ##!#    '/home/bsorenson/data/OMI/H5_files/OMI-A*_%Ym%m%d*.he5'))

        ##!## If no OMI data, 
        ##!#if(len(omi_files) == 0):
        ##!#
        ##!#    # Download the OMI files for this date
        ##!#    # ------------------------------------
        ##!#    download_OMI_all_HDF(dstr)

        ##!#    omi_files = glob(dt_date_str.strftime(\
        ##!#        '/home/bsorenson/data/OMI/H5_files/OMI-A*_%Ym%m%d*.he5'))
           
        # Screen the OMI swaths to determine which contain significant
        # aerosol (remove bad files if desired)
        # ------------------------------------------------------------
        good_omi_times = identify_OMI_HDF_swaths(dstr, \
            minlat = 70., min_AI = 2.0, remove_bad = remove_bad_OMI)

        # Download matching TROPOMI files for each of the remaining 
        # good swaths
        # ---------------------------------------------------------
        for omi_date in good_omi_times:
            trop_name = download_TROPOMI_file(omi_date)

            # Retrieve the OMI filename for this date
            # ---------------------------------------
            omi_dt_date = datetime.strptime(omi_date, '%Y%m%d%H%M')
            omi_file_name = glob(omi_dt_date.strftime(home_dir + \
                '/data/OMI/H5_files/OMI-*V_%Ym%m%dt%H%M*.he5'\
                ))[0].strip().split('/')[-1]

            trop_dtime = datetime.strptime(\
                trop_name.strip().split('/')[-1][33:49], '%Ym%m%dt%H%M%S')
            local_trop_time = trop_dtime.strftime('%Y%m%d%H%M')

            print(omi_file_name.strip().split('/')[-1], trop_name)
            #print(omi_file_name, trop_dtime.strftime('%Y%m%d%H%M'))

            data_dir = home_dir + '/Research/TROPOMI/prep_data/'

            # Make a data storage directory, if one does not already exist
            # ------------------------------------------------------------
            save_dir = data_dir + omi_dt_date.strftime('%Y%m%d%H%M') + '/'
            short_save_dir = omi_dt_date.strftime('%Y%m%d%H%M/')
            print(save_dir)

            if(not os.path.exists(save_dir)):
                print('Making ', save_dir)
                os.system('mkdir ' +  save_dir)
    
            # Run the TROPOMI converter to prep the files for movement
            # to raindrop
            # --------------------------------------------------------
            trop_file_name = convert_TROPOMI_to_HDF5(local_trop_time, \
                omi_name = omi_file_name, save_path = save_dir)
    
            # Now, write the name of the OMI file to a text file, which
            # will be zipped up with the TROPOMI file.
            # ---------------------------------------------------------
            omi_text_out = save_dir + 'omi_filename.txt'
            with open(omi_text_out, 'w') as fout:
                fout.write(omi_file_name)
    
            # Finally, gzip the data
            # ---------------------
            os.chdir(data_dir)
            cmnd = omi_dt_date.strftime('tar -cvzf ' + data_dir + \
                '/prepped_trop_data_%Y%m%d%H%M.tar.gz ' + short_save_dir)
            print(cmnd)
            os.system(cmnd)

        #trop_files = glob(dt_date_str.strftime(data_dir + \
        #    'TROPOMI*_%Ym%m%dt*-o*.nc'))
        #print(trop_files)
        
        # ENDIF
         
        #omi_files = glob(dt_date_str.strftime(\
        #    '/home/bsorenson/data/OMI/H5_files/OMI-A*_%Ym%m%d*.he5'))
    
        # For each TROPOMI file / date, 
        #for trop_file in trop_files:

            #trop_dtime = datetime.strptime(\
            #    trop_file.strip().split('/')[-1][33:49], '%Ym%m%dt%H%M%S')

            #local_trop_time = trop_dtime.strftime('%Y%m%d%H%M')
            #print(local_trop_time)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Downloading functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# This downloads the TROPOMI l1b HDF5 file that is closest to the passed
# date string from the GES DISC archive. 
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

    return local_file

def download_TROPOMI_match_OMI(date_str, \
        save_path = home_dir + '/Research/TROPOMI'):

    if(len(date_str) == 4):
        str_fmt = '%Y'
        out_fmt = '%Ym'
    elif(len(date_str) == 6):
        str_fmt = '%Y%m'
        out_fmt = '%Ym%m'
    elif(len(date_str) == 8):
        str_fmt = '%Y%m%d'
        out_fmt = '%Ym%m%d'
    elif(len(date_str) == 10):
        str_fmt = '%Y%m%d%H'
        out_fmt = '%Ym%m%dt%H'
    else:
        print("INVALID DATE STRING. RETURNING")
        sys.exit()
    
    dt_date_str = datetime.strptime(date_str, str_fmt)
    
    files = glob(dt_date_str.strftime(\
        '/home/bsorenson/data/OMI/H5_files/OMI-A*_' + out_fmt + '*.he5'))

    for ii, fname in enumerate(files):
        dt_name = datetime.strptime(fname.strip().split('/')[-1][20:34],\
            '%Ym%m%dt%H%M')

        local_date_str = dt_name.strftime('%Y%m%d%H%M')
        #print(dt_name.strftime('%Y%m%d%H%M'))
        found_file = download_TROPOMI_file(local_date_str, \
            dest_dir = '/home/bsorenson/data/TROPOMI/')
        trop_time = found_file.strip().split('/')[-1][33:49]
        print(trop_time)
        #convert_TROPOMI_to_HDF5(local_date_str)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Reading functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def convert_TROPOMI_to_HDF5(date_str, omi_name = None, \
        save_path = home_dir + '/Research/TROPOMI/'):

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
    dset.create_dataset('ssa0',      data = TROP_data['ssa0'][~TROP_data['uvai'].mask])
    dset.create_dataset('ssa1',      data = TROP_data['ssa1'][~TROP_data['uvai'].mask])
    dset.create_dataset('ssa2',      data = TROP_data['ssa2'][~TROP_data['uvai'].mask])

    if(omi_name is not None):
        dset.attrs['omi_filename'] = omi_name

    # Save, write, and close the HDF5 file
    # --------------------------------------
    dset.close()

    print("Saved file ",outfile)
    return outfile

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
    uvai   = data['SCIDATA/UVAerosolIndex'][:,:]
    ssa0   = data['SCIDATA/FinalAerosolSingleScattAlb'][:,:,0]
    ssa1   = data['SCIDATA/FinalAerosolSingleScattAlb'][:,:,1]
    ssa2   = data['SCIDATA/FinalAerosolSingleScattAlb'][:,:,2]
    lat    = data['GEODATA/latitude'][:,:]
    lon    = data['GEODATA/longitude'][:,:]
    time   = data['GEODATA/delta_time'][:]

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
    ssa0 = np.ma.masked_where(lat < minlat, ssa0)
    ssa1 = np.ma.masked_where(lat < minlat, ssa1)
    ssa2 = np.ma.masked_where(lat < minlat, ssa2)

    # Prep output dictionary
    # ----------------------
    TROP_data = {}
    TROP_data['date_str'] = date_str
    TROP_data['lat']  = lat
    TROP_data['lon']  = lon
    TROP_data['uvai'] = uvai
    TROP_data['ssa0'] = ssa0
    TROP_data['ssa1'] = ssa1
    TROP_data['ssa2'] = ssa2
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
    print(data.keys())
    for key in data.keys():
        TROP_coloc[key] = data[key][:,:]

    TROP_coloc['trop_ai'] = np.ma.masked_where(\
        TROP_coloc['trop_ai'] == -999.0, TROP_coloc['trop_ai'])
    TROP_coloc['trop_ssa0'] = np.ma.masked_where(\
        TROP_coloc['trop_ssa0'] == -999.0, TROP_coloc['trop_ssa0'])
    TROP_coloc['trop_ssa1'] = np.ma.masked_where(\
        TROP_coloc['trop_ssa1'] == -999.0, TROP_coloc['trop_ssa1'])
    TROP_coloc['trop_ssa2'] = np.ma.masked_where(\
        TROP_coloc['trop_ssa2'] == -999.0, TROP_coloc['trop_ssa2'])

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
    gl = ax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, \
        linewidth = 1, color = 'gray', alpha = 0.5, linestyle = '-',\
        y_inline = True, xlocs = range(-180, 180, 30), ylocs = range(70, 90, 10))
 
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

    plot_data = np.ma.masked_where(TROP_coloc['omi_lat'] < minlat, TROP_coloc[var])

    labelsize = 10
    mesh = ax.pcolormesh(TROP_coloc['omi_lon'], TROP_coloc['omi_lat'], plot_data, \
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

# slope: linear or thiel-sen
def plot_compare_OMI_TROPOMI(date_str, minlat = 65., slope = 'linear', \
        vmin = None, vmax = None, save = False):

    
    if(home_dir + '/Research/OMI/' not in sys.path):
        sys.path.append(home_dir + '/Research/OMI/')

    from OMILib import readOMI_swath_shawn, plotOMI_single_swath, \
        readOMI_swath_hdf

    # Read shawn (and raw) OMI data
    # -----------------------------
    omi_date = trop_to_omi_dict[date_str]
    #omi_date = date_str
    #OMI_base = readOMI_swath_shawn(omi_date, latmin = minlat, \
    #    shawn_path = home_dir + '/data/OMI/shawn_files/ltc3/')
    OMI_base = readOMI_swath_hdf(omi_date, 'control', latmin = minlat, \
        skiprows = [52])

    # Read high-res TROPOMI data
    # --------------------------
    TROP_full = read_TROPOMI(date_str)

    # Read colocated TROPOMI data (contains OMI data)
    # -----------------------------------------------
    TROP_coloc = read_TROPOMI_coloc(omi_date)
    #TROP_coloc = read_TROPOMI_coloc(date_str)

    if(TROP_coloc == -1):
        print("NO PLOTTING")
        return

    # Set up figure
    # -------------
    plt.close('all')
    fig = plt.figure(figsize = (11,7))
    ax1 = fig.add_subplot(2,3,1, projection = mapcrs) # Raw TROPOMI
    ax2 = fig.add_subplot(2,3,2, projection = mapcrs) # OMI raw
    ax3 = fig.add_subplot(2,3,3, projection = mapcrs) # OMI pert
    ax4 = fig.add_subplot(2,3,4, projection = mapcrs) # binned TROPOMI
    ax5 = fig.add_subplot(2,3,5)                      # TROP vs OMI raw
    ax6 = fig.add_subplot(2,3,6)                      # TROP vs OMI pert

    # Plot spatial OMI and TROPOMI data
    # ---------------------------------
    # Plot 1: raw TROPOMI
    plot_TROPOMI(TROP_full, ax = ax1, minlat = minlat, vmin = vmin, \
        vmax = vmax, circle_bound = True, zoom = True, save = False)
    ax1.set_title('TROPOMI UVAI')
    
    # Plot 2: raw OMI
    plotOMI_single_swath(ax2, OMI_base, pvar = 'UVAI', \
        title = 'OMI UVAI Raw', circle_bound = True, gridlines = False, \
        label = 'UVAI', vmin = vmin, vmax = vmax)

    ##!## Plot 3: perturbed OMI
    ##!#labelsize = 10
    ##!#mesh = ax3.pcolormesh(TROP_coloc['omi_lon'], TROP_coloc['omi_lat'], \
    ##!#    TROP_coloc['omi_uvai_pert'], transform = datacrs, shading = 'auto', \
    ##!#    cmap = 'jet', vmin = vmin, vmax = vmax)
    ##!#cbar = plt.colorbar(mesh, ax = ax3, orientation='vertical',\
    ##!#    pad=0.04, fraction = 0.040, extend = 'both')
    ##!#cbar.set_label('UV Aerosol Index', size = labelsize, weight = 'bold')
    ##!#ax3.coastlines()
    ##!#ax3.set_extent([-180, 180, minlat, 90], datacrs)
    ##!#ax3.set_boundary(circle, transform = ax3.transAxes)
    ##!#ax3.set_title('OMI UVAI Perturbation')

    # Plot 4: Plot binned TROPOMI
    plot_TROPOMI_coloc(TROP_coloc, 'trop_ai', ax = ax4, minlat = minlat, \
        vmin = vmin, vmax = vmax, circle_bound = True, zoom = True, \
        save = False)
    ax4.set_title('TROPOMI UVAI\nBinned to OMI Grid')
 
    # Plot scatter of OMI vs TROPOMI
    # ------------------------------
    # Plot 5: TROPOMI vs raw OMI   
    mask_OMI_raw  = np.ma.masked_where(OMI_base['UVAI'] < -2e5, OMI_base['UVAI'])
    mask_OMI_raw  = np.ma.masked_invalid(mask_OMI_raw)
    match_TROP    = TROP_coloc['trop_ai'][~mask_OMI_raw.mask]
    mask_OMI_raw  = mask_OMI_raw.compressed()

 
    # Get rid of missing TROPOMI pixels
    mask_OMI_raw   = np.ma.masked_where(match_TROP == -999.0, mask_OMI_raw)
    mask_TROP      = np.ma.masked_where(match_TROP == -999.0, match_TROP)
    final_OMI_raw  = mask_OMI_raw.compressed()
    final_TROP     = mask_TROP.compressed()

    print('HERE', np.nanmin(final_OMI_raw), np.nanmax(final_OMI_raw))
    print('HERE', np.nanmin(match_TROP), np.nanmax(match_TROP))
    print('HERE', np.nanmin(final_TROP), np.nanmax(final_TROP))

    xy = np.vstack([final_OMI_raw,final_TROP])
    z = stats.gaussian_kde(xy)(xy)
    ax5.scatter(final_OMI_raw, final_TROP, c = z, s = 6)
    zdata = plot_trend_line(ax5, final_OMI_raw, final_TROP, color='red', linestyle = '-', \
        linewidth = 1.5, slope = slope)
    ax5.set_title('TROPOMI = ' + str(np.round(zdata.slope,3)) + ' * OMI + ' + \
        str(np.round(zdata.intercept,3)) + '\n$r^{2}$ = ' + \
        str(np.round(zdata.rvalue**2., 3)))

    # Plot one to one line
    xlim = ax5.get_xlim()
    ylim = ax5.get_ylim()
    xvals = np.arange(xlim[0] - 1, xlim[-1] + 1)
    ax5.plot(xvals, xvals, '-', color = 'tab:cyan', linewidth = 1.5)
    ax5.set_xlim(xlim)   
    ax5.set_ylim(ylim)   
    ax5.grid()
    ax5.set_xlabel('OMI UVAI Raw')
    ax5.set_ylabel('TROPOMI UVAI')
    
    ##!## Plot 6: TROPOMI vs OMI pert
    ##!#mask_OMI = np.ma.masked_invalid(TROP_coloc['omi_uvai_pert'])
    ##!#match_TROP = TROP_coloc['trop_ai'][~mask_OMI.mask]
    ##!#mask_OMI   = mask_OMI.compressed()

    ##!## Get rid of missing TROPOMI pixels
    ##!#mask_OMI  = np.ma.masked_where(match_TROP == -999.0, mask_OMI)
    ##!#mask_TROP = np.ma.masked_where(match_TROP == -999.0, match_TROP)
    ##!#final_OMI  = mask_OMI.compressed()
    ##!#final_TROP = mask_TROP.compressed()

    ##!#xy = np.vstack([final_OMI,final_TROP])
    ##!#z = stats.gaussian_kde(xy)(xy)
    ##!#ax6.scatter(final_OMI, final_TROP, c = z, s = 6)
    ##!#zdata = plot_trend_line(ax6, final_OMI, final_TROP, color='red', linestyle = '-', \
    ##!#    linewidth = 1.5, slope = slope)
    ##!#ax6.set_title('TROPOMI = ' + str(np.round(zdata.slope,3)) + ' * OMI + ' + \
    ##!#    str(np.round(zdata.intercept,3)) + '\n$r^{2}$ = ' + \
    ##!#    str(np.round(zdata.rvalue**2., 3)))

    ##!## Plot one to one line
    ##!#xlim = ax6.get_xlim()
    ##!#ylim = ax6.get_ylim()
    ##!#xvals = np.arange(xlim[0] - 1, xlim[-1] + 1)
    ##!#ax6.plot(xvals, xvals, '-', color = 'tab:cyan', linewidth = 1.5)
    ##!#ax6.set_xlim(xlim)   
    ##!#ax6.set_ylim(ylim)   
    ##!#ax6.grid()
    ##!#ax6.set_xlabel('OMI UVAI Perturbed')
    ##!#ax6.set_ylabel('TROPOMI UVAI')

    plt.suptitle(date_str)
 
    # Polish the figure
    # ----------------- 
    fig.tight_layout()

    if(save):
        outname = 'tropomi_omi_combined_' + date_str + '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

def plot_TROPOMI_row_avg(date_str, plot_swath = False, minlat = 65., \
        save = True):
    
    if(len(date_str) == 4):
        str_fmt = '%Y'
        out_fmt = '%Ym'
    elif(len(date_str) == 6):
        str_fmt = '%Y%m'
        out_fmt = '%Ym%m'
    elif(len(date_str) == 8):
        str_fmt = '%Y%m%d'
        out_fmt = '%Ym%m%d'
    elif(len(date_str) == 10):
        str_fmt = '%Y%m%d%H'
        out_fmt = '%Ym%m%dt%H'
    else:
        print("INVALID DATE STRING. RETURNING")
        sys.exit()
    
    dt_date_str = datetime.strptime(date_str, str_fmt)
    
    files = glob(dt_date_str.strftime(\
        '/home/bsorenson/data/TROPOMI/TROPOMI-S*TROPOMAER_' + out_fmt + \
        '*.nc'))
    
    num_files = len(files)
    
    swath_ais = np.full((num_files, 4000, 450), np.nan)
    swath_lat = np.full((num_files, 4000, 450), np.nan)
    
    latmin = minlat
    for ii, fname in enumerate(files):
        print(fname.strip().split('/')[-1][33:47])
        dt_name = datetime.strptime(fname.strip().split('/')[-1][33:47],\
            '%Ym%m%dt%H%M')

        #if(download_match_TROP):    
        #    download_TROPOMI_file(dt_name.strftime('%Y%m%d%H%M'), \
        #        dest_dir = '/home/bsorenson/data/TROPOMI/')

        # Read and resample each file
        # ---------------------------
        local_date_str = dt_name.strftime('%Y%m%d%H%M')
        TROP_data = read_TROPOMI(local_date_str, minlat = minlat)
   
        if(plot_swath): 
            # Plot this swath
            # ---------------
            plot_TROPOMI_figure(local_date_str, minlat = 65., \
                vmin = -2, vmax = 4, circle_bound = True, \
                ptitle = '', zoom = True, save = True)
    
        swath_ais[ii, :TROP_data['uvai'].shape[0], :] = TROP_data['uvai'][:,:]
        swath_lat[ii, :TROP_data['lat'].shape[0], :]  = TROP_data['lat'][:,:]
    
    # Clean the data
    swath_ais = np.ma.masked_where(swath_ais < -50, swath_ais)
    swath_ais = np.ma.masked_invalid(swath_ais)
    swath_ais = np.ma.masked_where(swath_lat < latmin, swath_ais)
    
    # Average all the data for each day along each row, so there
    # is one set of row averages per day
    day_avgs = np.nanmean(swath_ais, axis = 1)
    # Average all the daily row averages
    row_avgs = np.nanmean(day_avgs, axis = 0)
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(row_avgs)
    ax.set_xlabel('Row #')
    ax.set_ylabel('Daily average AI')
    ax.set_title('TROPOMI daily average UVAI by row\nNorth of 65\n' + date_str)
    ax.grid()
    fig.tight_layout()
    
    outname = 'tropomi_row_avgs_' + date_str + '.png'
    fig.savefig(outname, dpi = 300)
    print("Saved image", outname)
    

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


