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
#from python_lib import plot_trend_line, plot_subplot_label, plot_figure_text, \
#    nearest_gridpoint, aerosol_event_dict, init_proj, \
#    convert_radiance_to_temp
from python_lib import *

datacrs = ccrs.PlateCarree()
mapcrs = ccrs.LambertConformal(central_longitude = -110, central_latitude = 40)

channel_dict = {
    'I1': {
        'Bandwidth': [0.6, 0.68]
    },\
    'I2': {
        'Bandwidth': [0.85, 0.88]
    },\
    'I3': {
        'Bandwidth': [1.58, 1.64]
    },\
    'I4': {
        'Bandwidth': [3.55, 3.93]
    },\
    'I5': {
        'Bandwidth': [10.5, 12.4]
    },\
    'M1': {
        'Bandwidth': [0.402, 0.422]
    },\
    'M2': {
        'Bandwidth': [0.436, 0.454]
    },\
    'M3': {
        'Bandwidth': [0.478, 0.488]
    },\
    'M4': {
        'Bandwidth': [0.545, 0.565]
    },\
    'M5': {
        'Bandwidth': [0.662, 0.682]
    },\
    'M6': {
        'Bandwidth': [0.739, 0.754]
    },\
    'M7': {
        'Bandwidth': [0.846, 0.885]
    },\
    'M8': {
        'Bandwidth': [1.23, 1.25]
    },\
    'M9': {
        'Bandwidth': [1.371, 1.386]
    },\
    'M10': {
        'Bandwidth': [1.58, 1.64]
    },\
    'M11': {
        'Bandwidth': [2.23, 2.28]
    },\
    'M12': {
        'Bandwidth': [3.61, 3.79]
    },\
    'M13': {
        'Bandwidth': [3.97, 4.13]
    },\
    'M14': {
        'Bandwidth': [8.4, 8.7]
    },\
    'M15': {
        'Bandwidth': [10.26, 11.26]
    },\
    'M16': {
        'Bandwidth': [11.54, 12.49]
    },\
    'DNB': {
        #'Name': 'EV_1KM_Emissive',\
        #'Index': 0,\
        'Bandwidth': [0.5, 0.9]
    }
}

def find_plume_VIIRS(filename):
    VIIRS_data5  = readVIIRS_granule(filename, band = 'M05', zoom = True)
    VIIRS_data7  = readVIIRS_granule(filename, band = 'M07', zoom = True)
    VIIRS_data11 = readVIIRS_granule(filename, band = 'M11', zoom = True)
    VIIRS_data16 = readVIIRS_granule(filename, band = 'M16', zoom = True)

    tmp_data = np.full(VIIRS_data5['data'].shape, np.nan)
    in_cloud = (((VIIRS_data5['data'] + VIIRS_data7['data'] > 1.2) |\
                 (VIIRS_data16['data'] < 265.)) |\
                ((VIIRS_data5['data'] + VIIRS_data7['data'] > 0.7) &\
                 (VIIRS_data16['data'] < 285.)))
    tmp_data[in_cloud] = 2.
    tmp_data[~in_cloud] = 0.
    tmp_data[~in_cloud & (VIIRS_data5['data'] - \
        VIIRS_data11['data'] > 0.05) & (VIIRS_data11['data'] > 0.05)] = 1
    mask_tmp_data = np.ma.masked_invalid(tmp_data)

    # For the hash data (smoke), remove any pixels that do not have
    # tmp data equal to 1
    # -------------------------------------------------------------
    hash_data   = np.ma.masked_where(mask_tmp_data != 1, mask_tmp_data)
    nohash_data = np.ma.masked_where(mask_tmp_data == 1, mask_tmp_data)

    # Free the memory
    # ---------------
    VIIRS_data5  = None
    VIIRS_data7  = None
    VIIRS_data11 = None
    VIIRS_data16 = None

    ##!#MODIS_ch1  = read_MODIS_channel(dt_date_str, 1,  zoom = True)
    ##!#MODIS_ch5  = read_MODIS_channel(dt_date_str, 5,  zoom = True)
    ##!#MODIS_ch31 = read_MODIS_channel(dt_date_str, 31, zoom = True)

    ##!#screen_limit = 0.05
    ##!#max_ch = 350.
    ##!#test_data = MODIS_ch1['data'] - MODIS_ch5['data']
    ##!#hash_data   = np.ma.masked_where(test_data <  \
    ##!#    (np.nanmean(test_data) * screen_limit), MODIS_ch1['data'])
    ##!#hash_data = np.ma.masked_where(MODIS_ch5['data'] < 0.05, hash_data)
    ##!#nohash_data = np.ma.masked_where(test_data >= \
    ##!#    (np.nanmean(test_data) * screen_limit), MODIS_ch1['data'])

    ##!#hash_data = np.ma.masked_where(MODIS_ch31['data'] > max_ch, hash_data)
    ##!#nohash_data = np.ma.masked_where(MODIS_ch31['data'] > max_ch, nohash_data)

    return hash_data, nohash_data

def readVIIRS_granule(filename, band = 'M15', zoom = True):
   
    VIIRS_data = {}
    VIIRS_data['filename'] = filename
    VIIRS_data['band']     = band
 
    # Read the filename to determine the file type
    # --------------------------------------------
    total_split = filename.split('/')
    name_split = total_split[-1].split('.')
    dataset_name = name_split[0]
    if(dataset_name == 'VNP46A1'):
        cmap = 'cividis'
        dtype = 'DNB'        

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

        viirs_ds.close()

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
                dtype = 'DNB'
                label = 'Radiance [nW cm$^{-2}$ Sr$^{-1}$]'
                cmap = 'Greys_r'
                dnb = viirs_ds['observation_data/DNB_observations'][:]

                dnb = np.ma.masked_where(dnb < 0, dnb)
                dnb = dnb * 1e9
                #dnb = np.log10(dnb)
                #vmax = np.log10(38.2)
                vmax = 10.0

            else:
                dtype = band
            
                dnb = viirs_ds['observation_data/' + band][:]
                dnb_scale = viirs_ds['observation_data/' + band].attrs['scale_factor']
                dnb_offset = viirs_ds['observation_data/' + band].attrs['add_offset']

                dnb = np.ma.masked_where(dnb < 0, dnb)
                dnb = (dnb - dnb_offset) * dnb_scale
                print(np.min(dnb), np.max(dnb))

                # Use the band number to determine the colormap
                # ---------------------------------------------
                if(int(band[1:]) > 11):
                    cmap = 'plasma'

                    dnb = convert_radiance_to_temp(\
                        (np.average(channel_dict[band]['Bandwidth'])), dnb)

                    print('wavelength avg = ', np.average(channel_dict[band]['Bandwidth']))
                    print('Radiances  ', np.min(dnb), np.max(dnb))
                    label = 'Brightness temperature [K]'
                    #label = 'Radiance [W m$^{-2}$ μm$^{-1}$ Sr$^{-1}$]'
    
                    if(np.max(dnb) > 330):
                        vmax = 330
                    else:
                        vmax = None
                else:
                    cmap = 'Greys_r'
                    label = 'Reflectance'
                    if(np.max(dnb) > 1.0):
                        vmax = 1.0
                    else:
                        vmax = None
                    


            # Extract the lat and lons
            # ------------------------
            x = viirs_ds_g['geolocation_data']['longitude'][:]
            y = viirs_ds_g['geolocation_data']['latitude'][:]

            viirs_ds.close()
            viirs_ds_g.close() 

    if(zoom):
        # Mask MODIS_data['data'] that are outside the desired range
        # --------------------------------------------
        dnb[(((y < 39.5 - 0.1) | \
              (y > 42.0 + 0.1)) | \
             ((x < -122.0 - 0.1) | \
              (x > -119.5 + 0.1)))] = np.nan
        ##!#MODIS_data['lat'][ (((MODIS_data['lat'] < aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][0]) | \
        ##!#                     (MODIS_data['lat'] > aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][1])) | \
        ##!#                    ((MODIS_data['lon'] < aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][0]) | \
        ##!#                     (MODIS_data['lon'] > aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][1])))] = -999.
        ##!#MODIS_data['lon'][ (((MODIS_data['lat'] < aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][0]) | \
        ##!#                     (MODIS_data['lat'] > aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][1])) | \
        ##!#                    ((MODIS_data['lon'] < aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][0]) | \
        ##!#                     (MODIS_data['lon'] > aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][1])))] = -999.

        dnb = np.ma.masked_invalid(dnb)
        #dnb = np.ma.masked_where(dnb == -999., dnb)
        ##!#MODIS_data['lat'] = np.ma.masked_where(MODIS_data['lat'] == -999., MODIS_data['lat'])
        ##!#MODIS_data['lon'] = np.ma.masked_where(MODIS_data['lon'] == -999., MODIS_data['lon'])
        

    VIIRS_data['filename'] = filename
    VIIRS_data['dtype']    = dtype
    VIIRS_data['label']    = label
    VIIRS_data['lon']      = x
    VIIRS_data['lat']      = y
    VIIRS_data['data']     = dnb
    VIIRS_data['vmax']     = vmax
    VIIRS_data['cmap']     = cmap

    return VIIRS_data

def plot_VIIRS_granule(VIIRS_data, ax = None, labelsize = 12, \
        labelticksize = 10, zoom = True, show_smoke = False, vmax = None, \
        save = False):

    in_ax = True
    if(ax is None):
        in_ax = False
        plt.close('all')
        mapcrs = init_proj('202107232155')
        fig = plt.figure(figsize=(6,6))
        ax = plt.axes(projection = mapcrs)

    if(vmax is None):
        vmax = VIIRS_data['vmax']

    pdata = np.copy(VIIRS_data['data'])


    mesh = ax.pcolormesh(VIIRS_data['lon'], VIIRS_data['lat'], pdata, \
        cmap = VIIRS_data['cmap'], vmin = None, vmax = vmax, \
        transform = ccrs.PlateCarree(), shading = 'auto')
    cbar = plt.colorbar(mesh, ax = ax, orientation='vertical',\
        pad=0.03, extend = 'both')

    cbar.set_label(VIIRS_data['label'], size = labelsize, weight = 'bold')
    cbar.ax.tick_params(labelsize = labelticksize)
    ax.set_title(VIIRS_data['dtype'])

    if(show_smoke):    
        hash_data, nohash_data = find_plume_VIIRS(VIIRS_data['filename'])
        #plt.rcParams.update({'hatch.color': 'r'})
        hatch_shape = '\\\\'
        #hash0 = ax2.pcolor(MODIS_data_ch1['lon'],MODIS_data_ch1['lat'],\
        #    hash_data[:MODIS_data_ch1['lat'].shape[0], :MODIS_data_ch1['lat'].shape[1]], hatch = hatch_shape, alpha=0., transform = datacrs)
        ax.pcolor(VIIRS_data['lon'],VIIRS_data['lat'],\
            hash_data, hatch = hatch_shape, alpha=0.40, transform = datacrs)
     

    if(zoom):
        ax.set_extent([-122.0,-119.5,39.5,42.0],\
                       ccrs.PlateCarree())
    if(not in_ax):
        plt.show()

def plot_VIIRS_figure(filename, band = 'M12', ax = None, show_smoke = False, \
        vmax = None, ptitle = '', zoom = True):

    in_ax = True
    if(ax is None):
        in_ax = False
        mapcrs = init_proj('202107232155')
        plt.close('all')
        fig = plt.figure(figsize = (6,6))
        ax = fig.add_subplot(1,1,1, projection = mapcrs)

    if(isinstance(filename, list)):

        for ff in filename:
            #plot_VIIRS_granule_DNB(ff,ax, band = band)

            # Read the data for this granule
            # ------------------------------
            VIIRS_data = readVIIRS_granule(ff, band = band, zoom = zoom)
        
            # Plot the data for this granule
            # ------------------------------
            plot_VIIRS_granule(VIIRS_data, ax = ax, zoom = zoom, \
                vmax = vmax, show_smoke = show_smoke)

    elif(isinstance(filename, str)):
        #plot_VIIRS_granule_DNB(filename,ax, band = band)

        # Read the data for this granule
        # ------------------------------
        VIIRS_data = readVIIRS_granule(filename, band = band, zoom = zoom)

        # Plot the data for this granule
        # ------------------------------
        plot_VIIRS_granule(VIIRS_data, ax = ax, zoom = zoom, \
                vmax = vmax, show_smoke = show_smoke)

    if(zoom):
        ax.set_extent([-122.0,-119.5,39.5,42.0],\
                       ccrs.PlateCarree())
       
    ax.set_title(ptitle)
 
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.STATES)

    if(not in_ax):
        plt.show()

# Nighttime
#   22 July (203): 1000
#   23 July (204): 0942
# Daytime
#   22 July (203): 2124
#   23 July (204): None yet
#
##!#viirs_event_dict = {
##!#    '20210722': {
##!#        'daytime':   'A2021203.1000.',
##!#        'nighttime': 'A2021203.2124.'
##!#    }
##!#    '20210723': {
##!#        #'daytime':   'A2021204.0942.',
##!#        'nighttime': 'A2021204.0942.'
##!#    }
##!#    
##!#}

# NOTE: MOD data first, then DNB
def plot_VIIRS_scatter(fname1, fname2, band1, band2, ax = None, save = False):
    in_ax = True
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure(figsize = (6,6))
        ax = fig.add_subplot(1,1,1)

    # Read the data for this granule
    # ------------------------------
    VIIRS_data1 = readVIIRS_granule(fname1, band = band1, zoom = True)
    VIIRS_data2 = readVIIRS_granule(fname2, band = band2, zoom = True)

    # Match the DNB data to the MOD data
    # ----------------------------------
    compress_dat1_lat  = VIIRS_data1['lat'][~VIIRS_data1['data'].mask]  
    compress_dat1_lon  = VIIRS_data1['lon'][~VIIRS_data1['data'].mask]  
    compress_dat1_data = VIIRS_data1['data'][~VIIRS_data1['data'].mask]  
    compress_dat2_lat  = VIIRS_data2['lat'][~VIIRS_data2['data'].mask]  
    compress_dat2_lon  = VIIRS_data2['lon'][~VIIRS_data2['data'].mask]  
    compress_dat2_data = VIIRS_data2['data'][~VIIRS_data2['data'].mask]  

    if(len(compress_dat1_data) == len(compress_dat2_data)):

        #ax.scatter(compress_dat1_data, compress_dat2_data, s = 6)
        ax.set_xlabel(VIIRS_data1['dtype'])
        ax.set_ylabel(VIIRS_data2['dtype']) 

        # Use find_plume to identify where the smoke is
        # ---------------------------------------------
        hash_data, nohash_data = find_plume_VIIRS(fname1)
        #plt.rcParams.update({'hatch.color': 'r'})
        hatch_shape = '\\\\'
        #hash0 = ax2.pcolor(MODIS_data_ch1['lon'],MODIS_data_ch1['lat'],\
        #    hash_data[:MODIS_data_ch1['lat'].shape[0], :MODIS_data_ch1['lat'].shape[1]], hatch = hatch_shape, alpha=0., transform = datacrs)
        #ax.pcolor(VIIRS_data1['lon'],VIIRS_data1['lat'],\
        #    hash_data, hatch = hatch_shape, alpha=0, transform = datacrs)

        ax.scatter(VIIRS_data1['data'][np.where(hash_data.mask)], \
            VIIRS_data2['data'][np.where(hash_data.mask)], s=6, \
            color='tab:orange', label='Outside Plume')
        ax.scatter(VIIRS_data1['data'][np.where(~hash_data.mask)], \
            VIIRS_data2['data'][np.where(~hash_data.mask)], s=6, \
            color='tab:blue', label='Inside Plume')
       
    else: 
        if(len(compress_dat1_data) > len(compress_dat2_data)):
            print(VIIRS_data1['dtype'], ' larger than ', VIIRS_data2['dtype'])
            big_label    = VIIRS_data1['dtype']
            small_label  = VIIRS_data2['dtype']
            big_data   = compress_dat1_data 
            big_lat    = compress_dat1_lat 
            big_lon    = compress_dat1_lon 
            small_data = compress_dat2_data 
            small_lat  = compress_dat2_lat 
            small_lon  = compress_dat2_lon 
        else:
            print(VIIRS_data2['dtype'], ' larger than ', VIIRS_data1['dtype'])
            big_label    = VIIRS_data2['dtype']
            small_label  = VIIRS_data1['dtype']
            big_data   = compress_dat2_data 
            big_lat    = compress_dat2_lat 
            big_lon    = compress_dat2_lon 
            small_data = compress_dat1_data 
            small_lat  = compress_dat1_lat 
            small_lon  = compress_dat1_lon 

        print("Before matching,")
        print("   big_data   = ", big_data.shape)
        print("   big_lat    = ", big_lat.shape)
        print("   big_lon    = ", big_lon.shape)
        print("   small_data = ", small_data.shape)
        print("   small_lat  = ", small_lat.shape)
        print("   small_lon  = ", small_lon.shape)
        big_match_to_small = np.array([big_data[nearest_gridpoint(tlat, tlon,\
             big_lat, big_lon)] for tlat, tlon in zip(small_lat, small_lon)])

        print("Before plotting,")
        print("   big_match_to_small = ", big_match_to_small.shape)
        print("   small_data         = ", small_data.shape)
        ax.scatter(small_data, big_match_to_small, s = 6)
        ax.set_ylim(0, 35)
        #ax.set_yscale('log')
        ax.set_xlabel(small_label)
        ax.set_ylabel(big_label) 
        #ax.scatter(VIIRS_data1['data'], VIIRS_data2['data'], s = 6)
        #ax.set_xlabel(VIIRS_data1['dtype'])
        #ax.set_ylabel(VIIRS_data2['dtype']) 
 
    if(not in_ax):
        plt.show() 


def plot_VIIRS_sixpanel(save = False):

    # Set up the figure
    # -----------------
    date_str = '20210723'
    mapcrs = init_proj('202107232155')
    fig = plt.figure(figsize = (11,8))
    ax1 = fig.add_subplot(2,3,1, projection = mapcrs) # Nighttime 0722 DNB
    ax2 = fig.add_subplot(2,3,2, projection = mapcrs) # Daytime 0722 VIS
    ax3 = fig.add_subplot(2,3,3, projection = mapcrs) # Nighttime 0723 DNB
    ax4 = fig.add_subplot(2,3,4, projection = mapcrs) # Nighttime 0722 IR
    ax5 = fig.add_subplot(2,3,5, projection = mapcrs) # Daytime 0722 IR
    ax6 = fig.add_subplot(2,3,6, projection = mapcrs) # Daytime 0723 IR

    # Read the VIIRS data for each time
    # ---------------------------------
    filename_n0722_dnb = glob('/home/bsorenson/data/VIIRS/DNB/VNP02DNB.A2021203.1000*.nc')
    filename_n0722_dat = glob('/home/bsorenson/data/VIIRS/DNB/VNP02MOD.A2021203.1000*.nc')
    filename_d0722_dat = glob('/home/bsorenson/data/VIIRS/DNB/VNP02MOD.A2021203.2124*.nc')
    filename_n0723_dnb = glob('/home/bsorenson/data/VIIRS/DNB/VNP02DNB.A2021204.0942*.nc')
    filename_n0723_dat = glob('/home/bsorenson/data/VIIRS/DNB/VNP02MOD.A2021204.0942*.nc')

    namesplit = filename_n0722_dat[0].split('/')[-1].split('.')
    date1 = datetime.strptime('.'.join(namesplit[1:3])[1:], '%Y%j.%H%M')
    namesplit = filename_d0722_dat[0].split('/')[-1].split('.')
    date2 = datetime.strptime('.'.join(namesplit[1:3])[1:], '%Y%j.%H%M')
    namesplit = filename_n0723_dnb[0].split('/')[-1].split('.')
    date3 = datetime.strptime('.'.join(namesplit[1:3])[1:], '%Y%j.%H%M')

    print(date1.strftime('%Y-%m-%d %H:%M'))
    print(date2.strftime('%Y-%m-%d %H:%M'))
    print(date3.strftime('%Y-%m-%d %H:%M'))

    plot_VIIRS_figure(filename_n0722_dnb, band = 'DNB', ax = ax1, vmax = 6,\
        ptitle = date1.strftime('%Y-%m-%d %H:%M'), zoom = True)
    plot_VIIRS_figure(filename_d0722_dat, band = 'M05', ax = ax2, \
        ptitle = date2.strftime('%Y-%m-%d %H:%M'), zoom = True)
    plot_VIIRS_figure(filename_n0723_dnb, band = 'DNB', ax = ax3, \
        ptitle = date3.strftime('%Y-%m-%d %H:%M'), zoom = True)

    plot_VIIRS_figure(filename_n0722_dat, band = 'M15', ax = ax4, zoom = True)
    plot_VIIRS_figure(filename_d0722_dat, band = 'M15', ax = ax5, zoom = True)
    plot_VIIRS_figure(filename_n0723_dat, band = 'M15', ax = ax6, zoom = True)

    plot_subplot_label(ax1, '(a)', backgroundcolor = 'white')
    plot_subplot_label(ax2, '(b)', backgroundcolor = 'white')
    plot_subplot_label(ax3, '(c)', backgroundcolor = 'white')
    plot_subplot_label(ax4, '(d)', backgroundcolor = 'white')
    plot_subplot_label(ax5, '(e)', backgroundcolor = 'white')
    plot_subplot_label(ax6, '(f)', backgroundcolor = 'white')

    plot_figure_text(ax1, 'VIIRS DNB', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax2, 'VIIRS 0.64 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax3, 'VIIRS DNB', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax4, 'VIIRS 10.76 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax5, 'VIIRS 10.76 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax6, 'VIIRS 10.76 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')

    fig.tight_layout()

    if(save):
        outname = 'viirs_sixpanel_' + date_str + '.png'
        plt.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else: 
        plt.show()
    
def plot_VIIRS_ninepanel(save = False):

    # Set up the figure
    # -----------------
    mapcrs = init_proj('202107232155')
    fig = plt.figure(figsize = (11,11))
    ax1 = fig.add_subplot(3,3,1, projection = mapcrs) # Nighttime 0722 DNB
    ax2 = fig.add_subplot(3,3,2, projection = mapcrs) # Daytime 0722 VIS
    ax3 = fig.add_subplot(3,3,3, projection = mapcrs) # Nighttime 0723 DNB
    ax4 = fig.add_subplot(3,3,4, projection = mapcrs) # Nighttime 0722 IR
    ax5 = fig.add_subplot(3,3,5, projection = mapcrs) # Daytime 0722 IR
    ax6 = fig.add_subplot(3,3,6, projection = mapcrs) # Daytime 0723 IR
    ax7 = fig.add_subplot(3,3,7)                      # Nighttime 0722 scatter
    ax8 = fig.add_subplot(3,3,8)                      # Daytime 0723 scatter
    ax9 = fig.add_subplot(3,3,9)                      # Daytime 0723 scatter

    # Read the VIIRS data for each time
    # ---------------------------------
    filename_n0722_dnb = glob('/home/bsorenson/data/VIIRS/DNB/VNP02DNB.A2021203.1000*.nc')
    filename_n0722_dat = glob('/home/bsorenson/data/VIIRS/DNB/VNP02MOD.A2021203.1000*.nc')
    filename_d0722_dat = glob('/home/bsorenson/data/VIIRS/DNB/VNP02MOD.A2021203.2124*.nc')
    filename_n0723_dnb = glob('/home/bsorenson/data/VIIRS/DNB/VNP02DNB.A2021204.0942*.nc')
    filename_n0723_dat = glob('/home/bsorenson/data/VIIRS/DNB/VNP02MOD.A2021204.0942*.nc')

    ##!#plot_VIIRS_figure(filename_n0722_dnb, band = 'DNB', ax = ax1, zoom = True)
    plot_VIIRS_figure(filename_d0722_dat, band = 'M05', ax = ax2, show_smoke = True, zoom = True)
    ##!#plot_VIIRS_figure(filename_n0723_dnb, band = 'DNB', ax = ax3, zoom = True)

    ##!#plot_VIIRS_figure(filename_n0722_dat, band = 'M15', ax = ax4, zoom = True)
    ##!#plot_VIIRS_figure(filename_d0722_dat, band = 'M15', ax = ax5, zoom = True)
    ##!#plot_VIIRS_figure(filename_n0723_dat, band = 'M15', ax = ax6, zoom = True)

    ##!## Plot the VIIRS data
    ##!## -------------------
    ##!#plot_VIIRS_scatter(filename_n0722_dat[0], filename_n0722_dnb[0], 'M15', 'DNB', \
    ##!#    ax = ax7, save = False)
    plot_VIIRS_scatter(filename_d0722_dat[0], filename_d0722_dat[0], 'M15', 'M05', \
        ax = ax8, save = False)
    ##!#plot_VIIRS_scatter(filename_n0723_dat[0], filename_n0723_dnb[0], 'M15', 'DNB', \
    ##!#    ax = ax9, save = False)

    plt.show()
    
