"""
  NAME:
    CrISLib

  PURPOSE:

"""

import numpy as np
import numpy.ma as ma
import sys
from netCDF4 import Dataset
from datetime import datetime,timedelta
from dateutil.relativedelta import relativedelta
import scipy.io
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
import time
import os

home_dir = os.environ['HOME']

sys.path.append(home_dir + '/')
#from python_lib import plot_trend_line, plot_subplot_label, plot_figure_text, \
#    nearest_gridpoint, aerosol_event_dict, init_proj, \
#    convert_radiance_to_temp
from python_lib import *

data_dir = home_dir + '/data/CrIS/'
datacrs = ccrs.PlateCarree()
mapcrs = ccrs.LambertConformal()
#mapcrs = ccrs.LambertConformal(central_longitude = -110, central_latitude = 40)

unit_dict = {
    'air_temp': 'K',\
    'dew_point': 'K',\
    'sfc_temp': 'K',\
    'spec_hum': 'g / kg', \
    'rel_hum': '%', \
    'gp_hgt': 'm',
    'wv': 'g/kg',
    'temp': 'K',
    'rh': '%'
}

label_dict = {
    'air_temp': 'Temperature',\
    'dew_point': 'Dew Point Temperature',\
    'sfc_temp': 'Temperature',\
    'spec_hum': 'Specific Humidity', \
    'rel_hum': 'Relative Humidity', \
    'gp_hgt': 'Geopotential Height',
    'wv': 'Mixing Ratio',
    'temp': 'Air Temperature', 
    'rh': 'Relative Humidity'
}

title_dict = {
    'air_temp': 'Air Temperature',\
    'dew_point': 'Dew Point Temperature',\
    'sfc_temp': 'Surface Skin Temperature',\
    'spec_hum': 'Specific Humidity', \
    'rel_hum': 'Relative Humidity', \
    'gp_hgt': 'Geopotential Height',
    'wv': 'Mixing Ratio',
    'temp': 'Air Temperature',
    'rh': 'Relative Humidity'
}

cris_loc_dict = {
    'mu': {
        'smoke_lat': 41.2328,
        'smoke_lon': -120.1568,
        'clear_lat1': 41.4879,
        'clear_lon1': -120.5508,
        'clear_lat2': 40.9871,
        'clear_lon2': -119.7220,
    },
    'md': {
        'smoke_lat': 41.0000,
        'smoke_lon': -120.4676,
        'clear_lat1': 41.3951,
        'clear_lon1': -121.0405,
        'clear_lat2': 40.6634,
        'clear_lon2': -120.1207,
    },
    'ml': {    
        'smoke_lat': 40.5672,
        'smoke_lon': -120.9731,
        'clear_lat1': 40.9128,
        'clear_lon1': -121.3236,
        'clear_lat2': 40.4173,
        'clear_lon2': -120.5538,
    },
    'low': {
        #smoke_lat = 40.1929
        #smoke_lon = -121.2456
        'smoke_lat': 40.3261,
        'smoke_lon': -121.0302,
        'clear_lat1': 40.6438,
        'clear_lon1': -121.5822,
        'clear_lat2': 39.9617,
        'clear_lon2': -120.4403,
    },
}

def find_plume_CrIS(filename):
    CrIS_data5  = readCrIS_granule(filename, band = 'M05', zoom = True)
    CrIS_data7  = readCrIS_granule(filename, band = 'M07', zoom = True)
    CrIS_data11 = readCrIS_granule(filename, band = 'M11', zoom = True)
    CrIS_data16 = readCrIS_granule(filename, band = 'M16', zoom = True)

    tmp_data = np.full(CrIS_data5['data'].shape, np.nan)
    in_cloud = (((CrIS_data5['data'] + CrIS_data7['data'] > 1.2) |\
                 (CrIS_data16['data'] < 265.)) |\
                ((CrIS_data5['data'] + CrIS_data7['data'] > 0.7) &\
                 (CrIS_data16['data'] < 285.)))
    tmp_data[in_cloud] = 2.
    tmp_data[~in_cloud] = 0.
    tmp_data[~in_cloud & (CrIS_data5['data'] - \
        CrIS_data11['data'] > 0.05) & (CrIS_data11['data'] > 0.05)] = 1
    mask_tmp_data = np.ma.masked_invalid(tmp_data)

    # For the hash data (smoke), remove any pixels that do not have
    # tmp data equal to 1
    # -------------------------------------------------------------
    hash_data   = np.ma.masked_where(mask_tmp_data != 1, mask_tmp_data)
    nohash_data = np.ma.masked_where(mask_tmp_data == 1, mask_tmp_data)

    # Free the memory
    # ---------------
    CrIS_data5  = None
    CrIS_data7  = None
    CrIS_data11 = None
    CrIS_data16 = None

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

# Send in the CrIS netCDF object, returns regridded lat, lon, and rad_sw, 
# rad_mw, and rad_lw arrays
# ----------------------------------------------------------------------
def regrid_CrIS(data, option = 0, sw_idx = 0, mw_idx = 0, lw_idx = 0):
  
    if(option == 0):
        lat    = np.full((data['rad_sw'].shape[0] * 3, data['rad_sw'].shape[1] * 3), np.nan)
        lon    = np.full((data['rad_sw'].shape[0] * 3, data['rad_sw'].shape[1] * 3), np.nan)
        rad_sw = np.full((data['rad_sw'].shape[0] * 3, data['rad_sw'].shape[1] * 3), np.nan)
        rad_mw = np.full((data['rad_sw'].shape[0] * 3, data['rad_sw'].shape[1] * 3), np.nan)
        rad_lw = np.full((data['rad_sw'].shape[0] * 3, data['rad_sw'].shape[1] * 3), np.nan)
        start = time.time() 

        for ii in range(data['rad_sw'].shape[0]):
            xi = (ii * 3)+ 1
            for jj in range(data['rad_sw'].shape[1]):
                yj = (jj * 3) + 1
                lat[xi-1:xi+2, yj-1:yj+2] = data['lat'][ii,jj,:].reshape((3,3))
                lon[xi-1:xi+2, yj-1:yj+2] = data['lon'][ii,jj,:].reshape((3,3))
                rad_sw[xi-1:xi+2, yj-1:yj+2] = data['rad_sw'][ii,jj,:,sw_idx].reshape((3,3))
                rad_mw[xi-1:xi+2, yj-1:yj+2] = data['rad_mw'][ii,jj,:,mw_idx].reshape((3,3))
                rad_lw[xi-1:xi+2, yj-1:yj+2] = data['rad_lw'][ii,jj,:,lw_idx].reshape((3,3))
        print(sw_idx, mw_idx, lw_idx)
        end = time.time()
        print("Time = ", end - start)

        out_dict = {}
        out_dict['lat'] = lat
        out_dict['lon'] = lon
        out_dict['rad_sw'] = rad_sw
        out_dict['rad_mw'] = rad_mw
        out_dict['rad_lw'] = rad_lw

        return out_dict

    else:
        print('here')
        start = time.time() 
        outdata = np.concatenate([ \
            np.concatenate([data[ii,jj,:,0].reshape((3,3)) \
            for jj in range(data.shape[1])], axis = 1) \
            for ii in range(data.shape[0])])
        end = time.time()
        print("Time = ", end - start)
    
        return outdata
   
def plot_CrIS_granule_test(CrIS_data, zoom = True, save = False):

    plt.close('all')
    mapcrs = init_proj('202107222110')

    fig = plt.figure()
    ax1 = fig.add_subplot(1,3,1, projection = mapcrs)
    ax2 = fig.add_subplot(1,3,2, projection = mapcrs)
    ax3 = fig.add_subplot(1,3,3, projection = mapcrs)

    ax1.pcolormesh(CrIS_data['lon'], CrIS_data['lat'], CrIS_data['rad_sw'], \
        transform = datacrs, shading = 'auto')
    ax2.pcolormesh(CrIS_data['lon'], CrIS_data['lat'], CrIS_data['rad_mw'], \
        transform = datacrs, shading = 'auto')
    ax3.pcolormesh(CrIS_data['lon'], CrIS_data['lat'], CrIS_data['rad_lw'], \
        transform = datacrs, shading = 'auto')

    ax1.coastlines()
    ax2.coastlines()
    ax3.coastlines()

    ax1.add_feature(cfeature.STATES)
    ax2.add_feature(cfeature.STATES)
    ax3.add_feature(cfeature.STATES)
    if(zoom):
        lat_lims1 = [39.5, 42.0]
        lon_lims1 = [-122.0, -119.5]
        ax1.set_extent([lon_lims1[0],lon_lims1[1],lat_lims1[0],\
            lat_lims1[1]],crs = datacrs)
        ax2.set_extent([lon_lims1[0],lon_lims1[1],lat_lims1[0],\
            lat_lims1[1]],crs = datacrs)
        ax3.set_extent([lon_lims1[0],lon_lims1[1],lat_lims1[0],\
            lat_lims1[1]],crs = datacrs)

    ax1.set_title(str(np.round(CrIS_data['lbda_sw'][CrIS_data['sw_idx']], 2)))
    ax2.set_title(str(np.round(CrIS_data['lbda_mw'][CrIS_data['mw_idx']], 2)))
    ax3.set_title(str(np.round(CrIS_data['lbda_lw'][CrIS_data['lw_idx']], 2)))

    plt.suptitle(CrIS_data['date'])

    plt.show()

def plotCrIS_layer(CrIS_volume, variable, layer_idx, pax = None, zoom = True, \
        save = False):

    in_ax = True
    if(pax is None):
        in_ax = False
        plt.close('all')
        mapcrs = init_proj('202107232155')
        fig = plt.figure(figsize=(6,6))
        pax = fig.add_subplot(1,1,1, projection = mapcrs)
   
    mesh = pax.pcolormesh(CrIS_volume['lon'], CrIS_volume['lat'], \
        CrIS_volume[variable][:,:,layer_idx], \
        transform = datacrs, shading = 'auto',\
        cmap = 'plasma')
    cbar = plt.colorbar(mesh, ax = pax, orientation='vertical',\
        pad=0.03, extend = 'both', label = label_dict[variable] + \
        ' [' + unit_dict[variable] + ']')
    pax.coastlines()
    pax.add_feature(cfeature.STATES)
    if((variable == 'spec_hum') | (variable == 'rel_hum')):
        z_var = 'air_pres_h2o'
    else:
        z_var = 'air_pres'
    pax.set_title(title_dict[variable] + '\n' + \
        str(int(CrIS_volume[z_var][layer_idx])) + ' mb')

    if(zoom):
        pax.set_extent([-122.0-0.1,-119.5+0.1,39.5-0.1,42.0+0.1],\
                       ccrs.PlateCarree())
    
    if(not in_ax): 
        if(save):
            print('name')
        else:
            plt.show()

def plotCrIS_profile(CrIS_volume, variable, slat, slon, pax = None, \
        save = False):

    in_ax = True
    if(pax is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure(figsize=(6,6))
        pax = fig.add_subplot(1,1,1)
 
    # Get the profile nearest to slat and slon
    # ----------------------------------------
    idxs = nearest_gridpoint(slat, slon, CrIS_volume['lat'], \
        CrIS_volume['lon'])
 
    if((variable == 'spec_hum') | (variable == 'rel_hum')):
        z_var = 'air_pres_h2o'
    else:
        z_var = 'air_pres'
    pax.plot(CrIS_volume[variable][idxs].flatten(), CrIS_volume[z_var][:])
    pax.set_xlabel(label_dict[variable] + ' [' + unit_dict[variable] + ']')
    
    pax.set_ylim(1000., 100.)
    pax.set_yscale('log')
 
    pax.set_title(title_dict[variable] + '\n' + \
        str(slat) + ' N, ' + str(slon) + ' E')

    if(not in_ax): 
        if(save):
            print('name')
        else:
            plt.show()

def plotCrIS_profile_2panel(date_str, variable, layer_idx, slat, slon, \
        zoom = True, save = False):
    
    # Read the data
    # -------------
    CrIS_volume = readCrIS_volume(date_str)

    # Set up the figure
    # -----------------
    plt.close('all')
    fig = plt.figure(figsize = (9, 4))
    mapcrs = init_proj('202107232155')
    ax1 = fig.add_subplot(1,2,1, projection = mapcrs)
    ax2 = fig.add_subplot(1,2,2)

    # Plot the map
    # ------------
    plotCrIS_layer(CrIS_volume, variable, layer_idx, pax = ax1, zoom = zoom, \
        save = False)
    
    # Add the point
    # -------------
    if(not isinstance(slat, list)):
        slat = [slat]
    if(not isinstance(slon, list)):
        slon = [slon]
    msize = 8 
    for tlat, tlon in zip(slat, slon):
        ax1.plot(tlon, tlat, linewidth=2, markersize = msize, marker='.',
                 color = 'black', transform=datacrs)
        ax1.plot(tlon, tlat, linewidth=2, markersize = msize - 2, marker='.',
                 transform=datacrs)

        # Plot the profile
        # ----------------
        plotCrIS_profile(CrIS_volume, variable, tlat, tlon, pax = ax2)
    
    plt.show()

# Add a VIIRS plot too
def plotCrIS_profile_3panel(date_str, variable, layer_idx, slat, slon, \
        viirs_channel = 'M15', zoom = True, save = False):
   
    if('/home/bsorenson/Research/VIIRS' not in sys.path):
        sys.path.append('/home/bsorenson/Research/VIIRS')
    from VIIRSLib import plot_VIIRS_figure

    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
 
    # Read the data
    # -------------
    CrIS_volume = readCrIS_volume(date_str)

    # Use the date string to get the VIIRS JPSS image
    # -----------------------------------------------
    if(viirs_channel == 'DNB'):
        filename = glob(dt_date_str.strftime(\
            '/home/bsorenson/data/VIIRS/DNB/JPSS/VJ102DNB.A%Y%j.%H%M*.nc'))
    else:
        filename = glob(dt_date_str.strftime(\
            '/home/bsorenson/data/VIIRS/DNB/JPSS/VJ102MOD.A%Y%j.%H%M*.nc'))
     

    # Set up the figure
    # -----------------
    plt.close('all')
    fig = plt.figure(figsize = (11, 4))
    mapcrs = init_proj('202107232155')
    ax1 = fig.add_subplot(1,3,1, projection = mapcrs)
    ax2 = fig.add_subplot(1,3,2, projection = mapcrs)
    ax3 = fig.add_subplot(1,3,3)

    plot_VIIRS_figure(filename, band = viirs_channel, ax = ax1, \
        ptitle = dt_date_str.strftime('%Y-%m-%d %H:%M'), zoom = zoom)

    # Plot the map
    # ------------
    plotCrIS_layer(CrIS_volume, variable, layer_idx, pax = ax2, zoom = zoom, \
        save = False)
    
    # Add the point
    # -------------
    if(not isinstance(slat, list)):
        slat = [slat]
    if(not isinstance(slon, list)):
        slon = [slon]
    msize = 8 
    for tlat, tlon in zip(slat, slon):
        ax1.plot(tlon, tlat, linewidth=2, markersize = msize, marker='.',
                 color = 'black', transform=datacrs)
        ax1.plot(tlon, tlat, linewidth=2, markersize = msize - 2, marker='.',
                 transform=datacrs)
        ax2.plot(tlon, tlat, linewidth=2, markersize = msize, marker='.',
                 color = 'black', transform=datacrs)
        ax2.plot(tlon, tlat, linewidth=2, markersize = msize - 2, marker='.',
                 transform=datacrs)

        # Plot the profile
        # ----------------
        plotCrIS_profile(CrIS_volume, variable, tlat, tlon, pax = ax3)
 
    ##!#print(ax3.get_yticks()) 
    ##!#print([f'{y:.2f}' for y in ax3.get_yticks()])
    ##!#ax3.set_yticks(ax3.get_yticks(), [f'{y:.2f}' for y in ax3.get_yticks()])
     
    fig.tight_layout()
    if(save):
        outname = 'cris_3panel_' + viirs_channel + '_' + variable + '_' + \
            date_str + '.png' 
        fig.savefig(outname, dpi = 300)
        print('Saved image', outname)
    else:
        plt.show()

# Reads one of Bill Smith Sr.'s CrIS retrieval files at a press level
# -------------------------------------------------------------------
def readCrIS_retrieval_level(date_str, press):

    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M%S')

    # Read in the file
    # ----------------
    data = scipy.io.loadmat(dt_date_str.strftime(home_dir + '/Research/CrIS/PHSS_npp_%Y%m%d_%H%M%S_rtv_m3F.mat'))
    #data = scipy.io.loadmat(dt_date_str.strftime('PHSS_npp_20210723_093719_rtv_m3F.mat')

    # To access all lat/lon
    lats  = np.squeeze(data['iovars']['plat'][0][0])
    lons  = np.squeeze(data['iovars']['plon'][0][0])
    plev  = np.squeeze(data['iovars']['plevs'][0][0])
    spres = np.squeeze(data['iovars']['spres'][0][0])

    # Find the pressure index closest to the desired level
    # ----------------------------------------------------
    p_idx = np.argmin(abs(plev - press)) 

    # Mask cloud-contaminated pixels (pixels with  no sfc temp)
    # ---------------------------------------------------------
    sfc_temp = np.squeeze(data['rtvdat']['sfct'][0][0])[:] 
    temps    = np.ma.masked_where(np.isnan(sfc_temp), \
        np.squeeze(data['rtvdat']['temp'][0][0])[:,p_idx])
    wv       = np.ma.masked_where(np.isnan(sfc_temp), \
        np.squeeze(data['rtvdat']['wv'][0][0])[:,p_idx])
    rh       = np.ma.masked_where(np.isnan(sfc_temp), \
        np.squeeze(data['rtvdat']['rh'][0][0])[:,p_idx])

    CrIS_data = {}
    CrIS_data['lat']       = lats
    CrIS_data['lon']       = lons
    CrIS_data['press']     = press
    CrIS_data['sfc_press'] = spres
    CrIS_data['sfc_temp']  = sfc_temp
    CrIS_data['temp']      = temps
    CrIS_data['wv']        = wv
    CrIS_data['rh']        = rh
  
    return CrIS_data

# Reads one of Bill Smith Sr.'s CrIS retrieval files at a press level
# -------------------------------------------------------------------
def readCrIS_retrieval_profile(date_str, plat, plon):

    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M%S')

    # Read in the file
    # ----------------
    data = scipy.io.loadmat(dt_date_str.strftime(home_dir + '/Research/CrIS/PHSS_npp_%Y%m%d_%H%M%S_rtv_m3F.mat'))
    #data = scipy.io.loadmat(dt_date_str.strftime('PHSS_npp_20210723_093719_rtv_m3F.mat')

    # To access all lat/lon
    lats  = np.squeeze(data['iovars']['plat'][0][0])
    lons  = np.squeeze(data['iovars']['plon'][0][0])
    plev  = np.squeeze(data['iovars']['plevs'][0][0])
    spres = np.squeeze(data['iovars']['spres'][0][0])

    # Find the pressure index closest to the desired level
    # ----------------------------------------------------
    loc_idx = nearest_gridpoint(plat, plon, lats, lons)
    #p_idx = np.argmin(abs(plev - press)) 

    CrIS_data = {}
    CrIS_data['plat']      = plat
    CrIS_data['plon']      = plon
    CrIS_data['lat']       = lats[loc_idx][0]
    CrIS_data['lon']       = lons[loc_idx][0]
    CrIS_data['sfc_press'] = spres[loc_idx][0]
    CrIS_data['press']     = plev
    CrIS_data['temp']      = np.squeeze(data['rtvdat']['temp'][0][0][loc_idx,:])
    CrIS_data['sfc_temp']  = np.squeeze(data['rtvdat']['sfct'][0][0])[loc_idx]
    CrIS_data['wv']        = np.squeeze(data['rtvdat']['wv'][0][0][loc_idx,:])
    CrIS_data['rh']        = np.squeeze(data['rtvdat']['rh'][0][0][loc_idx,:])
 
    return CrIS_data
 
def readCrIS_volume(date_str, zoom = True):

    data_dir = '/home/bsorenson/data/CrIS/CRIMSS/'
    
    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')

    filename = glob(dt_date_str.strftime(data_dir + 'SNDR.J1.CRIMSS.%Y%m%dT%H%M*.nc'))

    if(len(filename) == 0):
        print("ERROR: No files found for",date_str)
        return

    data = Dataset(filename[0],'r')

    CrIS_data = {}
    CrIS_data['air_temp']      = data['air_temp'][:,:,:]
    CrIS_data['sfc_temp']      = data['surf_air_temp'][:,:]
    CrIS_data['spec_hum']      = data['spec_hum'][:,:,:] * 1000.
    CrIS_data['rel_hum']       = data['rel_hum'][:,:,:] * 100.
    CrIS_data['gp_hgt']        = data['gp_hgt'][:,:,:]
    CrIS_data['air_pres']      = data['air_pres'][:] / 100.
    CrIS_data['air_pres_h2o']  = data['air_pres_h2o'][:] / 100.
    CrIS_data['lat']           = data['lat'][:,:]
    CrIS_data['lon']           = data['lon'][:,:]

    data.close()

    return CrIS_data

def readCrIS_granule(date_str, satellite = 'SNPP', resolution = 'NSR', \
        zoom = True):

    if(satellite == 'JPSS'):
        sat_add = 'J1'  
        print("WARNING: Only FSR available for JPSS")
        resolution = 'FSR'
    else:
        sat_add = 'SNPP'

    base_dir = data_dir + resolution + '/'

    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')

    # Use the date string (and glob) to find the matching file
    # --------------------------------------------------------
    base_filename = 'SNDR.' + sat_add + '.CRIS.'
    find_string = dt_date_str.strftime(base_dir + base_filename + '%Y%m%dT%H%M*')
    print(find_string)
    filename = glob(find_string)

    if(len(filename) == 0):
        print("ERROR: No CrIS files found for dtg",date_str)
        return
    else:
        print(filename)

   
    # Open and read the file
    # ----------------------
    data = Dataset(filename[0], 'r')

    if(resolution == 'NSR'):
        sw_idx = 150
        mw_idx = 100
        lw_idx = 420
    else:
        sw_idx = 350
        mw_idx = 400
        lw_idx = 422
    
    CrIS_data = regrid_CrIS(data, option = 0, sw_idx = sw_idx, mw_idx = mw_idx, \
        lw_idx = lw_idx)

    # data = Dataset(filename, 'r')
    wnum_sw = data['wnum_sw'][:].data
    lbda_sw = convert_wavenumber_to_wavelength(wnum_sw)
    wnum_mw = data['wnum_mw'][:].data
    lbda_mw = convert_wavenumber_to_wavelength(wnum_mw)
    wnum_lw = data['wnum_lw'][:].data
    lbda_lw = convert_wavenumber_to_wavelength(wnum_lw)
 
    CrIS_data['wnum_sw'] = wnum_sw
    CrIS_data['lbda_sw'] = lbda_sw
    CrIS_data['wnum_mw'] = wnum_mw
    CrIS_data['lbda_mw'] = lbda_mw
    CrIS_data['wnum_lw'] = wnum_lw
    CrIS_data['lbda_lw'] = lbda_lw

    CrIS_data['sw_idx'] = sw_idx
    CrIS_data['mw_idx'] = mw_idx
    CrIS_data['lw_idx'] = lw_idx

    CrIS_data['date'] = date_str
 
    return CrIS_data

 
    #CrIS_data = {}
    #CrIS_data['filename'] = filename
    #CrIS_data['band']     = band

    #return CrIS_data

def plot_CrIS_retrieval_level(CrIS_data, pvar, ax = None, labelsize = 12, \
        labelticksize = 10, markersize = 500, zoom = True, \
        show_smoke = False, vmin = None, vmax = None, alpha = 1.0, \
        plabel = '', save = False):

    in_ax = True
    if(ax is None):
        in_ax = False
        plt.close('all')
        mapcrs = init_proj('202107232155')
        fig = plt.figure(figsize=(6,6))
        ax = fig.add_subplot(1,1,1, projection = mapcrs)

    #if(vmax is None):
    #    vmax = CrIS_data['vmax']

    plot_data = np.ma.masked_invalid(CrIS_data[pvar])

    scat = ax.scatter(CrIS_data['lon'], CrIS_data['lat'], s = markersize, \
        c = plot_data, transform = ccrs.PlateCarree(), cmap = 'jet', \
        vmin = vmin, vmax = vmax, alpha = alpha)
    cbar = plt.colorbar(scat, ax = ax, orientation='vertical',\
        pad=0.03, extend = 'both')
    if(plabel == None):
        elabel = pvar
        
    cbar.set_label(plabel, fontsize = 12)

    ax.coastlines()
    ax.add_feature(cfeature.STATES)

    if(pvar == 'sfc_temp'):
        ax.set_title('SNPP CrIS\nSkin Temperature (K)')
    else:
        ax.set_title('SNPP CrIS\n' + str(int(CrIS_data['press'])) + ' hPa ' + \
            label_dict[pvar])

    #if(show_smoke):    
    #    hash_data, nohash_data = find_plume_CrIS(CrIS_data['filename'])
    #    #plt.rcParams.update({'hatch.color': 'r'})
    #    hatch_shape = '\\\\'
    #    #hash0 = ax2.pcolor(MODIS_data_ch1['lon'],MODIS_data_ch1['lat'],\
    #    #    hash_data[:MODIS_data_ch1['lat'].shape[0], :MODIS_data_ch1['lat'].shape[1]], hatch = hatch_shape, alpha=0., transform = datacrs)
    #    ax.pcolor(CrIS_data['lon'],CrIS_data['lat'],\
    #        hash_data, hatch = hatch_shape, alpha=0.40, transform = datacrs)
     

    if(zoom):
        ax.set_extent([-122.0,-119.5,39.5,42.0],\
                       ccrs.PlateCarree())
    if(not in_ax):
        plt.show()

def plot_CrIS_retrieval_profile(CrIS_data, pvar, lat, lon, ax = None, \
        labelsize = 12, labelticksize = 10, zoom = True, vmin = None, \
        vmax = None, ymin = 1000, ymax = 100, save = False):

    in_ax = True
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure(figsize=(6,6))
        ax = fig.add_subplot(1,1,1)

    ax.plot(np.ma.masked_invalid(CrIS_data[pvar]), CrIS_data['press'])
    ax.invert_yaxis()
    ax.set_ylim(ymax, ymin)
    #ax.set_xlim(vmin, vmax)
    ax.set_yscale('log')
    
    ax.set_title('SNPP CrIS\n' + title_dict[pvar] + ' Profile')

    if(not in_ax):
        plt.show()

def plot_CrIS_retrieval_combined(date_str, press = 500., pvar = 'wv',\
        file_adder = 'NP0', satellite = 'SNPP', row_str = 'md', \
        plot_skin_temp = False, alpha = 1.0, dot_size = 150, \
        save = False):

    if('/home/bsorenson/Research/VIIRS' not in sys.path):
        sys.path.append('/home/bsorenson/Research/VIIRS')
    from VIIRSLib import plot_VIIRS_figure
    
    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M%S')

    smoke_lat = cris_loc_dict[row_str]['smoke_lat']
    smoke_lon = cris_loc_dict[row_str]['smoke_lon']
    clear_lat1 = cris_loc_dict[row_str]['clear_lat1']
    clear_lon1 = cris_loc_dict[row_str]['clear_lon1']
    clear_lat2 = cris_loc_dict[row_str]['clear_lat2']
    clear_lon2 = cris_loc_dict[row_str]['clear_lon2']

    wv_max = 0.8
    wv_min = 0.0
    sfc_tmp_max = 330
    sfc_tmp_min = 280
    air_tmp_max = 305
    air_tmp_min = 285
    #pvar = 'sfc_temp'
  
    # Read the VIIRS data
    local_dt_date = dt_date_str + timedelta(minutes = 3)
    filename_d0722_dat = glob(local_dt_date.strftime(\
        home_dir + '/data/VIIRS/DNB/'+satellite+'/V'+file_adder + \
        '2MOD.A%Y%j.%H%M*.nc'))
    print(local_dt_date.strftime(home_dir + '/data/VIIRS/DNB/'+satellite+'/V'+file_adder +'2MOD.A%Y%j.%H%M*.nc'))
    namesplit = filename_d0722_dat[0].split('/')[-1].split('.')
    date1 = datetime.strptime('.'.join(namesplit[1:3])[1:], '%Y%j.%H%M')
 
    # Read the level data 
    CrIS_level_day   = readCrIS_retrieval_level(date_str,   press)
    CrIS_level_day2  = readCrIS_retrieval_level(date_str,   800.)
    
    # Read the daytime profiles
    CrIS_data_day_smoke  =  readCrIS_retrieval_profile(date_str, smoke_lat, smoke_lon)
    CrIS_data_day_clear1 =  readCrIS_retrieval_profile(date_str, clear_lat1, clear_lon1)
    CrIS_data_day_clear2 =  readCrIS_retrieval_profile(date_str, clear_lat2, clear_lon2)
   
    plt.close('all') 
    fig = plt.figure(figsize = (8, 12)) 
    mapcrs = init_proj('202107232155')
    ax1 = fig.add_subplot(3,2,1, projection = mapcrs)
    ax2 = fig.add_subplot(3,2,2, projection = mapcrs)
    ax3 = fig.add_subplot(3,2,3, projection = mapcrs)
    ax4 = fig.add_subplot(3,2,4, projection = mapcrs)
    ax5 = fig.add_subplot(3,2,5)
    ax6 = fig.add_subplot(3,2,6)
    
    # Plot the spatial data
    # ---------------------
    plot_VIIRS_figure(filename_d0722_dat, band = 'M05', ax = ax1, \
        ptitle = 'SNPP VIIRS\n0.67 μm Reflectance', zoom = True, vmax = 0.7)
    #plot_VIIRS_figure(filename_d0722_dat, band = 'M05', ax = ax2, \
    #    ptitle = date1.strftime('%Y-%m-%d %H:%M'), zoom = True, vmax = 0.7)

    point_size = 15
    #spatial_var = 'temp'
    spatial_var = pvar
    plot_CrIS_retrieval_level(CrIS_level_day, 'sfc_temp', ax = ax2, labelsize = 12, \
        labelticksize = 10, zoom = True, show_smoke = False, vmin = sfc_tmp_min, \
        vmax = sfc_tmp_max, save = False, markersize = dot_size, alpha = alpha, \
        plabel = 'Temperature (K)')
    plot_CrIS_retrieval_level(CrIS_level_day, spatial_var, ax = ax3, labelsize = 12, \
        labelticksize = 10, zoom = True, show_smoke = False, vmin = wv_min, \
        vmax = wv_max, save = False, markersize = dot_size, alpha = alpha, \
        plabel = 'Mixing Ratio (g/kg)')
    plot_CrIS_retrieval_level(CrIS_level_day2, 'temp', ax = ax4, labelsize = 12, \
        labelticksize = 10, zoom = True, show_smoke = False, vmin = air_tmp_min, \
        vmax = air_tmp_max, save = False, markersize = dot_size, alpha = alpha, \
        plabel = 'Temperature (K)')
    
    ax1.plot(smoke_lon, smoke_lat, linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax1.plot(smoke_lon, smoke_lat, linewidth=2, markersize = point_size, marker='.',
            transform=datacrs, color = 'tab:blue')
    ax1.plot(clear_lon1, clear_lat1, linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax1.plot(clear_lon1, clear_lat1, linewidth=2, markersize = point_size, marker='.',
            transform=datacrs, color = 'tab:orange')
    ax1.plot(clear_lon2, clear_lat2, linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax1.plot(clear_lon2, clear_lat2, linewidth=2, markersize = point_size, marker='.',
            transform=datacrs, color = 'tab:green')
    
    ax2.plot(smoke_lon, smoke_lat, linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax2.plot(smoke_lon, smoke_lat, linewidth=2, markersize = point_size, marker='.',
            transform=datacrs, color = 'tab:blue')
    ax2.plot(clear_lon1, clear_lat1, linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax2.plot(clear_lon1, clear_lat1, linewidth=2, markersize = point_size, marker='.',
            transform=datacrs, color = 'tab:orange')
    ax2.plot(clear_lon2, clear_lat2, linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax2.plot(clear_lon2, clear_lat2, linewidth=2, markersize = point_size, marker='.',
            transform=datacrs, color = 'tab:green')

    ax3.plot(smoke_lon, smoke_lat, linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax3.plot(smoke_lon, smoke_lat, linewidth=2, markersize = point_size, marker='.',
            transform=datacrs, color = 'tab:blue')
    ax3.plot(clear_lon1, clear_lat1, linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax3.plot(clear_lon1, clear_lat1, linewidth=2, markersize = point_size, marker='.',
            transform=datacrs, color = 'tab:orange')
    ax3.plot(clear_lon2, clear_lat2, linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax3.plot(clear_lon2, clear_lat2, linewidth=2, markersize = point_size, marker='.',
            transform=datacrs, color = 'tab:green')

    ax4.plot(smoke_lon, smoke_lat, linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax4.plot(smoke_lon, smoke_lat, linewidth=2, markersize = point_size, marker='.',
            transform=datacrs, color = 'tab:blue')
    ax4.plot(clear_lon1, clear_lat1, linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax4.plot(clear_lon1, clear_lat1, linewidth=2, markersize = point_size, marker='.',
            transform=datacrs, color = 'tab:orange')
    ax4.plot(clear_lon2, clear_lat2, linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax4.plot(clear_lon2, clear_lat2, linewidth=2, markersize = point_size, marker='.',
            transform=datacrs, color = 'tab:green')
    
    # Plot CrIS mixing ratio / relative humidity
    prof_var = 'wv'
    plot_CrIS_retrieval_profile(CrIS_data_day_smoke, prof_var, smoke_lat, smoke_lon, ax = ax5)
    plot_CrIS_retrieval_profile(CrIS_data_day_clear1, prof_var, clear_lat1, clear_lon1, ax = ax5)
    plot_CrIS_retrieval_profile(CrIS_data_day_clear2, prof_var, clear_lat2, clear_lon2, ax = ax5)
    
    prof_var = 'temp'
    plot_CrIS_retrieval_profile(CrIS_data_day_smoke, prof_var, smoke_lat, smoke_lon, ax = ax6)
    plot_CrIS_retrieval_profile(CrIS_data_day_clear1, prof_var, clear_lat1, clear_lon1, ax = ax6)
    plot_CrIS_retrieval_profile(CrIS_data_day_clear2, prof_var, clear_lat2, clear_lon2, ax = ax6)
   
    if(plot_skin_temp): 
        ax6.scatter(CrIS_data_day_smoke['sfc_temp'],CrIS_data_day_smoke['sfc_press'],  color = 'tab:blue')
        ax6.scatter(CrIS_data_day_clear1['sfc_temp'],CrIS_data_day_clear1['sfc_press'],  color = 'tab:orange')
        ax6.scatter(CrIS_data_day_clear2['sfc_temp'],CrIS_data_day_clear2['sfc_press'],  color = 'tab:green')
    
    ax5.axhline(press, linestyle = '--', color = 'k')
    ax6.axhline(press, linestyle = '--', color = 'k')
    
    # Highlight regions where blue wv is higher than green wv with cyan
    # Highlight regions where blue tp is lower than green tp with purple
    higher_wv = np.where(((CrIS_data_day_smoke['wv'] - CrIS_data_day_clear1['wv']) > 0.1) & \
                         ((CrIS_data_day_smoke['wv'] - CrIS_data_day_clear2['wv']) > 0.1))
    lower_tmp = np.where(((CrIS_data_day_clear1['temp'] - CrIS_data_day_smoke['temp']) > 0.1) & \
                         ((CrIS_data_day_clear2['temp'] - CrIS_data_day_smoke['temp']) > 0.1))
    
    tmp_idxs = np.split(lower_tmp[0],np.where(np.squeeze(np.diff(lower_tmp)) != 1)[0]+1)
   
    try: 
        ax5.axhspan(CrIS_data_day_smoke['press'][higher_wv[0][0]], CrIS_data_day_smoke['press'][higher_wv[0][-1]], color = 'cyan', alpha = 0.5)
        ax6.axhspan(CrIS_data_day_smoke['press'][higher_wv[0][0]], CrIS_data_day_smoke['press'][higher_wv[0][-1]], color = 'cyan', alpha = 0.5)

        for ii in range(len(tmp_idxs)):
            ax5.axhspan(CrIS_data_day_smoke['press'][tmp_idxs[ii][0]], CrIS_data_day_smoke['press'][tmp_idxs[ii][-1]], color = 'purple', alpha = 0.5)
            ax6.axhspan(CrIS_data_day_smoke['press'][tmp_idxs[ii][0]], CrIS_data_day_smoke['press'][tmp_idxs[ii][-1]], color = 'purple', alpha = 0.5)
    except IndexError:
        print("ERROR: unable to show shaded regions") 
    
    
    ax5.set_xlabel('Mixing Ratio (g/kg)')
    ax6.set_xlabel('Air Temperature (K)')
    ax5.set_ylabel('Pressure (hPa)')
    #ax6.set_ylabel('Pressure (hPa)')

    ax5.set_ylim(1000, 350)
    ax6.set_ylim(1000, 350)
     
    #ax3.set_xlim(250, 315)
    ax6.set_xlim(250, 315)

    plot_subplot_label(ax1, '(a)', backgroundcolor = 'white', location = 'upper_right', fontsize = 12)
    plot_subplot_label(ax2, '(b)', backgroundcolor = 'white', location = 'upper_right', fontsize = 12)
    plot_subplot_label(ax3, '(c)', backgroundcolor = 'white', location = 'upper_right', fontsize = 12)
    plot_subplot_label(ax4, '(d)', backgroundcolor = 'white', location = 'upper_right', fontsize = 12)
    plot_subplot_label(ax5, '(e)', backgroundcolor = 'white', location = 'upper_right', fontsize = 12)
    plot_subplot_label(ax6, '(f)', backgroundcolor = 'white', location = 'upper_right', fontsize = 12)
   
    plt.suptitle(dt_date_str.strftime('%d %B %Y %H:%M UTC'))
 
    fig.tight_layout()
   
    if(save):
        outname = dt_date_str.strftime('cris_rtv_' + row_str + '_%Y%m%d%H%M%S_v2.png')
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()


def plot_CrIS_granule(CrIS_data, ax = None, labelsize = 12, \
        labelticksize = 10, zoom = True, show_smoke = False, vmin = None, \
        vmax = None, save = False):

    in_ax = True
    if(ax is None):
        in_ax = False
        plt.close('all')
        mapcrs = init_proj('202107232155')
        fig = plt.figure(figsize=(6,6))
        ax = plt.axes(projection = mapcrs)

    if(vmax is None):
        vmax = CrIS_data['vmax']

    pdata = np.copy(CrIS_data['data'])


    mesh = ax.pcolormesh(CrIS_data['lon'], CrIS_data['lat'], pdata, \
        cmap = CrIS_data['cmap'], vmin = vmin, vmax = vmax, \
        transform = ccrs.PlateCarree(), shading = 'auto')
    cbar = plt.colorbar(mesh, ax = ax, orientation='vertical',\
        pad=0.03, extend = 'both')

    cbar.set_label(CrIS_data['label'], size = labelsize, weight = 'bold')
    cbar.ax.tick_params(labelsize = labelticksize)
    ax.set_title(CrIS_data['dtype'])

    if(show_smoke):    
        hash_data, nohash_data = find_plume_CrIS(CrIS_data['filename'])
        #plt.rcParams.update({'hatch.color': 'r'})
        hatch_shape = '\\\\'
        #hash0 = ax2.pcolor(MODIS_data_ch1['lon'],MODIS_data_ch1['lat'],\
        #    hash_data[:MODIS_data_ch1['lat'].shape[0], :MODIS_data_ch1['lat'].shape[1]], hatch = hatch_shape, alpha=0., transform = datacrs)
        ax.pcolor(CrIS_data['lon'],CrIS_data['lat'],\
            hash_data, hatch = hatch_shape, alpha=0.40, transform = datacrs)
     

    if(zoom):
        ax.set_extent([-122.0,-119.5,39.5,42.0],\
                       ccrs.PlateCarree())
    if(not in_ax):
        plt.show()

def plot_CrIS_figure(filename, band = 'M12', ax = None, show_smoke = False, \
        vmin = None, vmax = None, ptitle = '', zoom = True):

    in_ax = True
    if(ax is None):
        in_ax = False
        mapcrs = init_proj('202107232155')
        plt.close('all')
        fig = plt.figure(figsize = (6,6))
        ax = fig.add_subplot(1,1,1, projection = mapcrs)

    if(isinstance(filename, list)):

        for ff in filename:
            #plot_CrIS_granule_DNB(ff,ax, band = band)

            # Read the data for this granule
            # ------------------------------
            CrIS_data = readCrIS_granule(ff, band = band, zoom = zoom)
        
            # Plot the data for this granule
            # ------------------------------
            plot_CrIS_granule(CrIS_data, ax = ax, zoom = zoom, \
                vmin = vmin, vmax = vmax, show_smoke = show_smoke)

    elif(isinstance(filename, str)):
        #plot_CrIS_granule_DNB(filename,ax, band = band)

        # Read the data for this granule
        # ------------------------------
        CrIS_data = readCrIS_granule(filename, band = band, zoom = zoom)

        # Plot the data for this granule
        # ------------------------------
        plot_CrIS_granule(CrIS_data, ax = ax, zoom = zoom, \
                vmin = vmin, vmax = vmax, show_smoke = show_smoke)

    if(zoom):
        ax.set_extent([-122.0,-119.5,39.5,42.0],\
                       ccrs.PlateCarree())
       
    ax.set_title(ptitle, weight = 'bold')
 
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
def plot_CrIS_scatter(fname1, fname2, band1, band2, ax = None, save = False):
    in_ax = True
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure(figsize = (6,6))
        ax = fig.add_subplot(1,1,1)

    # Read the data for this granule
    # ------------------------------
    CrIS_data1 = readCrIS_granule(fname1, band = band1, zoom = True)
    CrIS_data2 = readCrIS_granule(fname2, band = band2, zoom = True)

    # Match the DNB data to the MOD data
    # ----------------------------------
    compress_dat1_lat  = CrIS_data1['lat'][~CrIS_data1['data'].mask]  
    compress_dat1_lon  = CrIS_data1['lon'][~CrIS_data1['data'].mask]  
    compress_dat1_data = CrIS_data1['data'][~CrIS_data1['data'].mask]  
    compress_dat2_lat  = CrIS_data2['lat'][~CrIS_data2['data'].mask]  
    compress_dat2_lon  = CrIS_data2['lon'][~CrIS_data2['data'].mask]  
    compress_dat2_data = CrIS_data2['data'][~CrIS_data2['data'].mask]  

    if(len(compress_dat1_data) == len(compress_dat2_data)):

        #ax.scatter(compress_dat1_data, compress_dat2_data, s = 6)
        ax.set_xlabel(CrIS_data1['dtype'])
        ax.set_ylabel(CrIS_data2['dtype']) 

        # Use find_plume to identify where the smoke is
        # ---------------------------------------------
        hash_data, nohash_data = find_plume_CrIS(fname1)
        #plt.rcParams.update({'hatch.color': 'r'})
        hatch_shape = '\\\\'
        #hash0 = ax2.pcolor(MODIS_data_ch1['lon'],MODIS_data_ch1['lat'],\
        #    hash_data[:MODIS_data_ch1['lat'].shape[0], :MODIS_data_ch1['lat'].shape[1]], hatch = hatch_shape, alpha=0., transform = datacrs)
        #ax.pcolor(CrIS_data1['lon'],CrIS_data1['lat'],\
        #    hash_data, hatch = hatch_shape, alpha=0, transform = datacrs)

        ax.scatter(CrIS_data1['data'][np.where(hash_data.mask)], \
            CrIS_data2['data'][np.where(hash_data.mask)], s=6, \
            color='tab:orange', label='Outside Plume')
        ax.scatter(CrIS_data1['data'][np.where(~hash_data.mask)], \
            CrIS_data2['data'][np.where(~hash_data.mask)], s=6, \
            color='tab:blue', label='Inside Plume')
       
    else: 
        if(len(compress_dat1_data) > len(compress_dat2_data)):
            print(CrIS_data1['dtype'], ' larger than ', CrIS_data2['dtype'])
            big_label    = CrIS_data1['dtype']
            small_label  = CrIS_data2['dtype']
            big_data   = compress_dat1_data 
            big_lat    = compress_dat1_lat 
            big_lon    = compress_dat1_lon 
            small_data = compress_dat2_data 
            small_lat  = compress_dat2_lat 
            small_lon  = compress_dat2_lon 
        else:
            print(CrIS_data2['dtype'], ' larger than ', CrIS_data1['dtype'])
            big_label    = CrIS_data2['dtype']
            small_label  = CrIS_data1['dtype']
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
        #ax.scatter(CrIS_data1['data'], CrIS_data2['data'], s = 6)
        #ax.set_xlabel(CrIS_data1['dtype'])
        #ax.set_ylabel(CrIS_data2['dtype']) 
 
    if(not in_ax):
        plt.show() 


def plot_CrIS_sixpanel(save = False):

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

    # Read the CrIS data for each time
    # ---------------------------------
    #filename_n0722_dnb = glob('/home/bsorenson/data/CrIS/DNB/VNP02DNB.A2021203.1000*.nc')
    #filename_n0722_dat = glob('/home/bsorenson/data/CrIS/DNB/VNP02MOD.A2021203.1000*.nc')
    #filename_n0722_dnb = glob('/home/bsorenson/data/CrIS/DNB/VNP02DNB.A2021203.1000*.nc')
    #filename_n0722_dat = glob('/home/bsorenson/data/CrIS/DNB/VNP02MOD.A2021203.1000*.nc')
    filename_d0722_dat = glob('/home/bsorenson/data/CrIS/DNB/VNP02MOD.A2021203.2124*.nc')
    filename_n0723_dnb = glob('/home/bsorenson/data/CrIS/DNB/VNP02DNB.A2021204.0942*.nc')
    filename_n0723_dat = glob('/home/bsorenson/data/CrIS/DNB/VNP02MOD.A2021204.0942*.nc')
    filename_d0723_dat = glob('/home/bsorenson/data/CrIS/DNB/VNP02MOD.A2021204.2100*.nc')

    namesplit = filename_d0722_dat[0].split('/')[-1].split('.')
    date1 = datetime.strptime('.'.join(namesplit[1:3])[1:], '%Y%j.%H%M')
    namesplit = filename_n0723_dat[0].split('/')[-1].split('.')
    date2 = datetime.strptime('.'.join(namesplit[1:3])[1:], '%Y%j.%H%M')
    namesplit = filename_d0723_dat[0].split('/')[-1].split('.')
    date3 = datetime.strptime('.'.join(namesplit[1:3])[1:], '%Y%j.%H%M')

    print(date1.strftime('%Y-%m-%d %H:%M'))
    print(date2.strftime('%Y-%m-%d %H:%M'))
    print(date3.strftime('%Y-%m-%d %H:%M'))

    plot_CrIS_figure(filename_d0722_dat, band = 'M05', ax = ax1, \
        ptitle = date1.strftime('%Y-%m-%d %H:%M'), zoom = True)
    plot_CrIS_figure(filename_n0723_dnb, band = 'DNB', ax = ax2, \
        ptitle = date2.strftime('%Y-%m-%d %H:%M'), zoom = True)
    plot_CrIS_figure(filename_d0723_dat, band = 'M05', ax = ax3,\
        ptitle = date3.strftime('%Y-%m-%d %H:%M'), zoom = True)

    plot_CrIS_figure(filename_d0722_dat, band = 'M15', ax = ax4, zoom = True, vmin = 260, vmax = 330)
    plot_CrIS_figure(filename_n0723_dat, band = 'M15', ax = ax5, zoom = True, vmin = 260, vmax = 330)
    plot_CrIS_figure(filename_d0723_dat, band = 'M15', ax = ax6, zoom = True, vmin = 260, vmax = 330)

    plot_subplot_label(ax1, '(a)', backgroundcolor = 'white')
    plot_subplot_label(ax2, '(b)', backgroundcolor = 'white')
    plot_subplot_label(ax3, '(c)', backgroundcolor = 'white')
    plot_subplot_label(ax4, '(d)', backgroundcolor = 'white')
    plot_subplot_label(ax5, '(e)', backgroundcolor = 'white')
    plot_subplot_label(ax6, '(f)', backgroundcolor = 'white')

    plot_figure_text(ax1, 'CrIS 0.67 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax2, 'CrIS DNB', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax3, 'CrIS 0.67 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax4, 'CrIS 10.76 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax5, 'CrIS 10.76 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax6, 'CrIS 10.76 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')

    fig.tight_layout()

    if(save):
        outname = 'viirs_sixpanel_' + date_str + '_v2.png'
        plt.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else: 
        plt.show()
    
def plot_CrIS_ninepanel(save = False):

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

    # Read the CrIS data for each time
    # ---------------------------------
    filename_n0722_dnb = glob('/home/bsorenson/data/CrIS/DNB/VNP02DNB.A2021203.1000*.nc')
    filename_n0722_dat = glob('/home/bsorenson/data/CrIS/DNB/VNP02MOD.A2021203.1000*.nc')
    filename_d0722_dat = glob('/home/bsorenson/data/CrIS/DNB/VNP02MOD.A2021203.2124*.nc')
    filename_n0723_dnb = glob('/home/bsorenson/data/CrIS/DNB/VNP02DNB.A2021204.0942*.nc')
    filename_n0723_dat = glob('/home/bsorenson/data/CrIS/DNB/VNP02MOD.A2021204.0942*.nc')

    ##!#plot_CrIS_figure(filename_n0722_dnb, band = 'DNB', ax = ax1, zoom = True)
    plot_CrIS_figure(filename_d0722_dat, band = 'M05', ax = ax2, show_smoke = True, zoom = True)
    ##!#plot_CrIS_figure(filename_n0723_dnb, band = 'DNB', ax = ax3, zoom = True)

    ##!#plot_CrIS_figure(filename_n0722_dat, band = 'M15', ax = ax4, zoom = True)
    ##!#plot_CrIS_figure(filename_d0722_dat, band = 'M15', ax = ax5, zoom = True)
    ##!#plot_CrIS_figure(filename_n0723_dat, band = 'M15', ax = ax6, zoom = True)

    ##!## Plot the CrIS data
    ##!## -------------------
    ##!#plot_CrIS_scatter(filename_n0722_dat[0], filename_n0722_dnb[0], 'M15', 'DNB', \
    ##!#    ax = ax7, save = False)
    plot_CrIS_scatter(filename_d0722_dat[0], filename_d0722_dat[0], 'M15', 'M05', \
        ax = ax8, save = False)
    ##!#plot_CrIS_scatter(filename_n0723_dat[0], filename_n0723_dnb[0], 'M15', 'DNB', \
    ##!#    ax = ax9, save = False)

    plt.show()
    
