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
import subprocess
from scipy import stats
from pyhdf import SD
import pandas as pd
import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from satpy.scene import Scene
from satpy import find_files_and_readers
from satpy.writers import get_enhanced_image
from glob import glob
import os

home_dir = os.environ['HOME']

sys.path.append(home_dir + '/')
viirs_dir = home_dir  + '/data/VIIRS/DNB/'
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

sat_changer = {
    'SNPP': 'NP0',\
    'JPSS': 'J10'
}

# param: MOD or DNB
def identify_VIIRS_file(date_str, param, viirs_dir, satellite = 'SNPP'):

    local_data_dir  = viirs_dir + satellite + '/'

    base_url = 'https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5200/V' + \
        sat_changer[satellite] + '2' + param + '/'
    base_url_geo = 'https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5200/V' + \
        sat_changer[satellite] + '3' + param + '/'

    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')

    # For each desired channel, figure out the closest file time
    # to the input date
    # ----------------------------------------------------------
    try:
        files = listFD(dt_date_str.strftime(base_url + '/%Y/%j/'), ext = '.nc')
        geo_files = listFD(dt_date_str.strftime(base_url_geo + '/%Y/%j/'), ext = '.nc')
    except subprocess.CalledProcessError:
        print("ERROR: No VIIRS files for the input DTG",date_str)
        return -2

    if(len(files) == 0):
        print("ERROR: No VIIRS files returned from the request. Exiting")
        return -1
    
    # Remove the timestamps from the file strings
    # -------------------------------------------
    files_only = [tfile.strip().split('/')[-1] for tfile in files]
    geo_files_only = [tfile.strip().split('/')[-1] for tfile in geo_files]

    # Use the actual file times to get timestamps
    # -------------------------------------------
    file_dates = [datetime.strptime(tfile[10:22],'%Y%j.%H%M') for tfile in files_only]
    geo_file_dates = [datetime.strptime(tfile[10:22],'%Y%j.%H%M') for tfile in geo_files_only]

    # Figure out where the closest files to the desired time are
    # ----------------------------------------------------------
    time_diffs = np.array([abs((dt_date_str - ddate).total_seconds()) \
        for ddate in file_dates])

    # Extract the index of the matching MODIS file
    # --------------------------------------------
    file_idx = np.argmin(time_diffs)
    found_file = files_only[file_idx]
    geo_file = geo_files_only[file_idx]

    # Check if the file is already downloaded
    # ---------------------------------------
    if(os.path.exists(local_data_dir + found_file)):
        print(found_file + ' already exists. Not downloading')
    else: 
        # Download the data file
        # ----------------------
        cmnd = dt_date_str.strftime("wget \"" + base_url + '/%Y/%j/' + found_file + \
            "\" --header \"Authorization: Bearer " + laads_daac_key + "\" -P .")
        print(cmnd)
        os.system(cmnd)

        # Move the file to the destination folder
        # ---------------------------------------
        cmnd = "mv " + found_file + " " + local_data_dir
        print(cmnd) 
        os.system(cmnd)

        # Now, grab the mathing MOD file 
        # ------------------------------
        #geo_file = 'V' + sat_changer[satellite] + '3' + param + '.' + found_file[10:] 
        print('\nHEERE\n',base_url_geo, geo_file)
        cmnd = dt_date_str.strftime("wget \"" + base_url_geo + '/%Y/%j/' + geo_file + \
            "\" --header \"Authorization: Bearer " + laads_daac_key + "\" -P .")
        print(cmnd)
        os.system(cmnd)

        # Move the file to the destination folder
        # ---------------------------------------
        cmnd = "mv " + geo_file + " " + local_data_dir
        print(cmnd) 
        os.system(cmnd)

    return found_file 


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

def readVIIRS_granule(filename, band = 'M15', \
        zoom = True):
  
    VIIRS_data = {}
    VIIRS_data['filename'] = filename
    VIIRS_data['band']     = band
    # Read the filename to determine the file type
    # --------------------------------------------
    total_split = filename.split('/')
    name_split = total_split[-1].split('.')
    dataset_name = name_split[0]
    sat_id = dataset_name[:5]
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

    #else:
    #elif(dataset_name[:5] == 'V'+file_adder + '2'):
    elif(dataset_name[4:5] == '2'):

        # Extract whether DNB or MOD
        # --------------------------
        data_type = dataset_name[5:]

        # Determine if geolocation data is present in same path
        # -----------------------------------------------------
        print("Searching for geolocation data for",'/'.join(total_split[:-1]) \
            + '/'+sat_id + data_type + '.' + '.'.join(name_split[1:4]))
        geoloc_list = glob('/'.join(total_split[:-1]) + '/'+ sat_id[:4] + \
            '3' + data_type + '.' + '.'.join(name_split[1:4]) + '*') 
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
                if(str(band) == 'true_color'):
                    dtype = 'true_color'
                    # Pull out the red, green, and blue reflectance data
                    # --------------------------------------------------
                    red = viirs_ds['observation_data/M05'][:]
                    red_scale = viirs_ds['observation_data/M05'].attrs['scale_factor']
                    red_offset = viirs_ds['observation_data/M05'].attrs['add_offset']
                    
                    green = viirs_ds['observation_data/M04'][:]
                    green_scale = viirs_ds['observation_data/M04'].attrs['scale_factor']
                    green_offset = viirs_ds['observation_data/M04'].attrs['add_offset']
                    
                    blue = viirs_ds['observation_data/M03'][:]
                    blue_scale = viirs_ds['observation_data/M03'].attrs['scale_factor']
                    blue_offset = viirs_ds['observation_data/M03'].attrs['add_offset']
                   
                    """ 
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
                    
                    ##!## Close the satellite file object
                    ##!## -------------------------------
                    ##!#modis.end()
                    """ 
                    # The RGB reflectances are on a much higher resolution than the lats and
                    
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
                    
                    dnb = np.ma.masked_where(np.isnan(image),image)
                    label = 'True color'
                    colors = 'True color'
                    vmax = None
                    cmap = None
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
        

    VIIRS_data['filename']  = filename
    VIIRS_data['satellite'] = filename
    VIIRS_data['dtype']     = dtype
    VIIRS_data['label']     = label
    VIIRS_data['lon']       = x
    VIIRS_data['lat']       = y
    VIIRS_data['data']      = dnb
    VIIRS_data['vmax']      = vmax
    VIIRS_data['cmap']      = cmap
    if(str(band) == 'true_color'):
        VIIRS_data['colortuple'] = colortuple
    else:
        VIIRS_data['colortuple'] = None

    return VIIRS_data

# dt_date_str is of format YYYYMMDDHHMM
def read_VIIRS_channel(date_str, channel, zoom = False, swath = False, \
        satellite = 'SNPP', add_time = 5, lat_lims = None, lon_lims = None):

    data_dir = viirs_dir + satellite + '/'
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    if(channel == 'DNB'):
        # DNB file reading doesn't work with Satpy file finder.
        # Just grab the file given the DTG passed
        # -----------------------------------------------------
        filename = glob(dt_date_str.strftime(data_dir + '*02' + channel + '.A%Y%j.%H%M*.nc'))[0]

        VIIRS_final = readVIIRS_granule(filename, channel, zoom = zoom)
        
        # If desired, clip the data to the lat/lon ranges
        # ----------------------------------------------- 
        if(lat_lims is not None):
            mask_lat = np.ma.masked_where((VIIRS_final['lat'] < lat_lims[0]) | \
                                          (VIIRS_final['lat'] > lat_lims[1]) | \
                                          (VIIRS_final['lon'] < lon_lims[0]) | \
                                          (VIIRS_final['lon'] > lon_lims[1]), \
                                           VIIRS_final['lat'])
            
            keep_long_idx = np.array([False in mask_lat.mask[ii,:] for ii in range(mask_lat.shape[0])])
            keep_short_idx = np.array([False in mask_lat.mask[:,ii] for ii in range(mask_lat.shape[1])])

            VIIRS_final['lat']  = VIIRS_final['lat'][keep_long_idx, :][:,keep_short_idx]
            VIIRS_final['lon']  = VIIRS_final['lon'][keep_long_idx, :][:,keep_short_idx]
            VIIRS_final['data'] = VIIRS_final['data'][keep_long_idx, :][:,keep_short_idx]
 
    else:
        # Extract the filename given the channel
        # --------------------------------------
        dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
        dt_date_str_beg = datetime.strptime(date_str,"%Y%m%d%H%M")
        dt_date_str_end = dt_date_str + timedelta(minutes = add_time)

        files = find_files_and_readers(start_time = dt_date_str_beg, \
            end_time = dt_date_str_end, base_dir = data_dir, \
            reader = 'viirs_l1b')

        # Grab only the 02MOD files
        # -------------------------
        is_02MOD = ['02MOD' in tfile for tfile in files['viirs_l1b']]
        filename = sorted(list(np.array(files['viirs_l1b'])[is_02MOD]))
 
        VIIRS_holder = {}
        counter = 0
        for ii, ifile in enumerate(filename):

            print("Reading VIIRS channel",channel," from ",ifile)

            VIIRS_data = readVIIRS_granule(ifile, channel, zoom = zoom)

            VIIRS_holder[str(counter)] = VIIRS_data 

            del VIIRS_data

            counter += 1

        print("HERE:", VIIRS_holder.keys())
        VIIRS_final = {}
        VIIRS_final['lat']  = np.concatenate([VIIRS_holder[key]['lat']  for key in VIIRS_holder.keys()], axis = 0)
        VIIRS_final['lon']  = np.concatenate([VIIRS_holder[key]['lon']  for key in VIIRS_holder.keys()], axis = 0)
        VIIRS_final['data'] = np.concatenate([VIIRS_holder[key]['data'] for key in VIIRS_holder.keys()], axis = 0)
        VIIRS_final['cmap'] = VIIRS_holder['0']['cmap']
        #VIIRS_final['vmax'] = np.max([VIIRS_holder[key]['vmax'] for key in VIIRS_holder.keys()])
        VIIRS_final['vmax'] = VIIRS_holder['0']['vmax']
        VIIRS_final['dtype'] = VIIRS_holder['0']['dtype']
        VIIRS_final['satellite'] = VIIRS_holder['0']['satellite']
        VIIRS_final['filename'] = filename

        # If desired, clip the data to the lat/lon ranges
        # ----------------------------------------------- 
        if(lat_lims is not None):
            mask_lat = np.ma.masked_where((VIIRS_final['lat'] < lat_lims[0]) | \
                                          (VIIRS_final['lat'] > lat_lims[1]) | \
                                          (VIIRS_final['lon'] < lon_lims[0]) | \
                                          (VIIRS_final['lon'] > lon_lims[1]), \
                                           VIIRS_final['lat'])
            
            keep_long_idx = np.array([False in mask_lat.mask[ii,:] for ii in range(mask_lat.shape[0])])
            keep_short_idx = np.array([False in mask_lat.mask[:,ii] for ii in range(mask_lat.shape[1])])

            VIIRS_final['lat']  = VIIRS_final['lat'][keep_long_idx, :][:,keep_short_idx]
            VIIRS_final['lon']  = VIIRS_final['lon'][keep_long_idx, :][:,keep_short_idx]
            VIIRS_final['data'] = VIIRS_final['data'][keep_long_idx, :, :][:,keep_short_idx]
        
        print(VIIRS_final['data'].shape, VIIRS_final['lat'].shape, VIIRS_final['lon'].shape)

        if(str(channel) == 'true_color'):
        
            # Convert the color values into a format usable for plotting
            # ----------------------------------------------------------
            colortuple = tuple(np.array([VIIRS_final['data'][:,:,0].flatten(), \
                VIIRS_final['data'][:,:,1].flatten(), \
                VIIRS_final['data'][:,:,2].flatten()]).transpose().tolist())
            VIIRS_final['colortuple'] = colortuple
        else: 
            VIIRS_final['colortuple'] = None

    print('VIIRS LAT SHAPE: ',VIIRS_final['lat'].shape)

    return VIIRS_final


# SLIGHTLY DIFFERENT THAN read_MODIS_satpy
# add_time: when wanting to grab a swath, specify the extra time desired, in minutes.
#           By default, this value is 5, which should make the code 
#           grab only single files
# ----------------------------------------------------------------------------
def read_VIIRS_satpy(date_str, channel, satellite = 'SNPP', \
        composite = False, swath = False, \
        zoom = True, return_xy = False, add_time = 5):
#def read_true_color(date_str,composite = False):

    data_dir = viirs_dir + satellite + '/'
    
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    dt_date_str_beg = datetime.strptime(date_str,"%Y%m%d%H%M")
    dt_date_str_end = dt_date_str + timedelta(minutes = add_time)

    files = find_files_and_readers(start_time = dt_date_str_beg, \
        end_time = dt_date_str_end, base_dir = data_dir, \
        reader = 'viirs_l1b')
   
    """ 
    if(swath):

        # Determine the correct GOES files associated with the date
        # ---------------------------------------------------------
        dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
        #dt_date_str_beg = dt_date_str - timedelta(minutes = 10) datetime.strptime(date_str,"%Y%m%d%H%M")
        dt_date_str_beg = datetime.strptime(date_str,"%Y%m%d%H%M")
        dt_date_str_end = dt_date_str + timedelta(minutes = add_time)

        cmpst_add = '_composite'
        #lat_lims = [60, 90]
        #lon_lims = [-180, 180]

        try:
            # Use the Satpy find_files_and_readers to grab the files
            # ------------------------------------------------------
            files = find_files_and_readers(start_time = dt_date_str_beg, \
                end_time = dt_date_str_end, base_dir = data_dir, \
                reader = 'viirs_l1b')
        except ValueError:
            print("ERROR: no files found for dtg",date_str)
            return

        day_filenames = files
        #lat_lims = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis_Lat']
        #lon_lims = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis_Lon']

    else:
        # Determine the correct MODIS file associated with the date
        # ---------------------------------------------------------
        dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
        filename = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis']
        print(filename)

        if(type(filename) is list):
            day_filenames = filename
            cmpst_add = ''
        else:
            if(composite):
                day_filenames = glob(filename[:50]+'*')
                cmpst_add = '_composite'
            else:
                day_filenames = glob(filename)
                cmpst_add = ''

        # Extract the modis true-color plot limits
        # ----------------------------------------
        lat_lims = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis_Lat']
        lon_lims = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis_Lon']
    """ 

    # Use satpy (Scene) to open the file
    # ----------------------------------
    scn = Scene(reader = 'viirs_l1b', filenames = files)

    # Load true-color data
    scn.load(['true_color'])
    if(channel != 'true_color'): 
        scn.load([str(channel)])
        # NOTE: Removing the "calibration" section for VIIRS data
        #scn.load([str(channel)], calibration = [channel_dict[str(channel)]['Unit_name']])
    ##!#else:
    ##!#    scn.load([str(channel)])

    ##!#if(channel == 'true_color'):
    ##!#    # Set the map projection and center the data
    ##!#    # ------------------------------------------
    #my_area = scn[str(channel)].attrs['area'].compute_optimal_bb_area({\
    #if(not swath):
    my_area = scn['true_color'].attrs['area'].compute_optimal_bb_area()
    #my_area = scn['true_color'].attrs['area'].compute_optimal_bb_area({\
    #    'proj':'lcc', 'lon_0': lon_lims[0], 'lat_0': lat_lims[0], \
    #    'lat_1': lat_lims[0], 'lat_2': lat_lims[0]})
    new_scn = scn.resample(my_area)
    #else:
    #    new_scn = scn

    ##!#if(zoom):
    ##!#    scn = scn.crop(ll_bbox = (lon_lims[0] + 0.65, lat_lims[0], \
    ##!#        lon_lims[1] - 0.65, lat_lims[1]))

    # Extract the lats, lons, and data
    # -----------------------------------------------------
    lons, lats = scn[str(channel)].attrs['area'].get_lonlats()
    #lons, lats = new_scn[str(channel)].attrs['area'].get_lonlats()

    if(channel == 'true_color'):
        # Enhance the image for plotting
        # ------------------------------
        var = get_enhanced_image(new_scn[str(channel)]).data
        var = var.transpose('y','x','bands')
    else:
        var = scn[str(channel)].data
    #var = new_scn[str(channel)].data

    if(not swath):
        # Extract the map projection from the data for plotting
        # -----------------------------------------------------
        crs = new_scn['true_color'].attrs['area'].to_cartopy_crs()
        #crs = new_scn[str(channel)].attrs['area'].to_cartopy_crs()
    else:
        #crs = ccrs.NorthPolarStereo()
        crs = new_scn['true_color'].attrs['area'].to_cartopy_crs()


    ##!#if(channel != 'true_color'):
    ##!#    var = var.data

    # Extract the appropriate units
    # -----------------------------
    if(channel == 'true_color'):
        plabel = ''
    else:
        plabel = channel
        #plabel = channel_dict[str(channel)]['Unit_name'].title() + ' [' + \
        #    channel_dict[str(channel)]['Unit'] + ']'

    if(return_xy):

        y = new_scn[channel].y 
        x = new_scn[channel].x 
        del scn 
        del new_scn
        return var, crs, lons, lats, plabel, x, y
        #return var, crs, lons, lats, lat_lims, lon_lims, plabel, x, y
    else:
        del scn 
        del new_scn
        return var, crs, lons, lats, plabel
        #return var, crs, lons, lats, lat_lims, lon_lims, plabel
    #return var, crs, lat_lims, lon_lims


def plot_VIIRS_granule(VIIRS_data, ax = None, labelsize = 12, \
        labelticksize = 10, zoom = True, show_smoke = False, vmin = None, \
        vmax = None, colorbar = True, show_title = True, save = False):

    in_ax = True
    if(ax is None):
        in_ax = False
        plt.close('all')
        mapcrs = init_proj('202107232155')
        fig = plt.figure(figsize=(6,6))
        ax = plt.axes(projection = mapcrs)

    if(vmax is None):
        vmax = VIIRS_data['vmax']

    if(VIIRS_data['colortuple'] is not None):
        ax.pcolormesh(VIIRS_data['lon'],VIIRS_data['lat'],\
            VIIRS_data['data'][:,:,0],color= VIIRS_data['colortuple'], \
            shading='auto', transform = ccrs.PlateCarree()) 
    else:
        pdata = np.copy(VIIRS_data['data'])
        pdata = np.ma.masked_where((pdata < 0) | (pdata > 500), pdata)

        print('MIN LON', np.min(VIIRS_data['lon']), 'MAX LON', np.max(VIIRS_data['lat']))
        mesh = ax.pcolormesh(VIIRS_data['lon'], VIIRS_data['lat'], pdata, \
            cmap = VIIRS_data['cmap'], vmin = vmin, vmax = vmax, \
            transform = ccrs.PlateCarree(), shading = 'auto')
        if(colorbar):
            cbar = plt.colorbar(mesh, ax = ax, orientation='vertical',\
                pad=0.03, shrink = 1.00, extend = 'max')

            cbar.set_label(VIIRS_data['label'], size = labelsize)
            cbar.ax.tick_params(labelsize = labelticksize)
        
        if(show_title):
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

def plot_VIIRS_figure(filename, satellite = 'SNPP', band = 'M12', \
        ax = None, show_smoke = False, \
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
            #plot_VIIRS_granule_DNB(ff,ax, band = band)

            # Read the data for this granule
            # ------------------------------
            VIIRS_data = readVIIRS_granule(ff, \
                band = band, zoom = zoom)
        
            # Plot the data for this granule
            # ------------------------------
            plot_VIIRS_granule(VIIRS_data, ax = ax, zoom = zoom, \
                vmin = vmin, vmax = vmax, show_smoke = show_smoke)

    elif(isinstance(filename, str)):
        #plot_VIIRS_granule_DNB(filename,ax, band = band)

        # Read the data for this granule
        # ------------------------------
        VIIRS_data = readVIIRS_granule(filename,\
            band = band, zoom = zoom)

        # Plot the data for this granule
        # ------------------------------
        plot_VIIRS_granule(VIIRS_data, ax = ax, zoom = zoom, \
                vmin = vmin, vmax = vmax, show_smoke = show_smoke)

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

# channel must be an integer between 1 and 16
# Am requiring lat_lims and lon_lims here 
def plot_VIIRS_satpy(date_str, channel, lat_lims, lon_lims, \
        ax = None, var = None, crs = None, lons = None, lats = None, \
        vmin = None, vmax = None, ptitle = None, plabel = None, \
        use_xy = False, satellite = 'SNPP', add_time = 5, \
        labelsize = 10, colorbar = True, swath = False, zoom=True,save=False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    if(var is None): 
        if(use_xy):
            var, crs, lons, lats, plabel, xx, yy = read_VIIRS_satpy(date_str, str(channel), \
                    satellite = satellite, zoom = True, return_xy = True, \
                    add_time = add_time)
        else:
            var, crs, lons, lats, plabel= read_VIIRS_satpy(date_str, str(channel), \
                    satellite = satellite, zoom = True, return_xy = False, \
                    add_time = add_time)
    
    # Plot the GOES data
    # ------------------
    in_ax = True 
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure()
        if(use_xy):
            mapcrs = init_proj(None)
            ax = fig.add_subplot(1,1,1, projection=mapcrs)
        else:
            ax = fig.add_subplot(1,1,1, projection=crs)


    ##!#ax.imshow(var.data, transform = crs, extent=(var.x[0], var.x[-1], \
    ##!#    var.y[-1], var.y[0]), vmin = vmin, vmax = vmax, origin='upper', \
    ##!#    cmap = 'Greys_r')
    if(channel == 'true_color'):
        ax.imshow(var.data, transform = crs, extent=(var.x[0], var.x[-1], \
            var.y[-1], var.y[0]), origin='upper')
    else:
        if(use_xy):
            im1 = ax.pcolormesh(xx, yy, var, transform = crs, \
                vmin = vmin, vmax = vmax, \
                cmap = channel_dict[str(channel)]['cmap'], \
                shading = 'auto')
        else:
            #im1 = ax.imshow(var, transform = crs, vmin = vmin, vmax = vmax, \
            im1 = ax.pcolormesh(lons, lats, var, transform = datacrs, \
                vmin = vmin, vmax = vmax, \
                cmap = channel_dict[str(channel)]['cmap'], \
                shading = 'auto')
        if(colorbar):
            cbar = plt.colorbar(im1, ax = ax, pad = 0.03, fraction = 0.052, \
                extend = 'both')
            cbar.set_label(plabel.replace('_',' '), size = labelsize, weight = 'bold')
    #ax.add_feature(cfeature.STATES)

    # Zoom in the figure if desired
    # -----------------------------
    if(zoom):
        ax.set_extent([lon_lims[0],lon_lims[1],lat_lims[0],lat_lims[1]],\
                       crs = ccrs.PlateCarree())
        zoom_add = '_zoom'
    else:
        zoom_add = ''

    #if(counties):
    #    ax.add_feature(USCOUNTIES.with_scale('5m'), alpha = 0.5)    

    # NOTE: commented out after removing the 'enhanced_image' code because
    #       it doesn't work now .
    ##!#ax.coastlines(resolution = '50m')
    ##!#ax.add_feature(cfeature.STATES)
    ##!#ax.add_feature(cfeature.BORDERS)
    if(ptitle is None):
        if(channel == 'true_color'):
            ax.set_title('VIIRS '+ dt_date_str.strftime('%Y-%m-%d %H:%M'))
        else:
            ax.set_title('VIIRS Ch ' + str(channel) + ' (' + \
                str(np.round(np.mean(channel_dict[str(channel)]['Bandwidth']),\
                 3)) + ' μm )\n'+\
                dt_date_str.strftime('%Y-%m-%d %H:%M'))
    ##!#axcs.set_xlabel('Ch. ' + str(MODIS_data0['channel']) +' [' + \
    ##!#    channel_dict[str(channel0)]['Bandwidth_label'] + \
    ##!#    MODIS_data0['variable'])
    else:
        ax.set_title(ptitle)

    if(not in_ax): 
        print('here')
        if(save):
            outname = 'viirs_satpy_' + date_str + zoom_add + '.png'
            plt.savefig(outname,dpi=300)
            print("Saved image",outname)
        else:
            plt.show()





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

# city: 'TriCities', 
def plot_VIIRS_twopanel(date_str1, date_str2, band, satellite = 'SNPP', \
        check_download = False, vmin = None, vmax = None, zoom = False, \
        city = 'TriCities', save = False):

    file_adder = sat_changer[satellite]

    mapcrs = init_proj('202107232155')
    plt.close('all')
    fig = plt.figure(figsize = (11,5))
    ax1 = fig.add_subplot(1,2,1, projection = mapcrs)
    ax2 = fig.add_subplot(1,2,2, projection = mapcrs)
    
   
    if(check_download): 
        identify_VIIRS_file(date_str1, band, viirs_dir, satellite = 'SNPP')
        identify_VIIRS_file(date_str2, band, viirs_dir, satellite = 'SNPP')
   
    # Grab the filenames
    # ------------------ 
    dt_date_str1 = datetime.strptime(date_str1, '%Y%m%d%H%M')
    dt_date_str2 = datetime.strptime(date_str2, '%Y%m%d%H%M')
    filename1 = glob(dt_date_str1.strftime(viirs_dir + satellite + \
        '/*02DNB.A%Y%j.%H%M*.nc'))[0]
    filename2 = glob(dt_date_str2.strftime(viirs_dir + satellite + \
        '/*02DNB.A%Y%j.%H%M*.nc'))[0]
    print(filename1)
    print(filename2)
  
    # --------------------------------------------------------------
    # 
    # Plot data for the first file
    #
    # --------------------------------------------------------------
    
    # Read the data for this granule
    # ------------------------------
    VIIRS_data = readVIIRS_granule(filename1,\
        band = band, zoom = zoom)
    
    # Plot the data for this granule
    # ------------------------------
    plot_VIIRS_granule(VIIRS_data, ax = ax1, zoom = zoom, \
            vmin = vmin, vmax = vmax, show_smoke = False)
    
    ax1.set_title(dt_date_str1.strftime('%Y-%m-%d %H:%M UTC\nSmoke-free'))
    
    ax1.coastlines()
    ax1.add_feature(cfeature.BORDERS)
    ax1.add_feature(cfeature.STATES)
    
    
    # --------------------------------------------------------------
    # 
    # Plot data for the second file
    #
    # --------------------------------------------------------------
    
    # Read the data for this granule
    # ------------------------------
    VIIRS_data = readVIIRS_granule(filename2,\
        band = band, zoom = zoom)
    
    # Plot the data for this granule
    # ------------------------------
    plot_VIIRS_granule(VIIRS_data, ax = ax2, zoom = zoom, \
            vmin = vmin, vmax = vmax, show_smoke = False)
    
    ax2.set_title(dt_date_str2.strftime('%Y-%m-%d %H:%M UTC\nSmoky'))
    
    ax2.coastlines()
    ax2.add_feature(cfeature.BORDERS)
    ax2.add_feature(cfeature.STATES)

    if(city == 'TriCities'):   
        ax1.set_extent([-119.5, -118.9, 46., 46.5], datacrs)
        ax2.set_extent([-119.5, -118.9, 46., 46.5], datacrs)
    elif(city == 'Yakima'):   
        ax1.set_extent([-120.9, -120.1, 46.3, 46.9], datacrs)
        ax2.set_extent([-120.9, -120.1, 46.3, 46.9], datacrs)
    else:
        ax1.set_extent([-125., -113., 41., 49.5], datacrs)
        ax2.set_extent([-125., -113., 41., 49.5], datacrs)
 
    title_str =  'Suomi-NPP VIIRS DNB'
    if(city == 'TriCities'):
        title_str = title_str + '\nTri-Cities, WA'
    if(city == 'Yakima'):
        title_str = title_str + '\nYakima, WA'
    plt.suptitle(title_str)
    plot_subplot_label(ax1, '(a)', backgroundcolor = 'white')
    plot_subplot_label(ax2, '(b)', backgroundcolor = 'white')
 
    fig.tight_layout()
    

    if(save):
        if(city is not None):
            outname = 'viirs_twopanel_' + satellite + '_' + date_str1 + \
                '_' + date_str2 + '_' + city + '.png'
        else:
            outname = 'viirs_twopanel_' + satellite + '_' + date_str1 + \
                '_' + date_str2 + '.png'
        plt.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else: 
        plt.show()
    

def plot_VIIRS_sixpanel(satellite = 'SNPP', save = False):

    file_adder = sat_changer[satellite]

    if(satellite == 'SNPP'):
        date_str1 = '202107222124'
        date_str2 = '202107230942'
        date_str3 = '202107232100'
    elif(satellite == 'JPSS'):
        date_str1 = '202107222030'
        date_str2 = '202107231030'
        date_str3 = '202107232012'

    dt_date_str1 = datetime.strptime(date_str1, '%Y%m%d%H%M')
    dt_date_str2 = datetime.strptime(date_str2, '%Y%m%d%H%M')
    dt_date_str3 = datetime.strptime(date_str3, '%Y%m%d%H%M')

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
    #filename_n0722_dnb = glob('home_dir + /data/VIIRS/DNB/VNP02DNB.A2021203.1000*.nc')
    #filename_n0722_dat = glob('home_dir + /data/VIIRS/DNB/VNP02MOD.A2021203.1000*.nc')
    #filename_n0722_dnb = glob('home_dir + /data/VIIRS/DNB/VNP02DNB.A2021203.1000*.nc')
    #filename_n0722_dat = glob('home_dir + /data/VIIRS/DNB/VNP02MOD.A2021203.1000*.nc')
    filename_d0722_dat = glob(dt_date_str1.strftime(\
        home_dir + '/data/VIIRS/DNB/'+satellite+'/V'+file_adder + \
        '2MOD.A%Y%j.%H%M*.nc'))
    filename_n0723_dnb = glob(dt_date_str2.strftime(\
        home_dir + '/data/VIIRS/DNB/'+satellite+'/V'+file_adder + \
        '2DNB.A%Y%j.%H%M*.nc'))
    filename_n0723_dat = glob(dt_date_str2.strftime(\
        home_dir + '/data/VIIRS/DNB/'+satellite+'/V' + file_adder + \
        '2MOD.A%Y%j.%H%M*.nc'))
    filename_d0723_dat = glob(dt_date_str3.strftime(\
        home_dir + '/data/VIIRS/DNB/'+satellite+'/V'+file_adder + \
        '2MOD.A%Y%j.%H%M*.nc'))

    print(filename_d0722_dat)
    print(filename_n0723_dnb)
    print(filename_n0723_dat)
    print(filename_d0723_dat)

    ##!#filename_d0722_dat = glob('home_dir + /data/VIIRS/DNB/V'+file_adder + \
    ##!#    '2MOD.A2021203.2124*.nc')
    ##!#filename_n0723_dnb = glob('home_dir + /data/VIIRS/DNB/V'+file_adder + \
    ##!#    '2DNB.A2021204.0942*.nc')
    ##!#filename_n0723_dat = glob('home_dir + /data/VIIRS/DNB/V'+file_adder + \
    ##!#    '2MOD.A2021204.0942*.nc')
    ##!#filename_d0723_dat = glob('home_dir + /data/VIIRS/DNB/V'+file_adder + \
    ##!#    '2MOD.A2021204.2100*.nc')

    namesplit = filename_d0722_dat[0].split('/')[-1].split('.')
    date1 = datetime.strptime('.'.join(namesplit[1:3])[1:], '%Y%j.%H%M')
    namesplit = filename_n0723_dat[0].split('/')[-1].split('.')
    date2 = datetime.strptime('.'.join(namesplit[1:3])[1:], '%Y%j.%H%M')
    namesplit = filename_d0723_dat[0].split('/')[-1].split('.')
    date3 = datetime.strptime('.'.join(namesplit[1:3])[1:], '%Y%j.%H%M')

    print(date1.strftime('%Y-%m-%d %H:%M'))
    print(date2.strftime('%Y-%m-%d %H:%M'))
    print(date3.strftime('%Y-%m-%d %H:%M'))

    plot_VIIRS_figure(filename_d0722_dat, band = 'M05', ax = ax1, \
        ptitle = date1.strftime('%d %B %Y %H:%M UTC'), vmax = 0.7, zoom = True)
    plot_VIIRS_figure(filename_n0723_dnb, band = 'DNB', ax = ax2, \
        ptitle = date2.strftime('%d %B %Y %H:%M UTC'), zoom = True)
    plot_VIIRS_figure(filename_d0723_dat, band = 'M05', ax = ax3,\
        ptitle = date3.strftime('%d %B %Y %H:%M UTC'), vmax = 0.7, zoom = True)

    plot_VIIRS_figure(filename_d0722_dat, band = 'M15', ax = ax4, \
        zoom = True, vmin = 260, vmax = 330)
    plot_VIIRS_figure(filename_n0723_dat, band = 'M15', ax = ax5, \
        zoom = True, vmin = 270, vmax = 290)
    plot_VIIRS_figure(filename_d0723_dat, band = 'M15', ax = ax6, \
        zoom = True, vmin = 260, vmax = 330)

    plot_subplot_label(ax1, '(a)', backgroundcolor = 'white')
    plot_subplot_label(ax2, '(b)', backgroundcolor = 'white')
    plot_subplot_label(ax3, '(c)', backgroundcolor = 'white')
    plot_subplot_label(ax4, '(d)', backgroundcolor = 'white')
    plot_subplot_label(ax5, '(e)', backgroundcolor = 'white')
    plot_subplot_label(ax6, '(f)', backgroundcolor = 'white')

    plot_figure_text(ax1, 'VIIRS 0.67 μm', xval = None, yval = None, \
        transform = None, color = 'red', fontsize = 12, \
        backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax2, 'VIIRS DNB', xval = None, yval = None, \
        transform = None, color = 'red', fontsize = 12, \
        backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax3, 'VIIRS 0.67 μm', xval = None, yval = None, \
        transform = None, color = 'red', fontsize = 12, \
        backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax4, 'VIIRS 10.76 μm', xval = None, yval = None, \
        transform = None, color = 'red', fontsize = 12, \
        backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax5, 'VIIRS 10.76 μm', xval = None, yval = None, \
        transform = None, color = 'red', fontsize = 12, \
        backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax6, 'VIIRS 10.76 μm', xval = None, yval = None, \
        transform = None, color = 'red', fontsize = 12, \
        backgroundcolor = 'white', halign = 'right')

    fig.tight_layout()

    if(save):
        outname = 'viirs_sixpanel_' + satellite + '_' + date_str + '_v2.png'
        plt.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else: 
        plt.show()
    
def plot_VIIRS_ninepanel(save = False):

    file_adder = sat_changer[satellite]

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
    filename_n0722_dnb = glob(home_dir + '/data/VIIRS/DNB/V' + \
        file_adder + '2DNB.A2021203.1000*.nc')
    filename_n0722_dat = glob(home_dir + '/data/VIIRS/DNB/V' + \
        file_adder + '2MOD.A2021203.1000*.nc')
    filename_d0722_dat = glob(home_dir + '/data/VIIRS/DNB/V' + \
        file_adder + '2MOD.A2021203.2124*.nc')
    filename_n0723_dnb = glob(home_dir + '/data/VIIRS/DNB/V' + \
        file_adder + '2DNB.A2021204.0942*.nc')
    filename_n0723_dat = glob(home_dir + '/data/VIIRS/DNB/V' + \
        file_adder + '2MOD.A2021204.0942*.nc')

    ##!#plot_VIIRS_figure(filename_n0722_dnb, band = 'DNB', ax = ax1, zoom = True)
    plot_VIIRS_figure(filename_d0722_dat, band = 'M05', ax = ax2, \
        show_smoke = True, zoom = True)
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
    
