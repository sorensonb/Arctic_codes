"""
  NAME:

  PURPOSE:
  
    NOTE: The AIRS channel and true color functions are designed to work with
    HDF AIRS files retriefed from 
    the data ordering website at this address:
    https://ladsweb.modaps.eosdis.nasa.gov/search/order/1/AIRS:Aqua


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
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as cm
from matplotlib.dates import DateFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.tri import Triangulation
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from satpy.scene import Scene
from satpy.writers import get_enhanced_image
from glob import glob

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# Set up global variables
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
mapcrs = ccrs.NorthPolarStereo(central_longitude = 45.)
datacrs = ccrs.PlateCarree()


plot_limits_dict = {
    "2021-07-20": {
        '2125': {
            'asos': 'asos_data_20210720.csv',
            'modis': '/home/bsorenson/data/AIRS/Aqua/MYD021KM.A2021201.2125.061.2021202154814.hdf',
            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072010-2021072021.nc',
            'Lat': [39.5, 42.0],
            'Lon': [-122.0, -119.5]
        }
    },
    "2021-07-21": {
        '2030': {
            'asos': 'asos_data_20210722.csv',
            'modis': '/home/bsorenson/data/AIRS/Aqua/MYD021KM.A2021202.2030.061.2021203174050.hdf',
            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072109-2021072122.nc',
            'Lat': [39.5, 42.0],
            'Lon': [-122.0, -119.5],
            'modis_Lat': [39.0, 42.5],
            'modis_Lon': [-123., -119.]
        }
    },
    "2021-07-22": {
        '2110': {
            'asos': 'asos_data_20210722.csv',
            'modis': '/home/bsorenson/data/AIRS/Aqua/MYD021KM.A2021203.2110.061.2021204155922.hdf',
            'mdswv': '/home/bsorenson/data/AIRS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072210-2021072221.nc',
            'Lat': [39.5, 42.0],
            'Lon': [-122.0, -119.5],
            'modis_Lat': [39.0, 42.5],
            'modis_Lon': [-123., -119.]
        }
    },
    "2021-08-05": {
        '2120': {
            'asos': 'asos_data_20210805.csv',
            'modis': '/home/bsorenson/data/AIRS/Aqua/MYD021KM.A2021217.2120.061.2021218164201.hdf',
            'Lat': [36.0, 39.0],
            'Lon': [-118.0, -114.0]
        },
        '2125': {
            'asos': 'asos_california_20210805.csv',
            'modis': '/home/bsorenson/data/AIRS/Aqua/MYD021KM.A2021217.2125.061.2021218161010.hdf',
            'Lat': [39.5, 42.5],
            'Lon': [-121.5, -119.5],
            'modis_Lat': [39.0, 42.5],
            'modis_Lon': [-123., -119.]
        }
    },
    "2021-08-06": {
        '2025': {
            'asos': 'asos_nevada_20210806.csv',
            'modis': '/home/bsorenson/data/AIRS/Aqua/MYD021KM.A2021218.2025.061.2021219151802.hdf',
            'mdswv': '/home/bsorenson/data/AIRS/Aqua/MYD05_L2.A2021218.2025.061.2021219152751.hdf',
            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
            'Lat': [36.0, 39.0],
            'Lon': [-118.0, -114.0],
            'modis_Lat': [35.0, 40.0],
            'modis_Lon': [-119., -113.]
        }
    },
    "2021-08-17": {
        '2145': {
            'asos': 'asos_data_20210817.csv',
            'Lat': [38.0, 42.0],
            'Lon': [-122.0, -117.0]
        }
    },
    "2021-08-30": {
        '2115': {
            'asos': 'asos_data_20210830.csv',
            'modis': '/home/bsorenson/data/AIRS/Aqua/MYD021KM.A2021242.2115.061.2021243183953.hdf',
            'Lat': [38.0, 40.0],
            'Lon': [-121.0, -118.5]
        }
    } 
}

# Find the gridpoint in the gridded lat/lon data that 
# corresponds to the station at slat and slon
# ---------------------------------------------------- 
def nearest_gridpoint(slat, slon, grid_lat, grid_lon):
    fun_c = np.maximum(np.abs(grid_lat - slat), \
        np.abs(grid_lon - slon))
    m_idx = np.where(fun_c == np.min(fun_c))
    return m_idx
 
# Extract the AIRS information from a given channel at each ob point
# -------------------------------------------------------------------
def nearest_grid_values(AIRS_data):
    # Read in the correct ASOS file 
    asos_file = plot_limits_dict[AIRS_data['cross_date']][AIRS_data['file_time']]['asos']
    df = pd.read_csv(asos_file)
    df['valid'] = pd.to_datetime(df['valid'])
    df = df.set_index('valid')

    # Pull the event time from the plot_limits_dict
    event_date = datetime.strptime(AIRS_data['cross_date'], "%Y-%m-%d")
    first_time = AIRS_data['file_time']
    event_dtime = event_date + timedelta(hours = int(first_time[:2]))
    # Test what happens if looking at the data near the hour
    #event_dtime = event_date + timedelta(hours = int(first_time[:2]), \
    #    minutes = int(first_time[2:4]))
    print("Overpass time",event_dtime)

    begin_range = event_dtime - timedelta(minutes = 15)
    end_range   = event_dtime + timedelta(minutes = 15)

    compare_dict = {}
    station_names = sorted(set(df['station'].values))
    compare_dict['modis_time'] = event_dtime.strftime('%Y%m%d%H%M')
    compare_dict['stations'] = station_names
    compare_dict['stn_data'] = np.zeros((len(station_names)))
    compare_dict['mds_data'] = np.zeros((len(station_names)))
 
    for ii, station in enumerate(station_names):
        # Get the correct ob data
        stn_df = df[df['station'] == station]
        lat_stn = stn_df['lat'].values[0]
        lon_stn = stn_df['lon'].values[0]

        #print(stn_df.index)
        ##if(station == 'SVE'):
        ##    for val in stn_df.index:
        ##        print(val) 
        #print(stn_df[event_dtime - timedelta(minutes = 30), event_dtime + timedelta(minutes = 30)])
        stn_df = stn_df[ begin_range : end_range]
        if(len(stn_df.index) == 0):
            # No data for this time.
            print("No valid data for station",station)
            compare_dict['stn_data'][ii] = np.nan 
        else:
            s_idx = np.argmin(np.abs((stn_df.index - event_dtime).total_seconds())) 
            if(stn_df['tmpc'][s_idx] == 'M'):
                stn_tmps = pd.to_numeric(stn_df['tmpc'], errors='coerce').values
                s_idx = np.where(~np.isnan(stn_tmps))[0][0]
            compare_dict['stn_data'][ii] = stn_df['tmpc'][s_idx]

        # Find the matching grid index
        m_idx = nearest_gridpoint(lat_stn, lon_stn, AIRS_data['lat'], \
            AIRS_data['lon'])
        m_data = AIRS_data['data'][m_idx][0]
      
        compare_dict['mds_data'][ii] = m_data
 
        ##print(station, lat_stn, AIRS_data['lat'][m_idx][0], \
        ##    lon_stn, AIRS_data['lon'][m_idx][0], compare_dict['stn_data'][ii], m_data )

    return compare_dict

 
# Plot the downloaded ASOS stations for each case
# ----------------------------------------------- 
def plot_ASOS_locs(pax,AIRS_data,color='red'):
    # Read in the correct ASOS file 
    asos_file = plot_limits_dict[AIRS_data['cross_date']][AIRS_data['file_time']]['asos']
    df = pd.read_csv(asos_file)

    station_names = set(df['station'].values)
    for ii, station in enumerate(station_names):
        lat_stn = df['lat'][df['station'] == station].values[0]
        lon_stn = df['lon'][df['station'] == station].values[0]
        pax.text(lon_stn, lat_stn, station, fontsize=10,weight='bold', \
            transform=datacrs, color=color)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Plotting functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def plot_AIRS_tlevel(filename,level,zoom=True,save=False):

    # Read in the AIRS file
    # ---------------------
    data = Dataset(filename,'r')

    # Extract orbit information
    # -------------------------
    idx = getattr(data,'coremetadata').split().index('EQUATORCROSSINGDATE')
    cross_date = getattr(data,'coremetadata').split()[idx+9]
    idx = getattr(data,'coremetadata').split().index('EQUATORCROSSINGTIME')
    cross_time = getattr(data,'coremetadata').split()[idx+9]

    # Extract lat, lon, and temperature data for level
    # ------------------------------------------------
    lon = data.variables['lonAIRS'][:,:,1,1]
    lat = data.variables['latAIRS'][:,:,1,1]

    tmp = data.variables['TAirStd'][:,:,level]

    data.close()

    # Plot the data
    # -------------
    plt.close()
    ax = plt.axes(projection = ccrs.LambertConformal())

    ax.coastlines()
    mesh = ax.pcolormesh(lon, lat, tmp, transform = datacrs, \
        cmap = 'plasma', shading = 'auto')
    cbar = plt.colorbar(mesh, ax = ax)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.STATES)
    ax.set_extent([np.nanmin(lon), np.nanmax(lon),\
        np.nanmin(lat), np.nanmax(lat)], datacrs)
    ax.set_title(cross_date + ' ' + cross_time + '\n' + str(level))
    plt.show() 

    """
    # Determine the correct AIRS file associated with the date
    # ---------------------------------------------------------
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    filename = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis']
    print(filename)
    if(composite):
        day_filenames = glob(filename[:50]+'*')
        cmpst_add = '_composite'
    else:
        day_filenames = glob(filename)
        cmpst_add = ''

    # Extract the modis true-color plot limits
    # ----------------------------------------
    lat_lims = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis_Lat']
    lon_lims = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis_Lon']

    # Use satpy (Scene) to open the file
    # ----------------------------------
    scn = Scene(reader = 'modis_l1b', filenames = day_filenames)

    # Load true-color data 
    scn.load(['true_color'])

    # Set the map projection and center the data
    # ------------------------------------------
    my_area = scn['true_color'].attrs['area'].compute_optimal_bb_area({\
        'proj':'lcc', 'lon_0': lon_lims[0], 'lat_0': lat_lims[0], 'lat_1': lat_lims[0], 'lat_2': lat_lims[0]})
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
    ax.imshow(var.data, transform = crs, extent=(var.x[0], var.x[-1], var.y[-1], var.y[0]), origin='upper')

    # Zoom in the figure if desired
    # -----------------------------
    if(zoom):
        ax.set_extent([lon_lims[0],lon_lims[1],lat_lims[0],lat_lims[1]],\
                       crs = ccrs.PlateCarree())
        zoom_add = '_zoom'
    else:
        zoom_add = ''

    ax.set_title('Aqua AIRS\n'+dt_date_str.strftime('%Y-%m-%d %H:%M'))

    if(save):
        outname = 'modis_true_color_' + date_str + zoom_add + cmpst_add + '.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

    """

def plot_AIRS_tprof(filename,slat,slon,zoom=True,save=False):

    '/INVENTORYMETADATA/MEASUREDPARAMETER/ORBITCALCULATEDSPATIALDOMAIN/ORBITALCALCULATEDSPATIALDOMAINCONTAINER/EQUATORCROSSINGTIME'
 
    # Read in the AIRS file
    # ---------------------
    data = Dataset(filename,'r')

    # Extract lat, lon, and temperature data for level
    # ------------------------------------------------
    lon = data.variables['lonAIRS'][:,:,1,1]
    lat = data.variables['latAIRS'][:,:,1,1]

    # Determine the indices where the desired lat and lon are located
    # ---------------------------------------------------------------
    # Find the matching grid index
    a_idx = nearest_gridpoint(slat, slon, lat, lon)

    print(a_idx[0][0], a_idx[1][0])
    tmp = data.variables['TAirStd'][a_idx[0][0], a_idx[1][0],:]
    ght = data.variables['GP_Height'][a_idx[0][0], a_idx[1][0],:]
    print(lat[a_idx], lon[a_idx])

    for tx, gx in zip(tmp, ght):
        print(tx, gx)

    data.close()

    # Plot the data
    # -------------
    plt.close()
    fig, ax = plt.subplots()
    ax.plot(tmp, ght)
    ax.set_ylabel('Geopotential Height [m]')
    ax.set_xlabel('Temperature [K]') 
    ax.set_title(str(lat[a_idx]) + ' N , ' + str(lon[a_idx]) + ' E')
    plt.show() 


# dt_date_str is of format YYYYMMDDHHMM
def read_AIRS_channel(date_str, channel, zoom = False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # Extract the filename given the channel
    # --------------------------------------
    if(str(channel)[:2] == 'wv'):
        filename = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['mdswv']
    else:
        filename = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis']
    
    print("Reading AIRS channel",channel," from ",filename)

    AIRS_data = {}

    modis = SD.SD(filename)

    dat = modis.attributes().get('CoreMetadata.0').split()
    indx = dat.index('EQUATORCROSSINGDATE')+9
    cross_date = dat[indx][1:len(dat[indx])-1]

    print(cross_date)

    lat5 = modis.select('Latitude').get()
    lon5 = modis.select('Longitude').get()

    data  = modis.select(channel_dict[str(channel)]['Name']).get()
    if(str(channel)[2:] != '_ir'):
        if(str(channel) != 'wv_nir'):
            data = data[channel_dict[str(channel)]['Index']]
        data  = data[::5,::5]

        if(data.shape != lat5.shape):
            data = data[:lat5.shape[0], :lat5.shape[1]]


    if(str(channel)[:2] == 'wv'):
        data_scale    = modis.select(channel_dict[str(channel)]['Name']\
            ).attributes().get('scale_factor')
        data_offset   = modis.select(channel_dict[str(channel)]['Name']\
            ).attributes().get('add_offset')
        # Extract the fill value and make sure any missing values are
        # removed.
        # -----------------------------------------------------------
        mask_val = modis.select(channel_dict[str(channel)]['Name']\
            ).attributes().get('_FillValue')
        bad_locations = np.where(data == mask_val)

        # Calculate reflectance using the scales and offsets
        # -------------------------------------------------
        data = ((data - data_offset) * data_scale)

        # Convert any missing data to nans
        # --------------------------------
        data[bad_locations] = np.nan
        data = np.ma.masked_invalid(data)

        colors = 'plasma'
        try:
            label = modis.select(channel_dict[str(channel)]['Name']\
                ).attributes().get('long_name') +  ' [' + \
                modis.select(channel_dict[str(channel)]['Name']\
                ).attributes().get('units') + ']'
        except:
            label = modis.select(channel_dict[str(channel)]['Name']\
                ).attributes().get('long_name') +  ' [' + \
                modis.select(channel_dict[str(channel)]['Name']\
                ).attributes().get('unit') + ']'
    else:
        channel = int(channel)
        # Thermal emission data
        if((channel >= 20) & (channel != 26)):
            data_scale    = modis.select(channel_dict[str(channel)]['Name']\
                ).attributes().get('radiance_scales')[channel_dict[str(channel)]['Index']]
            data_offset   = modis.select(channel_dict[str(channel)]['Name']\
                ).attributes().get('radiance_offsets')[channel_dict[str(channel)]['Index']]

            # Extract the fill value and make sure any missing values are
            # removed.
            # -----------------------------------------------------------
            mask_val = modis.select(channel_dict[str(channel)]['Name']\
                ).attributes().get('_FillValue')
            bad_locations = np.where(data == mask_val)

            data = ((data - data_offset) * data_scale)

            # Define constants for converting radiances to temperatures
            lmbda = (1e-6) * (np.average(channel_dict[str(channel)]['Bandwidth'])) # in m
            print("Average wavelength = ",np.average(channel_dict[str(channel)]['Bandwidth']))
            c_const = 3e8
            h_const = 6.626e-34 # J*s
            k_const = 1.381e-23 # J/K

            data = (h_const * c_const) / \
                (lmbda * k_const * np.log( ((2.0 * h_const * (c_const**2.0) ) / \
                ((lmbda**4.) * (lmbda / 1e-6) * data ) ) + 1.0 ) )
            #data = ((h_const * c_const)/(k_const * lmbda)) * (np.log((2.0 * h_const * (c_const ** 2.0) / \
            #    ((lmbda**5.0) * data)) + 1) ** -1.)

            # Convert any missing data to nans
            # --------------------------------
            data[bad_locations] = np.nan
            data = np.ma.masked_invalid(data)

            colors = 'plasma'
            label = 'Blackbody Temperature [K]'

        # Reflectances
        else:
            data_scale    = modis.select(channel_dict[str(channel)]['Name']\
                ).attributes().get('reflectance_scales')[channel_dict[str(channel)]['Index']]
            data_offset   = modis.select(channel_dict[str(channel)]['Name']\
                ).attributes().get('reflectance_offsets')[channel_dict[str(channel)]['Index']]

            # Extract the fill value and make sure any missing values are
            # removed.
            # -----------------------------------------------------------
            mask_val = modis.select(channel_dict[str(channel)]['Name']\
                ).attributes().get('_FillValue')
            bad_locations = np.where(data == mask_val)

            # Calculate reflectance using the scales and offsets
            # -------------------------------------------------
            data = ((data - data_offset) * data_scale)

            # Convert any missing data to nans
            # --------------------------------
            data[bad_locations] = np.nan
            data = np.ma.masked_invalid(data)

            colors = 'Greys_r'
            label = 'Reflectance'
 
    modis.end()

    # Mask any fire pixels
    # --------------------
    data = np.ma.masked_where(data > 340, data)

    AIRS_data['data'] = data
    AIRS_data['lat']  = lat5
    AIRS_data['lon']  = lon5
    AIRS_data['variable']  = label
    AIRS_data['cross_date']  = cross_date
    AIRS_data['channel']  = channel
    AIRS_data['colors']  = colors
    AIRS_data['file_time'] = filename.strip().split('/')[-1].split('.')[2]

    if(zoom):
        # Mask AIRS_data['data'] that are outside the desired range
        # --------------------------------------------
        AIRS_data['data'][(((AIRS_data['lat'] < plot_limits_dict[AIRS_data['cross_date']][AIRS_data['file_time']]['Lat'][0]) | \
                             (AIRS_data['lat'] > plot_limits_dict[AIRS_data['cross_date']][AIRS_data['file_time']]['Lat'][1])) | \
                            ((AIRS_data['lon'] < plot_limits_dict[AIRS_data['cross_date']][AIRS_data['file_time']]['Lon'][0]) | \
                             (AIRS_data['lon'] > plot_limits_dict[AIRS_data['cross_date']][AIRS_data['file_time']]['Lon'][1])))] = -999.
        AIRS_data['lat'][ (((AIRS_data['lat'] < plot_limits_dict[AIRS_data['cross_date']][AIRS_data['file_time']]['Lat'][0]) | \
                             (AIRS_data['lat'] > plot_limits_dict[AIRS_data['cross_date']][AIRS_data['file_time']]['Lat'][1])) | \
                            ((AIRS_data['lon'] < plot_limits_dict[AIRS_data['cross_date']][AIRS_data['file_time']]['Lon'][0]) | \
                             (AIRS_data['lon'] > plot_limits_dict[AIRS_data['cross_date']][AIRS_data['file_time']]['Lon'][1])))] = -999.
        AIRS_data['lon'][ (((AIRS_data['lat'] < plot_limits_dict[AIRS_data['cross_date']][AIRS_data['file_time']]['Lat'][0]) | \
                             (AIRS_data['lat'] > plot_limits_dict[AIRS_data['cross_date']][AIRS_data['file_time']]['Lat'][1])) | \
                            ((AIRS_data['lon'] < plot_limits_dict[AIRS_data['cross_date']][AIRS_data['file_time']]['Lon'][0]) | \
                             (AIRS_data['lon'] > plot_limits_dict[AIRS_data['cross_date']][AIRS_data['file_time']]['Lon'][1])))] = -999.

        AIRS_data['data'] = np.ma.masked_where(AIRS_data['data'] == -999., AIRS_data['data'])
        AIRS_data['lat'] = np.ma.masked_where(AIRS_data['lat'] == -999., AIRS_data['lat'])
        AIRS_data['lon'] = np.ma.masked_where(AIRS_data['lon'] == -999., AIRS_data['lon'])


    return AIRS_data

def plot_AIRS_channel(date_str,channel,zoom=True,show_smoke=False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    filename = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis']

    if(channel == 'red'):
        channel = 1
    elif(channel == 'green'):
        channel = 4
    elif(channel == 'blue'):
        channel = 3

    # Call read_AIRS_channel to read the desired AIRS data from the
    # file and put it in a dictionary
    # ---------------------------------------------------------------
    AIRS_data = read_AIRS_channel(dt_date_str, channel)

    print("Data max = ",np.max(AIRS_data['data']), "  Data min = ",np.min(AIRS_data['data']))

    plt.close('all')
    fig1 = plt.figure()
    datacrs = ccrs.PlateCarree()
    mapcrs = ccrs.LambertConformal()
    ax = plt.axes(projection = mapcrs)

    mesh = ax.pcolormesh(AIRS_data['lon'],AIRS_data['lat'],\
        AIRS_data['data'],cmap = AIRS_data['colors'], shading='auto', \
        transform = datacrs) 

    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(filename) 
        hash0 = ax0.pcolor(AIRS_data['lon'],AIRS_data['lat'],\
            hash_data1, hatch = '///', alpha=0., transform = datacrs)

    cbar = plt.colorbar(mesh,orientation='horizontal',pad=0,\
        aspect=50,shrink = 0.850)
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label(AIRS_data['variable'],fontsize=16,weight='bold')
    
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.STATES)
    ax.coastlines()
    if(zoom):
        ax.set_extent([plot_limits_dict[AIRS_data['cross_date']][AIRS_data['file_time']]['Lon'][0], \
                       plot_limits_dict[AIRS_data['cross_date']][AIRS_data['file_time']]['Lon'][1], \
                       plot_limits_dict[AIRS_data['cross_date']][AIRS_data['file_time']]['Lat'][0], \
                       plot_limits_dict[AIRS_data['cross_date']][AIRS_data['file_time']]['Lat'][1]],\
                       ccrs.PlateCarree())
    ax.set_title('Channel ' + str(channel) + '\n' + \
        channel_dict[str(channel)]['Bandwidth_label']) 

    plt.show()

# date_str of format "YYYYMMDDHHMM"
def compare_AIRS_3panel(date_str,channel1,channel2,channel3,zoom=True,save=False,\
        plot_ASOS_loc = False, show_smoke = True, compare_OMI = False, \
        compare_CERES = False, return_AIRS = False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    filename = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis']

    if(channel1== 'red'):
        channel1= 1
    elif(channel1== 'green'):
        channel1= 4
    elif(channel1== 'blue'):
        channel1= 3

    # Step 1: Call read_AIRS_channel to read the desired AIRS data from the
    # file and put it in a dictionary
    # -----------------------------------------------------------------------
    AIRS_data1 = read_AIRS_channel(dt_date_str.strftime("%Y%m%d%H%M"), channel1, zoom = True)
    AIRS_data2 = read_AIRS_channel(dt_date_str.strftime("%Y%m%d%H%M"), channel2, zoom = True)
    AIRS_data3 = read_AIRS_channel(dt_date_str.strftime("%Y%m%d%H%M"), channel3, zoom = True)

    ##!#print("Data1 max = ",np.max(AIRS_data1['data']), "  Data1 min = ",np.min(AIRS_data1['data']))
    ##!#print("Data2 max = ",np.max(AIRS_data2['data']), "  Data2 min = ",np.min(AIRS_data2['data']))
    ##!#print("Data3 max = ",np.max(AIRS_data3['data']), "  Data3 min = ",np.min(AIRS_data3['data']))

    # Create copies to set the minimum and maximums for each plot
    # -----------------------------------------------------------
    cpy_1 = np.copy(AIRS_data1['data'])
    cpy_2 = np.copy(AIRS_data2['data'])
    cpy_3 = np.copy(AIRS_data3['data'])

    cpy_1 = np.ma.masked_where((((AIRS_data1['lat'] < plot_limits_dict[AIRS_data1['cross_date']][AIRS_data1['file_time']]['Lat'][0]) | \
                         (AIRS_data1['lat'] > plot_limits_dict[AIRS_data1['cross_date']][AIRS_data1['file_time']]['Lat'][1])) | \
                        ((AIRS_data1['lon'] < plot_limits_dict[AIRS_data1['cross_date']][AIRS_data1['file_time']]['Lon'][0]) | \
                         (AIRS_data1['lon'] > plot_limits_dict[AIRS_data1['cross_date']][AIRS_data1['file_time']]['Lon'][1]))), cpy_1)
    cpy_2 = np.ma.masked_where((((AIRS_data2['lat'] < plot_limits_dict[AIRS_data2['cross_date']][AIRS_data2['file_time']]['Lat'][0]) | \
                         (AIRS_data2['lat'] > plot_limits_dict[AIRS_data2['cross_date']][AIRS_data2['file_time']]['Lat'][1])) | \
                        ((AIRS_data2['lon'] < plot_limits_dict[AIRS_data2['cross_date']][AIRS_data2['file_time']]['Lon'][0]) | \
                         (AIRS_data2['lon'] > plot_limits_dict[AIRS_data2['cross_date']][AIRS_data2['file_time']]['Lon'][1]))), cpy_2)
    cpy_3 = np.ma.masked_where((((AIRS_data3['lat'] < plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][0]) | \
                         (AIRS_data3['lat'] > plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][1])) | \
                        ((AIRS_data3['lon'] < plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][0]) | \
                         (AIRS_data3['lon'] > plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][1]))), cpy_3)


    # Step 2: Set up figure to have 3 panels
    # --------------------------------------
    datacrs = ccrs.PlateCarree() 
    mapcrs = ccrs.LambertConformal()

    plt.close('all')
    if(compare_OMI and compare_CERES):
        fig = plt.figure(figsize=(12,9))
        ax0 = fig.add_subplot(2,3,1,projection = mapcrs)
        ax1 = fig.add_subplot(2,3,2,projection = mapcrs)
        ax2 = fig.add_subplot(2,3,3,projection = mapcrs)
        axo = fig.add_subplot(2,3,4,projection = mapcrs)
        axcs = fig.add_subplot(2,3,5,projection = mapcrs)
        axcl = fig.add_subplot(2,3,6,projection = mapcrs)
    elif((compare_CERES and not compare_OMI)):
        fig = plt.figure(figsize=(12,9))
        ax0 = fig.add_subplot(2,3,1,projection = mapcrs)
        ax1 = fig.add_subplot(2,3,2,projection = mapcrs)
        ax2 = fig.add_subplot(2,3,3,projection = mapcrs)
        axcs = fig.add_subplot(2,3,4,projection = mapcrs)
        axcl = fig.add_subplot(2,3,5,projection = mapcrs)
    elif((compare_OMI and not compare_CERES)):
        fig = plt.figure(figsize=(9,9))
        ax0 = fig.add_subplot(2,2,1,projection = mapcrs)
        ax1 = fig.add_subplot(2,2,2,projection = mapcrs)
        ax2 = fig.add_subplot(2,2,3,projection = mapcrs)
        axo = fig.add_subplot(2,2,4,projection = mapcrs)
    else:
        fig = plt.figure(figsize=(14.5,5))
        ax0 = fig.add_subplot(1,3,1,projection = mapcrs)
        ax1 = fig.add_subplot(1,3,2,projection = mapcrs)
        ax2 = fig.add_subplot(1,3,3,projection = mapcrs)

    # Step 3: Plot the AIRS channel data in the first 2 panels
    # ---------------------------------------------------------
    mesh0 = ax0.pcolormesh(AIRS_data1['lon'],AIRS_data1['lat'],\
        cpy_1,cmap = AIRS_data1['colors'], shading='auto', \
        vmin = np.nanmin(AIRS_data1['data']), \
        vmax = np.nanmax(AIRS_data1['data']), transform = datacrs) 

    cbar0 = plt.colorbar(mesh0,ax=ax0,orientation='vertical',\
        pad=0.03,label=AIRS_data1['variable'])

    ax0.add_feature(cfeature.BORDERS)
    ax0.add_feature(cfeature.STATES)
    ax0.coastlines()
    if(zoom):
        ax0.set_extent([plot_limits_dict[AIRS_data1['cross_date']][AIRS_data1['file_time']]['Lon'][0], \
                        plot_limits_dict[AIRS_data1['cross_date']][AIRS_data1['file_time']]['Lon'][1], \
                        plot_limits_dict[AIRS_data1['cross_date']][AIRS_data1['file_time']]['Lat'][0], \
                        plot_limits_dict[AIRS_data1['cross_date']][AIRS_data1['file_time']]['Lat'][1]],\
                        datacrs)
    #ax0.set_title('AIRS Ch. ' + str(channel1) + '\n' + \
    #    str(channel_dict[str(channel1)]['Bandwidth'][0]) + ' μm - ' + \
    #    str(channel_dict[str(channel1)]['Bandwidth'][1]) + ' μm')
    ax0.set_title('Channel ' + str(channel1) + '\n' + \
        channel_dict[str(channel1)]['Bandwidth_label']) 


    # Plot channel 2
    mesh1 = ax1.pcolormesh(AIRS_data2['lon'],AIRS_data2['lat'],\
        AIRS_data2['data'],cmap = AIRS_data2['colors'], shading='auto', \
        vmin = np.nanmin(AIRS_data2['data']), \
        vmax = 2., transform = datacrs) 


    cbar1 = plt.colorbar(mesh1,ax=ax1,orientation='vertical',\
        pad=0.03,label=AIRS_data2['variable'])
    
    ax1.add_feature(cfeature.BORDERS)
    ax1.add_feature(cfeature.STATES)
    ax1.coastlines()
    if(zoom):
        ax1.set_extent([plot_limits_dict[AIRS_data2['cross_date']][AIRS_data2['file_time']]['Lon'][0], \
                        plot_limits_dict[AIRS_data2['cross_date']][AIRS_data2['file_time']]['Lon'][1], \
                        plot_limits_dict[AIRS_data2['cross_date']][AIRS_data2['file_time']]['Lat'][0], \
                        plot_limits_dict[AIRS_data2['cross_date']][AIRS_data2['file_time']]['Lat'][1]],\
                        datacrs)
    #ax1.set_title('AIRS Ch. ' + str(channel2) + '\n' + \
    #    str(channel_dict[str(channel2)]['Bandwidth'][0]) + ' μm - ' + \
    #    str(channel_dict[str(channel2)]['Bandwidth'][1]) + ' μm')
    ax1.set_title('Channel ' + str(channel2) + '\n' + \
        channel_dict[str(channel2)]['Bandwidth_label']) 

    # Plot channel 3
    mesh2 = ax2.pcolormesh(AIRS_data3['lon'],AIRS_data3['lat'],\
        AIRS_data3['data'],cmap = AIRS_data3['colors'], shading='auto', \
        vmin = np.nanmin(AIRS_data3['data']), \
        vmax = np.nanmax(AIRS_data3['data']), transform = datacrs) 
    cbar2 = plt.colorbar(mesh2,ax=ax2,orientation='vertical',\
        pad=0.03,label=AIRS_data3['variable'])
    
    ax2.add_feature(cfeature.BORDERS)
    ax2.add_feature(cfeature.STATES)
    ax2.coastlines()
    if(zoom):
        ax2.set_extent([plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][0], \
                        plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][1], \
                        plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][0], \
                        plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][1]],\
                        datacrs)
    #ax2.set_title('AIRS Ch. ' + str(channel3) + '\n' + \
    #    str(channel_dict[str(channel3)]['Bandwidth'][0]) + ' μm - ' + \
    #    str(channel_dict[str(channel3)]['Bandwidth'][1]) + ' μm')
    ax2.set_title('Channel ' + str(channel3) + '\n' + \
        channel_dict[str(channel3)]['Bandwidth_label']) 


    if(compare_OMI):
        print("Reading OMI data")
        data = h5py.File(plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['omi'],'r')
        LAT   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude'][:,:]
        LON   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude'][:,:]
        UVAI  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex'][:,:]
        XTRACK = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags'][:,:]
        mask_UVAI = np.ma.masked_where((XTRACK < -2e5) | (UVAI < -2e5), UVAI)
        mask_UVAI = np.ma.masked_where((((LAT < plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][0]) | \
                             (LAT > plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][1])) | \
                            ((LON < plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][0]) | \
                             (LON > plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][1]))), mask_UVAI)

        mesh3 = axo.pcolormesh(LON,LAT,mask_UVAI, cmap = 'plasma', shading='auto', \
            vmin = np.nanmin(mask_UVAI), vmax = np.nanmax(mask_UVAI), transform = datacrs) 
        axo.add_feature(cfeature.BORDERS)
        axo.add_feature(cfeature.STATES)
        axo.coastlines()
        if(zoom):
            axo.set_extent([plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][0], \
                            plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][1], \
                            plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][0], \
                            plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][1]],\
                            datacrs)
        cbar3 = plt.colorbar(mesh3,ax=axo,orientation='vertical',\
            pad=0.03,label='OMI UVAI')
        axo.set_title('OMI UVAI')
        data.close()

    if(compare_CERES):
        print("Reading CERES data")

        # NOTE: need to use the time variable to screen out any data
        # that are not close to the event time.
        base_date = datetime(year=1970,month=1,day=1)
        start_date = dt_date_str - timedelta(hours = 1)
        end_date   = dt_date_str + timedelta(hours = 2)

        print(start_date, end_date)

        data = Dataset(plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['ceres'],'r')
        LAT   = 90. - data.variables['Colatitude_of_CERES_FOV_at_surface'][:]
        LON   = data.variables['Longitude_of_CERES_FOV_at_surface'][:]
        LON[LON>179.99] = -360.+LON[LON>179.99]
        swflux  = data.variables['CERES_SW_TOA_flux___upwards'][:]
        lwflux  = data.variables['CERES_LW_TOA_flux___upwards'][:]
        time  = data.variables['time'][:]
        local_time = np.array([base_date + relativedelta(days = ttime) for ttime in time])

        mask_LAT = LAT[ \
            (LAT >= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][0]) & \
            (LAT <= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][1]) & \
            (LON >= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][0]) & \
            (LON <= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][1]) & \
            (swflux > 0) & (lwflux > 0)]
        mask_LON = LON[ \
            (LAT >= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][0]) & \
            (LAT <= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][1]) & \
            (LON >= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][0]) & \
            (LON <= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][1]) & \
            (swflux > 0) & (lwflux > 0)]
        mask_swf = swflux[ \
            (LAT >= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][0]) & \
            (LAT <= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][1]) & \
            (LON >= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][0]) & \
            (LON <= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][1]) & \
            (swflux > 0) & (lwflux > 0)]
        mask_lwf = lwflux[ \
            (LAT >= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][0]) & \
            (LAT <= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][1]) & \
            (LON >= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][0]) & \
            (LON <= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][1]) & \
            (swflux > 0) & (lwflux > 0)]
        mask_time = local_time[ \
            (LAT >= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][0]) & \
            (LAT <= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][1]) & \
            (LON >= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][0]) & \
            (LON <= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][1]) & \
            (swflux > 0) & (lwflux > 0)]

        # Removed masked data
        triMesh = Triangulation(mask_LON,mask_LAT)

        ##!## Reshape the data to make it work with pcolormesh
        ##!#first_size = int(np.sqrt(mask_swf.shape)) 
        ##!#max_size = int(first_size ** 2.)
        ##!#mask_swf = mask_swf[:max_size].reshape((first_size, first_size))
        ##!#mask_lwf = mask_lwf[:max_size].reshape((first_size, first_size))
        ##!#LAT = LAT[:max_size].reshape((first_size, first_size))
        ##!#LON = LON[:max_size].reshape((first_size, first_size))
      
        print(np.nanmax(mask_LAT.compressed()), np.min(mask_LAT.compressed()))
        print(np.nanmax(LAT), np.nanmin(LAT))
        print(plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'])
 
        #scat3 = axcs.scatter(mask_LON, mask_LAT,mask_swf, transform = datacrs)
        mesh3 = axcs.scatter(mask_LON.compressed(), mask_LAT.compressed(),s = 120,marker='s',c = mask_swf.compressed(),cmap='plasma', transform = datacrs)
        #mesh3 = axcs.tricontourf(triMesh,mask_swf, cmap = 'plasma', shading='auto', \
        #    vmin = np.nanmin(mask_swf), vmax = np.nanmax(mask_swf), transform = datacrs) 
#        mesh3 = axcs.pcolormesh(LON,LAT,mask_swf, cmap = 'plasma', shading='auto', \
#            vmin = np.nanmin(mask_swf), vmax = np.nanmax(mask_swf), transform = datacrs) 
        axcs.add_feature(cfeature.BORDERS)
        axcs.add_feature(cfeature.STATES)
        axcs.coastlines()
        if(zoom):
            axcs.set_extent([plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][0], \
                             plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][1], \
                             plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][0], \
                             plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][1]],\
                             datacrs)
        cbar3 = plt.colorbar(mesh3,ax=axcs,orientation='vertical',\
            pad=0.03,label='TOA SWF [W/m2]')
        axcs.set_title('CERES SWF')

        mesh4 = axcl.scatter(mask_LON.compressed(), mask_LAT.compressed(),s = 120,marker = 's',c = mask_lwf.compressed(),cmap='plasma', transform = datacrs)
        #mesh4 = axcl.tricontourf(triMesh,mask_lwf, cmap = 'plasma', shading='auto', \
        #    vmin = np.nanmin(mask_lwf), vmax = np.nanmax(mask_lwf), transform = datacrs) 
        ##!#mesh4 = axcl.tricontourf(LON,LAT,mask_lwf, cmap = 'plasma', shading='auto', \
        ##!#    vmin = np.nanmin(mask_lwf), vmax = np.nanmax(mask_lwf), transform = datacrs) 
        axcl.add_feature(cfeature.BORDERS)
        axcl.add_feature(cfeature.STATES)
        axcl.coastlines()
        if(zoom):
            axcl.set_extent([plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][0], \
                            plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][1], \
                            plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][0], \
                            plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][1]],\
                            datacrs)
            cbar4 = plt.colorbar(mesh4,ax=axcl,orientation='vertical',\
                pad=0.03,label='TOA LWF [W/m2]')
        axcl.set_title('CERES LWF')

        data.close()

    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(dt_date_str.strftime('%Y%m%d%H%M')) 
        hash0 = ax0.pcolor(AIRS_data1['lon'],AIRS_data1['lat'],\
            hash_data1[:AIRS_data1['lat'].shape[0], :AIRS_data1['lat'].shape[1]], hatch = '///', alpha=0., transform = datacrs)
        hash1 = ax1.pcolor(AIRS_data2['lon'],AIRS_data2['lat'],\
            hash_data1[:AIRS_data2['lat'].shape[0], :AIRS_data2['lat'].shape[1]], hatch = '///', alpha=0., transform = datacrs)
        hash2 = ax2.pcolor(AIRS_data3['lon'],AIRS_data3['lat'],\
            hash_data1[:AIRS_data3['lat'].shape[0], :AIRS_data3['lat'].shape[1]], hatch = '///', alpha=0., transform = datacrs)
        if(compare_OMI): 
            hash3 = axo.pcolor(AIRS_data3['lon'],AIRS_data3['lat'],\
                hash_data1[:AIRS_data3['lat'].shape[0], :AIRS_data3['lat'].shape[1]], hatch = '///', alpha=0., transform = datacrs)
        if(compare_CERES): 
            hash4 = axcs.pcolor(AIRS_data3['lon'],AIRS_data3['lat'],\
                hash_data1[:AIRS_data3['lat'].shape[0], :AIRS_data3['lat'].shape[1]], hatch = '///', alpha=0., transform = datacrs)
            hash5 = axcl.pcolor(AIRS_data3['lon'],AIRS_data3['lat'],\
                hash_data1[:AIRS_data3['lat'].shape[0], :AIRS_data3['lat'].shape[1]], hatch = '///', alpha=0., transform = datacrs)
    

    if(plot_ASOS_loc):
        print("Plotting ASOS site")
        plot_ASOS_locs(ax0,AIRS_data1,color='lime')
        plot_ASOS_locs(ax1,AIRS_data1,color='lime')
        plot_ASOS_locs(ax2,AIRS_data1,color='lime')
        if(compare_OMI): plot_ASOS_locs(axo,AIRS_data1,color='black')
        if(compare_CERES): 
            plot_ASOS_locs(axcs,AIRS_data1,color='black')
            plot_ASOS_locs(axcl,AIRS_data1,color='black')

    cross_date = AIRS_data1['cross_date']
    file_time  = AIRS_data1['file_time']
    plt.suptitle(cross_date + ' ' + file_time)
    if(save):
        pdate = cross_date[:4] + cross_date[5:7] + cross_date[8:10] + file_time
        outname = 'modis_compare_ch' + str(channel1) + '_vs_ch' + \
            str(channel2) + '_ch' + str(channel3) + '_' + pdate + '_3panel.png'
        if(compare_OMI):
            outname = 'modis_compare_ch' + str(channel1) + '_vs_ch' + \
                str(channel2) + '_ch' + str(channel3) + '_omi_' + pdate + '_3panel.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else: 
        plt.show()

    if(return_AIRS):
        return AIRS_data1
    

def compare_AIRS_channels(date_str,channel1,channel2,zoom=True,save=False,\
        plot_ASOS_loc = False,show_smoke = True):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    filename = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis']

    if(channel1 == 'red'):
        channel1 = 1
    elif(channel1 == 'green'):
        channel1 = 4
    elif(channel1 == 'blue'):
        channel1 = 3

    # -----------------------------------------------------------------------
    # Step 1: Call read_AIRS_channel to read the desired AIRS data from the
    # file and put it in a dictionary
    # -----------------------------------------------------------------------
    AIRS_data1 = read_AIRS_channel(dt_date_str.strftime('%Y%m%d%H%M'), channel1, zoom = zoom)
    AIRS_data2 = read_AIRS_channel(dt_date_str.strftime('%Y%m%d%H%M'), channel2, zoom = zoom)

    print("Data1 max = ",np.max(AIRS_data1['data']), "  Data1 min = ",np.min(AIRS_data1['data']))
    print("Data2 max = ",np.max(AIRS_data2['data']), "  Data2 min = ",np.min(AIRS_data2['data']))

    print(AIRS_data1['data'].shape, AIRS_data2['data'].shape)

    # --------------------------------------
    # Step 2: Set up figure to have 3 panels
    # --------------------------------------
    datacrs = ccrs.PlateCarree() 
    mapcrs = ccrs.LambertConformal()

    plt.close('all')
    fig = plt.figure(figsize=(14.5,5))
    ax0 = fig.add_subplot(1,3,2,projection = mapcrs)
    ax1 = fig.add_subplot(1,3,3,projection = mapcrs)
    ax2 = fig.add_subplot(1,3,1)

    # Step 3: Plot the AIRS channel data in the first 2 panels
    # ---------------------------------------------------------
    mesh0 = ax0.pcolormesh(AIRS_data1['lon'],AIRS_data1['lat'],\
        AIRS_data1['data'],cmap = AIRS_data1['colors'], shading='auto', \
        transform = datacrs)

    cbar0 = plt.colorbar(mesh0,ax=ax0,orientation='vertical',\
        pad=0.03,label=AIRS_data1['variable'])
    #cbar0 = plt.colorbar(mesh0,cax = ax0, orientation='horizontal',pad=0,\
    #    aspect=50,shrink = 0.850)
    #cbar0.ax.tick_params(labelsize=14)
    #cbar0.set_label(AIRS_data1['variable'],fontsize=16,weight='bold')
   
    ax0.add_feature(cfeature.BORDERS)
    ax0.add_feature(cfeature.STATES)
    ax0.coastlines()
    if(zoom):
        ax0.set_extent([plot_limits_dict[AIRS_data1['cross_date']][AIRS_data1['file_time']]['Lon'][0], \
                        plot_limits_dict[AIRS_data1['cross_date']][AIRS_data1['file_time']]['Lon'][1], \
                        plot_limits_dict[AIRS_data1['cross_date']][AIRS_data1['file_time']]['Lat'][0], \
                        plot_limits_dict[AIRS_data1['cross_date']][AIRS_data1['file_time']]['Lat'][1]],\
                        datacrs)
    #ax0.set_title('AIRS Ch. ' + str(channel1) + '\n' + \
    #    str(channel_dict[str(channel1)]['Bandwidth'][0]) + ' μm - ' + \
    #    str(channel_dict[str(channel1)]['Bandwidth'][1]) + ' μm')
    ax0.set_title('Channel ' + str(channel1) + '\n' + \
        channel_dict[str(channel1)]['Bandwidth_label']) 

    # Plot channel 2
    mesh1 = ax1.pcolormesh(AIRS_data2['lon'],AIRS_data2['lat'],\
        AIRS_data2['data'],cmap = AIRS_data2['colors'], shading='auto', \
        transform = datacrs) 

    cbar1 = plt.colorbar(mesh1,ax=ax1,orientation='vertical',\
        pad=0.03,label=AIRS_data2['variable'])
    #cbar1 = plt.colorbar(mesh1,cax = ax1,orientation='horizontal',pad=1,\
    #    aspect=51,shrink = 1.851)
    #cbar1.ax.tick_params(labelsize=14)
    #cbar1.set_label(AIRS_data2['variable'],fontsize=16,weight='bold')
    
    ax1.add_feature(cfeature.BORDERS)
    ax1.add_feature(cfeature.STATES)
    ax1.coastlines()
    if(zoom):
        ax1.set_extent([plot_limits_dict[AIRS_data2['cross_date']][AIRS_data2['file_time']]['Lon'][0], \
                        plot_limits_dict[AIRS_data2['cross_date']][AIRS_data2['file_time']]['Lon'][1], \
                        plot_limits_dict[AIRS_data2['cross_date']][AIRS_data2['file_time']]['Lat'][0], \
                        plot_limits_dict[AIRS_data2['cross_date']][AIRS_data2['file_time']]['Lat'][1]],\
                        datacrs)
    #ax1.set_title('AIRS Ch. ' + str(channel2) + '\n' + \
    #    str(channel_dict[str(channel2)]['Bandwidth'][0]) + ' μm - ' + \
    #    str(channel_dict[str(channel2)]['Bandwidth'][1]) + ' μm')
    ax1.set_title('Channel ' + str(channel2) + '\n' + \
        channel_dict[str(channel2)]['Bandwidth_label']) 

    ##!## Shade areas that are inside the plume, defined currently as areas
    ##!## with reflectance below the mean
    ##!## -----------------------------------------------------------------
    ##!#if(AIRS_data1['variable'] == 'Reflectance'):
    ##!#    hash_data1 = np.ma.masked_where(AIRS_data1['data'] < \
    ##!#        np.nanmean(AIRS_data1['data']),AIRS_data1['data'])
    ##!#    hash_data2 = np.ma.masked_where(AIRS_data1['data'] < \
    ##!#        np.nanmean(AIRS_data1['data']),AIRS_data2['data'])

    ##!#    nohash_data1 = np.ma.masked_where(AIRS_data1['data'] >= \
    ##!#        np.nanmean(AIRS_data1['data']),AIRS_data1['data'])
    ##!#    nohash_data2 = np.ma.masked_where(AIRS_data1['data'] >= \
    ##!#        np.nanmean(AIRS_data1['data']),AIRS_data2['data'])

    ##!#    hash1 = ax0.pcolor(AIRS_data1['lon'],AIRS_data1['lat'],\
    ##!#        hash_data1, hatch = '///', alpha=0., transform = datacrs),
    ##!#    hash2 = ax1.pcolor(AIRS_data2['lon'],AIRS_data2['lat'],\
    ##!#        hash_data2, hatch = '///', alpha=0., transform = datacrs),
    ##!#if(AIRS_data2['variable'] == 'Reflectance'):
    ##!#    hash_data2 = np.ma.masked_where(AIRS_data2['data'] < \
    ##!#        np.nanmean(AIRS_data2['data']),AIRS_data2['data'])
    ##!#    hash_data1 = np.ma.masked_where(AIRS_data2['data'] < \
    ##!#        np.nanmean(AIRS_data2['data']),AIRS_data1['data'])

    ##!#    nohash_data2 = np.ma.masked_where(AIRS_data2['data'] >= \
    ##!#        np.nanmean(AIRS_data2['data']),AIRS_data2['data'])
    ##!#    nohash_data1 = np.ma.masked_where(AIRS_data2['data'] >= \
    ##!#        np.nanmean(AIRS_data2['data']),AIRS_data1['data'])
    ##!#

    ##!#    hash2 = ax1.pcolor(AIRS_data2['lon'],AIRS_data2['lat'],\
    ##!#        hash_data2, hatch = '///', alpha=0., transform = datacrs),
    ##!#    hash1 = ax0.pcolor(AIRS_data1['lon'],AIRS_data1['lat'],\
    ##!#        hash_data1, hatch = '///', alpha=0., transform = datacrs),
    
    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(dt_date_str.strftime('%Y%m%d%H%M')) 
        hash2 = ax1.pcolor(AIRS_data2['lon'],AIRS_data2['lat'],\
            hash_data1, hatch = '///', alpha=0., transform = datacrs),
        hash1 = ax0.pcolor(AIRS_data1['lon'],AIRS_data1['lat'],\
            hash_data1, hatch = '///', alpha=0., transform = datacrs),
    
    # Step 4: Plot scatter AIRS channel comparison in third panel
    # ------------------------------------------------------------

    # Figure out which pixels are in the plume
    # ----------------------------------------
    max_ch = 350.

    tmp_data1 = np.copy(AIRS_data1['data'])
    tmp_data2 = np.copy(AIRS_data2['data'])

    tmp_data1 = np.ma.masked_where( (AIRS_data1['data'] > max_ch) | \
        (AIRS_data2['data'] > max_ch), tmp_data1)
    tmp_data2 = np.ma.masked_where( (AIRS_data1['data'] > max_ch) | \
        (AIRS_data2['data'] > max_ch), tmp_data2)

    ###hash_data1 = np.ma.masked_where((hash_data1 > max_ch) | \
    ###    (hash_data2 > max_ch), hash_data1)
    ###hash_data2 = np.ma.masked_where((hash_data1 > max_ch) | \
    ###    (hash_data2 > max_ch), hash_data2)
    ###nohash_data1 = np.ma.masked_where((nohash_data1 > max_ch) | \
    ###    (nohash_data2 > max_ch), nohash_data1)
    ###nohash_data2 = np.ma.masked_where((nohash_data1 > max_ch) | \
    ###    (nohash_data2 > max_ch), nohash_data2)

    plot_data1 = tmp_data1.compressed()
    plot_data2 = tmp_data2.compressed()

    ###print("Hash 1 size:", hash_data1.compressed().shape, "  Hash 2 size: ", hash_data2.compressed().shape)
    ###print("No Hash 1 size:", nohash_data1.compressed().shape, "  No Hash 2 size: ", nohash_data2.compressed().shape)

    rval_p = pearsonr(plot_data1,plot_data2)[0]
    #rval_s = spearmanr(tmp_data1,tmp_data2)[0]
    print("Pearson:  ",rval_p)
    #print("Spearman: ",rval_s)

    xy = np.vstack([plot_data1,plot_data2])
    z = stats.gaussian_kde(xy)(xy)

    #ax2.scatter(plot_data1,plot_data2,c=z,s=6)
    ax2.scatter(AIRS_data1['data'][np.where(~hash_data1.mask)], \
        AIRS_data2['data'][np.where(~hash_data1.mask)], s=6, \
        color='tab:blue', label='Inside Plume')
    ax2.scatter(AIRS_data1['data'][np.where(hash_data1.mask)], \
        AIRS_data2['data'][np.where(hash_data1.mask)], s=6, \
        color='tab:orange', label='Outside Plume')
    #ax2.scatter(nohash_data1.compressed(), nohash_data2.compressed(), s=6, \
    #    color='tab:orange', label='Outside Plume')
    ax2.set_xlabel('Ch. ' + str(AIRS_data1['channel']) +' ' + AIRS_data1['variable'])
    ax2.set_ylabel('Ch. ' + str(AIRS_data2['channel']) +' ' + AIRS_data2['variable'])
    ax2.set_title('Pearson correlation: '+str(np.round(rval_p, 3)))
    ax2.legend()

    if(plot_ASOS_loc):
        print("Plotting ASOS site")
        plot_ASOS_locs(ax0,AIRS_data1)
        plot_ASOS_locs(ax1,AIRS_data1)

    if(save):
        cross_date = AIRS_data1['cross_date']
        file_time  = AIRS_data1['file_time']
        pdate = cross_date[:4] + cross_date[5:7] + cross_date[8:10] + file_time
        outname = 'modis_compare_ch' + str(channel1) + '_ch' + str(channel2) + '_' + pdate + '_2pannelscatter.png'
        plt.savefig(outname,dpi=300)
        print('Saved image',outname)
    else:
        plt.show()

    #return hash_data1, nohash_data1, AIRS_data1, AIRS_data2

def compare_AIRS_3scatter(date_str,channel0,channel1,channel2,channel3,\
        save=False, compare_OMI = False, compare_CERES = False, \
        avg_pixel = False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    filename = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis']

    if(channel1 == 'red'):
        channel1 = 1
    elif(channel1 == 'green'):
        channel1 = 4
    elif(channel1 == 'blue'):
        channel1 = 3

    # Step 1: Call read_AIRS_channel to read the desired AIRS data from the
    # file and put it in a dictionary
    # -----------------------------------------------------------------------
    AIRS_data0 = read_AIRS_channel(dt_date_str.strftime('%Y%m%d%H%M'), channel0, zoom = True)
    AIRS_data1 = read_AIRS_channel(dt_date_str.strftime('%Y%m%d%H%M'), channel1, zoom = True)
    AIRS_data2 = read_AIRS_channel(dt_date_str.strftime('%Y%m%d%H%M'), channel2, zoom = True)
    AIRS_data3 = read_AIRS_channel(dt_date_str.strftime('%Y%m%d%H%M'), channel3, zoom = True)

    print("Data0 max = ",np.max(AIRS_data0['data']), "  Data0 min = ",np.min(AIRS_data0['data']))
    print("Data1 max = ",np.max(AIRS_data1['data']), "  Data1 min = ",np.min(AIRS_data1['data']))
    print("Data2 max = ",np.max(AIRS_data2['data']), "  Data2 min = ",np.min(AIRS_data2['data']))
    print("Data3 max = ",np.max(AIRS_data3['data']), "  Data3 min = ",np.min(AIRS_data3['data']))

    max_ch = 350.

    # Determine where the smoke is located
    # ------------------------------------
    hash_data1, nohash_data1 = find_plume(dt_date_str.strftime('%Y%m%d%H%M')) 

    tmp_data0 = np.copy(AIRS_data0['data'])
    tmp_data1 = np.copy(AIRS_data1['data'])
    tmp_data2 = np.copy(AIRS_data2['data'])
    tmp_data3 = np.copy(AIRS_data3['data'])

    tmp_data0 = np.ma.masked_where( (AIRS_data0['data'] > max_ch) | \
        (AIRS_data1['data'] > max_ch) | (AIRS_data2['data'] > max_ch) | \
        (AIRS_data3['data'] > max_ch), tmp_data0)
    tmp_data1 = np.ma.masked_where( (AIRS_data0['data'] > max_ch) | \
        (AIRS_data1['data'] > max_ch) | (AIRS_data2['data'] > max_ch) | \
        (AIRS_data3['data'] > max_ch), tmp_data1)
    tmp_data2 = np.ma.masked_where( (AIRS_data0['data'] > max_ch) | \
        (AIRS_data1['data'] > max_ch) | (AIRS_data2['data'] > max_ch) | \
        (AIRS_data3['data'] > max_ch), tmp_data2)
    tmp_data3 = np.ma.masked_where( (AIRS_data0['data'] > max_ch) | \
        (AIRS_data1['data'] > max_ch) | (AIRS_data2['data'] > max_ch) | \
        (AIRS_data3['data'] > max_ch), tmp_data3)
    tmp_lat0  = np.ma.masked_where( (AIRS_data0['data'] > max_ch) | \
        (AIRS_data1['data'] > max_ch) | (AIRS_data2['data'] > max_ch) | \
        (AIRS_data3['data'] > max_ch), AIRS_data0['lat'])
    tmp_lon0  = np.ma.masked_where( (AIRS_data0['data'] > max_ch) | \
        (AIRS_data1['data'] > max_ch) | (AIRS_data2['data'] > max_ch) | \
        (AIRS_data3['data'] > max_ch), AIRS_data0['lon'])

    ##!#plot_data0 = tmp_data0.compressed()
    ##!#plot_data1 = tmp_data1.compressed()
    ##!#plot_data2 = tmp_data2.compressed()
    ##!#plot_data3 = tmp_data3.compressed()
    ##!#plot_lat0  = tmp_lat0.compressed()
    ##!#plot_lon0  = tmp_lon0.compressed()

    cross_date = AIRS_data0['cross_date']
    file_time  = AIRS_data0['file_time']

    plt.close('all')
    if(compare_OMI and compare_CERES):
        fig = plt.figure(figsize=(12,9))
        ax0 = fig.add_subplot(2,3,1)
        ax1 = fig.add_subplot(2,3,2)
        ax2 = fig.add_subplot(2,3,3)
        axo = fig.add_subplot(2,3,4)
        axcs = fig.add_subplot(2,3,5)
        axcl = fig.add_subplot(2,3,6)
    elif((compare_CERES and not compare_OMI)):
        fig = plt.figure(figsize=(17,9))
        ax0 = fig.add_subplot(2,3,1)
        ax1 = fig.add_subplot(2,3,2)
        ax2 = fig.add_subplot(2,3,3)
        axcs = fig.add_subplot(2,3,4)
        axcl = fig.add_subplot(2,3,5)
    elif((compare_OMI and not compare_CERES)):
        fig = plt.figure(figsize=(9,9))
        ax0 = fig.add_subplot(2,2,1)
        ax1 = fig.add_subplot(2,2,2)
        ax2 = fig.add_subplot(2,2,3)
        axo = fig.add_subplot(2,2,4)
    else:
        fig = plt.figure(figsize=(17,5))
        ax0 = fig.add_subplot(1,3,1)
        ax1 = fig.add_subplot(1,3,2)
        ax2 = fig.add_subplot(1,3,3)

    #xy = np.vstack([plot_data0,plot_data1])
    #z = stats.gaussian_kde(xy)(xy)
    #ax0.scatter(plot_data0,plot_data1,c=z,s=6)

    ax0.scatter(tmp_data0[np.where(~hash_data1.mask)], \
        tmp_data1[np.where(~hash_data1.mask)], s=6, \
        color='tab:blue')
    ax0.scatter(tmp_data0[np.where(hash_data1.mask)], \
        tmp_data1[np.where(hash_data1.mask)], s=6, \
        color='tab:orange')

    ax0.set_xlabel('Ch. ' + str(AIRS_data0['channel']) +' [' + \
        channel_dict[str(channel0)]['Bandwidth_label'] + \
        AIRS_data0['variable'])
    ax0.set_ylabel('Ch. ' + str(AIRS_data1['channel']) +' [' + \
        channel_dict[str(channel1)]['Bandwidth_label'] + \
        AIRS_data1['variable'])
    plt.suptitle('Aqua AIRS ' + cross_date + ' ' + file_time)
    #ax0.set_title('Pearson correlation: '+str(np.round(rval_p, 3)))

    #xy = np.vstack([plot_data0,plot_data2])
    #z = stats.gaussian_kde(xy)(xy)
    #ax1.scatter(plot_data0,plot_data2,c=z,s=6)
    ax1.scatter(tmp_data0[np.where(~hash_data1.mask)], \
        tmp_data2[np.where(~hash_data1.mask)], s=6, \
        color='tab:blue')
    ax1.scatter(tmp_data0[np.where(hash_data1.mask)], \
        tmp_data2[np.where(hash_data1.mask)], s=6, \
        color='tab:orange')

    ax1.set_xlabel('Ch. ' + str(AIRS_data0['channel']) +' [' + \
        channel_dict[str(channel0)]['Bandwidth_label'] + \
        AIRS_data0['variable'])
    ax1.set_ylabel('Ch. ' + str(AIRS_data2['channel']) +' [' + \
        channel_dict[str(channel2)]['Bandwidth_label'] + \
        AIRS_data2['variable'])
    #ax0.set_title('Pearson correlation: '+str(np.round(rval_p, 3)))

    #xy = np.vstack([plot_data0,plot_data3])
    #z = stats.gaussian_kde(xy)(xy)
    #ax2.scatter(plot_data0,plot_data3,c=z,s=6)
    ax2.scatter(tmp_data0[np.where(~hash_data1.mask)], \
        tmp_data3[np.where(~hash_data1.mask)], s=6, \
        color='tab:blue',label = 'Inside plume')
    ax2.scatter(tmp_data0[np.where(hash_data1.mask)], \
        tmp_data3[np.where(hash_data1.mask)], s=6, \
        color='tab:orange', label = 'Outside plume')

    ax2.set_xlabel('Ch. ' + str(AIRS_data0['channel']) +' [' + \
        channel_dict[str(channel0)]['Bandwidth_label'] + \
        AIRS_data0['variable'])
    ax2.set_ylabel('Ch. ' + str(AIRS_data3['channel']) +' [' + \
        channel_dict[str(channel3)]['Bandwidth_label'] + \
        AIRS_data3['variable'])
    #ax0.set_title('Pearson correlation: '+str(np.round(rval_p, 3)))

    # Remove the nans from the AIRS data, lats, and lons for both the
    # in-the-plume (hashed) and outside-the-plume (nohashed) data
    hash_plot_data0   = tmp_data0[np.where(~hash_data1.mask)].compressed()
    nohash_plot_data0 = tmp_data0[np.where(hash_data1.mask)].compressed()
    hash_plot_lat0   = tmp_lat0[np.where(~hash_data1.mask)].compressed()
    nohash_plot_lat0 = tmp_lat0[np.where(hash_data1.mask)].compressed()
    hash_plot_lon0   = tmp_lon0[np.where(~hash_data1.mask)].compressed()
    nohash_plot_lon0 = tmp_lon0[np.where(hash_data1.mask)].compressed()

    hash_plot_data1   = tmp_data1[np.where(~hash_data1.mask)].compressed()
    nohash_plot_data1 = tmp_data1[np.where(hash_data1.mask)].compressed()
    hash_plot_lat1   = tmp_lat0[np.where(~hash_data1.mask)].compressed()
    nohash_plot_lat1 = tmp_lat0[np.where(hash_data1.mask)].compressed()
    hash_plot_lon1   = tmp_lon0[np.where(~hash_data1.mask)].compressed()
    nohash_plot_lon1 = tmp_lon0[np.where(hash_data1.mask)].compressed()


    if(compare_OMI):
        print("Reading OMI data")
        data = h5py.File(plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['omi'],'r')
        LAT   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude'][:,:]
        LON   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude'][:,:]
        UVAI  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex'][:,:]
        XTRACK = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags'][:,:]
        mask_LAT = np.ma.masked_where((XTRACK < -2e5) | (UVAI < -2e5), LAT)
        mask_LON = np.ma.masked_where((XTRACK < -2e5) | (UVAI < -2e5), LON)
        mask_UVAI = np.ma.masked_where((XTRACK < -2e5) | (UVAI < -2e5), UVAI)
        mask_LAT  = np.ma.masked_where((((LAT < plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][0]) | \
                             (LAT > plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][1])) | \
                            ((LON < plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][0]) | \
                             (LON > plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][1]))), mask_LAT)
        mask_LON  = np.ma.masked_where((((LAT < plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][0]) | \
                             (LAT > plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][1])) | \
                            ((LON < plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][0]) | \
                             (LON > plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][1]))), mask_LON)
        mask_UVAI = np.ma.masked_where((((LAT < plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][0]) | \
                             (LAT > plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][1])) | \
                            ((LON < plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][0]) | \
                             (LON > plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][1]))), mask_UVAI)


        #hash_data1, nohash_data1 = find_plume(filename) 
        # Colocate the masked OMI data with the AIRS channel0 data
        # ---------------------------------------------------------
        print("Colocating OMI data")

        if(avg_pixel):
            print('averaging AIRS pixels')

            ##!## Determine the extreme bounds of the current 
            ##!#max_AIRS_lat = np.nanmax(tmp_lat0)
            ##!#min_AIRS_lat = np.nanmin(tmp_lat0)
            ##!#max_AIRS_lon = np.nanmax(tmp_lon0)
            ##!#min_AIRS_lon = np.nanmin(tmp_lon0)

            # Extract the OMI corner latitudes and longitudes
            # ----------------------------------------------- 
            crnr_LAT   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/'+\
                'Geolocation Fields/FoV75CornerLatitude'][:,:,:]
            crnr_LON   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/'+\
                'Geolocation Fields/FoV75CornerLongitude'][:,:,:]

            # Pull out the OMI lat/lon corner pixels that are associated with
            # the masked OMI data
            # ---------------------------------------------------------------
            work_LAT   = mask_LAT.compressed()
            work_LON   = mask_LON.compressed()
            work_UVAI  = mask_UVAI.compressed()

            work_crnrLAT = crnr_LAT[np.where(~np.isnan(mask_UVAI)==True)[0],\
                np.where(~np.isnan(mask_UVAI)==True)[1],:] 
            work_crnrLON = crnr_LON[np.where(~np.isnan(mask_UVAI)==True)[0],\
                np.where(~np.isnan(mask_UVAI)==True)[1],:] 

            # Declare full arrays to hold the averaged AIRS data
            # Must be dimensioned to the shape of the OMI data
            # ---------------------------------------------------
            print("work UVAI shape = ",work_UVAI.shape)
            hash_avg_modis   = np.full(work_UVAI.shape,np.nan)
            hash_avg_mlat    = np.full(work_UVAI.shape,np.nan)
            hash_avg_mlon    = np.full(work_UVAI.shape,np.nan)
            nohash_avg_modis = np.full(work_UVAI.shape,np.nan)
            nohash_avg_mlat  = np.full(work_UVAI.shape,np.nan)
            nohash_avg_mlon  = np.full(work_UVAI.shape,np.nan)

            # Loop over screened out OMI pixels
            # ---------------------------------
            print("corners shape:",work_crnrLAT.shape)
            for ii in range(work_crnrLAT.shape[0]):

                # For each OMI pixel, find the in-plume and out-of-plume AIRS
                # pixels that are close to the OMI corner values
                # -------------------------------------------------------------
                hash_in = np.where((hash_plot_lat0 >= np.nanmin(work_crnrLAT[ii,:])) & \
                                   (hash_plot_lat0 <= np.nanmax(work_crnrLAT[ii,:])) & \
                                   (hash_plot_lon0 >= np.nanmin(work_crnrLON[ii,:])) & \
                                   (hash_plot_lon0 <= np.nanmax(work_crnrLON[ii,:])) )
                nohash_in = np.where((nohash_plot_lat0 >= np.nanmin(work_crnrLAT[ii,:])) & \
                                     (nohash_plot_lat0 <= np.nanmax(work_crnrLAT[ii,:])) & \
                                     (nohash_plot_lon0 >= np.nanmin(work_crnrLON[ii,:])) & \
                                     (nohash_plot_lon0 <= np.nanmax(work_crnrLON[ii,:])) )

                # Use the indices to pull out the AIRS pixels
                hash_modis_data   = hash_plot_data1[hash_in]
                hash_modis_lat    = hash_plot_lat1[hash_in]
                hash_modis_lon    = hash_plot_lon1[hash_in]
                nohash_modis_data = nohash_plot_data1[nohash_in]
                nohash_modis_lat  = nohash_plot_lat1[nohash_in]
                nohash_modis_lon  = nohash_plot_lon1[nohash_in]

                #print(work_crnrLAT[ii,:], work_crnrLON[ii,:], hash_modis_lat, \
                #    nohash_modis_lat)
            
                # Use the corner values to make a Polygon and determine if
                # each close AIRS pixel is within the OMI pixel
                # --------------------------------------------------------

                # Generate a Polygon using the OMI corner points
                points = []
                for x, y in zip(work_crnrLON[ii,:], work_crnrLAT[ii,:]):
                    points.append((x,y))
                omi_poly = Polygon(points)
              
                inside_hash = [omi_poly.contains(Point(hlon, hlat)) for \
                    hlon, hlat in zip(hash_modis_lon, hash_modis_lat)]
                inside_nohash = [omi_poly.contains(Point(nlon, nlat)) for \
                    nlon, nlat in zip(nohash_modis_lon, nohash_modis_lat)]
       
                if((len(inside_hash) > 0) & (True in inside_hash)):
                    hash_avg_modis[ii] = np.average(hash_modis_data[inside_hash])
                    hash_avg_mlat[ii]  = np.average(hash_modis_lat[inside_hash])
                    hash_avg_mlon[ii]  = np.average(hash_modis_lon[inside_hash])

                if((len(inside_nohash) > 0) & (True in inside_nohash)):
                    nohash_avg_modis[ii] = np.average(nohash_modis_data[inside_nohash])
                    nohash_avg_mlat[ii]  = np.average(nohash_modis_lat[inside_nohash])
                    nohash_avg_mlon[ii]  = np.average(nohash_modis_lon[inside_nohash])

                ##!## Loop over the close AIRS pixels and determine which ones
                ##!## are within the polygon
                ##!## ---------------------------------------------------------
                ##!#for jj in range(len(hash_modis_data)):

                ##!#    # Average each AIRS value that is within the polygon
                ##!#    # into a value, and insert the value into the avg
                ##!#    # modis array for hash and nohash.
                ##!#    # --------------------------------------------------- 
                ##!#    mpoint = Point(hash_modis_lon[jj], hash_modis_lat[jj])

                ##!#    inside_vals = [omi_poly.contains(mpoint)  ] 

            # end screened corner data loop

            # Plot the screened OMI data vs the averaged AIRS data
            # -----------------------------------------------------

            #print('hashed averaged AIRS data = ',hash_avg_modis)
            #print('nohashed averaged AIRS data = ',nohash_avg_modis)

            print(work_UVAI[~np.ma.masked_invalid(hash_avg_modis).mask])
            
            axo.scatter(np.ma.masked_invalid(hash_avg_modis).compressed(), \
                work_UVAI[~np.ma.masked_invalid(hash_avg_modis).mask], s=6, \
                color='tab:blue')
            #axo.scatter(np.ma.masked_invalid(nohash_avg_modis).compressed(), \
            #    work_UVAI[~np.ma.masked_invalid(nohash_avg_modis).mask], s=6, \
            #    color='tab:orange')
            #axo.scatter(nohash_plot_data0, nohash_match_OMI, s=6, \
            #    color='tab:orange', label='Outside Plume')

            axo.set_xlabel('Averaged Ch. ' + str(AIRS_data1['channel']) +' [' + \
                channel_dict[str(channel1)]['Bandwidth_label'] + \
                AIRS_data1['variable'])
            axo.set_ylabel('OMI UVAI')

            # Plot hash_avg_modis against mask_UVAI 

            # Plot nohash_avg_modis against mask_UVAI 
 
        else:
            hash_match_OMI   = np.full(hash_plot_data1.shape,-9.)
            hash_match_LAT   = np.full(hash_plot_lat1.shape,-9.)
            hash_match_LON   = np.full(hash_plot_lon1.shape,-9.)
            nohash_match_OMI = np.full(nohash_plot_data1.shape,-9.)
            nohash_match_LAT = np.full(nohash_plot_lat1.shape,-9.)
            nohash_match_LON = np.full(nohash_plot_lon1.shape,-9.)

            print(hash_plot_data0.shape)
            for ii in range(hash_match_OMI.shape[0]):
                # Find the gridpoint in the gridded lat/lon data that 
                # corresponds to the station at slat and slon
                # ---------------------------------------------------- 
                o_idx = nearest_gridpoint(hash_plot_lat1[ii], hash_plot_lon1[ii],\
                    mask_LAT, mask_LON)

                if(len(o_idx[0]) > 1):
                    o_idx = (np.array([o_idx[0][0]])), (np.array([o_idx[1][0]]))
                hash_match_OMI[ii] = mask_UVAI[o_idx]
                hash_match_LAT[ii] = mask_LAT[o_idx] 
                hash_match_LON[ii] = mask_LON[o_idx] 

            print(nohash_plot_data1.shape)
            for ii in range(nohash_match_OMI.shape[0]):
                # Find the gridpoint in the gridded lat/lon data that 
                # corresponds to the station at slat and slon
                # ---------------------------------------------------- 
                o_idx = nearest_gridpoint(nohash_plot_lat1[ii], nohash_plot_lon1[ii],\
                    mask_LAT, mask_LON)

                if(len(o_idx[0]) > 1):
                    o_idx = (np.array([o_idx[0][0]])), (np.array([o_idx[1][0]]))
                nohash_match_OMI[ii] = mask_UVAI[o_idx]
                nohash_match_LAT[ii] = mask_LAT[o_idx] 
                nohash_match_LON[ii] = mask_LON[o_idx] 
 
            #xy = np.vstack([plot_data0,match_OMI])
            #z = stats.gaussian_kde(xy)(xy)
            #axo.scatter(plot_data0,match_OMI,c=z,s=6)

            axo.scatter(hash_plot_data1, hash_match_OMI, s=6, \
                color='tab:blue')
            #axo.scatter(nohash_plot_data0, nohash_match_OMI, s=6, \
            #    color='tab:orange')

            axo.set_xlabel('Ch. ' + str(AIRS_data1['channel']) +' [' + \
                channel_dict[str(channel1)]['Bandwidth_label'] + \
                AIRS_data1['variable'])
            axo.set_ylabel('OMI UVAI')
        
        data.close()

    # End compare_OMI

    if(compare_CERES):
        print("Reading CERES data")

        # NOTE: need to use the time variable to screen out any data
        # that are not close to the event time.
        base_date = datetime(year=1970,month=1,day=1)
        start_date = dt_date_str - timedelta(hours = 1)
        end_date   = dt_date_str + timedelta(hours = 2)

        print(start_date, end_date)

        data = Dataset(plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['ceres'],'r')
        LAT   = 90. - data.variables['Colatitude_of_CERES_FOV_at_surface'][:]
        LON   = data.variables['Longitude_of_CERES_FOV_at_surface'][:]
        LON[LON>179.99] = -360.+LON[LON>179.99]
        swflux  = data.variables['CERES_SW_TOA_flux___upwards'][:]
        lwflux  = data.variables['CERES_LW_TOA_flux___upwards'][:]
        time  = data.variables['time'][:]
        local_time = np.array([base_date + relativedelta(days = ttime) for ttime in time])
        data.close()

        mask_LAT = LAT[ \
            (LAT >= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][0]) & \
            (LAT <= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][1]) & \
            (LON >= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][0]) & \
            (LON <= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][1]) & \
            (swflux > 0) & (lwflux > 0)]
        mask_LON = LON[ \
            (LAT >= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][0]) & \
            (LAT <= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][1]) & \
            (LON >= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][0]) & \
            (LON <= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][1]) & \
            (swflux > 0) & (lwflux > 0)]
        mask_swf = swflux[ \
            (LAT >= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][0]) & \
            (LAT <= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][1]) & \
            (LON >= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][0]) & \
            (LON <= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][1]) & \
            (swflux > 0) & (lwflux > 0)]
        mask_lwf = lwflux[ \
            (LAT >= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][0]) & \
            (LAT <= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][1]) & \
            (LON >= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][0]) & \
            (LON <= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][1]) & \
            (swflux > 0) & (lwflux > 0)]
        mask_time = local_time[ \
            (LAT >= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][0]) & \
            (LAT <= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lat'][1]) & \
            (LON >= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][0]) & \
            (LON <= plot_limits_dict[AIRS_data3['cross_date']][AIRS_data3['file_time']]['Lon'][1]) & \
            (swflux > 0) & (lwflux > 0)]


        if(avg_pixel):
            print('(NOT) averaging AIRS pixels')
        else:
            hash_match_SWF   = np.full(hash_plot_data0.shape,-9.)
            hash_match_LWF   = np.full(hash_plot_data0.shape,-9.)
            hash_match_LAT   = np.full(hash_plot_lat0.shape,-9.)
            hash_match_LON   = np.full(hash_plot_lon0.shape,-9.)
            nohash_match_SWF = np.full(nohash_plot_data0.shape,-9.)
            nohash_match_LWF = np.full(nohash_plot_data0.shape,-9.)
            nohash_match_LAT = np.full(nohash_plot_lat0.shape,-9.)
            nohash_match_LON = np.full(nohash_plot_lon0.shape,-9.)

            print(hash_plot_data0.shape)
            for ii in range(hash_match_SWF.shape[0]):
                # Find the gridpoint in the gridded lat/lon data that 
                # corresponds to the station at slat and slon
                # ---------------------------------------------------- 
                o_idx = nearest_gridpoint(hash_plot_lat0[ii], hash_plot_lon0[ii],\
                    mask_LAT, mask_LON)

                if(len(o_idx[0]) > 1):
                    o_idx = (np.array([o_idx[0][0]])), (np.array([o_idx[1][0]]))
                hash_match_SWF[ii] = mask_swf[o_idx]
                hash_match_LWF[ii] = mask_lwf[o_idx]
                hash_match_LAT[ii] = mask_LAT[o_idx] 
                hash_match_LON[ii] = mask_LON[o_idx] 

            print(nohash_plot_data0.shape)
            for ii in range(nohash_match_SWF.shape[0]):
                # Find the gridpoint in the gridded lat/lon data that 
                # corresponds to the station at slat and slon
                # ---------------------------------------------------- 
                o_idx = nearest_gridpoint(nohash_plot_lat0[ii], nohash_plot_lon0[ii],\
                    mask_LAT, mask_LON)

                if(len(o_idx[0]) > 1):
                    o_idx = (np.array([o_idx[0][0]])), (np.array([o_idx[1][0]]))
                nohash_match_SWF[ii] = mask_swf[o_idx]
                nohash_match_LWF[ii] = mask_lwf[o_idx]
                nohash_match_LAT[ii] = mask_LAT[o_idx] 
                nohash_match_LON[ii] = mask_LON[o_idx] 
 
            #xy = np.vstack([plot_data0,match_OMI])
            #z = stats.gaussian_kde(xy)(xy)
            #axo.scatter(plot_data0,match_OMI,c=z,s=6)

            axcl.scatter(hash_plot_data0, hash_match_LWF,\
                s = 6, color='tab:blue')
            axcs.scatter(hash_plot_data0, hash_match_SWF,\
                s = 6, color='tab:blue')
            ##!#axcl.scatter(nohash_plot_data0, nohash_match_LWF,\
            ##!#    s = 6, color='tab:orange')
            ##!#axcs.scatter(nohash_plot_data0, nohash_match_SWF,\
            ##!#    s = 6, color='tab:orange')

            axcl.set_xlabel('Ch. ' + str(AIRS_data0['channel']) +' [' + \
                channel_dict[str(channel0)]['Bandwidth_label'] + \
                AIRS_data0['variable'])
            axcs.set_xlabel('Ch. ' + str(AIRS_data0['channel']) +' [' + \
                channel_dict[str(channel0)]['Bandwidth_label'] + \
                AIRS_data0['variable'])
            axcl.set_ylabel('CERES LWF [W/m2]')
            axcs.set_ylabel('CERES SWF [W/m2]')

    # end compare_CERES
       
    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    if(compare_OMI):
        plt.legend(lines, labels, loc = 'lower center', bbox_to_anchor = (0, 0.01, 1, 1),\
            bbox_transform = plt.gcf().transFigure, ncol=2)
    else:
        plt.legend(lines, labels, loc = 'right', bbox_to_anchor = (0, 0., 1.0, 1),\
            bbox_transform = plt.gcf().transFigure, ncol=1)

    if(save):
        pdate = cross_date[:4] + cross_date[5:7] + cross_date[8:10] + file_time
        outname = 'modis_compare_ch' + str(channel0) + '_vs_ch' + \
            str(channel1) + '_ch' + str(channel2) + '_ch' + \
            str(channel3) + '_' + pdate + '_3scatter.png'
        if(compare_OMI):
            outname = 'modis_compare_ch' + str(channel0) + '_vs_ch' + \
                str(channel1) + '_ch' + str(channel2) + '_ch' + \
                str(channel3) + '_omi_' + pdate + '_3scatter.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else: 
        plt.show()


# Compare colocated AIRS and ob data for two dates
def colocate_comparison(date1, date2, channel = 31):
    # Read in AIRS data for both cases
    # ---------------------------------
    dt_date_str1 = datetime.strptime(date1,"%Y%m%d%H%M")
    dt_date_str2 = datetime.strptime(date2,"%Y%m%d%H%M")
    filename1 = plot_limits_dict[dt_date_str1.strftime('%Y-%m-%d')][dt_date_str1.strftime('%H%M')]['modis']
    filename2 = plot_limits_dict[dt_date_str2.strftime('%Y-%m-%d')][dt_date_str2.strftime('%H%M')]['modis']

    AIRS_data1 = read_AIRS_channel(dt_date_str1.strftime('%Y%m%d%H%M'), channel, zoom = True)
    AIRS_data2 = read_AIRS_channel(dt_date_str2.strftime('%Y%m%d%H%M'), channel, zoom = True)

    #print(AIRS_data1,AIRS_data2)

    # Use the colocation code to extract matching data
    # ------------------------------------------------
    compare_data1 = nearest_grid_values(AIRS_data1)
    compare_data2 = nearest_grid_values(AIRS_data2)

    # Loop over the data and print statistics
    # ---------------------------------------
    for ii in range(len(compare_data1['stn_data'])):
        stn_diff = compare_data2['stn_data'][ii] - compare_data1['stn_data'][ii]
        mds_diff = compare_data2['mds_data'][ii] - compare_data1['mds_data'][ii]
        print(compare_data1['stations'][ii], stn_diff, mds_diff)

    # Plot the data
    # -------------

    plt.close('all')
    fig, ax = plt.subplots(figsize=(8,5))
    print(compare_data1['stations'])
    ax.bar(np.arange(len(compare_data1['stations'])), \
        compare_data1['stn_data'], color='tab:blue', \
        label = 'Obs ' + dt_date_str1.strftime('%d/%m/%y'))
    ax.bar(np.arange(len(compare_data1['stations'])), \
        compare_data2['stn_data'], color = 'tab:blue', \
        linestyle = '--', label = 'Obs ' + dt_date_str2.strftime('%d/%m/%y'))

    ax2 = ax.twinx()
    ax2.bar(np.arange(len(compare_data1['stations'])),\
        compare_data1['mds_data'], color='tab:orange', label = '05 Modis')
    ax2.bar(np.arange(len(compare_data1['stations'])),\
        compare_data2['mds_data'], color='tab:orange', linestyle = '--', \
        label = '06 Modis')

    ##!#lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    ##!#lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    ##!#if(compare_OMI):
    ##!#    plt.legend(lines, labels, loc = 'lower center', bbox_to_anchor = (0, 0.01, 1, 1),\
    ##!#        bbox_transform = plt.gcf().transFigure, ncol=2)
    ##!#else:
    ##!#    plt.legend(lines, labels, loc = 'right', bbox_to_anchor = (0, 0., 1.0, 1),\
    ##!#        bbox_transform = plt.gcf().transFigure, ncol=1)
 
    ax.set_xticks(np.arange(len(compare_data1['stations']))) 
    ax.set_xticklabels(compare_data1['stations'])
    ax.set_ylabel('Observed 2m Temperature [degC]',color='tab:blue') 
    ax2.set_ylabel('Colocated AIRS Channel ' + str(channel) + ' Brightness Temp [K] ',color='tab:orange') 
    plt.legend() 
    plt.show()
