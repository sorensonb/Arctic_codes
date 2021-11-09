"""
  NAME:

  PURPOSE:
  
    NOTE: The MODIS channel and true color functions are designed to work with
    HDF MODIS files retriefed from 
    the data ordering website at this address:
    https://ladsweb.modaps.eosdis.nasa.gov/search/order/1/MODIS:Aqua


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
mapcrs = ccrs.LambertConformal()
datacrs = ccrs.PlateCarree()

zoom_dict = {
    'Finland': [10,55,65,80]
}

proj_dict = {
    'Finland': ccrs.NorthPolarStereo(central_longitude = 35.) 
}

channel_dict = {
    '1': {
        'Name': 'EV_250_Aggr1km_RefSB',\
        'Index': 0,\
        'Bandwidth': [0.620, 0.670]
    },\
    '2': {
        'Name': 'EV_250_Aggr1km_RefSB',\
        'Index': 1,\
        'Bandwidth': [0.841, 0.876]
    },\
    '3': {
        'Name': 'EV_500_Aggr1km_RefSB',\
        'Index': 0,\
        'Bandwidth': [0.459, 0.479]
    },\
    '4': {
        'Name': 'EV_500_Aggr1km_RefSB',\
        'Index': 1,\
        'Bandwidth': [0.545, 0.565]
    },\
    '5': {
        'Name': 'EV_500_Aggr1km_RefSB',\
        'Index': 2,\
        'Bandwidth': [1.230, 1.250]
    },\
    '6': {
        'Name': 'EV_500_Aggr1km_RefSB',\
        'Index': 3,\
        'Bandwidth': [1.628, 1.652]
    },\
    '7': {
        'Name': 'EV_500_Aggr1km_RefSB',\
        'Index': 4,\
        'Bandwidth': [2.105, 2.155]
    },\
    '8': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 0,\
        'Bandwidth': [0.405, 0.420]
    },\
    '9': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 1,\
        'Bandwidth': [0.438, 0.448]
    },\
    '10': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 2,\
        'Bandwidth': [0.483, 0.493]
    },\
    '11': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 3,\
        'Bandwidth': [0.526, 0.536]
    },\
    '12': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 4,\
        'Bandwidth': [0.546, 0.556]
    },\
    '13lo': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 5,\
        'Bandwidth': [0.662, 0.672]
    },\
    '13hi': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 6,\
        'Bandwidth': [0.662, 0.672]
    },\
    '14lo': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 7,\
        'Bandwidth': [0.673, 0.683]
    },\
    '14hi': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 8,\
        'Bandwidth': [0.673, 0.683]
    },\
    '15': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 9,\
        'Bandwidth': [0.743, 0.753]
    },\
    '16': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 10,\
        'Bandwidth': [0.862, 0.877]
    },\
    '17': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 11,\
        'Bandwidth': [0.890, 0.920]
    },\
    '18': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 12,\
        'Bandwidth': [0.931, 0.941]
    },\
    '19': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 13,\
        'Bandwidth': [0.915, 0.965]
    },\
    '20': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 0,\
        'Bandwidth': [3.660, 3.840]
    },\
    '21': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 1,\
        'Bandwidth': [3.929, 3.989]
    },\
    '22': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 2,\
        'Bandwidth': [3.929, 3.989]
    },\
    '23': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 3,\
        'Bandwidth': [4.020, 4.080]
    },\
    '24': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 4,\
        'Bandwidth': [4.433, 4.498]
    },\
    '25': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 5,\
        'Bandwidth': [4.482, 4.549]
    },\
    '26': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 14,\
        'Bandwidth': [1.360, 1.390]
    },\
    '27': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 6,\
        'Bandwidth': [6.535, 6.895]
    },\
    '28': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 7,\
        'Bandwidth': [7.175, 7.475]
    },\
    '29': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 8,\
        'Bandwidth': [8.400, 8.700]
    },\
    '30': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 9,\
        'Bandwidth': [9.580, 9.880]
    },\
    '31': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 10,\
        'Bandwidth': [10.780, 11.280]
    },\
    '32': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 11,\
        'Bandwidth': [11.770, 12.270]
    },\
    '33': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 12,\
        'Bandwidth': [13.185, 13.485]
    },\
    '34': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 13,\
        'Bandwidth': [13.485, 13.785]
    },\
    '35': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 14,\
        'Bandwidth': [13.785, 14.085]
    },\
    '36': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 15,\
        'Bandwidth': [14.085, 14.385]
    },\
    'wv_ir': {
        'Name': 'Water_Vapor_Infrared',\
        'Index': None,
        'Bandwidth': None 
    },\
    'wv_nir': {
        'Name': 'Water_Vapor_Near_Infrared',\
        'Index': None,
        'Bandwidth': None 
    }
}

for key in channel_dict.keys():
    if(channel_dict[key]['Bandwidth'] is not None):
        channel_dict[key]['Bandwidth_label'] = \
            str(channel_dict[key]['Bandwidth'][0]) + ' μm - ' + \
            str(channel_dict[key]['Bandwidth'][1]) + ' μm'
    else:
        channel_dict[key]['Bandwidth_label'] = ''

plot_limits_dict = {
    "2021-07-13": {
        '2110': {
            'asos': 'asos_data_20210713.csv',
            #'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021203.2110.061.2021204155922.hdf',
            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
            #'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072210-2021072221.nc',
            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
            'Lat': [39.5, 42.0],
            'Lon': [-122.0, -119.5],
            'modis_Lat': [39.0, 42.5],
            'modis_Lon': [-123., -119.]
        }
    },
    "2021-07-20": {
        '2125': {
            'asos': 'asos_data_20210720.csv',
            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021201.2125.061.2021202154814.hdf',
            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072010-2021072021.nc',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.20.214.L2.SUBS2RET.v6.0.32.0.G21202153435.hdf'],
            'Lat': [39.5, 42.0],
            'Lon': [-122.0, -119.5]
        }
    },
    "2021-07-21": {
        '2030': {
            'asos': 'asos_data_20210722.csv',
            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021202.2030.061.2021203174050.hdf',
            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072109-2021072122.nc',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.21.205.L2.SUBS2RET.v6.0.32.0.G21203185004.hdf'],
            'Lat': [39.5, 42.0],
            'Lon': [-122.0, -119.5],
            'modis_Lat': [39.0, 42.5],
            'modis_Lon': [-123., -119.]
        }
    },
    "2021-07-22": {
        '2110': {
            'asos': 'asos_data_20210722.csv',
            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021203.2110.061.2021204155922.hdf',
            'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
            'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072210-2021072221.nc',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
            'Lat': [39.5, 42.0],
            'Lon': [-122.0, -119.5],
            'modis_Lat': [39.0, 42.5],
            'modis_Lon': [-123., -119.]
        }
    },
    "2021-08-04": {
        '2110': {
            'asos': 'asos_data_20210806.csv',
            #'asos': 'asos_nevada_20210806.csv',
            #'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021218.2025.061.2021219151802.hdf',
            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021218.2025.061.2021219152751.hdf',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.04.206.L2.SUBS2RET.v6.0.32.0.G21217152448.hdf',\
                     '/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.04.207.L2.SUBS2RET.v6.0.32.0.G21217152904.hdf'],\
            #'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
            'Lat': [36.0, 39.0],
            'Lon': [-118.0, -114.0],
            'modis_Lat': [35.0, 40.0],
            'modis_Lon': [-119., -113.]
        }
    },
    "2021-08-05": {
        '2120': {
            'asos': 'asos_data_20210805.csv',
            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021217.2120.061.2021218164201.hdf',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.05.214.L2.SUBS2RET.v6.0.32.0.G21218175548.hdf'],
            'Lat': [36.0, 39.0],
            'Lon': [-118.0, -114.0],
            'modis_Lat': [35.0, 40.0],
            'modis_Lon': [-119., -113.]
        },
        '2125': {
            'asos': 'asos_california_20210805.csv',
            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021217.2125.061.2021218161010.hdf',
            'ceres': '/home/bsorenson/data/CERES/FLASHFlux/Aqua/FLASH_SSF_Aqua_Version4A_Subset_2021080508-2021080521.nc',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.05.214.L2.SUBS2RET.v6.0.32.0.G21218175548.hdf'],
            'Lat': [39.5, 42.5],
            'Lon': [-121.5, -119.5],
            'modis_Lat': [39.0, 42.5],
            'modis_Lon': [-123., -119.]
        }
    },
    "2021-08-06": {
        '2025': {
            'asos': 'asos_data_20210806.csv',
            #'asos': 'asos_nevada_20210806.csv',
            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021218.2025.061.2021219151802.hdf',
            'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021218.2025.061.2021219152751.hdf',
            'ceres': '/home/bsorenson/data/CERES/FLASHFlux/Aqua/FLASH_SSF_Aqua_Version4A_Subset_2021080609-2021080620.nc',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.06.204.L2.SUBS2RET.v6.0.32.0.G21219130523.hdf',\
                     '/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.06.205.L2.SUBS2RET.v6.0.32.0.G21219130455.hdf'],
            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
            'Lat': [36.0, 39.0],
            'Lon': [-118.0, -114.0],
            'modis_Lat': [35.0, 40.0],
            'modis_Lon': [-119., -113.]
        }
    },
    "2021-08-07": {
        '2110': {
            'asos': 'asos_data_20210806.csv',
            #'asos': 'asos_nevada_20210806.csv',
            #'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021218.2025.061.2021219151802.hdf',
            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021218.2025.061.2021219152751.hdf',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.07.212.L2.SUBS2RET.v6.0.32.0.G21220123225.hdf'],\
            #'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
            'Lat': [36.0, 39.0],
            'Lon': [-118.0, -114.0],
            'modis_Lat': [35.0, 40.0],
            'modis_Lon': [-119., -113.]
        }
    },
    "2021-08-08": {
        '2110': {
            'asos': 'asos_data_20210806.csv',
            #'asos': 'asos_nevada_20210806.csv',
            #'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021218.2025.061.2021219151802.hdf',
            #'mdswv': '/home/bsorenson/data/MODIS/Aqua/MYD05_L2.A2021218.2025.061.2021219152751.hdf',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.08.202.L2.SUBS2RET.v6.0.32.0.G21221124420.hdf',\
                     '/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.08.203.L2.SUBS2RET.v6.0.32.0.G21221185932.hdf'],\
            #'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
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
            'modis': '/home/bsorenson/data/MODIS/Aqua/MYD021KM.A2021242.2115.061.2021243183953.hdf',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.30.213.L2.SUBS2RET.v6.0.32.0.G21243151114.hdf'],
            'Lat': [38.0, 40.0],
            'Lon': [-121.0, -118.5]
        }
    } 
}

def getCorners(centers):
    one = centers[:-1,:]
    two = centers[1:,:]
    d1 = (two - one) / 2.0
    one = one - d1
    two = two + d1
    stepOne = np.zeros((centers.shape[0] + 1, centers.shape[1]))
    stepOne[:-2,:] = one
    stepOne[-2:,:] = two[-2:,:]
    one = stepOne[:,:-1]
    two = stepOne[:,1:]
    d2 = (two - one) / 2.
    one = one - d2
    two = two + d2
    stepTwo = np.zeros((centers.shape[0] + 1, centers.shape[1] + 1))
    stepTwo[:,:-2] = one
    stepTwo[:,-2:] = two[:,-2:]
    return stepTwo

# Find the gridpoint in the gridded lat/lon data that 
# corresponds to the station at slat and slon
# ---------------------------------------------------- 
def nearest_gridpoint(slat, slon, grid_lat, grid_lon):
    fun_c = np.maximum(np.abs(grid_lat - slat), \
        np.abs(grid_lon - slon))
    m_idx = np.where(fun_c == np.min(fun_c))
    return m_idx
 
# Extract the MODIS information from a given channel at each ob point
# -------------------------------------------------------------------
def nearest_grid_values(MODIS_data):
    # Read in the correct ASOS file 
    asos_file = plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['asos']
    df = pd.read_csv(asos_file)
    df['valid'] = pd.to_datetime(df['valid'])
    df = df.set_index('valid')

    # Pull the event time from the plot_limits_dict
    event_date = datetime.strptime(MODIS_data['cross_date'], "%Y-%m-%d")
    first_time = MODIS_data['file_time']
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
        m_idx = nearest_gridpoint(lat_stn, lon_stn, MODIS_data['lat'], \
            MODIS_data['lon'])
        m_data = MODIS_data['data'][m_idx][0]
      
        compare_dict['mds_data'][ii] = m_data
 
        ##print(station, lat_stn, MODIS_data['lat'][m_idx][0], \
        ##    lon_stn, MODIS_data['lon'][m_idx][0], compare_dict['stn_data'][ii], m_data )

    return compare_dict

 
# Plot the downloaded ASOS stations for each case
# ----------------------------------------------- 
def plot_ASOS_locs(pax,MODIS_data,color='red'):
    # Read in the correct ASOS file 
    asos_file = plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['asos']
    df = pd.read_csv(asos_file)

    station_names = set(df['station'].values)
    for ii, station in enumerate(station_names):
        lat_stn = df['lat'][df['station'] == station].values[0]
        lon_stn = df['lon'][df['station'] == station].values[0]
        pax.text(lon_stn, lat_stn, station, fontsize=10,weight='bold', \
            transform=datacrs, color=color)

# Determine areas of an image that are in smoke, defined by:
#    (ch1_refl - ch5_refl) < mean(ch1_refl - ch5_refl) * 0.25
def find_plume(dt_date_str):
    MODIS_ch1  = read_MODIS_channel(dt_date_str, 1,  zoom = True)
    MODIS_ch5  = read_MODIS_channel(dt_date_str, 5,  zoom = True)
    MODIS_ch31 = read_MODIS_channel(dt_date_str, 31, zoom = True)

    screen_limit = 0.20
    max_ch = 350.
    test_data = MODIS_ch1['data'] - MODIS_ch5['data']
    hash_data   = np.ma.masked_where(test_data <  \
        (np.nanmean(test_data) * screen_limit), MODIS_ch1['data'])
    nohash_data = np.ma.masked_where(test_data >= \
        (np.nanmean(test_data) * screen_limit), MODIS_ch1['data'])

    hash_data = np.ma.masked_where(MODIS_ch31['data'] > max_ch, hash_data)
    nohash_data = np.ma.masked_where(MODIS_ch31['data'] > max_ch, nohash_data)

    return hash_data, nohash_data

# Start_date and end_date must be formatted as 
# "YYYYMMDD"
# Variable must be either "CHL" or "PIC"
def read_MODIS(variable,start_date,latmin = 60):

    ##if(monthly == True):
    ##    dformat = '%Y%m' 
    ##    end_string_idx = -5
    ##else:
    ##    dformat = '%Y%m%d' 
    ##    end_string_idx = -3
    dformat = '%Y%m%d' 
    # Set up starting and ending datetime objects
    sdate = datetime.strptime(start_date,dformat)
    ##edate = datetime.strptime(end_date,dformat)
 
    
    # Grab all the modis files
    if(variable == 'CHL'):
        grabber = 'chlor-a'
        label = 'Chlorophyll Concentration, OCI Algorithm (mg m$^{-3}$)'
        ptitle = 'Chlorophyll-α'
        pticks = [0.01,0.02,0.05,0.1,0.2,0.5,1.0,2.0,5.0,10.0,20.0]
        ptick_labels = ['0.01','0.02','0.05','0.1','0.2','0.5','1','2','5',\
            '10','20']
    elif(variable == 'PIC'):
        grabber = 'pic'
        label = 'Calcite Concentration, Balch and Gordon (mol m$^{-3}$)'
        ptitle = 'Particulate Inorganic Carbon'
        pticks = [1e-5,2e-5,5e-5,1e-4,2e-4,5e-4,1e-3,2e-3,5e-3,0.01,0.02,0.05]
        ptick_labels = ['1e-5','2e-5','5e-5','1e-4','2e-4','5e-4','1e-3',\
            '2e-3','5e-3','0.01','0.02','0.05']
    else:
        print("ERROR: variable must be CHL or PIC")
        return
    #cmnd = "ls /home/bsorenson/data/MODIS/CHL/A*.nc"
    cmnd = "ls /home/bsorenson/data/MODIS/"+variable+"/A*.nc"
    
    
    file_initial = subprocess.check_output(cmnd,shell=True).decode('utf-8').strip().split('\n')
    file_names = []

    for fname in file_initial:
        splitter = fname.split('/')[-1]
        fdate = datetime.strptime(splitter[1:5],'%Y') + relativedelta(days = int(splitter[5:8])-1)
        if(fdate == sdate):
            file_names.append(fname)
   
    lat_ranges = np.arange(latmin,90.25,0.25)
    lon_ranges = np.arange(-180.,180.25,0.25)
 
    # Read in the latitude, longitude, and area data
    num_files = len(file_names)
    modis_data = {}
    modis_data['data'] = np.full((num_files,len(lat_ranges),len(lon_ranges)),-9.)
    modis_data['lat']  = lat_ranges
    modis_data['lon']  = lon_ranges
    modis_data['variable'] = variable
    modis_data['ptitle'] = ptitle
    modis_data['label'] = label
    modis_data['grabber'] = grabber
    modis_data['pticks'] = pticks
    modis_data['ptick_labels'] = ptick_labels
    modis_data['titles']  = []
    
    count = 0
    for fname in file_names:
    
        total_name = fname
        #total_name = data_loc+fname
        print(total_name)
        data = Dataset(total_name,'r')
        
        for xi in range(len(lat_ranges)-1):       
            print(lat_ranges[xi])
            for yj in range(len(lon_ranges)-1):
                lat_indices = np.where((data.variables['lat'][:] >= lat_ranges[xi]) & (data.variables['lat'][:] < lat_ranges[xi+1]))
                lon_indices = np.where((data.variables['lon'][:] >= lon_ranges[yj]) & (data.variables['lon'][:] < lon_ranges[yj+1]))
                modis_data['data'][count,xi,yj] =\
                    np.average(data.variables[grabber][lat_indices[0],lon_indices[0]])
        data.close() 
        modis_data['titles'].append(fname)
        #modis_data['dates'].append(fname[-11:end_string_idx])
        count+=1

    # Convert the longitude values from 0 - 360 to -180 - 180
    return modis_data

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Trend calculating functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def trend_calc(modis_data,x_ind,y_ind,variable,thielsen=False):
    temp_data = modis_data[variable][:,x_ind,y_ind]
    temp_data[temp_data < -999.] = np.nan
    #temp_data = np.ma.masked_where(temp_data < -999., temp_data)
    # Don't calculate trends if there are less than 2 valid values
    if(np.count_nonzero(~np.isnan(temp_data)) < 2):
        trend = pcnt_change = np.nan
    else:
        avg_modis = np.nanmean(temp_data)
        interpx = np.arange(len(temp_data))
        #interpx = years
        ##interper = np.poly1d(np.polyfit(interpx,temp_data,1)) 
        ### Normalize trend by dividing by number of years
        ##trend = (interper(interpx[-1])-interper(interpx[0]))

        slope, intercept, r_value, p_value, std_err = stats.linregress(interpx,temp_data)
        ##print(slope/len(test_dict[dictkey].keys())
        regress_y = interpx*slope+intercept
        trend = regress_y[-1] - regress_y[0]

        if(thielsen==True):
            #The slope
            S=0
            sm=0
            nx = len(temp_data)
            num_d=int(nx*(nx-1)/2)  # The number of elements in avgs
            Sn=np.zeros(num_d)
            for si in range(0,nx-1):
                for sj in range(si+1,nx):
                    # Find the slope between the two points
                    Sn[sm] = (temp_data[si]-temp_data[sj])/(si-sj) 
                    sm=sm+1
                # Endfor
            # Endfor
            Snsorted=sorted(Sn)
            sm=int(num_d/2.)
            if(2*sm    == num_d):
                tsslope=0.5*(Snsorted[sm]+Snsorted[sm+1])
            if(2*sm+1 == num_d): 
                tsslope=Snsorted[sm+1]
            regress_ts = interpx*tsslope+intercept
            trend = regress_ts[-1]-regress_ts[0]


        pcnt_change = (trend/avg_modis)*100.
    # Find the percent change per decade
    return trend,pcnt_change


# modis_trendCalc calculates the trends over the time period at each
# grid point on the 25x25 km grid.
def modis_trendCalc(modis_data,thielSen=False):
    # Loop over the data and calculate trends
    for i in range(448):
        print(i)
        max_trend = -99.
        min_trend = 99.
        for j in range(304):
            modis_data['thick_trends'][i,j],temp_pcnt = \
                trend_calc(modis_data,i,j,'ice_thick',thielsen=thielSen)
            #modis_data['thick_land_trends'][i,j],temp_pcnt = \
            #    trend_calc(modis_data,i,j,'ice_thick',thielsen=thielSen)
            modis_data['con_trends'][i,j],temp_pcnt = \
                trend_calc(modis_data,i,j,'ice_con',thielsen=thielSen)
            #modis_data['con_land_trends'][i,j],temp_pcnt = \
            #    trend_calc(modis_data,i,j,'ice_con',thielsen=thielSen)
            ##temp_trend = modis_data['trends'][i,j]
            ##if(temp_trend>max_trend):
            ##    max_trend = temp_trend
            ##if(temp_trend<min_trend):
            ##    min_trend = temp_trend
        ##print("  max trend = ",max_trend)
        ##print("  min trend = ",min_trend)
        # Deal with land masks
        #good_indices = np.where(modis_data['data'][0,i,:]<251.)
        #land_indices = np.where(modis_data['data'][0,i,:]>=251.)
        #modis_data['trends'][i,land_indices] = np.nan
        #modis_data['land_trends'][i,good_indices] = np.nan

    return modis_data

# Calculate trends on the 1x1 degree lat/lon grid
def modis_gridtrendCalc(modis_data,area=True,thielSen=False):
    if(area==True):
        print("\nCalculating area trends\n")
    else:
        print("\nCalculating % concentration trends\n")
    modis_data['month_fix'] = '_monthfix'

    #lat_ranges = np.arange(minlat,90.5,1.0)
    #lon_ranges = np.arange(0.5,360.5,1.0)
    lat_ranges = modis_data['grid_lat']
    lon_ranges = modis_data['grid_lon']
    ## Create array to hold monthly averages over region
    #initial_avgs = np.full(len(CERES_dict['dates']),-9999.)
    #initial_years  = np.zeros(len(initial_avgs))
    #initial_months = np.zeros(len(initial_avgs))
   
    # Loop over the lat and lon dimensions
    for xi in range(len(lat_ranges)):
        max_trend = -999
        min_trend = 999
        for yj in range(len(lon_ranges)):
            # Calculate the trend at the current box
            if(area==True):
                interp_data = (modis_data['grid_modis_conc'][:,xi,yj]/100.)*modis_data['grid_total_area'][xi,yj]
                good_indices = np.where(np.isnan(interp_data)==False)
            else:
                interp_data = modis_data['grid_modis_conc'][:,xi,yj]
                good_indices = np.where(interp_data!=-99.)
            if(len(good_indices[0])==0):
                total_trend = -9999.
            else:
                interpx = np.arange(len(interp_data))
                #print(CERES_dict['data'][:,xi,yj][good_indices])
                total_interper = np.poly1d(np.polyfit( \
                    interpx[good_indices],\
                    interp_data[good_indices],1)) 
                # Normalize trend by dividing by number of years
                total_trend = (total_interper(interpx[-1])-\
                               total_interper(interpx[0]))
                if(total_trend>max_trend):
                    max_trend = total_trend
                if(total_trend<min_trend):
                    min_trend = total_trend
            if(area==True):
                modis_data['grid_modis_area_trend'][xi,yj] = total_trend
        print(xi)
        #print("max trend = ",max_trend)
        #print("min trend = ",min_trend)

    return modis_data

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Gridding functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# grid_data_conc grids the 25x25 km gridded modis concentration data into
# a 1x1 degree lat/lon grid
def grid_data_values(modis_data):
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(30.,90.,1.0)
    grid_con    = np.full((len(modis_data['ice_con'][:,0,0]),len(lat_ranges),len(lon_ranges)),-99.)
    grid_thick  = np.full((len(modis_data['ice_thick'][:,0,0]),len(lat_ranges),len(lon_ranges)),-99.)
    grid_con_cc = np.full((len(modis_data['ice_con'][:,0,0]),len(lat_ranges),len(lon_ranges)),-999.)
    grid_thick_cc = np.full((len(modis_data['ice_thick'][:,0,0]),len(lat_ranges),len(lon_ranges)),-999.)
    print("Size of grid array: ",grid_con.shape)
    for nt in range(len(modis_data['ice_con'][:,0,0])):
        print(nt)
        for xi in range(448):
            # Don't grid the data if any portion of the lat/lon box is over land.
            # Don't include land data
            for yj in range(304):
                lat_index = np.where(np.floor(modis_data['lat'][xi,yj])>=lat_ranges)[-1][-1]
                lon_index = np.where(np.floor(modis_data['lon'][xi,yj])>=lon_ranges)[-1][-1]
                # Add the current pixel area into the correct grid box, no
                # matter if the current box is missing or not.
                #if((lat_index==20) & (lon_index==10)):
                #    print("Current grid area = ",grid_modis_area[lat_index,lon_index])
                #if((lat_index==20) & (lon_index==10)):
                #    print("    New grid area = ",grid_modis_area[lat_index,lon_index])
                if(modis_data['ice_con'][nt,xi,yj] != -9999.0):
                    #if(nt==0): grid_modis_area[lat_index,lon_index] += modis_data['area'][xi,yj]
                    if(grid_con_cc[nt,lat_index,lon_index]==-999.):
                        grid_con[nt,lat_index,lon_index] = modis_data['ice_con'][nt,xi,yj]
                        grid_con_cc[nt,lat_index,lon_index] = 1.
                    else:
                        grid_con[nt,lat_index,lon_index] = ((grid_con[nt,lat_index,\
                            lon_index]*grid_con_cc[nt,lat_index,lon_index])+\
                            modis_data['ice_con'][nt,xi,yj])/(grid_con_cc[nt,lat_index,\
                            lon_index]+1.)
                        grid_con_cc[nt,lat_index,lon_index]+=1
                if(modis_data['ice_thick'][nt,xi,yj] != -9999.0):
                    #if(nt==0): grid_modis_area[lat_index,lon_index] += modis_data['area'][xi,yj]
                    if(grid_thick_cc[nt,lat_index,lon_index]==-999.):
                        grid_thick[nt,lat_index,lon_index] = modis_data['ice_thick'][nt,xi,yj]
                        grid_thick_cc[nt,lat_index,lon_index] = 1.
                    else:
                        grid_thick[nt,lat_index,lon_index] = ((grid_thick[nt,lat_index,\
                            lon_index]*grid_thick_cc[nt,lat_index,lon_index])+\
                            modis_data['ice_thick'][nt,xi,yj])/(grid_thick_cc[nt,lat_index,\
                            lon_index]+1.)
                        grid_thick_cc[nt,lat_index,lon_index]+=1
                #else:
                #    if(nt==0): grid_modis_area[lat_index,lon_index] = np.nan

                    # end else
                # end if good modis check
            # end y grid loop
        # end x grid loop
    # end time loop 
    modis_data['grid_con']   = grid_con
    modis_data['grid_thick'] = grid_thick
    modis_data['grid_con_cc'] = grid_con_cc
    modis_data['grid_thick_cc'] = grid_thick_cc
    modis_data['grid_lat'] = lat_ranges
    modis_data['grid_lon'] = lon_ranges
    
    return modis_data

# grid_data averages the 25x25 km trends into a 1x1 degree grid
# The 'grid_modis' paramter in the dictionary therefore contains
# the gridded trends.
def grid_data_trends(modis_dict):
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(30.,90.,1.0)
    grid_modis    = np.full((len(lat_ranges),len(lon_ranges)),-999.)
    grid_modis_cc = np.full((len(lat_ranges),len(lon_ranges)),-999.)
    for xi in range(448):
        print(xi)
        for yj in range(304):
            if(not np.isnan(modis_dict['trends'][xi,yj])):
                lat_index = np.where(np.floor(modis_dict['lat'][xi,yj])>=lat_ranges)[-1][-1]
                lon_index = np.where(np.floor(modis_dict['lon'][xi,yj])>=lon_ranges)[-1][-1]
                if(grid_modis_cc[lat_index,lon_index]==-999.):
                    grid_modis[lat_index,lon_index] = modis_dict['trends'][xi,yj]
                    grid_modis_cc[lat_index,lon_index] = 1.
                else:
                    grid_modis[lat_index,lon_index] = ((grid_modis[lat_index,lon_index]*grid_modis_cc[lat_index,lon_index])+\
                                                     modis_dict['trends'][xi,yj])/(grid_modis_cc[lat_index,lon_index]+1.)
                    grid_modis_cc[lat_index,lon_index]+=1
    modis_dict['grid_modis'] = grid_modis
    modis_dict['grid_modis_cc'] = grid_modis_cc
    modis_dict['grid_lat'] = lat_ranges
    modis_dict['grid_lon'] = lon_ranges
    return modis_dict

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Plotting functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def read_true_color(date_str,composite = False):
    # Determine the correct MODIS file associated with the date
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

    return var, crs, lat_lims, lon_lims

def plot_true_color_satpy(date_str,ax = None,zoom=True,save=False,composite=False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    var, crs, lat_lims, lon_lims = read_true_color(date_str,composite=composite)

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

    ax.set_title('Aqua MODIS\n'+dt_date_str.strftime('%Y-%m-%d %H:%M'))

    if(save):
        outname = 'modis_true_color_' + date_str + zoom_add + cmpst_add + '.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

    """

    # Extract the true-color data
    # ---------------------------
    scn.load(['true_color'])

    scn.save_dataset('true_color','test_image_true'+date_str+'.png')
    scn.show('true_color')

    # Overlay the imagery on cartopy
    # ------------------------------
    ##my_area

    # Save the image
    # --------------
    """

def plot_true_color(filename,zoom=True):
    modis = SD.SD(filename)

    dat = modis.attributes().get('CoreMetadata.0').split()
    indx = dat.index('EQUATORCROSSINGDATE')+9
    cross_date = dat[indx][1:len(dat[indx])-1]
    cross_time = filename.strip().split('/')[-1].split('.')[2]

    lat5 = modis.select('Latitude').get()
    lon5 = modis.select('Longitude').get()

    red   = modis.select('EV_250_Aggr1km_RefSB').get()[0]
    green = modis.select('EV_500_Aggr1km_RefSB').get()[1]
    blue  = modis.select('EV_500_Aggr1km_RefSB').get()[0]

    red   = red[::5,::5]
    green = green[::5,::5]
    blue  = blue[::5,::5]

    red_scale    = modis.select('EV_250_Aggr1km_RefSB').attributes().get('reflectance_scales')[0]
    red_offset   = modis.select('EV_250_Aggr1km_RefSB').attributes().get('reflectance_offsets')[0]
    green_scale  = modis.select('EV_500_Aggr1km_RefSB').attributes().get('reflectance_scales')[1]
    green_offset = modis.select('EV_500_Aggr1km_RefSB').attributes().get('reflectance_offsets')[1]
    blue_scale   = modis.select('EV_500_Aggr1km_RefSB').attributes().get('reflectance_scales')[0]
    blue_offset  = modis.select('EV_500_Aggr1km_RefSB').attributes().get('reflectance_offsets')[0]

    modis.end()

    red   = (red - red_offset) * red_scale
    green = (green - green_offset) * green_scale
    blue  = (blue - blue_offset) * blue_scale

    red   = red*(255./1.1) 
    green = green*(255./1.1) 
    blue  = blue*(255./1.1) 

    old_val = np.array([0,30,60,120,190,255])
    ehn_val = np.array([0,110,160,210,240,255])
    red     = np.interp(red,old_val,ehn_val) / 255.
    old_val = np.array([0,30,60,120,190,255])
    ehn_val = np.array([0,110,160,200,230,240])
    green   = np.interp(green,old_val,ehn_val) / 255.
    old_val = np.array([0,30,60,120,190,255])
    ehn_val = np.array([0,100,150,210,240,255])
    blue    = np.interp(blue,old_val,ehn_val) / 255.
  
    image = np.zeros((red.shape[0],red.shape[1],3))
    image[:,:,0] = red
    image[:,:,1] = green
    image[:,:,2] = blue
   
    colortuple = tuple(np.array([image[:,:,0].flatten(), \
        image[:,:,1].flatten(), image[:,:,2].flatten()]).transpose().tolist())

    cornerLats = getCorners(lat5) ; cornerLons = getCorners(lon5)
    
    plt.close('all') 
    fig1 = plt.figure()
    ax = plt.axes(projection = mapcrs)

    image = np.ma.masked_where(np.isnan(image),image)

    print(lon5.shape, lat5.shape, image.shape)

    ax.pcolormesh(cornerLons,cornerLats,image[:,:,0],color= colortuple, shading='auto', \
        transform = ccrs.PlateCarree()) 
    
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.STATES)
    ax.coastlines()
    if(zoom):
        ax.set_extent([plot_limits_dict[cross_date][cross_time]['Lon'][0], \
                       plot_limits_dict[cross_date][cross_time]['Lon'][1], \
                       plot_limits_dict[cross_date][cross_time]['Lat'][0], \
                       plot_limits_dict[cross_date][cross_time]['Lat'][1]], \
            ccrs.PlateCarree())

    plt.show()

# dt_date_str is of format YYYYMMDDHHMM
def read_MODIS_channel(date_str, channel, zoom = False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # Extract the filename given the channel
    # --------------------------------------
    if(str(channel)[:2] == 'wv'):
        filename = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['mdswv']
    else:
        filename = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis']
    
    print("Reading MODIS channel",channel," from ",filename)

    MODIS_data = {}

    modis = SD.SD(filename)

    dat = modis.attributes().get('CoreMetadata.0').split()
    indx = dat.index('EQUATORCROSSINGDATE')+9
    cross_date = dat[indx][1:len(dat[indx])-1]
    indx = dat.index('EQUATORCROSSINGTIME')+9
    cross_time = dat[indx][1:len(dat[indx])-1]

    print('MODIS orbit info',cross_date,cross_time)

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

    MODIS_data['data'] = data
    MODIS_data['lat']  = lat5
    MODIS_data['lon']  = lon5
    MODIS_data['variable']  = label
    MODIS_data['cross_date']  = cross_date
    MODIS_data['channel']  = channel
    MODIS_data['colors']  = colors
    MODIS_data['file_time'] = filename.strip().split('/')[-1].split('.')[2]

    if(zoom):
        # Mask MODIS_data['data'] that are outside the desired range
        # --------------------------------------------
        MODIS_data['data'][(((MODIS_data['lat'] < plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][0]) | \
                             (MODIS_data['lat'] > plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][1])) | \
                            ((MODIS_data['lon'] < plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][0]) | \
                             (MODIS_data['lon'] > plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][1])))] = -999.
        MODIS_data['lat'][ (((MODIS_data['lat'] < plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][0]) | \
                             (MODIS_data['lat'] > plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][1])) | \
                            ((MODIS_data['lon'] < plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][0]) | \
                             (MODIS_data['lon'] > plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][1])))] = -999.
        MODIS_data['lon'][ (((MODIS_data['lat'] < plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][0]) | \
                             (MODIS_data['lat'] > plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][1])) | \
                            ((MODIS_data['lon'] < plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][0]) | \
                             (MODIS_data['lon'] > plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][1])))] = -999.

        MODIS_data['data'] = np.ma.masked_where(MODIS_data['data'] == -999., MODIS_data['data'])
        MODIS_data['lat'] = np.ma.masked_where(MODIS_data['lat'] == -999., MODIS_data['lat'])
        MODIS_data['lon'] = np.ma.masked_where(MODIS_data['lon'] == -999., MODIS_data['lon'])


    return MODIS_data

def plot_MODIS_channel(date_str,channel,zoom=True,show_smoke=False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    filename = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis']

    if(channel == 'red'):
        channel = 1
    elif(channel == 'green'):
        channel = 4
    elif(channel == 'blue'):
        channel = 3

    # Call read_MODIS_channel to read the desired MODIS data from the
    # file and put it in a dictionary
    # ---------------------------------------------------------------
    MODIS_data = read_MODIS_channel(date_str, channel)

    print("Data max = ",np.max(MODIS_data['data']), "  Data min = ",np.min(MODIS_data['data']))

    plt.close('all')
    fig1 = plt.figure()
    datacrs = ccrs.PlateCarree()
    #mapcrs = ccrs.Miller()
    ax = plt.axes(projection = mapcrs)

    plot_MODIS_spatial(MODIS_data, ax, zoom)
#def plot_MODIS_spatial(MODIS_data, pax, zoom):

    #mesh = ax.pcolormesh(MODIS_data['lon'],MODIS_data['lat'],\
    #    MODIS_data['data'],cmap = MODIS_data['colors'], shading='auto', \
    #    transform = datacrs) 

    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(date_str) 
        #hash_data1, nohash_data1 = find_plume(filename) 
        hash0 = ax.pcolor(MODIS_data['lon'],MODIS_data['lat'],\
            hash_data1, hatch = '///', alpha=0., transform = datacrs)

    ##!#cbar = plt.colorbar(mesh,orientation='horizontal',pad=0,\
    ##!#    aspect=50,shrink = 0.850)
    ##!#cbar.ax.tick_params(labelsize=14)
    ##!#cbar.set_label(MODIS_data['variable'],fontsize=16,weight='bold')
    ##!#
    ##!#ax.add_feature(cfeature.BORDERS)
    ##!#ax.add_feature(cfeature.STATES)
    ##!#ax.coastlines()
    ##!#if(zoom):
    ##!#    ax.set_extent([plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][0], \
    ##!#                   plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][1], \
    ##!#                   plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][0], \
    ##!#                   plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][1]],\
    ##!#                   ccrs.PlateCarree())
    ##!#ax.set_title('Channel ' + str(channel) + '\n' + \
    ##!#    channel_dict[str(channel)]['Bandwidth_label']) 

    plt.show()

def read_OMI_match_MODIS(date_str):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    modis_date = dt_date_str.strftime('%Y-%m-%d')

    data = h5py.File(plot_limits_dict[modis_date][date_str[8:]]['omi'],'r')
    LAT   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude'][:,:]
    LON   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude'][:,:]
    UVAI  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex'][:,:]
    XTRACK = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags'][:,:]
    mask_UVAI = np.ma.masked_where((XTRACK < -2e5) | (UVAI < -2e5), UVAI)
    mask_UVAI = np.ma.masked_where((\
        ((LAT < plot_limits_dict[modis_date][date_str[8:]]['Lat'][0]) | \
         (LAT > plot_limits_dict[modis_date][date_str[8:]]['Lat'][1])) | \
        ((LON < plot_limits_dict[modis_date][date_str[8:]]['Lon'][0]) | \
         (LON > plot_limits_dict[modis_date][date_str[8:]]['Lon'][1]))), mask_UVAI)
    data.close()

    return LAT, LON, mask_UVAI

def read_CERES_match_MODIS(date_str):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # NOTE: need to use the time variable to screen out any data
    # that are not close to the event time.
    ##!#base_date = datetime(year=1970,month=1,day=1)
    ##!#start_date = dt_date_str - timedelta(hours = 1)
    ##!#end_date   = dt_date_str + timedelta(hours = 2)

    ##!#print(start_date, end_date)

    modis_date = dt_date_str.strftime('%Y-%m-%d')

    data = Dataset(plot_limits_dict[modis_date][date_str[8:]]['ceres'],'r')
    LAT   = 90. - data.variables['Colatitude_of_CERES_FOV_at_surface'][:]
    LON   = data.variables['Longitude_of_CERES_FOV_at_surface'][:]
    LON[LON>179.99] = -360.+LON[LON>179.99]
    swflux  = data.variables['CERES_SW_TOA_flux___upwards'][:]
    lwflux  = data.variables['CERES_LW_TOA_flux___upwards'][:]
    time  = data.variables['time'][:]
    #local_time = np.array([base_date + relativedelta(days = ttime) for ttime in time])

    data.close()

    mask_LAT = LAT[ \
        (LAT >= plot_limits_dict[modis_date][date_str[8:]]['Lat'][0]) & \
        (LAT <= plot_limits_dict[modis_date][date_str[8:]]['Lat'][1]) & \
        (LON >= plot_limits_dict[modis_date][date_str[8:]]['Lon'][0]) & \
        (LON <= plot_limits_dict[modis_date][date_str[8:]]['Lon'][1]) & \
        (swflux > 0) & (lwflux > 0)]
    mask_LON = LON[ \
        (LAT >= plot_limits_dict[modis_date][date_str[8:]]['Lat'][0]) & \
        (LAT <= plot_limits_dict[modis_date][date_str[8:]]['Lat'][1]) & \
        (LON >= plot_limits_dict[modis_date][date_str[8:]]['Lon'][0]) & \
        (LON <= plot_limits_dict[modis_date][date_str[8:]]['Lon'][1]) & \
        (swflux > 0) & (lwflux > 0)]
    mask_swf = swflux[ \
        (LAT >= plot_limits_dict[modis_date][date_str[8:]]['Lat'][0]) & \
        (LAT <= plot_limits_dict[modis_date][date_str[8:]]['Lat'][1]) & \
        (LON >= plot_limits_dict[modis_date][date_str[8:]]['Lon'][0]) & \
        (LON <= plot_limits_dict[modis_date][date_str[8:]]['Lon'][1]) & \
        (swflux > 0) & (lwflux > 0)]
    mask_lwf = lwflux[ \
        (LAT >= plot_limits_dict[modis_date][date_str[8:]]['Lat'][0]) & \
        (LAT <= plot_limits_dict[modis_date][date_str[8:]]['Lat'][1]) & \
        (LON >= plot_limits_dict[modis_date][date_str[8:]]['Lon'][0]) & \
        (LON <= plot_limits_dict[modis_date][date_str[8:]]['Lon'][1]) & \
        (swflux > 0) & (lwflux > 0)]
    ##!#mask_time = local_time[ \
    ##!#    (LAT >= plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][0]) & \
    ##!#    (LAT <= plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][1]) & \
    ##!#    (LON >= plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][0]) & \
    ##!#    (LON <= plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][1]) & \
    ##!#    (swflux > 0) & (lwflux > 0)]

    return mask_LAT, mask_LON, mask_swf, mask_lwf

# date_str of format "YYYYMMDDHHMM"
def compare_MODIS_3panel(date_str,channel1,channel2,channel3,zoom=True,save=False,\
        plot_ASOS_loc = False, show_smoke = True, compare_OMI = False, \
        compare_CERES = False, return_MODIS = False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    filename = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis']

    if(channel1== 'red'):
        channel1= 1
    elif(channel1== 'green'):
        channel1= 4
    elif(channel1== 'blue'):
        channel1= 3

    # Step 1: Call read_MODIS_channel to read the desired MODIS data from the
    # file and put it in a dictionary
    # -----------------------------------------------------------------------
    MODIS_data1 = read_MODIS_channel(dt_date_str.strftime("%Y%m%d%H%M"), channel1, zoom = True)
    MODIS_data2 = read_MODIS_channel(dt_date_str.strftime("%Y%m%d%H%M"), channel2, zoom = True)
    MODIS_data3 = read_MODIS_channel(dt_date_str.strftime("%Y%m%d%H%M"), channel3, zoom = True)

    ##!#print("Data1 max = ",np.max(MODIS_data1['data']), "  Data1 min = ",np.min(MODIS_data1['data']))
    ##!#print("Data2 max = ",np.max(MODIS_data2['data']), "  Data2 min = ",np.min(MODIS_data2['data']))
    ##!#print("Data3 max = ",np.max(MODIS_data3['data']), "  Data3 min = ",np.min(MODIS_data3['data']))

    # Create copies to set the minimum and maximums for each plot
    # -----------------------------------------------------------
    cpy_1 = np.copy(MODIS_data1['data'])
    cpy_2 = np.copy(MODIS_data2['data'])
    cpy_3 = np.copy(MODIS_data3['data'])

    cpy_1 = np.ma.masked_where((((MODIS_data1['lat'] < plot_limits_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lat'][0]) | \
                         (MODIS_data1['lat'] > plot_limits_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lat'][1])) | \
                        ((MODIS_data1['lon'] < plot_limits_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lon'][0]) | \
                         (MODIS_data1['lon'] > plot_limits_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lon'][1]))), cpy_1)
    cpy_2 = np.ma.masked_where((((MODIS_data2['lat'] < plot_limits_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lat'][0]) | \
                         (MODIS_data2['lat'] > plot_limits_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lat'][1])) | \
                        ((MODIS_data2['lon'] < plot_limits_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lon'][0]) | \
                         (MODIS_data2['lon'] > plot_limits_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lon'][1]))), cpy_2)
    cpy_3 = np.ma.masked_where((((MODIS_data3['lat'] < plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) | \
                         (MODIS_data3['lat'] > plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1])) | \
                        ((MODIS_data3['lon'] < plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) | \
                         (MODIS_data3['lon'] > plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]))), cpy_3)


    # Step 2: Set up figure to have 3 panels
    # --------------------------------------
    datacrs = ccrs.PlateCarree() 
    #mapcrs = ccrs.LambertConformal()

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

    # Plot channels 1, 2, and 3
    plot_MODIS_spatial(MODIS_data1, ax0, zoom = zoom)
    plot_MODIS_spatial(MODIS_data2, ax1, zoom = zoom)
    plot_MODIS_spatial(MODIS_data3, ax2, zoom = zoom)

    ##!## Step 3: Plot the MODIS channel data in the first 2 panels
    ##!## ---------------------------------------------------------
    ##!#mesh0 = ax0.pcolormesh(MODIS_data1['lon'],MODIS_data1['lat'],\
    ##!#    cpy_1,cmap = MODIS_data1['colors'], shading='auto', \
    ##!#    vmin = np.nanmin(MODIS_data1['data']), \
    ##!#    vmax = np.nanmax(MODIS_data1['data']), transform = datacrs) 

    ##!#cbar0 = plt.colorbar(mesh0,ax=ax0,orientation='vertical',\
    ##!#    pad=0.03,label=MODIS_data1['variable'])

    ##!#ax0.add_feature(cfeature.BORDERS)
    ##!#ax0.add_feature(cfeature.STATES)
    ##!#ax0.coastlines()
    ##!#if(zoom):
    ##!#    ax0.set_extent([plot_limits_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lon'][0], \
    ##!#                    plot_limits_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lon'][1], \
    ##!#                    plot_limits_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lat'][0], \
    ##!#                    plot_limits_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lat'][1]],\
    ##!#                    datacrs)
    ##!##ax0.set_title('MODIS Ch. ' + str(channel1) + '\n' + \
    ##!##    str(channel_dict[str(channel1)]['Bandwidth'][0]) + ' μm - ' + \
    ##!##    str(channel_dict[str(channel1)]['Bandwidth'][1]) + ' μm')
    ##!#ax0.set_title('Channel ' + str(channel1) + '\n' + \
    ##!#    channel_dict[str(channel1)]['Bandwidth_label']) 


    ##!## Plot channel 2
    ##!#mesh1 = ax1.pcolormesh(MODIS_data2['lon'],MODIS_data2['lat'],\
    ##!#    MODIS_data2['data'],cmap = MODIS_data2['colors'], shading='auto', \
    ##!#    vmin = np.nanmin(MODIS_data2['data']), \
    ##!#    vmax = np.nanmax(MODIS_data2['data']), transform = datacrs) 


    ##!#cbar1 = plt.colorbar(mesh1,ax=ax1,orientation='vertical',\
    ##!#    pad=0.03,label=MODIS_data2['variable'])
    ##!#
    ##!#ax1.add_feature(cfeature.BORDERS)
    ##!#ax1.add_feature(cfeature.STATES)
    ##!#ax1.coastlines()
    ##!#if(zoom):
    ##!#    ax1.set_extent([plot_limits_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lon'][0], \
    ##!#                    plot_limits_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lon'][1], \
    ##!#                    plot_limits_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lat'][0], \
    ##!#                    plot_limits_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lat'][1]],\
    ##!#                    datacrs)
    ##!##ax1.set_title('MODIS Ch. ' + str(channel2) + '\n' + \
    ##!##    str(channel_dict[str(channel2)]['Bandwidth'][0]) + ' μm - ' + \
    ##!##    str(channel_dict[str(channel2)]['Bandwidth'][1]) + ' μm')
    ##!#ax1.set_title('Channel ' + str(channel2) + '\n' + \
    ##!#    channel_dict[str(channel2)]['Bandwidth_label']) 

    ##!## Plot channel 3
    ##!#mesh2 = ax2.pcolormesh(MODIS_data3['lon'],MODIS_data3['lat'],\
    ##!#    MODIS_data3['data'],cmap = MODIS_data3['colors'], shading='auto', \
    ##!#    vmin = np.nanmin(MODIS_data3['data']), \
    ##!#    vmax = np.nanmax(MODIS_data3['data']), transform = datacrs) 
    ##!#cbar2 = plt.colorbar(mesh2,ax=ax2,orientation='vertical',\
    ##!#    pad=0.03,label=MODIS_data3['variable'])
    ##!#
    ##!#ax2.add_feature(cfeature.BORDERS)
    ##!#ax2.add_feature(cfeature.STATES)
    ##!#ax2.coastlines()
    ##!#if(zoom):
    ##!#    ax2.set_extent([plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0], \
    ##!#                    plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1], \
    ##!#                    plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0], \
    ##!#                    plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1]],\
    ##!#                    datacrs)
    ##!##ax2.set_title('MODIS Ch. ' + str(channel3) + '\n' + \
    ##!##    str(channel_dict[str(channel3)]['Bandwidth'][0]) + ' μm - ' + \
    ##!##    str(channel_dict[str(channel3)]['Bandwidth'][1]) + ' μm')
    ##!#ax2.set_title('Channel ' + str(channel3) + '\n' + \
    ##!#    channel_dict[str(channel3)]['Bandwidth_label']) 


    if(compare_OMI):
        print("Reading OMI data")
        LAT, LON, mask_UVAI = read_OMI_match_MODIS(date_str)
        ##!#data = h5py.File(plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['omi'],'r')
        ##!#LAT   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude'][:,:]
        ##!#LON   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude'][:,:]
        ##!#UVAI  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex'][:,:]
        ##!#XTRACK = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags'][:,:]
        ##!#mask_UVAI = np.ma.masked_where((XTRACK < -2e5) | (UVAI < -2e5), UVAI)
        ##!#mask_UVAI = np.ma.masked_where((((LAT < plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) | \
        ##!#                     (LAT > plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1])) | \
        ##!#                    ((LON < plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) | \
        ##!#                     (LON > plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]))), mask_UVAI)

        plot_OMI_spatial(date_str, LAT, LON, mask_UVAI, axo, zoom = zoom)

        ##!#mesh3 = axo.pcolormesh(LON,LAT,mask_UVAI, cmap = 'plasma', shading='auto', \
        ##!#    vmin = np.nanmin(mask_UVAI), vmax = np.nanmax(mask_UVAI), transform = datacrs) 
        ##!#axo.add_feature(cfeature.BORDERS)
        ##!#axo.add_feature(cfeature.STATES)
        ##!#axo.coastlines()
        ##!#if(zoom):
        ##!#    axo.set_extent([plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0], \
        ##!#                    plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1], \
        ##!#                    plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0], \
        ##!#                    plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1]],\
        ##!#                    datacrs)
        ##!#cbar3 = plt.colorbar(mesh3,ax=axo,orientation='vertical',\
        ##!#    pad=0.03,label='OMI UVAI')
        ##!#axo.set_title('OMI UVAI')

    if(compare_CERES):
        print("Reading CERES data")

        mask_LAT, mask_LON, mask_swf, mask_lwf = read_CERES_match_MODIS(date_str)

        ##!#base_date = datetime(year=1970,month=1,day=1)
        ##!#start_date = dt_date_str - timedelta(hours = 1)
        ##!#end_date   = dt_date_str + timedelta(hours = 2)

        ##!#print(start_date, end_date)

        ##!#data = Dataset(plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['ceres'],'r')
        ##!#LAT   = 90. - data.variables['Colatitude_of_CERES_FOV_at_surface'][:]
        ##!#LON   = data.variables['Longitude_of_CERES_FOV_at_surface'][:]
        ##!#LON[LON>179.99] = -360.+LON[LON>179.99]
        ##!#swflux  = data.variables['CERES_SW_TOA_flux___upwards'][:]
        ##!#lwflux  = data.variables['CERES_LW_TOA_flux___upwards'][:]
        ##!#time  = data.variables['time'][:]
        ##!#local_time = np.array([base_date + relativedelta(days = ttime) for ttime in time])

        ##!#mask_LAT = LAT[ \
        ##!#    (LAT >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) & \
        ##!#    (LAT <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1]) & \
        ##!#    (LON >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) & \
        ##!#    (LON <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]) & \
        ##!#    (swflux > 0) & (lwflux > 0)]
        ##!#mask_LON = LON[ \
        ##!#    (LAT >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) & \
        ##!#    (LAT <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1]) & \
        ##!#    (LON >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) & \
        ##!#    (LON <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]) & \
        ##!#    (swflux > 0) & (lwflux > 0)]
        ##!#mask_swf = swflux[ \
        ##!#    (LAT >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) & \
        ##!#    (LAT <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1]) & \
        ##!#    (LON >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) & \
        ##!#    (LON <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]) & \
        ##!#    (swflux > 0) & (lwflux > 0)]
        ##!#mask_lwf = lwflux[ \
        ##!#    (LAT >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) & \
        ##!#    (LAT <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1]) & \
        ##!#    (LON >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) & \
        ##!#    (LON <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]) & \
        ##!#    (swflux > 0) & (lwflux > 0)]
        ##!###!#mask_time = local_time[ \
        ##!###!#    (LAT >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) & \
        ##!###!#    (LAT <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1]) & \
        ##!###!#    (LON >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) & \
        ##!###!#    (LON <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]) & \
        ##!###!#    (swflux > 0) & (lwflux > 0)]

        # Removed masked data
        ##!#triMesh = Triangulation(mask_LON,mask_LAT)

        ##!## Reshape the data to make it work with pcolormesh
        ##!#first_size = int(np.sqrt(mask_swf.shape)) 
        ##!#max_size = int(first_size ** 2.)
        ##!#mask_swf = mask_swf[:max_size].reshape((first_size, first_size))
        ##!#mask_lwf = mask_lwf[:max_size].reshape((first_size, first_size))
        ##!#LAT = LAT[:max_size].reshape((first_size, first_size))
        ##!#LON = LON[:max_size].reshape((first_size, first_size))
      
        plot_CERES_spatial(date_str,mask_LAT, mask_LON, mask_swf, 'SWF', axcs, zoom = zoom)
        plot_CERES_spatial(date_str,mask_LAT, mask_LON, mask_lwf, 'LWF', axcl, zoom = zoom)

        ##!#print(np.nanmax(mask_LAT.compressed()), np.min(mask_LAT.compressed()))
        ##!#print(np.nanmax(LAT), np.nanmin(LAT))
        ##!#print(plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'])
 
        ##!##scat3 = axcs.scatter(mask_LON, mask_LAT,mask_swf, transform = datacrs)
        ##!#mesh3 = axcs.scatter(mask_LON.compressed(), mask_LAT.compressed(),\
        ##!#    s = 120,marker='s',c = mask_swf.compressed(),cmap='plasma', \
        ##!#    transform = datacrs)
        ##!##mesh3 = axcs.tricontourf(triMesh,mask_swf, cmap = 'plasma', shading='auto', \
        ##!##    vmin = np.nanmin(mask_swf), vmax = np.nanmax(mask_swf), transform = datacrs) 
#       ##!# mesh3 = axcs.pcolormesh(LON,LAT,mask_swf, cmap = 'plasma', shading='auto', \
#       ##!#     vmin = np.nanmin(mask_swf), vmax = np.nanmax(mask_swf), transform = datacrs) 
        ##!#axcs.add_feature(cfeature.BORDERS)
        ##!#axcs.add_feature(cfeature.STATES)
        ##!#axcs.coastlines()
        ##!#if(zoom):
        ##!#    axcs.set_extent([plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0], \
        ##!#                     plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1], \
        ##!#                     plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0], \
        ##!#                     plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1]],\
        ##!#                     datacrs)
        ##!#cbar3 = plt.colorbar(mesh3,ax=axcs,orientation='vertical',\
        ##!#    pad=0.03,label='TOA SWF [W/m2]')
        ##!#axcs.set_title('CERES SWF')



        ##!#print(mask_LON.compressed().shape, mask_LAT.compressed().shape, mask_lwf.compressed().shape)
        ##!#mesh4 = axcl.scatter(mask_LON, mask_LAT,\
        ##!#    s = 120,marker = 's',c = mask_lwf,cmap='plasma', \
        ##!#    transform = datacrs)
        ##!##mesh4 = axcl.tricontourf(triMesh,mask_lwf, cmap = 'plasma', shading='auto', \
        ##!##    vmin = np.nanmin(mask_lwf), vmax = np.nanmax(mask_lwf), transform = datacrs) 
        ##!###!#mesh4 = axcl.tricontourf(LON,LAT,mask_lwf, cmap = 'plasma', shading='auto', \
        ##!###!#    vmin = np.nanmin(mask_lwf), vmax = np.nanmax(mask_lwf), transform = datacrs) 
        ##!#axcl.add_feature(cfeature.BORDERS)
        ##!#axcl.add_feature(cfeature.STATES)
        ##!#axcl.coastlines()
        ##!#if(zoom):
        ##!#    axcl.set_extent([plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0], \
        ##!#                    plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1], \
        ##!#                    plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0], \
        ##!#                    plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1]],\
        ##!#                    datacrs)
        ##!#    cbar4 = plt.colorbar(mesh4,ax=axcl,orientation='vertical',\
        ##!#        pad=0.03,label='TOA LWF [W/m2]')
        ##!#axcl.set_title('CERES LWF')

        ##!#data.close()

    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(dt_date_str.strftime('%Y%m%d%H%M')) 
        hash0 = ax0.pcolor(MODIS_data1['lon'],MODIS_data1['lat'],\
            hash_data1[:MODIS_data1['lat'].shape[0], :MODIS_data1['lat'].shape[1]], hatch = '///', alpha=0., transform = datacrs)
        hash1 = ax1.pcolor(MODIS_data2['lon'],MODIS_data2['lat'],\
            hash_data1[:MODIS_data2['lat'].shape[0], :MODIS_data2['lat'].shape[1]], hatch = '///', alpha=0., transform = datacrs)
        hash2 = ax2.pcolor(MODIS_data3['lon'],MODIS_data3['lat'],\
            hash_data1[:MODIS_data3['lat'].shape[0], :MODIS_data3['lat'].shape[1]], hatch = '///', alpha=0., transform = datacrs)
        if(compare_OMI): 
            hash3 = axo.pcolor(MODIS_data3['lon'],MODIS_data3['lat'],\
                hash_data1[:MODIS_data3['lat'].shape[0], :MODIS_data3['lat'].shape[1]], hatch = '///', alpha=0., transform = datacrs)
        if(compare_CERES): 
            hash4 = axcs.pcolor(MODIS_data3['lon'],MODIS_data3['lat'],\
                hash_data1[:MODIS_data3['lat'].shape[0], :MODIS_data3['lat'].shape[1]], hatch = '///', alpha=0., transform = datacrs)
            hash5 = axcl.pcolor(MODIS_data3['lon'],MODIS_data3['lat'],\
                hash_data1[:MODIS_data3['lat'].shape[0], :MODIS_data3['lat'].shape[1]], hatch = '///', alpha=0., transform = datacrs)
    

    if(plot_ASOS_loc):
        print("Plotting ASOS site")
        plot_ASOS_locs(ax0,MODIS_data1,color='lime')
        plot_ASOS_locs(ax1,MODIS_data1,color='lime')
        plot_ASOS_locs(ax2,MODIS_data1,color='lime')
        if(compare_OMI): plot_ASOS_locs(axo,MODIS_data1,color='black')
        if(compare_CERES): 
            plot_ASOS_locs(axcs,MODIS_data1,color='black')
            plot_ASOS_locs(axcl,MODIS_data1,color='black')

    cross_date = MODIS_data1['cross_date']
    file_time  = MODIS_data1['file_time']
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

    if(return_MODIS):
        return MODIS_data1
    

def compare_MODIS_channels(date_str,channel1,channel2,zoom=True,save=False,\
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
    # Step 1: Call read_MODIS_channel to read the desired MODIS data from the
    # file and put it in a dictionary
    # -----------------------------------------------------------------------
    MODIS_data1 = read_MODIS_channel(dt_date_str.strftime('%Y%m%d%H%M'), channel1, zoom = zoom)
    MODIS_data2 = read_MODIS_channel(dt_date_str.strftime('%Y%m%d%H%M'), channel2, zoom = zoom)

    print("Data1 max = ",np.max(MODIS_data1['data']), "  Data1 min = ",np.min(MODIS_data1['data']))
    print("Data2 max = ",np.max(MODIS_data2['data']), "  Data2 min = ",np.min(MODIS_data2['data']))

    print(MODIS_data1['data'].shape, MODIS_data2['data'].shape)

    # --------------------------------------
    # Step 2: Set up figure to have 3 panels
    # --------------------------------------
    datacrs = ccrs.PlateCarree() 
    #mapcrs = ccrs.LambertConformal()

    plt.close('all')
    fig = plt.figure(figsize=(14.5,5))
    ax0 = fig.add_subplot(1,3,2,projection = mapcrs)
    ax1 = fig.add_subplot(1,3,3,projection = mapcrs)
    ax2 = fig.add_subplot(1,3,1)

    # Step 3: Plot the MODIS channel data in the first 2 panels
    # ---------------------------------------------------------
    mesh0 = ax0.pcolormesh(MODIS_data1['lon'],MODIS_data1['lat'],\
        MODIS_data1['data'],cmap = MODIS_data1['colors'], shading='auto', \
        transform = datacrs)

    cbar0 = plt.colorbar(mesh0,ax=ax0,orientation='vertical',\
        pad=0.03,label=MODIS_data1['variable'])
    #cbar0 = plt.colorbar(mesh0,cax = ax0, orientation='horizontal',pad=0,\
    #    aspect=50,shrink = 0.850)
    #cbar0.ax.tick_params(labelsize=14)
    #cbar0.set_label(MODIS_data1['variable'],fontsize=16,weight='bold')
   
    ax0.add_feature(cfeature.BORDERS)
    ax0.add_feature(cfeature.STATES)
    ax0.coastlines()
    if(zoom):
        ax0.set_extent([plot_limits_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lon'][0], \
                        plot_limits_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lon'][1], \
                        plot_limits_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lat'][0], \
                        plot_limits_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lat'][1]],\
                        datacrs)
    #ax0.set_title('MODIS Ch. ' + str(channel1) + '\n' + \
    #    str(channel_dict[str(channel1)]['Bandwidth'][0]) + ' μm - ' + \
    #    str(channel_dict[str(channel1)]['Bandwidth'][1]) + ' μm')
    ax0.set_title('Channel ' + str(channel1) + '\n' + \
        channel_dict[str(channel1)]['Bandwidth_label']) 

    # Plot channel 2
    mesh1 = ax1.pcolormesh(MODIS_data2['lon'],MODIS_data2['lat'],\
        MODIS_data2['data'],cmap = MODIS_data2['colors'], shading='auto', \
        transform = datacrs) 

    cbar1 = plt.colorbar(mesh1,ax=ax1,orientation='vertical',\
        pad=0.03,label=MODIS_data2['variable'])
    #cbar1 = plt.colorbar(mesh1,cax = ax1,orientation='horizontal',pad=1,\
    #    aspect=51,shrink = 1.851)
    #cbar1.ax.tick_params(labelsize=14)
    #cbar1.set_label(MODIS_data2['variable'],fontsize=16,weight='bold')
    
    ax1.add_feature(cfeature.BORDERS)
    ax1.add_feature(cfeature.STATES)
    ax1.coastlines()
    if(zoom):
        ax1.set_extent([plot_limits_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lon'][0], \
                        plot_limits_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lon'][1], \
                        plot_limits_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lat'][0], \
                        plot_limits_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lat'][1]],\
                        datacrs)
    #ax1.set_title('MODIS Ch. ' + str(channel2) + '\n' + \
    #    str(channel_dict[str(channel2)]['Bandwidth'][0]) + ' μm - ' + \
    #    str(channel_dict[str(channel2)]['Bandwidth'][1]) + ' μm')
    ax1.set_title('Channel ' + str(channel2) + '\n' + \
        channel_dict[str(channel2)]['Bandwidth_label']) 

    ##!## Shade areas that are inside the plume, defined currently as areas
    ##!## with reflectance below the mean
    ##!## -----------------------------------------------------------------
    ##!#if(MODIS_data1['variable'] == 'Reflectance'):
    ##!#    hash_data1 = np.ma.masked_where(MODIS_data1['data'] < \
    ##!#        np.nanmean(MODIS_data1['data']),MODIS_data1['data'])
    ##!#    hash_data2 = np.ma.masked_where(MODIS_data1['data'] < \
    ##!#        np.nanmean(MODIS_data1['data']),MODIS_data2['data'])

    ##!#    nohash_data1 = np.ma.masked_where(MODIS_data1['data'] >= \
    ##!#        np.nanmean(MODIS_data1['data']),MODIS_data1['data'])
    ##!#    nohash_data2 = np.ma.masked_where(MODIS_data1['data'] >= \
    ##!#        np.nanmean(MODIS_data1['data']),MODIS_data2['data'])

    ##!#    hash1 = ax0.pcolor(MODIS_data1['lon'],MODIS_data1['lat'],\
    ##!#        hash_data1, hatch = '///', alpha=0., transform = datacrs),
    ##!#    hash2 = ax1.pcolor(MODIS_data2['lon'],MODIS_data2['lat'],\
    ##!#        hash_data2, hatch = '///', alpha=0., transform = datacrs),
    ##!#if(MODIS_data2['variable'] == 'Reflectance'):
    ##!#    hash_data2 = np.ma.masked_where(MODIS_data2['data'] < \
    ##!#        np.nanmean(MODIS_data2['data']),MODIS_data2['data'])
    ##!#    hash_data1 = np.ma.masked_where(MODIS_data2['data'] < \
    ##!#        np.nanmean(MODIS_data2['data']),MODIS_data1['data'])

    ##!#    nohash_data2 = np.ma.masked_where(MODIS_data2['data'] >= \
    ##!#        np.nanmean(MODIS_data2['data']),MODIS_data2['data'])
    ##!#    nohash_data1 = np.ma.masked_where(MODIS_data2['data'] >= \
    ##!#        np.nanmean(MODIS_data2['data']),MODIS_data1['data'])
    ##!#

    ##!#    hash2 = ax1.pcolor(MODIS_data2['lon'],MODIS_data2['lat'],\
    ##!#        hash_data2, hatch = '///', alpha=0., transform = datacrs),
    ##!#    hash1 = ax0.pcolor(MODIS_data1['lon'],MODIS_data1['lat'],\
    ##!#        hash_data1, hatch = '///', alpha=0., transform = datacrs),
    
    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(dt_date_str.strftime('%Y%m%d%H%M')) 
        hash2 = ax1.pcolor(MODIS_data2['lon'],MODIS_data2['lat'],\
            hash_data1, hatch = '///', alpha=0., transform = datacrs),
        hash1 = ax0.pcolor(MODIS_data1['lon'],MODIS_data1['lat'],\
            hash_data1, hatch = '///', alpha=0., transform = datacrs),
    
    # Step 4: Plot scatter MODIS channel comparison in third panel
    # ------------------------------------------------------------

    # Figure out which pixels are in the plume
    # ----------------------------------------
    max_ch = 350.

    tmp_data1 = np.copy(MODIS_data1['data'])
    tmp_data2 = np.copy(MODIS_data2['data'])

    tmp_data1 = np.ma.masked_where( (MODIS_data1['data'] > max_ch) | \
        (MODIS_data2['data'] > max_ch), tmp_data1)
    tmp_data2 = np.ma.masked_where( (MODIS_data1['data'] > max_ch) | \
        (MODIS_data2['data'] > max_ch), tmp_data2)

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
    ax2.scatter(MODIS_data1['data'][np.where(~hash_data1.mask)], \
        MODIS_data2['data'][np.where(~hash_data1.mask)], s=6, \
        color='tab:blue', label='Inside Plume')
    ax2.scatter(MODIS_data1['data'][np.where(hash_data1.mask)], \
        MODIS_data2['data'][np.where(hash_data1.mask)], s=6, \
        color='tab:orange', label='Outside Plume')
    #ax2.scatter(nohash_data1.compressed(), nohash_data2.compressed(), s=6, \
    #    color='tab:orange', label='Outside Plume')
    ax2.set_xlabel('Ch. ' + str(MODIS_data1['channel']) +' ' + MODIS_data1['variable'])
    ax2.set_ylabel('Ch. ' + str(MODIS_data2['channel']) +' ' + MODIS_data2['variable'])
    ax2.set_title('Pearson correlation: '+str(np.round(rval_p, 3)))
    ax2.legend()

    if(plot_ASOS_loc):
        print("Plotting ASOS site")
        plot_ASOS_locs(ax0,MODIS_data1)
        plot_ASOS_locs(ax1,MODIS_data1)

    if(save):
        cross_date = MODIS_data1['cross_date']
        file_time  = MODIS_data1['file_time']
        pdate = cross_date[:4] + cross_date[5:7] + cross_date[8:10] + file_time
        outname = 'modis_compare_ch' + str(channel1) + '_ch' + str(channel2) + '_' + pdate + '_2pannelscatter.png'
        plt.savefig(outname,dpi=300)
        print('Saved image',outname)
    else:
        plt.show()

    #return hash_data1, nohash_data1, MODIS_data1, MODIS_data2

def compare_MODIS_3scatter(date_str,channel0,channel1,channel2,channel3,\
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

    # Step 1: Call read_MODIS_channel to read the desired MODIS data from the
    # file and put it in a dictionary
    # -----------------------------------------------------------------------
    MODIS_data0 = read_MODIS_channel(dt_date_str.strftime('%Y%m%d%H%M'), channel0, zoom = True)
    MODIS_data1 = read_MODIS_channel(dt_date_str.strftime('%Y%m%d%H%M'), channel1, zoom = True)
    MODIS_data2 = read_MODIS_channel(dt_date_str.strftime('%Y%m%d%H%M'), channel2, zoom = True)
    MODIS_data3 = read_MODIS_channel(dt_date_str.strftime('%Y%m%d%H%M'), channel3, zoom = True)

    print("Data0 max = ",np.max(MODIS_data0['data']), "  Data0 min = ",np.min(MODIS_data0['data']))
    print("Data1 max = ",np.max(MODIS_data1['data']), "  Data1 min = ",np.min(MODIS_data1['data']))
    print("Data2 max = ",np.max(MODIS_data2['data']), "  Data2 min = ",np.min(MODIS_data2['data']))
    print("Data3 max = ",np.max(MODIS_data3['data']), "  Data3 min = ",np.min(MODIS_data3['data']))

    max_ch = 350.

    # Determine where the smoke is located
    # ------------------------------------
    hash_data1, nohash_data1 = find_plume(dt_date_str.strftime('%Y%m%d%H%M')) 

    tmp_data0 = np.copy(MODIS_data0['data'])
    tmp_data1 = np.copy(MODIS_data1['data'])
    tmp_data2 = np.copy(MODIS_data2['data'])
    tmp_data3 = np.copy(MODIS_data3['data'])
    tmp_lat0  = np.copy(MODIS_data0['lat'])
    tmp_lon0  = np.copy(MODIS_data0['lon'])

    if(not (tmp_data0.shape == tmp_data1.shape == tmp_data2.shape == tmp_data3.shape)):
        print("shape mismatch")
        shapes = []
        shapes.append(tmp_data0.shape)
        shapes.append(tmp_data1.shape)
        shapes.append(tmp_data2.shape)
        shapes.append(tmp_data3.shape)

        min_shape = min(shapes)
        print(min_shape)

        tmp_data0  = tmp_data0[:min_shape[0],:min_shape[1]]
        tmp_data1  = tmp_data1[:min_shape[0],:min_shape[1]]
        tmp_data2  = tmp_data2[:min_shape[0],:min_shape[1]]
        tmp_data3  = tmp_data3[:min_shape[0],:min_shape[1]]
        tmp_lat0   = tmp_lat0[:min_shape[0],:min_shape[1]]
        tmp_lon0   = tmp_lon0[:min_shape[0],:min_shape[1]]
        hash_data1 = hash_data1[:min_shape[0],:min_shape[1]]

    tmp_data0 = np.ma.masked_where( (abs(tmp_data0) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data2) > max_ch) | \
        (abs(tmp_data3) > max_ch), tmp_data0)
    tmp_data1 = np.ma.masked_where( (abs(tmp_data0) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data2) > max_ch) | \
        (abs(tmp_data3) > max_ch), tmp_data1)
    tmp_data2 = np.ma.masked_where( (abs(tmp_data0) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data2) > max_ch) | \
        (abs(tmp_data3) > max_ch), tmp_data2)
    tmp_data3 = np.ma.masked_where( (abs(tmp_data0) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data2) > max_ch) | \
        (abs(tmp_data3) > max_ch), tmp_data3)
    tmp_lat0  = np.ma.masked_where( (abs(tmp_data0) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data2) > max_ch) | \
        (abs(tmp_data3) > max_ch), tmp_lat0)
    tmp_lon0  = np.ma.masked_where( (abs(tmp_data0) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data2) > max_ch) | \
        (abs(tmp_data3) > max_ch), tmp_lon0)

    ##!#plot_data0 = tmp_data0.compressed()
    ##!#plot_data1 = tmp_data1.compressed()
    ##!#plot_data2 = tmp_data2.compressed()
    ##!#plot_data3 = tmp_data3.compressed()
    ##!#plot_lat0  = tmp_lat0.compressed()
    ##!#plot_lon0  = tmp_lon0.compressed()

    cross_date = MODIS_data0['cross_date']
    file_time  = MODIS_data0['file_time']

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

    plot_scatter(ax0, tmp_data0, tmp_data1, MODIS_data0, MODIS_data1, \
        channel0, channel1, hash_data1)
    plot_scatter(ax1, tmp_data0, tmp_data2, MODIS_data0, MODIS_data2, \
        channel0, channel2, hash_data1)
    plot_scatter(ax2, tmp_data0, tmp_data3, MODIS_data0, MODIS_data3, \
        channel0, channel3, hash_data1)


    # Remove the nans from the MODIS data, lats, and lons for both the
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
        data = h5py.File(plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['omi'],'r')
        LAT   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude'][:,:]
        LON   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude'][:,:]
        UVAI  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex'][:,:]
        XTRACK = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags'][:,:]
        mask_LAT = np.ma.masked_where( (XTRACK < -2e5) | (UVAI < 2.), LAT)
        mask_LON = np.ma.masked_where( (XTRACK < -2e5) | (UVAI < 2.), LON)
        mask_UVAI = np.ma.masked_where((XTRACK < -2e5) | (UVAI < 2.), UVAI)
        mask_LAT  = np.ma.masked_where((((LAT < plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) | \
                             (LAT > plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1])) | \
                            ((LON < plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) | \
                             (LON > plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]))), mask_LAT)
        mask_LON  = np.ma.masked_where((((LAT < plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) | \
                             (LAT > plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1])) | \
                            ((LON < plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) | \
                             (LON > plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]))), mask_LON)
        mask_UVAI = np.ma.masked_where((((LAT < plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) | \
                             (LAT > plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1])) | \
                            ((LON < plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) | \
                             (LON > plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]))), mask_UVAI)


        #hash_data1, nohash_data1 = find_plume(filename) 
        # Colocate the masked OMI data with the MODIS channel0 data
        # ---------------------------------------------------------
        print("Colocating OMI data")

        if(avg_pixel):
            print('averaging MODIS pixels')

            ##!## Determine the extreme bounds of the current 
            ##!#max_MODIS_lat = np.nanmax(tmp_lat0)
            ##!#min_MODIS_lat = np.nanmin(tmp_lat0)
            ##!#max_MODIS_lon = np.nanmax(tmp_lon0)
            ##!#min_MODIS_lon = np.nanmin(tmp_lon0)

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

            # Declare full arrays to hold the averaged MODIS data
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

                # For each OMI pixel, find the in-plume and out-of-plume MODIS
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

                # Use the indices to pull out the MODIS pixels
                hash_modis_data   = hash_plot_data1[hash_in]
                hash_modis_lat    = hash_plot_lat1[hash_in]
                hash_modis_lon    = hash_plot_lon1[hash_in]
                nohash_modis_data = nohash_plot_data1[nohash_in]
                nohash_modis_lat  = nohash_plot_lat1[nohash_in]
                nohash_modis_lon  = nohash_plot_lon1[nohash_in]

                #print(work_crnrLAT[ii,:], work_crnrLON[ii,:], hash_modis_lat, \
                #    nohash_modis_lat)
            
                # Use the corner values to make a Polygon and determine if
                # each close MODIS pixel is within the OMI pixel
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

                ##!## Loop over the close MODIS pixels and determine which ones
                ##!## are within the polygon
                ##!## ---------------------------------------------------------
                ##!#for jj in range(len(hash_modis_data)):

                ##!#    # Average each MODIS value that is within the polygon
                ##!#    # into a value, and insert the value into the avg
                ##!#    # modis array for hash and nohash.
                ##!#    # --------------------------------------------------- 
                ##!#    mpoint = Point(hash_modis_lon[jj], hash_modis_lat[jj])

                ##!#    inside_vals = [omi_poly.contains(mpoint)  ] 

            # end screened corner data loop

            # Plot the screened OMI data vs the averaged MODIS data
            # -----------------------------------------------------

            #print('hashed averaged MODIS data = ',hash_avg_modis)
            #print('nohashed averaged MODIS data = ',nohash_avg_modis)

            print(work_UVAI[~np.ma.masked_invalid(hash_avg_modis).mask])
            
            hrval_p = spearmanr(np.ma.masked_invalid(hash_avg_modis).compressed(), \
                work_UVAI[~np.ma.masked_invalid(hash_avg_modis).mask])[0]
            axo.scatter(np.ma.masked_invalid(hash_avg_modis).compressed(), \
                work_UVAI[~np.ma.masked_invalid(hash_avg_modis).mask], s=6, \
                color='tab:blue')
            #axo.scatter(np.ma.masked_invalid(nohash_avg_modis).compressed(), \
            #    work_UVAI[~np.ma.masked_invalid(nohash_avg_modis).mask], s=6, \
            #    color='tab:orange')
            #axo.scatter(nohash_plot_data0, nohash_match_OMI, s=6, \
            #    color='tab:orange', label='Outside Plume')

            axo.set_xlabel('Averaged Ch. ' + str(MODIS_data1['channel']) +' [' + \
                channel_dict[str(channel1)]['Bandwidth_label'] + \
                MODIS_data1['variable'])
            axo.set_title('Smoke correlation: '+str(np.round(hrval_p, 3)))
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

            hrval_p = spearmanr(hash_plot_data1, \
                hash_match_OMI)[0]
            axo.scatter(hash_plot_data1, hash_match_OMI, s=6, \
                color='tab:blue')
            #axo.scatter(nohash_plot_data0, nohash_match_OMI, s=6, \
            #    color='tab:orange')

            axo.set_xlabel('Ch. ' + str(MODIS_data1['channel']) +' [' + \
                channel_dict[str(channel1)]['Bandwidth_label'] + \
                MODIS_data1['variable'])
            axo.set_title('Smoke correlation: '+str(np.round(hrval_p, 3)))
            axo.set_ylabel('OMI UVAI')
        
        data.close()

    # End compare_OMI

    if(compare_CERES):
        print("Reading CERES data")

        mask_LAT, mask_LOn, mask_swf, mask_lwf = read_CERES_match_MODIS(date_str)

        ##!## NOTE: need to use the time variable to screen out any data
        ##!## that are not close to the event time.
        ##!#base_date = datetime(year=1970,month=1,day=1)
        ##!#start_date = dt_date_str - timedelta(hours = 1)
        ##!#end_date   = dt_date_str + timedelta(hours = 2)

        ##!#print(start_date, end_date)

        ##!#data = Dataset(plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['ceres'],'r')
        ##!#LAT   = 90. - data.variables['Colatitude_of_CERES_FOV_at_surface'][:]
        ##!#LON   = data.variables['Longitude_of_CERES_FOV_at_surface'][:]
        ##!#LON[LON>179.99] = -360.+LON[LON>179.99]
        ##!#swflux  = data.variables['CERES_SW_TOA_flux___upwards'][:]
        ##!#lwflux  = data.variables['CERES_LW_TOA_flux___upwards'][:]
        ##!#time  = data.variables['time'][:]
        ##!#local_time = np.array([base_date + relativedelta(days = ttime) for ttime in time])
        ##!#data.close()

        ##!#mask_LAT = LAT[ \
        ##!#    (LAT >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) & \
        ##!#    (LAT <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1]) & \
        ##!#    (LON >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) & \
        ##!#    (LON <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]) & \
        ##!#    (swflux > 0) & (lwflux > 0)]
        ##!#mask_LON = LON[ \
        ##!#    (LAT >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) & \
        ##!#    (LAT <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1]) & \
        ##!#    (LON >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) & \
        ##!#    (LON <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]) & \
        ##!#    (swflux > 0) & (lwflux > 0)]
        ##!#mask_swf = swflux[ \
        ##!#    (LAT >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) & \
        ##!#    (LAT <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1]) & \
        ##!#    (LON >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) & \
        ##!#    (LON <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]) & \
        ##!#    (swflux > 0) & (lwflux > 0)]
        ##!#mask_lwf = lwflux[ \
        ##!#    (LAT >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) & \
        ##!#    (LAT <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1]) & \
        ##!#    (LON >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) & \
        ##!#    (LON <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]) & \
        ##!#    (swflux > 0) & (lwflux > 0)]
        ##!###!#mask_time = local_time[ \
        ##!###!#    (LAT >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) & \
        ##!###!#    (LAT <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1]) & \
        ##!###!#    (LON >= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) & \
        ##!###!#    (LON <= plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]) & \
        ##!###!#    (swflux > 0) & (lwflux > 0)]


        if(avg_pixel):
            print('Averaging MODIS pixels')
            
            # Average a certain number of MODIS pixels around each CERES
            # pixel.
            # ----------------------------------------------------------
            hash_match_SWF    = np.full(mask_swf.shape, np.nan)
            hash_match_LWF    = np.full(mask_lwf.shape, np.nan)
            hash_match_LAT    = np.full(mask_LAT.shape, np.nan)
            hash_match_LON    = np.full(mask_LON.shape, np.nan)

            # Loop over the lower-resolution AIRS data.
            for ii in range(mask_swf.shape[0]):

                # Find the gridpoint in the MODIS lat/lon data that 
                # corresponds to the AIRS pixel
                # ---------------------------------------------------- 
                ho_idx = nearest_gridpoint(mask_LAT[ii], mask_LON[ii],\
                    hash_plot_lat0, hash_plot_lon0)

                if(len(ho_idx[0]) > 1):
                    ho_idx = (np.array([ho_idx[0][0]])), \
                        (np.array([ho_idx[1][0]]))
        
                # Average the n MODIS pixels around the current CERES
                # pixel.
                # --------------------------------------------------- 
                hash_match_SWF[ii] = np.nanmean(hash_plot_data0[ho_idx[0][0]-1:ho_idx[0][0]+2])
                hash_match_LWF[ii] = np.nanmean(hash_plot_data0[ho_idx[0][0]-1:ho_idx[0][0]+2])
                #hash_match_SWF[ii] = np.nanmean(hash_plot_data0[ho_idx[0]-1:ho_idx[0]+2,\
                #    ho_idx[1]-1:ho_idx[1]+2])
                #hash_match_LWF[ii] = np.nanmean(hash_plot_data0[ho_idx[0]-1:ho_idx[0]+2,\
                #    ho_idx[1]-1:ho_idx[1]+2])
                hash_match_LAT[ii] = hash_plot_lat0[ho_idx]
                hash_match_LON[ii] = hash_plot_lon0[ho_idx]
            #xy = np.vstack([plot_data0,match_OMI])
            #z = stats.gaussian_kde(xy)(xy)
            #axo.scatter(plot_data0,match_OMI,c=z,s=6)

            lrval_p = spearmanr(hash_match_LWF, \
                mask_lwf)[0]
            axcl.scatter(hash_match_LWF, mask_lwf,\
                s = 6, color='tab:blue')

            srval_p = spearmanr(hash_match_SWF, \
                mask_swf)[0]
            axcs.scatter(hash_match_SWF, mask_swf,\
                s = 6, color='tab:blue')
            ##!#axcl.scatter(nohash_plot_data0, nohash_match_LWF,\
            ##!#    s = 6, color='tab:orange')
            ##!#axcs.scatter(nohash_plot_data0, nohash_match_SWF,\
            ##!#    s = 6, color='tab:orange')

            axcl.set_xlabel('Ch. ' + str(MODIS_data0['channel']) +' [' + \
                channel_dict[str(channel0)]['Bandwidth_label'] + \
                MODIS_data0['variable'])
            axcl.set_title('Smoke correlation: '+str(np.round(lrval_p, 3)))
            axcs.set_xlabel('Ch. ' + str(MODIS_data0['channel']) +' [' + \
                channel_dict[str(channel0)]['Bandwidth_label'] + \
                MODIS_data0['variable'])
            axcs.set_title('Smoke correlation: '+str(np.round(srval_p, 3)))
            axcl.set_ylabel('CERES LWF [W/m2]')
            axcs.set_ylabel('CERES SWF [W/m2]')

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

            lrval_p = spearmanr(hash_plot_data0, \
                hash_match_LWF)[0]
            axcl.scatter(hash_plot_data0, hash_match_LWF,\
                s = 6, color='tab:blue')

            srval_p = spearmanr(hash_plot_data0, \
                hash_match_SWF)[0]
            axcs.scatter(hash_plot_data0, hash_match_SWF,\
                s = 6, color='tab:blue')
            ##!#axcl.scatter(nohash_plot_data0, nohash_match_LWF,\
            ##!#    s = 6, color='tab:orange')
            ##!#axcs.scatter(nohash_plot_data0, nohash_match_SWF,\
            ##!#    s = 6, color='tab:orange')

            axcl.set_xlabel('Ch. ' + str(MODIS_data0['channel']) +' [' + \
                channel_dict[str(channel0)]['Bandwidth_label'] + \
                MODIS_data0['variable'])
            axcl.set_title('Smoke correlation: '+str(np.round(lrval_p, 3)))
            axcs.set_xlabel('Ch. ' + str(MODIS_data0['channel']) +' [' + \
                channel_dict[str(channel0)]['Bandwidth_label'] + \
                MODIS_data0['variable'])
            axcs.set_title('Smoke correlation: '+str(np.round(srval_p, 3)))
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

def plot_MODIS_spatial(MODIS_data, pax, zoom):

    if(MODIS_data['channel'] == 'wv_ir'):
        vmax = 1.5
    else:
        vmax = np.nanmax(MODIS_data['data'])

    #mesh = ax.pcolormesh(MODIS_data['lon'],MODIS_data['lat'],\
    #    MODIS_data['data'],cmap = MODIS_data['colors'], shading='auto', \
    #    transform = datacrs) 

    # Plot channel 1
    mesh1 = pax.pcolormesh(MODIS_data['lon'],MODIS_data['lat'],\
        MODIS_data['data'],cmap = MODIS_data['colors'], shading='auto', \
        vmin = np.nanmin(MODIS_data['data']), \
        vmax = vmax, transform = datacrs) 

    cbar1 = plt.colorbar(mesh1,ax=pax,orientation='vertical',\
        pad=0.03,label=MODIS_data['variable'])

    pax.add_feature(cfeature.BORDERS)
    pax.add_feature(cfeature.STATES)
    pax.coastlines()
    if(zoom):
        pax.set_extent([plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][0], \
                        plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][1], \
                        plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][0], \
                        plot_limits_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][1]],\
                        datacrs)
    pax.set_title('Channel ' + str(MODIS_data['channel']) + '\n' + \
        channel_dict[str(MODIS_data['channel'])]['Bandwidth_label']) 


# dtype must be 'SWF' or 'LWF'
def plot_CERES_spatial(date_str, mask_LAT, mask_LON, mask_data, dtype, pax, \
        zoom = False):

    #print(plot_limits_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'])

    plabel = 'TOA ' + dtype + '[W/m2]'
    ptitle = 'CERES TOA ' + dtype

    dt_date_str = datetime.strptime(date_str, "%Y%m%d%H%M")

    #scat3 = axcs.scatter(mask_LON, mask_LAT,mask_swf, transform = datacrs)
    mesh3 = pax.scatter(mask_LON.compressed(), mask_LAT.compressed(),\
        s = 120,marker='s',c = mask_data.compressed(),cmap='plasma', \
        transform = datacrs)
    pax.add_feature(cfeature.BORDERS)
    pax.add_feature(cfeature.STATES)
    pax.coastlines()
    if(zoom):
        pax.set_extent([plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][0], \
                        plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][1], \
                        plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][0], \
                        plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][1]],\
                        datacrs)
    cbar3 = plt.colorbar(mesh3,ax=pax,orientation='vertical',\
        pad=0.03,label=plabel)
    pax.set_title(ptitle)

def plot_scatter(lax,data1, data2, MODIS_data1, MODIS_data2, channel1, \
        channel2, hash_data1):

    hrval_p = spearmanr(data1[np.where(~hash_data1.mask)].compressed(), \
        data2[np.where(~hash_data1.mask)].compressed())[0]
    nrval_p = spearmanr(data1[np.where(hash_data1.mask)].compressed(), \
        data2[np.where(hash_data1.mask)].compressed())[0]
    #print("0-1 Pearson:  ",hrval_p)
    lax.scatter(data1[np.where(hash_data1.mask)].compressed(), \
        data2[np.where(hash_data1.mask)].compressed(), s=6, \
        color='tab:orange')
    lax.scatter(data1[np.where(~hash_data1.mask)].compressed(), \
        data2[np.where(~hash_data1.mask)].compressed(), s=6, \
        color='tab:blue')

    lax.set_xlabel('Ch. ' + str(MODIS_data1['channel']) +' [' + \
        channel_dict[str(channel1)]['Bandwidth_label'] + '] '+ \
        MODIS_data1['variable'])
    lax.set_ylabel('Ch. ' + str(MODIS_data2['channel']) +' [' + \
        channel_dict[str(channel2)]['Bandwidth_label'] + ']' + \
        MODIS_data2['variable'])
    plt.suptitle('Aqua MODIS ' + MODIS_data1['cross_date'] + ' ' + \
        MODIS_data1['file_time'])
    lax.set_title('Smoke correlation: '+str(np.round(hrval_p, 3))+'\n'+\
                  'Clear correlation: '+str(np.round(nrval_p, 3)))
    

def plot_OMI_spatial(date_str, LAT, LON, mask_UVAI, pax, zoom = False):

    dt_date_str = datetime.strptime(date_str, "%Y%m%d%H%M")

    mesh3 = pax.pcolormesh(LON,LAT,mask_UVAI, cmap = 'plasma', shading='auto', \
        vmin = np.nanmin(mask_UVAI), vmax = np.nanmax(mask_UVAI), transform = datacrs) 
    pax.add_feature(cfeature.BORDERS)
    pax.add_feature(cfeature.STATES)
    pax.coastlines()
    if(zoom):
        pax.set_extent([plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][0], \
                        plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][1], \
                        plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][0], \
                        plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][1]],\
                        datacrs)
    cbar3 = plt.colorbar(mesh3,ax=pax,orientation='vertical',\
        pad=0.03,label='OMI UVAI')
    pax.set_title('OMI UVAI',fontsize=10)


def plot_combined_imagery(date_str,channel1 = 1, channel2 = 5, channel3 = 31,\
        zoom=True,save=False,composite=True):
    
    # ------------------------------------------------------------------------
    #
    # For 20210722:
    #   1. True color 20210721       2. True color 20210722
    #   3. Ch 1 20210722             4. Ch 5 20210722
    #   5. Ch 31 20210722            6. WV 20210722
    #   7. CERES SW 20210722         8. CERES LW 20210722
    #
    # ------------------------------------------------------------------------
    
    plt.close('all')
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # Read true color data for this date
    var, crs, lat_lims, lon_lims = read_true_color(date_str,composite=composite)

    if(date_str == '202107222110'):
        # Read true color data for the previous date
        prev_date_str = datetime.strptime('202107212030',"%Y%m%d%H%M")
        var1, crs1, lat_lims1, lon_lims1 = read_true_color('202107212030',\
            composite=composite)

        fig = plt.figure(figsize=(9,12))
        ax5 = fig.add_subplot(4,2,1,projection = crs1)   # Previous day true      
        ax0 = fig.add_subplot(4,2,2,projection = crs)    # true color
        ax1 = fig.add_subplot(4,2,3,projection = mapcrs) # Ch 1
        ax2 = fig.add_subplot(4,2,5,projection = mapcrs) # Ch 5
        ax3 = fig.add_subplot(4,2,7,projection = mapcrs) # Ch 31
        ax4 = fig.add_subplot(4,2,4,projection = mapcrs) # WV
        ax6 = fig.add_subplot(4,2,6,projection = mapcrs) # SW
        ax7 = fig.add_subplot(4,2,8,projection = mapcrs) # LW

        # Plot the true-color data for the previous date
        # ----------------------------------------------
        ax5.imshow(var1.data, transform = crs1, extent=(var1.x[0], var1.x[-1], \
            var1.y[-1], var1.y[0]), origin='upper')

        # Zoom in the figure if desired
        # -----------------------------
        if(zoom):
            ax5.set_extent([lon_lims1[0],lon_lims1[1],lat_lims1[0],\
                lat_lims1[1]],crs = datacrs)
        ax5.set_title('Aqua MODIS\n' + prev_date_str.strftime('%Y-%m-%d %H:%M'))

    elif(date_str == '202108062025'):

        fig = plt.figure(figsize=(9,12))
        ax0 = fig.add_subplot(4,2,1,projection = crs)    # True color
        ax1 = fig.add_subplot(4,2,3,projection = mapcrs) # Ch 1
        ax2 = fig.add_subplot(4,2,5,projection = mapcrs) # Ch 5
        ax3 = fig.add_subplot(4,2,7,projection = mapcrs) # Ch 31
        ax4 = fig.add_subplot(4,2,4,projection = mapcrs) # WV 
        ax6 = fig.add_subplot(4,2,6,projection = mapcrs) # SW
        ax7 = fig.add_subplot(4,2,8,projection = mapcrs) # LW
        ax5 = fig.add_subplot(4,2,2,projection = mapcrs) # OMI AI
   
        # ---------------------------------------------------------------------
        #
        # Read the OMI data
        #
        # ---------------------------------------------------------------------

        LAT, LON, mask_UVAI = read_OMI_match_MODIS(date_str)
 
        plot_OMI_spatial(date_str, LAT, LON, mask_UVAI, ax5, zoom = zoom)

    # ----------------------------------------------------------------------
    #
    # Plot the true-color data
    #
    # ----------------------------------------------------------------------
    ax0.imshow(var.data, transform = crs, extent=(var.x[0], var.x[-1], \
        var.y[-1], var.y[0]), origin='upper')

    # Zoom in the figure if desired
    # -----------------------------
    if(zoom):
        ax0.set_extent([lon_lims[0],lon_lims[1],lat_lims[0],lat_lims[1]],\
                       crs = datacrs)
        zoom_add = '_zoom'
    else:
        zoom_add = ''

    ax0.set_title('Aqua MODIS\n'+dt_date_str.strftime('%Y-%m-%d %H:%M'))

    # ----------------------------------------------------------------------
    #
    # Read the MODIS data
    #
    # ----------------------------------------------------------------------

    # Call read_MODIS_channel to read the desired MODIS data from the file 
    # and put it in a dictionary
    # ---------------------------------------------------------------------
    MODIS_data1 = read_MODIS_channel(dt_date_str.strftime("%Y%m%d%H%M"), \
        channel1, zoom = True)
    MODIS_data2 = read_MODIS_channel(dt_date_str.strftime("%Y%m%d%H%M"), \
        channel2, zoom = True)
    MODIS_data3 = read_MODIS_channel(dt_date_str.strftime("%Y%m%d%H%M"), \
        channel3, zoom = True)
    MODIS_dataw = read_MODIS_channel(dt_date_str.strftime("%Y%m%d%H%M"), \
        'wv_ir', zoom = True)

    # ----------------------------------------------------------------------
    #
    # Plot the MODIS data
    #
    # ----------------------------------------------------------------------
    
    # Plot channels 1, 5, 31, and WV data
    plot_MODIS_spatial(MODIS_data1, ax1, zoom = zoom)
    plot_MODIS_spatial(MODIS_data2, ax2, zoom = zoom)
    plot_MODIS_spatial(MODIS_data3, ax3, zoom = zoom)
    plot_MODIS_spatial(MODIS_dataw, ax4, zoom = zoom)

    # ----------------------------------------------------------------------
    #
    # Read the CERES data
    #
    # ----------------------------------------------------------------------

    mask_LAT, mask_LON, mask_swf, mask_lwf = \
        read_CERES_match_MODIS(date_str)

    # ----------------------------------------------------------------------
    #
    # Plot the CERES data
    #
    # ----------------------------------------------------------------------

    plot_CERES_spatial(date_str,mask_LAT, mask_LON, mask_swf, 'SWF', ax6, zoom = zoom)
    plot_CERES_spatial(date_str,mask_LAT, mask_LON, mask_lwf, 'LWF', ax7, zoom = zoom)

    plt.show()
    ##!#if(save):
    ##!#    outname = 'modis_true_color_' + date_str + zoom_add + cmpst_add + '.png'
    ##!#    plt.savefig(outname,dpi=300)
    ##!#    print("Saved image",outname)
    ##!#else:
    ##!#    plt.show()

# Compare colocated MODIS and ob data for two dates
def colocate_comparison(date1, date2, channel = 31):
    # Read in MODIS data for both cases
    # ---------------------------------
    dt_date_str1 = datetime.strptime(date1,"%Y%m%d%H%M")
    dt_date_str2 = datetime.strptime(date2,"%Y%m%d%H%M")
    filename1 = plot_limits_dict[dt_date_str1.strftime('%Y-%m-%d')][dt_date_str1.strftime('%H%M')]['modis']
    filename2 = plot_limits_dict[dt_date_str2.strftime('%Y-%m-%d')][dt_date_str2.strftime('%H%M')]['modis']

    MODIS_data1 = read_MODIS_channel(dt_date_str1.strftime('%Y%m%d%H%M'), channel, zoom = True)
    MODIS_data2 = read_MODIS_channel(dt_date_str2.strftime('%Y%m%d%H%M'), channel, zoom = True)

    #print(MODIS_data1,MODIS_data2)

    # Use the colocation code to extract matching data
    # ------------------------------------------------
    compare_data1 = nearest_grid_values(MODIS_data1)
    compare_data2 = nearest_grid_values(MODIS_data2)

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
    ax2.set_ylabel('Colocated MODIS Channel ' + str(channel) + ' Brightness Temp [K] ',color='tab:orange') 
    plt.legend() 
    plt.show()
 
