"""
  NAME:

  PURPOSE:
  
    NOTE: The GOES channel and true color functions are designed to work with
    HDF GOES files retriefed from 
    the data ordering website at this address:
    https://ladsweb.modaps.eosdis.nasa.gov/search/order/1/GOES:Aqua


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
from satpy import find_files_and_readers
from satpy.scene import Scene
from satpy.writers import get_enhanced_image
from glob import glob

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# Set up global variables
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
datacrs = ccrs.PlateCarree()
data_dir = '/home/bsorenson/data/GOES/goes17_abi/'
debug = False

zoom_dict = {
    'Finland': [10,55,65,80]
}

proj_dict = {
    'Finland': ccrs.NorthPolarStereo(central_longitude = 35.) 
}

# This contains the smokey station for each case
case_dict = {
    '202107132110': 'O05',
    '202107222110': 'O05',
    '202108062025': 'TPH',
}

# - 0.64  (band 2, Red)
# - 2.25  (band 6, Cloud particle size)
# - 3.90  (band 7, Shortwave Window)
# - 6.18  (band 8, Upper-Level Water Vapor)
# - 6.95  (band 9, Mid-Level Water Vapor)
# - 7.34  (band 10, Lower-Level Water Vapor)
# - 10.35 (band 13, Clean IR Longwave Window)
channel_dict = {
    '1': {
        'name': 'Blue',
        'wavelength': 0.47
    },\
    '2': {
        'name': 'Red',
        'wavelength': 0.64
    },\
    '3': {
        'name': 'Veggie',
        'wavelength': 0.86
    },\
    '4': {
        'name': 'Cirrus',
        'wavelength': 1.38
    },\
    '5': {
        'name': 'Snow/Ice',
        'wavelength': 1.61
    },\
    '6': {
        'name': 'Cloud Particle Size',
        'wavelength': 2.25
    },\
    '7': {
        'name': 'Shortwave Window',
        'wavelength': 3.90
    },\
    '8': {
        'name': 'Upper-Level Water Vapor',
        'wavelength': 6.18
    },\
    '9': {
        'name': 'Mid-Level Water Vapor',
        'wavelength': 6.95
    },\
    '10': {
        'name': 'Lower-Level Water Vapor',
        'wavelength': 7.34
    },\
    '11': {
        'name': 'Cloud-Top Phase',
        'wavelength': 8.50
    },\
    '12': {
        'name': 'Ozone',
        'wavelength': 9.61
    },\
    '13': {
        'name': 'Clean IR Longwave Window',
        'wavelength': 10.35
    },\
    '14': {
        'name': 'IR Longwave Window',
        'wavelength': 11.20
    },\
    '15': {
        'name': 'Dirty IR Longwave Window',
        'wavelength': 12.30
    },\
    '16': {
        'name': 'CO2 Longwave Infrared',
        'wavelength': 13.30
    }
}

for key in channel_dict.keys():
    if(channel_dict[key]['wavelength'] is not None):
        channel_dict[key]['wavelength_label'] = \
            str(channel_dict[key]['wavelength']) + ' μm'
    else:
        channel_dict[key]['wavelength_label'] = ''

plot_limits_dict = {
    "2021-07-13": {
        '2110': {
            'asos': 'asos_data_20210713.csv',
            #'goes': '/home/bsorenson/data/GOES/Aqua/MYD021KM.A2021203.2110.061.2021204155922.hdf',
            #'mdswv': '/home/bsorenson/data/GOES/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
            #'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072210-2021072221.nc',
            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
            'Lat': [39.5, 42.0],
            'Lon': [-122.0, -119.5],
            'goes_Lat': [39.0, 42.5],
            'goes_Lon': [-123., -119.]
        }
    },
    "2021-07-20": {
        'asos': 'asos_data_20210720.csv',
        #'goes': '/home/bsorenson/data/GOES/Aqua/MYD021KM.A2021201.2125.061.2021202154814.hdf',
        #'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072010-2021072021.nc',
        #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.20.214.L2.SUBS2RET.v6.0.32.0.G21202153435.hdf'],
        'Lat': [39.5, 42.0],
        'Lon': [-122.0, -119.5],
        'data_lim': {
            1:  [0.05, 0.5],
            31: [270., 330.],
        },
        'goes_Lat': [39.5, 42.0],
        'goes_Lon': [-122.0, -119.5]
    },
    "2021-07-21": {
        'asos': 'asos_data_20210722_4.csv',
        #'goes': '/home/bsorenson/data/GOES/Aqua/MYD021KM.A2021202.2030.061.2021203174050.hdf',
        #'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072109-2021072122.nc',
        #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.21.205.L2.SUBS2RET.v6.0.32.0.G21203185004.hdf'],
        'Lat': [39.5, 42.0],
        'Lon': [-122.0, -119.5],
        'Lat': [39.5, 42.0],
        'Lon': [-122.0, -119.5],
        'data_lim': {
            1:  [0.05, 0.5],
            31: [270., 330.],
        },
        'goes_Lat': [39.5, 42.0],
        'goes_Lon': [-122.0, -119.5]
    },
    "2021-07-22": {
        'asos': 'asos_data_20210722_4.csv',
        #'goes': '/home/bsorenson/data/GOES/Aqua/MYD021KM.A2021203.2110.061.2021204155922.hdf',
        #'mdswv': '/home/bsorenson/data/GOES/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
        #'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072210-2021072221.nc',
        'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
        'Lat': [39.5, 42.0],
        'Lon': [-122.0, -119.5],
        #'Lat': [39.5, 42.0],
        #'Lon': [-122.0, -119.5],
        'data_lim': {
            1:  [0.05, 0.5],
            5:  [None, None],
            31: [270., 330.],
            32: [270., 330.],
            'wv_ir': [0.2, 1.5],
        },
        #'goes_Lat': [39.5, 42.0],
        #'goes_Lon': [-122.0, -119.5]
        'goes_Lat': [39.5, 42.0],
        'goes_Lon': [-122.0, -119.5]
    },
    "2021-07-23": {
        'asos': 'asos_data_20210722_2.csv',
        'goes': '/home/bsorenson/data/GOES/Aqua/MYD021KM.A2021204.2155.061.2021205153516.hdf',
        #'mdswv': '/home/bsorenson/data/GOES/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
        #'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072210-2021072221.nc',
        #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
        'Lat': [39.5, 42.0],
        'Lon': [-122.0, -119.5],
        #'Lat': [39.5, 42.0],
        #'Lon': [-122.0, -119.5],
        'data_lim': {
            1:  [0.05, 0.5],
            31: [270., 330.],
        },
        #'goes_Lat': [39.5, 42.0],
        #'goes_Lon': [-122.0, -119.5]
        'goes_Lat': [39.5, 42.0],
        'goes_Lon': [-122.0, -119.5]
    },
    "2021-08-04": {
        'asos': 'asos_data_20210806.csv',
        #'asos': 'asos_nevada_20210806.csv',
        #'goes': '/home/bsorenson/data/GOES/Aqua/MYD021KM.A2021218.2025.061.2021219151802.hdf',
        #'mdswv': '/home/bsorenson/data/GOES/Aqua/MYD05_L2.A2021218.2025.061.2021219152751.hdf',
        'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.04.206.L2.SUBS2RET.v6.0.32.0.G21217152448.hdf',\
                 '/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.04.207.L2.SUBS2RET.v6.0.32.0.G21217152904.hdf'],\
        #'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
        'Lat': [36.0, 39.0],
        'Lon': [-118.0, -114.0],
        'goes_Lat': [35.0, 40.0],
        'goes_Lon': [-119., -113.]
    },
    "2021-08-05": {
        '2120': {
            'asos': 'asos_data_20210805.csv',
            'goes': '/home/bsorenson/data/GOES/Aqua/MYD021KM.A2021217.2120.061.2021218164201.hdf',
            'mdswv': '/home/bsorenson/data/GOES/Aqua/MYD05_L2.A2021217.2120.061.2021218165546.hdf',
            'ceres': '/home/bsorenson/data/CERES/FLASHFlux/Aqua/FLASH_SSF_Aqua_Version4A_Subset_2021080508-2021080521.nc',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.05.214.L2.SUBS2RET.v6.0.32.0.G21218175548.hdf'],
            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0805t2038-o90733_v003-2021m0807t014855.he5',
            'Lat': [36.5, 39.0],
            'Lon': [-118.0, -114.0],
            'goes_Lat': [36.0, 39.0],
            'goes_Lon': [-118., -114.]
        },
        '2125': {
            'asos': 'asos_california_20210805.csv',
            'goes': '/home/bsorenson/data/GOES/Aqua/MYD021KM.A2021217.2125.061.2021218161010.hdf',
            'ceres': '/home/bsorenson/data/CERES/FLASHFlux/Aqua/FLASH_SSF_Aqua_Version4A_Subset_2021080508-2021080521.nc',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.05.214.L2.SUBS2RET.v6.0.32.0.G21218175548.hdf'],
            'Lat': [39.5, 42.0],
            'Lon': [-122.0, -119.0],
            'goes_Lat': [39.5, 42.0],
            'goes_Lon': [-122., -119.]
        }
    },
    "2021-08-06": {
        '2025': {
            'asos': 'asos_data_20210806.csv',
            #'asos': 'asos_nevada_20210806.csv',
            'goes': '/home/bsorenson/data/GOES/Aqua/MYD021KM.A2021218.2025.061.2021219151802.hdf',
            'mdswv': '/home/bsorenson/data/GOES/Aqua/MYD05_L2.A2021218.2025.061.2021219152751.hdf',
            'ceres': '/home/bsorenson/data/CERES/FLASHFlux/Aqua/FLASH_SSF_Aqua_Version4A_Subset_2021080609-2021080620.nc',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.06.204.L2.SUBS2RET.v6.0.32.0.G21219130523.hdf',\
                     '/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.06.205.L2.SUBS2RET.v6.0.32.0.G21219130455.hdf'],
            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
            'Lat': [36.0, 39.0],
            'Lon': [-118.0, -114.0],
            'goes_Lat': [36.0, 39.0],
            'goes_Lon': [-118., -114.]
        }
    },
    "2021-08-07": {
        '2110': {
            'asos': 'asos_data_20210806.csv',
            #'asos': 'asos_nevada_20210806.csv',
            'goes': '/home/bsorenson/data/GOES/Aqua/MYD021KM.A2021219.2110.061.2021220151612.hdf',
            #'mdswv': '/home/bsorenson/data/GOES/Aqua/MYD05_L2.A2021218.2025.061.2021219152751.hdf',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.07.212.L2.SUBS2RET.v6.0.32.0.G21220123225.hdf'],\
            #'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
            'Lat': [36.0, 39.0],
            'Lon': [-118.0, -114.0],
            'goes_Lat': [35.0, 40.0],
            'goes_Lon': [-119., -113.]
        }
    },
    "2021-08-08": {
        '2110': {
            'asos': 'asos_data_20210806.csv',
            #'asos': 'asos_nevada_20210806.csv',
            #'goes': '/home/bsorenson/data/GOES/Aqua/MYD021KM.A2021218.2025.061.2021219151802.hdf',
            #'mdswv': '/home/bsorenson/data/GOES/Aqua/MYD05_L2.A2021218.2025.061.2021219152751.hdf',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.08.202.L2.SUBS2RET.v6.0.32.0.G21221124420.hdf',\
                     '/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.08.203.L2.SUBS2RET.v6.0.32.0.G21221185932.hdf'],\
            #'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
            'Lat': [36.0, 39.0],
            'Lon': [-118.0, -114.0],
            'goes_Lat': [35.0, 40.0],
            'goes_Lon': [-119., -113.]
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
            'goes': '/home/bsorenson/data/GOES/Aqua/MYD021KM.A2021242.2115.061.2021243183953.hdf',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.30.213.L2.SUBS2RET.v6.0.32.0.G21243151114.hdf'],
            'Lat': [38.0, 40.0],
            'Lon': [-121.0, -118.5]
        }
    },
    "2021-09-01": {
        '2105': {
            'asos': 'asos_data_20210830.csv',
            'goes': '/home/bsorenson/data/GOES/Aqua/MYD021KM.A2021244.2105.061.2021245152256.hdf',
            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.30.213.L2.SUBS2RET.v6.0.32.0.G21243151114.hdf'],
            'Lat': [38.0, 42.0],
            'Lon': [-121.5, -118.0],
            'goes_Lat': [38.0, 42.0],
            'goes_Lon': [-121.5, -118.]
        }
    } 
}

def init_proj(date_str):
    #mapcrs = Miller()
    if(date_str == None):
        mapcrs = ccrs.LambertConformal()
    else:
        dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

        mapcrs = ccrs.LambertConformal(central_longitude = \
            np.mean(plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['Lon']),\
            central_latitude = \
            np.mean(plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['Lat']))

    return mapcrs

def plot_trend_line(pax, xdata, ydata, color='black', linestyle = '-', \
        slope = 'thiel-sen'):

    if(slope == 'thiel-sen'):
        res = stats.theilslopes(ydata, xdata, 0.95)
        print("Theil-Sen: {0}x + {1}".format(res[0], res[1]))

        # Then, plot the trend line on the figure
        pax.plot(xdata, res[1] + res[0] * xdata, \
            color='k', linewidth = 2.5, linestyle = linestyle)
        # Then, plot the trend line on the figure
        pax.plot(xdata, res[1] + res[0] * xdata, \
            color=color, linestyle = linestyle)
    else:
        # First, calculate the trend
        zdata = np.polyfit(xdata, ydata, 1)

        print("{0}x + {1}".format(*zdata))

        # Then, plot the trend line on the figure
        pax.plot(np.unique(xdata), np.poly1d(zdata)(np.unique(xdata)), \
            color=color, linestyle = linestyle)

def plot_subplot_label(ax, label, xval = None, yval = None, transform = None, \
        color = 'black', backgroundcolor = None, fontsize = 14, \
        location = 'upper_left'):

    if(location == 'upper_left'):
        y_lim = 0.90
        x_lim = 0.05
    elif(location == 'lower_left'):
        y_lim = 0.05
        x_lim = 0.05
    elif(location == 'upper_right'):
        y_lim = 0.90
        x_lim = 0.90
    elif(location == 'lower_right'):
        y_lim = 0.05
        x_lim = 0.90

    if(xval is None):
        xval = ax.get_xlim()[0] + (ax.get_xlim()[1] - ax.get_xlim()[0]) * x_lim
    if(yval is None):
        yval = ax.get_ylim()[0] + (ax.get_ylim()[1] - ax.get_ylim()[0]) * y_lim
    print('Xval = ',xval, 'Yval = ',yval)

    if(transform is None):
        if(backgroundcolor is None):
            ax.text(xval,yval,label, \
                color=color, weight='bold', \
                fontsize=fontsize)
        else:
            ax.text(xval,yval,label, \
                color=color, weight='bold', \
                fontsize=fontsize, backgroundcolor = backgroundcolor)
    else:
        if(backgroundcolor is None):
            ax.text(xval,yval,label, \
                color=color, weight='bold', \
                transform = transform, fontsize=fontsize)
        else:
            ax.text(xval,yval,label, \
                color=color, weight='bold', \
                transform = transform, fontsize=fontsize, \
                backgroundcolor = backgroundcolor)

def plot_figure_text(ax, text, xval = None, yval = None, transform = None, \
        color = 'black', fontsize = 12, backgroundcolor = 'white',\
        halign = 'left'):

    if(xval is None):
        print(len(text))
        xval = ax.get_xlim()[0] + (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.95
    if(yval is None):
        yval = ax.get_ylim()[0] + (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.05
    print('Xval = ',xval, 'Yval = ',yval)

    if(transform is None):
        ax.text(xval,yval,text, \
            color=color, weight='bold', \
            fontsize=fontsize, backgroundcolor = backgroundcolor, \
            horizontalalignment = halign)
    else:
        ax.text(xval,yval,text, \
            color=color, weight='bold', \
            transform = transform, fontsize=fontsize, \
            backgroundcolor = backgroundcolor, \
            horizontalalignment = halign)

def getCorners_1d(centers):
    one = centers[:-1]
    two = centers[1:]
    d1 = (two - one) / 2.0
    one = one - d1
    two = two + d1
    stepOne = np.zeros((centers.shape[0] + 1))
    stepOne[:-2] = one
    stepOne[-2:] = two[-2:]
    return stepOne

    ##!#one = stepOne[:,:-1]
    ##!#two = stepOne[:,1:]
    ##!#d2 = (two - one) / 2.
    ##!#one = one - d2
    ##!#two = two + d2
    ##!#stepTwo = np.zeros((centers.shape[0] + 1, centers.shape[1] + 1))
    ##!#stepTwo[:,:-2] = one
    ##!#stepTwo[:,-2:] = two[:,-2:]
    ##!#return stepTwo

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
 
# Extract the GOES information from a given channel at each ob point
# -------------------------------------------------------------------
def nearest_grid_values(GOES_data):
    # Read in the correct ASOS file 
    asos_file = plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['asos']
    df = pd.read_csv(asos_file)
    df['valid'] = pd.to_datetime(df['valid'])
    df = df.set_index('valid')

    # Pull the event time from the plot_limits_dict
    event_date = datetime.strptime(GOES_data['cross_date'], "%Y-%m-%d")
    first_time = GOES_data['file_time']
    event_dtime = event_date + timedelta(hours = int(first_time[:2]))
    # Test what happens if looking at the data near the hour
    #event_dtime = event_date + timedelta(hours = int(first_time[:2]), \
    #    minutes = int(first_time[2:4]))
    if(debug):
        print("Overpass time",event_dtime)

    begin_range = event_dtime - timedelta(minutes = 15)
    end_range   = event_dtime + timedelta(minutes = 15)

    compare_dict = {}
    station_names = sorted(set(df['station'].values))
    compare_dict['goes_time'] = event_dtime.strftime('%Y%m%d%H%M')
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
        m_idx = nearest_gridpoint(lat_stn, lon_stn, GOES_data['lat'], \
            GOES_data['lon'])
        m_data = GOES_data['data'][m_idx][0]
      
        compare_dict['mds_data'][ii] = m_data
 
        ##print(station, lat_stn, GOES_data['lat'][m_idx][0], \
        ##    lon_stn, GOES_data['lon'][m_idx][0], compare_dict['stn_data'][ii], m_data )

    return compare_dict

 
# Plot the downloaded ASOS stations for each case
#
# NOTE: if you want to plot ASOS stations from a specific file, and not
# the default file associated with a YYYYMMDDHHMM case, pass in the
# ASOS file name in the 'date_str' argument
# ---------------------------------------------------------------------
def plot_ASOS_locs(pax,date_str,crs = datacrs, color='red', \
        sites = None):
    if(len(date_str.split('.')) > 1):
        asos_file = date_str
    else:
        dt_date_str = datetime.strptime(date_str,'%Y%m%d%H%M')
        cross_date = dt_date_str.strftime('%Y-%m-%d')
        file_date  = dt_date_str.strftime('%H%M')
        asos_file = plot_limits_dict[cross_date][file_date]['asos']

    # Read in the correct ASOS file 
    df = pd.read_csv(asos_file)

    if(color == 'default'):
        color = None

    station_names = set(df['station'].values)
    if(sites == None):
        sites = station_names
    for ii, station in enumerate(station_names):
        if(station in sites):
            lat_stn = df['lat'][df['station'] == station].values[0]
            lon_stn = df['lon'][df['station'] == station].values[0]

            pax.plot(lon_stn, lat_stn,
                     color=color, linewidth=2, marker='o',
                     transform=crs,
                     )

            pax.text(lon_stn + 0.1, lat_stn + 0.1, station, fontsize=10, \
                weight='bold', transform=crs, color=color, backgroundcolor = 'white')

# Determine areas of an image that are in smoke, defined by:
#    (ch1_refl - ch5_refl) < mean(ch1_refl - ch5_refl) * 0.25
def find_plume(dt_date_str):
    GOES_ch1  = read_GOES_channel(dt_date_str, 1,  zoom = True)
    GOES_ch5  = read_GOES_channel(dt_date_str, 5,  zoom = True)
    GOES_ch31 = read_GOES_channel(dt_date_str, 31, zoom = True)

    screen_limit = 0.05
    max_ch = 350.
    test_data = GOES_ch1['data'] - GOES_ch5['data']
    hash_data   = np.ma.masked_where(test_data <  \
        (np.nanmean(test_data) * screen_limit), GOES_ch1['data'])
    hash_data = np.ma.masked_where(GOES_ch5['data'] < 0.05, hash_data)
    nohash_data = np.ma.masked_where(test_data >= \
        (np.nanmean(test_data) * screen_limit), GOES_ch1['data'])

    hash_data = np.ma.masked_where(GOES_ch31['data'] > max_ch, hash_data)
    nohash_data = np.ma.masked_where(GOES_ch31['data'] > max_ch, nohash_data)

    return hash_data, nohash_data

# Start_date and end_date must be formatted as 
# "YYYYMMDD"
# Variable must be either "CHL" or "PIC"
def read_GOES(variable,start_date,latmin = 60):

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
 
    
    # Grab all the goes files
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
    #cmnd = "ls /home/bsorenson/data/GOES/CHL/A*.nc"
    cmnd = "ls /home/bsorenson/data/GOES/"+variable+"/A*.nc"
    
    
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
    goes_data = {}
    goes_data['data'] = np.full((num_files,len(lat_ranges),len(lon_ranges)),-9.)
    goes_data['lat']  = lat_ranges
    goes_data['lon']  = lon_ranges
    goes_data['variable'] = variable
    goes_data['ptitle'] = ptitle
    goes_data['label'] = label
    goes_data['grabber'] = grabber
    goes_data['pticks'] = pticks
    goes_data['ptick_labels'] = ptick_labels
    goes_data['titles']  = []
    
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
                goes_data['data'][count,xi,yj] =\
                    np.average(data.variables[grabber][lat_indices[0],lon_indices[0]])
        data.close() 
        goes_data['titles'].append(fname)
        #goes_data['dates'].append(fname[-11:end_string_idx])
        count+=1

    # Convert the longitude values from 0 - 360 to -180 - 180
    return goes_data

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Trend calculating functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def trend_calc(goes_data,x_ind,y_ind,variable,thielsen=False):
    temp_data = goes_data[variable][:,x_ind,y_ind]
    temp_data[temp_data < -999.] = np.nan
    #temp_data = np.ma.masked_where(temp_data < -999., temp_data)
    # Don't calculate trends if there are less than 2 valid values
    if(np.count_nonzero(~np.isnan(temp_data)) < 2):
        trend = pcnt_change = np.nan
    else:
        avg_goes = np.nanmean(temp_data)
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


        pcnt_change = (trend/avg_goes)*100.
    # Find the percent change per decade
    return trend,pcnt_change


# goes_trendCalc calculates the trends over the time period at each
# grid point on the 25x25 km grid.
def goes_trendCalc(goes_data,thielSen=False):
    # Loop over the data and calculate trends
    for i in range(448):
        print(i)
        max_trend = -99.
        min_trend = 99.
        for j in range(304):
            goes_data['thick_trends'][i,j],temp_pcnt = \
                trend_calc(goes_data,i,j,'ice_thick',thielsen=thielSen)
            #goes_data['thick_land_trends'][i,j],temp_pcnt = \
            #    trend_calc(goes_data,i,j,'ice_thick',thielsen=thielSen)
            goes_data['con_trends'][i,j],temp_pcnt = \
                trend_calc(goes_data,i,j,'ice_con',thielsen=thielSen)
            #goes_data['con_land_trends'][i,j],temp_pcnt = \
            #    trend_calc(goes_data,i,j,'ice_con',thielsen=thielSen)
            ##temp_trend = goes_data['trends'][i,j]
            ##if(temp_trend>max_trend):
            ##    max_trend = temp_trend
            ##if(temp_trend<min_trend):
            ##    min_trend = temp_trend
        ##print("  max trend = ",max_trend)
        ##print("  min trend = ",min_trend)
        # Deal with land masks
        #good_indices = np.where(goes_data['data'][0,i,:]<251.)
        #land_indices = np.where(goes_data['data'][0,i,:]>=251.)
        #goes_data['trends'][i,land_indices] = np.nan
        #goes_data['land_trends'][i,good_indices] = np.nan

    return goes_data

# Calculate trends on the 1x1 degree lat/lon grid
def goes_gridtrendCalc(goes_data,area=True,thielSen=False):
    if(area==True):
        print("\nCalculating area trends\n")
    else:
        print("\nCalculating % concentration trends\n")
    goes_data['month_fix'] = '_monthfix'

    #lat_ranges = np.arange(minlat,90.5,1.0)
    #lon_ranges = np.arange(0.5,360.5,1.0)
    lat_ranges = goes_data['grid_lat']
    lon_ranges = goes_data['grid_lon']
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
                interp_data = (goes_data['grid_goes_conc'][:,xi,yj]/100.)*goes_data['grid_total_area'][xi,yj]
                good_indices = np.where(np.isnan(interp_data)==False)
            else:
                interp_data = goes_data['grid_goes_conc'][:,xi,yj]
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
                goes_data['grid_goes_area_trend'][xi,yj] = total_trend
        print(xi)
        #print("max trend = ",max_trend)
        #print("min trend = ",min_trend)

    return goes_data

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Gridding functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# grid_data_conc grids the 25x25 km gridded goes concentration data into
# a 1x1 degree lat/lon grid
def grid_data_values(goes_data):
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(30.,90.,1.0)
    grid_con    = np.full((len(goes_data['ice_con'][:,0,0]),len(lat_ranges),len(lon_ranges)),-99.)
    grid_thick  = np.full((len(goes_data['ice_thick'][:,0,0]),len(lat_ranges),len(lon_ranges)),-99.)
    grid_con_cc = np.full((len(goes_data['ice_con'][:,0,0]),len(lat_ranges),len(lon_ranges)),-999.)
    grid_thick_cc = np.full((len(goes_data['ice_thick'][:,0,0]),len(lat_ranges),len(lon_ranges)),-999.)
    print("Size of grid array: ",grid_con.shape)
    for nt in range(len(goes_data['ice_con'][:,0,0])):
        print(nt)
        for xi in range(448):
            # Don't grid the data if any portion of the lat/lon box is over land.
            # Don't include land data
            for yj in range(304):
                lat_index = np.where(np.floor(goes_data['lat'][xi,yj])>=lat_ranges)[-1][-1]
                lon_index = np.where(np.floor(goes_data['lon'][xi,yj])>=lon_ranges)[-1][-1]
                # Add the current pixel area into the correct grid box, no
                # matter if the current box is missing or not.
                #if((lat_index==20) & (lon_index==10)):
                #    print("Current grid area = ",grid_goes_area[lat_index,lon_index])
                #if((lat_index==20) & (lon_index==10)):
                #    print("    New grid area = ",grid_goes_area[lat_index,lon_index])
                if(goes_data['ice_con'][nt,xi,yj] != -9999.0):
                    #if(nt==0): grid_goes_area[lat_index,lon_index] += goes_data['area'][xi,yj]
                    if(grid_con_cc[nt,lat_index,lon_index]==-999.):
                        grid_con[nt,lat_index,lon_index] = goes_data['ice_con'][nt,xi,yj]
                        grid_con_cc[nt,lat_index,lon_index] = 1.
                    else:
                        grid_con[nt,lat_index,lon_index] = ((grid_con[nt,lat_index,\
                            lon_index]*grid_con_cc[nt,lat_index,lon_index])+\
                            goes_data['ice_con'][nt,xi,yj])/(grid_con_cc[nt,lat_index,\
                            lon_index]+1.)
                        grid_con_cc[nt,lat_index,lon_index]+=1
                if(goes_data['ice_thick'][nt,xi,yj] != -9999.0):
                    #if(nt==0): grid_goes_area[lat_index,lon_index] += goes_data['area'][xi,yj]
                    if(grid_thick_cc[nt,lat_index,lon_index]==-999.):
                        grid_thick[nt,lat_index,lon_index] = goes_data['ice_thick'][nt,xi,yj]
                        grid_thick_cc[nt,lat_index,lon_index] = 1.
                    else:
                        grid_thick[nt,lat_index,lon_index] = ((grid_thick[nt,lat_index,\
                            lon_index]*grid_thick_cc[nt,lat_index,lon_index])+\
                            goes_data['ice_thick'][nt,xi,yj])/(grid_thick_cc[nt,lat_index,\
                            lon_index]+1.)
                        grid_thick_cc[nt,lat_index,lon_index]+=1
                #else:
                #    if(nt==0): grid_goes_area[lat_index,lon_index] = np.nan

                    # end else
                # end if good goes check
            # end y grid loop
        # end x grid loop
    # end time loop 
    goes_data['grid_con']   = grid_con
    goes_data['grid_thick'] = grid_thick
    goes_data['grid_con_cc'] = grid_con_cc
    goes_data['grid_thick_cc'] = grid_thick_cc
    goes_data['grid_lat'] = lat_ranges
    goes_data['grid_lon'] = lon_ranges
    
    return goes_data

# grid_data averages the 25x25 km trends into a 1x1 degree grid
# The 'grid_goes' paramter in the dictionary therefore contains
# the gridded trends.
def grid_data_trends(goes_dict):
    lon_ranges  = np.arange(-180.,180.,1.0)
    lat_ranges  = np.arange(30.,90.,1.0)
    grid_goes    = np.full((len(lat_ranges),len(lon_ranges)),-999.)
    grid_goes_cc = np.full((len(lat_ranges),len(lon_ranges)),-999.)
    for xi in range(448):
        print(xi)
        for yj in range(304):
            if(not np.isnan(goes_dict['trends'][xi,yj])):
                lat_index = np.where(np.floor(goes_dict['lat'][xi,yj])>=lat_ranges)[-1][-1]
                lon_index = np.where(np.floor(goes_dict['lon'][xi,yj])>=lon_ranges)[-1][-1]
                if(grid_goes_cc[lat_index,lon_index]==-999.):
                    grid_goes[lat_index,lon_index] = goes_dict['trends'][xi,yj]
                    grid_goes_cc[lat_index,lon_index] = 1.
                else:
                    grid_goes[lat_index,lon_index] = ((grid_goes[lat_index,lon_index]*grid_goes_cc[lat_index,lon_index])+\
                                                     goes_dict['trends'][xi,yj])/(grid_goes_cc[lat_index,lon_index]+1.)
                    grid_goes_cc[lat_index,lon_index]+=1
    goes_dict['grid_goes'] = grid_goes
    goes_dict['grid_goes_cc'] = grid_goes_cc
    goes_dict['grid_lat'] = lat_ranges
    goes_dict['grid_lon'] = lon_ranges
    return goes_dict

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Plotting functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Channel must be:
# - 0.64  (band 2, Red)
# - 2.25  (band 6, Cloud particle size)
# - 3.90  (band 7, Shortwave Window)
# - 6.18  (band 8, Upper-Level Water Vapor)
# - 6.95  (band 9, Mid-Level Water Vapor)
# - 7.34  (band 10, Lower-Level Water Vapor)
# - 10.35 (band 13, Clean IR Longwave Window)
def read_GOES_satpy(date_str, channel):

    # Extract the channel wavelength using the input string
    # -----------------------------------------------------
    channel = channel_dict[str(channel)]['wavelength']

    # Determine the correct GOES files associated with the date
    # ---------------------------------------------------------
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    dt_date_str_end = dt_date_str + timedelta(minutes = 10)

    # Use the Satpy find_files_and_readers to grab the files
    # ------------------------------------------------------
    files = find_files_and_readers(start_time = dt_date_str, \
        end_time = dt_date_str_end, base_dir = data_dir, reader = 'abi_l1b')

    print(files)

    # Extract the goes true-color plot limits
    # ----------------------------------------
    lat_lims = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]['goes_Lat']
    lon_lims = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]['goes_Lon']

    # Use satpy (Scene) to open the file
    # ----------------------------------
    scn = Scene(reader = 'abi_l1b', filenames = files)

    # Load the desired channel data
    # -----------------------------
    scn.load([channel])

    ## Set the map projection and center the data
    ## ------------------------------------------
    #my_area = scn[channel].attrs['area'].compute_optimal_bb_area({\
    #    'proj':'lcc', 'lon_0': lon_lims[0], 'lat_0': lat_lims[0], \
    #    'lat_1': lat_lims[0], 'lat_2': lat_lims[0]})
    #new_scn = scn.resample(my_area)

    # Enhance the image for plotting
    # ------------------------------
    var = get_enhanced_image(scn[channel]).data
    var = var.transpose('y','x','bands')

    # Extract the map projection from the data for plotting
    # -----------------------------------------------------
    crs = scn[channel].attrs['area'].to_cartopy_crs()

    return var, crs, lat_lims, lon_lims

# channel must be an integer between 1 and 16
def plot_GOES_satpy(date_str, channel, ax = None, zoom=True,save=False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    var, crs, lat_lims, lon_lims = read_GOES_satpy(date_str, channel)

    # Plot the GOES data
    # ------------------
    in_ax = True 
    if(ax is None):
        in_ax = False
        plt.close('all')
        ax = plt.axes(projection=crs)

    ax.imshow(var.data, transform = crs, extent=(var.x[0], var.x[-1], \
        var.y[-1], var.y[0]), origin='upper', cmap = 'Greys_r')

    # Zoom in the figure if desired
    # -----------------------------
    if(zoom):
        ax.set_extent([lon_lims[0],lon_lims[1],lat_lims[0],lat_lims[1]],\
                       crs = ccrs.PlateCarree())
        zoom_add = '_zoom'
    else:
        zoom_add = ''

    ax.coastlines(resolution = '50m')
    ax.add_feature(cfeature.STATES)
    ax.add_feature(cfeature.BORDERS)
    ax.set_title('GOES-17\n'+dt_date_str.strftime('%Y-%m-%d %H:%M'))

    if(not in_ax): 
        print('here')
        if(save):
            outname = 'goes__' + date_str + zoom_add + cmpst_add + '.png'
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

def plot_GOES_satpy_6panel(date_str, ch1, ch2, ch3, ch4, ch5, ch6, \
        zoom = True, save = False):
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    plt.close('all')
    fig1 = plt.figure(figsize = (13,8))
    var0, crs0, lat_lims, lon_lims = read_GOES_satpy(date_str, ch1)

    var1, crs1, lat_lims, lon_lims = read_GOES_satpy(date_str, ch2)
    var2, crs2, lat_lims, lon_lims = read_GOES_satpy(date_str, ch3)
    var3, crs3, lat_lims, lon_lims = read_GOES_satpy(date_str, ch4)
    var4, crs4, lat_lims, lon_lims = read_GOES_satpy(date_str, ch5)
    var5, crs5, lat_lims, lon_lims = read_GOES_satpy(date_str, ch6)

    ax0 = fig1.add_subplot(2,3,1, projection = crs0)
    ax1 = fig1.add_subplot(2,3,2, projection = crs1)
    ax2 = fig1.add_subplot(2,3,3, projection = crs2)
    ax3 = fig1.add_subplot(2,3,4, projection = crs3)
    ax4 = fig1.add_subplot(2,3,5, projection = crs4)
    ax5 = fig1.add_subplot(2,3,6, projection = crs5)

    max0 = 1.0
    max1 = 1.0
    max2 = 1.0 
    max3 = 1.0 
    max4 = 1.0 
    max5 = 1.0 
    min0 = 0.0
    min1 = 0.0
    min2 = 0.5
    min3 = None
    min4 = None
    min5 = None

    v3lon, v3lat = var3.attrs['area'].get_lonlats()
    v2data = np.squeeze(var2.values)
    v3data = np.squeeze(var3.values)
    v4data = np.squeeze(var4.values)
    v5data = np.squeeze(var5.values)
    v2data = np.ma.masked_where( ~((v3lat >= lat_lims[0] - 0.5) & (v3lat <= lat_lims[1] + 0.5) &\
           (v3lon >= lon_lims[0] - 1.) & (v3lon <= lon_lims[1] + 1.)), v2data)
    v3data = np.ma.masked_where( ~((v3lat >= lat_lims[0] - 0.5) & (v3lat <= lat_lims[1] + 0.5) &\
           (v3lon >= lon_lims[0] - 1.) & (v3lon <= lon_lims[1] + 1.)), v3data)
    v4data = np.ma.masked_where( ~((v3lat >= lat_lims[0] - 0.5) & (v3lat <= lat_lims[1] + 0.5) &\
           (v3lon >= lon_lims[0] - 1.) & (v3lon <= lon_lims[1] + 1.)), v4data)
    v5data = np.ma.masked_where( ~((v3lat >= lat_lims[0] - 0.5) & (v3lat <= lat_lims[1] + 0.5) &\
           (v3lon >= lon_lims[0] - 1.) & (v3lon <= lon_lims[1] + 1.)), v5data)
    #my_area = scn[channel].attrs['area'].compute_optimal_bb_area({\
    #    'proj':'lcc', 'lon_0': lon_lims[0], 'lat_0': lat_lims[0], \
    #    'lat_1': lat_lims[0], 'lat_2': lat_lims[0]})

    mesh0 = ax0.imshow(var0.data, transform = crs0, extent=(var0.x[0], var0.x[-1], \
        var0.y[-1], var0.y[0]), origin='upper', cmap = 'Greys_r', \
        vmin = min0, vmax = max0)
    mesh1 = ax1.imshow(var1.data, transform = crs1, extent=(var1.x[0], var1.x[-1], \
        var1.y[-1], var1.y[0]), origin='upper', cmap = 'Greys_r', \
        vmin = min1, vmax = max1)
    mesh2 = ax2.imshow(v2data, transform = crs2, extent=(var2.x[0], var2.x[-1], \
        var2.y[-1], var2.y[0]), origin='upper', cmap = 'Greys_r', \
        vmin = min2, vmax = max2)
    mesh3 = ax3.imshow(v3data, transform = crs3, extent=(var3.x[0], var3.x[-1], \
        var3.y[-1], var3.y[0]), origin='upper', cmap = 'Greys_r', \
        vmin = min3, vmax = max3)
    mesh4 = ax4.imshow(v4data, transform = crs4, extent=(var4.x[0], var4.x[-1], \
        var4.y[-1], var4.y[0]), origin='upper', cmap = 'Greys_r', \
        vmin = min4, vmax = max4)
    mesh5 = ax5.imshow(v5data, transform = crs5, extent=(var5.x[0], var5.x[-1], \
        var5.y[-1], var5.y[0]), origin='upper', cmap = 'Greys_r', \
        vmin = min5, vmax = max5)
    cbar0 = plt.colorbar(mesh0,ax=ax0,orientation='vertical',\
        pad=0.03, shrink = 0.67)
    cbar1 = plt.colorbar(mesh1,ax=ax1,orientation='vertical',\
        pad=0.03, shrink = 0.67)
    cbar2 = plt.colorbar(mesh2,ax=ax2,orientation='vertical',\
        pad=0.03, shrink = 0.67)
    cbar3 = plt.colorbar(mesh3,ax=ax3,orientation='vertical',\
        pad=0.03, shrink = 0.67)
    cbar4 = plt.colorbar(mesh4,ax=ax4,orientation='vertical',\
        pad=0.03, shrink = 0.67)
    cbar5 = plt.colorbar(mesh5,ax=ax5,orientation='vertical',\
        pad=0.03, shrink = 0.67)

    # Zoom in the figure if desired
    # -----------------------------
    if(zoom):
        ax0.set_extent([lon_lims[0],lon_lims[1],lat_lims[0],lat_lims[1]],\
                       crs = ccrs.PlateCarree())
        ax1.set_extent([lon_lims[0],lon_lims[1],lat_lims[0],lat_lims[1]],\
                       crs = ccrs.PlateCarree())
        ax2.set_extent([lon_lims[0],lon_lims[1],lat_lims[0],lat_lims[1]],\
                       crs = ccrs.PlateCarree())
        ax3.set_extent([lon_lims[0],lon_lims[1],lat_lims[0],lat_lims[1]],\
                       crs = ccrs.PlateCarree())
        ax4.set_extent([lon_lims[0],lon_lims[1],lat_lims[0],lat_lims[1]],\
                       crs = ccrs.PlateCarree())
        ax5.set_extent([lon_lims[0],lon_lims[1],lat_lims[0],lat_lims[1]],\
                       crs = ccrs.PlateCarree())
        zoom_add = '_zoom'
    else:
        zoom_add = ''

    ax0.coastlines(resolution = '50m')
    ax0.add_feature(cfeature.STATES)
    ax0.add_feature(cfeature.BORDERS)
    ax0.set_title('GOES-17 Band ' + str(ch1) + '\n' + \
        channel_dict[str(ch1)]['name'] + '\n' + \
        dt_date_str.strftime('%Y-%m-%d %H:%M'))

    ax1.coastlines(resolution = '50m')
    ax1.add_feature(cfeature.STATES)
    ax1.add_feature(cfeature.BORDERS)
    ax1.set_title('GOES-17 Band ' + str(ch2) + '\n' + \
        channel_dict[str(ch2)]['name'] + '\n' + \
        dt_date_str.strftime('%Y-%m-%d %H:%M'))

    ax2.coastlines(resolution = '50m')
    ax2.add_feature(cfeature.STATES)
    ax2.add_feature(cfeature.BORDERS)
    ax2.set_title('GOES-17 Band ' + str(ch3) + '\n' + \
        channel_dict[str(ch3)]['name'] + '\n' + \
        dt_date_str.strftime('%Y-%m-%d %H:%M'))

    ax3.coastlines(resolution = '50m')
    ax3.add_feature(cfeature.STATES)
    ax3.add_feature(cfeature.BORDERS)
    ax3.set_title('GOES-17 Band ' + str(ch4) + '\n' + \
        channel_dict[str(ch4)]['name'] + '\n' + \
        dt_date_str.strftime('%Y-%m-%d %H:%M'))

    ax4.coastlines(resolution = '50m')
    ax4.add_feature(cfeature.STATES)
    ax4.add_feature(cfeature.BORDERS)
    ax4.set_title('GOES-17 Band ' + str(ch5) + '\n' + \
        channel_dict[str(ch5)]['name'] + '\n' + \
        dt_date_str.strftime('%Y-%m-%d %H:%M'))

    ax5.coastlines(resolution = '50m')
    ax5.add_feature(cfeature.STATES)
    ax5.add_feature(cfeature.BORDERS)
    ax5.set_title('GOES-17 Band ' + str(ch6) + '\n' + \
        channel_dict[str(ch6)]['name'] + '\n' + \
        dt_date_str.strftime('%Y-%m-%d %H:%M'))

    fig1.tight_layout()

    if(save):
        outname = 'goes17_'+date_str+'_6panel.png'
        fig1.savefig(outname, dpi = 300)
        print('Saved image', outname)
    else:
        plt.show()
    

def plot_true_color(filename,zoom=True):
    goes = SD.SD(filename)

    dat = goes.attributes().get('CoreMetadata.0').split()
    indx = dat.index('EQUATORCROSSINGDATE')+9
    cross_date = dat[indx][1:len(dat[indx])-1]
    cross_time = filename.strip().split('/')[-1].split('.')[2]

    lat5 = goes.select('Latitude').get()
    lon5 = goes.select('Longitude').get()

    red   = goes.select('EV_250_Aggr1km_RefSB').get()[0]
    green = goes.select('EV_500_Aggr1km_RefSB').get()[1]
    blue  = goes.select('EV_500_Aggr1km_RefSB').get()[0]

    red   = red[::5,::5]
    green = green[::5,::5]
    blue  = blue[::5,::5]

    red_scale    = goes.select('EV_250_Aggr1km_RefSB').attributes().get('reflectance_scales')[0]
    red_offset   = goes.select('EV_250_Aggr1km_RefSB').attributes().get('reflectance_offsets')[0]
    green_scale  = goes.select('EV_500_Aggr1km_RefSB').attributes().get('reflectance_scales')[1]
    green_offset = goes.select('EV_500_Aggr1km_RefSB').attributes().get('reflectance_offsets')[1]
    blue_scale   = goes.select('EV_500_Aggr1km_RefSB').attributes().get('reflectance_scales')[0]
    blue_offset  = goes.select('EV_500_Aggr1km_RefSB').attributes().get('reflectance_offsets')[0]

    goes.end()

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
    mapcrs = init_proj()
    ax = plt.axes(projection = mapcrs)

    image = np.ma.masked_where(np.isnan(image),image)
    
    if(debug):
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
def read_GOES_channel(date_str, channel, zoom = False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # Extract the filename given the channel
    # --------------------------------------
    if(str(channel)[:2] == 'wv'):
        filename = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['mdswv']
    else:
        filename = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['goes']
    
    print("Reading GOES channel",channel," from ",filename)

    GOES_data = {}

    goes = SD.SD(filename)

    dat = goes.attributes().get('CoreMetadata.0').split()
    indx = dat.index('EQUATORCROSSINGDATE')+9
    cross_date = dat[indx][1:len(dat[indx])-1]
    indx = dat.index('EQUATORCROSSINGTIME')+9
    cross_time = dat[indx][1:len(dat[indx])-1]

    if(debug):
        print('GOES orbit info',cross_date,cross_time)

    lat5 = goes.select('Latitude').get()
    lon5 = goes.select('Longitude').get()

    if(str(channel)[:2] != 'wv'):
        sza = goes.select('SolarZenith').get() * goes.select('SolarZenith').attributes().get('scale_factor')
        vza = goes.select('SensorZenith').get() * goes.select('SensorZenith').attributes().get('scale_factor')

    #print('Sensor zenith: ', goes.select('SensorZenith').get().shape)
    #print('Sensor zenith scale: ', goes.select('SensorZenith').attributes().get('scale_factor'))
    #print('Solar zenith: ', goes.select('SolarZenith').get().shape)
    #print('Solar zenith scale: ', goes.select('SolarZenith').attributes().get('scale_factor'))

    data  = goes.select(channel_dict[str(channel)]['Name']).get()
    if(str(channel)[2:] != '_ir'):
        if(str(channel) != 'wv_nir'):
            data = data[channel_dict[str(channel)]['Index']]
        data  = data[::5,::5]

        if(data.shape != lat5.shape):
            data = data[:lat5.shape[0], :lat5.shape[1]]


    if(str(channel)[:2] == 'wv'):
        data_scale    = goes.select(channel_dict[str(channel)]['Name']\
            ).attributes().get('scale_factor')
        data_offset   = goes.select(channel_dict[str(channel)]['Name']\
            ).attributes().get('add_offset')
        # Extract the fill value and make sure any missing values are
        # removed.
        # -----------------------------------------------------------
        mask_val = goes.select(channel_dict[str(channel)]['Name']\
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
        label = 'IR Precipitable water vapor [cm]'
        ##!#try:
        ##!#    label = goes.select(channel_dict[str(channel)]['Name']\
        ##!#        ).attributes().get('long_name') +  ' [' + \
        ##!#        goes.select(channel_dict[str(channel)]['Name']\
        ##!#        ).attributes().get('units') + ']'
        ##!#except:
        ##!#    label = goes.select(channel_dict[str(channel)]['Name']\
        ##!#        ).attributes().get('long_name') +  ' [' + \
        ##!#        goes.select(channel_dict[str(channel)]['Name']\
        ##!#        ).attributes().get('unit') + ']'
    else:
        channel = int(channel)
        # Thermal emission data
        if((channel >= 20) & (channel != 26)):
            data_scale    = goes.select(channel_dict[str(channel)]['Name']\
                ).attributes().get('radiance_scales')[channel_dict[str(channel)]['Index']]
            data_offset   = goes.select(channel_dict[str(channel)]['Name']\
                ).attributes().get('radiance_offsets')[channel_dict[str(channel)]['Index']]

            # Extract the fill value and make sure any missing values are
            # removed.
            # -----------------------------------------------------------
            mask_val = goes.select(channel_dict[str(channel)]['Name']\
                ).attributes().get('_FillValue')
            bad_locations = np.where(data == mask_val)

            data = ((data - data_offset) * data_scale)

            # Define constants for converting radiances to temperatures
            lmbda = (1e-6) * (np.average(channel_dict[str(channel)]['wavelength'])) # in m
            if(debug):
                print("Average wavelength = ",np.average(channel_dict[str(channel)]['wavelength']))
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
            label = 'Brightness Temperature [K]'

        # Reflectances
        else:
            data_scale    = goes.select(channel_dict[str(channel)]['Name']\
                ).attributes().get('reflectance_scales')[channel_dict[str(channel)]['Index']]
            data_offset   = goes.select(channel_dict[str(channel)]['Name']\
                ).attributes().get('reflectance_offsets')[channel_dict[str(channel)]['Index']]

            # Extract the fill value and make sure any missing values are
            # removed.
            # -----------------------------------------------------------
            mask_val = goes.select(channel_dict[str(channel)]['Name']\
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
 
    goes.end()

    # Mask any fire pixels
    # --------------------
    data = np.ma.masked_where(data > 340, data)

    GOES_data['data'] = data
    GOES_data['lat']  = lat5
    GOES_data['lon']  = lon5
    if(str(channel)[:2] != 'wv'):
        GOES_data['sza']  = sza
        GOES_data['vza']  = vza
    GOES_data['variable']  = label
    GOES_data['cross_date']  = cross_date
    GOES_data['channel']  = channel
    GOES_data['colors']  = colors
    GOES_data['file_time'] = filename.strip().split('/')[-1].split('.')[2]

    if(zoom):
        # Mask GOES_data['data'] that are outside the desired range
        # --------------------------------------------
        GOES_data['data'][(((GOES_data['lat'] < plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['Lat'][0]) | \
                             (GOES_data['lat'] > plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['Lat'][1])) | \
                            ((GOES_data['lon'] < plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['Lon'][0]) | \
                             (GOES_data['lon'] > plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['Lon'][1])))] = -999.
        GOES_data['lat'][ (((GOES_data['lat'] < plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['Lat'][0]) | \
                             (GOES_data['lat'] > plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['Lat'][1])) | \
                            ((GOES_data['lon'] < plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['Lon'][0]) | \
                             (GOES_data['lon'] > plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['Lon'][1])))] = -999.
        GOES_data['lon'][ (((GOES_data['lat'] < plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['Lat'][0]) | \
                             (GOES_data['lat'] > plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['Lat'][1])) | \
                            ((GOES_data['lon'] < plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['Lon'][0]) | \
                             (GOES_data['lon'] > plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['Lon'][1])))] = -999.

        GOES_data['data'] = np.ma.masked_where(GOES_data['data'] == -999., GOES_data['data'])
        GOES_data['lat'] = np.ma.masked_where(GOES_data['lat'] == -999., GOES_data['lat'])
        GOES_data['lon'] = np.ma.masked_where(GOES_data['lon'] == -999., GOES_data['lon'])


    return GOES_data

def plot_GOES_channel(date_str,channel,zoom=True,show_smoke=False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    filename = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['goes']

    if(channel == 'red'):
        channel = 1
    elif(channel == 'green'):
        channel = 4
    elif(channel == 'blue'):
        channel = 3

    # Call read_GOES_channel to read the desired GOES data from the
    # file and put it in a dictionary
    # ---------------------------------------------------------------
    GOES_data = read_GOES_channel(date_str, channel)

    if(debug):
        print("Data max = ",np.max(GOES_data['data']), "  Data min = ",\
        np.min(GOES_data['data']))

    plt.close('all')
    fig1 = plt.figure()
    datacrs = ccrs.PlateCarree()
    #mapcrs = ccrs.Miller()
    mapcrs = init_proj(date_str)
    ax = plt.axes(projection = mapcrs)

    plot_GOES_spatial(GOES_data, ax, zoom)

    #mesh = ax.pcolormesh(GOES_data['lon'],GOES_data['lat'],\
    #    GOES_data['data'],cmap = GOES_data['colors'], shading='auto', \
    #    transform = datacrs) 

    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(date_str) 
        #hash_data1, nohash_data1 = find_plume(filename) 
        hash0 = ax.pcolor(GOES_data['lon'],GOES_data['lat'],\
            hash_data1, hatch = '///', alpha=0., transform = datacrs)

    ##!#cbar = plt.colorbar(mesh,orientation='horizontal',pad=0,\
    ##!#    aspect=50,shrink = 0.850)
    ##!#cbar.ax.tick_params(labelsize=14)
    ##!#cbar.set_label(GOES_data['variable'],fontsize=16,weight='bold')
    ##!#
    ##!#ax.add_feature(cfeature.BORDERS)
    ##!#ax.add_feature(cfeature.STATES)
    ##!#ax.coastlines()
    ##!#if(zoom):
    ##!#    ax.set_extent([plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['Lon'][0], \
    ##!#                   plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['Lon'][1], \
    ##!#                   plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['Lat'][0], \
    ##!#                   plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['Lat'][1]],\
    ##!#                   ccrs.PlateCarree())
    ##!#ax.set_title('Channel ' + str(channel) + '\n' + \
    ##!#    channel_dict[str(channel)]['wavelength_label']) 

    plt.show()

def read_OMI_match_GOES(date_str, min_AI = -2e5, corners = False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    goes_date = dt_date_str.strftime('%Y-%m-%d')

    data = h5py.File(plot_limits_dict[goes_date][date_str[8:]]['omi'],'r')
    LAT   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/'+\
        'Latitude'][:,:]
    LON   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/'+\
        'Longitude'][:,:]
    UVAI  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/'+\
        'UVAerosolIndex'][:,:]
    XTRACK = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/'+\
        'XTrackQualityFlags'][:,:]
    if(corners):
        crnr_LAT   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/'+\
            'Geolocation Fields/FoV75CornerLatitude'][:,:,:]
        crnr_LON   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/'+\
            'Geolocation Fields/FoV75CornerLongitude'][:,:,:]

    mask_UVAI = np.ma.masked_where((XTRACK < -2e5) | (UVAI < min_AI), UVAI)
    mask_UVAI = np.ma.masked_where((\
        ((LAT < plot_limits_dict[goes_date][date_str[8:]]['Lat'][0]) | \
         (LAT > plot_limits_dict[goes_date][date_str[8:]]['Lat'][1])) | \
        ((LON < plot_limits_dict[goes_date][date_str[8:]]['Lon'][0]) | \
         (LON > plot_limits_dict[goes_date][date_str[8:]]['Lon'][1]))), \
        mask_UVAI)

    data.close()

    if(corners):
        return LAT, LON, mask_UVAI, crnr_LAT, crnr_LON
    else:
        return LAT, LON, mask_UVAI

def read_CERES_match_GOES(date_str):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # NOTE: need to use the time variable to screen out any data
    # that are not close to the event time.
    ##!#base_date = datetime(year=1970,month=1,day=1)
    ##!#start_date = dt_date_str - timedelta(hours = 1)
    ##!#end_date   = dt_date_str + timedelta(hours = 2)

    ##!#print(start_date, end_date)

    goes_date = dt_date_str.strftime('%Y-%m-%d')

    data = Dataset(plot_limits_dict[goes_date][date_str[8:]]['ceres'],'r')
    LAT   = 90. - data.variables['Colatitude_of_CERES_FOV_at_surface'][:]
    LON   = data.variables['Longitude_of_CERES_FOV_at_surface'][:]
    LON[LON>179.99] = -360.+LON[LON>179.99]
    swflux  = data.variables['CERES_SW_TOA_flux___upwards'][:]
    lwflux  = data.variables['CERES_LW_TOA_flux___upwards'][:]
    time  = data.variables['time'][:]
    #local_time = np.array([base_date + relativedelta(days = ttime) for ttime in time])

    data.close()

    mask_LAT = LAT[ \
        (LAT >= plot_limits_dict[goes_date][date_str[8:]]['Lat'][0]) & \
        (LAT <= plot_limits_dict[goes_date][date_str[8:]]['Lat'][1]) & \
        (LON >= plot_limits_dict[goes_date][date_str[8:]]['Lon'][0]) & \
        (LON <= plot_limits_dict[goes_date][date_str[8:]]['Lon'][1]) & \
        (swflux > 0) & (lwflux > 0)]
    mask_LON = LON[ \
        (LAT >= plot_limits_dict[goes_date][date_str[8:]]['Lat'][0]) & \
        (LAT <= plot_limits_dict[goes_date][date_str[8:]]['Lat'][1]) & \
        (LON >= plot_limits_dict[goes_date][date_str[8:]]['Lon'][0]) & \
        (LON <= plot_limits_dict[goes_date][date_str[8:]]['Lon'][1]) & \
        (swflux > 0) & (lwflux > 0)]
    mask_swf = swflux[ \
        (LAT >= plot_limits_dict[goes_date][date_str[8:]]['Lat'][0]) & \
        (LAT <= plot_limits_dict[goes_date][date_str[8:]]['Lat'][1]) & \
        (LON >= plot_limits_dict[goes_date][date_str[8:]]['Lon'][0]) & \
        (LON <= plot_limits_dict[goes_date][date_str[8:]]['Lon'][1]) & \
        (swflux > 0) & (lwflux > 0)]
    mask_lwf = lwflux[ \
        (LAT >= plot_limits_dict[goes_date][date_str[8:]]['Lat'][0]) & \
        (LAT <= plot_limits_dict[goes_date][date_str[8:]]['Lat'][1]) & \
        (LON >= plot_limits_dict[goes_date][date_str[8:]]['Lon'][0]) & \
        (LON <= plot_limits_dict[goes_date][date_str[8:]]['Lon'][1]) & \
        (swflux > 0) & (lwflux > 0)]
    ##!#mask_time = local_time[ \
    ##!#    (LAT >= plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['Lat'][0]) & \
    ##!#    (LAT <= plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['Lat'][1]) & \
    ##!#    (LON >= plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['Lon'][0]) & \
    ##!#    (LON <= plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['Lon'][1]) & \
    ##!#    (swflux > 0) & (lwflux > 0)]

    return mask_LAT, mask_LON, mask_swf, mask_lwf

# date_str of format "YYYYMMDDHHMM"
def compare_GOES_3panel(date_str,channel1,channel2,channel3,zoom=True,save=False,\
        plot_ASOS_loc = False, show_smoke = True, compare_OMI = False, \
        compare_CERES = False, return_GOES = False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    filename = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['goes']

    if(channel1== 'red'):
        channel1= 1
    elif(channel1== 'green'):
        channel1= 4
    elif(channel1== 'blue'):
        channel1= 3

    # Step 1: Call read_GOES_channel to read the desired GOES data from the
    # file and put it in a dictionary
    # -----------------------------------------------------------------------
    GOES_data1 = read_GOES_channel(dt_date_str.strftime("%Y%m%d%H%M"), channel1, zoom = True)
    GOES_data2 = read_GOES_channel(dt_date_str.strftime("%Y%m%d%H%M"), channel2, zoom = True)
    GOES_data3 = read_GOES_channel(dt_date_str.strftime("%Y%m%d%H%M"), channel3, zoom = True)

    ##!#print("Data1 max = ",np.max(GOES_data1['data']), "  Data1 min = ",np.min(GOES_data1['data']))
    ##!#print("Data2 max = ",np.max(GOES_data2['data']), "  Data2 min = ",np.min(GOES_data2['data']))
    ##!#print("Data3 max = ",np.max(GOES_data3['data']), "  Data3 min = ",np.min(GOES_data3['data']))

    # Create copies to set the minimum and maximums for each plot
    # -----------------------------------------------------------
    cpy_1 = np.copy(GOES_data1['data'])
    cpy_2 = np.copy(GOES_data2['data'])
    cpy_3 = np.copy(GOES_data3['data'])

    cpy_1 = np.ma.masked_where((((GOES_data1['lat'] < plot_limits_dict[GOES_data1['cross_date']][GOES_data1['file_time']]['Lat'][0]) | \
                         (GOES_data1['lat'] > plot_limits_dict[GOES_data1['cross_date']][GOES_data1['file_time']]['Lat'][1])) | \
                        ((GOES_data1['lon'] < plot_limits_dict[GOES_data1['cross_date']][GOES_data1['file_time']]['Lon'][0]) | \
                         (GOES_data1['lon'] > plot_limits_dict[GOES_data1['cross_date']][GOES_data1['file_time']]['Lon'][1]))), cpy_1)
    cpy_2 = np.ma.masked_where((((GOES_data2['lat'] < plot_limits_dict[GOES_data2['cross_date']][GOES_data2['file_time']]['Lat'][0]) | \
                         (GOES_data2['lat'] > plot_limits_dict[GOES_data2['cross_date']][GOES_data2['file_time']]['Lat'][1])) | \
                        ((GOES_data2['lon'] < plot_limits_dict[GOES_data2['cross_date']][GOES_data2['file_time']]['Lon'][0]) | \
                         (GOES_data2['lon'] > plot_limits_dict[GOES_data2['cross_date']][GOES_data2['file_time']]['Lon'][1]))), cpy_2)
    cpy_3 = np.ma.masked_where((((GOES_data3['lat'] < plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][0]) | \
                         (GOES_data3['lat'] > plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][1])) | \
                        ((GOES_data3['lon'] < plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][0]) | \
                         (GOES_data3['lon'] > plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][1]))), cpy_3)


    # Step 2: Set up figure to have 3 panels
    # --------------------------------------
    datacrs = ccrs.PlateCarree() 
    #mapcrs = ccrs.LambertConformal()
    mapcrs = init_proj(date_str)

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
    plot_GOES_spatial(GOES_data1, ax0, ptitle = '', zoom = zoom)
    plot_GOES_spatial(GOES_data2, ax1, ptitle = '', zoom = zoom)
    plot_GOES_spatial(GOES_data3, ax2, ptitle = '', zoom = zoom)

    plot_figure_text(ax2, 'GOES 11 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 15, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax0, 'GOES 0.64 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 15, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax1, 'GOES 1.24 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 15, backgroundcolor = 'white', halign = 'right')
    ##!## Step 3: Plot the GOES channel data in the first 2 panels
    ##!## ---------------------------------------------------------
    ##!#mesh0 = ax0.pcolormesh(GOES_data1['lon'],GOES_data1['lat'],\
    ##!#    cpy_1,cmap = GOES_data1['colors'], shading='auto', \
    ##!#    vmin = np.nanmin(GOES_data1['data']), \
    ##!#    vmax = np.nanmax(GOES_data1['data']), transform = datacrs) 

    ##!#cbar0 = plt.colorbar(mesh0,ax=ax0,orientation='vertical',\
    ##!#    pad=0.03,label=GOES_data1['variable'])

    ##!#ax0.add_feature(cfeature.BORDERS)
    ##!#ax0.add_feature(cfeature.STATES)
    ##!#ax0.coastlines()
    ##!#if(zoom):
    ##!#    ax0.set_extent([plot_limits_dict[GOES_data1['cross_date']][GOES_data1['file_time']]['Lon'][0], \
    ##!#                    plot_limits_dict[GOES_data1['cross_date']][GOES_data1['file_time']]['Lon'][1], \
    ##!#                    plot_limits_dict[GOES_data1['cross_date']][GOES_data1['file_time']]['Lat'][0], \
    ##!#                    plot_limits_dict[GOES_data1['cross_date']][GOES_data1['file_time']]['Lat'][1]],\
    ##!#                    datacrs)
    ##!##ax0.set_title('GOES Ch. ' + str(channel1) + '\n' + \
    ##!##    str(channel_dict[str(channel1)]['wavelength'][0]) + ' μm - ' + \
    ##!##    str(channel_dict[str(channel1)]['wavelength'][1]) + ' μm')
    ##!#ax0.set_title('Channel ' + str(channel1) + '\n' + \
    ##!#    channel_dict[str(channel1)]['wavelength_label']) 


    ##!## Plot channel 2
    ##!#mesh1 = ax1.pcolormesh(GOES_data2['lon'],GOES_data2['lat'],\
    ##!#    GOES_data2['data'],cmap = GOES_data2['colors'], shading='auto', \
    ##!#    vmin = np.nanmin(GOES_data2['data']), \
    ##!#    vmax = np.nanmax(GOES_data2['data']), transform = datacrs) 


    ##!#cbar1 = plt.colorbar(mesh1,ax=ax1,orientation='vertical',\
    ##!#    pad=0.03,label=GOES_data2['variable'])
    ##!#
    ##!#ax1.add_feature(cfeature.BORDERS)
    ##!#ax1.add_feature(cfeature.STATES)
    ##!#ax1.coastlines()
    ##!#if(zoom):
    ##!#    ax1.set_extent([plot_limits_dict[GOES_data2['cross_date']][GOES_data2['file_time']]['Lon'][0], \
    ##!#                    plot_limits_dict[GOES_data2['cross_date']][GOES_data2['file_time']]['Lon'][1], \
    ##!#                    plot_limits_dict[GOES_data2['cross_date']][GOES_data2['file_time']]['Lat'][0], \
    ##!#                    plot_limits_dict[GOES_data2['cross_date']][GOES_data2['file_time']]['Lat'][1]],\
    ##!#                    datacrs)
    ##!##ax1.set_title('GOES Ch. ' + str(channel2) + '\n' + \
    ##!##    str(channel_dict[str(channel2)]['wavelength'][0]) + ' μm - ' + \
    ##!##    str(channel_dict[str(channel2)]['wavelength'][1]) + ' μm')
    ##!#ax1.set_title('Channel ' + str(channel2) + '\n' + \
    ##!#    channel_dict[str(channel2)]['wavelength_label']) 

    ##!## Plot channel 3
    ##!#mesh2 = ax2.pcolormesh(GOES_data3['lon'],GOES_data3['lat'],\
    ##!#    GOES_data3['data'],cmap = GOES_data3['colors'], shading='auto', \
    ##!#    vmin = np.nanmin(GOES_data3['data']), \
    ##!#    vmax = np.nanmax(GOES_data3['data']), transform = datacrs) 
    ##!#cbar2 = plt.colorbar(mesh2,ax=ax2,orientation='vertical',\
    ##!#    pad=0.03,label=GOES_data3['variable'])
    ##!#
    ##!#ax2.add_feature(cfeature.BORDERS)
    ##!#ax2.add_feature(cfeature.STATES)
    ##!#ax2.coastlines()
    ##!#if(zoom):
    ##!#    ax2.set_extent([plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][0], \
    ##!#                    plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][1], \
    ##!#                    plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][0], \
    ##!#                    plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][1]],\
    ##!#                    datacrs)
    ##!##ax2.set_title('GOES Ch. ' + str(channel3) + '\n' + \
    ##!##    str(channel_dict[str(channel3)]['wavelength'][0]) + ' μm - ' + \
    ##!##    str(channel_dict[str(channel3)]['wavelength'][1]) + ' μm')
    ##!#ax2.set_title('Channel ' + str(channel3) + '\n' + \
    ##!#    channel_dict[str(channel3)]['wavelength_label']) 


    if(compare_OMI):
        print("Reading OMI data")
        LAT, LON, mask_UVAI = read_OMI_match_GOES(date_str)
        ##!#data = h5py.File(plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['omi'],'r')
        ##!#LAT   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude'][:,:]
        ##!#LON   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude'][:,:]
        ##!#UVAI  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex'][:,:]
        ##!#XTRACK = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags'][:,:]
        ##!#mask_UVAI = np.ma.masked_where((XTRACK < -2e5) | (UVAI < -2e5), UVAI)
        ##!#mask_UVAI = np.ma.masked_where((((LAT < plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][0]) | \
        ##!#                     (LAT > plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][1])) | \
        ##!#                    ((LON < plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][0]) | \
        ##!#                     (LON > plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][1]))), mask_UVAI)

        plot_OMI_spatial(date_str, LAT, LON, mask_UVAI, axo, zoom = zoom)

        ##!#mesh3 = axo.pcolormesh(LON,LAT,mask_UVAI, cmap = 'plasma', shading='auto', \
        ##!#    vmin = np.nanmin(mask_UVAI), vmax = np.nanmax(mask_UVAI), transform = datacrs) 
        ##!#axo.add_feature(cfeature.BORDERS)
        ##!#axo.add_feature(cfeature.STATES)
        ##!#axo.coastlines()
        ##!#if(zoom):
        ##!#    axo.set_extent([plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][0], \
        ##!#                    plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][1], \
        ##!#                    plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][0], \
        ##!#                    plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][1]],\
        ##!#                    datacrs)
        ##!#cbar3 = plt.colorbar(mesh3,ax=axo,orientation='vertical',\
        ##!#    pad=0.03,label='OMI UVAI')
        ##!#axo.set_title('OMI UVAI')

    if(compare_CERES):
        print("Reading CERES data")

        mask_LAT, mask_LON, mask_swf, mask_lwf = read_CERES_match_GOES(date_str)

        ##!#base_date = datetime(year=1970,month=1,day=1)
        ##!#start_date = dt_date_str - timedelta(hours = 1)
        ##!#end_date   = dt_date_str + timedelta(hours = 2)

        ##!#print(start_date, end_date)

        ##!#data = Dataset(plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['ceres'],'r')
        ##!#LAT   = 90. - data.variables['Colatitude_of_CERES_FOV_at_surface'][:]
        ##!#LON   = data.variables['Longitude_of_CERES_FOV_at_surface'][:]
        ##!#LON[LON>179.99] = -360.+LON[LON>179.99]
        ##!#swflux  = data.variables['CERES_SW_TOA_flux___upwards'][:]
        ##!#lwflux  = data.variables['CERES_LW_TOA_flux___upwards'][:]
        ##!#time  = data.variables['time'][:]
        ##!#local_time = np.array([base_date + relativedelta(days = ttime) for ttime in time])

        ##!#mask_LAT = LAT[ \
        ##!#    (LAT >= plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][0]) & \
        ##!#    (LAT <= plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][1]) & \
        ##!#    (LON >= plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][0]) & \
        ##!#    (LON <= plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][1]) & \
        ##!#    (swflux > 0) & (lwflux > 0)]
        ##!#mask_LON = LON[ \
        ##!#    (LAT >= plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][0]) & \
        ##!#    (LAT <= plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][1]) & \
        ##!#    (LON >= plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][0]) & \
        ##!#    (LON <= plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][1]) & \
        ##!#    (swflux > 0) & (lwflux > 0)]
        ##!#mask_swf = swflux[ \
        ##!#    (LAT >= plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][0]) & \
        ##!#    (LAT <= plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][1]) & \
        ##!#    (LON >= plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][0]) & \
        ##!#    (LON <= plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][1]) & \
        ##!#    (swflux > 0) & (lwflux > 0)]
        ##!#mask_lwf = lwflux[ \
        ##!#    (LAT >= plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][0]) & \
        ##!#    (LAT <= plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][1]) & \
        ##!#    (LON >= plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][0]) & \
        ##!#    (LON <= plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][1]) & \
        ##!#    (swflux > 0) & (lwflux > 0)]
        ##!###!#mask_time = local_time[ \
        ##!###!#    (LAT >= plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][0]) & \
        ##!###!#    (LAT <= plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][1]) & \
        ##!###!#    (LON >= plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][0]) & \
        ##!###!#    (LON <= plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][1]) & \
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
        ##!#print(plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'])
 
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
        ##!#    axcs.set_extent([plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][0], \
        ##!#                     plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][1], \
        ##!#                     plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][0], \
        ##!#                     plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][1]],\
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
        ##!#    axcl.set_extent([plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][0], \
        ##!#                    plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][1], \
        ##!#                    plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][0], \
        ##!#                    plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][1]],\
        ##!#                    datacrs)
        ##!#    cbar4 = plt.colorbar(mesh4,ax=axcl,orientation='vertical',\
        ##!#        pad=0.03,label='TOA LWF [W/m2]')
        ##!#axcl.set_title('CERES LWF')

        ##!#data.close()

    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(dt_date_str.strftime('%Y%m%d%H%M')) 

        plt.rcParams.update({'hatch.color': 'r'})
        #ax3.pcolor(GOES_data_ch1['lon'],GOES_data_ch1['lat'],\
        #    hash_data1, hatch = '\\\\', alpha=0., transform = datacrs,\
        #    cmap = 'plasma')

        hatch_shape = '\\\\'
        hash0 = ax0.pcolor(GOES_data1['lon'],GOES_data1['lat'],\
            hash_data1[:GOES_data1['lat'].shape[0], :GOES_data1['lat'].shape[1]], hatch = hatch_shape, alpha=0., transform = datacrs)
        hash1 = ax1.pcolor(GOES_data2['lon'],GOES_data2['lat'],\
            hash_data1[:GOES_data2['lat'].shape[0], :GOES_data2['lat'].shape[1]], hatch = hatch_shape, alpha=0., transform = datacrs)
        hash2 = ax2.pcolor(GOES_data3['lon'],GOES_data3['lat'],\
            hash_data1[:GOES_data3['lat'].shape[0], :GOES_data3['lat'].shape[1]], hatch = hatch_shape, alpha=0., transform = datacrs)
        if(compare_OMI): 
            hash3 = axo.pcolor(GOES_data3['lon'],GOES_data3['lat'],\
                hash_data1[:GOES_data3['lat'].shape[0], :GOES_data3['lat'].shape[1]], hatch = '///', alpha=0., transform = datacrs)
        if(compare_CERES): 
            hash4 = axcs.pcolor(GOES_data3['lon'],GOES_data3['lat'],\
                hash_data1[:GOES_data3['lat'].shape[0], :GOES_data3['lat'].shape[1]], hatch = '///', alpha=0., transform = datacrs)
            hash5 = axcl.pcolor(GOES_data3['lon'],GOES_data3['lat'],\
                hash_data1[:GOES_data3['lat'].shape[0], :GOES_data3['lat'].shape[1]], hatch = '///', alpha=0., transform = datacrs)
    

    if(plot_ASOS_loc):
        print("Plotting ASOS site")
        plot_ASOS_locs(ax0,date_str,color='lime')
        plot_ASOS_locs(ax1,date_str,color='lime')
        plot_ASOS_locs(ax2,date_str,color='lime')
        if(compare_OMI): plot_ASOS_locs(axo,date_str,color='black')
        if(compare_CERES): 
            plot_ASOS_locs(axcs,date_str,color='black')
            plot_ASOS_locs(axcl,date_str,color='black')

    # Add subplot labels
    # ------------------
    plot_subplot_label(ax0, '(a)', backgroundcolor = 'white')
    plot_subplot_label(ax1, '(b)', backgroundcolor = 'white')
    plot_subplot_label(ax2, '(c)', backgroundcolor = 'white')


    cross_date = GOES_data1['cross_date']
    file_time  = GOES_data1['file_time']

    fig.tight_layout()

    #plt.suptitle(cross_date + ' ' + file_time)
    if(save):
        pdate = cross_date[:4] + cross_date[5:7] + cross_date[8:10] + file_time
        outname = 'goes_compare_ch' + str(channel1) + '_vs_ch' + \
            str(channel2) + '_ch' + str(channel3) + '_' + pdate + '_3panel.png'
        if(compare_OMI):
            outname = 'goes_compare_ch' + str(channel1) + '_vs_ch' + \
                str(channel2) + '_ch' + str(channel3) + '_omi_' + pdate + '_3panel.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else: 
        plt.show()

    if(return_GOES):
        return GOES_data1
    

def compare_GOES_channels(date_str,channel1,channel2,zoom=True,save=False,\
        plot_ASOS_loc = False,show_smoke = True):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    filename = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['goes']

    if(channel1 == 'red'):
        channel1 = 1
    elif(channel1 == 'green'):
        channel1 = 4
    elif(channel1 == 'blue'):
        channel1 = 3

    # -----------------------------------------------------------------------
    # Step 1: Call read_GOES_channel to read the desired GOES data from the
    # file and put it in a dictionary
    # -----------------------------------------------------------------------
    GOES_data1 = read_GOES_channel(dt_date_str.strftime('%Y%m%d%H%M'), channel1, zoom = zoom)
    GOES_data2 = read_GOES_channel(dt_date_str.strftime('%Y%m%d%H%M'), channel2, zoom = zoom)

    if(debug):
        print("Data1 max = ",np.max(GOES_data1['data']), "  Data1 min = ",\
        np.min(GOES_data1['data']))
        print("Data2 max = ",np.max(GOES_data2['data']), "  Data2 min = ",np.min(GOES_data2['data']))
    
        print(GOES_data1['data'].shape, GOES_data2['data'].shape)

    # --------------------------------------
    # Step 2: Set up figure to have 3 panels
    # --------------------------------------
    datacrs = ccrs.PlateCarree() 
    #mapcrs = ccrs.LambertConformal()
    mapcrs = init_proj(date_str)
    #mapcrs = ccrs.LambertConformal(central_longitude = \
    #    np.mean(plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['Lon']),\
    #    central_latitude = \
    #    np.mean(plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['Lat']))

    plt.close('all')
    fig = plt.figure(figsize=(14.5,5))
    ax0 = fig.add_subplot(1,3,2,projection = mapcrs)
    ax1 = fig.add_subplot(1,3,3,projection = mapcrs)
    ax2 = fig.add_subplot(1,3,1)

    # Step 3: Plot the GOES channel data in the first 2 panels
    # ---------------------------------------------------------
    mesh0 = ax0.pcolormesh(GOES_data1['lon'],GOES_data1['lat'],\
        GOES_data1['data'],cmap = GOES_data1['colors'], shading='auto', \
        transform = datacrs)

    cbar0 = plt.colorbar(mesh0,ax=ax0,orientation='vertical',\
        pad=0.03,label=GOES_data1['variable'])
    #cbar0 = plt.colorbar(mesh0,cax = ax0, orientation='horizontal',pad=0,\
    #    aspect=50,shrink = 0.850)
    #cbar0.ax.tick_params(labelsize=14)
    #cbar0.set_label(GOES_data1['variable'],fontsize=16,weight='bold')
   
    ax0.add_feature(cfeature.BORDERS)
    ax0.add_feature(cfeature.STATES)
    ax0.coastlines()
    if(zoom):
        ax0.set_extent([plot_limits_dict[GOES_data1['cross_date']][GOES_data1['file_time']]['Lon'][0], \
                        plot_limits_dict[GOES_data1['cross_date']][GOES_data1['file_time']]['Lon'][1], \
                        plot_limits_dict[GOES_data1['cross_date']][GOES_data1['file_time']]['Lat'][0], \
                        plot_limits_dict[GOES_data1['cross_date']][GOES_data1['file_time']]['Lat'][1]],\
                        datacrs)
    #ax0.set_title('GOES Ch. ' + str(channel1) + '\n' + \
    #    str(channel_dict[str(channel1)]['wavelength'][0]) + ' μm - ' + \
    #    str(channel_dict[str(channel1)]['wavelength'][1]) + ' μm')
    ax0.set_title('Channel ' + str(channel1) + '\n' + \
        channel_dict[str(channel1)]['wavelength_label']) 

    # Plot channel 2
    mesh1 = ax1.pcolormesh(GOES_data2['lon'],GOES_data2['lat'],\
        GOES_data2['data'],cmap = GOES_data2['colors'], shading='auto', \
        transform = datacrs) 

    cbar1 = plt.colorbar(mesh1,ax=ax1,orientation='vertical',\
        pad=0.03,label=GOES_data2['variable'])
    #cbar1 = plt.colorbar(mesh1,cax = ax1,orientation='horizontal',pad=1,\
    #    aspect=51,shrink = 1.851)
    #cbar1.ax.tick_params(labelsize=14)
    #cbar1.set_label(GOES_data2['variable'],fontsize=16,weight='bold')
    
    ax1.add_feature(cfeature.BORDERS)
    ax1.add_feature(cfeature.STATES)
    ax1.coastlines()
    if(zoom):
        ax1.set_extent([plot_limits_dict[GOES_data2['cross_date']][GOES_data2['file_time']]['Lon'][0], \
                        plot_limits_dict[GOES_data2['cross_date']][GOES_data2['file_time']]['Lon'][1], \
                        plot_limits_dict[GOES_data2['cross_date']][GOES_data2['file_time']]['Lat'][0], \
                        plot_limits_dict[GOES_data2['cross_date']][GOES_data2['file_time']]['Lat'][1]],\
                        datacrs)
    #ax1.set_title('GOES Ch. ' + str(channel2) + '\n' + \
    #    str(channel_dict[str(channel2)]['wavelength'][0]) + ' μm - ' + \
    #    str(channel_dict[str(channel2)]['wavelength'][1]) + ' μm')
    ax1.set_title('Channel ' + str(channel2) + '\n' + \
        channel_dict[str(channel2)]['wavelength_label']) 

    ##!## Shade areas that are inside the plume, defined currently as areas
    ##!## with reflectance below the mean
    ##!## -----------------------------------------------------------------
    ##!#if(GOES_data1['variable'] == 'Reflectance'):
    ##!#    hash_data1 = np.ma.masked_where(GOES_data1['data'] < \
    ##!#        np.nanmean(GOES_data1['data']),GOES_data1['data'])
    ##!#    hash_data2 = np.ma.masked_where(GOES_data1['data'] < \
    ##!#        np.nanmean(GOES_data1['data']),GOES_data2['data'])

    ##!#    nohash_data1 = np.ma.masked_where(GOES_data1['data'] >= \
    ##!#        np.nanmean(GOES_data1['data']),GOES_data1['data'])
    ##!#    nohash_data2 = np.ma.masked_where(GOES_data1['data'] >= \
    ##!#        np.nanmean(GOES_data1['data']),GOES_data2['data'])

    ##!#    hash1 = ax0.pcolor(GOES_data1['lon'],GOES_data1['lat'],\
    ##!#        hash_data1, hatch = '///', alpha=0., transform = datacrs),
    ##!#    hash2 = ax1.pcolor(GOES_data2['lon'],GOES_data2['lat'],\
    ##!#        hash_data2, hatch = '///', alpha=0., transform = datacrs),
    ##!#if(GOES_data2['variable'] == 'Reflectance'):
    ##!#    hash_data2 = np.ma.masked_where(GOES_data2['data'] < \
    ##!#        np.nanmean(GOES_data2['data']),GOES_data2['data'])
    ##!#    hash_data1 = np.ma.masked_where(GOES_data2['data'] < \
    ##!#        np.nanmean(GOES_data2['data']),GOES_data1['data'])

    ##!#    nohash_data2 = np.ma.masked_where(GOES_data2['data'] >= \
    ##!#        np.nanmean(GOES_data2['data']),GOES_data2['data'])
    ##!#    nohash_data1 = np.ma.masked_where(GOES_data2['data'] >= \
    ##!#        np.nanmean(GOES_data2['data']),GOES_data1['data'])
    ##!#

    ##!#    hash2 = ax1.pcolor(GOES_data2['lon'],GOES_data2['lat'],\
    ##!#        hash_data2, hatch = '///', alpha=0., transform = datacrs),
    ##!#    hash1 = ax0.pcolor(GOES_data1['lon'],GOES_data1['lat'],\
    ##!#        hash_data1, hatch = '///', alpha=0., transform = datacrs),
    
    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(dt_date_str.strftime('%Y%m%d%H%M')) 
        hash2 = ax1.pcolor(GOES_data2['lon'],GOES_data2['lat'],\
            hash_data1, hatch = '///', alpha=0., transform = datacrs),
        hash1 = ax0.pcolor(GOES_data1['lon'],GOES_data1['lat'],\
            hash_data1, hatch = '///', alpha=0., transform = datacrs),
    
    # Step 4: Plot scatter GOES channel comparison in third panel
    # ------------------------------------------------------------

    # Figure out which pixels are in the plume
    # ----------------------------------------
    max_ch = 350.

    tmp_data1 = np.copy(GOES_data1['data'])
    tmp_data2 = np.copy(GOES_data2['data'])

    tmp_data1 = np.ma.masked_where( (GOES_data1['data'] > max_ch) | \
        (GOES_data2['data'] > max_ch), tmp_data1)
    tmp_data2 = np.ma.masked_where( (GOES_data1['data'] > max_ch) | \
        (GOES_data2['data'] > max_ch), tmp_data2)

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
    ax2.scatter(GOES_data1['data'][np.where(~hash_data1.mask)], \
        GOES_data2['data'][np.where(~hash_data1.mask)], s=6, \
        color='tab:blue', label='Inside Plume')
    ax2.scatter(GOES_data1['data'][np.where(hash_data1.mask)], \
        GOES_data2['data'][np.where(hash_data1.mask)], s=6, \
        color='tab:orange', label='Outside Plume')
    #ax2.scatter(nohash_data1.compressed(), nohash_data2.compressed(), s=6, \
    #    color='tab:orange', label='Outside Plume')
    ax2.set_xlabel('Ch. ' + str(GOES_data1['channel']) +' ' + GOES_data1['variable'])
    ax2.set_ylabel('Ch. ' + str(GOES_data2['channel']) +' ' + GOES_data2['variable'])
    ax2.set_title('Pearson correlation: '+str(np.round(rval_p, 3)))
    ax2.legend()

    fig.tight_layout()

    if(plot_ASOS_loc):
        print("Plotting ASOS site")
        plot_ASOS_locs(ax0,date_str)
        plot_ASOS_locs(ax1,date_str)

    if(save):
        cross_date = GOES_data1['cross_date']
        file_time  = GOES_data1['file_time']
        pdate = cross_date[:4] + cross_date[5:7] + cross_date[8:10] + file_time
        outname = 'goes_compare_ch' + str(channel1) + '_ch' + str(channel2) + '_' + pdate + '_2pannelscatter.png'
        plt.savefig(outname,dpi=300)
        print('Saved image',outname)
    else:
        plt.show()

    #return hash_data1, nohash_data1, GOES_data1, GOES_data2

def compare_GOES_3scatter(date_str,channel0,channel1,channel2,channel3,\
        save=False, compare_OMI = False, compare_CERES = False, \
        avg_pixel = False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    filename = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['goes']

    if(channel1 == 'red'):
        channel1 = 1
    elif(channel1 == 'green'):
        channel1 = 4
    elif(channel1 == 'blue'):
        channel1 = 3

    # Step 1: Call read_GOES_channel to read the desired GOES data from the
    # file and put it in a dictionary
    # -----------------------------------------------------------------------
    GOES_data0 = read_GOES_channel(dt_date_str.strftime('%Y%m%d%H%M'), channel0, zoom = True)
    GOES_data1 = read_GOES_channel(dt_date_str.strftime('%Y%m%d%H%M'), channel1, zoom = True)
    GOES_data2 = read_GOES_channel(dt_date_str.strftime('%Y%m%d%H%M'), channel2, zoom = True)
    GOES_data3 = read_GOES_channel(dt_date_str.strftime('%Y%m%d%H%M'), channel3, zoom = True)

    if(debug):
        print("Data0 max = ",np.max(GOES_data0['data']), "  Data0 min = ",\
        np.min(GOES_data0['data']))
        print("Data1 max = ",np.max(GOES_data1['data']), "  Data1 min = ",\
            np.min(GOES_data1['data']))
        print("Data2 max = ",np.max(GOES_data2['data']), "  Data2 min = ",\
            np.min(GOES_data2['data']))
        print("Data3 max = ",np.max(GOES_data3['data']), "  Data3 min = ",\
            np.min(GOES_data3['data']))

    max_ch = 350.

    # Determine where the smoke is located
    # ------------------------------------
    hash_data1, nohash_data1 = find_plume(dt_date_str.strftime('%Y%m%d%H%M')) 

    tmp_data0 = np.copy(GOES_data0['data'])
    tmp_data1 = np.copy(GOES_data1['data'])
    tmp_data2 = np.copy(GOES_data2['data'])
    tmp_data3 = np.copy(GOES_data3['data'])
    tmp_lat0  = np.copy(GOES_data0['lat'])
    tmp_lon0  = np.copy(GOES_data0['lon'])

    if(not (tmp_data0.shape == tmp_data1.shape == tmp_data2.shape == tmp_data3.shape)):
        if(debug):
            print("shape mismatch")
        shapes = []
        shapes.append(tmp_data0.shape)
        shapes.append(tmp_data1.shape)
        shapes.append(tmp_data2.shape)
        shapes.append(tmp_data3.shape)

        min_shape = min(shapes)
        if(debug):
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

    cross_date = GOES_data0['cross_date']
    file_time  = GOES_data0['file_time']

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

    plot_scatter(ax0, tmp_data0, tmp_data1, GOES_data0, GOES_data1, \
        hash_data1)
    plot_scatter(ax1, tmp_data0, tmp_data2, GOES_data0, GOES_data2, \
        hash_data1)
    plot_scatter(ax2, tmp_data0, tmp_data3, GOES_data0, GOES_data3, \
        hash_data1)


    if(compare_OMI):
        print("Reading OMI data")

        #colocate_OMI(date_str, axo, avg_pixel=True)
        plot_scatter_OMI(date_str, GOES_data0, axo, avg_pixel = avg_pixel)
        ##!#colocate_OMI(date_str, tmp_data0, tmp_lat0, tmp_lon0, hash_data1,\
        ##!#             axo, avg_pixel = avg_pixel)

        ##!#data = h5py.File(plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['omi'],'r')
        ##!#LAT   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Latitude'][:,:]
        ##!#LON   = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/Longitude'][:,:]
        ##!#UVAI  = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Data Fields/UVAerosolIndex'][:,:]
        ##!#XTRACK = data['HDFEOS/SWATHS/Aerosol NearUV Swath/Geolocation Fields/XTrackQualityFlags'][:,:]
        ##!#mask_LAT = np.ma.masked_where( (XTRACK < -2e5) | (UVAI < 2.), LAT)
        ##!#mask_LON = np.ma.masked_where( (XTRACK < -2e5) | (UVAI < 2.), LON)
        ##!#mask_UVAI = np.ma.masked_where((XTRACK < -2e5) | (UVAI < 2.), UVAI)
        ##!#mask_LAT  = np.ma.masked_where((((LAT < plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][0]) | \
        ##!#                     (LAT > plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][1])) | \
        ##!#                    ((LON < plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][0]) | \
        ##!#                     (LON > plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][1]))), mask_LAT)
        ##!#mask_LON  = np.ma.masked_where((((LAT < plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][0]) | \
        ##!#                     (LAT > plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][1])) | \
        ##!#                    ((LON < plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][0]) | \
        ##!#                     (LON > plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][1]))), mask_LON)
        ##!#mask_UVAI = np.ma.masked_where((((LAT < plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][0]) | \
        ##!#                     (LAT > plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'][1])) | \
        ##!#                    ((LON < plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][0]) | \
        ##!#                     (LON > plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lon'][1]))), mask_UVAI)


        #hash_data1, nohash_data1 = find_plume(filename) 
        # Colocate the masked OMI data with the GOES channel0 data
        # ---------------------------------------------------------
        print("Colocating OMI data")


            # Plot hash_avg_goes against mask_UVAI 

            # Plot nohash_avg_goes against mask_UVAI 
 
    # End compare_OMI

    if(compare_CERES):
        print("Reading CERES data")

        plot_scatter_CERES(date_str, GOES_data0, axcs, avg_pixel = avg_pixel)

    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    if(compare_OMI):
        plt.legend(lines, labels, loc = 'lower center', bbox_to_anchor = (0, 0.01, 1, 1),\
            bbox_transform = plt.gcf().transFigure, ncol=2)
    else:
        plt.legend(lines, labels, loc = 'right', bbox_to_anchor = (0, 0., 1.0, 1),\
            bbox_transform = plt.gcf().transFigure, ncol=1)

    fig.tight_layout()

    if(save):
        pdate = cross_date[:4] + cross_date[5:7] + cross_date[8:10] + file_time
        outname = 'goes_compare_ch' + str(channel0) + '_vs_ch' + \
            str(channel1) + '_ch' + str(channel2) + '_ch' + \
            str(channel3) + '_' + pdate + '_3scatter.png'
        if(compare_OMI):
            outname = 'goes_compare_ch' + str(channel0) + '_vs_ch' + \
                str(channel1) + '_ch' + str(channel2) + '_ch' + \
                str(channel3) + '_omi_' + pdate + '_3scatter.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else: 
        plt.show()

def plot_GOES_spatial(GOES_data, pax, zoom, vmin = None, vmax = None, \
        ptitle = None):

    if(vmin is None):
        if('data_lim' in plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']].keys()):
            vmin = plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['data_lim'][GOES_data['channel']][0]
        else:
            vmin = np.nanmin(GOES_data['data'])
    if(vmax is None):
        if('data_lim' in plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']].keys()):
            vmax = plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['data_lim'][GOES_data['channel']][1]
        else:
            if(GOES_data['channel'] == 'wv_ir'):
                vmax = 1.5
            else:
                vmax = np.nanmax(GOES_data['data'])

    #mesh = ax.pcolormesh(GOES_data['lon'],GOES_data['lat'],\
    #    GOES_data['data'],cmap = GOES_data['colors'], shading='auto', \
    #    transform = datacrs) 

    # Plot channel 1
    mesh1 = pax.pcolormesh(GOES_data['lon'],GOES_data['lat'],\
        GOES_data['data'],cmap = GOES_data['colors'], shading='auto', \
        vmin = vmin, \
        vmax = vmax, transform = datacrs) 

    cbar1 = plt.colorbar(mesh1,ax=pax,orientation='vertical',\
        pad=0.03)
        #shrink = 0.870, pad=0.03,label=GOES_data['variable'])
    cbar1.set_label(GOES_data['variable'], size = 13, weight = 'bold')
    cbar1.ax.tick_params(labelsize = 11)

    pax.add_feature(cfeature.BORDERS)
    pax.add_feature(cfeature.STATES)
    pax.coastlines()
    if(zoom):
        pax.set_extent([plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['Lon'][0], \
                        plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['Lon'][1], \
                        plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['Lat'][0], \
                        plot_limits_dict[GOES_data['cross_date']][GOES_data['file_time']]['Lat'][1]],\
                        datacrs)
    if(ptitle == None):
        pax.set_title('Channel ' + str(GOES_data['channel']) + '\n' + \
            channel_dict[str(GOES_data['channel'])]['wavelength_label']) 
    else:
        pax.set_title(ptitle)


# dtype must be 'SWF' or 'LWF'
def plot_CERES_spatial(date_str, mask_LAT, mask_LON, mask_data, dtype, pax, \
        vmin = None, vmax = None, markersize = 170, ptitle = None, zoom = False):

    #print(plot_limits_dict[GOES_data3['cross_date']][GOES_data3['file_time']]['Lat'])

    plabel = 'TOA Flux [Wm$^{-2}$]'
    #plabel = 'TOA ' + dtype + ' [W/m2]'

    dt_date_str = datetime.strptime(date_str, "%Y%m%d%H%M")

    #scat3 = axcs.scatter(mask_LON, mask_LAT,mask_swf, transform = datacrs)
    mesh3 = pax.scatter(mask_LON.compressed(), mask_LAT.compressed(),\
        #s =  170, marker = 's', c = mask_data.compressed(), cmap = 'jet', \
        #s =  170, marker = 's', c = mask_data.compressed(), cmap = 'plasma', \
        s =  170, marker = 's', c = mask_data.compressed(), cmap = 'plasma', \
        vmin = vmin, vmax = vmax, transform = datacrs)
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
        pad=0.03)
    cbar3.set_label(plabel, size = 14, weight = 'bold')
    cbar3.ax.tick_params(labelsize = 12)
    if(ptitle == None):
        pax.set_title('CERES TOA ' + dtype)
    else:
        pax.set_title(ptitle)

def plot_scatter(lax,data1, data2, GOES_data1, GOES_data2, hash_data1,\
        xlabel = None, ylabel = None, slope = 'thiel-sen', \
        plot_legend = False):

    hrval_p = spearmanr(data1[np.where(~hash_data1.mask)].compressed(), \
        data2[np.where(~hash_data1.mask)].compressed())[0]
    nrval_p = spearmanr(data1[np.where(hash_data1.mask)].compressed(), \
        data2[np.where(hash_data1.mask)].compressed())[0]
    #print("0-1 Pearson:  ",hrval_p)
    lax.scatter(data1[np.where(hash_data1.mask)].compressed(), \
        data2[np.where(hash_data1.mask)].compressed(), s=6, \
        color='tab:orange', label = 'Outside plume')
    lax.scatter(data1[np.where(~hash_data1.mask)].compressed(), \
        data2[np.where(~hash_data1.mask)].compressed(), s=6, \
        color='tab:blue', label = 'Inside plume')

    # Plot trend lines for each set
    # -----------------------------
    print("Calculating GOES " + str(GOES_data1['channel']) + ' vs ' + \
        str(GOES_data2['channel']) + ' smoke trend')
    plot_trend_line(lax, \
        data1[np.where(~hash_data1.mask)].compressed(), \
        data2[np.where(~hash_data1.mask)].compressed(), \
        color='tab:blue', slope = slope)

    print("Calculating GOES " + str(GOES_data1['channel']) + ' vs ' + \
        str(GOES_data2['channel']) + ' clear trend')
    plot_trend_line(lax, \
        data1[np.where(hash_data1.mask)].compressed(), \
        data2[np.where(hash_data1.mask)].compressed(), \
        color='tab:orange')

    if(GOES_data2['channel'] == 'wv_ir'):
        lax.set_ylim(lax.get_ylim()[0], 2.0)
        #vmax = 1.5
    #else:
    #    vmax = np.nanmax(GOES_data['data'])


    if(xlabel == None):
        xlabel = 'Ch. ' + str(GOES_data1['channel']) +' [' + \
        channel_dict[str(GOES_data1['channel'])]['wavelength_label'] + '] '+ \
        GOES_data1['variable']
    if(ylabel == None):
        ylabel ='Ch. ' + str(GOES_data2['channel']) +' [' + \
        channel_dict[str(GOES_data2['channel'])]['wavelength_label'] + ']' + \
        GOES_data2['variable'] 

    lax.set_xlabel(xlabel, fontsize = 14)
    lax.set_ylabel(ylabel, fontsize = 14)

    if(plot_legend):
        lax.legend()

    #plt.suptitle('Aqua GOES ' + GOES_data1['cross_date'] + ' ' + \
    #    GOES_data1['file_time'])
    #lax.set_title('Smoke correlation: '+str(np.round(hrval_p, 3))+'\n'+\
    #              'Clear correlation: '+str(np.round(nrval_p, 3)))
   
 
#def plot_scatter_OMI(date_str, tmp_data0, tmp_lat0, tmp_lon0, channel, \
def plot_scatter_OMI(date_str, GOES_data, pax, avg_pixel = False, \
        xlabel = None, ylabel = None, ptitle = None):

    # Determine where the smoke is located
    # ------------------------------------
    hash_data, nohash_data = find_plume(date_str) 

    tmp_data = np.copy(GOES_data['data'])
    tmp_lat  = np.copy(GOES_data['lat'])
    tmp_lon  = np.copy(GOES_data['lon'])

    if(tmp_data.shape != hash_data.shape):
        print("shape mismatch")
        shapes = []
        shapes.append(tmp_data.shape)
        shapes.append(hash_data.shape)

        min_shape = min(shapes)
        print(min_shape)

        tmp_data  = tmp_data[:min_shape[0],:min_shape[1]]
        tmp_lat   = tmp_lat[:min_shape[0],:min_shape[1]]
        tmp_lon   = tmp_lon[:min_shape[0],:min_shape[1]]
        hash_data = hash_data[:min_shape[0],:min_shape[1]]

    print(hash_data.shape, tmp_data.shape)

    max_ch = 350.

    tmp_data = np.ma.masked_where((abs(tmp_data) > max_ch), tmp_data)
    tmp_lat  = np.ma.masked_where((abs(tmp_data) > max_ch), tmp_lat)
    tmp_lon  = np.ma.masked_where((abs(tmp_data) > max_ch), tmp_lon)

    # Remove the nans from the GOES data, lats, and lons for both the
    # in-the-plume (hashed) and outside-the-plume (nohashed) data
    hash_plot_data   = tmp_data[np.where(~hash_data.mask)].compressed()
    nohash_plot_data = tmp_data[np.where(hash_data.mask)].compressed()
    hash_plot_lat    = tmp_lat[np.where(~hash_data.mask)].compressed()
    nohash_plot_lat  = tmp_lat[np.where(hash_data.mask)].compressed()
    hash_plot_lon    = tmp_lon[np.where(~hash_data.mask)].compressed()
    nohash_plot_lon  = tmp_lon[np.where(hash_data.mask)].compressed()

    if(avg_pixel):

        LAT, LON, mask_UVAI, crnr_LAT, crnr_LON = \
            read_OMI_match_GOES(date_str, corners = True, min_AI = 1.0)

        print('averaging GOES pixels')

        # Pull out the OMI lat/lon corner pixels that are associated with
        # the masked OMI data
        # ---------------------------------------------------------------
        work_UVAI  = mask_UVAI.compressed()

        work_crnrLAT = crnr_LAT[np.where(~np.isnan(mask_UVAI)==True)[0],\
            np.where(~np.isnan(mask_UVAI)==True)[1],:] 
        work_crnrLON = crnr_LON[np.where(~np.isnan(mask_UVAI)==True)[0],\
            np.where(~np.isnan(mask_UVAI)==True)[1],:] 

        # Declare full arrays to hold the averaged GOES data
        # Must be dimensioned to the shape of the OMI data
        # ---------------------------------------------------
        if(debug):
            print("work UVAI shape = ",work_UVAI.shape)
        hash_avg_goes   = np.full(work_UVAI.shape,np.nan)
        hash_avg_mlat    = np.full(work_UVAI.shape,np.nan)
        hash_avg_mlon    = np.full(work_UVAI.shape,np.nan)
        nohash_avg_goes = np.full(work_UVAI.shape,np.nan)
        nohash_avg_mlat  = np.full(work_UVAI.shape,np.nan)
        nohash_avg_mlon  = np.full(work_UVAI.shape,np.nan)

        # Loop over screened out OMI pixels
        # ---------------------------------
        if(debug):
            print("corners shape:",work_crnrLAT.shape)
        for ii in range(work_crnrLAT.shape[0]):

            # For each OMI pixel, find the in-plume and out-of-plume GOES
            # pixels that are within the the OMI corner values
            # -------------------------------------------------------------
            hash_in = np.where((hash_plot_lat >= np.nanmin(work_crnrLAT[ii,:])) & \
                               (hash_plot_lat <= np.nanmax(work_crnrLAT[ii,:])) & \
                               (hash_plot_lon >= np.nanmin(work_crnrLON[ii,:])) & \
                               (hash_plot_lon <= np.nanmax(work_crnrLON[ii,:])) )
            nohash_in = np.where((nohash_plot_lat >= np.nanmin(work_crnrLAT[ii,:])) & \
                                 (nohash_plot_lat <= np.nanmax(work_crnrLAT[ii,:])) & \
                                 (nohash_plot_lon >= np.nanmin(work_crnrLON[ii,:])) & \
                                 (nohash_plot_lon <= np.nanmax(work_crnrLON[ii,:])) )

            # Use the indices to pull out the GOES pixels
            hash_goes_data   = hash_plot_data[hash_in]
            hash_goes_lat    = hash_plot_lat[hash_in]
            hash_goes_lon    = hash_plot_lon[hash_in]
            nohash_goes_data = nohash_plot_data[nohash_in]
            nohash_goes_lat  = nohash_plot_lat[nohash_in]
            nohash_goes_lon  = nohash_plot_lon[nohash_in]

            # Use the corner values to make a Polygon and determine if
            # each close GOES pixel is within the OMI pixel
            # --------------------------------------------------------
            # Generate a Polygon using the OMI corner points
            points = []
            for x, y in zip(work_crnrLON[ii,:], work_crnrLAT[ii,:]):
                points.append((x,y))
            omi_poly = Polygon(points)
          
            inside_hash = [omi_poly.contains(Point(hlon, hlat)) for \
                hlon, hlat in zip(hash_goes_lon, hash_goes_lat)]
            inside_nohash = [omi_poly.contains(Point(nlon, nlat)) for \
                nlon, nlat in zip(nohash_goes_lon, nohash_goes_lat)]
  
            # If at least 1 smokey GOES pixel is within the OMI pixel bounds, 
            # classify the OMI pixel as smokey
            # ----------------------------------------------------------------
            if((len(inside_hash) > 0) & (True in inside_hash)):
                hash_avg_goes[ii] = np.average(hash_goes_data[inside_hash])
                hash_avg_mlat[ii]  = np.average(hash_goes_lat[inside_hash])
                hash_avg_mlon[ii]  = np.average(hash_goes_lon[inside_hash])

            if((len(inside_nohash) > 0) & (True in inside_nohash)):
                nohash_avg_goes[ii] = np.average(nohash_goes_data[inside_nohash])
                nohash_avg_mlat[ii]  = np.average(nohash_goes_lat[inside_nohash])
                nohash_avg_mlon[ii]  = np.average(nohash_goes_lon[inside_nohash])

        # end screened corner data loop

        # Plot the screened OMI data vs the averaged GOES data
        # -----------------------------------------------------

        #print(work_UVAI[~np.ma.masked_invalid(hash_avg_goes).mask])
        
        hrval_p = spearmanr(np.ma.masked_invalid(hash_avg_goes).compressed(), \
            work_UVAI[~np.ma.masked_invalid(hash_avg_goes).mask])[0]
        pax.scatter(np.ma.masked_invalid(hash_avg_goes).compressed(), \
            work_UVAI[~np.ma.masked_invalid(hash_avg_goes).mask], s=6, \
            color='tab:blue')
        #axo.scatter(np.ma.masked_invalid(nohash_avg_goes).compressed(), \
        #    work_UVAI[~np.ma.masked_invalid(nohash_avg_goes).mask], s=6, \
        #    color='tab:orange')
        #axo.scatter(nohash_plot_data0, nohash_match_OMI, s=6, \
        #    color='tab:orange', label='Outside Plume')

        # Plot trend lines for each set
        # -----------------------------
        print("Calculating GOES/OMI trend")
        plot_trend_line(pax, \
            np.ma.masked_invalid(hash_avg_goes).compressed(), \
            work_UVAI[~np.ma.masked_invalid(hash_avg_goes).mask], \
            color='tab:blue')

        if(xlabel == None):
            xlabel = 'Averaged Ch. ' + str(GOES_data['channel']) +' [' + \
            channel_dict[str(GOES_data['channel'])]['wavelength_label'] + \
            GOES_data['variable']
        pax.set_xlabel(xlabel, fontsize = 14)
        if(ptitle is None):
            ptitle = 'Smoke correlation: '+str(np.round(hrval_p, 3))
        pax.set_title(ptitle)
        pax.set_ylabel('OMI UVAI', fontsize = 14)

    else:
        LAT, LON, mask_UVAI, crnr_LAT, crnr_LON = \
            read_OMI_match_GOES(date_str, corners = True)

        hash_match_OMI   = np.full(hash_plot_data.shape,-9.)
        hash_match_LAT   = np.full(hash_plot_lat.shape,-9.)
        hash_match_LON   = np.full(hash_plot_lon.shape,-9.)
        nohash_match_OMI = np.full(nohash_plot_data.shape,-9.)
        nohash_match_LAT = np.full(nohash_plot_lat.shape,-9.)
        nohash_match_LON = np.full(nohash_plot_lon.shape,-9.)

        #print(hash_plot_data.shape)
        for ii in range(hash_match_OMI.shape[0]):
            # Find the gridpoint in the gridded lat/lon data that 
            # corresponds to the station at slat and slon
            # ---------------------------------------------------- 
            o_idx = nearest_gridpoint(hash_plot_lat[ii], hash_plot_lon[ii],\
                LAT, LON)
                #mask_LAT, mask_LON)

            if(len(o_idx[0]) > 1):
                o_idx = (np.array([o_idx[0][0]])), (np.array([o_idx[1][0]]))
            hash_match_OMI[ii] = mask_UVAI[o_idx]
            hash_match_LAT[ii] = LAT[o_idx] 
            hash_match_LON[ii] = LON[o_idx] 
            #hash_match_LAT[ii] = mask_LAT[o_idx] 
            #hash_match_LON[ii] = mask_LON[o_idx] 

        #print(nohash_plot_data.shape)
        for ii in range(nohash_match_OMI.shape[0]):
            # Find the gridpoint in the gridded lat/lon data that 
            # corresponds to the station at slat and slon
            # ---------------------------------------------------- 
            o_idx = nearest_gridpoint(nohash_plot_lat[ii], nohash_plot_lon[ii],\
                LAT, LON)
                #mask_LAT, mask_LON)

            if(len(o_idx[0]) > 1):
                o_idx = (np.array([o_idx[0][0]])), (np.array([o_idx[1][0]]))
            nohash_match_OMI[ii] = mask_UVAI[o_idx]
            nohash_match_LAT[ii] = LAT[o_idx] 
            nohash_match_LON[ii] = LON[o_idx] 
            #nohash_match_LAT[ii] = mask_LAT[o_idx] 
            #nohash_match_LON[ii] = mask_LON[o_idx] 

        #xy = np.vstack([plot_data0,match_OMI])
        #z = stats.gaussian_kde(xy)(xy)
        #axo.scatter(plot_data0,match_OMI,c=z,s=6)

        hrval_p = spearmanr(hash_plot_data, \
            hash_match_OMI)[0]
        pax.scatter(hash_plot_data, hash_match_OMI, s=6, \
            color='tab:blue')

        # Plot trend lines for each set
        # -----------------------------
        print("Calculating GOES/OMI trend")
        plot_trend_line(pax, \
            hash_plot_data, \
            hash_match_OMI, \
            color='tab:blue')
        #axo.scatter(nohash_plot_data0, nohash_match_OMI, s=6, \
        #    color='tab:orange')

        if(xlabel == None):
            xlabel = 'Ch. ' + str(GOES_data['channel']) +' [' + \
            channel_dict[str(GOES_data['channel'])]['wavelength_label'] + \
            GOES_data['variable']
        pax.set_xlabel(xlabel, fontsize = 14)
        if(ptitle is None):
            ptitle = 'Smoke correlation: '+str(np.round(hrval_p, 3))
        pax.set_title(ptitle)
        pax.set_ylabel('OMI UVAI', fontsize = 14)
        

def plot_scatter_CERES(date_str, GOES_data, pax, avg_pixel = False,\
        plume_only = True, plot_total_flux = True):

    # Determine where the smoke is located
    # ------------------------------------
    hash_data, nohash_data = find_plume(date_str) 

    tmp_data = np.copy(GOES_data['data'])
    tmp_lat  = np.copy(GOES_data['lat'])
    tmp_lon  = np.copy(GOES_data['lon'])

    if(tmp_data.shape != hash_data.shape):
        print("shape mismatch")
        shapes = []
        shapes.append(tmp_data.shape)
        shapes.append(hash_data.shape)

        min_shape = min(shapes)
        #print(min_shape)

        tmp_data  = tmp_data[:min_shape[0],:min_shape[1]]
        tmp_lat   = tmp_lat[:min_shape[0],:min_shape[1]]
        tmp_lon   = tmp_lon[:min_shape[0],:min_shape[1]]
        hash_data = hash_data[:min_shape[0],:min_shape[1]]

    #print(hash_data.shape, tmp_data.shape)

    max_ch = 350.

    tmp_data = np.ma.masked_where((abs(tmp_data) > max_ch), tmp_data)
    tmp_lat  = np.ma.masked_where((abs(tmp_data) > max_ch), tmp_lat)
    tmp_lon  = np.ma.masked_where((abs(tmp_data) > max_ch), tmp_lon)

    hash_check_lat    = tmp_lat[np.where(~hash_data.mask)].compressed()
    hash_check_lon    = tmp_lon[np.where(~hash_data.mask)].compressed()
    if(plume_only):
        # Remove the nans from the GOES data, lats, and lons for both the
        # in-the-plume (hashed) and outside-the-plume (nohashed) data
        hash_plot_data   = tmp_data[np.where(~hash_data.mask)].compressed()
        nohash_plot_data = tmp_data[np.where(hash_data.mask)].compressed()
        hash_plot_lat    = tmp_lat[np.where(~hash_data.mask)].compressed()
        nohash_plot_lat  = tmp_lat[np.where(hash_data.mask)].compressed()
        hash_plot_lon    = tmp_lon[np.where(~hash_data.mask)].compressed()
        nohash_plot_lon  = tmp_lon[np.where(hash_data.mask)].compressed()
    else:
        hash_plot_data = tmp_data    
        hash_plot_lat  = tmp_lat    
        hash_plot_lon  = tmp_lon    
    mask_LAT, mask_LON, mask_swf, mask_lwf = read_CERES_match_GOES(date_str)

    crnr_LAT = getCorners_1d(mask_LAT.data) ; crnr_LON = getCorners_1d(mask_LON)


    #if(avg_pixel):
    
    # Average a certain number of GOES pixels around each CERES
    # pixel.
    # ----------------------------------------------------------

    # total arrays are dimensioned to the size of the CERES
    # arrays. Contain all matched GOES data (in and out)
    total_match_SWF    = np.full(mask_swf.shape, np.nan)
    total_match_LWF    = np.full(mask_lwf.shape, np.nan)
    total_match_LAT    = np.full(mask_LAT.shape, np.nan)
    total_match_LON    = np.full(mask_LON.shape, np.nan)
    # hash arrays are dimensioned to the size of the CERES
    # arrays. Contain only in-plume matched GOES data
    hash_match_SWF    = np.full(mask_swf.shape, np.nan)
    hash_match_LWF    = np.full(mask_lwf.shape, np.nan)
    hash_match_LAT    = np.full(mask_LAT.shape, np.nan)
    hash_match_LON    = np.full(mask_LON.shape, np.nan)
    # nohash arrays are dimensioned to the size of the CERES
    # arrays. Contain only out=of-plume matched GOES data
    nohash_match_SWF    = np.full(mask_swf.shape, np.nan)
    nohash_match_LWF    = np.full(mask_lwf.shape, np.nan)
    nohash_match_LAT    = np.full(mask_LAT.shape, np.nan)
    nohash_match_LON    = np.full(mask_LON.shape, np.nan)

    #print(hash_plot_lat.shape)

    # Loop over the lower-resolution CERES data.    
    # Either grab the nearest GOES pixel to the CERES pixel
    # OR average the GOES pixels within a certain range of the
    # CERES pixel.
    # ---------------------------------------------------------
    for ii in range(mask_swf.shape[0]):

        lat_list = sorted([crnr_LAT[ii], crnr_LAT[ii+1]])
        lon_list = sorted([crnr_LON[ii], crnr_LON[ii+1]])

        # NEW APPROACH: use calculated pixel corners and use where 
        # statement to find all GOES pixels within the CERES bounds
        # ----------------------------------------------------------
        #in_idx = np.where( ((hash_plot_lat >= lat_list[0]) & \
        #    (hash_plot_lat <= lat_list[1])) & \
        #   ((hash_plot_lon >= lon_list[0]) & \
        #    (hash_plot_lon <= lon_list[1])))


        ##!##fun_c = np.maximum(np.abs(grid_lat - slat), \
        ##!##    np.abs(grid_lon - slon))
        ##!##m_idx = np.where(fun_c == np.min(fun_c))

        if(avg_pixel):
            ho_idx = np.where( (abs((mask_LAT[ii] - hash_plot_lat)) < 0.10) & \
                (abs((mask_LON[ii] - hash_plot_lon)) < 0.10) )
            
            #hash_plot_lat    = tmp_lat[np.where(~hash_data.mask)].compressed()
            #nohash_plot_lat  = tmp_lat[np.where(hash_data.mask)].compressed()
            #hash_plot_lon    = tmp_lon[np.where(~hash_data.mask)].compressed()
        else:
            # Find the gridpoint in the GOES lat/lon data that 
            # corresponds to the AIRS pixel
            # ---------------------------------------------------- 
            ho_idx = nearest_gridpoint(mask_LAT[ii], mask_LON[ii],\
                hash_plot_lat, hash_plot_lon)

            if(len(ho_idx[0]) > 1):
                ho_idx = (np.array([ho_idx[0][0]])), \
                    (np.array([ho_idx[1][0]]))
   
        if(len(hash_plot_data[ho_idx]) > 0):
            # if any smoke-classified GOES pixel is within the CERES pixel,
            # classify the CERES pixel as "in-plume"
            # -------------------------------------------------------------- 
            if(True in ~hash_data.mask[ho_idx]): 
                hash_match_SWF[ii] = np.nanmean(hash_plot_data[ho_idx])
                hash_match_LWF[ii] = np.nanmean(hash_plot_data[ho_idx])
                hash_match_LAT[ii] = mask_LAT[ii]
                hash_match_LON[ii] = mask_LON[ii]
            else: 
                nohash_match_SWF[ii] = np.nanmean(hash_plot_data[ho_idx])
                nohash_match_LWF[ii] = np.nanmean(hash_plot_data[ho_idx])
                nohash_match_LAT[ii] = mask_LAT[ii]
                nohash_match_LON[ii] = mask_LON[ii]
            #print(mask_LAT[ii], ~hash_data.mask[ho_idx])

            # Average the n GOES pixels around the current CERES
            # pixel. 
            # total_match arrays contain all colocated GOES data
            # both in-plume and out-of-plume.
            # --------------------------------------------------- 
            total_match_SWF[ii] = np.nanmean(hash_plot_data[ho_idx])
            total_match_LWF[ii] = np.nanmean(hash_plot_data[ho_idx])
            total_match_LAT[ii] = mask_LAT[ii]
            total_match_LON[ii] = mask_LON[ii]


    #plt.close('all')
    fig = plt.figure(figsize=(8,5))
    #dt_date_str = datetime.strptime(date_str, "%Y%m%d%H%M")
    mapcrs = init_proj(date_str)
    #mapcrs = ccrs.LambertConformal(central_longitude = \
    #    np.mean(plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['Lon']),\
    #    central_latitude = \
    #np.mean(plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['Lat']))

    tax1 = fig.add_subplot(2,3,1,projection = mapcrs) # smoke SW
    tax2 = fig.add_subplot(2,3,2,projection = mapcrs) # smoke LW
    tax3 = fig.add_subplot(2,3,3,projection = mapcrs) # smoke CH31
    tax4 = fig.add_subplot(2,3,4,projection = mapcrs) # total SW
    tax5 = fig.add_subplot(2,3,5,projection = mapcrs) # total LW
    tax6 = fig.add_subplot(2,3,6,projection = mapcrs) # total CH31

    #print('old: ',hash_match_LWF.shape)
    total_match_LAT = np.ma.masked_invalid(total_match_LAT)
    total_match_LON = np.ma.masked_invalid(total_match_LON)
    total_match_LWF = np.ma.masked_invalid(total_match_LWF)
    total_match_SWF = np.ma.masked_invalid(total_match_SWF)

    hash_match_LAT = np.ma.masked_invalid(hash_match_LAT)
    hash_match_LON = np.ma.masked_invalid(hash_match_LON)
    hash_match_LWF = np.ma.masked_invalid(hash_match_LWF)
    hash_match_SWF = np.ma.masked_invalid(hash_match_SWF)

    nohash_match_LAT = np.ma.masked_invalid(nohash_match_LAT)
    nohash_match_LON = np.ma.masked_invalid(nohash_match_LON)
    nohash_match_LWF = np.ma.masked_invalid(nohash_match_LWF)
    nohash_match_SWF = np.ma.masked_invalid(nohash_match_SWF)

    mask_swf = np.ma.masked_invalid(mask_swf)
    mask_lwf = np.ma.masked_invalid(mask_lwf)
    mask_LAT = np.ma.masked_invalid(mask_LAT)
    mask_LON = np.ma.masked_invalid(mask_LON)
  
    ##!#total_match_SWF = np.ma.masked_invalid(total_match_SWF[~total_match_SWF.mask])
    ##!#total_match_LWF = np.ma.masked_invalid(total_match_LWF[~total_match_LWF.mask])
    ##!#total_match_LAT = np.ma.masked_invalid(total_match_LAT[~total_match_LAT.mask])
    ##!#total_match_LON = np.ma.masked_invalid(total_match_LON[~total_match_LON.mask])
    total_mask_swf  = np.ma.masked_invalid(mask_swf[~total_match_SWF.mask])
    total_mask_lwf  = np.ma.masked_invalid(mask_lwf[~total_match_LWF.mask])

    ##!#hash_match_SWF = np.ma.masked_invalid(hash_match_SWF[~hash_match_SWF.mask])
    ##!#hash_match_LWF = np.ma.masked_invalid(hash_match_LWF[~hash_match_LWF.mask])
    ##!#hash_match_LAT = np.ma.masked_invalid(hash_match_LAT[~hash_match_LAT.mask])
    ##!#hash_match_LON = np.ma.masked_invalid(hash_match_LON[~hash_match_LON.mask])
    hash_mask_swf  = np.ma.masked_invalid(mask_swf[~hash_match_SWF.mask])
    hash_mask_lwf  = np.ma.masked_invalid(mask_lwf[~hash_match_LWF.mask])
 
    ##!#nohash_match_SWF = np.ma.masked_invalid(nohash_match_SWF[~nohash_match_SWF.mask])
    ##!#nohash_match_LWF = np.ma.masked_invalid(nohash_match_LWF[~nohash_match_LWF.mask])
    ##!#nohash_match_LAT = np.ma.masked_invalid(nohash_match_LAT[~nohash_match_LAT.mask])
    ##!#nohash_match_LON = np.ma.masked_invalid(nohash_match_LON[~nohash_match_LON.mask])
    nohash_mask_swf  = np.ma.masked_invalid(mask_swf[~nohash_match_SWF.mask])
    nohash_mask_lwf  = np.ma.masked_invalid(mask_lwf[~nohash_match_LWF.mask])

    # ----------------------------------------------------------------------- 
    #
    # Plot the CERES and gridded GOES data on a temporary figure
    # Currently not saved anywhere.
    #
    # ----------------------------------------------------------------------- 
    plot_CERES_spatial(date_str, total_match_LAT, total_match_LON, \
        total_mask_swf, 'SWF', tax1, zoom=True)
    plot_CERES_spatial(date_str, total_match_LAT, total_match_LON, \
        total_mask_lwf, 'LWF', tax2, zoom=True)
    plot_CERES_spatial(date_str, total_match_LAT, total_match_LON, \
        total_match_LWF, 'SWF', tax3, vmin = np.nanmin(GOES_data['data']),\
        vmax = np.nanmax(GOES_data['data']), zoom=True)
    plot_CERES_spatial(date_str, mask_LAT, mask_LON, \
        mask_swf, 'SWF', tax4, zoom=True)
    plot_CERES_spatial(date_str, mask_LAT, mask_LON, \
        mask_lwf, 'LWF', tax5, zoom=True)
    plot_GOES_spatial(GOES_data, tax6, zoom = True)
 
    mask_total = mask_swf + mask_lwf

    lhrval_p = spearmanr(hash_match_LWF.compressed(), \
        hash_mask_lwf)[0]
    lnrval_p = spearmanr(nohash_match_LWF.compressed(), \
        nohash_mask_lwf)[0]
    pax.scatter(hash_match_LWF.compressed(), hash_mask_lwf,\
        s = 18, color='tab:blue', marker='x',label = 'LWF smoke')
    pax.scatter(nohash_match_LWF.compressed(), nohash_mask_lwf,\
        s = 18, color='tab:orange', marker='x',label = 'LWF clear')

    shrval_p = spearmanr(hash_match_SWF.compressed(), \
        hash_mask_swf)[0]
    snrval_p = spearmanr(nohash_match_SWF.compressed(), \
        nohash_mask_swf)[0]
    pax.scatter(hash_match_SWF.compressed(), hash_mask_swf,\
        s = 18, color='tab:blue', marker='o',label = 'SWF smoke')
    pax.scatter(nohash_match_SWF.compressed(), nohash_mask_swf,\
        s = 18, color='tab:orange', marker='o',label = 'SWF clear')

    # Plot total flux data
    # --------------------
    hash_match_TF   = hash_match_SWF + hash_match_LWF
    nohash_match_TF = nohash_match_SWF + nohash_match_LWF
    hash_mask_TF   = hash_mask_swf + hash_mask_lwf
    nohash_mask_TF = nohash_mask_swf + nohash_mask_lwf

    thrval_p = spearmanr(hash_match_TF.compressed(), \
        hash_mask_swf)[0]
    tnrval_p = spearmanr(nohash_match_TF.compressed(), \
        nohash_mask_swf)[0]
    if(plot_total_flux):
        pax.scatter(hash_match_SWF.compressed(), hash_mask_TF,\
            s = 18, color='tab:blue', marker='*',label = 'Total smoke')
        pax.scatter(nohash_match_SWF.compressed(), nohash_mask_TF,\
            s = 18, color='tab:orange', marker='*',label = 'Total clear')

    # Plot trend lines for each set
    # -----------------------------
    print("Calculating GOES/SWF smoke trend")
    plot_trend_line(pax, hash_match_SWF.compressed(), hash_mask_swf, \
        color='tab:blue')
    print("Calculating GOES/SWF clear trend")
    plot_trend_line(pax, nohash_match_SWF.compressed(), nohash_mask_swf, \
        color='tab:orange')

    print("Calculating GOES/LWF smoke trend")
    plot_trend_line(pax, hash_match_LWF.compressed(), hash_mask_lwf, \
        color='tab:blue',  linestyle = '--')
    print("Calculating GOES/LWF clear trend")
    plot_trend_line(pax, nohash_match_LWF.compressed(), nohash_mask_lwf, \
        color='tab:orange', linestyle = '--')

    if(plot_total_flux):
        print("Calculating GOES/Total smoke trend")
        plot_trend_line(pax, hash_match_SWF.compressed(), hash_mask_TF, \
            color='tab:blue',  linestyle = '--')
        print("Calculating GOES/Total clear trend")
        plot_trend_line(pax, nohash_match_SWF.compressed(), nohash_mask_TF, \
            color='tab:orange', linestyle = '--')

    pax.set_title('SWF smoke: '+str(np.round(shrval_p,3)) + '  SWF clear: ' + \
         str(np.round(snrval_p,3)) + \
        '\nLWF smoke: '+str(np.round(lhrval_p,3)) + '  LWF clear: ' + \
        str(np.round(lnrval_p,3))) 
    pax.set_xlabel('Ch. ' + str(GOES_data['channel']) +' [' + \
        channel_dict[str(GOES_data['channel'])]['wavelength_label'] + \
        GOES_data['variable'])
    pax.set_ylabel('CERES TOA Flux [Wm$^{-2}$]', fontsize = 14)
    pax.legend()
    #pax.set_title('Smoke correlation: '+str(np.round(lrval_p, 3)))


    ##!#axcl.scatter(nohash_plot_data0, nohash_match_LWF,\
    ##!#    s = 6, color='tab:orange')
    ##!#axcs.scatter(nohash_plot_data0, nohash_match_SWF,\
    ##!#    s = 6, color='tab:orange')

    ##!#axcl.set_xlabel('Ch. ' + str(GOES_data0['channel']) +' [' + \
    ##!#    channel_dict[str(channel0)]['wavelength_label'] + \
    ##!#    GOES_data0['variable'])
    ##!#axcl.set_title('Smoke correlation: '+str(np.round(lrval_p, 3)))
    ##!#axcs.set_xlabel('Ch. ' + str(GOES_data0['channel']) +' [' + \
    ##!#    channel_dict[str(channel0)]['wavelength_label'] + \
    ##!#    GOES_data0['variable'])
    ##!#axcs.set_title('Smoke correlation: '+str(np.round(srval_p, 3)))
    ##!#axcl.set_ylabel('CERES LWF [W/m2]')
    ##!#axcs.set_ylabel('CERES SWF [W/m2]')

    ##!#else:
    ##!#    hash_match_SWF   = np.full(hash_plot_data.shape,np.nan)
    ##!#    hash_match_LWF   = np.full(hash_plot_data.shape,np.nan)
    ##!#    hash_match_LAT   = np.full(hash_plot_lat.shape,np.nan)
    ##!#    hash_match_LON   = np.full(hash_plot_lon.shape,np.nan)
    ##!#    nohash_match_SWF = np.full(nohash_plot_data.shape,np.nan)
    ##!#    nohash_match_LWF = np.full(nohash_plot_data.shape,np.nan)
    ##!#    nohash_match_LAT = np.full(nohash_plot_lat.shape,np.nan)
    ##!#    nohash_match_LON = np.full(nohash_plot_lon.shape,np.nan)

    ##!#    for ii in range(hash_match_SWF.shape[0]):
    ##!#        # Find the gridpoint in the gridded lat/lon data that 
    ##!#        # corresponds to the station at slat and slon
    ##!#        # ---------------------------------------------------- 
    ##!#        o_idx = nearest_gridpoint(hash_plot_lat[ii], hash_plot_lon[ii],\
    ##!#            mask_LAT, mask_LON)

    ##!#        if(len(o_idx[0]) > 1):
    ##!#            o_idx = (np.array([o_idx[0][0]])), (np.array([o_idx[1][0]]))
    ##!#        hash_match_SWF[ii] = mask_swf[o_idx]
    ##!#        hash_match_LWF[ii] = mask_lwf[o_idx]
    ##!#        hash_match_LAT[ii] = mask_LAT[o_idx] 
    ##!#        hash_match_LON[ii] = mask_LON[o_idx] 

    ##!#    print(nohash_plot_data.shape)
    ##!#    for ii in range(nohash_match_SWF.shape[0]):
    ##!#        # Find the gridpoint in the gridded lat/lon data that 
    ##!#        # corresponds to the station at slat and slon
    ##!#        # ---------------------------------------------------- 
    ##!#        o_idx = nearest_gridpoint(nohash_plot_lat[ii], nohash_plot_lon[ii],\
    ##!#            mask_LAT, mask_LON)

    ##!#        #print(o_idx, o_idx[0], o_idx[0].shape)
    ##!#        if(len(o_idx[0].shape) > 1):
    ##!#            o_idx = (np.array([o_idx[0][0]])), (np.array([o_idx[1][0]]))
    ##!#        nohash_match_SWF[ii] = mask_swf[o_idx]
    ##!#        nohash_match_LWF[ii] = mask_lwf[o_idx]
    ##!#        nohash_match_LAT[ii] = mask_LAT[o_idx] 
    ##!#        nohash_match_LON[ii] = mask_LON[o_idx] 

    ##!#    #xy = np.vstack([plot_data0,match_OMI])
    ##!#    #z = stats.gaussian_kde(xy)(xy)
    ##!#    #axo.scatter(plot_data0,match_OMI,c=z,s=6)

    ##!#    lrval_p = spearmanr(hash_plot_data, \
    ##!#        hash_match_LWF)[0]
    ##!#    ##plt.scatter(hash_plot_data, hash_match_SWF,\
    ##!#    ##    s = 6, color='tab:blue')
    ##!#    ##plt.scatter(hash_plot_data, hash_match_LWF,\
    ##!#    ##    s = 6, color='tab:orange')
    ##!#    #axcl.scatter(hash_plot_data, hash_match_LWF,\
    ##!#    #    s = 6, color='tab:blue')
    ##!#    mask_total = mask_swf + mask_lwf
    ##!#    srval_p = spearmanr(hash_plot_data, \
    ##!#        hash_match_SWF)[0]
    ##!#    pax.scatter(hash_plot_data, hash_match_SWF,\
    ##!#        s = 6, color='tab:blue', label = 'SWF')
    ##!#    pax.scatter(hash_plot_data, hash_match_LWF,\
    ##!#        s = 6, color='tab:orange', label = 'LWF')
    ##!#    ##!#axcl.scatter(nohash_plot_data0, nohash_match_LWF,\
    ##!#    ##!#    s = 6, color='tab:orange')
    ##!#    ##!#axcs.scatter(nohash_plot_data0, nohash_match_SWF,\
    ##!#    ##!#    s = 6, color='tab:orange')

    ##!###axcl.set_xlabel('Ch. ' + str(GOES_data['channel']) +' [' + \
    ##!###    channel_dict[str(GOES_data['channel'])]['wavelength_label'] + \
    ##!###    GOES_data0['variable'])
    ##!###axcl.set_title('Smoke correlation: '+str(np.round(lrval_p, 3)))
    ##!#pax.set_xlabel('Ch. ' + str(GOES_data['channel']) +' [' + \
    ##!#    channel_dict[str(GOES_data['channel'])]['wavelength_label'] + \
    ##!#    GOES_data['variable'])
    ##!#pax.set_title('Smoke correlation: '+str(np.round(srval_p, 3)))
    ##!##axcl.set_ylabel('CERES LWF [W/m2]')
    ##!#pax.set_ylabel('CERES Flux [W/m2]')
    ##!#pax.legend()
    ##!## end compare_CERES

def plot_scatter_OMI_CERES(date_str, GOES_data, pax, avg_pixel = False,\
        plume_only = True, xlabel = None, ylabel = None, ptitle = None):

    # Determine where the smoke is located
    # ------------------------------------
    hash_data, nohash_data = find_plume(date_str) 

    tmp_data = np.copy(GOES_data['data'])
    tmp_lat  = np.copy(GOES_data['lat'])
    tmp_lon  = np.copy(GOES_data['lon'])

    if(tmp_data.shape != hash_data.shape):
        print("shape mismatch")
        shapes = []
        shapes.append(tmp_data.shape)
        shapes.append(hash_data.shape)

        min_shape = min(shapes)
        #print(min_shape)

        tmp_data  = tmp_data[:min_shape[0],:min_shape[1]]
        tmp_lat   = tmp_lat[:min_shape[0],:min_shape[1]]
        tmp_lon   = tmp_lon[:min_shape[0],:min_shape[1]]
        hash_data = hash_data[:min_shape[0],:min_shape[1]]

    max_ch = 350.

    tmp_data = np.ma.masked_where((abs(tmp_data) > max_ch), tmp_data)
    tmp_lat  = np.ma.masked_where((abs(tmp_data) > max_ch), tmp_lat)
    tmp_lon  = np.ma.masked_where((abs(tmp_data) > max_ch), tmp_lon)

    hash_check_lat    = tmp_lat[np.where(~hash_data.mask)].compressed()
    hash_check_lon    = tmp_lon[np.where(~hash_data.mask)].compressed()
    if(plume_only):
        # Remove the nans from the GOES data, lats, and lons for both the
        # in-the-plume (hashed) and outside-the-plume (nohashed) data
        hash_plot_data   = tmp_data[np.where(~hash_data.mask)].compressed()
        nohash_plot_data = tmp_data[np.where(hash_data.mask)].compressed()
        hash_plot_lat    = tmp_lat[np.where(~hash_data.mask)].compressed()
        nohash_plot_lat  = tmp_lat[np.where(hash_data.mask)].compressed()
        hash_plot_lon    = tmp_lon[np.where(~hash_data.mask)].compressed()
        nohash_plot_lon  = tmp_lon[np.where(hash_data.mask)].compressed()
    else:
        hash_plot_data = tmp_data    
        hash_plot_lat  = tmp_lat    
        hash_plot_lon  = tmp_lon    

    # ------------------------------------------------------------------------
    # 
    #  Read in OMI and CERES data
    #
    # ------------------------------------------------------------------------

    # Read in the corresponding CERES data
    mask_LAT, mask_LON, mask_swf, mask_lwf = read_CERES_match_GOES(date_str)
    
    # Calculate the corners of the CERES pixels
    crnr_LAT = getCorners_1d(mask_LAT.data) ; crnr_LON = getCorners_1d(mask_LON)

    LAT, LON, mask_UVAI = read_OMI_match_GOES(date_str)
    #LAT, LON, mask_UVAI, crnr_LAT, crnr_LON = \
    #    read_OMI_match_GOES(date_str, corners = True)

    #if(avg_pixel):
   
    # Loop over the lower resolution data (CERES) and find the 
    # OMI pixel that is closest to each CERES pixel 
    # ----------------------------------------------------------

    # total arrays are dimensioned to the size of the CERES
    # arrays. Contain all matched GOES data (in and out)
    total_omi_match_SWF    = np.full(mask_swf.shape, np.nan)
    total_omi_match_LWF    = np.full(mask_lwf.shape, np.nan)
    total_omi_match_LAT    = np.full(mask_LAT.shape, np.nan)
    total_omi_match_LON    = np.full(mask_LON.shape, np.nan)
    ##!## hash arrays are dimensioned to the size of the CERES
    ##!## arrays. Contain only in-plume matched GOES data
    ##!#hash_match_SWF    = np.full(mask_swf.shape, np.nan)
    ##!#hash_match_LWF    = np.full(mask_lwf.shape, np.nan)
    ##!#hash_match_LAT    = np.full(mask_LAT.shape, np.nan)
    ##!#hash_match_LON    = np.full(mask_LON.shape, np.nan)
    ##!## nohash arrays are dimensioned to the size of the CERES
    ##!## arrays. Contain only out=of-plume matched GOES data
    ##!#nohash_match_SWF    = np.full(mask_swf.shape, np.nan)
    ##!#nohash_match_LWF    = np.full(mask_lwf.shape, np.nan)
    ##!#nohash_match_LAT    = np.full(mask_LAT.shape, np.nan)
    ##!#nohash_match_LON    = np.full(mask_LON.shape, np.nan)

    if(debug):
        print(hash_plot_lat.shape)

    # Loop over the lower-resolution CERES data.    
    # Either grab the nearest GOES pixel to the CERES pixel
    # OR average the GOES pixels within a certain range of the
    # CERES pixel.
    # ---------------------------------------------------------
    for ii in range(mask_swf.shape[0]):

        #lat_list = sorted([crnr_LAT[ii], crnr_LAT[ii+1]])
        #lon_list = sorted([crnr_LON[ii], crnr_LON[ii+1]])

        # NEW APPROACH: use calculated pixel corners and use where 
        # statement to find all GOES pixels within the CERES bounds
        # ----------------------------------------------------------
        #in_idx = np.where( ((hash_plot_lat >= lat_list[0]) & \
        #    (hash_plot_lat <= lat_list[1])) & \
        #   ((hash_plot_lon >= lon_list[0]) & \
        #    (hash_plot_lon <= lon_list[1])))


        ##!##fun_c = np.maximum(np.abs(grid_lat - slat), \
        ##!##    np.abs(grid_lon - slon))
        ##!##m_idx = np.where(fun_c == np.min(fun_c))

        #if(avg_pixel):
        #    ho_idx = np.where( (abs((mask_LAT[ii] - hash_plot_lat)) < 0.10) & \
        #        (abs((mask_LON[ii] - hash_plot_lon)) < 0.10) )
            
            #hash_plot_lat    = tmp_lat[np.where(~hash_data.mask)].compressed()
            #nohash_plot_lat  = tmp_lat[np.where(hash_data.mask)].compressed()
            #hash_plot_lon    = tmp_lon[np.where(~hash_data.mask)].compressed()
        #else:
        # Find the gridpoint in the GOES lat/lon data that 
        # corresponds to the AIRS pixel
        # ---------------------------------------------------- 
        ho_idx = nearest_gridpoint(mask_LAT[ii], mask_LON[ii],\
            LAT, LON)

        if(len(ho_idx[0]) > 1):
            ho_idx = (np.array([ho_idx[0][0]])), \
                (np.array([ho_idx[1][0]]))
   
        if(len(mask_UVAI[ho_idx]) > 0):
            ##!## if any smoke-classified GOES pixel is within the CERES pixel,
            ##!## classify the CERES pixel as "in-plume"
            ##!## -------------------------------------------------------------- 
            ##!#if(True in ~hash_data.mask[ho_idx]): 
            ##!#    hash_match_SWF[ii] = np.nanmean(hash_plot_data[ho_idx])
            ##!#    hash_match_LWF[ii] = np.nanmean(hash_plot_data[ho_idx])
            ##!#    hash_match_LAT[ii] = mask_LAT[ii]
            ##!#    hash_match_LON[ii] = mask_LON[ii]
            ##!#else: 
            ##!#    nohash_match_SWF[ii] = np.nanmean(hash_plot_data[ho_idx])
            ##!#    nohash_match_LWF[ii] = np.nanmean(hash_plot_data[ho_idx])
            ##!#    nohash_match_LAT[ii] = mask_LAT[ii]
            ##!#    nohash_match_LON[ii] = mask_LON[ii]
            #print(mask_LAT[ii], mask_LON[ii], LAT[ho_idx], LON[ho_idx])
            # ~hash_data.mask[ho_idx])

            # Average the n GOES pixels around the current CERES
            # pixel. 
            # total_match arrays contain all colocated GOES data
            # both in-plume and out-of-plume.
            # --------------------------------------------------- 
            #print(ho_idx, mask_UVAI[ho_idx])
            total_omi_match_SWF[ii] = mask_UVAI[ho_idx]
            total_omi_match_LWF[ii] = mask_UVAI[ho_idx]
            total_omi_match_LAT[ii] = mask_LAT[ii]
            total_omi_match_LON[ii] = mask_LON[ii]


    ##!##plt.close('all')
    ##!#fig = plt.figure(figsize=(8,5))
    ##!#dt_date_str = datetime.strptime(date_str, "%Y%m%d%H%M")
    ##!#mapcrs = ccrs.LambertConformal(central_longitude = \
    ##!#    np.mean(plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['Lon']),\
    ##!#    central_latitude = \
    ##!#np.mean(plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['Lat']))

    ##!#tax1 = fig.add_subplot(2,3,1,projection = mapcrs) # smoke SW
    ##!#tax2 = fig.add_subplot(2,3,2,projection = mapcrs) # smoke LW
    ##!#tax3 = fig.add_subplot(2,3,3,projection = mapcrs) # smoke CH31
    ##!#tax4 = fig.add_subplot(2,3,4,projection = mapcrs) # total SW
    ##!#tax5 = fig.add_subplot(2,3,5,projection = mapcrs) # total LW
    ##!#tax6 = fig.add_subplot(2,3,6,projection = mapcrs) # total CH31

    #print('old: ',hash_match_LWF.shape)
    total_omi_match_LAT = np.ma.masked_invalid(total_omi_match_LAT)
    total_omi_match_LON = np.ma.masked_invalid(total_omi_match_LON)
    total_omi_match_LWF = np.ma.masked_invalid(total_omi_match_LWF)
    total_omi_match_SWF = np.ma.masked_invalid(total_omi_match_SWF)

    ##!#hash_match_LAT = np.ma.masked_invalid(hash_match_LAT)
    ##!#hash_match_LON = np.ma.masked_invalid(hash_match_LON)
    ##!#hash_match_LWF = np.ma.masked_invalid(hash_match_LWF)
    ##!#hash_match_SWF = np.ma.masked_invalid(hash_match_SWF)

    ##!#nohash_match_LAT = np.ma.masked_invalid(nohash_match_LAT)
    ##!#nohash_match_LON = np.ma.masked_invalid(nohash_match_LON)
    ##!#nohash_match_LWF = np.ma.masked_invalid(nohash_match_LWF)
    ##!#nohash_match_SWF = np.ma.masked_invalid(nohash_match_SWF)

    mask_swf = np.ma.masked_invalid(mask_swf)
    mask_lwf = np.ma.masked_invalid(mask_lwf)
    mask_LAT = np.ma.masked_invalid(mask_LAT)
    mask_LON = np.ma.masked_invalid(mask_LON)
  
    ##!#total_match_SWF = np.ma.masked_invalid(total_match_SWF[~total_match_SWF.mask])
    ##!#total_match_LWF = np.ma.masked_invalid(total_match_LWF[~total_match_LWF.mask])
    ##!#total_match_LAT = np.ma.masked_invalid(total_match_LAT[~total_match_LAT.mask])
    ##!#total_match_LON = np.ma.masked_invalid(total_match_LON[~total_match_LON.mask])
    total_mask_swf  = np.ma.masked_invalid(mask_swf[~total_omi_match_SWF.mask])
    total_mask_lwf  = np.ma.masked_invalid(mask_lwf[~total_omi_match_LWF.mask])

    ##!#hash_match_SWF = np.ma.masked_invalid(hash_match_SWF[~hash_match_SWF.mask])
    ##!#hash_match_LWF = np.ma.masked_invalid(hash_match_LWF[~hash_match_LWF.mask])
    ##!#hash_match_LAT = np.ma.masked_invalid(hash_match_LAT[~hash_match_LAT.mask])
    ##!#hash_match_LON = np.ma.masked_invalid(hash_match_LON[~hash_match_LON.mask])
    ##!#hash_mask_swf  = np.ma.masked_invalid(mask_swf[~hash_match_SWF.mask])
    ##!#hash_mask_lwf  = np.ma.masked_invalid(mask_lwf[~hash_match_LWF.mask])
 
    ##!#nohash_match_SWF = np.ma.masked_invalid(nohash_match_SWF[~nohash_match_SWF.mask])
    ##!#nohash_match_LWF = np.ma.masked_invalid(nohash_match_LWF[~nohash_match_LWF.mask])
    ##!#nohash_match_LAT = np.ma.masked_invalid(nohash_match_LAT[~nohash_match_LAT.mask])
    ##!#nohash_match_LON = np.ma.masked_invalid(nohash_match_LON[~nohash_match_LON.mask])
    ##!#nohash_mask_swf  = np.ma.masked_invalid(mask_swf[~nohash_match_SWF.mask])
    ##!#nohash_mask_lwf  = np.ma.masked_invalid(mask_lwf[~nohash_match_LWF.mask])

    # ----------------------------------------------------------------------- 
    #
    # Plot the CERES and gridded OMI data on a temporary figure
    # Currently not saved anywhere.
    #
    # ----------------------------------------------------------------------- 
    ##!#plot_CERES_spatial(date_str, total_omi_match_LAT, total_omi_match_LON, \
    ##!#    total_mask_swf, 'SWF', tax1, zoom=True)
    ##!#plot_CERES_spatial(date_str, total_omi_match_LAT, total_omi_match_LON, \
    ##!#    total_mask_lwf, 'LWF', tax2, zoom=True)
    ##!#plot_CERES_spatial(date_str, total_omi_match_LAT, total_omi_match_LON, \
    ##!#    total_omi_match_LWF, 'SWF', tax3, vmin = np.nanmin(mask_UVAI),\
    ##!#    vmax = np.nanmax(mask_UVAI), zoom=True)
    ##!###!#plot_CERES_spatial(date_str, mask_LAT, mask_LON, \
    ##!###!#    mask_swf, 'SWF', tax4, zoom=True)
    ##!###!#plot_CERES_spatial(date_str, mask_LAT, mask_LON, \
    ##!###!#    mask_lwf, 'LWF', tax5, zoom=True)
    ##!#plot_OMI_spatial(date_str, LAT, LON, mask_UVAI, tax6, zoom = True)
 

    lrval_p = spearmanr(total_omi_match_LWF.compressed(), \
        total_mask_lwf)[0]

    #print(hash_match_LWF.compressed(), hash_mask_lwf)
    pax.scatter(total_omi_match_LWF.compressed(), total_mask_lwf,\
        s = 18, color='tab:orange', marker='o',label = 'LW')
    print("Calculating OMI/LWF trend")
    plot_trend_line(pax, total_omi_match_LWF.compressed(), total_mask_lwf, \
        color='tab:orange')
    ##!#pax.scatter(nohash_match_LWF.compressed(), nohash_mask_lwf,\
    ##!#    s = 18, color='tab:orange', marker='$L$',label = 'LWF clear')

    srval_p = spearmanr(total_omi_match_SWF.compressed(), \
        total_mask_swf)[0]
    pax.scatter(total_omi_match_SWF.compressed(), total_mask_swf,\
        s = 18, color='tab:blue', marker='o',label = 'SW')
    print("Calculating OMI/SWF trend")
    plot_trend_line(pax, total_omi_match_SWF.compressed(), total_mask_swf, \
        color='tab:blue')

    # ------------------------------------------------------------------------
    #
    # Calculate net forcing, which is found by
    #   Tot = SW_N + LW_N
    #   SW_N = SW_C - SW_A
    #   LW_N = LW_C - LW_A
    #   LW_C = average of all LW values for AI below 0.7
    #   SW_C = average of all SW values for AI below 0.7
    #
    # ------------------------------------------------------------------------
    clr_SW = np.nanmean(total_mask_swf[np.where(total_omi_match_SWF.compressed() < 0.7)])
    clr_LW = np.nanmean(total_mask_lwf[np.where(total_omi_match_SWF.compressed() < 0.7)])
    print("clear_SW = ",clr_SW," clear_LW = ",clr_LW)

    #net_SW = total_mask_swf - clr_SW
    #net_LW = total_mask_lwf - clr_LW
    net_SW = total_mask_swf
    net_LW = total_mask_lwf

    total_net = net_SW + net_LW

    mask_total = total_mask_swf + total_mask_lwf
    trval_p = spearmanr(total_omi_match_SWF.compressed(), \
        total_net)[0]
    ##!#mask_total = total_mask_swf + total_mask_lwf
    ##!#trval_p = spearmanr(total_omi_match_SWF.compressed(), \
    ##!#    mask_total)[0]


    #pax2 = pax.twinx()
    pax.scatter(total_omi_match_SWF.compressed(),total_net, s = 18, \
        color='tab:olive',marker='o', label = 'Total')
    print("Calculating OMI/Total Net trend")
    plot_trend_line(pax, total_omi_match_SWF.compressed(), total_net, \
        color='tab:olive')

    if(xlabel is None):
        xlabel = 'OMI AI'
    if(ylabel is None):
        ylabel = 'CERES TOA Flux [W/m2]'
    pax.set_xlabel(xlabel, fontsize = 14)
    pax.set_ylabel(ylabel, fontsize = 14)

    ##!#lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    ##!#lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    ##!#if(compare_OMI):
    ##!#    plt.legend(lines, labels, loc = 'lower center', bbox_to_anchor = (0, 0.01, 1, 1),\
    ##!#        bbox_transform = plt.gcf().transFigure, ncol=2)
    ##!#else:
    ##!#    plt.legend(lines, labels, loc = 'right', bbox_to_anchor = (0, 0., 1.0, 1),\
    ##!#        bbox_transform = plt.gcf().transFigure, ncol=1)

    pax.legend()
    #pax.set_title('Smoke correlation: '+str(np.round(lrval_p, 3)))
    if(ptitle is None):
        ptitle = 'OMI/SWF: '+str(np.round(srval_p,3)) + '  OMI/LWF: ' + \
        str(np.round(lrval_p,3)) + '\nOMI/Total:  ' + str(np.round(trval_p,3))
    pax.set_title(ptitle)

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
        pad=0.03)
    cbar3.set_label('OMI UVAI', size = 14, weight = 'bold')
    cbar3.ax.tick_params(labelsize = 12)
    pax.set_title('OMI UVAI',fontsize=10)


def plot_combined_imagery(date_str,channel1 = 1, channel2 = 5, channel3 = 31,\
        zoom=True,save=False,composite=True, show_smoke = False):
    
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
    mapcrs = init_proj(date_str)
    #mapcrs = ccrs.LambertConformal(central_longitude = \
    #    np.mean(plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['Lon']),\
    #    central_latitude = \
    #    np.mean(plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['Lat']))

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
        ax5.set_title('Aqua GOES\n' + prev_date_str.strftime('%Y-%m-%d %H:%M'))

    elif((date_str == '202108062025') | (date_str == '202108052120')):

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

        LAT, LON, mask_UVAI = read_OMI_match_GOES(date_str)
 
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

    ax0.set_title('Aqua GOES\n'+dt_date_str.strftime('%Y-%m-%d %H:%M'))

    # ----------------------------------------------------------------------
    #
    # Read the GOES data
    #
    # ----------------------------------------------------------------------

    # Call read_GOES_channel to read the desired GOES data from the file 
    # and put it in a dictionary
    # ---------------------------------------------------------------------
    GOES_data1 = read_GOES_channel(dt_date_str.strftime("%Y%m%d%H%M"), \
        channel1, zoom = True)
    GOES_data2 = read_GOES_channel(dt_date_str.strftime("%Y%m%d%H%M"), \
        channel2, zoom = True)
    GOES_data3 = read_GOES_channel(dt_date_str.strftime("%Y%m%d%H%M"), \
        channel3, zoom = True)
    GOES_dataw = read_GOES_channel(dt_date_str.strftime("%Y%m%d%H%M"), \
        'wv_ir', zoom = True)

    # ----------------------------------------------------------------------
    #
    # Plot the GOES data
    #
    # ----------------------------------------------------------------------
    
    # Plot channels 1, 5, 31, and WV data
    plot_GOES_spatial(GOES_data1, ax1, zoom = zoom)
    plot_GOES_spatial(GOES_data2, ax2, zoom = zoom)
    plot_GOES_spatial(GOES_data3, ax3, zoom = zoom)
    plot_GOES_spatial(GOES_dataw, ax4, zoom = zoom)

    # ----------------------------------------------------------------------
    #
    # Read the CERES data
    #
    # ----------------------------------------------------------------------

    mask_LAT, mask_LON, mask_swf, mask_lwf = \
        read_CERES_match_GOES(date_str)

    # ----------------------------------------------------------------------
    #
    # Plot the CERES data
    #
    # ----------------------------------------------------------------------

    plot_CERES_spatial(date_str,mask_LAT, mask_LON, mask_swf, 'SWF', ax6, zoom = zoom)
    plot_CERES_spatial(date_str,mask_LAT, mask_LON, mask_lwf, 'LWF', ax7, zoom = zoom)

    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(dt_date_str.strftime('%Y%m%d%H%M')) 

        ax0.pcolor(GOES_data1['lon'],GOES_data1['lat'],\
            hash_data1, hatch = '////', alpha=0., transform = datacrs)
        ax1.pcolor(GOES_data1['lon'],GOES_data1['lat'],\
            hash_data1, hatch = '////', alpha=0., transform = datacrs)
        ax2.pcolor(GOES_data1['lon'],GOES_data1['lat'],\
            hash_data1, hatch = '////', alpha=0., transform = datacrs)
        ax3.pcolor(GOES_data1['lon'],GOES_data1['lat'],\
            hash_data1, hatch = '////', alpha=0., transform = datacrs)
        ax4.pcolor(GOES_data1['lon'],GOES_data1['lat'],\
            hash_data1, hatch = '////', alpha=0., transform = datacrs)
        ax5.pcolor(GOES_data1['lon'],GOES_data1['lat'],\
            hash_data1, hatch = '////', alpha=0., transform = datacrs)
        ax6.pcolor(GOES_data1['lon'],GOES_data1['lat'],\
            hash_data1, hatch = '////', alpha=0., transform = datacrs)
        ax7.pcolor(GOES_data1['lon'],GOES_data1['lat'],\
            hash_data1, hatch = '////', alpha=0., transform = datacrs)

    fig.tight_layout()

    if(save):
        outname = 'goes_combined_imagery_' + date_str + zoom_add + '.png'
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

def plot_combined_scatter(date_str,channel0 = 31, channel1 = 1, channel2 = 5,\
        zoom=True,save=False,composite=True,avg_pixel=False,plume_only=False):

    
    plt.close('all')
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # ----------------------------------------------------------------------
    #
    # Read the GOES data
    #
    # ----------------------------------------------------------------------

    # Call read_GOES_channel to read the desired GOES data from the file 
    # and put it in a dictionary
    # ---------------------------------------------------------------------
    GOES_data0 = read_GOES_channel(dt_date_str.strftime("%Y%m%d%H%M"), \
        channel0, zoom = zoom)
    GOES_data1 = read_GOES_channel(dt_date_str.strftime("%Y%m%d%H%M"), \
        channel1, zoom = zoom)
    GOES_data2 = read_GOES_channel(dt_date_str.strftime("%Y%m%d%H%M"), \
        channel2, zoom = zoom)
    GOES_dataw = read_GOES_channel(dt_date_str.strftime("%Y%m%d%H%M"), \
        'wv_ir', zoom = zoom)
    
    # ----------------------------------------------------------------------
    #
    # Read the CERES data
    #
    # ----------------------------------------------------------------------
    mask_LAT, mask_LON, mask_swf, mask_lwf = \
        read_CERES_match_GOES(date_str)

    # ------------------------------------------------------------------------
    #
    # For 20210722:
    #
    # ------------------------------------------------------------------------
    if((date_str == '202107222110')):
        # Read true color data for the previous date

        fig = plt.figure(figsize=(9,9))
        ax0 = fig.add_subplot(2,2,1) # 31 vs 1
        ax1 = fig.add_subplot(2,2,2) # 31 vs 5
        ax2 = fig.add_subplot(2,2,3) # 31 vs WV
        ax3 = fig.add_subplot(2,2,4) # 31 vs SW/LW/Total

    elif((date_str == '202108062025') | (date_str == '202108052120')):

        fig = plt.figure(figsize=(12,14))
        ax0 = fig.add_subplot(3,2,1) # 31 vs 1
        ax1 = fig.add_subplot(3,2,2) # 31 vs 5
        ax2 = fig.add_subplot(3,2,3) # 31 vs WV
        ax3 = fig.add_subplot(3,2,4) # 31 vs SW/LW/Total
        ax4 = fig.add_subplot(3,2,5) # 31 vs AI
        ax5 = fig.add_subplot(3,2,6) # AI vs Sw/LW/Total
   
        # ---------------------------------------------------------------------
        #
        # Read OMI data and plot Ch 31 vs OMI AI
        #
        # ---------------------------------------------------------------------
        plot_scatter_OMI(date_str, GOES_data0, ax4, avg_pixel = avg_pixel)

        # ----------------------------------------------------------------------
        #
        # Colocate the OMI and CERES data, then plot
        #
        # ----------------------------------------------------------------------
        plot_scatter_OMI_CERES(date_str, GOES_data0, ax5, \
            avg_pixel = avg_pixel, plume_only = plume_only)
    
        plot_subplot_label(ax4, '(e)')
        plot_subplot_label(ax5, '(f)')

    # Determine where the smoke is located
    # ------------------------------------
    hash_data, nohash_data = find_plume(dt_date_str.strftime('%Y%m%d%H%M')) 

    tmp_data0 = np.copy(GOES_data0['data'])
    tmp_data1 = np.copy(GOES_data1['data'])
    tmp_data2 = np.copy(GOES_data2['data'])
    tmp_dataw = np.copy(GOES_dataw['data'])
    tmp_lat0  = np.copy(GOES_data0['lat'])
    tmp_lon0  = np.copy(GOES_data0['lon'])

    if(not (tmp_data0.shape == tmp_data1.shape == tmp_data2.shape == \
            tmp_dataw.shape == hash_data.shape)):
        print("shape mismatch")
        shapes = []
        shapes.append(tmp_data0.shape)
        shapes.append(tmp_data1.shape)
        shapes.append(tmp_data2.shape)
        shapes.append(tmp_dataw.shape)
        shapes.append(hash_data.shape)

        min_shape = min(shapes)
        print(min_shape)

        tmp_data0  = tmp_data0[:min_shape[0],:min_shape[1]]
        tmp_data1  = tmp_data1[:min_shape[0],:min_shape[1]]
        tmp_data2  = tmp_data2[:min_shape[0],:min_shape[1]]
        tmp_dataw  = tmp_dataw[:min_shape[0],:min_shape[1]]
        tmp_lat0   = tmp_lat0[:min_shape[0],:min_shape[1]]
        tmp_lon0   = tmp_lon0[:min_shape[0],:min_shape[1]]
        hash_data  = hash_data[:min_shape[0],:min_shape[1]]

    max_ch = 350.

    tmp_data0 = np.ma.masked_where( (abs(tmp_data1) > max_ch) | \
        (abs(tmp_data0) > max_ch) | (abs(tmp_data2) > max_ch) | \
        (abs(tmp_dataw) > max_ch), tmp_data0)
    tmp_data1 = np.ma.masked_where( (abs(tmp_data1) > max_ch) | \
        (abs(tmp_data0) > max_ch) | (abs(tmp_data2) > max_ch) | \
        (abs(tmp_dataw) > max_ch), tmp_data1)
    tmp_data2 = np.ma.masked_where( (abs(tmp_data1) > max_ch) | \
        (abs(tmp_data0) > max_ch) | (abs(tmp_data2) > max_ch) | \
        (abs(tmp_dataw) > max_ch), tmp_data2)
    tmp_dataw = np.ma.masked_where( (abs(tmp_data1) > max_ch) | \
        (abs(tmp_data0) > max_ch) | (abs(tmp_data2) > max_ch) | \
        (abs(tmp_dataw) > max_ch), tmp_dataw)
    ##!#hash_data = np.ma.masked_where( (abs(tmp_data1) > max_ch) | \
    ##!#    (abs(tmp_data0) > max_ch) | \
    ##!#    (abs(tmp_data2) > max_ch) | (abs(tmp_dataw) > max_ch) | \
    ##!#    (abs(hash_data) > max_ch), hash_data)
    tmp_lat0 = np.ma.masked_where( (abs(tmp_data1) > max_ch) | \
        (abs(tmp_data0) > max_ch) | (abs(tmp_data2) > max_ch) | \
        (abs(tmp_dataw) > max_ch), tmp_lat0)
    tmp_lon0 = np.ma.masked_where( (abs(tmp_data1) > max_ch) | \
        (abs(tmp_data0) > max_ch) | (abs(tmp_data2) > max_ch) | \
        (abs(tmp_dataw) > max_ch), tmp_lon0)


    # ----------------------------------------------------------------------
    #
    # Plot the GOES data
    #
    # ----------------------------------------------------------------------

    plot_scatter(ax0, tmp_data0, tmp_data1, GOES_data0, GOES_data1, \
        hash_data)
    plot_scatter(ax1, tmp_data0, tmp_data2, GOES_data0, GOES_data2, \
        hash_data)
    plot_scatter(ax2, tmp_data0, tmp_dataw, GOES_data0, GOES_dataw, \
        hash_data)

    # ----------------------------------------------------------------------
    #
    # Plot the CERES data
    #
    # ----------------------------------------------------------------------
    plot_scatter_CERES(date_str, GOES_data0, ax3, avg_pixel = avg_pixel,\
        plume_only = plume_only)

    # Add labels to all the subplots
    # ------------------------------
    plot_subplot_label(ax0, '(a)')
    plot_subplot_label(ax1, '(b)')
    plot_subplot_label(ax2, '(c)')
    plot_subplot_label(ax3, '(d)')

    fig.tight_layout()

    ##!#plot_CERES_spatial(date_str,mask_LAT, mask_LON, mask_swf, 'SWF', ax6, zoom = zoom)
    ##!#plot_CERES_spatial(date_str,mask_LAT, mask_LON, mask_lwf, 'LWF', ax7, zoom = zoom)
    if(save):
        plume_add = '_onlyplume'
        pixel_add = ''
        if(plume_only == False):
            plume_add = '_onlyplume'
        if(avg_pixel):
            pixel_add = '_avgpixel' 
        outname = 'goes_combined_scatter_' + date_str + plume_add + \
            pixel_add + '.png'
        fig.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

def plot_scatter_select(composite = True, avg_pixel = True, \
        plume_only = False,zoom = True, save=False):
    
    plt.close('all')
    date_str_July = '202107222110'
    date_str_Aug  = '202108062025'
    dt_date_July = datetime.strptime(date_str_July,"%Y%m%d%H%M")
    dt_date_Aug  = datetime.strptime(date_str_Aug,"%Y%m%d%H%M")

    # ----------------------------------------------------------------------
    #
    # Read the GOES data
    #
    # ----------------------------------------------------------------------

    # Call read_GOES_channel to read the desired GOES data from the file 
    # and put it in a dictionary
    # ---------------------------------------------------------------------
    GOES_dataJul31 = read_GOES_channel(date_str_July,31, zoom = zoom)
    GOES_dataAug31 = read_GOES_channel(date_str_Aug, 31, zoom = zoom)

    ##!## Determine where the smoke is located
    ##!## ------------------------------------
    ##!#hash_data, nohash_data = find_plume(dt_date_str.strftime('%Y%m%d%H%M')) 

    ##!#tmp_data0 = np.copy(GOES_data0['data'])
    ##!#tmp_data1 = np.copy(GOES_data1['data'])
    ##!#tmp_data2 = np.copy(GOES_data2['data'])
    ##!#tmp_dataw = np.copy(GOES_dataw['data'])
    ##!#tmp_lat0  = np.copy(GOES_data0['lat'])
    ##!#tmp_lon0  = np.copy(GOES_data0['lon'])

    ##!#if(not (tmp_data0.shape == tmp_data1.shape == tmp_data2.shape == \
    ##!#        tmp_dataw.shape == hash_data.shape)):
    ##!#    print("shape mismatch")
    ##!#    shapes = []
    ##!#    shapes.append(tmp_data0.shape)
    ##!#    shapes.append(tmp_data1.shape)
    ##!#    shapes.append(tmp_data2.shape)
    ##!#    shapes.append(tmp_dataw.shape)
    ##!#    shapes.append(hash_data.shape)

    ##!#    min_shape = min(shapes)
    ##!#    print(min_shape)

    ##!#    tmp_data0  = tmp_data0[:min_shape[0],:min_shape[1]]
    ##!#    tmp_data1  = tmp_data1[:min_shape[0],:min_shape[1]]
    ##!#    tmp_data2  = tmp_data2[:min_shape[0],:min_shape[1]]
    ##!#    tmp_dataw  = tmp_dataw[:min_shape[0],:min_shape[1]]
    ##!#    tmp_lat0   = tmp_lat0[:min_shape[0],:min_shape[1]]
    ##!#    tmp_lon0   = tmp_lon0[:min_shape[0],:min_shape[1]]
    ##!#    hash_data  = hash_data[:min_shape[0],:min_shape[1]]

    ##!#max_ch = 350.

    ##!#tmp_data0 = np.ma.masked_where( (abs(tmp_data1) > max_ch) | \
    ##!#    (abs(tmp_data0) > max_ch) | (abs(tmp_data2) > max_ch) | \
    ##!#    (abs(tmp_dataw) > max_ch), tmp_data0)
    ##!#tmp_data1 = np.ma.masked_where( (abs(tmp_data1) > max_ch) | \
    ##!#    (abs(tmp_data0) > max_ch) | (abs(tmp_data2) > max_ch) | \
    ##!#    (abs(tmp_dataw) > max_ch), tmp_data1)
    ##!#tmp_data2 = np.ma.masked_where( (abs(tmp_data1) > max_ch) | \
    ##!#    (abs(tmp_data0) > max_ch) | (abs(tmp_data2) > max_ch) | \
    ##!#    (abs(tmp_dataw) > max_ch), tmp_data2)
    ##!#tmp_dataw = np.ma.masked_where( (abs(tmp_data1) > max_ch) | \
    ##!#    (abs(tmp_data0) > max_ch) | (abs(tmp_data2) > max_ch) | \
    ##!#    (abs(tmp_dataw) > max_ch), tmp_dataw)
    ##!###!#hash_data = np.ma.masked_where( (abs(tmp_data1) > max_ch) | \
    ##!###!#    (abs(tmp_data0) > max_ch) | \
    ##!###!#    (abs(tmp_data2) > max_ch) | (abs(tmp_dataw) > max_ch) | \
    ##!###!#    (abs(hash_data) > max_ch), hash_data)
    ##!#tmp_lat0 = np.ma.masked_where( (abs(tmp_data1) > max_ch) | \
    ##!#    (abs(tmp_data0) > max_ch) | (abs(tmp_data2) > max_ch) | \
    ##!#    (abs(tmp_dataw) > max_ch), tmp_lat0)
    ##!#tmp_lon0 = np.ma.masked_where( (abs(tmp_data1) > max_ch) | \
    ##!#    (abs(tmp_data0) > max_ch) | (abs(tmp_data2) > max_ch) | \
    ##!#    (abs(tmp_dataw) > max_ch), tmp_lon0)
    
    # ----------------------------------------------------------------------
    #
    # Read the CERES data
    #
    # ----------------------------------------------------------------------
    mask_LAT_July, mask_LON_July, mask_swf_July, mask_lwf_July = \
        read_CERES_match_GOES(date_str_July)
    mask_LAT_Aug, mask_LON_Aug, mask_swf_Aug, mask_lwf_Aug = \
        read_CERES_match_GOES(date_str_Aug)

    # ------------------------------------------------------------------------
    #
    # Set up the figure panels
    #
    # ------------------------------------------------------------------------
    fig = plt.figure(figsize=(9,9))
    ax0 = fig.add_subplot(2,2,1) # 0722 Ch 31 vs SW/LW
    ax1 = fig.add_subplot(2,2,2) # 0806 Ch 31 vs SW/LW
    ax2 = fig.add_subplot(2,2,3) # 0806 Ch 31 vs AI
    ax3 = fig.add_subplot(2,2,4) # 0806 AI vs SW/LW

    # ---------------------------------------------------------------------
    #
    # Read OMI data and plot Ch 31 vs OMI AI
    #
    # ---------------------------------------------------------------------
    plot_scatter_OMI(date_str_Aug, GOES_dataAug31, ax2, avg_pixel = avg_pixel)

    # ----------------------------------------------------------------------
    #
    # Colocate the OMI and CERES data, then plot
    #
    # ----------------------------------------------------------------------
    plot_scatter_OMI_CERES(date_str_Aug, GOES_dataAug31, ax3, \
        avg_pixel = avg_pixel, plume_only = plume_only)

    # ----------------------------------------------------------------------
    #
    # Plot the CERES data
    #
    # ----------------------------------------------------------------------
    plot_scatter_CERES(date_str_July, GOES_dataJul31, ax0, \
        avg_pixel = avg_pixel,plume_only = plume_only)
    plot_scatter_CERES(date_str_Aug, GOES_dataAug31, ax1, \
        avg_pixel = avg_pixel,plume_only = plume_only)

    fig.tight_layout()

    ##!#plot_CERES_spatial(date_str,mask_LAT, mask_LON, mask_swf, 'SWF', ax6, zoom = zoom)
    ##!#plot_CERES_spatial(date_str,mask_LAT, mask_LON, mask_lwf, 'LWF', ax7, zoom = zoom)
    if(save):
        plume_add = '_onlyplume'
        pixel_add = ''
        if(plume_only == False):
            plume_add = '_onlyplume'
        if(avg_pixel):
            pixel_add = '_avgpixel' 
        outname = 'goes_combined_scatter_select' + plume_add + \
            pixel_add + '.png'
        fig.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

def plot_combined_figure1(zoom = True, show_smoke = True, composite = True, \
        save=False):

    date_str = '202107222110'
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # ----------------------------------------------------------------------
    #
    # Read the GOES and CERES data
    #
    # ----------------------------------------------------------------------

    # Call read_GOES_channel to read the desired GOES data from the file 
    # and put it in a dictionary
    # ---------------------------------------------------------------------
    GOES_data_ch1  = read_GOES_channel(date_str, 1, zoom = zoom)
    GOES_data_ch5  = read_GOES_channel(date_str, 5, zoom = zoom)
    GOES_data_ch31 = read_GOES_channel(date_str, 31, zoom = zoom)
    GOES_data_wv   = read_GOES_channel(date_str, 'wv_ir', zoom = zoom)

    # Determine where the smoke is located
    # ------------------------------------
    hash_data, nohash_data = find_plume(dt_date_str.strftime('%Y%m%d%H%M')) 

    tmp_data1  = np.copy(GOES_data_ch1['data'])
    tmp_data5  = np.copy(GOES_data_ch5['data'])
    tmp_data31 = np.copy(GOES_data_ch31['data'])
    tmp_datawv = np.copy(GOES_data_wv['data'])
    tmp_lat0   = np.copy(GOES_data_ch1['lat'])
    tmp_lon0   = np.copy(GOES_data_ch1['lon'])

    if(not (tmp_data1.shape == tmp_data5.shape == tmp_data31.shape == \
            tmp_datawv.shape == hash_data.shape)):
        print("shape mismatch")
        shapes = []
        shapes.append(tmp_data1.shape)
        shapes.append(tmp_data5.shape)
        shapes.append(tmp_data31.shape)
        shapes.append(tmp_datawv.shape)
        shapes.append(hash_data.shape)

        min_shape = min(shapes)
        print(min_shape)

        tmp_data1  = tmp_data1[:min_shape[0],:min_shape[1]]
        tmp_data5  = tmp_data5[:min_shape[0],:min_shape[1]]
        tmp_data31 = tmp_data31[:min_shape[0],:min_shape[1]]
        tmp_datawv = tmp_datawv[:min_shape[0],:min_shape[1]]
        tmp_lat0   = tmp_lat0[:min_shape[0],:min_shape[1]]
        tmp_lon0   = tmp_lon0[:min_shape[0],:min_shape[1]]
        hash_data  = hash_data[:min_shape[0],:min_shape[1]]

    max_ch = 350.

    tmp_data1 = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) | \
        (abs(tmp_datawv) > max_ch), tmp_data1)
    tmp_data5 = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) | \
        (abs(tmp_datawv) > max_ch), tmp_data5)
    tmp_data31 = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) | \
        (abs(tmp_datawv) > max_ch), tmp_data31)
    tmp_datawv = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) | \
        (abs(tmp_datawv) > max_ch), tmp_datawv)
    tmp_lat0 = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) | \
        (abs(tmp_datawv) > max_ch), tmp_lat0)
    tmp_lon0 = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) | \
        (abs(tmp_datawv) > max_ch), tmp_lon0)

    # Read the true color data
    # ------------------------
    var1, crs1, lat_lims1, lon_lims1 = read_true_color(date_str,\
        composite=composite)

    # Read in the CERES data
    # ----------------------
    mask_LAT, mask_LON, mask_swf, mask_lwf = \
        read_CERES_match_GOES(date_str)
   
    # ----------------------------------------------------------------------
    #
    #  Set up the 9-panel figure
    #
    # ----------------------------------------------------------------------
 
    mapcrs = init_proj(date_str)
    plt.close('all')
    fig = plt.figure(figsize=(12,10))
    ax1 = fig.add_subplot(3,3,1,projection = crs1)   # true color    
    ax2 = fig.add_subplot(3,3,2,projection = mapcrs) # Ch 31
    ax3 = fig.add_subplot(3,3,3,projection = mapcrs) # Ch 1
    ax4 = fig.add_subplot(3,3,4,projection = mapcrs) # Ch 5
    ax5 = fig.add_subplot(3,3,5)                     # Ch 31 vs 1 scatter
    ax6 = fig.add_subplot(3,3,6)                     # Ch 31 vs WV scatter
    ax7 = fig.add_subplot(3,3,7,projection = mapcrs) # WV
    ax8 = fig.add_subplot(3,3,8,projection = mapcrs) # LW
    ax9 = fig.add_subplot(3,3,9,projection = mapcrs) # SW

    # Plot the true-color data for the previous date
    # ----------------------------------------------
    ax1.imshow(var1.data, transform = crs1, extent=(var1.x[0], var1.x[-1], \
        var1.y[-1], var1.y[0]), origin='upper')

    # Zoom in the figure if desired
    # -----------------------------
    if(zoom):
        ax1.set_extent([lon_lims1[0],lon_lims1[1],lat_lims1[0],\
            lat_lims1[1]],crs = datacrs)

    # ----------------------------------------------------------------------
    #
    # Plot the data in the figure
    #
    # ----------------------------------------------------------------------

    # Plot channel 1, 5, 31, and WV data spatial data
    # -----------------------------------------------
    plot_GOES_spatial(GOES_data_ch31, ax2, zoom = zoom, ptitle = '')
    plot_GOES_spatial(GOES_data_ch1,  ax3, zoom = zoom, ptitle = '')
    plot_GOES_spatial(GOES_data_ch5,  ax4, zoom = zoom, ptitle = '')
    plot_GOES_spatial(GOES_data_wv,   ax7, zoom = zoom, ptitle = '')

    # Plot the GOES scatter data
    # ---------------------------
    plot_scatter(ax5, tmp_data31, tmp_data1, GOES_data_ch31, GOES_data_ch1, \
        hash_data, xlabel = '11 μm brightness temperature', \
        ylabel = '0.64 μm reflectance', plot_legend = True)
    plot_scatter(ax6, tmp_data31, tmp_datawv, GOES_data_ch31, GOES_data_wv, \
        hash_data, xlabel = '11 μm brightness temperature', \
        ylabel = '1.24 μm reflectance')

    # Plot the CERES SWF and LWF data
    # -------------------------------
    plot_CERES_spatial(date_str, mask_LAT, mask_LON, mask_swf, 'SWF', ax8, \
        ptitle = '', zoom = zoom)
    plot_CERES_spatial(date_str, mask_LAT, mask_LON, mask_lwf, 'LWF', ax9, \
        ptitle = '', zoom = zoom)

    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(date_str) 

        plt.rcParams.update({'hatch.color': 'r'})
        ax3.pcolor(GOES_data_ch1['lon'],GOES_data_ch1['lat'],\
            hash_data1, hatch = '\\\\', alpha=0., transform = datacrs,\
            cmap = 'plasma')
    # Plot the ASOS locations on the map
    # ----------------------------------
    plot_ASOS_locs(ax1,date_str,color='lime', sites = ['O05','AAT'])

    # Add subplot labels
    # ------------------
    plot_subplot_label(ax1, '(a)', backgroundcolor = 'white')
    plot_subplot_label(ax2, '(b)', backgroundcolor = 'white')
    plot_subplot_label(ax3, '(c)', backgroundcolor = 'white')
    plot_subplot_label(ax4, '(d)', backgroundcolor = 'white')
    plot_subplot_label(ax5, '(e)')
    plot_subplot_label(ax6, '(f)')
    plot_subplot_label(ax7, '(g)', backgroundcolor = 'white')
    plot_subplot_label(ax8, '(h)', backgroundcolor = 'white')
    plot_subplot_label(ax9, '(i)', backgroundcolor = 'white')

    # Add plot text
    # -------------
    plot_figure_text(ax2, 'GOES 11 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax3, 'GOES 0.64 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax4, 'GOES 1.24 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax7, 'GOES IR WV', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax8, 'CERES SW', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax9, 'CERES LW', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')

    fig.tight_layout()

    GOES_data_ch1.clear()
    GOES_data_ch5.clear()
    GOES_data_ch31.clear()
    GOES_data_wv.clear()

    if(save):
        outname = 'goes_total_combined_' + date_str + '_fig1.png'
        fig.savefig(outname, dpi=300)
        print("Saved",outname)
    else:
        plt.show()
  
def plot_figure2(save=False, composite = True, calc_radiance = True, \
        satellite = 'goes_ch31'):

    if('/home/bsorenson/Research/SBDART' not in sys.path):
        sys.path.append('/home/bsorenson/Research/SBDART')
    from SBDART_Lib import run_sbdart, plot_bright_vza

    # Run SBDART for the different channel
    # ------------------------------------
    goes31 = run_sbdart(satellite, calc_radiance, run = True)

    # Read true color data for the date
    # ---------------------------------
    date_str22 = '202107222110'
    dt_date_str22 = datetime.strptime(date_str22,"%Y%m%d%H%M")
    var3, crs3, lat_lims3, lon_lims3 = read_true_color(date_str22,\
        composite=composite)

    # Set up the figure
    # -----------------
    fig = plt.figure(figsize=(8,12))
    gs1 = fig.add_gridspec(nrows = 3, ncols = 1, wspace = 0.10, hspace = 0.40)
    ax1 = fig.add_subplot(gs1[0], projection = crs3) # true color 7/22
    ax2 = fig.add_subplot(gs1[1]) # meteo 7/22
    ax3 = fig.add_subplot(gs1[2])
    #ax1 = fig.add_subplot(1,2,1,projection = crs3) # true color 7/22
    #ax2 = fig.add_subplot(1,2,2) # meteo 7/22

    # ----------------------------------------------------------------------
    #
    # Panel 1: Map
    #
    # ----------------------------------------------------------------------

    # Plot the true-color data for the case date
    # ------------------------------------------
    ax1.imshow(var3.data, transform = crs3, extent=(var3.x[0], var3.x[-1], \
        var3.y[-1], var3.y[0]), origin='upper')

    ax1.set_extent([lon_lims3[0],lon_lims3[1],lat_lims3[0],\
        lat_lims3[1]],crs = datacrs)

    # Plot the ASOS locations on the map
    # ----------------------------------
    plot_ASOS_locs(ax1,date_str22,color='lime', sites = ['O05','AAT'])

    # ----------------------------------------------------------------------
    #
    # Panel 2: Meteogram
    #
    # ----------------------------------------------------------------------
    plot_asos_diurnal(ax2, date_str22, 'O05', 'AAT')

    def plot_goes_line(dt_date, pax):
        local_goes_time = dt_date - timedelta(hours = 7)
        goes_diff = (local_goes_time - datetime(year=2021,month=7,\
            day=22,hour=0,minute=0)).seconds
        print(local_goes_time, goes_diff)
        pax.axvline(goes_diff,color='black',linestyle = '--', lw=2,alpha=0.75,\
            label='GOES')

    plot_goes_line(dt_date_str22, ax2)
    ax2.legend() 

    # ----------------------------------------------------------------------
    #
    # Panel 3: SBDART stuff
    #
    # ----------------------------------------------------------------------
    plot_bright_vza(goes31, pax = ax3)
    
   
    # Add subplot labels
    # ------------------
    dt_date_str = dt_date_str22
    xval = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lon'][0] + \
        (plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lon'][1] - \
        plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lon'][0])*0.05
    yval = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lat'][0] + \
        (plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lat'][1] - \
        plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lat'][0])*0.90

    plot_subplot_label(ax1, '(a)', xval = xval, yval = yval, \
        transform = datacrs, backgroundcolor = 'white')
    plot_subplot_label(ax2, '(b)', xval = 1200, yval = 5)
    plot_subplot_label(ax3, '(c)', xval = 17, location = 'lower_left')

    fig.tight_layout()
 
    if(save):
        outname = 'goes_asos_meteo_combined_20210722_fig2_v2.png'
        fig.savefig(outname, dpi=300)
        print("Saved image",outname)
    else: 
        plt.show() 

 
def plot_combined_figure3(zoom = True, show_smoke = False, composite = True, \
        plume_only = False, avg_pixel = True, save=False):

    date_str = '202108062025'
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # ----------------------------------------------------------------------
    #
    # Read the GOES and CERES data
    #
    # ----------------------------------------------------------------------

    # Call read_GOES_channel to read the desired GOES data from the file 
    # and put it in a dictionary
    # ---------------------------------------------------------------------
    GOES_data_ch5  = read_GOES_channel(date_str, 5, zoom = zoom)
    GOES_data_ch31 = read_GOES_channel(date_str, 31, zoom = zoom)
    GOES_data_wv   = read_GOES_channel(date_str, 'wv_ir', zoom = zoom)

    # Determine where the smoke is located MAYBE???
    # ---------------------------------------------

    # Read the true color data
    # ------------------------
    var1, crs1, lat_lims1, lon_lims1 = read_true_color(date_str,\
        composite=composite)

    # Read in the CERES data
    # ----------------------
    mask_LAT, mask_LON, mask_swf, mask_lwf = \
        read_CERES_match_GOES(date_str)
  
    # Read in the OMI data
    # --------------------
    LAT, LON, mask_UVAI = read_OMI_match_GOES(date_str)
 
 
    # ----------------------------------------------------------------------
    #
    #  Set up the 9-panel figure
    #
    # ----------------------------------------------------------------------
 
    mapcrs = init_proj(date_str)
    plt.close('all')
    fig = plt.figure(figsize=(13,10))
    ax1 = fig.add_subplot(3,3,1,projection = crs1)   # true color    
    ax2 = fig.add_subplot(3,3,2,projection = mapcrs) # Ch 31
    ax3 = fig.add_subplot(3,3,3,projection = mapcrs) # OMI AI
    ax4 = fig.add_subplot(3,3,4,projection = mapcrs) # Ch 5
    ax5 = fig.add_subplot(3,3,5,projection = mapcrs) # CERES SW
    ax6 = fig.add_subplot(3,3,6,projection = mapcrs) # CERES LW
    ax7 = fig.add_subplot(3,3,7,projection = mapcrs) # WV
    ax8 = fig.add_subplot(3,3,8)                     # AI scatter
    ax9 = fig.add_subplot(3,3,9)                     # AI vs flux scatter

    # Plot the true-color data for the previous date
    # ----------------------------------------------
    ax1.imshow(var1.data, transform = crs1, extent=(var1.x[0], var1.x[-1], \
        var1.y[-1], var1.y[0]), origin='upper')

    # Zoom in the figure if desired
    # -----------------------------
    if(zoom):
        ax1.set_extent([lon_lims1[0],lon_lims1[1],lat_lims1[0],\
            lat_lims1[1]],crs = datacrs)

    # ----------------------------------------------------------------------
    #
    # Plot the data in the figure
    #
    # ----------------------------------------------------------------------

    # Plot channel 5, 31, and WV data spatial data
    # -----------------------------------------------
    plot_GOES_spatial(GOES_data_ch31, ax2, zoom = zoom, ptitle = '')
    plot_GOES_spatial(GOES_data_ch5,  ax4, zoom = zoom, ptitle = '')
    plot_GOES_spatial(GOES_data_wv,   ax7, zoom = zoom, ptitle = '')

    # Plot the OMI spatial data
    # -------------------------
    plot_OMI_spatial(date_str, LAT, LON, mask_UVAI, ax3, zoom = zoom)

    # Plot the Ch 31 v OMI and OMI vs flux scatter data
    # -------------------------------------------------
    plot_scatter_OMI(date_str, GOES_data_ch31, ax8, avg_pixel = avg_pixel, \
        xlabel = '11 μm brightness temperature', \
        ylabel = 'OMI UVAI', ptitle = '')
    plot_scatter_OMI_CERES(date_str, GOES_data_ch31, ax9, \
        avg_pixel = avg_pixel, plume_only = plume_only, ptitle = '')

    ##!#plot_scatter(ax5, tmp_data31, tmp_data1, GOES_data_ch31, GOES_data_ch1, \
    ##!#    hash_data, xlabel = '11 μm brightness temperature', \
    ##!#    ylabel = '0.64 μm reflectance', plot_legend = True)
    ##!#plot_scatter(ax6, tmp_data31, tmp_datawv, GOES_data_ch31, GOES_data_wv, \
    ##!#    hash_data, xlabel = '11 μm brightness temperature', \
    ##!#    ylabel = '1.24 μm reflectance')

    # Plot the CERES SWF and LWF data
    # -------------------------------
    plot_CERES_spatial(date_str, mask_LAT, mask_LON, mask_swf, 'SWF', ax5, \
        ptitle = '', zoom = zoom)
    plot_CERES_spatial(date_str, mask_LAT, mask_LON, mask_lwf, 'LWF', ax6, \
        ptitle = '', zoom = zoom)

    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(date_str) 

        ax3.pcolor(GOES_data_ch1['lon'],GOES_data_ch1['lat'],\
            hash_data1, hatch = '////', alpha=0., transform = datacrs,\
            color = 'lime')

    # Add subplot labels
    # ------------------
    plot_subplot_label(ax1, '(a)', backgroundcolor = 'white')
    plot_subplot_label(ax2, '(b)', backgroundcolor = 'white')
    plot_subplot_label(ax3, '(c)', backgroundcolor = 'white')
    plot_subplot_label(ax4, '(d)', backgroundcolor = 'white')
    plot_subplot_label(ax5, '(e)', backgroundcolor = 'white')
    plot_subplot_label(ax6, '(f)', backgroundcolor = 'white')
    plot_subplot_label(ax7, '(g)', backgroundcolor = 'white')
    plot_subplot_label(ax8, '(h)')
    plot_subplot_label(ax9, '(i)')

    # Add plot text
    # -------------
    plot_figure_text(ax2, 'GOES 11 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax3, 'OMI AI', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax4, 'GOES 1.24 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax5, 'CERES SW', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax6, 'CERES LW', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax7, 'GOES IR WV', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')

    fig.tight_layout()

    GOES_data_ch5.clear()
    GOES_data_ch31.clear()
    GOES_data_wv.clear()

    if(save):
        outname = 'goes_total_combined_' + date_str + '_fig3.png'
        fig.savefig(outname, dpi=300)
        print("Saved",outname)
    else:
        plt.show()
    
def plot_figureS1(save=False, composite = True):

    # Read true color data for the previous date
    date_str20 = '202107202125'
    date_str21 = '202107212030'
    date_str22 = '202107222110'
    date_str23 = '202107232155'
    dt_date_str20 = datetime.strptime(date_str20,"%Y%m%d%H%M")
    var1, crs1, lat_lims1, lon_lims1 = read_true_color('202107202125',\
        composite=composite)
    dt_date_str21 = datetime.strptime(date_str21,"%Y%m%d%H%M")
    var2, crs2, lat_lims2, lon_lims2 = read_true_color('202107212030',\
        composite=composite)
    dt_date_str22 = datetime.strptime(date_str22,"%Y%m%d%H%M")
    dt_date_str23 = datetime.strptime(date_str23,"%Y%m%d%H%M")
    var4, crs4, lat_lims4, lon_lims4 = read_true_color('202107232155',\
        composite=composite)

    fig = plt.figure(figsize=(9,12))
    ax0 = fig.add_subplot(3,2,1,projection = crs1) # true color 7/20
    ax1 = fig.add_subplot(3,2,3,projection = crs2) # true color 7/21
    ax3 = fig.add_subplot(3,2,5,projection = crs4) # true color 7/23
    ax4 = fig.add_subplot(3,2,2) # meteo 7/20
    ax5 = fig.add_subplot(3,2,4)   # meteo 7/21   
    ax7 = fig.add_subplot(3,2,6) # meteo 7/23

    # Plot the true-color data for the case date
    # ----------------------------------------------
    ax0.imshow(var1.data, transform = crs1, extent=(var1.x[0], var1.x[-1], \
        var1.y[-1], var1.y[0]), origin='upper')
    ax1.imshow(var2.data, transform = crs2, extent=(var2.x[0], var2.x[-1], \
        var2.y[-1], var2.y[0]), origin='upper')
    ax3.imshow(var4.data, transform = crs4, extent=(var4.x[0], var4.x[-1], \
        var4.y[-1], var4.y[0]), origin='upper')

    ax0.set_extent([lon_lims1[0],lon_lims1[1],lat_lims1[0],\
        lat_lims1[1]],crs = datacrs)
    ax1.set_extent([lon_lims2[0],lon_lims2[1],lat_lims2[0],\
        lat_lims2[1]],crs = datacrs)
    ax3.set_extent([lon_lims4[0],lon_lims4[1],lat_lims4[0],\
        lat_lims4[1]],crs = datacrs)

    # Plot the time series of ASOS data from O05 and AAT
    # --------------------------------------------------
    plot_asos_diurnal(ax4, date_str20, 'O05', 'AAT')
    plot_asos_diurnal(ax5, date_str21, 'O05', 'AAT')
    plot_asos_diurnal(ax7, date_str23, 'O05', 'AAT')

    # Plot the locations of ASOS sites O05 and AAT
    # --------------------------------------------------
    plot_ASOS_locs(ax0,date_str20,color='lime', sites = ['O05','AAT'])
    plot_ASOS_locs(ax1,date_str21,color='lime', sites = ['O05','AAT'])
    plot_ASOS_locs(ax3,date_str23,color='lime', sites = ['O05','AAT'])

    dt_date_str = dt_date_str22
    xval = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lon'][0] + \
        (plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lon'][1] - \
        plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lon'][0])*0.05
    yval = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lat'][0] + \
        (plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lat'][1] - \
        plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lat'][0])*0.90
    plot_subplot_label(ax0, '(a)',color = 'white', xval = xval, yval = yval, \
        transform = datacrs)
    plot_subplot_label(ax1, '(b)',color = 'white', xval = xval, yval = yval, \
        transform = datacrs)
    plot_subplot_label(ax3, '(c)',color = 'white', xval = xval, yval = yval, \
        transform = datacrs)

    plot_subplot_label(ax4, '(d)', xval = 1200, location = 'lower_left')
    plot_subplot_label(ax5, '(e)', xval = 1200, location = 'lower_left')
    plot_subplot_label(ax7, '(f)', xval = 1200, location = 'lower_left')

    def plot_goes_line(dt_date, pax):
        local_goes_time = dt_date - timedelta(hours = 7)
        goes_diff = (local_goes_time - datetime(year=2021,month=7,\
            day=22,hour=0,minute=0)).seconds
        print(local_goes_time, goes_diff)
        pax.axvline(goes_diff,color='black',linestyle = '--', lw=2,alpha=0.75,\
            label='GOES')

    plot_goes_line(dt_date_str20, ax4)
    plot_goes_line(dt_date_str21, ax5)
    plot_goes_line(dt_date_str23, ax7)
 
    fig.tight_layout()
 
    if(save):
        outname = 'goes_asos_meteo_combined_20210722_figureS1.png'
        fig.savefig(outname, dpi=300)
        print("Saved image",outname)
    else: 
        plt.show() 

def plot_combined_figureS2(zoom = True, show_smoke = False, composite = True, \
        plume_only = False, avg_pixel = True, save=False):

    date_str = '202108062025'
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # ----------------------------------------------------------------------
    #
    # Read the GOES data
    #
    # ----------------------------------------------------------------------

    # Call read_GOES_channel to read the desired GOES data from the file 
    # and put it in a dictionary
    # ---------------------------------------------------------------------
    GOES_data_ch1  = read_GOES_channel(date_str, 1, zoom = zoom)
    GOES_data_ch5  = read_GOES_channel(date_str, 5, zoom = zoom)
    GOES_data_ch31 = read_GOES_channel(date_str, 31, zoom = zoom)
    GOES_data_wv   = read_GOES_channel(date_str, 'wv_ir', zoom = zoom)

    # Determine where the smoke is located
    # ------------------------------------
    hash_data, nohash_data = find_plume(dt_date_str.strftime('%Y%m%d%H%M')) 

    tmp_data1  = np.copy(GOES_data_ch1['data'])
    tmp_data5  = np.copy(GOES_data_ch5['data'])
    tmp_data31 = np.copy(GOES_data_ch31['data'])
    tmp_datawv = np.copy(GOES_data_wv['data'])
    tmp_lat0   = np.copy(GOES_data_ch1['lat'])
    tmp_lon0   = np.copy(GOES_data_ch1['lon'])

    if(not (tmp_data1.shape == tmp_data5.shape == tmp_data31.shape == \
            tmp_datawv.shape == hash_data.shape)):
        print("shape mismatch")
        shapes = []
        shapes.append(tmp_data1.shape)
        shapes.append(tmp_data5.shape)
        shapes.append(tmp_data31.shape)
        shapes.append(tmp_datawv.shape)
        shapes.append(hash_data.shape)

        min_shape = min(shapes)
        print(min_shape)

        tmp_data1  = tmp_data1[:min_shape[0],:min_shape[1]]
        tmp_data5  = tmp_data5[:min_shape[0],:min_shape[1]]
        tmp_data31 = tmp_data31[:min_shape[0],:min_shape[1]]
        tmp_datawv = tmp_datawv[:min_shape[0],:min_shape[1]]
        tmp_lat0   = tmp_lat0[:min_shape[0],:min_shape[1]]
        tmp_lon0   = tmp_lon0[:min_shape[0],:min_shape[1]]
        hash_data  = hash_data[:min_shape[0],:min_shape[1]]

    max_ch = 350.

    tmp_data1 = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) | \
        (abs(tmp_datawv) > max_ch), tmp_data1)
    tmp_data5 = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) | \
        (abs(tmp_datawv) > max_ch), tmp_data5)
    tmp_data31 = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) | \
        (abs(tmp_datawv) > max_ch), tmp_data31)
    tmp_datawv = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) | \
        (abs(tmp_datawv) > max_ch), tmp_datawv)
    tmp_lat0 = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) | \
        (abs(tmp_datawv) > max_ch), tmp_lat0)
    tmp_lon0 = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) | \
        (abs(tmp_datawv) > max_ch), tmp_lon0)


    # ----------------------------------------------------------------------
    #
    #  Set up the 4-panel figure
    #
    # ----------------------------------------------------------------------
 
    mapcrs = init_proj(date_str)
    plt.close('all')
    fig = plt.figure(figsize=(11,9))
    ax1 = fig.add_subplot(2,2,1, projection = mapcrs) # Ch 1
    ax2 = fig.add_subplot(2,2,2)                      # AI scatter
    ax3 = fig.add_subplot(2,2,3)                      # AI vs flux scatter
    ax4 = fig.add_subplot(2,2,4)                      # AI vs flux scatter

    # ----------------------------------------------------------------------
    #
    # Plot the data in the figure
    #
    # ----------------------------------------------------------------------

    # Plot channel 1 spatial data
    # -----------------------------------------------
    plot_GOES_spatial(GOES_data_ch1, ax1, zoom = zoom, ptitle = '')

    plot_scatter(ax2, tmp_data31, tmp_data1, GOES_data_ch31, GOES_data_ch1, \
        hash_data, xlabel = '11 μm brightness temperature', \
        ylabel = '0.64 μm reflectance', plot_legend = True)
    plot_scatter(ax3, tmp_data31, tmp_data5, GOES_data_ch31, GOES_data_ch5, \
        hash_data, xlabel = '11 μm brightness temperature', \
        ylabel = '1.24 μm reflectance')
    plot_scatter(ax4, tmp_data31, tmp_datawv, GOES_data_ch31, GOES_data_wv, \
        hash_data, xlabel = '11 μm brightness temperature', \
        ylabel = 'GOES IR Precipitable Water')

    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(date_str) 

        ax3.pcolor(GOES_data_ch1['lon'],GOES_data_ch1['lat'],\
            hash_data1, hatch = '////', alpha=0., transform = datacrs,\
            color = 'lime')

    # Add subplot labels
    # ------------------
    plot_subplot_label(ax1, '(a)', backgroundcolor = 'white')
    plot_subplot_label(ax2, '(b)')
    plot_subplot_label(ax3, '(c)')
    plot_subplot_label(ax4, '(d)')

    # Add plot text
    # -------------
    plot_figure_text(ax1, 'GOES 0.64 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    
    fig.tight_layout()

    GOES_data_ch1.clear()
    GOES_data_ch5.clear()
    GOES_data_ch31.clear()
    GOES_data_wv.clear()

    if(save):
        outname = 'goes_total_combined_' + date_str + '_figS2.png'
        fig.savefig(outname, dpi=300)
        print("Saved",outname)
    else:
        plt.show()


def plot_meteogram_compare(date_str, true_color = True, zoom=True, \
        save=False, composite=True, show_smoke = False, \
        plot_ASOS_loc = True, ir_data = True):
    
    plt.close('all')
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    mapcrs = init_proj(date_str)

    # -----------------------------------------------------------------------     
    #
    # Figure format:
    #   
    #       1. True color / ch1 clear day      2. True color / ch1 smoke day
    #
    #                       3. Meteogram
    #                       4. Meteogram
    #       
    # -----------------------------------------------------------------------     
    if(ir_data):
        fig = plt.figure(figsize=(12,15))
        gs  = plt.GridSpec(4, 3, wspace=0.4, hspace=0.05, figure = fig)
        ax3 = plt.subplot(gs[2,:])
        ax4 = plt.subplot(gs[3,:])
    else:
        fig = plt.figure(figsize=(12,12))
        gs  = plt.GridSpec(3, 3, wspace=0.4, hspace=0.05, figure = fig)
        ax3 = plt.subplot(gs[1,:])
        ax4 = plt.subplot(gs[2,:])

    # -----------------------------------------------------------------------     
    #
    # Axis 0: True color image of previous day
    #
    # -----------------------------------------------------------------------     

    if(date_str == '202107222110'):
        prev_str = '202107212030'
        post_str = '202107232155'
    elif(date_str == '202108062025'):
        prev_str = '202108052120'
        post_str = '202108072110'

    dt_prev_str = datetime.strptime(prev_str, '%Y%m%d%H%M')
    dt_post_str = datetime.strptime(post_str, '%Y%m%d%H%M')

    if(true_color):
        # Read true color data for both dates
        # -----------------------------------
        var1, crs1, lat_lims1, lon_lims1 = \
            read_true_color(prev_str,composite=composite)
        ax0 = plt.subplot(gs[0,0],projection=crs1)

        # Plot the true-color data for the previous date
        # ----------------------------------------------
        ax0.imshow(var1.data, transform = crs1, extent=(var1.x[0], var1.x[-1], \
            var1.y[-1], var1.y[0]), origin='upper')

        # Zoom in the figure if desired
        # -----------------------------
        if(zoom):
            ax0.set_extent([lon_lims1[0],lon_lims1[1],lat_lims1[0],\
                lat_lims1[1]],crs = datacrs)
    else:
        # Read GOES Ch 31 data
        # ---------------------
        GOES_data1 = read_GOES_channel(prev_str, 1, zoom = zoom)

        # Plot the GOES ch 31 data
        # -------------------------
        mapcrs = init_proj(date_str)
        ax0 = plt.subplot(gs[0,0], projection = mapcrs)
        plot_GOES_spatial(GOES_data1, ax0, zoom = zoom)

        if(ir_data):
            # Read GOES Ch 31 data
            # ---------------------
            GOES_data4 = read_GOES_channel(prev_str, 31, zoom = zoom)

            # Plot the GOES ch 31 data
            # -------------------------
            mapcrs = init_proj(date_str)
            ax5 = plt.subplot(gs[1,0], projection = mapcrs)
            plot_GOES_spatial(GOES_data4, ax5, zoom = zoom)


    ax0.set_title('Aqua GOES\n' + dt_prev_str.strftime('%Y-%m-%d %H:%M'))

    # -----------------------------------------------------------------------     
    #
    # Axis 1: True color image of smokey day
    #
    # -----------------------------------------------------------------------     
    if(true_color):
        var, crs, lat_lims, lon_lims     = \
            read_true_color(date_str,composite=composite)
        ax1 = plt.subplot(gs[0,1],projection=crs)

        # Plot the true-color data for the case date
        # ----------------------------------------------
        ax1.imshow(var.data, transform = crs, extent=(var.x[0], var.x[-1], \
            var.y[-1], var.y[0]), origin='upper')

        # Zoom in the figure if desired
        # -----------------------------
        if(zoom):
            ax1.set_extent([lon_lims[0],lon_lims[1],lat_lims[0],\
                lat_lims[1]],crs = datacrs)

    else:
        # Read GOES Ch 31 data
        # ---------------------
        GOES_data2 = read_GOES_channel(date_str, 1, zoom = zoom)

        # Plot the GOES ch 31 data
        # -------------------------
        mapcrs = init_proj(date_str)
        ax1 = plt.subplot(gs[0,1], projection = mapcrs)
        plot_GOES_spatial(GOES_data2, ax1, zoom = zoom)

        if(ir_data):
            # Read GOES Ch 31 data
            # ---------------------
            GOES_data5 = read_GOES_channel(date_str, 31, zoom = zoom)

            # Plot the GOES ch 31 data
            # -------------------------
            mapcrs = init_proj(date_str)
            ax6 = plt.subplot(gs[1,1], projection = mapcrs)
            plot_GOES_spatial(GOES_data5, ax6, zoom = zoom)

    ax1.set_title('Aqua GOES\n' + dt_date_str.strftime('%Y-%m-%d %H:%M'))

    # -----------------------------------------------------------------------     
    #
    # Axis 2: True color image of day after 
    #
    # -----------------------------------------------------------------------     
    if(true_color):
        var2, crs2, lat_lims2, lon_lims2     = \
            read_true_color(post_str,composite=composite)
        ax2 = plt.subplot(gs[0,2],projection=crs2)

        # Plot the true-color data for the case date
        # ----------------------------------------------
        ax2.imshow(var2.data, transform = crs2, extent=(var2.x[0], var2.x[-1], \
            var2.y[-1], var2.y[0]), origin='upper')

        # Zoom in the figure if desired
        # -----------------------------
        if(zoom):
            ax2.set_extent([lon_lims2[0],lon_lims2[1],lat_lims2[0],\
                lat_lims2[1]],crs = datacrs)

    else:
        # Read GOES Ch 31 data
        # ---------------------
        GOES_data3 = read_GOES_channel(post_str, 1, zoom = zoom)

        # Plot the GOES ch 31 data
        # -------------------------
        mapcrs = init_proj(post_str)
        ax2 = plt.subplot(gs[0,2], projection = mapcrs)
        plot_GOES_spatial(GOES_data3, ax2, zoom = zoom)

        if(ir_data):
            # Read GOES Ch 31 data
            # ---------------------
            GOES_data6 = read_GOES_channel(post_str, 31, zoom = zoom)

            # Plot the GOES ch 31 data
            # -------------------------
            mapcrs = init_proj(date_str)
            ax7 = plt.subplot(gs[1,2], projection = mapcrs)
            plot_GOES_spatial(GOES_data6, ax7, zoom = zoom)

    if(plot_ASOS_loc):

        # Plot the ASOS locations on the map
        # ----------------------------------
        plot_ASOS_locs(ax0,prev_str,color='lime')
        plot_ASOS_locs(ax1,date_str,color='lime')
        plot_ASOS_locs(ax2,post_str,color='lime')

        if(ir_data):
            plot_ASOS_locs(ax5,prev_str,color='lime')
            plot_ASOS_locs(ax6,date_str,color='lime')
            plot_ASOS_locs(ax7,post_str,color='lime')


    ax2.set_title('Aqua GOES\n' + dt_post_str.strftime('%Y-%m-%d %H:%M'))

    # -----------------------------------------------------------------------     
    #
    # Axis 3: Meteogram of ASOS data     
    # Axis 4: Difference in temperatures between smokey station and clear
    # stations
    #
    # -----------------------------------------------------------------------     
    colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple',\
        'tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']

    # Open the ASOS data file
    # -----------------------
    df = pd.read_csv(plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')][\
        dt_date_str.strftime('%H%M')]['asos'])

    # Read ASOS data for the current date
    # -----------------------------------
    station_names = sorted(set(df['station'].values))
    for ii, station in enumerate(station_names):
        time_stn = df['valid'][df['station'] == station]
        tmpc_stn = pd.to_numeric(df['tmpc'][df['station'] == station], errors='coerce').values
        drct_stn = pd.to_numeric(df['drct'][df['station'] == station], errors='coerce').values
        pwxc_stn = df['wxcodes'][df['station'] == station].values
        lat_stn = df['lat'][df['station'] == station].values[0]
        lon_stn = df['lon'][df['station'] == station].values[0]

        # Remove masked data
        # ------------------
        indices = np.where(~np.isnan(tmpc_stn))
        time_stn = time_stn.values[indices] 
        tmpc_stn = np.ma.masked_invalid(tmpc_stn).compressed()
        drct_stn = np.ma.masked_invalid(drct_stn).compressed()
        pwxc_stn = pwxc_stn[indices]

        time_stn = np.array([datetime.strptime(ttstn,"%Y-%m-%d %H:%M") \
            for ttstn in time_stn])

        dtime_stn = [ttime - timedelta(hours = 5) for ttime in time_stn]

        if(station != case_dict[dt_date_str.strftime('%Y%m%d%H%M')]):
            time_ref = df['valid'][df['station'] == \
                case_dict[dt_date_str.strftime('%Y%m%d%H%M')]]
            tmpc_ref = pd.to_numeric(df['tmpc'][df['station'] == \
                case_dict[dt_date_str.strftime('%Y%m%d%H%M')]], errors='coerce').values

            indices_ref = np.where(~np.isnan(tmpc_ref))
            time_ref = time_ref.values[indices_ref] 
            tmpc_ref = np.ma.masked_invalid(tmpc_ref).compressed()
        
            time_ref = np.array([datetime.strptime(ttref,"%Y-%m-%d %H:%M") \
                for ttref in time_ref])

            if(len(tmpc_ref) > len(tmpc_stn)):
                # Reference data are larger than station data. Make an array
                # of size of station data, loop over station times, find where
                # the reference time matches the current station time, insert
                # the reference data into the copy array.
                # ----------------------------------------------------------

                work_time = np.copy(time_stn)
                match_tmpc_ref = np.full((tmpc_stn.size), np.nan)
   
                # Loop over the station times and find the station time
                # that matches the current reference time
                for jj in range(len(time_stn)):
                    time_offsets = np.array([(ttref - \
                        time_stn[jj]).total_seconds()/60. \
                        for ttref in time_ref])
                    close_vals = np.where(abs(time_offsets) < 10.)[0]
                    if(len(close_vals) > 0):
                        locator = np.argmin(abs(time_offsets[close_vals])) 
                        match_tmpc_ref[jj] = tmpc_ref[close_vals[locator]]

                tmpc_dif = match_tmpc_ref - tmpc_stn

            else:
                # Station data are larger than reference data. Make an array
                # of size of reference data, loop over reference times, find 
                # where the station time matches the current reference time, 
                # insert the station data into the copy array.
                # ----------------------------------------------------------
                work_time = np.copy(time_ref)
                match_tmpc_stn = np.full((tmpc_ref.size), np.nan)
    
                # Loop over the reference times and find the station time
                # that matches the current reference time
                for jj in range(len(time_ref)):
                    time_offsets = np.array([(ttref - \
                        time_stn[jj]).total_seconds()/60. \
                        for ttref in time_ref])
                    close_vals = np.where(abs(time_offsets) < 10.)[0]
                    if(len(close_vals) > 0):
                        locator = np.argmin(abs(time_offsets[close_vals])) 
                        match_tmpc_stn[jj] = tmpc_stn[close_vals[locator]]

                ## Loop over the reference times and grab  
                tmpc_dif = tmpc_ref - match_tmpc_stn 
 
            ax4.plot(dtime_stn,tmpc_dif, color=colors[ii])
            #ax3.plot(dtime_stn,tmpc_dif,label=station, color=colors[ii])
    
        # Plot the temp data differently for when it reports smoke/haze
        # -------------------------------------------------------------
        tmpc_stn_hz = np.copy(tmpc_stn)
        tmpc_stn_nohz = np.copy(tmpc_stn)
        tmpc_stn_hz   = np.ma.masked_where((pwxc_stn != 'HZ') & \
            (pwxc_stn != 'FU'), tmpc_stn_hz)
        tmpc_stn_nohz = np.ma.masked_where((pwxc_stn == 'HZ') | \
            (pwxc_stn == 'FU'), tmpc_stn_nohz)
    
        ax3.plot(dtime_stn,tmpc_stn,label=station, color=colors[ii])
        #ax.plot(dtime_stn,tmpc_stn_nohz,label=station, color=colors[ii])
        #ax.plot(dtime_stn,tmpc_stn_hz,'--', label=station, color=colors[ii])
        print(station, lat_stn, lon_stn)
        #ax2.text(lon_stn, lat_stn, station, transform=datacrs, color=colors[ii])
    
        #ax3.plot(dtime_stn,drct_stn_nohz,label=station, color=colors[ii])
        #ax3.plot(dtime_stn,drct_stn_hz,'--', label=station, color=colors[ii])
    
    # Convert the file time to a plot_limits_dict format
    #event_date = datetime.strptime(infile.split('/')[-1].split('_')[-1][:8], "%Y%m%d")
    #grabber_date = event_date.strftime('%Y-%m-%d')
    #first_time = list(plot_limits_dict[grabber_date].keys())[0]
    
    # Pull the event time from the plot_limits_dict
    #event_dtime = event_date + timedelta(hours = int(first_time[:2]), \
    #    minutes = int(first_time[2:4]))

    # Plot a vertical line at the time of the clear GOES overpass
    # ------------------------------------------------------------
    ax3.axvline(dt_prev_str,color='black',linestyle = '--', lw=2,alpha=0.75,\
        label='' + dt_prev_str.strftime('%m/%d %H:%M'))
    ax4.axvline(dt_prev_str,color='black',linestyle = '--',lw=2,alpha=0.75)
    
    # Plot a vertical line at the time of the GOES overpass
    # ------------------------------------------------------
    ax3.axvline(dt_date_str,color='black',lw=2,alpha=0.75,label='' \
        + dt_date_str.strftime('%m/%d %H:%M'))
    ax4.axvline(dt_date_str,color='black',lw=2,alpha=0.75)

    # Plot a vertical line at the time of the post GOES overpass
    # -----------------------------------------------------------
    ax3.axvline(dt_post_str,color='black',linestyle = ':', lw=2,alpha=0.75,label='' \
        + dt_post_str.strftime('%m/%d %H:%M'))
    ax4.axvline(dt_post_str,color='black',linestyle = ':', lw=2,alpha=0.75)

    # Plot a horizontal line at 0
    # ---------------------------
    ax4.axhline(0,color='grey',linestyle='--',lw=2,alpha=0.75)
    #ax3.axvline(event_dtime,color='tab:red',lw=2,alpha=0.75)

    ax4.sharex(ax3)
    #ax1.sharey(ax3) 

 
    #ax3.set_xlabel('Time [UTC]')
    ax3.set_ylabel('2-m Temperature [$^{o}$C]')
    ax4.set_ylabel('T(' + case_dict[dt_date_str.strftime('%Y%m%d%H%M')] + \
        ') - T(other stations)')
    #ax2.legend()
    ax3.tick_params(axis="x", \
        labelbottom = False)

    # -----------------------------------------------------------------------     
    #
    # Extract GOES values at each pixel and plot on the meteograms
    #
    # -----------------------------------------------------------------------     
   

 
    # Add subplot letter labels
    # -------------------------
    xval = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lon'][0] + \
        (plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lon'][1] - \
        plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lon'][0])*0.05
    yval = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lat'][0] + \
        (plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lat'][1] - \
        plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lat'][0])*0.90
    plot_subplot_label(ax0, '(a)', xval = xval, yval = yval, \
        transform = datacrs)
    plot_subplot_label(ax1, '(b)', xval = xval, yval = yval, \
        transform = datacrs)
    plot_subplot_label(ax2, '(c)', xval = xval, yval = yval, \
        transform = datacrs)
    plot_subplot_label(ax3, '(d)', xval = dtime_stn[0])
    plot_subplot_label(ax4, '(e)', xval = dtime_stn[0])

    # Add a legend
    # ------------
    lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    plt.legend(lines, labels, loc = 'lower center', bbox_to_anchor = (0, 0.05, 1, 1),\
        bbox_transform = plt.gcf().transFigure, ncol=5)
    #ax3.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
    #          fancybox=True, shadow=True, ncol=5)
 
    fig.tight_layout()

    if(save):
        outname = 'goes_combined_station_' + date_str + '.png'
        fig.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

# Calculate the average diurnal cycle for each day in the provided ASOS file
# --------------------------------------------------------------------------
def calc_meteogram_climo(training_asos_file, work_stn, save=False):
    
    # Read in the ASOS data using pandas
    # ----------------------------------
    df = pd.read_csv(training_asos_file)

    # Modify the dataframe to use date as index
    # -----------------------------------------
    # Convert from UTC to PDT
    df['valid'] = pd.to_datetime(df['valid'], format = '%Y-%m-%d %H:%M') - \
        timedelta(hours = 7)
    df['tmpc'] = pd.to_numeric(df['tmpc'], errors = 'coerce').values
    time_df = df.set_index('valid')

    start_date = datetime(year = 2021, month = 7, day = 1, hour = 0, \
        minute = 0)

    # Determine all the unique local times
    # ------------------------------------
    unique_dtimes = np.array([start_date + timedelta(minutes = 15) + \
        timedelta(minutes = 20 * ii) for ii in range(72)])
    unique_times = np.array([dtime.strftime('%H:%M:%S') for dtime in \
        unique_dtimes])

    # Create arrays to hold the average and standard deviation of the
    # temperatures at each 20-minute interval
    # ---------------------------------------------------------------
    avg_curve_test = np.full(unique_times.shape, np.nan)
    std_curve_test = np.full(unique_times.shape, np.nan)
  
    # Convert the time index to a datetime string format for comparisons
    # ------------------------------------------------------------------ 
    time_idx = time_df[time_df['station'] == work_stn].index.strftime("%H:%M:%S")

    # Loop over each of the climatology times
    # ---------------------------------------
    for ii in range(len(unique_dtimes)):
        # Find the obs that are +/- 10 minutes around this ob
        # ---------------------------------------------------
        prev_time = (unique_dtimes[ii] - timedelta(minutes = 10)).strftime("%H:%M:%S")
        post_time = (unique_dtimes[ii] + timedelta(minutes = 10)).strftime("%H:%M:%S")
        locate_tmpc = time_df[time_df['station'] == work_stn]['tmpc'].values[\
            np.where((time_idx >= prev_time) & (time_idx < post_time))]

        # Insert the average and standard deviation of the selected obs
        # within this time range for all diurnal cycles into the arrays
        # -------------------------------------------------------------
        if(len(locate_tmpc) > 0):
            avg_curve_test[ii] = np.nanmean(locate_tmpc)
            std_curve_test[ii] = np.nanstd(locate_tmpc)
        """
        avg_curve_test[ii] = np.nanmean(\
            pd.to_numeric(time_df['tmpc'][time_df['station'] == work_stn], \
            errors = 'coerce').values[np.where(time_df[\
            time_df['station'] == work_stn].index.strftime('%H:%M:%S') == \
            unique_times[ii])])
        std_curve_test[ii] = np.nanstd(\
            pd.to_numeric(time_df['tmpc'][time_df['station'] == work_stn], \
            errors = 'coerce').values[np.where(time_df[\
            time_df['station'] == work_stn].index.strftime('%H:%M:%S') == \
            unique_times[ii])])
        #std_curve_test[ii] = np.nanstd(pd.to_numeric(time_df['tmpc'][time_df['station'] == work_stn], \
        #    errors = 'coerce').values[np.where(time_df[time_df['station'] == work_stn].index.strftime('%H:%M:%S') == \
        #    unique_times[ii])])
        """
  
    # Return the local times (in datetime format), the averages,
    # and the standard deviations to the main function
    # ----------------------------------------------------------

    return unique_dtimes, avg_curve_test, std_curve_test 
    #plus_vals = avg_curve_test + std_curve_test
    #minus_vals = avg_curve_test - std_curve_test

    ##fig1 = plt.figure()
    ##ax2 = fig1.add_subplot(1,1,1)
    ##print(avg_curve_test)
    ##ax2.plot(np.ma.masked_invalid(avg_curve_test).compressed())
    ##ax2.plot(np.ma.masked_invalid(plus_vals).compressed())
    ##ax2.plot(np.ma.masked_invalid(minus_vals).compressed())

    """
    # For now, just plot each of the diurnal cycles
    # ---------------------------------------------
    fig2 = plt.figure()
    ax = fig2.add_subplot(1,1,1)

    time_df = time_df[time_df['station'] == work_stn]
    while(start_date.day < 14):
        subset_df = time_df[(time_df.index > start_date) & \
            (time_df.index < start_date + timedelta(days = 1))]
        
        tmpc_stn = pd.to_numeric(subset_df['tmpc'], errors='coerce').values
        local_times = subset_df.index - start_date
        #print(local_times)
        tmpc_stn = np.ma.masked_invalid(tmpc_stn)
    
        ax.plot(local_times, tmpc_stn,label=str(start_date.day))
        start_date = start_date + timedelta(days = 1)

    ax.legend()
    plt.show() 
    """

# Plot a comparison of the diurnal cycle from the inputted day
# with the "climatological" average from the training period
# ------------------------------------------------------------
def plot_asos_diurnal(pax, date_str, work_stn, base_stn, plot_map = False, \
        save = False):
    
    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
    case_str = date_str[:8]

    # Read the ASOS data for the desired date
    # ---------------------------------------
    work_data = 'asos_data_20210722_4.csv'
    df = pd.read_csv(work_data)

    # Convert from UTC to PDT
    df['dt_valid'] = pd.to_datetime(df['valid'], format = '%Y-%m-%d %H:%M') - \
        timedelta(hours = 7)
    df['tmpc'] = pd.to_numeric(df['tmpc'], errors = 'coerce').values
    time_df = df.set_index('dt_valid')

    # Extract only the data for the desired station at the desired
    # date
    # ------------------------------------------------------------
    base_df = time_df[time_df['station'] == base_stn]
    time_df = time_df[time_df['station'] == work_stn]
    #time_df = df

    # Calculate the climatological cycle for the desired station
    # -----------------------------------------------------------
    base_dtimes, base_climo, base_std = \
        calc_meteogram_climo('asos_data_20210722_5.csv', base_stn)
    plt_dtimes, stn_climo, stn_std = \
        calc_meteogram_climo('asos_data_20210722_5.csv', work_stn)

    climo_xval = np.array([(plt_dtime - plt_dtimes[0]).seconds for \
             plt_dtime in plt_dtimes])

    base_climo_xval = np.array([(base_dtime - base_dtimes[0]).seconds for \
             base_dtime in base_dtimes])
    base_mask_climo       = np.ma.masked_invalid(base_climo) 
    base_mask_climo_times = base_climo_xval[~base_mask_climo.mask]
    base_mask_climo_std   = base_std[~base_mask_climo.mask]
    base_mask_climo       = base_mask_climo[~base_mask_climo.mask]

    mask_climo       = np.ma.masked_invalid(stn_climo) 
    mask_climo_times = climo_xval[~mask_climo.mask]
    mask_climo_std   = stn_std[~mask_climo.mask]
    mask_climo       = mask_climo[~mask_climo.mask]

    # Calculate the differences between the climatological profiles
    # -------------------------------------------------------------
    diff_climo = np.ma.masked_invalid(stn_climo) - np.ma.masked_invalid(base_climo)
    diff_climo       = np.ma.masked_invalid(diff_climo) 
    diff_climo_times = climo_xval[~diff_climo.mask]
    diff_climo       = diff_climo[~diff_climo.mask]
    
    # Temporally colocate the two data
    # --------------------------------
    start_date = datetime(year = 2021, month = 7, day = 1, hour = 0, \
        minute = 0)
    unique_dtimes = np.array([start_date + timedelta(minutes = 15) + \
        timedelta(minutes = 20 * ii) for ii in range(72)])
    unique_times = np.array([dtime.strftime('%H:%M:%S') for dtime in \
        unique_dtimes])

    base_avg = np.full(unique_times.shape, np.nan)
    time_avg = np.full(unique_times.shape, np.nan)

    # Convert the time index to a datetime string format for comparisons
    # ------------------------------------------------------------------ 
    work_uday = datetime.strptime(case_str, '%Y%m%d')
    #work_uday = datetime.strptime(date_str, '%Y/%m/%d')

    # Now pull out only the data for the desired day
    next_day = work_uday + timedelta(days = 1)

    base_df = base_df[work_uday : next_day]
    time_df = time_df[work_uday : next_day]

    base_time_idx = base_df.index.strftime("%H:%M:%S")
    time_time_idx = time_df.index.strftime("%H:%M:%S")

    for ii in range(len(unique_dtimes)):
        # Find the obs that are +/- 10 minutes around this ob
        # ---------------------------------------------------
        prev_time = (unique_dtimes[ii] - timedelta(minutes = 10)).strftime("%H:%M:%S")
        post_time = (unique_dtimes[ii] + timedelta(minutes = 10)).strftime("%H:%M:%S")
        locate_base = base_df['tmpc'].values[\
            np.where((base_time_idx >= prev_time) & (base_time_idx < post_time))]
        locate_time = time_df['tmpc'].values[\
            np.where((time_time_idx >= prev_time) & (time_time_idx < post_time))]


        # Insert the average and standard deviation of the selected obs
        # within this time range for all diurnal cycles into the arrays
        # -------------------------------------------------------------
        if(len(locate_base) > 0):
            base_avg[ii] = np.nanmean(locate_base)
        if(len(locate_time) > 0):
            time_avg[ii] = np.nanmean(locate_time)
    

    base_avg_xval = np.array([(unique_dtime - unique_dtimes[0]).seconds for \
             unique_dtime in unique_dtimes])
    base_avg_xval = np.array([(unique_dtime - unique_dtimes[0]).seconds for \
             unique_dtime in unique_dtimes])

    diff_local = np.ma.masked_invalid(time_avg) - np.ma.masked_invalid(base_avg)
    diff_local       = np.ma.masked_invalid(diff_local) 
    diff_local_times = base_avg_xval[~diff_local.mask]
    diff_local       = diff_local[~diff_local.mask]
    #base_avg_climo       = np.ma.masked_invalid(base_climo) 
    #base_avg_climo_times = base_climo_xval[~base_mask_climo.mask]

    mask_base_avg = np.ma.masked_invalid(base_avg)
    mask_time_avg = np.ma.masked_invalid(time_avg)

    mask_time_times = base_avg_xval[~mask_time_avg.mask]
    mask_time_data = time_avg[~mask_time_avg.mask]
    mask_base_times = base_avg_xval[~mask_base_avg.mask]
    mask_base_data  = base_avg[~mask_base_avg.mask]

    ##!#work_uday = datetime.strptime('2021/07/22', '%Y/%m/%d')

    ##!## Now pull out only the data for the desired day
    ##!#next_day = work_uday + timedelta(days = 1)

    ##!#base_time_df = base_df[work_uday : next_day]
    ##!#time_time_df = time_df[work_uday : next_day]

    ##!#base_local_times = np.array([datetime.strptime(tst, '%Y-%m-%d %H:%M') for \
    ##!#    tst in base_time_df['valid'].values])
    ##!#time_local_times = np.array([datetime.strptime(tst, '%Y-%m-%d %H:%M') for \
    ##!#    tst in time_time_df['valid'].values])

    ##!#stn_xval = np.array([(base_local_time - base_local_times[0]).seconds \
    ##!#    for base_local_time in base_local_times])



    ##!#base_mask_tmp = np.ma.masked_invalid(base_time_df['tmpc']) 
    ##!#time_mask_tmp = np.ma.masked_invalid(time_time_df['tmpc']) 

    ##!##diff_day = time_mask_tmp - base_mask_tmp

    ##!##print(base_mask_tmp.compressed().shape, time_mask_tmp.compressed().shape)

    ##!#mask_times = stn_xval[~mask_tmp.mask]
    ##!#mask_tmp   = mask_tmp[~mask_tmp.mask]



    # Plot the climatological data and the single day on a meteogram
    # ---------------------------------------------------------------
    test_times = np.array([start_date + timedelta(seconds = float(val)) for val in diff_local_times])
    test_strtimes = np.array([ttime.strftime('%H:%M') for ttime in test_times])

    pax.plot(base_mask_climo_times, base_mask_climo, color = 'tab:blue', label = 'μ$_{ ' + base_stn + '}$')
    pax.plot(mask_climo_times, mask_climo, color = 'tab:orange', label = 'μ$_{ ' + work_stn + '}$')
    pax.plot(mask_base_times, mask_base_data, color = 'tab:blue', linestyle = '--', label = 'T$_{ ' + base_stn + '}$')
    pax.plot(mask_time_times, mask_time_data, color = 'tab:orange', linestyle = '--', label = 'T$_{ ' + work_stn + '}$')
    pax.set_xticks(diff_local_times[::6])
    pax.set_xticklabels(test_strtimes[::6], fontsize = 12)
    pax.set_ylabel('2-m Temperature [$^{o}$C]', fontsize = 14, weight = 'bold')
    pax.tick_params(axis='y', labelsize = 12)
    pax.set_xlabel('Local time [PDT]', fontsize = 14, weight = 'bold')
    pax.set_title(dt_date_str.strftime('%Y-%m-%d'), fontsize = 14)
    pax.grid()
    pax.legend()
    #ax.plot(mask_climo_times, mask_climo, color = 'r', label = 'μ$_{T ' + work_stn + '}$')
    ##!#ax2.plot(diff_climo_times, diff_climo, color = 'tab:blue', label = 'μ$_{ ' + work_stn + '}$ - μ$_{ ' + base_stn + '}$')
    ##!#ax2.plot(diff_local_times, diff_local, color = 'tab:blue', linestyle = '--', label = 'T$_{' + work_stn + '}$ - T$_{' + base_stn+'}$')
    ##!#ax2.set_xticks(diff_local_times[::6])
    ##!#ax2.set_xticklabels(test_strtimes[::6])
    ##!#ax2.set_ylabel('2-m Temperature Difference [$^{o}$C]')
    ##!#ax2.set_xlabel('Local time [PDT]')
    ##!#ax2.set_title(dt_date_str.strftime('%Y-%m-%d'))
    ##!#ax2.grid()
    ##!#ax2.legend()


    """
    # Plot the diurnal cycle from each day in the ASOS file for the desired
    # station
    # ---------------------------------------------------------------------

    # Extract all the unique days
    unique_days = sorted(set([datetime.strptime(tdt,'%Y-%m-%d %H:%M'\
        ).strftime('%Y/%m/%d') for tdt in time_df['valid'].values]))

    for uday in unique_days:
        
        work_uday = datetime.strptime(uday, '%Y/%m/%d')

        # Now pull out only the data for the desired day
        next_day = work_uday + timedelta(days = 1)
        #stn_df = time_df[((time_df['dt_valid'] >= dt_date_str) & (time_df['dt_valid'] <= next_day))]
        stn_df = time_df[work_uday : next_day]

        
        local_times = np.array([datetime.strptime(tst, '%Y-%m-%d %H:%M') for \
            tst in stn_df['valid'].values]) - timedelta(hours = 8)

        stn_xval = np.array([(local_time - local_times[0]).seconds for local_time in local_times])

        mask_tmp = np.ma.masked_invalid(stn_df['tmpc']) 

        mask_times = stn_xval[~mask_tmp.mask]
        mask_tmp   = mask_tmp[~mask_tmp.mask]

        ax.plot(mask_times, mask_tmp, label = work_stn + ' ' + uday[5:])
    ax.set_ylabel('2-m temperature [$^{o}$C]')
    ax.set_title(work_stn)
    plt.legend()
    plt.show()

    """

def plot_total_asos_diurnal(save = False, composite = True):

    # Read true color data for the previous date
    date_str20 = '202107202125'
    date_str21 = '202107212030'
    date_str22 = '202107222110'
    date_str23 = '202107232155'
    dt_date_str20 = datetime.strptime(date_str20,"%Y%m%d%H%M")
    var1, crs1, lat_lims1, lon_lims1 = read_true_color('202107202125',\
        composite=composite)
    dt_date_str21 = datetime.strptime(date_str21,"%Y%m%d%H%M")
    var2, crs2, lat_lims2, lon_lims2 = read_true_color('202107212030',\
        composite=composite)
    dt_date_str22 = datetime.strptime(date_str22,"%Y%m%d%H%M")
    var3, crs3, lat_lims3, lon_lims3 = read_true_color('202107222110',\
        composite=composite)
    dt_date_str23 = datetime.strptime(date_str23,"%Y%m%d%H%M")
    var4, crs4, lat_lims4, lon_lims4 = read_true_color('202107232155',\
        composite=composite)

    fig = plt.figure(figsize=(9,15))
    ax0 = fig.add_subplot(4,2,1,projection = crs1) # true color 7/20
    ax1 = fig.add_subplot(4,2,3,projection = crs2) # true color 7/21
    ax2 = fig.add_subplot(4,2,5,projection = crs3) # true color 7/22
    ax3 = fig.add_subplot(4,2,7,projection = crs4) # true color 7/23
    ax4 = fig.add_subplot(4,2,2) # meteo 7/20
    ax5 = fig.add_subplot(4,2,4)   # meteo 7/21   
    ax6 = fig.add_subplot(4,2,6) # meteo 7/22
    ax7 = fig.add_subplot(4,2,8) # meteo 7/23

    # Plot the true-color data for the case date
    # ----------------------------------------------
    ax0.imshow(var1.data, transform = crs1, extent=(var1.x[0], var1.x[-1], \
        var1.y[-1], var1.y[0]), origin='upper')
    ax1.imshow(var2.data, transform = crs2, extent=(var2.x[0], var2.x[-1], \
        var2.y[-1], var2.y[0]), origin='upper')
    ax2.imshow(var3.data, transform = crs3, extent=(var3.x[0], var3.x[-1], \
        var3.y[-1], var3.y[0]), origin='upper')
    ax3.imshow(var4.data, transform = crs4, extent=(var4.x[0], var4.x[-1], \
        var4.y[-1], var4.y[0]), origin='upper')

    ax0.set_extent([lon_lims1[0],lon_lims1[1],lat_lims1[0],\
        lat_lims1[1]],crs = datacrs)
    ax1.set_extent([lon_lims2[0],lon_lims2[1],lat_lims2[0],\
        lat_lims2[1]],crs = datacrs)
    ax2.set_extent([lon_lims3[0],lon_lims3[1],lat_lims3[0],\
        lat_lims3[1]],crs = datacrs)
    ax3.set_extent([lon_lims4[0],lon_lims4[1],lat_lims4[0],\
        lat_lims4[1]],crs = datacrs)

    ax0.set_title('Aqua GOES\n'+dt_date_str20.strftime('%Y-%m-%d %H:%M'))
    ax1.set_title('Aqua GOES\n'+dt_date_str21.strftime('%Y-%m-%d %H:%M'))
    ax2.set_title('Aqua GOES\n'+dt_date_str22.strftime('%Y-%m-%d %H:%M'))
    ax3.set_title('Aqua GOES\n'+dt_date_str23.strftime('%Y-%m-%d %H:%M'))

    plot_asos_diurnal(ax4, date_str20, 'O05', 'AAT')
    plot_asos_diurnal(ax5, date_str21, 'O05', 'AAT')
    plot_asos_diurnal(ax6, date_str22, 'O05', 'AAT')
    plot_asos_diurnal(ax7, date_str23, 'O05', 'AAT')

    plot_ASOS_locs(ax0,'asos_data_case_locs.csv', color = 'red')
    plot_ASOS_locs(ax1,'asos_data_case_locs.csv', color = 'red')
    plot_ASOS_locs(ax2,'asos_data_case_locs.csv', color = 'red')
    plot_ASOS_locs(ax3,'asos_data_case_locs.csv', color = 'red')

    dt_date_str = dt_date_str22
    xval = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lon'][0] + \
        (plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lon'][1] - \
        plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lon'][0])*0.05
    yval = plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lat'][0] + \
        (plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lat'][1] - \
        plot_limits_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lat'][0])*0.90
    plot_subplot_label(ax0, '(a)',color = 'white', xval = xval, yval = yval, \
        transform = datacrs)
    plot_subplot_label(ax1, '(b)',color = 'white', xval = xval, yval = yval, \
        transform = datacrs)
    plot_subplot_label(ax2, '(c)',color = 'white', xval = xval, yval = yval, \
        transform = datacrs)
    plot_subplot_label(ax3, '(d)',color = 'white', xval = xval, yval = yval, \
        transform = datacrs)

    plot_subplot_label(ax4, '(e)', xval = 1200)
    plot_subplot_label(ax5, '(f)', xval = 1200)
    plot_subplot_label(ax6, '(g)', xval = 1200)
    plot_subplot_label(ax7, '(h)', xval = 1200)

    def plot_goes_line(dt_date, pax):
        local_goes_time = dt_date - timedelta(hours = 7)
        goes_diff = (local_goes_time - datetime(year=2021,month=7,\
            day=22,hour=0,minute=0)).seconds
        print(local_goes_time, goes_diff)
        pax.axvline(goes_diff,color='black',linestyle = '--', lw=2,alpha=0.75,\
            label='GOES')

    plot_goes_line(dt_date_str20, ax4)
    plot_goes_line(dt_date_str21, ax5)
    plot_goes_line(dt_date_str22, ax6)
    plot_goes_line(dt_date_str23, ax7)
 
    fig.tight_layout()
 
    if(save):
        outname = 'goes_asos_meteo_combined_20210722.png'
        fig.savefig(outname, dpi=300)
        print("Saved image",outname)
    else: 
        plt.show() 

# Compare colocated GOES and ob data for two dates
def colocate_comparison(date1, date2, channel = 31):
    # Read in GOES data for both cases
    # ---------------------------------
    dt_date_str1 = datetime.strptime(date1,"%Y%m%d%H%M")
    dt_date_str2 = datetime.strptime(date2,"%Y%m%d%H%M")
    filename1 = plot_limits_dict[dt_date_str1.strftime('%Y-%m-%d')][dt_date_str1.strftime('%H%M')]['goes']
    filename2 = plot_limits_dict[dt_date_str2.strftime('%Y-%m-%d')][dt_date_str2.strftime('%H%M')]['goes']

    GOES_data1 = read_GOES_channel(dt_date_str1.strftime('%Y%m%d%H%M'), channel, zoom = True)
    GOES_data2 = read_GOES_channel(dt_date_str2.strftime('%Y%m%d%H%M'), channel, zoom = True)

    #print(GOES_data1,GOES_data2)

    # Use the colocation code to extract matching data
    # ------------------------------------------------
    compare_data1 = nearest_grid_values(GOES_data1)
    compare_data2 = nearest_grid_values(GOES_data2)

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
    ax2.set_ylabel('Colocated GOES Channel ' + str(channel) + ' Brightness Temp [K] ',color='tab:orange') 
    plt.legend() 
    plt.show()

def split_window(date_str, zoom = True, save = False):
    
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # Read 0.64 μm data
    # --------------- 
    GOES_data1 = read_GOES_channel(date_str, 1, zoom = zoom)

    # Read 11 μm data
    # --------------- 
    GOES_data31 = read_GOES_channel(date_str, 31, zoom = zoom)

    # Read 12 μm data
    # --------------- 
    GOES_data32 = read_GOES_channel(date_str, 32, zoom = zoom)


    GOES_data_div = read_GOES_channel(date_str, 32, zoom = zoom)
    GOES_data_sub = read_GOES_channel(date_str, 32, zoom = zoom)

    # Calculate split window stuff
    # ----------------------------
    GOES_data_sub['data'] = GOES_data31['data'][:,:] - GOES_data32['data'][:,:]
    GOES_data_div['data'] = GOES_data31['data'][:,:] / GOES_data32['data'][:,:]

    GOES_data_sub['variable'] = '11 μm - 12 μm Split Window'
    GOES_data_div['variable'] = '11 μm / 12 μm Split Window'
    GOES_data_sub['colors'] = 'bwr'
    
    # Set up figure
    # ------------- 
    mapcrs = init_proj(date_str)
    plt.close('all')
    fig = plt.figure(figsize=(12,10))
    ax1 = fig.add_subplot(2,3,2,projection = mapcrs) # true color    
    ax2 = fig.add_subplot(2,3,3,projection = mapcrs) # Ch 31
    ax3 = fig.add_subplot(2,3,4,projection = mapcrs) # Ch 1
    ax4 = fig.add_subplot(2,3,5,projection = mapcrs) # Ch 5
    ax0 = fig.add_subplot(2,3,1,projection = mapcrs) # Ch 5
    
    plot_GOES_spatial(GOES_data1, ax0, zoom = zoom, ptitle = '0.64 μm reflectance')
    plot_GOES_spatial(GOES_data31, ax1, zoom = zoom, ptitle = '11 μm brightness temp')
    plot_GOES_spatial(GOES_data32, ax2, zoom = zoom, ptitle = '12 μm brightness temp')
    plot_GOES_spatial(GOES_data_sub,  ax3, zoom = zoom, vmin = -2, vmax = 2, ptitle = '')
    plot_GOES_spatial(GOES_data_div,  ax4, zoom = zoom, vmin = 0.995, vmax = 1, ptitle = '')

    if(save):
        if(zoom):
            zoom_add = '_zoom'
        else:
            zoom_add = ''
        outname = 'goes_split_window_' + date_str + zoom_add + '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

def plot_GOES_temporary(date_str, zoom = True, save = False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # Read 0.64 μm data
    # --------------- 
    GOES_data_ch1  = read_GOES_channel(date_str, 1, zoom = zoom)
    GOES_data_ch5  = read_GOES_channel(date_str, 5, zoom = zoom)
    GOES_data_ch31 = read_GOES_channel(date_str, 31, zoom = zoom)

    # Determine where the smoke is located
    # ------------------------------------
    hash_data, nohash_data = find_plume(date_str) 

    tmp_data1  = np.copy(GOES_data_ch1['data'])
    tmp_data5  = np.copy(GOES_data_ch5['data'])
    tmp_data31 = np.copy(GOES_data_ch31['data'])
    tmp_lat0   = np.copy(GOES_data_ch1['lat'])
    tmp_lon0   = np.copy(GOES_data_ch1['lon'])

    if(not (tmp_data1.shape == tmp_data5.shape == tmp_data31.shape == \
            hash_data.shape)):
        print("shape mismatch")
        shapes = []
        shapes.append(tmp_data1.shape)
        shapes.append(tmp_data5.shape)
        shapes.append(tmp_data31.shape)
        shapes.append(hash_data.shape)

        min_shape = min(shapes)
        print(min_shape)

        tmp_data1  = tmp_data1[:min_shape[0],:min_shape[1]]
        tmp_data5  = tmp_data5[:min_shape[0],:min_shape[1]]
        tmp_data31 = tmp_data31[:min_shape[0],:min_shape[1]]
        tmp_lat0   = tmp_lat0[:min_shape[0],:min_shape[1]]
        tmp_lon0   = tmp_lon0[:min_shape[0],:min_shape[1]]
        hash_data  = hash_data[:min_shape[0],:min_shape[1]]

    max_ch = 350.

    tmp_data1 = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) \
        , tmp_data1)
    tmp_data5 = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) \
        , tmp_data5)
    tmp_data31 = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) \
        , tmp_data31)
    tmp_lat0 = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) \
        , tmp_lat0)
    tmp_lon0 = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) \
        , tmp_lon0)

    # Set up figure
    # -------------
    plt.close('all')
    fig1 = plt.figure(figsize = (10,5))
    ax0 = fig1.add_subplot(1,2,1)
    ax1 = fig1.add_subplot(1,2,2)

    plot_scatter(ax0, tmp_data31, tmp_data1, GOES_data_ch31, GOES_data_ch1, \
        hash_data, xlabel = '11 μm  temperature', \
        ylabel = '0.64 μm reflectance', plot_legend = True)
    plot_scatter(ax1, tmp_data31, tmp_data5, GOES_data_ch31, GOES_data_ch5, \
        hash_data, xlabel = '11 μm brightness temperature', \
        ylabel = '1.24 μm reflectance')

    if(save):
        outname = 'goes_scatter_twopanel_' + date_str + '.png'
        fig1.savefig(outname, dpi = 300)
        print("Saved image",outname)
    else:
        plt.show()

def plot_GOES_temporary_4panel(date_str, zoom = True, composite = True, \
        show_smoke = True, save = False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # Read the true color data
    # ------------------------
    var1, crs1, lat_lims1, lon_lims1 = read_true_color(date_str,\
        composite=composite)


    # Read 0.64 μm data
    # --------------- 
    GOES_data_ch1  = read_GOES_channel(date_str, 1, zoom = zoom)
    GOES_data_ch31 = read_GOES_channel(date_str, 31, zoom = zoom)

    # Determine where the smoke is located
    # ------------------------------------
    hash_data, nohash_data = find_plume(date_str) 

    tmp_data1  = np.copy(GOES_data_ch1['data'])
    tmp_data31 = np.copy(GOES_data_ch31['data'])
    tmp_lat0   = np.copy(GOES_data_ch1['lat'])
    tmp_lon0   = np.copy(GOES_data_ch1['lon'])

    if(not (tmp_data1.shape == tmp_data31.shape == \
            hash_data.shape)):
        print("shape mismatch")
        shapes = []
        shapes.append(tmp_data1.shape)
        shapes.append(tmp_data31.shape)
        shapes.append(hash_data.shape)

        min_shape = min(shapes)
        print(min_shape)

        tmp_data1  = tmp_data1[:min_shape[0],:min_shape[1]]
        tmp_data31 = tmp_data31[:min_shape[0],:min_shape[1]]
        tmp_lat0   = tmp_lat0[:min_shape[0],:min_shape[1]]
        tmp_lon0   = tmp_lon0[:min_shape[0],:min_shape[1]]
        hash_data  = hash_data[:min_shape[0],:min_shape[1]]

    max_ch = 350.

    tmp_data1 = np.ma.masked_where( (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) \
        , tmp_data1)
    tmp_data31 = np.ma.masked_where( (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) \
        , tmp_data31)
    tmp_lat0 = np.ma.masked_where( (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) \
        , tmp_lat0)
    tmp_lon0 = np.ma.masked_where( (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) \
        , tmp_lon0)


    # Set up figure
    # -------------
    mapcrs = init_proj(date_str)
    plt.close('all')
    ##!#fig = plt.figure(figsize=(9,9))
    ##!##ax1 = fig.add_subplot(2,3,1,projection = mapcrs)   # true color    
    ##!#ax1 = fig.add_subplot(2,2,1,projection = crs1)   # true color    
    ##!#ax2 = fig.add_subplot(2,2,2,projection = mapcrs) # Ch 31
    ##!#ax3 = fig.add_subplot(2,2,3,projection = mapcrs) # Ch 1
    ##!#ax4 = fig.add_subplot(2,2,4) # Scatter 
    ##!##ax4 = fig.add_subplot(2,2,5) # Ch 5
    ##!##ax5 = fig.add_subplot(2,2,6) # Ch 5

    # Set up the figure
    # -----------------
    fig = plt.figure(figsize=(9,15))
    gs1 = fig.add_gridspec(nrows = 3, ncols = 2, wspace = 0.40, hspace = 0.30)
    ax1 = fig.add_subplot(gs1[0,0], projection = crs1) # true color 7/22
    ax2 = fig.add_subplot(gs1[0,1], projection = mapcrs) # vis 7/22
    ax3 = fig.add_subplot(gs1[1,0], projection = mapcrs) # IR 7/22
    ax4 = fig.add_subplot(gs1[1,1]) # IR / vis scatter
    ax5 = fig.add_subplot(gs1[2,:])

    # Plot the true-color data for the previous date
    # ----------------------------------------------
    ax1.imshow(var1.data, transform = crs1, extent=(var1.x[0], var1.x[-1], \
        var1.y[-1], var1.y[0]), origin='upper')

    # Zoom in the figure if desired
    # -----------------------------
    if(zoom):
        ax1.set_extent([lon_lims1[0],lon_lims1[1],lat_lims1[0],\
            lat_lims1[1]],crs = datacrs)

    # ----------------------------------------------------------------------
    #
    # Plot the data in the figure
    #
    # ----------------------------------------------------------------------

    # Plot channel 1, 5, 31, and WV data spatial data
    # -----------------------------------------------
    plot_GOES_spatial(GOES_data_ch1,  ax2, zoom = zoom, ptitle = '')
    plot_GOES_spatial(GOES_data_ch31, ax3, zoom = zoom, ptitle = '')

    plot_figure_text(ax1, 'True Color', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 15, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax2, 'GOES 0.64 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 15, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax3, 'GOES 11 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 15, backgroundcolor = 'white', halign = 'right')

    # Plot the scatter between IR and vis
    # -----------------------------------
    plot_scatter(ax4, tmp_data31, tmp_data1, GOES_data_ch31, GOES_data_ch1, \
        hash_data, xlabel = '11 μm brightness temperature', \
        ylabel = '0.64 μm reflectance', plot_legend = True)
    #plot_scatter(ax5, tmp_data31, tmp_data5, GOES_data_ch31, GOES_data_ch5, \
    #    hash_data, xlabel = '11 μm brightness temperature', \
    #    ylabel = '1.24 μm reflectance')

    # ----------------------------------------------------------------------
    #
    # Panel 5: Meteogram
    #
    # ----------------------------------------------------------------------
    plot_asos_diurnal(ax5, date_str, 'O05', 'AAT')

    def plot_goes_line(dt_date, pax):
        local_goes_time = dt_date - timedelta(hours = 7)
        goes_diff = (local_goes_time - datetime(year=2021,month=7,\
            day=22,hour=0,minute=0)).seconds
        print(local_goes_time, goes_diff)
        pax.axvline(goes_diff,color='black',linestyle = '--', lw=2,alpha=0.75,\
            label='GOES')

    plot_goes_line(dt_date_str, ax5)
    ax5.legend() 


    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(date_str) 

        plt.rcParams.update({'hatch.color': 'r'})
        ax2.pcolor(GOES_data_ch1['lon'],GOES_data_ch1['lat'],\
            hash_data1, hatch = '\\\\', alpha=0., transform = datacrs,\
            cmap = 'plasma')
        #ax4.pcolor(GOES_data_ch1['lon'],GOES_data_ch1['lat'],\
        #    hash_data1, hatch = '\\\\', alpha=0., transform = datacrs,\
        #    cmap = 'plasma')

    # Add ASOS site locations
    # -----------------------
    plot_ASOS_locs(ax3,'asos_data_case_locs.csv', color = 'red')

    # Add subplot letter labels
    # -------------------------
    plot_subplot_label(ax1, '(a)', backgroundcolor = 'white')
    plot_subplot_label(ax2, '(b)', backgroundcolor = 'white')
    plot_subplot_label(ax3, '(c)', backgroundcolor = 'white')
    plot_subplot_label(ax4, '(d)', backgroundcolor = 'white')
    plot_subplot_label(ax5, '(e)', xval = 1200, location = 'lower_left')

    #plt.suptitle(dt_date_str.strftime('%d %B %Y'), fontsize = 18)
    #fig.tight_layout()

    if(save):
        outname = 'goes_spatscat2_4panel_' + date_str + '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image",outname)
    else:
        plt.show()
def plot_GOES_temporary(date_str, zoom = True, save = False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # Read 0.64 μm data
    # --------------- 
    GOES_data_ch1  = read_GOES_channel(date_str, 1, zoom = zoom)
    GOES_data_ch5  = read_GOES_channel(date_str, 5, zoom = zoom)
    GOES_data_ch31 = read_GOES_channel(date_str, 31, zoom = zoom)

    # Determine where the smoke is located
    # ------------------------------------
    hash_data, nohash_data = find_plume(date_str) 

    tmp_data1  = np.copy(GOES_data_ch1['data'])
    tmp_data5  = np.copy(GOES_data_ch5['data'])
    tmp_data31 = np.copy(GOES_data_ch31['data'])
    tmp_lat0   = np.copy(GOES_data_ch1['lat'])
    tmp_lon0   = np.copy(GOES_data_ch1['lon'])

    if(not (tmp_data1.shape == tmp_data5.shape == tmp_data31.shape == \
            hash_data.shape)):
        print("shape mismatch")
        shapes = []
        shapes.append(tmp_data1.shape)
        shapes.append(tmp_data5.shape)
        shapes.append(tmp_data31.shape)
        shapes.append(hash_data.shape)

        min_shape = min(shapes)
        print(min_shape)

        tmp_data1  = tmp_data1[:min_shape[0],:min_shape[1]]
        tmp_data5  = tmp_data5[:min_shape[0],:min_shape[1]]
        tmp_data31 = tmp_data31[:min_shape[0],:min_shape[1]]
        tmp_lat0   = tmp_lat0[:min_shape[0],:min_shape[1]]
        tmp_lon0   = tmp_lon0[:min_shape[0],:min_shape[1]]
        hash_data  = hash_data[:min_shape[0],:min_shape[1]]

    max_ch = 350.

    tmp_data1 = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) \
        , tmp_data1)
    tmp_data5 = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) \
        , tmp_data5)
    tmp_data31 = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) \
        , tmp_data31)
    tmp_lat0 = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) \
        , tmp_lat0)
    tmp_lon0 = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) \
        , tmp_lon0)

    # Set up figure
    # -------------
    plt.close('all')
    fig1 = plt.figure(figsize = (10,5))
    ax0 = fig1.add_subplot(1,2,1)
    ax1 = fig1.add_subplot(1,2,2)

    plot_scatter(ax0, tmp_data31, tmp_data1, GOES_data_ch31, GOES_data_ch1, \
        hash_data, xlabel = '11 μm  temperature', \
        ylabel = '0.64 μm reflectance', plot_legend = True)
    plot_scatter(ax1, tmp_data31, tmp_data5, GOES_data_ch31, GOES_data_ch5, \
        hash_data, xlabel = '11 μm brightness temperature', \
        ylabel = '1.24 μm reflectance')

    if(save):
        outname = 'goes_scatter_twopanel_' + date_str + '.png'
        fig1.savefig(outname, dpi = 300)
        print("Saved image",outname)
    else:
        plt.show()