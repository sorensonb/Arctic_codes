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
from scipy.ndimage.filters import gaussian_filter
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

sys.path.append('/home/bsorenson/')
from python_lib import *

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
goes_channel_dict = {
    '1': {
        'limits': [None, None],\
        'name': 'Blue',
        'wavelength': 0.47
    },\
    '2': {
        'limits': [0, 80],\
        'name': 'Red',
        'wavelength': 0.64
    },\
    '3': {
        'limits': [None, None],\
        'name': 'Veggie',
        'wavelength': 0.86
    },\
    '4': {
        'limits': [None, None],\
        'name': 'Cirrus',
        'wavelength': 1.38
    },\
    '5': {
        'limits': [None, None],\
        'name': 'Snow/Ice',
        'wavelength': 1.61
    },\
    '6': {
        'limits': [0, 50],\
        'name': 'Cloud Particle Size',
        'wavelength': 2.25
    },\
    '7': {
        'limits': [None, None],\
        'name': 'Shortwave Window',
        'wavelength': 3.90
    },\
    '8': {
        'limits': [240, 250],\
        'name': 'Upper-Level Water Vapor',
        'wavelength': 6.18
    },\
    '9': {
        'limits': [250, 260],\
        'name': 'Mid-Level Water Vapor',
        'wavelength': 6.95
    },\
    '10': {
        'limits': [255, 270],\
        'name': 'Lower-Level Water Vapor',
        'wavelength': 7.34
    },\
    '11': {
        'limits': [None, None],\
        'name': 'Cloud-Top Phase',
        'wavelength': 8.50
    },\
    '12': {
        'limits': [None, None],\
        'name': 'Ozone',
        'wavelength': 9.61
    },\
    '13': {
        'limits': [270, 330],\
        'name': 'Clean IR Longwave Window',
        'wavelength': 10.35
    },\
    '14': {
        'limits': [None, None],\
        'name': 'IR Longwave Window',
        'wavelength': 11.20
    },\
    '15': {
        'limits': [None, None],\
        'name': 'Dirty IR Longwave Window',
        'wavelength': 12.30
    },\
    '16': {
        'limits': [None, None],\
        'name': 'CO2 Longwave Infrared',
        'wavelength': 13.30
    }
}

for key in goes_channel_dict.keys():
    if(goes_channel_dict[key]['wavelength'] is not None):
        goes_channel_dict[key]['wavelength_label'] = \
            str(goes_channel_dict[key]['wavelength']) + ' μm'
    else:
        goes_channel_dict[key]['wavelength_label'] = ''

goes_area_dict = {
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

##!#def init_proj(date_str):
##!#    #mapcrs = Miller()
##!#    if(date_str == None):
##!#        mapcrs = ccrs.LambertConformal()
##!#    else:
##!#        dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
##!#
##!#        mapcrs = ccrs.LambertConformal(central_longitude = \
##!#            np.mean(goes_area_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['Lon']),\
##!#            central_latitude = \
##!#            np.mean(goes_area_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['Lat']))
##!#
##!#    return mapcrs
##!#
##!#def plot_trend_line(pax, xdata, ydata, color='black', linestyle = '-', \
##!#        slope = 'thiel-sen'):
##!#
##!#    if(slope == 'thiel-sen'):
##!#        res = stats.theilslopes(ydata, xdata, 0.95)
##!#        print("Theil-Sen: {0}x + {1}".format(res[0], res[1]))
##!#
##!#        # Then, plot the trend line on the figure
##!#        pax.plot(xdata, res[1] + res[0] * xdata, \
##!#            color='k', linewidth = 2.5, linestyle = linestyle)
##!#        # Then, plot the trend line on the figure
##!#        pax.plot(xdata, res[1] + res[0] * xdata, \
##!#            color=color, linestyle = linestyle)
##!#    else:
##!#        # First, calculate the trend
##!#        zdata = np.polyfit(xdata, ydata, 1)
##!#
##!#        print("{0}x + {1}".format(*zdata))
##!#
##!#        # Then, plot the trend line on the figure
##!#        pax.plot(np.unique(xdata), np.poly1d(zdata)(np.unique(xdata)), \
##!#            color=color, linestyle = linestyle)
##!#
##!#def plot_subplot_label(ax, label, xval = None, yval = None, transform = None, \
##!#        color = 'black', backgroundcolor = None, fontsize = 14, \
##!#        location = 'upper_left'):
##!#
##!#    if(location == 'upper_left'):
##!#        y_lim = 0.90
##!#        x_lim = 0.05
##!#    elif(location == 'lower_left'):
##!#        y_lim = 0.05
##!#        x_lim = 0.05
##!#    elif(location == 'upper_right'):
##!#        y_lim = 0.90
##!#        x_lim = 0.90
##!#    elif(location == 'lower_right'):
##!#        y_lim = 0.05
##!#        x_lim = 0.90
##!#
##!#    if(xval is None):
##!#        xval = ax.get_xlim()[0] + (ax.get_xlim()[1] - ax.get_xlim()[0]) * x_lim
##!#    if(yval is None):
##!#        yval = ax.get_ylim()[0] + (ax.get_ylim()[1] - ax.get_ylim()[0]) * y_lim
##!#    print('Xval = ',xval, 'Yval = ',yval)
##!#
##!#    if(transform is None):
##!#        if(backgroundcolor is None):
##!#            ax.text(xval,yval,label, \
##!#                color=color, weight='bold', \
##!#                fontsize=fontsize)
##!#        else:
##!#            ax.text(xval,yval,label, \
##!#                color=color, weight='bold', \
##!#                fontsize=fontsize, backgroundcolor = backgroundcolor)
##!#    else:
##!#        if(backgroundcolor is None):
##!#            ax.text(xval,yval,label, \
##!#                color=color, weight='bold', \
##!#                transform = transform, fontsize=fontsize)
##!#        else:
##!#            ax.text(xval,yval,label, \
##!#                color=color, weight='bold', \
##!#                transform = transform, fontsize=fontsize, \
##!#                backgroundcolor = backgroundcolor)
##!#
##!#def plot_figure_text(ax, text, xval = None, yval = None, transform = None, \
##!#        color = 'black', fontsize = 12, backgroundcolor = 'white',\
##!#        halign = 'left'):
##!#
##!#    if(xval is None):
##!#        print(len(text))
##!#        xval = ax.get_xlim()[0] + (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.95
##!#    if(yval is None):
##!#        yval = ax.get_ylim()[0] + (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.05
##!#    print('Xval = ',xval, 'Yval = ',yval)
##!#
##!#    if(transform is None):
##!#        ax.text(xval,yval,text, \
##!#            color=color, weight='bold', \
##!#            fontsize=fontsize, backgroundcolor = backgroundcolor, \
##!#            horizontalalignment = halign)
##!#    else:
##!#        ax.text(xval,yval,text, \
##!#            color=color, weight='bold', \
##!#            transform = transform, fontsize=fontsize, \
##!#            backgroundcolor = backgroundcolor, \
##!#            horizontalalignment = halign)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Download data
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# date_str: YYYYMMDDHHMM
# channels is a list containing the numbers of the desired channels
def download_GOES_bucket(date_str, sat = 'goes17', \
        channels = [2,6,8,9,10,13], dest_dir = \
        '/home/bsorenson/data/GOES/goes17_abi/'):

    # Convert the input date_str to datetime
    # --------------------------------------
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    #if((dt_date_str.minute == 0)):
    #    dt_date_str = dt_date_str - timedelta(minutes = 5)

    # Pull the entire file list for this date and time
    # ------------------------------------------------
    request_add = dt_date_str.strftime('s3://noaa-'+sat+'/ABI-L1b-RadC/2021/%j/%H/')

    # For each desired channel, figure out the closest file time
    # to the input date
    # ----------------------------------------------------------
    try:
        files = subprocess.run(['aws','s3','ls','--no-sign-request',\
            request_add], check=True, \
            capture_output=True).stdout.decode('utf-8').strip().split('\n') 
    except subprocess.CalledProcessError:
        print("ERROR: No ",sat," files for the input DTG",date_str)
        return

    # Remove the timestamps from the file strings
    # -------------------------------------------
    files_only = [tfile.strip().split()[3] for tfile in files]

    # Use the actual file times to get timestamps
    # -------------------------------------------
    file_dates = [datetime.strptime(tfile[27:40],'%Y%j%H%M%S') for tfile in files_only] 

    # Extract the channels
    # --------------------
    file_channels = np.array([int(tfile[19:21]) for tfile in files_only])

    # Figure out where the closest files to the desired time are
    # ----------------------------------------------------------
    time_diffs = np.array([abs((dt_date_str - ddate).total_seconds()) \
        for ddate in file_dates])

    # Select those files. Should result in all 16 channels for a single
    # time
    # -----------------------------------------------------------------
    files_found = np.array(files_only)[np.where(\
        time_diffs == np.min(time_diffs))]

    # Use the input channels to extract only the desired channels
    # -----------------------------------------------------------
    np_channels = np.array(channels) - 1
    good_files = files_found[np_channels]

    # Loop over the selected files and download them
    # ----------------------------------------------
    for gf in good_files:
        request_add = dt_date_str.strftime('s3://noaa-'+sat+\
            '/ABI-L1b-RadC/2021/%j/%H/' + gf)
        #print(request_add)
        cmnd_list = ['aws','s3','cp','--no-sign-request',\
            request_add, dest_dir]
        print(' '.join(cmnd_list))
        subprocess.run(cmnd_list, check=True, \
            capture_output=False)

# begin_date : YYYYMMDDHHMM
# end_date   : YYYYMMDDHHMM
# interval   : minutes between desired images
def auto_GOES_download(begin_date, end_date, interval, sat = 'goes17', \
        channels = [2,6,8,9,10,13], dest_dir = \
        '/home/bsorenson/data/GOES/goes17_abi/'):

    # Convert the input date_str to datetime
    # --------------------------------------
    begin_dt_date = datetime.strptime(begin_date,"%Y%m%d%H%M")
    end_dt_date   = datetime.strptime(end_date,"%Y%m%d%H%M")

    # Using the interval, get the desired file times for each
    # GOES image
    # -------------------------------------------------------
    num_seconds_in_minute = 60
    num_minutes = (end_dt_date - begin_dt_date).total_seconds() / \
        num_seconds_in_minute
    num_imgs = num_minutes / interval


    image_times = [begin_dt_date + timedelta(minutes = interval * ii) \
        for ii in range(int(num_imgs))]

    for ttime in image_times:
        print(ttime)
        download_GOES_bucket(ttime.strftime('%Y%m%d%H%M'), \
            sat = sat, channels = channels, dest_dir = dest_dir)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Miscellaneous plotting functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def plot_zenith_angles():
    p_lats = np.arange(0, 60, 0.5)
    p_lons_137 = np.full(p_lats.shape, -137.2)
    p_lons_121 = np.full(p_lats.shape, -121.0570)
    p_lons_110 = np.full(p_lats.shape, -110.0570)

    # Calculate zenith angles for each point
    vzas_137 = calc_zenith_angle(p_lats, p_lons_137)
    vzas_122 = calc_zenith_angle(p_lats, p_lons_121)
    vzas_110 = calc_zenith_angle(p_lats, p_lons_110)

    # Plot the results
    plt.close('all')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(p_lats, vzas_137)
    ax.plot(p_lats, vzas_122)
    ax.plot(p_lats, vzas_110)
    ax.set_xlabel('Latitude')
    ax.set_ylabel('Viewing zenith angle [degrees]')
    plt.show()

def calc_zenith_angle(phi, lbda):
    # GOES-17 values
    sat_height = 35786.0234375 # km
    sat_lon    = -137.2
    r_e        = 6371.0 # km

    r_dist = r_e + sat_height
    
    # gamma_rad is the angle between the radius vectors of
    # point E (the satellite subpoint along the equator) and point
    # P (the point location)
    gamma_rad = np.arccos(np.cos(np.radians(phi)) * \
        np.cos(np.radians(sat_lon - lbda)))

    # d_dist is the distance from the point on the earth's surface 
    # to the satellite.
    d_dist = r_dist * (1. + (r_e / r_dist)**2. - 2. * \
        (r_e / r_dist) * np.cos(gamma_rad))**0.5 

    # Calculate the viewing zenith angle in degrees
    vza = np.degrees(np.arcsin((r_dist / d_dist) * np.sin(gamma_rad)))

    return vza
    

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
    asos_file = goes_area_dict[GOES_data['cross_date']][GOES_data['file_time']]['asos']
    df = pd.read_csv(asos_file)
    df['valid'] = pd.to_datetime(df['valid'])
    df = df.set_index('valid')

    # Pull the event time from the goes_area_dict
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
        asos_file = goes_area_dict[cross_date][file_date]['asos']

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

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Plotting functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

calib_dict = {
    0.64: 'reflectance',
    1.61: 'reflectance',
    2.25: 'reflectance',
    3.90: 'radiance',
    6.18: 'brightness_temperature',
    6.95: 'brightness_temperature',
    7.34: 'brightness_temperature',
    10.35: 'brightness_temperature'
}   
label_dict = {
    0.64: '%',
    1.61: '%',
    2.25: '%',
    3.90: 'mW m$^{-2}$ Sr$^{-1}$ (cm$^{-1}$)$^{-1}$',
    #6.18: 'mW m$^{-2}$ Sr$^{-1}$ (cm$^{-1}$)$^{-1}$',
    6.18: 'K',
    6.95: 'K',
    7.34: 'K',
    10.35: 'K'
}   
cmap_dict = {
    0.64:  'Greys_r',
    1.61:  'Greys_r',
    2.25:  'Greys_r',
    3.90:  'plasma',
    6.18:  'plasma',
    6.95:  'plasma',
    7.34:  'plasma',
    10.35: 'plasma'
}   
# Channel must be:
# - 0.64  (band 2, Red)
# - 1.61  (band 5, Snow / ice)
# - 2.25  (band 6, Cloud particle size)
# - 3.90  (band 7, Shortwave Window)
# - 6.18  (band 8, Upper-Level Water Vapor)
# - 6.95  (band 9, Mid-Level Water Vapor)
# - 7.34  (band 10, Lower-Level Water Vapor)
# - 10.35 (band 13, Clean IR Longwave Window)
def read_GOES_satpy(date_str, channel, scene_date = None, zoom = True):

    # Extract the channel wavelength using the input string
    # -----------------------------------------------------
    channel = goes_channel_dict[str(channel)]['wavelength']

    # Determine the correct GOES files associated with the date
    # ---------------------------------------------------------
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    dt_date_str_end = dt_date_str + timedelta(minutes = 10)

    try:
        # Use the Satpy find_files_and_readers to grab the files
        # ------------------------------------------------------
        files = find_files_and_readers(start_time = dt_date_str, \
            end_time = dt_date_str_end, base_dir = data_dir, reader = 'abi_l1b')
    except ValueError:
        print("ERROR: no files found for dtg",date_str)
        return

    # Extract the goes true-color plot limits
    # ----------------------------------------
    lat_lims = goes_area_dict[dt_date_str.strftime('%Y-%m-%d')]['goes_Lat']
    lon_lims = goes_area_dict[dt_date_str.strftime('%Y-%m-%d')]['goes_Lon']

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

    ##!## Enhance the image for plotting
    ##!## ------------------------------
    ##!#var = get_enhanced_image(scn[channel]).data
    ##!#var = var.transpose('y','x','bands')

    # Zoom the image on the desired area
    # ----------------------------------
    if(zoom):
        try:
            scn = scn.crop(ll_bbox = (lon_lims[0] + 0.65, lat_lims[0], \
                lon_lims[1] - 0.65, lat_lims[1]))
        except NotImplementedError:
            print("WARNING: zoom didn't work for dtg",date_str,\
                " and channel", channel)


    # Extract the lats, lons, and data
    # -----------------------------------------------------
    lons, lats = scn[channel].attrs['area'].get_lonlats()
    var = scn[channel].data

    # Extract the map projection from the data for plotting
    # -----------------------------------------------------
    crs = scn[channel].attrs['area'].to_cartopy_crs()

    # Extract the appropriate units
    # -----------------------------
    units = label_dict[channel]
    #units = scn[channel].units
    plabel = calib_dict[channel].title() + ' [' + units + ']'

    del scn 

    return var, crs, lons, lats, lat_lims, lon_lims, plabel
    #return var, crs, lat_lims, lon_lims

# channel must be an integer between 1 and 16
def plot_GOES_satpy(date_str, channel, ax = None, var = None, crs = None, \
        lons = None, lats = None, lat_lims = None, lon_lims = None, \
        vmin = None, vmax = None, ptitle = None, plabel = None, \
        labelsize = 10, colorbar = True, zoom=True,save=False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    if(var is None): 
        var, crs, lons, lats, lat_lims, lon_lims, plabel = read_GOES_satpy(date_str, channel)

    # Plot the GOES data
    # ------------------
    in_ax = True 
    if(ax is None):
        in_ax = False
        plt.close('all')
        ax = plt.axes(projection=crs)

    ##!#ax.imshow(var.data, transform = crs, extent=(var.x[0], var.x[-1], \
    ##!#    var.y[-1], var.y[0]), vmin = vmin, vmax = vmax, origin='upper', \
    ##!#    cmap = 'Greys_r')
    #im1 = ax.imshow(var, transform = crs, vmin = vmin, vmax = vmax, \
    im1 = ax.pcolormesh(lons, lats, var, transform = datacrs, \
        vmin = vmin, vmax = vmax, \
        cmap = cmap_dict[goes_channel_dict[str(channel)]['wavelength']], \
        shading = 'auto')
    ax.add_feature(cfeature.STATES)
    if(colorbar):
        cbar = plt.colorbar(im1, ax = ax, pad = 0.03, fraction = 0.052, \
            extend = 'both')
        cbar.set_label(plabel.replace('_',' '), size = labelsize, weight = 'bold')

    # Zoom in the figure if desired
    # -----------------------------
    if(zoom):
    ##!#    ax.set_extent([lon_lims[0]+0.55,lon_lims[1]-0.6,lat_lims[0],lat_lims[1]],\
    ##!#                   crs = ccrs.PlateCarree())
        zoom_add = '_zoom'
    else:
        zoom_add = ''

    # NOTE: commented out after removing the 'enhanced_image' code because
    #       it doesn't work now .
    ##!#ax.coastlines(resolution = '50m')
    ##!#ax.add_feature(cfeature.STATES)
    ##!#ax.add_feature(cfeature.BORDERS)
    if(ptitle is None):
        ax.set_title('GOES-17\n'+dt_date_str.strftime('%Y-%m-%d %H:%M'))
    else:
        ax.set_title(ptitle)

    if(not in_ax): 
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

def get_GOES_data_lat_lon(date_str, dlat, dlon, channel):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    var0, crs0, lons0, lats0, lat_lims, lon_lims, plabel0 = read_GOES_satpy(date_str, channel)

    #cd_idx = nearest_gridpoint(dlat, dlon,lats0, lons0)
    #goes_val = np.array(var0)[cd_idx]
    if(not isinstance(dlat, list)):
        dlat = [dlat]
        dlon = [dlon]

    #for tdlat, tdlon in zip(dlat, dlon):
    cd_idx = [nearest_gridpoint(tdlat, tdlon, \
        lats0, lons0) for tdlat, tdlon in zip(dlat, dlon)]
    goes_val = np.array([np.array(var0)[cidx] for cidx in cd_idx]).squeeze()
    goes_lat = np.array([np.array(lats0)[cidx] for cidx in cd_idx]).squeeze()
    goes_lon = np.array([np.array(lons0)[cidx] for cidx in cd_idx]).squeeze()

    return goes_val, goes_lat, goes_lon

def plot_GOES_satpy_6panel(date_str, ch1, ch2, ch3, ch4, ch5, ch6, \
        zoom = True, save_dir = './', save = False):
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    plt.close('all')
    fig1 = plt.figure(figsize = (10,6))
    var0, crs0, lons0, lats0, lat_lims, lon_lims, plabel0 = read_GOES_satpy(date_str, ch1)
    var1, crs1, lons1, lats1, lat_lims, lon_lims, plabel1 = read_GOES_satpy(date_str, ch2)
    var2, crs2, lons2, lats2, lat_lims, lon_lims, plabel2 = read_GOES_satpy(date_str, ch3)
    var3, crs3, lons3, lats3, lat_lims, lon_lims, plabel3 = read_GOES_satpy(date_str, ch4)
    var4, crs4, lons4, lats4, lat_lims, lon_lims, plabel4 = read_GOES_satpy(date_str, ch5)
    var5, crs5, lons5, lats5, lat_lims, lon_lims, plabel5 = read_GOES_satpy(date_str, ch6)

    ax0 = fig1.add_subplot(2,3,1, projection = crs0)
    ax1 = fig1.add_subplot(2,3,2, projection = crs1)
    ax2 = fig1.add_subplot(2,3,3, projection = crs2)
    ax3 = fig1.add_subplot(2,3,4, projection = crs3)
    ax4 = fig1.add_subplot(2,3,5, projection = crs4)
    ax5 = fig1.add_subplot(2,3,6, projection = crs5)

    ##!#ax1.set_title('GOES-17 Band ' + str(ch2) + '\n' + \
    ##!#    goes_channel_dict[str(ch2)]['name'] + '\n' + \
    labelsize = 11
    font_size = 10
    plot_GOES_satpy(date_str, ch1, ax = ax0, var = var0, crs = crs0, \
        lons = lons0, lats = lats0, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)
    plot_GOES_satpy(date_str, ch2, ax = ax1, var = var1, crs = crs0, \
        lons = lons1, lats = lats1, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = 5, vmax = 50, ptitle = '', plabel = plabel1, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)
    plot_GOES_satpy(date_str, ch3, ax = ax2, var = var2, crs = crs0, \
        lons = lons2, lats = lats2, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = 270, vmax = 330, ptitle = '', plabel = plabel2, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)
    plot_GOES_satpy(date_str, ch4, ax = ax3, var = var3, crs = crs0, \
        lons = lons3, lats = lats3, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = 240, vmax = 250, ptitle = '', plabel = plabel3, \
        colorbar = True, labelsize = labelsize, zoom=True,save=False)
    plot_GOES_satpy(date_str, ch5, ax = ax4, var = var4, crs = crs0, \
        lons = lons4, lats = lats4, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = 250, vmax = 260, ptitle = '', plabel = plabel4, \
        colorbar = True, labelsize = labelsize, zoom=True,save=False)
    plot_GOES_satpy(date_str, ch6, ax = ax5, var = var5, crs = crs0, \
        lons = lons5, lats = lats5, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = 255, vmax = 270, ptitle = '', plabel = plabel5, \
        colorbar = True, labelsize = labelsize, zoom=True,save=False)

    plot_figure_text(ax0, 'GOES-17 ' + \
        str(goes_channel_dict[str(ch1)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    plot_figure_text(ax1, 'GOES-17 ' + \
        str(goes_channel_dict[str(ch2)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    plot_figure_text(ax2, 'GOES-17 ' + \
        str(goes_channel_dict[str(ch3)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    plot_figure_text(ax3, 'GOES-17 ' + \
        str(goes_channel_dict[str(ch4)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    plot_figure_text(ax4, 'GOES-17 ' + \
        str(goes_channel_dict[str(ch5)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    plot_figure_text(ax5, 'GOES-17 ' + \
        str(goes_channel_dict[str(ch6)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')

    plot_subplot_label(ax0,  '(a)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax1,  '(b)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax2,  '(c)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax3,  '(d)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax4,  '(e)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax5,  '(f)', backgroundcolor = 'white', fontsize = font_size)

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    #
    # Extract the WV values in each position
    #
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
   
    # "Coldest" values
    # -------------
    c_lon_stn = -120.6530
    c_lat_stn = 40.7595
    cd_idx = nearest_gridpoint(c_lat_stn, c_lon_stn,\
        lats3, lons3)
    # "Cold" values
    # -------------
    lon_stn = -120.3877
    lat_stn = 41.2456
    c_idx = nearest_gridpoint(lat_stn, lon_stn,\
        lats3, lons3)
    # "Warm" values
    # -------------
    w_lon_stn = -120.9810
    w_lat_stn = 41.20980

    #w_lat_stn = 40.50900
    w_idx = nearest_gridpoint(w_lat_stn, w_lon_stn,\
        lats3, lons3)

    ax2.plot(c_lon_stn, c_lat_stn,
             color='tab:green', linewidth=2, marker='o',
             transform=datacrs)
    ax2.plot(lon_stn, lat_stn,
             color='tab:blue', linewidth=2, marker='o',
             transform=datacrs)
    ax2.plot(w_lon_stn, w_lat_stn,
             color='tab:purple', linewidth=2, marker='o',
             transform=datacrs)
    ax3.plot(lon_stn, lat_stn,
             color='tab:blue', linewidth=2, marker='o',
             transform=datacrs)
    ax3.plot(w_lon_stn, w_lat_stn,
             color='tab:purple', linewidth=2, marker='o',
             transform=datacrs)
    ax4.plot(lon_stn, lat_stn,
             color='tab:blue', linewidth=2, marker='o',
             transform=datacrs)
    ax4.plot(w_lon_stn, w_lat_stn,
             color='tab:purple', linewidth=2, marker='o',
             transform=datacrs)
    ax5.plot(lon_stn, lat_stn,
             color='tab:blue', linewidth=2, marker='o',
             transform=datacrs)
    ax5.plot(w_lon_stn, w_lat_stn,
             color='tab:purple', linewidth=2, marker='o',
             transform=datacrs)

    print("TIR")
    print("     Cold - ", np.array(var2)[cd_idx])
    print("     Warm - ", np.array(var2)[w_idx])
    print("Upper WV")
    print("     Cold - ", np.array(var3)[c_idx])
    print("     Warm - ", np.array(var3)[w_idx])
    print("Mid   WV")
    print("     Cold - ", np.array(var4)[c_idx])
    print("     Warm - ", np.array(var4)[w_idx])
    print("Lower WV")
    print("     Cold - ", np.array(var5)[c_idx])
    print("     Warm - ", np.array(var5)[w_idx])

    lon_stn = -120.7605
    lat_stn = 41.2098

    # Zoom in the figure if desired
    # -----------------------------
    if(zoom):
    ##!#    ax0.set_extent([lon_lims[0],lon_lims[1],lat_lims[0],lat_lims[1]],\
    ##!#                   crs = ccrs.PlateCarree())
    ##!#    ax1.set_extent([lon_lims[0],lon_lims[1],lat_lims[0],lat_lims[1]],\
    ##!#                   crs = ccrs.PlateCarree())
    ##!#    ax2.set_extent([lon_lims[0],lon_lims[1],lat_lims[0],lat_lims[1]],\
    ##!#                   crs = ccrs.PlateCarree())
    ##!#    ax3.set_extent([lon_lims[0],lon_lims[1],lat_lims[0],lat_lims[1]],\
    ##!#                   crs = ccrs.PlateCarree())
    ##!#    ax4.set_extent([lon_lims[0],lon_lims[1],lat_lims[0],lat_lims[1]],\
    ##!#                   crs = ccrs.PlateCarree())
    ##!#    ax5.set_extent([lon_lims[0],lon_lims[1],lat_lims[0],lat_lims[1]],\
    ##!#                   crs = ccrs.PlateCarree())
        zoom_add = '_zoom'
    else:
        zoom_add = ''

    fig1.tight_layout()

    if(save):
        outname = save_dir + 'goes17_'+date_str+'_6panel.png'
        fig1.savefig(outname, dpi = 300)
        print('Saved image', outname)
    else:
        plt.show()
   
def plot_GOES_6panel_auto(begin_date, end_date, ch1 = 2, ch2 = 6, \
        ch3 = 13, ch4 = 8, ch5 = 9, ch6 = 10, save_dir = './', \
        save = False):

    # Convert the input date_str to datetime
    # --------------------------------------
    begin_dt_date = datetime.strptime(begin_date,"%Y%m%d%H%M")
    end_dt_date   = datetime.strptime(end_date,"%Y%m%d%H%M")

    # Find all downloaded GOES filenames that are between these
    # two dates
    # ---------------------------------------------------------
    all_files = glob('/home/bsorenson/data/GOES/goes17_abi/*.nc')
    all_dates = [datetime.strptime(ffile.strip().split('/')[-1][27:40],\
        '%Y%j%H%M%S') for ffile in all_files]
   
    # Get just the file dates, with only 1 date per channel
    # ----------------------------------------------------- 
    unique_dates = sorted(list(set(all_dates)))

    in_all = np.array([((tdate > begin_dt_date) & (tdate < end_dt_date)) \
        for tdate in all_dates])
    in_unique = np.array([((tdate > begin_dt_date) & (tdate < end_dt_date)) \
        for tdate in unique_dates])
    
    in_all_idx = np.where(in_all == True)
    in_unique_idx = np.where(in_unique == True)

    if(len(in_all_idx[0]) == 0):
        print("ERROR: no data between",begin_date, end_date)
        print("     Run auto_GOES_download to get the data")
        return

    # good_all_files contains the filenames for all files (including
    # all channels) that have timestamps within the desired range
    #
    # good_unique_times contains just a single time for all channel
    # files that are within the range
    good_all_files    = np.array(all_files)[in_all_idx]
    good_all_channels = np.array([int(\
        tfile.strip().split('/')[-1][19:21]) for tfile in good_all_files])
    good_unique_times = np.array(unique_dates)[in_unique_idx[0]]
 
    # Check to make sure that there are channel files for all
    # desired times.
    # -------------------------------------------------------
    if(not ch1 in good_all_channels):
        print("ERROR: data for ",ch1," not downloaded")
        return
    if(not ch2 in good_all_channels):
        print("ERROR: data for ",ch2," not downloaded")
        return
    if(not ch3 in good_all_channels):
        print("ERROR: data for ",ch3," not downloaded")
        return
    if(not ch4 in good_all_channels):
        print("ERROR: data for ",ch4," not downloaded")
        return
    if(not ch5 in good_all_channels):
        print("ERROR: data for ",ch5," not downloaded")
        return
    if(not ch6 in good_all_channels):
        print("ERROR: data for ",ch6," not downloaded")
        return

    print("All data present")
    
    for ttime in good_unique_times:
        print(ttime.strftime('%Y%m%d%H%M'))
        date_str = ttime.strftime('%Y%m%d%H%M')
        if(date_str == '202107202126'):
            print("Not making image for this time")
        else:
            plot_GOES_satpy_6panel(date_str, ch1, ch2, ch3, ch4, \
                ch5, ch6, zoom = True, save_dir = save_dir, \
            save = True)

def plot_GOES_time_series_points(GOES_dict, time_idx = 20, \
        ch_idx = 0, save_dir = './', save = False):

    # Read GOES image data
    # --------------------
    print(GOES_dict['dt_dates'][time_idx])
    date_str = GOES_dict['dt_dates'][time_idx].strftime('%Y%m%d%H%M')
    var, crs, lons, lats, lat_lims, lon_lims, plabel = \
        read_GOES_satpy(date_str, GOES_dict['channels'][ch_idx])

    plt.close('all')
    fig = plt.figure(figsize = (11, 3.5))
    gs = fig.add_gridspec(nrows = 1, ncols = 4)
    ax1  = fig.add_subplot(gs[3], projection = crs)   # true color    
    ax2  = fig.add_subplot(gs[0:3]) # Ch 1
    #ax1 = fig.add_subplot(1,2,1, projection = crs)
    #ax2 = fig.add_subplot(1,2,2)

    labelsize = 10
    font_size = 8
    plot_GOES_satpy(date_str, GOES_dict['channels'][ch_idx], \
        ax = ax1, var = var, crs = crs, \
        lons = lons, lats = lats, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = goes_channel_dict[\
            str(GOES_dict['channels'][ch_idx])]['limits'][0], \
        vmax = goes_channel_dict[\
            str(GOES_dict['channels'][ch_idx])]['limits'][1], \
        ptitle = '', plabel = plabel, \
        #vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)

    # Add the locations of the points
    # -------------------------------
    for tlat, tlon in zip(GOES_dict['plat'], GOES_dict['plon']):
        print(tlat, tlon)
        ax1.plot(tlon, tlat, linewidth=2, markersize = 8, marker='.',
                 color = 'black', transform=datacrs)
        ax1.plot(tlon, tlat, linewidth=2, markersize = 5, marker='.',
                 transform=datacrs)
    plot_figure_text(ax1, 'GOES-17 ' + \
        str(goes_channel_dict[\
        str(GOES_dict['channels'][ch_idx])]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    ax1.set_title(GOES_dict['dt_dates'][time_idx].strftime(\
        '%Y-%m-%d %H:%MZ'), fontsize = 10)

    # Plot the first 2 channels
    # -------------------------
    for jj in range(GOES_dict['data'].shape[2]):
        ax2.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx,jj])
    #ax2.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx,1])
    #ax2.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx,2])
    ax2.axvline(GOES_dict['dt_dates'][time_idx], color = 'black',\
        linestyle = ':')
    ax2.set_ylabel(plabel.replace('_',' '), \
        size = labelsize, weight = 'bold')
    ax2.grid()
    ax2.xaxis.set_major_formatter(DateFormatter('%m/%d\n%H:%MZ'))
    ax2.tick_params(axis="x", labelsize = 9)
    ##!#    str(GOES_dict['channels'][ch_idx])]['wavelength']) + ' μm', \
    ts_title = 'GOES-17 ' + \
        goes_channel_dict[\
        str(GOES_dict['channels'][ch_idx])]['name'] + \
        ' (' + \
        str(goes_channel_dict[\
        str(GOES_dict['channels'][ch_idx])]['wavelength']) + \
        ' μm)'
        # + ' '.join(plabel.split()[0].split('_'))
    ax2.set_title(ts_title, fontsize = 10)
    fig.autofmt_xdate()

    fig.tight_layout()

    if(save):
        outname = save_dir + 'goes_time_series_points_ch' + \
            str(GOES_dict['channels'][ch_idx]) + '_' + \
            GOES_dict['dt_dates'][time_idx].strftime('%Y%m%d%H%M') + \
            '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

def plot_GOES_time_series_points_auto(GOES_dict, ch_idx, \
        save_dir = './'):

    # Loop over all the dates and save the images
    # -------------------------------------------
    for ii in range(len(GOES_dict['dt_dates'])):
        if(GOES_dict['dt_dates'][ii].strftime('%Y%m%d%H%M') == \
            '202107202126'):
            continue
        else:
            plot_GOES_time_series_points(GOES_dict, time_idx = ii, \
                ch_idx = ch_idx, save_dir = save_dir, save = True)

#def plot_GOES_time_series_points(GOES_dict, time_idx = 20, \
#        ch_idx = 0, save = False):

#def plot_GOES_time_series(dt_dates, goes_data, channels, save = False):
def plot_GOES_time_series_channels(GOES_dict, time_idx = 20, \
        idx = 0, save_dir = './', save = False):

    plt.close('all')
    fig = plt.figure(figsize = (10, 6))
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)

    # Plot the first 2 channels
    # -------------------------
    ln1 = ax1.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,0,idx], \
        label = str(goes_channel_dict[\
        str(GOES_dict['channels'][0])]['wavelength']) + \
        ' μm')
    ln2 = ax1.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,1,idx], \
        label = str(goes_channel_dict[\
        str(GOES_dict['channels'][1])]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:blue')
    ax12 = ax1.twinx()
    ln3 = ax12.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,2,idx], \
        label = str(goes_channel_dict[\
        str(GOES_dict['channels'][2])]['wavelength']) + \
        ' μm', color = 'tab:orange')
    ax1.set_ylabel('Reflectance [%]', color = 'tab:blue')
    ax12.set_ylabel('Brightness Temperature [K]', color = 'tab:orange')
    ax1.grid()
    ax1.xaxis.set_major_formatter(DateFormatter('%m/%d\n%H:%MZ'))
    ax1.tick_params(axis="x", labelsize = 9)
    ax1.tick_params(axis="y", color = 'tab:blue', labelcolor = 'tab:blue')
    ax12.tick_params(axis="y", color = 'tab:orange', labelcolor = 'tab:orange')

    lns = ln1+ln2+ln3
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc=0, fontsize = 9)

    ax2.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,3,idx], \
        #label = str(int(GOES_dict['channels'][3])))
        label = str(goes_channel_dict[\
        str(GOES_dict['channels'][3])]['wavelength']) + \
        ' μm')
    ax2.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,4,idx], \
        #label = str(int(GOES_dict['channels'][4])))
        label = str(goes_channel_dict[\
        str(GOES_dict['channels'][4])]['wavelength']) + \
        ' μm')
    #ax22 = ax2.twinx()
    ax2.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,5,idx], \
        #label = str(int(GOES_dict['channels'][5])))
        label = str(goes_channel_dict[\
        str(GOES_dict['channels'][5])]['wavelength']) + \
        ' μm')
    ax2.set_ylabel('Brightness Temperature [K]')
    ax2.grid()
    ax2.xaxis.set_major_formatter(DateFormatter('%m/%d\n%H:%MZ'))
    ax2.tick_params(axis="x", labelsize = 9)
    ax2.legend(fontsize = 9) 

    names = ['Blue','Orange','Green','Red','Purple','Brown','Pink','Grey']

    #if(idx == 0):
    #    point_name = 'Clearest Pixel'
    #elif(idx == 1):
    #    point_name = 'Smokiest pixel'
    #else:
    #    point_name = 'Lightly Smoky Pixel'
    #plt.title(point_name)
    plt.title(names[idx])

    fig.autofmt_xdate()
    fig.tight_layout()

    if(save):
        outname = save_dir + 'goes_time_series_channels_pt' + \
            str(idx) + '.png'
            #GOES_dict['dt_dates'][time_idx].strftime('%Y%m%d%H%M') + \
            #'.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

def plot_GOES_time_series_mesh(GOES_dict, ch_idx1 = 1, \
        ch_idx2 = 0, save_dir = './', \
        date_idx = 23, sigma = 0.8, save = False):

    channel = 2
    #dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
    dt_date_str = GOES_dict['dt_dates'][date_idx]
    date_str = dt_date_str.strftime('%Y%m%d%H%M')
    var, crs, lons, lats, lat_lims, lon_lims, plabel = \
        read_GOES_satpy(date_str, channel)
    var2, crs2, lons2, lats2, lat_lims2, lon_lims2, plabel2 = \
        read_GOES_satpy(date_str, 13)

    plt.close('all')
    fig = plt.figure(figsize = (11, 3.5))
    gs = fig.add_gridspec(nrows = 1, ncols = 4)
    ax2  = fig.add_subplot(gs[3], projection = crs)   # true color    
    ax1  = fig.add_subplot(gs[0:3]) # Ch 1

    if(GOES_dict['channels'][ch_idx1] > 6):
        cmap1 = 'plasma'
        cmap2 = 'Greys_r'
    else:
        cmap1 = 'Greys_r'
        cmap2 = 'plasma'

    mesh = ax1.pcolormesh(GOES_dict['dt_dates'], \
        GOES_dict['plat'][:],\
        GOES_dict['data'][:,ch_idx1,:].T, cmap = cmap1, \
        vmin = goes_channel_dict[str(GOES_dict['channels'][ch_idx1])\
            ]['limits'][0], \
        vmax = goes_channel_dict[str(GOES_dict['channels'][ch_idx1])\
            ]['limits'][1], \
        shading = 'auto')
   
    ax1.set_ylabel('Latitude') 
    ax1.xaxis.set_major_formatter(DateFormatter('%m/%d\n%H:%MZ'))
    ax1.tick_params(axis="x", labelsize = 9)

    cbar = plt.colorbar(mesh, ax = ax1, pad = 0.03, fraction = 0.052, \
        extend = 'both', label = 'Brightness Temperature [K]')

    ##!## Add 'x's for the location of the smoke in the time series
    ##!## ---------------------------------------------------------
    ##!#tmp_data = np.full(GOES_dict['data'][:,ch_idx2,:].T.shape, np.nan)
    ##!#in_smoke = (GOES_dict['data'][:,ch_idx2,:].T > 22.)
    ##!#tmp_data[in_smoke] = 1.
    ##!#tmp_data[~in_smoke] = 0.

    ##!#mask_tmp_data = np.ma.masked_invalid(tmp_data)

    ##!#hash_data   = np.ma.masked_where(mask_tmp_data != 1, mask_tmp_data)
    #nohash_data = np.ma.masked_where(mask_tmp_data == 1, mask_tmp_data)

    ##!#hash0 = ax1.pcolor(GOES_dict['dt_dates'], GOES_dict['plat'][:],\
    ##!#    hash_data, hatch = 'xx', alpha=0., shading = 'auto')
    smooth_data = gaussian_filter(GOES_dict['data'][:,ch_idx2,:].T, \
        sigma = sigma)
    cntr = ax1.contour(GOES_dict['dt_dates'], \
        GOES_dict['plat'][:],\
        smooth_data, cmap = cmap2, \
        vmin = goes_channel_dict[str(GOES_dict['channels'][ch_idx2])\
            ]['limits'][0], \
        vmax = goes_channel_dict[str(GOES_dict['channels'][ch_idx2])\
            ]['limits'][1])
    ax1.clabel(cntr, cntr.levels, inline=True, fontsize=8)

    # Add vertical line at the image time
    # -----------------------------------
    ax1.axvline(dt_date_str, color = 'black',\
        linestyle = ':')

    #cbar.set_label(plabel.replace('_',' '), size = labelsize, weight = 'bold')

    # Plot the GOES 0.64 micron image
    # -------------------------------
    labelsize = 10
    font_size = 8
    plot_GOES_satpy(date_str, 2, \
        ax = ax2, var = var, crs = crs, \
        lons = lons, lats = lats, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = goes_channel_dict['2']['limits'][0], \
        vmax = goes_channel_dict['2']['limits'][1], \
        ptitle = '', plabel = plabel, \
        #vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)

    # Smooth the other data too
    smooth_data2 = gaussian_filter(var2, \
        sigma = 1.0)
    cntr2 = ax2.contour(lons2, lats2, smooth_data2, \
        np.arange(280, 320, 5), transform = datacrs, cmap = 'plasma', \
        alpha = 0.5)
    ax2.clabel(cntr2, cntr2.levels, inline=True, fontsize=8)

    plot_figure_text(ax2, 'GOES-17 ' + \
        str(goes_channel_dict['2']['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    ax2.set_title(dt_date_str.strftime(\
        '%Y-%m-%d %H:%MZ'), fontsize = 10)
    # Add line showing the cross section
    # ----------------------------------
    max_lat = np.max(GOES_dict['plat'])
    max_lon = np.min(GOES_dict['plon'])
    min_lat = np.min(GOES_dict['plat'])
    min_lon = np.max(GOES_dict['plon'])

    ax2.plot([max_lon, min_lon], [max_lat, min_lat],
             color='black', linewidth=3, marker='o',
             transform=datacrs, markersize = 5, 
             )
    ax2.plot([max_lon, min_lon], [max_lat, min_lat],
             color='cyan', linewidth=1, marker='o',
             transform=datacrs, markersize = 3, 
             )


    fig.autofmt_xdate()
    fig.tight_layout()

    plt.show()

    ##!## Plot the first 2 channels
    ##!## -------------------------
    ##!#ln1 = ax1.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,0,idx], \
    ##!#    label = str(goes_channel_dict[\
    ##!#    str(GOES_dict['channels'][0])]['wavelength']) + \
    ##!#    ' μm')
    ##!#ln2 = ax1.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,1,idx], \
    ##!#    label = str(goes_channel_dict[\
    ##!#    str(GOES_dict['channels'][1])]['wavelength']) + \
    ##!#    ' μm', linestyle = '--', color = 'tab:blue')
    ##!#ax12 = ax1.twinx()
    ##!#ln3 = ax12.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,2,idx], \
    ##!#    label = str(goes_channel_dict[\
    ##!#    str(GOES_dict['channels'][2])]['wavelength']) + \
    ##!#    ' μm', color = 'tab:orange')
    ##!#ax1.set_ylabel('Reflectance [%]', color = 'tab:blue')
    ##!#ax12.set_ylabel('Brightness Temperature [K]', color = 'tab:orange')
    ##!#ax1.grid()
    ##!#ax1.xaxis.set_major_formatter(DateFormatter('%m/%d\n%H:%MZ'))
    ##!#ax1.tick_params(axis="x", labelsize = 9)
    ##!#ax1.tick_params(axis="y", color = 'tab:blue', labelcolor = 'tab:blue')
    ##!#ax12.tick_params(axis="y", color = 'tab:orange', labelcolor = 'tab:orange')

    ##!#lns = ln1+ln2+ln3
    ##!#labs = [l.get_label() for l in lns]
    ##!#ax1.legend(lns, labs, loc=0, fontsize = 9)

    ##!#ax2.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,3,idx], \
    ##!#    #label = str(int(GOES_dict['channels'][3])))
    ##!#    label = str(goes_channel_dict[\
    ##!#    str(GOES_dict['channels'][3])]['wavelength']) + \
    ##!#    ' μm')
    ##!#ax2.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,4,idx], \
    ##!#    #label = str(int(GOES_dict['channels'][4])))
    ##!#    label = str(goes_channel_dict[\
    ##!#    str(GOES_dict['channels'][4])]['wavelength']) + \
    ##!#    ' μm')
    ##!##ax22 = ax2.twinx()
    ##!#ax2.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,5,idx], \
    ##!#    #label = str(int(GOES_dict['channels'][5])))
    ##!#    label = str(goes_channel_dict[\
    ##!#    str(GOES_dict['channels'][5])]['wavelength']) + \
    ##!#    ' μm')
    ##!#ax2.set_ylabel('Brightness Temperature [K]')
    ##!#ax2.grid()
    ##!#ax2.xaxis.set_major_formatter(DateFormatter('%m/%d\n%H:%MZ'))
    ##!#ax2.tick_params(axis="x", labelsize = 9)
    ##!#ax2.legend(fontsize = 9) 

    ##!#names = ['Blue','Orange','Green','Red','Purple','Brown','Pink','Grey']

    ##!##if(idx == 0):
    ##!##    point_name = 'Clearest Pixel'
    ##!##elif(idx == 1):
    ##!##    point_name = 'Smokiest pixel'
    ##!##else:
    ##!##    point_name = 'Lightly Smoky Pixel'
    ##!##plt.title(point_name)
    ##!#plt.title(names[idx])

    ##!#fig.autofmt_xdate()
    ##!#fig.tight_layout()

    ##!#if(save):
    ##!#    outname = save_dir + 'goes_time_series_channels_pt' + \
    ##!#        str(idx) + '.png'
    ##!#        #GOES_dict['dt_dates'][time_idx].strftime('%Y%m%d%H%M') + \
    ##!#        #'.png'
    ##!#    fig.savefig(outname, dpi = 300)
    ##!#    print("Saved image", outname)
    ##!#else:
    ##!#    plt.show()

def plot_GOES_cross_channels(GOES_dict, time_idx = 20, \
        ch_idx1 = 0, ch_idx2 = 1, save_dir = './', save = False):

    channel1 = GOES_dict['channels'][ch_idx1]
    channel2 = GOES_dict['channels'][ch_idx2]
    #dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
    dt_date_str = GOES_dict['dt_dates'][time_idx]
    date_str = dt_date_str.strftime('%Y%m%d%H%M')
    var1, crs1, lons1, lats1, lat_lims1, lon_lims1, plabel1 = \
        read_GOES_satpy(date_str, channel1)
    var2, crs2, lons2, lats2, lat_lims2, lon_lims2, plabel2 = \
        read_GOES_satpy(date_str, channel2)

    plt.close('all')
    fig = plt.figure(figsize = (9, 6))
    gs = fig.add_gridspec(nrows = 2, ncols = 4)
    ax2  = fig.add_subplot(gs[0, 0:2], projection = crs1)   # true color    
    ax3  = fig.add_subplot(gs[0, 2:4], projection = crs2)   # true color    
    ax1  = fig.add_subplot(gs[1, 0:4]) # Ch 1


    ln1 = ax1.plot(GOES_dict['plat'], \
            GOES_dict['data'][time_idx,ch_idx1,:], \
            label = str(goes_channel_dict[\
            str(GOES_dict['channels'][ch_idx1])]['wavelength']) + \
            ' μm')
    ax12 = ax1.twinx()
    ln2 = ax12.plot(GOES_dict['plat'], \
        GOES_dict['data'][time_idx,ch_idx2,:], \
        label = str(goes_channel_dict[\
        str(GOES_dict['channels'][ch_idx2])]['wavelength']) + \
        ' μm', color = 'tab:orange')
    ##!#ln2 = ax1.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,1,idx], \
    ##!#    label = str(goes_channel_dict[\
    ##!#    str(GOES_dict['channels'][1])]['wavelength']) + \
    ##!#    ' μm', linestyle = '--', color = 'tab:blue')
    ax1.set_ylabel('Reflectance [%]', color = 'tab:blue')
    ax12.set_ylabel('Brightness Temperature [K]', color = 'tab:orange')
    ax1.grid()
    #ax1.xaxis.set_major_formatter(DateFormatter('%m/%d\n%H:%MZ'))
    ax1.tick_params(axis="x", labelsize = 9)
    ax1.tick_params(axis="y", color = 'tab:blue', labelcolor = 'tab:blue')
    ax12.tick_params(axis="y", color = 'tab:orange', labelcolor = 'tab:orange')

    # Plot the GOES 0.64 micron image
    # -------------------------------
    labelsize = 10
    font_size = 8
    plot_GOES_satpy(date_str, channel1, \
        ax = ax2, var = var1, crs = crs1, \
        lons = lons1, lats = lats1, lat_lims = lat_lims1, 
        lon_lims = lon_lims1, \
        vmin = goes_channel_dict[str(channel1)]['limits'][0], \
        vmax = goes_channel_dict[str(channel1)]['limits'][1], \
        ptitle = '', plabel = plabel1, \
        #vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)
    plot_figure_text(ax2, 'GOES-17 ' + \
        str(goes_channel_dict[str(channel1)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    ax2.set_title(dt_date_str.strftime(\
        '%Y-%m-%d %H:%MZ'), fontsize = 10)
    # Add line showing the cross section
    # ----------------------------------
    max_lat = np.max(GOES_dict['plat'])
    max_lon = np.min(GOES_dict['plon'])
    min_lat = np.min(GOES_dict['plat'])
    min_lon = np.max(GOES_dict['plon'])

    ax2.plot([max_lon, min_lon], [max_lat, min_lat],
             color='black', linewidth=3, marker='o',
             transform=datacrs, markersize = 5, 
             )
    ax2.plot([max_lon, min_lon], [max_lat, min_lat],
             color='cyan', linewidth=1, marker='o',
             transform=datacrs, markersize = 3, 
             )

    plot_GOES_satpy(date_str, channel2, \
        ax = ax3, var = var2, crs = crs2, \
        lons = lons2, lats = lats2, lat_lims = lat_lims2, \
        lon_lims = lon_lims2, \
        vmin = goes_channel_dict[str(channel2)]['limits'][0], \
        vmax = goes_channel_dict[str(channel2)]['limits'][1], \
        ptitle = '', plabel = plabel2, \
        #vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)
    plot_figure_text(ax3, 'GOES-17 ' + \
        str(goes_channel_dict[str(channel2)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    ax3.set_title(dt_date_str.strftime(\
        '%Y-%m-%d %H:%MZ'), fontsize = 10)

    # Add line showing the cross section
    # ----------------------------------
    ax3.plot([max_lon, min_lon], [max_lat, min_lat],
             color='black', linewidth=3, marker='o',
             transform=datacrs, markersize = 5, 
             )
    ax3.plot([max_lon, min_lon], [max_lat, min_lat],
             color='cyan', linewidth=1, marker='o',
             transform=datacrs, markersize = 3, 
             )

    fig.autofmt_xdate()
    fig.tight_layout()

    if(save):
        outname = save_dir + 'goes_cross_channels_' + \
            GOES_dict['dt_dates'][time_idx].strftime('%Y%m%d%H%M') + \
            '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

def read_GOES_time_series_auto(begin_date, end_date, \
        channels = [2, 6, 13, 8, 9, 10], save_dir = './', \
        dlat = [40.750520, \
                 40.672445,\
                 40.624068, \
                 40.575759, \
                 40.546941, \
                 40.487687, \
                 40.397445, \
                 40.339418], \
        dlon = [-121.040965, \
                 -120.906663, \
                 -120.812962, \
                 -120.719259, \
                 -120.659137, \
                 -120.588284, \
                 -120.526473, \
                 -120.426204]):

    # Convert the input date_str to datetime
    # --------------------------------------
    begin_dt_date = datetime.strptime(begin_date,"%Y%m%d%H%M")
    end_dt_date   = datetime.strptime(end_date,"%Y%m%d%H%M")

    # Find all downloaded GOES filenames that are between these
    # two dates
    # ---------------------------------------------------------
    all_files = np.array(glob('/home/bsorenson/data/GOES/goes17_abi/*.nc'))
    all_dates = np.array([datetime.strptime(ffile.strip().split('/')[-1][27:40],\
        '%Y%j%H%M%S') for ffile in all_files])
    all_date_strs = np.array([tdate.strftime('%Y%m%d%H%M') for tdate in \
        all_dates])

    # Get rid of 202107202126 because it doesn't work for some reason
    # ----------------------------------------------------------------
    all_files = list(all_files[all_date_strs != '202107202126']) 
    all_dates = list(all_dates[all_date_strs != '202107202126']) 
   
    # Get just the file dates, with only 1 date per channel
    # ----------------------------------------------------- 
    unique_dates = sorted(list(set(all_dates)))

    in_all = np.array([((tdate > begin_dt_date) & (tdate < end_dt_date)) \
        for tdate in all_dates])
    in_unique = np.array([((tdate > begin_dt_date) & (tdate < end_dt_date)) \
        for tdate in unique_dates])
    
    in_all_idx = np.where(in_all == True)
    in_unique_idx = np.where(in_unique == True)

    if(len(in_all_idx[0]) == 0):
        print("ERROR: no data between",begin_date, end_date)
        print("     Run auto_GOES_download to get the data")
        return

    # good_all_files contains the filenames for all files (including
    # all channels) that have timestamps within the desired range
    #
    # good_unique_times contains just a single time for all channel
    # files that are within the range
    good_all_files    = np.array(all_files)[in_all_idx]
    good_all_channels = np.array([int(\
        tfile.strip().split('/')[-1][19:21]) for tfile in good_all_files])
    good_unique_times = np.array(unique_dates)[in_unique_idx[0]]
 
    # Check to make sure that there are channel files for all
    # desired times.
    # -------------------------------------------------------
    for channel in channels:
        if(not channel in good_all_channels):
            print("ERROR: data for ",channel," not downloaded")
            return
    ##!#if(not ch1 in good_all_channels):
    ##!#    print("ERROR: data for ",ch1," not downloaded")
    ##!#    return
    ##!#if(not ch2 in good_all_channels):
    ##!#    print("ERROR: data for ",ch2," not downloaded")
    ##!#    return
    ##!#if(not ch3 in good_all_channels):
    ##!#    print("ERROR: data for ",ch3," not downloaded")
    ##!#    return
    ##!#if(not ch4 in good_all_channels):
    ##!#    print("ERROR: data for ",ch4," not downloaded")
    ##!#    return
    ##!#if(not ch5 in good_all_channels):
    ##!#    print("ERROR: data for ",ch5," not downloaded")
    ##!#    return
    ##!#if(not ch6 in good_all_channels):
    ##!#    print("ERROR: data for ",ch6," not downloaded")
    ##!#    return

    print("All data present")
    
    # Set up arrays to hold the data
    # ------------------------------
    channels = np.array(channels)
    #channels = np.array([ch1, ch2, ch3, ch4, ch5, ch6])

    # Idx 1 : "Clear" pixel
    # Idx 2 : coldest IR pixel
    # Idx 3 : hit by WV channel cold bubbles
    #dlat = [41.20980,40.7595, 41.2456]
    #dlon = [-120.9810, -120.6530, -120.3877]
    ##!#dlat = [40.750520, \
    ##!#         40.672445,\
    ##!#         40.624068, \
    ##!#         40.575759, \
    ##!#         40.546941, \
    ##!#         40.487687, \
    ##!#         40.397445, \
    ##!#         40.339418]
    ##!#dlon = [-121.040965, \
    ##!#         -120.906663, \
    ##!#         -120.812962, \
    ##!#         -120.719259, \
    ##!#         -120.659137, \
    ##!#         -120.588284, \
    ##!#         -120.526473, \
    ##!#         -120.426204]

    goes_data = np.full((len(good_unique_times), len(channels), \
        len(dlat)), \
        np.nan)
    goes_lats = np.full((len(good_unique_times), len(channels), \
        len(dlat)), \
        np.nan)
    goes_lons = np.full((len(good_unique_times), len(channels), \
        len(dlat)), \
        np.nan)

    for ii, ttime in enumerate(good_unique_times):
        print(ttime.strftime('%Y%m%d%H%M'))
        date_str = ttime.strftime('%Y%m%d%H%M')
        #if(date_str == '202107202126'):
        #    print("Not making image for this time")
        #else:
        for jj, tch in enumerate(channels):
            # Extract the GOES values for the current time
            # --------------------------------------------
            goes_vals, goes_lats_local, goes_lons_local  = \
                get_GOES_data_lat_lon(date_str, dlat, dlon, tch)
            goes_data[ii,jj,:] = goes_vals / 1.
            goes_lats[ii,jj,:] = goes_lats_local / 1.
            goes_lons[ii,jj,:] = goes_lons_local / 1.

    # Put the data in a dictionary
    # ----------------------------
    out_dict = {}
    out_dict['data'] = goes_data
    out_dict['goes_lats'] = goes_lats
    out_dict['goes_lons'] = goes_lons
    out_dict['dt_dates'] = good_unique_times
    out_dict['channels'] = channels
    out_dict['plat'] = dlat
    out_dict['plon'] = dlon

    return out_dict

def plot_GOES_satpy_point_test(date_str, channel = 2, \
        plats = [40.750520, \
                 40.672445,\
                 40.624068, \
                 40.575759, \
                 40.546941, \
                 40.487687, \
                 40.397445, \
                 40.339418], \
        plons = [-121.040965, \
                 -120.906663, \
                 -120.812962, \
                 -120.719259, \
                 -120.659137, \
                 -120.588284, \
                 -120.526473, \
                 -120.426204], \
        plats_2 = None, plons_2 = None, \
        save = False):

    # Read GOES image data
    # --------------------
    #print(GOES_dict['dt_dates'][time_idx])
    #date_str = GOES_dict['dt_dates'][time_idx].strftime('%Y%m%d%H%M')
    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
    var, crs, lons, lats, lat_lims, lon_lims, plabel = \
        read_GOES_satpy(date_str, channel)

    plt.close('all')
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1, projection = crs)

    labelsize = 10
    font_size = 8
    plot_GOES_satpy(date_str, channel, \
        ax = ax1, var = var, crs = crs, \
        lons = lons, lats = lats, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = goes_channel_dict[str(channel)]['limits'][0], \
        vmax = goes_channel_dict[str(channel)]['limits'][1], \
        ptitle = '', plabel = plabel, \
        #vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)

    # Add the locations of the points
    # -------------------------------
    for tlat, tlon in zip(plats, plons):
        print(tlat, tlon)
        ax1.plot(tlon, tlat, linewidth=2, markersize = 8, marker='.',
                 color = 'black', transform=datacrs)
        ax1.plot(tlon, tlat, linewidth=2, markersize = 5, marker='.',
                 transform=datacrs)
    if(plats_2 is not None):
        # Add the locations of the second points, if desired
        # --------------------------------------------------
        for tlat, tlon in zip(plats_2, plons_2):
            print(tlat, tlon)
            ax1.plot(tlon, tlat, linewidth=2, markersize = 5, marker='^',
                     color = 'black', transform=datacrs)
            ax1.plot(tlon, tlat, linewidth=2, markersize = 3, marker='^',
                     transform=datacrs)
        
    plot_figure_text(ax1, 'GOES-17 ' + \
        str(goes_channel_dict[str(channel)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    ax1.set_title(dt_date_str.strftime(\
        '%Y-%m-%d %H:%MZ'), fontsize = 10)

    fig.tight_layout()

    if(save):
        print("SAVE HERE")
    else:
        plt.show()

def plot_GOES_satpy_hatch(date_str, save = False):

    # Read GOES image data
    # --------------------
    #print(GOES_dict['dt_dates'][time_idx])
    #date_str = GOES_dict['dt_dates'][time_idx].strftime('%Y%m%d%H%M')
    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
    var1, crs1, lons1, lats1, lat_lims1, lon_lims1, plabel1 = \
        read_GOES_satpy(date_str, 2)
    var2, crs2, lons2, lats2, lat_lims2, lon_lims2, plabel2 = \
        read_GOES_satpy(date_str, 13)

    plt.close('all')
    fig = plt.figure()
    ax1 = fig.add_subplot(1,2,1, projection = crs1)
    ax2 = fig.add_subplot(1,2,2, projection = crs2)

    labelsize = 10
    font_size = 8
    plot_GOES_satpy(date_str, 2, \
        ax = ax1, var = var1, crs = crs1, \
        lons = lons1, lats = lats1, lat_lims = lat_lims1, \
        lon_lims = lon_lims1, \
        vmin = goes_channel_dict[str(2)]['limits'][0], \
        vmax = goes_channel_dict[str(2)]['limits'][1], \
        ptitle = '', plabel = plabel1, \
        #vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,\
        save=False)
    plot_GOES_satpy(date_str, 13, \
        ax = ax2, var = var2, crs = crs2, \
        lons = lons2, lats = lats2, lat_lims = lat_lims2, \
        lon_lims = lon_lims2, \
        vmin = goes_channel_dict[str(13)]['limits'][0], \
        vmax = goes_channel_dict[str(13)]['limits'][1], \
        ptitle = '', plabel = plabel2, \
        #vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,\
        save=False)

    # Play with hatching
    # ------------------
    tmp_data = np.full(var1.shape, np.nan)
    in_smoke = (var1 > 20.)
    tmp_data[in_smoke] = 1.
    tmp_data[~in_smoke] = 0.

    mask_tmp_data = np.ma.masked_invalid(tmp_data)

    hash_data   = np.ma.masked_where(mask_tmp_data != 1, mask_tmp_data)
    #nohash_data = np.ma.masked_where(mask_tmp_data == 1, mask_tmp_data)

    hash0 = ax2.pcolor(lons1, lats1,\
        hash_data, hatch = 'xxx', alpha=0., transform = datacrs, \
        shading = 'auto')


    # Add labels
    # ----------

    plot_figure_text(ax1, 'GOES-17 ' + \
        str(goes_channel_dict[str(2)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    ax1.set_title(dt_date_str.strftime(\
        '%Y-%m-%d %H:%MZ'), fontsize = 10)

    plot_figure_text(ax2, 'GOES-17 ' + \
        str(goes_channel_dict[str(13)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    ax2.set_title(dt_date_str.strftime(\
        '%Y-%m-%d %H:%MZ'), fontsize = 10)

    fig.tight_layout()

    if(save):
        print("SAVE HERE")
    else:
        plt.show()
