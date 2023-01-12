"""
  NAME:

  PURPOSE:
    Read and plot GOES data using satpy. Also contains functions for 
        automatically downloading GOES data from the Amazon bucket. 

    NOTE: must be running a Conda environment with satpy installed. 

  SYNTAX:
    # = = = = = = = = = = =
    # Notes for downloading
    # = = = = = = = = = = =
    Before running the downloading scripts, the main of which is 
        auto_GOES_download(), you MUST have the Amazon AWS CLI 
        interface installed on your Linux terminal. To do this, 
        you must have sudo permissions on your computer. If you
        have sudo permissions, run the following in the command
        line:
    
    $ sudo apt install awscli

    Once this is installed, and if you have a Conda environment
        installed that has satpy and datetime (and the other essential
        packages listed in the "import" section), you may download
        GOES data automatically using auto_GOES_download. To run
        the function, use the following syntax:

    >>> from GOESLib_lite import *
    >>> begin_date = '202107201200' # begin_date in YYYYMMDDHHMM format
    >>> end_date   = '202107210000' # end_date in YYYYMMDDHHMM format
    >>> interval   = 30  # this is the number of minutes between each file
    >>> auto_GOES_download(begin_date, end_date, interval)
  
    # = = = = = = = = = = 
    # Notes for plotting
    # = = = = = = = = = =
    When running in a separate Python script, use the following syntax:
    >>> from GOESLib_lite import *
    >>> date_str = '202107202126'  # date_str in YYYYMMDDHHMM format, for whatever the time of the GOES file 
    >>> plot_GOES_satpy(date_str, 2)   # Plots for channel 2 data
    
    If you want to save the image, do:
    >>> plot_GOES_satpy(date_str, 2, save = True)


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
from matplotlib.lines import Line2D
from metpy.plots import USCOUNTIES
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from satpy import find_files_and_readers
from satpy.scene import Scene
from satpy.writers import get_enhanced_image
from glob import glob
import os

home_dir = os.environ['HOME']

sys.path.append(home_dir + '/')
from python_lib import *

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# Set up global variables
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
datacrs = ccrs.PlateCarree()
data_dir = home_dir + '/data/GOES/goes17_abi/'
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
        'short_name': 'Blue',
        'wavelength': 0.47
    },\
    '2': {
        'limits': [0, 80],\
        'name': 'Red',
        'short_name': 'Red',
        'wavelength': 0.64
    },\
    '3': {
        'limits': [None, None],\
        'name': 'Veggie',
        'short_name': 'Veggie',
        'wavelength': 0.86
    },\
    '4': {
        'limits': [None, None],\
        'name': 'Cirrus',
        'short_name': 'Cirrus',
        'wavelength': 1.38
    },\
    '5': {
        'limits': [None, None],\
        'name': 'Snow/Ice',
        'short_name': 'Snow/Ice',
        'wavelength': 1.61
    },\
    '6': {
        'limits': [0, 50],\
        'name': 'Cloud Particle Size',
        'short_name': 'Cloud Particle Size',
        'wavelength': 2.25
    },\
    '7': {
        'limits': [None, None],\
        'name': 'Shortwave Window',
        'short_name': 'Shortwave Window',
        'wavelength': 3.90
    },\
    '8': {
        'limits': [240, 250],\
        'name': 'Upper-Level Water Vapor',
        'short_name': 'Upper-Level WV',
        'wavelength': 6.18
    },\
    '9': {
        'limits': [250, 260],\
        'name': 'Mid-Level Water Vapor',
        'short_name': 'Mid-Level WV',
        'wavelength': 6.95
    },\
    '10': {
        'limits': [255, 270],\
        'name': 'Lower-Level Water Vapor',
        'short_name': 'Lower-Level WV',
        'wavelength': 7.34
    },\
    '11': {
        'limits': [None, None],\
        'name': 'Cloud-Top Phase',
        'short_name': 'Cloud-Top Phase',
        'wavelength': 8.50
    },\
    '12': {
        'limits': [None, None],\
        'name': 'Ozone',
        'short_name': 'Ozone',
        'wavelength': 9.61
    },\
    '13': {
        'limits': [270, 330],\
        'name': 'Clean IR Longwave Window',
        'short_name': 'Clean TIR',
        'wavelength': 10.35
    },\
    '14': {
        'limits': [None, None],\
        'name': 'IR Longwave Window',
        'short_name': 'TIR',
        'wavelength': 11.20
    },\
    '15': {
        'limits': [None, None],\
        'name': 'Dirty IR Longwave Window',
        'short_name': 'Dirty TIR',
        'wavelength': 12.30
    },\
    '16': {
        'limits': [None, None],\
        'name': 'CO2 Longwave Infrared',
        'short_name': 'CO2 TIR',
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
        ##!#'2110': {
        'asos': 'asos_data_20210713.csv',
        #'goes': home_dir + '/data/GOES/Aqua/MYD021KM.A2021203.2110.061.2021204155922.hdf',
        #'mdswv': home_dir + '/data/GOES/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
        #'ceres': home_dir + '/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072210-2021072221.nc',
        #'airs': [home_dir + '/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
        'Lat': [39.5, 42.0],
        'Lon': [-122.0, -119.5],
        ##!#'goes_Lat': [39.0, 42.5],
        ##!#'goes_Lon': [-123., -119.]
        'goes_Lat': [39.5, 42.0],
        'goes_Lon': [-122.0, -119.5],
        ##!#}
    },
    "2021-07-14": {
        ##!#'2110': {
        'asos': 'asos_data_20210713.csv',
        #'goes': home_dir + '/data/GOES/Aqua/MYD021KM.A2021203.2110.061.2021204155922.hdf',
        #'mdswv': home_dir + '/data/GOES/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
        #'ceres': home_dir + '/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072210-2021072221.nc',
        #'airs': [home_dir + '/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
        'Lat': [39.5, 42.0],
        'Lon': [-122.0, -119.5],
        ##!#'goes_Lat': [39.0, 42.5],
        ##!#'goes_Lon': [-123., -119.]
        'goes_Lat': [39.5, 42.0],
        'goes_Lon': [-122.0, -119.5],
        ##!#}
    },
    "2021-07-20": {
        'asos': 'asos_data_20210720.csv',
        #'goes': home_dir + '/data/GOES/Aqua/MYD021KM.A2021201.2125.061.2021202154814.hdf',
        #'ceres': home_dir + '/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072010-2021072021.nc',
        #'airs': [home_dir + '/data/AIRS/Aqua/AIRS.2021.07.20.214.L2.SUBS2RET.v6.0.32.0.G21202153435.hdf'],
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
        #'goes': home_dir + '/data/GOES/Aqua/MYD021KM.A2021202.2030.061.2021203174050.hdf',
        #'ceres': home_dir + '/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072109-2021072122.nc',
        #'airs': [home_dir + '/data/AIRS/Aqua/AIRS.2021.07.21.205.L2.SUBS2RET.v6.0.32.0.G21203185004.hdf'],
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
        #'goes': home_dir + '/data/GOES/Aqua/MYD021KM.A2021203.2110.061.2021204155922.hdf',
        #'mdswv': home_dir + '/data/GOES/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
        #'ceres': home_dir + '/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072210-2021072221.nc',
        'airs': [home_dir + '/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
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
        'goes': home_dir + '/data/GOES/Aqua/MYD021KM.A2021204.2155.061.2021205153516.hdf',
        #'mdswv': home_dir + '/data/GOES/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
        #'ceres': home_dir + '/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072210-2021072221.nc',
        #'airs': [home_dir + '/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
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
        #'goes': home_dir + '/data/GOES/Aqua/MYD021KM.A2021218.2025.061.2021219151802.hdf',
        #'mdswv': home_dir + '/data/GOES/Aqua/MYD05_L2.A2021218.2025.061.2021219152751.hdf',
        'airs': [home_dir + '/data/AIRS/Aqua/AIRS.2021.08.04.206.L2.SUBS2RET.v6.0.32.0.G21217152448.hdf',\
                 home_dir + '/data/AIRS/Aqua/AIRS.2021.08.04.207.L2.SUBS2RET.v6.0.32.0.G21217152904.hdf'],\
        #'omi': home_dir + '/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
        'Lat': [36.0, 39.0],
        'Lon': [-118.0, -114.0],
        'goes_Lat': [35.0, 40.0],
        'goes_Lon': [-119., -113.]
    },
    "2021-08-05": {
        '2120': {
            'asos': 'asos_data_20210805.csv',
            'goes': home_dir + '/data/GOES/Aqua/MYD021KM.A2021217.2120.061.2021218164201.hdf',
            'mdswv': home_dir + '/data/GOES/Aqua/MYD05_L2.A2021217.2120.061.2021218165546.hdf',
            'ceres': home_dir + '/data/CERES/FLASHFlux/Aqua/FLASH_SSF_Aqua_Version4A_Subset_2021080508-2021080521.nc',
            'airs': [home_dir + '/data/AIRS/Aqua/AIRS.2021.08.05.214.L2.SUBS2RET.v6.0.32.0.G21218175548.hdf'],
            'omi': home_dir + '/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0805t2038-o90733_v003-2021m0807t014855.he5',
            'Lat': [36.5, 39.0],
            'Lon': [-118.0, -114.0],
            'goes_Lat': [36.0, 39.0],
            'goes_Lon': [-118., -114.]
        },
        '2125': {
            'asos': 'asos_california_20210805.csv',
            'goes': home_dir + '/data/GOES/Aqua/MYD021KM.A2021217.2125.061.2021218161010.hdf',
            'ceres': home_dir + '/data/CERES/FLASHFlux/Aqua/FLASH_SSF_Aqua_Version4A_Subset_2021080508-2021080521.nc',
            'airs': [home_dir + '/data/AIRS/Aqua/AIRS.2021.08.05.214.L2.SUBS2RET.v6.0.32.0.G21218175548.hdf'],
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
            'goes': home_dir + '/data/GOES/Aqua/MYD021KM.A2021218.2025.061.2021219151802.hdf',
            'mdswv': home_dir + '/data/GOES/Aqua/MYD05_L2.A2021218.2025.061.2021219152751.hdf',
            'ceres': home_dir + '/data/CERES/FLASHFlux/Aqua/FLASH_SSF_Aqua_Version4A_Subset_2021080609-2021080620.nc',
            'airs': [home_dir + '/data/AIRS/Aqua/AIRS.2021.08.06.204.L2.SUBS2RET.v6.0.32.0.G21219130523.hdf',\
                     home_dir + '/data/AIRS/Aqua/AIRS.2021.08.06.205.L2.SUBS2RET.v6.0.32.0.G21219130455.hdf'],
            'omi': home_dir + '/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
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
            'goes': home_dir + '/data/GOES/Aqua/MYD021KM.A2021219.2110.061.2021220151612.hdf',
            #'mdswv': home_dir + '/data/GOES/Aqua/MYD05_L2.A2021218.2025.061.2021219152751.hdf',
            'airs': [home_dir + '/data/AIRS/Aqua/AIRS.2021.08.07.212.L2.SUBS2RET.v6.0.32.0.G21220123225.hdf'],\
            #'omi': home_dir + '/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
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
            #'goes': home_dir + '/data/GOES/Aqua/MYD021KM.A2021218.2025.061.2021219151802.hdf',
            #'mdswv': home_dir + '/data/GOES/Aqua/MYD05_L2.A2021218.2025.061.2021219152751.hdf',
            'airs': [home_dir + '/data/AIRS/Aqua/AIRS.2021.08.08.202.L2.SUBS2RET.v6.0.32.0.G21221124420.hdf',\
                     home_dir + '/data/AIRS/Aqua/AIRS.2021.08.08.203.L2.SUBS2RET.v6.0.32.0.G21221185932.hdf'],\
            #'omi': home_dir + '/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
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
            'goes': home_dir + '/data/GOES/Aqua/MYD021KM.A2021242.2115.061.2021243183953.hdf',
            'airs': [home_dir + '/data/AIRS/Aqua/AIRS.2021.08.30.213.L2.SUBS2RET.v6.0.32.0.G21243151114.hdf'],
            'Lat': [38.0, 40.0],
            'Lon': [-121.0, -118.5]
        }
    },
    "2021-09-01": {
        '2105': {
            'asos': 'asos_data_20210830.csv',
            'goes': home_dir + '/data/GOES/Aqua/MYD021KM.A2021244.2105.061.2021245152256.hdf',
            #'airs': [home_dir + '/data/AIRS/Aqua/AIRS.2021.08.30.213.L2.SUBS2RET.v6.0.32.0.G21243151114.hdf'],
            'Lat': [38.0, 42.0],
            'Lon': [-121.5, -118.0],
            'goes_Lat': [38.0, 42.0],
            'goes_Lon': [-121.5, -118.]
        }
    } 
}

min_dict = {
    2: 5,
    6: 5,
    8: 240, 
    9: 250, 
    10: 255, 
    13: 270,
}
max_dict = {
    2: 80,
    6: 30, 
    8: 250, 
    9: 260, 
    10: 270, 
    13: 330,
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
        home_dir + '/data/GOES/goes17_abi/'):

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
        # Determine if the file is already downloaded
        # -------------------------------------------
        total_fname = dest_dir + gf
        if(os.path.exists(total_fname)):
            print("File",gf,"already exists. Not downloading")
        else:
            print("File",gf,"not found. Downloading")

            request_add = dt_date_str.strftime('s3://noaa-'+sat+\
                '/ABI-L1b-RadC/%Y/%j/%H/' + gf)
            #print(request_add)
            cmnd_list = ['aws','s3','cp','--no-sign-request',\
                request_add, dest_dir]
            print(' '.join(cmnd_list))
            subprocess.run(cmnd_list, check=True, \
                capture_output=False)

# begin_date : YYYYMMDDHHMM
# end_date   : YYYYMMDDHHMM
# interval   : minutes between desired images
# sat        : default is 'goes17', but can set to 'goes16' or 
#              'goes18' to download other data.
# channels   : a list containing the channel numbers you want
#              to download. For example, [2] would download
#              just data for the 0.64 micron channel (channel 2)
#              [2,6,13] downloads data for the 0.64 micron channel,
#              the 2.25 micron channel, and the 10.35 micron channel.
# dest_dir   : the location where files are going to be downloaded to.
#              NOTE: This MUST be changed when running on a computer
#              that does not belong to bsorenson.
def auto_GOES_download(begin_date, end_date, interval, sat = 'goes17', \
        channels = [2,6,8,9,10,13]):

    dest_dir = home_dir + '/data/GOES/' + sat + '_abi/'

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

def test_thermal_remove(date_str, colorbar = True):

    # Read the GOES data using the two calibrations
    # ---------------------------------------------
    var3_rad, crs3, lons3, lats3, lat_lims3, lon_lims3, plabel3 = \
        read_GOES_satpy(date_str, 7, calibration = 'radiance')
    var3_btmp,   _,     _,     _,         _,         _, _ = \
        read_GOES_satpy(date_str, 7, calibration = 'brightness_temperature')
    var10_btmp,   _,     _,     _,         _,         _, _ = \
        read_GOES_satpy(date_str, 13, calibration = 'brightness_temperature')

    # Convert the brightness temperatures to radiances
    # ------------------------------------------------
    var3_calc_rad = convert_temp_to_radiance(3.9, np.array(var3_btmp))
    var10_calc_rad = convert_temp_to_radiance(10.35, np.array(var3_btmp))

    var3_minus_bght_rad = var3_rad - var3_calc_rad

    # Set up a figure
    # ---------------
    plt.close('all')
    fig = plt.figure(figsize = (8, 8))
    ax1 = fig.add_subplot(2,2,1, projection = crs3)
    ax2 = fig.add_subplot(2,2,2, projection = crs3)
    ax3 = fig.add_subplot(2,2,3, projection = crs3)
    ax4 = fig.add_subplot(2,2,4, projection = crs3)

    ##!#im1 = ax1.pcolormesh(lons3, lats3, var3_rad, transform = datacrs, \
    ##!#    vmin = None, vmax = 3, \
    ##!#    cmap = 'plasma', \
    ##!#    shading = 'auto')
    im1 = ax1.pcolormesh(lons3, lats3, var10_btmp, transform = datacrs, \
        vmin = None, vmax = None, \
        cmap = 'plasma', \
        shading = 'auto')
    ax1.add_feature(cfeature.STATES)
    im2 = ax2.pcolormesh(lons3, lats3, var3_btmp, transform = datacrs, \
        vmin = None, vmax = None, \
        #vmin = None, vmax = 330, \
        cmap = 'plasma', \
        shading = 'auto')
    ax2.add_feature(cfeature.STATES)
    im3 = ax3.pcolormesh(lons3, lats3, var10_calc_rad, transform = datacrs, \
    #im3 = ax3.pcolormesh(lons3, lats3, var3_calc_rad, transform = datacrs, \
        vmin = None, vmax = None, \
        #vmin = None, vmax = 2, \
        cmap = 'plasma', \
        shading = 'auto')
    ax3.add_feature(cfeature.STATES)
    im4 = ax4.pcolormesh(lons3, lats3, var3_calc_rad, transform = datacrs, \
    #im4 = ax4.pcolormesh(lons3, lats3, var3_minus_bght_rad, transform = datacrs, \
        vmin = None, vmax = None, \
        #vmin = None, vmax = 1, \
        cmap = 'plasma', \
        shading = 'auto')
    ax4.add_feature(cfeature.STATES)
    if(colorbar):
        labelsize = 10
        cbar1 = plt.colorbar(im1, ax = ax1, pad = 0.03, fraction = 0.052, \
            extend = 'both')
        cbar1.set_label('10.35 micron bright temp', size = labelsize, weight = 'bold')
        #cbar1.set_label('Radiance', size = labelsize, weight = 'bold')
        cbar2 = plt.colorbar(im2, ax = ax2, pad = 0.03, fraction = 0.052, \
            extend = 'both')
        cbar2.set_label('3.9 micron bright temp', size = labelsize, weight = 'bold')
        #cbar2.set_label('Brightness Temperature', size = labelsize, weight = 'bold')
        cbar3 = plt.colorbar(im3, ax = ax3, pad = 0.03, fraction = 0.052, \
            extend = 'both')
        cbar3.set_label('10.35 micron radiance', size = labelsize, weight = 'bold')
        #cbar3.set_label('Radiance from Bright Temp', size = labelsize, weight = 'bold')
        cbar4 = plt.colorbar(im4, ax = ax4, pad = 0.03, fraction = 0.052, \
            extend = 'both')
        cbar4.set_label('3.9 micron radiance', size = labelsize, weight = 'bold')
        #cbar4.set_label('Radiance - Bright Rad', size = labelsize, weight = 'bold')

    # Set up another scatter plot
    # ---------------------------
    fig2 = plt.figure()
    axs1 = fig2.add_subplot(1,1,1)
    #scatter_data = [var3_rad.flatten(), var3_calc_rad.flatten(), \
    #    var3_minus_bght_rad.flatten()]
    axs1.scatter(var3_btmp.flatten(), var10_btmp.flatten(), s = 6)
    #axs1.scatter(var3_btmp.flatten(), var10_btmp.flatten(), s = 6)
    #axs1.boxplot(scatter_data)


    plt.show()

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

# phi  = latitude of the point
# lbda = longitude of the point
def calc_zenith_angle(phi, lbda, sat = 'goes17'):
    
    if(sat == 'goes17'):
        # GOES-17 values
        sat_height = 35786.0234375 # km
        sat_lon    = -137.2 # goes17
    elif(sat == 'goes16'):
        # GOES-16 values
        sat_height = 35786.023 # km
        sat_lon    = -75.2 # goes17


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

def select_plume_data(date_str):
    var_ch13, _, lons_ch13, lats_ch13, _, _, _ = \
        read_GOES_satpy(date_str, 13)

    ##!## Prep the data for co-locating
    ##!## -----------------------------
    ##!#flat_var_ch13 = var_ch13.flatten()
    ##!#flat_lats_ch13 = lats_ch13.flatten()
    ##!#flat_lons_ch13 = lons_ch13.flatten()

    # Co-locate the ch2 data with the ch13 data
    # -----------------------------------------
    flat_var_ch2, flat_lats_ch2, flat_lons_ch2  = \
        get_GOES_data_lat_lon(date_str, lats_ch13, \
        lons_ch13, 2, version = 1, verbose = True)
    flat_var_ch6, flat_lats_ch6, flat_lons_ch6  = \
        get_GOES_data_lat_lon(date_str, lats_ch13, \
        lons_ch13, 6, version = 1, verbose = True)

    return flat_var_ch2, flat_lats_ch2, flat_lons_ch2, \
           flat_var_ch6, flat_lats_ch6, flat_lons_ch6, \
           flat_var_ch13, flat_lats_ch13, flat_lons_ch13

def screen_plume_data(flat_var_ch2, flat_lats_ch2, flat_lons_ch2, \
           flat_var_ch6, flat_lats_ch6, flat_lons_ch6, \
           flat_var_ch13, flat_lats_ch13, flat_lons_ch13):

    screen_limit = 0.05
    max_ch6 = 10.
    max_ch = 350.
    test_data = flat_var_ch2 - flat_var_ch6
    #test_data = GOES_ch1['data'] - GOES_ch5['data']
    hash_data   = np.ma.masked_where(test_data <  \
        (np.nanmean(test_data) * screen_limit), flat_var_ch2)
    hash_data = np.ma.masked_where(flat_var_ch6 < max_ch6, hash_data)
    nohash_data = np.ma.masked_where(test_data >= \
        (np.nanmean(test_data) * screen_limit),flat_var_ch2) 

    hash_data = np.ma.masked_where(var_ch13 > max_ch, hash_data)
    nohash_data = np.ma.masked_where(var_ch13 > max_ch, nohash_data)

    return hash_data, nohash_data

# Determine areas of an image that are in smoke, defined by:
#    (ch1_refl - ch5_refl) < mean(ch1_refl - ch5_refl) * 0.25
def find_plume_GOES(date_str):
    ##!###!#var_ch2, _, lons_ch2, lats_ch2, _, _, _ = \
    ##!###!#    read_GOES_satpy(date_str, 2)
    ##!###!#var_ch6, _, lons_ch6, lats_ch6, _, _, _ = \
    ##!###!#    read_GOES_satpy(date_str, 6)
    var_ch13, _, lons_ch13, lats_ch13, _, _, _ = \
        read_GOES_satpy(date_str, 13)

    ##!## Prep the data for co-locating
    ##!## -----------------------------
    ##!###!#flat_var_ch2 = var_ch2.flatten()
    ##!###!#flat_lats_ch2 = lats_ch2.flatten()
    ##!###!#flat_lons_ch2 = lons_ch2.flatten()
    ##!###!#flat_var_ch6 = var_ch6.flatten()
    ##!###!#flat_lats_ch6 = lats_ch6.flatten()
    ##!###!#flat_lons_ch6 = lons_ch6.flatten()
    ##!#flat_var_ch13 = var_ch13.flatten()
    ##!#flat_lats_ch13 = lats_ch13.flatten()
    ##!#flat_lons_ch13 = lons_ch13.flatten()

    # Co-locate the ch2 data with the ch13 data
    # -----------------------------------------
    var_ch2, lats_ch2, lons_ch2  = \
        get_GOES_data_lat_lon(date_str, lats_ch13, \
        lons_ch13, 2, version = 1, verbose = True)
    var_ch6, lats_ch6, lons_ch6  = \
        get_GOES_data_lat_lon(date_str, lats_ch13, \
        lons_ch13, 6, version = 1, verbose = True)

    print(var_ch2.shape)

    ##!#flat_var_ch2, flat_lats_ch2, flat_lons_ch2, \
    ##!#flat_var_ch6, flat_lats_ch6, flat_lons_ch6, \
    ##!#flat_var_ch13, flat_lats_ch13, flat_lons_ch13 = \
    ##!#    select_plume_data(date_str)

    screen_limit = 0.05
    max_ch6 = 10.
    max_ch = 350.
    test_data = var_ch2 - var_ch6
    #test_data = GOES_ch1['data'] - GOES_ch5['data']
    hash_data   = np.ma.masked_where(test_data <  \
        (np.nanmean(test_data) * screen_limit), var_ch2)
    hash_data = np.ma.masked_where(var_ch6 < max_ch6, hash_data)
    nohash_data = np.ma.masked_where(test_data >= \
        (np.nanmean(test_data) * screen_limit),var_ch2) 

    hash_data = np.ma.masked_where(var_ch13 > max_ch, hash_data)
    nohash_data = np.ma.masked_where(var_ch13 > max_ch, nohash_data)

    ##!#hash_data, nohash_data = \
    ##!#    screen_plume_data(flat_var_ch2, flat_lats_ch2, flat_lons_ch2, \
    ##!#           flat_var_ch6, flat_lats_ch6, flat_lons_ch6, \
    ##!#           flat_var_ch13, flat_lats_ch13, flat_lons_ch13)

    ##!#screen_limit = 0.05
    ##!#max_ch6 = 10.
    ##!#max_ch = 350.
    ##!#test_data = flat_var_ch2 - flat_var_ch6
    ##!##test_data = GOES_ch1['data'] - GOES_ch5['data']
    ##!#hash_data   = np.ma.masked_where(test_data <  \
    ##!#    (np.nanmean(test_data) * screen_limit), flat_var_ch2)
    ##!#hash_data = np.ma.masked_where(flat_var_ch6 < max_ch6, hash_data)
    ##!#nohash_data = np.ma.masked_where(test_data >= \
    ##!#    (np.nanmean(test_data) * screen_limit),flat_var_ch2) 

    ##!#hash_data = np.ma.masked_where(var_ch13 > max_ch, hash_data)
    ##!#nohash_data = np.ma.masked_where(var_ch13 > max_ch, nohash_data)

    return hash_data, var_ch2, var_ch6

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Plotting functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

calib_dict = {
    0.47: 'reflectance',
    0.64: 'reflectance',
    0.86: 'reflectance',
    1.61: 'reflectance',
    2.25: 'reflectance',
    3.90: 'brightness_temperature',
    #3.90: 'radiance',
    6.18: 'brightness_temperature',
    6.95: 'brightness_temperature',
    7.34: 'brightness_temperature',
    10.35: 'brightness_temperature'
}   
label_dict = {
    0.47: '%',
    0.64: '%',
    0.86: '%',
    1.61: '%',
    2.25: '%',
    3.90: 'K',
    #3.90: 'mW m$^{-2}$ Sr$^{-1}$ (cm$^{-1}$)$^{-1}$',
    #6.18: 'mW m$^{-2}$ Sr$^{-1}$ (cm$^{-1}$)$^{-1}$',
    6.18: 'K',
    6.95: 'K',
    7.34: 'K',
    10.35: 'K'
}   
cmap_dict = {
    0.47:  'Greys_r',
    0.64:  'Greys_r',
    0.86:  'Greys_r',
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
def read_GOES_satpy(date_str, channel, scene_date = None, \
        calibration = None, sat = 'goes17', zoom = True):

    data_dir = home_dir + '/data/GOES/' + sat + '_abi/'

    # Extract the channel wavelength using the input string
    # -----------------------------------------------------
    if(channel != 'true_color'):
        channel = goes_channel_dict[str(channel)]['wavelength']

    # Determine the correct GOES files associated with the date
    # ---------------------------------------------------------
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    dt_date_str_end = dt_date_str + timedelta(minutes = 3)

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

    print('channel = ',channel)

    # Load the desired channel data
    # -----------------------------
    if(channel != 'true_color'):
        if(calibration == None):
            calibration = calib_dict[channel]
        scn.load([channel], calibration = calibration)

        # Zoom the image on the desired area
        # ----------------------------------
        if(zoom):
            try:
                scn = scn.crop(ll_bbox = (lon_lims[0] + 0.65, lat_lims[0], \
                    lon_lims[1] - 0.65, lat_lims[1]))
            except NotImplementedError:
                print("WARNING: zoom didn't work for dtg",date_str,\
                    " and channel", channel)
                scn = scn

    else:
        scn.load([0.47, 0.64, 0.86])
        scn = scn.resample(scn.min_area(), resampler = 'native')
        scn.load(['true_color'])

        try:
            scn = scn.crop(ll_bbox = (lon_lims[0] + 0.65, lat_lims[0], \
                lon_lims[1] - 0.65, lat_lims[1]))
        except NotImplementedError:
            print("WARNING: zoom didn't work for dtg",date_str,\
                " and channel", channel)


        var = get_enhanced_image(scn[channel]).data
        var = var.transpose('y','x','bands')

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


    # Extract the map projection from the data for plotting
    # -----------------------------------------------------
    crs = scn[channel].attrs['area'].to_cartopy_crs()

    # Extract the lats, lons, and data
    # -----------------------------------------------------
    lons, lats = scn[channel].attrs['area'].get_lonlats()
    if(channel != 'true_color'):
        var = scn[channel].data

        # Extract the appropriate units
        # -----------------------------
        units = label_dict[channel]
        #units = scn[channel].units
        plabel = calib_dict[channel].title() + ' [' + units + ']'
    else:
        units = ''
        plabel = ''

    del scn 

    return var, crs, lons, lats, lat_lims, lon_lims, plabel
    #return var, crs, lat_lims, lon_lims

# channel must be an integer between 1 and 16
def plot_GOES_satpy(date_str, channel, ax = None, var = None, crs = None, \
        lons = None, lats = None, lat_lims = None, lon_lims = None, \
        vmin = None, vmax = None, ptitle = None, plabel = None, \
        calibration = None, labelsize = 10, \
        colorbar = True, counties = False, zoom=True,\
        save=False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    if(var is None): 
        var, crs, lons, lats, lat_lims, lon_lims, plabel = \
            read_GOES_satpy(date_str, channel, calibration = calibration)

    # Plot the GOES data
    # ------------------
    in_ax = True 
    if(ax is None):
        in_ax = False
        plt.close('all')
        ax = plt.axes(projection=crs)

    if(channel == 'true_color'):
        ax.imshow(var.data, transform = crs, extent=(var.x[0], var.x[-1], \
            var.y[-1], var.y[0]), vmin = vmin, vmax = vmax, origin='upper')
    else:
        #im1 = ax.imshow(var, transform = crs, vmin = vmin, vmax = vmax, \
        im1 = ax.pcolormesh(lons, lats, var, transform = datacrs, \
            vmin = vmin, vmax = vmax, \
            cmap = cmap_dict[goes_channel_dict[str(channel)]['wavelength']], \
            shading = 'auto')
        if(colorbar):
            cbar = plt.colorbar(im1, ax = ax, pad = 0.03, fraction = 0.052, \
                extend = 'both')
            cbar.set_label(plabel.replace('_',' '), size = labelsize, weight = 'bold')
    ax.add_feature(cfeature.STATES)

    if(counties):
        ax.add_feature(USCOUNTIES.with_scale('5m'))    

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
        if(channel == 'true_color'):
            ax.set_title('GOES-17 True Color\n' + \
                dt_date_str.strftime('%Y-%m-%d %H:%M UTC'))
        else:
            ax.set_title('GOES-17 ' + \
                str(np.round(goes_channel_dict[str(channel)]['wavelength'],2)) \
                + ' μm\n' + dt_date_str.strftime('%Y-%m-%d %H:%M UTC'))
    else:
        ax.set_title(ptitle)

    if(not in_ax): 
        if(save):
            outname = 'goes_ch' + str(channel) +'_'+ date_str + zoom_add + '.png'
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

# This function gets the GOES data from a higher-resolution channel
# and co-locates it with data from a lower-resolution channel.
def get_GOES_data_lat_lon(date_str, dlat, dlon, channel, version = 0,\
        verbose = False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    var0, _, lons0, lats0, _, _, _ = read_GOES_satpy(date_str, channel)

    #cd_idx = nearest_gridpoint(dlat, dlon,lats0, lons0)
    #goes_val = np.array(var0)[cd_idx]

    #for tdlat, tdlon in zip(dlat, dlon):
    if(version == 0):
        if(not isinstance(dlat, list)):
            dlat = [dlat]
            dlon = [dlon]
        cd_idx = [nearest_gridpoint(tdlat, tdlon, \
            lats0, lons0) for tdlat, tdlon in zip(dlat, dlon)]
        goes_val = np.array([np.array(var0)[cidx] \
            for cidx in cd_idx]).squeeze()
        goes_lat = np.array([np.array(lats0)[cidx] \
            for cidx in cd_idx]).squeeze()
        goes_lon = np.array([np.array(lons0)[cidx] \
            for cidx in cd_idx]).squeeze()
    else:

        arr_var0 = np.array(var0)
        arr_lats0 = np.array(lats0)
        arr_lons0 = np.array(lons0)

        goes_val = np.array([[arr_var0[nearest_gridpoint(\
            tlat1, tlon1, lats0, lons0)] for \
            tlat1, tlon1 in zip(tlat2, tlon2)] for \
            tlat2, tlon2 in zip(dlat, dlon)]).squeeze()
        goes_lat = np.array([[arr_lats0[nearest_gridpoint(\
            tlat1, tlon1, lats0, lons0)] for \
            tlat1, tlon1 in zip(tlat2, tlon2)] for \
            tlat2, tlon2 in zip(dlat, dlon)]).squeeze()
        goes_lon = np.array([[arr_lons0[nearest_gridpoint(\
            tlat1, tlon1, lats0, lons0)] for \
            tlat1, tlon1 in zip(tlat2, tlon2)] for \
            tlat2, tlon2 in zip(dlat, dlon)]).squeeze()

        ##!#goes_val = np.full(len(cd_idx), np.nan)
        ##!#goes_lat = np.full(len(cd_idx), np.nan)
        ##!#goes_lon = np.full(len(cd_idx), np.nan)
       
        ##!#for ii, cidx in enumerate(cd_idx):
        ##!#    if(verbose):
        ##!#        print(ii)
        ##!#    goes_val[ii] = np.array(var0)[cd_idx[ii]][0] 
        ##!#    goes_lat[ii] = np.array(lats0)[cd_idx[ii]][0] 
        ##!#    goes_lon[ii] = np.array(lons0)[cd_idx[ii]][0] 

    return goes_val, goes_lat, goes_lon

# NOTE: This is for playing around with find_plume_GOES
def plot_GOES_scatter(date_str):

    # Read in the base channel, usually Ch13
    # --------------------------------------
    channel1 = 13
    var, crs, lons, lats, lat_lims, lon_lims, plabel = \
        read_GOES_satpy(date_str, channel1)

    # Use find_plume_GOES to both co-locate the high-res
    # visible (and low-res ch6 data) with the ch13 data
    # AND mask where the data are not in the plume
    # 
    # NOTE: As of 6/10/2022, find_plume_GOES does not work
    # as well as find_plume from MODISLib, mainly because
    # the 2.25 micron limit used with the  MODIS data does not
    # apply here.
    # --------------------------------------------------------
    hash_data, low_var2, low_var6 = find_plume_GOES(date_str)

    # Make a plot of something
    # ------------------------
    plt.close('all')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    data = [low_var2.flatten(), low_var6.flatten()]
    ax.boxplot(data)
    plt.show()

def plot_GOES_satpy_6panel(date_str, ch1, ch2, ch3, ch4, ch5, ch6, \
        zoom = True, save_dir = './', sat = 'goes17', save = False):
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    plt.close('all')
    fig1 = plt.figure(figsize = (10,6.5))
    var0, crs0, lons0, lats0, lat_lims, lon_lims, plabel0 = read_GOES_satpy(date_str, ch1, sat = sat)
    var1, crs1, lons1, lats1, lat_lims, lon_lims, plabel1 = read_GOES_satpy(date_str, ch2, sat = sat)
    var2, crs2, lons2, lats2, lat_lims, lon_lims, plabel2 = read_GOES_satpy(date_str, ch3, sat = sat)
    var3, crs3, lons3, lats3, lat_lims, lon_lims, plabel3 = read_GOES_satpy(date_str, ch4, sat = sat)
    var4, crs4, lons4, lats4, lat_lims, lon_lims, plabel4 = read_GOES_satpy(date_str, ch5, sat = sat)
    var5, crs5, lons5, lats5, lat_lims, lon_lims, plabel5 = read_GOES_satpy(date_str, ch6, sat = sat)

    ax0 = fig1.add_subplot(2,3,1, projection = crs0)
    ax1 = fig1.add_subplot(2,3,2, projection = crs1)
    ax2 = fig1.add_subplot(2,3,3, projection = crs2)
    ax3 = fig1.add_subplot(2,3,4, projection = crs3)
    ax4 = fig1.add_subplot(2,3,5, projection = crs4)
    ax5 = fig1.add_subplot(2,3,6, projection = crs5)

    min_dict = {
        2: 5,
        6: 0,
        8: 240, 
        9: 245, 
        10: 250, 
        13: 270,
    }
    max_dict = {
        2: 80,
        6: 40, 
        8: 250, 
        9: 260, 
        10: 270, 
        13: 330,
    }

    ##!#ax1.set_title('GOES-17 Band ' + str(ch2) + '\n' + \
    ##!#    goes_channel_dict[str(ch2)]['name'] + '\n' + \
    labelsize = 11
    font_size = 10
    if(ch1 == 'true_color'):
        plot_GOES_satpy(date_str, ch1, ax = ax0, var = var0, crs = crs0, \
            lons = lons0, lats = lats0, lat_lims = lat_lims, lon_lims = lon_lims, \
            ptitle = '', plabel = plabel0, \
            colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)
    else:
        plot_GOES_satpy(date_str, ch1, ax = ax0, var = var0, crs = crs0, \
            lons = lons0, lats = lats0, lat_lims = lat_lims, lon_lims = lon_lims, \
            vmin = min_dict[ch1], vmax = max_dict[ch1], ptitle = '', plabel = plabel0, \
            colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)
    plot_GOES_satpy(date_str, ch2, ax = ax1, var = var1, crs = crs0, \
        lons = lons1, lats = lats1, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = min_dict[ch2], vmax = max_dict[ch2], ptitle = '', plabel = plabel1, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)
    plot_GOES_satpy(date_str, ch3, ax = ax2, var = var2, crs = crs0, \
        lons = lons2, lats = lats2, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = min_dict[ch3], vmax = max_dict[ch3], ptitle = '', plabel = plabel2, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)
    plot_GOES_satpy(date_str, ch4, ax = ax3, var = var3, crs = crs0, \
        lons = lons3, lats = lats3, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = min_dict[ch4], vmax = max_dict[ch4], ptitle = '', plabel = plabel3, \
        colorbar = True, labelsize = labelsize, zoom=True,save=False)
    plot_GOES_satpy(date_str, ch5, ax = ax4, var = var4, crs = crs0, \
        lons = lons4, lats = lats4, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = min_dict[ch5], vmax = max_dict[ch5], ptitle = '', plabel = plabel4, \
        colorbar = True, labelsize = labelsize, zoom=True,save=False)
    plot_GOES_satpy(date_str, ch6, ax = ax5, var = var5, crs = crs0, \
        lons = lons5, lats = lats5, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = min_dict[ch6], vmax = max_dict[ch6], ptitle = '', plabel = plabel5, \
        colorbar = True, labelsize = labelsize, zoom=True,save=False)


    smoke_lat = 40.5672
    smoke_lon = -120.9731
    clear_lat1 = 40.9128
    clear_lon1 = -121.3236

    point_size = 5
    plot_point_on_map(ax2, smoke_lat, smoke_lon, markersize = point_size)
    plot_point_on_map(ax2, clear_lat1, clear_lon1, markersize = point_size)
    sw_idx_s = nearest_gridpoint(smoke_lat, smoke_lon, lats2, lons2)
    sw_idx_c = nearest_gridpoint(clear_lat1, clear_lon1, lats2, lons2)
    print("TIR")
    print("     Smoky  - ", np.array(var2)[sw_idx_s])
    print("     Clear1 - ", np.array(var2)[sw_idx_c])

    plot_point_on_map(ax3, smoke_lat, smoke_lon, markersize = point_size)
    plot_point_on_map(ax3, clear_lat1, clear_lon1, markersize = point_size)
    sw_idx_s = nearest_gridpoint(smoke_lat, smoke_lon, lats3, lons3)
    sw_idx_c = nearest_gridpoint(clear_lat1, clear_lon1, lats3, lons3)
    print("Upper WV")
    print("     Smoky  - ", np.array(var3)[sw_idx_s])
    print("     Clear1 - ", np.array(var3)[sw_idx_c])

    plot_point_on_map(ax4, smoke_lat, smoke_lon, markersize = point_size)
    plot_point_on_map(ax4, clear_lat1, clear_lon1, markersize = point_size)
    sw_idx_s = nearest_gridpoint(smoke_lat, smoke_lon, lats4, lons4)
    sw_idx_c = nearest_gridpoint(clear_lat1, clear_lon1, lats4, lons4)
    print("Mid WV")
    print("     Smoky  - ", np.array(var4)[sw_idx_s])
    print("     Clear1 - ", np.array(var4)[sw_idx_c])

    plot_point_on_map(ax5, smoke_lat, smoke_lon, markersize = point_size)
    plot_point_on_map(ax5, clear_lat1, clear_lon1, markersize = point_size)
    sw_idx_s = nearest_gridpoint(smoke_lat, smoke_lon, lats5, lons5)
    sw_idx_c = nearest_gridpoint(clear_lat1, clear_lon1, lats5, lons5)
    print("Lower WV")
    print("     Smoky  - ", np.array(var5)[sw_idx_s])
    print("     Clear1 - ", np.array(var5)[sw_idx_c])

    if(ch1 == 'true_color'):
        plot_figure_text(ax0, 'True Color', \
            xval = None, yval = None, transform = None, \
            color = 'red', fontsize = font_size, backgroundcolor = 'white', \
            halign = 'right')
    else: 
        plot_figure_text(ax0, \
            str(goes_channel_dict[str(ch1)]['wavelength']) + ' μm', \
            xval = None, yval = None, transform = None, \
            color = 'red', fontsize = font_size, backgroundcolor = 'white', \
            halign = 'right')
        plot_figure_text(ax0, \
            str(goes_channel_dict[str(ch1)]['short_name']), \
            xval = None, yval = None, transform = None, \
            color = 'red', fontsize = font_size - 1, backgroundcolor = 'white', \
            location = 'upper_right', halign = 'right')

    # 2nd channel
    # -----------
    plot_figure_text(ax1, \
        str(goes_channel_dict[str(ch2)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    plot_figure_text(ax1, \
        str(goes_channel_dict[str(ch2)]['short_name']), \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size - 1, backgroundcolor = 'white', \
        location = 'upper_right', halign = 'right')
    # 3rd channel
    # -----------
    plot_figure_text(ax2, \
        str(goes_channel_dict[str(ch3)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    plot_figure_text(ax2, \
        str(goes_channel_dict[str(ch3)]['short_name']), \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size - 1, backgroundcolor = 'white', \
        location = 'upper_right', halign = 'right')
    # 4th channel
    # -----------
    plot_figure_text(ax3, \
        str(goes_channel_dict[str(ch4)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    plot_figure_text(ax3, \
        str(goes_channel_dict[str(ch4)]['short_name']), \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size - 1, backgroundcolor = 'white', \
        location = 'upper_right', halign = 'right')
    # 5th channel
    # -----------
    plot_figure_text(ax4, \
        str(goes_channel_dict[str(ch5)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    plot_figure_text(ax4, \
        str(goes_channel_dict[str(ch5)]['short_name']), \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size - 1, backgroundcolor = 'white', \
        location = 'upper_right', halign = 'right')
    # 6th channel
    # -----------
    plot_figure_text(ax5, \
        str(goes_channel_dict[str(ch6)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    plot_figure_text(ax5, \
        str(goes_channel_dict[str(ch6)]['short_name']), \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size - 1, backgroundcolor = 'white', \
        location = 'upper_right', halign = 'right')

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
   
    ##!## "Coldest" values
    ##!## -------------
    ##!#c_lon_stn = -120.6530
    ##!#c_lat_stn = 40.7595
    ##!#cd_idx = nearest_gridpoint(c_lat_stn, c_lon_stn,\
    ##!#    lats3, lons3)
    ##!## "Cold" values
    ##!## -------------
    ##!#lon_stn = -120.3877
    ##!#lat_stn = 41.2456
    ##!#c_idx = nearest_gridpoint(lat_stn, lon_stn,\
    ##!#    lats3, lons3)
    ##!## "Warm" values
    ##!## -------------
    ##!#w_lon_stn = -120.9810
    ##!#w_lat_stn = 41.20980

    ##!##w_lat_stn = 40.50900
    ##!#w_idx = nearest_gridpoint(w_lat_stn, w_lon_stn,\
    ##!#    lats3, lons3)

    ##!#ax2.plot(c_lon_stn, c_lat_stn,
    ##!#         color='tab:green', linewidth=2, marker='o',
    ##!#         transform=datacrs)
    ##!#ax2.plot(lon_stn, lat_stn,
    ##!#         color='tab:blue', linewidth=2, marker='o',
    ##!#         transform=datacrs)
    ##!#ax2.plot(w_lon_stn, w_lat_stn,
    ##!#         color='tab:purple', linewidth=2, marker='o',
    ##!#         transform=datacrs)
    ##!#ax3.plot(lon_stn, lat_stn,
    ##!#         color='tab:blue', linewidth=2, marker='o',
    ##!#         transform=datacrs)
    ##!#ax3.plot(w_lon_stn, w_lat_stn,
    ##!#         color='tab:purple', linewidth=2, marker='o',
    ##!#         transform=datacrs)
    ##!#ax4.plot(lon_stn, lat_stn,
    ##!#         color='tab:blue', linewidth=2, marker='o',
    ##!#         transform=datacrs)
    ##!#ax4.plot(w_lon_stn, w_lat_stn,
    ##!#         color='tab:purple', linewidth=2, marker='o',
    ##!#         transform=datacrs)
    ##!#ax5.plot(lon_stn, lat_stn,
    ##!#         color='tab:blue', linewidth=2, marker='o',
    ##!#         transform=datacrs)
    ##!#ax5.plot(w_lon_stn, w_lat_stn,
    ##!#         color='tab:purple', linewidth=2, marker='o',
    ##!#         transform=datacrs)

    ##!#print("TIR")
    ##!#print("     Cold - ", np.array(var2)[cd_idx])
    ##!#print("     Warm - ", np.array(var2)[w_idx])
    ##!#print("Upper WV")
    ##!#print("     Cold - ", np.array(var3)[c_idx])
    ##!#print("     Warm - ", np.array(var3)[w_idx])
    ##!#print("Mid   WV")
    ##!#print("     Cold - ", np.array(var4)[c_idx])
    ##!#print("     Warm - ", np.array(var4)[w_idx])
    ##!#print("Lower WV")
    ##!#print("     Cold - ", np.array(var5)[c_idx])
    ##!#print("     Warm - ", np.array(var5)[w_idx])

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

    if(sat == 'goes17'):
        title_str = 'GOES-17\n'
    else:
        title_str = 'GOES-16\n'

    fig1.suptitle(title_str + \
        dt_date_str.strftime('%Y/%m/%d %H:%M UTC'))

    fig1.tight_layout()

    if(save):
        outname = save_dir + sat + '_'+date_str+'_6panel.png'
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
    all_files = glob(home_dir + '/data/GOES/goes17_abi/*.nc')
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
        if(ch1 != 'true_color'):
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

def plot_GOES_figure2_v2(date_str = '202107210000', \
        goes_ch1 = 'true_color', goes_ch2 = 6, goes_ch3 = 13, \
        goes_ch4 = 8, goes_ch5 = 9, goes_ch6 = 10, \
        ch_idx1 = 0, ch_idx2 = 1, ch_idx3 = 2,\
        ttype1 = 'low', ttype2 = 'ml', \
        idx1 = 3, idx2 = 8, idx3 = 5, idx4 = 15, idx5 = 20,\
        date_idx = 25, 
        show_smoke = False, composite = True, double_fig = False, \
        zoom = True, save=False):

    date_str2 = '202107210000'
    dt_date_str2 = datetime.strptime(date_str,"%Y%m%d%H%M")

    # Read the GOES data
    # ------------------------
    var2, crs0, lons, lats, lat_lims2, lon_lims2, plabel2   = \
        read_GOES_satpy(date_str2, goes_ch1)
    var3, crs0, lons3, lats3, lat_lims3, lon_lims3, plabel3 = \
        read_GOES_satpy(date_str2, goes_ch2)
    #var3, crs0, lons2, lats2, lat_lims0, lon_lims0, plabel3 = read_GOES_satpy(date_str2, 6)
    var4, crs0, lons2, lats2, lat_lims0, lon_lims0, plabel4 = \
        read_GOES_satpy(date_str2, goes_ch3)
    var5, crs0, lons2, lats2, lat_lims0, lon_lims2, plabel5 = \
        read_GOES_satpy(date_str2, goes_ch4)
    var6, crs0, lons2, lats2, lat_lims0, lon_lims0, plabel6 = \
        read_GOES_satpy(date_str2, goes_ch5)
    var7, crs0, lons2, lats2, lat_lims0, lon_lims0, plabel7 = \
        read_GOES_satpy(date_str2, goes_ch6)

    # Read the GOES time series data
    # ------------------------------
    file_name1 = home_dir + '/Research/GOES/goes_cross_data_' + \
        ttype1 + '_202107201201_202107210331_v3.nc'
    file_name2 = home_dir + '/Research/GOES/goes_cross_data_' + \
        ttype2 + '_202107201201_202107210331_v3.nc'
    GOES_dict  = read_GOES_time_series_NCDF(file_name1)
    GOES_dict2 = read_GOES_time_series_NCDF(file_name2)

    # ----------------------------------------------------------------------
    #
    #  Set up the 6-panel figure
    #
    # ----------------------------------------------------------------------

    plt.close('all')
    fig   = plt.figure(figsize=(9.5,9))
    gs    = fig.add_gridspec(nrows = 3, ncols = 3)
    ax4   = fig.add_subplot(gs[0,0],  projection = crs0)   # GOES True color
    ax5   = fig.add_subplot(gs[0,1],  projection = crs0)   # GOES SWIR
    ax6   = fig.add_subplot(gs[0,2],  projection = crs0)   # GOES TIR 
    ax7   = fig.add_subplot(gs[1,0],  projection = crs0)   # GOES upper WV
    ax8   = fig.add_subplot(gs[1,1],  projection = crs0)   # GOES midle WV
    ax9   = fig.add_subplot(gs[1,2],  projection = crs0)   # GOES lower WV
    ax10  = fig.add_subplot(gs[2,:])                       # GOES time series
    
    labelsize = 10
    plot_GOES_satpy(date_str, goes_ch1, ax = ax4, var = var2, crs = crs0, \
        lons = lons, lats = lats, lat_lims = lat_lims2, lon_lims = lon_lims2, \
        ptitle = '', plabel = plabel2, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)

    ##!## Plot channel 1, 5, 31, and WV data spatial data
    ##!## -----------------------------------------------
    ##!#plot_GOES_satpy(date_str, 2, ax = ax4, var = var2, crs = crs0, \
    ##!#    lons = lons, lats = lats, lat_lims = lat_lims2, lon_lims = lon_lims2, \
    ##!#    vmin = None, vmax = 60, \
    ##!#    ptitle = '', plabel = plabel2, colorbar = True, labelsize = labelsize, \
    ##!#    zoom=True,save=False)
    plot_GOES_satpy(date_str, 6, ax = ax5, var = var3, crs = crs0, \
        lons = lons3, lats = lats3, lat_lims = lat_lims3, lon_lims = lon_lims3, \
        vmin = None, vmax = 40, \
        ptitle = '', plabel = plabel3, colorbar = True, labelsize = labelsize, \
        zoom=True,save=False)
    plot_GOES_satpy(date_str, 13, ax = ax6, var = var4, crs = crs0, \
        lons = lons2, lats = lats2, lat_lims = lat_lims2, lon_lims = lon_lims2, \
        vmin = None, vmax = None, \
        ptitle = '', plabel = plabel4, colorbar = True, labelsize = labelsize, \
        zoom=True,save=False)
    plot_GOES_satpy(date_str, 8, ax = ax7, var = var5, crs = crs0, \
        lons = lons2, lats = lats2, lat_lims = lat_lims2, lon_lims = lon_lims2, \
        vmin = None, vmax = None, \
        ptitle = '', plabel = plabel5, colorbar = True, labelsize = labelsize, \
        zoom=True,save=False)
    plot_GOES_satpy(date_str, 9, ax = ax8, var = var6, crs = crs0, \
        lons = lons2, lats = lats2, lat_lims = lat_lims2, lon_lims = lon_lims2, \
        vmin = None, vmax = None, \
        ptitle = '', plabel = plabel6, colorbar = True, labelsize = labelsize, \
        zoom=True,save=False)
    plot_GOES_satpy(date_str, 10, ax = ax9, var = var7, crs = crs0, \
        lons = lons2, lats = lats2, lat_lims = lat_lims2, lon_lims = lon_lims2, \
        vmin = None, vmax = None, \
        ptitle = '', plabel = plabel7, colorbar = True, labelsize = labelsize, \
        zoom=True,save=False)

    point_size = 5
    plot_point_on_map(ax4, GOES_dict['plat'][idx1], GOES_dict['plon'][idx1],\
        markersize = point_size)
    plot_point_on_map(ax4, GOES_dict['plat'][idx2], GOES_dict['plon'][idx2],\
        markersize = point_size)
    plot_point_on_map(ax4, GOES_dict2['plat'][idx3], GOES_dict2['plon'][idx3],\
        markersize = point_size)
    sw_idx_b = nearest_gridpoint(GOES_dict['plat'][idx1], GOES_dict['plon'][idx1],\
        lats3, lons3)
    sw_idx_o = nearest_gridpoint(GOES_dict['plat'][idx2], GOES_dict['plon'][idx2],\
        lats3, lons3)
    sw_idx_r = nearest_gridpoint(GOES_dict['plat'][idx4], GOES_dict['plon'][idx4],\
        lats3, lons3)
    print("SWIR")
    print("     Blue   - ", np.array(var3)[sw_idx_b])
    print("     Orange - ", np.array(var3)[sw_idx_o])
    print("     Red    - ", np.array(var3)[sw_idx_r])

    plot_point_on_map(ax5, GOES_dict['plat'][idx1], GOES_dict['plon'][idx1],\
        markersize = point_size)
    plot_point_on_map(ax5, GOES_dict['plat'][idx2], GOES_dict['plon'][idx2],\
        markersize = point_size)
    plot_point_on_map(ax5, GOES_dict2['plat'][idx3], GOES_dict2['plon'][idx3],\
        markersize = point_size)
    lw_idx_b = nearest_gridpoint(GOES_dict['plat'][idx1], GOES_dict['plon'][idx1],\
        lats2, lons2)
    lw_idx_o = nearest_gridpoint(GOES_dict['plat'][idx2], GOES_dict['plon'][idx2],\
        lats2, lons2)
    lw_idx_r = nearest_gridpoint(GOES_dict['plat'][idx4], GOES_dict['plon'][idx4],\
        lats2, lons2)
    lw_idx_p = nearest_gridpoint(GOES_dict['plat'][idx5], GOES_dict['plon'][idx5],\
        lats2, lons2)
    print("TIR")
    print("     Blue   - ", np.array(var4)[lw_idx_b])
    print("     Orange - ", np.array(var4)[lw_idx_o])
    print("     Red    - ", np.array(var4)[lw_idx_r])
    print("     Purple - ", np.array(var4)[lw_idx_p])

    plot_point_on_map(ax6, GOES_dict['plat'][idx1], GOES_dict['plon'][idx1],\
        markersize = point_size)
    plot_point_on_map(ax6, GOES_dict['plat'][idx2], GOES_dict['plon'][idx2],\
        markersize = point_size)
    plot_point_on_map(ax6, GOES_dict2['plat'][idx3], GOES_dict2['plon'][idx3],\
        markersize = point_size)
    print("Upper WV")
    print("     Blue   - ", np.array(var5)[lw_idx_b])
    print("     Orange - ", np.array(var5)[lw_idx_o])
    print("     Red    - ", np.array(var5)[lw_idx_r])
    print("     Purple - ", np.array(var5)[lw_idx_p])

    plot_point_on_map(ax7, GOES_dict['plat'][idx1], GOES_dict['plon'][idx1],\
        markersize = point_size)
    plot_point_on_map(ax7, GOES_dict['plat'][idx2], GOES_dict['plon'][idx2],\
        markersize = point_size)
    plot_point_on_map(ax7, GOES_dict2['plat'][idx3], GOES_dict2['plon'][idx3],\
        markersize = point_size)
    print("Mid WV")
    print("     Blue   - ", np.array(var6)[lw_idx_b])
    print("     Orange - ", np.array(var6)[lw_idx_o])
    print("     Red    - ", np.array(var6)[lw_idx_r])
    print("     Purple - ", np.array(var6)[lw_idx_p])

    plot_point_on_map(ax8, GOES_dict['plat'][idx1], GOES_dict['plon'][idx1],\
        markersize = point_size)
    plot_point_on_map(ax8, GOES_dict['plat'][idx2], GOES_dict['plon'][idx2],\
        markersize = point_size)
    plot_point_on_map(ax8, GOES_dict2['plat'][idx3], GOES_dict2['plon'][idx3],\
        markersize = point_size)
    print("Lower WV")
    print("     Blue   - ", np.array(var7)[lw_idx_b])
    print("     Orange - ", np.array(var7)[lw_idx_o])
    print("     Red    - ", np.array(var7)[lw_idx_r])
    print("     Purple - ", np.array(var7)[lw_idx_p])
    plot_point_on_map(ax9, GOES_dict['plat'][idx1], GOES_dict['plon'][idx1],\
        markersize = point_size)
    plot_point_on_map(ax9, GOES_dict['plat'][idx2], GOES_dict['plon'][idx2],\
        markersize = point_size)
    plot_point_on_map(ax9, GOES_dict2['plat'][idx3], GOES_dict2['plon'][idx3],\
        markersize = point_size)

    # Plot the two channel data for the first point
    # ---------------------------------------------
    ln11 = ax10.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx1,idx1], \
        label = str(goes_channel_dict[\
        str(GOES_dict['channels'][ch_idx1])]['wavelength']) + \
        ' μm', color = 'tab:blue')
    ln21 = ax10.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx1,idx2], \
        label = str(goes_channel_dict[\
        str(GOES_dict['channels'][ch_idx1])]['wavelength']) + \
        ' μm', color = 'tab:orange')
    ln41 = ax10.plot(GOES_dict2['dt_dates'], GOES_dict2['data'][:,ch_idx1,idx3], \
        label = str(goes_channel_dict[\
        str(GOES_dict2['channels'][ch_idx1])]['wavelength']) + \
        ' μm', color = 'tab:green')

    # Plot the two channel data for the second point
    # ----------------------------------------------
    ln12 = ax10.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx2,idx1], \
        label = str(goes_channel_dict[\
        str(GOES_dict['channels'][ch_idx2])]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:blue')
    ln22 = ax10.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx2,idx2], \
        label = str(goes_channel_dict[\
        str(GOES_dict['channels'][ch_idx2])]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:orange')
    ln42 = ax10.plot(GOES_dict2['dt_dates'], GOES_dict2['data'][:,ch_idx2,idx3], \
        label = str(goes_channel_dict[\
        str(GOES_dict2['channels'][ch_idx2])]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:green')
    #ax10.axvline(GOES_dict['dt_dates'][date_idx], color = 'black',\
    #    linestyle = ':')

    ##!#lns = ln11 + ln21 + ln41 + ln12 + ln22 + ln42

    if(ch_idx3 is not None):
        ax102 = ax10.twinx()
        ln31 = ax102.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx3,idx1], \
            label = str(goes_channel_dict[\
            str(GOES_dict['channels'][ch_idx3])]['wavelength']) + \
            ' μm', linestyle = ':', color = 'tab:blue')
        ln32 = ax102.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx3,idx2], \
            label = str(goes_channel_dict[\
            str(GOES_dict['channels'][ch_idx3])]['wavelength']) + \
            ' μm', linestyle = ':', color = 'tab:orange')
        ln33 = ax102.plot(GOES_dict2['dt_dates'], GOES_dict2['data'][:,\
            ch_idx3,idx3], \
            label = str(goes_channel_dict[\
            str(GOES_dict2['channels'][ch_idx3])]['wavelength']) + \
            ' μm', linestyle = ':', color = 'tab:green')
        ax102.set_ylabel('Brightness Temperature [K]')

        ##!#lns = lns + ln31 + ln32 

    labelsize = 10
    font_size = 10
    ax10.set_ylabel('Reflectance [%]', weight = 'bold')
    ax102.set_ylabel('Brightness Temperature [K]', weight = 'bold')
    ax10.grid()
    ax10.xaxis.set_major_formatter(DateFormatter('%m/%d\n%H:%MZ'))
    ax10.tick_params(axis="x", labelsize = font_size + 1)

    ##!#labs = [l.get_label() for l in lns]

    custom_lines = [Line2D([0], [0], color='k'),
                    Line2D([0], [0], color='k', linestyle = '--'),
                    Line2D([0], [0], color='k', linestyle = ':')]

    ax10.legend(custom_lines, ['0.64 μm Refl.', '2.25 μm Refl.', '10.35 μm Bright. Temp.'],\
        fontsize = font_size, loc = 2)
    #ax10.legend(lns, labs, fontsize = font_size)

    # Add subplot labels
    # ------------------
    font_size = 10
    plabels = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)']
    plot_subplot_label(ax4, plabels[0], backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax5, plabels[1], backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax6, plabels[2], backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax7, plabels[3], backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax8, plabels[4], backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax9, plabels[5], backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax10, plabels[6], backgroundcolor = 'white', fontsize = font_size, location = 'upper_right')

    # Add plot text
    # -------------
    plot_figure_text(ax4, 'True Color' , color = 'red', \
        fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax5, \
        str(goes_channel_dict[str(goes_ch2)]['wavelength']) \
        + ' μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    #plot_figure_text(ax5, 'GOES-17 2.25 μm', xval = None, yval = None, transform = None, \
    #    color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax6, \
        str(goes_channel_dict[str(goes_ch3)]['wavelength']) \
        + ' μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax7, \
        str(goes_channel_dict[str(goes_ch4)]['wavelength']) \
        + ' μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax8, \
        str(goes_channel_dict[str(goes_ch5)]['wavelength']) \
        + ' μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax9, \
        str(goes_channel_dict[str(goes_ch6)]['wavelength']) \
        + ' μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')

    # Add text labels to the GOES-17 imagery
    # --------------------------------------
    plot_figure_text(ax5, \
        str(goes_channel_dict[str(goes_ch2)]['short_name']), \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size - 1, backgroundcolor = 'white', \
        location = 'upper_right', halign = 'right')
    plot_figure_text(ax6, \
        str(goes_channel_dict[str(goes_ch3)]['short_name']), \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size - 1, backgroundcolor = 'white', \
        location = 'upper_right', halign = 'right')
    plot_figure_text(ax7, \
        str(goes_channel_dict[str(goes_ch4)]['short_name']), \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size - 1, backgroundcolor = 'white', \
        location = 'upper_right', halign = 'right')
    plot_figure_text(ax8, \
        str(goes_channel_dict[str(goes_ch5)]['short_name']), \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size - 1, backgroundcolor = 'white', \
        location = 'upper_right', halign = 'right')
    plot_figure_text(ax9, \
        str(goes_channel_dict[str(goes_ch6)]['short_name']), \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size - 1, backgroundcolor = 'white', \
        location = 'upper_right', halign = 'right')

    fig.suptitle('GOES-17 Imagery of the Dixie Fire\n00:00 UTC 21 July 2021')
    #fig.suptitle('GOES-17\n00:00 UTC 21 July 2021')
    fig.tight_layout()

    if(save):
        outname = 'goes_combined_fig2_' + date_str + '_v2.png'
        fig.savefig(outname, dpi=300)
        print("Saved",outname)
    else:
        plt.show()

def plot_GOES_ASOS_comp(date_str = '202107210000', \
        goes_ch1 = 'true_color', goes_ch2 = 6, goes_ch3 = 13, \
        goes_ch4 = 8, goes_ch5 = 9, goes_ch6 = 10, \
        ch_idx1 = 0, ch_idx2 = 1, ch_idx3 = 2,\
        ttype1 = 'asos', ttype2 = 'ml', \
        idx1 = 3, idx2 = 8, idx3 = 5, idx4 = 15, idx5 = 20,\
        date_idx = 25, 
        show_smoke = False, composite = True, double_fig = False, \
        zoom = True, save=False):

    date_str2 = date_str
    #date_str2 = '202107202100'
    dt_date_str2 = datetime.strptime(date_str,"%Y%m%d%H%M")

    # Read the GOES data
    # ------------------------
    var2, crs0, lons, lats, lat_lims2, lon_lims2, plabel2   = \
        read_GOES_satpy(date_str2, goes_ch1)
    var3, crs0, lons3, lats3, lat_lims3, lon_lims3, plabel3 = \
        read_GOES_satpy(date_str2, goes_ch2)
    #var3, crs0, lons2, lats2, lat_lims0, lon_lims0, plabel3 = read_GOES_satpy(date_str2, 6)
    var4, crs0, lons2, lats2, lat_lims0, lon_lims0, plabel4 = \
        read_GOES_satpy(date_str2, goes_ch3)
    var5, crs0, lons2, lats2, lat_lims0, lon_lims2, plabel5 = \
        read_GOES_satpy(date_str2, goes_ch4)
    var6, crs0, lons2, lats2, lat_lims0, lon_lims0, plabel6 = \
        read_GOES_satpy(date_str2, goes_ch5)
    var7, crs0, lons2, lats2, lat_lims0, lon_lims0, plabel7 = \
        read_GOES_satpy(date_str2, goes_ch6)

    # Read the GOES time series data
    # ------------------------------
    file_name1 = home_dir + '/Research/GOES/goes_cross_data_' + \
        ttype1 + '_202107201201_202107210331.nc'
    GOES_dict  = read_GOES_time_series_NCDF(file_name1)

    # ----------------------------------------------------------------------
    #
    #  Set up the 6-panel figure
    #
    # ----------------------------------------------------------------------

    plt.close('all')
    fig   = plt.figure(figsize=(9.5,9))
    gs    = fig.add_gridspec(nrows = 3, ncols = 3)
    ax4   = fig.add_subplot(gs[0,0],  projection = crs0)   # GOES True color
    ax5   = fig.add_subplot(gs[0,1],  projection = crs0)   # GOES SWIR
    ax6   = fig.add_subplot(gs[0,2],  projection = crs0)   # GOES TIR 
    ax7   = fig.add_subplot(gs[1,0],  projection = crs0)   # GOES upper WV
    ax8   = fig.add_subplot(gs[1,1],  projection = crs0)   # GOES midle WV
    ax9   = fig.add_subplot(gs[1,2],  projection = crs0)   # GOES lower WV
    ax10  = fig.add_subplot(gs[2,:])                       # GOES time series
    
    labelsize = 10
    plot_GOES_satpy(date_str, goes_ch1, ax = ax4, var = var2, crs = crs0, \
        lons = lons, lats = lats, lat_lims = lat_lims2, lon_lims = lon_lims2, \
        ptitle = '', plabel = plabel2, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)

    ##!## Plot channel 1, 5, 31, and WV data spatial data
    ##!## -----------------------------------------------
    ##!#plot_GOES_satpy(date_str, 2, ax = ax4, var = var2, crs = crs0, \
    ##!#    lons = lons, lats = lats, lat_lims = lat_lims2, lon_lims = lon_lims2, \
    ##!#    vmin = None, vmax = 60, \
    ##!#    ptitle = '', plabel = plabel2, colorbar = True, labelsize = labelsize, \
    ##!#    zoom=True,save=False)
    plot_GOES_satpy(date_str, 6, ax = ax5, var = var3, crs = crs0, \
        lons = lons3, lats = lats3, lat_lims = lat_lims3, lon_lims = lon_lims3, \
        vmin = None, vmax = 40, \
        ptitle = '', plabel = plabel3, colorbar = True, labelsize = labelsize, \
        zoom=True,save=False)
    plot_GOES_satpy(date_str, 13, ax = ax6, var = var4, crs = crs0, \
        lons = lons2, lats = lats2, lat_lims = lat_lims2, lon_lims = lon_lims2, \
        vmin = None, vmax = None, \
        ptitle = '', plabel = plabel4, colorbar = True, labelsize = labelsize, \
        zoom=True,save=False)
    plot_GOES_satpy(date_str, 8, ax = ax7, var = var5, crs = crs0, \
        lons = lons2, lats = lats2, lat_lims = lat_lims2, lon_lims = lon_lims2, \
        vmin = None, vmax = None, \
        ptitle = '', plabel = plabel5, colorbar = True, labelsize = labelsize, \
        zoom=True,save=False)
    plot_GOES_satpy(date_str, 9, ax = ax8, var = var6, crs = crs0, \
        lons = lons2, lats = lats2, lat_lims = lat_lims2, lon_lims = lon_lims2, \
        vmin = None, vmax = None, \
        ptitle = '', plabel = plabel6, colorbar = True, labelsize = labelsize, \
        zoom=True,save=False)
    plot_GOES_satpy(date_str, 10, ax = ax9, var = var7, crs = crs0, \
        lons = lons2, lats = lats2, lat_lims = lat_lims2, lon_lims = lon_lims2, \
        vmin = None, vmax = None, \
        ptitle = '', plabel = plabel7, colorbar = True, labelsize = labelsize, \
        zoom=True,save=False)

    point_size = 5
    ax4.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax4.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)
    ax4.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax4.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)
    sw_idx_b = nearest_gridpoint(GOES_dict['plat'][idx1], GOES_dict['plon'][idx1],\
        lats3, lons3)
    sw_idx_o = nearest_gridpoint(GOES_dict['plat'][idx2], GOES_dict['plon'][idx2],\
        lats3, lons3)
    print("SWIR")
    print("     Blue   - ", np.array(var3)[sw_idx_b])
    print("     Orange - ", np.array(var3)[sw_idx_o])
    ax5.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax5.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)
    ax5.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax5.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)
    lw_idx_b = nearest_gridpoint(GOES_dict['plat'][idx1], GOES_dict['plon'][idx1],\
        lats2, lons2)
    lw_idx_o = nearest_gridpoint(GOES_dict['plat'][idx2], GOES_dict['plon'][idx2],\
        lats2, lons2)
    print("TIR")
    print("     Blue   - ", np.array(var4)[lw_idx_b])
    print("     Orange - ", np.array(var4)[lw_idx_o])
    ax6.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax6.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)
    ax6.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax6.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)
    print("Upper WV")
    print("     Blue   - ", np.array(var5)[lw_idx_b])
    print("     Orange - ", np.array(var5)[lw_idx_o])
    ax7.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax7.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)
    ax7.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax7.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)
    print("Mid WV")
    print("     Blue   - ", np.array(var6)[lw_idx_b])
    print("     Orange - ", np.array(var6)[lw_idx_o])
    ax8.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax8.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)
    ax8.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax8.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)
    print("Lower WV")
    print("     Blue   - ", np.array(var7)[lw_idx_b])
    print("     Orange - ", np.array(var7)[lw_idx_o])
    ax9.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax9.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)
    ax9.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax9.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)

    # Plot the two channel data for the first point
    # ---------------------------------------------
    ln11 = ax10.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx1,idx1], \
        label = str(goes_channel_dict[\
        str(GOES_dict['channels'][ch_idx1])]['wavelength']) + \
        ' μm', color = 'tab:blue')
    ln21 = ax10.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx1,idx2], \
        label = str(goes_channel_dict[\
        str(GOES_dict['channels'][ch_idx1])]['wavelength']) + \
        ' μm', color = 'tab:orange')

    # Plot the two channel data for the second point
    # ----------------------------------------------
    ln12 = ax10.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx2,idx1], \
        label = str(goes_channel_dict[\
        str(GOES_dict['channels'][ch_idx2])]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:blue')
    ln22 = ax10.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx2,idx2], \
        label = str(goes_channel_dict[\
        str(GOES_dict['channels'][ch_idx2])]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:orange')
    #ax10.axvline(GOES_dict['dt_dates'][date_idx], color = 'black',\
    #    linestyle = ':')

    ##!#lns = ln11 + ln21 + ln41 + ln12 + ln22 + ln42

    if(ch_idx3 is not None):
        ax102 = ax10.twinx()
        ln31 = ax102.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx3,idx1], \
            label = str(goes_channel_dict[\
            str(GOES_dict['channels'][ch_idx3])]['wavelength']) + \
            ' μm', linestyle = ':', color = 'tab:blue')
        ln32 = ax102.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx3,idx2], \
            label = str(goes_channel_dict[\
            str(GOES_dict['channels'][ch_idx3])]['wavelength']) + \
            ' μm', linestyle = ':', color = 'tab:orange')
        ax102.set_ylabel('Brightness Temperature [K]')

        ##!#lns = lns + ln31 + ln32 

    labelsize = 10
    font_size = 10
    ax10.set_ylabel('Reflectance [%]', weight = 'bold')
    ax102.set_ylabel('Brightness Temperature [K]', weight = 'bold')
    ax10.grid()
    ax10.xaxis.set_major_formatter(DateFormatter('%m/%d\n%H:%MZ'))
    ax10.tick_params(axis="x", labelsize = font_size + 1)

    ##!#labs = [l.get_label() for l in lns]

    custom_lines = [Line2D([0], [0], color='k'),
                    Line2D([0], [0], color='k', linestyle = '--'),
                    Line2D([0], [0], color='k', linestyle = ':')]

    ax10.legend(custom_lines, ['0.64 μm Refl.', '2.25 μm Refl.', '10.35 μm Bright. Temp.'],\
        fontsize = font_size, loc = 2)
    #ax10.legend(lns, labs, fontsize = font_size)

    # Add subplot labels
    # ------------------
    font_size = 10
    plabels = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)']
    plot_subplot_label(ax4, plabels[0], backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax5, plabels[1], backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax6, plabels[2], backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax7, plabels[3], backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax8, plabels[4], backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax9, plabels[5], backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax10, plabels[6], backgroundcolor = 'white', fontsize = font_size, location = 'upper_right')

    # Add plot text
    # -------------
    plot_figure_text(ax4, 'True Color' , color = 'red', \
        fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax5, \
        str(goes_channel_dict[str(goes_ch2)]['wavelength']) \
        + ' μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    #plot_figure_text(ax5, 'GOES-17 2.25 μm', xval = None, yval = None, transform = None, \
    #    color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax6, \
        str(goes_channel_dict[str(goes_ch3)]['wavelength']) \
        + ' μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax7, \
        str(goes_channel_dict[str(goes_ch4)]['wavelength']) \
        + ' μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax8, \
        str(goes_channel_dict[str(goes_ch5)]['wavelength']) \
        + ' μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax9, \
        str(goes_channel_dict[str(goes_ch6)]['wavelength']) \
        + ' μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')

    fig.suptitle('GOES-17\n00:00 UTC 21 July 2021')
    fig.tight_layout()

    if(save):
        outname = 'goes_combined_asos_comp_' + date_str + '_v2.png'
        fig.savefig(outname, dpi=300)
        print("Saved",outname)
    else:
        plt.show()
##for a in ax:
##    a.set_xticklabels([])
##    a.set_yticklabels([])
##    a.set_aspect('equal')
##
##fig.subplots_adjust(wspace=0, hspace=0)

def plot_GOES_figure2(date_str1 = '202107201330', date_str2 = '202107201830',\
        date_str3 = '202107202300', date_str4 = '202107210330', \
        ch1 = 'true_color', ch2 = 6, ch3 = 9, \
        ch4 = 10, ch5 = 13, \
        ch_idx1 = 0, ch_idx2 = 1, ch_idx3 = 2,\
        ttype1 = 'low', ttype2 = 'ml',\
        idx1 = 3, idx2 = 8, idx3 = 5, \
        show_smoke = False, composite = True, double_fig = False, \
        zoom = True, add_wv_time = True, save=False):

    dt_date_str1 = datetime.strptime(date_str1,"%Y%m%d%H%M")
    dt_date_str2 = datetime.strptime(date_str2,"%Y%m%d%H%M")
    dt_date_str3 = datetime.strptime(date_str3,"%Y%m%d%H%M")
    dt_date_str4 = datetime.strptime(date_str4,"%Y%m%d%H%M")

    ##!## ----------------------------------------------------------------------
    ##!##
    ##!## Read the GOES data
    ##!##
    ##!## ----------------------------------------------------------------------
    ##!## Read the GOES data for time 1
    ##!## -----------------------------
    ##!#var2, crs0, lons, lats, lat_lims2, lon_lims2, plabel2   = \
    ##!#    read_GOES_satpy(date_str2, goes_ch1)
    ##!#var3, crs0, lons3, lats3, lat_lims3, lon_lims3, plabel3 = \
    ##!#    read_GOES_satpy(date_str2, goes_ch2)
    ##!##var3, crs0, lons2, lats2, lat_lims0, lon_lims0, plabel3 = read_GOES_satpy(date_str2, 6)
    ##!#var4, crs0, lons2, lats2, lat_lims0, lon_lims0, plabel4 = \
    ##!#    read_GOES_satpy(date_str2, goes_ch3)
    ##!#var5, crs0, lons2, lats2, lat_lims0, lon_lims2, plabel5 = \
    ##!#    read_GOES_satpy(date_str2, goes_ch4)
    ##!#var6, crs0, lons2, lats2, lat_lims0, lon_lims0, plabel6 = \
    ##!#    read_GOES_satpy(date_str2, goes_ch5)
    ##!#var7, crs0, lons2, lats2, lat_lims0, lon_lims0, plabel7 = \
    ##!#    read_GOES_satpy(date_str2, goes_ch6)

    # Read the GOES time series data
    # ------------------------------
    file_name1 = home_dir + '/Research/GOES/goes_cross_data_' + \
        ttype1 + '_202107201201_202107210331.nc'
    file_name2 = home_dir + '/Research/GOES/goes_cross_data_' + \
        ttype2 + '_202107201201_202107210331.nc'
    GOES_dict  = read_GOES_time_series_NCDF(file_name1)
    GOES_dict2 = read_GOES_time_series_NCDF(file_name2)

    # ----------------------------------------------------------------------
    #
    #  Set up the 6-panel figure
    #
    # ----------------------------------------------------------------------
    plt.close('all')
    if(add_wv_time):
        figsize = (6.0, 15)
        numrows = 7
    else:
        figsize = (5.65, 12)
        numrows = 6
    fig = plt.figure(figsize = figsize)
    gs = fig.add_gridspec(nrows = numrows, ncols = 4)

    dt_date_strs = [dt_date_str1, dt_date_str2, dt_date_str3,\
        dt_date_str4]
    date_strs = [date_str1, date_str2, date_str3, date_str4]
    channels = [ch1, ch2, ch3, ch4, ch5]

    colorbar = False 
    labelsize = 8
    if(add_wv_time):
        point_size = 4
    else:
        point_size = 3
    for ii, dstr in enumerate(date_strs):
        #if(ii == len(date_strs) - 1):
        #    colorbar = True
        print(ii, dstr)
        for jj, channel in enumerate(channels):
            print('  ', jj, channel)
            # Read in the five data types
            # ---------------------------
            var1, crs1, lons1, lats1, lat_lims1, lon_lims1, plabel1 = \
                read_GOES_satpy(dstr, channel)
            tax = fig.add_subplot(gs[jj,ii], projection = crs1)

            if(channel == 'true_color'):
                plot_GOES_satpy(dstr, channel, ax = tax, var = var1, \
                    crs = crs1, lons = lons1, lats = lats1, \
                    lat_lims = lat_lims1, lon_lims = lon_lims1, \
                    ptitle = '', plabel = '', \
                    colorbar = False, labelsize = labelsize + 1, zoom=True,\
                    save=False)
            else:
                plot_GOES_satpy(dstr, channel, ax = tax, var = var1, crs = crs1, \
                     lons = lons1, lats = lats1, lat_lims = lat_lims1, \
                     lon_lims = lon_lims1, vmin = min_dict[channel], \
                     vmax = max_dict[channel], \
                     ptitle = '', plabel = plabel1, colorbar = colorbar, \
                     labelsize = labelsize, zoom=True,save=False)

            tax.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
                    linewidth=2, markersize = point_size + 2, marker='.',
                    color = 'black', transform=datacrs)
            tax.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
                    linewidth=2, markersize = point_size, marker='.',
                    transform=datacrs)
            tax.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
                    linewidth=2, markersize = point_size + 2, marker='.',
                    color = 'black', transform=datacrs)
            tax.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
                    linewidth=2, markersize = point_size, marker='.',
                    transform=datacrs)
            tax.plot(GOES_dict2['plon'][idx3], GOES_dict2['plat'][idx3], \
                    linewidth=2, markersize = point_size + 2, marker='.',
                    color = 'black', transform=datacrs)
            tax.plot(GOES_dict2['plon'][idx3], GOES_dict2['plat'][idx3], \
                    linewidth=2, markersize = point_size, marker='.',
                    transform=datacrs)

            if(jj == 0):
                tax.set_title(dt_date_strs[ii].strftime('%Y/%m/%d\n%H:%M UTC'),\
                    fontsize = 8)
   
    ax10  = fig.add_subplot(gs[5,:]) # time series of GOES data

    # Plot the two channel data for the first point
    # ---------------------------------------------
    ln11 = ax10.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx1,idx1], \
        label = str(goes_channel_dict[\
        str(GOES_dict['channels'][ch_idx1])]['wavelength']) + \
        ' μm', color = 'tab:blue')
    ln21 = ax10.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx1,idx2], \
        label = str(goes_channel_dict[\
        str(GOES_dict['channels'][ch_idx1])]['wavelength']) + \
        ' μm', color = 'tab:orange')
    ln41 = ax10.plot(GOES_dict2['dt_dates'], GOES_dict2['data'][:,ch_idx1,idx3], \
        label = str(goes_channel_dict[\
        str(GOES_dict2['channels'][ch_idx1])]['wavelength']) + \
        ' μm', color = 'tab:green')

    # Plot the two channel data for the second point
    # ----------------------------------------------
    ln12 = ax10.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx2,idx1], \
        label = str(goes_channel_dict[\
        str(GOES_dict['channels'][ch_idx2])]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:blue')
    ln22 = ax10.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx2,idx2], \
        label = str(goes_channel_dict[\
        str(GOES_dict['channels'][ch_idx2])]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:orange')
    ln42 = ax10.plot(GOES_dict2['dt_dates'], GOES_dict2['data'][:,ch_idx2,idx3], \
        label = str(goes_channel_dict[\
        str(GOES_dict2['channels'][ch_idx2])]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:green')

    if(ch_idx3 is not None):
        ax102 = ax10.twinx()
        ln31 = ax102.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx3,idx1], \
            label = str(goes_channel_dict[\
            str(GOES_dict['channels'][ch_idx3])]['wavelength']) + \
            ' μm', linestyle = ':', color = 'tab:blue')
        ln32 = ax102.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx3,idx2], \
            label = str(goes_channel_dict[\
            str(GOES_dict['channels'][ch_idx3])]['wavelength']) + \
            ' μm', linestyle = ':', color = 'tab:orange')
        ln33 = ax102.plot(GOES_dict2['dt_dates'], GOES_dict2['data'][:,\
            ch_idx3,idx3], \
            label = str(goes_channel_dict[\
            str(GOES_dict2['channels'][ch_idx3])]['wavelength']) + \
            ' μm', linestyle = ':', color = 'tab:green')
        ax102.set_ylabel('BT [K]')
        ax102.tick_params(axis="y", labelsize = 7)
    ax10.xaxis.set_major_formatter(DateFormatter('%m/%d\n%H:%MZ'))
    ax10.tick_params(axis="x", labelsize = labelsize)
    ax10.tick_params(axis="y", labelsize = 7)
    ax10.set_ylabel('Reflectance')
    ax10.grid()

    custom_lines = [Line2D([0], [0], color='k'),
                    Line2D([0], [0], color='k', linestyle = '--'),\
                    Line2D([0], [0], color='k', linestyle = ':')]

    ax10.legend(custom_lines, ['0.64 μm', '2.25 μm', '10.35 μm'],\
        fontsize = 8, loc = 2)

    if(add_wv_time):
        ax11  = fig.add_subplot(gs[6,:]) # time series of GOES data

        # Plot the two channel data for the first point
        # ---------------------------------------------
        ln111 = ax11.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,4,idx1], \
            label = str(goes_channel_dict[\
            str(GOES_dict['channels'][3])]['wavelength']) + \
            ' μm', color = 'tab:blue')
        ln211 = ax11.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,4,idx2], \
            label = str(goes_channel_dict[\
            str(GOES_dict['channels'][3])]['wavelength']) + \
            ' μm', color = 'tab:orange')
        ln411 = ax11.plot(GOES_dict2['dt_dates'], GOES_dict2['data'][:,4,idx3], \
            label = str(goes_channel_dict[\
            str(GOES_dict2['channels'][3])]['wavelength']) + \
            ' μm', color = 'tab:green')

        # Plot the two channel data for the second point
        # ----------------------------------------------
        ln121 = ax11.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,5,idx1], \
            label = str(goes_channel_dict[\
            str(GOES_dict['channels'][5])]['wavelength']) + \
            ' μm', linestyle = '--', color = 'tab:blue')
        ln221 = ax11.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,5,idx2], \
            label = str(goes_channel_dict[\
            str(GOES_dict['channels'][5])]['wavelength']) + \
            ' μm', linestyle = '--', color = 'tab:orange')
        ln421 = ax11.plot(GOES_dict2['dt_dates'], GOES_dict2['data'][:,5,idx3], \
            label = str(goes_channel_dict[\
            str(GOES_dict2['channels'][5])]['wavelength']) + \
            ' μm', linestyle = '--', color = 'tab:green')

        ##!#if(ch_idx3 is not None):
        ##!#    ax102 = ax10.twinx()
        ##!#    ln31 = ax102.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx3,idx1], \
        ##!#        label = str(goes_channel_dict[\
        ##!#        str(GOES_dict['channels'][ch_idx3])]['wavelength']) + \
        ##!#        ' μm', linestyle = ':', color = 'tab:blue')
        ##!#    ln32 = ax102.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx3,idx2], \
        ##!#        label = str(goes_channel_dict[\
        ##!#        str(GOES_dict['channels'][ch_idx3])]['wavelength']) + \
        ##!#        ' μm', linestyle = ':', color = 'tab:orange')
        ##!#    ln33 = ax102.plot(GOES_dict2['dt_dates'], GOES_dict2['data'][:,\
        ##!#        ch_idx3,idx3], \
        ##!#        label = str(goes_channel_dict[\
        ##!#        str(GOES_dict2['channels'][ch_idx3])]['wavelength']) + \
        ##!#        ' μm', linestyle = ':', color = 'tab:green')
        ##!#    ax102.set_ylabel('Brightness Temperature [K]')
        ax11.xaxis.set_major_formatter(DateFormatter('%m/%d\n%H:%MZ'))
        ax11.tick_params(axis="x", labelsize = labelsize + 1)
        ax11.tick_params(axis="y", labelsize = 7)
        ax11.set_ylabel('BT')
        ax11.grid()
        custom_lines = [Line2D([0], [0], color='k'),
                        Line2D([0], [0], color='k', linestyle = '--')]

        ax11.legend(custom_lines, ['6.95 μm', '7.34 μm'],\
            fontsize = 8, loc = 2)

    # Add plot text
    # -------------
    row_label_size = 10
    if(add_wv_time):
        fig.text(0.10, 0.83, 'True Color', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)
        fig.text(0.10, 0.74, '2.25 μm', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)
        fig.text(0.10, 0.64, '6.95 μm', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)
        fig.text(0.10, 0.54, '7.34 μm', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)
        fig.text(0.10, 0.45, '10.35 μm', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)
    else:
        fig.text(0.10, 0.82, 'True Color', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)
        fig.text(0.10, 0.72, '2.25 μm', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)
        fig.text(0.10, 0.60, '6.95 μm', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)
        fig.text(0.10, 0.49, '7.34 μm', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)
        fig.text(0.10, 0.37, '10.35 μm', ha='center', va='center', \
            rotation='vertical',weight='bold',fontsize=row_label_size)

    fig.autofmt_xdate()
    fig.subplots_adjust(wspace=0, hspace=0)

    if(save):
        if(add_wv_time):
            outname = 'goes_multi_big_file_wv_test.png'
        else:
            outname = 'goes_multi_big_file_test.png'
        fig.savefig(outname, dpi = 300)
        print('Saved image', outname)
    else:
        plt.show()
    return
    ##!#axi11  = fig.add_subplot(gs[0,0], projection = crs1) # true color    
    ##!#axi12  = fig.add_subplot(gs[1,0], projection = crs8) # MODIS Ch 7
    ##!#axi13  = fig.add_subplot(gs[2,0], projection = crs8) # MODIS Ch 31
    ##!#axi14  = fig.add_subplot(gs[3,0], projection = crs0) # GOES VIS 
    ##!#axi15  = fig.add_subplot(gs[4,0], projection = crs0)   # GOES SWIR
    ##!#axi21  = fig.add_subplot(gs[0,1], projection = crs1) # true color    
    ##!#axi22  = fig.add_subplot(gs[1,1], projection = crs8) # MODIS Ch 7
    ##!#axi23  = fig.add_subplot(gs[2,1], projection = crs8) # MODIS Ch 31
    ##!#axi24  = fig.add_subplot(gs[3,1], projection = crs0) # GOES VIS 
    ##!#axi25  = fig.add_subplot(gs[4,1], projection = crs0)   # GOES SWIR
    axt1  = fig.add_subplot(gs[5,:]) # time series of GOES data
    ##!#axt2  = fig.add_subplot(gs[6,:]) # time series of GOES data


    ##!## Plot the true-color data for the previous date
    ##!## ----------------------------------------------
    ##!#ax1.imshow(var1.data, transform = crs1, extent=(var1.x[0], var1.x[-1], \
    ##!#    var1.y[-1], var1.y[0]), origin='upper')

    ##!## Zoom in the figure if desired
    ##!## -----------------------------
    ##!#if(zoom):
    ##!#    ax1.set_extent([lon_lims1[0],lon_lims1[1],lat_lims1[0],\
    ##!#        lat_lims1[1]],crs = datacrs)

    # ----------------------------------------------------------------------
    #
    # Plot the data in the figure
    #
    # ----------------------------------------------------------------------

    # Plot channel 1, 5, 31, and WV data spatial data
    # -----------------------------------------------
    ##!#plot_MODIS_spatial(MODIS_data_ch7,  ax2, zoom = zoom, ptitle = '', labelsize = 10, labelticksize = 9)
    ##!#plot_MODIS_spatial(MODIS_data_ch31, ax3, zoom = zoom, ptitle = '', labelsize = 10, labelticksize = 9)
    labelsize = 10

    # Plot channel 1, 5, 31, and WV data spatial data
    # -----------------------------------------------
    plot_GOES_satpy(date_str, 2, ax = ax4, var = var2, crs = crs0, \
        lons = lons, lats = lats, lat_lims = lat_lims2, lon_lims = lon_lims2, \
        vmin = None, vmax = 60, \
        ptitle = '', plabel = plabel2, colorbar = True, labelsize = labelsize, \
        zoom=True,save=False)
    plot_GOES_satpy(date_str, 6, ax = ax5, var = var3, crs = crs0, \
        lons = lons3, lats = lats3, lat_lims = lat_lims3, lon_lims = lon_lims3, \
        vmin = None, vmax = 40, \
        ptitle = '', plabel = plabel3, colorbar = True, labelsize = labelsize, \
        zoom=True,save=False)
    plot_GOES_satpy(date_str, 13, ax = ax6, var = var4, crs = crs0, \
        lons = lons2, lats = lats2, lat_lims = lat_lims2, lon_lims = lon_lims2, \
        vmin = None, vmax = None, \
        ptitle = '', plabel = plabel4, colorbar = True, labelsize = labelsize, \
        zoom=True,save=False)
    plot_GOES_satpy(date_str, 8, ax = ax7, var = var5, crs = crs0, \
        lons = lons2, lats = lats2, lat_lims = lat_lims2, lon_lims = lon_lims2, \
        vmin = None, vmax = None, \
        ptitle = '', plabel = plabel5, colorbar = True, labelsize = labelsize, \
        zoom=True,save=False)
    plot_GOES_satpy(date_str, 9, ax = ax8, var = var6, crs = crs0, \
        lons = lons2, lats = lats2, lat_lims = lat_lims2, lon_lims = lon_lims2, \
        vmin = None, vmax = None, \
        ptitle = '', plabel = plabel6, colorbar = True, labelsize = labelsize, \
        zoom=True,save=False)
    plot_GOES_satpy(date_str, 10, ax = ax9, var = var7, crs = crs0, \
        lons = lons2, lats = lats2, lat_lims = lat_lims2, lon_lims = lon_lims2, \
        vmin = None, vmax = None, \
        ptitle = '', plabel = plabel7, colorbar = True, labelsize = labelsize, \
        zoom=True,save=False)

    point_size = 5
    ax4.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax4.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)
    ax4.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax4.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)
    ax4.plot(GOES_dict2['plon'][idx3], GOES_dict2['plat'][idx3], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax4.plot(GOES_dict2['plon'][idx3], GOES_dict2['plat'][idx3], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)
    ax5.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax5.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)
    ax5.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax5.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)
    ax5.plot(GOES_dict2['plon'][idx3], GOES_dict2['plat'][idx3], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax5.plot(GOES_dict2['plon'][idx3], GOES_dict2['plat'][idx3], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)
    ax6.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax6.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)
    ax6.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax6.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)
    ax6.plot(GOES_dict2['plon'][idx3], GOES_dict2['plat'][idx3], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax6.plot(GOES_dict2['plon'][idx3], GOES_dict2['plat'][idx3], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)

    # Plot the two channel data for the first point
    # ---------------------------------------------
    ln11 = ax10.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx1,idx1], \
        label = str(goes_channel_dict[\
        str(GOES_dict['channels'][ch_idx1])]['wavelength']) + \
        ' μm', color = 'tab:blue')
    ln21 = ax10.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx1,idx2], \
        label = str(goes_channel_dict[\
        str(GOES_dict['channels'][ch_idx1])]['wavelength']) + \
        ' μm', color = 'tab:orange')
    ln41 = ax10.plot(GOES_dict2['dt_dates'], GOES_dict2['data'][:,ch_idx1,idx3], \
        label = str(goes_channel_dict[\
        str(GOES_dict2['channels'][ch_idx1])]['wavelength']) + \
        ' μm', color = 'tab:green')

    # Plot the two channel data for the second point
    # ----------------------------------------------
    ln12 = ax10.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx2,idx1], \
        label = str(goes_channel_dict[\
        str(GOES_dict['channels'][ch_idx2])]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:blue')
    ln22 = ax10.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx2,idx2], \
        label = str(goes_channel_dict[\
        str(GOES_dict['channels'][ch_idx2])]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:orange')
    ln42 = ax10.plot(GOES_dict2['dt_dates'], GOES_dict2['data'][:,ch_idx2,idx3], \
        label = str(goes_channel_dict[\
        str(GOES_dict2['channels'][ch_idx2])]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:green')
    #ax10.axvline(GOES_dict['dt_dates'][date_idx], color = 'black',\
    #    linestyle = ':')

    ##!#lns = ln11 + ln21 + ln41 + ln12 + ln22 + ln42

    if(ch_idx3 is not None):
        ax102 = ax10.twinx()
        ln31 = ax102.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx3,idx1], \
            label = str(goes_channel_dict[\
            str(GOES_dict['channels'][ch_idx3])]['wavelength']) + \
            ' μm', linestyle = ':', color = 'tab:blue')
        ln32 = ax102.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx3,idx2], \
            label = str(goes_channel_dict[\
            str(GOES_dict['channels'][ch_idx3])]['wavelength']) + \
            ' μm', linestyle = ':', color = 'tab:orange')
        ln33 = ax102.plot(GOES_dict2['dt_dates'], GOES_dict2['data'][:,\
            ch_idx3,idx3], \
            label = str(goes_channel_dict[\
            str(GOES_dict2['channels'][ch_idx3])]['wavelength']) + \
            ' μm', linestyle = ':', color = 'tab:green')
        ax102.set_ylabel('Brightness Temperature [K]')

        ##!#lns = lns + ln31 + ln32 

    ax10.set_title('GOES-17')
    labelsize = 10
    font_size = 10
    ax10.set_ylabel(plabel1.replace('_',' '), weight = 'bold')
    ax102.set_ylabel(plabel7.replace('_',' '), weight = 'bold')
    ax10.grid()
    ax10.xaxis.set_major_formatter(DateFormatter('%m/%d\n%H:%MZ'))
    ax10.tick_params(axis="x", labelsize = font_size + 1)

    ##!#labs = [l.get_label() for l in lns]

    custom_lines = [Line2D([0], [0], color='k'),
                    Line2D([0], [0], color='k', linestyle = '--'),
                    Line2D([0], [0], color='k', linestyle = ':')]

    ax10.legend(custom_lines, ['0.64 μm', '2.25 μm', '10.35 μm'],\
        fontsize = font_size, loc = 2)
    #ax10.legend(lns, labs, fontsize = font_size)

    # Add subplot labels
    # ------------------
    font_size = 10
    if(double_fig):
        plabels = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)']
    else:
        plabels = ['(d)','(e)','(f)','(g)','(h)','(i)','(j)']
    plot_subplot_label(ax1, '(a)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax2, '(b)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax3, '(c)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax4, plabels[0], backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax5, plabels[1], backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax6, plabels[2], backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax7, plabels[3], backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax8, plabels[4], backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax9, plabels[5], backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax10, plabels[6], backgroundcolor = 'white', fontsize = font_size, location = 'upper_right')

    # Add plot text
    # -------------
    plot_figure_text(ax1, 'MODIS ' + \
        str(np.round(np.mean(channel_dict[str(1)]['Bandwidth']), 2))  \
        + ' μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax2, 'MODIS ' + \
        str(np.round(np.mean(channel_dict[str(modis_ch1)]['Bandwidth']), 2)) \
        + ' μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    #plot_figure_text(ax2, 'MODIS 2.2 μm', xval = None, yval = None, transform = None, \
    #    color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax3, 'MODIS ' + \
        str(np.round(np.mean(channel_dict[str(modis_ch2)]['Bandwidth']), 2)) \
        + ' μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax4, 'GOES-17 ' + \
        str(goes_channel_dict[str(goes_ch1)]['wavelength']) \
        + ' μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax5, 'GOES-17 ' + \
        str(goes_channel_dict[str(goes_ch2)]['wavelength']) \
        + ' μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    #plot_figure_text(ax5, 'GOES-17 2.25 μm', xval = None, yval = None, transform = None, \
    #    color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax6, 'GOES-17 ' + \
        str(goes_channel_dict[str(goes_ch3)]['wavelength']) \
        + ' μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax7, 'GOES-17 ' + \
        str(goes_channel_dict[str(goes_ch4)]['wavelength']) \
        + ' μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax8, 'GOES-17 ' + \
        str(goes_channel_dict[str(goes_ch5)]['wavelength']) \
        + ' μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax9, 'GOES-17 ' + \
        str(goes_channel_dict[str(goes_ch6)]['wavelength']) \
        + ' μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')

    if(double_fig):
        fig1.suptitle('21:10 UTC 22 July 2021')
        fig2.suptitle('00:00 UTC 21 July 2021')
        fig1.tight_layout()
        fig2.tight_layout()
    else:
        fig.text(0.09, 0.880, '21:10 UTC 2021/07/22', ha='center', va='center', \
            rotation='vertical', weight = 'bold', fontsize = labelsize + 2)
        fig.text(0.09, 0.520, '-------------------- 00:00 UTC 2021/07/21 --------------------', ha='center', va='center', \
            rotation='vertical', weight = 'bold', fontsize = labelsize + 2)
        #plt.suptitle(date_str)

        fig.tight_layout()
    #MODIS_data_ch7.clear()
    #MODIS_data_ch31.clear()

    if(save):
        if(double_fig):
            outname1 = 'modis_total_combined_' + date_str + '_fig1_v5_modis.png'
            outname2 = 'modis_total_combined_' + date_str + '_fig1_v5_goes.png'
            fig1.savefig(outname1, dpi=300)
            fig2.savefig(outname2, dpi=300)
            print("Saved",outname1)
            print("Saved",outname2)
        else:
            outname = 'modis_total_combined_' + date_str + '_fig1_v52.png'
            fig.savefig(outname, dpi=300)
            print("Saved",outname)
    else:
        plt.show()
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

# This is written to compare time series of any two desired channels
# (as long as they're in the GOES_dict structure)
def plot_GOES_time_series_channel_comp(GOES_dict, ch_idx1, ch_idx2, \
        idx1, idx2, ch_idx3 = None, date_idx = 20, save_dir = './', \
        save = False):

    channel1 = int(GOES_dict['channels'][ch_idx1])
    channel2 = int(GOES_dict['channels'][ch_idx2])
    if(ch_idx3 is not None):
        channel3 = int(GOES_dict['channels'][ch_idx3])
        figsize = (9.5, 3.5)
    else:
        figsize = (8.5, 3.5)
        

    dt_date_str = GOES_dict['dt_dates'][date_idx]
    date_str = dt_date_str.strftime('%Y%m%d%H%M')
    var, crs, lons, lats, lat_lims, lon_lims, plabel = \
        read_GOES_satpy(date_str, channel1)

    plt.close('all')
    fig = plt.figure(figsize = figsize)
    gs = fig.add_gridspec(nrows = 1, ncols = 3)
    ax2  = fig.add_subplot(gs[2], projection = crs)   # true color    
    ax1  = fig.add_subplot(gs[0:2]) # Ch 1

    # Plot the two channel data for the first point
    # ---------------------------------------------
    ln11 = ax1.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx1,idx1], \
        label = str(goes_channel_dict[\
        str(channel1)]['wavelength']) + \
        ' μm', color = 'tab:blue')
    ln21 = ax1.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx1,idx2], \
        label = str(goes_channel_dict[\
        str(channel1)]['wavelength']) + \
        ' μm', color = 'tab:orange')

    # Plot the two channel data for the second point
    # ----------------------------------------------
    ln12 = ax1.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx2,idx1], \
        label = str(goes_channel_dict[\
        str(channel2)]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:blue')
    ln22 = ax1.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx2,idx2], \
        label = str(goes_channel_dict[\
        str(channel2)]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:orange')
    ax1.axvline(dt_date_str, color = 'black',\
        linestyle = ':')

    lns = ln11 + ln21 + ln12 + ln22 

    if(ch_idx3 is not None):
        ax12 = ax1.twinx()
        ln31 = ax12.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx3,idx1], \
            label = str(goes_channel_dict[\
            str(channel3)]['wavelength']) + \
            ' μm', linestyle = ':', color = 'tab:blue')
        ln32 = ax12.plot(GOES_dict['dt_dates'], GOES_dict['data'][:,ch_idx3,idx2], \
            label = str(goes_channel_dict[\
            str(channel3)]['wavelength']) + \
            ' μm', linestyle = ':', color = 'tab:orange')
        ax12.set_ylabel('Brightness Temperature [K]')

        lns = lns + ln31 + ln32 
        print(len(lns))

    labelsize = 10
    font_size = 8
    ax1.set_ylabel(plabel.replace('_',' '))
    ax1.grid()
    ax1.xaxis.set_major_formatter(DateFormatter('%m/%d\n%H:%MZ'))
    ax1.tick_params(axis="x", labelsize = font_size)

    # Print out the values of each channel at each point
    # --------------------------------------------------
    print(str(goes_channel_dict[str(channel1)]['wavelength']) + \
        ' μm\n')
    print("    Point blue   - ", GOES_dict['data'][date_idx,ch_idx1,idx1])
    print("    Point orange - ", GOES_dict['data'][date_idx,ch_idx1,idx2])
    print(str(goes_channel_dict[str(channel2)]['wavelength']) + \
        ' μm\n')
    print("    Point blue   - ", GOES_dict['data'][date_idx,ch_idx2,idx1])
    print("    Point orange - ", GOES_dict['data'][date_idx,ch_idx2,idx2])
    if(ch_idx3 is not None):
        print(str(goes_channel_dict[str(channel3)]['wavelength']) + \
            ' μm\n')
        print("    Point blue   - ", GOES_dict['data'][date_idx,ch_idx3,idx1])
        print("    Point orange - ", GOES_dict['data'][date_idx,ch_idx3,idx2])

    ##!#print("Point Blue:\n")
    ##!#print("    " + str(goes_channel_dict[str(channel1)]['wavelength']) + \
    ##!#    ' μm - ', GOES_dict['data'][date_idx,ch_idx1,idx1])
    ##!#print("    " + str(goes_channel_dict[str(channel2)]['wavelength']) + \
    ##!#    ' μm - ', GOES_dict['data'][date_idx,ch_idx2,idx1])
    ##!#print("Point Orange:\n")
    ##!#print("    " + str(goes_channel_dict[str(channel1)]['wavelength']) + \
    ##!#    ' μm - ', GOES_dict['data'][date_idx,ch_idx1,idx2])
    ##!#print("    " + str(goes_channel_dict[str(channel2)]['wavelength']) + \
    ##!#    ' μm - ', GOES_dict['data'][date_idx,ch_idx2,idx2])

    #lns = ln1+ln2+ln3
    #labs = [l.get_label() for l in lns]
    #ax1.legend(fontsize = font_size)
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc=0, fontsize = font_size)


    # Plot the GOES 0.64 micron image
    # -------------------------------
    plot_GOES_satpy(date_str, channel1, \
        ax = ax2, var = var, crs = crs, \
        lons = lons, lats = lats, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = goes_channel_dict[str(channel1)]['limits'][0], \
        vmax = goes_channel_dict[str(channel1)]['limits'][1], \
        ptitle = '', plabel = plabel, \
        #vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)

    # Plot the currently-analyzed points
    # ----------------------------------
    point_size = 5
    ax2.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax2.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)
    ax2.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax2.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)

    plot_figure_text(ax2, 'GOES-17 ' + \
        str(goes_channel_dict[str(channel1)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    ax2.set_title(dt_date_str.strftime(\
        '%Y-%m-%d %H:%MZ'), fontsize = labelsize)

    plot_subplot_label(ax1, '(a)', backgroundcolor = 'white', \
        fontsize = font_size + 2, location = 'upper_right')
    plot_subplot_label(ax2, '(b)', backgroundcolor = 'white', \
        fontsize = font_size + 2)

    fig.autofmt_xdate()
    fig.tight_layout()

    if(save):
        outname = save_dir + 'goes_time_series_channel_comp_' + \
        GOES_dict['ptype'] + '_ch' + \
            str(channel1) + '_ch' + str(channel2)
        if(ch_idx3 is not None):
            outname = outname + '_ch' + str(channel3)
        outname = outname + '_pt'+ str(idx1) + '_pt' + str(idx2)\
             + '_' + date_str + \
            '.png'
            #GOES_dict['dt_dates'][time_idx].strftime('%Y%m%d%H%M') + \
            #'.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

# This is written to compare time series of any two desired channels
# (as long as they're in the GOES_dict structure)
def plot_GOES_time_series_channel_comp_2loc(GOES_dict1, GOES_dict2, \
        ch_idx1, ch_idx2, idx1, idx2, idx3, idx4, \
        ch_idx3 = None, date_idx = 20, save_dir = './', \
        save = False):

    # Make the plot for panels 1 and 2
    # --------------------------------

    channel1 = int(GOES_dict1['channels'][ch_idx1])
    channel2 = int(GOES_dict1['channels'][ch_idx2])
    if(ch_idx3 is not None):
        channel3 = int(GOES_dict1['channels'][ch_idx3])
        figsize = (11, 7)
    else:
        figsize = (9, 7)
        

    dt_date_str = GOES_dict1['dt_dates'][date_idx]
    date_str = dt_date_str.strftime('%Y%m%d%H%M')
    var, crs, lons, lats, lat_lims, lon_lims, plabel = \
        read_GOES_satpy(date_str, channel1)

    plt.close('all')
    fig = plt.figure(figsize = figsize)
    gs = fig.add_gridspec(nrows = 2, ncols = 3)
    ax2  = fig.add_subplot(gs[0,2], projection = crs)   # true color    
    ax1  = fig.add_subplot(gs[0,0:2]) # Ch 1
    ax4  = fig.add_subplot(gs[1,2], projection = crs)   # true color    
    ax3  = fig.add_subplot(gs[1,0:2]) # Ch 1

    # Plot the two channel data for the first point
    # ---------------------------------------------
    ln11 = ax1.plot(GOES_dict1['dt_dates'], GOES_dict1['data'][:,ch_idx1,idx1], \
        label = str(goes_channel_dict[\
        str(channel1)]['wavelength']) + \
        ' μm', color = 'tab:blue')
    ln21 = ax1.plot(GOES_dict1['dt_dates'], GOES_dict1['data'][:,ch_idx1,idx2], \
        label = str(goes_channel_dict[\
        str(channel1)]['wavelength']) + \
        ' μm', color = 'tab:orange')

    # Plot the two channel data for the second point
    # ----------------------------------------------
    ln12 = ax1.plot(GOES_dict1['dt_dates'], GOES_dict1['data'][:,ch_idx2,idx1], \
        label = str(goes_channel_dict[\
        str(channel2)]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:blue')
    ln22 = ax1.plot(GOES_dict1['dt_dates'], GOES_dict1['data'][:,ch_idx2,idx2], \
        label = str(goes_channel_dict[\
        str(channel2)]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:orange')
    ax1.axvline(dt_date_str, color = 'black',\
        linestyle = ':')

    lns = ln11 + ln21 + ln12 + ln22 

    if(ch_idx3 is not None):
        ax12 = ax1.twinx()
        ln31 = ax12.plot(GOES_dict1['dt_dates'], GOES_dict1['data'][:,ch_idx3,idx1], \
            label = str(goes_channel_dict[\
            str(channel3)]['wavelength']) + \
            ' μm', linestyle = ':', color = 'tab:blue')
        ln32 = ax12.plot(GOES_dict1['dt_dates'], GOES_dict1['data'][:,ch_idx3,idx2], \
            label = str(goes_channel_dict[\
            str(channel3)]['wavelength']) + \
            ' μm', linestyle = ':', color = 'tab:orange')
        ax12.set_ylabel('Brightness Temperature [K]')

        lns = lns + ln31 + ln32 
        print(len(lns))

    labelsize = 10
    font_size = 8
    ax1.set_ylabel(plabel.replace('_',' '))
    ax1.grid()
    ax1.xaxis.set_major_formatter(DateFormatter('%m/%d\n%H:%MZ'))
    ax1.tick_params(axis="x", labelsize = font_size)

    # Print out the values of each channel at each point
    # --------------------------------------------------
    print(str(goes_channel_dict[str(channel1)]['wavelength']) + \
        ' μm\n')
    print("    Point blue   - ", GOES_dict1['data'][date_idx,ch_idx1,idx1])
    print("    Point orange - ", GOES_dict1['data'][date_idx,ch_idx1,idx2])
    print(str(goes_channel_dict[str(channel2)]['wavelength']) + \
        ' μm\n')
    print("    Point blue   - ", GOES_dict1['data'][date_idx,ch_idx2,idx1])
    print("    Point orange - ", GOES_dict1['data'][date_idx,ch_idx2,idx2])
    if(ch_idx3 is not None):
        print(str(goes_channel_dict[str(channel3)]['wavelength']) + \
            ' μm\n')
        print("    Point blue   - ", GOES_dict1['data'][date_idx,ch_idx3,idx1])
        print("    Point orange - ", GOES_dict1['data'][date_idx,ch_idx3,idx2])

    ##!#print("Point Blue:\n")
    ##!#print("    " + str(goes_channel_dict[str(channel1)]['wavelength']) + \
    ##!#    ' μm - ', GOES_dict1['data'][date_idx,ch_idx1,idx1])
    ##!#print("    " + str(goes_channel_dict[str(channel2)]['wavelength']) + \
    ##!#    ' μm - ', GOES_dict1['data'][date_idx,ch_idx2,idx1])
    ##!#print("Point Orange:\n")
    ##!#print("    " + str(goes_channel_dict[str(channel1)]['wavelength']) + \
    ##!#    ' μm - ', GOES_dict1['data'][date_idx,ch_idx1,idx2])
    ##!#print("    " + str(goes_channel_dict[str(channel2)]['wavelength']) + \
    ##!#    ' μm - ', GOES_dict1['data'][date_idx,ch_idx2,idx2])

    #lns = ln1+ln2+ln3
    #labs = [l.get_label() for l in lns]
    #ax1.legend(fontsize = font_size)
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc=0, fontsize = font_size)


    # Plot the GOES 0.64 micron image
    # -------------------------------
    plot_GOES_satpy(date_str, channel1, \
        ax = ax2, var = var, crs = crs, \
        lons = lons, lats = lats, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = goes_channel_dict[str(channel1)]['limits'][0], \
        vmax = goes_channel_dict[str(channel1)]['limits'][1], \
        ptitle = '', plabel = plabel, \
        #vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)

    # Plot the currently-analyzed points
    # ----------------------------------
    point_size = 5
    ax2.plot(GOES_dict1['plon'][idx1], GOES_dict1['plat'][idx1], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax2.plot(GOES_dict1['plon'][idx1], GOES_dict1['plat'][idx1], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)
    ax2.plot(GOES_dict1['plon'][idx2], GOES_dict1['plat'][idx2], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax2.plot(GOES_dict1['plon'][idx2], GOES_dict1['plat'][idx2], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)

    plot_figure_text(ax2, 'GOES-17 ' + \
        str(goes_channel_dict[str(channel1)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    ax2.set_title(dt_date_str.strftime(\
        '%Y-%m-%d %H:%MZ'), fontsize = labelsize)

    # Make plot for panels 3 and 4
    # ----------------------------
    # Plot the two channel data for the first point
    # ---------------------------------------------
    ln31 = ax3.plot(GOES_dict2['dt_dates'], GOES_dict2['data'][:,ch_idx1,idx3], \
        label = str(goes_channel_dict[\
        str(channel1)]['wavelength']) + \
        ' μm', color = 'tab:blue')
    ln41 = ax3.plot(GOES_dict2['dt_dates'], GOES_dict2['data'][:,ch_idx1,idx4], \
        label = str(goes_channel_dict[\
        str(channel1)]['wavelength']) + \
        ' μm', color = 'tab:orange')

    # Plot the two channel data for the second point
    # ----------------------------------------------
    ln32 = ax3.plot(GOES_dict2['dt_dates'], GOES_dict2['data'][:,ch_idx2,idx3], \
        label = str(goes_channel_dict[\
        str(channel2)]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:blue')
    ln42 = ax3.plot(GOES_dict2['dt_dates'], GOES_dict2['data'][:,ch_idx2,idx4], \
        label = str(goes_channel_dict[\
        str(channel2)]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:orange')
    ax3.axvline(dt_date_str, color = 'black',\
        linestyle = ':')

    lns = ln11 + ln21 + ln12 + ln22 

    if(ch_idx3 is not None):
        ax32 = ax3.twinx()
        ln51 = ax32.plot(GOES_dict2['dt_dates'], GOES_dict2['data'][:,ch_idx3,idx3], \
            label = str(goes_channel_dict[\
            str(channel3)]['wavelength']) + \
            ' μm', linestyle = ':', color = 'tab:blue')
        ln52 = ax32.plot(GOES_dict2['dt_dates'], GOES_dict2['data'][:,ch_idx3,idx4], \
            label = str(goes_channel_dict[\
            str(channel3)]['wavelength']) + \
            ' μm', linestyle = ':', color = 'tab:orange')
        ax32.set_ylabel('Brightness Temperature [K]')

        lns = lns + ln31 + ln32 
        print(len(lns))

    labelsize = 10
    font_size = 8
    ax3.set_ylabel(plabel.replace('_',' '))
    ax3.grid()
    ax3.xaxis.set_major_formatter(DateFormatter('%m/%d\n%H:%MZ'))
    ax3.tick_params(axis="x", labelsize = font_size)

    # Print out the values of each channel at each point
    # --------------------------------------------------
    print(str(goes_channel_dict[str(channel1)]['wavelength']) + \
        ' μm\n')
    print("    Point blue   - ", GOES_dict2['data'][date_idx,ch_idx1,idx3])
    print("    Point orange - ", GOES_dict2['data'][date_idx,ch_idx1,idx4])
    print(str(goes_channel_dict[str(channel2)]['wavelength']) + \
        ' μm\n')
    print("    Point blue   - ", GOES_dict2['data'][date_idx,ch_idx2,idx3])
    print("    Point orange - ", GOES_dict2['data'][date_idx,ch_idx2,idx4])
    if(ch_idx3 is not None):
        print(str(goes_channel_dict[str(channel3)]['wavelength']) + \
            ' μm\n')
        print("    Point blue   - ", GOES_dict2['data'][date_idx,ch_idx3,idx3])
        print("    Point orange - ", GOES_dict2['data'][date_idx,ch_idx3,idx4])

    ##!#print("Point Blue:\n")
    ##!#print("    " + str(goes_channel_dict[str(channel1)]['wavelength']) + \
    ##!#    ' μm - ', GOES_dict2['data'][date_idx,ch_idx1,idx1])
    ##!#print("    " + str(goes_channel_dict[str(channel2)]['wavelength']) + \
    ##!#    ' μm - ', GOES_dict2['data'][date_idx,ch_idx2,idx1])
    ##!#print("Point Orange:\n")
    ##!#print("    " + str(goes_channel_dict[str(channel1)]['wavelength']) + \
    ##!#    ' μm - ', GOES_dict2['data'][date_idx,ch_idx1,idx2])
    ##!#print("    " + str(goes_channel_dict[str(channel2)]['wavelength']) + \
    ##!#    ' μm - ', GOES_dict2['data'][date_idx,ch_idx2,idx2])

    #lns = ln1+ln2+ln3
    #labs = [l.get_label() for l in lns]
    #ax1.legend(fontsize = font_size)
    labs = [l.get_label() for l in lns]
    ax3.legend(lns, labs, loc=0, fontsize = font_size)


    # Plot the GOES 0.64 micron image
    # -------------------------------
    plot_GOES_satpy(date_str, channel1, \
        ax = ax4, var = var, crs = crs, \
        lons = lons, lats = lats, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = goes_channel_dict[str(channel1)]['limits'][0], \
        vmax = goes_channel_dict[str(channel1)]['limits'][1], \
        ptitle = '', plabel = plabel, \
        #vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)

    # Plot the currently-analyzed points
    # ----------------------------------
    point_size = 5
    ax4.plot(GOES_dict2['plon'][idx3], GOES_dict2['plat'][idx3], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax4.plot(GOES_dict2['plon'][idx3], GOES_dict2['plat'][idx3], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)
    ax4.plot(GOES_dict2['plon'][idx4], GOES_dict2['plat'][idx4], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax4.plot(GOES_dict2['plon'][idx4], GOES_dict2['plat'][idx4], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)

    plot_figure_text(ax4, 'GOES-17 ' + \
        str(goes_channel_dict[str(channel1)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    ax4.set_title(dt_date_str.strftime(\
        '%Y-%m-%d %H:%MZ'), fontsize = labelsize)

    plot_subplot_label(ax1, '(a)', backgroundcolor = 'white', \
        fontsize = font_size + 2, location = 'upper_right')
    plot_subplot_label(ax2, '(b)', backgroundcolor = 'white', \
        fontsize = font_size + 2)
    plot_subplot_label(ax3, '(c)', backgroundcolor = 'white', \
        fontsize = font_size + 2, location = 'upper_right')
    plot_subplot_label(ax4, '(d)', backgroundcolor = 'white', \
        fontsize = font_size + 2)

    fig.autofmt_xdate()
    fig.tight_layout()

    if(save):
        outname = save_dir + 'goes_time_series_channel_comp_2loc_' + \
            GOES_dict1['ptype'] + '_' + GOES_dict2['ptype'] +  '_ch' + \
            str(channel1) + '_ch' + str(channel2)
        if(ch_idx3 is not None):
            outname = outname + '_ch' + str(channel3)
        outname = outname + '_pt'+ str(idx1) + '_pt' + str(idx2)\
             + '_' + date_str + \
            '.png'
            #GOES_dict['dt_dates'][time_idx].strftime('%Y%m%d%H%M') + \
            #'.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

# This is written to compare time series of any two desired channels
# (as long as they're in the GOES_dict structure). Compares data
# from the same point(s) along the same cross section but from different
# dates.
def plot_GOES_time_series_channel_comp_2loc(GOES_dict1, GOES_dict2, \
        ch_idx1, ch_idx2, idx1, ch_idx3 = None, \
        date_idx = 20, save_dir = './', save = False):

    # Make the plot for panels 1 and 2
    # --------------------------------

    channel1 = int(GOES_dict1['channels'][ch_idx1])
    channel2 = int(GOES_dict1['channels'][ch_idx2])
    if(ch_idx3 is not None):
        channel3 = int(GOES_dict1['channels'][ch_idx3])
        figsize = (8, 7)
    else:
        figsize = (9, 7)
       

    dt_date_str1 = GOES_dict1['dt_dates'][date_idx]
    dt_date_str2 = GOES_dict2['dt_dates'][date_idx + 1]
    date_str1 = dt_date_str1.strftime('%Y%m%d%H%M')
    date_str2 = dt_date_str2.strftime('%Y%m%d%H%M')
    var1, crs1, lons1, lats1, lat_lims1, lon_lims1, plabel1 = \
        read_GOES_satpy(date_str1, channel2)
    var2, crs2, lons2, lats2, lat_lims2, lon_lims2, plabel2 = \
        read_GOES_satpy(date_str2, channel2)

    # Figure out relative dates to plot the same plot
    # -----------------------------------------------
    base_date1 = datetime.strptime('200007130000','%Y%m%d%H%M')
    base_date2 = datetime.strptime('200007200000','%Y%m%d%H%M')
    rel_date1 = GOES_dict1['dt_dates'] - base_date1
    rel_date2 = GOES_dict2['dt_dates'] - base_date2

    plt.close('all')
    fig = plt.figure(figsize = figsize)
    gs = fig.add_gridspec(nrows = 2, ncols = 2)
    ax2  = fig.add_subplot(gs[0,0], projection = crs1)   # true color    
    ax3  = fig.add_subplot(gs[0,1], projection = crs2)   # true color    
    ax1  = fig.add_subplot(gs[1,:]) # Ch 1

    # Plot the two channel data for the first date 
    # ---------------------------------------------
    ln11 = ax1.plot(GOES_dict1['dt_dates'], GOES_dict1['data'][:,ch_idx1,idx1], \
    #ln11 = ax1.plot(rel_date1, GOES_dict1['data'][:,ch_idx1,idx1], \
        label = str(goes_channel_dict[\
        str(channel1)]['wavelength']) + \
        ' μm', color = 'tab:blue')
    ax112 = ax1.twiny()
    ln12 = ax112.plot(GOES_dict2['dt_dates'], GOES_dict2['data'][:,ch_idx1,idx1], \
    #ln12 = ax1.plot(rel_date2, GOES_dict2['data'][:,ch_idx1,idx1], \
        label = str(goes_channel_dict[\
        str(channel1)]['wavelength']) + \
        ' μm', color = 'tab:orange')


    # Plot the two channel data for the second point
    # ----------------------------------------------
    ln21 = ax1.plot(GOES_dict1['dt_dates'], GOES_dict1['data'][:,ch_idx2,idx1], \
        label = str(goes_channel_dict[\
        str(channel2)]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:blue')
    ln22 = ax112.plot(GOES_dict2['dt_dates'], GOES_dict2['data'][:,ch_idx2,idx1], \
        label = str(goes_channel_dict[\
        str(channel1)]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:orange')
    ax1.axvline(dt_date_str1, color = 'black',\
        linestyle = ':')

    lns = ln11 + ln12 + ln21 + ln22

    labelsize = 10
    font_size = 10
    if(ch_idx3 is not None):
        ax12 = ax1.twinx()
        ax212 = ax12.twiny()
        ln31 = ax12.plot(GOES_dict1['dt_dates'], GOES_dict1['data'][:,ch_idx3,idx1], \
            label = str(goes_channel_dict[\
            str(channel3)]['wavelength']) + \
            ' μm', linestyle = ':', color = 'tab:blue')
        ln32 = ax212.plot(GOES_dict2['dt_dates'], GOES_dict2['data'][:,ch_idx3,idx1], \
            label = str(goes_channel_dict[\
            str(channel3)]['wavelength']) + \
            ' μm', linestyle = ':', color = 'tab:orange')
        ax12.set_ylabel('Brightness Temperature [K]', weight = 'bold')
        ax212.xaxis.set_major_formatter(DateFormatter(''))
        ax212.tick_params(axis="x", labelsize = font_size)

        lns = lns + ln31 + ln32
        print(len(lns))

    ax1.set_ylabel(plabel1.replace('_',' '), weight = 'bold')
    ax1.grid()
    ax1.xaxis.set_major_formatter(DateFormatter('%H:%MZ'))
    ax1.tick_params(axis="x", labelsize = font_size)
    ax112.xaxis.set_major_formatter(DateFormatter(''))
    ax112.tick_params(axis="x", labelsize = font_size)

    # Print out the values of each channel at each point
    # --------------------------------------------------
    print(str(goes_channel_dict[str(channel1)]['wavelength']) + \
        ' μm\n')
    print("    Date blue   - ", GOES_dict1['data'][date_idx,ch_idx1,idx1])
    print("    Date orange - ", GOES_dict2['data'][date_idx,ch_idx1,idx1])
    print(str(goes_channel_dict[str(channel2)]['wavelength']) + \
        ' μm\n')
    print("    Date blue   - ", GOES_dict1['data'][date_idx,ch_idx2,idx1])
    print("    Date orange - ", GOES_dict2['data'][date_idx,ch_idx2,idx1])
    if(ch_idx3 is not None):
        print(str(goes_channel_dict[str(channel3)]['wavelength']) + \
            ' μm\n')
        print("    Date blue   - ", GOES_dict1['data'][date_idx,ch_idx3,idx1])
        print("    Date orange - ", GOES_dict2['data'][date_idx,ch_idx3,idx1])

    #lns = ln1+ln2+ln3
    #labs = [l.get_label() for l in lns]
    #ax1.legend(fontsize = font_size)
    ##!#labs = [l.get_label() for l in lns]
    ##!#ax1.legend(lns, labs, loc=2, fontsize = font_size)

    custom_lines = [Line2D([0], [0], color='k'),
                    Line2D([0], [0], color='k', linestyle = '--'),
                    Line2D([0], [0], color='k', linestyle = ':'), 
                    Line2D([0], [0], color = 'tab:blue'), 
                    Line2D([0], [0], color = 'tab:orange')]

    ax1.legend(custom_lines, ['0.64 μm', '2.25 μm', '10.35 μm', \
            'Before', 'During'],\
        fontsize = font_size, loc = 2)


    # Plot the GOES 0.64 micron image for date 1
    # ------------------------------------------
    plot_GOES_satpy(date_str1, channel1, \
        ax = ax2, var = var1, crs = crs1, \
        lons = lons1, lats = lats1, lat_lims = lat_lims1, lon_lims = lon_lims1, \
        vmin = goes_channel_dict[str(channel2)]['limits'][0], \
        vmax = goes_channel_dict[str(channel2)]['limits'][1] - 20, \
        ptitle = '', plabel = plabel1, \
        #vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)

    # Plot the point for date 1
    # -------------------------
    point_size = 5
    ax2.plot(GOES_dict1['plon'][idx1], GOES_dict1['plat'][idx1], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax2.plot(GOES_dict1['plon'][idx1], GOES_dict1['plat'][idx1], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)

    plot_figure_text(ax2, 'Before Dixie Fire', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right', location = 'upper_right')
    plot_figure_text(ax2, 'GOES-17 ' + \
        str(goes_channel_dict[str(channel2)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    ax2.set_title(dt_date_str1.strftime(\
        '%Y-%m-%d %H:%MZ'), fontsize = labelsize)

    # Plot the GOES 0.64 micron image for date 2
    # ------------------------------------------
    plot_GOES_satpy(date_str2, channel1, \
        ax = ax3, var = var2, crs = crs2, \
        lons = lons2, lats = lats2, lat_lims = lat_lims2, lon_lims = lon_lims2, \
        vmin = goes_channel_dict[str(channel2)]['limits'][0], \
        vmax = goes_channel_dict[str(channel2)]['limits'][1] - 20, \
        ptitle = '', plabel = plabel2, \
        #vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)

    # Plot the currently-analyzed points
    # ----------------------------------
    point_size = 5
    ax3.plot(GOES_dict2['plon'][idx1], GOES_dict2['plat'][idx1], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax3.plot(GOES_dict2['plon'][idx1], GOES_dict2['plat'][idx1], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs, color = 'tab:orange')

    plot_figure_text(ax3, 'During Dixie Fire', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right', location = 'upper_right')
    plot_figure_text(ax3, 'GOES-17 ' + \
        str(goes_channel_dict[str(channel2)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    ax3.set_title(dt_date_str2.strftime(\
        '%Y-%m-%d %H:%MZ'), fontsize = labelsize)

    plot_subplot_label(ax1, '(c)', backgroundcolor = 'white', \
        fontsize = font_size + 2, location = 'upper_right')
    plot_subplot_label(ax2, '(a)', backgroundcolor = 'white', \
        fontsize = font_size + 2)
    plot_subplot_label(ax3, '(b)', backgroundcolor = 'white', \
        fontsize = font_size + 2)

    fig.autofmt_xdate()
    fig.tight_layout()

    if(save):
        outname = save_dir + 'goes_time_series_channel_comp_2date_' + \
            GOES_dict1['ptype'] + '_ch' + \
            str(channel1) + '_ch' + str(channel2)
        if(ch_idx3 is not None):
            outname = outname + '_ch' + str(channel3)
        outname = outname + '_pt'+ str(idx1) + '_' + date_str1 + \
            '_' + date_str2 + '.png'
            #GOES_dict['dt_dates'][time_idx].strftime('%Y%m%d%H%M') + \
            #'.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

def plot_GOES_time_series_mesh(GOES_dict, ch_idx1 = 1, \
        ch_idx2 = 0, save_dir = './', \
        date_idx = 23, sigma = 0.8, map_cntr = True, \
        show_points = False, save = False):

    channel1 = int(GOES_dict['channels'][ch_idx1])
    channel2 = int(GOES_dict['channels'][ch_idx2])

    #dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
    dt_date_str = GOES_dict['dt_dates'][date_idx]
    date_str = dt_date_str.strftime('%Y%m%d%H%M')
    var, crs, lons, lats, lat_lims, lon_lims, plabel = \
        read_GOES_satpy(date_str, channel1)
    if(map_cntr):
        var2, crs2, lons2, lats2, lat_lims2, lon_lims2, plabel2 = \
            read_GOES_satpy(date_str, channel2)

    plt.close('all')
    fig = plt.figure(figsize = (11, 3.5))
    gs = fig.add_gridspec(nrows = 1, ncols = 4)
    ax2  = fig.add_subplot(gs[3], projection = crs)   # true color    
    ax1  = fig.add_subplot(gs[0:3]) # Ch 1

    if(channel1 > 6):
        cmap1 = 'plasma'
    else:
        cmap1 = 'Greys_r'
    if(channel2 > 6):
        cmap2 = 'plasma'
        levels = np.arange(283, 323, 5)
    else:
        cmap2 = 'Greys_r'
        levels = np.arange(0, goes_channel_dict[str(channel1)\
            ]['limits'][1], 2)

    mesh = ax1.pcolormesh(GOES_dict['dt_dates'], \
        GOES_dict['plat'][:],\
        GOES_dict['data'][:,ch_idx1,:].T, cmap = cmap1, \
        vmin = goes_channel_dict[str(channel1)\
            ]['limits'][0], \
        vmax = goes_channel_dict[str(channel1)\
            ]['limits'][1], \
        shading = 'auto')
   
    ax1.set_ylabel('Latitude') 
    ax1.xaxis.set_major_formatter(DateFormatter('%m/%d\n%H:%MZ'))
    ax1.tick_params(axis="x", labelsize = 9)

    cbar = plt.colorbar(mesh, ax = ax1, pad = 0.03, fraction = 0.052, \
        extend = 'both', label = plabel)

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
        vmin = goes_channel_dict[str(channel2)]['limits'][0], \
        vmax = goes_channel_dict[str(channel2)]['limits'][1], \
        levels = levels)
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
    plot_GOES_satpy(date_str, channel1, \
        ax = ax2, var = var, crs = crs, \
        lons = lons, lats = lats, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = goes_channel_dict[str(channel1)]['limits'][0], \
        vmax = goes_channel_dict[str(channel1)]['limits'][1], \
        ptitle = '', plabel = plabel, \
        #vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)

    if(map_cntr):
        # Smooth the other data too
        smooth_data2 = gaussian_filter(var2, \
            sigma = 1.0)
        cntr2 = ax2.contour(lons2, lats2, smooth_data2, \
            levels, transform = datacrs, cmap = 'plasma', \
            alpha = 0.5)
        ax2.clabel(cntr2, cntr2.levels, inline=True, fontsize=8)

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

    # Show the locations of the points on the map and the 
    # time series, if desired
    # ----------------------------------------------------
    if(show_points):
        for tlat, tlon in zip(GOES_dict['plat'][1::2], GOES_dict['plon'][1::2]):
            print(tlat, tlon)
            ax1.plot(GOES_dict['dt_dates'][date_idx], tlat, \
                    linewidth=2, markersize = 5, marker='^',
                    color = 'black')
            ax1.plot(GOES_dict['dt_dates'][date_idx], tlat, \
                    linewidth=2, markersize = 3, marker='^')
            ax2.plot(tlon, tlat, linewidth=2, markersize = 5, marker='^',
                    color = 'black', transform=datacrs)
            ax2.plot(tlon, tlat, linewidth=2, markersize = 3, marker='^',
                    transform=datacrs)
            

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
    all_files = np.array(glob(home_dir + '/data/GOES/goes17_abi/*.nc'))
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

def read_GOES_time_series_NCDF(file_name):

    # Set up the output dictionary
    # ----------------------------
    GOES_dict = {}

    # Open the netCDF file
    # --------------------
    nc = Dataset(file_name, 'r')

    # Convert the seconds in the file to datetime objects
    # ---------------------------------------------------
    base_date = datetime.strptime(nc.base_date, '%Y%m%d%H%M')
    dates = np.array([(base_date + timedelta(seconds = int(ttime))) \
        for ttime in nc['time'][:]])

    # Insert the data into the dictionary
    # -----------------------------------
    GOES_dict['data']      = nc['data'][:,:,:].data
    GOES_dict['goes_lats'] = nc['goes_lat'][:,:,:].data
    GOES_dict['goes_lons'] = nc['goes_lon'][:,:,:].data
    GOES_dict['dt_dates']  = dates
    GOES_dict['channels']  = nc['channel_num'][:].data
    GOES_dict['plat']      = nc['point_lat'][:].data
    GOES_dict['plon']      = nc['point_lon'][:].data
    GOES_dict['ptype']     = nc.ptype

    nc.close()

    return GOES_dict
 
def write_GOES_time_series_NCDF(GOES_dict, save_dir = './'):

    file_name_start = save_dir + 'goes_cross_data_' + GOES_dict['ptype'] + \
        '_' + \
        GOES_dict['dt_dates'][0].strftime('%Y%m%d%H%M') + '_' + \
        GOES_dict['dt_dates'][-1].strftime('%Y%m%d%H%M')
    file_name = file_name_start + '.nc'

    # Check if the file already exists
    # --------------------------------
    while(os.path.exists(file_name)):
        print("WARNING: file",file_name," already exists. Making new version")

        # Check if this file exists
        # -------------------------
        if(file_name.strip().split('_')[-1].split('.')[0][0] == 'v'):
            version_num = int(file_name.strip().split('_')[-1].split('.')[0][1:])
            file_name = file_name_start + '_v' + str(version_num + 1) + '.nc'

        else:
            file_name = save_dir + file_name_start + '_v2.nc'

    print(file_name) 

    # Create a new netCDF dataset to write to
    # --------------------------------------- 
    nc = Dataset(file_name,'w',format='NETCDF4')
  
    # Dimensions = datetime, channel, point
    # Create the sizes of each dimension in the file. 
    # -----------------------------------------------
    num_point = len(GOES_dict['plat'])
    num_ch    = len(GOES_dict['channels'])
    num_time  = len(GOES_dict['dt_dates'])
  
    # For the dates, calculate the number of seconds between 
    # each date and a reference date, which is January 1st, 2000
    # at 00:00 UTC.
    # ---------------------------------------------------------- 
    base_date = datetime(year=2000,month=1,day=1,hour = 0, minute = 0,\
        second = 0) 
    times = np.array([(tdate - base_date).total_seconds() \
        for tdate in GOES_dict['dt_dates']])
  
    # Add an attribute to the file to contain the cross section
    # location
    # ---------------------------------------------------------
    nc.ptype = GOES_dict['ptype']
    nc.base_date = base_date.strftime('%Y%m%d%H%M')
 
    # Use the dimension size variables to actually create dimensions in 
    # the file.
    # ----------------------------------------------------------------- 
    n_point = nc.createDimension('n_point',num_point)
    n_ch    = nc.createDimension('n_ch',   num_ch)
    n_time  = nc.createDimension('n_time', num_time)

    # Create variables for the three dimensions. Note that since these
    # are variables, they are still given 'dimensions' using 'n_point',
    # 'n_ch', and 'n_time'.
    # ----------------------------------------------------------------
    TIME = nc.createVariable('time','f8',('n_time'))
    TIME.description = 'Seconds since 00:00 UTC 1 January 2000'
    CHANNEL = nc.createVariable('channel_num','i2',('n_ch'))
    CHANNEL.description = 'GOES channel number'
    POINT = nc.createVariable('cross_point','i2',('n_point'))
    POINT.description = 'Cross section point index'
   
    # Create variables
    # ----------------
    DATA = nc.createVariable('data','f4',('n_time','n_ch','n_point'))
    DATA.description = 'GOES reflectance and/or brightness temperature at each point'
    GOES_LAT = nc.createVariable('goes_lat','f4',('n_time','n_ch','n_point'))
    GOES_LAT.description = 'GOES latitude of each point for each channel'
    GOES_LON = nc.createVariable('goes_lon','f4',('n_time','n_ch','n_point'))
    GOES_LON.description = 'GOES longitude of each point for each channel'
    POINT_LAT = nc.createVariable('point_lat','f4',('n_point'))
    POINT_LAT.description = 'Latitude of each cross section point'
    POINT_LON = nc.createVariable('point_lon','f4',('n_point'))
    POINT_LON.description = 'Longitude of each cross section point'

    # Fill in dimension variables one-by-one.
    # NOTE: not sure if you can insert the entire data array into the
    # dimension (for example, doing MONTH = times), so this could be
    # something for you to try. Might make this faster if it works
    # ---------------------------------------------------------------
    TIME[:]          = times[:]
    CHANNEL[:]       = GOES_dict['channels'][:]
    POINT[:]         = np.arange(num_point)
    DATA[:,:,:]      = GOES_dict['data'][:,:,:]
    GOES_LAT[:,:,:]  = GOES_dict['goes_lats'][:,:,:]
    GOES_LON[:,:,:]  = GOES_dict['goes_lons'][:,:,:]
    POINT_LAT[:]     = GOES_dict['plat'][:]
    POINT_LON[:]     = GOES_dict['plon'][:]
    ##!#for i in range(num_time):
    ##!#    MONTH[i] = times[i]
    ##!#for i in range(num_lat):
    ##!#    for j in range(num_lon):
    ##!#        LAT[i,j]=lat_ranges[i]
    ##!#        LON[i,j]=lon_ranges[j]

    ##!## Fill in actual variables
    ##!## I have a bunch of extra stuff in here for handling the dictionary
    ##!## keys, which you will likely not need for your data
    ##!## ------------------------------------------------------------------
    ##!#for i in range(num_lat):
    ##!#    print(lat_ranges[i])
    ##!#    for j in range(num_lon):
    ##!#        dictkey = (str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j])))
    ##!#        if(dictkey not in OMI_data):
    ##!#            # Insert missing values for AI and count
    ##!#            AI[:,i,j] = [-999.9 for m in range(num_time)]
    ##!#            OB_COUNT[:,i,j] = [-99 for m in range(num_time)]
    ##!#        else:
    ##!#            for m in range(num_time):
    ##!#                timekey = testkeys[m]
    ##!#                if(timekey not in OMI_data[dictkey]):
    ##!#                    AI[m,i,j] = -999.9
    ##!#                    OB_COUNT[m,i,j] = -99
    ##!#                else:
    ##!#                    AI[m,i,j] = OMI_data[dictkey][timekey]['avg']
    ##!#                    OB_COUNT[m,i,j] = OMI_data[dictkey][timekey]['#_obs']

    # Close the netCDF file, which actually saves it.
    # -----------------------------------------------
    nc.close()
    print("Saved file",file_name)

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
