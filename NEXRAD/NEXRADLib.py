"""
  NAME:

  PURPOSE:
    Read and plot NEXRAD data using pyart. Also contains functions for 
        automatically downloading NEXRAD data from the Amazon bucket. 

    NOTE: must be running a Conda environment with satpy installed. 

  SYNTAX:
    # = = = = = = = = = = =
    # Notes for downloading
    # = = = = = = = = = = =
    Before running the downloading scripts, the main of which is 
        auto_NEXRAD_download(), you MUST have the Amazon AWS CLI 
        interface installed on your Linux terminal. To do this, 
        you must have sudo permissions on your computer. If you
        have sudo permissions, run the following in the command
        line:
    
    $ sudo apt install awscli

    Once this is installed, and if you have a Conda environment
        installed that has satpy and datetime (and the other essential
        packages listed in the "import" section), you may download
        NEXRAD data automatically using auto_NEXRAD_download. To run
        the function, use the following syntax:

    >>> from NEXRADLib_lite import *
    >>> begin_date = '202107201200' # begin_date in YYYYMMDDHHMM format
    >>> end_date   = '202107210000' # end_date in YYYYMMDDHHMM format
    >>> interval   = 30  # this is the number of minutes between each file
    >>> auto_NEXRAD_download(begin_date, end_date, interval)
  
    # = = = = = = = = = = 
    # Notes for plotting
    # = = = = = = = = = =
    When running in a separate Python script, use the following syntax:
    >>> from NEXRADLib_lite import *
    >>> date_str = '202107202126'  # date_str in YYYYMMDDHHMM format, for whatever the time of the NEXRAD file 
    >>> plot_NEXRAD_satpy(date_str, 2)   # Plots for channel 2 data
    
    If you want to save the image, do:
    >>> plot_NEXRAD_satpy(date_str, 2, save = True)


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
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from metpy.plots import USCOUNTIES
import pyart
from glob import glob
import os
from pathlib import Path

import copy
from netCDF4 import num2date
from pandas import to_datetime
from pyart.core import Radar
from scipy.interpolate import interp2d

sys.path.append('/home/bsorenson/')
from python_lib import *

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# Set up global variables
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
datacrs = ccrs.PlateCarree()
data_dir = '/home/bsorenson/data/NEXRAD/'
debug = False

nexrad_area_dict = {
    "2021-07-13": {
        ##!#'2110': {
        'asos': 'asos_data_20210713.csv',
        #'nexrad': '/home/bsorenson/data/NEXRAD/Aqua/MYD021KM.A2021203.2110.061.2021204155922.hdf',
        #'mdswv': '/home/bsorenson/data/NEXRAD/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
        #'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072210-2021072221.nc',
        #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
        'Lat': [39.5, 42.0],
        'Lon': [-122.0, -119.5],
        ##!#'nexrad_Lat': [39.0, 42.5],
        ##!#'nexrad_Lon': [-123., -119.]
        'nexrad_Lat': [39.5, 42.0],
        'nexrad_Lon': [-122.0, -119.5],
        ##!#}
    },
    "2021-07-14": {
        ##!#'2110': {
        'asos': 'asos_data_20210713.csv',
        #'nexrad': '/home/bsorenson/data/NEXRAD/Aqua/MYD021KM.A2021203.2110.061.2021204155922.hdf',
        #'mdswv': '/home/bsorenson/data/NEXRAD/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
        #'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072210-2021072221.nc',
        #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.22.212.L2.SUBS2RET.v6.0.32.0.G21204140844.hdf'],
        'Lat': [39.5, 42.0],
        'Lon': [-122.0, -119.5],
        ##!#'nexrad_Lat': [39.0, 42.5],
        ##!#'nexrad_Lon': [-123., -119.]
        'nexrad_Lat': [39.5, 42.0],
        'nexrad_Lon': [-122.0, -119.5],
        ##!#}
    },
    "2021-07-20": {
        'asos': 'asos_data_20210720.csv',
        #'nexrad': '/home/bsorenson/data/NEXRAD/Aqua/MYD021KM.A2021201.2125.061.2021202154814.hdf',
        #'ceres': '/home/bsorenson/data/CERES/SSF_Level2/Aqua/CERES_SSF_Aqua-XTRK_Edition4A_Subset_2021072010-2021072021.nc',
        #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.07.20.214.L2.SUBS2RET.v6.0.32.0.G21202153435.hdf'],
        'Lat': [39.5, 42.0],
        'Lon': [-122.0, -119.5],
        'data_lim': {
            1:  [0.05, 0.5],
            31: [270., 330.],
        },
        'nexrad_Lat': [39.5, 42.0],
        'nexrad_Lon': [-122.0, -119.5]
    },
    "2021-07-21": {
        'asos': 'asos_data_20210722_4.csv',
        #'nexrad': '/home/bsorenson/data/NEXRAD/Aqua/MYD021KM.A2021202.2030.061.2021203174050.hdf',
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
        'nexrad_Lat': [39.5, 42.0],
        'nexrad_Lon': [-122.0, -119.5]
    },
    "2021-07-22": {
        'asos': 'asos_data_20210722_4.csv',
        #'nexrad': '/home/bsorenson/data/NEXRAD/Aqua/MYD021KM.A2021203.2110.061.2021204155922.hdf',
        #'mdswv': '/home/bsorenson/data/NEXRAD/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
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
        #'nexrad_Lat': [39.5, 42.0],
        #'nexrad_Lon': [-122.0, -119.5]
        'nexrad_Lat': [39.5, 42.0],
        'nexrad_Lon': [-122.0, -119.5]
    },
    "2021-07-23": {
        'asos': 'asos_data_20210722_2.csv',
        'nexrad': '/home/bsorenson/data/NEXRAD/Aqua/MYD021KM.A2021204.2155.061.2021205153516.hdf',
        #'mdswv': '/home/bsorenson/data/NEXRAD/Aqua/MYD05_L2.A2021203.2110.061.2021204163638.hdf',
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
        #'nexrad_Lat': [39.5, 42.0],
        #'nexrad_Lon': [-122.0, -119.5]
        'nexrad_Lat': [39.5, 42.0],
        'nexrad_Lon': [-122.0, -119.5]
    },
    "2021-08-04": {
        'asos': 'asos_data_20210806.csv',
        #'asos': 'asos_nevada_20210806.csv',
        #'nexrad': '/home/bsorenson/data/NEXRAD/Aqua/MYD021KM.A2021218.2025.061.2021219151802.hdf',
        #'mdswv': '/home/bsorenson/data/NEXRAD/Aqua/MYD05_L2.A2021218.2025.061.2021219152751.hdf',
        'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.04.206.L2.SUBS2RET.v6.0.32.0.G21217152448.hdf',\
                 '/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.04.207.L2.SUBS2RET.v6.0.32.0.G21217152904.hdf'],\
        #'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
        'Lat': [36.0, 39.0],
        'Lon': [-118.0, -114.0],
        'nexrad_Lat': [35.0, 40.0],
        'nexrad_Lon': [-119., -113.]
    },
    "2021-08-05": {
        '2120': {
            'asos': 'asos_data_20210805.csv',
            'nexrad': '/home/bsorenson/data/NEXRAD/Aqua/MYD021KM.A2021217.2120.061.2021218164201.hdf',
            'mdswv': '/home/bsorenson/data/NEXRAD/Aqua/MYD05_L2.A2021217.2120.061.2021218165546.hdf',
            'ceres': '/home/bsorenson/data/CERES/FLASHFlux/Aqua/FLASH_SSF_Aqua_Version4A_Subset_2021080508-2021080521.nc',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.05.214.L2.SUBS2RET.v6.0.32.0.G21218175548.hdf'],
            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0805t2038-o90733_v003-2021m0807t014855.he5',
            'Lat': [36.5, 39.0],
            'Lon': [-118.0, -114.0],
            'nexrad_Lat': [36.0, 39.0],
            'nexrad_Lon': [-118., -114.]
        },
        '2125': {
            'asos': 'asos_california_20210805.csv',
            'nexrad': '/home/bsorenson/data/NEXRAD/Aqua/MYD021KM.A2021217.2125.061.2021218161010.hdf',
            'ceres': '/home/bsorenson/data/CERES/FLASHFlux/Aqua/FLASH_SSF_Aqua_Version4A_Subset_2021080508-2021080521.nc',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.05.214.L2.SUBS2RET.v6.0.32.0.G21218175548.hdf'],
            'Lat': [39.5, 42.0],
            'Lon': [-122.0, -119.0],
            'nexrad_Lat': [39.5, 42.0],
            'nexrad_Lon': [-122., -119.]
        }
    },
    "2021-08-06": {
        '2025': {
            'asos': 'asos_data_20210806.csv',
            #'asos': 'asos_nevada_20210806.csv',
            'nexrad': '/home/bsorenson/data/NEXRAD/Aqua/MYD021KM.A2021218.2025.061.2021219151802.hdf',
            'mdswv': '/home/bsorenson/data/NEXRAD/Aqua/MYD05_L2.A2021218.2025.061.2021219152751.hdf',
            'ceres': '/home/bsorenson/data/CERES/FLASHFlux/Aqua/FLASH_SSF_Aqua_Version4A_Subset_2021080609-2021080620.nc',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.06.204.L2.SUBS2RET.v6.0.32.0.G21219130523.hdf',\
                     '/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.06.205.L2.SUBS2RET.v6.0.32.0.G21219130455.hdf'],
            'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
            'Lat': [36.0, 39.0],
            'Lon': [-118.0, -114.0],
            'nexrad_Lat': [36.0, 39.0],
            'nexrad_Lon': [-118., -114.]
        }
    },
    "2021-08-07": {
        '2110': {
            'asos': 'asos_data_20210806.csv',
            #'asos': 'asos_nevada_20210806.csv',
            'nexrad': '/home/bsorenson/data/NEXRAD/Aqua/MYD021KM.A2021219.2110.061.2021220151612.hdf',
            #'mdswv': '/home/bsorenson/data/NEXRAD/Aqua/MYD05_L2.A2021218.2025.061.2021219152751.hdf',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.07.212.L2.SUBS2RET.v6.0.32.0.G21220123225.hdf'],\
            #'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
            'Lat': [36.0, 39.0],
            'Lon': [-118.0, -114.0],
            'nexrad_Lat': [35.0, 40.0],
            'nexrad_Lon': [-119., -113.]
        }
    },
    "2021-08-08": {
        '2110': {
            'asos': 'asos_data_20210806.csv',
            #'asos': 'asos_nevada_20210806.csv',
            #'nexrad': '/home/bsorenson/data/NEXRAD/Aqua/MYD021KM.A2021218.2025.061.2021219151802.hdf',
            #'mdswv': '/home/bsorenson/data/NEXRAD/Aqua/MYD05_L2.A2021218.2025.061.2021219152751.hdf',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.08.202.L2.SUBS2RET.v6.0.32.0.G21221124420.hdf',\
                     '/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.08.203.L2.SUBS2RET.v6.0.32.0.G21221185932.hdf'],\
            #'omi': '/home/bsorenson/data/OMI/H5_files/OMI-Aura_L2-OMAERUV_2021m0806t1943-o90747_v003-2021m0808t031152.he5',
            'Lat': [36.0, 39.0],
            'Lon': [-118.0, -114.0],
            'nexrad_Lat': [35.0, 40.0],
            'nexrad_Lon': [-119., -113.]
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
            'nexrad': '/home/bsorenson/data/NEXRAD/Aqua/MYD021KM.A2021242.2115.061.2021243183953.hdf',
            'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.30.213.L2.SUBS2RET.v6.0.32.0.G21243151114.hdf'],
            'Lat': [38.0, 40.0],
            'Lon': [-121.0, -118.5]
        }
    },
    "2021-09-01": {
        '2105': {
            'asos': 'asos_data_20210830.csv',
            'nexrad': '/home/bsorenson/data/NEXRAD/Aqua/MYD021KM.A2021244.2105.061.2021245152256.hdf',
            #'airs': ['/home/bsorenson/data/AIRS/Aqua/AIRS.2021.08.30.213.L2.SUBS2RET.v6.0.32.0.G21243151114.hdf'],
            'Lat': [38.0, 42.0],
            'Lon': [-121.5, -118.0],
            'nexrad_Lat': [38.0, 42.0],
            'nexrad_Lon': [-121.5, -118.]
        }
    } 
}

min_dict = {
    'composite_reflectivity': {\
        'ppi': -5, \
        'rhi': -5, \
    },
    'reflectivity': {\
        'ppi': -10,\
        'rhi': -10,\
    },
    'velocity':  {\
        'ppi': -30,\
        'rhi': -30,\
    },
    'differential_reflectivity': {\
        'ppi':-2.5, \
        'rhi':-2.5, \
    },
    'cross_correlation_ratio': {\
        'ppi': 0.5,\
        'rhi': 0.5,\
    }
}
max_dict = {
    'composite_reflectivity': {\
        'ppi': 90, \
        'rhi': 64, \
    },
    'reflectivity': {\
        'ppi': 90,\
        'rhi': 70,\
    },
    'velocity':  {\
        'ppi': 30,\
        'rhi': 30,\
    },
    'differential_reflectivity': {\
        'ppi': 5.0, \
        'rhi': 5.0, \
    },
    'cross_correlation_ratio': {\
        'ppi': 1.,\
        'rhi': 1.,\
    }
}
cmap_dict = {
    'composite_reflectivity': 'gist_ncar',\
    'reflectivity': 'gist_ncar',\
    'velocity':    'bwr',\
    'differential_reflectivity': 'gist_ncar', \
    'cross_correlation_ratio': 'pyart_RefDiff',\
}



def init_proj(NEXRAD_dict):
    #mapcrs = Miller()
    if(NEXRAD_dict == None):
        mapcrs = ccrs.LambertConformal()
    else:
        mapcrs = ccrs.LambertConformal(central_longitude = NEXRAD_dict['center_lon'],\
            central_latitude = NEXRAD_dict['center_lat'])

    return mapcrs


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Composite Reflectivity function from PyART 
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def composite_reflectivity(radar, field="reflectivity", gatefilter=None):
    """
    Composite Reflectivity

    Often a legacy product, composite reflectivity is:
    "A display or mapping of the maximum radar reflectivity factor at any
    altitude as a function of position on the ground." - AMS Glossary
    This is more useful for the dry regions of the world, where maximum
    reflectivity values are found aloft, as opposed to the lowest scan.
    Alternatively this is useful for comparing to NWP since composite Z
    is a standard output of NWP.

    Why this is not as easy as one would think: Turns out the data are
    not natively stored with index 0 being azimuth 0. Likely due to the
    physical spinning of the radar antenna.

    Author: Randy J. Chase (@dopplerchase)

    Parameters
    ----------
    radar : Radar
        Radar object used.
    field : str
        Reflectivity field name to use to look up reflectivity data. In the
        radar object. Default field name is 'reflectivity'.
    gatefilter : GateFilter
        GateFilter instance. None will result in no gatefilter mask being
        applied to data.

    Returns
    -------
    radar : Radar
        The radar object containing the radar dimensions, metadata and
        composite field.

    References
    ----------
    American Meteorological Society, 2022: "Composite reflectivity". Glossary of Meteorology,
    http://glossary.ametsoc.org/wiki/Composite_reflectivity

    """

    # Determine the lowest sweep (used for metadata and such)
    minimum_sweep = np.min(radar.sweep_number["data"])

    # loop over all measured sweeps
    for sweep in sorted(radar.sweep_number["data"]):
        # get start and stop index numbers
        sweep_slice = radar.get_slice(sweep)

        # grab radar data
        z = radar.get_field(sweep, field)
        z_dtype = z.dtype

        # Use gatefilter
        if gatefilter is not None:
            mask_sweep = gatefilter.gate_excluded[sweep_slice, :]
            z = np.ma.masked_array(z, mask_sweep)

        # extract lat lons
        lon = radar.gate_longitude["data"][sweep_slice, :]
        lat = radar.gate_latitude["data"][sweep_slice, :]

        # get the range and time
        ranges = radar.range["data"]
        time = radar.time["data"]

        # get azimuth
        az = radar.azimuth["data"][sweep_slice]
        # get order of azimuths
        az_ids = np.argsort(az)

        # reorder azs so they are in order
        az = az[az_ids]
        z = z[az_ids]
        lon = lon[az_ids]
        lat = lat[az_ids]
        time = time[az_ids]

        # if the first sweep, store re-ordered lons/lats
        if sweep == minimum_sweep:
            azimuth_final = az
            time_final = time
            lon_0 = copy.deepcopy(lon)
            lon_0[-1, :] = lon_0[0, :]
            lat_0 = copy.deepcopy(lat)
            lat_0[-1, :] = lat_0[0, :]

        else:
            # Configure the intperpolator
            z_interpolator = interp2d(ranges, az, z, kind="linear")

            # Apply the interpolation
            z = z_interpolator(ranges, azimuth_final)

        # if first sweep, create new dim, otherwise concat them up
        if sweep == minimum_sweep:
            z_stack = copy.deepcopy(z[np.newaxis, :, :])
        else:
            z_stack = np.concatenate([z_stack, z[np.newaxis, :, :]])

    # now that the stack is made, take max across vertical
    compz = z_stack.max(axis=0).astype(z_dtype)

    # since we are using the whole volume scan, report mean time
    try:
        dtime = to_datetime(
            num2date(radar.time["data"], radar.time["units"]).astype(str),
            format="ISO8601",
        )
    except ValueError:
        dtime = to_datetime(
            num2date(radar.time["data"], radar.time["units"]).astype(str)
        )
    dtime = dtime.mean()

    # return dict, because this is was pyart does with lots of things
    fields = {}
    fields["composite_reflectivity"] = {
        "data": compz,
        "units": "dBZ",
        "long_name": "composite_reflectivity",
        "comment": "composite reflectivity computed from calculating the max radar value in each radar gate vertically after reordering",
    }

    time = radar.time.copy()
    time["data"] = time_final
    time["mean"] = dtime

    gate_longitude = radar.gate_longitude.copy()
    gate_longitude["data"] = lon_0
    gate_longitude["comment"] = "reordered longitude grid, [az,range]"
    gate_latitude = radar.gate_latitude.copy()
    gate_latitude["data"] = lat_0
    gate_latitude["comment"] = "reordered latitude grid, [az,range]"

    _range = radar.range.copy()
    metadata = radar.metadata.copy()
    scan_type = radar.scan_type
    latitude = radar.latitude.copy()
    longitude = radar.longitude.copy()
    altitude = radar.altitude.copy()
    instrument_parameters = radar.instrument_parameters

    sweep_number = radar.sweep_number.copy()
    sweep_number["data"] = np.array([0], dtype="int32")
    sweep_mode = radar.sweep_mode.copy()
    sweep_mode["data"] = np.array([radar.sweep_mode["data"][0]])
    ray_shape = compz.shape[0]
    azimuth = radar.azimuth.copy()
    azimuth["data"] = azimuth_final
    elevation = radar.elevation.copy()
    elevation["data"] = np.zeros(ray_shape, dtype="float32")
    fixed_angle = radar.fixed_angle.copy()
    fixed_angle["data"] = np.array([0.0], dtype="float32")
    sweep_start_ray_index = radar.sweep_start_ray_index.copy()
    sweep_start_ray_index["data"] = np.array([0], dtype="int32")
    sweep_end_ray_index = radar.sweep_end_ray_index.copy()
    sweep_end_ray_index["data"] = np.array([ray_shape - 1], dtype="int32")

    return Radar(
        time,
        _range,
        fields,
        metadata,
        scan_type,
        latitude,
        longitude,
        altitude,
        sweep_number,
        sweep_mode,
        fixed_angle,
        sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth,
        elevation,
        instrument_parameters=instrument_parameters,
    ), z_stack



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Download data
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# date_str: YYYYMMDDHHMM
# channels is a list containing the numbers of the desired channels
def download_NEXRAD_bucket(date_str, radar, \
        dest_dir = data_dir):

    # Convert the input date_str to datetime
    # --------------------------------------
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    #if((dt_date_str.minute == 0)):
    #    dt_date_str = dt_date_str - timedelta(minutes = 5)

    # Pull the entire file list for this date and time
    # ------------------------------------------------
    request_add = dt_date_str.strftime('s3://noaa-nexrad-level2/%Y/%m/%d/' + \
        radar + '/')

    # For each desired channel, figure out the closest file time
    # to the input date
    # ----------------------------------------------------------
    try:
        files = subprocess.run(['aws','s3','ls','--no-sign-request',\
            request_add], check=True, \
            capture_output=True).stdout.decode('utf-8').strip().split('\n') 
    except subprocess.CalledProcessError:
        print("ERROR: No ",radar," files for the input DTG",date_str)
        return

    # Remove the timestamps from the file strings
    # -------------------------------------------
    files_only = [tfile.strip().split()[3] for tfile in files]

    # Use the actual file times to get timestamps
    # -------------------------------------------
    file_dates = [datetime.strptime(tfile[4:19],'%Y%m%d_%H%M%S') for tfile in files_only] 

    # Figure out where the closest files to the desired time are
    # ----------------------------------------------------------
    time_diffs = np.array([abs((dt_date_str - ddate).total_seconds()) \
        for ddate in file_dates])

    # Select those files. Should result in all 16 channels for a single
    # time
    # -----------------------------------------------------------------
    files_found = np.array(files_only)[np.where(\
        time_diffs == np.min(time_diffs))]
    good_file = files_found[0]

    # Download the file 
    # -----------------
    request_add = dt_date_str.strftime('s3://noaa-nexrad-level2/%Y/%m/%d/'\
        + radar + '/' + good_file)
    #print(request_add)
    cmnd_list = ['aws','s3','cp','--no-sign-request',\
        request_add, dest_dir + radar + '/']
    print(' '.join(cmnd_list))
    subprocess.run(cmnd_list, check=True, \
        capture_output=False)

# begin_date : YYYYMMDDHHMM
# end_date   : YYYYMMDDHHMM
# interval   : minutes between desired images
# radar      : WSR-88D radar to download from
# dest_dir   : the location where files are going to be downloaded to.
#              NOTE: This MUST be changed when running on a computer
#              that does not belong to bsorenson.
def auto_NEXRAD_download(begin_date, end_date, interval, radar, \
        dest_dir = data_dir):

    # Convert the input date_str to datetime
    # --------------------------------------
    begin_dt_date = datetime.strptime(begin_date,"%Y%m%d%H%M")
    end_dt_date   = datetime.strptime(end_date,"%Y%m%d%H%M")

    # Using the interval, get the desired file times for each
    # NEXRAD image
    # -------------------------------------------------------
    num_seconds_in_minute = 60
    num_minutes = (end_dt_date - begin_dt_date).total_seconds() / \
        num_seconds_in_minute
    num_imgs = num_minutes / interval


    image_times = [begin_dt_date + timedelta(minutes = interval * ii) \
        for ii in range(int(num_imgs))]

    for ttime in image_times:
        print(ttime)
        download_NEXRAD_bucket(ttime.strftime('%Y%m%d%H%M'), \
            radar = radar, dest_dir = dest_dir)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Miscellaneous plotting functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

def find_NEXRAD_value_point(radar, field, xidx, yidx, angle):

    #sweep_slice = radar.get_slice(angle)
   
    #data = radar.fields[field]['data'][sweep_slice]

    sweep_slice = radar.get_slice(angle)

    # grab radar data
    z = radar.get_field(angle, field)

    print(z.shape)

    point_value = z[xidx, yidx]
    #point_value = data[xidx, yidx]

    return point_value
 
    ##!## grab radar data
    ##!#z = radar.get_field(sweep, field)
    ##!#z_dtype = z.dtype

    ##!### Use gatefilter
    ##!##if gatefilter is not None:
    ##!##    mask_sweep = gatefilter.gate_excluded[sweep_slice, :]
    ##!##    z = np.ma.masked_array(z, mask_sweep)

    ##!## extract lat lons
    ##!#lon = radar.gate_longitude["data"][sweep_slice, :]
    ##!#lat = radar.gate_latitude["data"][sweep_slice, :]

    ##!## get the range and time
    ##!#ranges = radar.range["data"]
    ##!#time = radar.time["data"]

    ##!## get azimuth
    ##!#az = radar.azimuth["data"][sweep_slice]
    ##!## get order of azimuths
    ##!#az_ids = np.argsort(az)

    ##!## reorder azs so they are in order
    ##!#az = az[az_ids]
    ##!#z = z[az_ids]
    ##!#lon = lon[az_ids]
    ##!#lat = lat[az_ids]
    ##!#time = time[az_ids]

    ##!#azimuth_final = az
    ##!#time_final = time
    ##!#lon_0 = copy.deepcopy(lon)
    ##!#lon_0[-1, :] = lon_0[0, :]
    ##!#lat_0 = copy.deepcopy(lat)
    ##!#lat_0[-1, :] = lat_0[0, :]

    ##!## Configure the intperpolator
    ##!#z_interpolator = interp2d(ranges, az, z, kind="linear")

    ##!## Apply the interpolation
    ##!#z = z_interpolator(ranges, azimuth_final)
   
    ##!#return z[xidx, yidx]



def find_NEXRAD_profile_point(radar, field, xidx, yidx):

    # Determine the lowest sweep (used for metadata and such)
    minimum_sweep = np.min(radar.sweep_number["data"])

    # loop over all measured sweeps
    for sweep in sorted(radar.sweep_number["data"]):
        # get start and stop index numbers
        sweep_slice = radar.get_slice(sweep)

        # grab radar data
        z = radar.get_field(sweep, field)
        z_dtype = z.dtype

        ## Use gatefilter
        #if gatefilter is not None:
        #    mask_sweep = gatefilter.gate_excluded[sweep_slice, :]
        #    z = np.ma.masked_array(z, mask_sweep)

        # extract lat lons
        lon = radar.gate_longitude["data"][sweep_slice, :]
        lat = radar.gate_latitude["data"][sweep_slice, :]

        # get the range and time
        ranges = radar.range["data"]
        time = radar.time["data"]

        # get azimuth
        az = radar.azimuth["data"][sweep_slice]
        # get order of azimuths
        az_ids = np.argsort(az)

        # reorder azs so they are in order
        az = az[az_ids]
        z = z[az_ids]
        lon = lon[az_ids]
        lat = lat[az_ids]
        time = time[az_ids]

        # if the first sweep, store re-ordered lons/lats
        if sweep == minimum_sweep:
            azimuth_final = az
            time_final = time
            lon_0 = copy.deepcopy(lon)
            lon_0[-1, :] = lon_0[0, :]
            lat_0 = copy.deepcopy(lat)
            lat_0[-1, :] = lat_0[0, :]

        else:
            # Configure the intperpolator
            z_interpolator = interp2d(ranges, az, z, kind="linear")

            # Apply the interpolation
            z = z_interpolator(ranges, azimuth_final)

        # if first sweep, create new dim, otherwise concat them up
        if sweep == minimum_sweep:
            z_stack = copy.deepcopy(z[np.newaxis, :, :])
        else:
            z_stack = np.concatenate([z_stack, z[np.newaxis, :, :]])

    profile_val = z_stack[:, xidx, yidx]

    # now that the stack is made, take max across vertical
    #compz = z_stack.max(axis=0).astype(z_dtype)

    return profile_val 

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Plotting functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Adapted from Alec Sczepanski's ASGSA PyArt Tools of the Trade
# -------------------------------------------------------------
def read_NEXRAD(date_str, radar, angle = 0):

    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')

    # Use glob to find the radar files for this radar
    # -----------------------------------------------
    radar_files = glob(data_dir + radar + '/*')

    # Use the actual file times to get timestamps
    # -------------------------------------------
    file_dates = [datetime.strptime(tfile.strip().split('/')[-1][4:19],\
        '%Y%m%d_%H%M%S') for tfile in radar_files] 

    # Figure out where the closest files to the desired time are
    # ----------------------------------------------------------
    time_diffs = np.array([abs((dt_date_str - ddate).total_seconds()) \
        for ddate in file_dates])

    # If none of the time differences are less than 15 minutes, abort
    # ---------------------------------------------------------------
    if(np.min(time_diffs) > (15 * 60)):
        print("ERROR: No ",radar," swaths matching time",date_str)
        return -1

    # Select those files. Should result in 2 files, with one being MDM
    # -----------------------------------------------------------------
    files_found = np.array(radar_files)[np.where(\
        time_diffs == np.min(time_diffs))]
    good_file = files_found[0]
    print(good_file)

    # - - - - - - - - - - - - - -
    # Begin Alec's section
    # - - - - - - - - - - - - - -
    # Open the file:
    data  = pyart.io.read(good_file)
    
    # Establish radar center latitude and longitude:
    center_lon = data.longitude['data'][0]
    center_lat = data.latitude['data'][0]
    
    # Establish min and max latitude and longitude to be plotted:
    min_lat = center_lat - 2 #subtracting 2 degrees lat
    max_lat = center_lat + 2 #adding 2 degrees lat
    min_lon = center_lon - 2 #subtracting 2 degrees lon
    max_lon = center_lon + 2 #adding 2 degrees lon

    # Pull desired elevation angle from data:
    # (Can change the [0] in the next line to another number for another
    # elevation angle. [0] is the lowest elevation angle, typically 0.5 deg)
    angles = data.fixed_angle['data'][:]
    fixed_angle = angles[angle]
    
    # Pull time data:
    radar_date = pyart.graph.common.generate_radar_time_begin(data)
    radar_date_str =radar_date.strftime('%Y%m%d%H%M')
    
    # Generate figure title:
    figure_title = radar + ' - Elevation = %.3f' % (fixed_angle) + \
        ' - Time: ' + radar_date.strftime('%Y%m%d %H:%M')

    # - - - - - - - - - - - - - -
    # End Alec's section
    # - - - - - - - - - - - - - -

    NEXRAD_dict = {}
    NEXRAD_dict['radar']         = data 
    NEXRAD_dict['radar_name']    = radar
    NEXRAD_dict['center_lat']    = center_lat
    NEXRAD_dict['center_lon']    = center_lon
    NEXRAD_dict['min_lat']       = min_lat
    NEXRAD_dict['min_lon']       = min_lon
    NEXRAD_dict['max_lat']       = max_lat
    NEXRAD_dict['max_lon']       = max_lon
    NEXRAD_dict['angles']        = angles
    NEXRAD_dict['angle_idx']     = angle
    NEXRAD_dict['fixed_angle']   = fixed_angle
    NEXRAD_dict['radar_date']    = radar_date_str
    NEXRAD_dict['figure_title']  = figure_title
    
    return NEXRAD_dict

title_dict = {
    'composite_reflectivity': 'Composite Equiv. Refl. Fact. (dBZ)',\
    'reflectivity': 'Equiv. Refl. Fact. (dBZ)',\
    'velocity':     'Radial Velocity (m/s)',\
    'differential_reflectivity': 'Log Diff. Refl. Factor (dB)', \
    'cross_correlation_ratio': 'Cross Corr. Ratio',\
}

#def plot_NEXRAD_GOES(NEXRAD_dict, variable, channel, ax = None, \
def plot_NEXRAD_GOES_2panel(date_str, radar, variable, channel, ax = None, \
        angle = 0, ptitle = None, plabel = None, vmin = -5, vmax = 90, \
        labelsize = 10, colorbar = True, counties = True, save_dir = './',\
        alpha = 1.0, mask_outside = True, zoom=True, save=False):

    if('/home/bsorenson/Research/GOES' not in sys.path):
        sys.path.append('/home/bsorenson/Research/GOES')
    from GOESLib import read_GOES_satpy, plot_GOES_satpy
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    #dt_date_str = datetime.strptime(NEXRAD_dict['radar_date'],"%Y%m%d%H%M")

    # Read the NEXRAD data
    # --------------------
    NEXRAD_dict = read_NEXRAD(date_str, radar, angle = angle)

    # Read the GOES data
    # ------------------
    var0, crs0, lons0, lats0, lat_lims, lon_lims, plabel_goes = \
        read_GOES_satpy(date_str, channel)

    labelsize = 10
    # Set up the figure
    # -----------------
    plt.close('all')
    fig = plt.figure(figsize = (9, 4))
    ax1 = fig.add_subplot(1,2,1, projection = crs0)
    ax2 = fig.add_subplot(1,2,2, projection = crs0)
    plot_GOES_satpy(NEXRAD_dict['radar_date'], channel, ax = ax1, var = var0, \
        crs = crs0, lons = lons0, lats = lats0, lat_lims = lat_lims, \
        lon_lims = lon_lims, vmin = 5, vmax = 40, ptitle = '', \
        plabel = plabel_goes, colorbar = True, labelsize = labelsize + 1, \
        counties = counties, zoom=True,save=False) 
    plot_NEXRAD_ppi(NEXRAD_dict, variable, ax = ax2, angle = angle, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = None)
        #mask_outside = mask_outside, crs = crs0)

    fig.tight_layout()

    plt.show()

# Plots GOES data, KBBX and KRGX radar data
def plot_NEXRAD_GOES_3panel(date_str, radar1, radar2, variable, channel, ax = None, \
        angle = 0, ptitle = None, plabel = None, vmin = -5, vmax = 90, \
        labelsize = 10, colorbar = True, counties = True, save_dir = './',\
        alpha = 1.0, mask_outside = True, zoom=True, save=False):

    if('/home/bsorenson/Research/GOES' not in sys.path):
        sys.path.append('/home/bsorenson/Research/GOES')
    from GOESLib import read_GOES_satpy, plot_GOES_satpy
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    #dt_date_str = datetime.strptime(NEXRAD_dict['radar_date'],"%Y%m%d%H%M")

    # Read the NEXRAD data
    # --------------------
    NEXRAD_dict1 = read_NEXRAD(date_str, radar1, angle = angle)
    NEXRAD_dict2 = read_NEXRAD(date_str, radar2, angle = angle)

    # Read the GOES data
    # ------------------
    var0, crs0, lons0, lats0, lat_lims, lon_lims, plabel_goes = \
        read_GOES_satpy(date_str, channel)

    labelsize = 10
    # Set up the figure
    # -----------------
    plt.close('all')
    mapcrs = init_proj(NEXRAD_dict1)
    fig = plt.figure(figsize = (11, 4))
    ax1 = fig.add_subplot(1,3,1, projection = crs0)
    ax2 = fig.add_subplot(1,3,2, projection = mapcrs)
    ax3 = fig.add_subplot(1,3,3, projection = mapcrs)
    plot_GOES_satpy(NEXRAD_dict1['radar_date'], channel, ax = ax1, var = var0, \
        crs = crs0, lons = lons0, lats = lats0, lat_lims = lat_lims, \
        lon_lims = lon_lims, vmin = 5, vmax = 40, ptitle = '', \
        plabel = plabel_goes, colorbar = True, labelsize = labelsize + 1, \
        counties = counties, zoom=True,save=False) 
    plot_NEXRAD_ppi(NEXRAD_dict1, variable, ax = ax2, angle = angle, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = None)
    plot_NEXRAD_ppi(NEXRAD_dict2, variable, ax = ax3, angle = angle, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = None)

    plt.suptitle(date_str)

    fig.tight_layout()

    del(NEXRAD_dict1)
    del(NEXRAD_dict2)
    del(var0)

    if(save):
        outname = 'nexrad_goes_3panel_' + radar1 + '_' + radar2 + '_' + date_str + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()



#def plot_NEXRAD_GOES(NEXRAD_dict, variable, channel, ax = None, \
def plot_NEXRAD_GOES_4panel(date_str, radar, variable, channel, ax = None, \
        angle1 = 0, angle2 = 2, angle3 = 4, ptitle = None, plabel = None, \
        vmin = -5, vmax = 90, \
        labelsize = 10, colorbar = True, counties = True, save_dir = './',\
        alpha = 1.0, mask_outside = True, zoom=True, save=False):

    if('/home/bsorenson/Research/GOES' not in sys.path):
        sys.path.append('/home/bsorenson/Research/GOES')
    from GOESLib import read_GOES_satpy, plot_GOES_satpy
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    #dt_date_str = datetime.strptime(NEXRAD_dict['radar_date'],"%Y%m%d%H%M")

    # Read the NEXRAD data
    # --------------------
    NEXRAD_dict1 = read_NEXRAD(date_str, radar, angle = angle1)
    NEXRAD_dict2 = read_NEXRAD(date_str, radar, angle = angle2)
    NEXRAD_dict3 = read_NEXRAD(date_str, radar, angle = angle3)

    # Read the GOES data
    # ------------------
    var0, crs0, lons0, lats0, lat_lims, lon_lims, plabel_goes = \
        read_GOES_satpy(date_str, channel)

    labelsize = 10
    # Set up the figure
    # -----------------
    plt.close('all')
    fig = plt.figure(figsize = (7, 7))
    ax1 = fig.add_subplot(2,2,1, projection = crs0)
    ax2 = fig.add_subplot(2,2,2, projection = crs0)
    ax3 = fig.add_subplot(2,2,3, projection = crs0)
    ax4 = fig.add_subplot(2,2,4, projection = crs0)
    plot_GOES_satpy(NEXRAD_dict1['radar_date'], channel, ax = ax1, var = var0, \
        crs = crs0, lons = lons0, lats = lats0, lat_lims = lat_lims, \
        lon_lims = lon_lims, vmin = 2.5, vmax = 40, ptitle = '', \
        plabel = plabel_goes, colorbar = True, labelsize = labelsize + 1, \
        counties = counties, zoom=True,save=False) 
    plot_NEXRAD_ppi(NEXRAD_dict1, variable, ax = ax2, angle = angle1, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = None)
    plot_NEXRAD_ppi(NEXRAD_dict2, variable, ax = ax3, angle = angle2, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = None)
    plot_NEXRAD_ppi(NEXRAD_dict3, variable, ax = ax4, angle = angle3, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = None)
        #mask_outside = mask_outside, crs = crs0)

    fig.tight_layout()

    if(save):
        outname = 'nexrad_goes_4panel_' + radar + '_' + date_str + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

#def plot_NEXRAD_GOES(NEXRAD_dict, variable, channel, ax = None, \
def plot_NEXRAD_GOES_5panel(date_str, ch1, ch2, ch3, \
        variable = 'composite_reflectivity', \
        ax = None, ptitle = None, plabel = None, \
        vmin = -5, vmax = 90, angle1 = 0, angle2 = 0, \
        labelsize = 10, colorbar = True, counties = True, save_dir = './',\
        alpha = 1.0, mask_outside = True, zoom=True, save=False):

    if('/home/bsorenson/Research/GOES' not in sys.path):
        sys.path.append('/home/bsorenson/Research/GOES')
    from GOESLib import read_GOES_satpy, plot_GOES_satpy, goes_channel_dict
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    #dt_date_str = datetime.strptime(NEXRAD_dict['radar_date'],"%Y%m%d%H%M")

    # Read the NEXRAD data
    # --------------------
   
    #try: 
    NEXRAD_dict1 = read_NEXRAD(date_str, 'KBBX', angle = angle1)
    if(NEXRAD_dict1 == -1): 
        plot_radar_1 = False
        print("WARNING: Error reading KBBX data for ", date_str)
        print("         Continuing without the data") 
    else:
        plot_radar_1 = True

    NEXRAD_dict2 = read_NEXRAD(date_str, 'KRGX', angle = angle2) 
    if(NEXRAD_dict2 == -1): 
        print("WARNING: Error reading KRGX data for ", date_str)
        print("         Continuing without the data") 
        plot_radar_2 = False
    else: 
        plot_radar_2 = True 

    # Read the GOES data
    # ------------------
    var0, crs0, lons0, lats0, lat_lims, lon_lims, plabel_goes1 = \
        read_GOES_satpy(date_str, ch1)
    var1, crs1, lons1, lats1, lat_lims, lon_lims, plabel_goes2 = \
        read_GOES_satpy(date_str, ch2)
    var2, crs2, lons2, lats2, lat_lims, lon_lims, plabel_goes3 = \
        read_GOES_satpy(date_str, ch3)

    labelsize = 10
    # Set up the figure
    # -----------------
    plt.close('all')
    fig = plt.figure(figsize = (11, 7))
    ax1 = fig.add_subplot(2,3,1, projection = crs0)
    ax2 = fig.add_subplot(2,3,2, projection = crs0)
    ax3 = fig.add_subplot(2,3,3, projection = crs0)
    ax4 = fig.add_subplot(2,3,5, projection = crs0)
    ax5 = fig.add_subplot(2,3,6, projection = crs0)
    plot_GOES_satpy(NEXRAD_dict1['radar_date'], ch1, ax = ax1, var = var0, \
        crs = crs0, lons = lons0, lats = lats0, lat_lims = lat_lims, \
        lon_lims = lon_lims, vmin = 0, vmax = 60, ptitle = '', \
        plabel = plabel_goes1, colorbar = True, labelsize = labelsize + 1, \
        counties = counties, zoom=True,save=False) 
    plot_GOES_satpy(NEXRAD_dict1['radar_date'], ch2, ax = ax2, var = var1, \
        crs = crs1, lons = lons1, lats = lats1, lat_lims = lat_lims, \
        lon_lims = lon_lims, vmin = 2.5, vmax = 40, ptitle = '', \
        plabel = plabel_goes2, colorbar = True, labelsize = labelsize + 1, \
        counties = counties, zoom=True,save=False) 
    plot_GOES_satpy(NEXRAD_dict1['radar_date'], ch3, ax = ax3, var = var2, \
        crs = crs2, lons = lons2, lats = lats2, lat_lims = lat_lims, \
        lon_lims = lon_lims, vmin = 270, vmax = 320, ptitle = '', \
        plabel = plabel_goes3, colorbar = True, labelsize = labelsize + 1, \
        counties = counties, zoom=True,save=False) 
    if(plot_radar_1): 
        plot_NEXRAD_ppi(NEXRAD_dict1, variable, ax = ax4, angle = angle1, \
            counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
            mask_outside = mask_outside, crs = None)
        del(NEXRAD_dict1)
    if(plot_radar_2): 
        plot_NEXRAD_ppi(NEXRAD_dict2, variable, ax = ax5, angle = angle2, \
            counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
            mask_outside = mask_outside, crs = None)
            #mask_outside = mask_outside, crs = crs0)
        del(NEXRAD_dict2)

    font_size = 10
    plot_figure_text(ax1, \
        str(goes_channel_dict[str(ch1)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    plot_figure_text(ax2, \
        str(goes_channel_dict[str(ch2)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    plot_figure_text(ax3, \
        str(goes_channel_dict[str(ch3)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')

    plt.suptitle(date_str)
   
    plot_subplot_label(ax1, '(a)', backgroundcolor = 'white', \
        fontsize = font_size, location = 'upper_right')
    plot_subplot_label(ax2, '(b)', backgroundcolor = 'white', \
        fontsize = font_size)

    fig.tight_layout()

    del(var0)
    del(var1)
    del(var2)

    if(save):
        if(variable == 'composite_reflectivity'):
            outname = 'nexrad_goes_5panel_KBBX_KRGX_ang1cr_ang2cr_' + \
                date_str + '.png'
        else:   
            outname = 'nexrad_goes_5panel_KBBX_KRGX_ang1' + str(angle1).zfill(2) + \
                '_ang2' + str(angle2).zfill(2) + '_' + date_str + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

#def plot_NEXRAD_GOES(NEXRAD_dict, variable, channel, ax = None, \
def plot_NEXRAD_GOES_12panel(date_str, ch1, ch2, ch3, \
        azm1, azm2, azm3, 
        variable = 'composite_reflectivity', \
        az_tol = 0.8, 
        ax = None, ptitle = None, plabel = None, \
        vmin = -5, vmax = 90, \
        sat = 'goes17',
        labelsize = 10, colorbar = True, counties = True, save_dir = './',\
        alpha = 1.0, mask_outside = True, zoom=True, save=False):

    if('/home/bsorenson/Research/GOES' not in sys.path):
        sys.path.append('/home/bsorenson/Research/GOES')
    from GOESLib import read_GOES_satpy, plot_GOES_satpy, goes_channel_dict
    if('/home/bsorenson/Research/MODIS' not in sys.path):
        sys.path.append('/home/bsorenson/Research/MODIS/obs_smoke_forcing/')
    from MODISLib import read_MODIS_satpy
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    #dt_date_str = datetime.strptime(NEXRAD_dict['radar_date'],"%Y%m%d%H%M")

    # Read the NEXRAD data
    # --------------------
    NEXRAD_dict1 = read_NEXRAD(date_str, 'KBBX', angle = 0)
    NEXRAD_dict2 = read_NEXRAD(date_str, 'KRGX', angle = 0) 

    # Read the GOES data
    # ------------------
    if(sat == 'goes16'):
        var0, crs0, lons0, lats0, lat_lims, lon_lims, plabel_goes1, xx0, yy0 = \
            read_GOES_satpy(date_str, ch1, sat = sat, return_xy = True, zoom = False)
        var1, crs1, lons1, lats1, lat_lims, lon_lims, plabel_goes2, xx1, yy1 = \
            read_GOES_satpy(date_str, ch2, sat = sat, return_xy = True, zoom = False)
        var2, crs2, lons2, lats2, lat_lims, lon_lims, plabel_goes3, xx2, yy2 = \
            read_GOES_satpy(date_str, ch3, sat = sat, return_xy = True, zoom = False)
    else:
        var0, crs0, lons0, lats0, lat_lims, lon_lims, plabel_goes1 = \
            read_GOES_satpy(date_str, ch1, sat = sat)
        var1, crs1, lons1, lats1, lat_lims, lon_lims, plabel_goes2 = \
            read_GOES_satpy(date_str, ch2, sat = sat)
        var2, crs2, lons2, lats2, lat_lims, lon_lims, plabel_goes3 = \
            read_GOES_satpy(date_str, ch3, sat = sat)

    labelsize = 10
    # Set up the figure
    # -----------------
    plt.close('all')
    fig  = plt.figure(figsize = (11, 13))
    if(sat == 'goes16'):
        mapcrs = init_proj(NEXRAD_dict1)
        ax1  = fig.add_subplot(4,3,1, projection = mapcrs) # GOES VIS
        ax2  = fig.add_subplot(4,3,2, projection = mapcrs) # GOES SWIR
        ax3  = fig.add_subplot(4,3,3, projection = mapcrs) # GOES TIR
        ax4  = fig.add_subplot(4,3,4, projection = mapcrs) # KBBX REFL PPI AZM1
        ax5  = fig.add_subplot(4,3,5)                    # KBBX REFL RHI AZM1
        ax6  = fig.add_subplot(4,3,6)                    # KBBX  CC  RHI AZM1
        ax7  = fig.add_subplot(4,3,7, projection = mapcrs) # KRGX REFL PPI AZM1
        ax8  = fig.add_subplot(4,3,8)                    # KRGX REFL RHI AZM1
        ax9  = fig.add_subplot(4,3,9)                    # KRGX  CC  RHI AZM1
        ax10 = fig.add_subplot(4,3,10, projection = mapcrs) # KRGX REFL PPI AZM2
        ax11 = fig.add_subplot(4,3,11)                    # KRGX REFL RHI AZM2
        ax12 = fig.add_subplot(4,3,12)                    # KRGX  CC  RHI AZM
    else:
        ax1  = fig.add_subplot(4,3,1, projection = crs0) # GOES VIS
        ax2  = fig.add_subplot(4,3,2, projection = crs0) # GOES SWIR
        ax3  = fig.add_subplot(4,3,3, projection = crs0) # GOES TIR
        ax4  = fig.add_subplot(4,3,4, projection = crs0) # KBBX REFL PPI AZM1
        ax5  = fig.add_subplot(4,3,5)                    # KBBX REFL RHI AZM1
        ax6  = fig.add_subplot(4,3,6)                    # KBBX  CC  RHI AZM1
        ax7  = fig.add_subplot(4,3,7, projection = crs0) # KRGX REFL PPI AZM1
        ax8  = fig.add_subplot(4,3,8)                    # KRGX REFL RHI AZM1
        ax9  = fig.add_subplot(4,3,9)                    # KRGX  CC  RHI AZM1
        ax10 = fig.add_subplot(4,3,10, projection = crs0) # KRGX REFL PPI AZM2
        ax11 = fig.add_subplot(4,3,11)                    # KRGX REFL RHI AZM2
        ax12 = fig.add_subplot(4,3,12)                    # KRGX  CC  RHI AZM

    if(sat == 'goes16'):
        plot_GOES_satpy(NEXRAD_dict1['radar_date'], ch1, ax = ax1, var = var0, \
            crs = crs0, lons = lons0, lats = lats0, lat_lims = lat_lims, \
            xx = xx0, yy = yy0, 
            lon_lims = lon_lims, vmin = 0, vmax = 60, ptitle = '', \
            plabel = plabel_goes1, colorbar = True, labelsize = labelsize + 1, \
            counties = counties, zoom=True,save=False, use_xy = True) 
        plot_GOES_satpy(NEXRAD_dict1['radar_date'], ch2, ax = ax2, var = var1, \
            crs = crs1, lons = lons1, lats = lats1, lat_lims = lat_lims, \
            xx = xx1, yy = yy1, 
            lon_lims = lon_lims, vmin = 2.5, vmax = 40, ptitle = '', \
            plabel = plabel_goes2, colorbar = True, labelsize = labelsize + 1, \
            counties = counties, zoom=True,save=False, use_xy = True) 
        plot_GOES_satpy(NEXRAD_dict1['radar_date'], ch3, ax = ax3, var = var2, \
            crs = crs2, lons = lons2, lats = lats2, lat_lims = lat_lims, \
            xx = xx2, yy = yy2, 
            lon_lims = lon_lims, vmin = 270, vmax = 320, ptitle = '', \
            plabel = plabel_goes3, colorbar = True, labelsize = labelsize + 1, \
            counties = counties, zoom=True,save=False, use_xy = True) 
    else:
        plot_GOES_satpy(NEXRAD_dict1['radar_date'], ch1, ax = ax1, var = var0, \
            crs = crs0, lons = lons0, lats = lats0, lat_lims = lat_lims, \
            lon_lims = lon_lims, vmin = 0, vmax = 60, ptitle = '', \
            plabel = plabel_goes1, colorbar = True, labelsize = labelsize + 1, \
            counties = counties, zoom=True,save=False) 
        plot_GOES_satpy(NEXRAD_dict1['radar_date'], ch2, ax = ax2, var = var1, \
            crs = crs1, lons = lons1, lats = lats1, lat_lims = lat_lims, \
            lon_lims = lon_lims, vmin = 2.5, vmax = 40, ptitle = '', \
            plabel = plabel_goes2, colorbar = True, labelsize = labelsize + 1, \
            counties = counties, zoom=True,save=False) 
        plot_GOES_satpy(NEXRAD_dict1['radar_date'], ch3, ax = ax3, var = var2, \
            crs = crs2, lons = lons2, lats = lats2, lat_lims = lat_lims, \
            lon_lims = lon_lims, vmin = 270, vmax = 320, ptitle = '', \
            plabel = plabel_goes3, colorbar = True, labelsize = labelsize + 1, \
            counties = counties, zoom=True,save=False) 

    rmin2 = 110
    rmax2 = 180
    rmin3 = 110
    rmax3 = 180
    plot_azimuth_line(NEXRAD_dict1, azm1, 45, 120, ax1)
    plot_azimuth_line(NEXRAD_dict2, azm2, rmin2, rmax2, ax1)
    plot_azimuth_line(NEXRAD_dict2, azm3, rmin3, rmax3, ax1)
    plot_azimuth_line(NEXRAD_dict1, azm1, 45, 120, ax2)
    plot_azimuth_line(NEXRAD_dict2, azm2, rmin2, rmax2, ax2)
    plot_azimuth_line(NEXRAD_dict2, azm3, rmin3, rmax3, ax2)
    plot_azimuth_line(NEXRAD_dict1, azm1, 45, 120, ax3)
    plot_azimuth_line(NEXRAD_dict2, azm2, rmin2, rmax2, ax3)
    plot_azimuth_line(NEXRAD_dict2, azm3, rmin3, rmax3, ax3)

    plot_NEXRAD_ppi(NEXRAD_dict1, variable, ax = ax4, angle = 0, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = None)
    plot_azimuth_line(NEXRAD_dict1, azm1, 45, 120, ax4)
    plot_my_rhi(NEXRAD_dict1, 'reflectivity', azm1, \
            az_tol = az_tol, vmin = None, vmax = None, \
            range_min = 45, range_max = 120,\
            ax = ax5)
    plot_my_rhi(NEXRAD_dict1, 'cross_correlation_ratio', azm1, \
            az_tol = az_tol, vmin = None, vmax = None, \
            range_min = 45, range_max = 120,\
            ax = ax6)


    plot_NEXRAD_ppi(NEXRAD_dict2, variable, ax = ax7, angle = 0, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = None)
        #mask_outside = mask_outside, crs = crs0)
    plot_azimuth_line(NEXRAD_dict2, azm2, rmin2, rmax2, ax7)
    plot_my_rhi(NEXRAD_dict2, 'reflectivity', azm2, \
            az_tol = az_tol, vmin = None, vmax = None, \
            range_min = rmin2, range_max = rmax2,\
            ax = ax8)
    plot_my_rhi(NEXRAD_dict2, 'cross_correlation_ratio', azm2, \
            az_tol = az_tol, vmin = None, vmax = None, \
            range_min = rmin2, range_max = rmax2,\
            ax = ax9)


    plot_NEXRAD_ppi(NEXRAD_dict2, variable, ax = ax10, angle = 0, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = None)
    plot_azimuth_line(NEXRAD_dict2, azm3, rmin3, rmax3, ax10)
    plot_my_rhi(NEXRAD_dict2, 'reflectivity', azm3, \
            az_tol = az_tol, vmin = None, vmax = None, \
            range_min = rmin3, range_max = rmax3,\
            ax = ax11)
    plot_my_rhi(NEXRAD_dict2, 'cross_correlation_ratio', azm3, \
            az_tol = az_tol, vmin = None, vmax = None, \
            range_min = rmin3, range_max = rmax3,\
            ax = ax12)


    # Find the distances between the smoky point and KBBX
    # ---------------------------------------------------   
    if(sat == 'goes16'):
        radar_lat = NEXRAD_dict1['radar'].latitude['data'][0]
        radar_lon = NEXRAD_dict1['radar'].longitude['data'][0]
        dist1 = find_distance_between_points(40.2824, -121.2412, radar_lat, radar_lon)
        radar_lat = NEXRAD_dict2['radar'].latitude['data'][0]
        radar_lon = NEXRAD_dict2['radar'].longitude['data'][0]
        dist2 = find_distance_between_points(40.2824, -121.2412, radar_lat, radar_lon)
        msize = 8
        plot_point_on_map(ax1, 40.2824, -121.2412, markersize = msize, color = 'tab:blue')
        plot_point_on_map(ax1, 41.4914, -120.5644, markersize = msize, color = 'tab:orange')
        plot_point_on_map(ax2, 40.2824, -121.2412, markersize = msize, color = 'tab:blue')
        plot_point_on_map(ax2, 41.4914, -120.5644, markersize = msize, color = 'tab:orange')
        plot_point_on_map(ax3, 40.2824, -121.2412, markersize = msize, color = 'tab:blue')
        plot_point_on_map(ax3, 41.4914, -120.5644, markersize = msize, color = 'tab:orange')
        plot_point_on_map(ax4, 40.2824, -121.2412, markersize = msize, color = 'tab:blue')
        plot_point_on_map(ax4, 41.4914, -120.5644, markersize = msize, color = 'tab:orange')
        plot_point_on_map(ax7, 40.2824, -121.2412, markersize = msize, color = 'tab:blue')
        plot_point_on_map(ax7, 41.4914, -120.5644, markersize = msize, color = 'tab:orange')
        plot_point_on_map(ax10, 40.2824, -121.2412, markersize = msize, color = 'tab:blue')
        plot_point_on_map(ax10, 41.4914, -120.5644, markersize = msize, color = 'tab:orange')
        asos_height = 1379.67 / 1000.
        radar_height = NEXRAD_dict1['radar'].altitude['data'][0] / 1000.
        delta_height = asos_height
        ax5.plot([dist1],[delta_height], marker = 'o', color = 'tab:red')
        ax6.plot([dist1],[delta_height], marker = 'o', color = 'tab:red')
        radar_height = NEXRAD_dict2['radar'].altitude['data'][0] / 1000.
        delta_height = asos_height
        ax8.plot([dist2],[delta_height], marker = 'o', color = 'tab:red')
        ax9.plot([dist2],[delta_height], marker = 'o', color = 'tab:red')

    font_size = 10
    plot_figure_text(ax1, \
        str(goes_channel_dict[str(ch1)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right', weight = 'bold')
    plot_figure_text(ax2, \
        str(goes_channel_dict[str(ch2)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right', weight = 'bold')
    plot_figure_text(ax3, \
        str(goes_channel_dict[str(ch3)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right', weight = 'bold')

    plot_subplot_label(ax1, '(a)', backgroundcolor = 'white', \
        fontsize = font_size, location = 'upper_left')
    plot_subplot_label(ax2, '(b)', backgroundcolor = 'white', \
        fontsize = font_size, location = 'upper_left')
    plot_subplot_label(ax3, '(c)', backgroundcolor = 'white', \
        fontsize = font_size, location = 'upper_left')
    plot_subplot_label(ax4, '(d)', backgroundcolor = 'white', \
        fontsize = font_size, location = 'upper_left')
    plot_subplot_label(ax5, '(e)', backgroundcolor = 'white', \
        fontsize = font_size, location = 'upper_left')
    plot_subplot_label(ax6, '(f)', backgroundcolor = 'white', \
        fontsize = font_size, location = 'upper_left')
    plot_subplot_label(ax7, '(g)', backgroundcolor = 'white', \
        fontsize = font_size, location = 'upper_left')
    plot_subplot_label(ax8, '(h)', backgroundcolor = 'white', \
        fontsize = font_size, location = 'upper_left')
    plot_subplot_label(ax9, '(i)', backgroundcolor = 'white', \
        fontsize = font_size, location = 'upper_left')
    plot_subplot_label(ax10, '(j)', backgroundcolor = 'white', \
        fontsize = font_size, location = 'upper_left')
    plot_subplot_label(ax11, '(k)', backgroundcolor = 'white', \
        fontsize = font_size, location = 'upper_left')
    plot_subplot_label(ax12, '(l)', backgroundcolor = 'white', \
        fontsize = font_size, location = 'upper_left')

    plt.suptitle(dt_date_str.strftime('%d %B %Y %H:%M UTC'))
    
    fig.tight_layout()

    del(NEXRAD_dict1)
    del(NEXRAD_dict2)
    del(var0)
    del(var1)
    del(var2)

    if(save):
        outname = 'nexrad_goes_12panel_KBBX_KRGX_' + date_str + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

def plot_azimuth_line(NEXRAD_dict, azimuth, range_min, range_max, ax, \
        base_linewidth = 2.5, markersize = 2.0, marker = 'o', color = 'red'):

    xsect = pyart.util.cross_section_ppi(NEXRAD_dict['radar'], [azimuth])
    test_lats = xsect.gate_latitude['data'][0][::30]
    test_lons = xsect.gate_longitude['data'][0][::30]
    radar_lat = xsect.latitude['data'][0]
    radar_lon = xsect.longitude['data'][0]
    
    dists = np.array([find_distance_between_points(radar_lat, radar_lon, xx, yy) for xx, yy in zip(test_lats, test_lons)])
    keep_idxs = np.where( (dists >= range_min) & (dists <= range_max) )[0]
    lat1 = test_lats[keep_idxs[0]]
    lon1 = test_lons[keep_idxs[0]]
    lat2 = test_lats[keep_idxs[-1]]
    lon2 = test_lons[keep_idxs[-1]]
    
    ax.plot([lon1, lon2], [lat1, lat2], transform = datacrs, \
        linewidth = base_linewidth, markersize = markersize, \
        marker = marker, color = 'k')
    ax.plot([lon1, lon2], [lat1, lat2], transform = datacrs, \
        linewidth = base_linewidth - 1.5, markersize = markersize - 0.5, \
        marker = marker, color = color)

def plot_NEXRAD_rhi_multipanel(date_str, radar, azimuth, \
        angle_idx = 0, range_min = 0, range_max = 150, \
        height_lim = 8, save = True, save_dir = './'):

    NEXRAD_dict = read_NEXRAD(date_str, radar)
   
    plt.close('all') 
    mapcrs = init_proj(NEXRAD_dict)
    fig = plt.figure(figsize = (9, 10))
    ax1 = fig.add_subplot(3,3,1, projection = mapcrs)
    ax2 = fig.add_subplot(3,3,2, projection = mapcrs)
    ax3 = fig.add_subplot(3,3,3, projection = mapcrs)
    ax4 = fig.add_subplot(3,3,4)
    ax5 = fig.add_subplot(3,3,5)
    ax6 = fig.add_subplot(3,3,6)
    ax7 = fig.add_subplot(3,3,7)
    ax8 = fig.add_subplot(3,3,8)
    ax9 = fig.add_subplot(3,3,9)
    
    plot_NEXRAD_ppi(NEXRAD_dict, 'reflectivity', zoom = True, \
        angle = angle_idx, ax = ax1)
    plot_NEXRAD_ppi(NEXRAD_dict, 'cross_correlation_ratio', zoom = True, \
        angle = angle_idx, ax = ax2)
    plot_NEXRAD_ppi(NEXRAD_dict, 'differential_reflectivity', zoom = True, \
        angle = angle_idx, ax = ax3)
    
    az_tol = 0.8
    plot_my_rhi(NEXRAD_dict, 'reflectivity', azimuth, \
            az_tol = az_tol, vmin = None, vmax = None, \
            range_min = range_min, range_max = range_max,\
            ax = ax4, height_lim = height_lim)
    plot_my_rhi(NEXRAD_dict, 'cross_correlation_ratio', azimuth, \
            az_tol = az_tol, vmin = None, vmax = None, \
            range_min = range_min, range_max = range_max,\
            ax = ax5, height_lim = height_lim)
    plot_my_rhi(NEXRAD_dict, 'differential_reflectivity', azimuth, \
            az_tol = az_tol, vmin = None, vmax = None, \
            range_min = range_min, range_max = range_max,\
            ax = ax6, height_lim = height_lim)
    
    plot_NEXRAD_rhi(NEXRAD_dict, 'reflectivity', azimuth, \
        vmin = None, vmax = None, range_min = range_min, \
        range_max = range_max, ax = ax7, height_lim = height_lim)
    plot_NEXRAD_rhi(NEXRAD_dict, 'cross_correlation_ratio', azimuth, \
        vmin = None, vmax = None, range_min = range_min, \
        range_max = range_max, ax = ax8, height_lim = height_lim)
    plot_NEXRAD_rhi(NEXRAD_dict, 'differential_reflectivity', azimuth, \
        vmin = None, vmax = None, range_min = range_min, \
        range_max = range_max, ax = ax9, height_lim = height_lim)
    
    ##!#xsect = pyart.util.cross_section_ppi(NEXRAD_dict['radar'], [azimuth])
    ##!#test_lats = xsect.gate_latitude['data'][0][::30]
    ##!#test_lons = xsect.gate_longitude['data'][0][::30]
    ##!#radar_lat = xsect.latitude['data'][0]
    ##!#radar_lon = xsect.longitude['data'][0]
    ##!#
    ##!#dists = np.array([find_distance_between_points(radar_lat, radar_lon, xx, yy) for xx, yy in zip(test_lats, test_lons)])
    ##!#keep_idxs = np.where( (dists >= range_min) & (dists <= range_max) )[0]
    ##!#lat1 = test_lats[keep_idxs[0]]
    ##!#lon1 = test_lons[keep_idxs[0]]
    ##!#lat2 = test_lats[keep_idxs[-1]]
    ##!#lon2 = test_lons[keep_idxs[-1]]
    
    plot_azimuth_line(NEXRAD_dict, azimuth, range_min, range_max, ax1)
    plot_azimuth_line(NEXRAD_dict, azimuth, range_min, range_max, ax2)
    plot_azimuth_line(NEXRAD_dict, azimuth, range_min, range_max, ax3)
    
    fig.tight_layout()

    if(save):
        
        outname = save_dir + NEXRAD_dict['radar_name'] + \
            '_multi_myrhi_' + str(azimuth) + '_' + \
            NEXRAD_dict['radar_date'] + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

def plot_NEXRAD_rhi_multipanel_auto_varyazm(date_str, radar, min_azm, \
        max_azm, delta_azm = 3, angle_idx = 0, range_min = 0, range_max = 190, \
        height_lim = 8, save_dir = './', save = True):

    # Convert the input date_str to datetime
    # --------------------------------------
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # See if the output directory exists
    # ----------------------------------
    total_save_dir = save_dir + 'ppi_rhi/' + radar + '/vary_azm/' + date_str + '/'
        
    try:
        Path(total_save_dir).mkdir(parents = True, exist_ok = False)
    except FileExistsError:
        print("Save directory already exists. Overwriting images there")

    # Set up the range of azimuth angles
    # ----------------------------------
    azimuths = np.arange(min_azm, max_azm + 1, delta_azm)

    print(total_save_dir)

    for azm in azimuths:
        print(azm)
        plot_NEXRAD_rhi_multipanel(date_str, radar, azm, \
            angle_idx = angle_idx, range_min = range_min, range_max = range_max, \
            save_dir = total_save_dir, height_lim = height_lim, \
            save = True)

# delta_time : minutes
def plot_NEXRAD_rhi_multipanel_auto_varytime(radar, azm, \
        begin_str, end_str, delta_time = 30, \
        angle_idx = 0, range_min = 0, range_max = 190, \
        save_dir = './', save = True):

    # Convert the input date_str to datetime
    # --------------------------------------
    dt_begin_str = datetime.strptime(begin_str,"%Y%m%d%H%M")
    dt_end_str   = datetime.strptime(end_str,"%Y%m%d%H%M")

    # See if the output directory exists
    # ----------------------------------
    total_save_dir = save_dir + 'ppi_rhi/' + radar + '/vary_time/' + str(azm) + '/'
        
    try:
        Path(total_save_dir).mkdir(parents = True, exist_ok = False)
    except FileExistsError:
        print("Save directory already exists. Overwriting images there")

    # Set up the range of azimuth angles
    # ----------------------------------
    num_times = int(((dt_end_str - dt_begin_str).seconds / 60) / delta_time)

    print(total_save_dir, num_times)

    for ntime in range(num_times):
        local_date = dt_begin_str + timedelta(minutes = int(delta_time * ntime))
        print(local_date)
        local_str = local_date.strftime('%Y%m%d%H%M')      
  
        plot_NEXRAD_rhi_multipanel(local_str, radar, azm, \
            angle_idx = angle_idx, range_min = range_min, range_max = range_max, \
            save_dir = total_save_dir, save = True)


# Plots 
def plot_NEXRAD_multiangle(date_str, radar, variable, \
        ptitle = None, plabel = None, vmin = -5, vmax = 90, \
        labelsize = 10, colorbar = True, counties = True, save_dir = './',\
        alpha = 1.0, mask_outside = True, zoom=True, save=False):

    # Read the NEXRAD data
    # --------------------
    NEXRAD_dict = read_NEXRAD(date_str, radar, angle = 2)

    labelsize = 10
    # Set up the figure
    # -----------------
    plt.close('all')
    mapcrs = init_proj(NEXRAD_dict)
    fig = plt.figure(figsize = (13, 6))
    ax1  = fig.add_subplot(2,5,1, projection = mapcrs)
    ax2  = fig.add_subplot(2,5,2, projection = mapcrs)
    ax3  = fig.add_subplot(2,5,3, projection = mapcrs)
    ax4  = fig.add_subplot(2,5,4, projection = mapcrs)
    ax5  = fig.add_subplot(2,5,5, projection = mapcrs)
    ax6  = fig.add_subplot(2,5,6, projection = mapcrs)
    ax7  = fig.add_subplot(2,5,7, projection = mapcrs)
    ax8  = fig.add_subplot(2,5,8, projection = mapcrs)
    ax9  = fig.add_subplot(2,5,9, projection = mapcrs)
    ax10 = fig.add_subplot(2,5,10, projection = mapcrs)

    plot_NEXRAD_ppi(NEXRAD_dict, 'composite_reflectivity', ax = ax1, angle = None, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = mapcrs)
    plot_NEXRAD_ppi(NEXRAD_dict, 'reflectivity', ax = ax2, angle = 0, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = mapcrs)

    plot_NEXRAD_ppi(NEXRAD_dict, 'reflectivity', ax = ax3, angle = 1, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = mapcrs)
    plot_NEXRAD_ppi(NEXRAD_dict, 'reflectivity', ax = ax4, angle = 2, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = mapcrs)

    plot_NEXRAD_ppi(NEXRAD_dict, 'reflectivity', ax = ax5, angle = 3, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = mapcrs)
    plot_NEXRAD_ppi(NEXRAD_dict, 'reflectivity', ax = ax6, angle = 4, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = mapcrs)

    plot_NEXRAD_ppi(NEXRAD_dict, 'reflectivity', ax = ax7, angle = 5, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = mapcrs)
    plot_NEXRAD_ppi(NEXRAD_dict, 'reflectivity', ax = ax8, angle = 6, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = mapcrs)

    plot_NEXRAD_ppi(NEXRAD_dict, 'reflectivity', ax = ax9, angle = 7, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = mapcrs)
    plot_NEXRAD_ppi(NEXRAD_dict, 'reflectivity', ax = ax10, angle = 8, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = mapcrs)

    plt.suptitle(date_str)

    fig.tight_layout()

    del(NEXRAD_dict)
    #del(var0)

    if(save):
        outname = 'nexrad_multiangle_' + radar + '_' + date_str + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

# Plots GOES data, KBBX and KRGX radar data
def plot_NEXRAD_multiradar_multiangle(date_str, radar1, radar2, variable, \
        ptitle = None, plabel = None, vmin = -5, vmax = 90, \
        labelsize = 10, colorbar = True, counties = True, save_dir = './',\
        alpha = 1.0, mask_outside = True, zoom=True, save=False):

    # Read the NEXRAD data
    # --------------------
    NEXRAD_dict1 = read_NEXRAD(date_str, radar1, angle = 2)
    NEXRAD_dict2 = read_NEXRAD(date_str, radar2, angle = 2)

    labelsize = 10
    # Set up the figure
    # -----------------
    plt.close('all')
    mapcrs = init_proj(NEXRAD_dict1)
    fig = plt.figure(figsize = (13, 6))
    ax1  = fig.add_subplot(2,5,1, projection = mapcrs)
    ax2  = fig.add_subplot(2,5,6, projection = mapcrs)
    ax3  = fig.add_subplot(2,5,2, projection = mapcrs)
    ax4  = fig.add_subplot(2,5,7, projection = mapcrs)
    ax5  = fig.add_subplot(2,5,3, projection = mapcrs)
    ax6  = fig.add_subplot(2,5,8, projection = mapcrs)
    ax7  = fig.add_subplot(2,5,4, projection = mapcrs)
    ax8  = fig.add_subplot(2,5,9, projection = mapcrs)
    ax9  = fig.add_subplot(2,5,5, projection = mapcrs)
    ax10 = fig.add_subplot(2,5,10, projection = mapcrs)

    plot_NEXRAD_ppi(NEXRAD_dict1, 'composite_reflectivity', ax = ax1, angle = None, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = mapcrs)
    plot_NEXRAD_ppi(NEXRAD_dict2, 'composite_reflectivity', ax = ax2, angle = None, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = mapcrs)

    plot_NEXRAD_ppi(NEXRAD_dict1, 'reflectivity', ax = ax3, angle = 0, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = mapcrs)
    plot_NEXRAD_ppi(NEXRAD_dict2, 'reflectivity', ax = ax4, angle = 0, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = mapcrs)

    plot_NEXRAD_ppi(NEXRAD_dict1, 'reflectivity', ax = ax5, angle = 3, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = mapcrs)
    plot_NEXRAD_ppi(NEXRAD_dict2, 'reflectivity', ax = ax6, angle = 3, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = mapcrs)

    plot_NEXRAD_ppi(NEXRAD_dict1, 'reflectivity', ax = ax7, angle = 6, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = mapcrs)
    plot_NEXRAD_ppi(NEXRAD_dict2, 'reflectivity', ax = ax8, angle = 6, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = mapcrs)

    plot_NEXRAD_ppi(NEXRAD_dict1, 'reflectivity', ax = ax9, angle = 9, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = mapcrs)
    plot_NEXRAD_ppi(NEXRAD_dict2, 'reflectivity', ax = ax10, angle = 9, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside, crs = mapcrs)

    plt.suptitle(date_str)

    fig.tight_layout()

    del(NEXRAD_dict1)
    del(NEXRAD_dict2)
    #del(var0)

    if(save):
        outname = 'nexrad_multipanel_' + radar1 + '_' + radar2 + '_' + date_str + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()


#def plot_NEXRAD_GOES(NEXRAD_dict, variable, channel, ax = None, \
def plot_NEXRAD_GOES(date_str, radar, variable, channel, ax = None, \
        angle = 0, ptitle = None, plabel = None, vmin = -5, vmax = 90, \
        labelsize = 10, colorbar = True, counties = True, save_dir = './',\
        alpha = 1.0, mask_outside = True, zoom=True, save=False):

    if('/home/bsorenson/Research/GOES' not in sys.path):
        sys.path.append('/home/bsorenson/Research/GOES')
    from GOESLib import read_GOES_satpy, plot_GOES_satpy
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    #dt_date_str = datetime.strptime(NEXRAD_dict['radar_date'],"%Y%m%d%H%M")

    # Read the NEXRAD data
    # --------------------
    NEXRAD_dict = read_NEXRAD(date_str, radar, angle = angle)

    # Read the GOES data
    # ------------------
    var0, crs0, lons0, lats0, lat_lims, lon_lims, plabel_goes = \
        read_GOES_satpy(NEXRAD_dict['radar_date'], channel)

    labelsize = 10
    # Set up the figure
    # -----------------
    plt.close('all')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection = crs0)
    plot_GOES_satpy(NEXRAD_dict['radar_date'], channel, ax = ax, var = var0, \
        crs = crs0, lons = lons0, lats = lats0, lat_lims = lat_lims, \
        lon_lims = lon_lims, vmin = 5, vmax = 80, ptitle = '', \
        plabel = plabel_goes, colorbar = True, labelsize = labelsize + 1, \
        zoom=True,save=False) 
    plot_NEXRAD_ppi(NEXRAD_dict, variable, ax = ax, angle = angle, \
        counties = counties, vmin = vmin, vmax = vmax, alpha = alpha, \
        mask_outside = mask_outside)

    plt.show()

def plot_NEXRAD_ppi_figure(date_str, radar, variable, angle = 0, \
    save_dir = './', vmin = None, vmax = None,\
    mask_outside = True, zoom = True, save = False):

    # Read in the NEXRAD data
    # -----------------------
    NEXRAD_dict = read_NEXRAD(date_str, radar, angle = angle)

    # Plot the NEXRAD data
    # --------------------
    plot_NEXRAD_ppi(NEXRAD_dict, variable, angle = angle, save_dir = save_dir,\
        vmin = vmin, vmax = vmax, \
        mask_outside = mask_outside, zoom = zoom, save = save) 
    

# variable must be 'reflectivity', 'velocity', 'differential_reflectivity',
# or 'cross_correlation_ratio'
def plot_NEXRAD_ppi(NEXRAD_dict, variable, ax = None, angle = None, ptitle = None, \
        plabel = None, vmin = None, vmax = None, \
        labelsize = 10, colorbar = True, \
        counties = True, save_dir = './',\
        alpha = 1.0, mask_outside = True, crs = None, zoom=True, save=False):

    dt_date_str = datetime.strptime(NEXRAD_dict['radar_date'],"%Y%m%d%H%M")

    if(vmin is None):
        vmin = min_dict[variable]['ppi']
    if(vmax is None):
        vmax = max_dict[variable]['ppi']

    # Plot the NEXRAD data
    # ------------------
    in_ax = True 
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure()
        mapcrs = init_proj(NEXRAD_dict)
        crs = mapcrs
        #ax = fig.add_subplot(1,1,1)
        ax = fig.add_subplot(1,1,1, projection = mapcrs)

    # Begin section adapted from Alec Sczepanski's ASGSA PyArt Tools of the 
    # Trade
    # ---------------------------------------------------------------------
    gatefilter = pyart.filters.GateFilter(NEXRAD_dict['radar'])
    gatefilter.exclude_below('reflectivity', -20)
    if(variable == 'composite_reflectivity'):
        compz, z_stack = composite_reflectivity(NEXRAD_dict['radar'], \
            field = 'reflectivity', gatefilter = gatefilter)
        gatefilter = pyart.filters.GateFilter(compz)
        gatefilter.exclude_below('composite_reflectivity', -20)
        #compz = pyart.retrieve.composite_reflectivity(\
        #    NEXRAD_dict['radar'], field = 'reflectivity_horizontal', \
        #    gatefilter = gatefilter)
        display = pyart.graph.RadarMapDisplay(compz)
    else:
        display = pyart.graph.RadarMapDisplay(NEXRAD_dict['radar'])

    if(counties):
        ax.add_feature(USCOUNTIES.with_scale('5m'), alpha = 0.25)

    if(angle is None):
        angle = NEXRAD_dict['angle_idx']

    if(colorbar):
        colorbar_flag = 1
    else:
        colorbar_flag = 0

    if(ptitle is None):
        ptitle = NEXRAD_dict['figure_title']
        #ptitle = variable.title() + '\n' + \
        #    dt_date_str.strftime('%Y%m%d %H:%M')
        ptitle = NEXRAD_dict['radar_name'] + \
            ' - Elevation = %.3f' % (NEXRAD_dict['angles'][angle]) + \
            '\n' + dt_date_str.strftime('%Y%m%d %H:%M')

    #if(variable != 'composite_reflectivity'):
    gatefilter.exclude_below(variable, vmin)

    if(crs is None):
        crs = datacrs

    if(variable == 'composite_reflectivity'):
        #display.plot_ppi_map('composite_reflectivity', ax = ax, vmin = None, vmax = None, \
        #    cmap = 'pyart_HomeyerRainbow', gatefilter = gatefilter, projection = crs)
        display.plot_ppi_map(variable, # Variable name
            vmin = vmin,  # min value
            vmax = vmax,  # max value
            min_lon = NEXRAD_dict['min_lon'], # Min longitude of plot
            max_lon = NEXRAD_dict['max_lon'], # Max longitude of plot
            min_lat = NEXRAD_dict['min_lat'], # Min latitude of plot
            max_lat = NEXRAD_dict['max_lat'], # Max latitude of plot
            resolution = '10m', # Projection resolution
            projection = crs, # Set projection
            #projection = mapcrs, # Set projection
            colorbar_flag = colorbar_flag, # Turn on colorbar?
            #fig = fig, # Where to plot data
            ax = ax, 
            lat_0 = NEXRAD_dict['center_lat'], 
            lon_0 = NEXRAD_dict['center_lon'],
            cmap = cmap_dict[variable], 
            title = ptitle,
            colorbar_label  = title_dict[variable],\
            gatefilter = gatefilter, \
            alpha = alpha, \
            mask_outside = mask_outside)
    else:
        display.plot_ppi_map(variable, # Variable name
            angle, # Elevation angle
            vmin = vmin,  # min value
            vmax = vmax,  # max value
            min_lon = NEXRAD_dict['min_lon'], # Min longitude of plot
            max_lon = NEXRAD_dict['max_lon'], # Max longitude of plot
            min_lat = NEXRAD_dict['min_lat'], # Min latitude of plot
            max_lat = NEXRAD_dict['max_lat'], # Max latitude of plot
            resolution = '10m', # Projection resolution
            projection = crs, # Set projection
            #projection = mapcrs, # Set projection
            colorbar_flag = colorbar_flag, # Turn on colorbar?
            #fig = fig, # Where to plot data
            ax = ax, 
            lat_0 = NEXRAD_dict['center_lat'], 
            lon_0 = NEXRAD_dict['center_lon'],
            cmap = cmap_dict[variable], 
            title = ptitle,
            colorbar_label  = title_dict[variable],\
            gatefilter = gatefilter, \
            alpha = alpha, \
            mask_outside = mask_outside)

    #ax.add_feature(cfeature.STATES)

    # Zoom in the figure if desired
    # -----------------------------
    if(zoom):
        lat_lims = [39.5, 42.0]
        lon_lims = [-122.0, -119.5]
        ax.set_extent([lon_lims[0],lon_lims[1],lat_lims[0],lat_lims[1]],\
                       crs = ccrs.PlateCarree())
        zoom_add = '_zoom'
    else:
        zoom_add = ''

    if(not in_ax): 
        if(save):
            outname = save_dir + NEXRAD_dict['radar_name'] + '_' + variable + '_angle' \
                + str(angle) + '_' +  NEXRAD_dict['radar_date'] + zoom_add + \
                '.png'
            plt.savefig(outname,dpi=200)
            print("Saved image",outname)
        else:
            plt.show()

def calc_my_rhi(NEXRAD_dict, azimuth, variable, az_tol = 0.5):
    avg_rhi = np.full((NEXRAD_dict['radar'].nsweeps, NEXRAD_dict['radar'].ngates), np.nan)
    rng_rhi = np.full((NEXRAD_dict['radar'].nsweeps, NEXRAD_dict['radar'].ngates), np.nan)
    elv_rhi = np.full((NEXRAD_dict['radar'].nsweeps, NEXRAD_dict['radar'].ngates), np.nan)

    radar_height = NEXRAD_dict['radar'].altitude['data'][0] / 1e3
    for ii, sweep_slice in enumerate(NEXRAD_dict['radar'].iter_slice()):
        start_idx = sweep_slice.start
        end_idx   = start_idx+NEXRAD_dict['radar'].get_azimuth(ii).shape[0]
        az_idxs_swp = np.where( np.abs(NEXRAD_dict['radar'].get_azimuth(ii) - azimuth) < az_tol)[0]
        avg_rhi[ii,:] = np.nanmean(NEXRAD_dict['radar'].fields[variable]['data'][\
           start_idx:end_idx][az_idxs_swp,:], axis = 0)
        elv_rhi[ii,:] = (np.nanmean(NEXRAD_dict['radar'].gate_altitude['data'][start_idx:end_idx,:][\
           az_idxs_swp,:], axis = 0) / 1e3)
        rng_rhi[ii,:] = NEXRAD_dict['radar'].range['data'] / 1e3

    return avg_rhi, rng_rhi, elv_rhi

    #NEXRAD_dict['avg_rhi'] = avg_rhi
    #NEXRAD_dict['rng_rhi'] = rng_rhi
    #NEXRAD_dict['elv_rhi'] = elv_rhi
    #return NEXRAD_dict

def plot_my_rhi(NEXRAD_dict, variable, azimuth,  az_tol = 0.2, \
        range_min = 0, range_max = 100, height_lim = 8,\
        vmin = None, vmax = None, \
        ax = None
        ):
    
    #NEXRAD_dict = calc_my_rhi(NEXRAD_dict, azimuth, az_tol = az_tol)
    avg_rhi, rng_rhi, elv_rhi = calc_my_rhi(NEXRAD_dict, azimuth, variable, az_tol = az_tol)

    if(vmin is None):
        vmin = min_dict[variable]['rhi']
    if(vmax is None):
        vmax = max_dict[variable]['rhi']

    # Plot the NEXRAD data
    # ------------------
    in_ax = True 
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

    mask_rhi = np.ma.masked_invalid(avg_rhi)
    mesh = ax.pcolormesh(rng_rhi, elv_rhi, mask_rhi, vmin = vmin, vmax = vmax, \
        cmap = cmap_dict[variable], shading = 'auto') 
    cbar = plt.colorbar(mesh, ax = ax, label = title_dict[variable])
    ax.set_xlim(range_min, range_max)
    ax.set_ylim(0, height_lim)
    ax.set_xlabel('Distance from radar (km)')
    ax.set_ylabel('Height above sea level (km)')
    ax.set_title(NEXRAD_dict['radar_name'] + ' ' + str(azimuth) + '$^{o}$ RHI')

    if(not in_ax): 
        if(save):
            outname = save_dir + NEXRAD_dict['radar_name'] + '_' + variable + \
                '_myrhi' + str(int(azimuth[0])) + '_' + \
                NEXRAD_dict['radar_date'] + '.png'
            plt.savefig(outname,dpi=200)
            print("Saved image",outname)
        else:
            plt.show()


# plot_type: 'default' for the normal pyart version
#            anything else for my version
def plot_NEXRAD_rhi_figure(date_str, radar, variable, azimuth = 0, \
    save_dir = './', vmin = None, vmax = None,\
    range_min = 0, range_max = 100, \
    az_tol = 0.5, \
    mask_outside = True, zoom = True, \
    plot_type = 'default', save = False):

    # Read in the NEXRAD data
    # -----------------------
    NEXRAD_dict = read_NEXRAD(date_str, radar)

    # Plot the NEXRAD data
    # --------------------
    if(plot_type == 'default'):
        plot_NEXRAD_rhi(NEXRAD_dict, variable, azimuth, save_dir = save_dir,\
            vmin = vmin, vmax = vmax, range_min = range_min, \
            range_max = range_max, zoom = zoom, save = save)
    else:
        plot_my_rhi(NEXRAD_dict, azimuth, variable, az_tol = az_tol, \
            range_min = range_min, range_max = range_max, save = save)
 
# variable must be 'reflectivity', 'velocity', 'differential_reflectivity',
# or 'cross_correlation_ratio'
# azimuth may be a single number or a list of numbers
# -------------------------------------------------------------------------
def plot_NEXRAD_rhi(NEXRAD_dict, variable, azimuth, ax = None, angle = None, \
        ptitle = None, plabel = None, vmin = None, vmax = None, \
        az_tol = 1.00, \
        labelsize = 10, range_min = 0, range_max = 100, height_lim = 8,\
        colorbar = True, counties = True, save_dir = './',\
        alpha = 1.0, zoom=True, save=False):

    dt_date_str = datetime.strptime(NEXRAD_dict['radar_date'],"%Y%m%d%H%M")

    if(not isinstance(azimuth, list)):
        azimuth = [azimuth]

    if(vmin is None):
        vmin = min_dict[variable]['rhi']
    if(vmax is None):
        vmax = max_dict[variable]['rhi']

    # Plot the NEXRAD data
    # ------------------
    in_ax = True 
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure()
        num_row = len(azimuth)
        #axs = fig.subplots(nrows = num_row, ncols = 1)

    xsect = pyart.util.cross_section_ppi(NEXRAD_dict['radar'], azimuth, \
        az_tol = az_tol)

    if(colorbar):
        colorbar_flag = 1
    else:
        colorbar_flag = 0

    #gatefilter = pyart.filters.GateFilter(NEXRAD_dict['radar'])
    #gatefilter.exclude_below(variable, vmin)

    #display = pyart.graph.RadarDisplay(xsect)
    #display.plot(variable, 0, vmin = vmin, vmax = vmax)
    #ax.set_xlim(range_min, range_max)
    #ax.set_ylim(0, height_lim)

    display = pyart.graph.RadarDisplay(xsect)
    if(not in_ax):
        for ii, azm in enumerate(azimuth):
            ax = fig.add_subplot(num_row, 1, ii + 1) 
            display.plot(variable, ii, vmin = vmin, vmax = vmax, cmap = cmap_dict[variable])
            ax.set_xlim(range_min, range_max)
            ax.set_ylim(0, height_lim)
    else:
        display.plot(variable, 0, vmin = vmin, vmax = vmax, ax = ax, cmap = cmap_dict[variable])
        ax.set_xlim(range_min, range_max)
        ax.set_ylim(0, height_lim)
    
    #fig.tight_layout()

    if(not in_ax): 
        if(save):
            outname = save_dir + NEXRAD_dict['radar_name'] + '_' + variable + \
                '_cross' + str(int(azimuth[0])) + '_' + \
                NEXRAD_dict['radar_date'] + '.png'
            plt.savefig(outname,dpi=300)
            print("Saved image",outname)
        else:
            plt.show()

# This function gets the NEXRAD data from a higher-resolution channel
# and co-locates it with data from a lower-resolution channel.
def get_NEXRAD_data_lat_lon(date_str, dlat, dlon, channel, version = 0,\
        verbose = False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    var0, _, lons0, lats0, _, _, _ = read_NEXRAD_satpy(date_str, channel)

    #cd_idx = nearest_gridpoint(dlat, dlon,lats0, lons0)
    #nexrad_val = np.array(var0)[cd_idx]

    #for tdlat, tdlon in zip(dlat, dlon):
    if(version == 0):
        if(not isinstance(dlat, list)):
            dlat = [dlat]
            dlon = [dlon]
        cd_idx = [nearest_gridpoint(tdlat, tdlon, \
            lats0, lons0) for tdlat, tdlon in zip(dlat, dlon)]
        nexrad_val = np.array([np.array(var0)[cidx] \
            for cidx in cd_idx]).squeeze()
        nexrad_lat = np.array([np.array(lats0)[cidx] \
            for cidx in cd_idx]).squeeze()
        nexrad_lon = np.array([np.array(lons0)[cidx] \
            for cidx in cd_idx]).squeeze()
    else:

        arr_var0 = np.array(var0)
        arr_lats0 = np.array(lats0)
        arr_lons0 = np.array(lons0)

        nexrad_val = np.array([[arr_var0[nearest_gridpoint(\
            tlat1, tlon1, lats0, lons0)] for \
            tlat1, tlon1 in zip(tlat2, tlon2)] for \
            tlat2, tlon2 in zip(dlat, dlon)]).squeeze()
        nexrad_lat = np.array([[arr_lats0[nearest_gridpoint(\
            tlat1, tlon1, lats0, lons0)] for \
            tlat1, tlon1 in zip(tlat2, tlon2)] for \
            tlat2, tlon2 in zip(dlat, dlon)]).squeeze()
        nexrad_lon = np.array([[arr_lons0[nearest_gridpoint(\
            tlat1, tlon1, lats0, lons0)] for \
            tlat1, tlon1 in zip(tlat2, tlon2)] for \
            tlat2, tlon2 in zip(dlat, dlon)]).squeeze()

        ##!#nexrad_val = np.full(len(cd_idx), np.nan)
        ##!#nexrad_lat = np.full(len(cd_idx), np.nan)
        ##!#nexrad_lon = np.full(len(cd_idx), np.nan)
       
        ##!#for ii, cidx in enumerate(cd_idx):
        ##!#    if(verbose):
        ##!#        print(ii)
        ##!#    nexrad_val[ii] = np.array(var0)[cd_idx[ii]][0] 
        ##!#    nexrad_lat[ii] = np.array(lats0)[cd_idx[ii]][0] 
        ##!#    nexrad_lon[ii] = np.array(lons0)[cd_idx[ii]][0] 

    return nexrad_val, nexrad_lat, nexrad_lon

# NOTE: This is for playing around with find_plume_NEXRAD
def plot_NEXRAD_scatter(date_str):

    # Read in the base channel, usually Ch13
    # --------------------------------------
    channel1 = 13
    var, crs, lons, lats, lat_lims, lon_lims, plabel = \
        read_NEXRAD_satpy(date_str, channel1)

    # Use find_plume_NEXRAD to both co-locate the high-res
    # visible (and low-res ch6 data) with the ch13 data
    # AND mask where the data are not in the plume
    # 
    # NOTE: As of 6/10/2022, find_plume_NEXRAD does not work
    # as well as find_plume from MODISLib, mainly because
    # the 2.25 micron limit used with the  MODIS data does not
    # apply here.
    # --------------------------------------------------------
    hash_data, low_var2, low_var6 = find_plume_NEXRAD(date_str)

    # Make a plot of something
    # ------------------------
    plt.close('all')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    data = [low_var2.flatten(), low_var6.flatten()]
    ax.boxplot(data)
    plt.show()

   
def plot_NEXRAD_ppi_auto(begin_date, end_date, radar, variable, save_dir = './', \
        angle_idx = 0, zoom = True, save = False):

    # Convert the input date_str to datetime
    # --------------------------------------
    begin_dt_date = datetime.strptime(begin_date,"%Y%m%d%H%M")
    end_dt_date   = datetime.strptime(end_date,"%Y%m%d%H%M")

    # Find all downloaded NEXRAD filenames that are between these
    # two dates
    # ---------------------------------------------------------
    all_files = sorted(glob('/home/bsorenson/data/NEXRAD/' + radar + '/*'))
    all_dates = [datetime.strptime(ffile.strip().split('/')[-1][4:19],\
        '%Y%m%d_%H%M%S') for ffile in all_files]
  
    # Get just the file dates, with only 1 date per channel
    # ----------------------------------------------------- 
    in_all = np.array([((tdate > begin_dt_date) & (tdate < end_dt_date)) \
        for tdate in all_dates])
    
    in_all_idx = np.where(in_all == True)

    if(len(in_all_idx[0]) == 0):
        print("ERROR: no data between",begin_date, end_date)
        print("     Run auto_NEXRAD_download to get the data")
        return

    # good_all_files contains the filenames for all files 
    # that have timestamps within the desired range
    # ---------------------------------------------------
    good_all_files    = np.array(all_files)[in_all_idx]
    good_dates        = np.array(all_dates)[in_all_idx]

    # See if the output directory exists
    # ----------------------------------
    total_save_dir = save_dir + radar + '/' + variable + \
        '/angle' + str(angle_idx) + '/'
    try:
        Path(total_save_dir).mkdir(parents = True, exist_ok = False)
    except FileExistsError:
        print("Save directory already exists. Overwriting images there")

    for ttime in good_dates:
        print(ttime.strftime('%Y%m%d%H%M'))
        date_str = ttime.strftime('%Y%m%d%H%M')
        if(date_str == '202107202126'):
            print("Not making image for this time")
        else:
            NEXRAD_dict = read_NEXRAD(date_str, radar, angle = angle_idx)
            plot_NEXRAD_ppi(NEXRAD_dict, variable, zoom = True, \
                save_dir = total_save_dir, \
                angle = angle_idx, save = True)

# Plots azimuths between a min and max value for a single date
# (Can be modified to do multiple dates)
# ------------------------------------------------------------
def plot_NEXRAD_rhi_auto(date_str, radar, variable, \
        save_dir = './', min_azm = 0, max_azm = 50, range_lim = 100, \
        height_lim = 8, zoom = True, save = False):

    # Convert the input date_str to datetime
    # --------------------------------------
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # See if the output directory exists
    # ----------------------------------
    total_save_dir = save_dir + 'rhi/' + radar + '/' + variable + '/'
        
    try:
        Path(total_save_dir).mkdir(parents = True, exist_ok = False)
    except FileExistsError:
        print("Save directory already exists. Overwriting images there")

    # Set up the range of azimuth angles
    # ----------------------------------
    azimuths = np.arange(min_azm, max_azm + 1)

    for azm in azimuths:
        print(azm)
        NEXRAD_dict = read_NEXRAD(date_str, radar)
        plot_NEXRAD_rhi(NEXRAD_dict, variable, azm, \
            ptitle = None, plabel = None, vmin = -5, vmax = 90, \
            range_lim = range_lim, height_lim = height_lim, \
            labelsize = 10, colorbar = True, \
            save_dir = total_save_dir,save= True)

def plot_NEXRAD_time_series_points(NEXRAD_dict, time_idx = 20, \
        ch_idx = 0, save_dir = './', save = False):

    # Read NEXRAD image data
    # --------------------
    print(NEXRAD_dict['dt_dates'][time_idx])
    date_str = NEXRAD_dict['dt_dates'][time_idx].strftime('%Y%m%d%H%M')
    var, crs, lons, lats, lat_lims, lon_lims, plabel = \
        read_NEXRAD_satpy(date_str, NEXRAD_dict['channels'][ch_idx])

    plt.close('all')
    fig = plt.figure(figsize = (11, 3.5))
    gs = fig.add_gridspec(nrows = 1, ncols = 4)
    ax1  = fig.add_subplot(gs[3], projection = crs)   # true color    
    ax2  = fig.add_subplot(gs[0:3]) # Ch 1
    #ax1 = fig.add_subplot(1,2,1, projection = crs)
    #ax2 = fig.add_subplot(1,2,2)

    labelsize = 10
    font_size = 8
    plot_NEXRAD_satpy(date_str, NEXRAD_dict['channels'][ch_idx], \
        ax = ax1, var = var, crs = crs, \
        lons = lons, lats = lats, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = nexrad_channel_dict[\
            str(NEXRAD_dict['channels'][ch_idx])]['limits'][0], \
        vmax = nexrad_channel_dict[\
            str(NEXRAD_dict['channels'][ch_idx])]['limits'][1], \
        ptitle = '', plabel = plabel, \
        #vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)

    # Add the locations of the points
    # -------------------------------
    for tlat, tlon in zip(NEXRAD_dict['plat'], NEXRAD_dict['plon']):
        print(tlat, tlon)
        ax1.plot(tlon, tlat, linewidth=2, markersize = 8, marker='.',
                 color = 'black', transform=datacrs)
        ax1.plot(tlon, tlat, linewidth=2, markersize = 5, marker='.',
                 transform=datacrs)
    plot_figure_text(ax1, 'NEXRAD-17 ' + \
        str(nexrad_channel_dict[\
        str(NEXRAD_dict['channels'][ch_idx])]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    ax1.set_title(NEXRAD_dict['dt_dates'][time_idx].strftime(\
        '%Y-%m-%d %H:%MZ'), fontsize = 10)

    # Plot the first 2 channels
    # -------------------------
    for jj in range(NEXRAD_dict['data'].shape[2]):
        ax2.plot(NEXRAD_dict['dt_dates'], NEXRAD_dict['data'][:,ch_idx,jj])
    #ax2.plot(NEXRAD_dict['dt_dates'], NEXRAD_dict['data'][:,ch_idx,1])
    #ax2.plot(NEXRAD_dict['dt_dates'], NEXRAD_dict['data'][:,ch_idx,2])
    ax2.axvline(NEXRAD_dict['dt_dates'][time_idx], color = 'black',\
        linestyle = ':')
    ax2.set_ylabel(plabel.replace('_',' '), \
        size = labelsize, weight = 'bold')
    ax2.grid()
    ax2.xaxis.set_major_formatter(DateFormatter('%m/%d\n%H:%MZ'))
    ax2.tick_params(axis="x", labelsize = 9)
    ##!#    str(NEXRAD_dict['channels'][ch_idx])]['wavelength']) + ' μm', \
    ts_title = 'NEXRAD-17 ' + \
        nexrad_channel_dict[\
        str(NEXRAD_dict['channels'][ch_idx])]['name'] + \
        ' (' + \
        str(nexrad_channel_dict[\
        str(NEXRAD_dict['channels'][ch_idx])]['wavelength']) + \
        ' μm)'
        # + ' '.join(plabel.split()[0].split('_'))
    ax2.set_title(ts_title, fontsize = 10)
    fig.autofmt_xdate()

    fig.tight_layout()

    if(save):
        outname = save_dir + 'nexrad_time_series_points_ch' + \
            str(NEXRAD_dict['channels'][ch_idx]) + '_' + \
            NEXRAD_dict['dt_dates'][time_idx].strftime('%Y%m%d%H%M') + \
            '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

def plot_NEXRAD_time_series_points_auto(NEXRAD_dict, ch_idx, \
        save_dir = './'):

    # Loop over all the dates and save the images
    # -------------------------------------------
    for ii in range(len(NEXRAD_dict['dt_dates'])):
        if(NEXRAD_dict['dt_dates'][ii].strftime('%Y%m%d%H%M') == \
            '202107202126'):
            continue
        else:
            plot_NEXRAD_time_series_points(NEXRAD_dict, time_idx = ii, \
                ch_idx = ch_idx, save_dir = save_dir, save = True)

#def plot_NEXRAD_time_series_points(NEXRAD_dict, time_idx = 20, \
#        ch_idx = 0, save = False):

#def plot_NEXRAD_time_series(dt_dates, nexrad_data, channels, save = False):
def plot_NEXRAD_time_series_channels(NEXRAD_dict, time_idx = 20, \
        idx = 0, save_dir = './', save = False):

    plt.close('all')
    fig = plt.figure(figsize = (10, 6))
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)

    # Plot the first 2 channels
    # -------------------------
    ln1 = ax1.plot(NEXRAD_dict['dt_dates'], NEXRAD_dict['data'][:,0,idx], \
        label = str(nexrad_channel_dict[\
        str(NEXRAD_dict['channels'][0])]['wavelength']) + \
        ' μm')
    ln2 = ax1.plot(NEXRAD_dict['dt_dates'], NEXRAD_dict['data'][:,1,idx], \
        label = str(nexrad_channel_dict[\
        str(NEXRAD_dict['channels'][1])]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:blue')
    ax12 = ax1.twinx()
    ln3 = ax12.plot(NEXRAD_dict['dt_dates'], NEXRAD_dict['data'][:,2,idx], \
        label = str(nexrad_channel_dict[\
        str(NEXRAD_dict['channels'][2])]['wavelength']) + \
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

    ax2.plot(NEXRAD_dict['dt_dates'], NEXRAD_dict['data'][:,3,idx], \
        #label = str(int(NEXRAD_dict['channels'][3])))
        label = str(nexrad_channel_dict[\
        str(NEXRAD_dict['channels'][3])]['wavelength']) + \
        ' μm')
    ax2.plot(NEXRAD_dict['dt_dates'], NEXRAD_dict['data'][:,4,idx], \
        #label = str(int(NEXRAD_dict['channels'][4])))
        label = str(nexrad_channel_dict[\
        str(NEXRAD_dict['channels'][4])]['wavelength']) + \
        ' μm')
    #ax22 = ax2.twinx()
    ax2.plot(NEXRAD_dict['dt_dates'], NEXRAD_dict['data'][:,5,idx], \
        #label = str(int(NEXRAD_dict['channels'][5])))
        label = str(nexrad_channel_dict[\
        str(NEXRAD_dict['channels'][5])]['wavelength']) + \
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
        outname = save_dir + 'nexrad_time_series_channels_pt' + \
            str(idx) + '.png'
            #NEXRAD_dict['dt_dates'][time_idx].strftime('%Y%m%d%H%M') + \
            #'.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

# This is written to compare time series of any two desired channels
# (as long as they're in the NEXRAD_dict structure)
def plot_NEXRAD_time_series_channel_comp(NEXRAD_dict, ch_idx1, ch_idx2, \
        idx1, idx2, ch_idx3 = None, date_idx = 20, save_dir = './', \
        save = False):

    channel1 = int(NEXRAD_dict['channels'][ch_idx1])
    channel2 = int(NEXRAD_dict['channels'][ch_idx2])
    if(ch_idx3 is not None):
        channel3 = int(NEXRAD_dict['channels'][ch_idx3])
        figsize = (9.5, 3.5)
    else:
        figsize = (8.5, 3.5)
        

    dt_date_str = NEXRAD_dict['dt_dates'][date_idx]
    date_str = dt_date_str.strftime('%Y%m%d%H%M')
    var, crs, lons, lats, lat_lims, lon_lims, plabel = \
        read_NEXRAD_satpy(date_str, channel1)

    plt.close('all')
    fig = plt.figure(figsize = figsize)
    gs = fig.add_gridspec(nrows = 1, ncols = 3)
    ax2  = fig.add_subplot(gs[2], projection = crs)   # true color    
    ax1  = fig.add_subplot(gs[0:2]) # Ch 1

    # Plot the two channel data for the first point
    # ---------------------------------------------
    ln11 = ax1.plot(NEXRAD_dict['dt_dates'], NEXRAD_dict['data'][:,ch_idx1,idx1], \
        label = str(nexrad_channel_dict[\
        str(channel1)]['wavelength']) + \
        ' μm', color = 'tab:blue')
    ln21 = ax1.plot(NEXRAD_dict['dt_dates'], NEXRAD_dict['data'][:,ch_idx1,idx2], \
        label = str(nexrad_channel_dict[\
        str(channel1)]['wavelength']) + \
        ' μm', color = 'tab:orange')

    # Plot the two channel data for the second point
    # ----------------------------------------------
    ln12 = ax1.plot(NEXRAD_dict['dt_dates'], NEXRAD_dict['data'][:,ch_idx2,idx1], \
        label = str(nexrad_channel_dict[\
        str(channel2)]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:blue')
    ln22 = ax1.plot(NEXRAD_dict['dt_dates'], NEXRAD_dict['data'][:,ch_idx2,idx2], \
        label = str(nexrad_channel_dict[\
        str(channel2)]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:orange')
    ax1.axvline(dt_date_str, color = 'black',\
        linestyle = ':')

    lns = ln11 + ln21 + ln12 + ln22 

    if(ch_idx3 is not None):
        ax12 = ax1.twinx()
        ln31 = ax12.plot(NEXRAD_dict['dt_dates'], NEXRAD_dict['data'][:,ch_idx3,idx1], \
            label = str(nexrad_channel_dict[\
            str(channel3)]['wavelength']) + \
            ' μm', linestyle = ':', color = 'tab:blue')
        ln32 = ax12.plot(NEXRAD_dict['dt_dates'], NEXRAD_dict['data'][:,ch_idx3,idx2], \
            label = str(nexrad_channel_dict[\
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
    print(str(nexrad_channel_dict[str(channel1)]['wavelength']) + \
        ' μm\n')
    print("    Point blue   - ", NEXRAD_dict['data'][date_idx,ch_idx1,idx1])
    print("    Point orange - ", NEXRAD_dict['data'][date_idx,ch_idx1,idx2])
    print(str(nexrad_channel_dict[str(channel2)]['wavelength']) + \
        ' μm\n')
    print("    Point blue   - ", NEXRAD_dict['data'][date_idx,ch_idx2,idx1])
    print("    Point orange - ", NEXRAD_dict['data'][date_idx,ch_idx2,idx2])
    if(ch_idx3 is not None):
        print(str(nexrad_channel_dict[str(channel3)]['wavelength']) + \
            ' μm\n')
        print("    Point blue   - ", NEXRAD_dict['data'][date_idx,ch_idx3,idx1])
        print("    Point orange - ", NEXRAD_dict['data'][date_idx,ch_idx3,idx2])

    ##!#print("Point Blue:\n")
    ##!#print("    " + str(nexrad_channel_dict[str(channel1)]['wavelength']) + \
    ##!#    ' μm - ', NEXRAD_dict['data'][date_idx,ch_idx1,idx1])
    ##!#print("    " + str(nexrad_channel_dict[str(channel2)]['wavelength']) + \
    ##!#    ' μm - ', NEXRAD_dict['data'][date_idx,ch_idx2,idx1])
    ##!#print("Point Orange:\n")
    ##!#print("    " + str(nexrad_channel_dict[str(channel1)]['wavelength']) + \
    ##!#    ' μm - ', NEXRAD_dict['data'][date_idx,ch_idx1,idx2])
    ##!#print("    " + str(nexrad_channel_dict[str(channel2)]['wavelength']) + \
    ##!#    ' μm - ', NEXRAD_dict['data'][date_idx,ch_idx2,idx2])

    #lns = ln1+ln2+ln3
    #labs = [l.get_label() for l in lns]
    #ax1.legend(fontsize = font_size)
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc=0, fontsize = font_size)


    # Plot the NEXRAD 0.64 micron image
    # -------------------------------
    plot_NEXRAD_satpy(date_str, channel1, \
        ax = ax2, var = var, crs = crs, \
        lons = lons, lats = lats, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = nexrad_channel_dict[str(channel1)]['limits'][0], \
        vmax = nexrad_channel_dict[str(channel1)]['limits'][1], \
        ptitle = '', plabel = plabel, \
        #vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)

    # Plot the currently-analyzed points
    # ----------------------------------
    point_size = 5
    ax2.plot(NEXRAD_dict['plon'][idx1], NEXRAD_dict['plat'][idx1], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax2.plot(NEXRAD_dict['plon'][idx1], NEXRAD_dict['plat'][idx1], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)
    ax2.plot(NEXRAD_dict['plon'][idx2], NEXRAD_dict['plat'][idx2], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax2.plot(NEXRAD_dict['plon'][idx2], NEXRAD_dict['plat'][idx2], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)

    plot_figure_text(ax2, 'NEXRAD-17 ' + \
        str(nexrad_channel_dict[str(channel1)]['wavelength']) + ' μm', \
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
        outname = save_dir + 'nexrad_time_series_channel_comp_' + \
        NEXRAD_dict['ptype'] + '_ch' + \
            str(channel1) + '_ch' + str(channel2)
        if(ch_idx3 is not None):
            outname = outname + '_ch' + str(channel3)
        outname = outname + '_pt'+ str(idx1) + '_pt' + str(idx2)\
             + '_' + date_str + \
            '.png'
            #NEXRAD_dict['dt_dates'][time_idx].strftime('%Y%m%d%H%M') + \
            #'.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

# This is written to compare time series of any two desired channels
# (as long as they're in the NEXRAD_dict structure)
def plot_NEXRAD_time_series_channel_comp_2loc(NEXRAD_dict1, NEXRAD_dict2, \
        ch_idx1, ch_idx2, idx1, idx2, idx3, idx4, \
        ch_idx3 = None, date_idx = 20, save_dir = './', \
        save = False):

    # Make the plot for panels 1 and 2
    # --------------------------------

    channel1 = int(NEXRAD_dict1['channels'][ch_idx1])
    channel2 = int(NEXRAD_dict1['channels'][ch_idx2])
    if(ch_idx3 is not None):
        channel3 = int(NEXRAD_dict1['channels'][ch_idx3])
        figsize = (11, 7)
    else:
        figsize = (9, 7)
        

    dt_date_str = NEXRAD_dict1['dt_dates'][date_idx]
    date_str = dt_date_str.strftime('%Y%m%d%H%M')
    var, crs, lons, lats, lat_lims, lon_lims, plabel = \
        read_NEXRAD_satpy(date_str, channel1)

    plt.close('all')
    fig = plt.figure(figsize = figsize)
    gs = fig.add_gridspec(nrows = 2, ncols = 3)
    ax2  = fig.add_subplot(gs[0,2], projection = crs)   # true color    
    ax1  = fig.add_subplot(gs[0,0:2]) # Ch 1
    ax4  = fig.add_subplot(gs[1,2], projection = crs)   # true color    
    ax3  = fig.add_subplot(gs[1,0:2]) # Ch 1

    # Plot the two channel data for the first point
    # ---------------------------------------------
    ln11 = ax1.plot(NEXRAD_dict1['dt_dates'], NEXRAD_dict1['data'][:,ch_idx1,idx1], \
        label = str(nexrad_channel_dict[\
        str(channel1)]['wavelength']) + \
        ' μm', color = 'tab:blue')
    ln21 = ax1.plot(NEXRAD_dict1['dt_dates'], NEXRAD_dict1['data'][:,ch_idx1,idx2], \
        label = str(nexrad_channel_dict[\
        str(channel1)]['wavelength']) + \
        ' μm', color = 'tab:orange')

    # Plot the two channel data for the second point
    # ----------------------------------------------
    ln12 = ax1.plot(NEXRAD_dict1['dt_dates'], NEXRAD_dict1['data'][:,ch_idx2,idx1], \
        label = str(nexrad_channel_dict[\
        str(channel2)]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:blue')
    ln22 = ax1.plot(NEXRAD_dict1['dt_dates'], NEXRAD_dict1['data'][:,ch_idx2,idx2], \
        label = str(nexrad_channel_dict[\
        str(channel2)]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:orange')
    ax1.axvline(dt_date_str, color = 'black',\
        linestyle = ':')

    lns = ln11 + ln21 + ln12 + ln22 

    if(ch_idx3 is not None):
        ax12 = ax1.twinx()
        ln31 = ax12.plot(NEXRAD_dict1['dt_dates'], NEXRAD_dict1['data'][:,ch_idx3,idx1], \
            label = str(nexrad_channel_dict[\
            str(channel3)]['wavelength']) + \
            ' μm', linestyle = ':', color = 'tab:blue')
        ln32 = ax12.plot(NEXRAD_dict1['dt_dates'], NEXRAD_dict1['data'][:,ch_idx3,idx2], \
            label = str(nexrad_channel_dict[\
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
    print(str(nexrad_channel_dict[str(channel1)]['wavelength']) + \
        ' μm\n')
    print("    Point blue   - ", NEXRAD_dict1['data'][date_idx,ch_idx1,idx1])
    print("    Point orange - ", NEXRAD_dict1['data'][date_idx,ch_idx1,idx2])
    print(str(nexrad_channel_dict[str(channel2)]['wavelength']) + \
        ' μm\n')
    print("    Point blue   - ", NEXRAD_dict1['data'][date_idx,ch_idx2,idx1])
    print("    Point orange - ", NEXRAD_dict1['data'][date_idx,ch_idx2,idx2])
    if(ch_idx3 is not None):
        print(str(nexrad_channel_dict[str(channel3)]['wavelength']) + \
            ' μm\n')
        print("    Point blue   - ", NEXRAD_dict1['data'][date_idx,ch_idx3,idx1])
        print("    Point orange - ", NEXRAD_dict1['data'][date_idx,ch_idx3,idx2])

    ##!#print("Point Blue:\n")
    ##!#print("    " + str(nexrad_channel_dict[str(channel1)]['wavelength']) + \
    ##!#    ' μm - ', NEXRAD_dict1['data'][date_idx,ch_idx1,idx1])
    ##!#print("    " + str(nexrad_channel_dict[str(channel2)]['wavelength']) + \
    ##!#    ' μm - ', NEXRAD_dict1['data'][date_idx,ch_idx2,idx1])
    ##!#print("Point Orange:\n")
    ##!#print("    " + str(nexrad_channel_dict[str(channel1)]['wavelength']) + \
    ##!#    ' μm - ', NEXRAD_dict1['data'][date_idx,ch_idx1,idx2])
    ##!#print("    " + str(nexrad_channel_dict[str(channel2)]['wavelength']) + \
    ##!#    ' μm - ', NEXRAD_dict1['data'][date_idx,ch_idx2,idx2])

    #lns = ln1+ln2+ln3
    #labs = [l.get_label() for l in lns]
    #ax1.legend(fontsize = font_size)
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc=0, fontsize = font_size)


    # Plot the NEXRAD 0.64 micron image
    # -------------------------------
    plot_NEXRAD_satpy(date_str, channel1, \
        ax = ax2, var = var, crs = crs, \
        lons = lons, lats = lats, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = nexrad_channel_dict[str(channel1)]['limits'][0], \
        vmax = nexrad_channel_dict[str(channel1)]['limits'][1], \
        ptitle = '', plabel = plabel, \
        #vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)

    # Plot the currently-analyzed points
    # ----------------------------------
    point_size = 5
    ax2.plot(NEXRAD_dict1['plon'][idx1], NEXRAD_dict1['plat'][idx1], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax2.plot(NEXRAD_dict1['plon'][idx1], NEXRAD_dict1['plat'][idx1], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)
    ax2.plot(NEXRAD_dict1['plon'][idx2], NEXRAD_dict1['plat'][idx2], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax2.plot(NEXRAD_dict1['plon'][idx2], NEXRAD_dict1['plat'][idx2], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)

    plot_figure_text(ax2, 'NEXRAD-17 ' + \
        str(nexrad_channel_dict[str(channel1)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    ax2.set_title(dt_date_str.strftime(\
        '%Y-%m-%d %H:%MZ'), fontsize = labelsize)

    # Make plot for panels 3 and 4
    # ----------------------------
    # Plot the two channel data for the first point
    # ---------------------------------------------
    ln31 = ax3.plot(NEXRAD_dict2['dt_dates'], NEXRAD_dict2['data'][:,ch_idx1,idx3], \
        label = str(nexrad_channel_dict[\
        str(channel1)]['wavelength']) + \
        ' μm', color = 'tab:blue')
    ln41 = ax3.plot(NEXRAD_dict2['dt_dates'], NEXRAD_dict2['data'][:,ch_idx1,idx4], \
        label = str(nexrad_channel_dict[\
        str(channel1)]['wavelength']) + \
        ' μm', color = 'tab:orange')

    # Plot the two channel data for the second point
    # ----------------------------------------------
    ln32 = ax3.plot(NEXRAD_dict2['dt_dates'], NEXRAD_dict2['data'][:,ch_idx2,idx3], \
        label = str(nexrad_channel_dict[\
        str(channel2)]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:blue')
    ln42 = ax3.plot(NEXRAD_dict2['dt_dates'], NEXRAD_dict2['data'][:,ch_idx2,idx4], \
        label = str(nexrad_channel_dict[\
        str(channel2)]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:orange')
    ax3.axvline(dt_date_str, color = 'black',\
        linestyle = ':')

    lns = ln11 + ln21 + ln12 + ln22 

    if(ch_idx3 is not None):
        ax32 = ax3.twinx()
        ln51 = ax32.plot(NEXRAD_dict2['dt_dates'], NEXRAD_dict2['data'][:,ch_idx3,idx3], \
            label = str(nexrad_channel_dict[\
            str(channel3)]['wavelength']) + \
            ' μm', linestyle = ':', color = 'tab:blue')
        ln52 = ax32.plot(NEXRAD_dict2['dt_dates'], NEXRAD_dict2['data'][:,ch_idx3,idx4], \
            label = str(nexrad_channel_dict[\
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
    print(str(nexrad_channel_dict[str(channel1)]['wavelength']) + \
        ' μm\n')
    print("    Point blue   - ", NEXRAD_dict2['data'][date_idx,ch_idx1,idx3])
    print("    Point orange - ", NEXRAD_dict2['data'][date_idx,ch_idx1,idx4])
    print(str(nexrad_channel_dict[str(channel2)]['wavelength']) + \
        ' μm\n')
    print("    Point blue   - ", NEXRAD_dict2['data'][date_idx,ch_idx2,idx3])
    print("    Point orange - ", NEXRAD_dict2['data'][date_idx,ch_idx2,idx4])
    if(ch_idx3 is not None):
        print(str(nexrad_channel_dict[str(channel3)]['wavelength']) + \
            ' μm\n')
        print("    Point blue   - ", NEXRAD_dict2['data'][date_idx,ch_idx3,idx3])
        print("    Point orange - ", NEXRAD_dict2['data'][date_idx,ch_idx3,idx4])

    ##!#print("Point Blue:\n")
    ##!#print("    " + str(nexrad_channel_dict[str(channel1)]['wavelength']) + \
    ##!#    ' μm - ', NEXRAD_dict2['data'][date_idx,ch_idx1,idx1])
    ##!#print("    " + str(nexrad_channel_dict[str(channel2)]['wavelength']) + \
    ##!#    ' μm - ', NEXRAD_dict2['data'][date_idx,ch_idx2,idx1])
    ##!#print("Point Orange:\n")
    ##!#print("    " + str(nexrad_channel_dict[str(channel1)]['wavelength']) + \
    ##!#    ' μm - ', NEXRAD_dict2['data'][date_idx,ch_idx1,idx2])
    ##!#print("    " + str(nexrad_channel_dict[str(channel2)]['wavelength']) + \
    ##!#    ' μm - ', NEXRAD_dict2['data'][date_idx,ch_idx2,idx2])

    #lns = ln1+ln2+ln3
    #labs = [l.get_label() for l in lns]
    #ax1.legend(fontsize = font_size)
    labs = [l.get_label() for l in lns]
    ax3.legend(lns, labs, loc=0, fontsize = font_size)


    # Plot the NEXRAD 0.64 micron image
    # -------------------------------
    plot_NEXRAD_satpy(date_str, channel1, \
        ax = ax4, var = var, crs = crs, \
        lons = lons, lats = lats, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = nexrad_channel_dict[str(channel1)]['limits'][0], \
        vmax = nexrad_channel_dict[str(channel1)]['limits'][1], \
        ptitle = '', plabel = plabel, \
        #vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)

    # Plot the currently-analyzed points
    # ----------------------------------
    point_size = 5
    ax4.plot(NEXRAD_dict2['plon'][idx3], NEXRAD_dict2['plat'][idx3], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax4.plot(NEXRAD_dict2['plon'][idx3], NEXRAD_dict2['plat'][idx3], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)
    ax4.plot(NEXRAD_dict2['plon'][idx4], NEXRAD_dict2['plat'][idx4], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax4.plot(NEXRAD_dict2['plon'][idx4], NEXRAD_dict2['plat'][idx4], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)

    plot_figure_text(ax4, 'NEXRAD-17 ' + \
        str(nexrad_channel_dict[str(channel1)]['wavelength']) + ' μm', \
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
        outname = save_dir + 'nexrad_time_series_channel_comp_2loc_' + \
            NEXRAD_dict1['ptype'] + '_' + NEXRAD_dict2['ptype'] +  '_ch' + \
            str(channel1) + '_ch' + str(channel2)
        if(ch_idx3 is not None):
            outname = outname + '_ch' + str(channel3)
        outname = outname + '_pt'+ str(idx1) + '_pt' + str(idx2)\
             + '_' + date_str + \
            '.png'
            #NEXRAD_dict['dt_dates'][time_idx].strftime('%Y%m%d%H%M') + \
            #'.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

# This is written to compare time series of any two desired channels
# (as long as they're in the NEXRAD_dict structure). Compares data
# from the same point(s) along the same cross section but from different
# dates.
def plot_NEXRAD_time_series_channel_comp_2loc(NEXRAD_dict1, NEXRAD_dict2, \
        ch_idx1, ch_idx2, idx1, ch_idx3 = None, \
        date_idx = 20, save_dir = './', save = False):

    # Make the plot for panels 1 and 2
    # --------------------------------

    channel1 = int(NEXRAD_dict1['channels'][ch_idx1])
    channel2 = int(NEXRAD_dict1['channels'][ch_idx2])
    if(ch_idx3 is not None):
        channel3 = int(NEXRAD_dict1['channels'][ch_idx3])
        figsize = (8, 7)
    else:
        figsize = (9, 7)
       

    dt_date_str1 = NEXRAD_dict1['dt_dates'][date_idx]
    dt_date_str2 = NEXRAD_dict2['dt_dates'][date_idx + 1]
    date_str1 = dt_date_str1.strftime('%Y%m%d%H%M')
    date_str2 = dt_date_str2.strftime('%Y%m%d%H%M')
    var1, crs1, lons1, lats1, lat_lims1, lon_lims1, plabel1 = \
        read_NEXRAD_satpy(date_str1, channel2)
    var2, crs2, lons2, lats2, lat_lims2, lon_lims2, plabel2 = \
        read_NEXRAD_satpy(date_str2, channel2)

    # Figure out relative dates to plot the same plot
    # -----------------------------------------------
    base_date1 = datetime.strptime('200007130000','%Y%m%d%H%M')
    base_date2 = datetime.strptime('200007200000','%Y%m%d%H%M')
    rel_date1 = NEXRAD_dict1['dt_dates'] - base_date1
    rel_date2 = NEXRAD_dict2['dt_dates'] - base_date2

    plt.close('all')
    fig = plt.figure(figsize = figsize)
    gs = fig.add_gridspec(nrows = 2, ncols = 2)
    ax2  = fig.add_subplot(gs[0,0], projection = crs1)   # true color    
    ax3  = fig.add_subplot(gs[0,1], projection = crs2)   # true color    
    ax1  = fig.add_subplot(gs[1,:]) # Ch 1

    # Plot the two channel data for the first date 
    # ---------------------------------------------
    ln11 = ax1.plot(NEXRAD_dict1['dt_dates'], NEXRAD_dict1['data'][:,ch_idx1,idx1], \
    #ln11 = ax1.plot(rel_date1, NEXRAD_dict1['data'][:,ch_idx1,idx1], \
        label = str(nexrad_channel_dict[\
        str(channel1)]['wavelength']) + \
        ' μm', color = 'tab:blue')
    ax112 = ax1.twiny()
    ln12 = ax112.plot(NEXRAD_dict2['dt_dates'], NEXRAD_dict2['data'][:,ch_idx1,idx1], \
    #ln12 = ax1.plot(rel_date2, NEXRAD_dict2['data'][:,ch_idx1,idx1], \
        label = str(nexrad_channel_dict[\
        str(channel1)]['wavelength']) + \
        ' μm', color = 'tab:orange')


    # Plot the two channel data for the second point
    # ----------------------------------------------
    ln21 = ax1.plot(NEXRAD_dict1['dt_dates'], NEXRAD_dict1['data'][:,ch_idx2,idx1], \
        label = str(nexrad_channel_dict[\
        str(channel2)]['wavelength']) + \
        ' μm', linestyle = '--', color = 'tab:blue')
    ln22 = ax112.plot(NEXRAD_dict2['dt_dates'], NEXRAD_dict2['data'][:,ch_idx2,idx1], \
        label = str(nexrad_channel_dict[\
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
        ln31 = ax12.plot(NEXRAD_dict1['dt_dates'], NEXRAD_dict1['data'][:,ch_idx3,idx1], \
            label = str(nexrad_channel_dict[\
            str(channel3)]['wavelength']) + \
            ' μm', linestyle = ':', color = 'tab:blue')
        ln32 = ax212.plot(NEXRAD_dict2['dt_dates'], NEXRAD_dict2['data'][:,ch_idx3,idx1], \
            label = str(nexrad_channel_dict[\
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
    print(str(nexrad_channel_dict[str(channel1)]['wavelength']) + \
        ' μm\n')
    print("    Date blue   - ", NEXRAD_dict1['data'][date_idx,ch_idx1,idx1])
    print("    Date orange - ", NEXRAD_dict2['data'][date_idx,ch_idx1,idx1])
    print(str(nexrad_channel_dict[str(channel2)]['wavelength']) + \
        ' μm\n')
    print("    Date blue   - ", NEXRAD_dict1['data'][date_idx,ch_idx2,idx1])
    print("    Date orange - ", NEXRAD_dict2['data'][date_idx,ch_idx2,idx1])
    if(ch_idx3 is not None):
        print(str(nexrad_channel_dict[str(channel3)]['wavelength']) + \
            ' μm\n')
        print("    Date blue   - ", NEXRAD_dict1['data'][date_idx,ch_idx3,idx1])
        print("    Date orange - ", NEXRAD_dict2['data'][date_idx,ch_idx3,idx1])

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


    # Plot the NEXRAD 0.64 micron image for date 1
    # ------------------------------------------
    plot_NEXRAD_satpy(date_str1, channel1, \
        ax = ax2, var = var1, crs = crs1, \
        lons = lons1, lats = lats1, lat_lims = lat_lims1, lon_lims = lon_lims1, \
        vmin = nexrad_channel_dict[str(channel2)]['limits'][0], \
        vmax = nexrad_channel_dict[str(channel2)]['limits'][1] - 20, \
        ptitle = '', plabel = plabel1, \
        #vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)

    # Plot the point for date 1
    # -------------------------
    point_size = 5
    ax2.plot(NEXRAD_dict1['plon'][idx1], NEXRAD_dict1['plat'][idx1], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax2.plot(NEXRAD_dict1['plon'][idx1], NEXRAD_dict1['plat'][idx1], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs)

    plot_figure_text(ax2, 'Before Dixie Fire', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right', location = 'upper_right')
    plot_figure_text(ax2, 'NEXRAD-17 ' + \
        str(nexrad_channel_dict[str(channel2)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    ax2.set_title(dt_date_str1.strftime(\
        '%Y-%m-%d %H:%MZ'), fontsize = labelsize)

    # Plot the NEXRAD 0.64 micron image for date 2
    # ------------------------------------------
    plot_NEXRAD_satpy(date_str2, channel1, \
        ax = ax3, var = var2, crs = crs2, \
        lons = lons2, lats = lats2, lat_lims = lat_lims2, lon_lims = lon_lims2, \
        vmin = nexrad_channel_dict[str(channel2)]['limits'][0], \
        vmax = nexrad_channel_dict[str(channel2)]['limits'][1] - 20, \
        ptitle = '', plabel = plabel2, \
        #vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)

    # Plot the currently-analyzed points
    # ----------------------------------
    point_size = 5
    ax3.plot(NEXRAD_dict2['plon'][idx1], NEXRAD_dict2['plat'][idx1], \
            linewidth=2, markersize = point_size + 2, marker='.',
            color = 'black', transform=datacrs)
    ax3.plot(NEXRAD_dict2['plon'][idx1], NEXRAD_dict2['plat'][idx1], \
            linewidth=2, markersize = point_size, marker='.',
            transform=datacrs, color = 'tab:orange')

    plot_figure_text(ax3, 'During Dixie Fire', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right', location = 'upper_right')
    plot_figure_text(ax3, 'NEXRAD-17 ' + \
        str(nexrad_channel_dict[str(channel2)]['wavelength']) + ' μm', \
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
        outname = save_dir + 'nexrad_time_series_channel_comp_2date_' + \
            NEXRAD_dict1['ptype'] + '_ch' + \
            str(channel1) + '_ch' + str(channel2)
        if(ch_idx3 is not None):
            outname = outname + '_ch' + str(channel3)
        outname = outname + '_pt'+ str(idx1) + '_' + date_str1 + \
            '_' + date_str2 + '.png'
            #NEXRAD_dict['dt_dates'][time_idx].strftime('%Y%m%d%H%M') + \
            #'.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

def plot_NEXRAD_time_series_mesh(NEXRAD_dict, ch_idx1 = 1, \
        ch_idx2 = 0, save_dir = './', \
        date_idx = 23, sigma = 0.8, map_cntr = True, \
        show_points = False, save = False):

    channel1 = int(NEXRAD_dict['channels'][ch_idx1])
    channel2 = int(NEXRAD_dict['channels'][ch_idx2])

    #dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
    dt_date_str = NEXRAD_dict['dt_dates'][date_idx]
    date_str = dt_date_str.strftime('%Y%m%d%H%M')
    var, crs, lons, lats, lat_lims, lon_lims, plabel = \
        read_NEXRAD_satpy(date_str, channel1)
    if(map_cntr):
        var2, crs2, lons2, lats2, lat_lims2, lon_lims2, plabel2 = \
            read_NEXRAD_satpy(date_str, channel2)

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
        levels = np.arange(0, nexrad_channel_dict[str(channel1)\
            ]['limits'][1], 2)

    mesh = ax1.pcolormesh(NEXRAD_dict['dt_dates'], \
        NEXRAD_dict['plat'][:],\
        NEXRAD_dict['data'][:,ch_idx1,:].T, cmap = cmap1, \
        vmin = nexrad_channel_dict[str(channel1)\
            ]['limits'][0], \
        vmax = nexrad_channel_dict[str(channel1)\
            ]['limits'][1], \
        shading = 'auto')
   
    ax1.set_ylabel('Latitude') 
    ax1.xaxis.set_major_formatter(DateFormatter('%m/%d\n%H:%MZ'))
    ax1.tick_params(axis="x", labelsize = 9)

    cbar = plt.colorbar(mesh, ax = ax1, pad = 0.03, fraction = 0.052, \
        extend = 'both', label = plabel)

    ##!## Add 'x's for the location of the smoke in the time series
    ##!## ---------------------------------------------------------
    ##!#tmp_data = np.full(NEXRAD_dict['data'][:,ch_idx2,:].T.shape, np.nan)
    ##!#in_smoke = (NEXRAD_dict['data'][:,ch_idx2,:].T > 22.)
    ##!#tmp_data[in_smoke] = 1.
    ##!#tmp_data[~in_smoke] = 0.

    ##!#mask_tmp_data = np.ma.masked_invalid(tmp_data)

    ##!#hash_data   = np.ma.masked_where(mask_tmp_data != 1, mask_tmp_data)
    #nohash_data = np.ma.masked_where(mask_tmp_data == 1, mask_tmp_data)

    ##!#hash0 = ax1.pcolor(NEXRAD_dict['dt_dates'], NEXRAD_dict['plat'][:],\
    ##!#    hash_data, hatch = 'xx', alpha=0., shading = 'auto')
    smooth_data = gaussian_filter(NEXRAD_dict['data'][:,ch_idx2,:].T, \
        sigma = sigma)
    cntr = ax1.contour(NEXRAD_dict['dt_dates'], \
        NEXRAD_dict['plat'][:],\
        smooth_data, cmap = cmap2, \
        vmin = nexrad_channel_dict[str(channel2)]['limits'][0], \
        vmax = nexrad_channel_dict[str(channel2)]['limits'][1], \
        levels = levels)
    ax1.clabel(cntr, cntr.levels, inline=True, fontsize=8)

    # Add vertical line at the image time
    # -----------------------------------
    ax1.axvline(dt_date_str, color = 'black',\
        linestyle = ':')

    #cbar.set_label(plabel.replace('_',' '), size = labelsize, weight = 'bold')

    # Plot the NEXRAD 0.64 micron image
    # -------------------------------
    labelsize = 10
    font_size = 8
    plot_NEXRAD_satpy(date_str, channel1, \
        ax = ax2, var = var, crs = crs, \
        lons = lons, lats = lats, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = nexrad_channel_dict[str(channel1)]['limits'][0], \
        vmax = nexrad_channel_dict[str(channel1)]['limits'][1], \
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

    plot_figure_text(ax2, 'NEXRAD-17 ' + \
        str(nexrad_channel_dict[str(channel1)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    ax2.set_title(dt_date_str.strftime(\
        '%Y-%m-%d %H:%MZ'), fontsize = 10)

    # Add line showing the cross section
    # ----------------------------------
    max_lat = np.max(NEXRAD_dict['plat'])
    max_lon = np.min(NEXRAD_dict['plon'])
    min_lat = np.min(NEXRAD_dict['plat'])
    min_lon = np.max(NEXRAD_dict['plon'])

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
        for tlat, tlon in zip(NEXRAD_dict['plat'][1::2], NEXRAD_dict['plon'][1::2]):
            print(tlat, tlon)
            ax1.plot(NEXRAD_dict['dt_dates'][date_idx], tlat, \
                    linewidth=2, markersize = 5, marker='^',
                    color = 'black')
            ax1.plot(NEXRAD_dict['dt_dates'][date_idx], tlat, \
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
    ##!#ln1 = ax1.plot(NEXRAD_dict['dt_dates'], NEXRAD_dict['data'][:,0,idx], \
    ##!#    label = str(nexrad_channel_dict[\
    ##!#    str(NEXRAD_dict['channels'][0])]['wavelength']) + \
    ##!#    ' μm')
    ##!#ln2 = ax1.plot(NEXRAD_dict['dt_dates'], NEXRAD_dict['data'][:,1,idx], \
    ##!#    label = str(nexrad_channel_dict[\
    ##!#    str(NEXRAD_dict['channels'][1])]['wavelength']) + \
    ##!#    ' μm', linestyle = '--', color = 'tab:blue')
    ##!#ax12 = ax1.twinx()
    ##!#ln3 = ax12.plot(NEXRAD_dict['dt_dates'], NEXRAD_dict['data'][:,2,idx], \
    ##!#    label = str(nexrad_channel_dict[\
    ##!#    str(NEXRAD_dict['channels'][2])]['wavelength']) + \
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

    ##!#ax2.plot(NEXRAD_dict['dt_dates'], NEXRAD_dict['data'][:,3,idx], \
    ##!#    #label = str(int(NEXRAD_dict['channels'][3])))
    ##!#    label = str(nexrad_channel_dict[\
    ##!#    str(NEXRAD_dict['channels'][3])]['wavelength']) + \
    ##!#    ' μm')
    ##!#ax2.plot(NEXRAD_dict['dt_dates'], NEXRAD_dict['data'][:,4,idx], \
    ##!#    #label = str(int(NEXRAD_dict['channels'][4])))
    ##!#    label = str(nexrad_channel_dict[\
    ##!#    str(NEXRAD_dict['channels'][4])]['wavelength']) + \
    ##!#    ' μm')
    ##!##ax22 = ax2.twinx()
    ##!#ax2.plot(NEXRAD_dict['dt_dates'], NEXRAD_dict['data'][:,5,idx], \
    ##!#    #label = str(int(NEXRAD_dict['channels'][5])))
    ##!#    label = str(nexrad_channel_dict[\
    ##!#    str(NEXRAD_dict['channels'][5])]['wavelength']) + \
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
    ##!#    outname = save_dir + 'nexrad_time_series_channels_pt' + \
    ##!#        str(idx) + '.png'
    ##!#        #NEXRAD_dict['dt_dates'][time_idx].strftime('%Y%m%d%H%M') + \
    ##!#        #'.png'
    ##!#    fig.savefig(outname, dpi = 300)
    ##!#    print("Saved image", outname)
    ##!#else:
    ##!#    plt.show()

def plot_NEXRAD_cross_channels(NEXRAD_dict, time_idx = 20, \
        ch_idx1 = 0, ch_idx2 = 1, save_dir = './', save = False):

    channel1 = NEXRAD_dict['channels'][ch_idx1]
    channel2 = NEXRAD_dict['channels'][ch_idx2]
    #dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
    dt_date_str = NEXRAD_dict['dt_dates'][time_idx]
    date_str = dt_date_str.strftime('%Y%m%d%H%M')
    var1, crs1, lons1, lats1, lat_lims1, lon_lims1, plabel1 = \
        read_NEXRAD_satpy(date_str, channel1)
    var2, crs2, lons2, lats2, lat_lims2, lon_lims2, plabel2 = \
        read_NEXRAD_satpy(date_str, channel2)

    plt.close('all')
    fig = plt.figure(figsize = (9, 6))
    gs = fig.add_gridspec(nrows = 2, ncols = 4)
    ax2  = fig.add_subplot(gs[0, 0:2], projection = crs1)   # true color    
    ax3  = fig.add_subplot(gs[0, 2:4], projection = crs2)   # true color    
    ax1  = fig.add_subplot(gs[1, 0:4]) # Ch 1


    ln1 = ax1.plot(NEXRAD_dict['plat'], \
            NEXRAD_dict['data'][time_idx,ch_idx1,:], \
            label = str(nexrad_channel_dict[\
            str(NEXRAD_dict['channels'][ch_idx1])]['wavelength']) + \
            ' μm')
    ax12 = ax1.twinx()
    ln2 = ax12.plot(NEXRAD_dict['plat'], \
        NEXRAD_dict['data'][time_idx,ch_idx2,:], \
        label = str(nexrad_channel_dict[\
        str(NEXRAD_dict['channels'][ch_idx2])]['wavelength']) + \
        ' μm', color = 'tab:orange')
    ##!#ln2 = ax1.plot(NEXRAD_dict['dt_dates'], NEXRAD_dict['data'][:,1,idx], \
    ##!#    label = str(nexrad_channel_dict[\
    ##!#    str(NEXRAD_dict['channels'][1])]['wavelength']) + \
    ##!#    ' μm', linestyle = '--', color = 'tab:blue')
    ax1.set_ylabel('Reflectance [%]', color = 'tab:blue')
    ax12.set_ylabel('Brightness Temperature [K]', color = 'tab:orange')
    ax1.grid()
    #ax1.xaxis.set_major_formatter(DateFormatter('%m/%d\n%H:%MZ'))
    ax1.tick_params(axis="x", labelsize = 9)
    ax1.tick_params(axis="y", color = 'tab:blue', labelcolor = 'tab:blue')
    ax12.tick_params(axis="y", color = 'tab:orange', labelcolor = 'tab:orange')

    # Plot the NEXRAD 0.64 micron image
    # -------------------------------
    labelsize = 10
    font_size = 8
    plot_NEXRAD_satpy(date_str, channel1, \
        ax = ax2, var = var1, crs = crs1, \
        lons = lons1, lats = lats1, lat_lims = lat_lims1, 
        lon_lims = lon_lims1, \
        vmin = nexrad_channel_dict[str(channel1)]['limits'][0], \
        vmax = nexrad_channel_dict[str(channel1)]['limits'][1], \
        ptitle = '', plabel = plabel1, \
        #vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)
    plot_figure_text(ax2, 'NEXRAD-17 ' + \
        str(nexrad_channel_dict[str(channel1)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    ax2.set_title(dt_date_str.strftime(\
        '%Y-%m-%d %H:%MZ'), fontsize = 10)
    # Add line showing the cross section
    # ----------------------------------
    max_lat = np.max(NEXRAD_dict['plat'])
    max_lon = np.min(NEXRAD_dict['plon'])
    min_lat = np.min(NEXRAD_dict['plat'])
    min_lon = np.max(NEXRAD_dict['plon'])

    ax2.plot([max_lon, min_lon], [max_lat, min_lat],
             color='black', linewidth=3, marker='o',
             transform=datacrs, markersize = 5, 
             )
    ax2.plot([max_lon, min_lon], [max_lat, min_lat],
             color='cyan', linewidth=1, marker='o',
             transform=datacrs, markersize = 3, 
             )

    plot_NEXRAD_satpy(date_str, channel2, \
        ax = ax3, var = var2, crs = crs2, \
        lons = lons2, lats = lats2, lat_lims = lat_lims2, \
        lon_lims = lon_lims2, \
        vmin = nexrad_channel_dict[str(channel2)]['limits'][0], \
        vmax = nexrad_channel_dict[str(channel2)]['limits'][1], \
        ptitle = '', plabel = plabel2, \
        #vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)
    plot_figure_text(ax3, 'NEXRAD-17 ' + \
        str(nexrad_channel_dict[str(channel2)]['wavelength']) + ' μm', \
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
        outname = save_dir + 'nexrad_cross_channels_' + \
            NEXRAD_dict['dt_dates'][time_idx].strftime('%Y%m%d%H%M') + \
            '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

def read_NEXRAD_time_series_auto(begin_date, end_date, \
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

    # Find all downloaded NEXRAD filenames that are between these
    # two dates
    # ---------------------------------------------------------
    all_files = np.array(glob('/home/bsorenson/data/NEXRAD/nexrad17_abi/*.nc'))
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
        print("     Run auto_NEXRAD_download to get the data")
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

    nexrad_data = np.full((len(good_unique_times), len(channels), \
        len(dlat)), \
        np.nan)
    nexrad_lats = np.full((len(good_unique_times), len(channels), \
        len(dlat)), \
        np.nan)
    nexrad_lons = np.full((len(good_unique_times), len(channels), \
        len(dlat)), \
        np.nan)

    for ii, ttime in enumerate(good_unique_times):
        print(ttime.strftime('%Y%m%d%H%M'))
        date_str = ttime.strftime('%Y%m%d%H%M')
        #if(date_str == '202107202126'):
        #    print("Not making image for this time")
        #else:
        for jj, tch in enumerate(channels):
            # Extract the NEXRAD values for the current time
            # --------------------------------------------
            nexrad_vals, nexrad_lats_local, nexrad_lons_local  = \
                get_NEXRAD_data_lat_lon(date_str, dlat, dlon, tch)
            nexrad_data[ii,jj,:] = nexrad_vals / 1.
            nexrad_lats[ii,jj,:] = nexrad_lats_local / 1.
            nexrad_lons[ii,jj,:] = nexrad_lons_local / 1.

    # Put the data in a dictionary
    # ----------------------------
    out_dict = {}
    out_dict['data'] = nexrad_data
    out_dict['nexrad_lats'] = nexrad_lats
    out_dict['nexrad_lons'] = nexrad_lons
    out_dict['dt_dates'] = good_unique_times
    out_dict['channels'] = channels
    out_dict['plat'] = dlat
    out_dict['plon'] = dlon

    return out_dict

def read_NEXRAD_time_series_NCDF(file_name):

    # Set up the output dictionary
    # ----------------------------
    NEXRAD_dict = {}

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
    NEXRAD_dict['data']      = nc['data'][:,:,:].data
    NEXRAD_dict['nexrad_lats'] = nc['nexrad_lat'][:,:,:].data
    NEXRAD_dict['nexrad_lons'] = nc['nexrad_lon'][:,:,:].data
    NEXRAD_dict['dt_dates']  = dates
    NEXRAD_dict['channels']  = nc['channel_num'][:].data
    NEXRAD_dict['plat']      = nc['point_lat'][:].data
    NEXRAD_dict['plon']      = nc['point_lon'][:].data
    NEXRAD_dict['ptype']     = nc.ptype

    nc.close()

    return NEXRAD_dict
 
def write_NEXRAD_time_series_NCDF(NEXRAD_dict, save_dir = './'):

    file_name_start = save_dir + 'nexrad_cross_data_' + NEXRAD_dict['ptype'] + \
        '_' + \
        NEXRAD_dict['dt_dates'][0].strftime('%Y%m%d%H%M') + '_' + \
        NEXRAD_dict['dt_dates'][-1].strftime('%Y%m%d%H%M')
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
    num_point = len(NEXRAD_dict['plat'])
    num_ch    = len(NEXRAD_dict['channels'])
    num_time  = len(NEXRAD_dict['dt_dates'])
  
    # For the dates, calculate the number of seconds between 
    # each date and a reference date, which is January 1st, 2000
    # at 00:00 UTC.
    # ---------------------------------------------------------- 
    base_date = datetime(year=2000,month=1,day=1,hour = 0, minute = 0,\
        second = 0) 
    times = np.array([(tdate - base_date).total_seconds() \
        for tdate in NEXRAD_dict['dt_dates']])
  
    # Add an attribute to the file to contain the cross section
    # location
    # ---------------------------------------------------------
    nc.ptype = NEXRAD_dict['ptype']
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
    CHANNEL.description = 'NEXRAD channel number'
    POINT = nc.createVariable('cross_point','i2',('n_point'))
    POINT.description = 'Cross section point index'
   
    # Create variables
    # ----------------
    DATA = nc.createVariable('data','f4',('n_time','n_ch','n_point'))
    DATA.description = 'NEXRAD reflectance and/or brightness temperature at each point'
    NEXRAD_LAT = nc.createVariable('nexrad_lat','f4',('n_time','n_ch','n_point'))
    NEXRAD_LAT.description = 'NEXRAD latitude of each point for each channel'
    NEXRAD_LON = nc.createVariable('nexrad_lon','f4',('n_time','n_ch','n_point'))
    NEXRAD_LON.description = 'NEXRAD longitude of each point for each channel'
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
    CHANNEL[:]       = NEXRAD_dict['channels'][:]
    POINT[:]         = np.arange(num_point)
    DATA[:,:,:]      = NEXRAD_dict['data'][:,:,:]
    NEXRAD_LAT[:,:,:]  = NEXRAD_dict['nexrad_lats'][:,:,:]
    NEXRAD_LON[:,:,:]  = NEXRAD_dict['nexrad_lons'][:,:,:]
    POINT_LAT[:]     = NEXRAD_dict['plat'][:]
    POINT_LON[:]     = NEXRAD_dict['plon'][:]
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

def plot_NEXRAD_satpy_point_test(date_str, channel = 2, \
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

    # Read NEXRAD image data
    # --------------------
    #print(NEXRAD_dict['dt_dates'][time_idx])
    #date_str = NEXRAD_dict['dt_dates'][time_idx].strftime('%Y%m%d%H%M')
    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
    var, crs, lons, lats, lat_lims, lon_lims, plabel = \
        read_NEXRAD_satpy(date_str, channel)

    plt.close('all')
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1, projection = crs)

    labelsize = 10
    font_size = 8
    plot_NEXRAD_satpy(date_str, channel, \
        ax = ax1, var = var, crs = crs, \
        lons = lons, lats = lats, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = nexrad_channel_dict[str(channel)]['limits'][0], \
        vmax = nexrad_channel_dict[str(channel)]['limits'][1], \
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
        
    plot_figure_text(ax1, 'NEXRAD-17 ' + \
        str(nexrad_channel_dict[str(channel)]['wavelength']) + ' μm', \
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

def plot_NEXRAD_satpy_hatch(date_str, save = False):

    # Read NEXRAD image data
    # --------------------
    #print(NEXRAD_dict['dt_dates'][time_idx])
    #date_str = NEXRAD_dict['dt_dates'][time_idx].strftime('%Y%m%d%H%M')
    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
    var1, crs1, lons1, lats1, lat_lims1, lon_lims1, plabel1 = \
        read_NEXRAD_satpy(date_str, 2)
    var2, crs2, lons2, lats2, lat_lims2, lon_lims2, plabel2 = \
        read_NEXRAD_satpy(date_str, 13)

    plt.close('all')
    fig = plt.figure()
    ax1 = fig.add_subplot(1,2,1, projection = crs1)
    ax2 = fig.add_subplot(1,2,2, projection = crs2)

    labelsize = 10
    font_size = 8
    plot_NEXRAD_satpy(date_str, 2, \
        ax = ax1, var = var1, crs = crs1, \
        lons = lons1, lats = lats1, lat_lims = lat_lims1, \
        lon_lims = lon_lims1, \
        vmin = nexrad_channel_dict[str(2)]['limits'][0], \
        vmax = nexrad_channel_dict[str(2)]['limits'][1], \
        ptitle = '', plabel = plabel1, \
        #vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,\
        save=False)
    plot_NEXRAD_satpy(date_str, 13, \
        ax = ax2, var = var2, crs = crs2, \
        lons = lons2, lats = lats2, lat_lims = lat_lims2, \
        lon_lims = lon_lims2, \
        vmin = nexrad_channel_dict[str(13)]['limits'][0], \
        vmax = nexrad_channel_dict[str(13)]['limits'][1], \
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

    plot_figure_text(ax1, 'NEXRAD-17 ' + \
        str(nexrad_channel_dict[str(2)]['wavelength']) + ' μm', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    ax1.set_title(dt_date_str.strftime(\
        '%Y-%m-%d %H:%MZ'), fontsize = 10)

    plot_figure_text(ax2, 'NEXRAD-17 ' + \
        str(nexrad_channel_dict[str(13)]['wavelength']) + ' μm', \
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
