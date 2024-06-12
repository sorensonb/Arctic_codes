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
from matplotlib.cm import ScalarMappable
from matplotlib.dates import DateFormatter
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as mpatches
from metpy.plots import USCOUNTIES
#from matplotlib.tri import Triangulation
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from satpy.scene import Scene
from satpy import find_files_and_readers
from satpy.writers import get_enhanced_image
from glob import glob
import os

home_dir = os.environ['HOME']

sys.path.append(home_dir)
from python_lib import plot_trend_line, plot_subplot_label, plot_figure_text, \
    nearest_gridpoint, aerosol_event_dict, init_proj, \
    convert_radiance_to_temp, format_coord, circle, listFD, \
    laads_daac_key, plot_point_on_map

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# Set up global variables
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
modis_dir = home_dir  + '/data/MODIS/Aqua/'
data_dir  = modis_dir + 'MYD/'
myd06_dir  = modis_dir + 'MYD06/'
myd08_dir  = modis_dir + 'MYD08/'
cloud_dir = modis_dir + 'CLDMSK/'
cloudL3_dir = home_dir + '/data/MODIS/combined/'
cloudL3_daily_dir = home_dir + '/data/MODIS/combined/daily/'
cloudL3_monthly_dir = home_dir + '/data/MODIS/combined/monthly/'
datacrs = ccrs.PlateCarree()
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

channel_dict = {
    '1': {
        'Name': 'EV_250_Aggr1km_RefSB',\
        'Index': 0,\
        'cmap': 'Greys_r',
        'Unit': '%',
        'Unit_name': 'reflectance',
        'Bandwidth': [0.620, 0.670]
    },\
    '2': {
        'Name': 'EV_250_Aggr1km_RefSB',\
        'Index': 1,\
        'cmap': 'Greys_r',
        'Unit': '%',
        'Unit_name': 'reflectance',
        'Bandwidth': [0.841, 0.876]
    },\
    '3': {
        'Name': 'EV_500_Aggr1km_RefSB',\
        'Index': 0,\
        'cmap': 'Greys_r',
        'Unit': '%',
        'Unit_name': 'reflectance',
        'Bandwidth': [0.459, 0.479]
    },\
    '4': {
        'Name': 'EV_500_Aggr1km_RefSB',\
        'Index': 1,\
        'cmap': 'Greys_r',
        'Unit': '%',
        'Unit_name': 'reflectance',
        'Bandwidth': [0.545, 0.565]
    },\
    '5': {
        'Name': 'EV_500_Aggr1km_RefSB',\
        'Index': 2,\
        'cmap': 'Greys_r',
        'Unit': '%',
        'Unit_name': 'reflectance',
        'Bandwidth': [1.230, 1.250]
    },\
    '6': {
        'Name': 'EV_500_Aggr1km_RefSB',\
        'Index': 3,\
        'cmap': 'Greys_r',
        'Unit': '%',
        'Unit_name': 'reflectance',
        'Bandwidth': [1.628, 1.652]
    },\
    '7': {
        'Name': 'EV_500_Aggr1km_RefSB',\
        'Index': 4,\
        'cmap': 'Greys_r',
        'Unit': '%',
        'Unit_name': 'reflectance',
        'Bandwidth': [2.105, 2.155]
    },\
    '8': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 0,\
        'cmap': 'Greys_r',
        'Unit': '%',
        'Unit_name': 'reflectance',
        'Bandwidth': [0.405, 0.420]
    },\
    '9': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 1,\
        'cmap': 'Greys_r',
        'Unit': '%',
        'Unit_name': 'reflectance',
        'Bandwidth': [0.438, 0.448]
    },\
    '10': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 2,\
        'cmap': 'Greys_r',
        'Unit': '%',
        'Unit_name': 'reflectance',
        'Bandwidth': [0.483, 0.493]
    },\
    '11': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 3,\
        'cmap': 'Greys_r',
        'Unit': '%',
        'Unit_name': 'reflectance',
        'Bandwidth': [0.526, 0.536]
    },\
    '12': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 4,\
        'cmap': 'Greys_r',
        'Unit': '%',
        'Unit_name': 'reflectance',
        'Bandwidth': [0.546, 0.556]
    },\
    '13lo': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 5,\
        'cmap': 'Greys_r',
        'Unit': '%',
        'Unit_name': 'reflectance',
        'Bandwidth': [0.662, 0.672]
    },\
    '13hi': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 6,\
        'cmap': 'Greys_r',
        'Unit': '%',
        'Unit_name': 'reflectance',
        'Bandwidth': [0.662, 0.672]
    },\
    '14lo': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 7,\
        'cmap': 'Greys_r',
        'Unit': '%',
        'Unit_name': 'reflectance',
        'Bandwidth': [0.673, 0.683]
    },\
    '14hi': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 8,\
        'cmap': 'Greys_r',
        'Unit': '%',
        'Unit_name': 'reflectance',
        'Bandwidth': [0.673, 0.683]
    },\
    '15': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 9,\
        'cmap': 'Greys_r',
        'Unit': '%',
        'Unit_name': 'reflectance',
        'Bandwidth': [0.743, 0.753]
    },\
    '16': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 10,\
        'cmap': 'Greys_r',
        'Unit': '%',
        'Unit_name': 'reflectance',
        'Bandwidth': [0.862, 0.877]
    },\
    '17': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 11,\
        'cmap': 'Greys_r',
        'Unit': '%',
        'Unit_name': 'reflectance',
        'Bandwidth': [0.890, 0.920]
    },\
    '18': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 12,\
        'cmap': 'Greys_r',
        'Unit': '%',
        'Unit_name': 'reflectance',
        'Bandwidth': [0.931, 0.941]
    },\
    '19': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 13,\
        'cmap': 'Greys_r',
        'Unit': '%',
        'Unit_name': 'reflectance',
        'Bandwidth': [0.915, 0.965]
    },\
    '20': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 0,\
        'cmap': 'plasma',
        'Unit': 'K',
        'Unit_name': 'brightness_temperature',
        'Bandwidth': [3.660, 3.840]
    },\
    '21': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 1,\
        'cmap': 'plasma',
        'Unit': 'K',
        'Unit_name': 'brightness_temperature',
        'Bandwidth': [3.929, 3.989]
    },\
    '22': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 2,\
        'cmap': 'plasma',
        'Unit': 'K',
        'Unit_name': 'brightness_temperature',
        'Bandwidth': [3.929, 3.989]
    },\
    '23': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 3,\
        'cmap': 'plasma',
        'Unit': 'K',
        'Unit_name': 'brightness_temperature',
        'Bandwidth': [4.020, 4.080]
    },\
    '24': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 4,\
        'cmap': 'plasma',
        'Unit': 'K',
        'Unit_name': 'brightness_temperature',
        'Bandwidth': [4.433, 4.498]
    },\
    '25': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 5,\
        'cmap': 'plasma',
        'Unit': 'K',
        'Unit_name': 'brightness_temperature',
        'Bandwidth': [4.482, 4.549]
    },\
    '26': {
        'Name': 'EV_1KM_RefSB',\
        'Index': 14,\
        'cmap': 'plasma',
        'Unit': 'K',
        'Unit_name': 'brightness_temperature',
        'Bandwidth': [1.360, 1.390]
    },\
    '27': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 6,\
        'cmap': 'plasma',
        'Unit': 'K',
        'Unit_name': 'brightness_temperature',
        'Bandwidth': [6.535, 6.895]
    },\
    '28': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 7,\
        'cmap': 'plasma',
        'Unit': 'K',
        'Unit_name': 'brightness_temperature',
        'Bandwidth': [7.175, 7.475]
    },\
    '29': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 8,\
        'cmap': 'plasma',
        'Unit': 'K',
        'Unit_name': 'brightness_temperature',
        'Bandwidth': [8.400, 8.700]
    },\
    '30': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 9,\
        'cmap': 'plasma',
        'Unit': 'K',
        'Unit_name': 'brightness_temperature',
        'Bandwidth': [9.580, 9.880]
    },\
    '31': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 10,\
        'cmap': 'plasma',
        'Unit': 'K',
        'Unit_name': 'brightness_temperature',
        'Bandwidth': [10.780, 11.280]
    },\
    '32': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 11,\
        'cmap': 'plasma',
        'Unit': 'K',
        'Unit_name': 'brightness_temperature',
        'Bandwidth': [11.770, 12.270]
    },\
    '33': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 12,\
        'cmap': 'plasma',
        'Unit': 'K',
        'Unit_name': 'brightness_temperature',
        'Bandwidth': [13.185, 13.485]
    },\
    '34': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 13,\
        'cmap': 'plasma',
        'Unit': 'K',
        'Unit_name': 'brightness_temperature',
        'Bandwidth': [13.485, 13.785]
    },\
    '35': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 14,\
        'cmap': 'plasma',
        'Unit': 'K',
        'Unit_name': 'brightness_temperature',
        'Bandwidth': [13.785, 14.085]
    },\
    '36': {
        'Name': 'EV_1KM_Emissive',\
        'Index': 15,\
        'cmap': 'plasma',
        'Unit': 'K',
        'Unit_name': 'brightness_temperature',
        'Bandwidth': [14.085, 14.385]
    },\
    'wv_ir': {
        'Name': 'Water_Vapor_Infrared',\
        'Index': None,
        'cmap': 'plasma',
        'Unit': 'K',
        'Unit_name': 'brightness_temperature',
        'Bandwidth': None 
    },\
    'wv_nir': {
        'Name': 'Water_Vapor_Near_Infrared',\
        'Index': None,
        'cmap': 'plasma',
        'Unit': 'K',
        'Unit_name': 'brightness_temperature',
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

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Downloading  and writing functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Using a CERES swath as a base, download the matching MODIS granules
# -------------------------------------------------------------------
def download_MODIS_swath(CERES_date_str, \
        dest_dir = modis_dir, download = True, \
        download_cloud_mask = True, \
        download_myd06 = True):

    if(home_dir + '/Research/CERES/' not in sys.path):
        sys.path.append(home_dir + '/Research/CERES/')
    from gridCERESLib import readgridCERES_hrly_grid
    
    dt_date_str = datetime.strptime(CERES_date_str, '%Y%m%d%H')

    # Read in the CERES swath
    # -----------------------
    CERES_data_hrly = readgridCERES_hrly_grid(CERES_date_str, 'SWF', \
        satellite = 'Aqua', minlat = 65.0, season='all')
    
    # Use the pixel times to determine the granules to grab
    # -----------------------------------------------------
    min_dt = np.min(CERES_data_hrly['time_dt'])   
    max_dt = np.max(CERES_data_hrly['time_dt'])   

    min_modis = min_dt - timedelta(minutes = min_dt.minute % 5, \
                                   seconds = min_dt.second, \
                                   microseconds = min_dt.microsecond)

    max_modis = max_dt - timedelta(minutes = max_dt.minute % 5, \
                                   seconds = max_dt.second, \
                                   microseconds = max_dt.microsecond)

    local_data_dir  = modis_dir + 'MYD/'
    local_cloud_dir = modis_dir + 'CLDMSK/'

    modis_time = min_modis
    return_str = ''
    modis_file_list = ''
    cldmk_file_list = ''
    myd06_file_list = ''
    download_return_dict = {}
    while(modis_time <= max_modis):

        local_time = modis_time.strftime('%Y%m%d%H%M')
        print(local_time) 

        # First, check for the MYD and CLDMSK files
        # -----------------------------------------
        modis_file = glob(modis_time.strftime(local_data_dir + \
            '*%Y%j.%H%M.*.hdf'))

        cldmk_file = glob(modis_time.strftime(local_cloud_dir + \
            '*A%Y%j.%H%M.*.nc'))

        myd06_file = glob(modis_time.strftime(myd06_dir + \
            '*A%Y%j.%H%M.*.nc'))

        if((len(modis_file) == 0) | \
                ((len(cldmk_file) == 0) & (download_cloud_mask)) | \
                ((len(myd06_file) == 0) & (download_myd06))):
            print("MYD or CLDMSK or MYD06 file not found. Must download")

            if(download):
                #modis_file = download_MODIS_file(local_time, \
                return_dict = download_MODIS_file(local_time, \
                    dest_dir = dest_dir, \
                    download_cloud_mask = download_cloud_mask, \
                    download_myd06 = download_myd06)
                modis_file = return_dict['modis_file']
                cldmk_file = return_dict['cldmk_file']
                myd06_file = return_dict['myd06_file']

                # Check for download errors
                # -------------------------
                if(modis_file == -1):
                    print("No MODIS file returned. Continuing.")
                    continue
        else:
            print("MODIS data downloaded",modis_file[0], cldmk_file[0], \
                myd06_file[0])
            modis_file = modis_file[0]
            cldmk_file = cldmk_file[0]
            myd06_file = cldmk_file[0]
        
        ##!#try:
        ##!#    modis_file = subprocess.check_output('ls ' + \
        ##!#        modis_time.strftime(local_data_dir + '*%Y%j.%H%M.*.hdf'), \
        ##!#        shell = True).decode('utf-8').strip().split('\n')[0]
        ##!#    print(modis_file)
        ##!#except subprocess.CalledProcessError:
        ##!#    print("MYD file not found. Must download")

        ##!#    if(download):
        ##!#        found_file = download_MODIS_file(local_time, \
        ##!#            dest_dir = dest_dir, \
        ##!#            download_cloud_mask = download_cloud_mask)

        ##!#        # Check for download errors
        ##!#        # -------------------------
        ##!#        if(found_file == -1):
        ##!#            print("No MODIS file returned. Continuing.")
        ##!#            continue

        #if(os.path.exists(dest_dir + found_file)):
        #    print(found_file + ' already exists. Not downloading')

        return_str = return_str + local_time + ' '
        modis_file_list = modis_file_list + modis_file + ' '
        cldmk_file_list = cldmk_file_list + cldmk_file + ' '
        myd06_file_list = myd06_file_list + myd06_file + ' '

        modis_time = modis_time + timedelta(minutes = 5)
   
    download_return_dict['modis_date_list'] = return_str.split()
    download_return_dict['modis_file_list'] = modis_file_list.split()
    download_return_dict['cldmk_file_list'] = cldmk_file_list.split()
    download_return_dict['myd06_file_list'] = myd06_file_list.split()

    return download_return_dict
    #return return_str.split(), modis_file_list.split(), cldmk_file_list

def identify_MODIS_MYD(date_str, modis_dir):

    local_data_dir  = modis_dir + 'MYD/'

    base_url = 'https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/61/MYD021KM'

    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')

    # For each desired channel, figure out the closest file time
    # to the input date
    # ----------------------------------------------------------
    try:
        files = listFD(dt_date_str.strftime(base_url + '/%Y/%j/'), ext = '.hdf')
    except subprocess.CalledProcessError:
        print("ERROR: No MODIS files for the input DTG",date_str)
        return -2

    if(len(files) == 0):
        print("ERROR: No MODIS files returned from the request. Exiting")
        return -1
    
    # Remove the timestamps from the file strings
    # -------------------------------------------
    files_only = [tfile.strip().split('/')[-1] for tfile in files]

    # Use the actual file times to get timestamps
    # -------------------------------------------
    file_dates = [datetime.strptime(tfile[10:22],'%Y%j.%H%M') for tfile in files_only]

    # Figure out where the closest files to the desired time are
    # ----------------------------------------------------------
    time_diffs = np.array([abs((dt_date_str - ddate).total_seconds()) \
        for ddate in file_dates])

    # Extract the index of the matching MODIS file
    # --------------------------------------------
    file_idx = np.argmin(time_diffs)
    found_file = files_only[file_idx]

    # Check if the file is already downloaded
    # ---------------------------------------
    if(os.path.exists(local_data_dir + found_file)):
        print(found_file + ' already exists. Not downloading')
    else: 
        # Download the file
        # -----------------
        cmnd = dt_date_str.strftime("wget \"" + base_url + '/%Y/%j/' + found_file + \
            "\" --header \"Authorization: Bearer " + laads_daac_key + "\" -P .")
        print(cmnd)
        os.system(cmnd)

        # Move the file to the destination folder
        # ---------------------------------------
        cmnd = "mv " + found_file + " " + local_data_dir
        print(cmnd) 
        os.system(cmnd)


    return found_file 

def identify_MODIS_CLDMSK(date_str, modis_dir):

    local_data_dir  = modis_dir + 'CLDMSK/'

    base_url = 'https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/' + \
        '5110/CLDMSK_L2_MODIS_Aqua'

    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')

    # For each desired channel, figure out the closest file time
    # to the input date
    # ----------------------------------------------------------
    try:
        files = listFD(dt_date_str.strftime(base_url + '/%Y/%j/'), ext = '.nc')
    except subprocess.CalledProcessError:
        print("ERROR: No MODIS CLDMSK files for the input DTG",date_str)
        return -2

    if(len(files) == 0):
        print("ERROR: No MODIS CLDMSK files returned from the request. Exiting")
        return -1
    
    # Remove the timestamps from the file strings
    # -------------------------------------------
    files_only = [tfile.strip().split('/')[-1] for tfile in files]

    # Use the actual file times to get timestamps
    # -------------------------------------------
    file_dates = [datetime.strptime(tfile[22:34],'%Y%j.%H%M') for tfile in files_only]

    # Figure out where the closest files to the desired time are
    # ----------------------------------------------------------
    time_diffs = np.array([abs((dt_date_str - ddate).total_seconds()) \
        for ddate in file_dates])

    # Extract the index of the matching MODIS file
    # --------------------------------------------
    file_idx = np.argmin(time_diffs)
    found_file = files_only[file_idx]

    # Check if the file is already downloaded
    # ---------------------------------------
    if(os.path.exists(local_data_dir + found_file)):
        print(found_file + ' already exists. Not downloading')
    else: 
        # Download the file
        # -----------------
        cmnd = dt_date_str.strftime("wget \"" + base_url + '/%Y/%j/' + found_file + \
            "\" --header \"Authorization: Bearer " + laads_daac_key + "\" -P .")
        print(cmnd)
        os.system(cmnd)

        # Move the file to the destination folder
        # ---------------------------------------
        cmnd = "mv " + found_file + " " + local_data_dir
        print(cmnd) 
        os.system(cmnd)

    return found_file 

def identify_MODIS_CLDL3(date_str, dest_dir = cloudL3_daily_dir):

    local_data_dir  = dest_dir

    base_url = 'https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/62/'+\
        'MCD06COSP_D3_MODIS/'

    dt_date_str = datetime.strptime(date_str, '%Y%m%d')

    print(dt_date_str)

    # For each desired channel, figure out the closest file time
    # to the input date
    # ----------------------------------------------------------
    try:
        files = listFD(dt_date_str.strftime(base_url + '/%Y/%j'), ext = '.nc')
    except subprocess.CalledProcessError:
        print("ERROR: No MODIS CLDL3 files for the input DTG",date_str)
        return -2

    if(len(files) == 0):
        print("ERROR: No MODIS CLDMSK files returned from the request. Exiting")
        return -1
    
    # Remove the timestamps from the file strings
    # -------------------------------------------
    found_file = files[0].strip().split('/')[-1]

    # Check if the file is already downloaded
    # ---------------------------------------
    if(os.path.exists(local_data_dir + found_file)):
        print(found_file + ' already exists. Not downloading')
    else: 
        # Download the file
        # -----------------
        cmnd = dt_date_str.strftime("wget \"" + base_url + '/%Y/%j/' + found_file + \
            "\" --header \"Authorization: Bearer " + laads_daac_key + "\" -P .")
        print(cmnd)
        os.system(cmnd)

        # Move the file to the destination folder
        # ---------------------------------------
        cmnd = "mv " + found_file + " " + local_data_dir
        print(cmnd) 
        os.system(cmnd)

    return found_file 


# date_str: YYYYMMDDHHMM
def identify_MODIS_MYD06(date_str, dest_dir = myd06_dir):

    local_data_dir  = dest_dir

    base_url = 'https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/61/'+\
        'MYD06_L2/'

    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')

    print(dt_date_str)

    # For each desired channel, figure out the closest file time
    # to the input date
    # ----------------------------------------------------------
    try:
        files = listFD(dt_date_str.strftime(base_url + '/%Y/%j'), ext = '.hdf')
    except subprocess.CalledProcessError:
        print("ERROR: No MODIS MYD06 files for the input DTG",date_str)
        return -2

    if(len(files) == 0):
        print("ERROR: No MODIS MYD06 files returned from the request. Exiting")
        return -1
    
    # Remove the timestamps from the file strings
    # -------------------------------------------
    files_only = [tfile.strip().split('/')[-1] for tfile in files]

    # Use the actual file times to get timestamps
    # -------------------------------------------
    file_dates = [datetime.strptime(tfile[10:22],'%Y%j.%H%M') for tfile in files_only]

    # Figure out where the closest files to the desired time are
    # ----------------------------------------------------------
    time_diffs = np.array([abs((dt_date_str - ddate).total_seconds()) \
        for ddate in file_dates])

    # Extract the index of the matching MODIS file
    # --------------------------------------------
    file_idx = np.argmin(time_diffs)
    found_file = files_only[file_idx]

    ## Remove the timestamps from the file strings
    ## -------------------------------------------
    #found_file = files[0].strip().split('/')[-1]

    # Check if the file is already downloaded
    # ---------------------------------------
    if(os.path.exists(local_data_dir + found_file)):
        print(found_file + ' already exists. Not downloading')
    else: 
        # Download the file
        # -----------------
        cmnd = dt_date_str.strftime("wget \"" + base_url + '/%Y/%j/' + found_file + \
            "\" --header \"Authorization: Bearer " + laads_daac_key + "\" -P .")
        print(cmnd)
        os.system(cmnd)

        # Move the file to the destination folder
        # ---------------------------------------
        cmnd = "mv " + found_file + " " + local_data_dir
        print(cmnd) 
        os.system(cmnd)

    return found_file 



# date_str: %Y%m or YYYYMMDD
def identify_MODIS_MYD08(date_str, dest_dir = myd08_dir):

    if(len(date_str) == 6):
        dt_date_str = datetime.strptime(date_str, '%Y%m')
        dtype = 'monthly'
    elif(len(date_str) == 8):
        dt_date_str = datetime.strptime(date_str, '%Y%m%d')
        dtype = 'daily'
    else:
        print("ERROR: Invalid date format. Must be YYYYMM or YYYYMMDD")
        return -3 

    local_data_dir  = dest_dir + dtype + '/'

    base_url = 'https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/61/'+\
        'MYD08_' + dtype[0].upper() + '3/'


    print(dt_date_str)

    # For each desired channel, figure out the closest file time
    # to the input date
    # ----------------------------------------------------------
    try:
        files = listFD(dt_date_str.strftime(base_url + '/%Y/%j'), ext = '.hdf')
    except subprocess.CalledProcessError:
        print("ERROR: No MODIS MYD08 files for the input DTG",date_str)
        return -2

    if(len(files) == 0):
        print("ERROR: No MODIS MYD08 files returned from the request. Exiting")
        return -1
    
    # Remove the timestamps from the file strings
    # -------------------------------------------
    found_file = files[0].strip().split('/')[-1]

    # Check if the file is already downloaded
    # ---------------------------------------
    if(os.path.exists(local_data_dir + found_file)):
        print(found_file + ' already exists. Not downloading')
    else: 
        # Download the file
        # -----------------
        cmnd = dt_date_str.strftime("wget \"" + base_url + '/%Y/%j/' + found_file + \
            "\" --header \"Authorization: Bearer " + laads_daac_key + "\" -P .")
        print(cmnd)
        os.system(cmnd)

        # Move the file to the destination folder
        # ---------------------------------------
        cmnd = "mv " + found_file + " " + local_data_dir
        print(cmnd) 
        os.system(cmnd)

    return found_file 


# This downloads the MODIS l1b HDF5 file that is closest to the passed
# date string from the LAADS DAAC archive. 
# --------------------------------------------------------------------
def download_MODIS_file(date_str, dest_dir = modis_dir, \
        download_cloud_mask = True, \
        download_myd06 = True):

    local_cloud_dir = modis_dir + 'CLDMSK/'

    found_file = identify_MODIS_MYD(date_str, modis_dir)
 
    if(download_cloud_mask):
        cloud_file = identify_MODIS_CLDMSK(date_str, modis_dir)
    else:
        cloud_file = ''

    if(download_myd06):
        myd06_file = identify_MODIS_MYD06(date_str, myd06_dir)
    else:
        myd06_file = ''   
 
    return_dict = {}
    return_dict['modis_file'] = found_file
    return_dict['cldmk_file'] = cloud_file
    return_dict['myd06_file'] = myd06_file
    
    return return_dict
    #return found_file

# date_str: YYYYMM
def download_MODIS_CLDL3_monthly(date_str, dest_dir = cloudL3_daily_dir):

    dt_date_str = datetime.strptime(date_str, '%Y%m')

    local_date_str = dt_date_str
    while(local_date_str.month == dt_date_str.month):

        found_file = identify_MODIS_CLDL3(local_date_str.strftime('%Y%m%d'), \
            dest_dir)

        print(found_file)

        local_date_str  = local_date_str + timedelta(days = 1)

# Writes a MODIS channel dictionary to HDF5 for Fortran colocation
# NOTE: can take either a date string OR a dictionary
# ----------------------------------------------------------------
def write_MODIS_to_HDF5(MODIS_data, channel = 1, swath = True, \
        save_path = './', minlat = 20., remove_empty_scans = False, \
        include_cloud_mask = True, include_myd06 = True, \
        remove_large_files = False):

    if(isinstance(MODIS_data, str)):
        MODIS_data = read_MODIS_channel(MODIS_data, channel, swath = swath, \
            include_cloud_mask = include_cloud_mask, \
            include_myd06 = include_myd06)

    if(remove_empty_scans):
        mask_data = np.ma.masked_where(MODIS_data['lat'] < minlat, \
            MODIS_data['data'])
        mask_dims = np.array([ (False in mask_data[ii,:].mask) for ii in \
            range(mask_data.shape[0])])
        keep_idxs = np.where(mask_dims == True)[0]
    else:
        keep_idxs = None 

    # Convert the filename object to datetime
    # ---------------------------------------
    file_date = MODIS_data['date']
    dt_date_str = datetime.strptime(file_date, '%Y%m%d%H%M')

    # Create a new netCDF dataset to write to the file
    # ------------------------------------------------
    outfile = save_path + 'modis_ch' + str(MODIS_data['channel']) + \
        '_subset_'+ file_date + '.hdf5'
    dset = h5py.File(outfile,'w')
 
    dset.create_dataset('latitude',  data = \
        MODIS_data['lat'][keep_idxs,:].squeeze())
    dset.create_dataset('longitude', data = \
        MODIS_data['lon'][keep_idxs,:].squeeze())
    dset.create_dataset('data', data = \
        MODIS_data['data'][keep_idxs,:].squeeze())
    if(include_cloud_mask):
        dset.create_dataset('cld', data = \
            MODIS_data['cloud_mask'][keep_idxs,:].squeeze())
    if(include_myd06):
        dset.create_dataset('cod', data = \
            MODIS_data['cloud_optical_depth'][keep_idxs,:].squeeze())
        dset.create_dataset('cld_top_pres', data = \
            MODIS_data['cloud_top_pressure'][keep_idxs,:].squeeze())

    # Save, write, and close the HDF5 file
    # --------------------------------------
    dset.close()

    print("Saved file ",outfile)  

    if(remove_large_files):
       print("REMOVE MODIS CLOUD FILE HERE") 

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Miscellaneous functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
 
# Plot the downloaded ASOS stations for each case
#
# NOTE: if you want to plot ASOS stations from a specific file, and not
# the default file associated with a YYYYMMDDHHMM case, pass in the
# ASOS file name in the 'date_str' argument
# ---------------------------------------------------------------------
def plot_ASOS_locs(pax,date_str,crs = datacrs, color='red', \
        sites = None, no_text = False):
    if(len(date_str.split('.')) > 1):
        asos_file = date_str
    else:
        dt_date_str = datetime.strptime(date_str,'%Y%m%d%H%M')
        cross_date = dt_date_str.strftime('%Y-%m-%d')
        file_date  = dt_date_str.strftime('%H%M')
        asos_file = aerosol_event_dict[cross_date][file_date]['asos']

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
            
            if(not no_text):
                pax.text(lon_stn + 0.1, lat_stn + 0.1, station, fontsize=10, \
                    weight='bold', transform=crs, color=color, backgroundcolor = 'white')

# Determine areas of an image that are in smoke, defined by:
#    (ch1_refl - ch5_refl) < mean(ch1_refl - ch5_refl) * 0.25
# NEW VERSION: 
#
#    cbar_labels = ['Clear','Smoke','Cloud']
#
def find_plume(dt_date_str, new_method = False):
    if(new_method):
        MODIS_ch1  = read_MODIS_channel(dt_date_str, 1,  zoom = True)
        MODIS_ch2  = read_MODIS_channel(dt_date_str, 2,  zoom = True)
        MODIS_ch5  = read_MODIS_channel(dt_date_str, 5,  zoom = True)
        MODIS_ch7  = read_MODIS_channel(dt_date_str, 7,  zoom = True)
        MODIS_ch32 = read_MODIS_channel(dt_date_str, 32, zoom = True)

        tmp_data = np.full(MODIS_ch1['data'].shape, np.nan)
        in_cloud = (((MODIS_ch1['data'] + MODIS_ch2['data'] > 1.2) |\
                     (MODIS_ch32['data'] < 265.)) |\
                    ((MODIS_ch1['data'] + MODIS_ch2['data'] > 0.7) &\
                     (MODIS_ch32['data'] < 285.)))
        tmp_data[in_cloud] = 2.
        tmp_data[~in_cloud] = 0.
        tmp_data[~in_cloud & (MODIS_ch1['data'] - \
            MODIS_ch7['data'] > 0.05) & (MODIS_ch7['data'] > 0.05)] = 1
        mask_tmp_data = np.ma.masked_invalid(tmp_data)

        # For the hash data (smoke), remove any pixels that do not have
        # tmp data equal to 1
        # -------------------------------------------------------------
        hash_data   = np.ma.masked_where(mask_tmp_data != 1, mask_tmp_data)
        nohash_data = np.ma.masked_where(mask_tmp_data == 1, mask_tmp_data)

    else:
        MODIS_ch1  = read_MODIS_channel(dt_date_str, 1,  zoom = True)
        MODIS_ch5  = read_MODIS_channel(dt_date_str, 5,  zoom = True)
        MODIS_ch31 = read_MODIS_channel(dt_date_str, 31, zoom = True)

        screen_limit = 0.20
        max_ch = 350.
        test_data = MODIS_ch1['data'] - MODIS_ch5['data']
        hash_data   = np.ma.masked_where(test_data <  \
            (np.nanmean(test_data) * screen_limit), MODIS_ch1['data'])
        hash_data = np.ma.masked_where(MODIS_ch5['data'] < 0.05, hash_data)
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
    cmnd = "ls " + home_dir + "/data/MODIS/"+variable+"/A*.nc"
    
    
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

def plot_compare_MODIS_COD(date_str, swath = True, save = False):

    MODIS_ch1  = read_MODIS_channel(date_str, 1,  \
        zoom = False, swath = swath, include_cloud_mask = True, \
        include_myd06 = True)
    MODIS_ch7  = read_MODIS_channel(date_str, 7,  \
        zoom = False, swath = swath, include_cloud_mask = True, \
        include_myd06 = True)
    
    fig = plt.figure(figsize = (10,6))
    ax1 = fig.add_subplot(2,3,1, projection = ccrs.NorthPolarStereo())
    ax2 = fig.add_subplot(2,3,2, projection = ccrs.NorthPolarStereo())
    ax3 = fig.add_subplot(2,3,3, projection = ccrs.NorthPolarStereo())
    ax4 = fig.add_subplot(2,3,4, projection = ccrs.NorthPolarStereo())
    ax5 = fig.add_subplot(2,3,5, projection = ccrs.NorthPolarStereo())
    
    mesh = ax1.pcolormesh(MODIS_ch1['lon'], MODIS_ch1['lat'], MODIS_ch1['data'], \
        transform = datacrs, shading = 'auto', cmap = MODIS_ch1['colors'], vmax = 0.7)
    cbar = plt.colorbar(mesh, ax = ax1, label = 'Reflectance')
    mesh = ax2.pcolormesh(MODIS_ch1['lon'], MODIS_ch1['lat'], MODIS_ch1['cloud_mask'], \
        transform = datacrs, \
        shading = 'auto', cmap = 'jet')
    cbar = plt.colorbar(mesh, ax = ax2, label = 'Cloud Mask')
    #ax2.pcolormesh(cldmsk_lon, cldmsk_lat, cldmsk_mask, transform = datacrs, \
    #    shading = 'auto', cmap = 'jet')
    mesh = ax3.pcolormesh(MODIS_ch1['lon'], MODIS_ch1['lat'], \
        MODIS_ch1['cloud_optical_depth'], transform = datacrs, \
        shading = 'auto', vmax = 40)
    cbar = plt.colorbar(mesh, ax = ax3, label = 'Cloud Optical Depth')
    mesh = ax4.pcolormesh(MODIS_ch7['lon'], MODIS_ch7['lat'], MODIS_ch7['data'], \
        transform = datacrs, vmax = 0.4, shading = 'auto', cmap = MODIS_ch1['colors'])
    cbar = plt.colorbar(mesh, ax = ax4, label = 'Reflectance')
    mesh = ax5.pcolormesh(MODIS_ch7['lon'], MODIS_ch7['lat'], \
        MODIS_ch1['cloud_top_pressure'], \
        transform = datacrs, shading = 'auto')
    cbar = plt.colorbar(mesh, ax = ax5, label = 'Cloud Top Pressure')
    #ax4.pcolormesh(myd06_lon, myd06_lat, myd06_dat, transform = datacrs, \
    #    shading = 'auto')
    
    ax1.coastlines()
    ax2.coastlines()
    ax3.coastlines()
    ax4.coastlines()
    ax5.coastlines()
    
    ax1.set_title('0.64 μm')
    ax2.set_title('Cloud Mask\nBlue = Cloud, Red = Clear')
    ax3.set_title('MODIS COD')
    ax4.set_title('MODIS 2.1 μm Refl.')
    ax5.set_title('MODIS CLD TOP PRES')
    
    ax1.set_boundary(circle, transform=ax1.transAxes)
    ax2.set_boundary(circle, transform=ax2.transAxes)
    ax3.set_boundary(circle, transform=ax3.transAxes)
    ax4.set_boundary(circle, transform=ax4.transAxes)
    ax5.set_boundary(circle, transform=ax5.transAxes)
    
    ax1.set_extent([-180,180,65,90], datacrs)
    ax2.set_extent([-180,180,65,90], datacrs)
    ax3.set_extent([-180,180,65,90], datacrs)
    ax4.set_extent([-180,180,65,90], datacrs)
    ax5.set_extent([-180,180,65,90], datacrs)
    
    plt.suptitle(date_str)
    
    fig.tight_layout()
   
    if(save): 
        outname = 'modis_cod_compare_' + dt_date_str + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()
        


def plot_compare_MODIS_cloud(date_str, swath = True, save = False):

    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')

    plt.close('all')
    mapcrs = ccrs.NorthPolarStereo()
    fig = plt.figure(figsize = (9, 9))
    ax1 = fig.add_subplot(2,2,1, projection = mapcrs)
    ax2 = fig.add_subplot(2,2,2, projection = mapcrs)
    ax3 = fig.add_subplot(2,2,3, projection = mapcrs)
    ax4 = fig.add_subplot(2,2,4, projection = mapcrs)
    
    zoom = False
    plot_MODIS_channel(date_str, 'true_color', swath = swath, \
        zoom = zoom, ax = ax1)
    plot_MODIS_channel(date_str, 1, swath = swath, \
        zoom = zoom, ax = ax2, vmax = 0.7)
    plot_MODIS_channel(date_str, 7, swath = swath, \
        zoom = zoom, ax = ax3, vmax = 0.4)
    plot_MODIS_channel(date_str, 'cloud_mask', swath = swath, \
        zoom = zoom, ax = ax4, vmax = None)

    ax1.set_extent([-180, 180, 65, 90], ccrs.PlateCarree())
    ax2.set_extent([-180, 180, 65, 90], ccrs.PlateCarree())
    ax3.set_extent([-180, 180, 65, 90], ccrs.PlateCarree())
    ax4.set_extent([-180, 180, 65, 90], ccrs.PlateCarree())
    ax1.coastlines()
    ax2.coastlines()
    ax3.coastlines()
    ax4.coastlines()

    plt.suptitle(dt_date_str.strftime('Aqua MODIS %Y-%m-%d %H:%M'))

    fig.tight_layout()

    if(save):
        outname = 'modis_comp_cloud_' + date_str + '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image",outname)
    else:
        plt.show()


def read_MODIS_satpy(date_str, channel,  composite = False, swath = False, \
        zoom = True, return_xy = False):
#def read_true_color(date_str,composite = False):

    if(swath):

        # Determine the correct GOES files associated with the date
        # ---------------------------------------------------------
        dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
        #dt_date_str_beg = dt_date_str - timedelta(minutes = 10) datetime.strptime(date_str,"%Y%m%d%H%M")
        dt_date_str_beg = datetime.strptime(date_str,"%Y%m%d%H%M")
        dt_date_str_end = dt_date_str + timedelta(minutes = 10)

        cmpst_add = '_composite'
        #lat_lims = [60, 90]
        #lon_lims = [-180, 180]

        try:
            # Use the Satpy find_files_and_readers to grab the files
            # ------------------------------------------------------
            files = find_files_and_readers(start_time = dt_date_str_beg, \
                end_time = dt_date_str_end, base_dir = data_dir, \
                reader = 'modis_l1b')
        except ValueError:
            print("ERROR: no files found for dtg",date_str)
            return

        day_filenames = files
        lat_lims = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis_Lat']
        lon_lims = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis_Lon']

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

    # Use satpy (Scene) to open the file
    # ----------------------------------
    scn = Scene(reader = 'modis_l1b', filenames = day_filenames)

    # Load true-color data
    scn.load(['true_color'])
    if(channel != 'true_color'): 
        scn.load([str(channel)], calibration = [channel_dict[str(channel)]['Unit_name']])
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
        plabel = channel_dict[str(channel)]['Unit_name'].title() + ' [' + \
            channel_dict[str(channel)]['Unit'] + ']'

    if(return_xy):

        y = new_scn[channel].y 
        x = new_scn[channel].x 
        del scn 
        del new_scn
        return var, crs, lons, lats, lat_lims, lon_lims, plabel, x, y
    else:
        del scn 
        del new_scn
        return var, crs, lons, lats, lat_lims, lon_lims, plabel
    #return var, crs, lat_lims, lon_lims

# channel must be an integer between 1 and 16
def plot_MODIS_satpy(date_str, channel, ax = None, var = None, crs = None, \
        lons = None, lats = None, lat_lims = None, lon_lims = None, \
        vmin = None, vmax = None, ptitle = None, plabel = None, \
        use_xy = False, 
        labelsize = 10, colorbar = True, swath = False, zoom=True,save=False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    if(var is None): 
        var, crs, lons, lats, lat_lims, lon_lims, plabel = read_MODIS_satpy(\
            date_str, str(channel), swath = swath)

    if(var is None): 
        if(use_xy):
            var, crs, lons, lats, lat_lims, lon_lims, plabel, xx, yy = \
                read_MODIS_satpy(date_str, str(channel), \
                swath = swath, return_xy = True)
        else:
            var, crs, lons, lats, lat_lims, lon_lims, plabel = read_MODIS_satpy(\
                date_str, str(channel), swath = swath)
    
    # Plot the GOES data
    # ------------------
    in_ax = True 
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure()
        if(use_xy):
            mapcrs = init_proj('202107222110')
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
            ax.set_title('MODIS '+ dt_date_str.strftime('%Y-%m-%d %H:%M'))
        else:
            ax.set_title('MODIS Ch ' + str(channel) + ' (' + \
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
            outname = 'modis__' + date_str + zoom_add + cmpst_add + '.png'
            plt.savefig(outname,dpi=300)
            print("Saved image",outname)
        else:
            plt.show()

def plot_MODIS_satpy_point_test(date_str, channel1 = 1, \
        channel2 = 7, channel3 = 31, \
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
        vmin = None, vmax = None, save = False):

    # Read MODIS image data
    # --------------------
    #print(MODIS_dict['dt_dates'][time_idx])
    #date_str = MODIS_dict['dt_dates'][time_idx].strftime('%Y%m%d%H%M')
    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H%M')
    var1, crs, lons1, lats1, lat_lims, lon_lims, plabel1 = \
        read_MODIS_satpy(date_str, channel1)
    var2, crs, lons2, lats2, lat_lims, lon_lims, plabel2 = \
        read_MODIS_satpy(date_str, channel2)
    var3, crs, lons3, lats3, lat_lims, lon_lims, plabel3 = \
        read_MODIS_satpy(date_str, channel3)

    plt.close('all')
    fig = plt.figure(figsize = (11, 4))
    ax1 = fig.add_subplot(1,3,1, projection = crs)
    ax2 = fig.add_subplot(1,3,2, projection = crs)
    ax3 = fig.add_subplot(1,3,3, projection = crs)

    labelsize = 10
    font_size = 8
    print('lim = ',aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
        dt_date_str.strftime('%H%M')]['data_lim'][channel1])
    plot_MODIS_satpy(date_str, channel1, \
        ax = ax1, var = var1, crs = crs, \
        lons = lons1, lats = lats1, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
            dt_date_str.strftime('%H%M')]['data_lim'][channel1][0], \
        vmax = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
            dt_date_str.strftime('%H%M')]['data_lim'][channel1][1], \
        #vmin = channel_dict[str(channel)]['limits'][0], \
        #vmax = channel_dict[str(channel)]['limits'][1], \
        ptitle = '', plabel = plabel2, \
        #vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)
    plot_MODIS_satpy(date_str, channel2, \
        ax = ax2, var = var2, crs = crs, \
        lons = lons2, lats = lats2, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
            dt_date_str.strftime('%H%M')]['data_lim'][channel2][0], \
        vmax = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
            dt_date_str.strftime('%H%M')]['data_lim'][channel2][1], \
        #vmin = channel_dict[str(channel)]['limits'][0], \
        #vmax = channel_dict[str(channel)]['limits'][1], \
        ptitle = '', plabel = plabel2, \
        #vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)
    plot_MODIS_satpy(date_str, channel3, \
        ax = ax3, var = var3, crs = crs, \
        lons = lons3, lats = lats3, lat_lims = lat_lims, lon_lims = lon_lims, \
        vmin = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
            dt_date_str.strftime('%H%M')]['data_lim'][channel3][0], \
        vmax = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
            dt_date_str.strftime('%H%M')]['data_lim'][channel3][1], \
        #vmin = channel_dict[str(channel)]['limits'][0], \
        #vmax = channel_dict[str(channel)]['limits'][1], \
        ptitle = '', plabel = plabel3, \
        #vmin = 5, vmax = 80, ptitle = '', plabel = plabel0, \
        colorbar = True, labelsize = labelsize + 1, zoom=True,save=False)


    # Add the locations of the points
    # -------------------------------
    width = height = 0.1
    for tlat, tlon in zip(plats, plons):
        ax1.add_patch(mpatches.Rectangle(\
            xy=[tlon - width/2, tlat - height/2], width=width, height=height,
            edgecolor = 'black', facecolor=None, alpha=0.8,
            transform=datacrs))
        ax2.add_patch(mpatches.Rectangle(\
            xy=[tlon - width/2, tlat - height/2], width=width, height=height,
            edgecolor = 'black', facecolor=None, alpha=0.8,
            transform=datacrs))
        ax3.add_patch(mpatches.Rectangle(\
            xy=[tlon - width/2, tlat - height/2], width=width, height=height,
            edgecolor = 'black', facecolor=None, alpha=0.8,
            transform=datacrs))
        ax1.plot(tlon, tlat, linewidth=2, markersize = 8, marker='.',
                 color = 'black', transform=datacrs)
        ax1.plot(tlon, tlat, linewidth=2, markersize = 5, marker='.',
                 transform=datacrs)
        ax2.plot(tlon, tlat, linewidth=2, markersize = 8, marker='.',
                 color = 'black', transform=datacrs)
        ax2.plot(tlon, tlat, linewidth=2, markersize = 5, marker='.',
                 transform=datacrs)
        ax3.plot(tlon, tlat, linewidth=2, markersize = 8, marker='.',
                 color = 'black', transform=datacrs)
        ax3.plot(tlon, tlat, linewidth=2, markersize = 5, marker='.',
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
       
    # Using the lat and lons of the point, pull out the data
    # for each channel 
    # -------------------------------------------------------
    arr_var1 = np.array(var1)
    arr_var2 = np.array(var2)
    arr_var3 = np.array(var3)

    points_var1 = np.array([arr_var1[\
        (lats1 > (tplat - (width / 2))) & \
        (lats1 <= (tplat + (width / 2))) & \
        (lons1 > (tplon - (height / 2))) & \
        (lons1 <= (tplon + (height / 2)))] \
        for tplat, tplon in zip(plats, plons)])
    points_var2 = np.array([arr_var2[\
        (lats2 > (tplat - (width / 2))) & \
        (lats2 <= (tplat + (width / 2))) & \
        (lons2 > (tplon - (height / 2))) & \
        (lons2 <= (tplon + (height / 2)))] \
        for tplat, tplon in zip(plats, plons)])
    points_var3 = np.array([arr_var3[\
        (lats3 > (tplat - (width / 2))) & \
        (lats3 <= (tplat + (width / 2))) & \
        (lons3 > (tplon - (height / 2))) & \
        (lons3 <= (tplon + (height / 2)))] \
        for tplat, tplon in zip(plats, plons)])

    # Print out the mean and std dev point values for each channel
    # ------------------------------------------------------------
    avg_diff1 = np.nanmean(points_var1, axis = 1)
    avg_diff2 = np.nanmean(points_var2, axis = 1)
    avg_diff3 = np.nanmean(points_var3, axis = 1)
    print("Channel 1")
    print(avg_diff1, avg_diff1[0] - avg_diff1[1])
    print("Channel 2")
    print(avg_diff2, avg_diff2[0] - avg_diff2[1])
    print("Channel 3")
    print(avg_diff3, avg_diff3[0] - avg_diff3[1])
    ##!#print("Point 1:",plats[0], plons[0])
    ##!#print("   channel 1", np.mean(arr_var1[0]))
    ##!#print("   channel 2", np.mean(arr_var2[0]))
    ##!#print("   channel 3", np.mean(arr_var3[0]))
    ##!#print("Point 2:",plats[1], plons[1])
    ##!#print("   channel 1", np.mean(arr_var1[1]))
    ##!#print("   channel 2", np.mean(arr_var2[1]))
    ##!#print("   channel 3", np.mean(arr_var3[1]))
 
    # Add plot labels
    # ---------------
    plot_figure_text(ax1, 'MODIS ' + \
        str(np.round(np.mean(channel_dict[str(channel1)]['Bandwidth']), 3)\
        ) + ' μm', \
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
        ax.set_extent([aerosol_event_dict[cross_date][cross_time]['Lon'][0], \
                       aerosol_event_dict[cross_date][cross_time]['Lon'][1], \
                       aerosol_event_dict[cross_date][cross_time]['Lat'][0], \
                       aerosol_event_dict[cross_date][cross_time]['Lat'][1]], \
            ccrs.PlateCarree())

    plt.show()

def read_MODIS_granule(filename, channel, zoom = False):

    modis = SD.SD(filename)

    dat = modis.attributes().get('CoreMetadata.0').split()
    indx = dat.index('EQUATORCROSSINGDATE')+9
    cross_date = dat[indx][1:len(dat[indx])-1]
    indx = dat.index('EQUATORCROSSINGTIME')+9
    cross_time = dat[indx][1:len(dat[indx])-1]

    if(debug):
        print('MODIS orbit info',cross_date,cross_time)

    lat5 = modis.select('Latitude').get()
    lon5 = modis.select('Longitude').get()

    if(str(channel)[:2] != 'wv'):
        sza = modis.select('SolarZenith').get() * modis.select('SolarZenith').attributes().get('scale_factor')
        vza = modis.select('SensorZenith').get() * modis.select('SensorZenith').attributes().get('scale_factor')

    if((str(channel)[2:] != '_ir') & (str(channel) != 'true_color')):
        data  = modis.select(channel_dict[str(channel)]['Name']).get()
        if(str(channel) != 'wv_nir'):
            data = data[channel_dict[str(channel)]['Index']]
        data  = data[::5,::5]
        # NOTE: Added on 20230627 to match these with the MYD06 files
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

        colors = 'Greys_r'
        label = 'IR Precipitable water vapor [cm]'
        ##!#try:
        ##!#    label = modis.select(channel_dict[str(channel)]['Name']\
        ##!#        ).attributes().get('long_name') +  ' [' + \
        ##!#        modis.select(channel_dict[str(channel)]['Name']\
        ##!#        ).attributes().get('units') + ']'
        ##!#except:
        ##!#    label = modis.select(channel_dict[str(channel)]['Name']\
        ##!#        ).attributes().get('long_name') +  ' [' + \
        ##!#        modis.select(channel_dict[str(channel)]['Name']\
        ##!#        ).attributes().get('unit') + ']'
    elif(str(channel) == 'true_color'):

        # Pull out the red, green, and blue reflectance data
        # --------------------------------------------------
        red   = modis.select('EV_250_Aggr1km_RefSB').get()[0]
        green = modis.select('EV_500_Aggr1km_RefSB').get()[1]
        blue  = modis.select('EV_500_Aggr1km_RefSB').get()[0]
        
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
        
        data = np.ma.masked_where(np.isnan(image),image)
        label = 'True color'
        colors = 'True color'
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

            if(debug):
                print("Average wavelength = ",np.average(channel_dict[str(channel)]['Bandwidth']))

            data = convert_radiance_to_temp(\
                (np.average(channel_dict[str(channel)]['Bandwidth'])), data)

            ##!## BEGIN FUNCTION

            ##!## Define constants for converting radiances to temperatures
            ##!#lmbda = (1e-6) * (np.average(channel_dict[str(channel)]['Bandwidth'])) # in m
            ##!#c_const = 3e8
            ##!#h_const = 6.626e-34 # J*s
            ##!#k_const = 1.381e-23 # J/K

            ##!#data = (h_const * c_const) / \
            ##!#    (lmbda * k_const * np.log( ((2.0 * h_const * (c_const**2.0) ) / \
            ##!#    ((lmbda**4.) * (lmbda / 1e-6) * data ) ) + 1.0 ) )
            ##!##data = ((h_const * c_const)/(k_const * lmbda)) * (np.log((2.0 * h_const * (c_const ** 2.0) / \
            ##!##    ((lmbda**5.0) * data)) + 1) ** -1.)

            ##!## END FUNCTION

            # Convert any missing data to nans
            # --------------------------------
            data[bad_locations] = np.nan
            data = np.ma.masked_invalid(data)

            colors = 'plasma'
            label = 'Brightness Temperature [K]'

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

    MODIS_data = {}

    if(data.shape[1] == 271):
        data = np.ma.masked_invalid(data[:,:-1])
    if(lon5.shape[1] == 271):
        lon5 = np.ma.masked_invalid(lon5[:,:-1])
    if(lat5.shape[1] == 271):
        lat5 = np.ma.masked_invalid(lat5[:,:-1])
    if(sza.shape[1] == 271):
        sza = np.ma.masked_invalid(sza[:,:-1])
    if(vza.shape[1] == 271):
        vza = np.ma.masked_invalid(vza[:,:-1])

    MODIS_data['data'] = data
    MODIS_data['lat']  = lat5
    MODIS_data['lon']  = lon5
    if(str(channel)[:2] != 'wv'):
        MODIS_data['sza']  = sza
        MODIS_data['vza']  = vza
    MODIS_data['variable']  = label
    MODIS_data['cross_date']  = cross_date
    MODIS_data['channel']  = channel
    MODIS_data['colors']  = colors
    MODIS_data['file_time'] = filename.strip().split('/')[-1].split('.')[2]
    if(str(channel) == 'true_color'):
        MODIS_data['colortuple'] = colortuple
    else:
        MODIS_data['colortuple'] = None

    if(zoom):
        # Mask MODIS_data['data'] that are outside the desired range
        # --------------------------------------------
        MODIS_data['data'][(((MODIS_data['lat'] < aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][0]) | \
                             (MODIS_data['lat'] > aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][1])) | \
                            ((MODIS_data['lon'] < aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][0]) | \
                             (MODIS_data['lon'] > aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][1])))] = -999.
        MODIS_data['lat'][ (((MODIS_data['lat'] < aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][0]) | \
                             (MODIS_data['lat'] > aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][1])) | \
                            ((MODIS_data['lon'] < aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][0]) | \
                             (MODIS_data['lon'] > aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][1])))] = -999.
        MODIS_data['lon'][ (((MODIS_data['lat'] < aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][0]) | \
                             (MODIS_data['lat'] > aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][1])) | \
                            ((MODIS_data['lon'] < aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][0]) | \
                             (MODIS_data['lon'] > aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][1])))] = -999.

        MODIS_data['data'] = np.ma.masked_where(MODIS_data['data'] == -999., MODIS_data['data'])
        MODIS_data['lat'] = np.ma.masked_where(MODIS_data['lat'] == -999., MODIS_data['lat'])
        MODIS_data['lon'] = np.ma.masked_where(MODIS_data['lon'] == -999., MODIS_data['lon'])

    return MODIS_data

# dt_date_str is of format YYYYMMDDHHMM
def read_MODIS_channel(date_str, channel, zoom = False, swath = False, \
        include_cloud_mask = False, include_myd06 = False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # Extract the filename given the channel
    # --------------------------------------
    if(str(channel)[:2] == 'wv'):
        filename = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['mdswv']
    else:

        if(swath):
            times = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['swath']
            dt_times = [datetime.strptime(time, '%Y%m%d%H%M') for time in times]
            filename = [aerosol_event_dict[ttime.strftime('%Y-%m-%d')][\
                ttime.strftime('%H%M')]['modis'] for ttime in dt_times]
           
            if(include_cloud_mask):
                cloud_name = [aerosol_event_dict[ttime.strftime('%Y-%m-%d')][\
                    ttime.strftime('%H%M')]['modis_cloud'] for ttime in dt_times]
            if(include_myd06):
                myd06_name = [aerosol_event_dict[ttime.strftime('%Y-%m-%d')][\
                    ttime.strftime('%H%M')]['modis_myd06'] for ttime in dt_times]
 
        else:
            filename = [aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
                dt_date_str.strftime('%H%M')]['modis']]
            if(include_cloud_mask):
                cloud_name = [aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
                    dt_date_str.strftime('%H%M')]['modis_cloud']]
            if(include_myd06):
                myd06_name = [aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
                    dt_date_str.strftime('%H%M')]['modis_myd06']]
   
    MODIS_holder = {}
    counter = 0
    for ii, ifile in enumerate(filename):

        print("Reading MODIS channel",channel," from ",ifile)

        # Set the channel to an arbitrary "1" in this case
        # ------------------------------------------------ 
        if(channel == 'cloud_mask'):
            include_cloud_mask = True
            local_channel = 1
        else:
            local_channel = channel

            
        MODIS_data = read_MODIS_granule(ifile, local_channel, zoom = zoom)

        if(include_cloud_mask):
            cloud_data = Dataset(cloud_name[ii])
            MODIS_data['cloud_mask'] = cloud_data[\
                'geophysical_data/Integer_Cloud_Mask'][::5,::5]
            if(MODIS_data['cloud_mask'].shape[1] == 271):
                MODIS_data['cloud_mask'] = MODIS_data['cloud_mask'][:,:-1]
            cloud_data.close()
        if(include_myd06):
            myd06_data = Dataset(myd06_name[ii])
        
            # Load cloud optical depth data
            # -----------------------------
            local_data = myd06_data['Cloud_Optical_Thickness'][::5,::5]
            local_data[local_data < 0] = 0.0
            local_data[local_data > 10000] = np.nan
            MODIS_data['cloud_optical_depth'] = local_data
            #MODIS_data['cloud_optical_depth'] = local_data[::5,::5] * \
            #    myd06_data['Cloud_Optical_Thickness'].scale_factor
            if(MODIS_data['cloud_optical_depth'].shape[1] == 271):
                MODIS_data['cloud_optical_depth'] = MODIS_data['cloud_optical_depth'][:,:-1]
            MODIS_data['cloud_optical_depth'] = np.ma.masked_invalid(\
                MODIS_data['cloud_optical_depth'])

            # Load cloud top pressure data
            # ----------------------------
            local_data = myd06_data['Cloud_Top_Pressure'][:,:]
            local_data[local_data < 0] = 0.0
            local_data[local_data > 10000] = np.nan
            MODIS_data['cloud_top_pressure'] = local_data
            #MODIS_data['cloud_optical_depth'] = local_data[::5,::5] * \
            #    myd06_data['Cloud_Optical_Thickness'].scale_factor
            #if(MODIS_data['cloud_optical_depth'].shape[1] == 271):
            #    MODIS_data['cloud_optical_depth'] = MODIS_data['cloud_optical_depth'][:,:-1]
            MODIS_data['cloud_top_pressure'] = np.ma.masked_invalid(\
                MODIS_data['cloud_top_pressure'])


            myd06_data.close()
 
        MODIS_holder[str(counter)] = MODIS_data 

        counter += 1

    MODIS_final = {}
    MODIS_final['data'] = np.concatenate([MODIS_holder[key]['data'] for key in MODIS_holder.keys()], axis = 0)
    MODIS_final['lat']  = np.concatenate([MODIS_holder[key]['lat']  for key in MODIS_holder.keys()], axis = 0)
    MODIS_final['lon']  = np.concatenate([MODIS_holder[key]['lon']  for key in MODIS_holder.keys()], axis = 0)
    MODIS_final['sza']  = np.concatenate([MODIS_holder[key]['sza']  for key in MODIS_holder.keys()], axis = 0)
    MODIS_final['vza']  = np.concatenate([MODIS_holder[key]['vza']  for key in MODIS_holder.keys()], axis = 0)
    if(include_cloud_mask):
        MODIS_final['cloud_mask'] = \
            np.concatenate([MODIS_holder[key]['cloud_mask']  \
            for key in MODIS_holder.keys()], axis = 0)
    if(include_myd06):
        MODIS_final['cloud_optical_depth'] = \
            np.concatenate([MODIS_holder[key]['cloud_optical_depth']  \
            for key in MODIS_holder.keys()], axis = 0)
        MODIS_final['cloud_top_pressure'] = \
            np.concatenate([MODIS_holder[key]['cloud_top_pressure']  \
            for key in MODIS_holder.keys()], axis = 0)
    MODIS_final['variable'] = MODIS_holder['0']['variable']
    MODIS_final['cross_date'] = MODIS_holder['0']['cross_date']
    MODIS_final['channel'] = channel
    #MODIS_final['channel'] = MODIS_holder['0']['channel']
    MODIS_final['colors'] = MODIS_holder['0']['colors']
    MODIS_final['file_time'] = MODIS_holder['0']['file_time']
    MODIS_final['date'] = date_str
    
    print(MODIS_final['data'].shape, MODIS_final['lat'].shape, MODIS_final['lon'].shape)

    if(str(channel) == 'true_color'):
    
        # Convert the color values into a format usable for plotting
        # ----------------------------------------------------------
        colortuple = tuple(np.array([MODIS_final['data'][:,:,0].flatten(), \
            MODIS_final['data'][:,:,1].flatten(), \
            MODIS_final['data'][:,:,2].flatten()]).transpose().tolist())
        MODIS_final['colortuple'] = colortuple
    else: 
        MODIS_final['colortuple'] = None
    if(swath):
        MODIS_final['swath_times'] = times
    else:
        MODIS_final['swath_times'] = 'NONE'

    return MODIS_final

# YYYYMMDD
def read_MODIS_CLDL3_daily(date_str, minlat = 65.5):

    dt_date_str = datetime.strptime(date_str, '%Y%m%d')

    filename = glob(dt_date_str.strftime(cloudL3_dir + \
        'daily/MCD06COSP*.A%Y%j*.nc'))
    if(len(filename) == 0):
        print("ERROR: No MODIS CLDL3 file found for date:",date_str)
        return

    print(filename[0])
    data = Dataset(filename[0], 'r')

    # Check the lat ranges
    # --------------------
    lat_idxs = np.where(data['latitude'][:] >= minlat)[0]

    out_dict = {}
    out_dict['lat']           = data['latitude'][lat_idxs]
    out_dict['lon']           = data['longitude'][:]
    out_dict['cld_frac_mean'] = data['Cloud_Mask_Fraction/Mean'][:,lat_idxs].T
    out_dict['cld_frac_std']  = data['Cloud_Mask_Fraction/Standard_Deviation'][:,lat_idxs].T
    out_dict['ob_count']      = data['Cloud_Mask_Fraction/Pixel_Counts'][:,lat_idxs].T
    out_dict['date_str']      = date_str
   
    data.close()
 
    return out_dict
   
def read_MODIS_CLDL3_daily_allmonth(date_str, minlat = 65.5):

    dt_date_str = datetime.strptime(date_str, '%Y%m')

    # Grab all the filenames for the matching year
    # --------------------------------------------
    allfiles = glob(dt_date_str.strftime(cloudL3_dir + 'daily/MCD*.A%Y*.nc'))

    # Grab just the dates
    # -------------------
    datetimes  = np.array([datetime.strptime(\
        tname.strip().split('/')[-1].split('.')[1],'A%Y%j') \
        for tname in allfiles])

    # Figure out which files match with the desired month
    # ---------------------------------------------------
    just_months = np.array([dtime.month for dtime in datetimes])
    good_idxs = np.where(just_months == dt_date_str.month)[0]

    good_times = datetimes[good_idxs]
    good_dstrs = [gtime.strftime('%Y%m%d') for gtime in good_times]

    # Figure out which lat indices
    # ----------------------------
    lat_ranges = np.arange(-89.5, 90.5, 1.0)
    keep_lats = np.where(lat_ranges >= minlat)[0]

    # Make an array to hold all the data for the month
    # ------------------------------------------------
    all_cldfrac_data = np.full((len(good_dstrs),len(keep_lats),360), np.nan)
    all_std_data     = np.full((len(good_dstrs),len(keep_lats),360), np.nan)
    all_counts_data  = np.full((len(good_dstrs),len(keep_lats),360), np.nan)
    all_lats = lat_ranges[keep_lats]
    all_lons = np.arange(-179.5, 180.5, 1.0)

    for ii, dstr in enumerate(good_dstrs):
        print(dstr)
        cldfrac_data = read_MODIS_CLDL3_daily(dstr, minlat = minlat)
    
        all_cldfrac_data[ii,:,:] = cldfrac_data['cld_frac_mean']
        all_std_data[ii,:,:]     = cldfrac_data['cld_frac_std']
        all_counts_data[ii,:,:]  = cldfrac_data['ob_count']

    out_dict = {}
    out_dict['date_str']       = date_str
    out_dict['day_strs']       = good_dstrs
    out_dict['lat']            = all_lats
    out_dict['lon']            = all_lons
    out_dict['cld_frac_mean']  = np.nanmean(all_cldfrac_data, axis = 0)

    # Calculate the average of each of the daily standard deviations,
    # following the methodology I found on :
    # https://www.statology.org/averaging-standard-deviations/
    # -----------------------------------------------------------------
    out_dict['cld_frac_std'] = \
        np.sqrt((np.sum((all_counts_data - 1) * (all_std_data**2.), \
        axis = 0)) / \
        (np.sum(all_counts_data, axis = 0) - all_counts_data.shape[0]))

    out_dict['ob_count'] = np.sum(all_counts_data, axis = 0)
    
    del all_cldfrac_data
    del all_std_data
    del all_counts_data

    return out_dict

def read_MODIS_CLDL3_single_month(date_str, minlat = 65.5, \
        local_data_dir = cloudL3_monthly_dir,):

    dt_date_str = datetime.strptime(date_str, '%Y%m')

    # Grab all the filenames for the matching year
    # --------------------------------------------
    found_file = glob(dt_date_str.strftime(\
        local_data_dir + 'modis_cldL3_monthly_%Y%m.nc'))[0]

    print(found_file)
    data = Dataset(found_file,'r')

    out_dict = {}
    out_dict['date_str']        = date_str
    #out_dict['lat']             = data['Latitude'][:,:]
    #out_dict['lon']             = data['Longitude'][:,:]
    out_dict['cld_frac_mean']   = data['Cloud_Fraction_Mean'][:,:]
    out_dict['cld_frac_std']    = data['Cloud_Fraction_StDev'][:,:]
    out_dict['ob_count']        = data['Ob_Counts'][:,:]

    out_dict['cld_frac_mean'] = np.ma.masked_where(\
        out_dict['cld_frac_mean'] < 0, out_dict['cld_frac_mean'])
    
    # I forgot to add the latitudes to the NETCDF monthly writer, so
    # none of the processed monthly files have latitude or longitude.
    # Back fill the lat and lon here
    # ---------------------------------------------------------------
    #if(out_dict['lat'].compressed().size == 0):
    temp_lats = np.arange(-89.5,90.5,1)
    temp_lons = np.arange(-179.5,180.5,1)
    temp_lats = temp_lats[temp_lats >= minlat]
    out_dict['lat'] = temp_lats
    out_dict['lon'] = temp_lons
        
    data.close()

    return out_dict

def read_MODIS_MYD08_single(date_str, minlat = 65.5, \
        maxlat = 89.5, local_data_dir = myd08_dir, \
        transform_lat = True):

    if(len(date_str) == 6):
        date_fmt = '%Y%m'
        dtype = 'monthly'
    elif(len(date_str) == 8):
        date_fmt = '%Y%m%d'
        dtype = 'daily'
    else:
        print("ERROR: INVALID DATE FORMAT")
        return

    dt_date_str = datetime.strptime(date_str, date_fmt)

    local_data_dir = local_data_dir + dtype + '/'

    # Grab all the filenames for the matching year
    # --------------------------------------------
    found_file = glob(dt_date_str.strftime(\
        local_data_dir + 'modis_MYD08_subset_' + date_fmt + \
            '.nc'))[0]

    print(found_file)
    data = Dataset(found_file,'r')

    out_dict = {}
    out_dict['date_str']        = date_str
    out_dict['lat']             = data['Latitude'][:]

    keep_idxs = np.where((out_dict['lat'] >= minlat) & \
        (out_dict['lat'] <= maxlat))[0]

    out_dict['lat']               = data['Latitude'][keep_idxs]
    out_dict['lon']               = data['Longitude'][:]
    out_dict['cld_frac_mean']     = data['Cloud_Fraction_Mean'][keep_idxs,:]
    out_dict['cld_frac_std']      = data['Cloud_Fraction_StDev'][keep_idxs,:]
    out_dict['ob_count']          = data['Ob_Counts'][keep_idxs,:]
    out_dict['day_cld_frac_mean'] = data['Cloud_Fraction_Day_Mean'][keep_idxs,:]
    out_dict['day_cld_frac_std']  = data['Cloud_Fraction_Day_StDev'][keep_idxs,:]
    out_dict['day_ob_count']      = data['Day_Ob_Counts'][keep_idxs,:]


    # Since the MYD08 trends start at 89.5 and work down, while the CLDL3
    # trends start at the min lat and work up, the MYD08 trends need to 
    # be switched around
    # -------------------------------------------------------------------
    if(transform_lat):
        comp_vars = ['cld_frac_mean', 'cld_frac_std', 'ob_count',\
            'day_cld_frac_mean','day_cld_frac_std','day_ob_count']
        for cvar in comp_vars:
            new_data = np.full(out_dict[cvar].shape, np.nan)
            for ii in range(new_data.shape[0]-1, -1, -1):
                new_data[new_data.shape[0]-1 - ii,:] = out_dict[cvar][ii,:]
            out_dict[cvar] = new_data
        out_dict['lat'] = out_dict['lat'][::-1] 
        
    out_dict['cld_frac_mean'] = np.ma.masked_where(\
        (out_dict['cld_frac_mean'] < 0) | (out_dict['cld_frac_mean'] > 1.0) | \
        (out_dict['cld_frac_mean'] == data['Cloud_Fraction_Mean'][:].fill_value), \
        out_dict['cld_frac_mean'])
    out_dict['day_cld_frac_mean'] = np.ma.masked_where(\
        (out_dict['day_cld_frac_mean'] < 0) | \
        (out_dict['day_cld_frac_mean'] > 1.0), out_dict['day_cld_frac_mean'])


    data.close()

    return out_dict

# 2021/01/28: modified function to allow for Aqua data
def read_MODIS_CLDL3_monthrange(start_date,end_date,\
        minlat=65.5, calc_month = True,season = 'sunlight'):
    CERES_data = {}
  
    spring=False
    summer=False
    autumn=False
    winter=False
    sunlight=False
    if(season=='spring'):
        keep_months = [3,4,5]
        spring=True
    elif(season=='summer'):
        keep_months = [6,7,8]
        summer=True
    elif(season=='autumn'):
        keep_months = [9,10,11]
        autumn=True
    elif(season=='winter'):
        keep_months = [12,1,2]
        winter=True
    elif(season=='sunlight'):
        keep_months = [4,5,6,7,8,9]
        sunlight=True
  
    dt_start_date = datetime.strptime(start_date, '%Y%m')
    dt_end_date   = datetime.strptime(end_date, '%Y%m')
 
    # Grab all the filenames for the matching year
    # --------------------------------------------
    allfiles = glob(cloudL3_monthly_dir + 'modis*.nc')

    # Grab just the dates
    # -------------------
    datetimes  = np.array([datetime.strptime(\
        tname.strip().split('/')[-1].split('_')[3],'%Y%m.nc') \
        for tname in allfiles])

    # Figure out which files match with the desired month
    # ---------------------------------------------------
    good_idxs = np.where([tmonth.month in keep_months \
        for tmonth in datetimes])
    good_times = datetimes[good_idxs]

    # Make sure the files are within the desired time range
    # -----------------------------------------------------
    within_range = np.where((good_times >= dt_start_date) & \
                            (good_times <= dt_end_date))

    good_times = np.array(good_times[within_range])
    good_dstrs = np.array([gtime.strftime('%Y%m') for gtime in good_times])


    # Declare arrays to hold the data
    # -------------------------------
    cloud_data = read_MODIS_CLDL3_single_month(good_dstrs[0], \
        minlat = minlat)
    lat_ranges = cloud_data['lat'][:]
    lon_ranges = cloud_data['lon'][:]

    lat_indices = np.where(lat_ranges>=minlat)[0]
    print(lat_indices)

    # Grid the lat and data
    grid_lon, grid_lat = np.meshgrid(lon_ranges,lat_ranges[lat_indices])

    cloudL3_data = {}
    cloudL3_data['cld_frac_mean']  = np.full((len(good_dstrs),\
        len(lat_indices),len(lon_ranges)), np.nan)
    cloudL3_data['cld_frac_std']   = np.full((len(good_dstrs),\
        len(lat_indices),len(lon_ranges)), np.nan)
    cloudL3_data['ob_count']       = np.full((len(good_dstrs),\
        len(lat_indices),len(lon_ranges)), np.nan)
    cloudL3_data['dates']          = good_dstrs
    cloudL3_data['lat']            = grid_lat
    cloudL3_data['lon']            = grid_lon
    cloudL3_data['season'] =season

    del cloud_data

    # Loop over good files and insert data into dictionary
    for ii, dstr in enumerate(good_dstrs):
        print(dstr)
        cloud_data = read_MODIS_CLDL3_single_month(dstr, minlat = minlat)

        cloudL3_data['cld_frac_mean'][ii,:,:] = \
            cloud_data['cld_frac_mean'][lat_indices,:]
        cloudL3_data['cld_frac_std'][ii,:,:] = \
            cloud_data['cld_frac_std'][lat_indices,:]
        cloudL3_data['ob_count'][ii,:,:] = \
            cloud_data['ob_count'][lat_indices,:]

    cloudL3_data['cld_frac_mean'] = np.ma.masked_where(\
        cloudL3_data['cld_frac_mean'] < 0, cloudL3_data['cld_frac_mean'])

    ##!#if(calc_month == True):
    ##!#    CERES_data = calcCERES_MonthClimo(CERES_data)

    return cloudL3_data

# 2021/01/28: modified function to allow for Aqua data
def read_MODIS_MYD08_monthrange(start_date,end_date,\
        minlat=65.5, maxlat = 89.5, calc_month = True,season = 'sunlight'):
    CERES_data = {}
  
    spring=False
    summer=False
    autumn=False
    winter=False
    sunlight=False
    if(season=='spring'):
        keep_months = [3,4,5]
        spring=True
    elif(season=='summer'):
        keep_months = [6,7,8]
        summer=True
    elif(season=='autumn'):
        keep_months = [9,10,11]
        autumn=True
    elif(season=='winter'):
        keep_months = [12,1,2]
        winter=True
    elif(season=='sunlight'):
        keep_months = [4,5,6,7,8,9]
        sunlight=True
  
    dt_start_date = datetime.strptime(start_date, '%Y%m')
    dt_end_date   = datetime.strptime(end_date, '%Y%m')
 
    # Grab all the filenames for the matching year
    # --------------------------------------------
    allfiles = glob(myd08_dir + 'monthly/modis*subset*.nc')

    # Grab just the dates
    # -------------------
    datetimes  = np.array([datetime.strptime(\
        tname.strip().split('/')[-1].split('_')[3],'%Y%m.nc') \
        for tname in allfiles])

    # Figure out which files match with the desired month
    # ---------------------------------------------------
    good_idxs = np.where([tmonth.month in keep_months \
        for tmonth in datetimes])
    good_times = datetimes[good_idxs]

    # Make sure the files are within the desired time range
    # -----------------------------------------------------
    within_range = np.where((good_times >= dt_start_date) & \
                            (good_times <= dt_end_date))

    good_times = np.array(good_times[within_range])
    good_dstrs = np.array([gtime.strftime('%Y%m') for gtime in good_times])

    # Declare arrays to hold the data
    # -------------------------------
    cloud_data = read_MODIS_MYD08_single(good_dstrs[0], \
        minlat = minlat)
    lat_ranges = cloud_data['lat'][:]
    lon_ranges = cloud_data['lon'][:]

    lat_indices = np.where((lat_ranges>=minlat) & (lat_ranges <= maxlat))[0]
    print(lat_indices, lat_ranges[lat_indices])

    # Grid the lat and data
    grid_lon, grid_lat = np.meshgrid(lon_ranges,lat_ranges[lat_indices])

    myd08_data = {}
    myd08_data['cld_frac_mean']  = np.full((len(good_dstrs),\
        len(lat_indices),len(lon_ranges)), np.nan)
    myd08_data['cld_frac_std']   = np.full((len(good_dstrs),\
        len(lat_indices),len(lon_ranges)), np.nan)
    myd08_data['ob_count']       = np.full((len(good_dstrs),\
        len(lat_indices),len(lon_ranges)), np.nan)
    myd08_data['day_cld_frac_mean']  = np.full((len(good_dstrs),\
        len(lat_indices),len(lon_ranges)), np.nan)
    myd08_data['day_cld_frac_std']   = np.full((len(good_dstrs),\
        len(lat_indices),len(lon_ranges)), np.nan)
    myd08_data['day_ob_count']       = np.full((len(good_dstrs),\
        len(lat_indices),len(lon_ranges)), np.nan)
    myd08_data['dates']          = good_dstrs
    myd08_data['lat']            = grid_lat
    myd08_data['lon']            = grid_lon
    myd08_data['season'] =season

    del cloud_data

    # Loop over good files and insert data into dictionary
    for ii, dstr in enumerate(good_dstrs):
        print(dstr)
        cloud_data = read_MODIS_MYD08_single(dstr, minlat = minlat, \
            maxlat = maxlat)

        myd08_data['cld_frac_mean'][ii,:,:] = \
            cloud_data['cld_frac_mean'][:,:]
        myd08_data['cld_frac_std'][ii,:,:] = \
            cloud_data['cld_frac_std'][:,:]
        myd08_data['ob_count'][ii,:,:] = \
            cloud_data['ob_count'][:,:]
        myd08_data['day_cld_frac_mean'][ii,:,:] = \
            cloud_data['day_cld_frac_mean'][:,:]
        myd08_data['day_cld_frac_std'][ii,:,:] = \
            cloud_data['day_cld_frac_std'][:,:]
        myd08_data['day_ob_count'][ii,:,:] = \
            cloud_data['day_ob_count'][:,:]

    myd08_data['cld_frac_mean'] = np.ma.masked_where(\
        myd08_data['cld_frac_mean'] < 0, myd08_data['cld_frac_mean'])
    myd08_data['day_cld_frac_mean'] = np.ma.masked_where(\
        myd08_data['day_cld_frac_mean'] < 0, myd08_data['day_cld_frac_mean'])

    ##!#if(calc_month == True):
    ##!#    CERES_data = calcCERES_MonthClimo(CERES_data)

    return myd08_data


def calcMODIS_CLDL3_grid_trend(cloud_data, month_idx, trend_type, minlat,\
        norm_to_decade = False):

    if(month_idx == None):
        month_idx = 0
        index_jumper = 1
    else:
        month_adder = '_month'
        if(cloud_data['season'] == 'sunlight'):
            index_jumper = 6
        else:   
            index_jumper = 12

    lat_ranges = cloud_data['lat'][:,0]
    #lat_ranges = np.arange(minlat,90,1.0)
    lon_ranges = cloud_data['lon'][0,:]
    #lon_ranges = np.arange(-180,180,1.0)

    # Make copy of cloud_data array
    print(cloud_data['dates'][month_idx::index_jumper])
    local_data   = np.copy(cloud_data['cld_frac_mean'][month_idx::index_jumper,:,:])
    #local_counts = np.copy(cloud_data['OB_COUNT'][month_idx::index_jumper,:,:])
    #local_mask = np.ma.masked_where(local_counts == 0, local_data)
    local_mask = np.ma.masked_where(((local_data < 0) | \
        (cloud_data['lat'] < minlat)), local_data)
    ceres_trends = np.full(local_data.shape[1:], np.nan)
    ceres_pvals  = np.full(local_data.shape[1:], np.nan)
    ceres_uncert = np.full(local_data.shape[1:], np.nan)

    # Loop over all the keys and print the regression slopes 
    # Grab the averages for the key
    for i in range(0,len(lat_ranges)):
        for j in range(0,len(lon_ranges)):
            # Check the current max and min
            #print(local_mask[:,i,j])
            work_mask = local_mask[:,i,j]
            #work_mask = local_mask[:,i,j][~local_mask[:,i,j].mask][0]
            if(len(work_mask.compressed()) > 1):
                x_vals = np.arange(0,len(work_mask.compressed()))
                # Find the slope of the line of best fit for the time series of
                # average data
                if((trend_type=='standard') | (trend_type == 'linregress')): 
                    result = stats.linregress(x_vals, work_mask.compressed())
                    #slope, intercept, r_value, p_value, std_err = \
                    #    stats.linregress(x_vals,work_mask.compressed())
                    ceres_trends[i,j] = result.slope * len(x_vals)
                    ceres_pvals[i,j]  = result.pvalue
                    ceres_uncert[i,j] = result.stderr * len(x_vals)
                else:
                    res = stats.theilslopes(work_mask.compressed(), x_vals, 0.90)
                    ceres_trends[i,j] = res[0]*len(x_vals)
            else:
                print('no data')

    #ceres_trends = np.ma.masked_where(((cloud_data['lat'] < minlat) | \
    #    (ceres_trends == -999.)), ceres_trends)
    ceres_trends = np.ma.masked_where(cloud_data['lat'] < minlat, ceres_trends)
    ceres_pvals  = np.ma.masked_where(cloud_data['lat'] < minlat, ceres_pvals)
    ceres_uncert = np.ma.masked_where(cloud_data['lat'] < minlat, ceres_uncert)

    # If the user wants to normalize the trends so that they reflect the
    # trend per decade...
    # ------------------------------------------------------------------
    if(norm_to_decade):
        local_dates = cloud_data['dates'][month_idx::index_jumper]
        dt_begin_date = datetime.strptime(local_dates[0],'%Y%m')
        dt_end_date   = datetime.strptime(local_dates[-1],'%Y%m')
        num_decade = (dt_end_date.year - dt_begin_date.year) / 10
        print("Number of decades:",num_decade)
        ceres_trends = ceres_trends / num_decade

    return ceres_trends, ceres_pvals, ceres_uncert

def calcMODIS_MYD08_grid_trend(cloud_data, month_idx, trend_type, minlat,\
        norm_to_decade = False, dtype = 'all'):


    if(dtype == 'day'):
        dtype_adder = dtype + '_'
        print("Calculating daylight only")
    else:
        dtype_adder = ''

    if(month_idx == None):
        month_idx = 0
        index_jumper = 1
    else:
        month_adder = '_month'
        if(cloud_data['season'] == 'sunlight'):
            index_jumper = 6
        else:   
            index_jumper = 12

    lat_ranges = cloud_data['lat'][:,0]
    #lat_ranges = np.arange(minlat,90,1.0)
    lon_ranges = cloud_data['lon'][0,:]
    #lon_ranges = np.arange(-180,180,1.0)

    # Make copy of cloud_data array
    print(cloud_data['dates'][month_idx::index_jumper])
    local_data   = np.copy(cloud_data[dtype_adder + \
        'cld_frac_mean'][month_idx::index_jumper,:,:])
    #local_counts = np.copy(cloud_data['OB_COUNT'][month_idx::index_jumper,:,:])
    #local_mask = np.ma.masked_where(local_counts == 0, local_data)
    local_mask = np.ma.masked_where(((local_data < 0) | \
        (cloud_data['lat'] < minlat)), local_data)
    cloud_trends = np.full(local_data.shape[1:], np.nan)
    cloud_pvals  = np.full(local_data.shape[1:], np.nan)
    cloud_uncert = np.full(local_data.shape[1:], np.nan)

    # Loop over all the keys and print the regression slopes 
    # Grab the averages for the key
    for i in range(0,len(lat_ranges)):
        for j in range(0,len(lon_ranges)):
            # Check the current max and min
            #print(local_mask[:,i,j])
            work_mask = local_mask[:,i,j]
            #work_mask = local_mask[:,i,j][~local_mask[:,i,j].mask][0]
            if(len(work_mask.compressed()) > 1):
                x_vals = np.arange(0,len(work_mask.compressed()))
                # Find the slope of the line of best fit for the time series of
                # average data
                if((trend_type=='standard') | (trend_type == 'linregress')): 
                    result = stats.linregress(x_vals, work_mask.compressed())
                    #slope, intercept, r_value, p_value, std_err = \
                    #    stats.linregress(x_vals,work_mask.compressed())
                    cloud_trends[i,j] = result.slope * len(x_vals)
                    cloud_pvals[i,j]  = result.pvalue
                    cloud_uncert[i,j] = result.stderr * len(x_vals)
                else:
                    res = stats.theilslopes(work_mask.compressed(), x_vals, 0.90)
                    cloud_trends[i,j] = res[0]*len(x_vals)
            else:
                print('no data')

    #ceres_trends = np.ma.masked_where(((cloud_data['lat'] < minlat) | \
    #    (ceres_trends == -999.)), ceres_trends)
    cloud_trends = np.ma.masked_where(cloud_data['lat'] < minlat, cloud_trends)
    cloud_pvals  = np.ma.masked_where(cloud_data['lat'] < minlat, cloud_pvals)
    cloud_uncert = np.ma.masked_where(cloud_data['lat'] < minlat, cloud_uncert)

    # If the user wants to normalize the trends so that they reflect the
    # trend per decade...
    # ------------------------------------------------------------------
    if(norm_to_decade):
        local_dates = cloud_data['dates'][month_idx::index_jumper]
        dt_begin_date = datetime.strptime(local_dates[0],'%Y%m')
        dt_end_date   = datetime.strptime(local_dates[-1],'%Y%m')
        num_decade = (dt_end_date.year - dt_begin_date.year) / 10
        print("Number of decades:",num_decade)
        cloud_trends = cloud_trends / num_decade

    return cloud_trends, cloud_pvals, cloud_uncert



##!#def calcOMI_MonthClimo(cloud_data):
##!#
##!#    # Set up arrays to hold monthly climatologies
##!#    month_climo = np.zeros((6,cloud_data['cld_frac_mean'].shape[1],\
##!#        cloud_data['cld_frac_mean'].shape[2]))
##!#
##!#    # Mask the monthly averages
##!#    local_mean   = np.copy(cloud_data['cld_frac_mean'][:,:,:])
##!#    local_stdev  = np.copy(cloud_data['cld_frac_mean'][:,:,:])
##!#    local_counts = np.copy(cloud_data['ob_counts'][:,:,:])
##!# 
##!#    # Calculate monthly climatologies
##!#    for m_i in range(6):
##!#        month_climo[m_i,:,:] = np.nanmean(local_mask[m_i::6,:,:],axis=0)
##!#        print("Month: ",OMI_data['DATES'][m_i][4:]) 
##!#
##!#    # Insert data into dictionary
##!#    cloud_data['MONTH_CLIMO'] = month_climo
##!#
##!#    return cloud_data

def write_MODIS_CLDL3_monthly(cloud_data, minlat = 65.5, \
        save_path = cloudL3_monthly_dir):

    if(isinstance(cloud_data, str)):
        dt_date_str = datetime.strptime(cloud_data, '%Y%m')
        cloud_data = read_MODIS_CLDL3_daily_allmonth(cloud_data, minlat = minlat)
    else:
        dt_date_str = datetime.strptime(cloud_data['date_str'], '%Y%m')

    # Create a new netCDF dataset to write to the file
    # ------------------------------------------------
    outfile = dt_date_str.strftime(save_path + \
        'modis_cldL3_monthly_%Y%m.nc')

    nc = Dataset(outfile,'w',format='NETCDF4')
  
    # Dimensions = lon, lat
    # Create the sizes of each dimension in the file. In this case,
    # the dimensions are "# of latitude" and "# of longitude"
    # -------------------------------------------------------------
    num_y = cloud_data['cld_frac_mean'].shape[0]
    num_x = cloud_data['cld_frac_mean'].shape[1]
    
    # Use the dimension size variables to actually create dimensions in 
    # the file.
    # ----------------------------------------------------------------- 
    n_y  = nc.createDimension('ny',num_y)
    n_x  = nc.createDimension('nx',num_x)

    # Create variables for the three dimensions. Note that since these
    # are variables, they are still given 'dimensions' using 'dlat', 
    # and 'dlon'. Latitude and longitude are each 2-d grids, so they are 
    # given 2 dimensions (dlat, dlon).
    # ------------------------------------------------------------------
    LON=nc.createVariable('Longitude','f4',('nx'))
    LON.description='Longitude (-180 to 180)'
    LON.units='Degrees'
    LAT=nc.createVariable('Latitude','f4',('ny'))
    LAT.description='Latitude (-90 to 90)'
    LAT.units='Degrees'

    # Create a variable for the AI data, dimensioned using both 
    # dimensions.
    # ---------------------------------------------------------  
    CLD_FRAC_MEAN = nc.createVariable('Cloud_Fraction_Mean','f4',('ny','nx'))
    CLD_FRAC_MEAN.description='Monthly averaged MODIS CLDL3 Cloud Fraction'

    CLD_FRAC_STD = nc.createVariable('Cloud_Fraction_StDev','f4',('ny','nx'))
    CLD_FRAC_STD.description='Average standard deviation of MODIS CLDL3 '+\
        'cloud fraction'

    OB_COUNT = nc.createVariable('Ob_Counts','f4',('ny','nx'))
    OB_COUNT.description='Total Observation Counts'
    
    # Fill in dimension variables one-by-one.
    # ---------------------------------------------------------------
    LAT[:]             = cloud_data['lat'][:]
    LON[:]             = cloud_data['lon'][:]
    CLD_FRAC_MEAN[:,:] = cloud_data['cld_frac_mean'][:,:]
    CLD_FRAC_STD[:,:]  = cloud_data['cld_frac_std'][:,:]
    OB_COUNT[:,:]      = cloud_data['ob_count'][:,:]

    # Save, write, and close the netCDF file
    # --------------------------------------
    nc.close()

    print("Saved file ",outfile)  
   
# NOTE: This is not set up to read in the original MYD08 files. This
#       just finds the corresponding file in the data directory and
#       writes the desired files and desired lat ranges to a subset
#       file. 
# date_str: 'YYYYMM' or 'YYYYMMDD'
def write_MODIS_MYD08(date_str, minlat = 65.5, \
        save_path = myd08_dir):

    # Create a new netCDF dataset to write to the file
    # ------------------------------------------------
    if(len(date_str) == 6):
        out_fmt = '%Y%m'
        dtype = 'monthly'
    elif(len(date_str) == 8):
        out_fmt = '%Y%m%d'
        dtype = 'daily'
    else:
        print("ERROR: Invalid data type. Must be YYYYMM or YYYYMMDD")
        return -3   

    dt_date_str = datetime.strptime(date_str, out_fmt)
    data = identify_MODIS_MYD08(date_str)
    data = Dataset(myd08_dir + dtype + '/' + data, 'r')

 
    outfile = dt_date_str.strftime(save_path + \
        dtype + '/modis_MYD08_subset_' + out_fmt + '.nc')

    nc = Dataset(outfile,'w',format='NETCDF4')
 
    # Figure out the latitude dimensions to use
    # -----------------------------------------
    good_idx = np.where(data['YDim'][:] >= minlat)[0]
 
    # Dimensions = lon, lat
    # Create the sizes of each dimension in the file. In this case,
    # the dimensions are "# of latitude" and "# of longitude"
    # -------------------------------------------------------------
    num_y = data['YDim'][good_idx].shape[0]
    num_x = data['XDim'].shape[0]
    
    # Use the dimension size variables to actually create dimensions in 
    # the file.
    # ----------------------------------------------------------------- 
    n_y  = nc.createDimension('ny',num_y)
    n_x  = nc.createDimension('nx',num_x)
    n_c  = nc.createDimension('nc',3)

    # Create variables for the three dimensions. Note that since these
    # are variables, they are still given 'dimensions' using 'dlat', 
    # and 'dlon'. Latitude and longitude are each 2-d grids, so they are 
    # given 2 dimensions (dlat, dlon).
    # ------------------------------------------------------------------
    LON=nc.createVariable('Longitude','f4',('nx'))
    LON.description='Longitude (-180 to 180)'
    LON.units='Degrees'
    LAT=nc.createVariable('Latitude','f4',('ny'))
    LAT.description='Latitude (-90 to 90)'
    LAT.units='Degrees'

    # Create a variable for the AI data, dimensioned using both 
    # dimensions.
    # ---------------------------------------------------------  
    CLD_FRAC_MEAN = nc.createVariable('Cloud_Fraction_Mean','f4',('ny','nx'))
    temp_desc = 'MODIS MYD08 daily mean '+\
        'cloud fraction'
    if(dtype == 'monthly'):
        temp_desc = 'Monthly averaged ' + temp_desc
    CLD_FRAC_MEAN.description=temp_desc

    if(dtype == 'monthly'):
        var_desc = 'Mean standard deviation of MODIS MYD08 '+\
        'cloud fraction daily means'
    else:
        var_desc = 'Standard deviation of MODIS MYD08 '+\
        'cloud fraction'       
 
    CLD_FRAC_STD = nc.createVariable('Cloud_Fraction_StDev','f4',('ny','nx'))
    CLD_FRAC_STD.description=var_desc

    COD_MEAN = nc.createVariable('Cloud_Optical_Depth_Mean','f4',('ny','nx'))
    temp_desc = 'MODIS MYD08 daily mean '+\
        'cloud optical depth'
    if(dtype == 'monthly'):
        temp_desc = 'Monthly averaged ' + temp_desc
    COD_MEAN.description=temp_desc

    if(dtype == 'monthly'):
        var_desc = 'Mean standard deviation of MODIS MYD08 '+\
        'cloud optical depth daily means'
    else:
        var_desc = 'Standard deviation of MODIS MYD08 '+\
        'cloud optical depth'       
 
    COD_STD = nc.createVariable('Cloud_Optical_Depth_StDev','f4',('ny','nx'))
    COD_STD.description=var_desc


    OB_COUNT = nc.createVariable('Ob_Counts','f4',('ny','nx'))
    OB_COUNT.description='Total Observation Counts'
    
    DAY_CLD_FRAC_MEAN = nc.createVariable('Cloud_Fraction_Day_Mean','f4',('ny','nx'))
    DAY_CLD_FRAC_MEAN.description='Monthly averaged MODIS MYD08 daily mean '+\
        'cloud fraction for day only'

    DAY_CLD_FRAC_STD = nc.createVariable('Cloud_Fraction_Day_StDev','f4',('ny','nx'))
    DAY_CLD_FRAC_STD.description='Mean standard deviation of MODIS MYD08 '+\
        'cloud fraction daily means for day only'

    DAY_OB_COUNT = nc.createVariable('Day_Ob_Counts','f4',('ny','nx'))
    DAY_OB_COUNT.description='Total observation counts for day only'
   
    COP_CLOUDMASKCLEAR_COUNT = nc.createVariable('CloudMaskClear_Counts','i4', ('ny','nx'))
    COP_CLOUDMASKCLEAR_COUNT.description = \
        'Cloud Optical Properties Cloud Phase (Clear, CloudMaskClear(CSR=0))'

    COP_RESTORED_CLEAR_COUNT = nc.createVariable('RestoredToClear_Counts','i4', ('ny','nx'))
    COP_RESTORED_CLEAR_COUNT.description = \
        'Cloud Optical Properties Cloud Phase (Clear, RestoredToClear(CSR=2))'

    COP_PARTLY_CLOUDY_COUNT = nc.createVariable('Partly_Cloudy_Counts','i4', ('nc','ny','nx'))
    COP_PARTLY_CLOUDY_COUNT.description = \
        'Cloud Optical Properties Cloud Phase ' + \
        '(Partly Cloudy, CSR=1 or CSR=3): Pixel Counts ' + \
        'in 3 cats: Liquid(1st) Ice(2nd) Undet(3rd)'

    COP_CLOUDMASKCLOUDY_COUNT = nc.createVariable('CloudMaskCloudy_Counts','i4', ('nc','ny','nx'))
    COP_CLOUDMASKCLOUDY_COUNT.description = \
        'Cloud Optical Properties Cloud Phase ' + \
        '(Cloudy, Not Restored to Clear, CSR=0): Pixel Counts in 3 cats: ' + \
        'Liquid(1st) Ice(2nd) Undet(3rd)'


    # Fill in dimension variables one-by-one.
    # ---------------------------------------------------------------
    LAT[:]                 = data['YDim'][good_idx]
    LON[:]                 = data['XDim'][:]
    if(dtype == 'monthly'):
        CLD_FRAC_MEAN[:,:]     = data['Cloud_Fraction_Mean_Mean'][good_idx,:]
        CLD_FRAC_STD[:,:]      = data['Cloud_Fraction_Mean_Std'][good_idx,:]
        COD_FRAC[:,:]          = data['Cloud_Optical_Thickness_Combined_Mean_Mean'][good_idx,:]
        COD_STD[:,:]           = data['Cloud_Optical_Thickness_Combined_Mean_Std'][good_idx,:]
        OB_COUNT[:,:]          = data['Cloud_Fraction_Pixel_Counts'][good_idx,:]
        DAY_CLD_FRAC_MEAN[:,:] = data['Cloud_Fraction_Day_Mean_Mean'][good_idx,:]
        DAY_CLD_FRAC_STD[:,:]  = data['Cloud_Fraction_Day_Mean_Std'][good_idx,:]
        DAY_OB_COUNT[:,:]      = data['Cloud_Fraction_Day_Pixel_Counts'][good_idx,:]
        #NOTE: NEED TO FIGURE OUT THE NAMES OF THE HISTOGRAM COUNT VARIABLES IN 
        #      MONTHLY FILES
    else:
        CLD_FRAC_MEAN[:,:]     = data['Cloud_Fraction_Mean'][good_idx,:]
        CLD_FRAC_STD[:,:]      = data['Cloud_Fraction_Standard_Deviation'][good_idx,:]
        COD_MEAN[:,:]          = data['Cloud_Optical_Thickness_Combined_Mean'][good_idx,:]
        COD_STD[:,:]           = data['Cloud_Optical_Thickness_Combined_Standard_Deviation'][good_idx,:]
        OB_COUNT[:,:]          = data['Cloud_Fraction_Pixel_Counts'][good_idx,:]
        DAY_CLD_FRAC_MEAN[:,:] = data['Cloud_Fraction_Day_Mean'][good_idx,:]
        DAY_CLD_FRAC_STD[:,:]  = data['Cloud_Fraction_Day_Standard_Deviation'][good_idx,:]
        DAY_OB_COUNT[:,:]      = data['Cloud_Fraction_Day_Pixel_Counts'][good_idx,:]
        COP_CLOUDMASKCLEAR_COUNT[:,:]  = data['COP_Phase_CloudMaskClear_Histogram_Counts'][0, good_idx,:]
        COP_RESTORED_CLEAR_COUNT[:,:]  = data['COP_Phase_RestoredToClear_Histogram_Counts'][0, good_idx,:]
        COP_PARTLY_CLOUDY_COUNT[:,:]   = data['COP_Phase_Partly_Cloudy_Histogram_Counts'][:,good_idx,:]
        COP_CLOUDMASKCLOUDY_COUNT[:,:] = data['COP_Phase_Cloudy_Histogram_Counts'][:,good_idx,:]

    data.close()

    # Save, write, and close the netCDF file
    # --------------------------------------
    nc.close()

    print("Saved file ",outfile)  

 
def plot_MODIS_channel(date_str,channel,zoom=True,show_smoke=False, \
        ax = None, swath = False, vmin = None, vmax = None, \
        ptitle = '', \
        circle_bound = True, plot_borders = True, labelsize = 12):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    #filename = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis']

    include_cloud_mask = False
    if(channel == 'red'):
        channel = 1
    elif(channel == 'green'):
        channel = 4
    elif(channel == 'blue'):
        channel = 3
    elif(channel == 'cloud_mask'):
        include_cloud_mask = True

    # Call read_MODIS_channel to read the desired MODIS data from the
    # file and put it in a dictionary
    # ---------------------------------------------------------------
    MODIS_data = read_MODIS_channel(date_str, channel, swath = swath, \
        include_cloud_mask = include_cloud_mask)

    if(debug):
        print("Data max = ",np.max(MODIS_data['data']), "  Data min = ",\
        np.min(MODIS_data['data']))

    in_ax = True 
    if(ax is None): 
        in_ax = False
        plt.close('all')
        fig1 = plt.figure()
        #fig1 = plt.figure(figsize = (6,6))
        mapcrs = ccrs.NorthPolarStereo()
        #mapcrs = ccrs.Robinson()
        ax = fig1.add_subplot(1,1,1, projection = mapcrs)

    plot_MODIS_spatial(MODIS_data, ax, zoom, circle_bound = circle_bound, \
        vmin = vmin, vmax = vmax, ptitle = ptitle, \
        plot_borders = plot_borders, labelsize = labelsize)

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

    if(swath and not zoom):
        ax.set_extent([-180, 180, 65, 90], datacrs)

    ##!#cbar = plt.colorbar(mesh,orientation='horizontal',pad=0,\
    ##!#    aspect=50,shrink = 0.850)
    ##!#cbar.ax.tick_params(labelsize=14)
    ##!#cbar.set_label(MODIS_data['variable'],fontsize=16,weight='bold')
    ##!#
    ##!#ax.add_feature(cfeature.BORDERS)
    ##!#ax.add_feature(cfeature.STATES)
    ##!#ax.coastlines()
    ##!#if(zoom):
    ##!#    ax.set_extent([aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][0], \
    ##!#                   aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][1], \
    ##!#                   aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][0], \
    ##!#                   aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][1]],\
    ##!#                   ccrs.PlateCarree())
    ##!#ax.set_title('Channel ' + str(channel) + '\n' + \
    ##!#    channel_dict[str(channel)]['Bandwidth_label']) 

    if(not in_ax):
        if(save):
            print("SAVE HERE")
        else:
            plt.show()

def plot_MODIS_channel_time_diff(date_str1,date_str2,channel,zoom=True,\
        show_smoke=False, ax = None, swath = False, vmin = None, vmax = None,\
        circle_bound = True, plot_borders = True):

    dt_date_str = datetime.strptime(date_str2,"%Y%m%d%H%M")
    #filename = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis']

    if(channel == 'red'):
        channel = 1
    elif(channel == 'green'):
        channel = 4
    elif(channel == 'blue'):
        channel = 3

    # Call read_MODIS_channel to read the desired MODIS data from the
    # file and put it in a dictionary
    # ---------------------------------------------------------------
    MODIS_data1 = read_MODIS_channel(date_str1, channel, swath = swath)
    MODIS_data2 = read_MODIS_channel(date_str2, channel, swath = swath)

    if(debug):
        print("Data max = ",np.max(MODIS_data1['data']), "  Data min = ",\
        np.min(MODIS_data1['data']))
        print("Data max = ",np.max(MODIS_data2['data']), "  Data min = ",\
        np.min(MODIS_data2['data']))

    in_ax = True 
    if(ax is None): 
        in_ax = False
        plt.close('all')
        fig1 = plt.figure()
        #fig1 = plt.figure(figsize = (6,6))
        mapcrs = ccrs.NorthPolarStereo()
        #mapcrs = ccrs.Robinson()
        ax = fig1.add_subplot(1,1,1, projection = mapcrs)

    MODIS_data2['data'] = MODIS_data2['data'] - MODIS_data1['data']
    MODIS_data2['colors'] = 'bwr'

    plot_MODIS_spatial(MODIS_data2, ax, zoom, circle_bound = circle_bound, \
        vmin = vmin, vmax = vmax, plot_borders = plot_borders)

    #mesh = ax.pcolormesh(MODIS_data['lon'],MODIS_data['lat'],\
    #    MODIS_data['data'],cmap = MODIS_data['colors'], shading='auto', \
    #    transform = datacrs) 

    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(date_str) 
        #hash_data1, nohash_data1 = find_plume(filename) 
        hash0 = ax.pcolor(MODIS_data2['lon'],MODIS_data2['lat'],\
            hash_data1, hatch = '///', alpha=0., transform = datacrs)

    if(swath and not zoom):
        ax.set_extent([-180, 180, 65, 90], datacrs)

    ##!#cbar = plt.colorbar(mesh,orientation='horizontal',pad=0,\
    ##!#    aspect=50,shrink = 0.850)
    ##!#cbar.ax.tick_params(labelsize=14)
    ##!#cbar.set_label(MODIS_data['variable'],fontsize=16,weight='bold')
    ##!#
    ##!#ax.add_feature(cfeature.BORDERS)
    ##!#ax.add_feature(cfeature.STATES)
    ##!#ax.coastlines()
    ##!#if(zoom):
    ##!#    ax.set_extent([aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][0], \
    ##!#                   aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][1], \
    ##!#                   aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][0], \
    ##!#                   aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][1]],\
    ##!#                   ccrs.PlateCarree())
    ##!#ax.set_title('Channel ' + str(channel) + '\n' + \
    ##!#    channel_dict[str(channel)]['Bandwidth_label']) 

    if(not in_ax):
        if(save):
            print("SAVE HERE")
        else:
            plt.show()

def read_OMI_match_MODIS(date_str, min_AI = -2e5, corners = False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    modis_date = dt_date_str.strftime('%Y-%m-%d')

    data = h5py.File(aerosol_event_dict[modis_date][date_str[8:]]['omi'],'r')
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
        ((LAT < aerosol_event_dict[modis_date][date_str[8:]]['Lat'][0]) | \
         (LAT > aerosol_event_dict[modis_date][date_str[8:]]['Lat'][1])) | \
        ((LON < aerosol_event_dict[modis_date][date_str[8:]]['Lon'][0]) | \
         (LON > aerosol_event_dict[modis_date][date_str[8:]]['Lon'][1]))), \
        mask_UVAI)

    data.close()

    if(corners):
        return LAT, LON, mask_UVAI, crnr_LAT, crnr_LON
    else:
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

    data = Dataset(aerosol_event_dict[modis_date][date_str[8:]]['ceres'],'r')
    LAT   = 90. - data.variables['Colatitude_of_CERES_FOV_at_surface'][:]
    LON   = data.variables['Longitude_of_CERES_FOV_at_surface'][:]
    LON[LON>179.99] = -360.+LON[LON>179.99]
    swflux  = data.variables['CERES_SW_TOA_flux___upwards'][:]
    lwflux  = data.variables['CERES_LW_TOA_flux___upwards'][:]
    time  = data.variables['time'][:]
    #local_time = np.array([base_date + relativedelta(days = ttime) for ttime in time])

    data.close()

    mask_LAT = LAT[ \
        (LAT >= aerosol_event_dict[modis_date][date_str[8:]]['Lat'][0]) & \
        (LAT <= aerosol_event_dict[modis_date][date_str[8:]]['Lat'][1]) & \
        (LON >= aerosol_event_dict[modis_date][date_str[8:]]['Lon'][0]) & \
        (LON <= aerosol_event_dict[modis_date][date_str[8:]]['Lon'][1]) & \
        (swflux > 0) & (lwflux > 0)]
    mask_LON = LON[ \
        (LAT >= aerosol_event_dict[modis_date][date_str[8:]]['Lat'][0]) & \
        (LAT <= aerosol_event_dict[modis_date][date_str[8:]]['Lat'][1]) & \
        (LON >= aerosol_event_dict[modis_date][date_str[8:]]['Lon'][0]) & \
        (LON <= aerosol_event_dict[modis_date][date_str[8:]]['Lon'][1]) & \
        (swflux > 0) & (lwflux > 0)]
    mask_swf = swflux[ \
        (LAT >= aerosol_event_dict[modis_date][date_str[8:]]['Lat'][0]) & \
        (LAT <= aerosol_event_dict[modis_date][date_str[8:]]['Lat'][1]) & \
        (LON >= aerosol_event_dict[modis_date][date_str[8:]]['Lon'][0]) & \
        (LON <= aerosol_event_dict[modis_date][date_str[8:]]['Lon'][1]) & \
        (swflux > 0) & (lwflux > 0)]
    mask_lwf = lwflux[ \
        (LAT >= aerosol_event_dict[modis_date][date_str[8:]]['Lat'][0]) & \
        (LAT <= aerosol_event_dict[modis_date][date_str[8:]]['Lat'][1]) & \
        (LON >= aerosol_event_dict[modis_date][date_str[8:]]['Lon'][0]) & \
        (LON <= aerosol_event_dict[modis_date][date_str[8:]]['Lon'][1]) & \
        (swflux > 0) & (lwflux > 0)]
    ##!#mask_time = local_time[ \
    ##!#    (LAT >= aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][0]) & \
    ##!#    (LAT <= aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lat'][1]) & \
    ##!#    (LON >= aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][0]) & \
    ##!#    (LON <= aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['Lon'][1]) & \
    ##!#    (swflux > 0) & (lwflux > 0)]

    return mask_LAT, mask_LON, mask_swf, mask_lwf

# date_str of format "YYYYMMDDHHMM"
def compare_MODIS_3panel(date_str,channel1,channel2,channel3,zoom=True,save=False,\
        plot_ASOS_loc = False, show_smoke = True, compare_OMI = False, \
        compare_CERES = False, return_MODIS = False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    filename = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis']

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

    cpy_1 = np.ma.masked_where((((MODIS_data1['lat'] < aerosol_event_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lat'][0]) | \
                         (MODIS_data1['lat'] > aerosol_event_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lat'][1])) | \
                        ((MODIS_data1['lon'] < aerosol_event_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lon'][0]) | \
                         (MODIS_data1['lon'] > aerosol_event_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lon'][1]))), cpy_1)
    cpy_2 = np.ma.masked_where((((MODIS_data2['lat'] < aerosol_event_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lat'][0]) | \
                         (MODIS_data2['lat'] > aerosol_event_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lat'][1])) | \
                        ((MODIS_data2['lon'] < aerosol_event_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lon'][0]) | \
                         (MODIS_data2['lon'] > aerosol_event_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lon'][1]))), cpy_2)
    cpy_3 = np.ma.masked_where((((MODIS_data3['lat'] < aerosol_event_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][0]) | \
                         (MODIS_data3['lat'] > aerosol_event_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'][1])) | \
                        ((MODIS_data3['lon'] < aerosol_event_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][0]) | \
                         (MODIS_data3['lon'] > aerosol_event_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lon'][1]))), cpy_3)


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
    plot_MODIS_spatial(MODIS_data1, ax0, ptitle = '', zoom = zoom)
    plot_MODIS_spatial(MODIS_data2, ax1, ptitle = '', zoom = zoom)
    plot_MODIS_spatial(MODIS_data3, ax2, ptitle = '', zoom = zoom)

    plot_figure_text(ax2, 'MODIS 11 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 15, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax0, 'MODIS 0.64 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 15, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax1, 'MODIS 1.24 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 15, backgroundcolor = 'white', halign = 'right')

    if(compare_OMI):
        print("Reading OMI data")
        LAT, LON, mask_UVAI = read_OMI_match_MODIS(date_str)
        plot_OMI_spatial(date_str, LAT, LON, mask_UVAI, axo, zoom = zoom)

    if(compare_CERES):
        print("Reading CERES data")

        mask_LAT, mask_LON, mask_swf, mask_lwf = read_CERES_match_MODIS(date_str)

        plot_CERES_spatial(date_str,mask_LAT, mask_LON, mask_swf, 'SWF', axcs, zoom = zoom)
        plot_CERES_spatial(date_str,mask_LAT, mask_LON, mask_lwf, 'LWF', axcl, zoom = zoom)

    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(dt_date_str.strftime('%Y%m%d%H%M')) 

        plt.rcParams.update({'hatch.color': 'r'})
        #ax3.pcolor(MODIS_data_ch1['lon'],MODIS_data_ch1['lat'],\
        #    hash_data1, hatch = '\\\\', alpha=0., transform = datacrs,\
        #    cmap = 'plasma')

        hatch_shape = '\\\\'
        hash0 = ax0.pcolor(MODIS_data1['lon'],MODIS_data1['lat'],\
            hash_data1[:MODIS_data1['lat'].shape[0], :MODIS_data1['lat'].shape[1]], hatch = hatch_shape, alpha=0., transform = datacrs)
        hash1 = ax1.pcolor(MODIS_data2['lon'],MODIS_data2['lat'],\
            hash_data1[:MODIS_data2['lat'].shape[0], :MODIS_data2['lat'].shape[1]], hatch = hatch_shape, alpha=0., transform = datacrs)
        hash2 = ax2.pcolor(MODIS_data3['lon'],MODIS_data3['lat'],\
            hash_data1[:MODIS_data3['lat'].shape[0], :MODIS_data3['lat'].shape[1]], hatch = hatch_shape, alpha=0., transform = datacrs)
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


    cross_date = MODIS_data1['cross_date']
    file_time  = MODIS_data1['file_time']

    fig.tight_layout()

    #plt.suptitle(cross_date + ' ' + file_time)
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
    filename = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis']

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

    if(debug):
        print("Data1 max = ",np.max(MODIS_data1['data']), "  Data1 min = ",\
        np.min(MODIS_data1['data']))
        print("Data2 max = ",np.max(MODIS_data2['data']), "  Data2 min = ",np.min(MODIS_data2['data']))
    
        print(MODIS_data1['data'].shape, MODIS_data2['data'].shape)

    # --------------------------------------
    # Step 2: Set up figure to have 3 panels
    # --------------------------------------
    datacrs = ccrs.PlateCarree() 
    #mapcrs = ccrs.LambertConformal()
    mapcrs = init_proj(date_str)
    #mapcrs = ccrs.LambertConformal(central_longitude = \
    #    np.mean(aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['Lon']),\
    #    central_latitude = \
    #    np.mean(aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['Lat']))

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
        ax0.set_extent([aerosol_event_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lon'][0], \
                        aerosol_event_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lon'][1], \
                        aerosol_event_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lat'][0], \
                        aerosol_event_dict[MODIS_data1['cross_date']][MODIS_data1['file_time']]['Lat'][1]],\
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
        ax1.set_extent([aerosol_event_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lon'][0], \
                        aerosol_event_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lon'][1], \
                        aerosol_event_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lat'][0], \
                        aerosol_event_dict[MODIS_data2['cross_date']][MODIS_data2['file_time']]['Lat'][1]],\
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

    fig.tight_layout()

    if(plot_ASOS_loc):
        print("Plotting ASOS site")
        plot_ASOS_locs(ax0,date_str)
        plot_ASOS_locs(ax1,date_str)

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
    filename = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['modis']

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

    if(debug):
        print("Data0 max = ",np.max(MODIS_data0['data']), "  Data0 min = ",\
        np.min(MODIS_data0['data']))
        print("Data1 max = ",np.max(MODIS_data1['data']), "  Data1 min = ",\
            np.min(MODIS_data1['data']))
        print("Data2 max = ",np.max(MODIS_data2['data']), "  Data2 min = ",\
            np.min(MODIS_data2['data']))
        print("Data3 max = ",np.max(MODIS_data3['data']), "  Data3 min = ",\
            np.min(MODIS_data3['data']))

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
        hash_data1)
    plot_scatter(ax1, tmp_data0, tmp_data2, MODIS_data0, MODIS_data2, \
        hash_data1)
    plot_scatter(ax2, tmp_data0, tmp_data3, MODIS_data0, MODIS_data3, \
        hash_data1)


    if(compare_OMI):
        print("Reading OMI data")

        #colocate_OMI(date_str, axo, avg_pixel=True)
        plot_scatter_OMI(date_str, MODIS_data0, axo, avg_pixel = avg_pixel)

        #hash_data1, nohash_data1 = find_plume(filename) 
        # Colocate the masked OMI data with the MODIS channel0 data
        # ---------------------------------------------------------
        print("Colocating OMI data")


            # Plot hash_avg_modis against mask_UVAI 

            # Plot nohash_avg_modis against mask_UVAI 
 
    # End compare_OMI

    if(compare_CERES):
        print("Reading CERES data")

        plot_scatter_CERES(date_str, MODIS_data0, axcs, avg_pixel = avg_pixel)

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

def plot_MODIS_spatial(MODIS_data, pax, zoom, vmin = None, vmax = None, \
        #pvar = 'data', ptitle = None, labelsize = 13, labelticksize = 11, \
        pvar = 'data', ptitle = '', labelsize = 13, labelticksize = 11, \
        circle_bound = False, plot_borders = True):

    print(MODIS_data['date'])
    local_dstr = datetime.strptime(MODIS_data['date'], '%Y%m%d%H%M')
    local_date = local_dstr.strftime('%Y-%m-%d')
    local_time = local_dstr.strftime('%H%M')

    if(vmin is None):
        if('data_lim' in aerosol_event_dict[local_date][local_time].keys()):
            vmin = aerosol_event_dict[local_date][local_time]['data_lim'][MODIS_data['channel']][0]
        else:
            vmin = np.nanmin(MODIS_data[pvar])
    if(vmax is None):
        if('data_lim' in aerosol_event_dict[local_date][local_time].keys()):
            vmax = aerosol_event_dict[local_date][local_time]['data_lim'][MODIS_data['channel']][1]
        else:
            if(MODIS_data['channel'] == 'wv_ir'):
                vmax = 1.5
            else:
                vmax = np.nanmax(MODIS_data[pvar])

    #mesh = ax.pcolormesh(MODIS_data['lon'],MODIS_data['lat'],\
    #    MODIS_data['data'],cmap = MODIS_data['colors'], shading='auto', \
    #    transform = datacrs) 

    if(str(MODIS_data['channel']) == 'true_color'):
        pax.pcolormesh(MODIS_data['lon'],MODIS_data['lat'],\
            MODIS_data['data'][:,:,0],color= MODIS_data['colortuple'], \
            shading='auto', transform = ccrs.PlateCarree()) 
    elif(MODIS_data['channel'] == 'cloud_mask'):

        """
        cmap2 = plt.get_cmap('tab10')
        colorvals = np.arange(0, 11, 1)
        norm = cm.BoundaryNorm(colorvals, cmap2.N)
        """

        col_dict = {
            0: 'tab:blue', \
            1: 'tab:orange', \
            2: 'tab:green', \
            3: 'tab:red', \
        }
        labels = np.array([
            'Cloudy',
            'Prob.\nCloudy',
            'Prob.\nClear',
            'Clear',
        ])

        ccmm = cm.ListedColormap([col_dict[x] for x in col_dict.keys()]) 
        #bounds = np.array([k for k in col_dict.keys()])
        bounds = np.array([-0.5,0.5,1.5,2.5,3.5])
        #bounds = bounds[:-1] + np.diff(bounds) / 2
        #bounds = np.concatenate((np.array([0]), bounds, np.array([24])))
        norm = cm.BoundaryNorm(bounds, ccmm.N)
    
        print("UNIQUE:", np.unique(MODIS_data['cloud_mask'][:,:]))

        mesh = pax.pcolormesh(MODIS_data['lon'],MODIS_data['lat'],\
                       #MODIS_data['cloud_mask'][:,:], cmap = 'tab10', \
                       MODIS_data['cloud_mask'][:,:], cmap = ccmm, \
                       norm = norm, shading='auto', \
                       transform = ccrs.PlateCarree()) 

        cbar = plt.colorbar(ScalarMappable(norm = norm, cmap = ccmm),\
            ticks = [0,1,2,3], 
            ax = pax, orientation='vertical', pad = 0.03, fraction = 0.052)
        cbar.ax.set_yticklabels(labels)
    else:
        # Plot channel 1
        mesh1 = pax.pcolormesh(MODIS_data['lon'],MODIS_data['lat'],\
            MODIS_data[pvar],cmap = MODIS_data['colors'], shading='auto', \
            vmin = vmin, \
            vmax = vmax, transform = datacrs) 
    
        #cbar1 = plt.colorbar(mesh1,ax=pax,orientation='vertical',\
        #    pad=0.03, fraction = 0.046)
        cbar1 = plt.colorbar(mesh1, ax = pax, pad = 0.03, fraction = 0.052, \
            extend = 'both')
            #shrink = 0.870, pad=0.03,label=MODIS_data['variable'])
        cbar1.set_label(MODIS_data['variable'], size = labelsize, weight = 'bold')
        cbar1.ax.tick_params(labelsize = labelticksize)

    if(plot_borders):
        pax.add_feature(cfeature.BORDERS)
        pax.add_feature(cfeature.STATES)
        pax.coastlines()
    if(zoom):
        pax.set_extent([aerosol_event_dict[local_date][local_time]['Lon'][0], \
                        aerosol_event_dict[local_date][local_time]['Lon'][1], \
                        aerosol_event_dict[local_date][local_time]['Lat'][0], \
                        aerosol_event_dict[local_date][local_time]['Lat'][1]],\
                        datacrs)
    else:
        if(circle_bound):
            pax.set_boundary(circle, transform=pax.transAxes)

    #if(ptitle == None):
    if(ptitle == ''):
        if(str(MODIS_data['channel']) == 'true_color'):
            pax.set_title('True color') 
        elif(MODIS_data['channel'] == 'cloud_mask'):
            pax.set_title('MODIS Cloud Mask')
        else:
            pax.set_title('Channel ' + str(MODIS_data['channel']) + '\n' + \
                channel_dict[str(MODIS_data['channel'])]['Bandwidth_label']) 
    else:
        pax.set_title(ptitle)


# dtype must be 'SWF' or 'LWF'
def plot_CERES_spatial(date_str, mask_LAT, mask_LON, mask_data, dtype, pax, \
        vmin = None, vmax = None, markersize = 170, ptitle = None, zoom = False):

    #print(aerosol_event_dict[MODIS_data3['cross_date']][MODIS_data3['file_time']]['Lat'])

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
        pax.set_extent([aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][1], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][1]],\
                        datacrs)
    cbar3 = plt.colorbar(mesh3,ax=pax,orientation='vertical',\
        pad=0.03)
    cbar3.set_label(plabel, size = 14, weight = 'bold')
    cbar3.ax.tick_params(labelsize = 12)
    if(ptitle == None):
        pax.set_title('CERES TOA ' + dtype)
    else:
        pax.set_title(ptitle)

def plot_scatter(lax,data1, data2, MODIS_data1, MODIS_data2, hash_data1,\
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
    print("Calculating MODIS " + str(MODIS_data1['channel']) + ' vs ' + \
        str(MODIS_data2['channel']) + ' smoke trend')
    plot_trend_line(lax, \
        data1[np.where(~hash_data1.mask)].compressed(), \
        data2[np.where(~hash_data1.mask)].compressed(), \
        color='tab:blue', slope = slope)

    print("Calculating MODIS " + str(MODIS_data1['channel']) + ' vs ' + \
        str(MODIS_data2['channel']) + ' clear trend')
    plot_trend_line(lax, \
        data1[np.where(hash_data1.mask)].compressed(), \
        data2[np.where(hash_data1.mask)].compressed(), \
        color='tab:orange')

    if(MODIS_data2['channel'] == 'wv_ir'):
        lax.set_ylim(lax.get_ylim()[0], 2.0)
        #vmax = 1.5
    #else:
    #    vmax = np.nanmax(MODIS_data['data'])


    if(xlabel == None):
        xlabel = 'Ch. ' + str(MODIS_data1['channel']) +' [' + \
        channel_dict[str(MODIS_data1['channel'])]['Bandwidth_label'] + '] '+ \
        MODIS_data1['variable']
    if(ylabel == None):
        ylabel ='Ch. ' + str(MODIS_data2['channel']) +' [' + \
        channel_dict[str(MODIS_data2['channel'])]['Bandwidth_label'] + ']' + \
        MODIS_data2['variable'] 

    lax.set_xlabel(xlabel, fontsize = 12)
    lax.set_ylabel(ylabel, fontsize = 12)

    if(plot_legend):
        lax.legend()

    #plt.suptitle('Aqua MODIS ' + MODIS_data1['cross_date'] + ' ' + \
    #    MODIS_data1['file_time'])
    #lax.set_title('Smoke correlation: '+str(np.round(hrval_p, 3))+'\n'+\
    #              'Clear correlation: '+str(np.round(nrval_p, 3)))
   
 
#def plot_scatter_OMI(date_str, tmp_data0, tmp_lat0, tmp_lon0, channel, \
def plot_scatter_OMI(date_str, MODIS_data, pax, avg_pixel = False, \
        xlabel = None, ylabel = None, ptitle = None, labelsize = 13,\
        labelticksize = 11):

    # Determine where the smoke is located
    # ------------------------------------
    hash_data, nohash_data = find_plume(date_str) 

    tmp_data = np.copy(MODIS_data['data'])
    tmp_lat  = np.copy(MODIS_data['lat'])
    tmp_lon  = np.copy(MODIS_data['lon'])

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

    # Remove the nans from the MODIS data, lats, and lons for both the
    # in-the-plume (hashed) and outside-the-plume (nohashed) data
    hash_plot_data   = tmp_data[np.where(~hash_data.mask)].compressed()
    nohash_plot_data = tmp_data[np.where(hash_data.mask)].compressed()
    hash_plot_lat    = tmp_lat[np.where(~hash_data.mask)].compressed()
    nohash_plot_lat  = tmp_lat[np.where(hash_data.mask)].compressed()
    hash_plot_lon    = tmp_lon[np.where(~hash_data.mask)].compressed()
    nohash_plot_lon  = tmp_lon[np.where(hash_data.mask)].compressed()

    if(avg_pixel):

        LAT, LON, mask_UVAI, crnr_LAT, crnr_LON = \
            read_OMI_match_MODIS(date_str, corners = True, min_AI = 1.0)

        print('averaging MODIS pixels')

        # Pull out the OMI lat/lon corner pixels that are associated with
        # the masked OMI data
        # ---------------------------------------------------------------
        work_UVAI  = mask_UVAI.compressed()

        work_crnrLAT = crnr_LAT[np.where(~np.isnan(mask_UVAI)==True)[0],\
            np.where(~np.isnan(mask_UVAI)==True)[1],:] 
        work_crnrLON = crnr_LON[np.where(~np.isnan(mask_UVAI)==True)[0],\
            np.where(~np.isnan(mask_UVAI)==True)[1],:] 

        # Declare full arrays to hold the averaged MODIS data
        # Must be dimensioned to the shape of the OMI data
        # ---------------------------------------------------
        if(debug):
            print("work UVAI shape = ",work_UVAI.shape)
        hash_avg_modis   = np.full(work_UVAI.shape,np.nan)
        hash_avg_mlat    = np.full(work_UVAI.shape,np.nan)
        hash_avg_mlon    = np.full(work_UVAI.shape,np.nan)
        nohash_avg_modis = np.full(work_UVAI.shape,np.nan)
        nohash_avg_mlat  = np.full(work_UVAI.shape,np.nan)
        nohash_avg_mlon  = np.full(work_UVAI.shape,np.nan)

        # Loop over screened out OMI pixels
        # ---------------------------------
        if(debug):
            print("corners shape:",work_crnrLAT.shape)
        for ii in range(work_crnrLAT.shape[0]):

            # For each OMI pixel, find the in-plume and out-of-plume MODIS
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

            # Use the indices to pull out the MODIS pixels
            hash_modis_data   = hash_plot_data[hash_in]
            hash_modis_lat    = hash_plot_lat[hash_in]
            hash_modis_lon    = hash_plot_lon[hash_in]
            nohash_modis_data = nohash_plot_data[nohash_in]
            nohash_modis_lat  = nohash_plot_lat[nohash_in]
            nohash_modis_lon  = nohash_plot_lon[nohash_in]

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
  
            # If at least 1 smokey MODIS pixel is within the OMI pixel bounds, 
            # classify the OMI pixel as smokey
            # ----------------------------------------------------------------
            if((len(inside_hash) > 0) & (True in inside_hash)):
                hash_avg_modis[ii] = np.average(hash_modis_data[inside_hash])
                hash_avg_mlat[ii]  = np.average(hash_modis_lat[inside_hash])
                hash_avg_mlon[ii]  = np.average(hash_modis_lon[inside_hash])

            if((len(inside_nohash) > 0) & (True in inside_nohash)):
                nohash_avg_modis[ii] = np.average(nohash_modis_data[inside_nohash])
                nohash_avg_mlat[ii]  = np.average(nohash_modis_lat[inside_nohash])
                nohash_avg_mlon[ii]  = np.average(nohash_modis_lon[inside_nohash])

        # end screened corner data loop

        # Plot the screened OMI data vs the averaged MODIS data
        # -----------------------------------------------------

        #print(work_UVAI[~np.ma.masked_invalid(hash_avg_modis).mask])
        
        hrval_p = spearmanr(np.ma.masked_invalid(hash_avg_modis).compressed(), \
            work_UVAI[~np.ma.masked_invalid(hash_avg_modis).mask])[0]
        pax.scatter(np.ma.masked_invalid(hash_avg_modis).compressed(), \
            work_UVAI[~np.ma.masked_invalid(hash_avg_modis).mask], s=6, \
            color='tab:blue')
        #axo.scatter(np.ma.masked_invalid(nohash_avg_modis).compressed(), \
        #    work_UVAI[~np.ma.masked_invalid(nohash_avg_modis).mask], s=6, \
        #    color='tab:orange')
        #axo.scatter(nohash_plot_data0, nohash_match_OMI, s=6, \
        #    color='tab:orange', label='Outside Plume')

        # Plot trend lines for each set
        # -----------------------------
        print("Calculating MODIS/OMI trend")
        plot_trend_line(pax, \
            np.ma.masked_invalid(hash_avg_modis).compressed(), \
            work_UVAI[~np.ma.masked_invalid(hash_avg_modis).mask], \
            color='tab:blue')

        if(xlabel == None):
            xlabel = 'Averaged Ch. ' + str(MODIS_data['channel']) +' [' + \
            channel_dict[str(MODIS_data['channel'])]['Bandwidth_label'] + \
            MODIS_data['variable']
        pax.set_xlabel(xlabel, fontsize = labelsize, weight = 'bold')
        if(ptitle is None):
            ptitle = 'Smoke correlation: '+str(np.round(hrval_p, 3))
        pax.set_title(ptitle)
        pax.set_ylabel('OMI UVAI', fontsize = labelsize, weight = 'bold')

    else:
        LAT, LON, mask_UVAI, crnr_LAT, crnr_LON = \
            read_OMI_match_MODIS(date_str, corners = True)

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
        print("Calculating MODIS/OMI trend")
        plot_trend_line(pax, \
            hash_plot_data, \
            hash_match_OMI, \
            color='tab:blue')
        #axo.scatter(nohash_plot_data0, nohash_match_OMI, s=6, \
        #    color='tab:orange')

        if(xlabel == None):
            xlabel = 'Ch. ' + str(MODIS_data['channel']) +' [' + \
            channel_dict[str(MODIS_data['channel'])]['Bandwidth_label'] + \
            MODIS_data['variable']
        pax.set_xlabel(xlabel, fontsize = labelsize, weight = 'bold')
        if(ptitle is None):
            ptitle = 'Smoke correlation: '+str(np.round(hrval_p, 3))
        pax.set_title(ptitle)
        pax.set_ylabel('OMI UVAI', fontsize = labelsize, weight = 'bold')
        

def plot_scatter_CERES(date_str, MODIS_data, pax, avg_pixel = False,\
        plume_only = True, plot_total_flux = True, labelsize = 13,\
        labelticksize = 11, lw_pax = None, total_pax = None):

    # Determine where the smoke is located
    # ------------------------------------
    hash_data, nohash_data = find_plume(date_str) 

    tmp_data = np.copy(MODIS_data['data'])
    tmp_lat  = np.copy(MODIS_data['lat'])
    tmp_lon  = np.copy(MODIS_data['lon'])

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
        # Remove the nans from the MODIS data, lats, and lons for both the
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
    mask_LAT, mask_LON, mask_swf, mask_lwf = read_CERES_match_MODIS(date_str)

    crnr_LAT = getCorners_1d(mask_LAT.data) ; crnr_LON = getCorners_1d(mask_LON)


    #if(avg_pixel):
    
    # Average a certain number of MODIS pixels around each CERES
    # pixel.
    # ----------------------------------------------------------

    # total arrays are dimensioned to the size of the CERES
    # arrays. Contain all matched MODIS data (in and out)
    total_match_SWF    = np.full(mask_swf.shape, np.nan)
    total_match_LWF    = np.full(mask_lwf.shape, np.nan)
    total_match_LAT    = np.full(mask_LAT.shape, np.nan)
    total_match_LON    = np.full(mask_LON.shape, np.nan)
    # hash arrays are dimensioned to the size of the CERES
    # arrays. Contain only in-plume matched MODIS data
    hash_match_SWF    = np.full(mask_swf.shape, np.nan)
    hash_match_LWF    = np.full(mask_lwf.shape, np.nan)
    hash_match_LAT    = np.full(mask_LAT.shape, np.nan)
    hash_match_LON    = np.full(mask_LON.shape, np.nan)
    # nohash arrays are dimensioned to the size of the CERES
    # arrays. Contain only out=of-plume matched MODIS data
    nohash_match_SWF    = np.full(mask_swf.shape, np.nan)
    nohash_match_LWF    = np.full(mask_lwf.shape, np.nan)
    nohash_match_LAT    = np.full(mask_LAT.shape, np.nan)
    nohash_match_LON    = np.full(mask_LON.shape, np.nan)

    #print(hash_plot_lat.shape)

    # Loop over the lower-resolution CERES data.    
    # Either grab the nearest MODIS pixel to the CERES pixel
    # OR average the MODIS pixels within a certain range of the
    # CERES pixel.
    # ---------------------------------------------------------
    for ii in range(mask_swf.shape[0]):

        lat_list = sorted([crnr_LAT[ii], crnr_LAT[ii+1]])
        lon_list = sorted([crnr_LON[ii], crnr_LON[ii+1]])

        # NEW APPROACH: use calculated pixel corners and use where 
        # statement to find all MODIS pixels within the CERES bounds
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
            # Find the gridpoint in the MODIS lat/lon data that 
            # corresponds to the AIRS pixel
            # ---------------------------------------------------- 
            ho_idx = nearest_gridpoint(mask_LAT[ii], mask_LON[ii],\
                hash_plot_lat, hash_plot_lon)

            if(len(ho_idx[0]) > 1):
                ho_idx = (np.array([ho_idx[0][0]])), \
                    (np.array([ho_idx[1][0]]))
   
        if(len(hash_plot_data[ho_idx]) > 0):
            # if any smoke-classified MODIS pixel is within the CERES pixel,
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

            # Average the n MODIS pixels around the current CERES
            # pixel. 
            # total_match arrays contain all colocated MODIS data
            # both in-plume and out-of-plume.
            # --------------------------------------------------- 
            total_match_SWF[ii] = np.nanmean(hash_plot_data[ho_idx])
            total_match_LWF[ii] = np.nanmean(hash_plot_data[ho_idx])
            total_match_LAT[ii] = mask_LAT[ii]
            total_match_LON[ii] = mask_LON[ii]

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
  
    total_mask_swf  = np.ma.masked_invalid(mask_swf[~total_match_SWF.mask])
    total_mask_lwf  = np.ma.masked_invalid(mask_lwf[~total_match_LWF.mask])

    hash_mask_swf  = np.ma.masked_invalid(mask_swf[~hash_match_SWF.mask])
    hash_mask_lwf  = np.ma.masked_invalid(mask_lwf[~hash_match_LWF.mask])
 
    nohash_mask_swf  = np.ma.masked_invalid(mask_swf[~nohash_match_SWF.mask])
    nohash_mask_lwf  = np.ma.masked_invalid(mask_lwf[~nohash_match_LWF.mask])

    mask_total = mask_swf + mask_lwf

    markersize = 14
    lhrval_p = spearmanr(hash_match_LWF.compressed(), \
        hash_mask_lwf)[0]
    lnrval_p = spearmanr(nohash_match_LWF.compressed(), \
        nohash_mask_lwf)[0]
    ln1 = pax.scatter(hash_match_SWF.compressed(), hash_mask_swf,\
        s = markersize, color='tab:blue', marker='o',label = 'Smoke')
    ln2 = pax.scatter(nohash_match_SWF.compressed(), nohash_mask_swf,\
        s = markersize, color='tab:orange', marker='o',label = 'Clear')
    if(lw_pax is None):
        ln3 = pax.scatter(hash_match_LWF.compressed(), hash_mask_lwf,\
            s = markersize, color='tab:green', marker='o',label = 'Smoke')
        ln4 = pax.scatter(nohash_match_LWF.compressed(), nohash_mask_lwf,\
            s = markersize, color='tab:red', marker='o',label = 'Clear')
    else:
        ln3 = lw_pax.scatter(hash_match_LWF.compressed(), hash_mask_lwf,\
            s = markersize, color='tab:blue', marker='o',label = 'Smoke')
        ln4 = lw_pax.scatter(nohash_match_LWF.compressed(), nohash_mask_lwf,\
            s = markersize, color='tab:orange', marker='o',label = 'Clear')

    ##!#shrval_p = spearmanr(hash_match_SWF.compressed(), \
    ##!#    hash_mask_swf)[0]
    ##!#snrval_p = spearmanr(nohash_match_SWF.compressed(), \
    ##!#    nohash_mask_swf)[0]

    # Plot trend lines for each set
    # -----------------------------
    print("Calculating MODIS/SWF smoke trend")
    plot_trend_line(pax, hash_match_SWF.compressed(), hash_mask_swf, \
        color='tab:blue')
    print("Calculating MODIS/SWF clear trend")
    plot_trend_line(pax, nohash_match_SWF.compressed(), nohash_mask_swf, \
        color='tab:orange')

    print("Calculating MODIS/LWF smoke trend")
    print("Calculating MODIS/LWF clear trend")
    if(lw_pax is None):
        plot_trend_line(pax, hash_match_LWF.compressed(), hash_mask_lwf, \
            color='tab:green')
        plot_trend_line(pax, nohash_match_LWF.compressed(), nohash_mask_lwf, \
            color='tab:red')
    else:
        plot_trend_line(lw_pax, hash_match_LWF.compressed(), hash_mask_lwf, \
            color='tab:blue')
        plot_trend_line(lw_pax, nohash_match_LWF.compressed(), nohash_mask_lwf, \
            color='tab:orange')

    #lns = ln1 + ln2 + ln3 + ln4
    custom_lines = [Line2D([0], [0], color='tab:blue', label = 'SWF Smoke',),
                    Line2D([0], [0], color='tab:orange', label = 'SWF Clear',),
                    Line2D([0], [0], color='tab:green', label = 'LWF Smoke',),\
                    Line2D([0], [0], color='tab:red', label = 'LWF Clear')]

    if(plot_total_flux):
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

        if(total_pax is None):
            pax2 = pax.twinx()
            ln5 = pax2.scatter(hash_match_SWF.compressed(), hash_mask_TF,\
                s = markersize, color='tab:cyan', marker='o',label = 'Smoke')
            ln6 = pax2.scatter(nohash_match_SWF.compressed(), nohash_mask_TF,\
                s = markersize, color='tab:purple', marker='o',label = 'Clear')
            print("Calculating MODIS/Total smoke trend")
            plot_trend_line(pax2, hash_match_SWF.compressed(), hash_mask_TF, \
                color='tab:cyan')
            print("Calculating MODIS/Total clear trend")
            plot_trend_line(pax2, nohash_match_SWF.compressed(), nohash_mask_TF, \
                color='tab:purple')
        else:
            ln5 = total_pax.scatter(hash_match_SWF.compressed(), hash_mask_TF,\
                s = markersize, color='tab:blue', marker='o',label = 'Smoke')
            ln6 = total_pax.scatter(nohash_match_SWF.compressed(), nohash_mask_TF,\
                s = markersize, color='tab:orange', marker='o',label = 'Clear')
            print("Calculating MODIS/Total smoke trend")
            plot_trend_line(total_pax, hash_match_SWF.compressed(), hash_mask_TF, \
                color='tab:blue')
            print("Calculating MODIS/Total clear trend")
            plot_trend_line(total_pax, nohash_match_SWF.compressed(), nohash_mask_TF, \
                color='tab:orange')


        #lns = lns + ln5 + ln6
        custom_lines = [Line2D([0], [0], color='tab:blue', label = 'SWF Smoke',),
                        Line2D([0], [0], color='tab:orange', label = 'SWF Clear',),
                        Line2D([0], [0], color='tab:green', label = 'LWF Smoke',),\
                        Line2D([0], [0], color='tab:red', label = 'LWF Clear'), \
                        Line2D([0], [0], color='tab:cyan', label = 'Total Smoke'),
                        Line2D([0], [0], color='tab:purple', label = 'Total Clear')]


    if(lw_pax is None):
        labs = [l.get_label() for l in custom_lines]
        pax.legend(custom_lines, labs, bbox_to_anchor = (1.05, 1.0), loc = 0, fontsize = 9)
    else:
        ##!#pax.legend(fontsize = 9, bbox_to_anchor = (1.05, 1.0), loc = 0)
        ##!#lw_pax.legend(fontsize = 9, bbox_to_anchor = (1.05, 1.0), loc = 0)
        total_pax.legend(fontsize = 9, bbox_to_anchor = (1.05, 1.0), loc = 0)
    ##!#pax.set_title('SWF smoke: '+str(np.round(shrval_p,3)) + '  SWF clear: ' + \
    ##!#     str(np.round(snrval_p,3)) + \
    ##!#    '\nLWF smoke: '+str(np.round(lhrval_p,3)) + '  LWF clear: ' + \
    ##!#    str(np.round(lnrval_p,3))) 
    pax.set_xlabel('MODIS ' + str(np.round(np.average(channel_dict[\
        str(MODIS_data['channel'])]['Bandwidth']),1)) + ' μm T$_{B}$', \
        fontsize = labelsize, weight = 'bold')
    if(lw_pax is None):
        pax.set_ylabel('TOA LW Flux [Wm$^{-2}$]', fontsize = labelsize, weight = 'bold')
    else:
        pax.set_ylabel('CERES SW Flux [Wm$^{-2}$]', fontsize = labelsize, weight = 'bold')
        lw_pax.set_xlabel('MODIS ' + str(np.round(np.average(channel_dict[\
            str(MODIS_data['channel'])]['Bandwidth']),1)) + ' μm T$_{B}$', \
            fontsize = labelsize, weight = 'bold')
        lw_pax.set_ylabel('CERES LW Flux [Wm$^{-2}$]', fontsize = labelsize, weight = 'bold')
        total_pax.set_xlabel('MODIS ' + str(np.round(np.average(channel_dict[\
            str(MODIS_data['channel'])]['Bandwidth']),1)) + ' μm T$_{B}$', \
            fontsize = labelsize, weight = 'bold')
        total_pax.set_ylabel('CERES Total Flux [Wm$^{-2}$]', fontsize = labelsize, weight = 'bold')

def plot_scatter_OMI_CERES(date_str, MODIS_data, pax, avg_pixel = False,\
        plume_only = True, xlabel = None, ylabel = None, ptitle = None,\
        labelsize = 14, labelticksize = 14):

    # Determine where the smoke is located
    # ------------------------------------
    hash_data, nohash_data = find_plume(date_str) 

    tmp_data = np.copy(MODIS_data['data'])
    tmp_lat  = np.copy(MODIS_data['lat'])
    tmp_lon  = np.copy(MODIS_data['lon'])

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
        # Remove the nans from the MODIS data, lats, and lons for both the
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
    mask_LAT, mask_LON, mask_swf, mask_lwf = read_CERES_match_MODIS(date_str)
    
    # Calculate the corners of the CERES pixels
    crnr_LAT = getCorners_1d(mask_LAT.data) ; crnr_LON = getCorners_1d(mask_LON)

    LAT, LON, mask_UVAI = read_OMI_match_MODIS(date_str)
    #LAT, LON, mask_UVAI, crnr_LAT, crnr_LON = \
    #    read_OMI_match_MODIS(date_str, corners = True)

    #if(avg_pixel):
   
    # Loop over the lower resolution data (CERES) and find the 
    # OMI pixel that is closest to each CERES pixel 
    # ----------------------------------------------------------

    # total arrays are dimensioned to the size of the CERES
    # arrays. Contain all matched MODIS data (in and out)
    total_omi_match_SWF    = np.full(mask_swf.shape, np.nan)
    total_omi_match_LWF    = np.full(mask_lwf.shape, np.nan)
    total_omi_match_LAT    = np.full(mask_LAT.shape, np.nan)
    total_omi_match_LON    = np.full(mask_LON.shape, np.nan)
    ##!## hash arrays are dimensioned to the size of the CERES
    ##!## arrays. Contain only in-plume matched MODIS data
    ##!#hash_match_SWF    = np.full(mask_swf.shape, np.nan)
    ##!#hash_match_LWF    = np.full(mask_lwf.shape, np.nan)
    ##!#hash_match_LAT    = np.full(mask_LAT.shape, np.nan)
    ##!#hash_match_LON    = np.full(mask_LON.shape, np.nan)
    ##!## nohash arrays are dimensioned to the size of the CERES
    ##!## arrays. Contain only out=of-plume matched MODIS data
    ##!#nohash_match_SWF    = np.full(mask_swf.shape, np.nan)
    ##!#nohash_match_LWF    = np.full(mask_lwf.shape, np.nan)
    ##!#nohash_match_LAT    = np.full(mask_LAT.shape, np.nan)
    ##!#nohash_match_LON    = np.full(mask_LON.shape, np.nan)

    if(debug):
        print(hash_plot_lat.shape)

    # Loop over the lower-resolution CERES data.    
    # Either grab the nearest MODIS pixel to the CERES pixel
    # OR average the MODIS pixels within a certain range of the
    # CERES pixel.
    # ---------------------------------------------------------
    for ii in range(mask_swf.shape[0]):

        #lat_list = sorted([crnr_LAT[ii], crnr_LAT[ii+1]])
        #lon_list = sorted([crnr_LON[ii], crnr_LON[ii+1]])

        # NEW APPROACH: use calculated pixel corners and use where 
        # statement to find all MODIS pixels within the CERES bounds
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
        # Find the gridpoint in the MODIS lat/lon data that 
        # corresponds to the AIRS pixel
        # ---------------------------------------------------- 
        ho_idx = nearest_gridpoint(mask_LAT[ii], mask_LON[ii],\
            LAT, LON)

        if(len(ho_idx[0]) > 1):
            ho_idx = (np.array([ho_idx[0][0]])), \
                (np.array([ho_idx[1][0]]))
   
        if(len(mask_UVAI[ho_idx]) > 0):
            ##!## if any smoke-classified MODIS pixel is within the CERES pixel,
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

            # Average the n MODIS pixels around the current CERES
            # pixel. 
            # total_match arrays contain all colocated MODIS data
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
    ##!#    np.mean(aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['Lon']),\
    ##!#    central_latitude = \
    ##!#np.mean(aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['Lat']))

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
    print("Calculating OMI/Total trend")
    plot_trend_line(pax, total_omi_match_SWF.compressed(), total_net, \
        color='tab:olive')

    if(xlabel is None):
        xlabel = 'OMI AI'
    if(ylabel is None):
        ylabel = 'TOA Flux [Wm$^{-2}$]'
    pax.set_xlabel(xlabel, fontsize = labelsize, weight = 'bold')
    pax.set_ylabel(ylabel, fontsize = labelsize, weight = 'bold')

    ##!#lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
    ##!#lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    ##!#if(compare_OMI):
    ##!#    plt.legend(lines, labels, loc = 'lower center', bbox_to_anchor = (0, 0.01, 1, 1),\
    ##!#        bbox_transform = plt.gcf().transFigure, ncol=2)
    ##!#else:
    ##!#    plt.legend(lines, labels, loc = 'right', bbox_to_anchor = (0, 0., 1.0, 1),\
    ##!#        bbox_transform = plt.gcf().transFigure, ncol=1)

    if(labelsize is not None):
        pax.legend(bbox_to_anchor = (1.05, 1.0), loc = 'upper left', fontsize = labelsize - 2)
    else:
        pax.legend(bbox_to_anchor = (1.05, 1.0), loc = 'upper left')
    #pax.set_title('Smoke correlation: '+str(np.round(lrval_p, 3)))
    if(ptitle is None):
        ptitle = 'OMI/SWF: '+str(np.round(srval_p,3)) + '  OMI/LWF: ' + \
        str(np.round(lrval_p,3)) + '\nOMI/Total:  ' + str(np.round(trval_p,3))
    pax.set_title(ptitle)

def plot_OMI_spatial(date_str, LAT, LON, mask_UVAI, pax, labelsize = 13, \
        labelticksize = 11, ptitle = '', zoom = False):

    dt_date_str = datetime.strptime(date_str, "%Y%m%d%H%M")

    mesh3 = pax.pcolormesh(LON,LAT,mask_UVAI, cmap = 'plasma', shading='auto', \
        vmin = np.nanmin(mask_UVAI), vmax = np.nanmax(mask_UVAI), transform = datacrs) 
    pax.add_feature(cfeature.BORDERS)
    pax.add_feature(cfeature.STATES)
    pax.coastlines()
    if(zoom):
        pax.set_extent([aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][1], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][1]],\
                        datacrs)
    cbar3 = plt.colorbar(mesh3,ax=pax,orientation='vertical',\
        pad=0.03, fraction = 0.046)
    cbar3.set_label('OMI UVAI', size = labelsize, weight = 'bold')
    cbar3.ax.tick_params(labelsize = labelticksize)
    pax.set_title(ptitle,fontsize=10)


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
    #    np.mean(aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['Lon']),\
    #    central_latitude = \
    #    np.mean(aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][dt_date_str.strftime('%H%M')]['Lat']))

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

    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(dt_date_str.strftime('%Y%m%d%H%M')) 

        ax0.pcolor(MODIS_data1['lon'],MODIS_data1['lat'],\
            hash_data1, hatch = '////', alpha=0., transform = datacrs)
        ax1.pcolor(MODIS_data1['lon'],MODIS_data1['lat'],\
            hash_data1, hatch = '////', alpha=0., transform = datacrs)
        ax2.pcolor(MODIS_data1['lon'],MODIS_data1['lat'],\
            hash_data1, hatch = '////', alpha=0., transform = datacrs)
        ax3.pcolor(MODIS_data1['lon'],MODIS_data1['lat'],\
            hash_data1, hatch = '////', alpha=0., transform = datacrs)
        ax4.pcolor(MODIS_data1['lon'],MODIS_data1['lat'],\
            hash_data1, hatch = '////', alpha=0., transform = datacrs)
        ax5.pcolor(MODIS_data1['lon'],MODIS_data1['lat'],\
            hash_data1, hatch = '////', alpha=0., transform = datacrs)
        ax6.pcolor(MODIS_data1['lon'],MODIS_data1['lat'],\
            hash_data1, hatch = '////', alpha=0., transform = datacrs)
        ax7.pcolor(MODIS_data1['lon'],MODIS_data1['lat'],\
            hash_data1, hatch = '////', alpha=0., transform = datacrs)

    fig.tight_layout()

    if(save):
        outname = 'modis_combined_imagery_' + date_str + zoom_add + '.png'
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
    # Read the MODIS data
    #
    # ----------------------------------------------------------------------

    # Call read_MODIS_channel to read the desired MODIS data from the file 
    # and put it in a dictionary
    # ---------------------------------------------------------------------
    MODIS_data0 = read_MODIS_channel(dt_date_str.strftime("%Y%m%d%H%M"), \
        channel0, zoom = zoom)
    MODIS_data1 = read_MODIS_channel(dt_date_str.strftime("%Y%m%d%H%M"), \
        channel1, zoom = zoom)
    MODIS_data2 = read_MODIS_channel(dt_date_str.strftime("%Y%m%d%H%M"), \
        channel2, zoom = zoom)
    MODIS_dataw = read_MODIS_channel(dt_date_str.strftime("%Y%m%d%H%M"), \
        'wv_ir', zoom = zoom)
    
    # ----------------------------------------------------------------------
    #
    # Read the CERES data
    #
    # ----------------------------------------------------------------------
    mask_LAT, mask_LON, mask_swf, mask_lwf = \
        read_CERES_match_MODIS(date_str)

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
        plot_scatter_OMI(date_str, MODIS_data0, ax4, avg_pixel = avg_pixel)

        # ----------------------------------------------------------------------
        #
        # Colocate the OMI and CERES data, then plot
        #
        # ----------------------------------------------------------------------
        plot_scatter_OMI_CERES(date_str, MODIS_data0, ax5, \
            avg_pixel = avg_pixel, plume_only = plume_only)
    
        plot_subplot_label(ax4, '(e)')
        plot_subplot_label(ax5, '(f)')

    # Determine where the smoke is located
    # ------------------------------------
    hash_data, nohash_data = find_plume(dt_date_str.strftime('%Y%m%d%H%M')) 

    tmp_data0 = np.copy(MODIS_data0['data'])
    tmp_data1 = np.copy(MODIS_data1['data'])
    tmp_data2 = np.copy(MODIS_data2['data'])
    tmp_dataw = np.copy(MODIS_dataw['data'])
    tmp_lat0  = np.copy(MODIS_data0['lat'])
    tmp_lon0  = np.copy(MODIS_data0['lon'])

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
    # Plot the MODIS data
    #
    # ----------------------------------------------------------------------

    plot_scatter(ax0, tmp_data0, tmp_data1, MODIS_data0, MODIS_data1, \
        hash_data)
    plot_scatter(ax1, tmp_data0, tmp_data2, MODIS_data0, MODIS_data2, \
        hash_data)
    plot_scatter(ax2, tmp_data0, tmp_dataw, MODIS_data0, MODIS_dataw, \
        hash_data)

    # ----------------------------------------------------------------------
    #
    # Plot the CERES data
    #
    # ----------------------------------------------------------------------
    plot_scatter_CERES(date_str, MODIS_data0, ax3, avg_pixel = avg_pixel,\
        plume_only = plume_only, plot_total_flux = False)

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
        outname = 'modis_combined_scatter_' + date_str + plume_add + \
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
    # Read the MODIS data
    #
    # ----------------------------------------------------------------------

    # Call read_MODIS_channel to read the desired MODIS data from the file 
    # and put it in a dictionary
    # ---------------------------------------------------------------------
    MODIS_dataJul31 = read_MODIS_channel(date_str_July,31, zoom = zoom)
    MODIS_dataAug31 = read_MODIS_channel(date_str_Aug, 31, zoom = zoom)

    ##!## Determine where the smoke is located
    ##!## ------------------------------------
    ##!#hash_data, nohash_data = find_plume(dt_date_str.strftime('%Y%m%d%H%M')) 

    ##!#tmp_data0 = np.copy(MODIS_data0['data'])
    ##!#tmp_data1 = np.copy(MODIS_data1['data'])
    ##!#tmp_data2 = np.copy(MODIS_data2['data'])
    ##!#tmp_dataw = np.copy(MODIS_dataw['data'])
    ##!#tmp_lat0  = np.copy(MODIS_data0['lat'])
    ##!#tmp_lon0  = np.copy(MODIS_data0['lon'])

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
        read_CERES_match_MODIS(date_str_July)
    mask_LAT_Aug, mask_LON_Aug, mask_swf_Aug, mask_lwf_Aug = \
        read_CERES_match_MODIS(date_str_Aug)

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
    plot_scatter_OMI(date_str_Aug, MODIS_dataAug31, ax2, avg_pixel = avg_pixel)

    # ----------------------------------------------------------------------
    #
    # Colocate the OMI and CERES data, then plot
    #
    # ----------------------------------------------------------------------
    plot_scatter_OMI_CERES(date_str_Aug, MODIS_dataAug31, ax3, \
        avg_pixel = avg_pixel, plume_only = plume_only)

    # ----------------------------------------------------------------------
    #
    # Plot the CERES data
    #
    # ----------------------------------------------------------------------
    plot_scatter_CERES(date_str_July, MODIS_dataJul31, ax0, \
        avg_pixel = avg_pixel,plume_only = plume_only)
    plot_scatter_CERES(date_str_Aug, MODIS_dataAug31, ax1, \
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
        outname = 'modis_combined_scatter_select' + plume_add + \
            pixel_add + '.png'
        fig.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

def plot_spatial_scatter(date_str, zoom=True,save=False,composite=True,\
        avg_pixel=False,plume_only=False):

    if(home_dir + '/Research/CERES' not in sys.path):
        sys.path.append(home_dir + '/Research/CERES')
    from gridCERESLib import readgridCERES_hrly_grid, plotCERES_hrly
    
    plt.close('all')
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # ----------------------------------------------------------------------
    #
    # Read the MODIS data
    #
    # ----------------------------------------------------------------------

    # Call read_MODIS_channel to read the desired MODIS data from the file 
    # and put it in a dictionary
    # ---------------------------------------------------------------------
    MODIS_data0 = read_MODIS_channel(dt_date_str.strftime("%Y%m%d%H%M"), \
        31, zoom = zoom)
    var0, crs0, lons0, lats0, lat_lims0, lon_lims0, plabel0 = \
        read_MODIS_satpy(date_str, 31)
    
    # ----------------------------------------------------------------------
    #
    # Read the CERES data
    #
    # ----------------------------------------------------------------------
    CERES_data_hrly_swf = readgridCERES_hrly_grid(date_str[:10], 'SWF', minlat = 20.)

    ##!#mask_LAT, mask_LON, mask_swf, mask_lwf = \
    ##!#    read_CERES_match_MODIS(date_str)

    # ------------------------------------------------------------------------
    #
    # For 20210722:
    #
    # ------------------------------------------------------------------------
    mapcrs = init_proj(date_str)
    if(date_str == '202108062025'):
        fig = plt.figure(figsize=(11,9))
    else:
        fig = plt.figure(figsize=(8.0,10))
    ax0 = fig.add_subplot(3,2,1, projection = mapcrs) # CERES SW
    ax1 = fig.add_subplot(3,2,3, projection = mapcrs) # CERES LW
    ax2 = fig.add_subplot(3,2,5, projection = mapcrs) # CERES total 
    #ax2 = fig.add_subplot(2,2,3, projection = mapcrs) # MODIS ch 31
    ax3 = fig.add_subplot(3,2,2)                      # 31 vs SW
    ax4 = fig.add_subplot(3,2,4)                      # 31 vs LW
    ax5 = fig.add_subplot(3,2,6)                      # 31 vs Total

    ##!## Determine where the smoke is located
    ##!## ------------------------------------
    ##!#hash_data, nohash_data = find_plume(dt_date_str.strftime('%Y%m%d%H%M')) 

    ##!#tmp_data0 = np.copy(MODIS_data0['data'])
    ##!#tmp_lat0  = np.copy(MODIS_data0['lat'])
    ##!#tmp_lon0  = np.copy(MODIS_data0['lon'])

    ##!#if(not (tmp_data0.shape ==  hash_data.shape)):
    ##!#    print("shape mismatch")
    ##!#    shapes = []
    ##!#    shapes.append(tmp_data0.shape)
    ##!#    shapes.append(hash_data.shape)

    ##!#    min_shape = min(shapes)
    ##!#    print(min_shape)

    ##!#    tmp_data0  = tmp_data0[:min_shape[0],:min_shape[1]]
    ##!#    hash_data  = hash_data[:min_shape[0],:min_shape[1]]

    ##!#max_ch = 350.

    ##!#tmp_data0 = np.ma.masked_where( (abs(tmp_data0) > max_ch), tmp_data0)
    ##!#tmp_lat0  = np.ma.masked_where( (abs(tmp_data0) > max_ch), tmp_lat0)
    ##!#tmp_lon0  = np.ma.masked_where( (abs(tmp_data0) > max_ch), tmp_lon0)

    # ----------------------------------------------------------------------
    #
    # Plot the MODIS data
    #
    # ----------------------------------------------------------------------
    ##!#plot_MODIS_spatial(MODIS_data0, ax2, zoom = zoom, ptitle = '')
    ##!#plot_MODIS_satpy(date_str, 31, ax = ax2, var = var0, crs = crs0, \
    ##!#    lons = lons0, lats = lats0, lat_lims = lat_lims0, lon_lims = lon_lims0, \
    ##!#    vmin = 270, vmax = 330, ptitle = '', plabel = plabel0, \
    ##!#    labelsize = 12, colorbar = True, zoom=True,save=False)

    # ----------------------------------------------------------------------
    #
    # Plot the CERES data
    #
    # ----------------------------------------------------------------------
    vmax = None
    vmin = None
    if(date_str == '202107222110'):
        sw_vmax = 250
        sw_vmin = 130
        lw_vmax = 370
        lw_vmin = 300
        tot_vmax = 580
        tot_vmin = 460
    elif(date_str == '202108062025'):
        sw_vmax = 330
        sw_vmin = 190
        lw_vmax = 370
        lw_vmin = 300
   
    plotCERES_hrly(ax0, CERES_data_hrly_swf, 'swf', \
        vmin = sw_vmin, vmax = sw_vmax, title = '', label = 'TOA Flux [W/m$^{2}$]', \
        circle_bound = False, gridlines = False, grid_data = True, \
        zoom = True)
    plotCERES_hrly(ax1, CERES_data_hrly_swf, 'lwf', \
        vmin = lw_vmin, vmax = lw_vmax, title = '', label = 'TOA Flux [W/m$^{2}$]', \
        circle_bound = False, gridlines = False, grid_data = True, \
        zoom = True)
    plotCERES_hrly(ax2, CERES_data_hrly_swf, 'total', \
        vmin = tot_vmin, vmax = tot_vmax, title = '', label = 'TOA Flux [W/m$^{2}$]', \
        circle_bound = False, gridlines = False, grid_data = True, \
        zoom = True)
    if(zoom):
        ax0.set_extent([aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][1], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][1]],\
                        datacrs)
        ax1.set_extent([aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][1], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][1]],\
                        datacrs)
        ax2.set_extent([aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][1], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][1]],\
                        datacrs)

    plot_scatter_CERES(date_str, MODIS_data0, ax3, avg_pixel = avg_pixel,\
        plume_only = plume_only,  plot_total_flux = True, lw_pax = ax4, \
        total_pax = ax5, labelsize = 11)

    ax0.set_extent([lon_lims0[0],lon_lims0[1],lat_lims0[0],lat_lims0[1]],\
                   crs = ccrs.PlateCarree())
    ax1.set_extent([lon_lims0[0],lon_lims0[1],lat_lims0[0],lat_lims0[1]],\
                   crs = ccrs.PlateCarree())
    ax2.set_extent([lon_lims0[0],lon_lims0[1],lat_lims0[0],lat_lims0[1]],\
                   crs = ccrs.PlateCarree())

    # Add labels to all the subplots
    # ------------------------------
    plot_subplot_label(ax0, '(a)', backgroundcolor = 'white')
    plot_subplot_label(ax1, '(c)', backgroundcolor = 'white')
    plot_subplot_label(ax3, '(b)')
    plot_subplot_label(ax4, '(d)')
    plot_subplot_label(ax2, '(e)', backgroundcolor = 'white')
    plot_subplot_label(ax5, '(f)')

    font_size = 13
    plot_figure_text(ax0, 'CERES SW', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax1, 'CERES LW', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax2, 'CERES Total', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size - 2, backgroundcolor = 'white', halign = 'right')

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
        outname = 'modis_spatial_scatter_' + date_str + plume_add + \
            pixel_add + '_v3.png'
            #pixel_add + '.png'
        fig.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

def plot_ceres_scatter(date_str, zoom=True,save=False,composite=True,\
        avg_pixel=False,plume_only=False):

    if(home_dir + '/Research/CERES' not in sys.path):
        sys.path.append(home_dir + '/Research/CERES')
    from gridCERESLib import readgridCERES_hrly_grid, plotCERES_hrly
    
    plt.close('all')
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # ----------------------------------------------------------------------
    #
    # Read the MODIS data
    #
    # ----------------------------------------------------------------------

    # Call read_MODIS_channel to read the desired MODIS data from the file 
    # and put it in a dictionary
    # ---------------------------------------------------------------------
    MODIS_data0 = read_MODIS_channel(dt_date_str.strftime("%Y%m%d%H%M"), \
        31, zoom = zoom)
    var0, crs0, lons0, lats0, lat_lims0, lon_lims0, plabel0 = \
        read_MODIS_satpy(date_str, 31)
    
    # ----------------------------------------------------------------------
    #
    # Read the CERES data
    #
    # ----------------------------------------------------------------------
    CERES_data_hrly_swf = readgridCERES_hrly_grid(date_str[:10], 'SWF', minlat = 20.)

    # ------------------------------------------------------------------------
    #
    # For 20210722:
    #
    # ------------------------------------------------------------------------
    mapcrs = init_proj(date_str)
    fig = plt.figure(figsize = (11.5, 4))
    ax3 = fig.add_subplot(1,3,1)                      # 31 vs SW
    ax4 = fig.add_subplot(1,3,2)                      # 31 vs LW
    ax5 = fig.add_subplot(1,3,3)                      # 31 vs Total

    plot_scatter_CERES(date_str, MODIS_data0, ax3, avg_pixel = avg_pixel,\
        plume_only = plume_only,  plot_total_flux = True, lw_pax = ax4, \
        total_pax = ax5, labelsize = 11)

    # Add labels to all the subplots
    # ------------------------------
    plot_subplot_label(ax3, '(a)')
    plot_subplot_label(ax4, '(b)')
    plot_subplot_label(ax5, '(c)')

    plt.suptitle(dt_date_str.strftime(\
        'Aqua CERES vs. MODIS 11 μm Brightness Temperature\n%d %B %Y %H:%M UTC'))

    fig.tight_layout()

    if(save):
        plume_add = '_onlyplume'
        pixel_add = ''
        if(plume_only == False):
            plume_add = '_onlyplume'
        if(avg_pixel):
            pixel_add = '_avgpixel' 
        outname = 'modis_ceres_scatter_' + date_str + plume_add + \
            pixel_add + '_v4.png'
            #pixel_add + '.png'
        fig.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()


def plot_spatial_scatter_wAI(date_str, zoom=True,save=False,composite=True,\
        avg_pixel=False,plume_only=False):

    if(home_dir + '/Research/CERES' not in sys.path):
        sys.path.append(home_dir + '/Research/CERES')
    from gridCERESLib import readgridCERES_hrly_grid, plotCERES_hrly
    
    plt.close('all')
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # ----------------------------------------------------------------------
    #
    # Read the MODIS data
    #
    # ----------------------------------------------------------------------
    # Read the true color data
    # ------------------------
    var1, crs1, lat_lims1, lon_lims1 = read_true_color(date_str,\
        composite=composite)

    # Call read_MODIS_channel to read the desired MODIS data from the file 
    # and put it in a dictionary
    # ---------------------------------------------------------------------
    MODIS_data0 = read_MODIS_channel(dt_date_str.strftime("%Y%m%d%H%M"), \
        31, zoom = zoom)
    
    # ----------------------------------------------------------------------
    #
    # Read the CERES data
    #
    # ----------------------------------------------------------------------
    CERES_data_hrly_swf = readgridCERES_hrly_grid(date_str[:10], 'SWF')

    # ---------------------------------------------------------------------
    #
    # Read the OMI data
    #
    # ---------------------------------------------------------------------
    LAT, LON, mask_UVAI = read_OMI_match_MODIS(date_str)
 

    # ------------------------------------------------------------------------
    #
    # For 20210722:
    #
    # ------------------------------------------------------------------------
    plt.close('all')
    mapcrs = init_proj(date_str)
    if(date_str == '202108062025'):
        fig  = plt.figure(figsize=(9.5,4.5))
        fig2 = plt.figure(figsize=(5.5,7))
    else:
        fig = plt.figure(figsize=(9,9))
    ax1  = fig.add_subplot(2,3,1,projection = crs1)   # true color    
    ax2  = fig.add_subplot(2,3,2,projection = mapcrs) # MODIS Ch 31
    ax3  = fig.add_subplot(2,3,3,projection = mapcrs) # CERES SWF
    ax4  = fig.add_subplot(2,3,4,projection = mapcrs) # CERES LWF
    ax5  = fig.add_subplot(2,3,5,projection = mapcrs) # OMI AI
    
    ax6  = fig2.add_subplot(2,1,1)                     # IR vs AI 
    ax7  = fig2.add_subplot(2,1,2)                     # IR vs CERES
    #ax8  = fig2.add_subplot(3,1,3)                     # AI vs CERES
    ##!#gs = fig.add_gridspec(nrows = 8, ncols = 16)
    ##!#ax1  = fig.add_subplot(gs[0:4  ,0:4],projection = crs1)   # true color    
    ##!#ax2  = fig.add_subplot(gs[0:4  ,4:8],projection = mapcrs) # MODIS Ch 31
    ##!#ax3  = fig.add_subplot(gs[0:4  ,8:12],projection = mapcrs) # CERES SWF
    ##!#ax4  = fig.add_subplot(gs[0:4  ,12:16],projection = mapcrs) # CERES LWF
    ##!#ax5  = fig.add_subplot(gs[4:8  ,1:5],projection = mapcrs) # OMI AI
    ##!#ax6  = fig.add_subplot(gs[4:8  ,5:10])                     # IR vs CERES
    ##!#ax7  = fig.add_subplot(gs[4:8  ,10:15])                     # AI vs CERES

    # ----------------------------------------------------------------------
    #
    # Plot the MODIS data
    #
    # ----------------------------------------------------------------------
    # Plot the true-color data for the previous date
    # ----------------------------------------------
    ax1.imshow(var1.data, transform = crs1, extent=(var1.x[0], var1.x[-1], \
        var1.y[-1], var1.y[0]), origin='upper')

    # Zoom in the figure if desired
    # -----------------------------
    if(zoom):
        ax1.set_extent([lon_lims1[0],lon_lims1[1],lat_lims1[0],\
            lat_lims1[1]],crs = datacrs)

    # Plot the MODIS IR data
    # ----------------------
    plot_MODIS_spatial(MODIS_data0, ax2, zoom = zoom, ptitle = '', \
        labelsize = None, labelticksize = None)

    # ----------------------------------------------------------------------
    #
    # Plot the CERES data
    #
    # ----------------------------------------------------------------------
    vmax = None
    vmin = None
    if(date_str == '202107222110'):
        sw_vmax = 250
        sw_vmin = 130
        lw_vmax = 370
        lw_vmin = 300
    elif(date_str == '202108062025'):
        sw_vmax = 330
        sw_vmin = 190
        lw_vmax = 370
        lw_vmin = 300
    
    plotCERES_hrly(ax3, CERES_data_hrly_swf, 'swf', \
        vmin = sw_vmin, vmax = sw_vmax, title = '', label = 'TOA Flux [Wm$^{-2}$]', \
        circle_bound = False, gridlines = False, grid_data = True, \
        labelsize = None, labelticksize = None, zoom = True)
    plotCERES_hrly(ax4, CERES_data_hrly_swf, 'lwf', \
        vmin = lw_vmin, vmax = lw_vmax, title = '', label = 'TOA Flux [Wm$^{-2}$]', \
        circle_bound = False, gridlines = False, grid_data = True, \
        labelsize = None, labelticksize = None, zoom = True)
    if(zoom):
        ax3.set_extent([aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][1], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][1]],\
                        datacrs)
        ax4.set_extent([aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][1], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][1]],\
                        datacrs)

    ##!#plot_scatter_CERES(date_str, MODIS_data0, ax6, avg_pixel = avg_pixel,\
    ##!#    plume_only = plume_only,  plot_total_flux = False, labelsize = None)

    # ----------------------------------------------------------------------
    #
    # Plot the OMI data
    #
    # ----------------------------------------------------------------------
    plot_OMI_spatial(date_str, LAT, LON, mask_UVAI, ax5, labelsize = None, \
        labelticksize = None, zoom = zoom)

    plot_scatter_OMI(date_str, MODIS_data0, ax6, avg_pixel = avg_pixel, \
        xlabel = 'MODIS 11.0 μm T$_{B}$', \
        ylabel = 'OMI UVAI', ptitle = '', labelsize = None)

    plot_scatter_OMI_CERES(date_str, MODIS_data0, ax7, \
        avg_pixel = avg_pixel, plume_only = plume_only, ptitle = '', labelsize = None)

    # Add labels to all the subplots
    # ------------------------------
    plot_subplot_label(ax1, '(a)', fontsize = 11, backgroundcolor = 'white')
    plot_subplot_label(ax2, '(b)', fontsize = 11, backgroundcolor = 'white')
    plot_subplot_label(ax3, '(c)', fontsize = 11, backgroundcolor = 'white')
    plot_subplot_label(ax4, '(d)', fontsize = 11, backgroundcolor = 'white')
    plot_subplot_label(ax5, '(e)', fontsize = 11, backgroundcolor = 'white')
    plot_subplot_label(ax6, '(a)', fontsize = 12, location = 'lower_left')
    plot_subplot_label(ax7, '(b)', fontsize = 12, )
    #plot_subplot_label(ax8, '(c)', fontsize = 12, )

    font_size = 11
    plot_figure_text(ax1, 'MODIS True Color', xval = None, yval = None, \
        transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    plot_figure_text(ax2, 'MODIS 11.0 μm', xval = None, yval = None, \
        transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    plot_figure_text(ax3, 'CERES SW', xval = None, yval = None, \
        transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    plot_figure_text(ax4, 'CERES LW', xval = None, yval = None, \
        transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')
    plot_figure_text(ax5, 'OMI AI', xval = None, yval = None, \
        transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')

    fig.tight_layout()
    fig2.tight_layout()

    if(save):
        plume_add = '_onlyplume'
        pixel_add = ''
        if(plume_only == False):
            plume_add = '_onlyplume'
        if(avg_pixel):
            pixel_add = '_avgpixel' 
        outname = 'modis_spatial_wAI_' + date_str + plume_add + \
            pixel_add + '.png'
        fig.savefig(outname,dpi=300)
        print("Saved image",outname)

        outname = 'modis_scatter_wAI_' + date_str + plume_add + \
            pixel_add + '_v2.png'
        fig2.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Figure plotting functions for the smoke thermal IR paper
# (Sorenson et al., 2023b)
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


def plot_combined_figure1(date_str = '202107222110', zoom = True, \
        show_smoke = True, composite = True, save=False):

    if(home_dir + '/Research/CERES' not in sys.path):
        sys.path.append(home_dir + '/Research/CERES')
    from gridCERESLib import readgridCERES_hrly_grid, plotCERES_hrly

    #date_str = '202107202125'
    #date_str = '202107222110'
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # ----------------------------------------------------------------------
    #
    # Read the MODIS and CERES data
    #
    # ----------------------------------------------------------------------

    # Call read_MODIS_channel to read the desired MODIS data from the file 
    # and put it in a dictionary
    # ---------------------------------------------------------------------
    MODIS_data_ch1  = read_MODIS_channel(date_str, 1, zoom = zoom)
    MODIS_data_ch7  = read_MODIS_channel(date_str, 7, zoom = zoom)
    MODIS_data_ch31 = read_MODIS_channel(date_str, 31, zoom = zoom)
    MODIS_data_wv   = read_MODIS_channel(date_str, 'wv_ir', zoom = zoom)

    # Determine where the smoke is located
    # ------------------------------------
    hash_data, nohash_data = find_plume(dt_date_str.strftime('%Y%m%d%H%M')) 

    tmp_data1  = np.copy(MODIS_data_ch1['data'])
    tmp_data7  = np.copy(MODIS_data_ch7['data'])
    tmp_data31 = np.copy(MODIS_data_ch31['data'])
    tmp_datawv = np.copy(MODIS_data_wv['data'])
    tmp_lat0   = np.copy(MODIS_data_ch1['lat'])
    tmp_lon0   = np.copy(MODIS_data_ch1['lon'])

    if(not (tmp_data1.shape == tmp_data7.shape == tmp_data31.shape == \
            tmp_datawv.shape == hash_data.shape)):
        print("shape mismatch")
        shapes = []
        shapes.append(tmp_data1.shape)
        shapes.append(tmp_data7.shape)
        shapes.append(tmp_data31.shape)
        shapes.append(tmp_datawv.shape)
        shapes.append(hash_data.shape)

        min_shape = min(shapes)
        print(min_shape)

        tmp_data1  = tmp_data1[:min_shape[0],:min_shape[1]]
        tmp_data7  = tmp_data7[:min_shape[0],:min_shape[1]]
        tmp_data31 = tmp_data31[:min_shape[0],:min_shape[1]]
        tmp_datawv = tmp_datawv[:min_shape[0],:min_shape[1]]
        tmp_lat0   = tmp_lat0[:min_shape[0],:min_shape[1]]
        tmp_lon0   = tmp_lon0[:min_shape[0],:min_shape[1]]
        hash_data  = hash_data[:min_shape[0],:min_shape[1]]

    max_ch = 350.

    tmp_data1 = np.ma.masked_where( (abs(tmp_data7) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) | \
        (abs(tmp_datawv) > max_ch), tmp_data1)
    tmp_data7 = np.ma.masked_where( (abs(tmp_data7) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) | \
        (abs(tmp_datawv) > max_ch), tmp_data7)
    tmp_data31 = np.ma.masked_where( (abs(tmp_data7) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) | \
        (abs(tmp_datawv) > max_ch), tmp_data31)
    tmp_datawv = np.ma.masked_where( (abs(tmp_data7) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) | \
        (abs(tmp_datawv) > max_ch), tmp_datawv)
    tmp_lat0 = np.ma.masked_where( (abs(tmp_data7) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) | \
        (abs(tmp_datawv) > max_ch), tmp_lat0)
    tmp_lon0 = np.ma.masked_where( (abs(tmp_data7) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch) | \
        (abs(tmp_datawv) > max_ch), tmp_lon0)

    # Read the true color data
    # ------------------------
    var1, crs1, lat_lims1, lon_lims1 = read_true_color(date_str,\
        composite=composite)

    # Read in the CERES data
    # ----------------------
    ##!#mask_LAT, mask_LON, mask_swf, mask_lwf = \
    ##!#    read_CERES_match_MODIS(date_str)
    CERES_data_hrly_swf = readgridCERES_hrly_grid(date_str[:10], 'SWF')

    if(CERES_data_hrly_swf is None):
        print("ERROR: no data returned from readgridCERES_hrly_grid")
        print("Quitting")
        return
    
   
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
    plot_MODIS_spatial(MODIS_data_ch31, ax2, zoom = zoom, ptitle = '')
    plot_MODIS_spatial(MODIS_data_ch1,  ax3, zoom = zoom, ptitle = '')
    plot_MODIS_spatial(MODIS_data_ch7,  ax4, zoom = zoom, ptitle = '')
    plot_MODIS_spatial(MODIS_data_wv,   ax7, zoom = zoom, ptitle = '')

    # Plot the MODIS scatter data
    # ---------------------------
    plot_scatter(ax5, tmp_data31, tmp_data1, MODIS_data_ch31, MODIS_data_ch1, \
        hash_data, xlabel = '11 μm brightness temperature', \
        ylabel = '0.64 μm reflectance', plot_legend = True)
    plot_scatter(ax6, tmp_data31, tmp_datawv, MODIS_data_ch31, MODIS_data_wv, \
        hash_data, xlabel = '11 μm brightness temperature', \
        ylabel = '2.1 μm reflectance')

    # Plot the CERES SWF and LWF data
    # -------------------------------
    ##!#plot_CERES_spatial(date_str, mask_LAT, mask_LON, mask_swf, 'SWF', ax8, \
    ##!#    ptitle = '', zoom = zoom)
    ##!#plot_CERES_spatial(date_str, mask_LAT, mask_LON, mask_lwf, 'LWF', ax9, \
    ##!#    ptitle = '', zoom = zoom)
    plotCERES_hrly(ax8, CERES_data_hrly_swf, 'swf', \
        vmin = 130, vmax = 250, title = '', label = '', \
        circle_bound = False, gridlines = False, grid_data = True, \
        zoom = True)
    if(zoom):
        ax8.set_extent([aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
                        date_str[8:]]['Lon'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
                        date_str[8:]]['Lon'][1], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
                        date_str[8:]]['Lat'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
                        date_str[8:]]['Lat'][1]],\
                        datacrs)

    plotCERES_hrly(ax9, CERES_data_hrly_swf, 'lwf', \
        vmin = None, vmax = None, title = '', label = '', \
        circle_bound = False, gridlines = False, grid_data = True, \
        zoom = True)
    if(zoom):
        ax9.set_extent([aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
                        date_str[8:]]['Lon'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
                        date_str[8:]]['Lon'][1], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
                        date_str[8:]]['Lat'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
                        date_str[8:]]['Lat'][1]],\
                        datacrs)


    #plotCERES_hrly(ax8, CERES_data_hrly, minlat=65, \
    #    vmin = None, vmax = None, title = '', label = '', \
    #    circle_bound = True, gridlines = True, grid_data = True, \
    #    zoom = True):

    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(date_str) 

        plt.rcParams.update({'hatch.color': 'r'})
        ax3.pcolor(MODIS_data_ch1['lon'],MODIS_data_ch1['lat'],\
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
    plot_figure_text(ax2, 'MODIS 11 μm', xval = None, yval = None, \
        transform = None, color = 'red', fontsize = 12, \
        backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax3, 'MODIS 0.64 μm', xval = None, yval = None, \
        transform = None, color = 'red', fontsize = 12, \
        backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax4, 'MODIS 2.1 μm', xval = None, yval = None, \
        transform = None, color = 'red', fontsize = 12, \
        backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax7, 'MODIS IR WV', xval = None, yval = None, \
        transform = None, color = 'red', fontsize = 12, \
        backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax8, 'CERES SW', xval = None, yval = None, \
        transform = None, color = 'red', fontsize = 12, \
        backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax9, 'CERES LW', xval = None, yval = None, \
        transform = None, color = 'red', fontsize = 12, \
        backgroundcolor = 'white', halign = 'right')

    fig.tight_layout()

    MODIS_data_ch1.clear()
    MODIS_data_ch7.clear()
    MODIS_data_ch31.clear()
    MODIS_data_wv.clear()

    if(save):
        outname = 'modis_total_combined_' + date_str + '_fig1.png'
        fig.savefig(outname, dpi=300)
        print("Saved",outname)
    else:
        plt.show()

def plot_combined_figure1_v2(date_str = '202107202125', zoom = True, show_smoke = True, composite = True, \
        save=False):

    if(home_dir + '/Research/CERES' not in sys.path):
        sys.path.append(home_dir + '/Research/CERES')
    from gridCERESLib import readgridCERES_hrly_grid, plotCERES_hrly
    if(home_dir + '/Research/GOES' not in sys.path):
        sys.path.append(home_dir + '/Research/GOES')
    from GOESLib import read_GOES_satpy, plot_GOES_satpy

    #date_str = '202107202125'
    #date_str = '202107222110'
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # ----------------------------------------------------------------------
    #
    # Read the MODIS and CERES data
    #
    # ----------------------------------------------------------------------

    # Call read_MODIS_channel to read the desired MODIS data from the file 
    # and put it in a dictionary
    # ---------------------------------------------------------------------
    MODIS_data_ch1  = read_MODIS_channel(date_str, 1, zoom = zoom)
    MODIS_data_ch5  = read_MODIS_channel(date_str, 7, zoom = zoom)
    MODIS_data_ch31 = read_MODIS_channel(date_str, 31, zoom = zoom)

    # Determine where the smoke is located
    # ------------------------------------
    hash_data, nohash_data = find_plume(dt_date_str.strftime('%Y%m%d%H%M')) 

    tmp_data1  = np.copy(MODIS_data_ch1['data'])
    tmp_data5  = np.copy(MODIS_data_ch5['data'])
    tmp_data31 = np.copy(MODIS_data_ch31['data'])
    tmp_lat0   = np.copy(MODIS_data_ch1['lat'])
    tmp_lon0   = np.copy(MODIS_data_ch1['lon'])

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
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch), tmp_data1)
    tmp_data5 = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch), tmp_data5)
    tmp_data31 = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch), tmp_data31)
    tmp_lat0 = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch), tmp_lat0)
    tmp_lon0 = np.ma.masked_where( (abs(tmp_data5) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch), tmp_lon0)

    # Read the true color data
    # ------------------------
    var1, crs1, lat_lims1, lon_lims1 = read_true_color(date_str,\
        composite=composite)

    # Read the GOES data
    # ------------------------
    var2, crs0, lons, lats, lat_lims2, lon_lims2, plabel = read_GOES_satpy(date_str, 8)
    var3, crs0, lons, lats, lat_lims0, lon_lims0, plabel = read_GOES_satpy(date_str, 9)
    var4, crs0, lons, lats, lat_lims0, lon_lims0, plabel = read_GOES_satpy(date_str, 10)
    var0, crs0, lons, lats, lat_lims0, lon_lims0, plabel = read_GOES_satpy(date_str, 13)

    # Read in the CERES data
    # ----------------------
    ##!#mask_LAT, mask_LON, mask_swf, mask_lwf = \
    ##!#    read_CERES_match_MODIS(date_str)
    CERES_data_hrly_swf = readgridCERES_hrly_grid(date_str[:10], 'SWF')

    if(CERES_data_hrly_swf is None):
        print("ERROR: no data returned from readgridCERES_hrly_grid")
        print("Quitting")
        return
    
   
    # ----------------------------------------------------------------------
    #
    #  Set up the 11-panel figure
    #
    # ----------------------------------------------------------------------

    mosaic = \
        """
        ABCD
        ABCD
        AEFG
        HEFG
        HIJK
        HIJK
        """ 
    mapcrs = init_proj(date_str)
    plt.close('all')
    fig = plt.figure(figsize=(11,7))
    gs = fig.add_gridspec(nrows = 6, ncols = 4)
    ax1  = fig.add_subplot(gs[0:3,0],projection = crs1)   # true color    
    ax2  = fig.add_subplot(gs[0:2,1],projection = mapcrs) # Ch 1
    ax3  = fig.add_subplot(gs[0:2,2],projection = mapcrs) # Ch 5
    ax4  = fig.add_subplot(gs[0:2,3],projection = mapcrs) # Ch 31
    ax5  = fig.add_subplot(gs[2:4,1],projection = mapcrs) # CERES SWF
    ax6  = fig.add_subplot(gs[2:4,2],projection = mapcrs) # CERES LWF
    ax7  = fig.add_subplot(gs[2:4,3],projection = mapcrs) # CERES total
    ax8  = fig.add_subplot(gs[3:6,0],projection = crs0)   # GOES 11 μm
    ax9  = fig.add_subplot(gs[4:6,1],projection = crs0)   # GOES lower WV
    ax10 = fig.add_subplot(gs[4:6,2],projection = crs0)   # GOES mid WF
    ax11 = fig.add_subplot(gs[4:6,3],projection = crs0)   # GOES upper WV
    ##!#ax1 = fig.add_subplot(3,3,1,projection = crs1)   # true color    
    ##!#ax2 = fig.add_subplot(3,3,2,projection = mapcrs) # Ch 31
    ##!#ax3 = fig.add_subplot(3,3,3,projection = mapcrs) # Ch 1
    ##!#ax4 = fig.add_subplot(3,3,4,projection = mapcrs) # Ch 5
    ##!#ax5 = fig.add_subplot(3,3,5)                     # Ch 31 vs 1 scatter
    ##!#ax6 = fig.add_subplot(3,3,6)                     # Ch 31 vs WV scatter
    ##!#ax7 = fig.add_subplot(3,3,7,projection = mapcrs) # WV
    ##!#ax8 = fig.add_subplot(3,3,8,projection = mapcrs) # LW
    ##!#ax9 = fig.add_subplot(3,3,9,projection = mapcrs) # SW
    var2, crs0, lons, lats, lat_lims2, lon_lims2, plabel2 = read_GOES_satpy(date_str, 8)
    var3, crs0, lons, lats, lat_lims0, lon_lims0, plabel3 = read_GOES_satpy(date_str, 9)
    var4, crs0, lons, lats, lat_lims0, lon_lims0, plabel4 = read_GOES_satpy(date_str, 10)
    var0, crs0, lons, lats, lat_lims0, lon_lims0, plabel0 = read_GOES_satpy(date_str, 13)

    # Plot the true-color data for the previous date
    # ----------------------------------------------
    ax1.imshow(var1.data, transform = crs1, extent=(var1.x[0], var1.x[-1], \
        var1.y[-1], var1.y[0]), origin='upper')

    labelsize = 8
    plot_GOES_satpy(date_str, 13, ax = ax8, var = var0, crs = crs0, \
        lat_lims = lat_lims2, lon_lims = lon_lims2, vmin = None, vmax = None, \
        ptitle = '', plabel = plabel0, colorbar = False, labelsize = labelsize + 1, \
        zoom=True,save=False)
    plot_GOES_satpy(date_str, 8, ax = ax9, var = var2, crs = crs0, \
        lat_lims = lat_lims2, lon_lims = lon_lims2, vmin = None, vmax = None, \
        ptitle = '', plabel = plabel2, colorbar = True, labelsize = labelsize, \
        zoom=True,save=False)
    plot_GOES_satpy(date_str, 9, ax = ax10, var = var3, crs = crs0, \
        lat_lims = lat_lims2, lon_lims = lon_lims2, vmin = None, vmax = None, \
        ptitle = '', plabel = plabel3, colorbar = True, labelsize = labelsize, \
        zoom=True,save=False)
    plot_GOES_satpy(date_str, 10, ax = ax11, var = var4, crs = crs0, \
        lat_lims = lat_lims2, lon_lims = lon_lims2, vmin = None, vmax = None, \
        ptitle = '', plabel = plabel4, colorbar = True, labelsize = labelsize, \
        zoom=True,save=False)

    ##!## Plot the GOES data 
    ##!## ------------------
    ##!#max0 = 1.0
    ##!#max1 = 1.0
    ##!#max2 = 1.0 
    ##!#max3 = 1.0 
    ##!#max4 = 1.0 
    ##!#max5 = 1.0 
    ##!#min0 = 0.0
    ##!#min1 = 0.0
    ##!#min2 = 0.5
    ##!#min3 = 0.6
    ##!#min4 = 0.6
    ##!#min5 = 0.6
    ##!#mesh8 = ax8.imshow(var0.data, transform = crs0, extent=(var0.x[0], var0.x[-1], \
    ##!#    var0.y[-1], var0.y[0]), origin='upper', cmap = 'Greys_r', \
    ##!#    vmin = min2, vmax = max2)
    ##!#mesh9 = ax9.imshow(var2.data, transform = crs0, extent=(var0.x[0], var0.x[-1], \
    ##!#    var0.y[-1], var0.y[0]), origin='upper', cmap = 'Greys_r', \
    ##!#    vmin = min3, vmax = max3)
    ##!#mesh10 = ax10.imshow(var3.data, transform = crs0, extent=(var0.x[0], var0.x[-1], \
    ##!#    var0.y[-1], var0.y[0]), origin='upper', cmap = 'Greys_r', \
    ##!#    vmin = min4, vmax = max4)
    ##!#mesh11 = ax11.imshow(var4.data, transform = crs0, extent=(var0.x[0], var0.x[-1], \
    ##!#    var0.y[-1], var0.y[0]), origin='upper', cmap = 'Greys_r', \
    ##!#    vmin = min5, vmax = max5)
    ##!#cbar8 = plt.colorbar(mesh8,ax=ax8,orientation='vertical',\
    ##!#    pad=0.03, shrink = 0.67)
    ##!#cbar9 = plt.colorbar(mesh9,ax=ax9,orientation='vertical',\
    ##!#    pad=0.03, shrink = 0.67)
    ##!#cbar10 = plt.colorbar(mesh10,ax=ax10,orientation='vertical',\
    ##!#    pad=0.03, shrink = 0.67)
    ##!#cbar11 = plt.colorbar(mesh11,ax=ax11,orientation='vertical',\
    ##!#    pad=0.03, shrink = 0.67)

    # Zoom in the figure if desired
    # -----------------------------
    if(zoom):
        ax1.set_extent([lon_lims1[0],lon_lims1[1],lat_lims1[0],\
            lat_lims1[1]],crs = datacrs)
    ##!#    ax8.set_extent([lon_lims0[0],lon_lims0[1],lat_lims0[0],\
    ##!#        lat_lims0[1]],crs = datacrs)
    ##!#    ax9.set_extent([lon_lims0[0],lon_lims0[1],lat_lims0[0],\
    ##!#        lat_lims0[1]],crs = datacrs)
    ##!#    ax10.set_extent([lon_lims0[0],lon_lims0[1],lat_lims0[0],\
    ##!#        lat_lims0[1]],crs = datacrs)
    ##!#    ax11.set_extent([lon_lims0[0],lon_lims0[1],lat_lims0[0],\
    ##!#        lat_lims0[1]],crs = datacrs)

    # ----------------------------------------------------------------------
    #
    # Plot the data in the figure
    #
    # ----------------------------------------------------------------------

    # Plot channel 1, 5, 31, and WV data spatial data
    # -----------------------------------------------
    plot_MODIS_spatial(MODIS_data_ch1, ax2, zoom = zoom, ptitle = '', labelsize = 10, labelticksize = 9)
    plot_MODIS_spatial(MODIS_data_ch5,  ax3, zoom = zoom, ptitle = '', labelsize = 10, labelticksize = 9)
    plot_MODIS_spatial(MODIS_data_ch31,  ax4, zoom = zoom, ptitle = '', labelsize = 10, labelticksize = 9)
    #plot_MODIS_spatial(MODIS_data_wv,   ax7, zoom = zoom, ptitle = '')


    # Plot the CERES SWF and LWF data
    # -------------------------------
    ##!#plot_CERES_spatial(date_str, mask_LAT, mask_LON, mask_swf, 'SWF', ax8, \
    ##!#    ptitle = '', zoom = zoom)
    ##!#plot_CERES_spatial(date_str, mask_LAT, mask_LON, mask_lwf, 'LWF', ax9, \
    ##!#    ptitle = '', zoom = zoom)
    plotCERES_hrly(ax5, CERES_data_hrly_swf, 'swf', \
        vmin = 100, vmax = 300, title = '', label = 'TOA Flux [W/m$^{2}$]', \
        labelsize = 10, labelticksize = 9, circle_bound = False, \
        gridlines = False, grid_data = True, \
        zoom = True)
    plotCERES_hrly(ax6, CERES_data_hrly_swf, 'lwf', \
        vmin = None, vmax = None, title = '', label = 'TOA Flux [W/m$^{2}$]', \
        labelsize = 10, labelticksize = 9, circle_bound = False, \
        gridlines = False, grid_data = True, \
        zoom = True)
    plotCERES_hrly(ax7, CERES_data_hrly_swf, 'total', \
        vmin = None, vmax = None, title = '', label = 'TOA Flux [W/m$^{2}$]', \
        labelsize = 10, labelticksize = 9, circle_bound = False, \
        gridlines = False, grid_data = True, \
        zoom = True)
    if(zoom):
        ax5.set_extent([aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][1], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][1]],\
                        datacrs)
        ax6.set_extent([aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][1], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][1]],\
                        datacrs)
        ax7.set_extent([aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][1], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][1]],\
                        datacrs)


    #plotCERES_hrly(ax8, CERES_data_hrly, minlat=65, \
    #    vmin = None, vmax = None, title = '', label = '', \
    #    circle_bound = True, gridlines = True, grid_data = True, \
    #    zoom = True):

    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(date_str) 

        plt.rcParams.update({'hatch.color': 'r'})
        ax3.pcolor(MODIS_data_ch1['lon'],MODIS_data_ch1['lat'],\
            hash_data1, hatch = '\\\\', alpha=0., transform = datacrs,\
            cmap = 'plasma')


    ##!## Plot the ASOS locations on the map
    ##!## ----------------------------------
    ##!#plot_ASOS_locs(ax1,date_str,color='lime', sites = ['O05','AAT'])

    # Add subplot labels
    # ------------------
    font_size = 10
    plot_subplot_label(ax1,  '(a)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax2,  '(b)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax3,  '(c)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax4,  '(d)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax5,  '(e)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax6,  '(f)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax7,  '(g)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax8,  '(h)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax9,  '(i)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax10, '(j)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax11, '(k)', backgroundcolor = 'white', fontsize = font_size)

    # Add plot text
    # -------------
    plot_figure_text(ax1, 'MODIS True Color', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax2, 'MODIS 0.64 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax3, 'MODIS 2.2 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax4, 'MODIS 11.0 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax5, 'CERES SW', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax6, 'CERES LW', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax7, 'CERES Total', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax8, 'GOES-17 11.0 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax9, 'GOES-17 Upper WV', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size-1, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax10, 'GOES-17 Mid WV', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax11, 'GOES-17 Lower WV', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size-1, backgroundcolor = 'white', halign = 'right')

    #plt.suptitle(date_str)

    fig.tight_layout()

    MODIS_data_ch1.clear()
    MODIS_data_ch5.clear()
    MODIS_data_ch31.clear()

    if(save):
        outname = 'modis_total_combined_' + date_str + '_fig1_v2.png'
        fig.savefig(outname, dpi=300)
        print("Saved",outname)
    else:
        plt.show()

def plot_combined_figure1_v3(date_str = '202107222110', zoom = True, show_smoke = True, composite = True, \
        save=False):

    #date_str = '202107202125'
    #date_str = '202107222110'
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    if(home_dir + '/Research/GOES' not in sys.path):
        sys.path.append(home_dir + '/Research/GOES')
    from GOESLib import read_GOES_satpy, plot_GOES_satpy

    # ----------------------------------------------------------------------
    #
    # Read the MODIS and CERES data
    #
    # ----------------------------------------------------------------------

    # Call read_MODIS_channel to read the desired MODIS data from the file 
    # and put it in a dictionary
    # ---------------------------------------------------------------------
    MODIS_data_ch1  = read_MODIS_channel(date_str, 1, zoom = zoom)
    MODIS_data_ch7  = read_MODIS_channel(date_str, 7, zoom = zoom)
    MODIS_data_ch31 = read_MODIS_channel(date_str, 31, zoom = zoom)

    # Determine where the smoke is located
    # ------------------------------------
    hash_data, nohash_data = find_plume(dt_date_str.strftime('%Y%m%d%H%M')) 

    tmp_data1  = np.copy(MODIS_data_ch1['data'])
    tmp_data7  = np.copy(MODIS_data_ch7['data'])
    tmp_data31 = np.copy(MODIS_data_ch31['data'])
    tmp_lat0   = np.copy(MODIS_data_ch1['lat'])
    tmp_lon0   = np.copy(MODIS_data_ch1['lon'])

    if(not (tmp_data1.shape == tmp_data7.shape == tmp_data31.shape == \
            hash_data.shape)):
        print("shape mismatch")
        shapes = []
        shapes.append(tmp_data1.shape)
        shapes.append(tmp_data7.shape)
        shapes.append(tmp_data31.shape)
        shapes.append(hash_data.shape)

        min_shape = min(shapes)
        print(min_shape)

        tmp_data1  = tmp_data1[:min_shape[0],:min_shape[1]]
        tmp_data7  = tmp_data7[:min_shape[0],:min_shape[1]]
        tmp_data31 = tmp_data31[:min_shape[0],:min_shape[1]]
        tmp_lat0   = tmp_lat0[:min_shape[0],:min_shape[1]]
        tmp_lon0   = tmp_lon0[:min_shape[0],:min_shape[1]]
        hash_data  = hash_data[:min_shape[0],:min_shape[1]]

    max_ch = 350.

    tmp_data1 = np.ma.masked_where( (abs(tmp_data7) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch), tmp_data1)
    tmp_data7 = np.ma.masked_where( (abs(tmp_data7) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch), tmp_data7)
    tmp_data31 = np.ma.masked_where( (abs(tmp_data7) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch), tmp_data31)
    tmp_lat0 = np.ma.masked_where( (abs(tmp_data7) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch), tmp_lat0)
    tmp_lon0 = np.ma.masked_where( (abs(tmp_data7) > max_ch) | \
        (abs(tmp_data1) > max_ch) | (abs(tmp_data31) > max_ch), tmp_lon0)

    # Read the true color data
    # ------------------------
    var1, crs1, lat_lims1, lon_lims1 = read_true_color(date_str,\
        composite=composite)

    # Read the GOES data
    # ------------------------
    var2, crs0, lons, lats, lat_lims2, lon_lims2, plabel2 = read_GOES_satpy(date_str, 8)
    var3, crs0, lons, lats, lat_lims0, lon_lims0, plabel3 = read_GOES_satpy(date_str, 9)
    var4, crs0, lons, lats, lat_lims0, lon_lims0, plabel4 = read_GOES_satpy(date_str, 10)

    # ----------------------------------------------------------------------
    #
    #  Set up the 6-panel figure
    #
    # ----------------------------------------------------------------------

    mapcrs = init_proj(date_str)
    plt.close('all')
    fig = plt.figure(figsize=(9,9))
    #gs = fig.add_gridspec(nrows = 2, ncols = 8)
    ax1  = fig.add_subplot(3,3,1,  projection = crs1)   # true color    
    ax2  = fig.add_subplot(3,3,2,  projection = mapcrs) # Ch 1
    ax3  = fig.add_subplot(3,3,3,  projection = mapcrs) # Ch 7
    ax4  = fig.add_subplot(3,3,4,  projection = mapcrs) # Ch 31
    ax5  = fig.add_subplot(3,3,5)                       # IR vs VIS
    ax6  = fig.add_subplot(3,3,6)                       # IR vs SWIR
    ax7  = fig.add_subplot(3,3,7,  projection = crs0) # Ch 1
    ax8  = fig.add_subplot(3,3,8,  projection = crs0) # Ch 7
    ax9  = fig.add_subplot(3,3,9,  projection = crs0) # Ch 31
    ##!#ax1  = fig.add_subplot(gs[0:2,0:2],projection = crs1)   # true color    
    ##!#ax2  = fig.add_subplot(gs[0,2:4],  projection = mapcrs) # Ch 1
    ##!#ax3  = fig.add_subplot(gs[0,4:6],  projection = mapcrs) # Ch 7
    ##!#ax4  = fig.add_subplot(gs[0,6:8],  projection = mapcrs) # Ch 31
    ##!#ax5  = fig.add_subplot(gs[1,2:5])                       # IR vs VIS
    ##!#ax6  = fig.add_subplot(gs[1,5:8])                       # IR vs SWIR



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
    plot_MODIS_spatial(MODIS_data_ch1,  ax2, zoom = zoom, ptitle = '', labelsize = 10, labelticksize = 9)
    plot_MODIS_spatial(MODIS_data_ch7,  ax3, zoom = zoom, ptitle = '', labelsize = 10, labelticksize = 9)
    plot_MODIS_spatial(MODIS_data_ch31, ax4, zoom = zoom, ptitle = '', labelsize = 10, labelticksize = 9)

    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(date_str) 

        plt.rcParams.update({'hatch.color': 'r'})
        ax2.pcolor(MODIS_data_ch1['lon'],MODIS_data_ch1['lat'],\
            hash_data1, hatch = '\\\\', alpha=0., transform = datacrs,\
            cmap = 'plasma')

    # Plot the scatter data
    # ---------------------
    plot_scatter(ax5, tmp_data31, tmp_data1, MODIS_data_ch31, MODIS_data_ch1, \
        hash_data, xlabel = '11 μm brightness temperature', \
        ylabel = '0.64 μm reflectance', plot_legend = True)
    plot_scatter(ax6, tmp_data31, tmp_data7, MODIS_data_ch31, MODIS_data_ch7, \
        hash_data, xlabel = '11 μm brightness temperature', \
        ylabel = '2.1 μm reflectance', plot_legend = True)

    labelsize = 10
    # Plot channel 1, 5, 31, and WV data spatial data
    # -----------------------------------------------
    plot_GOES_satpy(date_str, 8, ax = ax7, var = var2, crs = crs0, \
        lons = lons, lats = lats, lat_lims = lat_lims2, lon_lims = lon_lims2, \
        vmin = None, vmax = None, \
        ptitle = '', plabel = plabel2, colorbar = True, labelsize = labelsize, \
        zoom=True,save=False)
    plot_GOES_satpy(date_str, 9, ax = ax8, var = var3, crs = crs0, \
        lons = lons, lats = lats, lat_lims = lat_lims2, lon_lims = lon_lims2, \
        vmin = None, vmax = None, \
        ptitle = '', plabel = plabel3, colorbar = True, labelsize = labelsize, \
        zoom=True,save=False)
    plot_GOES_satpy(date_str, 10, ax = ax9, var = var4, crs = crs0, \
        lons = lons, lats = lats, lat_lims = lat_lims2, lon_lims = lon_lims2, \
        vmin = None, vmax = None, \
        ptitle = '', plabel = plabel4, colorbar = True, labelsize = labelsize, \
        zoom=True,save=False)

    # Add subplot labels
    # ------------------
    font_size = 10
    plot_subplot_label(ax1,  '(a)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax2,  '(b)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax3,  '(c)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax4,  '(d)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax5,  '(e)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax6,  '(f)', backgroundcolor = 'white', fontsize = font_size)

    # Add plot text
    # -------------
    plot_figure_text(ax1, 'MODIS True Color', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax2, 'MODIS 0.64 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax3, 'MODIS 2.2 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax4, 'MODIS 11.0 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    ##!#plot_figure_text(ax5, 'CERES SW', xval = None, yval = None, transform = None, \
    ##!#    color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    ##!#plot_figure_text(ax6, 'CERES LW', xval = None, yval = None, transform = None, \
    ##!#    color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')

    #plt.suptitle(date_str)

    fig.tight_layout()

    MODIS_data_ch1.clear()
    MODIS_data_ch7.clear()
    MODIS_data_ch31.clear()

    if(save):
        outname = 'modis_total_combined_' + date_str + '_fig1_v3.png'
        fig.savefig(outname, dpi=300)
        print("Saved",outname)
    else:
        plt.show()

def plot_combined_figure1_v4(date_str = '202107222110', \
        zoom = True, \
        modis_ch1 = 7, \
        modis_ch2 = 31, \
        goes_ch1 = 2, \
        goes_ch2 = 6, \
        goes_ch3 = 13, \
        goes_ch4 = 8, \
        goes_ch5 = 9, \
        goes_ch6 = 10, \
        show_smoke = False, composite = True, double_fig = False, \
        save=False):

    #date_str = '202107202125'
    date_str = '202107222110'
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    date_str2 = '202107210000'
    dt_date_str2 = datetime.strptime(date_str,"%Y%m%d%H%M")

    if(home_dir + '/Research/GOES' not in sys.path):
        sys.path.append(home_dir + '/Research/GOES')
    from GOESLib import read_GOES_satpy, plot_GOES_satpy,\
        goes_channel_dict

    # ----------------------------------------------------------------------
    #
    # Read the MODIS and CERES data
    #
    # ----------------------------------------------------------------------

    # Call read_MODIS_channel to read the desired MODIS data from the file 
    # and put it in a dictionary
    # ---------------------------------------------------------------------
    ##!#MODIS_data_ch7  = read_MODIS_channel(date_str, 6, zoom = zoom)
    ##!##MODIS_data_ch7  = read_MODIS_channel(date_str, 7, zoom = zoom)
    ##!#MODIS_data_ch31 = read_MODIS_channel(date_str, 31, zoom = zoom)
    var8, crs8, lons8, lats8, lat_lims8, lon_lims8, plabel8 = \
        read_MODIS_satpy(date_str, modis_ch1)
    var9, crs8, lons8, lats8, lat_lims9, lon_lims9, plabel9 = \
        read_MODIS_satpy(date_str, modis_ch2)

    # Read the true color data
    # ------------------------
    var1, crs1, lons1, lats1, lat_lims1, lon_lims1, plabel = \
        read_MODIS_satpy(date_str,'true_color',\
    #var1, crs1, lat_lims1, lon_lims1 = read_true_color(date_str,\
        composite=composite)

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

    # ----------------------------------------------------------------------
    #
    #  Set up the 6-panel figure
    #
    # ----------------------------------------------------------------------

    mapcrs = init_proj(date_str)
    plt.close('all')
    if(double_fig):
        fig1 = plt.figure(figsize=(8.2,3))
        fig2 = plt.figure(figsize=(9.5,6))
        #gs = fig.add_gridspec(nrows = 2, ncols = 8)
        ax1  = fig1.add_subplot(1,3,1,  projection = crs1)   # true color    
        ax2  = fig1.add_subplot(1,3,2,  projection = crs8) # Ch 7
        ax3  = fig1.add_subplot(1,3,3,  projection = crs8) # Ch 31
        #ax2  = fig.add_subplot(3,3,2,  projection = mapcrs) # Ch 7
        #ax3  = fig.add_subplot(3,3,3,  projection = mapcrs) # Ch 31
        ax4  = fig2.add_subplot(2,3,1,  projection = crs0)   # GOES vis 
        ax5  = fig2.add_subplot(2,3,2,  projection = crs0)   # GOES SWIR
        ax6  = fig2.add_subplot(2,3,3,  projection = crs0)   # GOES TIR 
        ax7  = fig2.add_subplot(2,3,4,  projection = crs0)   # GOES upper WV
        ax8  = fig2.add_subplot(2,3,5,  projection = crs0)   # GOES midle WV
        ax9  = fig2.add_subplot(2,3,6,  projection = crs0)   # GOES lower WV
    else:
        fig = plt.figure(figsize=(10,9))
        #gs = fig.add_gridspec(nrows = 2, ncols = 8)
        ax1  = fig.add_subplot(3,3,1,  projection = crs1)   # true color    
        ax2  = fig.add_subplot(3,3,2,  projection = crs8) # Ch 7
        ax3  = fig.add_subplot(3,3,3,  projection = crs8) # Ch 31
        #ax2  = fig.add_subplot(3,3,2,  projection = mapcrs) # Ch 7
        #ax3  = fig.add_subplot(3,3,3,  projection = mapcrs) # Ch 31
        ax4  = fig.add_subplot(3,3,4,  projection = crs0)   # GOES vis 
        ax5  = fig.add_subplot(3,3,5,  projection = crs0)   # GOES SWIR
        ax6  = fig.add_subplot(3,3,6,  projection = crs0)   # GOES TIR 
        ax7  = fig.add_subplot(3,3,7,  projection = crs0)   # GOES upper WV
        ax8  = fig.add_subplot(3,3,8,  projection = crs0)   # GOES midle WV
        ax9  = fig.add_subplot(3,3,9,  projection = crs0)   # GOES lower WV


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
    ##!#plot_MODIS_spatial(MODIS_data_ch7,  ax2, zoom = zoom, ptitle = '', labelsize = 10, labelticksize = 9)
    ##!#plot_MODIS_spatial(MODIS_data_ch31, ax3, zoom = zoom, ptitle = '', labelsize = 10, labelticksize = 9)
    labelsize = 10
    plot_MODIS_satpy(date_str, 5, ax = ax2, var = var8, crs = crs8, \
        lons = lons8, lats = lats8, lat_lims = lat_lims8, lon_lims = lon_lims8, \
        vmin = None, vmax = None, ptitle = '', plabel = plabel8, \
        labelsize = 10, colorbar = True, zoom=True,save=False)
    plot_MODIS_satpy(date_str, 31, ax = ax3, var = var9, crs = crs8, \
        lons = lons8, lats = lats8, lat_lims = lat_lims9, lon_lims = lon_lims9, \
        vmin = 270, vmax = 330, ptitle = '', plabel = plabel9, \
        labelsize = 10, colorbar = True, zoom=True,save=False)

    # Plot channel 1, 5, 31, and WV data spatial data
    # -----------------------------------------------
    plot_GOES_satpy(date_str, 2, ax = ax4, var = var2, crs = crs0, \
        lons = lons, lats = lats, lat_lims = lat_lims2, lon_lims = lon_lims2, \
        vmin = None, vmax = 60, \
        ptitle = '', plabel = plabel2, colorbar = True, labelsize = labelsize, \
        zoom=True,save=False)
    print(var2.shape, lons.shape, lats.shape, lons2.shape, lats2.shape, var3.shape)
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

    # Add subplot labels
    # ------------------
    font_size = 10
    if(double_fig):
        plabels = ['(a)','(b)','(c)','(d)','(e)','(f)']
    else:
        plabels = ['(d)','(e)','(f)','(g)','(h)','(i)']
    plot_subplot_label(ax1, '(a)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax2, '(b)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax3, '(c)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax4, plabels[0], backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax5, plabels[1], backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax6, plabels[2], backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax7, plabels[3], backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax8, plabels[4], backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax9, plabels[5], backgroundcolor = 'white', fontsize = font_size)

    # Add plot text
    # -------------
    plot_figure_text(ax1, 'MODIS True Color', xval = None, yval = None, transform = None, \
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
        fig.text(0.03, 0.83, '---- 21:10 UTC 22 July 2021 ----', ha='center', va='center', \
            rotation='vertical', weight = 'bold', fontsize = labelsize + 2)
        fig.text(0.03, 0.335, '------------------------- 00:00 UTC 21 July 2021 -------------------------', ha='center', va='center', \
            rotation='vertical', weight = 'bold', fontsize = labelsize + 2)
        #plt.suptitle(date_str)

        fig.tight_layout()
    #MODIS_data_ch7.clear()
    #MODIS_data_ch31.clear()

    if(save):
        if(double_fig):
            outname1 = 'modis_total_combined_' + date_str + '_fig1_v4_modis.png'
            outname2 = 'modis_total_combined_' + date_str + '_fig1_v4_goes.png'
            fig1.savefig(outname1, dpi=300)
            fig2.savefig(outname2, dpi=300)
            print("Saved",outname1)
            print("Saved",outname2)
        else:
            outname = 'modis_total_combined_' + date_str + '_fig1_v42.png'
            fig.savefig(outname, dpi=300)
            print("Saved",outname)
    else:
        plt.show()

def plot_combined_figure1_v5(date_str = '202107222110', \
        modis_ch1 = 7, modis_ch2 = 31, \
        goes_ch1 = 2, goes_ch2 = 6, goes_ch3 = 13, \
        goes_ch4 = 8, goes_ch5 = 9, goes_ch6 = 10, \
        ch_idx1 = 0, ch_idx2 = 1, ch_idx3 = 2,\
        ttype1 = 'low', ttype2 = 'ml', \
        idx1 = 3, idx2 = 8, idx3 = 5, \
        date_idx = 25, 
        show_smoke = False, composite = True, double_fig = False, \
        zoom = True, save=False):

    #date_str = '202107202125'
    date_str = '202107222110'
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    date_str2 = '202107210000'
    dt_date_str2 = datetime.strptime(date_str,"%Y%m%d%H%M")

    if(home_dir + '/Research/GOES' not in sys.path):
        sys.path.append(home_dir + '/Research/GOES')
    from GOESLib import read_GOES_satpy, plot_GOES_satpy,\
        goes_channel_dict, read_GOES_time_series_NCDF

    # ----------------------------------------------------------------------
    #
    # Read the MODIS and CERES data
    #
    # ----------------------------------------------------------------------

    # Call read_MODIS_channel to read the desired MODIS data from the file 
    # and put it in a dictionary
    # ---------------------------------------------------------------------
    ##!#MODIS_data_ch7  = read_MODIS_channel(date_str, 6, zoom = zoom)
    ##!##MODIS_data_ch7  = read_MODIS_channel(date_str, 7, zoom = zoom)
    ##!#MODIS_data_ch31 = read_MODIS_channel(date_str, 31, zoom = zoom)
    var1, crs1, lons1, lats1, lat_lims1, lon_lims1, plabel1 = \
        read_MODIS_satpy(date_str, 1)
    var8, crs8, lons8, lats8, lat_lims8, lon_lims8, plabel8 = \
        read_MODIS_satpy(date_str, modis_ch1)
    var9, crs8, lons8, lats8, lat_lims9, lon_lims9, plabel9 = \
        read_MODIS_satpy(date_str, modis_ch2)

    # Read the true color data
    # ------------------------
    ##var1, crs1, lons1, lats1, lat_lims1, lon_lims1, plabel = \
    ##    read_MODIS_satpy(date_str,'true_color',\
    ###var1, crs1, lat_lims1, lon_lims1 = read_true_color(date_str,\
    ##    composite=composite)

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
        ttype1 + '_202107201201_202107210231.nc'
    file_name2 = home_dir + '/Research/GOES/goes_cross_data_' + \
        ttype2 + '_202107201201_202107210231.nc'
    GOES_dict  = read_GOES_time_series_NCDF(file_name1)
    GOES_dict2 = read_GOES_time_series_NCDF(file_name2)

    # ----------------------------------------------------------------------
    #
    #  Set up the 6-panel figure
    #
    # ----------------------------------------------------------------------

    mapcrs = init_proj(date_str)
    plt.close('all')
    if(double_fig):
        fig1 = plt.figure(figsize=(8.2,3))
        fig2 = plt.figure(figsize=(9.5,6))
        #gs = fig.add_gridspec(nrows = 2, ncols = 8)
        ax1  = fig1.add_subplot(1,3,1,  projection = crs1)   # true color    
        ax2  = fig1.add_subplot(1,3,2,  projection = crs8) # Ch 7
        ax3  = fig1.add_subplot(1,3,3,  projection = crs8) # Ch 31
        #ax2  = fig.add_subplot(3,3,2,  projection = mapcrs) # Ch 7
        #ax3  = fig.add_subplot(3,3,3,  projection = mapcrs) # Ch 31
        ax4  = fig2.add_subplot(2,3,1,  projection = crs0)   # GOES vis 
        ax5  = fig2.add_subplot(2,3,2,  projection = crs0)   # GOES SWIR
        ax6  = fig2.add_subplot(2,3,3,  projection = crs0)   # GOES TIR 
        ax7  = fig2.add_subplot(2,3,4,  projection = crs0)   # GOES upper WV
        ax8  = fig2.add_subplot(2,3,5,  projection = crs0)   # GOES midle WV
        ax9  = fig2.add_subplot(2,3,6,  projection = crs0)   # GOES lower WV
    else:
        fig = plt.figure(figsize=(9.5,11))
        gs = fig.add_gridspec(nrows = 4, ncols = 3)
        ax1  = fig.add_subplot(gs[0,0], projection = crs1) # true color    
        ax2  = fig.add_subplot(gs[0,1], projection = crs8) # MODIS Ch 7
        ax3  = fig.add_subplot(gs[0,2], projection = crs8) # MODIS Ch 31
        ax4  = fig.add_subplot(gs[1,0], projection = crs0) # GOES VIS 
        ax5  = fig.add_subplot(gs[1,1], projection = crs0)   # GOES SWIR
        ax6  = fig.add_subplot(gs[1,2], projection = crs0)   # GOES TIR 
        ax7  = fig.add_subplot(gs[2,0], projection = crs0)   # GOES upper WV
        ax8  = fig.add_subplot(gs[2,1], projection = crs0)   # GOES midle WV
        ax9  = fig.add_subplot(gs[2,2], projection = crs0)   # GOES lower WV
        ax10 = fig.add_subplot(gs[3,:]) # time series of GOES data


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
    plot_MODIS_satpy(date_str, 1, ax = ax1, var = var1, crs = crs1, \
        lons = lons1, lats = lats1, lat_lims = lat_lims1, lon_lims = lon_lims1, \
        vmin = 0, vmax = 70, ptitle = '', plabel = plabel1, \
        labelsize = 10, colorbar = True, zoom=True,save=False)
    plot_MODIS_satpy(date_str, 5, ax = ax2, var = var8, crs = crs8, \
        lons = lons8, lats = lats8, lat_lims = lat_lims8, lon_lims = lon_lims8, \
        vmin = None, vmax = None, ptitle = '', plabel = plabel8, \
        labelsize = 10, colorbar = True, zoom=True,save=False)
    plot_MODIS_satpy(date_str, 31, ax = ax3, var = var9, crs = crs8, \
        lons = lons8, lats = lats8, lat_lims = lat_lims9, lon_lims = lon_lims9, \
        vmin = 270, vmax = 330, ptitle = '', plabel = plabel9, \
        labelsize = 10, colorbar = True, zoom=True,save=False)

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

def plot_combined_figure1_v6(date_str = '202107222110', \
        modis_ch1 = 'true_color', modis_ch2 = 5, \
        modis_ch3 = 7, modis_ch4 = 20, \
        modis_ch5 = 28, modis_ch6 = 31, \
        show_smoke = False, composite = True, \
        zoom = True, save=False):

    if(home_dir + '/Research/CERES' not in sys.path):
        sys.path.append(home_dir + '/Research/CERES')
    from gridCERESLib import readgridCERES_hrly_grid, plotCERES_hrly

    date_str = '202107222110'
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # ----------------------------------------------------------------------
    #
    # Read the MODIS and CERES data
    #
    # ----------------------------------------------------------------------

    # Call read_MODIS_channel to read the desired MODIS data from the file 
    # and put it in a dictionary
    # ---------------------------------------------------------------------
    var1, crs1, lons1, lats1, lat_lims1, lon_lims1, plabel1 = \
        read_MODIS_satpy(date_str, modis_ch1)
    var2, crs2, lons2, lats2, _, _, plabel2 = \
        read_MODIS_satpy(date_str, modis_ch2)
    var3, crs3, lons3, lats3, _, _, plabel3 = \
        read_MODIS_satpy(date_str, modis_ch3)
    var4, crs4, lons4, lats4, _, _, plabel4 = \
        read_MODIS_satpy(date_str, modis_ch4)
    var5, crs5, lons5, lats5, _, _, plabel5 = \
        read_MODIS_satpy(date_str, modis_ch5)
    var6, crs6, lons6, lats6, _, _, plabel6 = \
        read_MODIS_satpy(date_str, modis_ch6)

    # ----------------------------------------------------------------------
    #
    # Read the CERES data
    #
    # ----------------------------------------------------------------------
    CERES_data_hrly_swf = readgridCERES_hrly_grid(date_str[:10], 'SWF', \
        minlat = 20., modis_comp = True)

    # Read the true color data
    # ------------------------
    var1, crs1, lons1, lats1, lat_lims1, lon_lims1, plabel = \
        read_MODIS_satpy(date_str,'true_color',\
    #var1, crs1, lat_lims1, lon_lims1 = read_true_color(date_str,\
        composite=composite)

    # ----------------------------------------------------------------------
    #
    #  Set up the 6-panel figure
    #
    # ----------------------------------------------------------------------

    mapcrs = init_proj(date_str)
    plt.close('all')
    fig = plt.figure(figsize=(10.5,9.5))
    gs = fig.add_gridspec(nrows = 3, ncols = 3)
    ax1  = fig.add_subplot(gs[0,0], projection = crs1) # true color    
    ax2  = fig.add_subplot(gs[0,1], projection = crs2) # MODIS Ch 6
    ax3  = fig.add_subplot(gs[0,2], projection = crs3) # MODIS Ch 7
    ax4  = fig.add_subplot(gs[1,0], projection = crs4) # MODIS Ch 20
    ax5  = fig.add_subplot(gs[1,1], projection = crs5) # MODIS Ch 27
    ax6  = fig.add_subplot(gs[1,2], projection = crs6) # MODIS Ch 31
    ax7  = fig.add_subplot(gs[2,0], projection = crs1) # CERES SW 
    ax8  = fig.add_subplot(gs[2,1], projection = crs1) # CERES LW 
    ax9  = fig.add_subplot(gs[2,2], projection = crs1) # CERES Total


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
    plot_MODIS_satpy(date_str, modis_ch1, ax = ax1, var = var1, crs = crs1, \
        lons = lons1, lats = lats1, lat_lims = lat_lims1, lon_lims = lon_lims1, \
        ptitle = '', plabel = plabel1, \
        labelsize = 10, zoom=True, save=False)
    plot_MODIS_satpy(date_str, modis_ch2, ax = ax2, var = var2, crs = crs2, \
        lons = lons2, lats = lats2, lat_lims = lat_lims1, lon_lims = lon_lims1, \
        vmin = None, vmax = None, ptitle = '', plabel = plabel2, \
        labelsize = 10, colorbar = True, zoom=True,save=False)
    plot_MODIS_satpy(date_str, modis_ch3, ax = ax3, var = var3, crs = crs3, \
        lons = lons3, lats = lats3, lat_lims = lat_lims1, lon_lims = lon_lims1, \
        vmin = None, vmax = None, ptitle = '', plabel = plabel3, \
        labelsize = 10, colorbar = True, zoom=True,save=False)
    plot_MODIS_satpy(date_str, modis_ch4, ax = ax4, var = var4, crs = crs4, \
        lons = lons4, lats = lats4, lat_lims = lat_lims1, lon_lims = lon_lims1, \
        vmin = 275, vmax = 340, ptitle = '', plabel = plabel4, \
        labelsize = 10, colorbar = True, zoom=True,save=False)
    plot_MODIS_satpy(date_str, modis_ch5, ax = ax5, var = var5, crs = crs5, \
        lons = lons5, lats = lats5, lat_lims = lat_lims1, lon_lims = lon_lims1, \
        vmin = 255, vmax = 275, ptitle = '', plabel = plabel5, \
        labelsize = 10, colorbar = True, zoom=True,save=False)
    plot_MODIS_satpy(date_str, modis_ch6, ax = ax6, var = var6, crs = crs6, \
        lons = lons6, lats = lats6, lat_lims = lat_lims1, lon_lims = lon_lims1, \
        vmin = 275, vmax = 330, ptitle = '', plabel = plabel6, \
        labelsize = 10, colorbar = True, zoom=True,save=False)

    # Plot the CERES  data
    # --------------------
    if(date_str == '202107222110'):
        sw_vmax = 250
        sw_vmin = 130
        lw_vmax = 370
        lw_vmin = 300
        tot_vmax = 580
        tot_vmin = 460
    elif(date_str == '202108062025'):
        sw_vmax = 330
        sw_vmin = 190
        lw_vmax = 370
        lw_vmin = 300
   
    plotCERES_hrly(ax7, CERES_data_hrly_swf, 'swf', \
        vmin = sw_vmin, vmax = sw_vmax, title = '', label = 'TOA Flux [W/m$^{2}$]', \
        circle_bound = False, gridlines = False, grid_data = True, \
        zoom = True)
    plotCERES_hrly(ax8, CERES_data_hrly_swf, 'lwf', \
        vmin = lw_vmin, vmax = lw_vmax, title = '', label = 'TOA Flux [W/m$^{2}$]', \
        circle_bound = False, gridlines = False, grid_data = True, \
        zoom = True)
    plotCERES_hrly(ax9, CERES_data_hrly_swf, 'total', \
        vmin = tot_vmin, vmax = tot_vmax, title = '', label = 'TOA Flux [W/m$^{2}$]', \
        circle_bound = False, gridlines = False, grid_data = True, \
        zoom = True)
    if(zoom):
        ax7.set_extent([lon_lims1[0],lon_lims1[1],lat_lims1[0],lat_lims1[1]],\
                       datacrs)
        ax8.set_extent([lon_lims1[0],lon_lims1[1],lat_lims1[0],lat_lims1[1]],\
                       datacrs)
        ax9.set_extent([lon_lims1[0],lon_lims1[1],lat_lims1[0],lat_lims1[1]],\
                       datacrs)

    ##!#point_size = 5
    ##!#ax4.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
    ##!#        linewidth=2, markersize = point_size + 2, marker='.',
    ##!#        color = 'black', transform=datacrs)
    ##!#ax4.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
    ##!#        linewidth=2, markersize = point_size, marker='.',
    ##!#        transform=datacrs)
    ##!#ax4.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
    ##!#        linewidth=2, markersize = point_size + 2, marker='.',
    ##!#        color = 'black', transform=datacrs)
    ##!#ax4.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
    ##!#        linewidth=2, markersize = point_size, marker='.',
    ##!#        transform=datacrs)
    ##!#ax4.plot(GOES_dict2['plon'][idx3], GOES_dict2['plat'][idx3], \
    ##!#        linewidth=2, markersize = point_size + 2, marker='.',
    ##!#        color = 'black', transform=datacrs)
    ##!#ax4.plot(GOES_dict2['plon'][idx3], GOES_dict2['plat'][idx3], \
    ##!#        linewidth=2, markersize = point_size, marker='.',
    ##!#        transform=datacrs)
    ##!#ax5.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
    ##!#        linewidth=2, markersize = point_size + 2, marker='.',
    ##!#        color = 'black', transform=datacrs)
    ##!#ax5.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
    ##!#        linewidth=2, markersize = point_size, marker='.',
    ##!#        transform=datacrs)
    ##!#ax5.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
    ##!#        linewidth=2, markersize = point_size + 2, marker='.',
    ##!#        color = 'black', transform=datacrs)
    ##!#ax5.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
    ##!#        linewidth=2, markersize = point_size, marker='.',
    ##!#        transform=datacrs)
    ##!#ax5.plot(GOES_dict2['plon'][idx3], GOES_dict2['plat'][idx3], \
    ##!#        linewidth=2, markersize = point_size + 2, marker='.',
    ##!#        color = 'black', transform=datacrs)
    ##!#ax5.plot(GOES_dict2['plon'][idx3], GOES_dict2['plat'][idx3], \
    ##!#        linewidth=2, markersize = point_size, marker='.',
    ##!#        transform=datacrs)
    ##!#ax6.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
    ##!#        linewidth=2, markersize = point_size + 2, marker='.',
    ##!#        color = 'black', transform=datacrs)
    ##!#ax6.plot(GOES_dict['plon'][idx1], GOES_dict['plat'][idx1], \
    ##!#        linewidth=2, markersize = point_size, marker='.',
    ##!#        transform=datacrs)
    ##!#ax6.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
    ##!#        linewidth=2, markersize = point_size + 2, marker='.',
    ##!#        color = 'black', transform=datacrs)
    ##!#ax6.plot(GOES_dict['plon'][idx2], GOES_dict['plat'][idx2], \
    ##!#        linewidth=2, markersize = point_size, marker='.',
    ##!#        transform=datacrs)
    ##!#ax6.plot(GOES_dict2['plon'][idx3], GOES_dict2['plat'][idx3], \
    ##!#        linewidth=2, markersize = point_size + 2, marker='.',
    ##!#        color = 'black', transform=datacrs)
    ##!#ax6.plot(GOES_dict2['plon'][idx3], GOES_dict2['plat'][idx3], \
    ##!#        linewidth=2, markersize = point_size, marker='.',
    ##!#        transform=datacrs)

    # Add subplot labels
    # ------------------
    font_size = 10
    plot_subplot_label(ax1, '(a)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax2, '(b)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax3, '(c)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax4, '(d)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax5, '(e)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax6, '(f)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax7, '(g)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax8, '(h)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax9, '(i)', backgroundcolor = 'white', fontsize = font_size)

    # Add plot text
    # -------------
    plot_figure_text(ax1, 'MODIS True Color', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax2, 'MODIS ' + \
        str(np.round(np.mean(channel_dict[str(modis_ch2)]['Bandwidth']), 2)) \
        + ' μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    #plot_figure_text(ax2, 'MODIS 2.2 μm', xval = None, yval = None, transform = None, \
    #    color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax3, 'MODIS ' + \
        str(np.round(np.mean(channel_dict[str(modis_ch3)]['Bandwidth']), 2)) \
        + ' μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax4, 'MODIS ' + \
        str(np.round(np.mean(channel_dict[str(modis_ch4)]['Bandwidth']), 2)) \
        + ' μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax5, 'MODIS ' + \
        str(np.round(np.mean(channel_dict[str(modis_ch5)]['Bandwidth']), 2)) \
        + ' μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax6, 'MODIS ' + \
        str(np.round(np.mean(channel_dict[str(modis_ch6)]['Bandwidth']), 2)) \
        + ' μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax7, 'CERES SW', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax8, 'CERES LW', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax9, 'CERES Total', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')

    fig.suptitle(dt_date_str.strftime('Aqua MODIS and CERES Imagery of the Dixie Fire\n%d %B %Y %H:%M UTC'))

    fig.tight_layout()

    if(save):
            outname = 'modis_total_combined_' + date_str + '_fig1_v6.png'
            fig.savefig(outname, dpi=300)
            print("Saved",outname)
    else:
        plt.show()

  
def plot_figure2(modis_ch1 = 'true_color', save=False, composite = True):
        

    # Read true color data for the date
    # ---------------------------------
    date_str = '202107222110'
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")
    var1, crs1, lons1, lats1, lat_lims1, lon_lims1, plabel = \
        read_MODIS_satpy(date_str,'true_color', composite=composite)

    # Set up the figure
    # -----------------
    fig = plt.figure(figsize=(10.5,4.5))
    mapcrs = init_proj(date_str)
    ax1 = fig.add_subplot(1,2,1,projection = crs1) # true color 7/22
    ax2 = fig.add_subplot(1,2,2) # meteo 7/22

    # ----------------------------------------------------------------------
    #
    # Panel 1: Visible image
    #
    # ----------------------------------------------------------------------

    # Plot channel 31 spatial data
    # -----------------------------------------------
    plot_MODIS_satpy(date_str, modis_ch1, ax = ax1, var = var1, crs = crs1, \
        lons = lons1, lats = lats1, lat_lims = lat_lims1, lon_lims = lon_lims1, \
        ptitle = '', plabel = plabel, \
        labelsize = 10, zoom=True, save=False)

    font_size = 10
    plot_figure_text(ax1, 'MODIS True Color', \
        xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', \
        halign = 'right')

    ##!#ax1.set_extent([lon_lims3[0],lon_lims3[1],lat_lims3[0],\
    ##!#    lat_lims3[1]],crs = datacrs)

    # Plot the ASOS locations on the map
    # ----------------------------------
    plot_ASOS_locs(ax1, date_str, color='tab:blue', sites = ['AAT'])
    plot_ASOS_locs(ax1, date_str, color='tab:orange', sites = ['O05'])


    # ----------------------------------------------------------------------
    #
    # Panel 2: Meteogram
    #
    # ----------------------------------------------------------------------
    plot_asos_diurnal(ax2, date_str, 'O05', 'AAT')

    def plot_modis_line(dt_date, pax):
        local_modis_time = dt_date - timedelta(hours = 7)
        modis_diff = (local_modis_time - datetime(year=2021,month=7,\
            day=22,hour=0,minute=0)).seconds
        print(local_modis_time, modis_diff)
        pax.axvline(modis_diff,color='black',linestyle = '--', lw=2,alpha=0.75,\
            )

    plot_modis_line(dt_date_str, ax2)
    ax2.legend() 
   
    # Add subplot labels
    # ------------------
    dt_date_str = dt_date_str
    ##!#xval = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]\
    ##!#    [dt_date_str.strftime('%H%M')]['Lon'][0] + \
    ##!#    (aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]\
    ##!#    [dt_date_str.strftime('%H%M')]['Lon'][1] - \
    ##!#    aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]\
    ##!#    [dt_date_str.strftime('%H%M')]['Lon'][0])*0.05
    ##!#yval = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]\
    ##!#    [dt_date_str.strftime('%H%M')]['Lat'][0] + \
    ##!#    (aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]\
    ##!#    [dt_date_str.strftime('%H%M')]['Lat'][1] - \
    ##!#    aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]\
    ##!#    [dt_date_str.strftime('%H%M')]['Lat'][0])*0.90
    #plot_subplot_label(ax1, '(a)', xval = xval, yval = yval, \
    #    transform = datacrs, backgroundcolor = 'white')
    plot_subplot_label(ax1, '(a)', backgroundcolor = 'white', location = 'upper_left')
    plot_subplot_label(ax2, '(b)', backgroundcolor = 'white', location = 'upper_right')

    fig.tight_layout()
 
    if(save):
        outname = 'modis_asos_meteo_combined_20210722_fig2.png'
        fig.savefig(outname, dpi=300)
        print("Saved image",outname)
    else: 
        plt.show() 

def plot_MODIS_GOES_SBDART(save=False, composite = True, calc_radiance = True):

    if(home_dir + '/Research/SBDART' not in sys.path):
        sys.path.append(home_dir + '/Research/SBDART')
    from SBDART_Lib import run_sbdart, plot_bright_vza, read_atmos_profile

    # Run SBDART for the different channel
    # ------------------------------------
    ##!#goes17_ch08 = run_sbdart('goes17_ch08', calc_radiance, run = True, atms_file = '/home/bsorenson/Research/SBDART/data/model/210722_220000_XXX_HRRR.txt')
    goes17_ch08 = run_sbdart('goes17_ch08', calc_radiance, run = True)
    goes17_ch09 = run_sbdart('goes17_ch09', calc_radiance, run = True)
    goes17_ch10 = run_sbdart('goes17_ch10', calc_radiance, run = True)
    modis31     = run_sbdart('modis_ch31',  calc_radiance, run = True)

    # Set up the figure
    # -----------------
    fig = plt.figure(figsize=(10,7))
    axs = fig.subplots(nrows = 2, ncols = 2)
    #ax1 = fig.add_subplot(2,2,1)
    #ax2 = fig.add_subplot(2,2,2)
    #ax3 = fig.add_subplot(2,2,3)
    #ax4 = fig.add_subplot(2,2,4)
    
##!#    fig2 = plt.figure()
##!#    ax5  = fig2.add_subplot(1,1,1)
##!#
##!#    atms_data = read_atmos_profile(infile = '/home/bsorenson/Research/SBDART/data/model/210722_220000_XXX_HRRR.txt')
##!#    atms_data_new = atms_data.copy(deep = True)
##!#    ax5.plot(atms_data['wv'], atms_data['z'], label = 'original')
##!#    # Add 2
##!#    atms_data_new['wv'][atms_data['z'] <= 5.] = atms_data['wv'][atms_data['z'] <= 5.] + 2.
##!#    ax5.plot(atms_data_new['wv'], atms_data['z'], label = '+ 2 g/m3')
##!#    # Add 4
##!#    atms_data_new['wv'][atms_data['z'] <= 5.] = atms_data['wv'][atms_data['z'] <= 5.] + 4.
##!#    ax5.plot(atms_data_new['wv'], atms_data['z'], label = '+ 4 g/m3')
##!#    # Add 8
##!#    atms_data_new['wv'][atms_data['z'] <= 5.] = atms_data['wv'][atms_data['z'] <= 5.] + 8.
##!#    ax5.plot(atms_data_new['wv'], atms_data['z'], label = '+ 8 g/m3')
    
    
    # ----------------------------------------------------------------------
    #
    # Run the GOES-17 and MODIS simulations
    #
    # ----------------------------------------------------------------------
    plot_bright_vza(goes17_ch08, pax = axs[0,0])
    plot_bright_vza(goes17_ch09, pax = axs[0,1])
    plot_bright_vza(goes17_ch10, pax = axs[1,0])
    plot_bright_vza(modis31, pax = axs[1,1])
  
    #ax3.legend(loc = 'upper center', bbox_to_anchor = (0.8, -0.05), \
    #    ncol = 4)
 
    # Add subplot labels
    # ------------------
    plot_subplot_label(axs[0,0], '(a)', location = 'upper_right')
    plot_subplot_label(axs[0,1], '(b)', location = 'upper_right')
    plot_subplot_label(axs[1,0], '(c)', location = 'upper_right')
    plot_subplot_label(axs[1,1], '(d)', location = 'upper_right')

    ## Add a legend
    ## ------------
    #axLine, axLabel = axs[0,0].get_legend_handles_labels()
    #fig.legend(axLine, axLabel, loc = 'upper center', bbox_to_anchor = (0.5, 0.05), ncol = 4)

    fig.tight_layout()
 
    if(save):
        outname = 'modis_goes_sbdart_comps.png'
        fig.savefig(outname, dpi=300)
        print("Saved image",outname)
    else: 
        plt.show() 

def plot_MODIS_VIIRS_SBDART(save=False, composite = True, calc_radiance = True):

    if(home_dir + '/Research/SBDART' not in sys.path):
        sys.path.append(home_dir + '/Research/SBDART')
    from SBDART_Lib import run_sbdart, plot_bright_vza

    # Run SBDART for the different channel
    # ------------------------------------
    modis31    = run_sbdart('modis_ch31', calc_radiance, run = True)
    viirs_ch20 = run_sbdart('viirs_ch20', calc_radiance, run = True)

    # Set up the figure
    # -----------------
    fig = plt.figure(figsize=(11,4))
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)

    # ----------------------------------------------------------------------
    #
    # Run the GOES-17 and MODIS simulations
    #
    # ----------------------------------------------------------------------
    plot_bright_vza(modis31, pax = ax1)
    plot_bright_vza(viirs_ch20, pax = ax2)
   
    # Add subplot labels
    # ------------------
    plot_subplot_label(ax1, '(a)', location = 'upper_right')
    plot_subplot_label(ax2, '(b)', location = 'upper_right')

    fig.tight_layout()
 
    if(save):
        outname = 'modis_viirs_sbdart_comps.png'
        fig.savefig(outname, dpi=300)
        print("Saved image",outname)
    else: 
        plt.show() 
 
def plot_combined_figure3(zoom = True, show_smoke = False, composite = True, \
        plume_only = False, avg_pixel = True, save=False):

    if(home_dir + '/Research/CERES' not in sys.path):
        sys.path.append(home_dir + '/Research/CERES')
    from gridCERESLib import readgridCERES_hrly_grid, plotCERES_hrly

    date_str = '202108062025'
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # ----------------------------------------------------------------------
    #
    # Read the MODIS and CERES data
    #
    # ----------------------------------------------------------------------

    # Call read_MODIS_channel to read the desired MODIS data from the file 
    # and put it in a dictionary
    # ---------------------------------------------------------------------
    MODIS_data_ch5  = read_MODIS_channel(date_str, 5, zoom = zoom)
    MODIS_data_ch31 = read_MODIS_channel(date_str, 31, zoom = zoom)
    MODIS_data_wv   = read_MODIS_channel(date_str, 'wv_ir', zoom = zoom)

    # Determine where the smoke is located MAYBE???
    # ---------------------------------------------

    # Read the true color data
    # ------------------------
    var1, crs1, lat_lims1, lon_lims1 = read_true_color(date_str,\
        composite=composite)

    # Read in the CERES data
    # ----------------------
    #mask_LAT, mask_LON, mask_swf, mask_lwf = \
    #    read_CERES_match_MODIS(date_str)
    CERES_data_hrly_swf = readgridCERES_hrly_grid(date_str[:10], 'SWF')
  
    # Read in the OMI data
    # --------------------
    LAT, LON, mask_UVAI = read_OMI_match_MODIS(date_str)
 
 
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
    plot_MODIS_spatial(MODIS_data_ch31, ax2, zoom = zoom, ptitle = '')
    plot_MODIS_spatial(MODIS_data_ch5,  ax4, zoom = zoom, ptitle = '')
    plot_MODIS_spatial(MODIS_data_wv,   ax7, zoom = zoom, ptitle = '')

    # Plot the OMI spatial data
    # -------------------------
    plot_OMI_spatial(date_str, LAT, LON, mask_UVAI, ax3, zoom = zoom)

    # Plot the Ch 31 v OMI and OMI vs flux scatter data
    # -------------------------------------------------
    plot_scatter_OMI(date_str, MODIS_data_ch31, ax8, avg_pixel = avg_pixel, \
        xlabel = '11 μm brightness temperature', \
        ylabel = 'OMI UVAI', ptitle = '')
    plot_scatter_OMI_CERES(date_str, MODIS_data_ch31, ax9, \
        avg_pixel = avg_pixel, plume_only = plume_only, ptitle = '')

    ##!#plot_scatter(ax5, tmp_data31, tmp_data1, MODIS_data_ch31, MODIS_data_ch1, \
    ##!#    hash_data, xlabel = '11 μm brightness temperature', \
    ##!#    ylabel = '0.64 μm reflectance', plot_legend = True)
    ##!#plot_scatter(ax6, tmp_data31, tmp_datawv, MODIS_data_ch31, MODIS_data_wv, \
    ##!#    hash_data, xlabel = '11 μm brightness temperature', \
    ##!#    ylabel = '1.24 μm reflectance')

    # Plot the CERES SWF and LWF data
    # -------------------------------
    ##!#plot_CERES_spatial(date_str, mask_LAT, mask_LON, mask_swf, 'SWF', ax5, \
    ##!#    ptitle = '', zoom = zoom)
    ##!#plot_CERES_spatial(date_str, mask_LAT, mask_LON, mask_lwf, 'LWF', ax6, \
    ##!#    ptitle = '', zoom = zoom)
    plotCERES_hrly(ax5, CERES_data_hrly_swf, 'swf', \
        vmin = 175, vmax = 325, title = '', label = 'TOA Flux [W/m$^{2}$]', \
        labelsize = 13, labelticksize = 11, circle_bound = False, \
        gridlines = False, grid_data = True, \
        zoom = True)
    plotCERES_hrly(ax6, CERES_data_hrly_swf, 'lwf', \
        vmin = None, vmax = None, title = '', label = 'TOA Flux [W/m$^{2}$]', \
        labelsize = 13, labelticksize = 11, circle_bound = False, \
        gridlines = False, grid_data = True, \
        zoom = True)
    if(zoom):
        ax5.set_extent([aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][1], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][1]],\
                        datacrs)
        ax6.set_extent([aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][1], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][1]],\
                        datacrs)

    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(date_str) 

        ax3.pcolor(MODIS_data_ch1['lon'],MODIS_data_ch1['lat'],\
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
    plot_figure_text(ax2, 'MODIS 11 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax3, 'OMI AI', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax4, 'MODIS 1.24 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax5, 'CERES SW', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax6, 'CERES LW', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax7, 'MODIS IR WV', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')

    fig.tight_layout()

    MODIS_data_ch5.clear()
    MODIS_data_ch31.clear()
    MODIS_data_wv.clear()

    if(save):
        outname = 'modis_total_combined_' + date_str + '_fig3.png'
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
    xval = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lon'][0] + \
        (aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lon'][1] - \
        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lon'][0])*0.05
    yval = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lat'][0] + \
        (aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lat'][1] - \
        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]\
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

    def plot_modis_line(dt_date, pax):
        local_modis_time = dt_date - timedelta(hours = 7)
        modis_diff = (local_modis_time - datetime(year=2021,month=7,\
            day=22,hour=0,minute=0)).seconds
        print(local_modis_time, modis_diff)
        pax.axvline(modis_diff,color='black',linestyle = '--', lw=2,alpha=0.75,\
            label='MODIS')

    plot_modis_line(dt_date_str20, ax4)
    plot_modis_line(dt_date_str21, ax5)
    plot_modis_line(dt_date_str23, ax7)
 
    fig.tight_layout()
 
    if(save):
        outname = 'modis_asos_meteo_combined_20210722_figureS1.png'
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
    # Read the MODIS data
    #
    # ----------------------------------------------------------------------

    # Call read_MODIS_channel to read the desired MODIS data from the file 
    # and put it in a dictionary
    # ---------------------------------------------------------------------
    MODIS_data_ch1  = read_MODIS_channel(date_str, 1, zoom = zoom)
    MODIS_data_ch5  = read_MODIS_channel(date_str, 5, zoom = zoom)
    MODIS_data_ch31 = read_MODIS_channel(date_str, 31, zoom = zoom)
    MODIS_data_wv   = read_MODIS_channel(date_str, 'wv_ir', zoom = zoom)

    # Determine where the smoke is located
    # ------------------------------------
    hash_data, nohash_data = find_plume(dt_date_str.strftime('%Y%m%d%H%M')) 

    tmp_data1  = np.copy(MODIS_data_ch1['data'])
    tmp_data5  = np.copy(MODIS_data_ch5['data'])
    tmp_data31 = np.copy(MODIS_data_ch31['data'])
    tmp_datawv = np.copy(MODIS_data_wv['data'])
    tmp_lat0   = np.copy(MODIS_data_ch1['lat'])
    tmp_lon0   = np.copy(MODIS_data_ch1['lon'])

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
    plot_MODIS_spatial(MODIS_data_ch1, ax1, zoom = zoom, ptitle = '')

    plot_scatter(ax2, tmp_data31, tmp_data1, MODIS_data_ch31, MODIS_data_ch1, \
        hash_data, xlabel = '11 μm brightness temperature', \
        ylabel = '0.64 μm reflectance', plot_legend = True)
    plot_scatter(ax3, tmp_data31, tmp_data5, MODIS_data_ch31, MODIS_data_ch5, \
        hash_data, xlabel = '11 μm brightness temperature', \
        ylabel = '1.24 μm reflectance')
    plot_scatter(ax4, tmp_data31, tmp_datawv, MODIS_data_ch31, MODIS_data_wv, \
        hash_data, xlabel = '11 μm brightness temperature', \
        ylabel = 'MODIS IR Precipitable Water')

    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(date_str) 

        ax3.pcolor(MODIS_data_ch1['lon'],MODIS_data_ch1['lat'],\
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
    plot_figure_text(ax1, 'MODIS 0.64 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    
    fig.tight_layout()

    MODIS_data_ch1.clear()
    MODIS_data_ch5.clear()
    MODIS_data_ch31.clear()
    MODIS_data_wv.clear()

    if(save):
        outname = 'modis_total_combined_' + date_str + '_figS2.png'
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
        # Read MODIS Ch 31 data
        # ---------------------
        MODIS_data1 = read_MODIS_channel(prev_str, 1, zoom = zoom)

        # Plot the MODIS ch 31 data
        # -------------------------
        mapcrs = init_proj(date_str)
        ax0 = plt.subplot(gs[0,0], projection = mapcrs)
        plot_MODIS_spatial(MODIS_data1, ax0, zoom = zoom)

        if(ir_data):
            # Read MODIS Ch 31 data
            # ---------------------
            MODIS_data4 = read_MODIS_channel(prev_str, 31, zoom = zoom)

            # Plot the MODIS ch 31 data
            # -------------------------
            mapcrs = init_proj(date_str)
            ax5 = plt.subplot(gs[1,0], projection = mapcrs)
            plot_MODIS_spatial(MODIS_data4, ax5, zoom = zoom)


    ax0.set_title('Aqua MODIS\n' + dt_prev_str.strftime('%Y-%m-%d %H:%M'))

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
        # Read MODIS Ch 31 data
        # ---------------------
        MODIS_data2 = read_MODIS_channel(date_str, 1, zoom = zoom)

        # Plot the MODIS ch 31 data
        # -------------------------
        mapcrs = init_proj(date_str)
        ax1 = plt.subplot(gs[0,1], projection = mapcrs)
        plot_MODIS_spatial(MODIS_data2, ax1, zoom = zoom)

        if(ir_data):
            # Read MODIS Ch 31 data
            # ---------------------
            MODIS_data5 = read_MODIS_channel(date_str, 31, zoom = zoom)

            # Plot the MODIS ch 31 data
            # -------------------------
            mapcrs = init_proj(date_str)
            ax6 = plt.subplot(gs[1,1], projection = mapcrs)
            plot_MODIS_spatial(MODIS_data5, ax6, zoom = zoom)

    ax1.set_title('Aqua MODIS\n' + dt_date_str.strftime('%Y-%m-%d %H:%M'))

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
        # Read MODIS Ch 31 data
        # ---------------------
        MODIS_data3 = read_MODIS_channel(post_str, 1, zoom = zoom)

        # Plot the MODIS ch 31 data
        # -------------------------
        mapcrs = init_proj(post_str)
        ax2 = plt.subplot(gs[0,2], projection = mapcrs)
        plot_MODIS_spatial(MODIS_data3, ax2, zoom = zoom)

        if(ir_data):
            # Read MODIS Ch 31 data
            # ---------------------
            MODIS_data6 = read_MODIS_channel(post_str, 31, zoom = zoom)

            # Plot the MODIS ch 31 data
            # -------------------------
            mapcrs = init_proj(date_str)
            ax7 = plt.subplot(gs[1,2], projection = mapcrs)
            plot_MODIS_spatial(MODIS_data6, ax7, zoom = zoom)

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


    ax2.set_title('Aqua MODIS\n' + dt_post_str.strftime('%Y-%m-%d %H:%M'))

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
    df = pd.read_csv(aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][\
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
    
    # Convert the file time to a aerosol_event_dict format
    #event_date = datetime.strptime(infile.split('/')[-1].split('_')[-1][:8], "%Y%m%d")
    #grabber_date = event_date.strftime('%Y-%m-%d')
    #first_time = list(aerosol_event_dict[grabber_date].keys())[0]
    
    # Pull the event time from the aerosol_event_dict
    #event_dtime = event_date + timedelta(hours = int(first_time[:2]), \
    #    minutes = int(first_time[2:4]))

    # Plot a vertical line at the time of the clear MODIS overpass
    # ------------------------------------------------------------
    ax3.axvline(dt_prev_str,color='black',linestyle = '--', lw=2,alpha=0.75,\
        label='' + dt_prev_str.strftime('%m/%d %H:%M'))
    ax4.axvline(dt_prev_str,color='black',linestyle = '--',lw=2,alpha=0.75)
    
    # Plot a vertical line at the time of the MODIS overpass
    # ------------------------------------------------------
    ax3.axvline(dt_date_str,color='black',lw=2,alpha=0.75,label='' \
        + dt_date_str.strftime('%m/%d %H:%M'))
    ax4.axvline(dt_date_str,color='black',lw=2,alpha=0.75)

    # Plot a vertical line at the time of the post MODIS overpass
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
    # Extract MODIS values at each pixel and plot on the meteograms
    #
    # -----------------------------------------------------------------------     
   

 
    # Add subplot letter labels
    # -------------------------
    xval = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lon'][0] + \
        (aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lon'][1] - \
        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lon'][0])*0.05
    yval = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lat'][0] + \
        (aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lat'][1] - \
        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]\
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
        outname = 'modis_combined_station_' + date_str + '.png'
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
        second_xaxis = True, save = False):
    
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
    if(second_xaxis):
        #local = pytz.timezone('America/Los_Angeles')
        
        second_dates = np.array([ ttime + timedelta(hours = 7) for ttime in test_times])
        second_date_strs = np.array([ttime.strftime('%H:%M') for ttime in second_dates])
        pax2 = pax.twiny()
        pax2.plot(mask_time_times, mask_time_data, color = 'tab:orange', \
            linestyle = '--')
        pax2.set_xticks(diff_local_times[::6])
        pax2.set_xticklabels(second_date_strs[::6], fontsize = 12)
        pax2.set_xlabel('Time [UTC]', fontsize = 10, weight = 'bold')
    pax.set_ylabel('2-m Temperature [$^{o}$C]', fontsize = 12, weight = 'bold')
    pax.tick_params(axis='y', labelsize = 12)
    pax.set_xlabel('Time [PDT]', fontsize = 10, weight = 'bold')
    pax.set_title(dt_date_str.strftime('ASOS 2-m Temperatures\n%d %B %Y'), fontsize = 12)
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

    ax0.set_title('Aqua MODIS\n'+dt_date_str20.strftime('%Y-%m-%d %H:%M'))
    ax1.set_title('Aqua MODIS\n'+dt_date_str21.strftime('%Y-%m-%d %H:%M'))
    ax2.set_title('Aqua MODIS\n'+dt_date_str22.strftime('%Y-%m-%d %H:%M'))
    ax3.set_title('Aqua MODIS\n'+dt_date_str23.strftime('%Y-%m-%d %H:%M'))

    plot_asos_diurnal(ax4, date_str20, 'O05', 'AAT')
    plot_asos_diurnal(ax5, date_str21, 'O05', 'AAT')
    plot_asos_diurnal(ax6, date_str22, 'O05', 'AAT')
    plot_asos_diurnal(ax7, date_str23, 'O05', 'AAT')

    plot_ASOS_locs(ax0,'asos_data_case_locs.csv', color = 'red')
    plot_ASOS_locs(ax1,'asos_data_case_locs.csv', color = 'red')
    plot_ASOS_locs(ax2,'asos_data_case_locs.csv', color = 'red')
    plot_ASOS_locs(ax3,'asos_data_case_locs.csv', color = 'red')

    dt_date_str = dt_date_str22
    xval = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lon'][0] + \
        (aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lon'][1] - \
        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lon'][0])*0.05
    yval = aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lat'][0] + \
        (aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]\
        [dt_date_str.strftime('%H%M')]['Lat'][1] - \
        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')]\
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

    def plot_modis_line(dt_date, pax):
        local_modis_time = dt_date - timedelta(hours = 7)
        modis_diff = (local_modis_time - datetime(year=2021,month=7,\
            day=22,hour=0,minute=0)).seconds
        print(local_modis_time, modis_diff)
        pax.axvline(modis_diff,color='black',linestyle = '--', lw=2,alpha=0.75,\
            label='MODIS')

    plot_modis_line(dt_date_str20, ax4)
    plot_modis_line(dt_date_str21, ax5)
    plot_modis_line(dt_date_str22, ax6)
    plot_modis_line(dt_date_str23, ax7)
 
    fig.tight_layout()
 
    if(save):
        outname = 'modis_asos_meteo_combined_20210722.png'
        fig.savefig(outname, dpi=300)
        print("Saved image",outname)
    else: 
        plt.show() 

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# End known functions for the smoke thermal IR paper
# (Sorenson et al., 2023b)
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Extract the MODIS information from a given channel at each ob point
# -------------------------------------------------------------------
def nearest_grid_values(MODIS_data):
    # Read in the correct ASOS file 
    asos_file = aerosol_event_dict[MODIS_data['cross_date']][MODIS_data['file_time']]['asos']
    df = pd.read_csv(asos_file)
    df['valid'] = pd.to_datetime(df['valid'])
    df = df.set_index('valid')

    # Pull the event time from the aerosol_event_dict
    event_date = datetime.strptime(MODIS_data['cross_date'], "%Y-%m-%d")
    first_time = MODIS_data['file_time']
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

 


# Compare colocated MODIS and ob data for two dates
def colocate_comparison(date1, date2, channel = 31):
    # Read in MODIS data for both cases
    # ---------------------------------
    dt_date_str1 = datetime.strptime(date1,"%Y%m%d%H%M")
    #dt_date_str2 = datetime.strptime(date2,"%Y%m%d%H%M")
    filename1 = aerosol_event_dict[dt_date_str1.strftime('%Y-%m-%d')][dt_date_str1.strftime('%H%M')]['modis']
    #filename2 = aerosol_event_dict[dt_date_str2.strftime('%Y-%m-%d')][dt_date_str2.strftime('%H%M')]['modis']

    MODIS_data1 = read_MODIS_channel(dt_date_str1.strftime('%Y%m%d%H%M'), channel, zoom = True)
    #MODIS_data2 = read_MODIS_channel(dt_date_str2.strftime('%Y%m%d%H%M'), channel, zoom = True)

    #print(MODIS_data1,MODIS_data2)

    # Use the colocation code to extract matching data
    # ------------------------------------------------
    compare_data1 = nearest_grid_values(MODIS_data1)
    #compare_data2 = nearest_grid_values(MODIS_data2)

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

def split_window(date_str, zoom = True, save = False):
    
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # Read 0.64 μm data
    # --------------- 
    MODIS_data1 = read_MODIS_channel(date_str, 1, zoom = zoom)

    # Read 11 μm data
    # --------------- 
    MODIS_data31 = read_MODIS_channel(date_str, 31, zoom = zoom)

    # Read 12 μm data
    # --------------- 
    MODIS_data32 = read_MODIS_channel(date_str, 32, zoom = zoom)


    MODIS_data_div = read_MODIS_channel(date_str, 32, zoom = zoom)
    MODIS_data_sub = read_MODIS_channel(date_str, 32, zoom = zoom)

    # Calculate split window stuff
    # ----------------------------
    MODIS_data_sub['data'] = MODIS_data31['data'][:,:] - MODIS_data32['data'][:,:]
    MODIS_data_div['data'] = MODIS_data31['data'][:,:] / MODIS_data32['data'][:,:]

    MODIS_data_sub['variable'] = '11 μm - 12 μm Split Window'
    MODIS_data_div['variable'] = '11 μm / 12 μm Split Window'
    MODIS_data_sub['colors'] = 'bwr'
    
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
    
    plot_MODIS_spatial(MODIS_data1, ax0, zoom = zoom, ptitle = '0.64 μm reflectance')
    plot_MODIS_spatial(MODIS_data31, ax1, zoom = zoom, ptitle = '11 μm brightness temp')
    plot_MODIS_spatial(MODIS_data32, ax2, zoom = zoom, ptitle = '12 μm brightness temp')
    plot_MODIS_spatial(MODIS_data_sub,  ax3, zoom = zoom, vmin = -2, vmax = 2, ptitle = '')
    plot_MODIS_spatial(MODIS_data_div,  ax4, zoom = zoom, vmin = 0.995, vmax = 1, ptitle = '')

    if(save):
        if(zoom):
            zoom_add = '_zoom'
        else:
            zoom_add = ''
        outname = 'modis_split_window_' + date_str + zoom_add + '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

def plot_MODIS_temporary(date_str, zoom = True, save = False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # Read 0.64 μm data
    # --------------- 
    MODIS_data_ch1  = read_MODIS_channel(date_str, 1, zoom = zoom)
    MODIS_data_ch5  = read_MODIS_channel(date_str, 5, zoom = zoom)
    MODIS_data_ch31 = read_MODIS_channel(date_str, 31, zoom = zoom)

    # Determine where the smoke is located
    # ------------------------------------
    hash_data, nohash_data = find_plume(date_str) 

    tmp_data1  = np.copy(MODIS_data_ch1['data'])
    tmp_data5  = np.copy(MODIS_data_ch5['data'])
    tmp_data31 = np.copy(MODIS_data_ch31['data'])
    tmp_lat0   = np.copy(MODIS_data_ch1['lat'])
    tmp_lon0   = np.copy(MODIS_data_ch1['lon'])

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

    plot_scatter(ax0, tmp_data31, tmp_data1, MODIS_data_ch31, MODIS_data_ch1, \
        hash_data, xlabel = '11 μm  temperature', \
        ylabel = '0.64 μm reflectance', plot_legend = True)
    plot_scatter(ax1, tmp_data31, tmp_data5, MODIS_data_ch31, MODIS_data_ch5, \
        hash_data, xlabel = '11 μm brightness temperature', \
        ylabel = '1.24 μm reflectance')

    if(save):
        outname = 'modis_scatter_twopanel_' + date_str + '.png'
        fig1.savefig(outname, dpi = 300)
        print("Saved image",outname)
    else:
        plt.show()

def plot_MODIS_temporary_4panel(date_str, zoom = True, composite = True, \
        show_smoke = True, save = False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # Read the true color data
    # ------------------------
    var1, crs1, lat_lims1, lon_lims1 = read_true_color(date_str,\
        composite=composite)


    # Read 0.64 μm data
    # --------------- 
    MODIS_data_ch1  = read_MODIS_channel(date_str, 1, zoom = zoom)
    MODIS_data_ch31 = read_MODIS_channel(date_str, 31, zoom = zoom)

    # Determine where the smoke is located
    # ------------------------------------
    hash_data, nohash_data = find_plume(date_str) 

    tmp_data1  = np.copy(MODIS_data_ch1['data'])
    tmp_data31 = np.copy(MODIS_data_ch31['data'])
    tmp_lat0   = np.copy(MODIS_data_ch1['lat'])
    tmp_lon0   = np.copy(MODIS_data_ch1['lon'])

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
    plot_MODIS_spatial(MODIS_data_ch1,  ax2, zoom = zoom, ptitle = '')
    plot_MODIS_spatial(MODIS_data_ch31, ax3, zoom = zoom, ptitle = '')

    plot_figure_text(ax1, 'True Color', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 15, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax2, 'MODIS 0.64 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 15, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax3, 'MODIS 11 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 15, backgroundcolor = 'white', halign = 'right')

    # Plot the scatter between IR and vis
    # -----------------------------------
    plot_scatter(ax4, tmp_data31, tmp_data1, MODIS_data_ch31, MODIS_data_ch1, \
        hash_data, xlabel = '11 μm brightness temperature', \
        ylabel = '0.64 μm reflectance', plot_legend = True)
    #plot_scatter(ax5, tmp_data31, tmp_data5, MODIS_data_ch31, MODIS_data_ch5, \
    #    hash_data, xlabel = '11 μm brightness temperature', \
    #    ylabel = '1.24 μm reflectance')

    # ----------------------------------------------------------------------
    #
    # Panel 5: Meteogram
    #
    # ----------------------------------------------------------------------
    plot_asos_diurnal(ax5, date_str, 'O05', 'AAT')

    def plot_modis_line(dt_date, pax):
        local_modis_time = dt_date - timedelta(hours = 7)
        modis_diff = (local_modis_time - datetime(year=2021,month=7,\
            day=22,hour=0,minute=0)).seconds
        print(local_modis_time, modis_diff)
        pax.axvline(modis_diff,color='black',linestyle = '--', lw=2,alpha=0.75,\
            label='MODIS')

    plot_modis_line(dt_date_str, ax5)
    ax5.legend() 


    if(show_smoke):
        # Determine where the smoke is located
        # ------------------------------------
        hash_data1, nohash_data1 = find_plume(date_str) 

        plt.rcParams.update({'hatch.color': 'r'})
        ax2.pcolor(MODIS_data_ch1['lon'],MODIS_data_ch1['lat'],\
            hash_data1, hatch = '\\\\', alpha=0., transform = datacrs,\
            cmap = 'plasma')
        #ax4.pcolor(MODIS_data_ch1['lon'],MODIS_data_ch1['lat'],\
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
        outname = 'modis_spatscat2_4panel_' + date_str + '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image",outname)
    else:
        plt.show()

def plot_MODIS_temporary(date_str, zoom = True, save = False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # Read 0.64 μm data
    # --------------- 
    MODIS_data_ch1  = read_MODIS_channel(date_str, 1, zoom = zoom)
    MODIS_data_ch5  = read_MODIS_channel(date_str, 5, zoom = zoom)
    MODIS_data_ch31 = read_MODIS_channel(date_str, 31, zoom = zoom)

    # Determine where the smoke is located
    # ------------------------------------
    hash_data, nohash_data = find_plume(date_str) 

    tmp_data1  = np.copy(MODIS_data_ch1['data'])
    tmp_data5  = np.copy(MODIS_data_ch5['data'])
    tmp_data31 = np.copy(MODIS_data_ch31['data'])
    tmp_lat0   = np.copy(MODIS_data_ch1['lat'])
    tmp_lon0   = np.copy(MODIS_data_ch1['lon'])

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

    plot_scatter(ax0, tmp_data31, tmp_data1, MODIS_data_ch31, MODIS_data_ch1, \
        hash_data, xlabel = '11 μm  temperature', \
        ylabel = '0.64 μm reflectance', plot_legend = True)
    plot_scatter(ax1, tmp_data31, tmp_data5, MODIS_data_ch31, MODIS_data_ch5, \
        hash_data, xlabel = '11 μm brightness temperature', \
        ylabel = '1.24 μm reflectance')

    if(save):
        outname = 'modis_scatter_twopanel_' + date_str + '.png'
        fig1.savefig(outname, dpi = 300)
        print("Saved image",outname)
    else:
        plt.show()

def plot_MODIS_CERES_3panel(zoom = True, show_smoke = True, composite = True, \
        save=False):

    if(home_dir + '/Research/CERES' not in sys.path):
        sys.path.append(home_dir + '/Research/CERES')
    from gridCERESLib import readgridCERES_hrly_grid, plotCERES_hrly
    
    date_str = '202108062025'
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # ----------------------------------------------------------------------
    #
    # Read the CERES data
    #
    # ----------------------------------------------------------------------
    CERES_data_hrly = readgridCERES_hrly_grid(date_str[:10], 'SWF', \
        satellite = 'Aqua', minlat = 20.)

    if(CERES_data_hrly is None):
        print("ERROR: no data returned from readgridCERES_hrly_grid")
        print("Quitting")
        return

    ##!## Read in the CERES data
    ##!## ----------------------
    ##!#mask_LAT, mask_LON, mask_swf, mask_lwf = \
    ##!#    read_CERES_match_MODIS(date_str)
   
    # ----------------------------------------------------------------------
    #
    #  Set up the 9-panel figure
    #
    # ----------------------------------------------------------------------
 
    mapcrs = init_proj(date_str)
    plt.close('all')
    fig = plt.figure(figsize=(11,3.5))
    ax1 = fig.add_subplot(1, 3, 1, projection = mapcrs) # SW   
    ax2 = fig.add_subplot(1, 3, 2, projection = mapcrs) # LW
    ax3 = fig.add_subplot(1, 3, 3, projection = mapcrs) # Total

    # ----------------------------------------------------------------------
    #
    # Plot the data in the figure
    #
    # ----------------------------------------------------------------------
    # Plot the CERES SWF and LWF data
    # -------------------------------
    ##!#plot_CERES_spatial(date_str, mask_LAT, mask_LON, mask_swf, 'SWF', ax1, \
    ##!#    ptitle = '', zoom = zoom)
    ##!#plot_CERES_spatial(date_str, mask_LAT, mask_LON, mask_lwf, 'LWF', ax2, \
    ##!#    ptitle = '', zoom = zoom)
    ##!#plot_CERES_spatial(date_str, mask_LAT, mask_LON, mask_lwf, 'total', ax3, \
    ##!#    ptitle = '', zoom = zoom)

    plotCERES_hrly(ax1, CERES_data_hrly, 'SWF', minlat = 20., \
        vmin = 160, vmax = 325, title = ' ', label = '', \
        circle_bound = False, gridlines = False, grid_data = True)
    plotCERES_hrly(ax2, CERES_data_hrly, 'LWF', minlat = 20., \
        vmin = 300, vmax = 370, title = ' ', label = '', \
        circle_bound = False, gridlines = False, grid_data = True)
    plotCERES_hrly(ax3, CERES_data_hrly, 'total', minlat = 20., \
        vmin = 520, vmax = 625, title = ' ', label = '', \
        circle_bound = False, gridlines = False, grid_data = True)
    if(zoom):
        ax1.set_extent([aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][1], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][1]],\
                        datacrs)
        ax2.set_extent([aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][1], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][1]],\
                        datacrs)
        ax3.set_extent([aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lon'][1], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][0], \
                        aerosol_event_dict[dt_date_str.strftime('%Y-%m-%d')][date_str[8:]]['Lat'][1]],\
                        datacrs)

    ##!## Plot the ASOS locations on the map
    ##!## ----------------------------------
    ##!#plot_ASOS_locs(ax1,date_str,color='lime', sites = ['O05','AAT'])

    # Add subplot labels
    # ------------------
    plot_subplot_label(ax1, '(a)', backgroundcolor = 'white')
    plot_subplot_label(ax2, '(b)', backgroundcolor = 'white')
    plot_subplot_label(ax3, '(c)', backgroundcolor = 'white')

    # Add plot text
    # -------------
    plot_figure_text(ax1, 'CERES SW', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax2, 'CERES LW', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax3, 'CERES SW+LW', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = 12, backgroundcolor = 'white', halign = 'right')

    fig.tight_layout()

    if(save):
        outname = 'modis_ceres_combined_' + date_str + '_3panel.png'
        fig.savefig(outname, dpi=300)
        print("Saved",outname)
    else:
        plt.show()

def plot_MODIS_detection(date_str, zoom = True, save = False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # Set up the figure
    # -----------------
    fig = plt.figure(figsize=(14,8))
    mapcrs = init_proj(date_str)
    #gs1 = fig.add_gridspec(nrows = 3, ncols = 2, wspace = 0.40, hspace = 0.30)
    ax1 = fig.add_subplot(2,3,1, projection = mapcrs) # true color 7/22
    ax2 = fig.add_subplot(2,3,2, projection = mapcrs) # vis 7/22
    ax3 = fig.add_subplot(2,3,3, projection = mapcrs) # IR / vis scatter
    ax4 = fig.add_subplot(2,3,4, projection = mapcrs) # IR / vis scatter
    ax5 = fig.add_subplot(2,3,5, projection = mapcrs) # IR / vis scatter
    ax6 = fig.add_subplot(2,3,6, projection = mapcrs) # IR / vis scatter

    # Read 0.64 μm data
    # --------------- 
    MODIS_data_ch1  = read_MODIS_channel(date_str, 1, zoom = zoom)
    MODIS_data_ch2  = read_MODIS_channel(date_str, 2, zoom = zoom)
    MODIS_data_ch5  = read_MODIS_channel(date_str, 5, zoom = zoom)
    MODIS_data_ch7  = read_MODIS_channel(date_str, 7, zoom = zoom)
    MODIS_data_ch32 = read_MODIS_channel(date_str, 32, zoom = zoom)

    plot_MODIS_spatial(MODIS_data_ch1,  ax1, zoom = zoom)
    plot_MODIS_spatial(MODIS_data_ch2,  ax2, zoom = zoom)
    plot_MODIS_spatial(MODIS_data_ch5,  ax3, zoom = zoom)
    plot_MODIS_spatial(MODIS_data_ch7,  ax4, zoom = zoom)
    plot_MODIS_spatial(MODIS_data_ch32, ax5, zoom = zoom)

    # Play around with cloud and smoke detection
    # ------------------------------------------
    tmp_data = np.full(MODIS_data_ch1['data'].shape, np.nan)
    in_cloud = (((MODIS_data_ch1['data'] + MODIS_data_ch2['data'] > 1.2) |\
            (MODIS_data_ch32['data'] < 265.)) |\
        ((MODIS_data_ch1['data'] + MODIS_data_ch2['data'] > 0.7) &\
            (MODIS_data_ch32['data'] < 285.)))
    tmp_data[in_cloud] = 2.
    tmp_data[~in_cloud] = 0.
    tmp_data[~in_cloud & (MODIS_data_ch1['data'] - \
        MODIS_data_ch7['data'] > 0.05) & (MODIS_data_ch7['data'] > 0.05)] = 1
    mask_tmp_data = np.ma.masked_invalid(tmp_data)

    hash_data = np.ma.masked_where(mask_tmp_data != 1, mask_tmp_data)

    print('numbers = ',int(np.nanmax(mask_tmp_data)) - \
        int(np.nanmin(mask_tmp_data))+2)
    cmap = plt.get_cmap('plasma', int(np.nanmax(mask_tmp_data)) - \
        int(np.nanmin(mask_tmp_data))+1)
    cbar_labels = ['Clear','Smoke','Cloud']

    #tmp_data = MODIS_data_ch1['data'] - MODIS_data_ch7['data'] 
    mesh1 = ax6.pcolormesh(MODIS_data_ch1['lon'],MODIS_data_ch1['lat'],\
        tmp_data,cmap = cmap, vmin = 0,\
        vmax = 2,shading='auto', \
        transform = datacrs) 
        #tmp_data,cmap = 'plasma', shading='auto', transform = datacrs) 
    #mesh1 = ax8.pcolormesh(MODIS_data_ch1['lon'],MODIS_data_ch1['lat'],\
    #    tmp_data,cmap = 'bwr', vmin = -0.25, vmax = 0.25, shading='auto', transform = datacrs) 
    #cbar1 = plt.colorbar(mesh1,ax=ax6,orientation='vertical',\
    #    pad=0.03)
    cbar = plt.colorbar(mesh1,ax = ax6, orientation='vertical',\
        ticks = np.arange(len(cbar_labels)))
    print(cbar_labels[int(np.nanmin(mask_tmp_data)):int(np.nanmax(mask_tmp_data))+1])
    cbar.ax.set_yticklabels(cbar_labels,\
        fontsize=10,weight = 'bold', rotation=0)

    plt.rcParams.update({'hatch.color': 'r'})
    hatch_shape = '\\\\'
    #hash0 = ax2.pcolor(MODIS_data_ch1['lon'],MODIS_data_ch1['lat'],\
    #    hash_data[:MODIS_data_ch1['lat'].shape[0], :MODIS_data_ch1['lat'].shape[1]], hatch = hatch_shape, alpha=0., transform = datacrs)
    ax1.pcolor(MODIS_data_ch1['lon'],MODIS_data_ch1['lat'],\
        hash_data, hatch = hatch_shape, alpha=0.40, transform = datacrs)

    ax6.add_feature(cfeature.BORDERS)
    ax6.add_feature(cfeature.STATES)
    ax6.coastlines()
    if(zoom):
        ax6.set_extent([aerosol_event_dict[MODIS_data_ch1['cross_date']][MODIS_data_ch1['file_time']]['Lon'][0], \
                        aerosol_event_dict[MODIS_data_ch1['cross_date']][MODIS_data_ch1['file_time']]['Lon'][1], \
                        aerosol_event_dict[MODIS_data_ch1['cross_date']][MODIS_data_ch1['file_time']]['Lat'][0], \
                        aerosol_event_dict[MODIS_data_ch1['cross_date']][MODIS_data_ch1['file_time']]['Lat'][1]],\
                        datacrs)

    if(save):
        outname = 'modis_detection_' + date_str + '_6panel.png'
        fig.savefig(outname, dpi=300)
        print("Saved",outname)
    else:
        plt.show()

def plot_scatter_OMI_CERES_figure(zoom = True, show_smoke = False, composite = True, \
        plume_only = False, avg_pixel = True, save=False):

    date_str = '202108062025'
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # ----------------------------------------------------------------------
    #
    # Read the MODIS and CERES data
    #
    # ----------------------------------------------------------------------

    # Call read_MODIS_channel to read the desired MODIS data from the file 
    # and put it in a dictionary
    # ---------------------------------------------------------------------
    MODIS_data_ch31 = read_MODIS_channel(date_str, 31, zoom = zoom)

    # Determine where the smoke is located MAYBE???
    # ---------------------------------------------

    # Read in the CERES data
    # ----------------------
    mask_LAT, mask_LON, mask_swf, mask_lwf = \
        read_CERES_match_MODIS(date_str)
  
    # Read in the OMI data
    # --------------------
    LAT, LON, mask_UVAI = read_OMI_match_MODIS(date_str)
 
 
    # ----------------------------------------------------------------------
    #
    #  Set up the figure
    #
    # ----------------------------------------------------------------------
 
    mapcrs = init_proj(date_str)
    plt.close('all')
    fig = plt.figure(figsize=(6,6))
    ax1 = fig.add_subplot(1,1,1) # Ch 31

    # ----------------------------------------------------------------------
    #
    # Plot the data in the figure
    #
    # ----------------------------------------------------------------------
    plot_scatter_OMI_CERES(date_str, MODIS_data_ch31, ax1, \
        avg_pixel = avg_pixel, plume_only = plume_only, ptitle = '')

    #fig.tight_layout()

    MODIS_data_ch31.clear()

    if(save):
        outname = 'modis_scat_ceres_omi_' + date_str + '.png'
        fig.savefig(outname, dpi=300)
        print("Saved",outname)
    else:
        plt.show()


def plot_CERES_swaths(date_str = '202107222110', \
        show_smoke = False, composite = True, \
        zoom = True, save=False):

    if(home_dir + '/Research/CERES' not in sys.path):
        sys.path.append(home_dir + '/Research/CERES')
    from gridCERESLib import readgridCERES_hrly_grid, plotCERES_hrly

    date_str = '202107222110'
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # ----------------------------------------------------------------------
    #
    # Read the CERES data
    #
    # ----------------------------------------------------------------------
    CERES_data_hrly_swf = readgridCERES_hrly_grid(date_str[:10], 'SWF', \
        minlat = 20., modis_comp = True)

    var1, crs1, lons1, lats1, lat_lims1, lon_lims1, plabel = \
        read_MODIS_satpy(date_str,'true_color',\
    #var1, crs1, lat_lims1, lon_lims1 = read_true_color(date_str,\
        composite=True)

    # ----------------------------------------------------------------------
    #
    #  Set up the 6-panel figure
    #
    # ----------------------------------------------------------------------

    mapcrs = init_proj(date_str)
    plt.close('all')
    fig = plt.figure(figsize=(9.0,3.3))
    gs = fig.add_gridspec(nrows = 1, ncols = 3)
    ax1  = fig.add_subplot(gs[0,0], projection = mapcrs) # true color    
    ax2  = fig.add_subplot(gs[0,1], projection = mapcrs) # MODIS Ch 6
    ax3  = fig.add_subplot(gs[0,2], projection = mapcrs) # MODIS Ch 7

    # Plot the CERES  data
    # --------------------
    if(date_str == '202107222110'):
        sw_vmax = 250
        sw_vmin = 130
        lw_vmax = 370
        lw_vmin = 300
        tot_vmax = 580
        tot_vmin = 460
    elif(date_str == '202108062025'):
        sw_vmax = 330
        sw_vmin = 190
        lw_vmax = 370
        lw_vmin = 300
  
    plotCERES_hrly(ax1, CERES_data_hrly_swf, 'swf', \
        vmin = sw_vmin, vmax = sw_vmax, title = '', label = 'TOA Flux [W/m$^{2}$]', \
        circle_bound = False, gridlines = False, grid_data = True, \
        zoom = True)
    plotCERES_hrly(ax2, CERES_data_hrly_swf, 'lwf', \
        vmin = lw_vmin, vmax = lw_vmax, title = '', label = 'TOA Flux [W/m$^{2}$]', \
        circle_bound = False, gridlines = False, grid_data = True, \
        zoom = True)
    plotCERES_hrly(ax3, CERES_data_hrly_swf, 'total', \
        vmin = tot_vmin, vmax = tot_vmax, title = '', label = 'TOA Flux [W/m$^{2}$]', \
        circle_bound = False, gridlines = False, grid_data = True, \
        zoom = True)
    if(zoom):
        ax1.set_extent([lon_lims1[0],lon_lims1[1],lat_lims1[0],lat_lims1[1]],\
                       datacrs)
        ax2.set_extent([lon_lims1[0],lon_lims1[1],lat_lims1[0],lat_lims1[1]],\
                       datacrs)
        ax3.set_extent([lon_lims1[0],lon_lims1[1],lat_lims1[0],lat_lims1[1]],\
                       datacrs)

    # Add subplot labels
    # ------------------
    font_size = 10
    plot_subplot_label(ax1, '(a)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax2, '(b)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax3, '(c)', backgroundcolor = 'white', fontsize = font_size)

    # Add plot text
    # -------------
    plot_figure_text(ax1, 'CERES SW', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax2, 'CERES LW', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax3, 'CERES Total', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')

    fig.suptitle(dt_date_str.strftime('Aqua CERES imagery of the Dixie Fire\n%Y/%m/%d %H:%M UTC'))

    fig.tight_layout()

    if(save):
            outname = 'ceres_' + date_str + '_fig1.png'
            fig.savefig(outname, dpi=300)
            print("Saved",outname)
    else:
        plt.show()

def plot_viewing_geometry(date_str = '202107222110', zoom = True, show_smoke = False, composite = True, \
        save=False):

    #date_str = '202107202125'
    #date_str = '202107222110'
    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    ##!#if('/home/bsorenson/Research/GOES' not in sys.path):
    ##!#    sys.path.append('/home/bsorenson/Research/GOES')
    ##!#from GOESLib import read_GOES_satpy, plot_GOES_satpy

    # ----------------------------------------------------------------------
    #
    # Read the MODIS and CERES data
    #
    # ----------------------------------------------------------------------

    # Call read_MODIS_channel to read the desired MODIS data from the file 
    # and put it in a dictionary
    # ---------------------------------------------------------------------
    MODIS_data_ch7  = read_MODIS_channel(date_str, 7, zoom = zoom)
    MODIS_data_ch31 = read_MODIS_channel(date_str, 31, zoom = zoom)

    # Read the true color data
    # ------------------------
    var1, crs1, lat_lims1, lon_lims1 = read_true_color(date_str,\
        composite=composite)

    # ----------------------------------------------------------------------
    #
    #  Set up the 6-panel figure
    #
    # ----------------------------------------------------------------------

    mapcrs = init_proj(date_str)
    plt.close('all')
    fig = plt.figure(figsize=(7,7))
    #gs = fig.add_gridspec(nrows = 2, ncols = 8)
    ax1  = fig.add_subplot(2,2,1,  projection = crs1)   # true color    
    ax2  = fig.add_subplot(2,2,2,  projection = mapcrs) # Ch 7
    ax3  = fig.add_subplot(2,2,3,  projection = mapcrs) # Ch 31
    ax4  = fig.add_subplot(2,2,4,  projection = mapcrs)   # GOES vis 


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
    plot_MODIS_spatial(MODIS_data_ch7,  ax2, zoom = zoom, ptitle = '', labelsize = 10, labelticksize = 9)
    plot_MODIS_spatial(MODIS_data_ch31, ax3, zoom = zoom, ptitle = '', labelsize = 10, labelticksize = 9)
    pvar = 'sza'
    plot_MODIS_spatial(MODIS_data_ch31, ax4, zoom = zoom, pvar = pvar, ptitle = '', \
        vmin = np.min(MODIS_data_ch7[pvar]), vmax = np.max(MODIS_data_ch7[pvar]), \
        labelsize = 10, labelticksize = 9)

    labelsize = 10

    # Add subplot labels
    # ------------------
    font_size = 10
    plot_subplot_label(ax1, '(a)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax2, '(b)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax3, '(c)', backgroundcolor = 'white', fontsize = font_size)
    plot_subplot_label(ax4, '(d)', backgroundcolor = 'white', fontsize = font_size)

    # Add plot text
    # -------------
    plot_figure_text(ax1, 'MODIS True Color', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax2, 'MODIS 2.2 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax3, 'MODIS 11.0 μm', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')
    plot_figure_text(ax4, 'MODIS VZA', xval = None, yval = None, transform = None, \
        color = 'red', fontsize = font_size, backgroundcolor = 'white', halign = 'right')

    # Print the viewing geometry values
    # ---------------------------------
    lon_stn = -121.0570
    lat_stn = 40.4186
    c_idx = nearest_gridpoint(lat_stn, lon_stn,\
        MODIS_data_ch7['lat'], MODIS_data_ch7['lon'])

    ax1.plot(lon_stn, lat_stn,
             color='tab:blue', linewidth=2, marker='o',
             transform=datacrs)
    ax2.plot(lon_stn, lat_stn,
             color='tab:blue', linewidth=2, marker='o',
             transform=datacrs)
    ax3.plot(lon_stn, lat_stn,
             color='tab:blue', linewidth=2, marker='o',
             transform=datacrs)
    ax4.plot(lon_stn, lat_stn,
             color='tab:blue', linewidth=2, marker='o',
             transform=datacrs)

    print("SZA")
    print("     - ", MODIS_data_ch7['sza'][c_idx])
    print("VZA")
    print("     - ", MODIS_data_ch7['vza'][c_idx])

    # Calculate the GOES viewing geometry as well
    sys.path.append(home_dir + '/Research/GOES')
    from GOESLib import calc_zenith_angle
    goes_vza = calc_zenith_angle(lat_stn, lon_stn)
    print("")
    print("GOES 17 VZA - ",goes_vza)

    fig.tight_layout()

    MODIS_data_ch7.clear()
    MODIS_data_ch31.clear()

    if(save):
        outname = 'modis_view_geom_' + date_str + '.png'
        fig.savefig(outname, dpi=300)
        print("Saved",outname)
    else:
        plt.show()

# This function shows the locations of ASOS stations around
# the Dixie Fire
# ---------------------------------------------------------
def plot_MODIS_asos_sites(date_str, sites = None, save = False):

    dt_date_str = datetime.strptime(date_str,"%Y%m%d%H%M")

    # ----------------------------------------------------------------------
    #
    # Read the MODIS  data
    #
    # ----------------------------------------------------------------------

    # Call read_MODIS_channel to read the desired MODIS data from the file 
    # and put it in a dictionary
    # ---------------------------------------------------------------------
    var1, crs1, lons1, lats1, lat_lims1, lon_lims1, plabel1 = \
        read_MODIS_satpy(date_str, 'true_color')
    
    # ----------------------------------------------------------------------
    #
    #  Set up the figure
    #
    # ----------------------------------------------------------------------
    plt.close('all')
    fig = plt.figure()
    ax1  = fig.add_subplot(1,1,1, projection = crs1) # true color    

    # ----------------------------------------------------------------------
    #
    # Plot the data in the figure
    #
    # ----------------------------------------------------------------------
    labelsize = 10
    plot_MODIS_satpy(date_str, 'true_color', ax = ax1, var = var1, crs = crs1, \
        lons = lons1, lats = lats1, lat_lims = lat_lims1, lon_lims = lon_lims1, \
        ptitle = '', plabel = plabel1, \
        labelsize = 10, zoom=True, save=False)

    # Plot the ASOS locations on the map
    # ----------------------------------
    plot_ASOS_locs(ax1,date_str,color='lime')

    plt.show()

# pvar = 'cld_frac_mean', 'cld_frac_std', 'ob_count'      
def plot_MODIS_CLDL3(cloud_data, pvar = 'cld_frac_mean', \
        minlat = 65.5, ax = None, colorbar = True, save = False):

    if(isinstance(cloud_data, str)):
        if(len(cloud_data) == 6):
            dt_date_str = datetime.strptime(cloud_data, '%Y%m')
            cloud_data = read_MODIS_CLDL3_daily_allmonth(cloud_data, minlat = minlat)
        elif(len(cloud_data) == 8):
            dt_date_str = datetime.strptime(cloud_data, '%Y%m%d')
            cloud_data = read_MODIS_CLDL3_daily(cloud_data, minlat = minlat)
    else:
        if(len(cloud_data['date_str']) == 6):
            dt_date_str = datetime.strptime(cloud_data['date_str'], '%Y%m')
        elif(len(cloud_data['date_str']) == 8):
            dt_date_str = datetime.strptime(cloud_data['date_str'], '%Y%m%d')
    
    in_ax = True 
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, projection = ccrs.NorthPolarStereo())
    
    mesh = ax.pcolormesh(cloud_data['lon'], cloud_data['lat'], \
        cloud_data[pvar],transform = datacrs, shading = 'auto')
    if(colorbar):
        cbar = plt.colorbar(mesh, ax = ax, \
            extend = 'both')
        cbar.set_label(pvar.replace('_',' '))
    ax.coastlines()
    ax.set_boundary(circle, transform=ax.transAxes)
    ax.set_extent([-180,180,minlat,90], datacrs)

    ax.set_title('MODIS CLDL3 ' + pvar + '\n' + cloud_data['date_str'])

    if(not in_ax):
        if(save):
            outname = 'modis_cldl3_' + pvar + '_' + cloud_data['date_str'] + \
                '.png'
            fig.savefig(outname, dpi = 300)
            print("Saved image", outname)
           
        else:
            plt.show() 

def plot_MODIS_CLDL3_2panel(cloud_data, var1 = 'cld_frac_mean', \
        var2 = 'cld_frac_std',minlat = 65.5, ax = None, colorbar = True, \
        save = False):

    if(isinstance(cloud_data, str)):
        if(len(cloud_data) == 6):
            strfmt = '%Y%m'
            dt_date_str = datetime.strptime(cloud_data, strfmt)
            cloud_data = read_MODIS_CLDL3_single_month(cloud_data)
        elif(len(cloud_data) == 8):
            strfmt = '%Y%m%d'
            dt_date_str = datetime.strptime(cloud_data, strfmt)
            cloud_data = read_MODIS_CLDL3_daily(cloud_data)
    else:
        if(len(cloud_data) == 6):
            strfmt = '%Y%m'
            dt_date_str = datetimme.strptime(cloud_data['date_str'], strfmt)
        elif(len(cloud_data) == 8):
            strfmt = '%Y%m%d'
            dt_date_str = datetimme.strptime(cloud_data['date_str'], strfmt)
    
    fig = plt.figure(figsize = (7,3))
    ax1 = fig.add_subplot(1,2,1, projection = ccrs.NorthPolarStereo())
    ax2 = fig.add_subplot(1,2,2, projection = ccrs.NorthPolarStereo())

    plot_MODIS_CLDL3(cloud_data, pvar = var1, \
        minlat = 65.5, ax = ax1, save = False)
    plot_MODIS_CLDL3(cloud_data, pvar = var2, \
        minlat = 65.5, ax = ax2, save = False)
    outstr = dt_date_str.strftime(strfmt)
    #plt.suptitle(outstr)
    fig.tight_layout()

    if(save):
        outname = '_'.join(['modis','cldl3','2panel',var1, var2,outstr]) + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

def plot_compare_CLDL3_MYD08_trend(cloud_data, myd08_data, month_idx = None, \
        minlat = 65.5, season = '', myd08_var = 'all', ax = None, \
        colorbar = True, trend_type = 'standard', comp_xidx = None, comp_yidx = None, \
        norm_to_decade = True, show_pval = False, save = False):

    fig = plt.figure(figsize = (9.5,7.5))
    ax1 = fig.add_subplot(3,3,1, projection = ccrs.NorthPolarStereo())
    ax2 = fig.add_subplot(3,3,2, projection = ccrs.NorthPolarStereo())
    ax5 = fig.add_subplot(3,3,3, projection = ccrs.NorthPolarStereo())
    ax6 = fig.add_subplot(3,3,4, projection = ccrs.NorthPolarStereo())
    ax7 = fig.add_subplot(3,3,5, projection = ccrs.NorthPolarStereo())
    ax3 = fig.add_subplot(3,3,7)
    ax4 = fig.add_subplot(3,3,8)

    myd08_trend = plotMODIS_MYD08_MonthTrend(myd08_data,month_idx=month_idx,\
        trend_type=trend_type,pvar = myd08_var, season=season,minlat=65.5,\
        colorbar = True, title = 'MYD08', \
        ax = ax1, show_pval = show_pval, uncert_ax = ax6, \
        norm_to_decade = norm_to_decade,vmin = -0.1, vmax = 0.1, \
        return_trend = True)

    cldl3_trend = plotMODIS_CLDL3_MonthTrend(cloud_data,month_idx=month_idx,\
        trend_type=trend_type,season=season,minlat=minlat, \
        colorbar = colorbar, title = 'CLDL3', \
        ax = ax2, show_pval = show_pval, uncert_ax = ax7, \
        norm_to_decade = norm_to_decade,vmin = -0.1, vmax = 0.1, \
        return_trend = True)

    trend_diffs = myd08_trend - cldl3_trend
    mesh = ax5.pcolormesh(cloud_data['lon'], cloud_data['lat'], trend_diffs, \
        transform = datacrs, shading = 'auto', cmap = 'seismic', \
        vmin = -0.05, vmax = 0.05)
    cbar = plt.colorbar(mesh, ax = ax5, pad = 0.03, extend = 'both')
    cbar.set_label('Cloud Fraction Trend Difference')
    ax5.coastlines()
    ax5.set_extent([-180, 180, minlat, 90], datacrs) 
    ax5.set_boundary(circle, transform=ax5.transAxes)
    ax5.set_title('MYD08 - CLDL3')

    xy = np.vstack([cldl3_trend.flatten(), myd08_trend.flatten()])
    z = stats.gaussian_kde(xy)(xy)       
    
    ax3.scatter(cldl3_trend.flatten(), myd08_trend.flatten(), s = 6, c = z)
    ax3.set_xlabel('CLDL3 trend')
    ax3.set_ylabel('MYD08 trend')
    ax3_ylims = ax3.get_ylim()
    ax3_xlims = ax3.get_xlim()
    lims = [
        np.min([ax3_xlims, ax3_ylims]),\
        np.max([ax3_xlims, ax3_ylims])
    ]
    ax3.plot(lims, lims, 'r-', alpha = 0.5)
    ax3.set_xlim(lims)
    ax3.set_ylim(lims)

    r2_val = np.corrcoef(cldl3_trend.flatten(), myd08_trend.flatten())[0,1]**2.
    ax3.set_title('r$^{2}$ = ' + str(np.round(r2_val,3)))

    plt.suptitle(cloud_data['dates'][0] + ' - ' + \
                 cloud_data['dates'][-1])
    

    # Plot time series
    # ----------------
    site_adder = ''
    if((comp_xidx is not None) and (comp_yidx is not None)):
        site_adder = '_' + str(comp_xidx) + 'x' + str(comp_yidx)
        if(month_idx == None):
            idx_jumper = 1
        else:
            idx_jumper = 6
        ax4.plot(np.arange(len(cloud_data['dates'][month_idx::idx_jumper])), \
            myd08_data['cld_frac_mean'][month_idx::idx_jumper,comp_xidx,\
            comp_yidx], label = 'MYD08')
        ax4.plot(np.arange(len(cloud_data['dates'][month_idx::idx_jumper])), \
            cloud_data['cld_frac_mean'][month_idx::idx_jumper,comp_xidx,\
            comp_yidx], label = 'CLDL3')
        ax4.legend()
        plot_point_on_map(ax1, cloud_data['lat'][comp_xidx,comp_yidx], \
            cloud_data['lon'][comp_xidx,comp_yidx], markersize = 10)
        plot_point_on_map(ax2, cloud_data['lat'][comp_xidx,comp_yidx], \
            cloud_data['lon'][comp_xidx,comp_yidx], markersize = 10)
        plot_point_on_map(ax5, cloud_data['lat'][comp_xidx,comp_yidx], \
            cloud_data['lon'][comp_xidx,comp_yidx], markersize = 10)
        plot_point_on_map(ax6, cloud_data['lat'][comp_xidx,comp_yidx], \
            cloud_data['lon'][comp_xidx,comp_yidx], markersize = 10)
        plot_point_on_map(ax7, cloud_data['lat'][comp_xidx,comp_yidx], \
            cloud_data['lon'][comp_xidx,comp_yidx], markersize = 10)

        plot_trend_line(ax4, \
            np.arange(len(myd08_data['dates'][month_idx::idx_jumper])), \
            myd08_data['cld_frac_mean'][\
                month_idx::idx_jumper,comp_xidx,comp_yidx], \
            color='tab:blue')
        plot_trend_line(ax4, \
            np.arange(len(cloud_data['dates'][month_idx::idx_jumper])), \
            cloud_data['cld_frac_mean'][\
                month_idx::idx_jumper,comp_xidx,comp_yidx], \
            color='tab:orange')

        r2_val = np.corrcoef(cloud_data['cld_frac_mean'][\
                month_idx::idx_jumper,comp_xidx,comp_yidx], \
                myd08_data['cld_frac_mean'][\
                month_idx::idx_jumper,comp_xidx,comp_yidx])[0,1]**2.
        ax4.set_title(str(cloud_data['lat'][comp_xidx,comp_yidx]) + 'N, ' + \
                      str(cloud_data['lon'][comp_xidx,comp_yidx]) + 'E, \n'\
                      'r$^{2}$ = ' + str(np.round(r2_val,3)))


        ax4.set_ylabel('Monthly Cloud Fraction')
                    

    #plt.suptitle(outstr)
    fig.tight_layout()

    if(save):
        out_date = cloud_data['dates'][month_idx::idx_jumper][0] + '_' + \
            cloud_data['dates'][month_idx::idx_jumper][-1]
        outname = '_'.join(['modis','cldl3','trend','comp',out_date]) + \
            site_adder + '.png'
        fig.savefig(outname, dpi = 200)
        print("Saved image", outname)
    else:
        plt.show()

# Designed to work with the netCDF data
def plotMODIS_CLDL3_MonthTrend(cloud_data,month_idx=None,save=False,\
        trend_type='standard',season='',minlat=65.,return_trend=False, \
        colorbar = True, colorbar_label_size = None,title = None, \
        ax = None, show_pval = False, uncert_ax = None, \
        norm_to_decade = True, vmin = None, vmax = None):

    trend_label=''
    if(trend_type=='thiel-sen'):
        trend_label='_thielSen'

    if(month_idx == None):
        month_adder = ''
        month_idx = None
        index_jumper = 1
        do_month = False
        #vmax = 0.2
        #vmin = -0.2
    else:
        month_adder = '_month'
        if(cloud_data['season'] == 'sunlight'):
            index_jumper = 6
        else:   
            index_jumper = 12
        do_month = True
        #vmax = 0.20
        #vmin = -0.20

    #lat_ranges = np.arange(minlat,90.5,1.0)
    #lon_ranges = np.arange(-179.5,180.5,1.0)

    # --------------------------------------------------------------
    #
    # Use calcCERES_grid_trend to calculate the trends in the AI data
    #
    # --------------------------------------------------------------
    cloud_trends, cloud_pvals, cloud_uncert = \
        calcMODIS_CLDL3_grid_trend(cloud_data, month_idx, trend_type, \
        minlat, norm_to_decade = norm_to_decade)

    print('Max trend',np.max(cloud_trends))

    if(not show_pval):
        cloud_pvals = None
    else:
        print('month_idx = ',month_idx,' PVAL nanmean = ', \
            np.nanmean(cloud_pvals))

    if(uncert_ax is None):
        cloud_uncert = None
    # --------------------------------------------------------------
    #
    # Plot the calculated trends on a figure
    #
    # --------------------------------------------------------------

    # Set up mapping variables 
    colormap = plt.cm.seismic

    # Pull the beginning and ending dates into datetime objects
    start_date = datetime.strptime(cloud_data['dates'][month_idx::index_jumper][0],'%Y%m')
    end_date   = datetime.strptime(cloud_data['dates'][month_idx::index_jumper][-1],'%Y%m')
    
    # Make figure title
    #date_month = datetime(year = 1,month = month_idx+1, day = 1).strftime('%B')
    month_string = ''
    if(do_month == True):
        month_string = start_date.strftime('%B') + ' '

    if(title is None):
        title = 'MODIS CLDL3 Cloud Fraction\n' + month_string + 'Trends'\
            '\n'+start_date.strftime('%b. %Y') + ' - ' +\
            end_date.strftime('%b. %Y')

    # Call plotCERES_spatial to add the data to the figure

    in_ax = True 
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig1 = plt.figure()
        mapcrs = ccrs.NorthPolarStereo()
        ax = fig1.add_subplot(1,1,1, projection = mapcrs)

    if((vmin is None) and (vmax is None)):
        vmax = np.max(cloud_trends)
        vmin = -vmax

    colormap = plt.cm.seismic
    norm = cm.BoundaryNorm(np.arange(vmin - 0.025, vmax + 0.05, 0.025), colormap.N)
    mesh = ax.pcolormesh(cloud_data['lon'], cloud_data['lat'], cloud_trends, \
        shading = 'auto', transform = datacrs,\
        cmap = colormap,\
        norm = norm)
        #cmap = 'RdYlBu_r')
    if(colorbar):
        cbar = plt.colorbar(mesh, ax = ax, pad = 0.03, extend = 'both')
        cbar.set_label('Cloud Fraction Trend')
    ax.coastlines()
    ax.set_extent([-180, 180, minlat, 90], datacrs)   
    ax.set_boundary(circle, transform=ax.transAxes)
    ax.set_title(title)
    
    if(not in_ax):
        fig1.tight_layout()
        if(save == True):
            month_adder = ''
            if(do_month == True):
                month_adder = '_' + start_date.strftime('%b') 
            out_name = 'cloud_frac_trend'+ month_adder + '_' + \
                start_date.strftime('%Y%m') + '_' + end_date.strftime('%Y%m') + \
                '_min' + str(int(minlat)) + '.png'
            fig1.savefig(out_name,dpi=300)
            print("Saved image",out_name)
        else:
            plt.show()
    ##!#else:
    ##!#    plotCERES_spatial(pax, cloud_data['lat'], cloud_data['lon'], \
    ##!#        cloud_trends, 'trend', ptitle = title, plabel = 'W m$^{-2}$ per study period', \
    ##!#        vmin = v_min, vmax = v_max, colorbar_label_size = colorbar_label_size, \
    ##!#        minlat = minlat)

    if(uncert_ax is not None):
        vmin = 0.0
        vmax = 0.08
        colormap = plt.cm.viridis
        norm = cm.BoundaryNorm(np.arange(vmin, vmax + 0.02, 0.02), colormap.N)
        mesh = uncert_ax.pcolormesh(cloud_data['lon'], cloud_data['lat'], \
            cloud_uncert, \
            shading = 'auto', transform = datacrs,\
            cmap = colormap,\
            norm = norm)
            #cmap = 'RdYlBu_r')
        if(colorbar):
            cbar = plt.colorbar(mesh, ax = uncert_ax, pad = 0.03, extend = 'both')
            cbar.set_label('Cld. Frac. Trend Uncert.')
        uncert_ax.coastlines()
        uncert_ax.set_extent([-180, 180, minlat, 90], datacrs)   
        uncert_ax.set_boundary(circle, transform=uncert_ax.transAxes)
        uncert_ax.set_title(title)

    if(return_trend == True):
        return cloud_trends

# Designed to work with the netCDF data
# pvar = 'day' --> only plot day data
def plotMODIS_MYD08_MonthTrend(cloud_data,month_idx=None,save=False,\
        trend_type='standard',pvar = '', season='',minlat=65.,\
        return_trend=False, colorbar = True, colorbar_label_size = None,\
        title = None, ax = None, show_pval = False, uncert_ax = None, \
        norm_to_decade = True, vmin = None, vmax = None):

    trend_label=''
    if(trend_type=='thiel-sen'):
        trend_label='_thielSen'

    if(month_idx == None):
        month_adder = ''
        month_idx = None
        index_jumper = 1
        do_month = False
        #vmax = 0.2
        #vmin = -0.2
    else:
        month_adder = '_month'
        if(cloud_data['season'] == 'sunlight'):
            index_jumper = 6
        else:   
            index_jumper = 12
        do_month = True
        #vmax = 0.20
        #vmin = -0.20

    #lat_ranges = np.arange(minlat,90.5,1.0)
    #lon_ranges = np.arange(-179.5,180.5,1.0)

    # --------------------------------------------------------------
    #
    # Use calcCERES_grid_trend to calculate the trends in the AI data
    #
    # --------------------------------------------------------------
    cloud_trends, cloud_pvals, cloud_uncert = \
        calcMODIS_MYD08_grid_trend(cloud_data, month_idx, trend_type, \
        minlat, norm_to_decade = norm_to_decade, dtype = pvar)

    print('Max trend',np.max(cloud_trends))

    if(not show_pval):
        cloud_pvals = None
    else:
        print('month_idx = ',month_idx,' PVAL nanmean = ', \
            np.nanmean(cloud_pvals))

    if(uncert_ax is None):
        cloud_uncert = None
    # --------------------------------------------------------------
    #
    # Plot the calculated trends on a figure
    #
    # --------------------------------------------------------------

    # Set up mapping variables 
    colormap = plt.cm.bwr

    # Pull the beginning and ending dates into datetime objects
    start_date = datetime.strptime(cloud_data['dates'][month_idx::index_jumper][0],'%Y%m')
    end_date   = datetime.strptime(cloud_data['dates'][month_idx::index_jumper][-1],'%Y%m')
    
    # Make figure title
    #date_month = datetime(year = 1,month = month_idx+1, day = 1).strftime('%B')
    month_string = ''
    if(do_month == True):
        month_string = start_date.strftime('%B') + ' '

    if(title is None):
        title = 'MODIS MYD08 Cloud Fraction\n' + month_string + 'Trends'\
            '\n'+start_date.strftime('%b. %Y') + ' - ' +\
            end_date.strftime('%b. %Y')

    # Call plotCERES_spatial to add the data to the figure

    in_ax = True 
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig1 = plt.figure()
        mapcrs = ccrs.NorthPolarStereo()
        ax = fig1.add_subplot(1,1,1, projection = mapcrs)

    if((vmin is None) and (vmax is None)):
        vmax = np.max(cloud_trends)
        vmin = -vmax

    colormap = plt.cm.seismic
    norm = cm.BoundaryNorm(np.arange(vmin - 0.025, vmax + 0.05, 0.025), colormap.N)
    mesh = ax.pcolormesh(cloud_data['lon'], cloud_data['lat'], cloud_trends, \
        shading = 'auto', transform = datacrs,\
        cmap = colormap,\
        norm = norm)
        #cmap = 'RdYlBu_r')
    if(colorbar):
        cbar = plt.colorbar(mesh, ax = ax, pad = 0.03, extend = 'both')
        cbar.set_label('Cloud Fraction Trend')
    ax.coastlines()
    ax.set_extent([-180, 180, minlat, 90], datacrs)   
    ax.set_boundary(circle, transform=ax.transAxes)
    ax.set_title(title)

    if(not in_ax):
        fig1.tight_layout()
        if(save == True):
            month_adder = ''
            if(do_month == True):
                month_adder = '_' + start_date.strftime('%b') 
            out_name = 'cloud_frac_trend'+ month_adder + '_' + \
                start_date.strftime('%Y%m') + '_' + end_date.strftime('%Y%m') + \
                '_min' + str(int(minlat)) + '.png'
            fig1.savefig(out_name,dpi=300)
            print("Saved image",out_name)
        else:
            plt.show()
    ##!#else:
    ##!#    plotCERES_spatial(pax, cloud_data['lat'], cloud_data['lon'], \
    ##!#        cloud_trends, 'trend', ptitle = title, plabel = 'W m$^{-2}$ per study period', \
    ##!#        vmin = v_min, vmax = v_max, colorbar_label_size = colorbar_label_size, \
    ##!#        minlat = minlat)

    if(uncert_ax is not None):
        vmin = 0.0
        vmax = 0.08
        colormap = plt.cm.viridis
        norm = cm.BoundaryNorm(np.arange(vmin, vmax + 0.02, 0.02), colormap.N)
        mesh = uncert_ax.pcolormesh(cloud_data['lon'], cloud_data['lat'], \
            cloud_uncert, \
            shading = 'auto', transform = datacrs,\
            cmap = colormap,\
            norm = norm)
        if(colorbar):
            cbar = plt.colorbar(mesh, ax = uncert_ax, pad = 0.03, extend = 'both')
            cbar.set_label('Cld. Frac. Trend Uncert.')
        uncert_ax.coastlines()
        uncert_ax.set_extent([-180, 180, minlat, 90], datacrs)   
        uncert_ax.set_boundary(circle, transform=uncert_ax.transAxes)
        uncert_ax.set_title(title)

    if(return_trend == True):
        return cloud_trends
