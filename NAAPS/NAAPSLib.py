"""
  NAME:
    NAAPSLib

  PURPOSE:

"""

import numpy as np
import numpy.ma as ma
import sys
from netCDF4 import Dataset
from datetime import datetime,timedelta
from dateutil.relativedelta import relativedelta
from scipy.stats import pearsonr,spearmanr, mannwhitneyu
import subprocess
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as cm
from matplotlib.dates import DateFormatter
import matplotlib.gridspec as gridspec
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from satpy import find_files_and_readers
from satpy.scene import Scene
from satpy.writers import get_enhanced_image
from glob import glob
import os

home_dir = os.environ['HOME']
sys.path.append(home_dir)
if(home_dir + '/Research/CERES' not in sys.path):
    sys.path.append(home_dir + '/Research/CERES')
if(home_dir + '/Research/NCEP' not in sys.path):
    sys.path.append(home_dir + '/Research/NCEP')
if(home_dir + '/Research/OMI' not in sys.path):
    sys.path.append(home_dir + '/Research/OMI')
if(home_dir + '/Research/MODIS/obs_smoke_forcing/' not in sys.path):
    sys.path.append(home_dir + '/Research/MODIS/obs_smoke_forcing')

#sys.path.append(os.environ['RESEARCH_PATH'] + '/CERES')
from gridCERESLib import *
from NCEP_Lib import *
from OMILib import *
from MODISLib import *
#sys.path.append('/home/bsorenson/')
#from python_lib import plot_trend_line, plot_subplot_label, plot_figure_text, \
#    nearest_gridpoint, aerosol_event_dict, init_proj, \
#    convert_radiance_to_temp
from python_lib import *

data_dir = home_dir + '/data/NAAPS/'

datacrs = ccrs.PlateCarree()
mapcrs = ccrs.NorthPolarStereo()

event_dict = {
    '20080422': {
        'before_start': '20080407',
        'before_end':   '20080420',
        'start': '2008042100',
        'end':   '2008042500',
        'end_start': '20080426',
        'end_end':   '20080510',
    },
    '20120615': {
        'before_start': '20120609',
        'before_end':   '20120613',
        #'before_start': '20120601',
        #'before_end':   '20120613',
        'start': '2012061400',
        'end':   '2012061900',
        'end_start': '20120620',
        'end_end':   '20120624',
        #'end_start': '20120620',
        #'end_end':   '20120630',
    },
    '20140726': {
        'before_start': '20140721',
        'before_end':   '20140725',
        'start': '2014072600',
        'end':   '2014072800',
        'end_start': '20140729',
        'end_end':   '20140803',
    },
    '20140802': {
        'before_start': '20140727',
        'before_end':   '20140801',
        'start': '2014080200',
        'end':   '2014080800',
        'end_start': '20140808',
        'end_end':   '20140812',
    },
    '20170816': {
        'before_start': '20170811',
        'before_end':   '20170815',
        #'before_start': '20170804',
        #'before_end':   '20170815',
        'start': '2017081600',
        'end':   '2017082200',
        'end_start': '20170823',
        'end_end':   '20170827',
        #'end_start': '20170820',
        #'end_end':   '20170831',
    },
    '20200824': {
        'before_start': '20200810',
        'before_end':   '20200823',
        'start': '2020082400',
        'end':   '2020090800',
        'end_start': '20200909',
        'end_end':   '20200922',
    },

#    '200804221841'
#    '200804222020'
#    '200804222159'
#    '201708161504'
#    '201708161643'
#    '201708161821'
#    '201708171408'
#    '201708171547'
#    '201708171726'
#    '201708171905'
#    '201708172043'
#    '201708181312'
#    '201708181451'
#    '201708181630'
#    '201708181809'
#    '201708181948'
#    '201708191355'
#    '201708191534'
#    '201708191713'
#    '201807051856'
#    '201808241343'
#    '201908102115'
#    '201908102254'
#    '201908110033'
#    '201908110351'
}

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Reading functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

def read_NAAPS(date_str, minlat = 65., dtype = 'no_AI'):

    print(date_str) 
    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H')

    if(dtype == 'no_AI'):
        local_type = 'CLIMO1misr'
    else:
        local_type = 'OMI' 
    filename = dt_date_str.strftime(data_dir + dtype + '/' + \
        '%Y/%Y%m/NVA_' + local_type + '_sfc_conc_sink_%Y%m%d%H.nc')

    data = Dataset(filename, 'r')

    idx1 = np.where(data['lat'][:] >= minlat)[0][0]
    idx2 = np.where(data['lat'][:] >= minlat)[0][-1] + 1

    xx, yy = np.meshgrid(data['lon'][:], data['lat'][idx1:idx2])

    NAAPS_data = {}
    NAAPS_data['filename']       = filename
    NAAPS_data['data_type']      = dtype
    NAAPS_data['date']           = date_str
    NAAPS_data['smoke_conc_sfc'] = data['smoke_conc_sfc'][idx1:idx2,:]
    NAAPS_data['smoke_wetsink']  = data['smoke_wetsink'][idx1:idx2,:]
    NAAPS_data['smoke_drysink']  = data['smoke_drysink'][idx1:idx2,:]
    NAAPS_data['lon']            = xx
    NAAPS_data['lat']            = yy

    over_180    = np.where(NAAPS_data['lon'][0,:] < 0. )
    under_180   = np.where(NAAPS_data['lon'][0,:] > 0.)

    for ii in range(yy.shape[0]):
        NAAPS_data['smoke_conc_sfc'][ii,:] = \
            np.concatenate([NAAPS_data['smoke_conc_sfc'][ii,:][under_180],\
            NAAPS_data['smoke_conc_sfc'][ii,:][over_180]])
        NAAPS_data['smoke_wetsink'][ii,:] = \
            np.concatenate([NAAPS_data['smoke_wetsink'][ii,:][under_180],\
            NAAPS_data['smoke_wetsink'][ii,:][over_180]])
        NAAPS_data['smoke_drysink'][ii,:] = \
            np.concatenate([NAAPS_data['smoke_drysink'][ii,:][under_180],\
            NAAPS_data['smoke_drysink'][ii,:][over_180]])
        NAAPS_data['lat'][ii,:] = \
            np.concatenate([NAAPS_data['lat'][ii,:][under_180],\
            NAAPS_data['lat'][ii,:][over_180]])
        NAAPS_data['lon'][ii,:] = \
            np.concatenate([NAAPS_data['lon'][ii,:][under_180],\
            NAAPS_data['lon'][ii,:][over_180] + 360.])

    data.close()
 
    return NAAPS_data

def read_NAAPS_event(date_str, minlat = 65., dtype = 'no_AI'):
 
    print("HERE2:",date_str)
    dt_date_str = datetime.strptime(date_str, '%Y%m%d')

    # Get the event starting and ending times
    # ---------------------------------------
    begin_str = event_dict[date_str]['start']
    end_str   = event_dict[date_str]['end']

    dt_begin_str = datetime.strptime(begin_str, '%Y%m%d%H')
    dt_end_str   = datetime.strptime(end_str, '%Y%m%d%H')

    # Grab the files from the data directory
    # --------------------------------------
    files = subprocess.check_output(dt_begin_str.strftime('ls ' + data_dir + \
        dtype + '/%Y/%Y%m/*.nc'), \
        shell = True).decode('utf-8').strip().split('\n')

    # Account for events that cross month boundaries
    # ----------------------------------------------
    if(dt_begin_str.month != dt_end_str.month):
        files = files + subprocess.check_output(dt_end_str.strftime('ls ' + \
            data_dir + dtype + '/%Y/%Y%m/*.nc'), \
            shell = True).decode('utf-8').strip().split('\n')

    # Use the actual file times to get timestamps
    # -------------------------------------------
    file_dates = np.array([datetime.strptime(tfile[-13:-3],'%Y%m%d%H') \
        for tfile in files]) 

    # Figure out where the closest files to the desired time are
    # ----------------------------------------------------------
    in_time = np.where( (file_dates >= dt_begin_str) & \
        (file_dates <= dt_end_str))

    # Select those files. Should result in all 16 channels for a single
    # time
    # -----------------------------------------------------------------
    dates_found = file_dates[in_time]

    # Set up the first data object
    # ----------------------------
    NAAPS_data = read_NAAPS(dates_found[0].strftime('%Y%m%d%H'), \
        minlat = minlat, dtype = dtype) 
    NAAPS_data['dates'] = [NAAPS_data['date']]

    # Loop over the remaining objects
    # -------------------------------
    for ttime in dates_found[1:]:
        local_data = read_NAAPS(ttime.strftime('%Y%m%d%H'), \
            minlat = minlat, dtype = dtype)
        NAAPS_data['smoke_conc_sfc'] += local_data['smoke_conc_sfc']
        NAAPS_data['smoke_wetsink']  += local_data['smoke_wetsink']
        NAAPS_data['smoke_drysink']  += local_data['smoke_drysink']

        NAAPS_data['dates'].append(ttime.strftime('%Y%m%d%H'))

    NAAPS_data['dt_begin_date'] = dt_begin_str
    NAAPS_data['dt_end_date']   = dt_end_str
    
    return NAAPS_data

# NAAPS data is only passed and returned so that, if desired,
# the NAAPS data can be masked in the regions where the CERES data
# are masked for not having high enough ice or smoke.
def read_zoom_CERES_data_single(date1, date2, ceres_var, min_ice, \
        min_smoke, max_smoke, vmin2, NAAPS_data, minlat, lat_bounds, \
        lon_bounds, satellite = 'All', mask_NAAPS = False):

    CERES_data1 = readgridCERES_daily(date1,end_str = date2, \
        satellite = satellite, minlat = minlat)
    # Test masking the data
    # ---------------------
    avg_CERES1     = np.nanmean(CERES_data1[ceres_var], axis = 0)
    avg_CERES1_alb = np.nanmean(CERES_data1['alb_clr'], axis = 0)
    avg_CERES1_ice = np.nanmean(CERES_data1['ice_conc'], axis = 0)
    #avg_CERES12 = np.nanmean(CERES_data12[ceres_var], axis = 0)
    #avg_CERES22 = np.nanmean(CERES_data22[ceres_var], axis = 0)
    #avg_CERES12_ice = np.nanmean(CERES_data12['ice_conc'], axis = 0)

    # Test pulling out the data that only exists in both zones
    # --------------------------------------------------------
    #avg_CERES1 = np.ma.masked_where((avg_CERES1_ice <= min_ice) | \
    if(lon_bounds[0] < lon_bounds[1]):

        # Contains the smoky data
        CERES_data1[ceres_var] = np.ma.masked_where((\
            avg_CERES1_ice <= min_ice) | \
            (NAAPS_data['smoke_conc_sfc'] < min_smoke) | \
            (NAAPS_data['smoke_conc_sfc'] > max_smoke), avg_CERES1)
        CERES_data1[ceres_var] = np.ma.masked_where(~ \
            ((CERES_data1['lon'] >= lon_bounds[0]) & \
            (CERES_data1['lon'] <= lon_bounds[1]) & \
            (CERES_data1['lat'] >= lat_bounds[0]) & \
            (CERES_data1['lat'] <= lat_bounds[1])) , CERES_data1[ceres_var])
        CERES_data1[ceres_var] = np.ma.masked_where(\
            avg_CERES1_alb < vmin2, CERES_data1[ceres_var])
        # Contains the non-smoky data
        CERES_data1[ceres_var+'_nosmoke'] = np.ma.masked_where((\
            avg_CERES1_ice <= min_ice) | \
            (NAAPS_data['smoke_conc_sfc'] > min_smoke), avg_CERES1)
        CERES_data1[ceres_var+'_nosmoke'] = np.ma.masked_where(~ \
            ((CERES_data1['lon'] >= lon_bounds[0]) & \
            (CERES_data1['lon'] <= lon_bounds[1]) & \
            (CERES_data1['lat'] >= lat_bounds[0]) & \
            (CERES_data1['lat'] <= lat_bounds[1])) , \
            CERES_data1[ceres_var+'_nosmoke'])
        CERES_data1[ceres_var+'_nosmoke'] = np.ma.masked_where(\
            avg_CERES1_alb < vmin2, CERES_data1[ceres_var+'_nosmoke'])

        if(mask_NAAPS):
           # Smoky data
           NAAPS_data['smoke_conc_sfc'] = np.ma.masked_where(\
               (avg_CERES1_ice <= min_ice) | \
               (NAAPS_data['smoke_conc_sfc'] < min_smoke) | \
               (NAAPS_data['smoke_conc_sfc'] > max_smoke), \
               NAAPS_data['smoke_conc_sfc'])
           NAAPS_data['smoke_conc_sfc'] = np.ma.masked_where(~ \
               ((CERES_data1['lon'] >= lon_bounds[0]) & \
               (CERES_data1['lon'] <= lon_bounds[1]) & \
               (CERES_data1['lat'] >= lat_bounds[0]) & \
               (CERES_data1['lat'] <= lat_bounds[1])), \
               NAAPS_data['smoke_conc_sfc'])
           # Non-Smoky data

        #NAAPS_data['smoke_conc_sfc'] = np.ma.masked_where(CERES_data1['alb_clr'] < vmin2, NAAPS_data['smoke_conc_sfc'])

    else:
        # Contains the smoky data
        CERES_data1[ceres_var] = np.ma.masked_where((\
            avg_CERES1_ice <= min_ice) | \
            (NAAPS_data['smoke_conc_sfc'] < min_smoke) | \
            (NAAPS_data['smoke_conc_sfc'] > max_smoke), avg_CERES1)
        CERES_data1[ceres_var] = np.ma.masked_where(~\
            (((CERES_data1['lon'] >= lon_bounds[0]) | \
            (CERES_data1['lon'] <= lon_bounds[1])) & \
            (CERES_data1['lat'] >= lat_bounds[0]) & \
            (CERES_data1['lat'] <= lat_bounds[1])) , CERES_data1[ceres_var])
        CERES_data1[ceres_var] = np.ma.masked_where(\
            CERES_data1['alb_clr'] < vmin2, CERES_data1[ceres_var])
        # Contains the non-smoky data
        CERES_data1[ceres_var+'_nosmoke'] = np.ma.masked_where((\
            avg_CERES1_ice <= min_ice) | \
            (NAAPS_data['smoke_conc_sfc'] > min_smoke), avg_CERES1)
        CERES_data1[ceres_var+'_nosmoke'] = np.ma.masked_where(~\
            (((CERES_data1['lon'] >= lon_bounds[0]) | \
            (CERES_data1['lon'] <= lon_bounds[1])) & \
            (CERES_data1['lat'] >= lat_bounds[0]) & \
            (CERES_data1['lat'] <= lat_bounds[1])) , CERES_data1[ceres_var+'_nosmoke'])
        CERES_data1[ceres_var+'_nosmoke'] = np.ma.masked_where(\
            CERES_data1['alb_clr'] < vmin2, CERES_data1[ceres_var+'_nosmoke'])

        if(mask_NAAPS):
            NAAPS_data['smoke_conc_sfc'] = np.ma.masked_where(\
                (avg_CERES1_ice <= min_ice) | \
                (NAAPS_data['smoke_conc_sfc'] < min_smoke) | \
                (NAAPS_data['smoke_conc_sfc'] > max_smoke), \
                NAAPS_data['smoke_conc_sfc'])
            NAAPS_data['smoke_conc_sfc'] = np.ma.masked_where(~ (((CERES_data1['lon'] >= lon_bounds[0]) | \
                (CERES_data1['lon'] <= lon_bounds[1])) & (CERES_data1['lat'] >= lat_bounds[0]) & \
                (CERES_data1['lat'] <= lat_bounds[1])), NAAPS_data['smoke_conc_sfc'])
            #NAAPS_data['smoke_conc_sfc'] = np.ma.masked_where(CERES_data1['alb_clr'] < vmin2, NAAPS_data['smoke_conc_sfc'])

    return CERES_data1, NAAPS_data


# NAAPS data is only passed and returned so that, if desired,
# the NAAPS data can be masked in the regions where the CERES data
# are masked for not having high enough ice or smoke.
def read_zoom_CERES_data(date1, date2, date3, date4, ceres_var, min_ice, \
        min_smoke, max_smoke, vmin2, NAAPS_data, minlat, lat_bounds, \
        lon_bounds, satellite = 'All', mask_NAAPS = False):

    CERES_data1 = readgridCERES_daily(date1,end_str = date2, \
        satellite = satellite, minlat = minlat)
    CERES_data2 = readgridCERES_daily(date3,end_str = date4, \
        satellite = satellite, minlat = minlat)

    # Test masking the data
    # ---------------------
    avg_CERES1     = np.nanmean(CERES_data1[ceres_var], axis = 0)
    avg_CERES2     = np.nanmean(CERES_data2[ceres_var], axis = 0)
    avg_CERES1_alb = np.nanmean(CERES_data1['alb_clr'], axis = 0)
    avg_CERES2_alb = np.nanmean(CERES_data2['alb_clr'], axis = 0)
    avg_CERES1_ice = np.nanmean(CERES_data1['ice_conc'], axis = 0)
    #avg_CERES12 = np.nanmean(CERES_data12[ceres_var], axis = 0)
    #avg_CERES22 = np.nanmean(CERES_data22[ceres_var], axis = 0)
    #avg_CERES12_ice = np.nanmean(CERES_data12['ice_conc'], axis = 0)

    # Test pulling out the data that only exists in both zones
    # --------------------------------------------------------
    #avg_CERES1 = np.ma.masked_where((avg_CERES1_ice <= min_ice) | \
    if(lon_bounds[0] < lon_bounds[1]):

        # Contains the smoky data
        CERES_data1[ceres_var] = np.ma.masked_where((\
            avg_CERES1_ice <= min_ice) | \
            (NAAPS_data['smoke_conc_sfc'] < min_smoke) | \
            (NAAPS_data['smoke_conc_sfc'] > max_smoke), avg_CERES1)
        CERES_data1[ceres_var] = np.ma.masked_where(~ \
            ((CERES_data1['lon'] >= lon_bounds[0]) & \
            (CERES_data1['lon'] <= lon_bounds[1]) & \
            (CERES_data1['lat'] >= lat_bounds[0]) & \
            (CERES_data1['lat'] <= lat_bounds[1])) , CERES_data1[ceres_var])
        CERES_data1[ceres_var] = np.ma.masked_where(\
            avg_CERES1_alb < vmin2, CERES_data1[ceres_var])
        # Contains the non-smoky data
        CERES_data1[ceres_var+'_nosmoke'] = np.ma.masked_where((\
            avg_CERES1_ice <= min_ice) | \
            (NAAPS_data['smoke_conc_sfc'] > min_smoke), avg_CERES1)
        CERES_data1[ceres_var+'_nosmoke'] = np.ma.masked_where(~ \
            ((CERES_data1['lon'] >= lon_bounds[0]) & \
            (CERES_data1['lon'] <= lon_bounds[1]) & \
            (CERES_data1['lat'] >= lat_bounds[0]) & \
            (CERES_data1['lat'] <= lat_bounds[1])) , \
            CERES_data1[ceres_var+'_nosmoke'])
        CERES_data1[ceres_var+'_nosmoke'] = np.ma.masked_where(\
            avg_CERES1_alb < vmin2, CERES_data1[ceres_var+'_nosmoke'])

        # Contains the smoky data
        CERES_data2[ceres_var] = np.ma.masked_where((\
            avg_CERES1_ice <= min_ice) | \
            (NAAPS_data['smoke_conc_sfc'] < min_smoke) | \
            (NAAPS_data['smoke_conc_sfc'] > max_smoke), avg_CERES2)
        CERES_data2[ceres_var] = np.ma.masked_where(~\
            ((CERES_data2['lon'] >= lon_bounds[0]) & \
            (CERES_data2['lon'] <= lon_bounds[1]) & \
            (CERES_data2['lat'] >= lat_bounds[0]) & \
            (CERES_data2['lat'] <= lat_bounds[1])) , CERES_data2[ceres_var])
        CERES_data2[ceres_var] = np.ma.masked_where(\
            avg_CERES2_alb < vmin2, CERES_data2[ceres_var])
        # Contains the non-smoky data
        CERES_data2[ceres_var+'_nosmoke'] = np.ma.masked_where((\
            avg_CERES1_ice <= min_ice) | \
            (NAAPS_data['smoke_conc_sfc'] > min_smoke), avg_CERES2)
        CERES_data2[ceres_var+'_nosmoke'] = np.ma.masked_where(~\
            ((CERES_data2['lon'] >= lon_bounds[0]) & \
            (CERES_data2['lon'] <= lon_bounds[1]) & \
            (CERES_data2['lat'] >= lat_bounds[0]) & \
            (CERES_data2['lat'] <= lat_bounds[1])) , \
            CERES_data2[ceres_var+'_nosmoke'])
        CERES_data2[ceres_var+'_nosmoke'] = np.ma.masked_where(\
            avg_CERES2_alb < vmin2, CERES_data2[ceres_var+'_nosmoke'])

        if(mask_NAAPS):
           # Smoky data
           NAAPS_data['smoke_conc_sfc'] = np.ma.masked_where(\
               (avg_CERES1_ice <= min_ice) | \
               (NAAPS_data['smoke_conc_sfc'] < min_smoke) | \
               (NAAPS_data['smoke_conc_sfc'] > max_smoke), \
               NAAPS_data['smoke_conc_sfc'])
           NAAPS_data['smoke_conc_sfc'] = np.ma.masked_where(~ \
               ((CERES_data1['lon'] >= lon_bounds[0]) & \
               (CERES_data1['lon'] <= lon_bounds[1]) & \
               (CERES_data1['lat'] >= lat_bounds[0]) & \
               (CERES_data1['lat'] <= lat_bounds[1])), \
               NAAPS_data['smoke_conc_sfc'])
           # Non-Smoky data

        #NAAPS_data['smoke_conc_sfc'] = np.ma.masked_where(CERES_data1['alb_clr'] < vmin2, NAAPS_data['smoke_conc_sfc'])

    else:
        # Contains the smoky data
        CERES_data1[ceres_var] = np.ma.masked_where((\
            avg_CERES1_ice <= min_ice) | \
            (NAAPS_data['smoke_conc_sfc'] < min_smoke) | \
            (NAAPS_data['smoke_conc_sfc'] > max_smoke), avg_CERES1)
        CERES_data1[ceres_var] = np.ma.masked_where(~\
            (((CERES_data1['lon'] >= lon_bounds[0]) | \
            (CERES_data1['lon'] <= lon_bounds[1])) & \
            (CERES_data1['lat'] >= lat_bounds[0]) & \
            (CERES_data1['lat'] <= lat_bounds[1])) , CERES_data1[ceres_var])
        CERES_data1[ceres_var] = np.ma.masked_where(\
            CERES_data1['alb_clr'] < vmin2, CERES_data2[ceres_var])
        # Contains the non-smoky data
        CERES_data1[ceres_var+'_nosmoke'] = np.ma.masked_where((\
            avg_CERES1_ice <= min_ice) | \
            (NAAPS_data['smoke_conc_sfc'] > min_smoke), avg_CERES1)
        CERES_data1[ceres_var+'_nosmoke'] = np.ma.masked_where(~\
            (((CERES_data1['lon'] >= lon_bounds[0]) | \
            (CERES_data1['lon'] <= lon_bounds[1])) & \
            (CERES_data1['lat'] >= lat_bounds[0]) & \
            (CERES_data1['lat'] <= lat_bounds[1])) , CERES_data1[ceres_var+'_nosmoke'])
        CERES_data1[ceres_var+'_nosmoke'] = np.ma.masked_where(\
            CERES_data1['alb_clr'] < vmin2, CERES_data2[ceres_var+'_nosmoke'])

        # Contains the smoky data
        CERES_data2[ceres_var] = np.ma.masked_where((\
            avg_CERES1_ice <= min_ice) | \
            (NAAPS_data['smoke_conc_sfc'] < min_smoke) | \
            (NAAPS_data['smoke_conc_sfc'] > max_smoke), avg_CERES2)
        CERES_data2[ceres_var] = np.ma.masked_where(~\
            (((CERES_data2['lon'] >= lon_bounds[0]) | \
            (CERES_data2['lon'] <= lon_bounds[1])) & \
            (CERES_data2['lat'] >= lat_bounds[0]) & \
            (CERES_data2['lat'] <= lat_bounds[1])) , CERES_data2[ceres_var])
        CERES_data2[ceres_var] = np.ma.masked_where(\
            CERES_data2['alb_clr'] < vmin2, CERES_data2[ceres_var])
        # Contains the non-smoky data
        CERES_data2[ceres_var+'_nosmoke'] = np.ma.masked_where((\
            avg_CERES1_ice <= min_ice) | \
            (NAAPS_data['smoke_conc_sfc'] > min_smoke), avg_CERES2)
        CERES_data2[ceres_var+'_nosmoke'] = np.ma.masked_where(~\
            (((CERES_data2['lon'] >= lon_bounds[0]) | \
            (CERES_data2['lon'] <= lon_bounds[1])) & \
            (CERES_data2['lat'] >= lat_bounds[0]) & \
            (CERES_data2['lat'] <= lat_bounds[1])) , \
            CERES_data2[ceres_var+'_nosmoke'])
        CERES_data2[ceres_var+'_nosmoke'] = np.ma.masked_where(\
            CERES_data2['alb_clr'] < vmin2, CERES_data2[ceres_var+'_nosmoke'])


        if(mask_NAAPS):
            NAAPS_data['smoke_conc_sfc'] = np.ma.masked_where(\
                (avg_CERES1_ice <= min_ice) | \
                (NAAPS_data['smoke_conc_sfc'] < min_smoke) | \
                (NAAPS_data['smoke_conc_sfc'] > max_smoke), \
                NAAPS_data['smoke_conc_sfc'])
            NAAPS_data['smoke_conc_sfc'] = np.ma.masked_where(~ (((CERES_data1['lon'] >= lon_bounds[0]) | \
                (CERES_data1['lon'] <= lon_bounds[1])) & (CERES_data1['lat'] >= lat_bounds[0]) & \
                (CERES_data1['lat'] <= lat_bounds[1])), NAAPS_data['smoke_conc_sfc'])
            #NAAPS_data['smoke_conc_sfc'] = np.ma.masked_where(CERES_data1['alb_clr'] < vmin2, NAAPS_data['smoke_conc_sfc'])

    return CERES_data1, CERES_data2, NAAPS_data

# plot_daily_data: plots a 5-panel figure showing the total
# smoke deposition, the daily CERES variable values within
# the smoky region for the before and after period for that
# year, and the smoky and non-smoky region values for the
# before period.
def read_CERES_historic(dt_begin_str1, dt_end_str1, dt_begin_str2, \
        dt_end_str2, ceres_var, min_ice, min_smoke, max_smoke, vmin2, \
        satellite, NAAPS_data, minlat, lat_bounds, lon_bounds, \
        plot_daily_data = False):

    combined_data = {}
    combined_data_nosmoke = {}
    average_data = np.zeros(len(np.arange(-4, 5)) * 2)
    cnt = 0

    # Check the year ranges
    # ---------------------
    final_end_date = datetime(year = 2021, month = 9, day = 30)
    check_end_date = dt_end_str2 + relativedelta(years = 4)
    if(check_end_date > final_end_date):
        end_year_offset = check_end_date.year - final_end_date.year
        end_idx = 5 - end_year_offset
        beg_idx = -4 - end_year_offset
    else:
        beg_idx = -4
        end_idx = 5
        

    for ii in np.arange(beg_idx, end_idx):
        dt_begin_local1 = dt_begin_str1 + relativedelta(years = ii)
        dt_end_local1   = dt_end_str1 + relativedelta(years = ii)
        dt_begin_local2 = dt_begin_str2 + relativedelta(years = ii)
        dt_end_local2   = dt_end_str2 + relativedelta(years = ii)

        combined_data[dt_begin_local1.strftime('%Y')] = {}
        combined_data_nosmoke[dt_begin_local1.strftime('%Y')] = {}

        print(dt_begin_local1, dt_end_local1, dt_begin_local2, dt_end_local2)

        # Used for the historical purposes
        CERES_data1, CERES_data2, NAAPS_data = read_zoom_CERES_data(\
            dt_begin_local1.strftime('%Y%m%d'), \
            dt_end_local1.strftime('%Y%m%d'), \
            dt_begin_local2.strftime('%Y%m%d'), \
            dt_end_local2.strftime('%Y%m%d'), \
            ceres_var, min_ice, min_smoke, max_smoke, vmin2, NAAPS_data, minlat, \
            lat_bounds, lon_bounds, \
            satellite = satellite, mask_NAAPS = False)

        # Used for the single day plotting purposes
        CERES_data12 = readgridCERES_daily(dt_begin_local1.strftime('%Y%m%d'),\
            end_str = dt_end_local1.strftime('%Y%m%d'), \
            satellite = satellite, minlat = minlat)
        CERES_data22 = readgridCERES_daily(dt_begin_local2.strftime('%Y%m%d'),\
            end_str = dt_end_local2.strftime('%Y%m%d'), \
            satellite = satellite, minlat = minlat)
        avg_CERES1 = np.nanmean(CERES_data12[ceres_var], axis = 0)
        avg_CERES2 = np.nanmean(CERES_data22[ceres_var], axis = 0)
        CERES_data12[ceres_var] = np.ma.masked_where(~ (\
            (CERES_data12['lon'] >= lon_bounds[0]) & \
            (CERES_data12['lon'] <= lon_bounds[1]) & \
            (CERES_data12['lat'] >= lat_bounds[0]) & \
            (CERES_data12['lat'] <= lat_bounds[1])) , avg_CERES1)
        CERES_data22[ceres_var] = np.ma.masked_where(~ (\
            (CERES_data22['lon'] >= lon_bounds[0]) & \
            (CERES_data22['lon'] <= lon_bounds[1]) & \
            (CERES_data22['lat'] >= lat_bounds[0]) & \
            (CERES_data22['lat'] <= lat_bounds[1])) , avg_CERES2)

        print(dt_begin_local1.strftime('%Y%m%d'), \
            CERES_data1[ceres_var].compressed().shape, \
            CERES_data1[ceres_var + '_nosmoke'].compressed().shape)

        bdata1 = CERES_data1[ceres_var].flatten().compressed()
        bdata2 = CERES_data2[ceres_var].flatten().compressed()
        bdata3 = CERES_data1[ceres_var+'_nosmoke'].flatten().compressed()
        bdata4 = CERES_data2[ceres_var+'_nosmoke'].flatten().compressed()
  
        combined_data[dt_begin_local1.strftime('%Y')][\
            dt_begin_local1.strftime('%m/%d') + ' - ' + \
            dt_end_local1.strftime('%m/%d')] = bdata1
        combined_data[dt_begin_local2.strftime('%Y')][\
            dt_begin_local2.strftime('%m/%d') + ' - ' + \
            dt_end_local2.strftime('%m/%d')] = bdata2

        combined_data_nosmoke[dt_begin_local1.strftime('%Y')][\
            dt_begin_local1.strftime('%m/%d') + ' - ' + \
            dt_end_local1.strftime('%m/%d')] = bdata3
        combined_data_nosmoke[dt_begin_local2.strftime('%Y')][\
            dt_begin_local2.strftime('%m/%d') + ' - ' + \
            dt_end_local2.strftime('%m/%d')] = bdata4

        # Add the average data
        # --------------------
        average_data[cnt * 2]     = np.nanmean(bdata1)
        average_data[cnt * 2 + 1] = np.nanmean(bdata2)

        cnt += 1

        if(plot_daily_data):
            # Make a figure of the local data
            # -------------------------------
            fig3 = plt.figure(figsize = (9,7))
            ax10 = fig3.add_subplot(2,3,1, projection = mapcrs)
            ax11 = fig3.add_subplot(2,3,2, projection = mapcrs)
            ax12 = fig3.add_subplot(2,3,3, projection = mapcrs)
            ax21 = fig3.add_subplot(2,3,5, projection = mapcrs)
            ax22 = fig3.add_subplot(2,3,6, projection = mapcrs)

            plot_NAAPS(NAAPS_data, 'smoke_conc_sfc', ax = ax10, zoom = True, \
                    minlat = minlat, vmin = 0, vmax = 50, plot_log = False)

            # Plot the CERES data
            #plotCERES_daily(CERES_data1, ceres_var, end_str = dt_end_local1.strftime('%Y%m%d'), \
            plotCERES_daily(CERES_data12, ceres_var, end_str = \
                dt_end_local1.strftime('%Y%m%d'), satellite = satellite,  \
                only_sea_ice = False, minlat = minlat, \
                vmin = vmin2, vmax = 0.7, \
                avg_data = True, ax = ax11, save = False, min_ice = min_ice, \
                circle_bound = True, colorbar = True)
            #plotCERES_daily(CERES_data2, ceres_var, end_str = dt_end_local2.strftime('%Y%m%d'), \
            plotCERES_daily(CERES_data22, ceres_var, end_str = \
                dt_end_local2.strftime('%Y%m%d'), satellite = satellite,  \
                only_sea_ice = False, minlat = minlat, vmin = vmin2, \
                vmax = 0.7, \
                avg_data = True, ax = ax12, save = False, min_ice = min_ice, \
                circle_bound = True, colorbar = True)
            plotCERES_daily(CERES_data1, ceres_var, \
                end_str = dt_end_local1.strftime('%Y%m%d'), \
                satellite = satellite,  only_sea_ice = False, minlat = minlat, \
                vmin = vmin2, vmax = 0.7, \
                avg_data = True, ax = ax21, save = False, min_ice = min_ice, \
                circle_bound = True, colorbar = True)
            #plotCERES_daily(CERES_data2, ceres_var, end_str = dt_end_local2.strftime('%Y%m%d'), \
            plotCERES_daily(CERES_data1, ceres_var + '_nosmoke', \
                end_str = dt_end_local1.strftime('%Y%m%d'), \
                satellite = satellite,  only_sea_ice = False, minlat = minlat, \
                vmin = vmin2, vmax = 0.7, \
                avg_data = True, ax = ax22, save = False, min_ice = min_ice, \
                circle_bound = True, colorbar = True)

            ax10.coastlines()
            ax10.set_title('NAAPS-RA smoke_conc_sfc\n' + \
                NAAPS_data['dt_begin_date'].strftime('%Y-%m-%d') + ' - ' + \
                NAAPS_data['dt_end_date'].strftime('%Y-%m-%d'))
            ax11.set_title('CERES Average Clear-sky Albedo\n' + \
                dt_begin_local1.strftime('%Y%m%d') + ' - ' + \
                dt_end_local1.strftime('%Y%m%d'))
            ax12.set_title('CERES Average Clear-sky Albedo\n' + \
                dt_begin_local2.strftime('%Y%m%d') + ' - ' + \
                dt_end_local2.strftime('%Y%m%d'))
            ax21.set_title('CERES Average Clear-sky Albedo\n' + \
                dt_begin_local1.strftime('%Y%m%d') + ' - ' + \
                dt_end_local1.strftime('%Y%m%d'))
            ax22.set_title('CERES Average Clear-sky Albedo\n(No Smoke)\n' + \
                dt_begin_local1.strftime('%Y%m%d') + ' - ' + \
                dt_end_local1.strftime('%Y%m%d'))
            fig3.tight_layout()
            outname = 'naaps_ceres_historic_' + \
                dt_end_local1.strftime('%Y%m%d') + '.png'
            fig3.savefig(outname, dpi = 300)
            print("Saved image", outname)


    second_labels = np.array([[bkey for akey in \
        combined_data[bkey].keys()] for bkey in combined_data.keys()])
    labels = second_labels.flatten()
    second_arrays = np.array([[combined_data[bkey][akey] for akey in \
        combined_data[bkey].keys()] for bkey in combined_data.keys()])
    second_arrays_nosmoke = np.array([[\
        combined_data_nosmoke[bkey][akey] for akey in \
        combined_data_nosmoke[bkey].keys()] for bkey in \
        combined_data_nosmoke.keys()])

    return second_labels, second_arrays, average_data, second_arrays_nosmoke

# dt_begin_str = beginning of time window for traveling averaging
# dt_end_str   = end of time window for traveling averaging
# interval     = length of time window for averaging (in days)
# ceres_var    = 'alb_clr', 'swf_clr', 'lwf_clr'
# min_ice      = minimum NSIDC ice threshold for analysis
# plats        = point lats for pulling data at a single grid point
# plons        = point lons for pulling data at a single grid point
def read_CERES_region_time_series(dt_begin_str, dt_end_str, interval, \
        ceres_var, min_ice, min_smoke, max_smoke, vmin2, satellite, \
        NAAPS_data, minlat, lat_bounds, lon_bounds, \
        plot_daily_data = False):

    # Find the number of intervals between the two dates
    # --------------------------------------------------
    num_times = (dt_end_str - dt_begin_str).days - (interval - 1)

    print('num_times = ',num_times)

    combined_data = {}
    combined_data_nosmoke = {}
    average_data = np.zeros(num_times)
    cnt = 0

    for ii in np.arange(num_times):
        print(ii)
        dt_begin_local1 = dt_begin_str + timedelta(days = int(ii))
        dt_end_local1   = dt_begin_local1 + timedelta(days = interval)

        ##!#combined_data[dt_begin_local1.strftime('%Y')] = {}
        ##!#combined_data_nosmoke[dt_begin_local1.strftime('%Y')] = {}

        print(dt_begin_local1, dt_end_local1)

        CERES_data1, NAAPS_data = read_zoom_CERES_data_single(\
            dt_begin_local1.strftime('%Y%m%d'), \
            dt_end_local1.strftime('%Y%m%d'), \
            ceres_var, min_ice, \
            min_smoke, max_smoke, vmin2, NAAPS_data, minlat, lat_bounds, \
            lon_bounds, satellite = 'All', mask_NAAPS = False)

        bdata1 = CERES_data1[ceres_var].flatten().compressed()
        bdata3 = CERES_data1[ceres_var+'_nosmoke'].flatten().compressed()
  
        combined_data[dt_begin_local1.strftime('%m/%d') + ' - ' + \
            dt_end_local1.strftime('%m/%d')] = bdata1

        combined_data_nosmoke[dt_begin_local1.strftime('%m/%d') + ' - ' + \
            dt_end_local1.strftime('%m/%d')] = bdata3

        # Add the average data
        # --------------------
        average_data[cnt]     = np.nanmean(bdata1)

        cnt += 1

    second_labels = np.array([bkey for bkey in combined_data.keys()])
    labels = second_labels.flatten()
    second_arrays = np.array([combined_data[bkey] for bkey in \
        combined_data.keys()], dtype = object)
    second_arrays_nosmoke = np.array([combined_data_nosmoke[bkey] for \
        bkey in combined_data_nosmoke.keys()], dtype = object)
    
    #print(combined_data['2012']['06/25 - 06/30'])

    return second_labels, second_arrays, average_data, second_arrays_nosmoke

def read_CERES_region_time_series_multiyear(dt_begin_str, \
        dt_end_str, begin_year, end_year, interval, ceres_var, min_ice, \
        min_smoke, max_smoke, vmin2, satellite, NAAPS_data, minlat, \
        lat_bounds, lon_bounds):

    # Find the number of intervals between the two dates
    # --------------------------------------------------
    num_times = (dt_end_str - dt_begin_str).days - (interval - 1)
    years     = np.arange(begin_year,end_year+1)
     
    total_smoke_arrays   = np.full((len(years), num_times), np.nan)
    total_nosmoke_arrays = np.full((len(years), num_times), np.nan)

    for ii, year in enumerate(years):
        dt_local_begin = dt_begin_str.replace(year = years[ii])
        dt_local_end   = dt_end_str.replace(year = years[ii])

        second_labels, second_arrays, average_CERES, second_arrays_nosmoke = \
            read_CERES_region_time_series(dt_local_begin, dt_local_end, \
            interval, ceres_var, min_ice, min_smoke, max_smoke, vmin2, \
            satellite, NAAPS_data, minlat, lat_bounds, lon_bounds, \
            plot_daily_data = False)

        smoke_means   = np.array([np.nanmean(tdata) for tdata in second_arrays])
        nosmoke_means = np.array([np.nanmean(tdata) for tdata in \
            second_arrays_nosmoke])

        total_smoke_arrays[ii,:] = smoke_means
        total_nosmoke_arrays[ii,:] = nosmoke_means
   

    return total_smoke_arrays, total_nosmoke_arrays, second_labels
    #return second_arrays, second_arrays_nosmoke

# dt_begin_str = beginning of time window for traveling averaging
# dt_end_str   = end of time window for traveling averaging
# interval     = length of time window for averaging (in days)
# ceres_var    = 'alb_clr', 'swf_clr', 'lwf_clr'
# min_ice      = minimum NSIDC ice threshold for analysis
# plats        = point lats for pulling data at a single grid point
# plons        = point lons for pulling data at a single grid point
def read_CERES_points_time_series(dt_begin_str, dt_end_str, interval, \
        plats, plons, ceres_var, min_ice, min_smoke, max_smoke, \
        vmin2, satellite, NAAPS_data, minlat, lat_bounds, lon_bounds, \
        plot_daily_data = False):

    # Find the number of intervals between the two dates
    # --------------------------------------------------
    num_times = (dt_end_str - dt_begin_str).days - (interval - 1)

    print('num_times = ',num_times)

    combined_data = {}
    cnt = 0
    
    combined_data['plats'] = plats
    combined_data['plons'] = plats

    # Grab single point NAAPS data here
    # ---------------------------------
    combined_data['smoke_conc_sfc'] = {}
    combined_data['ceres_data'] = {}

    for jj in range(len(plats)):
        c_idx = nearest_gridpoint(plats[jj], plons[jj],\
            NAAPS_data['lat'], NAAPS_data['lon'])
        if(len(c_idx[0]) > 1):
            c_idx = (np.array([c_idx[0][0]])), (np.array([c_idx[1][0]]))

        combined_data['smoke_conc_sfc'][str(jj)] = \
            NAAPS_data['smoke_conc_sfc'][c_idx]

    for ii in np.arange(num_times):
        print(ii)
        dt_begin_local1 = dt_begin_str + timedelta(days = int(ii))
        dt_end_local1   = dt_begin_local1 + timedelta(days = interval)

        print(dt_begin_local1, dt_end_local1)

        CERES_data1, NAAPS_data = read_zoom_CERES_data_single(\
            dt_begin_local1.strftime('%Y%m%d'), \
            dt_end_local1.strftime('%Y%m%d'), \
            ceres_var, min_ice, \
            min_smoke, max_smoke, vmin2, NAAPS_data, minlat, lat_bounds, \
            lon_bounds, satellite = 'All', mask_NAAPS = False)

        combined_data['ceres_data'][dt_begin_local1.strftime('%m/%d') + ' - ' + \
            dt_end_local1.strftime('%m/%d')] = {}

        # Grab single point data here
        # --------------------------------------------------
        for jj in range(len(plats)):
             
            c_idx = nearest_gridpoint(plats[jj], plons[jj],\
                CERES_data1['lat'], CERES_data1['lon'])
            if(len(c_idx[0]) > 1):
                c_idx = (np.array([c_idx[0][0]])), (np.array([c_idx[1][0]]))

            pval = CERES_data1[ceres_var][c_idx]

            combined_data['ceres_data'][\
                dt_begin_local1.strftime('%m/%d') + ' - ' + \
                dt_end_local1.strftime('%m/%d')][str(jj)] = pval

    return combined_data

# dt_begin_str = beginning of time window for traveling averaging
# dt_end_str   = end of time window for traveling averaging
# interval     = length of time window for averaging (in days)
# ceres_var    = 'alb_clr', 'swf_clr', 'lwf_clr'
# min_ice      = minimum NSIDC ice threshold for analysis
# plats        = point lats for pulling data at a single grid point
# plons        = point lons for pulling data at a single grid point
def read_CERES_all_region(dt_begin_str1, dt_end_str1, \
        dt_begin_str2, dt_end_str2, ceres_var, \
        min_ice, min_smoke, max_smoke, vmin2, satellite, \
        NAAPS_data, minlat, lat_bounds, lon_bounds, \
        plot_daily_data = False):

    # Read a temporary CERES dataset to figure out the data dimensions
    # ----------------------------------------------------------------
    CERES_data1, CERES_data2, NAAPS_data = read_zoom_CERES_data(\
        dt_begin_str1.strftime('%Y%m%d'), \
        dt_end_str1.strftime('%Y%m%d'), \
        dt_begin_str2.strftime('%Y%m%d'), \
        dt_end_str2.strftime('%Y%m%d'), \
        ceres_var, min_ice, min_smoke, max_smoke, vmin2, NAAPS_data, minlat, \
        lat_bounds, lon_bounds, \
        satellite = satellite, mask_NAAPS = False)

    # Mask both CERES data so that only gridpoints within the lat/lon
    # bounds that contain data for both time periods are reserved
    # ---------------------------------------------------------------
    mask_CERES_data1 = np.ma.masked_where(CERES_data2[ceres_var].mask, \
        CERES_data1[ceres_var])
    mask_CERES_data2 = np.ma.masked_where(CERES_data1[ceres_var].mask, \
        CERES_data2[ceres_var])

    mask_lats = CERES_data1['lat'][~mask_CERES_data1.mask]
    mask_lons = CERES_data1['lon'][~mask_CERES_data1.mask]
    mask_NAAPS = NAAPS_data['smoke_conc_sfc'][~mask_CERES_data1.mask]
    mask_CERES_data1 = mask_CERES_data1.compressed()
    mask_CERES_data2 = mask_CERES_data2.compressed()
   
    combined_data = {}
    combined_data['lats']          = mask_lats
    combined_data['lons']          = mask_lons
    combined_data['NAAPS_data']    = mask_NAAPS
    combined_data['CERES_data1']   = mask_CERES_data1
    combined_data['CERES_data2']   = mask_CERES_data2
    combined_data['dt_begin_str1'] = dt_begin_str1
    combined_data['dt_end_str1']   = dt_end_str1
    combined_data['dt_begin_str2'] = dt_begin_str2
    combined_data['dt_end_str2']   = dt_end_str2
 
    return combined_data



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Calculation functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

# Could modify to allow user to calculate either daily or monthly
# average.
def calc_NAAPS_month_avg(date_str, minlat = 65., mask_zero = False):

    # Make sure date string is of correct format
    # ------------------------------------------
    if(len(date_str) != 6):
        print("ERROR: invalid date_str value")
        print("       Must be of format YYYYMM")
        return -1
    
    dt_date_str = datetime.strptime(date_str, '%Y%m')

    # Determine if the associated directory is untarred
    # -------------------------------------------------
    if(not os.path.isdir(dt_date_str.strftime(data_dir + '%Y/%Y%m/'))):
        print("ERROR: data directory for ",date_str," is not untarred")
        print("       Must untar associated NAAPS gzip tarfile before")
        print("       continuing")
        return -2 

    # Find all relevant files in the associated directory
    # ---------------------------------------------------
    files = glob(dt_date_str.strftime(data_dir + '%Y/%Y%m/NVA_*%Y%m*.nc'))
    num_files = len(files)

    # Allocate array to hold all the data north of minlat before averaging
    # --------------------------------------------------------------------
    num_lats = int(90 - minlat)
    test_lats = np.arange(-89.5, 90.5, 1)
    all_data = np.full((num_files,  num_lats, 360), np.nan)

    # Insert all the data for this month into the array
    # -------------------------------------------------
    beg_idx = np.where(test_lats > minlat)[0][0]
    print(beg_idx)
    for ii in range(num_files):
        print(files[ii])
        data = Dataset(files[ii], 'r')
 
        all_data[ii,:,:] = data['smoke_conc_sfc'][beg_idx:,:]

        data.close()

    # Average the data down into a single monthly average
    # ---------------------------------------------------
    if(mask_zero):
        all_data = np.ma.masked_where(all_data == 0., all_data)

    all_data = np.nanmean(all_data, axis = 0)

    # Return the averaged data
    return all_data

def calc_NAAPS_all_avgs(begin_date, end_date, minlat = 65., \
        mask_zero = False):

    begin_date_str = datetime.strptime(begin_date, '%Y%m')
    end_date_str   = datetime.strptime(end_date, '%Y%m')

    # Find the number of months to analyze
    # ------------------------------------
    local_date_str = begin_date_str
    num_months = 0
    date_strs = []
    while(local_date_str <= end_date_str):
        if((local_date_str.month >= 4 ) & (local_date_str.month <= 9)):
            num_months += 1
            date_strs.append(local_date_str.strftime('%Y%m'))
    
        local_date_str = local_date_str + relativedelta(months = 1)

    # Prep the data file
    num_lats = int(90 - minlat)
    local_lats = np.arange(minlat + 0.5, 90.5, 1)
    local_lons = np.arange(-179.5, 180.5, 1)
    all_data = np.full((num_months, num_lats, 360), np.nan)

    for ii, dstr in enumerate(date_strs):
        print(dstr)
        all_data[ii,:,:] = calc_NAAPS_month_avg(dstr, minlat = minlat, \
            mask_zero = mask_zero)
    
    grid_lons, grid_lats = np.meshgrid(local_lons, local_lats)
     
    NAAPS_avgs = {}
    NAAPS_avgs['begin_str'] = begin_date
    NAAPS_avgs['end_str']   = end_date
    NAAPS_avgs['dates']     = date_strs
    NAAPS_avgs['data']      = all_data   
    NAAPS_avgs['lats']      = grid_lats
    NAAPS_avgs['lons']      = grid_lons
 
    return NAAPS_avgs

def calcNAAPS_grid_trend(NAAPS_data, month_idx, trend_type, minlat):

    if(month_idx == None):
        month_idx = 0
        index_jumper = 1
    else:
        month_adder = '_month'
        index_jumper = 6 

    lat_ranges = np.arange(minlat+0.5,90.5,1.0)
    lon_ranges = np.arange(-179.5,180.5,1.0)

    # Make copy of NAAPS_data array
    print('HERE:',NAAPS_data['dates'][month_idx::index_jumper])
    local_data   = np.copy(NAAPS_data['data'][month_idx::index_jumper,:,:])
    local_mask = np.ma.masked_where((local_data == -999.9) & \
        (NAAPS_data['lats'] < minlat), local_data)

    smoke_trends = np.full(local_data.shape[1:], np.nan)
    smoke_pvals  = np.full(local_data.shape[1:], np.nan)
    smoke_uncert = np.full(local_data.shape[1:], np.nan)

    print(local_data.shape)

    # Loop over all the keys and print the regression slopes 
    # Grab the averages for the key
    for i in range(0,len(np.arange(np.min(NAAPS_data['lats']),90))):
        for j in range(0,len(lon_ranges)):
            # Check the current max and min
            x_vals = np.arange(0,len(local_mask[:,i,j]))
            # Find the slope of the line of best fit for the time series of
            # average data
            if(trend_type=='standard'): 
                #slope, intercept, r_value, p_value, std_err, intcpt_stderr = \
                result = stats.linregress(x_vals,local_mask[:,i,j])
                smoke_trends[i,j] = result.slope * len(x_vals)
                smoke_pvals[i,j]  = result.pvalue
                smoke_uncert[i,j] = result.stderr * len(x_vals)
            else:
                res = stats.theilslopes(local_mask[:,i,j], x_vals, 0.90)
                smoke_trends[i,j] = res[0]*len(x_vals)

    smoke_trends = np.ma.masked_where(NAAPS_data['lats'] < minlat, smoke_trends)
    smoke_pvals  = np.ma.masked_where(NAAPS_data['lats'] < minlat, smoke_pvals)
    smoke_uncert = np.ma.masked_where(NAAPS_data['lats'] < minlat, smoke_uncert)

    print('in trend calc')
    for x, y in zip(NAAPS_data['lats'][:,10], smoke_trends[:,10]):
        print(x,y)

    return smoke_trends, smoke_pvals, smoke_uncert

# This function assumes the data is being read from the netCDF file
# NOTE: Assume user is using new NAAPS climo file which starts in January
def calcNAAPS_MonthClimo(NAAPS_data):

    # Set up arrays to hold monthly climatologies
    month_climo = np.zeros((6, NAAPS_data['data'].shape[1], \
        NAAPS_data['data'].shape[2]))

    # Mask the monthly averages
    local_data   = np.copy(NAAPS_data['data'][:,:,:])
    local_mask = np.ma.masked_where(local_data == -999.9, local_data)
 
    # Calculate monthly climatologies
    for m_i in range(6):
        month_climo[m_i,:,:] = np.nanmean(local_data[m_i::6,:,:],axis=0)
        print("Month: ",NAAPS_data['dates'][m_i][4:]) 
        ##!#month_climo[m_i,:,:] = np.nanmean(local_mask[m_i::12,:,:],axis=0)
        ##!#print("Month: ",NAAPS_data['DATES'][m_i][4:], NAAPS_data['DATES'][::12]) 
    #month_climo[0,:,:]  = np.nanmean(local_mask[0::12,:,:],axis=0)  # January 
    #month_climo[1,:,:]  = np.nanmean(local_mask[1::12,:,:],axis=0)  # February
    #month_climo[2,:,:]  = np.nanmean(local_mask[2::12,:,:],axis=0)  # March 
    #month_climo[3,:,:]  = np.nanmean(local_mask[3::12,:,:],axis=0)  # April 
    #month_climo[4,:,:]  = np.nanmean(local_mask[4::12,:,:],axis=0)  # May 
    #month_climo[5,:,:]  = np.nanmean(local_mask[5::12,:,:],axis=0)  # June 
    #month_climo[6,:,:]  = np.nanmean(local_mask[6::12,:,:],axis=0)  # July 
    #month_climo[7,:,:]  = np.nanmean(local_mask[7::12,:,:],axis=0)  # August 
    #month_climo[8,:,:]  = np.nanmean(local_mask[8::12,:,:],axis=0)  # September 
    #month_climo[9,:,:]  = np.nanmean(local_mask[9::12,:,:],axis=0)  # October 
    #month_climo[10,:,:] = np.nanmean(local_mask[10::12,:,:],axis=0) # November
    #month_climo[11,:,:] = np.nanmean(local_mask[11::12,:,:],axis=0) # December

    # Insert data into dictionary
    NAAPS_data['MONTH_CLIMO'] = month_climo

    return NAAPS_data


def writeNAAPS_toNCDF(NAAPS_data,file_name, minlat = 65.):

    lat_ranges = np.arange(minlat + 0.5,90.5,1.0)
    lon_ranges = np.arange(-179.5,180.5,1.0)
   
    # Create a new netCDF dataset to write to
    # --------------------------------------- 
    nc = Dataset(file_name,'w',format='NETCDF4')
  
    # Dimensions = lat, lon, time
    # Create the sizes of each dimension in the file. In this case,
    # the dimensions are "# of latitude", "# of longitude", and 
    # "# of times in the file". 
    # -------------------------------------------------------------
    num_lat  = len(lat_ranges)
    num_lon  = len(lon_ranges)
    num_time = len(NAAPS_data['dates'])
    times = np.arange(num_time)
  
    # Instead of using a simple 'arange' function to define the 'times'
    # variable, actually calculate the number of months between each date
    # variable and a reference date, which is January 2005. This will have
    # no effect on any total NAAPS processes but will allow compatibility with
    # testing processes where just a few months between January 2005 and
    # July 2019 are selected. 
    base_date = datetime(year=2005,month=1,day=1) 
    times = np.array([(datetime.strptime(tmpx,'%Y%m').year - base_date.year) \
        * 12 + datetime.strptime(tmpx,'%Y%m').month - base_date.month \
        for tmpx in NAAPS_data['dates']])
   
    # Use the dimension size variables to actually create dimensions in 
    # the file.
    # ----------------------------------------------------------------- 
    n_time = nc.createDimension('nmth',num_time)
    n_lat  = nc.createDimension('dlat',num_lat)
    n_lon  = nc.createDimension('dlon',num_lon)

    # Create variables for the three dimensions. Note that since these
    # are variables, they are still given 'dimensions' using 'nmth',
    # 'dlat', and 'dlon'. Latitude and longitude are each 2-d grids, 
    # so they are given 2 dimensions (dlat, dlon).
    # ----------------------------------------------------------------
    MONTH = nc.createVariable('MONTH','i2',('nmth'))
    MONTH.description = 'Months since January 2005'
    LAT = nc.createVariable('Latitude','f4',('dlat','dlon'))
    LAT.description = 'Latitude (degrees North)'
    LAT.units = 'Degrees'
    LON = nc.createVariable('Longitude','f4',('dlat','dlon'))
    LON.description = 'Longitude (-180 - 180)'
    LON.units = 'Degrees'
   
    # Create a variable for the AI data, dimensioned using all three
    # dimensions.
    # --------------------------------------------------------------  
    CONC = nc.createVariable('smoke_conc_sfc','f4',('nmth','dlat','dlon'))
    CONC.description = 'Monthly averaged NAAPS smoke aerosol concentration '+\
        'in the lowest NAAPS model layer'

    # Fill in dimension variables one-by-one.
    # NOTE: not sure if you can insert the entire data array into the
    # dimension (for example, doing MONTH = times), so this could be
    # something for you to try. Might make this faster if it works
    # ---------------------------------------------------------------
    MONTH[:] = times[:]
    LAT[:,:] = NAAPS_data['lats']
    LON[:,:] = NAAPS_data['lons']
    ##for i in range(num_time):
    #    MONTH[i] = times[i]
    #for i in range(num_lat):
    #    for j in range(num_lon):
    #        LAT[i,j]=lat_ranges[i]
    #        LON[i,j]=lon_ranges[j]

    # Fill in actual variables
    # I have a bunch of extra stuff in here for handling the dictionary
    # keys, which you will likely not need for your data
    # ------------------------------------------------------------------
    CONC[:,:,:] = NAAPS_data['data'][:,:,:]
    #for ii in range(num_time):
    ##!#for i in range(num_lat):
    ##!#    print(lat_ranges[i])
    ##!#    for j in range(num_lon):
    ##!#        dictkey = (str(int(lat_ranges[i]))+'x'+str(int(lon_ranges[j])))
    ##!#        if(dictkey not in NAAPS_data):
    ##!#            # Insert missing values for AI and count
    ##!#            AI[:,i,j] = [-999.9 for m in range(num_time)]
    ##!#            OB_COUNT[:,i,j] = [-99 for m in range(num_time)]
    ##!#        else:
    ##!#            for m in range(num_time):
    ##!#                timekey = testkeys[m]
    ##!#                if(timekey not in NAAPS_data[dictkey]):
    ##!#                    AI[m,i,j] = -999.9
    ##!#                    OB_COUNT[m,i,j] = -99
    ##!#                else:
    ##!#                    AI[m,i,j] = NAAPS_data[dictkey][timekey]['avg']
    ##!#                    OB_COUNT[m,i,j] = NAAPS_data[dictkey][timekey]['#_obs']

    # Close the netCDF file, which actually saves it.
    # -----------------------------------------------
    nc.close()
    print("Saved file", file_name)

def readgridNAAPS_NCDF(infile=home_dir + \
        '/Research/NAAPS/naaps_grid_smoke_conc_sfc_2005_2020_noAI.nc',\
        start_date = 200504, end_date = 202009, calc_month = True, \
        minlat = 65.5):

    # Read in data to netCDF object
    in_data = Dataset(infile,'r')

    # Set up dictionary to hold data
    NAAPS_data = {}

    # Set up date strings in the file
    base_date = datetime(year = 2005, month = 1, day = 1)
    str_dates = \
        np.array([(base_date + relativedelta(months=mi)).strftime('%Y%m') for mi in \
            in_data['MONTH']])

    # Use minlat to restrict the data to only the desired minimum
    # latitude
    # ------------------------------------------------------------
    min_idx =  np.where(in_data['Latitude'][:,0] > minlat)[0][0]

    # Use the file dates to select just the data within the user-specified
    # timeframe
    # --------------------------------------------------------------------
    int_str_dates = np.array([int(tdate) for tdate in str_dates])

    time_indices = np.where((int_str_dates >= start_date) & \
        (int_str_dates <= end_date))[0]

    NAAPS_data['begin_str'] = str_dates[0]
    NAAPS_data['end_str']   = str_dates[-1]
    NAAPS_data['dates']     = str_dates
    NAAPS_data['data']      = in_data['smoke_conc_sfc'][:,min_idx:,:]
    NAAPS_data['lats']      = in_data['Latitude'][min_idx:,:]
    NAAPS_data['lons']      = in_data['Longitude'][min_idx:,:]
 

    if(calc_month == True):
        NAAPS_data = calcNAAPS_MonthClimo(NAAPS_data)

    # to add months to datetime object, do
    ###from dateutil.relativedelta import relativedelta
    ###datetime.datetime(year=2004,month=10,day=1) + relativedelta(months=1)
    
    in_data.close()
   
    return NAAPS_data


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#
# Plottng functions
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

# ptype is either 'trend', 'climo', or 'monthclimo'
# title is the plot title
def plotNAAPS_spatial(pax, plat, plon, pdata, ptype, ptitle = '', plabel = '', \
        vmin = None, vmax = None, colorbar_label_size = 14, colorbar = True, \
        minlat = 65., pvals = None):

    if(vmin == None):
        vmin = np.nanmin(pdata)
    if(vmax == None):
        vmax = np.nanmax(pdata)

    if(ptype == 'trend'):
        colormap = plt.cm.bwr
        #colormap = plt.cm.get_cmap('bwr', 5)
    elif(ptype == 'uncert'):
        #colormap = plt.cm.plasma
        colormap = plt.cm.get_cmap('jet', 6)
    else:
        colormap = plt.cm.jet

    # Make copy of NAAPS_data array
    local_data  = np.copy(pdata)
    # Grid the data, fill in white space
    cyclic_data,cyclic_lons = add_cyclic_point(local_data,plon[0,:])
    plat2,plon2 = np.meshgrid(plat[:,0],cyclic_lons)   
  
    # Mask any missing values
    mask_smoke = np.ma.masked_where(cyclic_data < -998.9, cyclic_data)
    mask_smoke = np.ma.masked_where(plat2.T < minlat, mask_smoke)

    # Plot lat/lon lines
    gl = pax.gridlines(crs = ccrs.PlateCarree(), draw_labels = False, \
        linewidth = 1, color = 'gray', alpha = 0.5, linestyle = '-',\
        y_inline = True, xlocs = range(-180, 180, 30), ylocs = range(70, 90, 10))
    #gl.top_labels = False
    #gl.bottom_labels = False
    #gl.left_labels = False
    #gl.right_labels = False
    #gl.xlines = False
    #gl.xlocator = mticker.FixedLocator(np.arange(-180, 180, 30))
    #gl.ylocator = mticker.FixedLocator(np.arange(70, 90, 10))
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    #gl.xlabel_style = {'size': 5, 'color': 'gray'}
    gl.xlabel_style = {'color': 'gray', 'weight': 'bold'}
    gl.ylabel_style = {'size': 15, 'color': 'gray'}
    gl.ylabel_style = {'color': 'black', 'weight': 'bold'}
    #pax.gridlines()

    pax.coastlines(resolution='50m')
    mesh = pax.pcolormesh(plon2, plat2,\
            mask_smoke.T,transform = datacrs,\
            cmap = colormap,vmin=vmin,vmax=vmax, shading = 'auto')
    if(pvals is not None):
        cyclic_pvals, cyclic_lons = add_cyclic_point(pvals, plon[0,:])
        print(pvals.shape, plat.shape, plat2.shape)
        mask_pvals = np.ma.masked_where((plat < minlat) | \
            (pvals > 0.05), pvals)
        pax.pcolor(plon, plat, mask_pvals, hatch = '...', alpha = 0.0, \
            shading = 'auto', transform = datacrs)

    pax.set_extent([-180,180,minlat,90],datacrs)
    pax.set_boundary(circle, transform=pax.transAxes)
    #ax.set_xlim(-3430748.535086173,3430748.438879491)
    #ax.set_ylim(-3413488.8763307533,3443353.899053069)
    #cbar = plt.colorbar(mesh,ticks = np.arange(-2.0,4.1,0.2),\
    if(colorbar):
        cbar = plt.colorbar(mesh,\
            ax = pax, orientation='vertical',shrink = 0.8, extend = 'both')
        #cbar10.set_label('UV Aerosol Index',weight='bold',fontsize=colorbar_label_size)
        #cbar.ax.tick_params(labelsize=14)
        cbar.set_label(plabel,fontsize=colorbar_label_size,weight='bold')
    pax.set_title(ptitle)

# Designed to work with the netCDF data
def plotNAAPS_MonthTrend(NAAPS_data,month_idx=None,save=False,\
        trend_type='standard',minlat=65.,return_trend=False, colorbar = True, \
        title = '', label = '', colorbar_label_size = 14, pax = None, \
        show_pval = False, uncert_ax = None, vmax = None):

    trend_label=''
    if(trend_type=='theil-sen'):
        trend_label='_theilSen'

    fillval = 1.0
    max_vals = [fillval, fillval, fillval, fillval, fillval, fillval]
    #max_vals = [0.5, 0.7, 3.0, 7.0, 12.0, 10.0]

    if(month_idx == None):
        month_adder = ''
        do_month = False
        v_max = 0.5
        v_min = -0.5
    else:
        month_adder = '_month'
        do_month = True
        index_jumper = 6 
        if(vmax is None):
            vmax = max_vals[month_idx]
            vmin = -max_vals[month_idx]
        else:
            vmin = -vmax
    # --------------------------------------------------------------
    #
    # Use calcNAAPS_grid_trend to calculate the trends in the AI data
    #
    # --------------------------------------------------------------
    smoke_trends, smoke_pvals, smoke_uncert = \
        calcNAAPS_grid_trend(NAAPS_data, month_idx, trend_type, \
    #ai_trends, ai_pvals = calcNAAPS_grid_trend(NAAPS_data, month_idx, trend_type, \
        minlat)

    if(not show_pval):
        smoke_pvals = None
    else:
        print('month_idx = ',month_idx,' PVAL nanmean = ', \
            np.nanmean(smoke_pvals))

    if(uncert_ax is None):
        smoke_uncert = None
    # --------------------------------------------------------------
    #
    # Plot the calculated trends on a figure
    #
    # --------------------------------------------------------------
    colormap = plt.cm.bwr

    # Pull the beginning and ending dates into datetime objects
    start_date = datetime.strptime(NAAPS_data['dates'][\
        month_idx::index_jumper][0],'%Y%m')
    end_date   = datetime.strptime(NAAPS_data['dates'][\
        month_idx::index_jumper][-1],'%Y%m')
    
    # Make figure title
    month_string = ''
    if(do_month == True):
        month_string = start_date.strftime('%B') + ' '

    if(title == ''):
        title = 'NAAPS Smoke Conc Sfc ' + month_string + 'Trends\n'+\
            start_date.strftime('%b. %Y') + ' - ' + \
            end_date.strftime('%b. %Y')
    if(label == ''):
        label = 'NAAPS Smoke Conc Sfc Trend'

    # Call plotNAAPS_spatial to add the data to the figure
    # Make figure
    if(pax is None):
        # Make figure
        plt.close('all')
        fig1 = plt.figure(figsize = (8,8))
        ax = plt.axes(projection = mapcrs)

        plotNAAPS_spatial(ax, NAAPS_data['lats'], NAAPS_data['lons'], \
            smoke_trends, 'trend', ptitle = title, plabel = label, \
            colorbar = colorbar, colorbar_label_size = colorbar_label_size, \
            vmin = vmin, vmax = vmax, minlat = minlat, pvals = smoke_pvals)

        fig1.tight_layout()

        if(save == True):
            month_adder = ''
            if(do_month == True):
                month_adder = '_' + start_date.strftime('%b') 
            out_name = 'naaps_grid_trend' + month_adder + '_' + \
                start_date.strftime('%Y%m') + '_' + end_date.strftime('%Y%m') + \
                '.png'
            plt.savefig(out_name,dpi=300)
            print("Saved image",out_name)
        else:
            plt.show()
    else:
        plotNAAPS_spatial(pax, NAAPS_data['lats'], NAAPS_data['lons'], \
            smoke_trends, 'trend', ptitle = title, plabel = label, \
            colorbar = colorbar, colorbar_label_size = colorbar_label_size, \
            vmin = vmin, vmax = vmax, minlat = minlat, pvals = smoke_pvals)

    if(uncert_ax is not None):
        plotNAAPS_spatial(uncert_ax, NAAPS_data['lats'], NAAPS_data['lons'], \
            smoke_uncert, 'uncert', ptitle = title, plabel = label, \
            colorbar = colorbar, colorbar_label_size = colorbar_label_size, \
            vmin = 0, vmax = 4.0, minlat = minlat)

    if(return_trend == True):
        return smoke_trends

# Generate a 15-panel figure comparing the climatology and trend between 3
# versions of the NAAPS data for all months
def plotNAAPS_ClimoTrend_all(NAAPS_data,\
        trend_type = 'standard', minlat=65,save=False):

    colormap = plt.cm.jet

    lat_ranges = np.arange(minlat,90,1.0)
    lon_ranges = np.arange(-180,180,1.0)

    index_jumper = 6 

    colorbar_label_size = 7
    axis_title_size = 8
    row_label_size = 10 
    #colorbar_label_size = 13
    #axis_title_size = 14.5
    #row_label_size = 14.5

    #fig = plt.figure()
    plt.close('all')
    fig = plt.figure(figsize=(9,13))
    #plt.suptitle('NAAPS Comparisons: '+start_date.strftime("%B"),y=0.95,\
    #    fontsize=18,fontweight=4,weight='bold')
    gs = gridspec.GridSpec(nrows=6, ncols=3, hspace = 0.001, wspace = 0.15)

    # - - - - - - - - - - - - - - - - - - - - -
    # Plot the climatologies along the top row
    # - - - - - - - - - - - - - - - - - - - - -
       
    # Plot DATA1 climos
    # -----------------
    ##!## Make copy of NAAPS_data array
    local_data1_Apr  = np.copy(NAAPS_data['MONTH_CLIMO'][0,:,:])
    local_data1_May  = np.copy(NAAPS_data['MONTH_CLIMO'][1,:,:])
    local_data1_Jun  = np.copy(NAAPS_data['MONTH_CLIMO'][2,:,:])
    local_data1_Jul  = np.copy(NAAPS_data['MONTH_CLIMO'][3,:,:])
    local_data1_Aug  = np.copy(NAAPS_data['MONTH_CLIMO'][4,:,:])
    local_data1_Sep  = np.copy(NAAPS_data['MONTH_CLIMO'][5,:,:])

    mask_AI1_Apr = np.ma.masked_where(local_data1_Apr == -999.9, local_data1_Apr)
    mask_AI1_May = np.ma.masked_where(local_data1_May == -999.9, local_data1_May)
    mask_AI1_Jun = np.ma.masked_where(local_data1_Jun == -999.9, local_data1_Jun)
    mask_AI1_Jul = np.ma.masked_where(local_data1_Jul == -999.9, local_data1_Jul)
    mask_AI1_Aug = np.ma.masked_where(local_data1_Aug == -999.9, local_data1_Aug)
    mask_AI1_Sep = np.ma.masked_where(local_data1_Sep == -999.9, local_data1_Sep)
    ##!## Grid the data, fill in white space
    ##!#cyclic_data,cyclic_lons = add_cyclic_point(local_data,NAAPS_data1['LON'][0,:])
    ##!#plat,plon = np.meshgrid(NAAPS_data1['LAT'][:,0],cyclic_lons)   
  
    # Mask any missing values
    #mask_AI = np.ma.masked_where(plat.T < minlat, mask_AI)
    ax00 = plt.subplot(gs[0,0], projection=mapcrs)   # April climo original
    ax01 = plt.subplot(gs[0,1], projection=mapcrs)   # April climo screened
    ax02 = plt.subplot(gs[0,2], projection=mapcrs)   # April trend original
    ax10 = plt.subplot(gs[1,0], projection=mapcrs)   # May climo original
    ax11 = plt.subplot(gs[1,1], projection=mapcrs)   # May climo screened
    ax12 = plt.subplot(gs[1,2], projection=mapcrs)   # May trend original
    ax20 = plt.subplot(gs[2,0], projection=mapcrs)   # June climo original
    ax21 = plt.subplot(gs[2,1], projection=mapcrs)   # June climo screened
    ax22 = plt.subplot(gs[2,2], projection=mapcrs)   # June trend original
    ax30 = plt.subplot(gs[3,0], projection=mapcrs)   # July climo original
    ax31 = plt.subplot(gs[3,1], projection=mapcrs)   # July climo screened
    ax32 = plt.subplot(gs[3,2], projection=mapcrs)   # July trend original
    ax40 = plt.subplot(gs[4,0], projection=mapcrs)   # August climo original
    ax41 = plt.subplot(gs[4,1], projection=mapcrs)   # August climo screened
    ax42 = plt.subplot(gs[4,2], projection=mapcrs)   # August trend original
    ax50 = plt.subplot(gs[5,0], projection=mapcrs)   # September climo original
    ax51 = plt.subplot(gs[5,1], projection=mapcrs)   # September climo screened
    ax52 = plt.subplot(gs[5,2], projection=mapcrs)   # September trend original

    # Plot the figures in the first row: April
    # ---------------------------------------
    cbar_switch = True
    plotNAAPS_spatial(ax00, NAAPS_data['lats'], NAAPS_data['lons'], mask_AI1_Apr, \
        'climo', ptitle = ' ', plabel = '', vmin = 0, vmax = 1.0, \
        minlat = minlat, colorbar = cbar_switch)
    plotNAAPS_MonthTrend(NAAPS_data,month_idx=0,trend_type=trend_type,label = ' ',\
        minlat=minlat,title = ' ', pax = ax01, colorbar = cbar_switch, \
        colorbar_label_size = colorbar_label_size, show_pval = True, \
        uncert_ax = ax02)

    # Plot the figures in the first row: May
    # ---------------------------------------
    plotNAAPS_spatial(ax10, NAAPS_data['lats'], NAAPS_data['lons'], mask_AI1_May, \
        'climo', ptitle = ' ', plabel = '', vmin = 0, vmax = 2.0, \
        minlat = minlat, colorbar = cbar_switch)
    plotNAAPS_MonthTrend(NAAPS_data,month_idx=1,trend_type=trend_type,label = ' ',\
        minlat=minlat,title = ' ', pax = ax11, colorbar = cbar_switch, \
        colorbar_label_size = colorbar_label_size, show_pval = True, \
        uncert_ax = ax12)

    # Plot the figures in the first row: June
    # ---------------------------------------
    plotNAAPS_spatial(ax20, NAAPS_data['lats'], NAAPS_data['lons'], mask_AI1_Jun, \
        'climo', ptitle = ' ', plabel = '', vmin = 0, vmax = 3.0, \
        minlat = minlat, colorbar = cbar_switch)
    plotNAAPS_MonthTrend(NAAPS_data,month_idx=2,trend_type=trend_type,label = ' ',\
        minlat=minlat,title = ' ', pax = ax21, colorbar = cbar_switch, \
        colorbar_label_size = colorbar_label_size, show_pval = True, \
        uncert_ax = ax22)

    # Plot the figures in the second row: July
    # ----------------------------------------
    plotNAAPS_spatial(ax30, NAAPS_data['lats'], NAAPS_data['lons'], mask_AI1_Jul, \
        'climo', ptitle = ' ', plabel = '', vmin = 0, vmax = 3.0, \
        minlat = minlat, colorbar = cbar_switch)
    plotNAAPS_MonthTrend(NAAPS_data,month_idx=3,trend_type=trend_type,label = ' ',\
        minlat=minlat,title = ' ', pax = ax31, colorbar = cbar_switch, \
        colorbar_label_size = colorbar_label_size, show_pval = True, \
        uncert_ax = ax32)

    # Plot the figures in the third row: August
    # -----------------------------------------
    plotNAAPS_spatial(ax40, NAAPS_data['lats'], NAAPS_data['lons'], mask_AI1_Aug, \
        'climo', ptitle = ' ', plabel = '', vmin = 0, vmax = 3.0, \
        minlat = minlat, colorbar = cbar_switch)
    plotNAAPS_MonthTrend(NAAPS_data,month_idx=4,trend_type=trend_type,label = ' ',\
        minlat=minlat,title = ' ', pax = ax41, colorbar = cbar_switch, \
        colorbar_label_size = colorbar_label_size, show_pval = True, \
        uncert_ax = ax42)

    # Plot the figures in the third row: September
    # --------------------------------------------
    plotNAAPS_spatial(ax50, NAAPS_data['lats'], NAAPS_data['lons'], mask_AI1_Sep, \
        'climo', ptitle = ' ', plabel = '', vmin = 0, vmax = 1.0, \
        minlat = minlat, colorbar = cbar_switch)
    plotNAAPS_MonthTrend(NAAPS_data,month_idx=5,trend_type=trend_type,label = ' ',\
        minlat=minlat,title = ' ', pax = ax51, colorbar = cbar_switch, \
        colorbar_label_size = colorbar_label_size, show_pval = True, \
        uncert_ax = ax52)

    fig.text(0.10, 0.82, 'April', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.692, 'May', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.565, 'June', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.435, 'July', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.305, 'August', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)
    fig.text(0.10, 0.18, 'September', ha='center', va='center', \
        rotation='vertical',weight='bold',fontsize=row_label_size + 1)

    fig.text(0.250, 0.90, 'Climatology', ha='center', va='center', \
        rotation='horizontal',weight='bold',fontsize=row_label_size)
    fig.text(0.500, 0.90, 'Trend', ha='center', va='center', \
        rotation='horizontal',weight='bold',fontsize=row_label_size)
    fig.text(0.75, 0.90, 'Standard Error of\nScreened Trend (Slope)', \
        ha='center', va='center', rotation='horizontal',weight='bold',\
        fontsize=row_label_size)


    #cax = fig.add_axes([0.15, 0.09, 0.35, 0.01])
    #norm = mpl.colors.Normalize(vmin = -0.5, vmax = 0.5)
    #cb1 = mpl.colorbar.ColorbarBase(cax, cmap = plt.cm.bwr, norm = norm, \
    #    orientation = 'horizontal', extend = 'both')
    #cb1.set_label('Sfc Smoke Trend (Smoke / Study Period)', \
    #    weight = 'bold')

    #cax2 = fig.add_axes([0.530, 0.09, 0.35, 0.01])
    #norm2 = mpl.colors.Normalize(vmin = 0.0, vmax = 0.3)
    #cmap = plt.cm.get_cmap('jet', 6)
    #cb2 = mpl.colorbar.ColorbarBase(cax2, cmap = cmap, norm = norm2, \
    #    orientation = 'horizontal', extend = 'both')
    #cb2.set_label('Standard Error of Smoke Trend (Slope)', weight = 'bold')

    #fig.tight_layout()

    outname = 'naaps_grid_comps_all6_2.png'
    if(save == True):
        plt.savefig(outname,dpi=300)
        print("Saved image",outname)
    else:
        plt.show()



def plot_NAAPS(NAAPS_data, var, ax = None, labelsize = 12, \
        plot_log = True, labelticksize = 10, zoom = True, vmin = None, \
        min_val = 0.0, \
        minlat = 65., circle_bound = True, vmax = None, save = False):

    in_ax = True
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, projection = mapcrs)

    if(plot_log):
        plot_data = np.ma.masked_where(NAAPS_data[var] < min_val, NAAPS_data[var])
        plot_data = np.log10(plot_data)
    else:
        plot_data = np.ma.masked_where(NAAPS_data[var] < min_val, NAAPS_data[var])
        plot_data = plot_data

    mesh = ax.pcolormesh(NAAPS_data['lon'], NAAPS_data['lat'], \
        plot_data, \
        cmap = 'viridis', vmin = vmin, vmax = vmax, \
        transform = ccrs.PlateCarree(), shading = 'auto')
    cbar = plt.colorbar(mesh, ax = ax, orientation='vertical',\
        pad=0.04, fraction = 0.040, extend = 'both')
    cbar.set_label(var, size = labelsize, weight = 'bold')
    #cbar.ax.tick_params(labelsize = labelticksize)
    ax.set_title(NAAPS_data['date'])

    if(circle_bound):
        ax.set_boundary(circle, transform = ax.transAxes)

    if(zoom):
        ax.set_extent([-180.0,180.0,minlat,90.0],\
                       ccrs.PlateCarree())
    else:
        ax.set_extent([-180.0,180.0,65.0,90.0],\
                       ccrs.PlateCarree())
    if(not in_ax):
        plt.show()

def plot_NAAPS_figure(date_str, var, minlat = 65., vmin = None, vmax = None, \
        circle_bound = True, plot_log = True, ptitle = '', zoom = True, \
        save_dir = home_dir + '/Research/NAAPS/single_time_images/', \
        save = False):

    plt.close('all')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection = mapcrs)

    # Read the data for this granule
    # ------------------------------
    NAAPS_data = read_NAAPS(date_str, minlat = minlat)

    # Plot the data for this granule
    # ------------------------------
    plot_NAAPS(NAAPS_data, var, ax = ax, zoom = zoom, \
            vmin = vmin, vmax = vmax, plot_log = plot_log, \
            circle_bound = circle_bound, minlat = minlat)

    ax.set_title('NAAPS-RA ' + var + '\n' + NAAPS_data['date'])
 
    ax.coastlines()

    if(save):
        outname = save_dir + 'naaps_' + var + '_' + date_str + '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

def plot_NAAPS_event(date_str, var, minlat = 65., ax = None, vmin = None, \
        vmax = None, plot_log = True, ptitle = '', circle_bound = True, \
        dtype = 'no_AI', zoom = True, save = False):

    plt.close('all')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection = mapcrs)

    # Read the data for this granule
    # ------------------------------
    NAAPS_data = read_NAAPS_event(date_str, minlat = minlat)

    # Plot the data for this granule
    # ------------------------------
    plot_NAAPS(NAAPS_data, var, ax = ax, zoom = zoom, \
            vmin = vmin, vmax = vmax, plot_log = plot_log, \
            circle_bound = circle_bound, minlat = minlat)

    ax.set_title('NAAPS-RA ' + var + '\n' + \
        NAAPS_data['dt_begin_date'].strftime('%Y-%m-%d') + ' - ' + \
        NAAPS_data['dt_end_date'].strftime('%Y-%m-%d'))
 
    ax.coastlines()

    if(save):
        outname = 'naaps_event_' + var + '_' + date_str + '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

# NOTE: Works for either total events or single times.
#       For a whole event, use YYYYMMDD
#       For a single time, use YYYYMMDDHH
def plot_NAAPS_compare_types(date_str, minlat = 65., ax = None, \
        vmin = None, vmax = None, plot_log = True, ptitle = '', \
        circle_bound = True, dtype = 'no_AI', zoom = True, save = False):
   
    plt.close('all')
    fig = plt.figure(figsize = (9,8.5))
    ax1 = fig.add_subplot(3,3,1, projection = mapcrs)
    ax2 = fig.add_subplot(3,3,2, projection = mapcrs)
    ax3 = fig.add_subplot(3,3,3, projection = mapcrs)
    ax4 = fig.add_subplot(3,3,4, projection = mapcrs)
    ax5 = fig.add_subplot(3,3,5, projection = mapcrs)
    ax6 = fig.add_subplot(3,3,6, projection = mapcrs)
    ax7 = fig.add_subplot(3,3,7, projection = mapcrs)
    ax8 = fig.add_subplot(3,3,8, projection = mapcrs)
    ax9 = fig.add_subplot(3,3,9, projection = mapcrs)

    # Read the data for this granule
    # ------------------------------
    if(len(date_str) == 10):
        max_dict = {
            'smoke_conc_sfc': 8,
            'smoke_drysink': 0.1,
            'smoke_wetsink': 3,
        }

        l_event = False
        NAAPS_noAI    = read_NAAPS(date_str, minlat = minlat, \
            dtype = 'no_AI')
        NAAPS_withAI  = read_NAAPS(date_str, minlat = minlat, \
            dtype = 'with_AI')
        dt_date_str = datetime.strptime(date_str, '%Y%m%d%H')
        plot_title = 'NAAPS-RA Comparisons\n' + \
            dt_date_str.strftime('%Y-%m-%d %H UTC')
        fileadd = ''

    elif(len(date_str) == 8):
        max_dict = {
            'smoke_conc_sfc': 120,
            'smoke_drysink': 1,
            'smoke_wetsink': 30,
        }
    
        l_event = True
        NAAPS_noAI    = read_NAAPS_event(date_str, minlat = minlat, \
            dtype = 'no_AI')
        NAAPS_withAI  = read_NAAPS_event(date_str, minlat = minlat, \
            dtype = 'with_AI')
        plot_title = 'Event Total NAAPS-RA Values\n' + \
            NAAPS_noAI['dt_begin_date'].strftime('%Y-%m-%d') + ' - ' + \
            NAAPS_noAI['dt_end_date'].strftime('%Y-%m-%d')
        fileadd = '_event'
    else:
        print("ERROR: Invalid date_str. Must be either YYYYMMDD or YYYYMMDDHH")
        return

    # Plot the smoke_conc_sfc data and comps
    # --------------------------------------
    var = 'smoke_conc_sfc'
    vmax = max_dict[var]
    plot_NAAPS(NAAPS_noAI, var, ax = ax1, zoom = zoom, \
            vmin = vmin, vmax = vmax, plot_log = plot_log, \
            circle_bound = circle_bound, minlat = minlat)
    plot_NAAPS(NAAPS_withAI, var, ax = ax2, zoom = zoom, \
            vmin = vmin, vmax = vmax, plot_log = plot_log, \
            circle_bound = circle_bound, minlat = minlat)

    diff_data = NAAPS_withAI[var] - NAAPS_noAI[var]

    labelsize = 12
    mesh = ax3.pcolormesh(NAAPS_withAI['lon'], NAAPS_noAI['lat'], \
        diff_data, \
        cmap = 'bwr', vmin = -(vmax/4), vmax = (vmax/4), \
        transform = ccrs.PlateCarree(), shading = 'auto')
    cbar = plt.colorbar(mesh, ax = ax3, orientation='vertical',\
        pad=0.04, fraction = 0.040, extend = 'both')
    cbar.set_label('Δ' + var, size = labelsize, weight = 'bold')
    ax3.set_boundary(circle, transform=ax3.transAxes)
    ax3.set_extent([-180.0,180.0,minlat,90.0],datacrs)
    ax1.set_title('AOD Only\n' + var)
    ax2.set_title('AOD + AI\n' + var)
    ax3.set_title('Difference (with_AI - no_AI)\n' + var)

    # Plot the smoke_conc_sfc data and comps
    # --------------------------------------
    var = 'smoke_drysink'
    vmax = max_dict[var]
    plot_NAAPS(NAAPS_noAI, var, ax = ax4, zoom = zoom, \
            vmin = vmin, vmax = vmax, plot_log = plot_log, \
            circle_bound = circle_bound, minlat = minlat)
    plot_NAAPS(NAAPS_withAI, var, ax = ax5, zoom = zoom, \
            vmin = vmin, vmax = vmax, plot_log = plot_log, \
            circle_bound = circle_bound, minlat = minlat)

    diff_data = NAAPS_withAI[var] - NAAPS_noAI[var]

    labelsize = 12
    mesh = ax6.pcolormesh(NAAPS_withAI['lon'], NAAPS_noAI['lat'], \
        diff_data, \
        cmap = 'bwr', vmin = -(vmax/4), vmax = (vmax/4), \
        transform = ccrs.PlateCarree(), shading = 'auto')
    cbar = plt.colorbar(mesh, ax = ax6, orientation='vertical',\
        pad=0.04, fraction = 0.040, extend = 'both')
    cbar.set_label('Δ' + var, size = labelsize, weight = 'bold')
    ax6.set_boundary(circle, transform=ax6.transAxes)
    ax6.set_extent([-180.0,180.0,minlat,90.0],datacrs)
    ax4.set_title('AOD Only\n' + var)
    ax5.set_title('AOD + AI\n' + var)
    ax6.set_title('Difference (with_AI - no_AI)\n' + var)

    # Plot the smoke_conc_sfc data and comps
    # --------------------------------------
    var = 'smoke_wetsink'
    vmax = max_dict[var]
    plot_NAAPS(NAAPS_noAI, var, ax = ax7, zoom = zoom, \
            vmin = vmin, vmax = vmax, plot_log = plot_log, \
            circle_bound = circle_bound, minlat = minlat)
    plot_NAAPS(NAAPS_withAI, var, ax = ax8, zoom = zoom, \
            vmin = vmin, vmax = vmax, plot_log = plot_log, \
            circle_bound = circle_bound, minlat = minlat)

    diff_data = NAAPS_withAI[var] - NAAPS_noAI[var]

    labelsize = 12
    mesh = ax9.pcolormesh(NAAPS_withAI['lon'], NAAPS_noAI['lat'], \
        diff_data, \
        cmap = 'bwr', vmin = -(vmax/4), vmax = (vmax/4), \
        transform = ccrs.PlateCarree(), shading = 'auto')
    cbar = plt.colorbar(mesh, ax = ax9, orientation='vertical',\
        pad=0.04, fraction = 0.040, extend = 'both')
    cbar.set_label('Δ' + var, size = labelsize, weight = 'bold')
    ax9.set_boundary(circle, transform=ax9.transAxes)
    ax9.set_extent([-180.0,180.0,minlat,90.0],datacrs)
    ax7.set_title('NAAPS-RA (No AI)\n' + var)
    ax8.set_title('NAAPS-RA (With AI)\n' + var)
    ax9.set_title('Difference (with_AI - no_AI)\n' + var)

    plt.suptitle(plot_title)

 
    ax1.coastlines()
    ax2.coastlines()
    ax3.coastlines()
    ax4.coastlines()
    ax5.coastlines()
    ax6.coastlines()
    ax7.coastlines()
    ax8.coastlines()
    ax9.coastlines()

    fig.tight_layout()
    
    if(save):
        outname = 'naaps' + fileadd + '_comparetype_allvar_' + \
            date_str + '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()






def plot_NAAPS_event_CERES(date_str, var, ceres_var = 'alb_clr', \
        minlat = 70.5, vmin = None, vmax = None, vmin2 = None, vmax2 = None, \
        min_ice = 80., min_smoke = 0, max_smoke = 2e5, plot_log = True, \
        satellite = 'Aqua', ptitle = '', lat_bounds = [-90, 90], \
        lon_bounds = [-180, 180], plot_daily_data = False, \
        zoom = True, save = False):

    plt.close('all')
    fig = plt.figure(figsize = (12,12))
    gs  = fig.add_gridspec(nrows = 3, ncols = 3)
    ax1 = plt.subplot(gs[0,0], projection = mapcrs)
    ax2 = plt.subplot(gs[0,1], projection = mapcrs)
    ax3 = plt.subplot(gs[0,2], projection = mapcrs)
    ax4 = plt.subplot(gs[1,:])
    ax5 = plt.subplot(gs[2,:])
    #ax4 = fig.add_subplot(2,2,4, projection = mapcrs)

    # Read the data for this granule
    # ------------------------------
    ##!#if(date_str[:6] == '201708'):
    ##!#    cdate_begin_str1 = '20170804'
    ##!#    cdate_end_str1   = '20170815'
    ##!#    cdate_begin_str2 = '20170820'
    ##!#    cdate_end_str2   = '20170831'
    ##!#elif(date_str[:6] == '201206'):
    ##!#    cdate_begin_str1 = '20120601'
    ##!#    cdate_end_str1   = '20120613'
    ##!#    cdate_begin_str2 = '20120620'
    ##!#    cdate_end_str2   = '20120630'
    ##!#elif(date_str[:6] == '200804'):
    ##!#    cdate_begin_str1 = '20080407'
    ##!#    cdate_end_str1   = '20080420'
    ##!#    cdate_begin_str2 = '20080426'
    ##!#    cdate_end_str2   = '20080510'
    ##!#elif(date_str[:6] == '202008'):
    ##!#    cdate_begin_str1 = '20200810'
    ##!#    cdate_end_str1   = '20200823'
    ##!#    cdate_begin_str2 = '20200909'
    ##!#    cdate_end_str2   = '20200922'

    cdate_begin_str1 = event_dict[date_str]['before_start']
    cdate_end_str1   = event_dict[date_str]['before_end']
    cdate_begin_str2 = event_dict[date_str]['end_start']
    cdate_end_str2   = event_dict[date_str]['end_end']

# beg: 2020082400
# end: 2020090800

    dt_begin_str1 = datetime.strptime(cdate_begin_str1, '%Y%m%d')
    dt_end_str1   = datetime.strptime(cdate_end_str1, '%Y%m%d')
    dt_begin_str2 = datetime.strptime(cdate_begin_str2, '%Y%m%d')
    dt_end_str2   = datetime.strptime(cdate_end_str2, '%Y%m%d')

    NAAPS_data = read_NAAPS_event(date_str, minlat = minlat)

    CERES_data1, CERES_data2, NAAPS_data = \
        read_zoom_CERES_data(cdate_begin_str1, cdate_end_str1, \
        cdate_begin_str2, cdate_end_str2, ceres_var, min_ice, \
        min_smoke, max_smoke, vmin2, NAAPS_data, minlat, lat_bounds, \
        lon_bounds, satellite = satellite, mask_NAAPS = False)

    print("HERE",CERES_data1[ceres_var].compressed().shape)

    # Plot the data for this granule
    # ------------------------------
    plot_NAAPS(NAAPS_data, var, ax = ax1, zoom = zoom, \
            minlat = minlat, vmin = vmin, vmax = vmax, plot_log = plot_log)

    # Plot the CERES data
    #plotCERES_daily(cdate_begin_str1, ceres_var, end_str = cdate_end_str1, \
    plotCERES_daily(CERES_data1, ceres_var, end_str = cdate_end_str1, \
        satellite = satellite,  only_sea_ice = False, minlat = minlat, \
        vmin = vmin2, vmax = vmax2, \
        avg_data = True, ax = ax2, save = False, min_ice = min_ice, \
        circle_bound = True, colorbar = True)
    plotCERES_daily(CERES_data2, ceres_var, end_str = cdate_end_str2, \
        satellite = satellite,  only_sea_ice = False, minlat = minlat, \
        vmin = vmin2, vmax = vmax2, \
        avg_data = True, ax = ax3, save = False, min_ice = min_ice, \
        circle_bound = True, colorbar = True)
    #plotCERES_daily(cdate_begin_str2, 'ice_conc', end_str = cdate_end_str2, \
    #    satellite = satellite,  only_sea_ice = False, minlat = minlat, \
    #    avg_data = True, ax = ax4, save = False, \
    #    circle_bound = True, colorbar = True)

    # Now, process the yearly data
    # ----------------------------
    second_labels, second_arrays, average_CERES, second_arrays_nosmoke = \
        read_CERES_historic(dt_begin_str1, dt_end_str1, dt_begin_str2, \
        dt_end_str2, ceres_var, min_ice, min_smoke, max_smoke, vmin2, \
        satellite, NAAPS_data, minlat, lat_bounds, lon_bounds, \
        plot_daily_data = plot_daily_data)

    # Read in the NCEP temp data
    # --------------------------
    second_labels_NCEP, second_arrays_NCEP, average_NCEP = \
        read_NCEP_historic(dt_begin_str1, dt_end_str1, dt_begin_str2, \
        dt_end_str2, minlat, lat_bounds, lon_bounds)
 
    ax6 = ax4.twinx()
 
    colors = ['tab:blue','tab:orange','limegreen','tab:red','tab:purple', 'cyan', 'tab:blue', 'tab:orange', 'limegreen'] 
    for ii in range(len(second_arrays)):
        # Plot the smoky data
        # -------------------
        boxs = ax4.boxplot([second_arrays[ii,0], second_arrays[ii,1]], \
            labels = [second_labels[ii,0], second_labels[ii,1]], \
            positions = [ii*2, ii*2 + 1], widths = [0.75, 0.75], \
            patch_artist = True)
        for item in ['boxes','whiskers','fliers','caps']:
            plt.setp(boxs[item], color = colors[ii])
        plt.setp(boxs['medians'], color = 'black')
        for patch in boxs['boxes']:
        #    patch.set_edgecolor(colors[ii])
            patch.set_facecolor('white')

        # Plot the non-smoky data
        # -----------------------
        boxs2 = ax5.boxplot([second_arrays_nosmoke[ii,0], \
            second_arrays_nosmoke[ii,1]], \
            labels = [second_labels[ii,0], second_labels[ii,1]], \
            positions = [ii*2, ii*2 + 1], widths = [0.75, 0.75], \
            patch_artist = True)
        for item in ['boxes','whiskers','fliers','caps']:
            plt.setp(boxs2[item], color = colors[ii])
        plt.setp(boxs2['medians'], color = 'black')
        for patch in boxs2['boxes']:
        #    patch.set_edgecolor(colors[ii])
            patch.set_facecolor('white')

    ax6.plot(average_NCEP)
    ax6.set_ylabel('NCEP Air Temperature')    

    label_dict = {
        'alb_clr': 'Clear-sky Albedo',
        'swf_clr': 'Clear-sky Shortwave Flux',
        'lwf_clr': 'Clear-sky Longwave Flux',
        'alb_all': 'All-sky Albedo',
        'swf_all': 'All-sky Shortwave Flux',
        'lwf_all': 'All-sky Longwave Flux',
    }

    ax4.set_ylabel(label_dict[ceres_var])
    ax5.set_ylabel(label_dict[ceres_var] + ' (No Smoke)')
    
    ax1.set_title('NAAPS-RA ' + var + '\n' + \
        NAAPS_data['dt_begin_date'].strftime('%Y-%m-%d') + ' - ' + \
        NAAPS_data['dt_end_date'].strftime('%Y-%m-%d'))
    ax2.set_title('CERES Average ' + label_dict[ceres_var] + '\n' + cdate_begin_str1 + ' - ' + \
        cdate_end_str1)
    ax3.set_title('CERES Average ' + label_dict[ceres_var] + '\n' + cdate_begin_str2 + ' - ' + \
        cdate_end_str2)
    #ax4.set_title('CERES Average Ice Concentration')
     
    ax1.coastlines()

    fig.tight_layout()

    # = = = = = = = = = = = = = = = = = = =
    #
    # Plot the smoky data histograms
    #
    # = = = = = = = = = = = = = = = = = = =
    fig2 = plt.figure(figsize = (9,9))
    axs = fig2.subplots(nrows = 3, ncols = 3)
    jj = 0
    for ii in range(second_arrays.shape[0]):
        if((ii > 2) & (ii < 6)):
            jj = 1
        elif(ii >= 6):
            jj = 2

        u_stat, p_val = mannwhitneyu(second_arrays[ii,0], second_arrays[ii,1],\
                         alternative = 'greater') 
        print(second_labels[ii,0], len(second_arrays[ii,0]), len(second_arrays[ii,1]))
        axs[jj,ii%3].hist(second_arrays[ii,0], alpha = 0.5, label = 'Before')
        axs[jj,ii%3].hist(second_arrays[ii,1], alpha = 0.5, label = 'After')
        #if(len(second_arrays[ii,0]) > len(second_arrays[ii,1])):
        #    axs[jj,ii%3].hist(second_arrays[ii,1], label = 'After')
        #    axs[jj,ii%3].hist(second_arrays[ii,0], label = 'Before')
        #else:
        #    axs[jj,ii%3].hist(second_arrays[ii,0], label = 'Before')
        #    axs[jj,ii%3].hist(second_arrays[ii,1], label = 'After')
        mean_before = np.nanmean(second_arrays[ii,0])
        mean_after  = np.nanmean(second_arrays[ii,1])
        alb_diff = np.round(mean_after - mean_before, 2)
        axs[jj,ii%3].set_title(second_labels[ii,0])
        #axs[jj,ii%3].set_title(second_labels[ii,0] + '\np_val = ' + \
        #    str(np.round(p_val, 2)) + '\nΔα = ' + str(alb_diff))
        axs[jj,ii%3].set_xlabel(ceres_var)
        axs[jj,ii%3].set_ylabel('Counts')
        axs[jj,ii%3].legend()

    plt.suptitle('Smoky Regions', weight = 'bold', fontsize = 14)

    fig2.tight_layout()

    # = = = = = = = = = = = = = = = = = = =
    #
    # Plot the non-smoky data histograms
    #
    # = = = = = = = = = = = = = = = = = = =
    fig3 = plt.figure(figsize = (9,9))
    axs = fig3.subplots(nrows = 3, ncols = 3)
    jj = 0
    for ii in range(second_arrays_nosmoke.shape[0]):
        if((ii > 2) & (ii < 6)):
            jj = 1
        elif(ii >= 6):
            jj = 2

        u_stat, p_val = mannwhitneyu(second_arrays_nosmoke[ii,0], \
                        second_arrays_nosmoke[ii,1],\
                        alternative = 'greater') 
        print(second_labels[ii,0], len(second_arrays_nosmoke[ii,0]), \
            len(second_arrays_nosmoke[ii,1]))
        axs[jj,ii%3].hist(second_arrays_nosmoke[ii,0], alpha = 0.5, label = 'Before')
        axs[jj,ii%3].hist(second_arrays_nosmoke[ii,1], alpha = 0.5, label = 'After')
        #if(len(second_arrays[ii,0]) > len(second_arrays[ii,1])):
        #    axs[jj,ii%3].hist(second_arrays[ii,1], label = 'After')
        #    axs[jj,ii%3].hist(second_arrays[ii,0], label = 'Before')
        #else:
        #    axs[jj,ii%3].hist(second_arrays[ii,0], label = 'Before')
        #    axs[jj,ii%3].hist(second_arrays[ii,1], label = 'After')
        mean_before = np.nanmean(second_arrays_nosmoke[ii,0])
        mean_after  = np.nanmean(second_arrays_nosmoke[ii,1])
        alb_diff = np.round(mean_after - mean_before, 2)
        axs[jj,ii%3].set_title(second_labels[ii,0])
        #axs[jj,ii%3].set_title(second_labels[ii,0] + '\np_val = ' + \
        #    str(np.round(p_val, 2)) + '\nΔα = ' + str(alb_diff))
        axs[jj,ii%3].set_xlabel(ceres_var)
        axs[jj,ii%3].set_ylabel('Counts')
        axs[jj,ii%3].legend()

    plt.suptitle('Non-smoky Regions', weight = 'bold', fontsize = 14)

    fig3.tight_layout()


    if(save):
        outname = '_'.join(['naaps','ceres','event',var,ceres_var,date_str,'v2.png'])
        #outname = 'naaps_ceres_event_' + var + '_' + ceres_var + '_' + date_str + '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
        #outname = 'naaps_ceres_event_boxplots_' + var + '_' + date_str + '.png'
        outname = '_'.join(['naaps','ceres','event','boxplots',var,ceres_var,date_str,'v2.png'])
        fig2.savefig(outname, dpi = 300)
        print("Saved image", outname)
        #outname = 'naaps_ceres_event_boxplots_' + var + '_nosmoke_' + date_str + '.png'
        outname = '_'.join(['naaps','ceres','event','boxplots',var,'nosmoke',ceres_var,date_str,'v2.png'])
        fig3.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

    return second_arrays, second_labels, second_arrays_nosmoke

#def plot_NAAPS_month_avg(date_str, ax = None, minlat = 65., mask_zero = False, \
#        plot_log = False, save = False):
#
#    dt_date_str = datetime.strptime(date_str, '%Y%m')
#    
#    in_ax = True 
#    if(ax is None):

def plot_NAAPS_event_CERES_region_time_series(date_str, begin_str, end_str, \
        interval, var, ceres_var = 'alb_clr', ax = None, \
        minlat = 70.5, vmin = None, vmax = None, vmin2 = None, vmax2 = None, \
        min_ice = 80., min_smoke = 0, max_smoke = 2e5, plot_log = True, \
        satellite = 'Aqua', ptitle = '', lat_bounds = [-90, 90], \
        lon_bounds = [-180, 180], plot_daily_data = False, \
        zoom = True, save = False):

    dt_date_str    = datetime.strptime(date_str, '%Y%m%d')
    dt_total_begin = datetime.strptime(begin_str, '%Y%m%d')
    dt_total_end   = datetime.strptime(end_str, '%Y%m%d')

# beg: 2020082400
# end: 2020090800

    smoke_begin_str = event_dict[date_str]['start']
    smoke_end_str   = event_dict[date_str]['end']

    dt_smoke_begin = datetime.strptime(smoke_begin_str, '%Y%m%d%H')
    dt_smoke_end   = datetime.strptime(smoke_end_str, '%Y%m%d%H')

    NAAPS_data = read_NAAPS_event(date_str, minlat = minlat)


    ##!## Plot the CERES data
    ##!##plotCERES_daily(cdate_begin_str1, ceres_var, end_str = cdate_end_str1, \
    ##!#plotCERES_daily(CERES_data1, ceres_var, end_str = cdate_end_str1, \
    ##!#    satellite = satellite,  only_sea_ice = False, minlat = minlat, \
    ##!#    vmin = vmin2, vmax = vmax2, \
    ##!#    avg_data = True, ax = ax2, save = False, min_ice = min_ice, \
    ##!#    circle_bound = True, colorbar = True)
    ##!#plotCERES_daily(CERES_data2, ceres_var, end_str = cdate_end_str2, \
    ##!#    satellite = satellite,  only_sea_ice = False, minlat = minlat, \
    ##!#    vmin = vmin2, vmax = vmax2, \
    ##!#    avg_data = True, ax = ax3, save = False, min_ice = min_ice, \
    ##!#    circle_bound = True, colorbar = True)
    #plotCERES_daily(cdate_begin_str2, 'ice_conc', end_str = cdate_end_str2, \
    #    satellite = satellite,  only_sea_ice = False, minlat = minlat, \
    #    avg_data = True, ax = ax4, save = False, \
    #    circle_bound = True, colorbar = True)

    # Now, process the yearly data
    # ----------------------------
    ##!#second_labels, second_arrays, average_CERES, second_arrays_nosmoke = \
    ##!#    read_CERES_historic(dt_begin_str1, dt_end_str1, dt_begin_str2, \
    ##!#    dt_end_str2, ceres_var, min_ice, min_smoke, max_smoke, vmin2, \
    ##!#    satellite, NAAPS_data, minlat, lat_bounds, lon_bounds, \
    ##!#    plot_daily_data = plot_daily_data)

    second_labels, second_arrays, average_CERES, second_arrays_nosmoke = \
        read_CERES_region_time_series(dt_total_begin, dt_total_end, interval, \
        ceres_var, min_ice, min_smoke, max_smoke, vmin2, satellite, \
        NAAPS_data, minlat, lat_bounds, lon_bounds, \
        plot_daily_data = False)

    final_labels = np.copy(second_labels)
    for ii, slbl in enumerate(second_labels):
        date1 = slbl.split()[0]
        date2 = slbl.split()[-1] 
        dt_date1 = datetime.strptime(date1, '%m/%d')
        dt_date2 = datetime.strptime(date2, '%m/%d')
        mid_time = int((dt_date2 - dt_date1).days / 2)
        final_date = dt_date1 + timedelta(days = int(mid_time))
        final_labels[ii] = final_date.strftime('%m/%d')
 
    # Make an array of the final labels in datetime format
    # ----------------------------------------------------
    final_dt_dates = np.array([datetime.strptime(flbl,'%m/%d').replace(\
        year = dt_date_str.year) for flbl in final_labels])

    within_dates = np.where( (final_dt_dates >= dt_smoke_begin) & \
        (final_dt_dates <= dt_smoke_end) )
 
    #str_final_labels = np.array([flbl.strftime('%m/%d') for flbl in final_labels]) 
    smoke_means   = np.array([np.nanmean(tdata) for tdata in second_arrays])
    nosmoke_means = np.array([np.nanmean(tdata) for tdata in \
        second_arrays_nosmoke])

    # Make the figure
    in_ax = True
    num_ticks = 4
    if(ax is None):
        num_ticks = 7
        in_ax = False
        plt.close('all')

        fig = plt.figure()
        #fig = plt.figure(figsize = (12,5))
        #ax1 = fig.add_subplot(1,2,2, projection = mapcrs)
        ax = fig.add_subplot(1,1,1)

    ##!## Plot the data for this granule
    ##!## ------------------------------
    ##!#print(np.nanmax(NAAPS_data[var]), vmin, vmax)
    ##!#plot_NAAPS(NAAPS_data, var, ax = ax1, zoom = zoom, \
    ##!#        minlat = minlat, vmin = vmin, vmax = vmax, plot_log = plot_log)
    ##!#ax1.coastlines()
    ##!#ax1.set_title('NAAPS-RA smoke_conc_sfc\n' + \
    ##!#    NAAPS_data['dt_begin_date'].strftime('%Y-%m-%d') + ' - ' + \
    ##!#    NAAPS_data['dt_end_date'].strftime('%Y-%m-%d'))


    xvals = np.arange(len(smoke_means))
    smoke_times = xvals[within_dates]    

    ax.plot(xvals, smoke_means, label = 'Smoke')
    ax.plot(xvals, nosmoke_means, label = 'No-smoke')
    ax.axvspan(smoke_times[0], smoke_times[-1], \
        color = 'cyan', alpha = 0.5)

    ax.set_xticks(xvals[::int(len(xvals)/num_ticks)])
    ax.set_xticklabels(final_labels[::int(len(xvals)/num_ticks)])

    ax.set_ylabel(ceres_var)
    if(ptitle == ''):
        ax.set_title(dt_date_str.strftime('CERES ' + ceres_var + ' ' + str(interval) + \
           '-day running average\n%d-%B-%Y'))
    else:
        ax.set_title(ptitle)

    ax.legend()

    if(not in_ax):
        if(save):
            outname = '_'.join(['naaps','ceres','time','series','region',ceres_var,\
                begin_str,end_str,'.png'])
            fig.savefig(outname, dpi = 300)
            print("Saved image", outname)
        else:
            plt.show()

def plot_NAAPS_event_CERES_region_time_series_allvars(date_str, begin_str, \
        end_str, interval, var, minlat = 70.5, vmin = None, vmax = None, \
        vmin2 = None, vmax2 = None, min_ice = 80., min_smoke = 0, \
        max_smoke = 2e5, plot_log = True, \
        satellite = 'Aqua', ptitle = '', lat_bounds = [-90, 90], \
        lon_bounds = [-180, 180], plot_daily_data = False, \
        zoom = True, save = False):

    dt_date_str    = datetime.strptime(date_str, '%Y%m%d')
    dt_total_begin = datetime.strptime(begin_str, '%Y%m%d')
    dt_total_end   = datetime.strptime(end_str, '%Y%m%d')

# beg: 2020082400
# end: 2020090800

    smoke_begin_str = event_dict[date_str]['start']
    smoke_end_str   = event_dict[date_str]['end']

    dt_smoke_begin = datetime.strptime(smoke_begin_str, '%Y%m%d%H')
    dt_smoke_end   = datetime.strptime(smoke_end_str, '%Y%m%d%H')

    NAAPS_data = read_NAAPS_event(date_str, minlat = minlat)


    ##!## Plot the CERES data
    ##!##plotCERES_daily(cdate_begin_str1, ceres_var, end_str = cdate_end_str1, \
    ##!#plotCERES_daily(CERES_data1, ceres_var, end_str = cdate_end_str1, \
    ##!#    satellite = satellite,  only_sea_ice = False, minlat = minlat, \
    ##!#    vmin = vmin2, vmax = vmax2, \
    ##!#    avg_data = True, ax = ax2, save = False, min_ice = min_ice, \
    ##!#    circle_bound = True, colorbar = True)
    ##!#plotCERES_daily(CERES_data2, ceres_var, end_str = cdate_end_str2, \
    ##!#    satellite = satellite,  only_sea_ice = False, minlat = minlat, \
    ##!#    vmin = vmin2, vmax = vmax2, \
    ##!#    avg_data = True, ax = ax3, save = False, min_ice = min_ice, \
    ##!#    circle_bound = True, colorbar = True)
    #plotCERES_daily(cdate_begin_str2, 'ice_conc', end_str = cdate_end_str2, \
    #    satellite = satellite,  only_sea_ice = False, minlat = minlat, \
    #    avg_data = True, ax = ax4, save = False, \
    #    circle_bound = True, colorbar = True)

    # Now, process the yearly data
    # ----------------------------
    ##!#second_labels, second_arrays, average_CERES, second_arrays_nosmoke = \
    ##!#    read_CERES_historic(dt_begin_str1, dt_end_str1, dt_begin_str2, \
    ##!#    dt_end_str2, ceres_var, min_ice, min_smoke, max_smoke, vmin2, \
    ##!#    satellite, NAAPS_data, minlat, lat_bounds, lon_bounds, \
    ##!#    plot_daily_data = plot_daily_data)

    # Read albedo data
    second_labels, second_arrays_alb, average_CERES, second_arrays_nosmoke_alb = \
        read_CERES_region_time_series(dt_total_begin, dt_total_end, interval, \
        'alb_clr', min_ice, min_smoke, max_smoke, vmin2, satellite, \
        NAAPS_data, minlat, lat_bounds, lon_bounds, \
        plot_daily_data = False)

    #str_final_labels = np.array([flbl.strftime('%m/%d') for flbl in final_labels]) 
    alb_smoke_means   = np.array([np.nanmean(tdata) for tdata in second_arrays_alb])
    alb_nosmoke_means = np.array([np.nanmean(tdata) for tdata in \
        second_arrays_nosmoke_alb])

    # Read SWF data
    second_labels, second_arrays_swf, average_CERES, second_arrays_nosmoke_swf = \
        read_CERES_region_time_series(dt_total_begin, dt_total_end, interval, \
        'swf_clr', min_ice, min_smoke, max_smoke, vmin2, satellite, \
        NAAPS_data, minlat, lat_bounds, lon_bounds, \
        plot_daily_data = False)

    #str_final_labels = np.array([flbl.strftime('%m/%d') for flbl in final_labels]) 
    swf_smoke_means   = np.array([np.nanmean(tdata) for tdata in second_arrays_swf])
    swf_nosmoke_means = np.array([np.nanmean(tdata) for tdata in \
        second_arrays_nosmoke_swf])

    # Read LWF data
    second_labels, second_arrays_lwf, average_CERES, second_arrays_nosmoke_lwf = \
        read_CERES_region_time_series(dt_total_begin, dt_total_end, interval, \
        'lwf_clr', min_ice, min_smoke, max_smoke, vmin2, satellite, \
        NAAPS_data, minlat, lat_bounds, lon_bounds, \
        plot_daily_data = False)

    #str_final_labels = np.array([flbl.strftime('%m/%d') for flbl in final_labels]) 
    lwf_smoke_means   = np.array([np.nanmean(tdata) for tdata in second_arrays_lwf])
    lwf_nosmoke_means = np.array([np.nanmean(tdata) for tdata in \
        second_arrays_nosmoke_lwf])

    final_labels = np.copy(second_labels)
    for ii, slbl in enumerate(second_labels):
        date1 = slbl.split()[0]
        date2 = slbl.split()[-1] 
        dt_date1 = datetime.strptime(date1, '%m/%d')
        dt_date2 = datetime.strptime(date2, '%m/%d')
        mid_time = int((dt_date2 - dt_date1).days / 2)
        final_date = dt_date1 + timedelta(days = int(mid_time))
        final_labels[ii] = final_date.strftime('%m/%d')
 
    # Make an array of the final labels in datetime format
    # ----------------------------------------------------
    final_dt_dates = np.array([datetime.strptime(flbl,'%m/%d').replace(\
        year = dt_date_str.year) for flbl in final_labels])

    within_dates = np.where( (final_dt_dates >= dt_smoke_begin) & \
        (final_dt_dates <= dt_smoke_end) )
 
    plt.close('all')
    fig = plt.figure(figsize = (12,4))
    ax1 = fig.add_subplot(1,3,1)
    ax2 = fig.add_subplot(1,3,2)
    ax3 = fig.add_subplot(1,3,3)

    xvals = np.arange(len(alb_smoke_means))
    smoke_times = xvals[within_dates]    

    num_ticks = 5

    # Plot ALB data
    ax1.plot(xvals, alb_smoke_means, label = 'Smoke')
    ax1.plot(xvals, alb_nosmoke_means, label = 'No-smoke')
    ax1.axvspan(smoke_times[0], smoke_times[-1], \
        color = 'cyan', alpha = 0.5)
    ax1.set_xticks(xvals[::int(len(xvals)/num_ticks)])
    ax1.set_xticklabels(final_labels[::int(len(xvals)/num_ticks)])
    ax1.set_ylabel('Albedo')
    ax1.set_title(dt_date_str.strftime('CERES Clear-sky Albedo\n' + str(interval) + \
        '-day running average'))
    ax1.legend()

    # Plot SWF data
    ax2.plot(xvals, swf_smoke_means, label = 'Smoke')
    ax2.plot(xvals, swf_nosmoke_means, label = 'No-smoke')
    ax2.axvspan(smoke_times[0], smoke_times[-1], \
        color = 'cyan', alpha = 0.5)
    ax2.set_xticks(xvals[::int(len(xvals)/num_ticks)])
    ax2.set_xticklabels(final_labels[::int(len(xvals)/num_ticks)])
    ax2.set_ylabel('Shortwave Flux [Wm$^{-2}$]')
    ax2.set_title(dt_date_str.strftime('CERES Clear-sky Shortwave Flux\n' + str(interval) + \
        '-day running average'))
    ax2.legend()

    # Plot LWF data
    ax3.plot(xvals, lwf_smoke_means, label = 'Smoke')
    ax3.plot(xvals, lwf_nosmoke_means, label = 'No-smoke')
    ax3.axvspan(smoke_times[0], smoke_times[-1], \
        color = 'cyan', alpha = 0.5)
    ax3.set_xticks(xvals[::int(len(xvals)/num_ticks)])
    ax3.set_xticklabels(final_labels[::int(len(xvals)/num_ticks)])
    ax3.set_ylabel('Longwave Flux [Wm$^{-2}$]')
    ax3.set_title(dt_date_str.strftime('CERES Clear-sky Longwave Flux\n' + str(interval) + \
        '-day running average'))
    ax3.legend()
   
    plt.suptitle(dt_date_str.strftime('%d-%B-%Y Smoke Event'))
 
    fig.tight_layout()

    if(save):
        outname = '_'.join(['naaps','ceres','time','series','region','all_vars',\
            begin_str,end_str + '.png'])
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()


    ##!## Read in the NCEP temp data
    ##!## --------------------------
    ##!#second_labels_NCEP, second_arrays_NCEP, average_NCEP = \
    ##!#    read_NCEP_historic(dt_begin_str1, dt_end_str1, dt_begin_str2, \
    ##!#    dt_end_str2, minlat, lat_bounds, lon_bounds)
 
    ##!#ax6 = ax4.twinx()
 
    ##!#colors = ['tab:blue','tab:orange','limegreen','tab:red','tab:purple', 'cyan', 'tab:blue', 'tab:orange', 'limegreen'] 
    ##!#for ii in range(len(second_arrays)):
    ##!#    # Plot the smoky data
    ##!#    # -------------------
    ##!#    boxs = ax4.boxplot([second_arrays[ii,0], second_arrays[ii,1]], \
    ##!#        labels = [second_labels[ii,0], second_labels[ii,1]], \
    ##!#        positions = [ii*2, ii*2 + 1], widths = [0.75, 0.75], \
    ##!#        patch_artist = True)
    ##!#    for item in ['boxes','whiskers','fliers','caps']:
    ##!#        plt.setp(boxs[item], color = colors[ii])
    ##!#    plt.setp(boxs['medians'], color = 'black')
    ##!#    for patch in boxs['boxes']:
    ##!#    #    patch.set_edgecolor(colors[ii])
    ##!#        patch.set_facecolor('white')

    ##!#    # Plot the non-smoky data
    ##!#    # -----------------------
    ##!#    boxs2 = ax5.boxplot([second_arrays_nosmoke[ii,0], \
    ##!#        second_arrays_nosmoke[ii,1]], \
    ##!#        labels = [second_labels[ii,0], second_labels[ii,1]], \
    ##!#        positions = [ii*2, ii*2 + 1], widths = [0.75, 0.75], \
    ##!#        patch_artist = True)
    ##!#    for item in ['boxes','whiskers','fliers','caps']:
    ##!#        plt.setp(boxs2[item], color = colors[ii])
    ##!#    plt.setp(boxs2['medians'], color = 'black')
    ##!#    for patch in boxs2['boxes']:
    ##!#    #    patch.set_edgecolor(colors[ii])
    ##!#        patch.set_facecolor('white')

    ##!#ax6.plot(average_NCEP)
    ##!#ax6.set_ylabel('NCEP Air Temperature')    

    ##!#label_dict = {
    ##!#    'alb_clr': 'Clear-sky Albedo',
    ##!#    'swf_clr': 'Clear-sky Shortwave Flux',
    ##!#    'lwf_clr': 'Clear-sky Longwave Flux',
    ##!#    'alb_all': 'All-sky Albedo',
    ##!#    'swf_all': 'All-sky Shortwave Flux',
    ##!#    'lwf_all': 'All-sky Longwave Flux',
    ##!#}

    ##!#ax4.set_ylabel(label_dict[ceres_var])
    ##!#ax5.set_ylabel(label_dict[ceres_var] + ' (No Smoke)')
    ##!#
    ##!#ax1.set_title('NAAPS-RA ' + var + '\n' + \
    ##!#    NAAPS_data['dt_begin_date'].strftime('%Y-%m-%d') + ' - ' + \
    ##!#    NAAPS_data['dt_end_date'].strftime('%Y-%m-%d'))
    ##!#ax2.set_title('CERES Average ' + label_dict[ceres_var] + '\n' + cdate_begin_str1 + ' - ' + \
    ##!#    cdate_end_str1)
    ##!#ax3.set_title('CERES Average ' + label_dict[ceres_var] + '\n' + cdate_begin_str2 + ' - ' + \
    ##!#    cdate_end_str2)
    ##!##ax4.set_title('CERES Average Ice Concentration')
    ##!# 
    ##!#ax1.coastlines()

    ##!#fig.tight_layout()

    ##!## = = = = = = = = = = = = = = = = = = =
    ##!##
    ##!## Plot the smoky data histograms
    ##!##
    ##!## = = = = = = = = = = = = = = = = = = =
    ##!#fig2 = plt.figure(figsize = (9,9))
    ##!#axs = fig2.subplots(nrows = 3, ncols = 3)
    ##!#jj = 0
    ##!#for ii in range(second_arrays.shape[0]):
    ##!#    if((ii > 2) & (ii < 6)):
    ##!#        jj = 1
    ##!#    elif(ii >= 6):
    ##!#        jj = 2

    ##!#    u_stat, p_val = mannwhitneyu(second_arrays[ii,0], second_arrays[ii,1],\
    ##!#                     alternative = 'greater') 
    ##!#    print(second_labels[ii,0], len(second_arrays[ii,0]), len(second_arrays[ii,1]))
    ##!#    axs[jj,ii%3].hist(second_arrays[ii,0], alpha = 0.5, label = 'Before')
    ##!#    axs[jj,ii%3].hist(second_arrays[ii,1], alpha = 0.5, label = 'After')
    ##!#    #if(len(second_arrays[ii,0]) > len(second_arrays[ii,1])):
    ##!#    #    axs[jj,ii%3].hist(second_arrays[ii,1], label = 'After')
    ##!#    #    axs[jj,ii%3].hist(second_arrays[ii,0], label = 'Before')
    ##!#    #else:
    ##!#    #    axs[jj,ii%3].hist(second_arrays[ii,0], label = 'Before')
    ##!#    #    axs[jj,ii%3].hist(second_arrays[ii,1], label = 'After')
    ##!#    mean_before = np.nanmean(second_arrays[ii,0])
    ##!#    mean_after  = np.nanmean(second_arrays[ii,1])
    ##!#    alb_diff = np.round(mean_after - mean_before, 2)
    ##!#    axs[jj,ii%3].set_title(second_labels[ii,0] + '\np_val = ' + \
    ##!#        str(np.round(p_val, 2)) + '\nΔα = ' + str(alb_diff))
    ##!#    axs[jj,ii%3].set_xlabel(ceres_var)
    ##!#    axs[jj,ii%3].set_ylabel('Counts')
    ##!#    axs[jj,ii%3].legend()

    ##!#plt.suptitle('Smoky Regions', weight = 'bold', fontsize = 14)

    ##!#fig2.tight_layout()

    ##!## = = = = = = = = = = = = = = = = = = =
    ##!##
    ##!## Plot the non-smoky data histograms
    ##!##
    ##!## = = = = = = = = = = = = = = = = = = =
    ##!#fig3 = plt.figure(figsize = (9,9))
    ##!#axs = fig3.subplots(nrows = 3, ncols = 3)
    ##!#jj = 0
    ##!#for ii in range(second_arrays_nosmoke.shape[0]):
    ##!#    if((ii > 2) & (ii < 6)):
    ##!#        jj = 1
    ##!#    elif(ii >= 6):
    ##!#        jj = 2

    ##!#    u_stat, p_val = mannwhitneyu(second_arrays_nosmoke[ii,0], \
    ##!#                    second_arrays_nosmoke[ii,1],\
    ##!#                    alternative = 'greater') 
    ##!#    print(second_labels[ii,0], len(second_arrays_nosmoke[ii,0]), \
    ##!#        len(second_arrays_nosmoke[ii,1]))
    ##!#    axs[jj,ii%3].hist(second_arrays_nosmoke[ii,0], alpha = 0.5, label = 'Before')
    ##!#    axs[jj,ii%3].hist(second_arrays_nosmoke[ii,1], alpha = 0.5, label = 'After')
    ##!#    #if(len(second_arrays[ii,0]) > len(second_arrays[ii,1])):
    ##!#    #    axs[jj,ii%3].hist(second_arrays[ii,1], label = 'After')
    ##!#    #    axs[jj,ii%3].hist(second_arrays[ii,0], label = 'Before')
    ##!#    #else:
    ##!#    #    axs[jj,ii%3].hist(second_arrays[ii,0], label = 'Before')
    ##!#    #    axs[jj,ii%3].hist(second_arrays[ii,1], label = 'After')
    ##!#    mean_before = np.nanmean(second_arrays_nosmoke[ii,0])
    ##!#    mean_after  = np.nanmean(second_arrays_nosmoke[ii,1])
    ##!#    alb_diff = np.round(mean_after - mean_before, 2)
    ##!#    axs[jj,ii%3].set_title(second_labels[ii,0] + '\np_val = ' + \
    ##!#        str(np.round(p_val, 2)) + '\nΔα = ' + str(alb_diff))
    ##!#    axs[jj,ii%3].set_xlabel(ceres_var)
    ##!#    axs[jj,ii%3].set_ylabel('Counts')
    ##!#    axs[jj,ii%3].legend()

    ##!#plt.suptitle('Non-smoky Regions', weight = 'bold', fontsize = 14)

    ##!#fig3.tight_layout()


    ##!#if(save):
    ##!#    outname = '_'.join(['naaps','ceres','event',var,ceres_var,date_str,'v2.png'])
    ##!#    #outname = 'naaps_ceres_event_' + var + '_' + ceres_var + '_' + date_str + '.png'
    ##!#    fig.savefig(outname, dpi = 300)
    ##!#    print("Saved image", outname)
    ##!#    #outname = 'naaps_ceres_event_boxplots_' + var + '_' + date_str + '.png'
    ##!#    outname = '_'.join(['naaps','ceres','event','boxplots',var,ceres_var,date_str,'v2.png'])
    ##!#    fig2.savefig(outname, dpi = 300)
    ##!#    print("Saved image", outname)
    ##!#    #outname = 'naaps_ceres_event_boxplots_' + var + '_nosmoke_' + date_str + '.png'
    ##!#    outname = '_'.join(['naaps','ceres','event','boxplots',var,'nosmoke',ceres_var,date_str,'v2.png'])
    ##!#    fig3.savefig(outname, dpi = 300)
    ##!#    print("Saved image", outname)
    ##!#else:
    ##!#    plt.show()

def plot_NAAPS_event_CERES_points_time_series(date_str, begin_str, end_str, \
        interval, var, plats, plons, ceres_var = 'alb_clr', \
        minlat = 70.5, vmin = None, vmax = None, vmin2 = None, vmax2 = None, \
        min_ice = 80., min_smoke = 0, max_smoke = 2e5, plot_log = True, \
        satellite = 'Aqua', ptitle = '', lat_bounds = [-90, 90], \
        lon_bounds = [-180, 180], plot_daily_data = False, \
        zoom = True, save = False):

    dt_date_str    = datetime.strptime(date_str, '%Y%m%d')
    dt_total_begin = datetime.strptime(begin_str, '%Y%m%d')
    dt_total_end   = datetime.strptime(end_str, '%Y%m%d')

# beg: 2020082400
# end: 2020090800

    smoke_begin_str = event_dict[date_str]['start']
    smoke_end_str   = event_dict[date_str]['end']

    dt_smoke_begin = datetime.strptime(smoke_begin_str, '%Y%m%d%H')
    dt_smoke_end   = datetime.strptime(smoke_end_str, '%Y%m%d%H')

    NAAPS_data = read_NAAPS_event(date_str, minlat = minlat)

    combined_data = read_CERES_points_time_series(dt_total_begin, \
        dt_total_end, interval, plats, plons, ceres_var, min_ice, \
        0, max_smoke, vmin2, satellite, NAAPS_data, minlat, \
        lat_bounds, lon_bounds, plot_daily_data = False)

    smoke_vals = np.array([combined_data[var][tkey][0] \
        for tkey in combined_data[var].keys()])

    ceres_vals = np.array([[combined_data['ceres_data'][dkey][tkey][0] \
        for tkey in combined_data['ceres_data'][dkey].keys()] \
        for dkey in combined_data['ceres_data'].keys()])
    ceres_vals = np.ma.masked_invalid(ceres_vals)
   
    second_labels = [dkey for dkey in combined_data['ceres_data'].keys()]
 
    #return combined_data

    final_labels = np.copy(second_labels)
    for ii, slbl in enumerate(second_labels):
        date1 = slbl.split()[0]
        date2 = slbl.split()[-1] 
        dt_date1 = datetime.strptime(date1, '%m/%d')
        dt_date2 = datetime.strptime(date2, '%m/%d')
        mid_time = int((dt_date2 - dt_date1).days / 2)
        final_date = dt_date1 + timedelta(days = int(mid_time))
        final_labels[ii] = final_date.strftime('%m/%d')
 
    # Make an array of the final labels in datetime format
    # ----------------------------------------------------
    final_dt_dates = np.array([datetime.strptime(flbl,'%m/%d').replace(\
        year = dt_date_str.year) for flbl in final_labels])

    within_dates = np.where( (final_dt_dates >= dt_smoke_begin) & \
        (final_dt_dates <= dt_smoke_end) )
 
    plt.close('all')
    fig = plt.figure(figsize = (12,5))
    ax1 = fig.add_subplot(1,2,2, projection = mapcrs)
    ax2 = fig.add_subplot(1,2,1)

    # Plot the data for this granule
    # ------------------------------
    print(np.nanmax(NAAPS_data[var]), vmin, vmax)
    plot_NAAPS(NAAPS_data, var, ax = ax1, zoom = zoom, \
            minlat = minlat, vmin = vmin, vmax = vmax, plot_log = plot_log)
    ax1.coastlines()
    ax1.set_title('NAAPS-RA smoke_conc_sfc\n' + \
        NAAPS_data['dt_begin_date'].strftime('%Y-%m-%d') + ' - ' + \
        NAAPS_data['dt_end_date'].strftime('%Y-%m-%d'))


    xvals = np.arange(ceres_vals.shape[0])
    smoke_times = xvals[within_dates]    

    for ii in range(ceres_vals.shape[1]):
        ax2.plot(xvals, ceres_vals[:,ii], \
            label = '{0} μg/m3'.format(int(smoke_vals[ii])))
        plot_point_on_map(ax1, plats[ii], plons[ii], markersize = 10)



    ax2.axvspan(smoke_times[0], smoke_times[-1], \
        color = 'cyan', alpha = 0.5)

    num_ticks = 7
    ax2.set_xticks(xvals[::int(len(xvals)/num_ticks)])
    ax2.set_xticklabels(final_labels[::int(len(xvals)/num_ticks)])

    ax2.set_ylabel(ceres_var)
    ax2.set_title(dt_date_str.strftime('CERES ' + ceres_var + ' ' + str(interval) + \
        '-day running average\n%d-%B-%Y'))

    ax2.legend()

    if(save):
        outname = '_'.join(['naaps','ceres','time','series','points',ceres_var,\
            begin_str,end_str,'.png'])
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

def plot_NAAPS_event_CERES_region_comp(date_str, var, ceres_var = 'alb_clr', \
        minlat = 70.5, vmin = None, vmax = None, vmin2 = None, vmax2 = None, \
        min_ice = 80., min_smoke = 0, max_smoke = 2e5, plot_log = True, \
        satellite = 'Aqua', ptitle = '', lat_bounds = [-90, 90], \
        lon_bounds = [-180, 180], ax = None, plot_daily_data = False, \
        in_year = None, markersize = None, zoom = True, save = False):

    cdate_begin_str1 = event_dict[date_str]['before_start']
    cdate_end_str1   = event_dict[date_str]['before_end']
    cdate_begin_str2 = event_dict[date_str]['end_start']
    cdate_end_str2   = event_dict[date_str]['end_end']

    if(in_year is not None):
        cdate_begin_str1 = str(in_year) + cdate_begin_str1[4:]
        cdate_end_str1   = str(in_year) + cdate_end_str1[4:]
        cdate_begin_str2 = str(in_year) + cdate_begin_str2[4:]
        cdate_end_str2   = str(in_year) + cdate_end_str2[4:]

# beg: 2020082400
# end: 2020090800

    dt_begin_str1 = datetime.strptime(cdate_begin_str1, '%Y%m%d')
    dt_end_str1   = datetime.strptime(cdate_end_str1, '%Y%m%d')
    dt_begin_str2 = datetime.strptime(cdate_begin_str2, '%Y%m%d')
    dt_end_str2   = datetime.strptime(cdate_end_str2, '%Y%m%d')

    NAAPS_data = read_NAAPS_event(date_str, minlat = minlat)

    combined_data = read_CERES_all_region(dt_begin_str1, dt_end_str1, \
        dt_begin_str2, dt_end_str2, ceres_var, \
        min_ice, 0, max_smoke, vmin2, satellite, \
        NAAPS_data, minlat, lat_bounds, lon_bounds, \
        plot_daily_data = False)

    # Calculate the difference between the CERES values before
    # and after the smoke
    # --------------------------------------------------------
    ceres_diff = combined_data['CERES_data2'] - combined_data['CERES_data1']
    
    # Make the figure
    in_ax = True
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

    scat = ax.scatter(combined_data['NAAPS_data'], ceres_diff, s = markersize, \
        #c = combined_data['lats'], \
        )
    #cbar = plt.colorbar(scat,ax=ax,label='Latitude')
    plot_trend_line(ax, combined_data['NAAPS_data'], ceres_diff, \
        color='k', linestyle = '-', slope = 'linregress')

    slope,intercept,r_value,p_value,std_err = \
        stats.linregress(combined_data['NAAPS_data'],ceres_diff)

    ax.set_xlabel('Event Total\nNAAPS ' + var)
    ax.set_ylabel('CERES Δ' + ceres_var)
    if(ptitle == ''):
        ax.set_title(cdate_begin_str1 + ' - ' + cdate_end_str1 + ' vs ' + \
            cdate_begin_str2 + ' - ' + cdate_end_str2 + '\nr$^{2}$ = ' + \
            str(np.round(r_value**2., 3)))
    elif(ptitle == 'RVAL'):
        ax.set_title('r$^{2}$ = ' + str(np.round(r_value**2., 3)))
    else:
        ax.set_title(ptitle + '\nr$^{2}$ = ' + str(np.round(r_value**2.,3)))    
    
    if(not in_ax):
        fig.tight_layout()
        if(save):
            outname = '_'.join(['naaps','ceres','region','diff','comp',ceres_var,\
                date_str, cdate_begin_str1,cdate_begin_str2,'.png'])
            fig.savefig(outname, dpi = 300)
            print("Saved image", outname)
        else:
            plt.show()

def plot_NAAPS_event_CERES_region_comp_twovars(date_str, var, ceres_var1 = 'swf_clr', \
        ceres_var2 = 'lwf_clr',\
        minlat = 70.5, vmin = None, vmax = None, vmin2 = None, vmax2 = None, \
        min_ice = 80., min_smoke = 0, max_smoke = 2e5, plot_log = True, \
        satellite = 'Aqua', ptitle = '', lat_bounds = [-90, 90], \
        lon_bounds = [-180, 180], plot_daily_data = False, \
        in_year = None, markersize = None, zoom = True, save = False):

    cdate_begin_str1 = event_dict[date_str]['before_start']
    cdate_end_str1   = event_dict[date_str]['before_end']
    cdate_begin_str2 = event_dict[date_str]['end_start']
    cdate_end_str2   = event_dict[date_str]['end_end']

    dt_begin_str1 = datetime.strptime(cdate_begin_str1, '%Y%m%d')
    dt_end_str1   = datetime.strptime(cdate_end_str1, '%Y%m%d')
    dt_begin_str2 = datetime.strptime(cdate_begin_str2, '%Y%m%d')
    dt_end_str2   = datetime.strptime(cdate_end_str2, '%Y%m%d')

    # Set up the figure
    plt.close('all')
    fig = plt.figure(figsize = (10, 5))
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
   
    plot_NAAPS_event_CERES_region_comp(date_str, var, ceres_var = ceres_var1, \
        minlat = minlat, vmin = vmin, vmax = vmax, vmin2 = vmin2, vmax2 = vmax2, \
        min_ice = min_ice, min_smoke = min_smoke, max_smoke = max_smoke, \
        plot_log = plot_log, satellite = satellite, ptitle = 'CERES Clear-sky SWF', \
        lat_bounds = lat_bounds, lon_bounds = lon_bounds, ax = ax1, \
        plot_daily_data = plot_daily_data, in_year = in_year, \
        markersize = markersize, save = False)

    plot_NAAPS_event_CERES_region_comp(date_str, var, ceres_var = ceres_var2, \
        minlat = minlat, vmin = vmin, vmax = vmax, vmin2 = vmin2, vmax2 = vmax2, \
        min_ice = min_ice, min_smoke = min_smoke, max_smoke = max_smoke, \
        plot_log = plot_log, satellite = satellite, ptitle = 'CERES Clear-sky LWF', \
        lat_bounds = lat_bounds, lon_bounds = lon_bounds, ax = ax2, \
        plot_daily_data = plot_daily_data, in_year = in_year, \
        markersize = markersize, save = False)

    plt.suptitle('CERES Flux Differences\nPost-smoke (' + \
        cdate_begin_str2[4:] + ' - ' + cdate_end_str2[4:] + \
        ') - Pre-smoke (' + cdate_begin_str1[4:] \
        + ' - ' + cdate_end_str1[4:] + ')')

    ax1.set_ylabel('ΔSWF [Wm$^{-2}$]')
    ax2.set_ylabel('ΔLWF [Wm$^{-2}$]')
    ax1.set_xlabel('Event Total\nNAAPS Near-surface Smoke Aerosol Conc. [μg m$^{-3}$]')
    ax2.set_xlabel('Event Total\nNAAPS Near-surface Smoke Aeroson Conc. [μg m$^{-3}$]')
    ax2.set_ylim(ax1.get_ylim())

    fig.tight_layout()
    if(save):
        outname = '_'.join(['naaps','ceres','region','diff','comp','twovar',\
            date_str, cdate_begin_str1,cdate_begin_str2 + '.png'])
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()


def plot_NAAPS_multi_CERES_region_comp(date_str, var, ceres_var = 'alb_clr', \
        minlat = 70.5, vmin = None, vmax = None, vmin2 = None, vmax2 = None, \
        min_ice = 80., min_smoke = 0, max_smoke = 2e5, plot_log = True, \
        satellite = 'Aqua', ptitle = '', lat_bounds = [-90, 90], \
        lon_bounds = [-180, 180], ax = None, plot_daily_data = False, \
        zoom = True, save = False):

    cdate_begin_str1 = event_dict[date_str]['before_start']
    cdate_end_str1   = event_dict[date_str]['before_end']
    cdate_begin_str2 = event_dict[date_str]['end_start']
    cdate_end_str2   = event_dict[date_str]['end_end']

    dt_begin_str1 = datetime.strptime(cdate_begin_str1, '%Y%m%d')
    dt_end_str1   = datetime.strptime(cdate_end_str1, '%Y%m%d')
    dt_begin_str2 = datetime.strptime(cdate_begin_str2, '%Y%m%d')
    dt_end_str2   = datetime.strptime(cdate_end_str2, '%Y%m%d')

    final_end_date = datetime(year = 2021, month = 9, day = 30)
    check_end_date = dt_end_str2 + relativedelta(years = 4)
    if(check_end_date > final_end_date):
        end_year_offset = check_end_date.year - final_end_date.year
        end_idx = 5 - end_year_offset
        beg_idx = -4 - end_year_offset
    else:
        beg_idx = -4
        end_idx = 5

    rangers = np.arange(beg_idx, end_idx)
    years = rangers + dt_begin_str1.year

    plt.close('all')
    fig2 = plt.figure(figsize = (9,9))
    axs = fig2.subplots(nrows = 3, ncols = 3)
    jj = 0
    for ii in range(len(years)):

        if((ii > 2) & (ii < 6)):
            jj = 1
        elif(ii >= 6):
            jj = 2

        # Call the plotter here
        plot_NAAPS_event_CERES_region_comp(date_str, var, ceres_var = ceres_var, \
            minlat = minlat, vmin = None, vmax = vmax, vmin2 = vmin2, vmax2 = None, \
            min_ice = 80., min_smoke = 0, max_smoke = max_smoke, plot_log = True, \
            satellite = 'All', ptitle = str(years[ii]), lat_bounds = lat_bounds, \
            lon_bounds = lon_bounds, ax = axs[jj,ii%3], plot_daily_data = False, \
            in_year = years[ii], markersize = 12, zoom = True, save = False)

    fig2.tight_layout()

    if(save):
        outname = '_'.join(['naaps','ceres','multi','region','comp',ceres_var,\
            date_str,'.png'])
        fig2.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

def plot_NAAPS_multi_CERES_region_time_series(date_str, begin_date, end_date, \
        interval, var, ceres_var = 'alb_clr', \
        minlat = 70.5, vmin = None, vmax = None, vmin2 = None, vmax2 = None, \
        min_ice = 80., min_smoke = 0, max_smoke = 2e5, plot_log = True, \
        satellite = 'Aqua', ptitle = '', lat_bounds = [-90, 90], \
        lon_bounds = [-180, 180], plot_daily_data = False, \
        zoom = True, save = False):

    #cdate_begin_str1 = event_dict[date_str]['before_start']
    #cdate_end_str1   = event_dict[date_str]['before_end']
    #cdate_begin_str2 = event_dict[date_str]['end_start']
    #cdate_end_str2   = event_dict[date_str]['end_end']

    #dt_begin_str1 = datetime.strptime(cdate_begin_str1, '%Y%m%d')
    #dt_end_str1   = datetime.strptime(cdate_end_str1, '%Y%m%d')
    #dt_begin_str2 = datetime.strptime(cdate_begin_str2, '%Y%m%d')
    #dt_end_str2   = datetime.strptime(cdate_end_str2, '%Y%m%d')

    dt_begin_str1 = datetime.strptime(begin_date, '%Y%m%d')
    dt_end_str1   = datetime.strptime(end_date, '%Y%m%d')

    final_end_date = datetime(year = 2021, month = 9, day = 30)
    check_end_date = dt_end_str1 + relativedelta(years = 4)
    if(check_end_date > final_end_date):
        end_year_offset = check_end_date.year - final_end_date.year
        end_idx = 5 - end_year_offset
        beg_idx = -4 - end_year_offset
    else:
        beg_idx = -4
        end_idx = 5

    rangers = np.arange(beg_idx, end_idx)
    years = rangers + dt_begin_str1.year

    plt.close('all')
    fig2 = plt.figure(figsize = (9,9))
    axs = fig2.subplots(nrows = 3, ncols = 3)
    jj = 0
    for ii in range(len(years)):

        local_begin_str = dt_begin_str1.replace(year = years[ii]).strftime('%Y%m%d')
        local_end_str   = dt_end_str1.replace(year = years[ii]).strftime('%Y%m%d')

        print(local_begin_str, local_end_str)
    
        if((ii > 2) & (ii < 6)):
            jj = 1
        elif(ii >= 6):
            jj = 2

        # Call the plotter here
        #plot_NAAPS_event_CERES_region_comp(date_str, var, ceres_var = ceres_var, \
        #    minlat = minlat, vmin = None, vmax = vmax, vmin2 = vmin2, vmax2 = None, \
        #    min_ice = 80., min_smoke = 0, max_smoke = max_smoke, plot_log = True, \
        #    satellite = 'All', ptitle = str(years[ii]), lat_bounds = lat_bounds, \
        #    lon_bounds = lon_bounds, ax = axs[jj,ii%3], plot_daily_data = False, \
        #    in_year = years[ii], markersize = 12, zoom = True, save = False)

        plot_NAAPS_event_CERES_region_time_series(date_str, local_begin_str, \
            local_end_str, interval, var, ceres_var = ceres_var, \
            ax = axs[jj,ii%3], minlat = minlat, vmin = vmin, vmax = vmax, \
            vmin2 = vmin2, vmax2 = vmax2, min_ice = min_ice, \
            min_smoke = min_smoke, max_smoke = max_smoke, \
            satellite = satellite, ptitle = str(years[ii]), \
            lat_bounds = lat_bounds, lon_bounds = lon_bounds, \
            plot_daily_data = False, zoom = True, save = False)

    fig2.tight_layout()

    if(save):
        outname = '_'.join(['naaps','ceres','multi','region','time','series',\
            ceres_var,date_str + '.png'])
        fig2.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

def plot_NAAPS_event_CERES_region_time_series_combined(date_str, begin_date, end_date, \
        interval, var, ceres_var = 'alb_clr', \
        minlat = 70.5, vmin = None, vmax = None, vmin2 = None, vmax2 = None, \
        min_ice = 80., min_smoke = 0, max_smoke = 2e5, plot_log = True, \
        satellite = 'Aqua', ptitle = '', lat_bounds = [-90, 90], \
        lon_bounds = [-180, 180], ax = None, plot_daily_data = False, \
        zoom = True, save = False):

    dt_date_str  = datetime.strptime(date_str, '%Y%m%d')
    dt_begin_str = datetime.strptime(begin_date, '%Y%m%d')
    dt_end_str   = datetime.strptime(end_date, '%Y%m%d')
    begin_year = 2008
    end_year = 2016

    smoke_begin_str = event_dict[date_str]['start']
    smoke_end_str   = event_dict[date_str]['end']

    dt_smoke_begin = datetime.strptime(smoke_begin_str, '%Y%m%d%H')
    dt_smoke_end   = datetime.strptime(smoke_end_str, '%Y%m%d%H')


    # Check the year ranges
    # ---------------------
    final_end_date = datetime(year = 2021, month = 9, day = 30)
    check_end_date = dt_end_str + relativedelta(years = 4)
    if(check_end_date > final_end_date):
        end_year_offset = check_end_date.year - final_end_date.year
        end_idx = 5 - end_year_offset
        beg_idx = -4 - end_year_offset
    else:
        beg_idx = -4
        end_idx = 5
       
    begin_year = dt_date_str.year + beg_idx 
    end_year   = dt_date_str.year + end_idx - 1
    
    NAAPS_data = read_NAAPS_event(date_str, minlat = minlat)
    
    smoke, nosmoke, labels = read_CERES_region_time_series_multiyear(dt_begin_str, \
        dt_end_str, begin_year, end_year, interval, ceres_var, min_ice, \
        min_smoke, max_smoke, vmin2, 'All', NAAPS_data, minlat, \
        lat_bounds, lon_bounds)
    
    final_labels = np.copy(labels)
    for ii, slbl in enumerate(labels):
        date1 = slbl.split()[0]
        date2 = slbl.split()[-1] 
        dt_date1 = datetime.strptime(date1, '%m/%d')
        dt_date2 = datetime.strptime(date2, '%m/%d')
        mid_time = int((dt_date2 - dt_date1).days / 2)
        final_date = dt_date1 + timedelta(days = int(mid_time))
        final_labels[ii] = final_date.strftime('%m/%d')

    # Extract just the 2012 data
    smoke_keep   = smoke[4,:]
    nosmoke_keep = nosmoke[4,:]
    
    # Get the climo, background data
    smoke_bkgd   = np.concatenate([smoke[:4,:], smoke[5:,:]])
    nosmoke_bkgd = np.concatenate([nosmoke[:4,:], nosmoke[5:,:]])
    
    smoke_bkgd_mean   = np.nanmean(smoke_bkgd, axis = 0)
    nosmoke_bkgd_mean = np.nanmean(nosmoke_bkgd, axis = 0)
    
    smoke_bkgd_stdv   = np.nanstd(smoke_bkgd, axis = 0)
    nosmoke_bkgd_stdv = np.nanstd(nosmoke_bkgd, axis = 0)
    
    smoke_bkgd_upper = smoke_bkgd_mean + smoke_bkgd_stdv
    smoke_bkgd_lower = smoke_bkgd_mean - smoke_bkgd_stdv
    nosmoke_bkgd_upper = nosmoke_bkgd_mean + nosmoke_bkgd_stdv
    nosmoke_bkgd_lower = nosmoke_bkgd_mean - nosmoke_bkgd_stdv

    # Make an array of the final labels in datetime format
    # ----------------------------------------------------
    final_dt_dates = np.array([datetime.strptime(flbl,'%m/%d').replace(\
        year = dt_date_str.year) for flbl in final_labels])

    within_dates = np.where( (final_dt_dates >= dt_smoke_begin) & \
        (final_dt_dates <= dt_smoke_end) )

    xvals = np.arange(len(smoke_keep))
    smoke_times = xvals[within_dates]    

    
    in_ax = True
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

    # Plot the smoky range
    xvals = np.arange(len(smoke_bkgd_upper))
    ax.plot(smoke_bkgd_mean, linestyle = ':', linewidth = 1.00, color = 'tab:blue')
    ax.plot(smoke_bkgd_upper, linestyle = '--', linewidth = 0.5, color = 'tab:blue')
    ax.plot(smoke_bkgd_lower, linestyle = '--', linewidth = 0.5, color = 'tab:blue')
    ax.fill_between(xvals, smoke_bkgd_upper, y2 = smoke_bkgd_lower, \
        color = 'tab:blue', alpha = 0.15)

    # Plot the nosmoky range
    ax.plot(nosmoke_bkgd_mean, linestyle = ':', linewidth = 1.00, color = 'tab:orange')
    ax.plot(nosmoke_bkgd_upper, linestyle = '--', linewidth = 0.5, color = 'tab:orange')
    ax.plot(nosmoke_bkgd_lower, linestyle = '--', linewidth = 0.5, color = 'tab:orange')
    ax.fill_between(xvals, nosmoke_bkgd_upper, y2 = nosmoke_bkgd_lower, \
        color = 'tab:orange', alpha = 0.15)

    ax.axvspan(smoke_times[0], smoke_times[-1], \
        color = 'cyan', alpha = 0.5)

    # Plot the main year data
    num_ticks = 7
    ax.plot(smoke_keep, color = 'tab:blue', label = str(dt_date_str.year) + \
        ' Smoke')
    ax.plot(nosmoke_keep, color = 'tab:orange', label = \
        str(dt_date_str.year) +  ' No-smoke')

    ax.legend()
    ax.set_xticks(xvals[::int(len(xvals)/num_ticks)])
    ax.set_xticklabels(final_labels[::int(len(xvals)/num_ticks)])
    ax.set_ylabel(ceres_var)
    ax.set_title(dt_date_str.strftime('%d-%b-%Y Smoke Event\n%Y CERES ' + \
        'Clear-sky ' + ceres_var + '\nvs. 2008 - 2016 Average'))


    if(not in_ax):
        fig.tight_layout()
        if(save):
            outname = '_'.join(['naaps','ceres','time','series','region',\
                'multiyear',ceres_var, date_str]) + '.png'
            fig.savefig(outname, dpi = 300)
            print("Saved image", outname)
        else:
            plt.show()

def plot_NCEP_region_time_series_combined(date_str, begin_date, end_date, \
        interval, var, \
        minlat = 70.5, vmin = None, vmax = None, vmin2 = None, vmax2 = None, \
        min_ice = 80., min_smoke = 0, max_smoke = 2e5, plot_log = True, \
        satellite = 'Aqua', ptitle = '', lat_bounds = [-90, 90], \
        lon_bounds = [-180, 180], ax = None, plot_daily_data = False, \
        zoom = True, save = False):

    dt_date_str  = datetime.strptime(date_str, '%Y%m%d')
    dt_begin_str = datetime.strptime(begin_date, '%Y%m%d')
    dt_end_str   = datetime.strptime(end_date, '%Y%m%d')
    begin_year = 2008
    end_year = 2016

    smoke_begin_str = event_dict[date_str]['start']
    smoke_end_str   = event_dict[date_str]['end']

    dt_smoke_begin = datetime.strptime(smoke_begin_str, '%Y%m%d%H')
    dt_smoke_end   = datetime.strptime(smoke_end_str, '%Y%m%d%H')


    # Check the year ranges
    # ---------------------
    final_end_date = datetime(year = 2021, month = 9, day = 30)
    check_end_date = dt_end_str + relativedelta(years = 4)
    if(check_end_date > final_end_date):
        end_year_offset = check_end_date.year - final_end_date.year
        end_idx = 5 - end_year_offset
        beg_idx = -4 - end_year_offset
    else:
        beg_idx = -4
        end_idx = 5
       
    begin_year = dt_date_str.year + beg_idx 
    end_year   = dt_date_str.year + end_idx - 1
    
    NAAPS_data = read_NAAPS_event(date_str, minlat = minlat)
   
    # READ NCEP DATA HERE
    # ------------------- 
    smoke, nosmoke, labels = read_NCEP_region_time_series_multiyear(dt_begin_str, \
        dt_end_str, begin_year, end_year, interval, min_ice, \
        min_smoke, max_smoke, vmin2, satellite, NAAPS_data, minlat, \
        lat_bounds, lon_bounds)
    
    final_labels = np.copy(labels)
    for ii, slbl in enumerate(labels):
        date1 = slbl.split()[0]
        date2 = slbl.split()[-1] 
        dt_date1 = datetime.strptime(date1, '%m/%d')
        dt_date2 = datetime.strptime(date2, '%m/%d')
        mid_time = int((dt_date2 - dt_date1).days / 2)
        final_date = dt_date1 + timedelta(days = int(mid_time))
        final_labels[ii] = final_date.strftime('%m/%d')

    # Extract just the 2012 data
    smoke_keep   = smoke[4,:]
    nosmoke_keep = nosmoke[4,:]
    
    # Get the climo, background data
    smoke_bkgd   = np.concatenate([smoke[:4,:], smoke[5:,:]])
    nosmoke_bkgd = np.concatenate([nosmoke[:4,:], nosmoke[5:,:]])
    
    smoke_bkgd_mean   = np.nanmean(smoke_bkgd, axis = 0)
    nosmoke_bkgd_mean = np.nanmean(nosmoke_bkgd, axis = 0)
    
    smoke_bkgd_stdv   = np.nanstd(smoke_bkgd, axis = 0)
    nosmoke_bkgd_stdv = np.nanstd(nosmoke_bkgd, axis = 0)
    
    smoke_bkgd_upper = smoke_bkgd_mean + smoke_bkgd_stdv
    smoke_bkgd_lower = smoke_bkgd_mean - smoke_bkgd_stdv
    nosmoke_bkgd_upper = nosmoke_bkgd_mean + nosmoke_bkgd_stdv
    nosmoke_bkgd_lower = nosmoke_bkgd_mean - nosmoke_bkgd_stdv

    # Make an array of the final labels in datetime format
    # ----------------------------------------------------
    final_dt_dates = np.array([datetime.strptime(flbl,'%m/%d').replace(\
        year = dt_date_str.year) for flbl in final_labels])

    within_dates = np.where( (final_dt_dates >= dt_smoke_begin) & \
        (final_dt_dates <= dt_smoke_end) )

    xvals = np.arange(len(smoke_keep))
    smoke_times = xvals[within_dates]    

    
    in_ax = True
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)

    # Plot the smoky range
    xvals = np.arange(len(smoke_bkgd_upper))
    ax.plot(smoke_bkgd_mean, linestyle = ':', linewidth = 1.00, color = 'tab:blue')
    ax.plot(smoke_bkgd_upper, linestyle = '--', linewidth = 0.5, color = 'tab:blue')
    ax.plot(smoke_bkgd_lower, linestyle = '--', linewidth = 0.5, color = 'tab:blue')
    ax.fill_between(xvals, smoke_bkgd_upper, y2 = smoke_bkgd_lower, \
        color = 'tab:blue', alpha = 0.15)

    # Plot the nosmoky range
    ax.plot(nosmoke_bkgd_mean, linestyle = ':', linewidth = 1.00, color = 'tab:orange')
    ax.plot(nosmoke_bkgd_upper, linestyle = '--', linewidth = 0.5, color = 'tab:orange')
    ax.plot(nosmoke_bkgd_lower, linestyle = '--', linewidth = 0.5, color = 'tab:orange')
    ax.fill_between(xvals, nosmoke_bkgd_upper, y2 = nosmoke_bkgd_lower, \
        color = 'tab:orange', alpha = 0.15)

    ax.axvspan(smoke_times[0], smoke_times[-1], \
        color = 'cyan', alpha = 0.5)

    ax.axhline(273.15, color = 'tab:red')

    # Plot the main year data
    num_ticks = 7
    ax.plot(smoke_keep, color = 'tab:blue', label = str(dt_date_str.year) + \
        ' Smoke Area')
    ax.plot(nosmoke_keep, color = 'tab:orange', label = \
        str(dt_date_str.year) +  ' No-smoke Area')

    ax.legend()
    ax.set_xticks(xvals[::int(len(xvals)/num_ticks)])
    ax.set_xticklabels(final_labels[::int(len(xvals)/num_ticks)])
    ax.set_ylabel('NCEP Reanalysis 2-m Air Temperature')
    ax.set_title(dt_date_str.strftime('%d-%b-%Y Smoke Event\n%Y NCEP ' + \
        'Reanalysis 2-m Air Temp\nvs. 2008 - 2016 Average'))


    if(not in_ax):
        fig.tight_layout()
        if(save):
            outname = '_'.join(['naaps','ncep','time','series','region',\
                'multiyear','airtmp', date_str]) + '.png'
            fig.savefig(outname, dpi = 300)
            print("Saved image", outname)
        else:
            plt.show()



def plot_NAAPS_CERES_flux_diffs(date_str, var, ceres_var = 'alb_clr', \
        minlat = 70.5, vmin = None, vmax = None, vmin2 = None, vmax2 = None, \
        min_ice = 80., min_smoke = 0, max_smoke = 2e5, plot_log = True, \
        satellite = 'Aqua', ptitle = '', lat_bounds = [-90, 90], \
        lon_bounds = [-180, 180], ax = None, plot_daily_data = False, \
        in_year = None, markersize = None, zoom = True, save = False):

    cdate_begin_str1 = event_dict[date_str]['before_start']
    cdate_end_str1   = event_dict[date_str]['before_end']
    cdate_begin_str2 = event_dict[date_str]['end_start']
    cdate_end_str2   = event_dict[date_str]['end_end']

# beg: 2020082400
# end: 2020090800

    dt_begin_str1 = datetime.strptime(cdate_begin_str1, '%Y%m%d')
    dt_end_str1   = datetime.strptime(cdate_end_str1, '%Y%m%d')
    dt_begin_str2 = datetime.strptime(cdate_begin_str2, '%Y%m%d')
    dt_end_str2   = datetime.strptime(cdate_end_str2, '%Y%m%d')

    NAAPS_data = read_NAAPS_event(date_str, minlat = minlat)

    # Used for the single day plotting purposes
    CERES_data1 = readgridCERES_daily(dt_begin_str1.strftime('%Y%m%d'),\
        end_str = dt_end_str1.strftime('%Y%m%d'), \
        satellite = satellite, minlat = minlat)
    CERES_data2 = readgridCERES_daily(dt_begin_str2.strftime('%Y%m%d'),\
        end_str = dt_end_str2.strftime('%Y%m%d'), \
        satellite = satellite, minlat = minlat)
    CERES_diffs = readgridCERES_daily(dt_begin_str2.strftime('%Y%m%d'),\
        end_str = dt_end_str2.strftime('%Y%m%d'), \
        satellite = satellite, minlat = minlat)

    # Find the differences
    CERES_diffs['swf_clr'] = CERES_data2['swf_clr'] - CERES_data1['swf_clr']
    CERES_diffs['lwf_clr'] = CERES_data2['lwf_clr'] - CERES_data1['lwf_clr']

    # Make a figure of the local data
    # -------------------------------
    plt.close('all')
    fig3 = plt.figure(figsize = (10,7))
    ax10 = fig3.add_subplot(3,3,2, projection = ccrs.NorthPolarStereo(\
        central_longitude = np.mean(lon_bounds))) # NAAPS
    ax11 = fig3.add_subplot(3,3,4, projection = ccrs.NorthPolarStereo(\
        central_longitude = np.mean(lon_bounds))) # CERES SWF before
    ax12 = fig3.add_subplot(3,3,5, projection = ccrs.NorthPolarStereo(\
        central_longitude = np.mean(lon_bounds))) # CERES SWF after
    ax13 = fig3.add_subplot(3,3,6, projection = ccrs.NorthPolarStereo(\
        central_longitude = np.mean(lon_bounds))) # CERES SWF diff
    ax21 = fig3.add_subplot(3,3,7, projection = ccrs.NorthPolarStereo(\
        central_longitude = np.mean(lon_bounds))) # CERES LWF before
    ax22 = fig3.add_subplot(3,3,8, projection = ccrs.NorthPolarStereo(\
        central_longitude = np.mean(lon_bounds))) # CERES LWF after
    ax23 = fig3.add_subplot(3,3,9, projection = ccrs.NorthPolarStereo(\
        central_longitude = np.mean(lon_bounds))) # CERES LWF diff

    plot_NAAPS(NAAPS_data, 'smoke_conc_sfc', ax = ax10, zoom = True, \
            minlat = minlat, vmin = 0, vmax = max_smoke, plot_log = False,\
            circle_bound = False)

    # Plot the CERES SWF data before
    plotCERES_daily(CERES_data1, 'swf_clr', end_str = \
        dt_end_str1.strftime('%Y%m%d'), satellite = satellite,  \
        only_sea_ice = False, minlat = minlat, \
        vmin = 220, vmax = 330, \
        avg_data = True, ax = ax11, save = False, min_ice = min_ice, \
        circle_bound = False, colorbar = True)
    # Plot the CERES SWF data after
    plotCERES_daily(CERES_data2, 'swf_clr', end_str = \
        dt_end_str1.strftime('%Y%m%d'), satellite = satellite,  \
        only_sea_ice = False, minlat = minlat, vmin = 220, \
        vmax = 330, \
        avg_data = True, ax = ax12, save = False, min_ice = min_ice, \
        circle_bound = False, colorbar = True)
    plotCERES_daily(CERES_diffs, 'swf_clr', end_str = \
        dt_end_str1.strftime('%Y%m%d'), satellite = satellite,  \
        only_sea_ice = False, minlat = minlat, vmin = -50, \
        vmax = 50, cmap = 'bwr',\
        avg_data = True, ax = ax13, save = False, min_ice = min_ice, \
        circle_bound = False, colorbar = True)
    # Plot the CERES LWF data before
    plotCERES_daily(CERES_data1, 'lwf_clr', end_str = \
        dt_end_str1.strftime('%Y%m%d'), satellite = satellite,  \
        only_sea_ice = False, minlat = minlat, \
        vmin = 210, vmax = 250, \
        avg_data = True, ax = ax21, save = False, min_ice = min_ice, \
        circle_bound = False, colorbar = True)
    # Plot the CERES LWF data after
    plotCERES_daily(CERES_data2, 'lwf_clr', end_str = \
        dt_end_str1.strftime('%Y%m%d'), satellite = satellite,  \
        only_sea_ice = False, minlat = minlat, vmin = 210, \
        vmax = 250, \
        avg_data = True, ax = ax22, save = False, min_ice = min_ice, \
        circle_bound = False, colorbar = True)
    plotCERES_daily(CERES_diffs, 'lwf_clr', end_str = \
        dt_end_str1.strftime('%Y%m%d'), satellite = satellite,  \
        only_sea_ice = False, minlat = minlat, vmin = -20, \
        vmax = 20, cmap = 'bwr', \
        avg_data = True, ax = ax23, save = False, min_ice = min_ice, \
        circle_bound = False, colorbar = True)

    ax10.set_extent([lon_bounds[0], lon_bounds[1], lat_bounds[0], \
        lat_bounds[1]], ccrs.PlateCarree()) 
    ax11.set_extent([lon_bounds[0], lon_bounds[1], lat_bounds[0], \
        lat_bounds[1]], ccrs.PlateCarree()) 
    ax12.set_extent([lon_bounds[0], lon_bounds[1], lat_bounds[0], \
        lat_bounds[1]], ccrs.PlateCarree()) 
    ax13.set_extent([lon_bounds[0], lon_bounds[1], lat_bounds[0], \
        lat_bounds[1]], ccrs.PlateCarree()) 
    ax21.set_extent([lon_bounds[0], lon_bounds[1], lat_bounds[0], \
        lat_bounds[1]], ccrs.PlateCarree()) 
    ax22.set_extent([lon_bounds[0], lon_bounds[1], lat_bounds[0], \
        lat_bounds[1]], ccrs.PlateCarree()) 
    ax23.set_extent([lon_bounds[0], lon_bounds[1], lat_bounds[0], \
        lat_bounds[1]], ccrs.PlateCarree()) 

    ax10.coastlines()
    ax10.set_title('NAAPS-RA smoke_conc_sfc\n' + \
        NAAPS_data['dt_begin_date'].strftime('%Y-%m-%d') + ' - ' + \
        NAAPS_data['dt_end_date'].strftime('%Y-%m-%d'))
    ax11.set_title('CERES Average Clear-sky SWF\n' + \
        dt_begin_str1.strftime('%Y%m%d') + ' - ' + \
        dt_end_str1.strftime('%Y%m%d'))
    ax12.set_title('CERES Average Clear-sky SWF\n' + \
        dt_begin_str2.strftime('%Y%m%d') + ' - ' + \
        dt_end_str2.strftime('%Y%m%d'))
    ax13.set_title('CERES Average Clear-sky SWF\nDifference')
    ax21.set_title('CERES Average Clear-sky LWF\n' + \
        dt_begin_str1.strftime('%Y%m%d') + ' - ' + \
        dt_end_str1.strftime('%Y%m%d'))
    ax22.set_title('CERES Average Clear-sky LWF\n' + \
        dt_begin_str2.strftime('%Y%m%d') + ' - ' + \
        dt_end_str2.strftime('%Y%m%d'))
    ax23.set_title('CERES Average Clear-sky LWF\nDifference')
    fig3.tight_layout()

    if(save):
        outname = '_'.join(['naaps','ceres','flux','diff','spatial',\
            date_str + '.png'])
        fig3.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

def plot_MODIS_data_before_after(date_str, save = False):

    if(date_str == '20120615'):
        before_str = '201206121525'
        after_str  = '201206231505'
    else:
        print("ERROR: No other dates allowed right now")
        return

    dt_before_str = datetime.strptime(before_str, '%Y%m%d%H%M')
    dt_after_str  = datetime.strptime(after_str, '%Y%m%d%H%M')

    plt.close('all')
    fig = plt.figure(figsize = (6, 10.5))
    ax1 = fig.add_subplot(4,2,1, projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ax2 = fig.add_subplot(4,2,2, projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ax3 = fig.add_subplot(4,2,3, projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ax4 = fig.add_subplot(4,2,4, projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ax5 = fig.add_subplot(4,2,5, projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ax6 = fig.add_subplot(4,2,6, projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ax7 = fig.add_subplot(4,2,7, projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ax8 = fig.add_subplot(4,2,8, projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ##!#gs  = fig.add_gridspec(nrows = 4, ncols = 4)
    ##!#ax1 = plt.subplot(gs[0,0], projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ##!#ax2 = plt.subplot(gs[0,1], projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ##!#ax3 = plt.subplot(gs[1,0], projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ##!#ax4 = plt.subplot(gs[1,1], projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ##!#ax5 = plt.subplot(gs[2,0], projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ##!#ax6 = plt.subplot(gs[2,1], projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ##!#ax7 = plt.subplot(gs[3,0], projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ##!#ax8 = plt.subplot(gs[3,1], projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ##!#ax01 = plt.subplot(gs[0:2,2:4])
    ##!#ax02 = plt.subplot(gs[2:4,2:4])
   
    channel = 3 
    lsize = 10
    # Time 1: first column
    plot_MODIS_channel(before_str, 'true_color', swath = True, \
        zoom = True, ax = ax1, plot_borders = True, vmax = None, labelsize = lsize)
    plot_MODIS_channel(before_str, 2, swath = True, \
        zoom = True, ax = ax3, plot_borders = True, vmax = 0.8, labelsize = lsize)
    plot_MODIS_channel(before_str, 5, swath = True, \
        zoom = True, ax = ax5, plot_borders = True, vmax = None, labelsize = lsize)
    plot_MODIS_channel(before_str, 31, swath = True, \
        zoom = True, ax = ax7, plot_borders = True, vmax = None, labelsize = lsize)
    # Time 2: second column
    plot_MODIS_channel(after_str, 'true_color', swath = True, \
        zoom = True, ax = ax2, plot_borders = True, vmax = None, labelsize = lsize)
    plot_MODIS_channel(after_str, 2, swath = True, \
        zoom = True, ax = ax4, plot_borders = True, vmax = 0.8, labelsize = lsize)
    plot_MODIS_channel(after_str, 5, swath = True, \
        zoom = True, ax = ax6, plot_borders = True, vmax = None, labelsize = lsize)
    plot_MODIS_channel(after_str, 31, swath = True, \
        zoom = True, ax = ax8, plot_borders = True, vmax = None, labelsize = lsize)
    #plot_MODIS_channel_time_diff('201206121525','201206231505',1,zoom=True,\
    #    ax = ax3, swath = True, vmin = -0.8, vmax = 0.8,\
    #    circle_bound = False, plot_borders = True)
    #ax1.coastlines()
    #ax2.coastlines()

    # Set up lat/lon points
    # ---------------------
    lons1 = np.linspace(-47.5000,-43.2181,3)
    lons2 = np.linspace(-46.50,-38.9105, 3)
    lons3 = np.linspace(-48.00,-37.0000, 3)
    lons4 = np.linspace(-47.00,-36.0000, 3)
    lats1 = np.full(lons1.shape, 65.0)
    lats2 = np.full(lons2.shape, 66.5)
    lats3 = np.full(lons3.shape, 68.0)
    lats4 = np.full(lons3.shape, 70.0)
    all_plons = np.concatenate([lons1, lons2, lons3, lons4])
    all_plats = np.concatenate([lats1, lats2, lats3, lats4])

    all_plats = np.array([all_plats[1], all_plats[3], all_plats[4], all_plats[6], all_plats[8], all_plats[10]])
    all_plons = np.array([all_plons[1], all_plons[3], all_plons[4], all_plons[6], all_plons[8], all_plons[10]])

    for ii in range(len(all_plons)):
        plot_point_on_map(ax5, all_plats[ii], all_plons[ii], markersize = 10)
    #for ii in range(len(lons2)):
    #    plot_point_on_map(ax5, lats2[ii], lons2[ii], markersize = 10)
    #for ii in range(len(lons3)):
    #    plot_point_on_map(ax5, lats3[ii], lons3[ii], markersize = 10)

    #ax.set_extent([-180, 180, 65, 90], ccrs.PlateCarree())
    ax1.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
    ax2.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
    ax3.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
    ax4.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
    ax5.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
    ax6.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
    ax7.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
    ax8.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
    #ax3.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
    
    #row_label_size = 10
    #fig.text(0.20, 0.99, '15:25 UTC\n12-June-2012', ha='center', va='center', \
    #    weight='bold',fontsize=row_label_size + 1)
    #fig.text(0.70, 0.98, '15:05 UTC\n23-June-2012', ha='center', va='center', \
    #    weight='bold',fontsize=row_label_size + 1)
    
    #plot_MODIS_satpy('201206121525', 'true_color', ax = ax1, var = None, crs = None, \
    #    lons = None, lats = None, lat_lims = None, lon_lims = None, \
    #    vmin = None, vmax = None, ptitle = None, plabel = None, \
    #    labelsize = 10, colorbar = True, swath = True, zoom=False,save=False)
    #plot_MODIS_satpy('201206231505', 'true_color', ax = ax2, var = None, crs = None, \
    #    lons = None, lats = None, lat_lims = None, lon_lims = None, \
    #    vmin = None, vmax = None, ptitle = None, plabel = None, \
    #    labelsize = 10, colorbar = True, swath = True, zoom=False,save=False)
    
    plt.suptitle(dt_before_str.strftime('%H:%M UTC %d-%b-%Y') + ' vs ' + \
        dt_after_str.strftime('%H:%M UTC %d-%b-%Y'), weight = 'bold')

    # Set up the third figure
    # ----------------------------------------------------------------------
    fig1  = plt.figure(figsize = (9, 10.5))
    ax11  = fig1.add_subplot(4,3,1,  projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ax12  = fig1.add_subplot(4,3,2,  projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ax13  = fig1.add_subplot(4,3,3,  projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ax14  = fig1.add_subplot(4,3,4,  projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ax15  = fig1.add_subplot(4,3,5,  projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ax16  = fig1.add_subplot(4,3,6,  projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ax17  = fig1.add_subplot(4,3,7,  projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ax18  = fig1.add_subplot(4,3,8,  projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ax19  = fig1.add_subplot(4,3,9,  projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ax110 = fig1.add_subplot(4,3,10, projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ax111 = fig1.add_subplot(4,3,11, projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ax112 = fig1.add_subplot(4,3,12, projection = ccrs.NorthPolarStereo(central_longitude = 320))

    plot_MODIS_channel(after_str, 1, swath = True, \
        zoom = True, ax = ax11, plot_borders = True, vmax = None, labelsize = lsize)
    plot_MODIS_channel(after_str, 2, swath = True, \
        zoom = True, ax = ax12, plot_borders = True, vmax = 0.8, labelsize = lsize)
    plot_MODIS_channel(after_str, 3, swath = True, \
        zoom = True, ax = ax13, plot_borders = True, vmax = None, labelsize = lsize)
    plot_MODIS_channel(after_str, 4, swath = True, \
        zoom = True, ax = ax14, plot_borders = True, vmax = None, labelsize = lsize)
    plot_MODIS_channel(after_str, 5, swath = True, \
        zoom = True, ax = ax15, plot_borders = True, vmax = None, labelsize = lsize)
    plot_MODIS_channel(after_str, 6, swath = True, \
        zoom = True, ax = ax16, plot_borders = True, vmax = None, labelsize = lsize)
    plot_MODIS_channel(after_str, 7, swath = True, \
        zoom = True, ax = ax17, plot_borders = True, vmax = None, labelsize = lsize)
    plot_MODIS_channel(after_str, 8, swath = True, \
        zoom = True, ax = ax18, plot_borders = True, vmax = None, labelsize = lsize)
    plot_MODIS_channel(after_str, 9, swath = True, \
        zoom = True, ax = ax19, plot_borders = True, vmax = None, labelsize = lsize)
    plot_MODIS_channel(after_str, 10, swath = True, \
        zoom = True, ax = ax110, plot_borders = True, vmax = None, labelsize = lsize)
    plot_MODIS_channel(after_str, 11, swath = True, \
        zoom = True, ax = ax111, plot_borders = True, vmax = None, labelsize = lsize)
    plot_MODIS_channel(after_str, 12, swath = True, \
        zoom = True, ax = ax112, plot_borders = True, vmax = None, labelsize = lsize)

    ax11.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
    ax12.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
    ax13.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
    ax14.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
    ax15.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
    ax16.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
    ax17.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
    ax18.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
    ax19.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
    ax110.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
    ax111.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
    ax112.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 

    # Set up the second figure: the spectral reflectances before
    # and after
    # ----------------------------------------------------------

    #plat1 = all_plats[0]
    #plon1 = all_plons[1]

    #channels = np.array([8,9,3,10,11,4,1,15,16,5,6,7])
    channels = np.array([3,4,1,2,5,6,7])
    before_refl = np.full((len(channels), len(all_plats)), np.nan)
    after_refl  = np.full((len(channels), len(all_plats)), np.nan)

    wavelengths = np.array([np.mean(channel_dict[str(ch)]['Bandwidth']) \
        for ch in channels])

    fig2 = plt.figure(figsize = (6,8))
    ax01 = fig2.add_subplot(2,1,1)
    ax02 = fig2.add_subplot(2,1,2)

    for ii in range(len(channels)):
        MODIS_data = read_MODIS_channel(before_str, channels[ii], swath = True)
        #var1, crs1, lons1, lats1, lat_lims1, lon_lims1, plabel1 = \
        #    read_MODIS_satpy(before_str, channels[ii])
        #m_idx = nearest_gridpoint(plat1, plon1, np.array(lats1), np.array(lons1))
        #before_refl[ii] = np.array(var1)[m_idx] 

        #m_idx = nearest_gridpoint(plat1, plon1, MODIS_data['lat'], MODIS_data['lon'])
        before_refl[ii, :] = np.array([MODIS_data['data'][nearest_gridpoint(\
            plat, plon, MODIS_data['lat'], MODIS_data['lon'])] for plat, \
            plon in zip(all_plats, all_plons)]).squeeze()
        #before_refl[ii] = MODIS_data['data'][m_idx]
       
        #var1, crs1, lons1, lats1, lat_lims1, lon_lims1, plabel1 = \
        #    read_MODIS_satpy(after_str, channels[ii])
        #m_idx = nearest_gridpoint(plat1, plon1, np.array(lats1), np.array(lons1))
        #after_refl[ii] = np.array(var1)[m_idx] 
        
        MODIS_data = read_MODIS_channel(after_str, channels[ii], swath = True)
        #m_idx = nearest_gridpoint(plat1, plon1, MODIS_data['lat'], MODIS_data['lon'])
        #after_refl[ii] = MODIS_data['data'][m_idx]
        after_refl[ii,:] = np.array([MODIS_data['data'][nearest_gridpoint(\
            plat, plon, MODIS_data['lat'], MODIS_data['lon'])] for plat, \
            plon in zip(all_plats, all_plons)]).squeeze()
   
    for jj in range(len(all_plats)): 
        line = ax01.plot(wavelengths, before_refl[:,jj]) 
        ax01.plot(wavelengths, after_refl[:,jj], linestyle = '--', \
            color = line[-1].get_color()) 

        ax02.plot(wavelengths, after_refl[:,jj] - before_refl[:,jj], \
            color = line[-1].get_color()) 

    ax01.set_xlabel('Wavelength [μm]')
    ax02.set_xlabel('Wavelength [μm]')
    ax01.set_ylabel('Reflectance')
    ax02.set_ylabel('ΔRefl')
    ax01.set_title(dt_before_str.strftime('Solid - Before (%d-%B-%Y)\n') + \
                   dt_after_str.strftime('Dashed - After  (%d-%B-%Y)\n'))
    fig2.tight_layout()
 
    fig.tight_layout()
    fig1.tight_layout()

    if(save):
        outname =  'modis_imagery_naaps_ceres_' + date_str + '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()
    
def plot_MODIS_OMI_data(date_str, save = False):

    if(date_str == '20120615'):
        work_str = '201206151420'
        omi_date = '201206151328'
        dt_work_str = datetime.strptime(work_str, '%Y%m%d%H%M')
    else:
        print("ERROR: No other dates allowed right now")
        return

    plt.close('all')
    fig = plt.figure(figsize = (9, 4.2))
    ax1 = fig.add_subplot(1,2,1, projection = ccrs.NorthPolarStereo(central_longitude = 320))
    ax2 = fig.add_subplot(1,2,2, projection = ccrs.NorthPolarStereo(central_longitude = 320))
    
    channel = 3 
    lsize = 10
    # Time 1: first column
    plot_MODIS_channel(work_str, 'true_color', swath = True, \
        zoom = True, ax = ax1, plot_borders = True, vmax = None, labelsize = lsize)

    plotOMI_single_swath_figure(omi_date, dtype = 'shawn',  \
        only_sea_ice = False, minlat = 60., skiprows = None, \
        lat_circles = None, save = False, zoom = False, \
        vmin = -1, vmax = 4, circle_bound = False, ax = ax2, \
        shawn_path = home_dir + '/data/OMI/shawn_files/ltc3/')

    #plot_MODIS_channel_time_diff('201206121525','201206231505',1,zoom=True,\
    #    ax = ax3, swath = True, vmin = -0.8, vmax = 0.8,\
    #    circle_bound = False, plot_borders = True)
    #ax1.coastlines()
    #ax2.coastlines()
    #ax.set_extent([-180, 180, 65, 90], ccrs.PlateCarree())
    ax1.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
    ax2.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
    #ax3.set_extent([305, 335, 60, 75], ccrs.PlateCarree()) 
    
    ax1.set_title('Aqua MODIS True Color')
    ax2.set_title('OMI UVAI')
     
    #plt.suptitle(dt_before_str.strftime('%H:%M UTC %d-%b-%Y') + ' vs ' + \
    #    dt_after_str.strftime('%H:%M U%C %d-%b-%Y'), weight = 'bold')
    
    fig.tight_layout()
   
    if(save):
        outname =  'modis_omi_imagery_naaps_smoke_' + date_str + '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()
    
