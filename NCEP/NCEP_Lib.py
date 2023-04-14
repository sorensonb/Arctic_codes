"""
  NAME:
    NCEPLib

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
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from glob import glob
import os
import importlib
from scipy.interpolate import RegularGridInterpolator

home_dir = os.environ['HOME']
sys.path.append(home_dir)
from python_lib import *

data_dir = home_dir + '/data/NCEP/'

if(home_dir + '/Research/NAAPS/' not in sys.path):
    sys.path.append(home_dir + '/Research/NAAPS/')
if(home_dir + '/Research/CERES/' not in sys.path):
    sys.path.append(home_dir + '/Research/CERES/')

from NAAPSLib import *
from gridCERESLib import *

datacrs = ccrs.PlateCarree()
mapcrs = ccrs.NorthPolarStereo()

def read_NCEP(begin_date, end_date, minlat = 65.):

    base_date = datetime(year = 1800, month = 1, day = 1)

    dt_begin_str = datetime.strptime(begin_date, '%Y%m%d')
    dt_end_str   = datetime.strptime(end_date, '%Y%m%d')

    # Calculate the number of hours between base and each date
    # --------------------------------------------------------
    begin_hours = divmod((dt_begin_str - base_date).total_seconds(), 3600)[0]
    end_hours = divmod((dt_end_str - base_date).total_seconds(), 3600)[0]

    print(dt_begin_str, dt_end_str, data_dir) 

    filename = dt_begin_str.strftime(home_dir + '/data/NCEP/air.sig995.%Y.nc')
    #filename = dt_begin_str.strftime(data_dir + 'air.sig995.%Y.nc')

    data = Dataset(filename, 'r')

    # Extract the time indices
    # ------------------------
    tidx = np.where((data['time'][:] >= begin_hours) & (data['time'][:] <= end_hours))

    # Pull out and average the data
    # -----------------------------
    arctic_air = np.nanmean(data['air'][tidx[0],:,:], axis = 0)

    # Only keep the data within the minlat range
    # ------------------------------------------
    idx1 = np.where(data['lat'][:] >= minlat)[0][0]
    idx2 = np.where(data['lat'][:] >= minlat)[0][-1] + 1
    xx, yy = np.meshgrid(data['lon'][:], data['lat'][idx1:idx2])

    arctic_air = arctic_air[idx1:idx2,:]

    data.close()

    NCEP_data = {}
    NCEP_data['filename']       = filename
    NCEP_data['begin_date']     = begin_date
    NCEP_data['end_date']       = end_date
    NCEP_data['data']           = arctic_air
    NCEP_data['lon']            = xx
    NCEP_data['lat']            = yy
    NCEP_data['dt_begin_date'] = dt_begin_str
    NCEP_data['dt_end_date']   = dt_end_str

    #over_180    = np.where(xx[0,:] >= 180.)
    #under_180   = np.where(xx[0,:] < 180.)


    #for ii in range(yy.shape[0]):
    #    NCEP_data['data'][ii,:] = \
    #        np.concatenate([NCEP_data['data'][ii,:][over_180],\
    #        NCEP_data['data'][ii,:][under_180]])
    #    NCEP_data['lat'][ii,:] = \
    #        np.concatenate([NCEP_data['lat'][ii,:][over_180],\
    #        NCEP_data['lat'][ii,:][under_180]])
    #    NCEP_data['lon'][ii,:] = \
    #        np.concatenate([NCEP_data['lon'][ii,:][over_180] - 360.,\
    #        NCEP_data['lon'][ii,:][under_180]])
 
    return NCEP_data

# NAAPS data is only passed and returned so that, if desired,
# the NAAPS data can be masked in the regions where the CERES data
# are masked for not having high enough ice or smoke.
def read_zoom_NCEP_data_single(date1, date2, min_ice, \
        min_smoke, max_smoke, vmin2, NAAPS_data, minlat, lat_bounds, \
        lon_bounds, satellite = 'All', mask_NAAPS = False, \
        plot_daily_data = False):

    NCEP_data = read_NCEP(date1, date2, minlat = minlat)



    # Read the corresponding CERES data
    CERES_data1 = readgridCERES_daily(date1,end_str = date2, \
        satellite = satellite, minlat = minlat)
    raw_CERES1_alb = np.nanmean(CERES_data1['alb_clr'], axis = 0)
    raw_CERES1_ice = np.nanmean(CERES_data1['ice_conc'], axis = 0)

    # Regrid the CERES and NAAPS data
    interper_NAAPS2NCEP = RegularGridInterpolator((NAAPS_data['lat'][:,0],\
        NAAPS_data['lon'][0,:]), NAAPS_data['smoke_conc_sfc'])
    interper_ice2NCEP = RegularGridInterpolator((CERES_data1['lat'][:,0],\
        CERES_data1['lon'][0,:]), raw_CERES1_ice)
    interper_alb2NCEP = RegularGridInterpolator((CERES_data1['lat'][:,0],\
        CERES_data1['lon'][0,:]), raw_CERES1_alb)

    # Trim the NCEP lats and lons to ensure that they do not go outside
    # the bounds of the NAAPS data or CERES data
    minlat_idx = np.where((NCEP_data['lat'][:,0] >= minlat) & \
        (NCEP_data['lat'][:,0] >= np.min(NAAPS_data['lat'])) & \
        (NCEP_data['lat'][:,0] <= np.max(NAAPS_data['lat'])))
    trim_NCEP_lat1d = NCEP_data['lat'][minlat_idx,0].squeeze()

    minlon_idx = np.where(\
        (NCEP_data['lon'][0,:] >= np.min(NAAPS_data['lon'])) & \
        (NCEP_data['lon'][0,:] <= np.max(NAAPS_data['lon'])))
    trim_NCEP_lon1d = NCEP_data['lon'][0,minlon_idx].squeeze()

    mesh_match_lons, mesh_match_lats = np.meshgrid(trim_NCEP_lon1d, \
        trim_NCEP_lat1d)

    NCEP_data['data'] = NCEP_data['data'][minlat_idx[0][0]:minlat_idx[0][-1]+1, \
        minlon_idx[0][0]:minlon_idx[0][-1]+1]
    NCEP_data['lon'] = NCEP_data['lon'][minlat_idx[0][0]:minlat_idx[0][-1]+1, \
        minlon_idx[0][0]:minlon_idx[0][-1]+1]
    NCEP_data['lat'] = NCEP_data['lat'][minlat_idx[0][0]:minlat_idx[0][-1]+1, \
        minlon_idx[0][0]:minlon_idx[0][-1]+1]
    work_NCEP     = NCEP_data['data']

    avg_NAAPS     = interper_NAAPS2NCEP((mesh_match_lats, mesh_match_lons))
    avg_CERES1_alb = interper_alb2NCEP((mesh_match_lats, mesh_match_lons))
    avg_CERES1_ice = interper_ice2NCEP((mesh_match_lats, mesh_match_lons))

    # Test masking the data
    # ---------------------

    # Test pulling out the data that only exists in both zones
    # --------------------------------------------------------
    #avg_CERES1 = np.ma.masked_where((avg_CERES1_ice <= min_ice) | \
    if(lon_bounds[0] < lon_bounds[1]):

        # Contains the smoky data
        NCEP_data['data'] = np.ma.masked_where(\
            (avg_CERES1_ice <= min_ice) | \
            (avg_NAAPS < min_smoke) | \
            (avg_NAAPS > max_smoke), work_NCEP)
        NCEP_data['data'] = np.ma.masked_where(~ \
            ((NCEP_data['lon'] >= lon_bounds[0]) & \
             (NCEP_data['lon'] <= lon_bounds[1]) & \
             (NCEP_data['lat'] >= lat_bounds[0]) & \
             (NCEP_data['lat'] <= lat_bounds[1])), NCEP_data['data'])
        NCEP_data['data'] = np.ma.masked_where(\
            avg_CERES1_alb < vmin2, NCEP_data['data'])

        # Contains the non-smoky data
        NCEP_data['data_nosmoke'] = np.ma.masked_where((\
            avg_CERES1_ice <= min_ice) | \
            (avg_NAAPS > min_smoke), work_NCEP)
        NCEP_data['data_nosmoke'] = np.ma.masked_where(~ \
            ((NCEP_data['lon'] >= lon_bounds[0]) & \
             (NCEP_data['lon'] <= lon_bounds[1]) & \
             (NCEP_data['lat'] >= lat_bounds[0]) & \
             (NCEP_data['lat'] <= lat_bounds[1])) , \
            NCEP_data['data_nosmoke'])
        NCEP_data['data_nosmoke'] = np.ma.masked_where(\
            avg_CERES1_alb < vmin2, NCEP_data['data_nosmoke'])

        avg_CERES1_ice = np.ma.masked_where(avg_CERES1_ice < min_ice, avg_CERES1_ice)
        avg_CERES1_alb = np.ma.masked_where(avg_CERES1_alb < vmin2, avg_CERES1_alb)

        ##!#if(mask_NAAPS):
        ##!#   # Smoky data
        ##!#   NAAPS_data['smoke_conc_sfc'] = np.ma.masked_where(\
        ##!#       (avg_CERES1_ice <= min_ice) | \
        ##!#       (NAAPS_data['smoke_conc_sfc'] < min_smoke) | \
        ##!#       (NAAPS_data['smoke_conc_sfc'] > max_smoke), \
        ##!#       NAAPS_data['smoke_conc_sfc'])
        ##!#   NAAPS_data['smoke_conc_sfc'] = np.ma.masked_where(~ \
        ##!#       ((CERES_data1['lon'] >= lon_bounds[0]) & \
        ##!#       (CERES_data1['lon'] <= lon_bounds[1]) & \
        ##!#       (CERES_data1['lat'] >= lat_bounds[0]) & \
        ##!#       (CERES_data1['lat'] <= lat_bounds[1])), \
        ##!#       NAAPS_data['smoke_conc_sfc'])
           # Non-Smoky data

        #NAAPS_data['smoke_conc_sfc'] = np.ma.masked_where(CERES_data1['alb_clr'] < vmin2, NAAPS_data['smoke_conc_sfc'])

    else:
        # Contains the smoky data
        NCEP_data['data'] = np.ma.masked_where((\
            avg_CERES1_ice <= min_ice) | \
            (avg_NAAPS < min_smoke) | \
            (avg_NAAPS > max_smoke), work_NCEP)
        NCEP_data['data'] = np.ma.masked_where(~\
            (((NCEP_data['lon'] >= lon_bounds[0]) | \
              (NCEP_data['lon'] <= lon_bounds[1])) & \
              (NCEP_data['lat'] >= lat_bounds[0]) & \
              (NCEP_data['lat'] <= lat_bounds[1])) , NCEP_data['data'])
        NCEP_data['data'] = np.ma.masked_where(\
            avg_CERES1_alb < vmin2, NCEP_data['data'])
        # Contains the non-smoky data
        NCEP_data['data_nosmoke'] = np.ma.masked_where((\
            avg_CERES1_ice <= min_ice) | \
            (avg_NAAPS > min_smoke), work_NCEP)
        NCEP_data['data_nosmoke'] = np.ma.masked_where(~\
            (((NCEP_data['lon'] >= lon_bounds[0]) | \
              (NCEP_data['lon'] <= lon_bounds[1])) & \
              (NCEP_data['lat'] >= lat_bounds[0]) & \
              (NCEP_data['lat'] <= lat_bounds[1])) , NCEP_data['data_nosmoke'])
        NCEP_data['data_nosmoke'] = np.ma.masked_where(\
            avg_CERES1_alb < vmin2, NCEP_data['data_nosmoke'])

        ##!#if(mask_NAAPS):
        ##!#    NAAPS_data['smoke_conc_sfc'] = np.ma.masked_where(\
        ##!#        (avg_CERES1_ice <= min_ice) | \
        ##!#        (NAAPS_data['smoke_conc_sfc'] < min_smoke) | \
        ##!#        (NAAPS_data['smoke_conc_sfc'] > max_smoke), \
        ##!#        NAAPS_data['smoke_conc_sfc'])
        ##!#    NAAPS_data['smoke_conc_sfc'] = np.ma.masked_where(~ (((CERES_data1['lon'] >= lon_bounds[0]) | \
        ##!#        (CERES_data1['lon'] <= lon_bounds[1])) & (CERES_data1['lat'] >= lat_bounds[0]) & \
        ##!#        (CERES_data1['lat'] <= lat_bounds[1])), NAAPS_data['smoke_conc_sfc'])
        ##!#    #NAAPS_data['smoke_conc_sfc'] = np.ma.masked_where(CERES_data1['alb_clr'] < vmin2, NAAPS_data['smoke_conc_sfc'])
    if(plot_daily_data):
        # Make a figure of the local data
        # -------------------------------
        plt.close('all')
        fig3 = plt.figure(figsize = (11,7))
        ax1 = fig3.add_subplot(2,4,1, projection = mapcrs)
        ax2 = fig3.add_subplot(2,4,2, projection = mapcrs)
        ax3 = fig3.add_subplot(2,4,3, projection = mapcrs)
        ax4 = fig3.add_subplot(2,4,4, projection = mapcrs)
        ax5 = fig3.add_subplot(2,4,5, projection = mapcrs)
        ax6 = fig3.add_subplot(2,4,6, projection = mapcrs)
        ax7 = fig3.add_subplot(2,4,7, projection = mapcrs)
        ax8 = fig3.add_subplot(2,4,8, projection = mapcrs)
   
        #      ICE   ALB   NAAPS   NCEP
        # raw   1     2      3      4
        #   
        # smth  5     6      7      8
        #   
        
        ax1.pcolormesh(CERES_data1['lon'], CERES_data1['lat'], raw_CERES1_ice, \
            transform = datacrs, shading = 'auto')
        ax1.set_extent([-180, 180, 60, 90], datacrs)
        ax1.coastlines()

        ax2.pcolormesh(CERES_data1['lon'], CERES_data1['lat'], raw_CERES1_alb, \
            transform = datacrs, shading = 'auto')
        ax2.set_extent([-180, 180, 60, 90], datacrs)
        ax2.coastlines()

        ax3.pcolormesh(NAAPS_data['lon'], NAAPS_data['lat'], NAAPS_data['smoke_conc_sfc'], \
            transform = datacrs, shading = 'auto', vmax = 50)
        ax3.set_extent([-180, 180, 60, 90], datacrs)
        ax3.coastlines()

        ax4.pcolormesh(NCEP_data['lon'], NCEP_data['lat'], NCEP_data['data'], \
            transform = datacrs, shading = 'auto')
        ax4.set_extent([-180, 180, 60, 90], datacrs)
        ax4.coastlines()

        ax5.pcolormesh(mesh_match_lons, mesh_match_lats, avg_CERES1_ice, \
            transform = datacrs, shading = 'auto')
        ax5.set_extent([-180, 180, 60, 90], datacrs)
        ax5.coastlines()

        ax6.pcolormesh(mesh_match_lons, mesh_match_lats, avg_CERES1_alb, \
            transform = datacrs, shading = 'auto')
        ax6.set_extent([-180, 180, 60, 90], datacrs)
        ax6.coastlines()

        ax7.pcolormesh(mesh_match_lons, mesh_match_lats, avg_NAAPS, \
            transform = datacrs, shading = 'auto', vmax = 50)
        ax7.set_extent([-180, 180, 60, 90], datacrs)
        ax7.coastlines()

        ax8.pcolormesh(NCEP_data['lon'], NCEP_data['lat'], NCEP_data['data_nosmoke'], \
            transform = datacrs, shading = 'auto')
        ax8.set_extent([-180, 180, 60, 90], datacrs)
        ax8.coastlines()

        fig3.tight_layout()
    
        plt.show()

    return NCEP_data



# dt_begin_str = beginning of time window for traveling averaging
# dt_end_str   = end of time window for traveling averaging
# interval     = length of time window for averaging (in days)
# ceres_var    = 'alb_clr', 'swf_clr', 'lwf_clr'
# min_ice      = minimum NSIDC ice threshold for analysis
# plats        = point lats for pulling data at a single grid point
# plons        = point lons for pulling data at a single grid point
def read_NCEP_region_time_series(dt_begin_str, dt_end_str, interval, \
        min_ice, min_smoke, max_smoke, vmin2, satellite, \
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

        NCEP_data = read_zoom_NCEP_data_single(\
            dt_begin_local1.strftime('%Y%m%d'), \
            dt_end_local1.strftime('%Y%m%d'), \
            min_ice, \
            min_smoke, max_smoke, vmin2, NAAPS_data, minlat, lat_bounds, \
            lon_bounds, satellite = 'All', mask_NAAPS = False, \
            plot_daily_data = False)

        bdata1 = NCEP_data['data'].flatten().compressed()
        bdata3 = NCEP_data['data_nosmoke'].flatten().compressed()
  
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



def read_NCEP_region_time_series_multiyear(dt_begin_str, \
        dt_end_str, begin_year, end_year, interval, min_ice, \
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

        second_labels, second_arrays, average_NCEP, second_arrays_nosmoke = \
            read_NCEP_region_time_series(dt_local_begin, dt_local_end, \
            interval, min_ice, min_smoke, max_smoke, vmin2, \
            satellite, NAAPS_data, minlat, lat_bounds, lon_bounds, \
            plot_daily_data = False)

        smoke_means   = np.array([np.nanmean(tdata) for tdata in second_arrays])
        nosmoke_means = np.array([np.nanmean(tdata) for tdata in \
            second_arrays_nosmoke])

        total_smoke_arrays[ii,:] = smoke_means
        total_nosmoke_arrays[ii,:] = nosmoke_means
   

    return total_smoke_arrays, total_nosmoke_arrays, second_labels
    #return second_arrays, second_arrays_nosmoke

##!#def read_NCEP_event(date_str, minlat = 65.):
##!# 
##!#    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H')
##!#
##!#    # Get the event starting and ending times
##!#    # ---------------------------------------
##!#    begin_str = event_dict[date_str]['start']
##!#    end_str   = event_dict[date_str]['end']
##!#
##!#    dt_begin_str = datetime.strptime(begin_str, '%Y%m%d%H')
##!#    dt_end_str   = datetime.strptime(end_str, '%Y%m%d%H')
##!#
##!#    # Grab the files from the data directory
##!#    # --------------------------------------
##!#    files = subprocess.check_output(dt_date_str.strftime('ls ' + data_dir + \
##!#        '%Y/%Y%m/*.nc'), \
##!#        shell = True).decode('utf-8').strip().split('\n')
##!#
##!#    # Use the actual file times to get timestamps
##!#    # -------------------------------------------
##!#    file_dates = np.array([datetime.strptime(tfile[-13:-3],'%Y%m%d%H') \
##!#        for tfile in files]) 
##!#
##!#    # Figure out where the closest files to the desired time are
##!#    # ----------------------------------------------------------
##!#    in_time = np.where( (file_dates >= dt_begin_str) & \
##!#        (file_dates <= dt_end_str))
##!#
##!#    # Select those files. Should result in all 16 channels for a single
##!#    # time
##!#    # -----------------------------------------------------------------
##!#    dates_found = file_dates[in_time]
##!#
##!#    # Set up the first data object
##!#    # ----------------------------
##!#    NCEP_data = read_NCEP(dates_found[0].strftime('%Y%m%d%H'), minlat = minlat) 
##!#    NCEP_data['dates'] = [NCEP_data['date']]
##!#
##!#    # Loop over the remaining objects
##!#    # -------------------------------
##!#    for ttime in dates_found[1:]:
##!#        local_data = read_NCEP(ttime.strftime('%Y%m%d%H'), minlat = minlat)
##!#        NCEP_data['smoke_conc_sfc'] += local_data['smoke_conc_sfc']
##!#        NCEP_data['smoke_wetsink']  += local_data['smoke_wetsink']
##!#        NCEP_data['smoke_drysink']  += local_data['smoke_drysink']
##!#
##!#        NCEP_data['dates'].append(ttime.strftime('%Y%m%d%H'))
##!#
##!#    NCEP_data['dt_begin_date'] = dt_begin_str
##!#    NCEP_data['dt_end_date']   = dt_end_str
##!#    
##!#    return NCEP_data


def plot_NCEP(NCEP_data, ax = None, labelsize = 12, \
        plot_log = True, labelticksize = 10, zoom = True, vmin = None, \
        minlat = 65., circle_bound = True, vmax = None, save = False):

    in_ax = True
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, projection = mapcrs)

    mesh = ax.pcolormesh(NCEP_data['lon'], NCEP_data['lat'], \
        NCEP_data['data'], \
        cmap = 'plasma', vmin = vmin, vmax = vmax, \
        transform = ccrs.PlateCarree(), shading = 'auto')

    cbar = plt.colorbar(mesh, ax = ax, orientation='vertical',\
        pad=0.04, fraction = 0.040, extend = 'both')
    cbar.set_label('Air Temperature', size = labelsize, weight = 'bold')
    #cbar.ax.tick_params(labelsize = labelticksize)
    ax.set_title('NCEP Air Temperature\n' + \
        NCEP_data['dt_begin_date'].strftime('%Y-%m-%d') + ' - ' + \
        NCEP_data['dt_end_date'].strftime('%Y-%m-%d'))

    ax.coastlines()
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

def plot_NCEP_figure(date_str, var, minlat = 65., vmin = None, vmax = None, \
        circle_bound = True, plot_log = True, ptitle = '', zoom = True, \
        save = False):

    plt.close('all')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection = mapcrs)

    # Read the data for this granule
    # ------------------------------
    NCEP_data = read_NCEP(date_str, minlat = minlat)

    # Plot the data for this granule
    # ------------------------------
    plot_NCEP(NCEP_data, var, ax = ax, zoom = zoom, \
            vmin = vmin, vmax = vmax, plot_log = plot_log, \
            circle_bound = circle_bound, minlat = minlat)

    ax.set_title('NCEP-RA ' + var + '\n' + NCEP_data['date'])
 
    ax.coastlines()

    if(save):
        outname = 'ncep_' + var + '_' + date_str + '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

def read_zoom_NCEP_data(date1, date2, minlat, lat_bounds, lon_bounds):

    NCEP_data1 = read_NCEP(date1, date2, minlat = minlat)

    if(lon_bounds[0] < lon_bounds[1]):
        NCEP_data1['data'] = np.ma.masked_where(~ ((NCEP_data1['lon'] >= lon_bounds[0]) & \
            (NCEP_data1['lon'] <= lon_bounds[1]) & (NCEP_data1['lat'] >= lat_bounds[0]) & \
            (NCEP_data1['lat'] <= lat_bounds[1])) , NCEP_data1['data'])
    else:
        NCEP_data1['data'] = np.ma.masked_where(~ (((NCEP_data1['lon'] >= lon_bounds[0]) | \
            (NCEP_data1['lon'] <= lon_bounds[1])) & (NCEP_data1['lat'] >= lat_bounds[0]) & \
            (NCEP_data1['lat'] <= lat_bounds[1])), NCEP_data1['data'])

    return NCEP_data1

def read_NCEP_historic(dt_begin_str1, dt_end_str1, dt_begin_str2, dt_end_str2,\
        minlat, lat_bounds, lon_bounds):

    final_end_date = datetime(year = 2021, month = 9, day = 30)
    check_end_date = dt_end_str2 + relativedelta(years = 4)
    if(check_end_date > final_end_date):
        end_year_offset = check_end_date.year - final_end_date.year
        end_idx = 5 - end_year_offset
        beg_idx = -4 - end_year_offset
    else:
        beg_idx = -4
        end_idx = 5

    combined_data = {}
    average_data = np.zeros(len(np.arange(-4, 5)) * 2)
    cnt = 0
    for ii in np.arange(beg_idx, end_idx):
        dt_begin_local1 = dt_begin_str1 + relativedelta(years = ii)
        dt_end_local1   = dt_end_str1 + relativedelta(years = ii)
        dt_begin_local2 = dt_begin_str2 + relativedelta(years = ii)
        dt_end_local2   = dt_end_str2 + relativedelta(years = ii)


        combined_data[dt_begin_local1.strftime('%Y')] = {}

        print(dt_begin_local1, dt_end_local1, dt_begin_local2, dt_end_local2)

        NCEP_data1 = read_zoom_NCEP_data(dt_begin_local1.strftime('%Y%m%d'), \
            dt_end_local1.strftime('%Y%m%d'), minlat, lat_bounds, lon_bounds)
        NCEP_data2 = read_zoom_NCEP_data(dt_begin_local2.strftime('%Y%m%d'), \
            dt_end_local2.strftime('%Y%m%d'), minlat, lat_bounds, lon_bounds)

        bdata1 = NCEP_data1['data'].flatten().compressed()
        bdata2 = NCEP_data2['data'].flatten().compressed()
  
        combined_data[dt_begin_local1.strftime('%Y')][\
            dt_begin_local1.strftime('%m/%d') + ' - ' + \
            dt_end_local1.strftime('%m/%d')] = bdata1
        combined_data[dt_begin_local2.strftime('%Y')][\
            dt_begin_local2.strftime('%m/%d') + ' - ' + \
            dt_end_local2.strftime('%m/%d')] = bdata2

        # Add the average data
        # --------------------
        average_data[cnt * 2]     = np.nanmean(bdata1)
        average_data[cnt * 2 + 1] = np.nanmean(bdata2)

        cnt += 1

    second_labels = np.array([[bkey for akey in \
        combined_data[bkey].keys()] for bkey in combined_data.keys()])
    labels = second_labels.flatten()
    second_arrays = np.array([[combined_data[bkey][akey] for akey in \
        combined_data[bkey].keys()] for bkey in combined_data.keys()])

    return second_labels, second_arrays, average_data

def plot_NCEP_event(date_str, minlat = 65., vmin = None, vmax = None, \
        ptitle = '', circle_bound = True, zoom = True, \
        lat_bounds = [-90, 90], lon_bounds = [-180, 180], \
        save = False):

    plt.close('all')
    fig = plt.figure(figsize = (12,10))
    gs  = fig.add_gridspec(nrows = 2, ncols = 3)
    ax1 = plt.subplot(gs[0,0], projection = mapcrs)
    ax2 = plt.subplot(gs[0,1], projection = mapcrs)
    ax3 = plt.subplot(gs[0,2], projection = mapcrs)
    ax4 = plt.subplot(gs[1,:])
    #ax4 = fig.add_subplot(2,2,4, projection = mapcrs)

    if(date_str[:6] == '201708'):
        cdate_begin_str1 = '20170804'
        cdate_end_str1   = '20170815'
        cdate_begin_str2 = '20170820'
        cdate_end_str2   = '20170831'
    elif(date_str[:6] == '201206'):
        cdate_begin_str1 = '20120601'
        cdate_end_str1   = '20120613'
        cdate_begin_str2 = '20120620'
        cdate_end_str2   = '20120630'
    elif(date_str[:6] == '200804'):
        cdate_begin_str1 = '20080407'
        cdate_end_str1   = '20080420'
        cdate_begin_str2 = '20080426'
        cdate_end_str2   = '20080510'

    dt_begin_str1 = datetime.strptime(cdate_begin_str1, '%Y%m%d')
    dt_end_str1   = datetime.strptime(cdate_end_str1, '%Y%m%d')
    dt_begin_str2 = datetime.strptime(cdate_begin_str2, '%Y%m%d')
    dt_end_str2   = datetime.strptime(cdate_end_str2, '%Y%m%d')

    # Read the data for this granule
    # ------------------------------
    NCEP_data1 = read_zoom_NCEP_data(cdate_begin_str1, \
        cdate_end_str1, minlat, lat_bounds, lon_bounds)
    NCEP_data2 = read_zoom_NCEP_data(cdate_begin_str2, \
        cdate_end_str2, minlat, lat_bounds, lon_bounds)

    ##!#NCEP_data1 = read_NCEP(cdate_begin_str1, cdate_end_str1, minlat = minlat)
    ##!#NCEP_data2 = read_NCEP(cdate_begin_str2, cdate_end_str2, minlat = minlat)

    ##!## Test masking the data
    ##!## ---------------------
    ##!##avg_CERES1 = np.nanmean(CERES_data1[ceres_var], axis = 0)
    ##!##avg_CERES2 = np.nanmean(CERES_data2[ceres_var], axis = 0)
    ##!##avg_CERES1_ice = np.nanmean(CERES_data1['ice_conc'], axis = 0)
    ##!##avg_CERES12 = np.nanmean(CERES_data12[ceres_var], axis = 0)
    ##!##avg_CERES22 = np.nanmean(CERES_data22[ceres_var], axis = 0)
    ##!##avg_CERES12_ice = np.nanmean(CERES_data12['ice_conc'], axis = 0)

    ##!## Test pulling out the data that only exists in both zones
    ##!## --------------------------------------------------------
    ##!##avg_CERES1 = np.ma.masked_where((avg_CERES1_ice <= min_ice) | \
    ##!#if(lon_bounds[0] < lon_bounds[1]):

    ##!#    NCEP_data1['data'] = np.ma.masked_where(~ ((NCEP_data1['lon'] >= lon_bounds[0]) & \
    ##!#        (NCEP_data1['lon'] <= lon_bounds[1]) & (NCEP_data1['lat'] >= lat_bounds[0]) & \
    ##!#        (NCEP_data1['lat'] <= lat_bounds[1])) , NCEP_data1['data'])
    ##!#    NCEP_data2['data'] = np.ma.masked_where(~ ((NCEP_data2['lon'] >= lon_bounds[0]) & \
    ##!#        (NCEP_data2['lon'] <= lon_bounds[1]) & (NCEP_data2['lat'] >= lat_bounds[0]) & \
    ##!#        (NCEP_data2['lat'] <= lat_bounds[1])) , NCEP_data2['data'])

    ##!#else:

    ##!#    NCEP_data1['data'] = np.ma.masked_where(~ (((NCEP_data1['lon'] >= lon_bounds[0]) | \
    ##!#        (NCEP_data1['lon'] <= lon_bounds[1])) & (NCEP_data1['lat'] >= lat_bounds[0]) & \
    ##!#        (NCEP_data1['lat'] <= lat_bounds[1])), NCEP_data1['data'])
    ##!#    NCEP_data2['data'] = np.ma.masked_where(~ (((NCEP_data2['lon'] >= lon_bounds[0]) | \
    ##!#        (NCEP_data2['lon'] <= lon_bounds[1])) & (NCEP_data2['lat'] >= lat_bounds[0]) & \
    ##!#        (NCEP_data2['lat'] <= lat_bounds[1])), NCEP_data2['data'])


    # Plot the data for this granule
    # ------------------------------
    plot_NCEP(NCEP_data1, ax = ax1, zoom = zoom, minlat = minlat, \
        vmin = vmin, vmax = vmax)
    plot_NCEP(NCEP_data2, ax = ax3, zoom = zoom, minlat = minlat, \
        vmin = vmin, vmax = vmax)

    # Now, process the yearly data
    # ----------------------------
    second_labels, second_arrays, average_data = \
        read_NCEP_historic(dt_begin_str1, dt_end_str1, dt_begin_str2, \
        dt_end_str2, minlat, lat_bounds, lon_bounds)

    colors = ['tab:blue','tab:orange','limegreen','tab:red','tab:purple', 'cyan', 'tab:blue', 'tab:orange', 'limegreen'] 
    for ii in range(len(second_arrays)):
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

    # Plot the line data
    # ------------------
    #ax6 = ax4.twinx()
    ax4.plot(average_data)

    ax4.set_ylabel('Surface air temperature')
    
    ax1.set_title('NCEP Air Temp\n' + \
        NCEP_data1['dt_begin_date'].strftime('%Y-%m-%d') + ' - ' + \
        NCEP_data1['dt_end_date'].strftime('%Y-%m-%d'))
    ax2.set_title('NCEP Air Temperature\n' + cdate_begin_str1 + ' - ' + \
        cdate_end_str1)
    ax3.set_title('NCEP Air Temperature\n' + cdate_begin_str2 + ' - ' + \
        cdate_end_str2)
    #ax4.set_title('CERES Average Ice Concentration')
     
    ax1.coastlines()
    ax3.coastlines()

    fig.tight_layout()

    fig2 = plt.figure(figsize = (9,9))
    axs = fig2.subplots(nrows = 3, ncols = 3)
    jj = 0
    for ii in range(second_arrays.shape[0]):
        if((ii > 2) & (ii < 6)):
            jj = 1
        elif(ii >= 6):
            jj = 2

        print(second_arrays[ii,0].shape, second_arrays[ii,1].shape)
        u_stat, p_val = mannwhitneyu(second_arrays[ii,0], second_arrays[ii,1],\
                         alternative = 'greater') 
        print(second_labels[ii,0], len(second_arrays[ii,0]), len(second_arrays[ii,1]))
        axs[jj,ii%3].hist(second_arrays[ii,0], alpha = 0.5, label = 'Before')
        axs[jj,ii%3].hist(second_arrays[ii,1], alpha = 0.5, label = 'After')

        axs[jj,ii%3].set_title(second_labels[ii,0] + '\np_val = ' + str(np.round(p_val, 2)))
        axs[jj,ii%3].set_xlabel('Air Temperature')
        axs[jj,ii%3].set_ylabel('Counts')
        axs[jj,ii%3].legend()


    fig2.tight_layout()

    if(save):
        outname = 'ncep_event_' + var + '_' + date_str + '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
        outname = 'ncep_event_boxplots_' + var + '_' + date_str + '.png'
        fig2.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

##!#def plot_NCEP_event_CERES(date_str, var, ceres_var = 'alb_clr', \
##!#        minlat = 70.5, vmin = None, vmax = None, vmin2 = None, vmax2 = None, \
##!#        min_ice = 80., min_smoke = 0, max_smoke = 2e5, plot_log = True, \
##!#        satellite = 'Aqua', ptitle = '', lat_bounds = [-90, 90], lon_bounds = [-180, 180], \
##!#        zoom = True, save = False):
##!#
##!#    plt.close('all')
##!#    fig = plt.figure(figsize = (12,10))
##!#    gs  = fig.add_gridspec(nrows = 2, ncols = 3)
##!#    ax1 = plt.subplot(gs[0,0], projection = mapcrs)
##!#    ax2 = plt.subplot(gs[0,1], projection = mapcrs)
##!#    ax3 = plt.subplot(gs[0,2], projection = mapcrs)
##!#    ax4 = plt.subplot(gs[1,:])
##!#    #ax4 = fig.add_subplot(2,2,4, projection = mapcrs)
##!#
##!#    # Read the data for this granule
##!#    # ------------------------------
##!#    NCEP_data = read_NCEP_event(date_str, minlat = minlat)
##!#
##!#    if(date_str[:6] == '201708'):
##!#        cdate_begin_str1 = '20170804'
##!#        cdate_end_str1   = '20170815'
##!#        cdate_begin_str2 = '20170820'
##!#        cdate_end_str2   = '20170831'
##!#    elif(date_str[:6] == '201206'):
##!#        cdate_begin_str1 = '20120601'
##!#        cdate_end_str1   = '20120613'
##!#        cdate_begin_str2 = '20120620'
##!#        cdate_end_str2   = '20120630'
##!#    elif(date_str[:6] == '200804'):
##!#        cdate_begin_str1 = '20080407'
##!#        cdate_end_str1   = '20080420'
##!#        cdate_begin_str2 = '20080426'
##!#        cdate_end_str2   = '20080510'
##!#
##!#    dt_begin_str1 = datetime.strptime(cdate_begin_str1, '%Y%m%d')
##!#    dt_end_str1   = datetime.strptime(cdate_end_str1, '%Y%m%d')
##!#    dt_begin_str2 = datetime.strptime(cdate_begin_str2, '%Y%m%d')
##!#    dt_end_str2   = datetime.strptime(cdate_end_str2, '%Y%m%d')
##!#
##!#    CERES_data1 = readgridCERES_daily(cdate_begin_str1,end_str = cdate_end_str1, \
##!#        satellite = satellite, minlat = minlat)
##!#    #CERES_data12 = readgridCERES_daily('20160804',end_str = '20160815', \
##!#    #    satellite = satellite, minlat = minlat)
##!#    CERES_data2 = readgridCERES_daily(cdate_begin_str2,end_str = cdate_end_str2, \
##!#        satellite = satellite, minlat = minlat)
##!#    #CERES_data22 = readgridCERES_daily('20160820',end_str = '20160831', \
##!#    #    satellite = satellite, minlat = minlat)
##!#
##!#    # Test masking the data
##!#    # ---------------------
##!#    avg_CERES1 = np.nanmean(CERES_data1[ceres_var], axis = 0)
##!#    avg_CERES2 = np.nanmean(CERES_data2[ceres_var], axis = 0)
##!#    avg_CERES1_ice = np.nanmean(CERES_data1['ice_conc'], axis = 0)
##!#    #avg_CERES12 = np.nanmean(CERES_data12[ceres_var], axis = 0)
##!#    #avg_CERES22 = np.nanmean(CERES_data22[ceres_var], axis = 0)
##!#    #avg_CERES12_ice = np.nanmean(CERES_data12['ice_conc'], axis = 0)
##!#
##!#    # Test pulling out the data that only exists in both zones
##!#    # --------------------------------------------------------
##!#    #avg_CERES1 = np.ma.masked_where((avg_CERES1_ice <= min_ice) | \
##!#    if(lon_bounds[0] < lon_bounds[1]):
##!#        CERES_data1[ceres_var] = np.ma.masked_where((avg_CERES1_ice <= min_ice) | \
##!#            (NCEP_data['smoke_conc_sfc'] < min_smoke) | \
##!#            (NCEP_data['smoke_conc_sfc'] > max_smoke), avg_CERES1)
##!#        CERES_data1[ceres_var] = np.ma.masked_where(~ ((CERES_data1['lon'] >= lon_bounds[0]) & \
##!#            (CERES_data1['lon'] <= lon_bounds[1]) & (CERES_data1['lat'] >= lat_bounds[0]) & \
##!#            (CERES_data1['lat'] <= lat_bounds[1])) , CERES_data1[ceres_var])
##!#        CERES_data1[ceres_var] = np.ma.masked_where(CERES_data1['alb_clr'] < vmin2, CERES_data1[ceres_var])
##!#
##!#        CERES_data2[ceres_var] = np.ma.masked_where((avg_CERES1_ice <= min_ice) | \
##!#            (NCEP_data['smoke_conc_sfc'] < min_smoke) | \
##!#            (NCEP_data['smoke_conc_sfc'] > max_smoke), avg_CERES2)
##!#        CERES_data2[ceres_var] = np.ma.masked_where(~ ((CERES_data2['lon'] >= lon_bounds[0]) & \
##!#            (CERES_data2['lon'] <= lon_bounds[1]) & (CERES_data2['lat'] >= lat_bounds[0]) & \
##!#            (CERES_data2['lat'] <= lat_bounds[1])) , CERES_data2[ceres_var])
##!#        CERES_data2[ceres_var] = np.ma.masked_where(CERES_data2['alb_clr'] < vmin2, CERES_data2[ceres_var])
##!#
##!#        NCEP_data['smoke_conc_sfc'] = np.ma.masked_where(\
##!#            (avg_CERES1_ice <= min_ice) | \
##!#            (NCEP_data['smoke_conc_sfc'] < min_smoke) | \
##!#            (NCEP_data['smoke_conc_sfc'] > max_smoke), \
##!#            NCEP_data['smoke_conc_sfc'])
##!#        NCEP_data['smoke_conc_sfc'] = np.ma.masked_where(~ ((CERES_data1['lon'] >= lon_bounds[0]) & \
##!#            (CERES_data1['lon'] <= lon_bounds[1]) & (CERES_data1['lat'] >= lat_bounds[0]) & \
##!#            (CERES_data1['lat'] <= lat_bounds[1])) , NCEP_data['smoke_conc_sfc'])
##!#        NCEP_data['smoke_conc_sfc'] = np.ma.masked_where(CERES_data1['alb_clr'] < vmin2, NCEP_data['smoke_conc_sfc'])
##!#
##!#    else:
##!#        CERES_data1[ceres_var] = np.ma.masked_where((avg_CERES1_ice <= min_ice) | \
##!#            (NCEP_data['smoke_conc_sfc'] < min_smoke) | \
##!#            (NCEP_data['smoke_conc_sfc'] > max_smoke), avg_CERES1)
##!#        CERES_data1[ceres_var] = np.ma.masked_where(~ (((CERES_data1['lon'] >= lon_bounds[0]) | \
##!#            (CERES_data1['lon'] <= lon_bounds[1])) & (CERES_data1['lat'] >= lat_bounds[0]) & \
##!#            (CERES_data1['lat'] <= lat_bounds[1])) , CERES_data1[ceres_var])
##!#        CERES_data1[ceres_var] = np.ma.masked_where(CERES_data1['alb_clr'] < vmin2, CERES_data1[ceres_var])
##!#
##!#        CERES_data2[ceres_var] = np.ma.masked_where((avg_CERES1_ice <= min_ice) | \
##!#            (NCEP_data['smoke_conc_sfc'] < min_smoke) | \
##!#            (NCEP_data['smoke_conc_sfc'] > max_smoke), avg_CERES2)
##!#        CERES_data2[ceres_var] = np.ma.masked_where(~ (((CERES_data2['lon'] >= lon_bounds[0]) | \
##!#            (CERES_data2['lon'] <= lon_bounds[1])) & (CERES_data2['lat'] >= lat_bounds[0]) & \
##!#            (CERES_data2['lat'] <= lat_bounds[1])) , CERES_data2[ceres_var])
##!#        CERES_data2[ceres_var] = np.ma.masked_where(CERES_data2['alb_clr'] < vmin2, CERES_data2[ceres_var])
##!#
##!#        NCEP_data['smoke_conc_sfc'] = np.ma.masked_where(\
##!#            (avg_CERES1_ice <= min_ice) | \
##!#            (NCEP_data['smoke_conc_sfc'] < min_smoke) | \
##!#            (NCEP_data['smoke_conc_sfc'] > max_smoke), \
##!#            NCEP_data['smoke_conc_sfc'])
##!#        NCEP_data['smoke_conc_sfc'] = np.ma.masked_where(~ (((CERES_data1['lon'] >= lon_bounds[0]) | \
##!#            (CERES_data1['lon'] <= lon_bounds[1])) & (CERES_data1['lat'] >= lat_bounds[0]) & \
##!#            (CERES_data1['lat'] <= lat_bounds[1])), NCEP_data['smoke_conc_sfc'])
##!#        NCEP_data['smoke_conc_sfc'] = np.ma.masked_where(CERES_data1['alb_clr'] < vmin2, NCEP_data['smoke_conc_sfc'])
##!#
##!#
##!#    # Plot the data for this granule
##!#    # ------------------------------
##!#    plot_NCEP(NCEP_data, var, ax = ax1, zoom = zoom, \
##!#            minlat = minlat, vmin = vmin, vmax = vmax, plot_log = plot_log)
##!#
##!#    # Plot the CERES data
##!#    #plotCERES_daily(cdate_begin_str1, ceres_var, end_str = cdate_end_str1, \
##!#    plotCERES_daily(CERES_data1, ceres_var, end_str = cdate_end_str1, \
##!#        satellite = satellite,  only_sea_ice = False, minlat = minlat, \
##!#        vmin = vmin2, vmax = vmax2, \
##!#        avg_data = True, ax = ax2, save = False, min_ice = min_ice, \
##!#        circle_bound = True, colorbar = True)
##!#    plotCERES_daily(CERES_data2, ceres_var, end_str = cdate_end_str2, \
##!#        satellite = satellite,  only_sea_ice = False, minlat = minlat, \
##!#        vmin = vmin2, vmax = vmax2, \
##!#        avg_data = True, ax = ax3, save = False, min_ice = min_ice, \
##!#        circle_bound = True, colorbar = True)
##!#    #plotCERES_daily(cdate_begin_str2, 'ice_conc', end_str = cdate_end_str2, \
##!#    #    satellite = satellite,  only_sea_ice = False, minlat = minlat, \
##!#    #    avg_data = True, ax = ax4, save = False, \
##!#    #    circle_bound = True, colorbar = True)
##!#
##!#    # Now, process the yearly data
##!#    # ----------------------------
##!#    combined_data = {}
##!#    for ii in np.arange(-4, 5):
##!#        dt_begin_local1 = dt_begin_str1 + relativedelta(years = ii)
##!#        dt_end_local1   = dt_end_str1 + relativedelta(years = ii)
##!#        dt_begin_local2 = dt_begin_str2 + relativedelta(years = ii)
##!#        dt_end_local2   = dt_end_str2 + relativedelta(years = ii)
##!#
##!#
##!#        combined_data[dt_begin_local1.strftime('%Y')] = {}
##!#
##!#        print(dt_begin_local1, dt_end_local1, dt_begin_local2, dt_end_local2)
##!#        CERES_data1 = readgridCERES_daily(dt_begin_local1.strftime('%Y%m%d'),\
##!#            end_str = dt_end_local1.strftime('%Y%m%d'), \
##!#            satellite = satellite, minlat = minlat)
##!#        CERES_data2 = readgridCERES_daily(dt_begin_local1.strftime('%Y%m%d'),\
##!#            end_str = dt_end_local2.strftime('%Y%m%d'), \
##!#            satellite = satellite, minlat = minlat)
##!#
##!#        avg_CERES1 = np.nanmean(CERES_data1[ceres_var], axis = 0)
##!#        avg_CERES2 = np.nanmean(CERES_data2[ceres_var], axis = 0)
##!#        avg_CERES1_ice = np.nanmean(CERES_data1['ice_conc'], axis = 0)
##!#
##!#        if(lon_bounds[0] < lon_bounds[1]):
##!#            CERES_data1[ceres_var] = np.ma.masked_where((avg_CERES1_ice <= min_ice) | \
##!#                (NCEP_data['smoke_conc_sfc'] < min_smoke) | \
##!#                (NCEP_data['smoke_conc_sfc'] > max_smoke), avg_CERES1)  
##!#            CERES_data1[ceres_var] = np.ma.masked_where(~ ((CERES_data1['lon'] >= lon_bounds[0]) & \
##!#                (CERES_data1['lon'] <= lon_bounds[1]) & (CERES_data1['lat'] >= lat_bounds[0]) & \
##!#                (CERES_data1['lat'] <= lat_bounds[1])) , CERES_data1[ceres_var])
##!#            CERES_data1[ceres_var] = np.ma.masked_where(CERES_data1['alb_clr'] < vmin2, CERES_data1[ceres_var])
##!#            CERES_data2[ceres_var] = np.ma.masked_where((avg_CERES1_ice <= min_ice) | \
##!#                (NCEP_data['smoke_conc_sfc'] < min_smoke) | \
##!#                (NCEP_data['smoke_conc_sfc'] > max_smoke), avg_CERES2)  
##!#            CERES_data2[ceres_var] = np.ma.masked_where(~ ((CERES_data2['lon'] >= lon_bounds[0]) & \
##!#                (CERES_data2['lon'] <= lon_bounds[1]) & (CERES_data2['lat'] >= lat_bounds[0]) & \
##!#                (CERES_data2['lat'] <= lat_bounds[1])) , CERES_data2[ceres_var])
##!#            CERES_data2[ceres_var] = np.ma.masked_where(CERES_data2['alb_clr'] < vmin2, CERES_data2[ceres_var])
##!#        else:
##!#            CERES_data1[ceres_var] = np.ma.masked_where((avg_CERES1_ice <= min_ice) | \
##!#                (NCEP_data['smoke_conc_sfc'] < min_smoke) | \
##!#                (NCEP_data['smoke_conc_sfc'] > max_smoke), avg_CERES1)
##!#            CERES_data1[ceres_var] = np.ma.masked_where(~ (((CERES_data1['lon'] >= lon_bounds[0]) | \
##!#                (CERES_data1['lon'] <= lon_bounds[1])) & (CERES_data1['lat'] >= lat_bounds[0]) & \
##!#                (CERES_data1['lat'] <= lat_bounds[1])) , CERES_data1[ceres_var])
##!#            CERES_data1[ceres_var] = np.ma.masked_where(CERES_data1['alb_clr'] < vmin2, CERES_data1[ceres_var])
##!#            CERES_data2[ceres_var] = np.ma.masked_where((avg_CERES1_ice <= min_ice) | \
##!#                (NCEP_data['smoke_conc_sfc'] < min_smoke) | \
##!#                (NCEP_data['smoke_conc_sfc'] > max_smoke), avg_CERES2)
##!#            CERES_data2[ceres_var] = np.ma.masked_where(~ (((CERES_data1['lon'] >= lon_bounds[0]) | \
##!#                (CERES_data1['lon'] <= lon_bounds[1])) & (CERES_data1['lat'] >= lat_bounds[0]) & \
##!#                (CERES_data1['lat'] <= lat_bounds[1])) , CERES_data2[ceres_var])
##!#            CERES_data2[ceres_var] = np.ma.masked_where(CERES_data2['alb_clr'] < vmin2, CERES_data2[ceres_var])
##!#
##!#        bdata1 = CERES_data1[ceres_var].flatten().compressed()
##!#        bdata2 = CERES_data2[ceres_var].flatten().compressed()
##!#  
##!#        combined_data[dt_begin_local1.strftime('%Y')][\
##!#            dt_begin_local1.strftime('%m/%d') + ' - ' + \
##!#            dt_end_local1.strftime('%m/%d')] = bdata1
##!#        combined_data[dt_begin_local2.strftime('%Y')][\
##!#            dt_begin_local2.strftime('%m/%d') + ' - ' + \
##!#            dt_end_local2.strftime('%m/%d')] = bdata2
##!#
##!#    second_labels = np.array([[bkey for akey in \
##!#        combined_data[bkey].keys()] for bkey in combined_data.keys()])
##!#    labels = second_labels.flatten()
##!#    second_arrays = np.array([[combined_data[bkey][akey] for akey in \
##!#        combined_data[bkey].keys()] for bkey in combined_data.keys()])
##!#    final_arrays = second_arrays.flatten()
##!#
##!#    #print(len(avg_CERES1.flatten().compressed()), len(avg_CERES2.flatten().compressed()))
##!#    #bdata1 = CERES_data1[ceres_var].flatten().compressed()
##!#    #bdata12 = CERES_data12[ceres_var].flatten().compressed()
##!#    #bdata2 = CERES_data2[ceres_var].flatten().compressed()
##!#    #bdata22 = CERES_data22[ceres_var].flatten().compressed()
##!#    #print(len(bdata1), len(bdata2))
##!#    #ax4.boxplot([bdata12, bdata22 , bdata1, bdata2])
##!#  
##!#    colors = ['tab:blue','tab:orange','limegreen','tab:red','tab:purple', 'cyan', 'tab:blue', 'tab:orange', 'limegreen'] 
##!#    for ii in range(len(second_arrays)):
##!#        boxs = ax4.boxplot([second_arrays[ii,0], second_arrays[ii,1]], \
##!#            labels = [second_labels[ii,0], second_labels[ii,1]], \
##!#            positions = [ii*2, ii*2 + 1], widths = [0.75, 0.75], \
##!#            patch_artist = True)
##!#        for item in ['boxes','whiskers','fliers','caps']:
##!#            plt.setp(boxs[item], color = colors[ii])
##!#        plt.setp(boxs['medians'], color = 'black')
##!#        for patch in boxs['boxes']:
##!#        #    patch.set_edgecolor(colors[ii])
##!#            patch.set_facecolor('white')
##!#
##!#    ax4.set_ylabel('Clear-sky Albedo')
##!#    
##!#    ax1.set_title('NCEP-RA ' + var + '\n' + \
##!#        NCEP_data['dt_begin_date'].strftime('%Y-%m-%d') + ' - ' + \
##!#        NCEP_data['dt_end_date'].strftime('%Y-%m-%d'))
##!#    ax2.set_title('CERES Average Clear-sky Albedo\n' + cdate_begin_str1 + ' - ' + \
##!#        cdate_end_str1)
##!#    ax3.set_title('CERES Average Clear-sky Albedo\n' + cdate_begin_str2 + ' - ' + \
##!#        cdate_end_str2)
##!#    #ax4.set_title('CERES Average Ice Concentration')
##!#     
##!#    ax1.coastlines()
##!#
##!#    fig.tight_layout()
##!#
##!#    fig2 = plt.figure(figsize = (9,9))
##!#    axs = fig2.subplots(nrows = 3, ncols = 3)
##!#    jj = 0
##!#    for ii in range(second_arrays.shape[0]):
##!#        if((ii > 2) & (ii < 6)):
##!#            jj = 1
##!#        elif(ii >= 6):
##!#            jj = 2
##!#
##!#        u_stat, p_val = mannwhitneyu(second_arrays[ii,0], second_arrays[ii,1],\
##!#                         alternative = 'greater') 
##!#        print(second_labels[ii,0], len(second_arrays[ii,0]), len(second_arrays[ii,1]))
##!#        axs[jj,ii%3].hist(second_arrays[ii,0], alpha = 0.5, label = 'Before')
##!#        axs[jj,ii%3].hist(second_arrays[ii,1], alpha = 0.5, label = 'After')
##!#        #if(len(second_arrays[ii,0]) > len(second_arrays[ii,1])):
##!#        #    axs[jj,ii%3].hist(second_arrays[ii,1], label = 'After')
##!#        #    axs[jj,ii%3].hist(second_arrays[ii,0], label = 'Before')
##!#        #else:
##!#        #    axs[jj,ii%3].hist(second_arrays[ii,0], label = 'Before')
##!#        #    axs[jj,ii%3].hist(second_arrays[ii,1], label = 'After')
##!#        axs[jj,ii%3].set_title(second_labels[ii,0] + '\np_val = ' + str(np.round(p_val, 2)))
##!#        axs[jj,ii%3].set_xlabel('Clear-sky albedo')
##!#        axs[jj,ii%3].set_ylabel('Counts')
##!#        axs[jj,ii%3].legend()
##!#
##!#
##!#    fig2.tight_layout()
##!#
##!#    if(save):
##!#        outname = 'ncep_event_' + var + '_' + date_str + '.png'
##!#        fig.savefig(outname, dpi = 300)
##!#        print("Saved image", outname)
##!#        outname = 'ncep_event_boxplots_' + var + '_' + date_str + '.png'
##!#        fig2.savefig(outname, dpi = 300)
##!#        print("Saved image", outname)
##!#    else:
##!#        plt.show()
##!#

def plot_NCEP_region_time_series_combined(date_str, begin_date, end_date, \
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


