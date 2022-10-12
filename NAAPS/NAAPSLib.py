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
from scipy.stats import pearsonr,spearmanr
import subprocess
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as cm
from matplotlib.dates import DateFormatter
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from satpy import find_files_and_readers
from satpy.scene import Scene
from satpy.writers import get_enhanced_image
from glob import glob

sys.path.append('/home/bsorenson/Research/CERES')
from gridCERESLib import *
sys.path.append('/home/bsorenson/')
#from python_lib import plot_trend_line, plot_subplot_label, plot_figure_text, \
#    nearest_gridpoint, aerosol_event_dict, init_proj, \
#    convert_radiance_to_temp
from python_lib import *

data_dir = '/home/bsorenson/data/NAAPS/'

datacrs = ccrs.PlateCarree()
mapcrs = ccrs.NorthPolarStereo()

event_dict = {
    '20080422': {
        'start': '2008042100',
        'end':   '2008042500'
    },
    '20170816': {
        'start': '2017081600',
        'end':   '2017082200'
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


def read_NAAPS(date_str, minlat = 65.):

    print(date_str) 
    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H')

    filename = dt_date_str.strftime(data_dir + '%Y/%Y%m/NVA_CLIMO1misr_sfc_conc_sink_%Y%m%d%H.nc')

    data = Dataset(filename, 'r')

    idx1 = np.where(data['lat'][:] >= minlat)[0][0]
    idx2 = np.where(data['lat'][:] >= minlat)[0][-1] + 1

    xx, yy = np.meshgrid(data['lon'][:], data['lat'][idx1:idx2])

    NAAPS_data = {}
    NAAPS_data['filename']       = filename
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
 
    return NAAPS_data

def read_NAAPS_event(date_str, minlat = 65.):
 
    dt_date_str = datetime.strptime(date_str, '%Y%m%d%H')

    # Get the event starting and ending times
    # ---------------------------------------
    begin_str = event_dict[date_str]['start']
    end_str   = event_dict[date_str]['end']

    dt_begin_str = datetime.strptime(begin_str, '%Y%m%d%H')
    dt_end_str   = datetime.strptime(end_str, '%Y%m%d%H')

    # Grab the files from the data directory
    # --------------------------------------
    files = subprocess.check_output(dt_date_str.strftime('ls ' + data_dir + \
        '%Y/%Y%m/*.nc'), \
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
    NAAPS_data = read_NAAPS(dates_found[0].strftime('%Y%m%d%H'), minlat = minlat) 
    NAAPS_data['dates'] = [NAAPS_data['date']]

    # Loop over the remaining objects
    # -------------------------------
    for ttime in dates_found[1:]:
        local_data = read_NAAPS(ttime.strftime('%Y%m%d%H'), minlat = minlat)
        NAAPS_data['smoke_conc_sfc'] += local_data['smoke_conc_sfc']
        NAAPS_data['smoke_wetsink']  += local_data['smoke_wetsink']
        NAAPS_data['smoke_drysink']  += local_data['smoke_drysink']

        NAAPS_data['dates'].append(ttime.strftime('%Y%m%d%H'))

    NAAPS_data['dt_begin_date'] = dt_begin_str
    NAAPS_data['dt_end_date']   = dt_end_str
    
    return NAAPS_data


def plot_NAAPS(NAAPS_data, var, ax = None, labelsize = 12, \
        plot_log = True, labelticksize = 10, zoom = True, vmin = None, \
        minlat = 65., circle_bound = True, vmax = None, save = False):

    in_ax = True
    if(ax is None):
        in_ax = False
        plt.close('all')
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, projection = mapcrs)

    min_val = 0.0
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
        pad=0.03, extend = 'both')

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
            circle_bound = circle_bound)

    ax.set_title('NAAPS-RA ' + var + '\n' + NAAPS_data['date'])
 
    ax.coastlines()

    if(save):
        outname = 'naaps_' + var + '_' + date_str + '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()

def plot_NAAPS_event(date_str, var, minlat = 65., vmin = None, vmax = None, \
        plot_log = True, ptitle = '', circle_bound = True, zoom = True, \
        save = False):

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
            circle_bound = circle_bound)

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

def plot_NAAPS_event_CERES(date_str, var, ceres_var = 'alb_clr', \
        minlat = 70.5, vmin = None, vmax = None, min_ice = 80., \
        min_smoke = 0, plot_log = True, \
        satellite = 'Aqua', ptitle = '', zoom = True, save = False):

    plt.close('all')
    fig = plt.figure(figsize = (12,10))
    ax1 = fig.add_subplot(2,2,1, projection = mapcrs)
    ax2 = fig.add_subplot(2,2,2, projection = mapcrs)
    ax3 = fig.add_subplot(2,2,3, projection = mapcrs)
    ax4 = fig.add_subplot(2,2,4)
    #ax4 = fig.add_subplot(2,2,4, projection = mapcrs)

    # Read the data for this granule
    # ------------------------------
    NAAPS_data = read_NAAPS_event(date_str, minlat = minlat)

    if(date_str[:6] == '201708'):
        cdate_begin_str1 = '20170801'
        cdate_end_str1   = '20170815'
        cdate_begin_str2 = '20170822'
        cdate_end_str2   = '20170831'

    CERES_data1 = readgridCERES_daily('20170801',end_str = '20170815', \
        satellite = satellite, minlat = minlat)
    CERES_data2 = readgridCERES_daily('20170820',end_str = '20170831', \
        satellite = satellite, minlat = minlat)

    # Test masking the data
    # ---------------------
    avg_CERES1 = np.nanmean(CERES_data1[ceres_var], axis = 0)
    avg_CERES2 = np.nanmean(CERES_data2[ceres_var], axis = 0)
    avg_CERES1_ice = np.nanmean(CERES_data1['ice_conc'], axis = 0)

    # Test pulling out the data that only exists in both zones
    # --------------------------------------------------------
    avg_CERES1 = np.ma.masked_where((avg_CERES1_ice <= min_ice) | \
        (NAAPS_data['smoke_conc_sfc'] < min_smoke), avg_CERES1)  
    avg_CERES2 = np.ma.masked_where((avg_CERES1_ice <= min_ice) | \
        (NAAPS_data['smoke_conc_sfc'] < min_smoke), avg_CERES2)  

    NAAPS_data['smoke_conc_sfc'] = np.ma.masked_where(\
        (avg_CERES1_ice <= min_ice) | \
        (NAAPS_data['smoke_conc_sfc'] < min_smoke), \
        NAAPS_data['smoke_conc_sfc'])


    # Plot the data for this granule
    # ------------------------------
    plot_NAAPS(NAAPS_data, var, ax = ax1, zoom = zoom, \
            minlat = minlat, vmin = vmin, vmax = vmax, plot_log = plot_log)

    # Plot the CERES data
    plotCERES_daily(cdate_begin_str1, ceres_var, end_str = cdate_end_str1, \
        satellite = satellite,  only_sea_ice = False, minlat = minlat, \
        avg_data = True, ax = ax2, save = False, min_ice = min_ice, \
        circle_bound = True, colorbar = True)
    plotCERES_daily(cdate_begin_str2, ceres_var, end_str = cdate_end_str2, \
        satellite = satellite,  only_sea_ice = False, minlat = minlat, \
        avg_data = True, ax = ax3, save = False, min_ice = min_ice, \
        circle_bound = True, colorbar = True)
    #plotCERES_daily(cdate_begin_str2, 'ice_conc', end_str = cdate_end_str2, \
    #    satellite = satellite,  only_sea_ice = False, minlat = minlat, \
    #    avg_data = True, ax = ax4, save = False, \
    #    circle_bound = True, colorbar = True)

    ax4.boxplot([avg_CERES1.flatten().compressed(), avg_CERES2.flatten().compressed()])

    ax1.set_title('NAAPS-RA ' + var + '\n' + \
        NAAPS_data['dt_begin_date'].strftime('%Y-%m-%d') + ' - ' + \
        NAAPS_data['dt_end_date'].strftime('%Y-%m-%d'))
    ax2.set_title('CERES Average Clear-sky Albedo\n' + cdate_begin_str1 + ' - ' + \
        cdate_end_str1)
    ax3.set_title('CERES Average Clear-sky Albedo\n' + cdate_begin_str2 + ' - ' + \
        cdate_end_str2)
    #ax4.set_title('CERES Average Ice Concentration')
     
    ax1.coastlines()

    fig.tight_layout()

    if(save):
        outname = 'naaps_ceres_event_' + var + '_' + date_str + '.png'
        fig.savefig(outname, dpi = 300)
        print("Saved image", outname)
    else:
        plt.show()


