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

home_dir = os.environ['HOME']
sys.path.append(home_dir)
from python_lib import *

data_dir = home_dir + '/data/NCEP/'

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

    print(dt_begin_str, dt_end_str) 

    filename = dt_begin_str.strftime(data_dir + 'air.sig995.%Y.nc')

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
